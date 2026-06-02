#!/usr/bin/env python3

"""
EdgeRDEFromNormalized — limma-trend DE on a pre-normalized expression matrix
SalmonStreamer Pipeline Module

Use this when you already have a TMM-normalized expression matrix (or any
matrix on a sample-comparable scale) and you cannot re-run the standard
edgeR quasi-likelihood pipeline because it requires raw integer counts.

Typical use case
----------------
1. Run EdgeRDE on raw counts with --export-normalized-expression to dump
   TMM_normalized_CPM.tsv / TMM_normalized_logCPM.tsv.
2. Post-process the normalized matrix outside the pipeline (e.g. collapse
   paralogs by summing TMM-normalized CPMs).
3. Feed the corrected matrix back into the pipeline via this subcommand to
   produce DE tables and the usual plots.

This wraps src/R_src/edger_de_from_normalized.R, which applies the limma-trend
workflow recommended by the limma user guide for already-normalized data:

    logCPM -> lmFit(design) -> contrasts.fit(contrast) -> eBayes(trend=TRUE)

Metadata format is identical to EdgeRDE: see that module's docstring for
details on --metadata-file vs --group-samples.

@Author: Luis Javier Madrigal-Roca
@Date:   2026-06-02
"""

import argparse
import os
import sys
import subprocess
import tempfile
import csv
import re


def validate_expression_file(filepath):
    """Validate the normalized expression matrix and return sample column names."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except IOError as e:
        raise ValueError(f"Cannot read expression file: {e}")

    if len(lines) < 2:
        raise ValueError(
            "Expression file must have at least 2 rows (header + 1 gene). "
            f"Found {len(lines)} rows."
        )

    header = lines[0].rstrip("\n").split('\t')
    if len(header) < 2:
        raise ValueError(
            "Expression file must have at least 2 columns "
            f"(gene_id + 1 sample). Found {len(header)} columns."
        )

    # Spot-check the first few data rows for numeric values
    for i, line in enumerate(lines[1:5], start=2):
        fields = line.rstrip("\n").split('\t')
        if len(fields) != len(header):
            raise ValueError(
                f"Row {i} has {len(fields)} columns but header has {len(header)}"
            )
        for v in fields[1:]:
            if v == "" or v == "NA":
                continue
            try:
                float(v)
            except ValueError:
                raise ValueError(
                    f"Row {i}, column {fields[0]}: expected numeric value, "
                    f"got '{v}'. (This subcommand expects an already-normalized "
                    "expression matrix, not raw counts.)"
                )

    return header[1:]


def detect_input_scale(filepath):
    """
    Auto-detect whether the matrix is on the log2-CPM scale.

    Heuristic: log-CPM matrices typically contain negative values
    (often around -4 to -8 for very lowly expressed rows after prior_count
    addition) while raw CPM matrices are non-negative. We sample the first
    ~5000 numeric values and return 'logcpm' if any are negative, else 'cpm'.
    """
    sampled = 0
    target = 5000
    with open(filepath, 'r') as f:
        next(f, None)  # header
        for line in f:
            for v in line.rstrip("\n").split('\t')[1:]:
                if v in ("", "NA"):
                    continue
                try:
                    x = float(v)
                except ValueError:
                    continue
                if x < 0:
                    return "logcpm"
                sampled += 1
                if sampled >= target:
                    return "cpm"
    return "cpm"  # default to cpm if file too small / all non-negative


def validate_metadata_file(filepath, sample_names, sample_suffix=None):
    """Same validation as EdgeRDE — check structure and sample overlap."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except IOError as e:
        raise ValueError(f"Cannot read metadata file: {e}")

    if len(lines) < 2:
        raise ValueError("Metadata file must have at least 2 rows (header + 1 sample)")

    header = lines[0].rstrip("\n").split('\t')
    header_lower = [h.lower() for h in header]

    if 'sample_name' not in header_lower:
        raise ValueError(
            "Metadata file must have 'sample_name' column. "
            f"Found columns: {', '.join(header)}"
        )

    has_group = 'group' in header_lower
    has_species_tissue = ('species' in header_lower and 'tissue' in header_lower)
    if not (has_group or has_species_tissue):
        raise ValueError(
            "Metadata file must have either 'group' OR both 'species' and "
            f"'tissue' columns. Found columns: {', '.join(header)}"
        )

    sample_col = header_lower.index('sample_name')
    metadata_samples = set()
    for i, line in enumerate(lines[1:], start=2):
        fields = line.rstrip("\n").split('\t')
        if len(fields) <= sample_col:
            raise ValueError(f"Metadata row {i} has too few columns")
        metadata_samples.add(fields[sample_col])

    stripped = set(sample_names)
    if sample_suffix:
        stripped = set(re.sub(sample_suffix, "", s) for s in sample_names)

    matching = metadata_samples & stripped
    if len(matching) == 0:
        raise ValueError(
            "No samples match between expression file and metadata. "
            f"Expression samples: {', '.join(sorted(list(stripped)[:5]))}... "
            f"Metadata samples: {', '.join(sorted(list(metadata_samples)[:5]))}..."
        )
    pct_matched = 100 * len(matching) / len(stripped)
    if pct_matched < 90:
        print(
            f"WARNING: Only {pct_matched:.1f}% of expression samples match metadata.",
            file=sys.stderr
        )


def build_metadata_from_groups(group_specs):
    rows = []
    for spec in group_specs:
        if ":" not in spec:
            raise ValueError(
                f"--group-samples entry '{spec}' must be in the form "
                "'GroupLabel:sample1,sample2,...'"
            )
        label, samples_str = spec.split(":", 1)
        for s in samples_str.split(","):
            s = s.strip()
            if s:
                rows.append((s, label.strip()))
    return rows


def write_temp_metadata(rows, tmp_dir):
    fh = tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=tmp_dir, delete=False
    )
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["sample_name", "group"])
    writer.writerows(rows)
    fh.close()
    return fh.name


def main(args):
    r_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "src", "R_src", "edger_de_from_normalized.R"
    )

    if not os.path.isfile(r_script):
        print(f"ERROR: R script not found at {r_script}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.rscript_executable):
        print(
            f"ERROR: Rscript executable not found at: {args.rscript_executable}\n"
            f"       Specify a different path with --rscript-executable",
            file=sys.stderr
        )
        sys.exit(1)

    if args.fdr_threshold < 0 or args.fdr_threshold > 1:
        print(
            f"ERROR: FDR threshold must be between 0 and 1 (got {args.fdr_threshold})",
            file=sys.stderr
        )
        sys.exit(1)
    if args.logfc_threshold < 0:
        print(
            f"ERROR: logFC threshold must be non-negative (got {args.logfc_threshold})",
            file=sys.stderr
        )
        sys.exit(1)
    if args.prior_count <= 0:
        print(
            f"ERROR: --prior-count must be > 0 (got {args.prior_count})",
            file=sys.stderr
        )
        sys.exit(1)

    print("Validating expression file...", file=sys.stderr)
    try:
        sample_names = validate_expression_file(args.expression_file)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Resolve input scale
    if args.input_scale == "auto":
        detected = detect_input_scale(args.expression_file)
        print(f"Detected input scale: {detected} (override with --input-scale)",
              file=sys.stderr)
        input_scale = detected
    else:
        input_scale = args.input_scale

    # Resolve metadata
    temp_meta_path = None
    if args.metadata_file:
        if not os.path.isfile(args.metadata_file):
            print(f"ERROR: Metadata file not found: {args.metadata_file}",
                  file=sys.stderr)
            sys.exit(1)
        metadata_path = args.metadata_file
    elif args.group_samples:
        try:
            rows = build_metadata_from_groups(args.group_samples)
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        os.makedirs(args.output_dir, exist_ok=True)
        temp_meta_path = write_temp_metadata(rows, args.output_dir)
        metadata_path = temp_meta_path
        print(f"Built temporary metadata file from --group-samples: {temp_meta_path}")
    else:
        print(
            "ERROR: Provide either --metadata-file or --group-samples.",
            file=sys.stderr
        )
        sys.exit(1)

    if args.metadata_file:
        print("Validating metadata file...", file=sys.stderr)
        try:
            validate_metadata_file(
                metadata_path, sample_names,
                sample_suffix=args.sample_suffix
            )
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    sample_suffix = args.sample_suffix if args.sample_suffix else "NULL"

    cmd = [
        args.rscript_executable,
        r_script,
        args.expression_file,
        args.output_dir,
        metadata_path,
        str(args.fdr_threshold),
        str(args.logfc_threshold),
        sample_suffix,
        input_scale,
        str(args.prior_count),
    ]

    print("=" * 70)
    print("SalmonStreamer EdgeRDEFromNormalized")
    print("=" * 70)
    print(f"Expression file : {args.expression_file}")
    print(f"Metadata file   : {metadata_path}")
    print(f"Output dir      : {args.output_dir}")
    print(f"Input scale     : {input_scale}")
    if input_scale == "cpm":
        print(f"Prior count     : {args.prior_count}")
    print(f"FDR threshold   : {args.fdr_threshold}")
    print(f"logFC threshold : {args.logfc_threshold}")
    if args.sample_suffix:
        print(f"Sample suffix   : {args.sample_suffix}")
    print()

    try:
        subprocess.run(cmd, check=True, capture_output=False)
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: R script failed with exit code {e.returncode}",
            file=sys.stderr
        )
        sys.exit(e.returncode)
    except FileNotFoundError as e:
        print(
            f"ERROR: Cannot find Rscript executable: {args.rscript_executable}\n"
            f"       {e}",
            file=sys.stderr
        )
        sys.exit(1)
    finally:
        if temp_meta_path and os.path.isfile(temp_meta_path):
            os.remove(temp_meta_path)

    print("\nEdgeRDEFromNormalized complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Run limma-trend differential expression on a pre-normalized "
            "expression matrix (SalmonStreamer EdgeRDEFromNormalized subcommand)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Typical workflow:
  1. SalmonStreamer.py EdgeRDE ... --export-normalized-expression
  2. (collapse paralogs / curate the TMM_normalized matrix outside the pipeline)
  3. SalmonStreamer.py EdgeRDEFromNormalized \\
       --expression-file TMM_normalized_logCPM.collapsed.tsv \\
       --metadata-file metadata.tsv \\
       --output-dir DE_results_collapsed/

Input scale:
  --input-scale logcpm  -> values used as-is (default if file contains
                            negative values, e.g. our TMM_normalized_logCPM.tsv)
  --input-scale cpm     -> log2(x + --prior-count) applied internally
                            (e.g. our TMM_normalized_CPM.tsv with prior 0.25)
  --input-scale auto    -> detect from the matrix (negative values -> logcpm)
        """
    )

    parser.add_argument(
        "--expression-file", required=True,
        help="Tab-separated normalized expression matrix (genes x samples, "
             "first column = TranscriptID)."
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory for all output files."
    )
    parser.add_argument(
        "--metadata-file", default=None,
        help=(
            "Tab-separated sample metadata. Required columns: sample_name, group "
            "(or species + tissue). Takes priority over --group-samples."
        )
    )
    parser.add_argument(
        "--group-samples", nargs="+", default=None, metavar="GROUP:s1,s2,...",
        help=(
            "Inline group specification when no metadata file is available. "
            "Format: 'GroupLabel:sample1,sample2,...' (one entry per group)."
        )
    )
    parser.add_argument(
        "--input-scale", choices=["auto", "logcpm", "cpm"], default="auto",
        help=(
            "Scale of the input matrix. 'logcpm' = use as-is; "
            "'cpm' = apply log2(x + prior_count); "
            "'auto' = detect from the data (default)."
        )
    )
    parser.add_argument(
        "--prior-count", type=float, default=0.25,
        help="Pseudo-count added before log2 when --input-scale is 'cpm' (default: 0.25)."
    )
    parser.add_argument(
        "--fdr-threshold", type=float, default=0.05,
        help="FDR significance cutoff (default: 0.05)."
    )
    parser.add_argument(
        "--logfc-threshold", type=float, default=1.0,
        help="|log2FC| threshold for significance categories (default: 1.0)."
    )
    parser.add_argument(
        "--sample-suffix", default=None,
        help=(
            "Regex pattern stripped from expression matrix column names before "
            "matching to metadata sample_name values (e.g. '_R1_filtered$'). Optional."
        )
    )
    parser.add_argument(
        "--rscript-executable",
        default=os.path.expanduser("~/.conda/envs/PyR/bin/Rscript"),
        help="Path to the Rscript executable (default: ~/.conda/envs/PyR/bin/Rscript)."
    )

    main(parser.parse_args())
