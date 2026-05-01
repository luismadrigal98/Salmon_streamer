#!/usr/bin/env python3

"""
EdgeRDE — edgeR-based Differential Expression Analysis
SalmonStreamer Pipeline Module

Wraps src/R_src/edger_de.R, which implements a quasi-likelihood edgeR pipeline
with TMM normalisation, adaptive filterByExpr filtering, all pairwise contrasts,
PCA, sample-correlation heatmap, volcano plots, and a summary report.

Metadata file format (tab-separated, --metadata-file):
  REQUIRED COLUMNS:
    sample_name    Must match count matrix column headers exactly
                   (after any --sample-suffix stripping).
    group          Experimental group label used for contrasts
                   (e.g. "SpeciesA_Leaf", "SpeciesB_Bud").
                   ALTERNATIVE: provide 'species' + 'tissue' columns instead;
                   group will be derived as <species>_<tissue>.

  OPTIONAL COLUMNS (ignored by edger_de.R, may be present):
    species, tissue, platform, library, flowcell, lane, replicate, ...

CLI groups (--group-samples, used when --metadata-file is absent):
  Format:  GroupLabel:sample1,sample2,sample3  (space-separated per group)
  Example: --group-samples SpeciesA:s1,s2,s3 SpeciesB:s4,s5,s6
  The Python wrapper converts this into a temporary metadata file and passes
  it to the R script; --metadata-file takes priority if both are supplied.

@Author: Luis Javier Madrigal-Roca
@Date:   2026-04-30
"""

import argparse
import os
import sys
import subprocess
import tempfile
import csv
import re


def validate_expression_file(filepath):
    """
    Validate expression file format and contents.
    
    Parameters
    ----------
    filepath : str
        Path to expression file.
    
    Raises
    ------
    ValueError
        If file format is invalid or contains errors.
    """
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
    
    # Parse header
    try:
        header = lines[0].strip().split('\t')
    except Exception as e:
        raise ValueError(f"Cannot parse expression file header: {e}")
    
    if len(header) < 2:
        raise ValueError(
            f"Expression file must have at least 2 columns (gene_id + 1 sample). "
            f"Found {len(header)} columns."
        )
    
    # Check that all value columns are numeric
    for i, line in enumerate(lines[1:5], start=2):  # Check first few data rows
        try:
            fields = line.strip().split('\t')
            if len(fields) != len(header):
                raise ValueError(
                    f"Row {i} has {len(fields)} columns but header has {len(header)}"
                )
            # Try to convert all count columns to numbers (accept int or float)
            for count_str in fields[1:]:
                try:
                    float(count_str)
                except ValueError:
                    raise ValueError(
                        f"Row {i}, column {fields[0]}: "
                        f"expected numeric count, got '{count_str}'"
                    )
        except Exception as e:
            raise ValueError(f"Error validating expression file row {i}: {e}")
    
    return header[1:]  # Return sample names


def validate_metadata_file(filepath, sample_names, sample_suffix=None):
    """
    Validate metadata file and check compatibility with expression file.
    
    Parameters
    ----------
    filepath : str
        Path to metadata file.
    sample_names : list[str]
        Sample names from expression file.
    sample_suffix : str, optional
        Regex pattern to strip from sample names.
    
    Raises
    ------
    ValueError
        If metadata file is invalid or incompatible with expression file.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except IOError as e:
        raise ValueError(f"Cannot read metadata file: {e}")
    
    if len(lines) < 2:
        raise ValueError("Metadata file must have at least 2 rows (header + 1 sample)")
    
    # Parse header (case-insensitive)
    header = lines[0].strip().split('\t')
    header_lower = [h.lower() for h in header]
    
    if 'sample_name' not in header_lower:
        raise ValueError(
            "Metadata file must have 'sample_name' column. "
            f"Found columns: {', '.join(header)}"
        )
    
    # Check for required group columns
    has_group = 'group' in header_lower
    has_species_tissue = ('species' in header_lower and 'tissue' in header_lower)
    
    if not (has_group or has_species_tissue):
        raise ValueError(
            "Metadata file must have either 'group' OR both 'species' and 'tissue' columns. "
            f"Found columns: {', '.join(header)}"
        )
    
    # Get sample_name column index
    sample_col = header_lower.index('sample_name')
    
    # Parse metadata samples and match to expression file
    metadata_samples = set()
    for i, line in enumerate(lines[1:], start=2):
        fields = line.strip().split('\t')
        if len(fields) <= sample_col:
            raise ValueError(f"Metadata row {i} has too few columns")
        metadata_samples.add(fields[sample_col])
    
    # Check sample matching (accounting for suffix stripping)
    stripped_sample_names = set(sample_names)
    if sample_suffix:
        stripped_sample_names = set(re.sub(sample_suffix, "", s) for s in sample_names)
    
    matching = metadata_samples & stripped_sample_names
    only_in_metadata = metadata_samples - stripped_sample_names
    only_in_expression = stripped_sample_names - metadata_samples
    
    if len(matching) == 0:
        raise ValueError(
            f"No samples match between expression file and metadata. "
            f"Expression samples: {', '.join(sorted(list(stripped_sample_names)[:5]))}... "
            f"Metadata samples: {', '.join(sorted(list(metadata_samples)[:5]))}..."
        )
    
    pct_matched = 100 * len(matching) / len(stripped_sample_names)
    if pct_matched < 90:
        print(
            f"WARNING: Only {pct_matched:.1f}% of expression samples match metadata. "
            f"In metadata but not expression: {', '.join(sorted(list(only_in_metadata)[:3]))}... "
            f"In expression but not metadata: {', '.join(sorted(list(only_in_expression)[:3]))}...",
            file=sys.stderr
        )


def build_metadata_from_groups(group_specs):
    """
    Parse --group-samples strings into a list of (sample_name, group) rows.

    Parameters
    ----------
    group_specs : list[str]
        Each element has the form "GroupLabel:s1,s2,s3".

    Returns
    -------
    list[tuple[str, str]]
        Rows of (sample_name, group) for every sample across all groups.
    """
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
    """Write (sample_name, group) rows to a temp TSV and return its path."""
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
        "src", "R_src", "edger_de.R"
    )

    if not os.path.isfile(r_script):
        print(f"ERROR: R script not found at {r_script}", file=sys.stderr)
        sys.exit(1)

    # Validate Rscript executable ----------------------------------------
    if not os.path.isfile(args.rscript_executable):
        print(
            f"ERROR: Rscript executable not found at: {args.rscript_executable}\n"
            f"       Specify a different path with --rscript-executable",
            file=sys.stderr
        )
        sys.exit(1)

    # Validate parameter ranges ------------------------------------------
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

    # Validate expression file -------------------------------------------
    print("Validating expression file...", file=sys.stderr)
    try:
        sample_names = validate_expression_file(args.expression_file)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Resolve metadata path -----------------------------------------------
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
        # For inline groups, skip metadata validation (metadata is auto-generated and correct)
        sample_names_to_check = set()
        for row in rows:
            sample_names_to_check.add(row[0])
        if not sample_names_to_check <= set(sample_names):
            missing = sample_names_to_check - set(sample_names)
            print(
                f"WARNING: Samples in --group-samples not in expression file: "
                f"{', '.join(sorted(missing))}",
                file=sys.stderr
            )
    else:
        print(
            "ERROR: Provide either --metadata-file or --group-samples.",
            file=sys.stderr
        )
        sys.exit(1)

    # Validate metadata file (if provided, not auto-generated) -----------
    if args.metadata_file:
        print("Validating metadata file...", file=sys.stderr)
        try:
            validate_metadata_file(
                metadata_path, 
                sample_names, 
                sample_suffix=args.sample_suffix
            )
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Build R command ------------------------------------------------------
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
    ]

    print("=" * 70)
    print("SalmonStreamer EdgeRDE")
    print("=" * 70)
    print(f"Expression file : {args.expression_file}")
    print(f"Metadata file   : {metadata_path}")
    print(f"Output dir      : {args.output_dir}")
    print(f"FDR threshold   : {args.fdr_threshold}")
    print(f"logFC threshold : {args.logfc_threshold}")
    if args.sample_suffix:
        print(f"Sample suffix   : {args.sample_suffix}")
    print()

    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: R script failed with exit code {e.returncode}\n"
            f"       Check the error messages above for details.",
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

    print("\nEdgeRDE complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Run edgeR-based differential expression analysis "
            "(SalmonStreamer EdgeRDE subcommand)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Metadata file required columns:
  sample_name   Matches count matrix column headers (after --sample-suffix).
  group         Experimental group label (or: species + tissue columns).

Examples:
  # Using a metadata file
  SalmonStreamer.py EdgeRDE \\
    --expression-file counts.tsv \\
    --output-dir DE_results/ \\
    --metadata-file metadata.tsv

  # Using inline group specification
  SalmonStreamer.py EdgeRDE \\
    --expression-file counts.tsv \\
    --output-dir DE_results/ \\
    --group-samples SpeciesA:s1,s2,s3 SpeciesB:s4,s5,s6

  # Strip a common suffix from sample names before matching metadata
  SalmonStreamer.py EdgeRDE \\
    --expression-file counts.tsv \\
    --metadata-file metadata.tsv \\
    --output-dir DE_results/ \\
    --sample-suffix "_R1_filtered"
        """
    )

    parser.add_argument(
        "--expression-file", required=True,
        help="Tab-separated count matrix (genes x samples, first column = gene IDs)."
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
            "Format: 'GroupLabel:sample1,sample2,...' (one entry per group). "
            "Example: --group-samples SpeciesA:s1,s2 SpeciesB:s3,s4"
        )
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
            "Regex pattern stripped from count matrix column names before "
            "matching to metadata sample_name values "
            "(e.g. '_R1_filtered$'). Optional."
        )
    )
    parser.add_argument(
        "--rscript-executable",
        default=os.path.expanduser("~/.conda/envs/PyR/bin/Rscript"),
        help="Path to the Rscript executable (default: ~/.conda/envs/PyR/bin/Rscript)."
    )

    main(parser.parse_args())
