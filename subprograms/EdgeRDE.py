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
    else:
        print(
            "ERROR: Provide either --metadata-file or --group-samples.",
            file=sys.stderr
        )
        sys.exit(1)

    # Validate expression file --------------------------------------------
    if not os.path.isfile(args.expression_file):
        print(f"ERROR: Expression file not found: {args.expression_file}",
              file=sys.stderr)
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
        result = subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: R script exited with code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)
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
