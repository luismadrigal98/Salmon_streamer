#!/usr/bin/env python3

"""
ASEIntegrate — Allele-Specific Expression + Differential Expression Integration
SalmonStreamer Pipeline Module

Wraps src/R_src/ase_de_integrate.R.  Takes the DE results produced by
EdgeRDE and per-gene ASE read counts produced by the EDGE_Penstemon
ase_read_counter.py script, then:

  1. (Always) Simple categorisation of every gene with both DE and ASE data:
       cis, trans_compensatory, trans_only, compensatory, no_divergence

  2. (Optional) Advanced Ad/Ed regulatory classification (Wittkopp et al.)
     when expression matrix + metadata + hybrid/parent group labels are given:
       Cis + Trans (Reinforcing), Cis x Trans (Opposing), Trans-only, Ambiguous

ASE counts file format (--ase-counts-file, tab-separated):
  REQUIRED COLUMNS:
    gene_id           Gene identifier.
    snp_count         Number of diagnostic SNPs in the gene.
    bias_ratio        Fraction of SNPs showing allelic bias.
    predominant_bias  'parent1_label' or 'parent2_label' (see --parent1-label).
    <p1>_ratio        Fraction of reads from parent 1 allele (column name is
                      <--parent1-label>_ratio, e.g. "barbatus_ratio").
    <p2>_ratio        Fraction of reads from parent 2 allele.

  This format matches the output of ase_read_counter.py from EDGE_Penstemon.

@Author: Luis Javier Madrigal-Roca
@Date:   2026-04-30
"""

import argparse
import os
import sys
import subprocess
import glob
import re


def validate_de_results_dir(de_dir):
    """
    Validate that DE results directory exists and contains DE result files.
    
    Parameters
    ----------
    de_dir : str
        Path to DE results directory.
    
    Returns
    -------
    list[str]
        List of found DE result files.
    
    Raises
    ------
    ValueError
        If directory doesn't exist or contains no DE results files.
    """
    if not os.path.isdir(de_dir):
        raise ValueError(f"DE results directory not found: {de_dir}")
    
    de_files = glob.glob(os.path.join(de_dir, "*_DE_results.tsv"))
    if not de_files:
        raise ValueError(
            f"No *_DE_results.tsv files found in: {de_dir}\n"
            f"       Run EdgeRDE first to generate DE results."
        )
    
    return de_files


def validate_ase_file(filepath, parent1_label, parent2_label):
    """
    Validate ASE counts file format and required columns.
    
    Parameters
    ----------
    filepath : str
        Path to ASE counts file.
    parent1_label : str
        Label for parent 1 (used to construct ratio column name).
    parent2_label : str
        Label for parent 2 (used to construct ratio column name).
    
    Raises
    ------
    ValueError
        If file format is invalid or required columns are missing.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except IOError as e:
        raise ValueError(f"Cannot read ASE counts file: {e}")
    
    if len(lines) < 2:
        raise ValueError("ASE counts file must have at least 2 rows (header + 1 gene)")
    
    # Parse header
    try:
        header = lines[0].strip().split('\t')
    except Exception as e:
        raise ValueError(f"Cannot parse ASE file header: {e}")
    
    header_lower = [h.lower() for h in header]
    
    # Check required base columns
    required_cols = ['gene_id', 'snp_count', 'bias_ratio', 'predominant_bias']
    missing_cols = [col for col in required_cols if col.lower() not in header_lower]
    if missing_cols:
        raise ValueError(
            f"ASE counts file is missing required columns: {', '.join(missing_cols)}\n"
            f"       Found columns: {', '.join(header)}"
        )
    
    # Check for ratio columns
    p1_ratio_col = f"{parent1_label}_ratio"
    p2_ratio_col = f"{parent2_label}_ratio"
    
    if p1_ratio_col not in header:
        raise ValueError(
            f"ASE counts file is missing column '{p1_ratio_col}'\n"
            f"       (based on --parent1-label '{parent1_label}')\n"
            f"       Found columns: {', '.join(header)}"
        )
    
    if p2_ratio_col not in header:
        raise ValueError(
            f"ASE counts file is missing column '{p2_ratio_col}'\n"
            f"       (based on --parent2-label '{parent2_label}')\n"
            f"       Found columns: {', '.join(header)}"
        )
    
    # Validate data rows
    p1_idx = header.index(p1_ratio_col)
    p2_idx = header.index(p2_ratio_col)
    snp_idx = [i for i, h in enumerate(header_lower) if h == 'snp_count'][0]
    bias_idx = [i for i, h in enumerate(header_lower) if h == 'bias_ratio'][0]
    
    for i, line in enumerate(lines[1:min(6, len(lines))], start=2):  # Check first few rows
        fields = line.strip().split('\t')
        if len(fields) < len(header):
            raise ValueError(
                f"ASE file row {i} has {len(fields)} columns but header has {len(header)}"
            )
        
        # Validate numeric fields
        try:
            snp_count = int(fields[snp_idx])
        except ValueError:
            raise ValueError(
                f"ASE file row {i}: 'snp_count' should be integer, got '{fields[snp_idx]}'"
            )
        
        try:
            bias_ratio = float(fields[bias_idx])
            if not (0.0 <= bias_ratio <= 1.0):
                raise ValueError(f"out of range [0, 1]: {bias_ratio}")
        except ValueError as e:
            raise ValueError(
                f"ASE file row {i}: 'bias_ratio' should be float in [0, 1], "
                f"got '{fields[bias_idx]}' ({e})"
            )
        
        try:
            p1_ratio = float(fields[p1_idx])
            if not (0.0 <= p1_ratio <= 1.0):
                raise ValueError(f"out of range [0, 1]: {p1_ratio}")
        except ValueError as e:
            raise ValueError(
                f"ASE file row {i}: '{p1_ratio_col}' should be float in [0, 1], "
                f"got '{fields[p1_idx]}' ({e})"
            )
        
        try:
            p2_ratio = float(fields[p2_idx])
            if not (0.0 <= p2_ratio <= 1.0):
                raise ValueError(f"out of range [0, 1]: {p2_ratio}")
        except ValueError as e:
            raise ValueError(
                f"ASE file row {i}: '{p2_ratio_col}' should be float in [0, 1], "
                f"got '{fields[p2_idx]}' ({e})"
            )


def main(args):
    r_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "src", "R_src", "ase_de_integrate.R"
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

    # Validate DE results directory ----------------------------------------
    print("Validating DE results directory...", file=sys.stderr)
    try:
        de_files = validate_de_results_dir(args.de_results_dir)
        print(f"Found {len(de_files)} DE results file(s)", file=sys.stderr)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate ASE counts file -------------------------------------------
    print("Validating ASE counts file...", file=sys.stderr)
    try:
        validate_ase_file(args.ase_counts_file, args.parent1_label, args.parent2_label)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate FDR threshold ------------------------------------------
    if args.fdr_threshold < 0 or args.fdr_threshold > 1:
        print(
            f"ERROR: FDR threshold must be between 0 and 1 (got {args.fdr_threshold})",
            file=sys.stderr
        )
        sys.exit(1)

    # Validate advanced-mode dependencies (all-or-nothing) ----------------
    advanced_args = [
        args.expression_file,
        args.metadata_file,
        args.hybrid_group,
        args.parent1_group,
        args.parent2_group,
    ]
    n_advanced = sum(a is not None for a in advanced_args)
    if 0 < n_advanced < len(advanced_args):
        print(
            "ERROR: Advanced Ad/Ed analysis requires ALL of: "
            "--expression-file, --metadata-file, --hybrid-group, "
            "--parent1-group, --parent2-group.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Validate advanced mode files if provided
    if n_advanced == len(advanced_args):
        if not os.path.isfile(args.expression_file):
            print(
                f"ERROR: Expression file not found: {args.expression_file}",
                file=sys.stderr
            )
            sys.exit(1)
        if not os.path.isfile(args.metadata_file):
            print(
                f"ERROR: Metadata file not found: {args.metadata_file}",
                file=sys.stderr
            )
            sys.exit(1)
        print("Advanced Ad/Ed analysis enabled", file=sys.stderr)

    def r_arg(val):
        return val if val is not None else "NULL"

    os.makedirs(args.output_dir, exist_ok=True)

    cmd = [
        args.rscript_executable,
        r_script,
        args.de_results_dir,
        args.ase_counts_file,
        args.output_dir,
        str(args.fdr_threshold),
        r_arg(args.expression_file),
        r_arg(args.metadata_file),
        r_arg(args.hybrid_group),
        r_arg(args.parent1_group),
        r_arg(args.parent2_group),
        args.parent1_label,
        args.parent2_label,
    ]

    print("=" * 70)
    print("SalmonStreamer ASEIntegrate")
    print("=" * 70)
    print(f"DE results dir  : {args.de_results_dir}")
    print(f"ASE counts file : {args.ase_counts_file}")
    print(f"Output dir      : {args.output_dir}")
    print(f"FDR threshold   : {args.fdr_threshold}")
    print(f"Parent 1 label  : {args.parent1_label}")
    print(f"Parent 2 label  : {args.parent2_label}")
    if args.expression_file:
        print(f"Expression file : {args.expression_file}")
        print(f"Metadata file   : {args.metadata_file}")
        print(f"Hybrid group    : {args.hybrid_group}")
        print(f"Parent 1 group  : {args.parent1_group}")
        print(f"Parent 2 group  : {args.parent2_group}")
    print()

    try:
        subprocess.run(cmd, check=True, capture_output=False)
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

    print("\nASEIntegrate complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Integrate ASE read counts with EdgeRDE differential expression results "
            "(SalmonStreamer ASEIntegrate subcommand)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
ASE counts file required columns:
  gene_id          Gene identifier.
  snp_count        Number of diagnostic SNPs.
  bias_ratio       Fraction of SNPs showing allelic bias.
  predominant_bias Which parent allele dominates (value = parent1_label or parent2_label).
  <p1>_ratio       Parent 1 allele read fraction (column = <--parent1-label>_ratio).
  <p2>_ratio       Parent 2 allele read fraction (column = <--parent2-label>_ratio).

Examples:
  # Simple categorisation only (no expression matrix needed)
  SalmonStreamer.py ASEIntegrate \\
    --de-results-dir DE_results/ \\
    --ase-counts-file ase_counts_per_gene.txt \\
    --output-dir ASE_results/ \\
    --parent1-label barbatus \\
    --parent2-label virgatus

  # Full Ad/Ed analysis (requires expression matrix and group labels)
  SalmonStreamer.py ASEIntegrate \\
    --de-results-dir DE_results/ \\
    --ase-counts-file ase_counts_per_gene.txt \\
    --output-dir ASE_results/ \\
    --parent1-label barbatus --parent2-label virgatus \\
    --expression-file counts.tsv \\
    --metadata-file metadata.tsv \\
    --hybrid-group F1_hybrid \\
    --parent1-group P_barbatus_Leaf \\
    --parent2-group P_virgatus_Leaf
        """
    )

    # Required arguments
    parser.add_argument(
        "--de-results-dir", required=True,
        help="Directory containing *_DE_results.tsv files from EdgeRDE."
    )
    parser.add_argument(
        "--ase-counts-file", required=True,
        help="Per-gene ASE read counts file (ase_counts_per_gene.txt format)."
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory for all output files."
    )

    # Labelling
    parser.add_argument(
        "--parent1-label", default="parent1",
        help=(
            "Value in the 'predominant_bias' column representing parent 1, "
            "and prefix for the ratio column (e.g. 'barbatus' → 'barbatus_ratio'). "
            "Default: 'parent1'."
        )
    )
    parser.add_argument(
        "--parent2-label", default="parent2",
        help=(
            "Value in the 'predominant_bias' column representing parent 2, "
            "and prefix for the ratio column. Default: 'parent2'."
        )
    )

    # Optional: thresholds
    parser.add_argument(
        "--fdr-threshold", type=float, default=0.05,
        help="FDR significance cutoff for DE genes (default: 0.05)."
    )

    # Optional: advanced Ad/Ed analysis (all must be supplied together)
    parser.add_argument(
        "--expression-file", default=None,
        help=(
            "Count matrix TSV (same format as --expression-file in EdgeRDE). "
            "Required for advanced Ad/Ed regulatory classification."
        )
    )
    parser.add_argument(
        "--metadata-file", default=None,
        help=(
            "Sample metadata TSV (same format as EdgeRDE). "
            "Required when --expression-file is provided."
        )
    )
    parser.add_argument(
        "--hybrid-group", default=None,
        help=(
            "Group label in metadata identifying hybrid/F1 samples "
            "(e.g. 'F1_hybrid'). Required with --expression-file."
        )
    )
    parser.add_argument(
        "--parent1-group", default=None,
        help=(
            "Group label for parent 1 samples in metadata "
            "(e.g. 'P_barbatus_Leaf'). Required with --expression-file."
        )
    )
    parser.add_argument(
        "--parent2-group", default=None,
        help=(
            "Group label for parent 2 samples in metadata "
            "(e.g. 'P_virgatus_Leaf'). Required with --expression-file."
        )
    )

    parser.add_argument(
        "--rscript-executable",
        default=os.path.expanduser("~/.conda/envs/PyR/bin/Rscript"),
        help="Path to the Rscript executable (default: ~/.conda/envs/PyR/bin/Rscript)."
    )

    main(parser.parse_args())
