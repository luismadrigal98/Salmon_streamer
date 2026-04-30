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


def main(args):
    r_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "src", "R_src", "ase_de_integrate.R"
    )

    if not os.path.isfile(r_script):
        print(f"ERROR: R script not found at {r_script}", file=sys.stderr)
        sys.exit(1)

    # Validate required inputs -------------------------------------------
    if not os.path.isdir(args.de_results_dir):
        print(f"ERROR: DE results directory not found: {args.de_results_dir}",
              file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.ase_counts_file):
        print(f"ERROR: ASE counts file not found: {args.ase_counts_file}",
              file=sys.stderr)
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
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: R script exited with code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)

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
