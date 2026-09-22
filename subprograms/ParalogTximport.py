#!/usr/bin/env python3

"""
ParalogTximport — aggregate Salmon quantifications to paralog-aware features.

Wraps src/R_src/tximport_features.R.

ParalogMerge can already emit a summed count matrix, and for many purposes that is
enough. What summing cannot produce is the **average transcript length** of a merged
feature, and that is exactly the quantity edgeR and DESeq2 use to correct for
composition differences between samples. tximport recomputes it from the Salmon
output, which is why it is the preferred route -- and why it should not require
anyone to open R and write the call by hand.

Two ways to use the result, selected with --counts-from-abundance:

  lengthScaledTPM (default)
      The length correction is folded into the counts. The matrix can be passed
      straight to EdgeRDE, or any other count-based DE tool, with no offset. This
      is the route that keeps the pipeline self-contained.
  no
      Raw counts plus a separate length matrix, for a downstream model that takes
      the length matrix as an offset. More faithful, but the consumer has to
      support offsets.

Outputs, written into --output-dir:
  paralog_feature_counts.tsv         features x samples
  paralog_feature_abundance_tpm.tsv  TPM
  paralog_feature_length.tsv         average transcript length per feature/sample

Requires the R package tximport (Bioconductor).

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2026-08-14

"""

import argparse
import csv
import os
import shutil
import re
import subprocess
import sys
import tempfile


def derive_sample_name(quant_path, strip_suffix=None):
    """Sample name for a quant.sf path: the directory that contains it.

    Salmon writes <sample_dir>/quant.sf, so the containing directory is the
    sample. Callers whose directories carry a decoration ('_quant', '_R1_filtered')
    can strip it with --strip-suffix.
    """
    name = os.path.basename(os.path.dirname(os.path.abspath(quant_path)))
    if strip_suffix:
        name = re.sub(strip_suffix, "", name)
    return name


def resolve_quant_files(args):
    """Turn --quant-dirs / --quant-files into [(sample_name, quant_path)]."""
    paths = []
    if args.quant_files:
        paths = list(args.quant_files)
    else:
        for d in args.quant_dirs:
            candidate = d if os.path.isfile(d) else os.path.join(d, "quant.sf")
            paths.append(candidate)

    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        print(f"ERROR: {len(missing)} quant.sf file(s) not found. First: {missing[0]}",
              file=sys.stderr)
        print("       --quant-dirs takes Salmon output directories (each holding "
              "quant.sf).", file=sys.stderr)
        sys.exit(1)

    if args.sample_names:
        if len(args.sample_names) != len(paths):
            print(f"ERROR: --sample-names has {len(args.sample_names)} entries but "
                  f"{len(paths)} quant file(s) were given.", file=sys.stderr)
            sys.exit(1)
        names = list(args.sample_names)
    else:
        names = [derive_sample_name(p, args.strip_suffix) for p in paths]

    duplicates = {n for n in names if names.count(n) > 1}
    if duplicates:
        print(f"ERROR: duplicate sample name(s) after derivation: "
              f"{sorted(duplicates)}", file=sys.stderr)
        print("       Pass --sample-names explicitly, or use --strip-suffix.",
              file=sys.stderr)
        sys.exit(1)

    return list(zip(names, paths))


def write_manifest(pairs, tmp_dir):
    """Write sample_name/quant_path rows to a temp TSV for the R script."""
    fh = tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=tmp_dir, delete=False, newline=""
    )
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample_name", "quant_path"])
    writer.writerows(pairs)
    fh.close()
    return fh.name


def main(args):
    if not getattr(args, "quant_dirs", None) and not getattr(args, "quant_files", None):
        print("ERROR: give --quant-dirs or --quant-files.", file=sys.stderr)
        sys.exit(1)

    r_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "src", "R_src", "tximport_features.R"
    )
    if not os.path.isfile(r_script):
        print(f"ERROR: R script not found at {r_script}", file=sys.stderr)
        sys.exit(1)

    rscript = shutil.which(os.path.expanduser(args.rscript_executable))
    if rscript is None:
        print(
            f"ERROR: Rscript executable not found: {args.rscript_executable}\n"
            f"       Put Rscript on your PATH or pass --rscript-executable",
            file=sys.stderr
        )
        sys.exit(1)
    args.rscript_executable = rscript

    if not os.path.isfile(args.tx2feature):
        print(f"ERROR: tx2feature map not found: {args.tx2feature}\n"
              f"       This is ParalogMerge's --out-tx2gene output.", file=sys.stderr)
        sys.exit(1)

    pairs = resolve_quant_files(args)
    os.makedirs(args.output_dir, exist_ok=True)
    manifest_path = write_manifest(pairs, args.output_dir)

    cmd = [
        args.rscript_executable,
        r_script,
        manifest_path,
        args.tx2feature,
        args.output_dir,
        args.counts_from_abundance,
        "TRUE" if args.ignore_tx_version else "FALSE",
    ]

    print("=" * 70)
    print("SalmonStreamer ParalogTximport")
    print("=" * 70)
    print(f"Samples          : {len(pairs)}")
    print(f"tx2feature       : {args.tx2feature}")
    print(f"Output dir       : {args.output_dir}")
    print(f"countsFromAbundance: {args.counts_from_abundance}")
    print()

    try:
        subprocess.run(cmd, check=True, capture_output=False)
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: R script failed with exit code {e.returncode}\n"
            f"       If the message mentions 'there is no package called tximport',\n"
            f"       install it with:\n"
            f'         R -e \'BiocManager::install("tximport")\'',
            file=sys.stderr
        )
        sys.exit(e.returncode)
    except FileNotFoundError as e:
        print(f"ERROR: Cannot find Rscript executable: "
              f"{args.rscript_executable}\n       {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        if os.path.isfile(manifest_path):
            os.remove(manifest_path)

    counts = os.path.join(args.output_dir, "paralog_feature_counts.tsv")
    print("\nParalogTximport complete.")
    print(f"  Feed to DE:  python SalmonStreamer.py EdgeRDE \\")
    print(f"                   --expression-file {counts} \\")
    print(f"                   --metadata-file metadata.tsv --output-dir DE_results/")
    return 0


def add_arguments(parser):
    source = parser.add_argument_group("quantification input (give one)")
    source.add_argument(
        "--quant-dirs", nargs="+", default=None,
        help="Salmon output directories, each containing quant.sf.",
    )
    source.add_argument(
        "--quant-files", nargs="+", default=None,
        help="quant.sf files directly, if the layout is non-standard.",
    )
    parser.add_argument(
        "--tx2feature", required=True,
        help="Transcript -> feature map from ParalogMerge --out-tx2gene. Using the "
             "plain txp2gene here instead would give an ordinary gene-level matrix "
             "with no paralog merging.",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory for the count, TPM and length matrices.",
    )
    parser.add_argument(
        "--counts-from-abundance", default="lengthScaledTPM",
        choices=["lengthScaledTPM", "scaledTPM", "no"],
        help="'lengthScaledTPM' (default) folds the transcript-length correction "
             "into the counts so the matrix goes straight to EdgeRDE with no "
             "offset. 'no' keeps raw counts and leaves the correction to a "
             "downstream tool that accepts the length matrix as an offset.",
    )
    parser.add_argument(
        "--sample-names", nargs="+", default=None,
        help="Sample names, in the same order as the quant inputs. Defaults to each "
             "quant.sf's containing directory name.",
    )
    parser.add_argument(
        "--strip-suffix", default=None,
        help=r"Regex stripped from derived sample names, e.g. '_quant$'. Ignored "
             r"when --sample-names is given.",
    )
    parser.add_argument(
        "--ignore-tx-version", action="store_true",
        help="Drop trailing transcript version suffixes ('.1') when matching "
             "transcripts to the map.",
    )
    parser.add_argument(
        "--rscript-executable",
        default="Rscript",
        help="Path to the Rscript executable (default: Rscript on your PATH).",
    )
    return parser


if __name__ == "__main__":
    standalone_parser = argparse.ArgumentParser(
        description="Aggregate Salmon quantifications to paralog-aware features "
                    "with tximport.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    add_arguments(standalone_parser)
    parsed = standalone_parser.parse_args()
    if not parsed.quant_dirs and not parsed.quant_files:
        standalone_parser.error("give --quant-dirs or --quant-files")
    sys.exit(main(parsed))
