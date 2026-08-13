#!/usr/bin/env python3

"""
Build IGV artefacts from a paralog group table so groups can be inspected by eye.

A paralog group is a claim about read ambiguity, not about the genome. Before acting
on one it is worth confirming it looks like a duplication in the alignments, and that
is a visual check: load the tracks written here alongside the BAMs and see whether the
shared reads really do pile over every member.

Consumes the table written by ParalogGroups (the `--output` TSV) and produces:

  <prefix>.bed        one feature per gene, named grp<NNN>_amb<X.XX>_<gene>, with the
                      BED score set to ambiguity * 1000 so IGV shades the strongly
                      ambiguous groups darker. Sorted by scaffold then start.
  <prefix>_loci.txt   one padded locus per group, ready to paste into the IGV location
                      box. Groups whose members sit on different scaffolds get one
                      locus per member instead of a single span.

Gene coordinates are read from the GFF3 used to build the transcriptome (`gene`
features, `ID=<gene>`); the gene IDs there must match those in the group table, which
is what `--txp2gene` already guarantees if the same annotation produced both.

When looking at the result in IGV, load the plain aligned BAM rather than a
secondary-filtered one -- the multimapping reads are the entire point, and a filtered
BAM has removed exactly the evidence you are trying to see. Aligners mark reads they
could not place uniquely with MAPQ 0 (STAR uses 0-3, reserving 255 for unique), and
IGV draws those white or transparent, so a genuine group reads as a pale pile-up
spanning every copy.
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from subprograms.ParalogGroups import _open_maybe_gzip          # noqa: E402

REQUIRED_COLUMNS = ("genes", "shared_reads", "n_samples_supporting",
                    "min_member_ambiguity")


def load_gene_coords(gff_path):
    """Map gene ID -> (scaffold, start, end, strand) from the `gene` rows of a GFF3.

    Coordinates stay 1-based inclusive here, exactly as the GFF states them; the BED
    conversion happens at the point of writing so the loci file can reuse them
    untouched.
    """
    coords = {}
    with _open_maybe_gzip(gff_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gene_id = fields[8]
            if gene_id.startswith("ID="):
                gene_id = gene_id[3:]
            gene_id = gene_id.split(";")[0]
            coords[gene_id] = (fields[0], int(fields[3]), int(fields[4]), fields[6])
    if not coords:
        raise ValueError(
            f"{gff_path}: no `gene` features found. ParalogTracks reads gene rows "
            f"with an ID= attribute; check this is the GFF3 the transcriptome was "
            f"built from and not a transcript-only or GTF-format file."
        )
    return coords


def select_groups(rows, min_ambiguity, min_samples):
    """Filter the group table, most-shared-reads first.

    Ordering by shared reads is what makes the group numbering meaningful: grp001 is
    the group carrying the most ambiguous signal, so scanning the loci file top-down
    walks the groups in descending order of how much they actually matter.
    """
    selected = [
        row for row in rows
        if float(row["min_member_ambiguity"]) >= min_ambiguity
        and int(row["n_samples_supporting"]) >= min_samples
    ]
    selected.sort(key=lambda row: -float(row["shared_reads"]))
    return selected


def build_bed_features(selected, coords):
    """Return (bed_rows, missing_genes) for the selected groups."""
    bed_rows = []
    missing = []
    for index, row in enumerate(selected, 1):
        ambiguity = float(row["min_member_ambiguity"])
        for gene in row["genes"].split(","):
            if gene not in coords:
                missing.append(gene)
                continue
            scaffold, start, end, strand = coords[gene]
            bed_rows.append((
                scaffold,
                start - 1,                                  # BED is 0-based half-open
                end,
                f"grp{index:03d}_amb{ambiguity:.2f}_{gene}",
                min(1000, int(ambiguity * 1000)),
                strand,
            ))
    bed_rows.sort(key=lambda feature: (feature[0], feature[1]))
    return bed_rows, missing


def group_locus(row, coords, pad_fraction=10, min_pad=2000):
    """Padded IGV locus for one group, or per-member loci if it spans scaffolds."""
    genes = [gene for gene in row["genes"].split(",") if gene in coords]
    if not genes:
        return ""
    scaffolds = {coords[gene][0] for gene in genes}
    if len(scaffolds) > 1:
        return " | ".join(
            f"{coords[gene][0]}:{coords[gene][1]}-{coords[gene][2]}" for gene in genes
        )
    scaffold = scaffolds.pop()
    low = min(coords[gene][1] for gene in genes)
    high = max(coords[gene][2] for gene in genes)
    pad = max(min_pad, (high - low) // pad_fraction)
    return f"{scaffold}:{max(1, low - pad)}-{high + pad}"


def add_arguments(parser):
    """Register this subcommand's options (shared by the standalone and CLI entry points)."""
    parser.add_argument(
        "--groups", required=True,
        help="Paralog group table written by ParalogGroups (its --output TSV).",
    )
    parser.add_argument(
        "--gff", required=True,
        help="GFF3 with gene features, used to look up coordinates. Should be the "
             "same annotation the transcriptome was built from.",
    )
    parser.add_argument(
        "--prefix", default="paralog_groups_tracks",
        help="Output prefix; writes <prefix>.bed and <prefix>_loci.txt "
             "(default: paralog_groups_tracks).",
    )
    parser.add_argument(
        "--min-ambiguity", type=float, default=0.0,
        help="Keep only groups whose least-ambiguous member has at least this "
             "fraction of its reads shared (default: 0, draw everything in the "
             "table). Calibrate this against gene trees with ParalogTreeCheck rather "
             "than guessing -- values near 0.9 look reasonable but in practice can "
             "exclude every tree-confirmed duplication in a dataset.",
    )
    parser.add_argument(
        "--min-samples", type=int, default=1,
        help="Keep only groups supported in at least this many samples (default: 1). "
             "This is a count, not a fraction, so it has to be set relative to how "
             "many samples went into the table.",
    )
    return parser


def main(args):
    coords = load_gene_coords(args.gff)

    with open(args.groups, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing_columns:
            print(
                f"error: {args.groups} lacks required column(s) "
                f"{', '.join(missing_columns)}. This should be the TSV written by "
                f"ParalogGroups.",
                file=sys.stderr,
            )
            return 1
        rows = list(reader)

    selected = select_groups(rows, args.min_ambiguity, args.min_samples)
    if not selected:
        print(
            f"warning: no groups passed --min-ambiguity {args.min_ambiguity} / "
            f"--min-samples {args.min_samples} ({len(rows)} in the table). Nothing "
            f"was written; loosen the cutoffs or calibrate them with ParalogTreeCheck.",
            file=sys.stderr,
        )
        return 1

    bed_rows, missing = build_bed_features(selected, coords)

    bed_path = f"{args.prefix}.bed"
    with open(bed_path, "w") as out:
        out.write(
            f'track name="{os.path.basename(args.prefix)}" '
            f'description="Salmon eq-class paralog groups, '
            f'ambiguity>={args.min_ambiguity}, >={args.min_samples} samples" '
            f'itemRgb="Off"\n'
        )
        for feature in bed_rows:
            out.write("\t".join(map(str, feature)) + "\n")

    loci_path = f"{args.prefix}_loci.txt"
    with open(loci_path, "w") as out:
        out.write("# paste a line's locus into the IGV location box\n")
        out.write("# load the plain aligned BAM, not a secondary-filtered one:\n")
        out.write("# the multimapping reads are the evidence, and filtering removes them\n")
        out.write("group\tambiguity\tlocus\tgenes\n")
        for index, row in enumerate(selected, 1):
            out.write(
                f"grp{index:03d}\t{row['min_member_ambiguity']}\t"
                f"{group_locus(row, coords)}\t{row['genes']}\n"
            )

    print(
        f"{len(selected)} group(s) -> {len(bed_rows)} BED features -> "
        f"{bed_path}, {loci_path}",
        file=sys.stderr,
    )
    if missing:
        # De-duplicated because one absent gene can appear in several groups, and a
        # systematic ID mismatch would otherwise print the same name hundreds of times.
        unique_missing = list(OrderedDict.fromkeys(missing))
        preview = ", ".join(unique_missing[:5])
        suffix = ", ..." if len(unique_missing) > 5 else ""
        print(
            f"warning: {len(unique_missing)} gene(s) had no coordinates in "
            f"{args.gff} and were skipped ({preview}{suffix}). If that is most of "
            f"them, the group table and the GFF are using different gene IDs.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    standalone_parser = argparse.ArgumentParser(
        description="Build IGV tracks and loci from a paralog group table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    add_arguments(standalone_parser)
    sys.exit(main(standalone_parser.parse_args()))
