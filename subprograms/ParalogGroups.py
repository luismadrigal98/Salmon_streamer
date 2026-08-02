#!/usr/bin/env python3

"""
Identify sets of paralogs from Salmon equivalence classes.

The criterion is the one Haylee and Lena want: group genes when RNAseq reads map
equally well to every copy. Salmon already computes exactly that structure and can
write it out directly -- there is no need to parse the SAM from --writeMappings.

Run the quantification with:

    salmon quant -i <index> -l A -1 r1.fq -2 r2.fq -o <outdir> \
        --dumpEq --dumpEqWeights --hardFilter

  --dumpEq        writes <outdir>/aux_info/eq_classes.txt
  --dumpEqWeights adds the per-transcript conditional probabilities to each line
  --hardFilter    keeps ONLY the equally-highest-scoring mappings per fragment,
                  which is what makes an equivalence class label mean "these reads
                  map equally well to exactly this set of transcripts"

Each line of eq_classes.txt is one equivalence class: the number of transcripts in
its label, the transcript IDs, (optionally the conditional probabilities), and the
number of fragments in the class. A class whose label holds >= 2 transcripts is a
direct observation of ambiguous reads; the label itself is the candidate group.

Note that with range factorization the same transcript set can appear on several
lines (different fragments induce different conditional probabilities), so counts
are aggregated by label here.
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from subprograms.ParalogTreeCheck import (          # noqa: E402
    add_tree_arguments, groups_from_rows, run_calibration,
)


def _open_maybe_gzip(path):
    """Salmon >=1.x writes aux_info/eq_classes.txt.gz; older versions wrote it plain.
    Sniff the magic number rather than trusting the extension."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path)


def parse_eq_classes(path):
    """Read a Salmon eq_classes.txt (plain or gzipped).

    Returns (transcript_names, {frozenset(label_indices): total_count}).
    Handles both the --dumpEq and the --dumpEq --dumpEqWeights layouts by
    inferring which is present from the field count of each line.
    """
    with _open_maybe_gzip(path) as handle:
        try:
            n_txp = int(next(handle).strip())
            n_eq = int(next(handle).strip())
        except (StopIteration, ValueError) as exc:
            raise ValueError(f"{path}: malformed header (expected two integers)") from exc

        names = [next(handle).strip() for _ in range(n_txp)]

        labels = defaultdict(float)
        seen = 0
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            seen += 1
            k = int(fields[0])
            # Without weights: 1 + k + 1 fields. With weights: 1 + k + k + 1.
            if len(fields) == k + 2:
                pass
            elif len(fields) == 2 * k + 2:
                pass
            else:
                raise ValueError(
                    f"{path} line {seen + n_txp + 2}: {len(fields)} fields is neither "
                    f"k+2 nor 2k+2 for k={k}"
                )
            idx = frozenset(int(x) for x in fields[1:k + 1])
            labels[idx] += float(fields[-1])

    if seen != n_eq:
        print(
            f"warning: {path} declared {n_eq} equivalence classes but {seen} lines were read",
            file=sys.stderr,
        )
    return names, labels


def load_txp2gene(path):
    """Two-column TSV: transcript_id <tab> gene_id."""
    mapping = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0].startswith("#"):
                continue
            mapping[fields[0]] = fields[1]
    return mapping


def summarize_sample(path, txp2gene):
    """Collapse one sample's equivalence classes to gene-set labels.

    Returns (shared, unique) where shared maps a frozenset of >=2 gene IDs to the
    read count carried by that exact set, and unique maps a gene ID to the count of
    reads that were unambiguous at the gene level.
    """
    names, labels = parse_eq_classes(path)

    unmapped_ids = set()
    shared = defaultdict(float)
    unique = defaultdict(float)

    for idx, count in labels.items():
        genes = set()
        for i in idx:
            txp = names[i]
            if txp2gene is None:
                genes.add(txp)
            elif txp in txp2gene:
                genes.add(txp2gene[txp])
            else:
                unmapped_ids.add(txp)
                genes.add(txp)
        # A label spanning several transcripts of ONE gene is isoform ambiguity,
        # not paralogy -- it collapses to a single gene and counts as unique.
        if len(genes) >= 2:
            shared[frozenset(genes)] += count
        else:
            unique[next(iter(genes))] += count

    if unmapped_ids:
        print(
            f"warning: {os.path.basename(path)}: {len(unmapped_ids)} transcript IDs "
            f"absent from the txp2gene map (e.g. {sorted(unmapped_ids)[:3]}); "
            f"used as their own gene",
            file=sys.stderr,
        )
    return shared, unique


def add_arguments(parser):
    """Register this subcommand's options (shared by the standalone and CLI entry points)."""
    parser.add_argument(
        "--eq-files", nargs="+", required=True,
        help="One or more aux_info/eq_classes.txt files (one per sample).",
    )
    parser.add_argument(
        "--txp2gene", default=None,
        help="Two-column TSV mapping transcript ID to gene ID. Without it, "
             "transcripts are treated as genes.",
    )
    parser.add_argument(
        "-o", "--output", default="paralog_groups.tsv",
        help="Output TSV (default: paralog_groups.tsv).",
    )
    parser.add_argument(
        "--min-reads", type=float, default=10.0,
        help="A group must carry at least this many shared reads in a sample for "
             "that sample to support it (default: 10).",
    )
    parser.add_argument(
        "--min-samples", type=int, default=1,
        help="A group must be supported in at least this many samples (default: 1).",
    )
    parser.add_argument(
        "--gene-table", default=None,
        help="Also write a per-gene TSV giving the raw read split (unique vs shared) "
             "behind each group, so the fractions can be audited directly.",
    )
    parser.add_argument(
        "--min-ambiguity", type=float, default=0.0,
        help="Keep only groups where the least-ambiguous member has at least this "
             "fraction of its reads shared with the rest of the group. Use ~0.9 for "
             "near-identical copies, 0 to report everything (default: 0). Calibrate "
             "this against gene trees with ParalogTreeCheck rather than guessing -- "
             "0.9 is far stricter than tree-confirmed duplicates require.",
    )
    # Optional inline gene-tree calibration. The same analysis is available as the
    # ParalogTreeCheck subcommand, which reads a finished table and so can be re-run
    # in seconds; use that when sweeping cutoffs, and these flags for a one-shot run.
    add_tree_arguments(parser, required_group=False)
    return parser


def main(args):
    per_sample_shared = []
    total_unique = defaultdict(float)
    total_shared_by_gene = defaultdict(float)

    txp2gene = load_txp2gene(args.txp2gene) if args.txp2gene else None

    for path in args.eq_files:
        shared, unique = summarize_sample(path, txp2gene)
        per_sample_shared.append(shared)
        for gene, count in unique.items():
            total_unique[gene] += count
        for label, count in shared.items():
            for gene in label:
                total_shared_by_gene[gene] += count

    # Pool the groups seen in any sample, then score them.
    all_labels = set()
    for shared in per_sample_shared:
        all_labels.update(shared)

    rows = []
    for label in all_labels:
        support = [s.get(label, 0.0) for s in per_sample_shared]
        n_supporting = sum(1 for c in support if c >= args.min_reads)
        if n_supporting < args.min_samples:
            continue
        group_reads = sum(support)

        # Per member: what share of that gene's reads is ambiguous within this group?
        member_ambiguity = {}
        for gene in label:
            denom = total_unique[gene] + total_shared_by_gene[gene]
            member_ambiguity[gene] = (group_reads / denom) if denom > 0 else 0.0
        weakest = min(member_ambiguity.values())
        if weakest < args.min_ambiguity:
            continue

        genes = sorted(label)
        rows.append({
            "n_genes": len(genes),
            "genes": ",".join(genes),
            "shared_reads": group_reads,
            "n_samples_supporting": n_supporting,
            "min_member_ambiguity": weakest,
            "max_member_ambiguity": max(member_ambiguity.values()),
            "per_member_ambiguity": ",".join(
                f"{g}:{member_ambiguity[g]:.3f}" for g in genes
            ),
            # Raw counts behind the fractions, so the split can be audited directly
            # rather than taken on faith.
            "per_member_unique_reads": ",".join(
                f"{g}:{total_unique[g]:.0f}" for g in genes
            ),
            "per_member_total_reads": ",".join(
                f"{g}:{total_unique[g] + total_shared_by_gene[g]:.0f}" for g in genes
            ),
        })

    rows.sort(key=lambda r: (-r["shared_reads"], r["genes"]))

    columns = [
        "n_genes", "genes", "shared_reads", "n_samples_supporting",
        "min_member_ambiguity", "max_member_ambiguity", "per_member_ambiguity",
        "per_member_unique_reads", "per_member_total_reads",
    ]
    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(
                f"{row[c]:.2f}" if isinstance(row[c], float) else str(row[c])
                for c in columns
            ) + "\n")

    # Optional per-gene view: the raw read split behind every group, one row per gene.
    if args.gene_table:
        in_groups = defaultdict(list)
        for r in rows:
            for g in r["genes"].split(","):
                in_groups[g].append(r)
        gcols = ["gene", "unique_reads", "shared_reads", "total_reads",
                 "fraction_shared", "n_groups", "largest_group_ambiguity", "partners"]
        with open(args.gene_table, "w") as out:
            out.write("\t".join(gcols) + "\n")
            for gene in sorted(in_groups):
                uniq = total_unique[gene]
                shar = total_shared_by_gene[gene]
                tot = uniq + shar
                grps = in_groups[gene]
                worst = max(float(r["per_member_ambiguity"].split(f"{gene}:")[1].split(",")[0])
                            for r in grps)
                partners = sorted({p for r in grps for p in r["genes"].split(",") if p != gene})
                out.write("\t".join([
                    gene, f"{uniq:.0f}", f"{shar:.0f}", f"{tot:.0f}",
                    f"{(shar / tot if tot else 0):.3f}", str(len(grps)),
                    f"{worst:.3f}", ",".join(partners),
                ]) + "\n")
        print(f"per-gene read split for {len(in_groups)} genes -> {args.gene_table}",
              file=sys.stderr)

    n_pairs = sum(1 for r in rows if r["n_genes"] == 2)
    print(
        f"{len(rows)} candidate paralog groups from {len(args.eq_files)} sample(s) "
        f"({n_pairs} pairs, {len(rows) - n_pairs} larger) -> {args.output}",
        file=sys.stderr,
    )

    # Optional inline calibration against gene trees.
    if getattr(args, "trees", None):
        if args.min_ambiguity > 0:
            print(
                f"warning: calibrating against trees with --min-ambiguity "
                f"{args.min_ambiguity} already applied. Groups below that cutoff are "
                f"absent from the table, so sensitivity will look artificially low. "
                f"Re-run with --min-ambiguity 0 to calibrate honestly.",
                file=sys.stderr,
            )
        run_calibration(groups_from_rows(rows), args.trees, args)

    return 0


if __name__ == "__main__":
    standalone_parser = argparse.ArgumentParser(
        description="Group genes into paralog sets from Salmon equivalence classes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    add_arguments(standalone_parser)
    sys.exit(main(standalone_parser.parse_args()))
