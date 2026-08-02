#!/usr/bin/env python3

"""
Collapse paralogous genes into single features so a DE analysis can be paralog-aware.

Where reads cannot be assigned between copies, the per-gene count is not a
measurement -- it is one arbitrary partition of a pool. The *sum* over the copies,
however, is well determined: every read in the ambiguous class came from somewhere
in the group. Summing therefore trades resolution you never had for a number that
is actually estimable, and this program prepares that input.

Three things it gets right that a hand-rolled sum does not:

1. GROUPS ARE NOT A PARTITION. ParalogGroups reports every ambiguous read class it
   sees, so a gene can appear in many groups (A+B, B+C, A+B+D...). Summing group by
   group would count the shared genes several times over. Features here are the
   connected components of the gene-sharing graph, so every gene lands in exactly
   one feature.

2. THE FEATURE SET MUST BE FIXED ACROSS SAMPLES. Groups are called per sample and
   pooled; if membership were allowed to vary by sample the columns of the count
   matrix would not be comparable. One global partition is built and applied to
   every sample.

3. AMBIGUITY IS NOT THE SAME AS UNIDENTIFIABILITY. The precision of a paralog's
   own estimate is governed by how many reads are unique to it, not by the fraction
   it shares: the coefficient of variation goes roughly as 1/sqrt(unique reads). Two
   genes that share half their reads but still hold thousands of unique ones are
   perfectly estimable apart, and merging them throws away real signal. Hence
   --min-unique-reads, which merges a group only when some member cannot stand on
   its own.

Counts are summed; TPM is not. An effective length for a merged feature is not
well defined, so for length-aware workflows use --out-tx2gene and let tximport do
the aggregation -- it computes the length offsets correctly, which hand-summing
cannot.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2026-08-02

"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# --------------------------------------------------------------- group input
def load_groups(path):
    """Read a ParalogGroups table, keeping what the merge policy needs."""
    groups = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for required in ("genes", "n_samples_supporting", "per_member_ambiguity"):
            if required not in idx:
                raise ValueError(f"{path}: missing column {required!r}. "
                                 f"Is this a ParalogGroups table?")
        has_unique = "per_member_unique_reads" in idx
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            genes = fields[idx["genes"]].split(",")
            ambiguity = {}
            for token in fields[idx["per_member_ambiguity"]].split(","):
                gene, value = token.rsplit(":", 1)
                ambiguity[gene] = float(value)
            unique = {}
            if has_unique:
                for token in fields[idx["per_member_unique_reads"]].split(","):
                    gene, value = token.rsplit(":", 1)
                    unique[gene] = float(value)
            groups.append({
                "genes": genes,
                "n_samples": int(fields[idx["n_samples_supporting"]]),
                "ambiguity": ambiguity,
                "unique": unique,
                "max_ambiguity": max(ambiguity.values()),
                "min_ambiguity": min(ambiguity.values()),
            })
    return groups


# --------------------------------------------------------------- count input
def read_quant_dirs(quant_dirs, value_column="NumReads"):
    """Build {gene: {sample: value}} from a list of Salmon output directories."""
    matrix = defaultdict(dict)
    samples = []
    for directory in quant_dirs:
        path = os.path.join(directory, "quant.sf")
        if not os.path.exists(path):
            # A Salmon output tree is created at startup, so its presence proves
            # nothing; the terminal artefact is quant.sf.
            print(f"warning: {path} not found; skipping {directory}", file=sys.stderr)
            continue
        sample = os.path.basename(directory.rstrip("/")).replace("_quant", "")
        samples.append(sample)
        with open(path) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            try:
                name_i = header.index("Name")
                value_i = header.index(value_column)
            except ValueError as exc:
                raise ValueError(f"{path}: expected 'Name' and {value_column!r} "
                                 f"columns, found {header}") from exc
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                matrix[fields[name_i]][sample] = float(fields[value_i])
    return matrix, samples


def read_count_matrix(path):
    """Read a combined counts table: first column feature, remaining columns samples."""
    matrix = defaultdict(dict)
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            row = {}
            for sample, value in zip(samples, fields[1:]):
                try:
                    row[sample] = float(value)
                except ValueError:
                    row[sample] = 0.0
            matrix[fields[0]] = row
    return matrix, samples


def collapse_to_genes(matrix, txp2gene):
    """Sum transcript rows onto their gene. A no-op when the annotation is 1:1."""
    out = defaultdict(lambda: defaultdict(float))
    unmapped = set()
    for name, row in matrix.items():
        gene = txp2gene.get(name)
        if gene is None:
            unmapped.add(name)
            gene = name
        for sample, value in row.items():
            out[gene][sample] += value
    if unmapped:
        print(f"note: {len(unmapped)} row ID(s) absent from --txp2gene "
              f"(e.g. {sorted(unmapped)[:3]}); used as their own gene",
              file=sys.stderr)
    return {gene: dict(row) for gene, row in out.items()}


def load_txp2gene(path):
    mapping = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0].startswith("#"):
                continue
            mapping[fields[0]] = fields[1]
    return mapping


# ------------------------------------------------------------------- policy
def infer_group_table_samples(groups):
    """How many samples went into the group table.

    per_member_unique_reads is a TOTAL over the samples ParalogGroups was run on,
    which is not necessarily the number of samples in the count matrix being merged
    -- calibrate on 35 libraries, merge a 4-sample subset, and dividing by the wrong
    denominator moves the estimability threshold by ~9x without any error. Taking the
    largest observed support count recovers the right divisor whenever at least one
    group is supported everywhere, which is the normal case.
    """
    return max((group["n_samples"] for group in groups), default=1)


def group_should_merge(group, args, n_group_table_samples):
    """Apply the merge policy to one called group."""
    if group["n_samples"] < args.min_samples:
        return False, "too few supporting samples"
    if group["max_ambiguity"] < args.min_ambiguity:
        return False, "below ambiguity cutoff"
    if args.require_all_members_ambiguous and group["min_ambiguity"] < args.min_ambiguity:
        return False, "not all members above ambiguity cutoff"
    if args.min_unique_reads > 0 and group["unique"]:
        per_sample = {gene: value / max(n_group_table_samples, 1)
                      for gene, value in group["unique"].items()}
        if min(per_sample.values()) >= args.min_unique_reads:
            return False, "every member individually estimable"
    return True, "merged"


def connected_components(mergeable_groups):
    """Union-find over the retained groups. Returns {gene: frozenset(component)}."""
    parent = {}

    def find(node):
        parent.setdefault(node, node)
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != root:                     # path compression
            parent[node], node = root, parent[node]
        return root

    def union(a, b):
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_a] = root_b

    for group in mergeable_groups:
        genes = group["genes"]
        for gene in genes[1:]:
            union(genes[0], gene)

    members = defaultdict(set)
    for gene in list(parent):
        members[find(gene)].add(gene)
    component_of = {}
    for group_members in members.values():
        frozen = frozenset(group_members)
        for gene in frozen:
            component_of[gene] = frozen
    return component_of


def feature_name(component, style):
    genes = sorted(component)
    if style == "concat":
        return "PARA_" + "+".join(genes)
    return f"PARA_{genes[0]}_n{len(genes)}"


# -------------------------------------------------------------------- output
def main(args):
    groups = load_groups(args.groups)
    print(f"{len(groups)} paralog group(s) loaded from {args.groups}", file=sys.stderr)

    # The feature definition is a property of the groups, not of any particular count
    # matrix, so it can be emitted on its own. That matters when the counts to merge
    # are the production quantification held by someone else: they need the map, not
    # a matrix built from this discovery-pass quant.
    map_only = not args.quant_dirs and not args.counts
    if map_only:
        if not (args.out_map or args.out_tx2gene):
            raise SystemExit(
                "give --counts or --quant-dirs to merge a matrix, or ask for "
                "--out-map / --out-tx2gene alone to emit just the feature definition"
            )
        if not args.txp2gene:
            raise SystemExit("--txp2gene is required when emitting the map alone, "
                             "so genes outside any group are still listed")
        txp2gene = load_txp2gene(args.txp2gene)
        matrix = {gene: {} for gene in set(txp2gene.values())}
        samples = []
        print(f"map-only mode: {len(matrix)} gene(s) from --txp2gene, no counts read",
              file=sys.stderr)
    else:
        if args.quant_dirs:
            matrix, samples = read_quant_dirs(args.quant_dirs, args.value_column)
        else:
            matrix, samples = read_count_matrix(args.counts)
        if not matrix:
            raise SystemExit("no counts were read; nothing to merge")
        if args.txp2gene:
            txp2gene = load_txp2gene(args.txp2gene)
            matrix = collapse_to_genes(matrix, txp2gene)
        print(f"{len(matrix)} gene row(s) x {len(samples)} sample(s)", file=sys.stderr)

    # ---- decide which groups to merge
    n_table_samples = (args.groups_n_samples if args.groups_n_samples > 0
                       else infer_group_table_samples(groups))
    if args.min_unique_reads > 0:
        if not any(group["unique"] for group in groups):
            raise SystemExit(
                "--min-unique-reads needs the per_member_unique_reads column, which "
                "this group table does not have. Re-run ParalogGroups to regenerate "
                "it, or drop the flag to use the ambiguity cutoff alone."
            )
        if args.groups_n_samples <= 0:
            print(f"note: unique-read totals in the group table are being divided by "
                  f"{n_table_samples} (inferred from the maximum support count). Pass "
                  f"--groups-n-samples if ParalogGroups was run on a different number "
                  f"of samples.", file=sys.stderr)
        if not map_only and n_table_samples != len(samples):
            print(f"note: the group table covers {n_table_samples} sample(s) but this "
                  f"matrix has {len(samples)}. That is fine -- the estimability "
                  f"threshold is per sample of the group table.", file=sys.stderr)

    merged_groups, reasons = [], defaultdict(int)
    for group in groups:
        keep, why = group_should_merge(group, args, n_table_samples)
        reasons[why] += 1
        if keep:
            merged_groups.append(group)
    print(f"{len(merged_groups)} group(s) pass the merge policy", file=sys.stderr)
    for why, count in sorted(reasons.items(), key=lambda kv: -kv[1]):
        if why != "merged":
            print(f"    {count:>6} skipped: {why}", file=sys.stderr)

    component_of = connected_components(merged_groups)
    components = {frozenset(c) for c in component_of.values()}
    if components:
        largest = max(len(c) for c in components)
        print(f"{len(component_of)} gene(s) fall into {len(components)} merged "
              f"feature(s); largest holds {largest} genes", file=sys.stderr)
        if args.max_component_size and largest > args.max_component_size:
            oversized = [c for c in components if len(c) > args.max_component_size]
            print(f"warning: {len(oversized)} component(s) exceed "
                  f"--max-component-size {args.max_component_size}. Single linkage "
                  f"chains A-B, B-C into one feature even when A and C share nothing. "
                  f"They are left UNMERGED and listed in the report.", file=sys.stderr)
            for component in oversized:
                for gene in component:
                    component_of.pop(gene, None)
            components = {frozenset(c) for c in component_of.values()}
    else:
        print("no groups passed the policy; the matrix will be written unchanged",
              file=sys.stderr)
        largest = 0

    # ---- build the merged matrix
    merged = {}
    feature_members = {}
    seen_components = {}
    for gene, row in matrix.items():
        component = component_of.get(gene)
        if component is None:
            if gene in merged:
                raise ValueError(f"duplicate gene row {gene!r} in the input matrix")
            merged[gene] = dict(row)
            feature_members[gene] = [gene]
            continue
        name = seen_components.get(component)
        if name is None:
            name = feature_name(component, args.feature_naming)
            seen_components[component] = name
            merged[name] = {sample: 0.0 for sample in samples}
            feature_members[name] = sorted(component)
        for sample, value in row.items():
            merged[name][sample] = merged[name].get(sample, 0.0) + value

    # Components whose members are all absent from the matrix never got created;
    # that is fine, but a component only partly present is worth flagging.
    partial = 0
    for component, name in seen_components.items():
        present = sum(1 for gene in component if gene in matrix)
        if present != len(component):
            partial += 1
    if partial:
        print(f"warning: {partial} merged feature(s) had members missing from the "
              f"count matrix; they were summed over the members present",
              file=sys.stderr)

    # ---- invariant: merging must not create or destroy reads
    before = sum(sum(row.values()) for row in matrix.values())
    after = sum(sum(row.values()) for row in merged.values())
    if abs(before - after) > max(1e-6, 1e-9 * max(before, 1.0)):
        raise AssertionError(
            f"read conservation violated: {before:.4f} before, {after:.4f} after. "
            f"This means a gene was counted into more than one feature."
        )
    print(f"read conservation checked: {before:.0f} counts in, {after:.0f} out",
          file=sys.stderr)

    # ---- write
    n_merged_features = len(seen_components)
    if map_only:
        print(f"{len(merged)} feature(s) ({n_merged_features} merged, "
              f"{len(merged) - n_merged_features} single-gene); no count matrix "
              f"written in map-only mode", file=sys.stderr)
    else:
        with open(args.out_counts, "w") as out:
            out.write("\t".join(["feature"] + samples) + "\n")
            for name in sorted(merged):
                row = merged[name]
                out.write("\t".join([name] + [f"{row.get(s, 0.0):.4f}" for s in samples])
                          + "\n")
        print(f"{len(merged)} feature(s) ({n_merged_features} merged, "
              f"{len(merged) - n_merged_features} single-gene) -> {args.out_counts}",
              file=sys.stderr)

    if args.out_map:
        with open(args.out_map, "w") as out:
            out.write("gene\tfeature\tn_genes_in_feature\tmerged\n")
            for name in sorted(feature_members):
                members = feature_members[name]
                for gene in members:
                    out.write(f"{gene}\t{name}\t{len(members)}\t"
                              f"{'yes' if len(members) > 1 else 'no'}\n")
        print(f"gene -> feature map -> {args.out_map}", file=sys.stderr)

    if args.out_tx2gene:
        if not args.txp2gene:
            print("warning: --out-tx2gene needs --txp2gene to know which transcripts "
                  "belong to which gene; skipped", file=sys.stderr)
        else:
            txp2gene = load_txp2gene(args.txp2gene)
            gene_to_feature = {gene: name
                               for name, members in feature_members.items()
                               for gene in members}
            with open(args.out_tx2gene, "w") as out:
                out.write("transcript\tfeature\n")
                for transcript, gene in sorted(txp2gene.items()):
                    out.write(f"{transcript}\t{gene_to_feature.get(gene, gene)}\n")
            print(f"tximport-ready transcript -> feature map -> {args.out_tx2gene}",
                  file=sys.stderr)

    if args.out_report:
        with open(args.out_report, "w") as out:
            out.write("feature\tn_genes\tgenes\ttotal_counts\n")
            for name in sorted(seen_components.values()):
                members = feature_members[name]
                total = sum(merged[name].values())
                out.write(f"{name}\t{len(members)}\t{','.join(members)}\t"
                          f"{total:.2f}\n")
        print(f"merged-feature report -> {args.out_report}", file=sys.stderr)

    return 0


def add_arguments(parser):
    parser.add_argument(
        "--groups", required=True,
        help="ParalogGroups output table defining the candidate groups.",
    )
    source = parser.add_argument_group("count input (give one)")
    source.add_argument(
        "--counts", default=None,
        help="Combined count matrix: first column the feature ID, one column per "
             "sample (the output of Salmon_output_processor).",
    )
    source.add_argument(
        "--quant-dirs", nargs="+", default=None,
        help="Salmon output directories to read quant.sf from directly.",
    )
    parser.add_argument(
        "--value-column", default="NumReads",
        help="Column to take from quant.sf (default: NumReads). Counts are what DE "
             "tools want; TPM cannot be summed across a merged feature because its "
             "effective length is undefined.",
    )
    parser.add_argument(
        "--txp2gene", default=None,
        help="Two-column TSV mapping transcript to gene, used to collapse a "
             "transcript-level matrix before merging and to emit --out-tx2gene.",
    )

    policy = parser.add_argument_group("merge policy")
    policy.add_argument(
        "--min-ambiguity", type=float, default=0.2,
        help="Merge a group only if its most ambiguous member shares at least this "
             "fraction of its reads within the group. The default of 0.2 comes from "
             "gene-tree calibration; verify it for your data with ParalogTreeCheck.",
    )
    policy.add_argument(
        "--min-samples", type=int, default=1,
        help="Merge a group only if it was supported in at least this many samples "
             "(default: 1).",
    )
    policy.add_argument(
        "--min-unique-reads", type=float, default=0.0,
        help="Merge a group only if some member averages fewer than this many "
             "unique reads per sample -- i.e. only when a member cannot be estimated "
             "on its own. Precision goes as 1/sqrt(unique reads), so ~20 corresponds "
             "to a coefficient of variation near 22 percent. 0 disables the test "
             "(default: 0).",
    )
    policy.add_argument(
        "--groups-n-samples", type=int, default=0,
        help="Number of samples ParalogGroups was run on, used to turn its unique-read "
             "totals into a per-sample rate. 0 infers it from the largest support "
             "count in the table (default: 0).",
    )
    policy.add_argument(
        "--require-all-members-ambiguous", action="store_true",
        help="Require every member to clear --min-ambiguity, not just one. Stricter, "
             "and it keeps a highly ambiguous small gene from dragging in a large "
             "well-determined partner.",
    )
    policy.add_argument(
        "--max-component-size", type=int, default=0,
        help="Refuse to merge components larger than this, and report them instead. "
             "Single linkage chains A-B and B-C into one feature even when A and C "
             "share no reads; this caps the damage. 0 means unlimited (default: 0).",
    )
    policy.add_argument(
        "--feature-naming", choices=["lead", "concat"], default="lead",
        help="'lead' names a feature after its first member plus a size suffix; "
             "'concat' joins every member ID (readable but long) (default: lead).",
    )

    out = parser.add_argument_group("output")
    out.add_argument(
        "-o", "--out-counts", default="paralog_merged_counts.tsv",
        help="Merged count matrix (default: paralog_merged_counts.tsv).",
    )
    out.add_argument(
        "--out-map", default=None,
        help="Write the gene -> feature map, so any DE hit can be traced back to the "
             "genes behind it.",
    )
    out.add_argument(
        "--out-tx2gene", default=None,
        help="Write a transcript -> feature map for tximport. Preferred over the "
             "summed matrix when the downstream model uses length offsets.",
    )
    out.add_argument(
        "--out-report", default=None,
        help="Write one row per merged feature listing its constituent genes.",
    )
    return parser


if __name__ == "__main__":
    standalone_parser = argparse.ArgumentParser(
        description="Collapse paralogous genes into single features for "
                    "paralog-aware differential expression.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    add_arguments(standalone_parser)
    sys.exit(main(standalone_parser.parse_args()))
