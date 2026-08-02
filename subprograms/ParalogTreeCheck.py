#!/usr/bin/env python3

"""
Ground-truth paralog groups against gene trees, and calibrate the ambiguity cutoff.

ParalogGroups calls paralogs from read behaviour: two genes are grouped when reads
map equally well to both. A gene tree calls paralogs from ancestry. These are not
the same question, and the difference is the point of this program:

  * a gene tree says "these two genes descend from a duplication"
  * an equivalence class says "a 150 bp read cannot tell these two genes apart"

A pair can be a true paralog and still be perfectly separable by reads -- that is a
*correct* non-detection, not a miss. So the trees are the right ground truth for
"is this a paralog?" and the wrong ground truth for "must I merge these counts?".
What the trees are genuinely good for is calibration: they tell you where on the
divergence axis your read-based detector stops firing, and whether that boundary
sits where sequence divergence says it should.

Given one or more Newick trees and a ParalogGroups table, this program:

  1. pulls the tips that belong to the quantified genome (the rest are outgroups),
  2. enumerates every within-genome pair, with its patristic distance and whether
     its MRCA is genome-specific (i.e. a duplication after the split from the other
     taxa in the tree -- the recent duplicates most likely to be unresolvable),
  3. reports which pairs the read-based method detected, at what ambiguity,
  4. sweeps the ambiguity cutoff and reports sensitivity against the tree-confirmed
     set alongside the genome-wide cost of that cutoff.

Tip labels are matched to gene IDs by regex, then resolved against the gene IDs
actually present in the group table, so the assorted label conventions that come
out of tree viewers ('Pkunth_PGA_scaffold5_000853.1 Copy', 'P. kunthii
(scaffold2_004562.1)') all land on the same gene without hand-editing.

Trees exported from a word processor are accepted directly: an RTF wrapper is
stripped before parsing, because that is how they usually arrive.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2026-08-02

"""

import argparse
import os
import re
import sys
from collections import defaultdict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# --------------------------------------------------------------------- input
def strip_rtf(raw):
    """Recover the text payload from an RTF file.

    Trees are routinely saved out of TextEdit/Word, which wraps the Newick string
    in RTF. Detecting and undoing that here is cheaper than asking every collaborator
    to re-export as plain text, and it is a no-op on files that are already plain.
    """
    if not raw.lstrip().startswith(r"{\rtf"):
        return raw
    for table in ("fonttbl", "colortbl", "expandedcolortbl", "stylesheet", "listtable"):
        raw = re.sub(r"\{\\\*?\\" + table + r"[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", "", raw)
    raw = re.sub(r"\\'([0-9a-fA-F]{2})", lambda m: chr(int(m.group(1), 16)), raw)
    raw = re.sub(r"\\[a-zA-Z]+-?\d* ?", " ", raw)
    raw = raw.replace("\\\n", "\n")
    return re.sub(r"[{}]", "", raw)


def extract_newick(text):
    """Return the first balanced parenthetical, terminated at its ';' if present."""
    start = text.find("(")
    if start < 0:
        return None
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                semi = text.find(";", i)
                return text[start:semi + 1] if semi >= 0 else text[start:i + 1]
    return None


# ------------------------------------------------------------------- parsing
class Node:
    """A Newick node. Deliberately minimal -- adding biopython or ete3 to this
    pipeline's dependency set for one recursive-descent parse is a bad trade."""

    __slots__ = ("name", "length", "children", "parent")

    def __init__(self):
        self.name = None
        self.length = 0.0
        self.children = []
        self.parent = None


def parse_newick(text):
    """Parse a Newick string into a tree of Node. Handles quoted labels, branch
    lengths, internal labels and support values."""
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    pos = [0]
    n = len(text)

    def parse_node():
        node = Node()
        if pos[0] < n and text[pos[0]] == "(":
            pos[0] += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if pos[0] < n and text[pos[0]] == ",":
                    pos[0] += 1
                    continue
                break
            if pos[0] < n and text[pos[0]] == ")":
                pos[0] += 1
        if pos[0] < n and text[pos[0]] in "'\"":
            # Gene names in these trees carry unescaped apostrophes -- F3'5'H, 3'UTR.
            # Strict Newick would double them, but hand-built trees do not, and
            # closing at the first inner quote silently truncates the whole tree to
            # one tip. Close only on a quote that sits in a structurally valid
            # position: followed by a delimiter, or doubled as the standard escape.
            quote = text[pos[0]]
            pos[0] += 1
            chars = []
            while pos[0] < n:
                char = text[pos[0]]
                if char != quote:
                    chars.append(char)
                    pos[0] += 1
                    continue
                if pos[0] + 1 < n and text[pos[0] + 1] == quote:
                    chars.append(quote)                  # doubled = literal quote
                    pos[0] += 2
                    continue
                look = pos[0] + 1
                while look < n and text[look] in " \t":
                    look += 1
                if look >= n or text[look] in ",():;":
                    pos[0] += 1                          # genuine closing quote
                    break
                chars.append(char)                       # embedded literal quote
                pos[0] += 1
            node.name = "".join(chars)
        else:
            start = pos[0]
            while pos[0] < n and text[pos[0]] not in ",():;":
                pos[0] += 1
            label = text[start:pos[0]].strip()
            if label:
                node.name = label
        if pos[0] < n and text[pos[0]] == ":":
            pos[0] += 1
            start = pos[0]
            while pos[0] < n and text[pos[0]] not in ",():;":
                pos[0] += 1
            try:
                node.length = float(text[start:pos[0]])
            except ValueError:
                node.length = 0.0
        return node

    root = parse_node()
    if pos[0] < n and text[pos[0]] not in ";":
        # Unconsumed structural text means the parse went wrong, and a partial parse
        # yields a plausible-looking tree with most of its tips missing. Fail loudly.
        raise ValueError(
            f"unparsed text remains at offset {pos[0]} of {n}: "
            f"{text[pos[0]:pos[0] + 60]!r} -- the tree was only partly read"
        )
    return root


def leaves(node):
    if not node.children:
        return [node]
    out = []
    for child in node.children:
        out.extend(leaves(child))
    return out


def _ancestors(node):
    chain, cursor = [], node
    while cursor is not None:
        chain.append(cursor)
        cursor = cursor.parent
    return chain


def patristic(a, b):
    """Summed branch length between two tips, plus their MRCA."""
    up_b = {id(x) for x in _ancestors(b)}
    distance, mrca = 0.0, None
    for node in _ancestors(a):
        if id(node) in up_b:
            mrca = node
            break
        distance += node.length
    cursor = b
    while cursor is not None and (mrca is None or id(cursor) != id(mrca)):
        distance += cursor.length
        cursor = cursor.parent
    return distance, mrca


# ------------------------------------------------------- tips -> our gene IDs
def build_resolver(gene_universe, gene_regex, taxon_regex, gene_prefix):
    """Return tip_label -> gene_id (or None if the tip is not one of our genes).

    Resolution goes through the set of gene IDs that actually appear in the group
    table, so label prefixes ('_PGA_') and transcript suffixes ('.1') do not have
    to be spelled out by the caller.
    """
    gene_re = re.compile(gene_regex)
    taxon_re = re.compile(taxon_regex) if taxon_regex else None
    universe = set(gene_universe)

    # index by the regex key so a tip resolves in one lookup
    by_key = defaultdict(set)
    for gene in universe:
        for m in gene_re.finditer(gene):
            by_key[m.group(0)].add(gene)

    unresolved, ambiguous = set(), {}

    def resolve(label):
        if not label:
            return None
        if taxon_re is not None and not taxon_re.search(label):
            return None
        match = gene_re.search(label)
        if not match:
            return None
        key = match.group(0)
        hits = by_key.get(key)
        if hits and len(hits) == 1:
            return next(iter(hits))
        if hits:
            ambiguous[key] = sorted(hits)
            return None
        candidate = f"{gene_prefix}{key}" if gene_prefix else key
        if candidate in universe:
            return candidate
        # The tip names a gene of ours that carries no shared reads at all, so it
        # never entered the group table. That is informative, not an error.
        unresolved.add(key)
        return candidate if gene_prefix else key

    resolve.unresolved = unresolved
    resolve.ambiguous = ambiguous
    return resolve


# --------------------------------------------------------------- group table
def load_group_table(path):
    """Read a ParalogGroups TSV into a list of dicts."""
    groups = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        needed = {"genes", "n_samples_supporting", "per_member_ambiguity"}
        missing = needed - set(header)
        if missing:
            raise ValueError(
                f"{path}: missing column(s) {sorted(missing)}. Is this a "
                f"ParalogGroups table?"
            )
        idx = {name: i for i, name in enumerate(header)}
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
                "set": frozenset(genes),
                "n_samples": int(fields[idx["n_samples_supporting"]]),
                "ambiguity": ambiguity,
                "unique": unique,
            })
    return groups


def groups_from_rows(rows):
    """Adapt ParalogGroups' in-memory rows to the shape used here, so the same
    calibration can run inline without a round trip through disk."""
    groups = []
    for row in rows:
        genes = row["genes"].split(",")
        ambiguity, unique = {}, {}
        for token in row["per_member_ambiguity"].split(","):
            gene, value = token.rsplit(":", 1)
            ambiguity[gene] = float(value)
        for token in row.get("per_member_unique_reads", "").split(","):
            if ":" in token:
                gene, value = token.rsplit(":", 1)
                unique[gene] = float(value)
        groups.append({
            "genes": genes, "set": frozenset(genes),
            "n_samples": row["n_samples_supporting"],
            "ambiguity": ambiguity, "unique": unique,
        })
    return groups


# ------------------------------------------------------------------ analysis
def collect_pairs(tree_files, resolve):
    """Every within-genome tip pair across all trees, with distance and clade type."""
    pairs = []
    for path in tree_files:
        family = os.path.splitext(os.path.basename(path))[0]
        with open(path, encoding="utf-8", errors="replace") as handle:
            text = strip_rtf(handle.read())
        newick = extract_newick(text)
        if not newick:
            print(f"warning: {path}: no Newick string found; skipped", file=sys.stderr)
            continue
        try:
            tree = parse_newick(newick)
        except Exception as exc:                                # noqa: BLE001
            print(f"warning: {path}: could not parse ({exc}); skipped", file=sys.stderr)
            continue

        tips = leaves(tree)
        ours = [(tip, resolve(tip.name)) for tip in tips]
        ours = [(tip, gene) for tip, gene in ours if gene]
        if len(ours) < 2:
            print(f"note: {family}: {len(ours)} tip(s) from this genome, no pairs to "
                  f"test (tree has {len(tips)} tips)", file=sys.stderr)
            continue

        for i in range(len(ours)):
            for j in range(i + 1, len(ours)):
                tip_a, gene_a = ours[i]
                tip_b, gene_b = ours[j]
                if gene_a == gene_b:
                    continue
                distance, mrca = patristic(tip_a, tip_b)
                clade_tips = leaves(mrca) if mrca is not None else []
                genome_specific = bool(clade_tips) and all(
                    resolve(t.name) for t in clade_tips
                )
                pairs.append({
                    "family": family,
                    "a": gene_a, "b": gene_b,
                    "distance": distance,
                    "genome_specific_clade": genome_specific,
                    "clade_size": len(clade_tips),
                })
    return pairs


def annotate_detection(pairs, groups):
    """Mark each tree pair with whether ParalogGroups put the two genes together."""
    by_gene = defaultdict(list)
    for group in groups:
        for gene in group["genes"]:
            by_gene[gene].append(group)

    for pair in pairs:
        shared = [g for g in by_gene.get(pair["a"], []) if pair["b"] in g["set"]]
        if shared:
            best = max(shared, key=lambda g: g["n_samples"])
            pair["detected"] = True
            pair["ambiguity"] = min(best["ambiguity"][pair["a"]],
                                    best["ambiguity"][pair["b"]])
            pair["n_samples"] = best["n_samples"]
            pair["unique_a"] = best["unique"].get(pair["a"])
            pair["unique_b"] = best["unique"].get(pair["b"])
        else:
            pair.update({"detected": False, "ambiguity": None, "n_samples": None,
                         "unique_a": None, "unique_b": None})
    return pairs


def expected_positive(pair, mode, max_distance):
    """Which tree pairs SHOULD the read-based method have grouped?

    Ancestry alone is too generous -- two genes can share a duplication 100 My ago
    and be trivially separable now. The defensible expectation is a recent
    duplication, identified either as a genome-specific clade or by a distance cut.
    """
    if mode == "distance":
        return pair["distance"] <= max_distance
    return pair["genome_specific_clade"]


def separation_diagnostic(pairs):
    """Is there a clean divergence boundary between detected and undetected pairs?

    If detection is really tracking sequence divergence, the two distance
    distributions should not interleave. Overlap is the signal that something other
    than divergence is driving calls.
    """
    detected = [p["distance"] for p in pairs if p["detected"]]
    missed = [p["distance"] for p in pairs if not p["detected"]]
    if not detected or not missed:
        return None
    return {
        "max_detected": max(detected),
        "min_missed": min(missed),
        "clean": max(detected) < min(missed),
        "n_detected": len(detected),
        "n_missed": len(missed),
        "overlap": [p for p in pairs
                    if (p["detected"] and p["distance"] > min(missed))
                    or (not p["detected"] and p["distance"] < max(detected))],
    }


def sweep_cutoffs(pairs, mode, max_distance, cutoffs, gene_costs=None):
    """Sensitivity against the tree-confirmed set at each ambiguity cutoff."""
    positives = [p for p in pairs if expected_positive(p, mode, max_distance)]
    rows = []
    for cutoff in cutoffs:
        kept = sum(1 for p in positives
                   if p["detected"] and p["ambiguity"] is not None
                   and p["ambiguity"] >= cutoff)
        row = {
            "cutoff": cutoff,
            "tree_pairs_retained": kept,
            "tree_pairs_total": len(positives),
            "sensitivity": (kept / len(positives)) if positives else float("nan"),
        }
        if gene_costs is not None:
            row["genes_merged"] = sum(1 for a in gene_costs if a >= cutoff)
        rows.append(row)
    return rows, positives


def gene_ambiguity_universe(groups):
    """Per gene, the highest ambiguity it reaches in any group -- the quantity a
    genome-wide cutoff is actually applied to."""
    worst = {}
    for group in groups:
        for gene in group["genes"]:
            value = group["ambiguity"][gene]
            if value > worst.get(gene, -1.0):
                worst[gene] = value
    return worst


# -------------------------------------------------------------------- report
def run_calibration(groups, tree_files, args, stream=sys.stderr):
    """Shared entry point: usable from this subcommand or inline from ParalogGroups."""
    universe = {gene for group in groups for gene in group["genes"]}
    resolve = build_resolver(universe, args.tree_gene_regex,
                             args.tree_taxon_regex, args.gene_prefix)

    pairs = collect_pairs(tree_files, resolve)
    if not pairs:
        print("no within-genome tip pairs found in any tree -- nothing to calibrate. "
              "Check --tree-gene-regex against your tip labels.", file=stream)
        return None
    annotate_detection(pairs, groups)

    if resolve.ambiguous:
        for key, hits in sorted(resolve.ambiguous.items())[:5]:
            print(f"warning: tip key {key!r} matches {len(hits)} gene IDs "
                  f"({hits[:3]}...); skipped", file=stream)
    if resolve.unresolved:
        sample = sorted(resolve.unresolved)[:5]
        print(f"note: {len(resolve.unresolved)} tip(s) resolved to genes absent from "
              f"the group table (e.g. {sample}). Those genes share no reads with "
              f"anything, which is itself a result: they are individually estimable.",
              file=stream)

    worst = gene_ambiguity_universe(groups)
    cutoffs = [0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    sweep, positives = sweep_cutoffs(pairs, args.expected_positives,
                                     args.tree_max_distance, cutoffs,
                                     gene_costs=list(worst.values()))
    diag = separation_diagnostic(pairs)

    if args.tree_report:
        with open(args.tree_report, "w") as out:
            columns = ["family", "gene_a", "gene_b", "patristic_distance",
                       "genome_specific_clade", "clade_size", "expected_positive",
                       "detected", "ambiguity", "n_samples_supporting",
                       "unique_reads_a", "unique_reads_b"]
            out.write("\t".join(columns) + "\n")
            for pair in sorted(pairs, key=lambda p: p["distance"]):
                out.write("\t".join([
                    pair["family"], pair["a"], pair["b"],
                    f"{pair['distance']:.6f}",
                    str(pair["genome_specific_clade"]),
                    str(pair["clade_size"]),
                    str(expected_positive(pair, args.expected_positives,
                                          args.tree_max_distance)),
                    str(pair["detected"]),
                    "" if pair["ambiguity"] is None else f"{pair['ambiguity']:.3f}",
                    "" if pair["n_samples"] is None else str(pair["n_samples"]),
                    "" if pair["unique_a"] is None else f"{pair['unique_a']:.0f}",
                    "" if pair["unique_b"] is None else f"{pair['unique_b']:.0f}",
                ]) + "\n")
        print(f"per-pair tree comparison -> {args.tree_report}", file=stream)

    if args.calibration_report:
        with open(args.calibration_report, "w") as out:
            out.write("ambiguity_cutoff\ttree_pairs_retained\ttree_pairs_total\t"
                      "sensitivity\tgenes_merged_genomewide\n")
            for row in sweep:
                out.write(f"{row['cutoff']:.2f}\t{row['tree_pairs_retained']}\t"
                          f"{row['tree_pairs_total']}\t{row['sensitivity']:.3f}\t"
                          f"{row.get('genes_merged', '')}\n")
        print(f"cutoff sweep -> {args.calibration_report}", file=stream)

    # ---- console summary
    print("", file=stream)
    print(f"{len(pairs)} within-genome pair(s) across "
          f"{len({p['family'] for p in pairs})} tree(s); "
          f"{len(positives)} treated as expected positives "
          f"({args.expected_positives} rule)", file=stream)

    if diag:
        if diag["clean"]:
            print(f"detection separates cleanly on divergence: every pair below "
                  f"{diag['max_detected']:.4f} was detected, every pair above "
                  f"{diag['min_missed']:.4f} was not "
                  f"({diag['n_detected']} vs {diag['n_missed']} pairs, no overlap)",
                  file=stream)
        else:
            print(f"detection does NOT separate cleanly on divergence: "
                  f"{len(diag['overlap'])} pair(s) sit on the wrong side of the "
                  f"boundary (detected up to {diag['max_detected']:.4f}, "
                  f"missed from {diag['min_missed']:.4f})", file=stream)

    detected_positives = [p for p in positives if p["detected"]]
    if detected_positives:
        floor = min(p["ambiguity"] for p in detected_positives)
        print(f"lowest ambiguity among detected expected-positives: {floor:.3f} "
              f"-- any cutoff above this drops a tree-confirmed pair", file=stream)
        n_at = sum(1 for a in worst.values() if a >= floor)
        print(f"a cutoff of {floor:.2f} flags {n_at} gene(s) genome-wide", file=stream)

    undetected = [p for p in positives if not p["detected"]]
    if undetected:
        print(f"{len(undetected)} expected-positive pair(s) NOT detected:", file=stream)
        for pair in sorted(undetected, key=lambda p: p["distance"])[:10]:
            print(f"    {pair['family']:<16}{pair['a']} + {pair['b']}  "
                  f"d={pair['distance']:.4f}", file=stream)
        print("  A tree pair that reads CAN separate is a correct non-detection for "
              "merging purposes; check the distance before treating it as a miss.",
              file=stream)

    return {"pairs": pairs, "sweep": sweep, "diagnostic": diag,
            "positives": positives}


def add_tree_arguments(parser, required_group=True):
    """Options shared between this subcommand and ParalogGroups' inline mode."""
    parser.add_argument(
        "--trees", nargs="+", required=required_group, metavar="TREE",
        help="Newick tree file(s). RTF-wrapped trees (saved from a word processor) "
             "are accepted directly.",
    )
    parser.add_argument(
        "--tree-gene-regex", default=r"scaffold[^\s'\").,]*_\d+",
        help="Regex pulling a gene key out of a tip label. The key is then resolved "
             "against the gene IDs present in the group table, so prefixes and "
             "transcript suffixes need not be spelled out "
             r"(default: scaffold[^\s'\").,]*_\d+).",
    )
    parser.add_argument(
        "--tree-taxon-regex", default=None,
        help="Optional extra filter: only tips whose label matches this regex are "
             "considered ours. Usually unnecessary, since a tip that does not "
             "resolve to a known gene is excluded anyway.",
    )
    parser.add_argument(
        "--gene-prefix", default="",
        help="Prefix to prepend to a tip key when it cannot be resolved against the "
             "group table (e.g. '_PGA_').",
    )
    parser.add_argument(
        "--expected-positives", choices=["clade", "distance"], default="clade",
        help="Which tree pairs the detector is expected to find: 'clade' = pairs "
             "whose MRCA contains only genes from this genome (a duplication after "
             "the split from the other taxa in the tree); 'distance' = pairs closer "
             "than --tree-max-distance (default: clade).",
    )
    parser.add_argument(
        "--tree-max-distance", type=float, default=0.1,
        help="Patristic distance defining a recent duplicate when "
             "--expected-positives distance is used (default: 0.1).",
    )
    parser.add_argument(
        "--tree-report", default=None,
        help="Write the per-pair comparison table here.",
    )
    parser.add_argument(
        "--calibration-report", default=None,
        help="Write the ambiguity-cutoff sweep here.",
    )
    return parser


def add_arguments(parser):
    parser.add_argument(
        "--groups", required=True,
        help="A ParalogGroups output table (run it with --min-ambiguity 0 so the "
             "full range is available to calibrate over).",
    )
    add_tree_arguments(parser, required_group=True)
    return parser


def main(args):
    groups = load_group_table(args.groups)
    print(f"{len(groups)} paralog group(s) loaded from {args.groups}", file=sys.stderr)
    if not args.tree_report and not args.calibration_report:
        args.tree_report = "paralog_tree_pairs.tsv"
        args.calibration_report = "paralog_cutoff_calibration.tsv"
        print("no report paths given; writing paralog_tree_pairs.tsv and "
              "paralog_cutoff_calibration.tsv", file=sys.stderr)
    result = run_calibration(groups, args.trees, args)
    return 0 if result else 1


if __name__ == "__main__":
    standalone_parser = argparse.ArgumentParser(
        description="Ground-truth paralog groups against gene trees and calibrate "
                    "the ambiguity cutoff.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    add_arguments(standalone_parser)
    sys.exit(main(standalone_parser.parse_args()))
