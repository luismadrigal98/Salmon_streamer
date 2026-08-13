#!/usr/bin/env python3

"""
Derive a transcript-to-gene map from a GFF3 annotation.

ParalogGroups and ParalogMerge need to know which transcripts belong to the same gene.
When the transcriptome is built here, `ExtractTranscriptome --txp2gene` emits that map
as a byproduct and the names are guaranteed to match. This subcommand covers the other
case: a transcriptome that already exists, indexed before the map was needed or built
by something else entirely.

The failure mode this guards against is silent. If the map's transcript names do not
match the names Salmon indexed, nothing errors -- the unmatched transcripts simply fall
through as their own genes, and the paralog groups come out empty or wrong for reasons
that look like biology. Pass `--transcriptome` (the FASTA that was indexed) and the
overlap is checked and reported before anything downstream consumes the map.

Salmon takes a transcript's name to be its FASTA header up to the first whitespace, so
that prefix is what is compared, not the full description line.
"""

import argparse
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from subprograms.ParalogGroups import _open_maybe_gzip          # noqa: E402


def parse_attributes(attribute_string):
    """Parse a GFF3 column 9 into a dict, tolerating trailing semicolons."""
    attributes = {}
    for pair in attribute_string.split(';'):
        pair = pair.strip()
        if not pair or '=' not in pair:
            continue
        key, value = pair.split('=', 1)
        attributes[key.strip()] = value.strip()
    return attributes


def read_txp2gene(gff_path, transcript_types):
    """Return (rows, multi_parent_count) of (transcript_id, gene_id) from a GFF3."""
    rows = []
    multi_parent = 0
    with _open_maybe_gzip(gff_path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] not in transcript_types:
                continue
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get('ID')
            parent = attributes.get('Parent')
            if not transcript_id:
                continue
            if not parent:
                # A transcript with no parent gene is its own gene, which is also what
                # happens to any name missing from the map downstream.
                rows.append((transcript_id, transcript_id))
                continue
            if ',' in parent:
                # GFF3 permits several parents; a transcript belonging to two genes has
                # no single answer, so take the first and say so rather than guessing
                # quietly.
                multi_parent += 1
                parent = parent.split(',')[0]
            rows.append((transcript_id, parent))
    return rows, multi_parent


def read_fasta_names(fasta_path):
    """Return the set of sequence names in a FASTA, truncated as Salmon truncates them."""
    names = set()
    with _open_maybe_gzip(fasta_path) as handle:
        for line in handle:
            if line.startswith('>'):
                names.add(line[1:].strip().split()[0])
    return names


def add_arguments(parser):
    """Register this subcommand's options (shared by the standalone and CLI entry points)."""
    parser.add_argument(
        "--gff", required=True,
        help="GFF3 annotation the transcriptome was built from.",
    )
    parser.add_argument(
        "-o", "--output", default="txp2gene.tsv",
        help="Output two-column TSV, no header (default: txp2gene.tsv).",
    )
    parser.add_argument(
        "--transcript-types", default="mRNA",
        help="Comma-separated GFF3 feature types to treat as transcripts "
             "(default: mRNA). Use 'transcript' for annotations that use that term.",
    )
    parser.add_argument(
        "--transcriptome", default=None,
        help="The transcriptome FASTA that was indexed. Strongly recommended: the map "
             "is checked against it and a name mismatch is reported here rather than "
             "silently producing empty paralog groups later.",
    )
    parser.add_argument(
        "--require-complete", action="store_true",
        help="Exit non-zero unless every sequence in --transcriptome is covered by the "
             "map. Use in scripted pipelines where a partial map is not acceptable.",
    )
    return parser


def main(args):
    transcript_types = [t.strip() for t in args.transcript_types.split(',') if t.strip()]
    rows, multi_parent = read_txp2gene(args.gff, transcript_types)

    if not rows:
        print(
            f"error: no features of type {'/'.join(transcript_types)} with an ID= "
            f"attribute found in {args.gff}. Check --transcript-types; annotations "
            f"vary between 'mRNA' and 'transcript'.",
            file=sys.stderr,
        )
        return 1

    exit_code = 0
    if args.transcriptome:
        fasta_names = read_fasta_names(args.transcriptome)
        mapped = {transcript for transcript, _ in rows}
        covered = fasta_names & mapped
        uncovered = fasta_names - mapped

        print(
            f"{len(fasta_names)} sequence(s) in {args.transcriptome}; "
            f"{len(covered)} covered by the map, {len(uncovered)} not.",
            file=sys.stderr,
        )
        if uncovered:
            preview = ", ".join(sorted(uncovered)[:5])
            suffix = ", ..." if len(uncovered) > 5 else ""
            print(
                f"warning: {len(uncovered)} transcriptome sequence(s) have no gene "
                f"mapping and would be treated as their own genes ({preview}{suffix}).",
                file=sys.stderr,
            )
            if not covered:
                print(
                    "error: NOTHING matched. The GFF3 IDs and the FASTA headers use "
                    "different conventions, so this map is unusable as written -- "
                    "check for added prefixes or version suffixes on one side.",
                    file=sys.stderr,
                )
                exit_code = 1
            elif args.require_complete:
                exit_code = 1

    # Written even on a coverage failure: seeing the map is what makes the mismatch
    # diagnosable. The non-zero exit is what stops a pipeline from consuming it.
    with open(args.output, 'w') as out:
        for transcript_id, gene_id in rows:
            out.write(f"{transcript_id}\t{gene_id}\n")

    n_genes = len({gene for _, gene in rows})
    print(
        f"{len(rows)} transcript(s) -> {n_genes} gene(s) -> {args.output}",
        file=sys.stderr,
    )
    if multi_parent:
        print(
            f"warning: {multi_parent} transcript(s) listed more than one Parent; "
            f"the first was used.",
            file=sys.stderr,
        )
    return exit_code


if __name__ == "__main__":
    standalone_parser = argparse.ArgumentParser(
        description="Derive a transcript-to-gene map from a GFF3 annotation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    add_arguments(standalone_parser)
    sys.exit(main(standalone_parser.parse_args()))
