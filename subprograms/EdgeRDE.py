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

Replication (--collapse-replicates-by):
  Every column of the count matrix is treated as an independent replicate of
  its group. That is only true when each column is a different biological
  unit. Several libraries from ONE individual (several tissue samples of one
  plant, re-sequenced libraries, lanes) are not independent, and leaving them
  as separate columns makes every dispersion too small and every p-value too
  optimistic. Name the metadata column that identifies the biological unit
  with --collapse-replicates-by and the libraries sharing a value are summed
  into one column before the analysis (what edgeR's sumTechReps does).

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
import shutil
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


def _find_column(header, *names):
    """Index of the first header entry matching one of `names`, case-insensitively."""
    lower = [h.strip().lower() for h in header]
    for name in names:
        if name.lower() in lower:
            return lower.index(name.lower())
    return None


def collapse_replicates(expression_file, metadata_file, unit_column, output_dir,
                        sample_suffix=None):
    """Sum the count columns of libraries that belong to one biological unit.

    Returns (collapsed_counts_path, collapsed_metadata_path). The collapsed
    metadata has one row per unit (sample_name = the unit's label, its group,
    and the libraries that went into it) so the provenance stays on disk.

    Raises ValueError when the unit column is missing, when a unit spans more
    than one group (it would have no well-defined group), or when no library
    in the metadata is found in the matrix.
    """
    with open(metadata_file, newline="") as handle:
        rows = [r for r in csv.reader(handle, delimiter="\t") if r and any(c.strip() for c in r)]
    if len(rows) < 2:
        raise ValueError("metadata file has no sample rows")
    header, body = rows[0], rows[1:]
    sample_i = _find_column(header, "sample_name", "sample")
    unit_i = _find_column(header, unit_column)
    group_i = _find_column(header, "group")
    species_i, tissue_i = _find_column(header, "species"), _find_column(header, "tissue")
    if sample_i is None:
        raise ValueError("metadata must contain a 'sample_name' or 'Sample' column")
    if unit_i is None:
        raise ValueError(
            f"--collapse-replicates-by '{unit_column}' is not a metadata column "
            f"(found: {', '.join(header)})")
    if group_i is None and (species_i is None or tissue_i is None):
        raise ValueError("metadata must contain 'group', or 'species' and 'tissue'")

    def group_of(row):
        if group_i is not None:
            return row[group_i].strip()
        return f"{row[species_i].strip()}_{row[tissue_i].strip()}"

    with open(expression_file) as handle:
        matrix_header = handle.readline().rstrip("\n").split("\t")
        matrix_rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    strip = (lambda name: re.sub(sample_suffix, "", name)) if sample_suffix else (lambda name: name)
    column_of = {strip(name): i for i, name in enumerate(matrix_header) if i > 0}

    units, unit_group, missing = {}, {}, []
    for row in body:
        sample, unit = row[sample_i].strip(), row[unit_i].strip()
        if not unit:
            raise ValueError(f"sample '{sample}' has an empty '{unit_column}' value")
        group = group_of(row)
        if unit_group.setdefault(unit, group) != group:
            raise ValueError(
                f"'{unit_column}' value '{unit}' appears in groups "
                f"'{unit_group[unit]}' and '{group}'. A biological unit must belong "
                f"to one group; if one individual was sampled for several groups "
                f"(e.g. tissues), use a column that is unique per individual and group.")
        if sample not in column_of:
            missing.append(sample)
            continue
        units.setdefault(unit, []).append(sample)
    if missing:
        print(f"WARNING: {len(missing)} metadata sample(s) not in the count matrix and "
              f"left out of the collapse: {', '.join(sorted(missing)[:5])}"
              f"{' ...' if len(missing) > 5 else ''}", file=sys.stderr)
    if not units:
        raise ValueError("no metadata sample matches a count-matrix column "
                         "(check sample_name values and --sample-suffix)")

    os.makedirs(output_dir, exist_ok=True)
    counts_path = os.path.join(output_dir, "collapsed_counts.tsv")
    meta_path = os.path.join(output_dir, "collapsed_metadata.tsv")
    unit_names = sorted(units)
    with open(counts_path, "w") as out:
        out.write("\t".join([matrix_header[0]] + unit_names) + "\n")
        for fields in matrix_rows:
            sums = []
            for unit in unit_names:
                total = sum(float(fields[column_of[s]]) for s in units[unit])
                sums.append(str(int(total)) if total == int(total) else repr(total))
            out.write("\t".join([fields[0]] + sums) + "\n")
    with open(meta_path, "w") as out:
        out.write("sample_name\tgroup\tn_libraries\tlibraries\n")
        for unit in unit_names:
            out.write(f"{unit}\t{unit_group[unit]}\t{len(units[unit])}\t"
                      f"{','.join(units[unit])}\n")

    per_group = {}
    for unit in unit_names:
        per_group[unit_group[unit]] = per_group.get(unit_group[unit], 0) + 1
    n_libraries = sum(len(v) for v in units.values())
    print(f"Collapsed {n_libraries} libraries into {len(unit_names)} biological "
          f"unit(s) by '{unit_column}': "
          + ", ".join(f"{g}={n}" for g, n in sorted(per_group.items())))
    single = sorted(g for g, n in per_group.items() if n < 2)
    if single:
        print(f"WARNING: group(s) with a single biological unit after collapsing: "
              f"{', '.join(single)}. They contribute nothing to the dispersion "
              f"estimate; contrasts involving them rest on the variability of the "
              f"other groups.", file=sys.stderr)
    return counts_path, meta_path


def main(args):
    r_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "src", "R_src", "edger_de.R"
    )

    if not os.path.isfile(r_script):
        print(f"ERROR: R script not found at {r_script}", file=sys.stderr)
        sys.exit(1)

    # Validate Rscript executable ----------------------------------------
    rscript = shutil.which(os.path.expanduser(args.rscript_executable))
    if rscript is None:
        print(
            f"ERROR: Rscript executable not found: {args.rscript_executable}\n"
            f"       Put Rscript on your PATH or pass --rscript-executable",
            file=sys.stderr
        )
        sys.exit(1)
    args.rscript_executable = rscript

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

    if getattr(args, "collapse_replicates_by", None) and not args.metadata_file:
        print("ERROR: --collapse-replicates-by needs --metadata-file (the column "
              "naming each library's biological unit lives there).", file=sys.stderr)
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

    # Collapse libraries of one biological unit, if asked ------------------
    expression_path = args.expression_file
    collapse_by = getattr(args, "collapse_replicates_by", None)
    if collapse_by:
        try:
            expression_path, metadata_path = collapse_replicates(
                args.expression_file, metadata_path, collapse_by,
                args.output_dir, sample_suffix=args.sample_suffix)
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        args.sample_suffix = None       # already applied while collapsing

    # Build R command ------------------------------------------------------
    sample_suffix = args.sample_suffix if args.sample_suffix else "NULL"
    export_normalized = getattr(args, "export_normalized_expression", False)
    normalized_format = getattr(args, "normalized_expression_format", "tsv")

    cmd = [
        args.rscript_executable,
        r_script,
        expression_path,
        args.output_dir,
        metadata_path,
        str(args.fdr_threshold),
        str(args.logfc_threshold),
        sample_suffix,
        "TRUE" if export_normalized else "FALSE",
        normalized_format,
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
    if export_normalized:
        print(f"Export normalized expression: yes ({normalized_format})")
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
        "--collapse-replicates-by", default=None, metavar="COLUMN",
        help=(
            "Metadata column naming the biological unit each library came from (e.g. 'individual'). Libraries sharing a value are summed into one column before the analysis, as edgeR's sumTechReps does. Use it whenever several libraries come from one individual: left as separate columns they are counted as independent replicates and the p-values are too small. Requires --metadata-file."
        )
    )
    parser.add_argument(
        "--export-normalized-expression", action="store_true",
        help=(
            "Export the TMM-normalized expression matrices (CPM and log2-CPM) "
            "plus per-sample normalization factors to the output directory "
            "for manual inspection."
        )
    )
    parser.add_argument(
        "--normalized-expression-format", choices=["tsv", "csv"], default="tsv",
        help=(
            "Output format for the TMM-normalized expression matrices "
            "when --export-normalized-expression is set (default: tsv)."
        )
    )
    parser.add_argument(
        "--rscript-executable",
        default="Rscript",
        help="Path to the Rscript executable (default: Rscript on your PATH)."
    )

    main(parser.parse_args())
