# Input Format Specification

## Overview

This document specifies the exact format requirements for input files used by `EdgeRDE` and `ASEIntegrate` subcommands.

---

## Count Matrix / Expression File

### Required Format

- **File type**: Tab-separated values (TSV)
- **First column**: Gene/transcript IDs
- **Remaining columns**: One per sample, raw read counts
- **Row names**: Gene IDs (no duplicates)
- **Column names**: Sample identifiers (must match metadata)

### Example

```
gene_id         sample_1        sample_2        sample_3        sample_4
GENE0001        150             142             168             155
GENE0002        0               5               12              8
GENE0003        1234            1456            1290            1512
GENE0004        0               0               0               0
GENE0005        567             489             601             523
```

### Rules

1. **No header row** for data values — the first row IS the header
2. **Gene IDs must be unique** across rows (no duplicates)
3. **Raw counts only** — do NOT pre-normalize, filter, or transform
4. **All counts must be integers** ≥ 0
5. **Do NOT include normalization factors** (edgeR handles this)
6. **Do NOT include log-transformed data** (edgeR expects raw counts)
7. **Column order**: Can be in any order; metadata determines grouping

### Important Notes

- Sample column names must match `sample_name` values in metadata file exactly
- Use `--sample-suffix` if column names have common patterns:
  ```bash
  --sample-suffix "_R1_filtered$"  # Remove '_R1_filtered' before matching metadata
  ```
- Gene IDs in this file will be matched against the ASE table
  - If the IDs carry a transcript suffix (e.g., `GENE0001.T1`, `GENE0001.T2`), ASEIntegrate removes `.T[digit]` before matching
  - The ASE table should contain gene-level IDs only

---

## Metadata File

### Required Format

- **File type**: Tab-separated values (TSV)
- **Encoding**: UTF-8
- **First row**: Column headers (case-insensitive matching in EdgeRDE)

### Required Columns

Choose ONE of these options:

**Option A** (Recommended):
- `sample_name`: Must match count matrix column headers exactly (after `--sample-suffix` stripping)
- `group`: Experimental group label (e.g., "Hybrid_Leaf", "Parent1_Root")
- *(optional)* a column naming the **biological unit** each library came from (for example `individual`), passed to `EdgeRDE --collapse-replicates-by`. Libraries sharing a value are summed into one column before the analysis. Needed whenever several libraries come from one individual; each value must belong to a single group

**Option B**:
- `sample_name`: Must match count matrix column headers exactly
- `species`: Species or line name (e.g., "LineA", "Parent1")
- `tissue`: Tissue type (e.g., "leaf", "root")
- EdgeRDE will derive `group` as `{species}_{tissue}`

### Optional Columns

You can include additional columns; they will be ignored by EdgeRDE:
- `platform`: Sequencing platform (e.g., "Illumina", "PacBio")
- `library`: Library prep method
- `flowcell`: Flowcell ID
- `lane`: Sequencing lane
- `replicate`: Replicate number
- `treatment`: Experimental treatment
- `timepoint`: Time point in time series
- Any other metadata

### Example (Option A: Direct `group`)

```
sample_name     group           platform        lane
sample_1        F1_hybrid       Illumina        lane1
sample_2        F1_hybrid       Illumina        lane1
sample_3        F1_hybrid       Illumina        lane2
sample_4        Parent1         Illumina        lane3
sample_5        Parent1         Illumina        lane3
sample_6        Parent2         Illumina        lane4
```

### Example (Option B: Derive `group` from `species` + `tissue`)

```
sample_name     species     tissue      replicate
sample_1        F1          leaf        1
sample_2        F1          leaf        2
sample_3        F1          leaf        3
sample_4        Parent1     leaf        1
sample_5        Parent1     leaf        2
sample_6        Parent2     leaf        1
```

**Result**: EdgeRDE derives groups as `F1_leaf`, `Parent1_leaf`, `Parent2_leaf` and runs all pairwise contrasts.

**For ASEIntegrate's advanced mode** the metadata file must have an explicit `group` column (Option A). ASEIntegrate does not derive groups from `species` + `tissue`, and it has no `--sample-suffix` option, so `sample_name` must equal the count-matrix column headers as they are. If EdgeRDE was run with Option B, add a `group` column holding the derived names (`Parent1_leaf`, ...) before running the advanced mode.

### Rules

1. **Sample names are case-sensitive** — must match count matrix column headers exactly
2. **Group names can be anything** — descriptive labels (no special characters recommended)
3. **At least 2 samples per group** recommended for statistical testing
4. **At least 2 groups** required to run contrasts
5. **No duplicate sample names** allowed
6. **Column names must be present** — will match case-insensitively

---

## ASE Counts File

`ASEIntegrate` does not compute allele-specific expression. `--ase-counts-file` takes a finished per-gene ASE table made beforehand with a separate tool. Two formats are accepted and told apart by their columns; lines starting with `#` are skipped in both.

In both formats `gene_id` is used as it is. It is matched against the DE table's `TranscriptID`, and in advanced mode against the first column of the expression matrix, after any `.T[digit]` suffix has been removed from those. The ASE table should therefore hold gene-level IDs without that suffix, and the match is exact and case-sensitive.

### Format A: ASErtain Gene Table (Preferred)

An [ASErtain](https://github.com/luismadrigal98/ASErtain) `<prefix>.gene_ase.tsv` can be passed to `--ase-counts-file` unchanged. It is recognised by its `gene_id`, `n_snps`, `direction` and `ase_call` columns plus a `mean_<label>_ratio` column, where `<label>` is one of the two values given as `--parent1-label` / `--parent2-label` (they must be the lineage labels of the ASErtain config). A gene is then treated as allelically biased when `ase_call` is true, which is a replicate-aware test.

| ASEIntegrate needs | Taken from the ASErtain table |
|--------------------|-------------------------------|
| Is the gene allelically biased? | `ase_call` |
| Toward which parent? | `direction`: the bias is toward parent 1 when the value equals `--parent1-label` |
| Parent 1 allele fraction | `mean_<label>_ratio`. ASErtain writes this column for one of the two lineages; when it is parent 2's, the parent 1 fraction is one minus it |
| SNP count | `n_snps` |

`q_value`, `log2_ratio` and `phase_discordance_cause` are copied to the output with an `ase_` prefix when present. No other column of the table is used.

If neither `mean_<parent1_label>_ratio` nor `mean_<parent2_label>_ratio` exists, the run stops and lists the columns it found. When the ASErtain run used its default labels, the column is `mean_variable_ratio` and the labels to pass are `variable` and `fixed`.

### Format B: Summary Table (Legacy)

The older six-column summary table is still accepted. It is the format written by `ase_read_counter.py` of the EDGE_Penstemon pipeline. Here a gene counts as allelically biased when `bias_ratio > 0.5`. That is a majority vote over SNPs, not a statistical test, and it takes no account of read depth or replication.

- **File type**: Tab-separated values (TSV)
- **First row**: Column headers
- **Numeric columns**: All ratios and counts must be numeric

#### Required Columns

1. **`gene_id`**: Gene identifier
   - Gene-level ID, matched as described above
   - One row per gene (duplicates are not checked for)

2. **`snp_count`**: Number of diagnostic SNPs in this gene
   - Integer ≥ 0
   - Copied to the output; the classification does not use it

3. **`bias_ratio`**: Fraction of diagnostic SNPs showing allelic bias
   - Numeric, range [0.0, 1.0], not a percentage. The first rows are checked before the run and a value above 1 stops it
   - The gene is treated as biased when the value is above 0.5

4. **`predominant_bias`**: Which parent's allele dominates
   - The value given as `--parent1-label` or `--parent2-label` (defaults: `parent1`, `parent2`)
   - Case-sensitive. Only equality with the parent 1 label is tested, so any other value is read as parent 2

5. **`{parent1_label}_ratio`**: Proportion of reads from parent 1 allele
   - Numeric, range [0.0, 1.0]
   - Default column name: `parent1_ratio`
   - With custom labels, for example `--parent1-label lineA`: `lineA_ratio`
   - This is the allele fraction used for the plot, `cis_log2` and `Ad`

6. **`{parent2_label}_ratio`**: Proportion of reads from parent 2 allele
   - Numeric, range [0.0, 1.0]
   - Default column name: `parent2_ratio`
   - Must be present; it is copied to the output and not used in the calculations

Any further columns are copied to the output unchanged.

#### Example (Default Labels)

```
gene_id     snp_count   bias_ratio   predominant_bias   parent1_ratio   parent2_ratio
GENE0001    8           0.75         parent1            0.78            0.22
GENE0002    4           0.50         parent1            0.51            0.49
GENE0003    12          0.92         parent2            0.08            0.92
GENE0004    2           0.50         parent1            0.45            0.55
GENE0005    15          0.87         parent1            0.85            0.15
```

#### Example (Custom Labels)

When running ASEIntegrate with `--parent1-label lineA --parent2-label lineB`:

```
gene_id     snp_count   bias_ratio   predominant_bias   lineA_ratio   lineB_ratio
GENE0001    8           0.75         lineA              0.78          0.22
GENE0002    4           0.50         lineA              0.51          0.49
GENE0003    12          0.92         lineB              0.08          0.92
```

#### Rules

1. **Parent labels must match** the ones provided to ASEIntegrate:
   - `--parent1-label parent1` expects column `parent1_ratio`
   - `--parent1-label lineA` expects column `lineA_ratio`

2. **Ratio columns should sum to about 1.0**
   - This is not checked. If they do not, look for other alleles or unassigned reads in the tool that made the table

3. **Genes without informative SNPs are still classified**
   - `snp_count` is not looked at. A gene listed with `snp_count = 0` and `bias_ratio = 0` counts as unbiased, so it comes out `trans_only` when it is DE between the parents and `no_divergence` otherwise; in advanced mode a ratio of 0.5 reads as balanced alleles
   - Leave such genes out of the table unless that is what you want

4. **Gene IDs must match the DE table**
   - Exact match (case-sensitive)
   - If the count matrix used `GENE0001.T1` and the ASE table has `GENE0001`, they are matched. The suffix is removed from the DE and expression IDs, never from the ASE table

5. **One row per gene**

#### Rules of Thumb for Data Quality

These are aids for reading the table. ASEIntegrate applies none of them.

- `snp_count ≥ 4` gives a more robust call than one or two SNPs
- `bias_ratio` is a fraction of SNPs, not an allele ratio: 0.5 means half of the SNPs show bias. The allele split itself is in the `{label}_ratio` columns

---

## File Naming Conventions

### Recommended

- Expression files: `{study}_counts.tsv` or `{cross}_raw_counts.tsv`
- Metadata files: `{study}_metadata.tsv` or `{cross}_sample_info.tsv`
- ASE files: `{study}_ase_counts.tsv` or `{cross}_allele_bias.tsv`
- Output directories: `{analysis}_DE_results/`, `{analysis}_ASE_DE/`

### Examples

```
# Named after the study
study1_counts.tsv
study1_metadata.tsv
study1_ase_counts.tsv

# Named after the cross
lineA_x_lineB_raw_counts.tsv
lineA_x_lineB_sample_info.tsv
lineA_x_lineB_allele_bias.tsv
```

---

## Validation Checklist

Before running EdgeRDE / ASEIntegrate, verify:

### Count Matrix
- [ ] Tab-separated (not spaces or commas)
- [ ] First column = gene IDs
- [ ] All numeric columns contain integers ≥ 0
- [ ] No duplicate gene IDs
- [ ] Column names match metadata `sample_name` values

### Metadata
- [ ] Tab-separated (not spaces or commas)
- [ ] Contains `sample_name` column
- [ ] Contains `group` OR both `species` AND `tissue`
- [ ] All `sample_name` values appear in count matrix column headers
- [ ] No duplicate sample names
- [ ] At least 2 samples per group
- [ ] At least 2 groups for contrasts

### ASE Counts File
- [ ] Tab-separated (not spaces or commas)
- [ ] Gene IDs are gene-level and match the count matrix IDs once any `.T[digit]` suffix is removed from the latter
- [ ] No duplicate gene IDs
- [ ] ASErtain gene table: `--parent1-label` / `--parent2-label` are the lineage labels of the ASErtain config
- [ ] Summary table: contains all six required columns
- [ ] Summary table: `snp_count` values are non-negative integers
- [ ] Summary table: `bias_ratio` values in range [0.0, 1.0]
- [ ] Summary table: `predominant_bias` values are the two parent labels (`parent1` / `parent2` by default)
- [ ] Summary table: ratio columns sum to approximately 1.0

---

## Common Formatting Errors

### Error: "No samples shared between metadata and count matrix"

**Cause**: Sample names don't match exactly (spaces, capitalization, suffixes, etc.)

**Fix**:
1. Compare count matrix headers with metadata `sample_name` values character-by-character
2. Use `--sample-suffix` if there are common suffixes to remove
3. Check for invisible characters (tabs, spaces at end of line)

**Example**:
```
Count matrix header: "sample_1"
Metadata sample_name: "sample_1 "  (trailing space!)
Result: NO MATCH

Fix: Remove trailing space in metadata
```

### Error: "Column 'parent1_ratio' not found in ASE file"

**Cause**: Custom parent labels used with a summary table, but column names don't match

**Fix**: Ensure column names match the labels you provide:
```bash
--parent1-label lineA --parent2-label lineB
# Requires columns: lineA_ratio, lineB_ratio

--parent1-label sp1 --parent2-label sp2
# Requires columns: sp1_ratio, sp2_ratio
```

### Error: "Cannot parse count as integer"

**Cause**: Expression file contains non-numeric values

**Fix**:
1. Remove any text or transformed values
2. Use raw read counts only
3. Check for accidental text in numeric columns (quotes, commas in numbers)

---

## Tools for Format Validation

### Quick Python check:

```python
import pandas as pd

# Check expression file
counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)
print(f"Shape: {counts.shape}")  # Should be (n_genes, n_samples)
print(f"All numeric: {counts.dtypes.apply(lambda x: x.kind == 'i' or x.kind == 'u').all()}")

# Check metadata
metadata = pd.read_csv('metadata.tsv', sep='\t')
print(f"Columns: {metadata.columns.tolist()}")
print(f"Samples in metadata: {len(metadata)}")
print(f"Samples in counts: {counts.shape[1]}")
print(f"Matching samples: {set(metadata['sample_name']) & set(counts.columns)}")

# Check ASE file (summary-table format, default labels)
ase = pd.read_csv('ase.tsv', sep='\t', index_col=0)
print(f"Genes in ASE: {len(ase)}")
print(f"Required columns present: {all(col in ase.columns for col in ['snp_count', 'bias_ratio', 'predominant_bias', 'parent1_ratio', 'parent2_ratio'])}")
```

---

## See Also

- [ASE + DE Integration Workflow](./ase_de_integration_workflow.md)
- [Output Files Reference](./output_files_reference.md)
- [Interpretation Guidelines](./interpretation_guidelines.md)
