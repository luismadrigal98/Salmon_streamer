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
AT1G01010       150             142             168             155
AT1G01020       0               5               12              8
AT1G01030       1234            1456            1290            1512
AT1G01040       0               0               0               0
AT1G01050       567             489             601             523
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
- Gene IDs in this file will be matched against ASE counts file
  - If transcript variants exist (e.g., `AT1G01010.T1`, `AT1G01010.T2`), the ASE script will automatically strip `.T[digit]` suffixes
  - ASE file should contain gene-level IDs only

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

**Option B**:
- `sample_name`: Must match count matrix column headers exactly
- `species`: Species or line name (e.g., "Arabidopsis", "Parent1")
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

### Rules

1. **Sample names are case-sensitive** — must match count matrix column headers exactly
2. **Group names can be anything** — descriptive labels (no special characters recommended)
3. **At least 2 samples per group** recommended for statistical testing
4. **At least 2 groups** required to run contrasts
5. **No duplicate sample names** allowed
6. **Column names must be present** — will match case-insensitively

---

## ASE Counts File

### Required Format

- **File type**: Tab-separated values (TSV)
- **First row**: Column headers
- **First column**: Gene IDs
- **Numeric columns**: All ratios and counts must be numeric

### Required Columns

1. **`gene_id`**: Gene identifier
   - Should match gene IDs in count matrix (with or without `.T[digit]` suffix)
   - ASE script automatically strips `.T1`, `.T2` suffixes
   - Must be unique (one row per gene)

2. **`snp_count`**: Number of diagnostic SNPs in this gene
   - Integer ≥ 0
   - Higher values = more robust bias estimates
   - Use for filtering: genes with ≥ 4 SNPs are more reliable

3. **`bias_ratio`**: Fraction of diagnostic SNPs showing allelic bias
   - Numeric, range [0.0, 1.0] (or percentage [0, 100])
   - Indicates consistency of allele bias across SNPs
   - High values (> 0.75) indicate strong, consistent bias

4. **`predominant_bias`**: Which parent's allele dominates
   - Must be exactly one of: `parent1` or `parent2` (or your custom labels if using `--parent1-label` / `--parent2-label`)
   - Case-sensitive
   - Indicates which parent's allele has more reads

5. **`{parent1_label}_ratio`**: Proportion of reads from parent 1 allele
   - Numeric, range [0.0, 1.0]
   - Default column name: `parent1_ratio`
   - If using custom labels: `bellpepper_ratio`, `chili_ratio`, etc.

6. **`{parent2_label}_ratio`**: Proportion of reads from parent 2 allele
   - Numeric, range [0.0, 1.0]
   - Default column name: `parent2_ratio`
   - Must sum (approximately) with corresponding `{parent1_label}_ratio` to 1.0

### Example (Standard Format)

```
gene_id     snp_count   bias_ratio   predominant_bias   parent1_ratio   parent2_ratio
AT1G01010   8           0.75         parent1            0.78            0.22
AT1G01020   4           0.50         parent1            0.51            0.49
AT1G01030   12          0.92         parent2            0.08            0.92
AT1G01040   2           0.50         parent1            0.45            0.55
AT1G01050   15          0.87         parent1            0.85            0.15
AT1G01060   0           0.00         parent1            0.50            0.50
```

### Example (Custom Labels)

When running ASEIntegrate with `--parent1-label BellPepper --parent2-label ChiliPepper`:

```
gene_id         snp_count   bias_ratio   predominant_bias   BellPepper_ratio   ChiliPepper_ratio
Caper_gene_001  8           0.75         BellPepper        0.78               0.22
Caper_gene_002  4           0.50         BellPepper        0.51               0.49
Caper_gene_003  12          0.92         ChiliPepper       0.08               0.92
```

### Rules

1. **Parent labels must match** the ones provided to ASEIntegrate:
   - `--parent1-label parent1` expects column `parent1_ratio`
   - `--parent1-label BellPepper` expects column `BellPepper_ratio`

2. **Ratio columns must sum ≈ 1.0**
   - Small rounding errors (±0.01) are OK
   - If they don't sum to ~1.0, check if there are other alleles or unmapped reads

3. **Zero SNPs are allowed**
   - Genes with `snp_count = 0` will not show allelic bias
   - Included in "no_divergence" category after integration

4. **Gene IDs must match count matrix**
   - Exact match (case-sensitive)
   - Transcript suffixes (`.T1`, `.T2`) are handled automatically
   - If genes in count matrix are `AT1G01010.T1` but ASE file has `AT1G01010`, they will be matched

5. **No duplicate gene IDs** allowed

### Data Quality Considerations

- **Reliable bias**: `snp_count ≥ 4` and `bias_ratio > 0.75`
- **Uncertain bias**: `snp_count < 4` or `bias_ratio < 0.5`
- **Balanced alleles**: `bias_ratio ≈ 0.5` indicates near-50:50 split

---

## File Naming Conventions

### Recommended

- Expression files: `{species}_counts.tsv` or `{cross}_raw_counts.tsv`
- Metadata files: `{species}_metadata.tsv` or `{cross}_sample_info.tsv`
- ASE files: `{species}_ase_counts.tsv` or `{cross}_allele_bias.tsv`
- Output directories: `{analysis}_DE_results/`, `{analysis}_ASE_DE/`

### Examples

```
# Arabidopsis × related species hybrid
arabidopsis_counts.tsv
arabidopsis_metadata.tsv
arabidopsis_ase_counts.tsv

# Bell Pepper × Chili Pepper cross
pepper_cross_counts.tsv
pepper_cross_metadata.tsv
pepper_cross_ase.tsv
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
- [ ] Contains all required columns
- [ ] Gene IDs approximately match count matrix IDs
- [ ] `snp_count` values are non-negative integers
- [ ] `bias_ratio` values in range [0.0, 1.0]
- [ ] `predominant_bias` values are `parent1` or `parent2`
- [ ] Ratio columns sum to approximately 1.0
- [ ] No duplicate gene IDs

---

## Creating ASE Counts File from EDGE_Penstemon

If using the EDGE_Penstemon pipeline, the `ase_read_counter.py` script produces output in the required format:

```bash
# From EDGE_Penstemon
python Python_scripts/ase_read_counter.py \
    --bam-file sample.sorted.bam \
    --vcf-file variants.vcf \
    --output ase_counts.tsv
```

Output will have the required columns:
```
gene_id     snp_count   bias_ratio   predominant_bias   parent1_ratio   parent2_ratio
```

Simply provide this file directly to ASEIntegrate via `--ase-counts-file`.

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

**Cause**: Custom parent labels used, but column names don't match

**Fix**: Ensure column names match the labels you provide:
```bash
--parent1-label BellPepper --parent2-label ChiliPepper
# Requires columns: BellPepper_ratio, ChiliPepper_ratio

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

# Check ASE file
ase = pd.read_csv('ase.tsv', sep='\t', index_col=0)
print(f"Genes in ASE: {len(ase)}")
print(f"Required columns present: {all(col in ase.columns for col in ['snp_count', 'bias_ratio', 'predominant_bias', 'parent1_ratio', 'parent2_ratio'])}")
```

---

## See Also

- [ASE + DE Integration Workflow](./ase_de_integration_workflow.md)
- [Output Files Reference](./output_files_reference.md)
- [Interpretation Guidelines](./interpretation_guidelines.md)
