# EdgeRDE Example: Arabidopsis Genotypes (Single-Genome DE)

This example demonstrates EdgeRDE for comparing expression between Arabidopsis genotypes.

## Files in This Example

1. `arabidopsis_counts.tsv` - Raw read count matrix
2. `arabidopsis_metadata.tsv` - Sample metadata with group assignments
3. `run_edgerde.sh` - Executable workflow script

## Data Description

**Experiment**: Expression comparison between WT and mutant Arabidopsis accessions in leaf tissue

**Samples**:
- 3 replicates of WT (Columbia ecotype)
- 3 replicates of Mutant (Ler ecotype)
- All grown in controlled conditions, leaf tissue harvested at 6 weeks

**Read type**: Paired-end, Illumina 100bp, ~20 million reads/sample

---

## Step-by-Step Walkthrough

### Step 1: Prepare Input Files

**File: `arabidopsis_counts.tsv`** (first 10 genes shown)

```
gene_id         WT_1        WT_2        WT_3        Mut_1       Mut_2       Mut_3
AT1G01010       150         142         168         145         125         130
AT1G01020       0           5           12          0           0           2
AT1G01030       1234        1456        1290        1512        1680        1545
AT1G01040       0           0           0           0           0           0
AT1G01050       567         489         601         523         612         498
AT1G01060       890         765         823         1245        1150        1189
AT1G01070       123         145         112         98          102         94
AT1G01080       45          52          48          61          58          64
AT1G01090       0           0           0           0           0           0
AT1G01100       2341        2156        2098        1876        1834        1923
```

**Rules**:
- Tab-separated (not spaces)
- First column: Gene IDs
- Remaining columns: Raw counts (integers ≥ 0)
- No normalization, transformation, or filtering
- Column names will be used to match metadata

---

**File: `arabidopsis_metadata.tsv`**

```
sample_name     group       genotype    tissue      replicate
WT_1            WT          Col-0       leaf        1
WT_2            WT          Col-0       leaf        2
WT_3            WT          Col-0       leaf        3
Mut_1           Mutant      Ler         leaf        1
Mut_2           Mutant      Ler         leaf        2
Mut_3           Mutant      Ler         leaf        3
```

**Rules**:
- Tab-separated
- Column `sample_name` must match count matrix column headers exactly
- Column `group` defines experimental groups for contrasts
- Optional columns (genotype, tissue, replicate) are ignored by EdgeRDE

**Alternative format** (using species + tissue):
```
sample_name     species     tissue      replicate
WT_1            Col-0       leaf        1
WT_2            Col-0       leaf        2
WT_3            Col-0       leaf        3
Mut_1           Ler         leaf        1
Mut_2           Ler         leaf        2
Mut_3           Ler         leaf        3
```
EdgeRDE will derive `group` as `{species}_{tissue}` → `Col-0_leaf`, `Ler_leaf`

---

### Step 2: Run EdgeRDE

**File: `run_edgerde.sh`**

```bash
#!/bin/bash
set -e

# Analysis parameters
PROJECT_DIR="/path/to/project"
COUNTS_FILE="${PROJECT_DIR}/arabidopsis_counts.tsv"
METADATA_FILE="${PROJECT_DIR}/arabidopsis_metadata.tsv"
OUTPUT_DIR="${PROJECT_DIR}/results/DE_analysis"

# Run EdgeRDE
python SalmonStreamer.py EdgeRDE \
    --expression-file "${COUNTS_FILE}" \
    --metadata-file "${METADATA_FILE}" \
    --output-dir "${OUTPUT_DIR}" \
    --fdr-threshold 0.05 \
    --logfc-threshold 1.0

echo "EdgeRDE analysis complete. Results in: ${OUTPUT_DIR}"
```

**Run it**:
```bash
bash run_edgerde.sh
```

---

### Step 3: Interpret Results

After EdgeRDE completes, you have:

**Main results file** (`WT_vs_Mutant_DE_results.tsv`):
```
TranscriptID    logFC       F-statistic p-value     FDR         Significance
AT1G01010       -0.12       0.45        0.51        0.78        Not Significant
AT1G01030       0.34        8.2         0.001       0.0034      Significant
AT1G01060       -0.55       15.3        2.1e-05     0.00012     Significant
```

**Significant genes only** (`WT_vs_Mutant_significant_genes.tsv`):
```
TranscriptID    logFC       F-statistic p-value     FDR
AT1G01030       0.34        8.2         0.001       0.0034
AT1G01060       -0.55       15.3        2.1e-05     0.00012
```

**QC plots**:
- `PCA_plot.pdf` - Check that replicates cluster; groups separate
- `sample_correlation_heatmap.pdf` - Verify within-group correlation
- `library_sizes.pdf` - Check sequencing depth consistency
- `WT_vs_Mutant_volcano.pdf` - Visual summary of DE results

**Interpretation**:
- `AT1G01030`: ~1.27-fold higher in WT (logFC=0.34); FDR=0.0034 (highly significant)
- `AT1G01060`: ~1.47-fold lower in WT (logFC=-0.55); FDR=0.00012 (highly significant)
- `AT1G01010`: Not significant (FDR=0.78)

---

## Common Next Steps

1. **Functional annotation**: What do the top DE genes do?
   - Search BY identifier on TAIR (arabidopsis.org)
   - GO enrichment analysis

2. **Validation**: qRT-PCR for top 5-10 genes
   ```bash
   # Design allele-specific primers for AT1G01030, AT1G01060, etc.
   ```

3. **If hybrid cross data available**: Run ASEIntegrate to categorize DE mechanism

4. **Pathway analysis**: 
   - Do DE genes cluster in specific pathways?
   - Example: "All DE genes are photosynthesis" suggests light-response changes

---

## Tips & Troubleshooting

**Issue**: "No samples shared between metadata and count matrix"
- **Fix**: Check that `sample_name` column matches count matrix headers exactly (case-sensitive)
- Verify no trailing spaces in either file

**Issue**: Very few significant genes
- **Check**: Is FDR threshold too strict? (try 0.1 instead of 0.05)
- May indicate low effect size or high biological noise

**Issue**: Groups don't separate on PCA
- **Possible causes**: 
  - Experiment didn't work (biological effect too small)
  - Technical noise dominates (batch effects from sequencing)
- **Fix**: Check metadata — are samples truly from different groups?

---

## Advanced Options

**Strip suffix from sample names** (if they contain patterns like `_R1`, `_filtered`):
```bash
python SalmonStreamer.py EdgeRDE \
    --expression-file arabidopsis_counts.tsv \
    --metadata-file arabidopsis_metadata.tsv \
    --output-dir results/ \
    --sample-suffix "_filtered$"  # Remove "_filtered" before matching metadata
```

**Inline group specification** (if no metadata file):
```bash
python SalmonStreamer.py EdgeRDE \
    --expression-file arabidopsis_counts.tsv \
    --group-samples WT:WT_1,WT_2,WT_3 Mutant:Mut_1,Mut_2,Mut_3 \
    --output-dir results/
```

---

## See Also

- [Input Format Specification](../docs/input_format_specification.md)
- [Output Files Reference](../docs/output_files_reference.md)
- [ASE + DE Integration Workflow](../docs/ase_de_integration_workflow.md) (if hybrid data available)
