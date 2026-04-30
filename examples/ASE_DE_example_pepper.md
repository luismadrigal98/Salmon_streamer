# ASE + DE Integration Example: Pepper Hybrid Cross

This example demonstrates the full ASE + DE integration workflow using pepper hybrid cross data.

## Experiment Overview

**Species**: Bell Pepper (*Capsicum annuum* var. grossum) × Chili Pepper (*Capsicum annuum* var. acuminatum)

**Design**:
- Parent 1: Bell Pepper (large fruited, mild)
- Parent 2: Chili Pepper (small fruited, hot)
- F1 Hybrid: Cross between the two

**Samples**:
- 3 replicates of Bell Pepper parent
- 3 replicates of Chili Pepper parent
- 3 replicates of F1 hybrid
- All fruit tissue at anthesis

**Total**: 9 samples, leaf tissue

---

## Files in This Example

1. `pepper_counts.tsv` - Raw counts matrix
2. `pepper_metadata.tsv` - Sample metadata
3. `pepper_ase.tsv` - Allele-specific expression counts
4. `run_analysis.sh` - Complete workflow script

---

## File Contents

### 1. Count Matrix: `pepper_counts.tsv`

```
gene_id         BellP_1     BellP_2     BellP_3     ChiliP_1    ChiliP_2    ChiliP_3    F1_1        F1_2        F1_3
Caper_00001     145         152         148         98          89          105         121         128         115
Caper_00002     2345        2156        2198        1876        1834        1923        2087        1945        2034
Caper_00003     0           0           0           0           0           0           0           0           0
Caper_00004     567         489         601         523         612         498         545         567         501
Caper_00005     890         765         823         1245        1150        1189        1023        945         1087
Caper_00006     123         145         112         98          102         94          110         124         108
Caper_00007     45          52          48          61          58          64          51          49          55
Caper_00008     234         267         245         198         212         201         221         234         215
Caper_00009     0           0           0           0           0           0           0           0           0
Caper_00010     1234        1456        1290        1512        1680        1545        1398        1567        1423
```

**Format**:
- Tab-separated
- First column: Gene IDs (pepper genes)
- Columns 2-10: Raw counts per sample
- No normalization
- NO `.T[digit]` suffixes (gene-level, not transcript-level)

---

### 2. Metadata: `pepper_metadata.tsv`

```
sample_name     group           species         tissue      replicate
BellP_1         BellPepper      bell_pepper     fruit       1
BellP_2         BellPepper      bell_pepper     fruit       2
BellP_3         BellPepper      bell_pepper     fruit       3
ChiliP_1        ChiliPepper     chili_pepper    fruit       1
ChiliP_2        ChiliPepper     chili_pepper    fruit       2
ChiliP_3        ChiliPepper     chili_pepper    fruit       3
F1_1            F1_hybrid       F1              fruit       1
F1_2            F1_hybrid       F1              fruit       2
F1_3            F1_hybrid       F1              fruit       3
```

**Key columns**:
- `sample_name`: Must match count matrix columns exactly
- `group`: Used for contrasts (BellPepper, ChiliPepper, F1_hybrid)
- Optional columns: species, tissue, replicate

**Contrasts that will be run**:
- BellPepper vs. ChiliPepper
- BellPepper vs. F1_hybrid
- ChiliPepper vs. F1_hybrid

---

### 3. ASE Counts: `pepper_ase.tsv`

```
gene_id         snp_count   bias_ratio   predominant_bias   BellPepper_ratio   ChiliPepper_ratio
Caper_00001     8           0.75         BellPepper        0.78               0.22
Caper_00002     4           0.50         BellPepper        0.51               0.49
Caper_00003     12          0.92         ChiliPepper       0.08               0.92
Caper_00004     6           0.67         BellPepper        0.65               0.35
Caper_00005     15          0.87         BellPepper        0.85               0.15
Caper_00006     2           0.50         ChiliPepper       0.48               0.52
Caper_00007     10          0.80         ChiliPepper       0.20               0.80
Caper_00008     5           0.60         BellPepper        0.62               0.38
Caper_00009     0           0.00         BellPepper        0.50               0.50
Caper_00010     9           0.78         BellPepper        0.76               0.24
```

**Column descriptions**:
- `gene_id`: Must match count matrix gene IDs (no suffixes)
- `snp_count`: Number of diagnostic SNPs distinguishing parents
- `bias_ratio`: Fraction of SNPs showing consistent bias
- `predominant_bias`: Which parent's allele is preferred
- `BellPepper_ratio`: Proportion of reads from Bell Pepper allele
- `ChiliPepper_ratio`: Proportion of reads from Chili Pepper allele

**Interpretation**:
- `Caper_00001`: 8 SNPs, 75% show bias, strongly favors Bell Pepper (78% vs 22%)
- `Caper_00002`: 4 SNPs, 50% show bias, essentially unbiased (51% vs 49%)
- `Caper_00003`: 12 SNPs, 92% show bias, strongly favors Chili Pepper (8% vs 92%)
- `Caper_00009`: No diagnostic SNPs, cannot determine bias

---

## Workflow Script

**File: `run_analysis.sh`**

```bash
#!/bin/bash
set -e

PROJECT_DIR="/path/to/pepper_project"
COUNTS="${PROJECT_DIR}/pepper_counts.tsv"
METADATA="${PROJECT_DIR}/pepper_metadata.tsv"
ASE_FILE="${PROJECT_DIR}/pepper_ase.tsv"
DE_OUT="${PROJECT_DIR}/results/DE"
ASE_OUT="${PROJECT_DIR}/results/ASE_DE"

echo "==============================================="
echo "Step 1: Running EdgeRDE (Differential Expression)"
echo "==============================================="

python SalmonStreamer.py EdgeRDE \
    --expression-file "${COUNTS}" \
    --metadata-file "${METADATA}" \
    --output-dir "${DE_OUT}" \
    --fdr-threshold 0.05 \
    --logfc-threshold 1.0

echo ""
echo "==============================================="
echo "Step 2: Running ASEIntegrate (Basic Analysis)"
echo "==============================================="

python SalmonStreamer.py ASEIntegrate \
    --de-results-dir "${DE_OUT}" \
    --ase-counts-file "${ASE_FILE}" \
    --output-dir "${ASE_OUT}/basic" \
    --fdr-threshold 0.05 \
    --parent1-label BellPepper \
    --parent2-label ChiliPepper

echo ""
echo "==============================================="
echo "Step 3: Running ASEIntegrate (Advanced Ad/Ed Analysis)"
echo "==============================================="

python SalmonStreamer.py ASEIntegrate \
    --de-results-dir "${DE_OUT}" \
    --ase-counts-file "${ASE_FILE}" \
    --output-dir "${ASE_OUT}/advanced" \
    --fdr-threshold 0.05 \
    --parent1-label BellPepper \
    --parent2-label ChiliPepper \
    --expression-file "${COUNTS}" \
    --metadata-file "${METADATA}" \
    --hybrid-group F1_hybrid \
    --parent1-group BellPepper \
    --parent2-group ChiliPepper

echo ""
echo "==============================================="
echo "Analysis Complete!"
echo "==============================================="
echo "DE results: ${DE_OUT}"
echo "ASE+DE results (basic): ${ASE_OUT}/basic"
echo "ASE+DE results (advanced): ${ASE_OUT}/advanced"
```

**Run it**:
```bash
bash run_analysis.sh
```

---

## Expected Results & Interpretation

### From EdgeRDE:

**DE Summary** (typical for this type of cross):
- ~2,000-3,000 significant genes (FDR < 0.05) between parents
- ~500-1,000 significant genes in each F1 comparison

**Example**: BellPepper vs. ChiliPepper contrast
```
Total genes: 15,000
Significant DE: 2,341 (15.6%)
Upregulated in BellPepper: 1,234
Upregulated in ChiliPepper: 1,107
```

**Quality checks**:
- PCA plot: Three clusters (one per group)
- Correlation heatmap: Replicates cluster together
- Library sizes: Similar across samples (range: 18-22 million reads)

---

### From ASEIntegrate (Basic):

**Gene category breakdown** (example):
```
Comparison: BellPepper_vs_ChiliPepper
Total DE genes: 2,341
├─ cis: 234 (10.0%)
├─ trans: 156 (6.7%)
├─ trans_compensatory: 89 (3.8%)
├─ compensatory: 45 (1.9%)
└─ no_divergence: 1,817 (77.6%)
```

**Interpretation**:
- ~10% of DE genes show cis regulation (strong candidates for enhancer work)
- ~7% are trans-regulated (coordinated expression changes)
- ~4% are trans-compensatory (interesting regulatory conflicts)
- ~78% lack diagnostic SNPs or show no ASE signal

**Top cis candidates** from `cis_regulatory_candidates_all.tsv`:
```
gene_id     logFC   FDR     snp_count   bias_ratio  predominant_bias
Caper_00001 1.23    0.0001  8           0.75        BellPepper
Caper_00005 -0.89   0.0015  15          0.87        BellPepper
Caper_00007 -1.45   0.0003  10          0.80        ChiliPepper
```

---

### From ASEIntegrate (Advanced):

**Ad/Ed classification** (subset of `ase_vs_de_ad_ed_data.tsv`):
```
gene_id     logFC   ad_est  ed_est  regulatory_class
Caper_00001 1.23    1.1     0.1     Cis_dominant
Caper_00002 0.89    0.2     0.7     Trans_dominant
Caper_00003 1.5     0.9     0.6     Cis_trans_reinforce
Caper_00007 -1.45   1.2     -2.4    Cis_trans_opposing
```

**Interpretation**:
- `Caper_00001`: Mostly cis (1.1 vs 0.1)
- `Caper_00002`: Mostly trans (0.2 vs 0.7)
- `Caper_00003`: Both cis and trans working together (same direction)
- `Caper_00007`: Cis increases expression, trans decreases (opposing)

---

## Next Steps from Example

### 1. Functional Annotation

```R
# In R, load cis candidates and search for enriched GO terms
library(topGO)
cis_genes <- read.csv("cis_regulatory_candidates_all.tsv", sep="\t")
# ... enrichment analysis
```

### 2. Validation (qRT-PCR)

```bash
# Design allele-specific primers for top cis genes
# Example primers for Caper_00001:
# Forward: GGCGATCTAGCTA
# Reverse (BellPepper-specific): GGTCGGATCCTA
# Reverse (ChiliPepper-specific): GGTCGGATCCTG
```

### 3. Sequence Alignment

Compare promoter/enhancer regions between parents for cis candidates:
- Find SNPs/indels in regulatory regions
- Look for transcription factor binding site differences

### 4. Phenotype Correlation

Do cis/trans regulatory differences correspond to phenotypic differences?
- Example: If fruit size genes are cis-regulated, size differences might be due to local regulation
- If color genes are trans-regulated, color might be controlled centrally (e.g., by TF master regulator)

---

## Tips for This Type of Analysis

1. **Check SNP distribution**: Do all genes have SNPs?
   - If `snp_count = 0` for most genes, you need more markers or re-sequencing

2. **Replicate quality**: Do F1 replicates cluster together on PCA?
   - If not, there might be contamination or heterogeneity in F1 population

3. **Parent divergence**: How many DE genes between parents?
   - If > 20% of genes show DE, parents are quite divergent (good for detecting divergence)
   - If < 1% of genes show DE, parents are very similar (harder to map effects)

4. **ASE in F1**: Do F1 samples show parent-of-origin effects?
   - Look at `parent1_ratio` distribution in F1 samples
   - If biased, might indicate imprinting or copy number variation

---

## Troubleshooting This Example

**Problem**: "No matching samples found"
- **Check**: Do sample_name values in metadata match count matrix column headers exactly?
- **Solution**: Verify spelling, capitalization, no trailing spaces

**Problem**: "Gene ID mismatch" or no cis genes found
- **Check**: Do gene IDs in count matrix match ASE file?
- **Solution**: Verify IDs are identical; ASE script auto-strips `.T[digit]` but not other variants

**Problem**: All genes classified as `no_divergence`
- **Check**: Are SNP counts zero for most genes?
- **Solution**: May need more diagnostic SNPs or parent resequencing to resolve alleles

---

## See Also

- [Input Format Specification](../docs/input_format_specification.md)
- [Output Files Reference](../docs/output_files_reference.md)
- [ASE + DE Integration Workflow](../docs/ase_de_integration_workflow.md)
- [Interpretation Guidelines](../docs/interpretation_guidelines.md)
