# ASE + DE Integration Example: Hybrid Cross (Illustrative)

This example walks through the full ASE + DE integration workflow for a cross between two parental lines. It is not tied to any organism: the lines are called `lineA` and `lineB`, the genes `GENE0001`, `GENE0002`, and so on.

**Every number in this file is invented.** The tables show the shape of the inputs and outputs, not results from a real dataset, and the four input files described below are not shipped with the repository.

## Experiment Overview

**Design**:
- Parent 1: `lineA`
- Parent 2: `lineB`
- F1 hybrid: cross between the two

**Samples**:
- 3 replicates of `lineA`
- 3 replicates of `lineB`
- 3 replicates of the F1 hybrid
- One tissue, sampled at the same stage in all three groups

**Total**: 9 samples

---

## Files in This Example

1. `counts.tsv` - Raw counts matrix
2. `metadata.tsv` - Sample metadata
3. `ase_counts.tsv` - Per-gene ASE table
4. `run_analysis.sh` - Complete workflow script

---

## File Contents

### 1. Count Matrix: `counts.tsv`

```
gene_id      lineA_1     lineA_2     lineA_3     lineB_1     lineB_2     lineB_3     F1_1        F1_2        F1_3
GENE0001     145         152         148         98          89          105         121         128         115
GENE0002     2345        2156        2198        1876        1834        1923        2087        1945        2034
GENE0003     0           0           0           0           0           0           0           0           0
GENE0004     567         489         601         523         612         498         545         567         501
GENE0005     890         765         823         1245        1150        1189        1023        945         1087
GENE0006     123         145         112         98          102         94          110         124         108
GENE0007     45          52          48          61          58          64          51          49          55
GENE0008     234         267         245         198         212         201         221         234         215
GENE0009     0           0           0           0           0           0           0           0           0
GENE0010     1234        1456        1290        1512        1680        1545        1398        1567        1423
```

**Format**:
- Tab-separated
- First column: Gene IDs
- Columns 2-10: Raw counts per sample
- No normalization
- No `.T[digit]` suffixes (gene-level, not transcript-level)

---

### 2. Metadata: `metadata.tsv`

```
sample_name     group           replicate
lineA_1         lineA           1
lineA_2         lineA           2
lineA_3         lineA           3
lineB_1         lineB           1
lineB_2         lineB           2
lineB_3         lineB           3
F1_1            F1_hybrid       1
F1_2            F1_hybrid       2
F1_3            F1_hybrid       3
```

**Key columns**:
- `sample_name`: Must match count matrix columns exactly
- `group`: Used for contrasts (`lineA`, `lineB`, `F1_hybrid`)
- Other columns (here `replicate`) are optional and ignored

**Contrasts that will be run** (EdgeRDE names each one `A_vs_B` with the group names in alphabetical order, and `logFC` is `A` minus `B`):
- `F1_hybrid_vs_lineA`
- `F1_hybrid_vs_lineB`
- `lineA_vs_lineB`

---

### 3. ASE Table: `ase_counts.tsv`

`ASEIntegrate` accepts two kinds of ASE table. If your ASE comes from [ASErtain](https://github.com/luismadrigal98/ASErtain), pass its `<prefix>.gene_ase.tsv` as it is, with `--parent1-label` / `--parent2-label` set to the lineage labels of the ASErtain config; that is the preferred input. This example uses the older summary table because it is small enough to read:

```
gene_id      snp_count   bias_ratio   predominant_bias   lineA_ratio   lineB_ratio
GENE0001     8           0.75         lineA              0.78          0.22
GENE0002     4           0.50         lineA              0.51          0.49
GENE0004     6           0.67         lineA              0.65          0.35
GENE0005     15          0.87         lineA              0.85          0.15
GENE0006     2           0.50         lineB              0.48          0.52
GENE0007     10          0.80         lineB              0.20          0.80
GENE0008     5           0.60         lineA              0.62          0.38
GENE0010     9           0.78         lineA              0.76          0.24
```

**Column descriptions**:
- `gene_id`: Must match count matrix gene IDs (no suffixes)
- `snp_count`: Number of diagnostic SNPs distinguishing parents
- `bias_ratio`: Fraction of SNPs showing bias; a gene counts as biased when this is above 0.5
- `predominant_bias`: Which parent's allele is preferred (`lineA` or `lineB`, the two labels)
- `lineA_ratio`: Proportion of F1 reads from the `lineA` allele
- `lineB_ratio`: Proportion of F1 reads from the `lineB` allele

**Reading the rows**:
- `GENE0001`: 8 SNPs, 75% show bias, favours `lineA` (78% vs 22%)
- `GENE0002`: 4 SNPs, 50% show bias, so the gene is not counted as biased (51% vs 49%)
- `GENE0007`: 10 SNPs, 80% show bias, favours `lineB` (20% vs 80%)
- `GENE0003` and `GENE0009` have no informative SNP and are left out of the table. A gene listed with `snp_count = 0` would still be classified, as unbiased

---

## Workflow Script

**File: `run_analysis.sh`**

```bash
#!/bin/bash
set -e

PROJECT_DIR="/path/to/project"
COUNTS="${PROJECT_DIR}/counts.tsv"
METADATA="${PROJECT_DIR}/metadata.tsv"
ASE_FILE="${PROJECT_DIR}/ase_counts.tsv"
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
    --parent1-label lineA \
    --parent2-label lineB \
    --parent1-group lineA \
    --parent2-group lineB

echo ""
echo "==============================================="
echo "Step 3: Running ASEIntegrate (Advanced Analysis)"
echo "==============================================="

python SalmonStreamer.py ASEIntegrate \
    --de-results-dir "${DE_OUT}" \
    --ase-counts-file "${ASE_FILE}" \
    --output-dir "${ASE_OUT}/advanced" \
    --fdr-threshold 0.05 \
    --parent1-label lineA \
    --parent2-label lineB \
    --expression-file "${COUNTS}" \
    --metadata-file "${METADATA}" \
    --hybrid-group F1_hybrid \
    --parent1-group lineA \
    --parent2-group lineB

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

Notes on the script:
- `--parent1-group` / `--parent2-group` name the parental comparison (`lineA_vs_lineB`) and its direction. In this example they could be omitted in Step 2, because the group names contain the labels and the comparison would be found from them.
- Step 3 writes everything Step 2 writes, plus the advanced files, so Step 2 can be skipped when you want both. The advanced analysis needs the R package `edgeR`.

---

## Expected Results & Interpretation

The file names below follow from the group names above. The numbers are illustrative; they are what a genome-scale run might look like, not the result of running the ten-gene matrix shown here.

### From EdgeRDE:

One set of files per contrast (`*_DE_results.tsv`, `*_significant_genes.tsv`, `*_volcano.pdf` / `.png`) for `F1_hybrid_vs_lineA`, `F1_hybrid_vs_lineB` and `lineA_vs_lineB`, plus the quality-control plots.

**Quality checks**:
- PCA plot: three clusters (one per group)
- Correlation heatmap: replicates cluster together
- Library sizes: similar across samples

---

### From ASEIntegrate (Basic):

Only the parental comparison is classified. The log says so for the other two:

```
Skipping: F1_hybrid_vs_lineA - not a parent 1 versus parent 2 comparison, so cis/trans classification does not apply.
Skipping: F1_hybrid_vs_lineB - not a parent 1 versus parent 2 comparison, so cis/trans classification does not apply.
Processing: lineA_vs_lineB (logFC used as named)
```

**Category summary** (`lineA_vs_lineB_ASE_category_summary.tsv`, illustrative):
```
ase_de_category      n_genes   avg_logFC_parent1_vs_parent2   avg_parent1_ratio
cis                  234       0.41                           0.61
compensatory         45        0.03                           0.58
no_divergence        2341      0.01                           0.50
trans_compensatory   89        -0.12                          0.57
trans_only           156       0.35                           0.50
```

**Interpretation**:
- The table covers the 2,865 genes that have both a DE row and an ASE row. Genes without an ASE row are not in it
- 479 of them are DE between the parents: 234 `cis`, 156 `trans_only`, 89 `trans_compensatory`
- 45 are allelically biased without being DE (`compensatory`), and 2,341 show neither (`no_divergence`)

**Top cis candidates** from `cis_regulatory_candidates_all.tsv` (selected columns, illustrative):
```
comparison       TranscriptID   logFC_parent1_vs_parent2   FDR      snp_count   bias_ratio   predominant_bias
lineA_vs_lineB   GENE0001       1.23                       0.0001   8           0.75         lineA
lineA_vs_lineB   GENE0007       -1.45                      0.0003   10          0.80         lineB
lineA_vs_lineB   GENE0005       0.89                       0.0015   15          0.87         lineA
```

In every `cis` row the bias points to the more highly expressed parent: toward `lineA` with a positive `logFC_parent1_vs_parent2`, toward `lineB` with a negative one.

---

### From ASEIntegrate (Advanced):

**Cis/trans classes** (selected columns of `ase_vs_de_ad_ed_data.tsv`, illustrative values):
```
gene_id    lineA_ratio   parents_differ   cis_log2   parental_log2   trans_log2   Ad      Ed      Regulation
GENE0001   0.78          TRUE             1.83       1.25            -0.58        -0.28   -0.05   Cis-only
GENE0002   0.51          TRUE             0.06       1.60            1.54         -0.01   -0.35   Trans-only
GENE0004   0.65          TRUE             0.89       2.10            1.21         -0.15   -0.10   Cis + Trans (Reinforcing)
GENE0007   0.20          TRUE             -2.00      1.45            3.45         0.30    0.30    Cis x Trans (Opposing)
GENE0010   0.76          FALSE            1.66       0.10            -1.56        -0.26   NA      Compensatory
```

**Interpretation**:
- `GENE0001`: the allelic ratio accounts for the parental difference; the trans component stays under its threshold (cis)
- `GENE0002`: the parents differ but the F1 alleles are balanced (trans)
- `GENE0004`: cis and trans push the same way
- `GENE0007`: the `lineB` allele is favoured in the F1 although `lineA` expresses more; cis and trans oppose
- `GENE0010`: the alleles are imbalanced but the parents are not differentially expressed, so `Ed` is not computed
- `Ed` describes dominance (where the F1's total expression sits between the parents), not the trans effect

The class thresholds are descriptive: "cis" is an allele fraction outside 0.5 ± 0.1, "trans" is `|trans_log2| ≥ 1`. See the [Output Files Reference](../docs/output_files_reference.md#ase_vs_de_ad_ed_datatsv).

---

## Next Steps from Example

### 1. Functional Annotation

```R
# In R, load the cis candidates as the gene set for an enrichment analysis
cis_genes <- read.delim("cis_regulatory_candidates_all.tsv")
# ... enrichment analysis with the tool of your choice
```

### 2. Validation (qRT-PCR)

Design allele-specific primers for the top cis candidates, with one primer placed on a diagnostic SNP so that each parental allele is amplified separately.

### 3. Sequence Alignment

Compare promoter/enhancer regions between parents for cis candidates:
- Find SNPs/indels in regulatory regions
- Look for transcription factor binding site differences

### 4. Phenotype Correlation

Do cis/trans regulatory differences correspond to phenotypic differences between the parents?
- A trait whose candidate genes are `cis` points to regulatory variation at those genes
- A trait whose candidate genes are `trans_only` points to a regulator elsewhere in the genome

---

## Tips for This Type of Analysis

1. **Check SNP distribution**: How many genes have an ASE row, and how many SNPs do they have?
   - Genes with few SNPs give uncertain calls; genes with none cannot be assessed

2. **Replicate quality**: Do F1 replicates cluster together on PCA?
   - If not, there might be contamination or heterogeneity in the F1 samples

3. **Parent divergence**: How many DE genes between parents?
   - Only parentally DE genes can be `cis`, `trans_only` or `trans_compensatory`

4. **Allele fraction across genes**: Look at the distribution of the parent 1 allele fraction
   - A distribution centred away from 0.5 across the whole genome should be understood (mapping bias is one cause) before single genes are interpreted

---

## Troubleshooting This Example

**Problem**: "No matching samples found"
- **Check**: Do sample_name values in metadata match count matrix column headers exactly?
- **Solution**: Verify spelling, capitalization, no trailing spaces

**Problem**: "No DE comparison could be identified as parent 1 versus parent 2"
- **Check**: Do `--parent1-group` / `--parent2-group` equal the two group names in the `*_DE_results.tsv` file name? Without them, do the group names contain the parent labels?
- **Check**: Do gene IDs in the DE table match the ASE table? The same error appears when the parental comparison was found but no gene has an ASE row
- **Solution**: Pass the two group names explicitly; verify the IDs are identical (`ASEIntegrate` removes `.T[digit]` from the DE side but no other variant)

**Problem**: All genes classified as `no_divergence`
- **Check**: Are any genes DE between the parents? Are any called biased (`is_biased`)?
- **Solution**: May need more diagnostic SNPs or deeper coverage to detect allelic bias

---

## See Also

- [Input Format Specification](../docs/input_format_specification.md)
- [Output Files Reference](../docs/output_files_reference.md)
- [ASE + DE Integration Workflow](../docs/ase_de_integration_workflow.md)
- [Interpretation Guidelines](../docs/interpretation_guidelines.md)
