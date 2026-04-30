# Output Files Reference

## Overview

This document describes all output files produced by `EdgeRDE` and `ASEIntegrate` subcommands, including file formats and how to interpret each output.

---

## EdgeRDE Outputs

Location: `--output-dir` specified during EdgeRDE run

### Per-Comparison Files

These files are generated for each pairwise contrast between groups.

#### `{comparison}_DE_results.tsv`

**Description**: Complete differential expression results for all genes in this contrast.

**Format**: Tab-separated values with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `TranscriptID` | string | Gene/transcript identifier |
| `logFC` | float | Log2-fold change (Group1 vs Group2) |
| `F-statistic` | float | Quasi-likelihood F-test statistic |
| `p-value` | float | Raw p-value from statistical test |
| `FDR` | float | Benjamini-Hochberg adjusted p-value (false discovery rate) |
| `Significance` | string | Category: "Significant" if FDR < threshold, "Not Significant" otherwise |

**Example**:
```
TranscriptID    logFC       F-statistic p-value     FDR         Significance
AT1G01010       2.156       45.23       1.2e-08     3.5e-07     Significant
AT1G01020       -1.345      12.34       0.0001      0.0015      Significant
AT1G01030       0.234       2.1         0.15        0.38        Not Significant
```

**Interpretation**:
- **Positive logFC**: Gene higher in Group1 (first group in comparison name)
- **Negative logFC**: Gene higher in Group2 (second group in comparison name)
- **FDR < 0.05**: Gene is significantly differentially expressed (default threshold)
- **|logFC| > 1**: Roughly 2-fold difference in expression

**Usage**: Use for downstream filtering, functional annotation, or literature searches.

---

#### `{comparison}_significant_genes.tsv`

**Description**: Filtered version of `*_DE_results.tsv` containing only significant genes (FDR < threshold).

**Columns**: Same as `*_DE_results.tsv`

**Number of rows**: Subset of `*_DE_results.tsv` (only genes with FDR < --fdr-threshold)

**Example**:
```
TranscriptID    logFC       F-statistic p-value     FDR         Significance
AT1G01010       2.156       45.23       1.2e-08     3.5e-07     Significant
AT1G01020       -1.345      12.34       0.0001      0.0015      Significant
```

**Usage**: Quick access to high-confidence DE genes; good for functional enrichment analysis.

---

#### `{comparison}_volcano.pdf` and `.png`

**Description**: Volcano plot visualizing DE results: logFC vs. -log10(p-value).

**Visual elements**:
- **X-axis**: log2-fold change (logFC)
- **Y-axis**: -log10(p-value)
- **Point colors**:
  - **Red points**: Significant genes (FDR < 0.05, typically above red horizontal line)
  - **Black points**: Non-significant genes
- **Horizontal dashed line**: -log10(0.05) threshold (FDR cutoff)
- **Vertical dashed lines**: logFC = ±1 (2-fold change boundaries)

**Interpretation**:
- **Upper left quadrant**: Significant DOWN-regulated genes
- **Upper right quadrant**: Significant UP-regulated genes
- **Center**: Not significantly different
- **Wide spread on right**: Many genes with large fold changes
- **Wide spread on top**: Many genes with strong statistical significance

**Usage**: Visual summary for presentations or papers; quick assessment of DE landscape.

---

### Global Quality Control Files

These files provide overall analysis quality and sample relationships.

#### `PCA_plot.pdf` and `.png`

**Description**: Principal Component Analysis plot showing sample relationships.

**Visual elements**:
- **X-axis**: PC1 (typically explains most variation)
- **Y-axis**: PC2
- **Points**: Individual samples, colored by group
- **Polygons**: Convex hulls around each group
- **Title**: Percentage of variance explained by each PC

**Interpretation**:
- **Clear separation**: Groups are expression-wise different (good for DE)
- **Overlap**: Groups have similar expression profiles
- **Outliers**: Points far from their group's cluster (potential QC issues)
- **Technical variation dominates**: Samples cluster by batch, not biology (potential confound)

**Usage**: Assess whether groups are distinct; identify sample outliers; verify metadata grouping is correct.

**Example**:
```
If PC1 = 45%, PC2 = 20%:
- First two PCs explain 65% of total variation
- Groups clearly separated on PC1: strong effect
```

---

#### `sample_correlation_heatmap.pdf`

**Description**: Hierarchical clustering of samples based on expression correlation.

**Visual elements**:
- **Color scale**: Red (high correlation) to White (low correlation)
- **Rows/Columns**: Sample names
- **Clustering**: Dendrogram showing which samples are most similar

**Interpretation**:
- **Within-group clustering**: Samples from same group cluster together (good)
- **Between-group mixing**: Groups don't separate cleanly (suggests low effect size or high noise)
- **Outlier samples**: Hang separately from their group's cluster

**Usage**: Identify samples with unusual expression profiles; verify replicate consistency.

---

#### `library_sizes.pdf` and `.png`

**Description**: Barplot of total reads per sample (library size / sequencing depth).

**X-axis**: Sample names  
**Y-axis**: Total reads (millions or billions, depending on sequencing depth)

**Interpretation**:
- **Similar heights**: Samples sequenced to similar depth (good)
- **Highly variable**: Some samples have 2-10x more reads than others
  - **Not necessarily bad**: edgeR's TMM normalization accounts for this
  - **Shows in output**: `normalization_factors.txt` shows how much each sample was adjusted

**Usage**: Identify severely under-sequenced samples; verify library prep consistency.

---

#### `dispersion_plot.pdf` and `.png`

**Description**: edgeR dispersion estimates vs. average read count.

**X-axis**: Average log10(count per million)  
**Y-axis**: Square root of dispersion (tagwise or trended)

**Points**: One per gene

**Red curves**: Trended and common dispersion estimates

**Interpretation**:
- **Points below curve**: Gene-specific dispersion lower than expected (more stable)
- **Points above curve**: Gene-specific dispersion higher than expected (more variable)
- **Spread decreases left-to-right**: Low-count genes more variable (expected)
- **High spread overall**: Might indicate overdispersion or strong biological heterogeneity

**Usage**: Assess whether statistical model assumptions are met; identify genes with unusual variance.

---

#### `QL_dispersion_plot.pdf` and `.png`

**Description**: Quasi-likelihood dispersion plot (from glmQLFTest, an alternative to GLM-LRT).

**Similar to `dispersion_plot.pdf`** but for quasi-likelihood fitting.

**Interpretation**: Similar to dispersion plot; shows model fit quality.

---

#### `DE_genes_heatmap.pdf`

**Description**: Heatmap of the top most-significant DE genes.

**Rows**: Top N significant genes (typically 20-50)  
**Columns**: Samples (colored by group)  
**Color scale**: Log CPM values (blue = low, red = high)

**Interpretation**:
- **Clear row/column patterns**: Strong consistent effects within groups
- **Mixed coloring within groups**: High within-group variation
- **Obvious block structure**: Good statistical separation between groups

**Usage**: Visually verify that top DE genes show expected patterns.

---

#### `analysis_summary.txt`

**Description**: Text summary of analysis parameters and results.

**Contents**:
```
=== SalmonStreamer EdgeRDE Analysis ===
Input file     : counts.tsv
Output dir     : DE_results/
Metadata file  : metadata.tsv
FDR threshold  : 0.05
logFC threshold: 1.0

Data dimensions: 15000 genes x 12 samples
Stripped suffix '' from sample names
Samples retained: 12
Groups          : F1_hybrid, Parent1, Parent2

Genes before filtering: 15000
Genes after filterByExpr: 12340

TMM factors: 1.002, 1.015, 0.998, ...

Contrasts: 
  F1_hybrid vs Parent1: 2344 significant genes
  F1_hybrid vs Parent2: 2156 significant genes
  Parent1 vs Parent2: 1823 significant genes
```

**Usage**: Quick reference for analysis parameters; record in lab notebook or methods.

---

#### `session_info.txt`

**Description**: R session information (R version, package versions, etc.).

**Contents**:
```
R version 4.x.x
Platform: x86_64-linux-gnu
Packages: edgeR 3.40.x, limma 3.54.x, ggplot2 3.4.x, pheatmap 1.0.x
```

**Usage**: Reproducibility; record versions for methods section.

---

### Summary of EdgeRDE Output Files

| File | Type | Count | Purpose |
|------|------|-------|---------|
| `{comparison}_DE_results.tsv` | Data | N = comparisons | Complete DE stats per contrast |
| `{comparison}_significant_genes.tsv` | Data | N = comparisons | Filtered significant genes |
| `{comparison}_volcano.pdf` | Plot | N = comparisons | Visual summary per contrast |
| `PCA_plot.pdf` | Plot | 1 | Sample relationships |
| `sample_correlation_heatmap.pdf` | Plot | 1 | Sample similarity |
| `library_sizes.pdf` | Plot | 1 | Sequencing depth |
| `dispersion_plot.pdf` | Plot | 1 | Variance model fit |
| `QL_dispersion_plot.pdf` | Plot | 1 | Alternative variance model |
| `DE_genes_heatmap.pdf` | Plot | 1 | Top DE genes expression |
| `analysis_summary.txt` | Text | 1 | Parameters & summary stats |
| `session_info.txt` | Text | 1 | R environment info |

**Total files**: ~11 + (2-3 × number of comparisons)

---

## ASEIntegrate Outputs

Location: `--output-dir` specified during ASEIntegrate run

### Per-Comparison Files

#### `{comparison}_with_ASE.tsv`

**Description**: Merged DE and ASE results with regulatory categorization.

**Columns**:

| Column | Type | Source | Description |
|--------|------|--------|-------------|
| `TranscriptID` | string | DE results | Gene/transcript identifier |
| `logFC` | float | DE results | Log2-fold change |
| `p_value` | float | DE results | Raw p-value |
| `FDR` | float | DE results | Adjusted p-value |
| `gene_category` | string | Integration | Regulatory classification (see below) |
| `snp_count` | int | ASE file | Number of diagnostic SNPs |
| `bias_ratio` | float | ASE file | Fraction of SNPs with bias |
| `predominant_bias` | string | ASE file | Which parent is biased ("parent1" or "parent2") |
| `parent1_ratio` | float | ASE file | Fraction reads from parent 1 allele |
| `parent2_ratio` | float | ASE file | Fraction reads from parent 2 allele |

**Gene Categories**:

| Category | DE? | ASE Bias? | Meaning |
|----------|-----|-----------|---------|
| `cis` | Yes (FDR<0.05) | Yes (biased) | Regulatory sequence differences |
| `trans` | Yes (FDR<0.05) | No (unbiased) | Trans-acting factor effects |
| `trans_compensatory` | Yes (FDR<0.05) | Yes (opposite to DE) | Trans factors overcorrect cis diff |
| `compensatory` | No (FDR≥0.05) | Yes (biased) | Cis differences masked by trans |
| `no_divergence` | No (FDR≥0.05) | No (unbiased) | No regulatory divergence |

**Example**:
```
TranscriptID    logFC   FDR     gene_category      snp_count   bias_ratio  predominant_bias   parent1_ratio   parent2_ratio
AT1G01010       2.15    0.001   cis                8           0.75        parent1            0.78            0.22
AT1G01020       1.89    0.008   trans              4           0.50        parent1            0.51            0.49
AT1G01030      -1.34    0.042   trans_compensatory 12          0.92        parent2            0.08            0.92
AT1G01040       0.12    0.89    compensatory       6           0.67        parent1            0.65            0.35
AT1G01050      -0.08    0.95    no_divergence      0           0.00        parent1            0.50            0.50
```

**Usage**: Primary output for analysis; filter by category for downstream work.

---

#### `{comparison}_ASE_category_summary.tsv`

**Description**: Summary counts of genes per category.

**Columns**:

| Column | Description |
|--------|-------------|
| `comparison` | Contrast name (e.g., "F1_hybrid_vs_Parent1") |
| `n_cis` | Number of cis-regulated genes |
| `n_trans` | Number of trans-regulated genes |
| `n_trans_compensatory` | Number of trans-compensatory genes |
| `n_compensatory` | Number of compensatory genes |
| `n_no_divergence` | Number with no divergence |
| `n_total_DE` | Total DE genes (cis + trans + trans_comp) |
| `n_total_ASE` | Total with ASE info (all except no_divergence) |
| `n_analyzed` | Total genes with both DE + ASE data |

**Example**:
```
comparison              n_cis   n_trans n_trans_comp    n_compensatory  n_no_divergence
F1_hybrid_vs_Parent1    234     156     89              45              2341
F1_hybrid_vs_Parent2    198     201     76              52              2288
Parent1_vs_Parent2      187     143     51              38              2301
```

**Usage**: Quick summary statistics for tables and figures.

---

#### `{comparison}_ASE_vs_DE.pdf` and `.png`

**Description**: Scatter plot of logFC (x-axis) vs. allele ratio (y-axis).

**Visual elements**:
- **X-axis**: logFC (expression difference)
- **Y-axis**: Allele ratio bias (0 = all parent2, 1 = all parent1)
- **Point colors**: Categories
  - Red: cis-regulated
  - Blue: trans-regulated
  - Green: trans-compensatory
  - Gray: compensatory
  - Black: no divergence
- **Reference lines**:
  - Horizontal at y = 0.5 (unbiased alleles)
  - Vertical at x = 0 (no expression difference)

**Interpretation**:
- **Upper right**: Cis + stronger in parent1 (expected concordance)
- **Upper left**: cis + stronger in parent2 (expected concordance)
- **Center + red dots**: Cis effects opposing expression difference
- **Center + blue dots**: Trans effects with ~balanced alleles
- **Spread of colors**: Mix of regulatory mechanisms

**Usage**: Visual summary of cis vs. trans regulatory landscape.

---

### Global Summary Files

#### `cis_regulatory_candidates_all.tsv`

**Description**: All genes classified as `cis` across all comparisons (useful for prioritizing candidates).

**Columns**: Same as `*_with_ASE.tsv`

**Rows**: Only genes with `gene_category == "cis"` in at least one comparison

**Additional column** (optional): `comparisons_with_cis`
- List of all comparisons where gene showed cis regulation

**Usage**: Create candidate gene list for functional validation or fine-mapping studies.

---

#### `ase_vs_de_ad_ed_data.tsv`

**Description**: Advanced regulatory classification using Ad/Ed framework (only if `--expression-file` provided).

**Columns**:

| Column | Type | Description |
|--------|------|-------------|
| `TranscriptID` | string | Gene identifier |
| `comparison` | string | Contrast name |
| `logFC_observed` | float | Observed log2-fold change in hybrid |
| `ad_estimate` | float | Estimated cis effect size |
| `ed_estimate` | float | Estimated trans effect size |
| `ad_plus_ed` | float | Sum of cis + trans (should ≈ logFC_observed) |
| `regulatory_class` | string | Classification (see below) |

**Regulatory Classifications**:

| Class | Meaning |
|-------|---------|
| `Cis_dominant` | Cis effect >> trans effect; mostly explains logFC |
| `Trans_dominant` | Trans effect >> cis effect; mostly explains logFC |
| `Cis_trans_reinforcing` | Both effects in same direction; additive |
| `Cis_trans_opposing` | Cis and trans effects in opposite directions |
| `Trans_only` | Cis effect ≈ 0; all logFC from trans |
| `Ambiguous` | Effects unclear or don't add up well |

**Example**:
```
TranscriptID    comparison          logFC_observed  ad_estimate  ed_estimate  ad_plus_ed  regulatory_class
AT1G01010       F1_hybrid_vs_P1     2.15            1.8          0.35         2.15        Cis_dominant
AT1G01020       F1_hybrid_vs_P1     1.89            0.3          1.6          1.9         Trans_dominant
AT1G01030       F1_hybrid_vs_P1    -1.34            0.8         -2.1         -1.3         Cis_trans_opposing
```

**Interpretation**:
- **Cis-dominant**: Regulatory mutations in this gene's promoter/enhancer likely explain expression difference
- **Trans-dominant**: Expression difference driven by factors acting on both alleles (TFs, chromatin modifiers)
- **Cis-trans opposing**: Cis changes would increase expression, but trans factors reduce it (and vice versa)
- **Reinforcing**: Cis and trans changes work together to create observed difference

**Usage**: Mechanistic understanding of expression evolution; prioritize genes for enhancer mapping or EMSA.

---

#### `regulatory_divergence_regressions.pdf`

**Description**: Plots showing Ad/Ed estimates with regression lines (only if `--expression-file` provided).

**Content**:
- Multiple panels: One per contrast
- **X-axis**: Cis effect estimate (ad)
- **Y-axis**: Trans effect estimate (ed)
- **Points**: Genes; color by regulatory class
- **Ellipses**: Confidence regions
- **Regression line**: Mean relationship between cis and trans

**Interpretation**:
- **Points in quadrant I (upper right)**: Cis and trans reinforcing
- **Points in quadrant III (lower left)**: Cis and trans opposing
- **Spread of points**: Heterogeneity in regulatory mechanisms across genes

**Usage**: Visualize overall regulatory landscape; publication figures.

---

#### `session_info.txt`

**Description**: R session information.

**Same format as EdgeRDE's `session_info.txt`**

---

### Summary of ASEIntegrate Output Files

| File | Type | Count | Purpose |
|------|------|-------|---------|
| `{comparison}_with_ASE.tsv` | Data | N = comparisons | Merged DE + ASE per contrast |
| `{comparison}_ASE_category_summary.tsv` | Data | N = comparisons | Category counts per contrast |
| `{comparison}_ASE_vs_DE.pdf` | Plot | N = comparisons | Scatter: logFC vs. allele bias |
| `cis_regulatory_candidates_all.tsv` | Data | 1 | All cis genes (all comparisons) |
| `ase_vs_de_ad_ed_data.tsv` | Data | 1 | Ad/Ed classification (if advanced mode) |
| `regulatory_divergence_regressions.pdf` | Plot | 1 | Ad/Ed visualization (if advanced mode) |
| `session_info.txt` | Text | 1 | R environment info |

**Total files**: ~7 + (2-3 × number of comparisons), or ~9-12 if advanced mode enabled

---

## File Organization Best Practices

### Recommended Directory Structure

```
project/
├── data/
│   ├── counts.tsv
│   ├── metadata.tsv
│   └── ase_counts.tsv
├── results/
│   ├── DE_results/
│   │   ├── F1_vs_P1_DE_results.tsv
│   │   ├── F1_vs_P1_volcano.pdf
│   │   ├── PCA_plot.pdf
│   │   ├── analysis_summary.txt
│   │   └── ...
│   ├── ASE_DE_results/
│   │   ├── F1_vs_P1_with_ASE.tsv
│   │   ├── F1_vs_P1_ASE_vs_DE.pdf
│   │   ├── cis_regulatory_candidates_all.tsv
│   │   └── ...
│   └── figures/  (manually curated plots for papers)
└── scripts/
    ├── 01_run_edgerde.sh
    └── 02_run_aseintegrate.sh
```

### Naming Convention for Analysis Runs

Use date or version numbers:
```
results_2026-04-30/
results_v1/
results_final/
```

This prevents accidentally overwriting results.

---

## Interpreting Combined Results

### Workflow Example: Pepper Cross

1. **EdgeRDE** identifies:
   - 2,344 DE genes (F1 vs Parent1, FDR<0.05)
   - Genes cluster well on PCA

2. **ASEIntegrate** categorizes those 2,344 into:
   - 234 cis-regulated (allele-biased, DE concordant)
   - 156 trans-regulated (unbiased alleles, DE still present)
   - 89 trans-compensatory (allele-biased opposite to DE)
   - 1,865 unmapped to ASE or ambiguous

3. **Summary interpretation**:
   - ~10% of DE genes show strong cis divergence
   - ~7% show interesting trans-compensatory pattern
   - ~80% lack diagnostic SNPs (cannot determine cis vs. trans)

4. **Next steps**:
   - Prioritize 234 cis candidates for enhancer mapping
   - Investigate 89 trans-compensatory genes for TF changes
   - Use Advanced mode to quantify effect sizes

---

## See Also

- [ASE + DE Integration Workflow](./ase_de_integration_workflow.md)
- [Input Format Specification](./input_format_specification.md)
- [Interpretation Guidelines](./interpretation_guidelines.md)
- [ASEIntegrate & EdgeRDE Audit Report](../AUDIT_ASEIntegrate_EdgeRDE.md)
