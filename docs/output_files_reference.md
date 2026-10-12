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
GENE0001        2.156       45.23       1.2e-08     3.5e-07     Significant
GENE0002        -1.345      12.34       0.0001      0.0015      Significant
GENE0003        0.234       2.1         0.15        0.38        Not Significant
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
GENE0001        2.156       45.23       1.2e-08     3.5e-07     Significant
GENE0002        -1.345      12.34       0.0001      0.0015      Significant
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

### Optional Normalized Expression Exports

Produced only when `--export-normalized-expression` is set. Format is selected by `--normalized-expression-format {tsv,csv}` (default `tsv`).

#### `TMM_normalized_CPM.{tsv,csv}`

**Description**: TMM-normalized counts-per-million matrix, computed as `edgeR::cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)` after `calcNormFactors(method = "TMM")` and `filterByExpr` filtering.

**Columns**:
- `TranscriptID`: gene / transcript identifier (same set as the `*_DE_results.tsv` files, i.e. post-filter)
- One column per retained sample, header = sample name (after any `--sample-suffix` stripping)

**Usage**: Manual inspection of expression on a sample-comparable scale; input for downstream visualization, clustering, or co-expression analyses outside the pipeline.

#### `TMM_normalized_logCPM.{tsv,csv}`

**Description**: log2-CPM matrix, computed as `edgeR::cpm(dge, log = TRUE, normalized.lib.sizes = TRUE)`. This is the exact matrix EdgeRDE uses internally for PCA and the sample-correlation heatmap, so values can be cross-checked against those plots.

**Columns**: same shape as `TMM_normalized_CPM` but on the log2 scale.

**Usage**: Recommended for plotting, heatmaps, and any downstream method that expects approximately variance-stabilized expression values.

#### `TMM_norm_factors.{tsv,csv}`

**Description**: Per-sample library size and TMM normalization factors used to derive the CPM matrices above.

**Columns**:
- `Sample`: sample identifier
- `Group`: experimental group from the metadata
- `LibSize`: raw library size (sum of counts after filtering)
- `NormFactor`: TMM normalization factor from `calcNormFactors`
- `EffectiveLibSize`: `LibSize * NormFactor` — the denominator effectively used by `cpm()`

**Usage**: Sanity-check normalization (factors should hover around 1.0); reproduce CPM values manually; flag samples with unusually low effective library sizes.

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
| `TMM_normalized_CPM.{tsv,csv}` | Data | 0 or 1 | TMM-normalized CPM matrix (optional) |
| `TMM_normalized_logCPM.{tsv,csv}` | Data | 0 or 1 | log2-CPM matrix used for PCA (optional) |
| `TMM_norm_factors.{tsv,csv}` | Data | 0 or 1 | Per-sample lib size & TMM factors (optional) |

**Total files**: ~11 + (2-3 × number of comparisons), plus 3 more when `--export-normalized-expression` is set

---

## EdgeRDEFromNormalized Outputs

Location: `--output-dir` specified during `EdgeRDEFromNormalized` run.

This subcommand runs limma-trend on a pre-normalized expression matrix (log-CPM or CPM with internal log-transform). Its output set is a strict subset of `EdgeRDE` — it drops the count-model-specific QC plots (`dispersion_plot.pdf`, `QL_dispersion_plot.pdf`, `library_sizes.pdf`) and the TMM export files (TMM is not re-run), but keeps everything ASEIntegrate consumes.

### File set

| File | Type | Count | Purpose |
|------|------|-------|---------|
| `{comparison}_DE_results.tsv` | Data | N = comparisons | Complete DE stats per contrast (limma-trend) |
| `{comparison}_significant_genes.tsv` | Data | N = comparisons | Filtered to FDR < `--fdr-threshold` |
| `{comparison}_volcano.{pdf,png}` | Plot | N = comparisons | Volcano plot per contrast |
| `PCA_plot.{pdf,png}` | Plot | 1 | Sample clustering on the input log-CPM matrix |
| `sample_correlation_heatmap.pdf` | Plot | 1 | Pearson correlations between samples |
| `DE_genes_heatmap.pdf` | Plot | 0 or 1 | Top DE genes heatmap (only when any genes pass `--fdr-threshold`) |
| `analysis_summary.txt` | Text | 1 | Parameters, dimensions, per-contrast significant counts |
| `session_info.txt` | Text | 1 | R environment / package versions |

**Total files**: ~5 + (2-3 × number of comparisons).

### DE table columns

`{comparison}_DE_results.tsv` keeps the same schema as `EdgeRDE` for downstream compatibility with `ASEIntegrate`:

| Column | Source | Meaning |
|--------|--------|---------|
| `TranscriptID` | rownames | Gene / transcript identifier |
| `logFC` | limma | log2 fold change for the contrast (group1 − group2) |
| `logCPM` | limma `AveExpr` | Mean log2-CPM across all samples (renamed) |
| `t` | limma | Moderated t-statistic |
| `PValue` | limma `P.Value` | Raw p-value (renamed) |
| `FDR` | limma `adj.P.Val` | Benjamini-Hochberg adjusted p-value (renamed) |
| `B` | limma | log-odds of differential expression |
| `Significance` | added | `Not Significant` / `Significant` / `Highly Significant` based on `--fdr-threshold` and `--logfc-threshold` |

Note that `EdgeRDE` instead produces an `F` column (quasi-likelihood F-statistic) where `EdgeRDEFromNormalized` produces `t` and `B`. ASEIntegrate uses `logFC` and `FDR` only, so the swap is harmless.

### Notes on interpretation

- **No `library_sizes.pdf` / `dispersion_plot.pdf`**: these are diagnostics of the count model; they don't apply when the input is already normalized.
- **PCA is on the input log-CPM matrix**, not a re-normalized one. If you've collapsed paralogs, the PCA reflects the collapsed structure — which is usually what you want.
- **`analysis_summary.txt`** records the resolved `input_scale` (`logcpm` vs `cpm`) and the `prior_count` actually used, so you can reproduce the log-transform later if needed.

---

## ASEIntegrate Outputs

Location: `--output-dir` specified during ASEIntegrate run

Files are written for the **parent 1 versus parent 2** comparison only. Cis/trans classification is not defined for parent-vs-F1 or other contrasts, so those `*_DE_results.tsv` files are skipped (the log names each one). `{comparison}` below is the name of the classified DE file without `_DE_results.tsv`, for example `Parent1_vs_Parent2`. How the parental comparison is found is described in the [workflow guide](./ase_de_integration_workflow.md#which-comparison-files-to-read). If none is found, or none of its genes has a row in the ASE table, the run stops with an error and writes no files.

Example values in this section are illustrative, not results from a real dataset.

### Per-Comparison Files

#### `{comparison}_with_ASE.tsv`

**Description**: The DE table joined to the ASE table, one row per gene present in both, in the order of the DE table, with the regulatory category. Genes without an ASE row are not written.

**Columns**, in order:

1. Every column of the DE table, unchanged (from `EdgeRDE`: `TranscriptID`, `logFC`, `logCPM`, `F`, `PValue`, `FDR`, `Significance`).
2. `gene_id_base` and `logFC_parent1_vs_parent2`.
3. The ASE columns, which depend on the kind of ASE table:
   - **ASErtain gene table**: `snp_count` (ASErtain's `n_snps`), `bias_ratio` (empty), `predominant_bias` (ASErtain's `direction`), `ase_source`, `<parent1_label>_ratio` and `<parent2_label>_ratio` (from `mean_<label>_ratio`; the other one is one minus it), and `ase_q_value`, `ase_log2_ratio`, `ase_phase_discordance_cause` (ASErtain's `q_value`, `log2_ratio` and `phase_discordance_cause`, copied as they are). No other ASErtain column is carried over.
   - **Summary (legacy) table**: every column of that table except `gene_id`, then `ase_source`.
4. `is_DE`, `is_biased`, `bias_p1`, `ase_de_category`, `cis_log2`, `trans_log2`.

**Columns added by ASEIntegrate**:

| Column | Type | Description |
|--------|------|-------------|
| `gene_id_base` | string | `TranscriptID` with any `.T[digit]` suffix removed; the join key, matched against `gene_id` of the ASE table |
| `logFC_parent1_vs_parent2` | float | The fold change oriented so that positive = higher in parent 1. Equals `logFC` when parent 1's group is named first in the comparison, `-logFC` otherwise |
| `ase_source` | string | `asertain` or `legacy` |
| `is_DE` | logical | `FDR` below the ASEIntegrate `--fdr-threshold`. The size of the fold change is not used, so this is not EdgeRDE's `Significance` column |
| `is_biased` | logical | Allelic bias. With an ASErtain table: ASErtain's `ase_call`. With a summary table: `bias_ratio > 0.5` (a majority vote over SNPs, not a test) |
| `bias_p1` | logical | `predominant_bias` equals `--parent1-label`. Every other value is read as "not parent 1", so a biased gene with a misspelled label is treated as biased toward parent 2 |
| `ase_de_category` | string | `cis`, `trans_only`, `trans_compensatory`, `compensatory` or `no_divergence` |
| `cis_log2` | float | log2 of the F1 allelic ratio, parent 1 over parent 2: `log2(p / (1 - p))` with `p` the parent 1 allele fraction, limited to 0.001–0.999 (descriptive) |
| `trans_log2` | float | `logFC_parent1_vs_parent2 - cis_log2` (descriptive) |

**Category rule** (on `logFC_parent1_vs_parent2`):

| `ase_de_category` | DE between parents (`is_DE`) | Allelic bias (`is_biased`) |
|-------------------|--------------------|--------------|
| `cis` | yes | yes, toward the more highly expressed parent |
| `trans_compensatory` | yes | yes, toward the less highly expressed parent |
| `trans_only` | yes | no |
| `compensatory` | no | yes |
| `no_divergence` | no | no |

"Toward the more highly expressed parent" means bias toward parent 1 with a positive `logFC_parent1_vs_parent2`, or toward parent 2 with a negative one.

`snp_count` is not part of the rule. In a summary table, a gene listed with `snp_count = 0` and `bias_ratio = 0` counts as unbiased, so it comes out `trans_only` when it is DE and `no_divergence` when it is not. Leave such genes out of the table if you do not want them classified.

**Example** (selected columns, summary-table input):
```
TranscriptID   logFC   logFC_parent1_vs_parent2   FDR     predominant_bias   parent1_ratio   is_DE   is_biased   ase_de_category      cis_log2   trans_log2
GENE0001       2.15    2.15                       0.001   parent1            0.78            TRUE    TRUE        cis                  1.83       0.32
GENE0002       1.89    1.89                       0.008   parent1            0.51            TRUE    FALSE       trans_only           0.06       1.83
GENE0003       1.34    1.34                       0.042   parent2            0.08            TRUE    TRUE        trans_compensatory   -3.52      4.86
GENE0004       0.12    0.12                       0.89    parent1            0.65            FALSE   TRUE        compensatory         0.89       -0.77
GENE0005       -0.08   -0.08                      0.95    parent2            0.49            FALSE   FALSE       no_divergence        -0.06      -0.02
GENE0006       -1.60   -1.60                      0.002   parent2            0.20            TRUE    TRUE        cis                  -2.00      0.40
```

#### `{comparison}_ASE_category_summary.tsv`

**Description**: One row per category that has at least one gene. The same table is printed in the log.

| Column | Description |
|--------|-------------|
| `ase_de_category` | Category name |
| `n_genes` | Number of genes |
| `avg_logFC_parent1_vs_parent2` | Mean oriented fold change of those genes |
| `avg_parent1_ratio` | Mean parent 1 allele fraction of those genes |

**Example**:
```
ase_de_category      n_genes   avg_logFC_parent1_vs_parent2   avg_parent1_ratio
cis                  234       0.41                           0.61
compensatory         45        0.03                           0.58
no_divergence        2341      0.01                           0.50
trans_compensatory   89        -0.12                          0.57
trans_only           156       0.35                           0.50
```

Only genes with both a DE row and an ASE row are counted. A DE gene that has no ASE row is in none of these rows.

#### `{comparison}_ASE_vs_DE.pdf` and `.png`

**Description**: Scatter plot of the oriented fold change (x) against the parent 1 allele fraction (y), coloured by `ase_de_category` (ColorBrewer Set2), with dashed reference lines at x = 0 and y = 0.5.

**Visual elements**:
- **Title**: `ASE bias vs Differential Expression: {comparison}`
- **X-axis**: `log2 Fold Change (<parent1_label> / <parent2_label>)`, the `logFC_parent1_vs_parent2` column
- **Y-axis**: `<parent1_label> allele ratio (ASE)`, the `<parent1_label>_ratio` column (0 = all parent 2, 1 = all parent 1)
- **Colours**: assigned to the categories present, in alphabetical order. A category can therefore have a different colour in another plot where one category is missing; read the legend under the plot.

**Interpretation**: the upper-right and lower-left quadrants are where expression and allelic bias agree (cis-regulatory divergence). The upper-left and lower-right quadrants are where they oppose (`trans_compensatory`).

The colour comes from the category, and the category uses the bias call (`is_biased`), not the position on the y-axis. A gene far from 0.5 that was not called biased is still drawn as `trans_only` or `no_divergence`.

### Global Summary Files

#### `cis_regulatory_candidates_all.tsv`

**Description**: Every gene labelled `cis`, with a leading `comparison` column followed by the columns of `{comparison}_with_ASE.tsv`, sorted by comparison and then by FDR. The file is written only when at least one gene is `cis`.

There is normally one classified comparison. When the parental comparison is inferred from the labels and more than one file qualifies (a design with several tissues, for instance), each one contributes its own rows, so a gene can appear once per comparison.

#### `ase_vs_de_ad_ed_data.tsv`

**Description**: Advanced mode only (`--expression-file`, `--metadata-file`, `--hybrid-group`, `--parent1-group`, `--parent2-group`). One row per gene present in both the ASE table and the expression matrix, matched after removing any `.T[digit]` suffix from the matrix IDs.

| Column | Type | Description |
|--------|------|-------------|
| `parent1_avg`, `parent2_avg`, `hybrid_avg` | float | Mean TMM-normalised CPM per group |
| `parental_log2` | float | `log2((parent1_avg + 0.5) / (parent2_avg + 0.5))` |
| `parental_FDR` | float | FDR of the gene in the parental DE comparison; empty for a gene that is not in that DE table |
| `parents_differ` | logical | `parental_FDR` is below `--fdr-threshold`. False for a gene that is not in the DE table |
| `Ad` | float | `0.5 -` parent 1 allele fraction. Negative = parent 1 allele over-represented |
| `parental_diff` | float | `parent1_avg - parent2_avg` |
| `position_ratio` | float | `(hybrid_avg - parent2_avg) / parental_diff`; empty unless `parents_differ` |
| `Ed` | float | `0.5 - position_ratio`: departure of the F1's total expression from mid-parent. Negative = parent 1-like; empty unless `parents_differ` |
| `cis_log2` | float | log2 of the F1 allelic ratio, parent 1 over parent 2, as in `{comparison}_with_ASE.tsv` |
| `trans_log2` | float | `parental_log2 - cis_log2`. Here the parental ratio comes from mean CPM, so the value differs from `trans_log2` in `{comparison}_with_ASE.tsv`, which uses the DE fold change |
| `Regulation` | string | Class from the standard decomposition (below) |

These columns follow the ASE columns (for an ASErtain table, the mapped ones listed under `{comparison}_with_ASE.tsv`), `is_biased_input` (the bias call that the basic categories use) and `ase_source`.

**`Regulation` classes** ("cis": `|Ad| > 0.1`, i.e. allele fraction outside 0.5 ± 0.1; "trans": `|trans_log2| ≥ 1`):

| Class | Condition |
|-------|-----------|
| `Cis-only` | Parents differ; cis, no trans |
| `Trans-only` | Parents differ; trans, no cis |
| `Cis + Trans (Reinforcing)` | Parents differ; cis and trans, `cis_log2` and `trans_log2` of the same sign |
| `Cis x Trans (Opposing)` | Parents differ; cis and trans, opposite signs |
| `Compensatory` | Parents do not differ; alleles imbalanced |
| `Conserved` | Parents do not differ; alleles balanced |
| `Ambiguous` | Parents differ; neither component clears its threshold |

The thresholds are descriptive; the trans component is not tested. `Ad` and `Ed` are deviations on a ratio scale, not log fold changes, and `Ed` measures dominance of expression rather than the trans effect.

Two cases to keep in mind when reading the classes:
- A gene that is in the expression matrix but not in the DE table (one removed by EdgeRDE's expression filter, for example) has `parents_differ` false, so it can only be `Compensatory` or `Conserved`.
- In a summary table, a gene listed with no informative SNP and a ratio of 0.5 has `Ad = 0` and is treated as having balanced alleles.

**Example** (selected columns):
```
gene_id    parent1_ratio   parents_differ   Ad      Ed      cis_log2   trans_log2   Regulation
GENE0001   0.78            TRUE             -0.28   -0.01   1.83       0.15         Cis-only
GENE0002   0.51            TRUE             -0.01   -0.39   0.06       1.79         Trans-only
GENE0003   0.08            TRUE             0.42    0.46    -3.52      4.77         Cis x Trans (Opposing)
GENE0004   0.65            FALSE            -0.15           0.89       -0.82        Compensatory
GENE0005   0.49            FALSE            0.01            -0.06      -0.07        Conserved
```

#### `regulatory_divergence_regressions.pdf` and `.png`

**Description**: Advanced mode only. `Ad` (x) against `Ed` (y) for parentally DE genes in a single panel, coloured by `Regulation`, with one linear regression line and confidence band per class and dotted lines at zero. `Ambiguous` genes and genes with `Ed` outside -1.5 to 1.5 are left out; `Compensatory` and `Conserved` genes have no `Ed` and do not appear. The file is written only when at least one gene remains.

**Visual elements**:
- **Colours**: `Cis + Trans (Reinforcing)` red, `Cis x Trans (Opposing)` blue, `Trans-only` grey, `Cis-only` green
- **Legend**: each class with its number of genes and the Pearson correlation (r, p) between `Ad` and `Ed` within the class
- **Axes**: negative `Ad` = bias toward the parent 1 allele; negative `Ed` = F1 expression is parent 1-like

The colour is the cis/trans class, which is assigned from `Ad` and `trans_log2`. The position on the y-axis (`Ed`) is not what separates the classes.

#### `session_info.txt`

**Description**: R `sessionInfo()` for the run.

### Summary of ASEIntegrate Output Files

| File | Type | Count | Purpose |
|------|------|-------|---------|
| `{comparison}_with_ASE.tsv` | Data | 1 per classified comparison | Joined DE + ASE with category |
| `{comparison}_ASE_category_summary.tsv` | Data | 1 per classified comparison | Genes per category |
| `{comparison}_ASE_vs_DE.pdf` / `.png` | Plot | 2 per classified comparison | Fold change vs allele fraction |
| `cis_regulatory_candidates_all.tsv` | Data | 0 or 1 | Genes labelled `cis` |
| `ase_vs_de_ad_ed_data.tsv` | Data | 1 (advanced mode) | Cis/trans classes, Ad and Ed |
| `regulatory_divergence_regressions.pdf` / `.png` | Plot | 0 or 2 (advanced mode) | Ad vs Ed |
| `session_info.txt` | Log | 1 | R session |

Only the parent 1 versus parent 2 comparison is classified, so "per classified comparison" is normally once per run.

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
│   │   ├── Parent1_vs_Parent2_DE_results.tsv
│   │   ├── Parent1_vs_Parent2_volcano.pdf
│   │   ├── F1_hybrid_vs_Parent1_DE_results.tsv
│   │   ├── PCA_plot.pdf
│   │   ├── analysis_summary.txt
│   │   └── ...
│   ├── ASE_DE_results/
│   │   ├── Parent1_vs_Parent2_with_ASE.tsv
│   │   ├── Parent1_vs_Parent2_ASE_vs_DE.pdf
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

### Worked Reading (Illustrative Numbers)

The counts below are invented to show how the two sets of outputs fit together. They are not results from a real dataset. The groups are `Parent1`, `Parent2` and `F1_hybrid`.

1. **EdgeRDE** writes three comparisons: `F1_hybrid_vs_Parent1`, `F1_hybrid_vs_Parent2` and `Parent1_vs_Parent2`. Suppose `Parent1_vs_Parent2_significant_genes.tsv` lists 1,823 genes (FDR < 0.05).

2. **ASEIntegrate** classifies `Parent1_vs_Parent2` only and skips the two parent-vs-F1 files. `Parent1_vs_Parent2_ASE_category_summary.tsv` covers the genes that have both a DE row and an ASE row, 2,865 of them here:
   - 234 `cis`, 156 `trans_only` and 89 `trans_compensatory`: the 479 DE genes that have ASE information
   - 45 `compensatory` and 2,341 `no_divergence`: genes that are not DE between the parents

3. **Reading it**:
   - The other 1,344 of the 1,823 DE genes have no row in the ASE table, for example because they carry no informative SNP. They are absent from every ASEIntegrate file, so nothing is said about cis versus trans for them. They are not `no_divergence`.
   - Proportions depend on the denominator, so state it. Over the 479 classified DE genes this is about 49% `cis`, 33% `trans_only` and 19% `trans_compensatory`.
   - `is_DE` uses only `FDR` below the ASEIntegrate `--fdr-threshold`. It agrees with `*_significant_genes.tsv` when both commands were run with the same threshold, and it ignores the `--logfc-threshold` that EdgeRDE uses for its `Significance` column.

4. **Next steps**:
   - Take candidate genes from `cis_regulatory_candidates_all.tsv`.
   - Run the advanced mode for the per-gene `cis_log2` and `trans_log2` values and the `Regulation` class.

---

## See Also

- [ASE + DE Integration Workflow](./ase_de_integration_workflow.md)
- [Input Format Specification](./input_format_specification.md)
- [Interpretation Guidelines](./interpretation_guidelines.md)
