# ASE + Differential Expression Integration Workflow

## Overview

This guide walks you through integrating allele-specific expression (ASE) analysis with differential expression (DE) analysis using SalmonStreamer's `EdgeRDE` and `ASEIntegrate` subcommands.

### When to Use This Workflow

- ✅ You have RNA-seq data from **hybrid crosses** (F1 crosses between two parental lines)
- ✅ You want to understand **which genes show cis vs. trans regulatory divergence**
- ✅ You have **diagnostic SNPs** that distinguish the two parental alleles
- ❌ You have no F1 hybrid samples, or no SNPs that separate the parental alleles (use `EdgeRDE` alone instead)

**A second reference genome is not required.** What this workflow needs is F1 hybrid samples plus SNPs that distinguish the two parental alleles. A single reference genome with a VCF covering the parents and the F1s is sufficient. `ASEIntegrate` itself never reads a genome: it joins a per-gene ASE table to the `EdgeRDE` results by gene ID. Two genomes are only needed if you choose to quantify with the dual-genome Salmon route (README, "Option B"), which is a separate decision.

---

## Workflow Steps

### Step 1: Prepare Your Data

Before running either command, ensure you have:

1. **Count matrix file** (tab-separated)
   - First column: Gene IDs
   - Remaining columns: Raw read counts per sample
   - No special formatting needed
   ```
   gene_id     sample_1   sample_2   sample_3
   GENE0001    150        142        168
   GENE0002    0          5          12
   ```

2. **Metadata file** (tab-separated)
   - Required columns: `sample_name`, and either `group` OR `species`+`tissue`
   - Example:
   ```
   sample_name    group           species    tissue
   sample_1       F1_hybrid       F1         leaf
   sample_2       F1_hybrid       F1         leaf
   sample_3       Parent1         parent1    leaf
   ```

3. **ASE counts file** (tab-separated, created separately)
   - See [Input Format Specification](./input_format_specification.md)

### Step 2: Run Differential Expression Analysis (EdgeRDE)

First, run `EdgeRDE` to identify genes with significant expression differences between groups.

```bash
python SalmonStreamer.py EdgeRDE \
    --expression-file counts.tsv \
    --metadata-file metadata.tsv \
    --output-dir DE_results/ \
    --fdr-threshold 0.05 \
    --logfc-threshold 1.0
```

**Output**: `DE_results/` directory with per-comparison results
- `Group1_vs_Group2_DE_results.tsv` — All genes with logFC, p-value, FDR
- `Group1_vs_Group2_significant_genes.tsv` — Filtered to FDR < 0.05
- `Group1_vs_Group2_volcano.pdf` — Volcano plot
- Quality control plots: PCA, correlation heatmap, library sizes

**Expected output files**: See [Output Files Reference](./output_files_reference.md#edgerde-outputs)

**⚠️ Important Notes**:
- Gene IDs in the count matrix should match your ASE counts file exactly (or with `.T1`, `.T2` suffixes removed)
- Sample names in the metadata file must match count matrix column headers exactly
- If using `--sample-suffix`, provide a **regex pattern** (e.g., `"_R1_filtered$"` to strip `_R1_filtered` suffix)

---

### Step 3: Prepare ASE Counts File

`ASEIntegrate` does not compute allele-specific expression. It consumes a finished per-gene ASE table and joins it to the DE results by gene ID, so the allele counting has to be done beforehand with a separate tool. Sources for that table:

- **EDGE_Penstemon pipeline**: `ase_read_counter.py` script produces this format
- **[ASErtain](https://github.com/luismadrigal98/ASErtain)**: computes ASE from a single reference genome, a VCF covering the parents and the F1s, and the F1 BAM files. Its output is **not yet in this format** — see [Using ASErtain output](#using-asertain-output) below
- **Custom analysis**: See template in [Input Format Specification](./input_format_specification.md#ase-counts-file-format)

**Required format** (tab-separated):
```
gene_id     snp_count   bias_ratio   predominant_bias   parent1_ratio   parent2_ratio
AT1G01010   8           0.75         parent1            0.78            0.22
AT1G01020   4           0.50         parent1            0.51            0.49
AT1G01030   12          0.92         parent2            0.08            0.92
```

**Column descriptions**:
- `gene_id`: Must match gene IDs in DE results (will auto-strip `.T[digit]` suffixes)
- `snp_count`: Number of diagnostic SNPs in this gene
- `bias_ratio`: Fraction of SNPs showing allelic bias (0.0 to 1.0)
- `predominant_bias`: Which parent's allele dominates ("parent1" or "parent2")
- `parent1_ratio`: Fraction of reads from parent 1 allele
- `parent2_ratio`: Fraction of reads from parent 2 allele

#### Using ASErtain output (preferred)

Pass ASErtain's per-gene table (`<prefix>.gene_ase.tsv`) to `--ase-counts-file` as it is, and set `--parent1-label` / `--parent2-label` to the two lineage labels of the ASErtain config. `ASEIntegrate` recognises the table from its columns and reads it like this:

| ASEIntegrate needs | Taken from the ASErtain table |
|--------------------|-------------------------------|
| is the gene allelically biased? | `ase_call` — ASErtain's replicate-aware call, with its effect floor and QC gates |
| toward which parent? | `direction` |
| parent 1 allele fraction | `mean_<label>_ratio` (ASErtain reports its variable lineage; the other is one minus it) |
| SNP count | `n_snps` |

`q_value`, `log2_ratio` and `phase_discordance_cause` are carried into the output with an `ase_` prefix. This is the input to prefer: with the legacy table above, "biased" means `bias_ratio > 0.5`, a majority vote over SNPs that takes no account of depth or replication.

ASErtain can also do the join itself (`asertain contrast --parental-de <EdgeRDE table> --de-gene-col TranscriptID --de-log2-col logFC --de-padj-col FDR`). There the fold change must already be oriented variable lineage over fixed lineage, and gene IDs must match exactly — it does not strip `.T[digit]` suffixes.

---

### Step 4: Run ASE + DE Integration (ASEIntegrate)

Now integrate the DE results with ASE data to categorize regulatory effects.

**Basic analysis** (categorize genes as cis, trans, etc.):

```bash
python SalmonStreamer.py ASEIntegrate \
    --de-results-dir DE_results/ \
    --ase-counts-file ase_counts.tsv \
    --output-dir ASE_DE_results/ \
    --fdr-threshold 0.05 \
    --parent1-label parent1 \
    --parent2-label parent2
```

**Advanced analysis** (adds the Ad/Ed classification described [below](#advanced-analysis-aded-classification)):

```bash
python SalmonStreamer.py ASEIntegrate \
    --de-results-dir DE_results/ \
    --ase-counts-file ase_counts.tsv \
    --output-dir ASE_DE_results/ \
    --fdr-threshold 0.05 \
    --parent1-label parent1 \
    --parent2-label parent2 \
    --expression-file counts.tsv \
    --metadata-file metadata.tsv \
    --hybrid-group F1_hybrid \
    --parent1-group Parent1 \
    --parent2-group Parent2
```

**Output**: `ASE_DE_results/` directory with:
- Per comparison file: `*_with_ASE.tsv`, `*_ASE_category_summary.tsv`, `*_ASE_vs_DE.pdf` / `.png`
- Across comparisons: `cis_regulatory_candidates_all.tsv`
- Advanced mode: `ase_vs_de_ad_ed_data.tsv`, `regulatory_divergence_regressions.pdf` / `.png`

#### Which comparison files to read

Cis/trans classification is defined for one comparison only, parent 1 versus parent 2, and its rule needs a positive fold change to mean "higher in parent 1". `EdgeRDE` writes all pairwise comparisons — with two parents and an F1, three files, two of them parent-vs-F1 — names each `A_vs_B` with the group names in alphabetical order, and reports `logFC` as `A` minus `B`. So neither the right file nor its sign can be taken for granted, and `ASEIntegrate` handles both:

- **It classifies the parental comparison only.** Every other `*_DE_results.tsv` is skipped, with a line in the log saying so.
- **It orients the fold change.** When parent 2's group is the one named first, the sign is flipped. The output keeps the original `logFC` and adds `logFC_parent1_vs_parent2`, which is what the categories and the plot use.

Tell it which comparison is the parental one with `--parent1-group` and `--parent2-group` (the two group names exactly as they appear in the file name). Without them it looks for the parent labels inside the group names: a comparison qualifies when one group name contains `--parent1-label` only and the other `--parent2-label` only. If no file qualifies, nothing is classified and the run stops with an error rather than guessing.

**Note on paralog-aware DE**: if the count matrix came from the paralog-aware route (`ParalogMerge` / `ParalogTximport`), genes folded into a merged `PARA_*` feature have no DE row of their own. The feature ID matches no `gene_id` in a per-gene ASE table, so those rows drop out of the ASE + DE join, and out of the Ad/Ed table, without a warning. Genes left unmerged keep their own ID and join as usual. To see which feature a gene went into, look it up in the gene-to-feature map written by `ParalogMerge --out-map` (`paralog_features_map.tsv` in the README example; columns `gene`, `feature`, `n_genes_in_feature`, `merged`).

---

## Interpreting Results

### Gene Categorization (Basic Analysis)

In the parental comparison, every gene that has both a DE row and an ASE row is assigned one category, written to the `ase_de_category` column:

| Category | DE Pattern | ASE Pattern | Biological Interpretation |
|----------|-----------|------------|--------------------------|
| `cis` | Yes | Biased, same direction as DE | Regulatory changes affecting one parental allele's expression |
| `trans_only` | Yes | Unbiased | Trans factors (transcription factors, chromatin) affect both alleles equally |
| `trans_compensatory` | Yes | Biased, opposite to DE | Trans factors act to partially compensate for cis differences |
| `compensatory` | No | Biased | Cis differences balanced by trans factors; net expression unchanged |
| `no_divergence` | No | Unbiased | No regulatory divergence detected |

The rule behind the table:
- **DE** means `FDR` below `--fdr-threshold`. The size of `logFC` is not used.
- **Biased** means `bias_ratio > 0.5`.
- **Direction** compares `predominant_bias` with the sign of `logFC`. Bias toward `--parent1-label` with a positive `logFC`, or bias toward the other parent with a negative `logFC`, counts as the same direction.

These labels are meaningful only for the parent-vs-parent comparison in the right orientation; see [Which comparison files to read](#which-comparison-files-to-read).

### Understanding Output Files

The example values in this section are illustrative, not results from a real dataset.

**Per-comparison summary** (`*_ASE_category_summary.tsv`): one file per comparison, one row per category.
```
ase_de_category      n_genes   avg_logFC_parent1_vs_parent2   avg_parent1_ratio
cis                  234       0.41                           0.61
compensatory         45        0.03                           0.58
no_divergence        2341      0.01                           0.50
trans_compensatory   89        -0.12                          0.57
trans_only           156       0.35                           0.50
```
Both averages are plain means over the genes in the category: the oriented fold change and the parent 1 allele fraction.

**Detailed results** (`*_with_ASE.tsv`), selected columns:
```
TranscriptID   logFC   FDR      snp_count   bias_ratio   predominant_bias   ase_de_category
AT1G01010      2.15    0.001    8           0.75         parent1            cis
AT1G01020      1.89    0.008    4           0.50         parent1            trans_only
AT1G01030      1.34    0.042    12          0.92         parent2            trans_compensatory
```
The file keeps every column of the DE table and of the ASE table, and adds `gene_id_base` (the ID used for the join), `logFC_parent1_vs_parent2`, `is_DE`, `is_biased`, `bias_p1`, `ase_de_category`, and the two descriptive effect sizes `cis_log2` (log2 of the F1 allelic ratio, parent 1 over parent 2) and `trans_log2` (`logFC_parent1_vs_parent2` minus `cis_log2`). With an ASErtain table the `snp_count`, `bias_ratio` and `predominant_bias` columns hold `n_snps`, nothing, and `direction`. Genes without an ASE row are not written.

### Cis vs. Trans Regulatory Effects

**Cis-regulatory** genes (`cis`) have:
- Allelic bias (`bias_ratio > 0.5`) in the same direction as the expression difference
- The biased allele typically comes from the species/line showing higher overall expression
- **Interpretation**: Regulatory sequences upstream/in the gene differ between parents

**Trans-regulatory** genes (`trans_only`) have:
- No allelic bias (`bias_ratio` of 0.5 or less)
- Significant expression difference between the two parents
- **Interpretation**: Diffusible factors (TFs, chromatin modifiers) regulate both alleles similarly

**Trans-compensatory** genes (`trans_compensatory`):
- Allelic bias (`bias_ratio > 0.5`), but opposite to the direction of DE
- Example: Parent1 has more expression, but reads map to Parent2 allele in F1
- **Interpretation**: Cis differences exist, but trans factors overcorrect expression levels

---

## Advanced Analysis: Ad/Ed Classification

When you provide the expression matrix, the metadata and the three group labels, `ASEIntegrate` adds a classification that uses the parents' and the F1's own expression levels, and reports two deviations, **Ad** and **Ed**.

Expression is TMM-normalised to CPM and averaged within each group, giving `parent1_avg`, `parent2_avg` and `hybrid_avg`. Then, per gene:

- `cis_log2` = log2 of the F1 allelic ratio, parent 1 over parent 2.
- `parental_log2` = log2 of `parent1_avg` over `parent2_avg`.
- `trans_log2 = parental_log2 - cis_log2`.
- `parents_differ`: the gene is differentially expressed between the parents (FDR below the threshold in the parental comparison; a two-fold difference in mean CPM if that comparison is not among the DE results).

`Regulation` is the standard decomposition (Wittkopp et al. 2004; McManus et al. 2010), with "cis" meaning the parent 1 allele fraction is outside 0.5 ± 0.1 and "trans" meaning `|trans_log2| ≥ 1`:

| `Regulation` | Condition |
|--------------|-----------|
| `Cis-only` | Parents differ; cis, no trans |
| `Trans-only` | Parents differ; trans, no cis |
| `Cis + Trans (Reinforcing)` | Parents differ; cis and trans with the same sign |
| `Cis x Trans (Opposing)` | Parents differ; cis and trans with opposite signs |
| `Compensatory` | Parents do not differ, but the alleles are imbalanced |
| `Conserved` | Parents do not differ, alleles balanced |
| `Ambiguous` | Parents differ, neither component clears its threshold |

The thresholds are descriptive. There is no statistical test of the trans component.

**Ad and Ed** are reported in the same table and drawn in the plot. They are not the cis and trans effects:

- `Ad = 0.5 - <parent1_label>_ratio`: the departure of the F1 allelic ratio from 50:50. Negative means the parent 1 allele is over-represented.
- `position_ratio = (hybrid_avg - parent2_avg) / (parent1_avg - parent2_avg)` places the F1's *total* expression between the parents: 1 is the parent 1 level, 0 the parent 2 level.
- `Ed = 0.5 - position_ratio`: the departure of F1 expression from the mid-parent value. Negative means the F1 is parent 1-like.

`Ed` measures dominance of expression, which is a different question from cis versus trans: a purely additive trans difference leaves `Ed` at zero, and an F1 can resemble the parent whose allele it expresses less. (Earlier versions classified on the signs of `Ad` and `Ed`; that labelled a pure-cis, additive gene `Ambiguous` and could label opposing cis and trans effects `Reinforcing`.) Because `Ed` divides by the parental difference, it is computed only where `parents_differ` is true and left empty elsewhere.

**Output file** (`ase_vs_de_ad_ed_data.tsv`), selected columns, illustrative values:
```
gene_id     cis_log2   parental_log2   trans_log2   Ad      Ed      Regulation
AT1G01010   2.00       1.98            -0.02        -0.30   -0.01   Cis-only
AT1G01020   0.00       1.50            1.50         0.00    -0.44   Trans-only
AT1G01030   1.91       -1.81           -3.72        -0.29   -0.47   Cis x Trans (Opposing)
AT1G01040   1.83       0.09            -1.73        -0.28           Compensatory
```
The file holds every column of the ASE table plus `parent1_avg`, `parent2_avg`, `hybrid_avg`, `parental_log2`, `parental_FDR`, `parents_differ`, `Ad`, `parental_diff`, `position_ratio`, `Ed`, `cis_log2`, `trans_log2` and `Regulation`.

**Things to know before interpreting it**:
- `Ad` and `Ed` are deviations on a ratio scale, not log fold changes. They do not sum to `logFC`; the effect sizes on the log scale are `cis_log2` and `trans_log2`.
- The table is built once per run from the ASE table and the expression matrix. The DE results enter only through `parents_differ`.
- Gene IDs are matched after stripping any `.T[digit]` suffix from the expression matrix, as in the basic join.
- `regulatory_divergence_regressions.pdf` / `.png` plots `Ad` against `Ed` for parentally DE genes, coloured by `Regulation` with one regression line per class, leaving out `Ambiguous` genes and genes with `Ed` outside -1.5 to 1.5.
- The Ad/Ed block needs the Bioconductor package `edgeR` for TMM normalisation; the basic analysis does not.

---

## Common Issues & Troubleshooting

### Issue: "No matching samples found"

**Problem**: Sample names in metadata don't match count matrix column headers.

**Solution**:
1. Check exact spelling (case-sensitive)
2. Use `--sample-suffix` if headers have common suffixes:
   ```bash
   python SalmonStreamer.py EdgeRDE \
       --expression-file counts.tsv \
       --metadata-file metadata.tsv \
       --sample-suffix "_R1_filtered$" \
       --output-dir DE_results/
   ```

### Issue: "Gene ID mismatch" or zero genes returned

**Problem**: Gene IDs in count matrix don't match ASE counts file.

**Solution**:
1. Check if count matrix IDs have `.T1`, `.T2` suffixes (transcript variants)
2. ASE file should just have gene IDs (the script strips suffixes automatically)
3. Verify IDs match exactly (TAIR IDs: `AT1G01010` not `at1g01010`)

### Issue: "Column '___' not found in ASE file"

**Problem**: ASE counts file is missing required columns.

**Solution**:
1. Verify column names exactly match:
   - `gene_id`, `snp_count`, `bias_ratio`, `predominant_bias`
   - `{parent1_label}_ratio` (e.g., `parent1_ratio`)
   - `{parent2_label}_ratio` (e.g., `parent2_ratio`)
2. Make sure parent labels match what you specified (check spelling)

### Issue: Advanced mode says "requires ALL of..."

**Problem**: Trying to run Ad/Ed analysis but missing required files.

**Solution**: For advanced analysis, provide all of:
```bash
--expression-file counts.tsv \
--metadata-file metadata.tsv \
--hybrid-group F1_hybrid \
--parent1-group Parent1 \
--parent2-group Parent2
```

Or remove all five arguments to skip advanced analysis.

---

## Worked Example (Illustrative): Hybrid Pepper Study

The file names and labels below are illustrative. This is a command walkthrough, not a report of results from a real dataset.

```bash
# Step 1: EdgeRDE on count data from Bell Pepper × Chili Pepper cross
python SalmonStreamer.py EdgeRDE \
    --expression-file pepper_counts.tsv \
    --metadata-file pepper_metadata.tsv \
    --output-dir pepper_DE/ \
    --fdr-threshold 0.05 \
    --logfc-threshold 1.0

# Step 2: ASEIntegrate with basic classification
python SalmonStreamer.py ASEIntegrate \
    --de-results-dir pepper_DE/ \
    --ase-counts-file pepper_ase.tsv \
    --output-dir pepper_ASE_DE/ \
    --parent1-label BellPepper \
    --parent2-label ChiliPepper

# Step 3: Advanced analysis
python SalmonStreamer.py ASEIntegrate \
    --de-results-dir pepper_DE/ \
    --ase-counts-file pepper_ase.tsv \
    --output-dir pepper_ASE_DE_advanced/ \
    --parent1-label BellPepper \
    --parent2-label ChiliPepper \
    --expression-file pepper_counts.tsv \
    --metadata-file pepper_metadata.tsv \
    --hybrid-group F1 \
    --parent1-group BellPepper \
    --parent2-group ChiliPepper
```

`ASEIntegrate` writes outputs for the parent-vs-parent comparison only (`BellPepper_vs_ChiliPepper_*`) and skips the two parent-vs-F1 files; see [Which comparison files to read](#which-comparison-files-to-read).

---

## Next Steps

After integrating ASE + DE data:

1. **Filter candidates**: Focus on cis-regulatory candidates from `cis_regulatory_candidates_all.tsv`
2. **Validate findings**: 
   - Check which genes have high SNP counts (more robust bias estimates)
   - Look for consistency across multiple comparisons
3. **Functional annotation**: 
   - Identify regulatory motifs in biased alleles
   - Check expression breadth (tissue-specific vs. constitutive)
4. **QTL mapping**: Use these genes as candidates for downstream QTL studies

---

## See Also

- [Input Format Specification](./input_format_specification.md)
- [Output Files Reference](./output_files_reference.md)
- [Interpretation Guidelines](./interpretation_guidelines.md)
