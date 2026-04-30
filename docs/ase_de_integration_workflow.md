# ASE + Differential Expression Integration Workflow

## Overview

This guide walks you through integrating allele-specific expression (ASE) analysis with differential expression (DE) analysis using SalmonStreamer's `EdgeRDE` and `ASEIntegrate` subcommands.

### When to Use This Workflow

- ✅ You have RNA-seq data from **hybrid crosses** (F1 crosses between two parental lines)
- ✅ You want to understand **which genes show cis vs. trans regulatory divergence**
- ✅ You have **diagnostic SNPs** that distinguish the two parental alleles
- ❌ You only have single-genome data (use `EdgeRDE` alone instead)

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
   AT1G01010   150        142        168
   AT1G01020   0          5          12
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

Create an ASE counts file with per-gene allele-specific expression information. This typically comes from:

- **EDGE_Penstemon pipeline**: `ase_read_counter.py` script produces this format
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

**Advanced analysis** (estimate cis and trans effects using Ad/Ed framework):

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
- Per-comparison: `*_with_ASE.tsv`, `*_ASE_category_summary.tsv`, `*_ASE_vs_DE.pdf`
- Global summary: `cis_regulatory_candidates_all.tsv`
- Advanced mode: `ase_vs_de_ad_ed_data.tsv`, regulatory divergence plots

---

## Interpreting Results

### Gene Categorization (Basic Analysis)

After ASEIntegrate runs, each gene is classified into one of these categories based on DE and ASE patterns:

| Category | DE Pattern | ASE Pattern | Biological Interpretation |
|----------|-----------|------------|--------------------------|
| **cis** | Yes (FDR < 0.05) | Biased | Regulatory changes affecting one parental allele's expression |
| **trans** | Yes (FDR < 0.05) | Unbiased | Trans factors (transcription factors, chromatin) affect both alleles equally |
| **trans_compensatory** | Yes (FDR < 0.05) | Biased opposite to DE | Trans factors act to partially compensate for cis differences |
| **compensatory** | No (FDR ≥ 0.05) | Biased | Cis differences balanced by trans factors; net expression unchanged |
| **no_divergence** | No (FDR ≥ 0.05) | Unbiased | No regulatory divergence detected |

### Understanding Output Files

**Per-comparison summary** (`*_ASE_category_summary.tsv`):
```
comparison         n_cis   n_trans   n_trans_comp   n_compensatory   n_no_div
P1_vs_F1           234     156       89             45               2341
P2_vs_F1           198     201       76             52               2288
```

**Detailed results** (`*_with_ASE.tsv`):
```
TranscriptID   logFC   FDR      gene_category          snp_count   bias_ratio   predominant_bias
AT1G01010      2.15    0.001    cis                    8           0.75         parent1
AT1G01020      1.89    0.008    trans                  4           0.50         parent1
AT1G01030     -1.34    0.042    trans_compensatory    12           0.92         parent2
```

### Cis vs. Trans Regulatory Effects

**Cis-regulatory** genes have:
- Significant allele bias (high `bias_ratio`)
- The biased allele typically comes from the species/line showing higher overall expression
- **Interpretation**: Regulatory sequences upstream/in the gene differ between parents

**Trans-regulatory** genes have:
- Little to no allele bias (low `bias_ratio`)
- Significant expression difference between hybrids and parents
- **Interpretation**: Diffusible factors (TFs, chromatin modifiers) regulate both alleles similarly

**Trans-compensatory** genes:
- High allele bias, but opposite to the direction of DE
- Example: Parent1 has more expression, but reads map to Parent2 allele in F1
- **Interpretation**: Cis differences exist, but trans factors overcorrect expression levels

---

## Advanced Analysis: Ad/Ed Framework

When you provide expression data, metadata, and group labels, ASEIntegrate runs an advanced analysis using the **Additivity/Dominance (Ad/Ed) framework** from Wittkopp et al.

This estimates:
- **Cis effect (ad)**: How much the regulatory divergence at this gene explains expression differences
- **Trans effect (ed)**: How much trans-acting factors explain expression differences
- **Classification**: Whether cis and trans effects reinforce or oppose each other

**Output file** (`ase_vs_de_ad_ed_data.tsv`):
```
TranscriptID   ad_estimate   ed_estimate   regulatory_class
AT1G01010      1.2           0.1           Cis_dominant
AT1G01020      0.3           1.8           Trans_dominant
AT1G01030      0.8           -0.7          Cis_trans_opposing
```

**Interpretation**:
- `ad > 0`: Cis effect acts in same direction as observed DE
- `ed > 0`: Trans effect acts in same direction as observed DE
- `ad + ed ≈ logFC`: The two effects roughly add up to the observed expression difference

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

Or remove all four arguments to skip advanced analysis.

---

## Real-World Example: Hybrid Pepper Study

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

Results show that 234 genes are cis-regulated, 156 are trans-regulated, and 89 show trans-compensatory effects in Bell Pepper vs F1 hybrids.

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
- [ASEIntegrate & EdgeRDE Audit Report](../AUDIT_ASEIntegrate_EdgeRDE.md)
