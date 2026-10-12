# Interpretation Guidelines for ASE + DE Analysis

## Overview

This guide helps you interpret the results of ASEIntegrate and EdgeRDE analyses, understand what the categories mean biologically, and troubleshoot common interpretive issues.

Column names, categories and thresholds are those written by `ASEIntegrate`; the exact definitions are in the [Output Files Reference](./output_files_reference.md#aseintegrate-outputs). The biological readings are possibilities to weigh against other evidence, not conclusions the pipeline reaches for you.

---

## Section 1: Understanding Gene Categories

### Before You Read the Categories

- **One comparison.** ASEIntegrate classifies the parent 1 versus parent 2 comparison only. Wherever this guide says "DE", it means differential expression between the two parents, read from `logFC_parent1_vs_parent2` (positive = higher in parent 1). Parent-vs-F1 comparisons are skipped and have no cis/trans reading.
- **"Biased" depends on the ASE table.** With an ASErtain gene table, a gene is biased when ASErtain called it (`ase_call`), which is a replicate-aware test. With the older summary table, a gene is biased when `bias_ratio > 0.5`, a majority vote over SNPs and not a test. The examples below use the summary table and show `bias_ratio`; with an ASErtain table that column holds `NA` and the call is in `is_biased`.
- **The examples are illustrative.** Gene IDs and numbers are made up to show how a row is read. They are not results from a real dataset.

### The Five Categories Explained

When ASEIntegrate categorizes genes, it combines two observations:
1. **Is the gene's expression different between the two parents?** (`is_DE`: `FDR` below `--fdr-threshold`, 0.05 by default; the size of the fold change is not used)
2. **Is one parental allele preferentially expressed in the F1?** (`is_biased`)

These two signals combine to suggest the *type* of regulatory change. The result is written to the `ase_de_category` column.

---

#### Category 1: `cis`

**Definition**: The gene is DE between the parents AND allelically biased in the F1, and the bias points to the more highly expressed parent.

**Key characteristics**:
- `is_DE` is true
- `is_biased` is true
- Bias toward parent 1 with a positive `logFC_parent1_vs_parent2`, or toward parent 2 with a negative one

**Biological interpretation**:
- In the F1 both alleles sit in the same nucleus and see the same trans factors. An allele that stays ahead there, in the same direction as the difference between the parents, points to a regulatory difference linked to the gene itself.
- Candidate causes include differences in promoter or enhancer sequence, methylation, or structural variants in regulatory regions. The category does not say which.

**Example** (illustrative):
```
Gene: GENE0001
- logFC_parent1_vs_parent2 = +2.15, FDR = 0.001 (parent 1 higher than parent 2)
- bias_ratio = 0.75 (75% of diagnostic SNPs show bias), predominant_bias = parent1
- parent1_ratio = 0.78 (78% of F1 reads come from the parent 1 allele)
- Category: cis. The parent 1 copy is expressed more both between the
  parents and within the F1.
```

**What to do next**:
1. **Sequence comparison**: Align promoter/enhancer regions between parents
2. **Motif analysis**: Search for differing transcription factor binding sites
3. **Epigenetics**: Check methylation/histone marks (if data available)
4. **Validation**: Real-time qPCR with allele-specific primers
5. **Functional assays**: Measure enhancer activity in both directions (if lab work planned)

**Tips for choosing cis candidates** (rules of thumb, not applied by ASEIntegrate):
- More SNPs behind a gene (`snp_count`) make the allelic estimate more reliable
- With an ASErtain table, `ase_q_value` holds ASErtain's q-value for the gene
- Genes with a large `logFC_parent1_vs_parent2` and a clear allelic ratio are easier to validate

---

#### Category 2: `trans_only`

**Definition**: The gene is DE between the parents but is NOT allelically biased in the F1.

**Key characteristics**:
- `is_DE` is true
- `is_biased` is false (summary table: `bias_ratio` of 0.5 or less)
- The allelic ratio is usually close to 50:50 (e.g., 0.51 vs 0.49), but the rule looks at the bias call, not at the ratio

**Biological interpretation**:
- The parents differ, yet the two alleles are expressed alike once they share one nucleus. The difference is then attributed to **trans-acting factors**, such as:
  - Diffusible transcription factors (TFs)
  - Chromatin modifiers
  - RNA-binding proteins
  - micro-RNAs
  - Signaling molecules
- These factors act on both parental alleles equally
- The two parents likely differ in the genes encoding these trans factors

**Caution**: "not biased" also covers "bias could not be shown". A gene with few SNPs or low coverage, and a summary-table gene listed with `snp_count = 0`, lands here when it is DE. `trans_only` therefore includes genes where a cis difference may exist but was not detected.

**Example** (illustrative):
```
Gene: GENE0002
- logFC_parent1_vs_parent2 = +1.89, FDR = 0.008 (parent 1 higher than parent 2)
- bias_ratio = 0.50, parent1_ratio = 0.51 (alleles close to 50:50 in the F1)
- Category: trans_only. The parents differ, but the two alleles are
  expressed alike in the F1.
```

**What to do next**:
1. **Causative trans factor**: Which gene encodes the regulator?
   - Check which known regulators are themselves DE between the parents
   - ChIP-seq if TF antibodies available
2. **Pathway analysis**: Do genes in this category share functional terms?
   - May point to a common TF or pathway
3. **Regulatory network**: Build networks of trans regulators and targets

**Challenge with trans-regulated genes**:
- Hard to map causative locus (could be on different chromosome)
- Usually requires additional experiments (ChIP, ATAC-seq, etc.)

---

#### Category 3: `trans_compensatory`

**Definition**: The gene is DE between the parents and allelically biased, but the bias points to the LESS highly expressed parent.

**Key characteristics**:
- `is_DE` is true
- `is_biased` is true
- Bias toward parent 1 with a negative `logFC_parent1_vs_parent2`, or toward parent 2 with a positive one

**Biological interpretation**:
- Two differences act in opposite directions. Within the F1, one allele is the more active one (the cis difference). Between the parents, the other parent is higher, so something outside the gene (trans) pushes the other way and outweighs the cis difference.
- In the advanced mode such a gene is labelled `Cis x Trans (Opposing)` when both components clear their thresholds.

**Example** (illustrative):
```
Gene: GENE0003
- logFC_parent1_vs_parent2 = +1.34, FDR = 0.042 (parent 1 higher than parent 2)
- bias_ratio = 0.92, predominant_bias = parent2
- parent1_ratio = 0.08 (F1 reads come mostly from the parent 2 allele)
- Category: trans_compensatory

Reading:
1. CIS: within the F1 the parent 2 allele is the more active one.
2. TRANS: between the parents, parent 1 is higher. A trans difference
   favours parent 1 and outweighs the cis difference.
```

**What to do next**:
1. **Look for the trans factor**: regulators that are DE between the parents are candidates
2. **Check the allelic estimate**: an opposing pattern resting on one or two SNPs deserves a second look before it is interpreted
3. **Phenotypic consequence**: does the phenotype follow the parental (trans-dominated) expression level?

---

#### Category 4: `compensatory`

**Definition**: The gene is NOT DE between the parents, but is allelically biased in the F1.

**Key characteristics**:
- `is_DE` is false (FDR at or above the threshold)
- `is_biased` is true

**Biological interpretation**:
- The alleles differ in the F1 although the parents do not differ overall. One explanation is a cis difference that is offset by a trans difference of opposite sign, leaving no net difference between the parents.
- "Not DE" is absence of evidence. A gene that just misses the FDR threshold lands here too, so this category is not proof of exact compensation.
- The basic analysis does not test the two components. `cis_log2` and `trans_log2` describe them for each gene.

**Example** (illustrative):
```
Gene: GENE0004
- logFC_parent1_vs_parent2 = +0.12, FDR = 0.89 (parents not different)
- bias_ratio = 0.67, predominant_bias = parent1, parent1_ratio = 0.65
- cis_log2 = 0.89, trans_log2 = -0.77
- Category: compensatory. The parent 1 allele is favoured in the F1,
  with no matching difference between the parents.
```

**What to do next**:
1. **Importance depends on context**:
   - If the phenotype does not differ between the parents, the allelic difference has no net effect on expression level
2. **Allele-specific RT-qPCR**: Confirm the allelic difference

---

#### Category 5: `no_divergence`

**Definition**: The gene is not DE between the parents and not allelically biased.

**Key characteristics**:
- `is_DE` is false
- `is_biased` is false

**Biological interpretation**:
- No regulatory divergence was detected at this gene: expression is similar between the parents and neither allele is favoured in the F1

**Two scenarios**:
1. **No divergence**: the gene is regulated alike in both parents
2. **No evidence**: too few SNPs or too little coverage to detect a bias. A summary-table gene listed with `snp_count = 0` is classified as unbiased, so it lands here when it is not DE

Genes with no row at all in the ASE table are a third case: they are not in the output and have no category.

**What to do next**:
- Usually treated as background
- Use as negative controls in validation experiments, choosing genes with enough SNPs for the unbiased call to mean something

---

## Section 2: Reading ASEIntegrate Output Files

### Quick Lookup: Key Metrics

When you open `*_with_ASE.tsv`, what matters?

| Column | What It Tells You | How to Use |
|--------|-------------------|-----------|
| `logFC_parent1_vs_parent2` | Direction + magnitude of the parental difference | Positive = parent 1 higher. `logFC` holds the same number in the orientation of the file name |
| `FDR` | Statistical confidence in DE | Below `--fdr-threshold` = DE (`is_DE`) |
| `snp_count` | Number of SNPs behind the ASE estimate | More SNPs = more robust; one or two = uncertain |
| `bias_ratio` | Summary table only: fraction of SNPs showing bias | Above 0.5 = biased; `NA` with an ASErtain table |
| `is_biased` | The bias call the category uses | ASErtain's `ase_call`, or `bias_ratio > 0.5` |
| `ase_de_category` | Regulatory category | See Section 1 |
| `<parent1_label>_ratio` | Fraction of F1 reads from the parent 1 allele | 0.5 = balanced; near 0 or 1 = strong imbalance |
| `cis_log2`, `trans_log2` | Descriptive sizes of the cis and trans components | `cis_log2` = log2 allelic ratio; `trans_log2` = parental difference minus `cis_log2` |

ASEIntegrate itself applies only the FDR threshold and the bias call. The remarks about SNP numbers are aids for reading the table.

### Example: Parsing a Single Gene Row

Selected columns, illustrative values, default labels:

```
TranscriptID: GENE0010
logFC: 2.45
logFC_parent1_vs_parent2: 2.45
FDR: 2.3e-05
snp_count: 7
bias_ratio: 0.86
predominant_bias: parent1
parent1_ratio: 0.82
parent2_ratio: 0.18
is_DE: TRUE
is_biased: TRUE
ase_de_category: cis
cis_log2: 2.19
trans_log2: 0.26
```

**What this means**:
1. **GENE0010** is about 5.5-fold higher in parent 1 than in parent 2 (2.45 log2 units)
2. **Significant** (FDR = 0.000023)
3. **Category is cis**: the F1 allelic bias points the same way as the parental difference
4. **ASE estimate**: 7 diagnostic SNPs, 86% of them showing bias
5. **Parent 1-biased**: the parent 1 allele supplies 82% of the F1 reads
6. **Sizes**: the allelic ratio (`cis_log2` = 2.19) accounts for most of the parental difference; the remainder (`trans_log2` = 0.26) is small

**Bottom line**: a cis candidate worth following up.

---

### Example: A Row That Carries No Information

```
TranscriptID: GENE0020
logFC_parent1_vs_parent2: 0.34
FDR: 0.67
snp_count: 0
bias_ratio: 0.00
predominant_bias: parent1
parent1_ratio: 0.50
parent2_ratio: 0.50
is_DE: FALSE
is_biased: FALSE
ase_de_category: no_divergence
```

**What this means**:
1. **Not DE** between the parents
2. **No SNPs**, so the "unbiased" call says nothing about the alleles
3. The label is `no_divergence` only because the summary table listed the gene with zero SNPs; had it been DE, the label would be `trans_only`

**Bottom line**: do not count rows like this as evidence for either category. Leaving such genes out of the ASE table avoids the problem.

---

## Section 3: Making Figures and Tables for Papers

The snippets below read the classified comparison's table. Replace `Parent1_vs_Parent2` with the name of your parental comparison, and `parent1_ratio` with `<parent1_label>_ratio` if you use custom labels.

### Table 1: Top Cis Candidates

```R
x <- read.delim("Parent1_vs_Parent2_with_ASE.tsv")
top <- subset(x, ase_de_category == "cis" & FDR < 0.01 & snp_count >= 4 &
                 abs(logFC_parent1_vs_parent2) > 1)
top <- top[order(top$FDR), c("TranscriptID", "logFC_parent1_vs_parent2", "FDR",
                             "snp_count", "cis_log2", "trans_log2")]
head(top, 15)
```

This gives you genes that:
- Are in the `cis` category
- Have at least 4 SNPs behind the ASE estimate
- Show a large parental fold change (easier to validate)

The cutoffs here (0.01, 4 SNPs, 1 log2 unit) are choices for this table, not values used by ASEIntegrate.

---

### Figure 1: Volcano Plot + Category Coloring

```R
palette_cat <- c(cis = "red", trans_only = "blue", trans_compensatory = "darkgreen",
                 compensatory = "orange", no_divergence = "grey60")

plot(x$logFC_parent1_vs_parent2, -log10(x$FDR),
     col = palette_cat[x$ase_de_category], pch = 19,
     xlab = "log2 fold change (parent 1 / parent 2)", ylab = "-log10(FDR)",
     main = "Parental DE coloured by ASE + DE category")
legend("topleft", legend = names(palette_cat), col = palette_cat, pch = 19, bty = "n")
```

Only genes with an ASE row are in the table, so this is not the full volcano plot of the comparison.

---

### Figure 2: Allele Fraction vs. Fold Change (Scatter)

ASEIntegrate already writes this plot as `{comparison}_ASE_vs_DE.pdf` / `.png`. To redraw it with your own styling:

```R
plot(x$logFC_parent1_vs_parent2, x$parent1_ratio,
     col = palette_cat[x$ase_de_category], pch = 19,
     xlab = "log2 fold change (parent 1 / parent 2)",
     ylab = "Parent 1 allele fraction in the F1")
abline(h = 0.5, v = 0, lty = 2, col = "grey50")  # Reference lines
```

Shows:
- **Upper right and lower left**: allelic bias and parental difference agree (`cis`)
- **Upper left and lower right**: they oppose (`trans_compensatory`)
- **Spread along x near y = 0.5**: parental difference without allelic imbalance (`trans_only`)

---

## Section 4: Troubleshooting Interpretation Issues

### Issue: "Most genes are in `no_divergence`; is analysis broken?"

**Not necessarily.** A large `no_divergence` class is expected when most genes are neither DE between the parents nor allelically biased.

The proportions depend on how diverged the parents are, on SNP density and read depth, on the thresholds, and on which ASE table was used (a tested call such as ASErtain's is stricter than `bias_ratio > 0.5`). This guide gives no reference percentages; compare categories within your own experiment.

---

### Issue: "All DE genes are `trans_only`; no cis effects detected"

**Possible causes**:
1. **Few genes are called biased**: with little SNP coverage or depth, bias is rarely detected and DE genes default to `trans_only`
   - Check: how many genes have `is_biased` true? What does `snp_count` look like across genes?
2. **Parents very similar**: few SNPs separate the alleles, so few genes can be assessed at all

**What to do**:
- Check the `snp_count` distribution across genes
- Count the genes with `is_biased` true, whatever their category
- Look at `cis_regulatory_candidates_all.tsv`: the file is only written when at least one gene is `cis`

---

### Issue: "`cis` and `trans_compensatory` split exactly by the sign of the fold change"

**Symptom**: every biased DE gene with a positive `logFC_parent1_vs_parent2` is `trans_compensatory` and every one with a negative value is `cis`.

**Possible cause**: the parent 1 label does not match the ASE table. The direction of bias is read by comparing `predominant_bias` with `--parent1-label`, and every other value counts as parent 2. If the values in the table are spelled differently from the label, every biased gene is read as biased toward parent 2.

**What to do**:
- Compare the values in `predominant_bias` (ASErtain: `direction`) with the two labels you passed; they are case-sensitive
- Check `bias_p1`: it should be true for some of the biased genes

---

### Issue: "A gene I expected is missing from the output"

**Possible reasons**:
1. **No ASE row**: only genes present in both the DE table and the ASE table are written
2. **Not in the DE table**: EdgeRDE filters out weakly expressed genes before testing
3. **Gene ID mismatch**: the ASE table must use the gene-level ID; `.T[digit]` suffixes are removed from the DE side only, and the match is case-sensitive
4. **Merged paralog**: a gene folded into a `PARA_*` feature has no DE row of its own

---

### Issue: "Results don't match my hypothesis about gene X"

**Possible reasons**:
1. **Low SNP count**: `snp_count` of 1 or 2 makes the ASE estimate noisy
2. **Different tissue**: cis and trans effects can differ between tissues
3. **Confounding factors**: if the parents differ in growth conditions as well as genotype, part of the parental difference is environmental and is counted as trans

**What to do**:
- Check the gene's row directly in the output file (search by ID)
- Look at `snp_count` and the allelic ratio for that gene
- Do independent validation (allele-specific qRT-PCR)

---

## Section 5: Advanced Analysis (cis/trans classes, Ad and Ed)

If you ran the advanced mode (`--expression-file`, `--metadata-file`, `--hybrid-group` and the two parent groups), `ase_vs_de_ad_ed_data.tsv` holds two different kinds of number. Keep them apart.

### The cis and trans components

- `cis_log2`: log2 of the F1 allelic ratio, parent 1 over parent 2. Both alleles share one nucleus in the F1, so this is the cis effect.
- `parental_log2`: log2 of the parental expression ratio. This is cis plus trans.
- `trans_log2 = parental_log2 - cis_log2`.

`Regulation` classifies each gene from these (`Cis-only`, `Trans-only`, `Cis + Trans (Reinforcing)`, `Cis x Trans (Opposing)`, `Compensatory`, `Conserved`, `Ambiguous`). "Cis" means the parent 1 allele fraction is outside 0.5 ± 0.1, "trans" means `|trans_log2| ≥ 1`, and "parents differ" means the gene is DE in the parental comparison (`parents_differ`). The thresholds are descriptive and the trans component is not tested, so treat a class as a description of where the gene sits, not as a significance call.

This classification is separate from `ase_de_category` of the basic analysis. It does not use the bias call at all: the cis side comes from the allele fraction alone. A gene can therefore be `trans_only` in the basic table and `Cis-only` here, or the reverse. A gene that is in the expression matrix but absent from the DE table counts as "parents do not differ" and can only be `Compensatory` or `Conserved`.

### Ad and Ed

- `Ad = 0.5 -` parent 1 allele fraction: the allelic imbalance on the fraction scale.
- `Ed = 0.5 - (F1 - P2) / (P1 - P2)`: where the F1's *total* expression sits between the parents. Zero is mid-parent (additive); negative is parent 1-like.

`Ed` is about dominance of expression, not about trans. They come apart in both directions: a purely additive trans difference gives `Ed = 0`, and an F1 can resemble the parent whose allele it expresses less. `Ed` is `NA` for genes whose parents are not differentially expressed, because there it would be a ratio of noise.

### Example readings

Illustrative values:

```
cis_log2 = 2.0, trans_log2 = 0.0, Ed = 0.0   -> Cis-only, additive inheritance
cis_log2 = 0.0, trans_log2 = 1.5, Ed = -0.4  -> Trans-only, F1 parent 1-like (dominant)
cis_log2 = 1.5, trans_log2 = 1.2             -> Cis + Trans (Reinforcing)
cis_log2 = 1.9, trans_log2 = -3.7            -> Cis x Trans (Opposing)
cis_log2 = 1.8, parents not different        -> Compensatory
```

---

## Section 6: Common Biological Interpretations

The scenarios below are hypotheses that a pattern of category counts is compatible with. The counts alone do not establish them: the proportions also depend on SNP density, read depth, the thresholds and the ASE table used.

### Scenario 1: Recent Speciation

**Pattern**:
- High proportion of `trans_only` genes
- Few `cis` genes
- Few `trans_compensatory`

**Interpretation**:
- Newly diverged species with mostly trans-acting changes
- Regulatory divergence is evolving rapidly
- Cis changes might be deleterious and purged

---

### Scenario 2: Ancient Divergence

**Pattern**:
- High proportion of `cis` genes
- Lower proportion of `trans_only`
- More `trans_compensatory` and `compensatory`

**Interpretation**:
- Older divergence allowed accumulation of cis mutations
- Many trans changes might have been deleterious (lost)
- Complex compensatory relationships evolved

---

### Scenario 3: Tissue-Specific Specialization

**Pattern**:
- Differences between tissues:
  - Tissue A: Mostly `cis` genes
  - Tissue B: Mostly `trans_only` genes

**Interpretation**:
- Different selection pressures in different tissues
- Cis regulation allows tissue-specific fine-tuning
- Trans regulation for coordinated global response

---

## Section 7: Downstream Analyses

### After you've categorized genes, what's next?

1. **Functional Enrichment** (cis genes):
   - Do cis genes share GO terms?
   - Example: Maybe all cis genes are "defense response"
   - Interpretation: This function has undergone regulatory divergence

2. **Phylogenetic Inference** (cis vs. trans):
   - Do cis changes accumulate at same rate as synonymous mutations?
   - If cis rate > syn rate: Under positive selection
   - If cis rate < syn rate: Under negative selection

3. **QTL Mapping**:
   - Use cis candidates as "targeted" markers
   - Co-localization with cis gene might indicate regulatory QTL
   - Trans genes less useful for QTL (effect is polygenic)

4. **Validation Experiments**:
   - Allele-specific RT-qPCR (small set of top genes)
   - ChIP-seq for trans-acting TFs
   - CRISPR knockouts for trans factors
   - Luciferase assays for promoter function

---

## Section 8: Literature References

### Key Papers on Cis vs. Trans Regulation

- **Wittkopp, Haerum & Clark (2004)** - Cis/trans decomposition from allelic expression in F1 hybrids
  - "Evolutionary changes in cis and trans gene regulation"
  - Nature 430:85-88

- **Tirosh & Barkai (2008)** - Yeast evolution
  - Shows ~30% cis, ~70% trans in transcriptional divergence
  - Molecular Systems Biology

- **McManus et al. (2010)** - Cis/trans classes from mRNA-seq of parents and F1 hybrids
  - "Regulatory divergence in Drosophila revealed by mRNA-seq"
  - Genome Research 20:816-825

- **Schaefke et al. (2015)** - Drosophila
  - Evidence for complex cis-trans interactions
  - Genome Biology

---

## See Also

- [Input Format Specification](./input_format_specification.md)
- [Output Files Reference](./output_files_reference.md)
- [ASE + DE Integration Workflow](./ase_de_integration_workflow.md)
