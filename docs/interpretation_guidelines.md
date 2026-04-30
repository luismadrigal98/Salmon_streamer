# Interpretation Guidelines for ASE + DE Analysis

## Overview

This guide helps you interpret the results of ASEIntegrate and EdgeRDE analyses, understand what the categories mean biologically, and troubleshoot common interpretive issues.

---

## Section 1: Understanding Gene Categories

### The Five Categories Explained

When ASEIntegrate categorizes genes, it combines two observations:
1. **Is the gene's expression different between groups?** (DE analysis)
2. **Is one parental allele preferentially expressed?** (ASE analysis)

These two signals combine to infer the *type* of regulatory change.

---

#### Category 1: `cis`-regulated

**Definition**: Gene shows BOTH significant DE AND allelic bias.

**Key characteristics**:
- FDR < 0.05 (statistically different expression between groups)
- High `bias_ratio` (typically > 0.75)
- `bias_ratio` column shows clear winner (e.g., 0.78 parent1 vs 0.22 parent2)
- The biased allele usually matches the direction of DE (e.g., parent1 biased AND parent1 higher in F1)

**Biological interpretation**:
- **The gene's regulatory sequence differs between parental lines**
- This could be:
  - Mutations in the promoter
  - Different enhancer elements
  - Different methylation patterns
  - Structural variants in regulatory regions
- One parent's version of the gene is "more active" than the other's

**Example scenario** (Arabidopsis × related species):
```
Gene: AT1G01010
- logFC = +2.15 (Arabidopsis > Lyrata in F1)
- bias_ratio = 0.75 (75% of diagnostic SNPs show bias)
- parent1_ratio = 0.78 (78% reads from Arabidopsis allele)
- Interpretation: Arabidopsis has evolved a MORE ACTIVE promoter/enhancer
  Result: Even in F1 hybrids, Arabidopsis allele used preferentially
```

**What to do next**:
1. **Sequence comparison**: Align promoter/enhancer regions between parents
2. **Motif analysis**: Search for differing transcription factor binding sites
3. **Epigenetics**: Check methylation/histone marks (if data available)
4. **Validation**: Real-time qPCR with allele-specific primers
5. **Functional assays**: Measure enhancer activity in both directions (if lab work planned)

**Tips for finding cis-regulatory variants**:
- `snp_count` ≥ 4 makes bias estimates more reliable
- Focus on genes with consistent bias (same SNPs showing bias repeatedly)
- Look for genes with extreme logFC and bias_ratio (easier to validate)

---

#### Category 2: `trans`-regulated

**Definition**: Gene shows significant DE BUT little/no allelic bias.

**Key characteristics**:
- FDR < 0.05 (statistically different expression)
- Low `bias_ratio` (typically < 0.5)
- Ratios approximately 50:50 (e.g., 0.51 parent1 vs 0.49 parent2)
- Both alleles expressed in similar proportions (unbiased)

**Biological interpretation**:
- **The gene's regulatory change is NOT in the gene itself**
- Instead, expression is controlled by **trans-acting factors** (DNA), such as:
  - Diffusible transcription factors (TFs)
  - Chromatin modifiers
  - RNA-binding proteins
  - micro-RNAs
  - Signaling molecules
- These factors act on BOTH parental alleles equally
- The two parents likely differ in the genes encoding these trans factors

**Example scenario** (pepper cross):
```
Gene: Caper_001030
- logFC = +1.89 (BellPepper > ChiliPepper in F1)
- bias_ratio = 0.50 (50% of SNPs show bias - total noise)
- parent1_ratio = 0.51 (essentially 50:50)
- Interpretation: BellPepper produces a trans factor (maybe a TF?)
  that activates BOTH the BellPepper AND ChiliPepper versions of this gene.
  Result: In F1, this gene is upregulated even though both alleles
  are equally used.
```

**What to do next**:
1. **Causative trans factor**: Which gene encodes the regulator?
   - Run DE analysis on known TFs; see which correlates
   - RNA-seq of TF candidates in parents vs F1
   - ChIP-seq if TF antibodies available
2. **Pathway analysis**: Do genes in this category share functional terms?
   - May point to a common TF or pathway
3. **Regulatory network**: Build networks of trans regulators and targets

**Challenge with trans-regulated genes**:
- Hard to map causative locus (could be on different chromosome)
- Usually requires additional experiments (ChIP, ATAC-seq, etc.)

---

#### Category 3: `trans_compensatory`

**Definition**: Gene shows significant DE, but allele bias OPPOSES the DE direction.

**Key characteristics**:
- FDR < 0.05 (expression is different)
- High `bias_ratio` (clear allele preference)
- BUT the biased allele is from the parent with LOWER overall expression
- Example: F1 shows more expression from Parent2, but allele bias favors Parent1

**Biological interpretation**:
- **Complex regulatory divergence**: Two opposing changes have occurred
- **Cis change**: Gene has regulatory mutations that favor one allele (call this direction the "cis effect")
- **Trans change**: But trans factors from the other parent override this, flipping the outcome
- Result: Despite cis regulation, the final expression direction is opposite

**Example scenario** (realistic):
```
Gene: AT1G01040
Scenario: 
- Parent1 (Barbara) vs Parent2 (Lena) in F1 cross
- logFC = -1.34 in F1 (F1 lower than Barbara, biased toward Lena)
- BUT allele bias shows Barbara bias (0.92 Barbara, 0.08 Lena)

Interpretation:
1. CIS: Barbara has evolved a STRONGER promoter 
   (intrinsically more active)
2. TRANS: But Lena produces a trans factor that REPRESSES both alleles
   (e.g., repressive TF or chromatin modifier)
3. RESULT: In F1 hybrids, this repressive trans factor dominates,
   flipping expression direction. Despite Barbara's allele being
   used 92% of the time, overall F1 expression is low.

Biological meaning: Lena has evolved a trans-acting repressor
that overcomes Barbara's cis-regulatory divergence.
```

**What to do next**:
1. **Identify the repressive trans factor**:
   - Look for DE repressor genes in Lena background
   - Test with knockdown/overexpression if possible
2. **Understand the conflict**: Why would these opposing changes exist?
   - Lena might have lost a TF (compensating for evolved promoter activity)
   - Or Lena evolved a new repressor as an adaptive response
3. **Phenotypic consequence**: 
   - Does the F1 phenotype reflect the final (trans-dominated) expression?
   - Or do cis effects matter in other tissues/conditions?

**Rarity**: Trans-compensatory genes are usually less common (10-15% of DE genes) because it requires two specific changes to have occurred.

---

#### Category 4: `compensatory`

**Definition**: Gene shows NO significant DE, but still has allelic bias.

**Key characteristics**:
- FDR ≥ 0.05 (expression NOT significantly different)
- High `bias_ratio` (clear allele preference)
- But net expression level is unchanged

**Biological interpretation**:
- **Cis difference is present, but masked by trans**
- One allele is genuinely more active (cis)
- BUT the less-active allele is compensated by trans factors
- Result: No detectable net difference in expression
- Example: Cis effect = +1.5 logFC, but trans effect = -1.5 logFC, sum = 0

**Example scenario**:
```
Gene: AT1G05000
- logFC = +0.12 (NOT significant; FDR = 0.89)
- bias_ratio = 0.67 (67% of SNPs favor parent1)
- parent1_ratio = 0.65

Interpretation:
- Parent1 allele is genuinely more active (allele bias shows this)
- BUT parent2 produces trans factors that specifically upregulate
  the parent2 allele
- Result: Both alleles contribute similarly to total expression
- Net effect: No DE detected

Biological meaning: Might be an example of "balancing selection" where
natural selection favors trans-acting compensatory mechanisms.
```

**What to do next**:
1. **Importance depends on context**:
   - If phenotype doesn't differ → compensatory is "working as designed"
   - If phenotype differs → may indicate incomplete compensation
2. **Candidate genes for stabilizing selection**: 
   - These genes show evidence of trans-acting compensation
   - May be under selection pressure for dosage balance
3. **Allele-specific RT-qPCR**: Validate cis effect is real despite no net DE

---

#### Category 5: `no_divergence`

**Definition**: No significant DE and no allelic bias.

**Key characteristics**:
- FDR ≥ 0.05 (expression NOT significantly different)
- Low `bias_ratio` (no allele preference)
- 50:50 allele ratio (or no diagnostic SNPs, `snp_count = 0`)

**Biological interpretation**:
- **No regulatory divergence between parents at this gene**
- Expression is similar between groups
- Both alleles used in similar proportions
- Most of the genome falls into this category (no divergence)

**Two scenarios**:
1. **True lack of divergence**: Gene identical between parents, no selection/mutation
2. **Lack of diagnostic SNPs** (`snp_count = 0`): Cannot determine if bias exists (insufficient data)

**What to do next**:
- Usually ignored (they're the background)
- Use as negative controls in validation experiments
- In QTL studies, these genes are unlikely to show linkage

---

## Section 2: Reading ASEIntegrate Output Files

### Quick Lookup: Key Metrics

When you open `*_with_ASE.tsv`, what matters?

| Metric | What It Tells You | How to Use |
|--------|-------------------|-----------|
| `logFC` | Direction + magnitude of expression difference | Positive = group1 higher; >2 or <-2 is "large" |
| `FDR` | Statistical confidence in DE | <0.05 usually means "believe it"; <0.01 more confident |
| `snp_count` | Reliability of ASE estimate | ≥4 means robust; <2 means uncertain |
| `bias_ratio` | Consistency of allele bias | >0.75 means very strong; 0.5 means no bias |
| `gene_category` | Regulatory mechanism | Tells you cis vs. trans |
| `parent1_ratio` | Which allele dominates | >0.75 = strong parent1 bias; ≈0.5 = no bias |

### Example: Parsing a Single Gene Row

```
TranscriptID: AT1G15234
logFC: 2.45
FDR: 2.3e-05
gene_category: cis
snp_count: 7
bias_ratio: 0.86
predominant_bias: parent1
parent1_ratio: 0.82
parent2_ratio: 0.18
```

**What this means**:
1. **AT1G15234** is upregulated ~5.5-fold (2.45 logFC) in group 1
2. **Very significant** (FDR = 0.000023; essentially certain)
3. **Mechanism is cis**: The gene's regulatory region evolved
4. **Robust ASE call**: 7 diagnostic SNPs, 86% show same bias pattern
5. **Parent1-biased**: Parent1's allele used 82% of the time
6. **Interpretation**: Parent1 has the more active version of this gene's promoter

**Bottom line**: High-confidence cis-regulated gene, worth investigating.

---

### Example: Genes to Deprioritize

```
TranscriptID: AT1G89456
logFC: 0.34
FDR: 0.67
gene_category: no_divergence
snp_count: 0
bias_ratio: 0.00
predominant_bias: parent1
parent1_ratio: 0.50
parent2_ratio: 0.50
```

**What this means**:
1. **Not significantly different** between groups
2. **No allele bias** (and no SNPs to measure it)
3. **Not diverged** between parents
4. **No interest for this analysis**

**Bottom line**: This is background; move on.

---

## Section 3: Making Figures and Tables for Papers

### Table 1: Top Cis-Regulated Genes

For your manuscript:

```sql
SELECT TranscriptID, logFC, FDR, bias_ratio, snp_count 
FROM *_with_ASE.tsv
WHERE gene_category = "cis" 
  AND FDR < 0.01
  AND snp_count >= 4
  AND abs(logFC) > 1.0
ORDER BY FDR
LIMIT 15;
```

This gives you:
- High-confidence genes
- Robust ASE estimates (≥4 SNPs)
- Large fold-changes (easier to validate)

---

### Figure 1: Volcano Plot + Category Coloring

```R
# Pseudocode for enhanced volcano plot
ase_data <- read.table("comparison_with_ASE.tsv")

colors <- ifelse(ase_data$gene_category == "cis", "red",
         ifelse(ase_data$gene_category == "trans", "blue",
         ifelse(ase_data$gene_category == "trans_compensatory", "green",
         ifelse(ase_data$gene_category == "compensatory", "orange", 
                "gray"))))

plot(ase_data$logFC, -log10(ase_data$FDR),
     col = colors, main = "DE Results Colored by Regulatory Mechanism",
     xlab = "log2FC", ylab = "-log10(FDR)")
```

**Interpretation for viewers**: One glance shows where regulatory mechanisms differ.

---

### Figure 2: Allele Bias vs. logFC (Scatter)

```R
# Plot bias against DE magnitude
plot(ase_data$parent1_ratio, ase_data$logFC,
     col = colors, pch = 19,
     xlab = "Parent1 Allele Ratio",
     ylab = "log2FC",
     main = "Allele Bias vs. Expression Difference")
abline(v=0.5, h=0, lty=2, col="gray")  # Reference lines
```

Shows:
- **Diagonal (lower-left to upper-right)**: Cis regulation (concordant bias and DE)
- **Horizontal spread**: Trans regulation (bias ≈ 0.5, but varying logFC)
- **Diagonal (upper-left to lower-right)**: Trans-compensatory (opposing bias and DE)

---

## Section 4: Troubleshooting Interpretation Issues

### Issue: "Most genes are in `no_divergence`; is analysis broken?"

**This is NORMAL.**

Most genes in any genome show no regulatory divergence. Typical breakdown:
- 60-80% `no_divergence` (background)
- 10-20% `cis`
- 10-15% `trans`
- 5-10% `trans_compensatory` + `compensatory`

If you have >50% DE genes, your groups might be too different (or same tissue not matched well).

---

### Issue: "All genes are trans-regulated; no cis effects detected"

**Possible causes**:
1. **Poor SNP coverage**: If `snp_count` is low across genes, ASE signal is weak
   - Check: Do most genes have `snp_count` = 0?
   - Fix: Need more genetic markers (whole-genome SNP set?)
2. **Parents too similar**: If parents are subspecies, regulatory differences might be small
3. **Hybrid imbalance**: Genome imbalance in hybrid might mask cis signals

**What to do**:
- Check `snp_count` distribution across genes
- Verify parents are genetically distinct (run DE on parent only samples)
- Look at `cis_regulatory_candidates_all.tsv` — are there ANY cis genes?

---

### Issue: "Results don't match my hypothesis about gene X"

**Possible reasons**:
1. **Gene ID mismatch**: Verify ID in output matches your reference
   - Example: Expected `GeneA` but got `GeneA.T1` or different nomenclature
2. **Low SNP count**: `snp_count = 1` or `2` = unreliable
   - ASE signal is noisy with few SNPs
3. **Different tissue**: Cis vs. trans varies by tissue
   - Maybe cis in leaf, trans in root
4. **Confounding factors**: If parents differ in multiple ways (location, growth conditions), effect size might be driven by environment

**What to do**:
- Check the gene's row directly in output file (search by ID)
- Look at `snp_count` specifically for your gene
- Do independent validation (qRT-PCR, allele-specific)

---

## Section 5: Advanced Analysis (Ad/Ed Framework)

If you ran the advanced mode with `--expression-file` and `--metadata-file`, you get the Ad/Ed regulatory framework output.

### Understanding Ad/Ed Estimates

**Ad = Additive effect (cis)**
- Positive `ad`: Cis change favors group1's allele
- Negative `ad`: Cis change favors group2's allele
- Large `|ad|`: Strong cis effect

**Ed = Dominance effect (trans)**
- Positive `ed`: Trans factors favor group1
- Negative `ed`: Trans factors favor group2
- Large `|ed|`: Strong trans effect

### Example Interpretations

```
logFC = 2.0, ad = 1.8, ed = 0.2
→ Cis DOMINATES (90% cis, 10% trans)

logFC = 2.0, ad = 0.2, ed = 1.8
→ Trans DOMINATES (10% cis, 90% trans)

logFC = 2.0, ad = 1.5, ed = 0.5
→ Cis and trans REINFORCE (same direction)

logFC = 2.0, ad = 1.3, ed = -0.7
→ Cis and trans OPPOSE (net effect = sum, but opposite signals)
```

---

## Section 6: Common Biological Interpretations

### Scenario 1: Recent Speciation

**Pattern**:
- High proportion of `trans` genes
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
- Lower proportion of `trans`
- More `trans_compensatory` and `compensatory`

**Interpretation**:
- Older divergence allowed accumulation of cis mutations
- Many trans changes might have been deleterious (lost)
- Complex compensatory relationships evolved

---

### Scenario 3: Tissue-Specific Specialization

**Pattern**:
- Differences between tissues:
  - Leaf tissue: Mostly `cis` genes (photosynthesis genes)
  - Root tissue: Mostly `trans` genes (nutrient uptake)

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

- **Wittkopp et al. (2004)** - Foundational Ad/Ed framework
  - "Evolutionary changes in cis and trans gene regulation"
  - Nature Reviews Genetics

- **Tirosh & Barkai (2008)** - Yeast evolution
  - Shows ~30% cis, ~70% trans in transcriptional divergence
  - Molecular Systems Biology

- **McManus et al. (2010)** - Great Apes
  - Tissue-specific cis vs. trans patterns
  - Nature Reviews Genetics

- **Schaefke et al. (2015)** - Drosophila
  - Evidence for complex cis-trans interactions
  - Genome Biology

---

## See Also

- [Input Format Specification](./input_format_specification.md)
- [Output Files Reference](./output_files_reference.md)
- [ASE + DE Integration Workflow](./ase_de_integration_workflow.md)
- [ASEIntegrate & EdgeRDE Audit Report](../AUDIT_ASEIntegrate_EdgeRDE.md)
