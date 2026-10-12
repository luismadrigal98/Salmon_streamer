# Salmon Streamer - Complete RNA-seq to QTL Analysis Pipeline

## Overview

Salmon Streamer is a comprehensive bioinformatics pipeline for RNA-seq data analysis, from initial quantification through genotype/phenotype processing and QTL input preparation. The pipeline integrates multiple tools and analysis steps into a unified command-line interface.

### Key Features

- **Complete RNA-seq Pipeline**: From raw reads to quantified expression
- **Paired-End Read Support**: Full support for both single-end and paired-end sequencing data
- **Allele-Specific Expression**: Handle dual-genome quantification for allele-specific mapping
- **Differential Expression Analysis**: Compare parental line expression with statistical testing
- **Genotype Calling**: Process transHouscript mapping data to call genotypes
- **QTL Preparation**: Generate all inputs needed for QTL analysis
- **Quality Control**: PCA-based sample filtering and outlier detection
- **Modular Design**: Run individual steps or complete workflows

### Available Commands

Salmon Streamer provides the following subcommands:

| Command | Description |
|---------|-------------|
| `AnnTransfer` | Transfer annotations from reference to alternative genome using Liftoff |
| `GenerateTranscriptome` | Generate combined transcriptome from liftover results (dual-genome) |
| `ExtractTranscriptome` | Extract transcriptome from a single reference genome using GFF3 annotations |
| `RunSalmonQuant` | Run Salmon pipeline for RNA-seq quantification (SE/PE) |
| `ProcessSalmonOut` | Process Salmon output into a combined table |
| `TranslateSalmon` | Translate Salmon outputs and organize read counts by allele |
| `CalculateRawReads` | Calculate raw reads per plant from Salmon outputs |
| `CalculateCPM` | Calculate CPM (counts per million) for PCA analysis |
| `ProcessPost` | Run complete post-processing pipeline (combines above three) |
| `PCA_QC` | Perform PCA-based quality control on expression data |
| `ParentalDE` | Differential expression analysis between parental lines |
| **`EdgeRDE`** | **edgeR-based differential expression analysis with PCA and quality control plots** |
| **`EdgeRDEFromNormalized`** | **limma-trend DE on a pre-normalized expression matrix (e.g. after collapsing paralogs in TMM-normalized data)** |
| **`ASEIntegrate`** | **Integrate allele-specific expression with DE results; classify genes as cis/trans-regulated** |
| **`MakeTxp2Gene`** | **Derive a transcript-to-gene map from a GFF3, validated against the transcriptome that was indexed** |
| **`ParalogGroups`** | **Group genes whose reads map equally well to every copy, from Salmon equivalence classes** |
| **`ParalogTreeCheck`** | **Ground-truth those groups against Newick gene trees and calibrate the ambiguity cutoff** |
| **`ParalogMerge`** | **Collapse unresolvable paralogs into single features to make a DE analysis paralog-aware** |
| **`ParalogTximport`** | **Aggregate Salmon quantifications to paralog-aware features with tximport, keeping transcript-length offsets** |
| **`ParalogTracks`** | **Write IGV tracks and loci so paralog groups can be checked against the alignments** |
| `ProcessGenotypes` | Process genotypes from transcript mapping data |
| `MakePhenotypes` | Generate phenotype files from expression data |
| `PrepareQTLInputs` | Prepare inputs for QTL analysis in R/qtl format |
| `RunQTL` | Run QTL analysis with automated job submission |
| `Voom` | Preprocess Salmon output for voom analysis (placeholder) |

### Pipeline Components

1. **Genome Preparation**
   - Annotation transfer between genomes (Liftoff)
   - Combined transcriptome generation (dual-genome)
   - Single-genome transcriptome extraction
   - Support for both single-genome and dual-genome quantification

2. **RNA-seq Quantification**
   - Salmon index building and quantification
   - **Single-end and paired-end read support**
   - Automatic R1/R2 file pairing for PE mode
   - Gzipped file handling
   - Output processing and combination

3. **Post-processing**
   - Allele-specific read count organization
   - CPM calculation for quality control
   - Alternative genome sorting for consecutive alleles

4. **Differential Expression Analysis**
   - Parental line comparison with statistical testing
   - Welch's t-test and Mann-Whitney U test options
   - Benjamini-Hochberg FDR correction
   - Fold change and effect size calculations

5. **Quality Control**
   - PCA-based sample filtering
   - IQR-based outlier detection
   - Read depth filtering

6. **Genotype Processing**
   - Transcript mapping analysis
   - Genotype calling from RNA-seq data
   - Error rate estimation and quality filtering
   - Chromosome-specific processing

7. **Phenotype Generation**
   - Expression data normalization (Box-Cox transformation)
   - Cross-specific phenotype file creation
   - Summary statistics per gene

8. **QTL Input Preparation**
   - R/qtl format file generation
   - Phenotype and genotype integration
   - Genetic map construction

9. **QTL Analysis**
   - Automated SLURM job submission
   - Permutation testing and significance analysis
   - LOD score calculation

## Requirements

### Software Dependencies
- Python 3.8+
- Salmon (for RNA-seq quantification)
- Liftoff (for annotation transfer)
- Minimap2 (for genome alignment)
- R with qtl package (for QTL analysis)

### Python Packages
- pandas
- numpy
- scipy
- tqdm
- argparse
- subprocess
- logging

### Optional for HPC
- SLURM (for cluster job submission)

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/luismadrigal98/Salmon_streamer
    cd Salmon_streamer
    ```

2. Create and activate the conda environment from the provided YAML file:
    ```bash
    conda env create -f conda_env.yml
    conda activate salmon
    ```

    This installs all required dependencies (Python 3.12, Salmon 1.10.3, Liftoff 1.6.2, Minimap2 2.28, pandas, numpy, scipy, tqdm, matplotlib, and more).

    > **Note:** If you prefer to install dependencies manually or need a different Python version, you can instead run:
    > ```bash
    > conda create --name salmon python=3.8
    > conda activate salmon
    > conda install -c bioconda salmon liftoff minimap2
    > pip install pandas numpy scipy tqdm
    > # R and qtl package (if running QTL analysis)
    > conda install -c conda-forge r-base
    > R -e "install.packages('qtl', repos='https://cran.r-project.org')"
    > ```

## Quick Start

> **Two quantification routes are available. Choose by the reference genomes you have:**
> - **Single-genome route** (one reference genome): Salmon quantifies reads against one transcriptome. Start at [Option A](#option-a-single-genome-workflow).
> - **Dual-genome route** (reference + alternative genome): Salmon quantifies reads against a combined transcriptome that holds one copy of each gene per genome. Start at [Option B](#option-b-dual-genome-workflow).
>
> This is a choice of how reads are quantified, not of whether allele-specific analysis is possible. The dual-genome route is not a requirement for it: with a single reference genome, ASE can be obtained from F1 hybrid samples and a VCF covering the parents and the F1s (see the note under [Option B](#option-b-dual-genome-workflow)).

---

### Option A: Single-Genome Workflow

Use this when you have **one reference genome**. Salmon quantifies reads against a single transcriptome and the result is one read-count table (e.g., as input for differential expression with `EdgeRDE`).

#### A1. Extract Transcriptome from Single Genome

Extract transcript sequences directly from your reference genome and its GFF3 annotation:

```bash
python SalmonStreamer.py ExtractTranscriptome \
    --genome-fasta reference_genome.fa \
    --gff-file reference_annotation.gff3 \
    --output-fasta transcriptome.fasta
```

To include gene IDs in FASTA headers (recommended for downstream tools):
```bash
python SalmonStreamer.py ExtractTranscriptome \
    --genome-fasta reference_genome.fa \
    --gff-file reference_annotation.gff3 \
    --output-fasta transcriptome.fasta \
    --include-gene-id
```

#### A2. Run RNA-seq Quantification (Single Genome)

**Single-End reads:**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome transcriptome.fasta \
    --reference reference_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type SE \
    --threads 8
```

**Paired-End reads:**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome transcriptome.fasta \
    --reference reference_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type PE \
    --r1-pattern _R1 \
    --r2-pattern _R2 \
    --threads 8
```

> **Note:** No `--alternative` argument is needed for single-genome quantification.

#### A3. Process Salmon Output

Combine per-sample Salmon outputs into a single read-count table:
```bash
python SalmonStreamer.py ProcessSalmonOut \
    --output ./salmon_output \
    --result_name combined_counts.txt \
    --mode python
```

#### A4. Quality Control with PCA (Optional)

```bash
python SalmonStreamer.py PCA_QC \
    --input-data combined_counts.txt \
    --label-rules-file examples/label_rules.json \
    --output-filtered-data-name filtered_counts.tsv \
    --output-dir ./qc_output \
    --iqr-multiplier 1.5
```

---

### Option B: Dual-Genome Workflow

Use this when you have **two genomes** (reference + alternative) and want allele-specific expression by competitive mapping: Salmon quantifies reads against a combined transcriptome built from both genomes (e.g., hybrid crosses).

> **The two-genome requirement applies to this quantification route only.** It is what Salmon competitive mapping needs, not a requirement of allele-specific analysis in general. With a single reference genome, ASE can be obtained from F1 hybrid samples and a VCF covering the parents and the F1s, and `ASEIntegrate` only joins tables by gene ID. See [Allele-Specific Expression + DE Integration](#allele-specific-expression--de-integration-with-aseintegrate).

### 1. Prepare Genomes and Transcriptome

Transfer annotations from reference to alternative genome:
```bash
python SalmonStreamer.py AnnTransfer \
    --target alternative_genome.fa \
    --reference reference_genome.fa \
    --annotation_gff3 reference_annotation.gff3 \
    --output transferred_annotation.gff3 \
    --intermediate_dir ./liftoff_temp
```

Generate combined transcriptome:
```bash
python SalmonStreamer.py GenerateTranscriptome \
    --alt-genome-id ALT \
    --ref-genome-id REF \
    --liftover-gff transferred_annotation.gff3 \
    --original-ref-gff reference_annotation.gff3 \
    --alt-genome-fasta alternative_genome.fa \
    --ref-genome-fasta reference_genome.fa \
    --output-dir ./transcriptome_output
```

`ALT` and `REF` are placeholders for your own genome identifiers; each sequence in the combined transcriptome is named `<id>_<gene ID>`.

> **Before choosing the identifiers:** the post-processing step ([step 4](#4-complete-post-processing-pipeline), `ProcessPost` / `TranslateSalmon`) is still site-specific. It has the reference label `IM767` hardcoded (`translate_salmon_outputs` in `src/postprocessing_utilities.py`) and does not accept any other reference identifier. If you intend to run that step, the reference identifier has to be `IM767` and each `--crosses` value has to be the `--alt-genome-id` of that cross, with no underscore in it.

### 2. Run RNA-seq Quantification

**Single-End Mode (default):**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome combined_transcriptome.fasta \
    --reference reference_genome.fa \
    --alternative alternative_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type SE \
    --threads 8
```

**Paired-End Mode:**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome combined_transcriptome.fasta \
    --reference reference_genome.fa \
    --alternative alternative_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type PE \
    --r1-pattern _R1 \
    --r2-pattern _R2 \
    --threads 8
```

### 3. Process Salmon Output

```bash
python SalmonStreamer.py ProcessSalmonOut \
    --output ./salmon_output \
    --result_name combined_counts.txt \
    --mode python
```

### 4. Complete Post-processing Pipeline

```bash
python SalmonStreamer.py ProcessPost \
    --crosses ALT1 ALT2 \
    --quant-results-files QUANT_RESULTS_REF_vs_ALT1 QUANT_RESULTS_REF_vs_ALT2 \
    --cpm-min 5.0 \
    --output-dir ./processed_output
```

One quant-results file per cross, in the same order as `--crosses`. This step currently works only when the reference identifier is `IM767`; see the note in [step 1](#1-prepare-genomes-and-transcriptome).

`--genes-file` is optional. If provided, only genes in that file are included in the output. If omitted, all genes are included.

### 5. Quality Control with PCA

```bash
python SalmonStreamer.py PCA_QC \
    --input-data RawSamples_forPCA \
    --label-rules-file examples/label_rules.json \
    --output-filtered-data-name filtered_counts.tsv \
    --output-dir ./qc_output \
    --iqr-multiplier 1.5
```

### 6. Differential Expression Analysis (Parental Lines)

Compare expression between parental lines to identify differentially expressed genes:

```bash
python SalmonStreamer.py ParentalDE \
    --expression-file combined_expression.txt \
    --parent1-samples parent1 \
    --parent2-samples parent2 \
    --parent1-name Parent1 \
    --parent2-name Parent2 \
    --output-dir ./DE_results \
    --test-method ttest \
    --fdr-threshold 0.05
```

Available statistical tests:
- `ttest`: Welch's t-test (default, assumes unequal variances)
- `mannwhitney`: Mann-Whitney U test (non-parametric alternative)

### 7. Process Genotypes

```bash
python SalmonStreamer.py ProcessGenotypes \
    --cross SWB \
    --allele-counts-file allele_counts.combined.updated767.SWB.txt \
    --samples-file SWB.combined.samples.txt \
    --genes-file Genes_to_updated_767_assembly.txt \
    --min-parental-lines 5 \
    --mapping-threshold 0.95 \
    --output-dir ./genotype_output
```

### 8. Generate Phenotype Files

```bash
python SalmonStreamer.py MakePhenotypes \
    --genes-by-cross-file Genes_by_cross.txt \
    --total-reads-files Total_reads_byplant.SF Total_reads_byplant.SWB \
    --readcounts-files Readcounts.updated767.f2_SF.txt Readcounts.updated767.f2_SWB.txt \
    --f2-lists SF.included.f2s.txt SWB.included.f2s.txt \
    --crosses SF SWB \
    --output-dir ./phenotype_files
```

### 9. Prepare QTL Inputs

```bash
python SalmonStreamer.py PrepareQTLInputs \
    --cross SWB \
    --genotype-pp-file SWB.F2_geno_PP.txt \
    --estimates-files estimates.SWB.Chr_01 estimates.SWB.Chr_02 \
    --genes-by-cross-file Genes_by_cross.txt \
    --phenotype-group gene.group1 \
    --output-dir ./qtl_inputs
```

### 10. Run QTL Analysis

```bash
python SalmonStreamer.py RunQTL \
    --phenofile-path phenotype_file.txt \
    --genfile-path genotype_file.txt \
    --outdir-base ./qtl_results \
    --permnum 1000 \
    --max-jobs 100
```

## Detailed Usage

### Individual Pipeline Steps

#### Annotation Transfer
Transfer gene annotations from a reference genome to an alternative genome using Liftoff.

```bash
python SalmonStreamer.py AnnTransfer --help
```

Key parameters:
- `--target`: Alternative genome FASTA
- `--reference`: Reference genome FASTA  
- `--annotation_gff3`: Reference annotation GFF3
- `--output`: Output annotation file

#### Transcriptome Generation (Dual-Genome)
Create a combined transcriptome containing genes from both reference and alternative genomes.

```bash
python SalmonStreamer.py GenerateTranscriptome --help
```

Key parameters:
- `--cov-threshold`: Minimum liftover coverage (default: 0.9)
- `--seqid-threshold`: Minimum sequence identity (default: 0.9)

#### Single-Genome Transcriptome Extraction

Extract transcript sequences directly from a single reference genome and its GFF3 annotation file. Use this when you only have one genome and want standard read-count quantification.

```bash
python SalmonStreamer.py ExtractTranscriptome \
    --genome-fasta reference_genome.fa \
    --gff-file reference_annotation.gff3 \
    --output-fasta transcriptome.fasta
```

Key parameters:
- `--genome-fasta` / `--genome`: Reference genome FASTA file (required)
- `--gff-file` / `--gff`: GFF3 annotation file (required)
- `--output-fasta` / `--output`: Output transcriptome FASTA file (required)
- `--transcript-types`: Comma-separated GFF3 feature types to treat as transcripts (default: `mRNA`)
- `--gene-types`: Comma-separated GFF3 feature types to treat as genes (default: `gene`)
- `--id-prefix`: Prefix to add to transcript IDs in output FASTA headers (optional)
- `--include-gene-id`: Include the parent gene ID in FASTA headers (flag, recommended for downstream tools)
- `--min-length`: Minimum transcript length to include in bp (default: 1)

Example with all optional parameters:
```bash
python SalmonStreamer.py ExtractTranscriptome \
    --genome-fasta reference_genome.fa \
    --gff-file reference_annotation.gff3 \
    --output-fasta transcriptome.fasta \
    --transcript-types mRNA,transcript \
    --include-gene-id \
    --min-length 100
```

#### RNA-seq Quantification
Run Salmon quantification on RNA-seq reads. Works with both single-genome and dual-genome transcriptomes. For single-genome workflows, omit the `--alternative` argument.

**Single-Genome, Single-End reads:**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome transcriptome.fasta \
    --reference reference_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type SE \
    --threads 8
```

**Single-Genome, Paired-End reads:**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome transcriptome.fasta \
    --reference reference_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type PE \
    --r1-pattern _R1 \
    --r2-pattern _R2 \
    --threads 8
```

**Dual-Genome, Single-End Mode:**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome combined_transcriptome.fasta \
    --reference reference_genome.fa \
    --alternative alternative_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type SE \
    --threads 8
```

**Dual-Genome, Paired-End Mode:**
```bash
python SalmonStreamer.py RunSalmonQuant \
    --input ./fastq_files \
    --transcriptome combined_transcriptome.fasta \
    --reference reference_genome.fa \
    --alternative alternative_genome.fa \
    --working_directory ./salmon_work \
    --output ./salmon_output \
    --library-type PE \
    --r1-pattern _R1 \
    --r2-pattern _R2 \
    --threads 8
```

Key parameters:
- `--library-type`: Library type - 'SE' for single-end (default) or 'PE' for paired-end
- `--alternative`: Alternative genome FASTA (omit for single-genome workflows)
- `--r1-pattern`: Pattern to identify R1 files in paired-end mode (default: '_R1')
- `--r2-pattern`: Pattern to identify R2 files in paired-end mode (default: '_R2')
- `--threads`: Number of parallel threads
- `--memory`: Memory per job (for cluster execution)
- `--quant_options`: Salmon quantification options

Supported R1/R2 patterns:
- `_R1` / `_R2` (default): sample_R1.fastq, sample_R2.fastq
- `_1` / `_2`: sample_1.fastq, sample_2.fastq
- `.R1` / `.R2`: sample.R1.fastq, sample.R2.fastq
- Custom patterns can be specified

#### Post-processing Components

**Translate Salmon Outputs**: Organize read counts by allele
```bash
python SalmonStreamer.py TranslateSalmon SWB \
    --quant-results-file QUANT_RESULTS.txt
```

`--genes-file` is optional. If provided, only genes listed in that file are included in the output. If omitted, all genes from the quant results are included.

**Calculate Raw Reads**: Sum total reads per sample
```bash
python SalmonStreamer.py CalculateRawReads \
    --salmon-files output1.txt output2.txt \
    --output-file raw_reads.txt
```

**Calculate CPM**: Normalize to counts per million for PCA
```bash
python SalmonStreamer.py CalculateCPM \
    --raw-reads-file raw_reads.txt \
    --salmon-files output1.txt output2.txt \
    --cpm-min 5.0
```

#### Differential Expression Analysis

Compare expression levels between two parental lines to identify differentially expressed genes:

```bash
python SalmonStreamer.py ParentalDE \
    --expression-file combined_expression.txt \
    --parent1-samples parent1 \
    --parent2-samples parent2 \
    --parent1-name Parent1 \
    --parent2-name Parent2 \
    --output-dir ./DE_results \
    --test-method ttest \
    --fdr-threshold 0.05
```

Key parameters:
- `--expression-file`: Tab-delimited expression matrix (genes as rows, samples as columns)
- `--parent1-samples`: Sample IDs or patterns for parent 1 (space-separated)
- `--parent2-samples`: Sample IDs or patterns for parent 2 (space-separated)
- `--parent1-name`: Name for parent 1 (e.g., Parent1)
- `--parent2-name`: Name for parent 2 (e.g., Parent2)
- `--test-method`: Statistical test to use ('ttest' or 'mannwhitney')
  - `ttest`: Welch's t-test (default, assumes unequal variances)
  - `mannwhitney`: Mann-Whitney U test (non-parametric, for non-normal distributions)
- `--fdr-threshold`: Significance level for FDR correction (default: 0.05)
- `--output-dir`: Output directory for results (default: ./DE_results)

Output includes:
- Gene-level statistics (mean, std, fold change)
- Test statistics (t-statistic/U-statistic, p-value)
- FDR-corrected q-values (Benjamini-Hochberg method)
- Significance flags

#### Differential Expression Analysis with EdgeRDE

**NEW**: Comprehensive edgeR-based differential expression analysis with publication-quality plots and quality control.

EdgeRDE implements a complete quasi-likelihood edgeR pipeline with TMM normalization, adaptive filtering, and all pairwise group contrasts. It produces:
- DE results for every gene across all group comparisons
- Publication-quality plots (volcano, PCA, sample correlation, dispersion)
- Quality control assessments
- Summary statistics

**Basic usage**:
```bash
python SalmonStreamer.py EdgeRDE \
    --expression-file counts.tsv \
    --metadata-file metadata.tsv \
    --output-dir DE_results/ \
    --fdr-threshold 0.05 \
    --logfc-threshold 1.0
```

**Key parameters**:
- `--expression-file`: Tab-separated count matrix (genes × samples, raw counts only)
- `--metadata-file`: Tab-separated metadata with `sample_name` and `group` columns (or `species`+`tissue`)
- `--group-samples`: Alternatively, specify groups inline: `Group1:s1,s2,s3 Group2:s4,s5,s6`
- `--sample-suffix`: Regex pattern to strip from column names before matching metadata (e.g., `"_R1_filtered$"`)
- `--collapse-replicates-by`: Metadata column naming the biological unit each library came from (e.g. `individual`). Libraries sharing a value are summed into one column before the analysis. **Use it whenever several libraries come from one individual** (several tissue samples of one plant, re-sequenced libraries): every column of the matrix is otherwise treated as an independent replicate, which makes the dispersions too small and the p-values too optimistic. The collapsed matrix and a metadata file listing which libraries went into each unit are written to the output directory
- `--fdr-threshold`: FDR cutoff for significance (default: 0.05)
- `--logfc-threshold`: Log2-fold-change threshold for categorization (default: 1.0)

**Example with inline groups**:
```bash
python SalmonStreamer.py EdgeRDE \
    --expression-file counts.tsv \
    --group-samples Parent1:P1_1,P1_2,P1_3 Parent2:P2_1,P2_2,P2_3 F1_hybrid:F1_1,F1_2,F1_3 \
    --output-dir DE_results/
```

**Output files**:
- `{group1}_vs_{group2}_DE_results.tsv`: Complete results (all genes, logFC, FDR, etc.)
- `{group1}_vs_{group2}_significant_genes.tsv`: Filtered to FDR < threshold
- `{group1}_vs_{group2}_volcano.pdf`: Volcano plots per contrast
- `PCA_plot.pdf`: Sample clustering by expression
- `sample_correlation_heatmap.pdf`: Sample similarity heatmap
- Quality control plots: library sizes, dispersion estimates, top DE genes heatmap
- `analysis_summary.txt`: Parameter summary and statistics

**Documentation**: See [ASE + DE Integration Workflow](docs/ase_de_integration_workflow.md) and [Output Files Reference](docs/output_files_reference.md).

#### Differential Expression on a Pre-Normalized Matrix with EdgeRDEFromNormalized

**NEW**: Run DE directly on an already-normalized expression matrix using the **limma-trend** workflow.

Use this when the standard EdgeRDE pipeline cannot be reused because you no longer have raw integer counts. Typical scenarios:
- You exported `TMM_normalized_CPM.tsv` / `TMM_normalized_logCPM.tsv` via EdgeRDE's `--export-normalized-expression` flag and then **collapsed paralogs** (or other grouped features) by summing TMM-normalized CPMs. Note that collapsing paralogs no longer has to happen outside the pipeline — see [Paralog-Aware Differential Expression](#paralog-aware-differential-expression) for the `ParalogGroups → ParalogTreeCheck → ParalogMerge` route, which derives the groups from the data and merges at the count level, before normalization.
- You received a normalized expression matrix from a collaborator and only need the DE step plus the standard plots.

`EdgeRDEFromNormalized` accepts either a log2-CPM matrix (used as-is) or a raw CPM matrix (log-transformed internally with a configurable prior). It fits per-group linear models with `limma::lmFit`, applies all pairwise contrasts with `contrasts.fit`, and shrinks variances with `eBayes(trend = TRUE)` — the workflow the limma user guide recommends for already-normalized expression data. Outputs mirror the EdgeRDE format (same DE table columns, same volcano / PCA / heatmap / summary files), so they slot directly into `ASEIntegrate` downstream.

**Basic usage** (log2-CPM input, scale auto-detected):
```bash
python SalmonStreamer.py EdgeRDEFromNormalized \
    --expression-file TMM_normalized_logCPM.collapsed.tsv \
    --metadata-file metadata.tsv \
    --output-dir DE_results_collapsed/
```

**With raw CPM input** (the script applies `log2(x + prior_count)` internally):
```bash
python SalmonStreamer.py EdgeRDEFromNormalized \
    --expression-file TMM_normalized_CPM.collapsed.tsv \
    --metadata-file metadata.tsv \
    --output-dir DE_results_collapsed/ \
    --input-scale cpm \
    --prior-count 0.25
```

**Key parameters**:
- `--expression-file`: Tab-separated normalized expression matrix (genes × samples, first column = `TranscriptID`)
- `--metadata-file`: Same format as EdgeRDE (`sample_name` + `group`, or `species` + `tissue`)
- `--group-samples`: Inline alternative to `--metadata-file` (e.g. `Group1:s1,s2,s3 Group2:s4,s5,s6`)
- `--input-scale {auto,logcpm,cpm}`: Scale of the input matrix. `auto` flags log-CPM whenever any negative values are present (default)
- `--prior-count`: Pseudo-count added before `log2` when `--input-scale cpm` (default `0.25`, matches edgeR's default)
- `--fdr-threshold`, `--logfc-threshold`, `--sample-suffix`, `--rscript-executable`: same semantics as `EdgeRDE`

**Output files**: identical layout to `EdgeRDE` — `{group1}_vs_{group2}_DE_results.tsv`, `{group1}_vs_{group2}_significant_genes.tsv`, `{group1}_vs_{group2}_volcano.{pdf,png}`, `PCA_plot.{pdf,png}`, `sample_correlation_heatmap.pdf`, `DE_genes_heatmap.pdf`, `analysis_summary.txt`, `session_info.txt`. The DE table columns are renamed from limma's defaults to match the EdgeRDE schema (`adj.P.Val` → `FDR`, `P.Value` → `PValue`, `AveExpr` → `logCPM`) so `ASEIntegrate` works without modification.

**Why limma-trend and not edgeR's QL test?** edgeR's quasi-likelihood pipeline assumes negative-binomial-distributed integer counts and uses per-gene dispersions estimated from the count-level model. Once data have been TMM-normalized to CPM (and especially once they've been summed across paralogs), those distributional assumptions no longer hold. limma-trend instead models log-CPM with a mean-variance trend, which is the standard recommendation for normalized data and is implemented as `lmFit → contrasts.fit → eBayes(trend = TRUE)`.

#### Paralog-Aware Differential Expression

Some genes cannot be told apart by RNAseq reads. When every read that hits one copy hits the others equally well, Salmon's per-gene count is one arbitrary split of a shared pool rather than a measurement, and a DE test on that number is testing the EM's initialisation. These four subcommands find those genes, verify the calls against gene trees, and collapse the unresolvable ones into single features.

The important design point: **ambiguity alone does not make a gene unusable.** The EM estimator is unbiased at every ambiguity level; what degrades is variance, roughly as `1/sqrt(unique reads)`. A gene sharing 90% of its reads but holding thousands of its own is estimated fine, while a gene sharing 40% with almost no unique reads is not. So the merge criterion is *"can this copy be estimated by itself?"* — the conjunction of `--min-ambiguity` and `--min-unique-reads` — not *"is it a paralog?"*.

**Step 0 — get a transcript-to-gene map.** Everything below works at the gene level, so it needs to know which transcripts belong to the same gene. Two routes, depending on whether the transcriptome already exists:

*Building the transcriptome here* — ask for the map at the same time, and it comes out consistent by construction:

```bash
python SalmonStreamer.py ExtractTranscriptome \
    --genome genome.fasta --gff annotation.gff3 \
    --output transcriptome.fasta \
    --txp2gene txp2gene.tsv
```

*Transcriptome already indexed* — derive the map from the GFF3 and check it against what was indexed:

```bash
python SalmonStreamer.py MakeTxp2Gene \
    --gff annotation.gff3 \
    --transcriptome transcriptome.fasta \
    -o txp2gene.tsv
```

Pass `--transcriptome` whenever you can. If the map's transcript names do not match the names Salmon indexed, nothing errors — the unmatched transcripts quietly become their own genes and the paralog groups come out empty or wrong in a way that looks like biology rather than plumbing. `MakeTxp2Gene` reports the overlap up front, fails outright when nothing matches, and takes `--require-complete` to make partial coverage fatal in a scripted pipeline.

Two details worth knowing. `ExtractTranscriptome` keys the map on the **FASTA header it actually writes**, so `--id-prefix` and `--include-gene-id` are accounted for automatically — a map derived separately from the GFF3 would key on the bare transcript ID and match nothing. And the dual-genome `GenerateTranscriptome` path extracts **gene** features directly, so its transcriptome is already gene-level and needs no map at all; omit `--txp2gene` there and each sequence is treated as its own gene, which is correct.

**Step 1 — quantify with equivalence classes.** The group detection reads Salmon's own equivalence-class output, so the quantification must be run with the relevant flags:

```bash
salmon quant -i <index> -l A -1 r1.fq -2 r2.fq -o <outdir> \
    --dumpEq --dumpEqWeights --hardFilter
```

`--hardFilter` is what makes a class label mean "these reads map equally well to exactly this set of transcripts"; without it the labels include sub-optimal mappings and the groups are not interpretable.

**Step 2 — call the groups.**

```bash
python SalmonStreamer.py ParalogGroups \
    --eq-files quant_dirs/*/aux_info/eq_classes.txt.gz \
    --txp2gene txp2gene.tsv \
    --min-samples 3 --min-reads 10 --min-ambiguity 0 \
    -o paralog_groups_all.tsv \
    --gene-table paralog_genes_readsplit.tsv
```

Keep `--min-ambiguity 0` here. Filtering at generation time bakes a cutoff into a file that takes hours to rebuild; the cutoff belongs downstream, where it can be re-derived in seconds.

**Step 3 — calibrate the cutoff against gene trees.** Do not guess this threshold.

```bash
python SalmonStreamer.py ParalogTreeCheck \
    --groups paralog_groups_all.tsv \
    --trees genetrees/*.rtf \
    --tree-report paralog_tree_pairs.tsv \
    --calibration-report paralog_cutoff_calibration.tsv \
    --conversion-report paralog_gene_conversion_candidates.tsv
```

This resolves tree tip labels to gene IDs, classifies every within-genome pair by patristic distance and whether its MRCA is a genome-specific clade, and sweeps the ambiguity cutoff reporting sensitivity against the tree-confirmed duplications. In practice the intuitive-looking `0.9` can retain *none* of them — read the `sensitivity` column rather than picking a round number.

Derive the cutoff from the report rather than reading a value off the grid: the decisive number is the lowest ambiguity among pairs the trees confirm *and* the reads detect, and the sweep includes that exact value as its own row. Re-run this step whenever the trees change, and re-derive the threshold from the new output — a tree edit moves the cutoff, and a threshold left over from the previous round matches neither set of trees.

**Gene-conversion screen.** `--conversion-report` collects pairs where read behaviour and tree position disagree: the reads cannot separate the copies, but the tree places them far apart, or places another taxon's gene closer to one copy than the copies are to each other. Under plain divergent evolution those two axes are redundant, which is why detection normally falls off cleanly with patristic distance; conversion breaks the redundancy by overwriting one copy with the other after the duplication. Tune with `--conversion-min-ambiguity` (default `0.2`) and `--conversion-min-distance` (default `0.05`, below the ~0.08 point where 150 bp reads normally start separating copies).

This is a screen, not a test. A misassembled duplicate, an under-sampled tree, or a badly aligned region produces the same pattern — and so, most commonly, does ordinary rate variation along a gene: a conserved domain keeps an old pair's reads ambiguous with no conversion at all. The screen only sees genes that appear in the supplied trees. Confirm a candidate with sequence-level evidence — a sliding-window identity scan across the pair, or a formal test such as GENECONV — before calling conversion. An empty report is a real result: it means read ambiguity tracks tree divergence with no outliers in the families you supplied.

**Step 4 — build the merged features.**

Read the cutoff out of step 3 rather than typing a literal — that is the whole point of
calibrating. The rule is the weakest pair your trees confirm *and* the reads detect:

```bash
MIN_AMB=$(awk -F'\t' 'NR>1 && $7=="True" && $8=="True" && $9!="" {
              if (!n++ || $9+0 < a) a = $9+0
          } END { if (n) printf "%.3f\n", a; else print "NA" }' paralog_tree_pairs.tsv)
echo "calibrated ambiguity cutoff: $MIN_AMB"
```

```bash
python SalmonStreamer.py ParalogMerge \
    --groups paralog_groups_all.tsv \
    --txp2gene txp2gene.tsv \
    --quant-dirs quant_dirs/*_quant \
    --min-ambiguity "$MIN_AMB" --min-unique-reads 20 \
    -o merged_counts.tsv \
    --out-map paralog_features_map.tsv \
    --out-tx2gene paralog_tx2feature.tsv \
    --out-report paralog_features_report.tsv
```

One thing to know about that number. The tree report's `ambiguity` column is the *less* ambiguous member of each pair, while `ParalogMerge --min-ambiguity` thresholds the *more* ambiguous member of a group (one copy that cannot stand alone is enough to merge). A cutoff read from `ambiguity` therefore keeps every tree-confirmed pair but admits at least as many groups as it was calibrated on. To threshold the same quantity that was calibrated, either derive the cutoff from the report's last column, `ambiguity_most_member`, or add `--require-all-members-ambiguous` to the merge.

`--min-unique-reads` is *not* tree-derived — it comes from the estimability argument above
(precision goes as `1/sqrt(unique reads)`, so 20 is roughly a 22% CV) and stays put when the
trees change.

Merged features are named `PARA_<first gene>_n<k>`. Omit `--quant-dirs`/`--counts` to emit only the feature definition — the map depends on the groups, not on any particular quantification, so one map applies to whichever counts the DE analysis uses.

Two routes downstream:
- **`ParalogTximport` (preferred)** — step 4b below. Aggregates with tximport, which recomputes the average transcript length of each merged feature. Summing counts by hand cannot produce that, and it is what edgeR/DESeq2 use to correct for a feature's effective length differing between samples when its members change in relative abundance (library composition is a separate matter, handled by TMM).
- **`-o merged_counts.tsv`** for a merged matrix directly. `ParalogMerge` verifies read conservation (total in = total out) and stops rather than writing a matrix that fails it. Simpler, and fine when the model does not use length offsets.

**Step 4b — aggregate with tximport (preferred route).**

```bash
python SalmonStreamer.py ParalogTximport --quant-dirs production_quant/*_quant --tx2feature paralog_tx2feature.tsv --output-dir tximport_out/ --strip-suffix '_quant$'
```

Writes `paralog_feature_counts.tsv`, `paralog_feature_abundance_tpm.tsv` and `paralog_feature_length.tsv` into `--output-dir`.

`--counts-from-abundance` decides how the length correction is delivered, and it is the only choice that matters here:

- **`lengthScaledTPM` (default)** folds the correction into the counts, so `paralog_feature_counts.tsv` goes straight into `EdgeRDE` — or any count-based DE tool — with no offset. This is what keeps the route inside the pipeline.
- **`no`** keeps raw counts and leaves the correction to the downstream model, which must then take `paralog_feature_length.tsv` as an offset. More faithful, but only worth it if your consumer supports offsets.

Requires the Bioconductor package `tximport`; install with `R -e 'BiocManager::install("tximport")'`. Sample names default to each `quant.sf`'s containing directory, so `--strip-suffix` (or `--sample-names`) is usually needed to match your metadata.

Groups are *not* a partition — a gene can appear in several called groups — so `ParalogMerge` resolves them to connected components, guaranteeing every gene lands in exactly one feature and nothing is double-counted.

**Step 5 — check the calls by eye (optional but recommended).**

```bash
python SalmonStreamer.py ParalogTracks \
    --groups paralog_groups_all.tsv \
    --gff annotation.gff \
    --min-ambiguity "$MIN_AMB" --min-samples 20 \
    --ambiguity-member most \
    --prefix paralog_groups_calibrated
```

`--ambiguity-member most` applies the cutoff to a group's most ambiguous member, which is what `ParalogMerge` does by default, so the tracks show the groups the merge considered. Without it the cutoff applies to the least ambiguous member (as in `ParalogGroups`) and the tracks show a stricter, smaller set.

Writes `<prefix>.bed` (BED score = ambiguity × 1000, so strong groups shade darker) and `<prefix>_loci.txt` (paste-able IGV coordinates). Load the **plain aligned BAM**, not a secondary-filtered one — the multimapping reads are the evidence, and filtering has removed exactly what you are looking for. Aligners give reads they could not place uniquely MAPQ 0 (STAR uses 0–3, reserving 255 for unique) and IGV draws those pale, so a genuine group reads as a washed-out pile-up spanning every copy.

**Interpreting a merged feature.** A merged feature is not a failed gene. The *sum* is often better measured than a typical single gene because it pools more reads — you lose which copy, not how much. Report the result at the feature level; the individual copy's expression is not separately measurable. Genes left unmerged kept enough unique reads to stand alone and can be read directly, though a gene sharing a substantial read fraction still deserves the caveat in the text.

**Merged features and ASE + DE integration.** `ASEIntegrate` joins DE results to a per-gene ASE table by gene ID. A gene folded into a `PARA_*` feature has no DE row of its own — the row belongs to the feature — and a feature ID matches no gene in the ASE table, so those rows drop out of the ASE + DE join without a warning. This follows from the point above: the copy's own expression is not separately measurable, so there is no per-gene DE value to pair with its ASE. Genes left unmerged keep their own ID and join as usual. To find which feature a gene went into, look it up in `paralog_features_map.tsv` (the `--out-map` table; columns `gene`, `feature`, `n_genes_in_feature`, `merged`).

#### Allele-Specific Expression + DE Integration with ASEIntegrate

**NEW**: Integrate differential expression with allele-specific expression to classify regulatory mechanisms (cis vs. trans).

ASEIntegrate combines DE results with per-gene allele bias to categorize genes into regulatory classes:
- **cis-regulated**: Gene has regulatory mutations (promoter/enhancer divergence)
- **trans-regulated**: Expression differences driven by trans factors (TFs, chromatin modifiers)
- **trans-compensatory**: Cis changes exist but trans factors override them
- **compensatory**: Cis differences masked by trans compensation
- **no divergence**: No regulatory divergence detected

This analysis is essential for understanding the evolutionary basis of expression differences in hybrid crosses.

**`ASEIntegrate` does not compute ASE.** It consumes a finished per-gene ASE table (`--ase-counts-file`) and joins it to the `*_DE_results.tsv` files from EdgeRDE by gene ID. It never reads a genome, so it needs neither a second reference nor the [dual-genome quantification route](#option-b-dual-genome-workflow). What the analysis does need is F1 hybrid samples and SNPs that distinguish the parental alleles, and the allele counting has to be done beforehand with a separate tool.

[ASErtain](https://github.com/luismadrigal98/ASErtain) is the single-reference way to obtain ASE: it works from one reference genome, a VCF covering the parents and the F1s, and the F1 BAM files. Pass its `<prefix>.gene_ase.tsv` straight to `--ase-counts-file`, with `--parent1-label` / `--parent2-label` set to the lineage labels of the ASErtain config. A gene then counts as allelically biased when ASErtain *called* it (a replicate-aware test with an effect floor), in ASErtain's direction. The older summary table (`snp_count`, `bias_ratio`, `predominant_bias`, `<label>_ratio`) is still accepted, but there "biased" means `bias_ratio > 0.5` — a majority vote over SNPs, not a test.

**Basic analysis** (cis/trans categorization):
```bash
python SalmonStreamer.py ASEIntegrate \
    --de-results-dir DE_results/ \
    --ase-counts-file ase_counts.tsv \
    --output-dir ASE_DE_results/ \
    --fdr-threshold 0.05 \
    --parent1-label parent1 \
    --parent2-label parent2
```

**Advanced analysis** (adds the cis/trans classification that uses the F1's own expression, plus the Ad/Ed deviations):
```bash
python SalmonStreamer.py ASEIntegrate \
    --de-results-dir DE_results/ \
    --ase-counts-file ase_counts.tsv \
    --output-dir ASE_DE_results_advanced/ \
    --fdr-threshold 0.05 \
    --parent1-label BellPepper \
    --parent2-label ChiliPepper \
    --expression-file counts.tsv \
    --metadata-file metadata.tsv \
    --hybrid-group F1_hybrid \
    --parent1-group BellPepper \
    --parent2-group ChiliPepper
```

**Key parameters**:
- `--de-results-dir`: Directory containing `*_DE_results.tsv` files from EdgeRDE
- `--ase-counts-file`: Tab-separated file with per-gene allele counts (see documentation)
- `--output-dir`: Output directory
- `--parent1-label`: Label for parent 1 allele (must match ASE file column name)
- `--parent2-label`: Label for parent 2 allele (must match ASE file column name)
- `--fdr-threshold`: FDR cutoff for DE significance
- `--parent1-group`, `--parent2-group`: The two parental group names as they appear in the DE file names. They say which comparison is the parental one and which way round it is; give them whenever the group names do not contain the parent labels
- `--expression-file`, `--metadata-file`, `--hybrid-group`: Add these (with the two parent groups) for the advanced analysis

**Output files**:
- `{comparison}_with_ASE.tsv`: Merged DE + ASE; the category is in the `ase_de_category` column
- `{comparison}_ASE_category_summary.tsv`: One row per category (`ase_de_category`, `n_genes`, `avg_logFC_parent1_vs_parent2`, `avg_parent1_ratio`)
- `{comparison}_ASE_vs_DE.pdf` / `.png`: Scatter plot (parent 1 / parent 2 logFC vs. parent 1 allele ratio)
- `cis_regulatory_candidates_all.tsv`: All genes labelled `cis`
- Advanced mode also produces:
  - `ase_vs_de_ad_ed_data.tsv`: Per-gene `cis_log2`, `parental_log2`, `trans_log2`, the `Regulation` class, and the `Ad` / `Ed` deviations
  - `regulatory_divergence_regressions.pdf` / `.png`: `Ad` vs. `Ed` for parentally DE genes, coloured by class

Only the parent 1 versus parent 2 comparison gets these files. `{comparison}_with_ASE.tsv` carries `logFC_parent1_vs_parent2` (the fold change re-oriented so that positive means higher in parent 1) next to the original `logFC`, plus `cis_log2` (log2 F1 allelic ratio) and `trans_log2` (parental log-ratio minus `cis_log2`).

The advanced `Regulation` class is the standard decomposition — `Cis-only`, `Trans-only`, `Cis + Trans (Reinforcing)`, `Cis x Trans (Opposing)`, `Compensatory`, `Conserved`, `Ambiguous` — with descriptive thresholds (allele fraction outside 0.5 ± 0.1 for cis, `|trans_log2| ≥ 1` for trans) and no test of the trans component. `Ad` and `Ed` are reported alongside: `Ed` is where the F1's total expression sits between the parents, a measure of dominance rather than of trans, and it is left empty where the parents are not differentially expressed. Definitions are in the [workflow guide](docs/ase_de_integration_workflow.md#advanced-analysis-aded-classification).

**Gene category interpretation** (values of `ase_de_category`):
| Category | DE? | Allele Bias? | Meaning |
|----------|-----|--------------|---------|
| cis | Yes | Yes (same direction) | Regulatory mutations in this gene |
| trans_only | Yes | No | Trans factors (e.g., TF) regulate expression |
| trans_compensatory | Yes | Yes (opposite) | Trans factors override cis effects |
| compensatory | No | Yes | Cis differences masked by trans |
| no_divergence | No | No | No regulatory divergence |

**Which comparison is classified.** Cis/trans is defined for parent 1 versus parent 2 only, and the rule needs a positive fold change to mean "higher in parent 1". EdgeRDE writes every pairwise comparison, names each `A_vs_B` with the groups in alphabetical order, and reports `logFC` as `A` minus `B`, so neither the right file nor its sign can be assumed. `ASEIntegrate` therefore finds the parental comparison, flips the sign when parent 2's group is named first, and skips every other file (it says so in the log). With `--parent1-group` / `--parent2-group` the match is exact; without them it is inferred from the parent labels occurring in the group names, which fails safely — nothing is classified and the run stops — when they do not.

**Example workflow** (illustrative names, Pepper hybrid cross):
```bash
# Step 1: Run EdgeRDE on expression counts
python SalmonStreamer.py EdgeRDE \
    --expression-file pepper_counts.tsv \
    --metadata-file pepper_metadata.tsv \
    --output-dir pepper_DE/

# Step 2: Run ASEIntegrate with ASE data
python SalmonStreamer.py ASEIntegrate \
    --de-results-dir pepper_DE/ \
    --ase-counts-file pepper_ase.tsv \
    --output-dir pepper_ASE_DE/ \
    --parent1-label BellPepper \
    --parent2-label ChiliPepper \
    --expression-file pepper_counts.tsv \
    --metadata-file pepper_metadata.tsv \
    --hybrid-group F1_hybrid \
    --parent1-group BellPepper \
    --parent2-group ChiliPepper
```

**Documentation & Examples**:
- [Complete workflow guide](docs/ase_de_integration_workflow.md)
- [Input format specification](docs/input_format_specification.md)
- [Output file reference](docs/output_files_reference.md)
- [Interpretation guidelines](docs/interpretation_guidelines.md)
- [EdgeRDE example](examples/EdgeRDE_example_arabidopsis.md)
- [ASE+DE example](examples/ASE_DE_example_pepper.md)

#### Genotype Processing Pipeline

The genotype processing pipeline performs comprehensive analysis of transcript mapping data:

1. **Statistical Analysis** (P1): Calculate allele-specific mapping statistics
2. **Gene Filtering** (P2): Filter genes based on mapping thresholds
3. **Genotype Calling** (P3): Call genotypes for F2 individuals
4. **Quality Control** (P4): Identify and remove problematic markers/samples  
5. **Error Estimation** (P5): Estimate genotyping error rates per chromosome
6. **Parameter Optimization** (P6): Optimize recombination and error parameters

Key parameters:
- `--min-parental-lines`: Minimum parental lines for a gene to be included
- `--mapping-threshold`: Threshold for allele-specific mapping accuracy
- `--min-reads-per-plant`: Minimum total reads per individual
- `--homozygous-threshold`: Threshold for homozygous genotype calls
- `--het-maf`: Minor allele frequency threshold for heterozygous calls

#### Phenotype Generation

Generate normalized phenotype files for QTL analysis:

- Applies Box-Cox transformation for normalization
- Creates cross-specific phenotype files
- Calculates summary statistics per gene
- Filters samples based on read depth

#### QTL Input Preparation

Prepare files in R/qtl format:

- Converts genotype data to R/qtl format
- Creates phenotype matrices for each gene group
- Generates genetic maps with marker positions
- Outputs ready-to-use R/qtl input files

### Advanced Configuration

#### SLURM Integration

For HPC environments, the pipeline supports SLURM job submission:

```bash
python SalmonStreamer.py RunQTL \
    --phenofile-path phenotypes.txt \
    --genfile-path genotypes.txt \
    --outdir-base ./results \
    --partition your_partition \
    --cpus-per-task 4 \
    --mem-per-cpu 8G \
    --max-concurrent-jobs 1000
```

`RunSalmonQuant` submits one SLURM job per sample. If your cluster caps how many jobs one user can have queued (sbatch fails with `QOSMaxSubmitJobPerUserLimit` or a similar `MaxSubmitJob` error), the pipeline waits and retries instead of dropping the sample. It retries every `--submit_retry_interval` seconds (default: 60) until your earlier jobs finish or `--submit_max_wait` seconds of total waiting have passed (default: 3600). Samples still unsubmitted at that point are listed in `<temporal_directory>/Jobs/unsubmitted_jobs.txt`, the run exits with an error, and `--clean` is skipped so the job scripts are kept. Submit them later with:

```bash
while read -r job; do sbatch "$job"; done < TEMP/Jobs/unsubmitted_jobs.txt
```

If you launch `RunSalmonQuant` itself as a SLURM job, give it enough `--time` to cover `--submit_max_wait` on top of building the index.

#### Quality Control Parameters

PCA-based quality control with customizable parameters:

```bash
python SalmonStreamer.py PCA_QC \
    --input-data expression_data.txt \
    --label-rules-file label_config.json \
    --iqr-multiplier 2.0 \
    --pc-to-retain 10 \
    --genes-as-rows
```

### Input File Formats

#### Label Rules JSON
For PCA quality control, create a JSON file defining sample groupings:

```json
{
    "source_patterns": {
        "IM767": ["767"],
        "Alternative": ["ALT", "SWB", "SF"]
    },
    "group_patterns": {
        "Parent": ["parent", "P"],
        "F1": ["f1", "F1"],
        "F2": ["f2", "F2"]
    }
}
```

#### Genes Mapping File
Tab-delimited file mapping genes between assemblies:
```
Chrom	stpos	endpos	old_name	new_name	62	155	444	502	541	664	909	1034	1192	scored_pops	new_chrom	new_start	new_end
Chr_01	13982	16715	MiIM7v11000002m.g	MgIM767.01G000100	yes	yes	no	no	yes	no	yes	yes	no	5	Chr_01	13982	16715
```

#### Sample Information File
Tab-delimited file defining sample types:
```
sample_id	type
J-P-001	SF
s4_62-107	f2
767-032	IM767
```

### Output Files

#### Genotype Processing
- `{cross}.marker.genes.txt`: List of high-quality marker genes
- `{cross}.F2_geno_PP.txt`: Genotype posterior probabilities
- `PC1.{cross}/CleanedCalls.stringent.{cross}.{chrom}`: Cleaned genotype calls
- `estimates.{cross}.{chrom}`: Error and recombination rate estimates

#### Phenotype Generation  
- `pfiles/f2_p_{gene_id}`: Individual gene phenotype files
- `f2summaries_by_gene.txt`: Summary statistics per gene

#### QTL Inputs
- `{cross}.rQTL.genotype.txt`: R/qtl format genotype file
- `{phenotype_group}_{cross}.txt`: Cross-specific phenotype matrices

## Troubleshooting

### Common Issues

1. **Memory Errors**: Increase `--memory` parameter for large datasets
2. **Missing Dependencies**: Ensure all required software is installed and in PATH
3. **File Format Issues**: Check that input files match expected formats
4. **Permission Errors**: Ensure write permissions in output directories

### Getting Help

For detailed help on any subcommand:
```bash
python SalmonStreamer.py <subcommand> --help
```

View all available subcommands:
```bash
python SalmonStreamer.py --help
```

### Performance Tips

1. **Use Multiple Threads**: Set `--threads` to number of available CPU cores
2. **Optimize Memory**: Adjust `--memory` based on dataset size and available RAM
3. **Batch Processing**: For many samples, consider splitting into smaller batches
4. **Storage**: Use fast storage (SSD) for intermediate files when possible

## Citation

If you use Salmon Streamer in your research, please cite:

### Software Citation
```
Madrigal-Roca, L.J., Veltsos, P., & Kelly, J.K. (2025). 
Salmon Streamer: A Comprehensive Pipeline for RNA-seq to QTL Analysis (Version 1.0.0) [Computer software]. 
Zenodo. 
```https://doi.org/10.5281/zenodo.20381998

### BibTeX Format
```bibtex
@software{madrigal_roca_2025_salmon_streamer,
  author       = {Madrigal-Roca, Luis Javier and
                  Veltsos, Paris and
                  Kelly, John K.},
  title        = {Salmon Streamer: A Comprehensive Pipeline for RNA-seq to QTL Analysis},
  version      = {1.2.0},
  publisher    = {Zenodo},
  year         = {2025},
  doi          = {10.5281/zenodo.20381998},
  url          = {https://github.com/luismadrigal98/Salmon_streamer}
}
```

### Author Information and ORCIDs
- **Luis Javier Madrigal-Roca** - University of Kansas  
  ORCID: [0000-0002-5485-6395](https://orcid.org/0000-0002-5485-6395)
  
- **Paris Veltsos** - Vrije Universiteit Brussel  
  ORCID: [0000-0002-8872-6281](https://orcid.org/0000-0002-8872-6281)
  
- **John K. Kelly** - University of Kansas  
  ORCID: [0000-0001-9480-1252](https://orcid.org/0000-0001-9480-1252)

### Related Publications
If you use specific components of the pipeline, please also consider citing the underlying tools:

- **Salmon**: Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature methods, 14(4), 417-419. DOI: [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)

- **Liftoff**: Shumate, A., & Salzberg, S. L. (2021). Liftoff: accurate mapping of gene annotations. Bioinformatics, 37(12), 1639-1643. DOI: [10.1093/bioinformatics/btaa1016](https://doi.org/10.1093/bioinformatics/btaa1016)

- **R/qtl**: Broman, K. W., Wu, H., Sen, Ś., & Churchill, G. A. (2003). R/qtl: QTL mapping in experimental crosses. Bioinformatics, 19(7), 889-890. DOI: [10.1093/bioinformatics/btg112](https://doi.org/10.1093/bioinformatics/btg112)

## License

This project is licensed under the terms specified in the LICENSE file.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Recent Updates

### Version 1.2.0 (May 2026)

**New Features:**
- ✨ **EdgeRDE module**: Full edgeR quasi-likelihood differential expression pipeline with TMM normalization, adaptive `filterByExpr` filtering, all pairwise group contrasts, PCA, sample-correlation heatmap, volcano plots, and a summary report
- 🔧 **Audit fixes (unreleased)**:
  - `ASEIntegrate` classifies the parent 1 versus parent 2 comparison only and orients its fold change (`logFC_parent1_vs_parent2`); before, every pairwise DE file was classified and `cis` / `trans_compensatory` came out exchanged whenever parent 2's group sorted first. `--parent1-group` / `--parent2-group` can be given without the advanced inputs to name the comparison
  - `ASEIntegrate` reads an ASErtain `gene_ase.tsv` directly (bias = ASErtain's call). The legacy summary table is still accepted, with a note that `bias_ratio > 0.5` is not a test
  - The advanced classification now follows the standard cis/trans decomposition (`cis_log2`, `trans_log2`; adds `Cis-only`, `Compensatory`, `Conserved`). `Ed` is reported as what it is, a dominance measure, and only for parentally DE genes. `edgeR` is needed for this block only
  - `EdgeRDE --collapse-replicates-by` sums libraries of one biological unit before the analysis
  - `ParalogTreeCheck` reports `ambiguity_most_member`; `ParalogTracks --ambiguity-member most` draws the groups `ParalogMerge` thresholds
- ✨ **EdgeRDEFromNormalized module**: limma-trend DE on a pre-normalized expression matrix. Designed for the "export TMM → collapse paralogs → re-analyze" workflow, but works on any matrix on a sample-comparable log-CPM or CPM scale. Outputs mirror EdgeRDE so `ASEIntegrate` plugs in unchanged.
- ✨ **ASEIntegrate module**: Combine EdgeRDE results with per-gene allele-specific expression to classify regulatory divergence (cis, trans, cis+trans, compensatory) and optionally estimate Ad/Ed effect sizes
- ✨ **TMM-normalized expression export**: New `--export-normalized-expression` flag on `EdgeRDE` writes the TMM-normalized CPM matrix, log2-CPM matrix, and per-sample normalization factors for manual inspection. Format selectable via `--normalized-expression-format {tsv,csv}` (default `tsv`)

**Improvements:**
- Expression file validation now accepts both integer and float counts
- Expanded examples directory with EdgeRDE (Arabidopsis) and ASE + DE (pepper) workflows
- Documentation: added ASE + DE integration workflow, input format spec, output files reference, and interpretation guidelines

### Version 1.1.0 (October 2025)

**New Features:**
- ✨ **Paired-End Read Support**: Full support for paired-end RNA-seq data
  - Automatic R1/R2 file pairing with customizable patterns
  - Support for gzipped paired-end files
  - Backward compatible with existing single-end workflows
  
- ✨ **Single-Genome Transcriptome Extraction**: New `ExtractTranscriptome` module
  - Extract transcriptomes directly from a single reference genome and GFF3 annotation
  - Enables standard RNA-seq read-count quantification without a second genome
  - Handles strand-specific extraction with reverse complement for minus-strand genes
  - Supports flexible feature type selection (mRNA, transcript, etc.)

- ✨ **Differential Expression Analysis**: New `ParentalDE` module for parental line comparisons
  - Welch's t-test for parametric comparisons
  - Mann-Whitney U test for non-parametric comparisons
  - Benjamini-Hochberg FDR correction
  - Comprehensive statistical output

**Bug Fixes:**
- Fixed IndexError in `translate_salmon_outputs` function with improved error handling
- Fixed malformed output in cmd mode for `combine_results` (duplicate headers issue)
- Improved file format validation in post-processing utilities

**Improvements:**
- Added alternative genome sorting for consecutive allele grouping in output tables
- Enhanced error messages and logging throughout the pipeline
- Improved file processing with better permission checking
- Refactored code for better maintainability

## Changelog

### v1.2.0 - May 2026
- Added `EdgeRDE` subcommand: edgeR quasi-likelihood DE pipeline with TMM normalization, adaptive filtering, all pairwise contrasts, PCA, sample-correlation heatmap, volcano plots, dispersion plots, and a summary report.
- Added `EdgeRDEFromNormalized` subcommand: limma-trend DE on a pre-normalized expression matrix (e.g. paralog-collapsed TMM CPM/log-CPM). Auto-detects log-CPM vs CPM input, produces the same DE table schema and plot set as `EdgeRDE`, and is compatible with `ASEIntegrate` downstream.
- Added `ASEIntegrate` subcommand: merges EdgeRDE results with per-gene ASE counts, classifies regulatory divergence, and optionally estimates Ad/Ed effect sizes.
- Added `--export-normalized-expression` and `--normalized-expression-format` to `EdgeRDE`, producing `TMM_normalized_CPM`, `TMM_normalized_logCPM`, and `TMM_norm_factors` files for manual inspection.
- Expression file validation now accepts integer or float counts.
- Expanded `docs/` and `examples/` with EdgeRDE / ASE + DE integration material.

### v1.1.1 - March 2026
- **Fixed** `--chrom_level` / `-c` flag now defaults to off (keep all sequences). Previously defaulted to `True`, which silently discarded all scaffold/contig-named sequences and caused an empty decoy file error for assemblies using scaffold-based naming.
- **Fixed** paired-end mode: `pair_fastq_files` was called with a directory path instead of a file list, causing no pairs to be found.
- **Fixed** paired-end mode: corrected iteration over `pair_fastq_files` return value (list of tuples, not a dict).
- Added validation after genome homogenization: clear error message if the homogenized file is empty.

### v1.1.0 - October 2025
- Added single-genome transcriptome extraction via `ExtractTranscriptome` module
- Added paired-end read support with automatic R1/R2 pairing
- Implemented ParentalDE module for differential expression analysis
- Fixed post-processing bugs in translation and combination functions
- Added sorting functionality for alternative genome outputs
- Enhanced error handling and validation

### v1.0.0 - Initial Release
- Complete RNA-seq to QTL pipeline
- Single-end read support
- Allele-specific expression quantification
- Genotype calling and processing
- QTL input preparation

## Legacy Usage Information

For reference, the original individual scripts can still be run separately:

### Running the Salmon Pipeline

The `salmon_runner.py` script takes several arguments to configure the pipeline. Below is an example of how to run it:

```sh
python salmon_runner.py -i <input_directory> --transcriptome <transcriptome.fasta> --reference <reference_genome.fasta> --alternative <alternative_genome.fasta> --working_directory <working_directory> --output <output_directory>
```

#### Arguments
- `-i, --input`: Directory with the fastq files (required)
- `--transcriptome, --tf`: Transcriptome in fasta format (required)
- `--reference, -r`: Reference genome in fasta format (required)
- `--alternative, -a`: Alternative genome in fasta format (required)
- `--working_directory, -w`: Working directory (required)
- `--output, -o`: Output directory (required)
- `--temporal_directory, --temp`: Temporal directory (optional, default: ./TEMP)
- `--reference_name, -rn`: Name of the reference genome (optional)
- `--alternative_name, -an`: Name of the alternative genome (optional)
- `--salmon_index_options, -sio`: Salmon index options as a single string (optional, default: --keepDuplicates -k 31)
- `--quant_options, -qo`: Salmon quant options as a single string (optional, default: --noLengthCorrection -l U -p 1)
- `--threads, -t`: Number of threads to use (optional, default: 1)
- `--chrom_level, -c`: If set, remove scaffolds and contigs from the genomes (optional flag, default: off — all sequences are kept)
- `--memory, -m`: Memory to use in the cluster for individual quantification jobs (optional, default: 2)
- `--clean, -cl`: Clean the temporal directory after the run (optional, default: False)

### Combining Quantification Results

The `salmon_output_processor.py` script combines the quantification results into a single table. Below is an example of how to run it:

```sh
python salmon_output_processor.py -o <output_directory> --result_name <result_file_name>
```

#### Arguments
- `-o, --output`: Directory with the output of the Salmon pipeline (required)
- `--result_name, --rn`: Name of the output file to store the combined results (optional, default: table.txt)
- `--mode, -m`: Mode to run the commands. Options are "cmd" for command line and "python" for Python code (optional, default: cmd)