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

> **Two workflows are available depending on your data:**
> - **Single-genome workflow** (one reference genome): Start at [Option A](#option-a-single-genome-workflow).
> - **Dual-genome workflow** (reference + alternative genome, for allele-specific expression): Start at [Option B](#option-b-dual-genome-workflow).

---

### Option A: Single-Genome Workflow

Use this when you only have **one reference genome** and want to quantify read counts (e.g., standard RNA-seq differential expression, no allele-specific analysis needed).

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

Use this when you have **two genomes** (reference + alternative) and want allele-specific expression analysis (e.g., hybrid crosses).

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
    --alt-genome-id SWB \
    --ref-genome-id IM767 \
    --liftover-gff transferred_annotation.gff3 \
    --original-ref-gff reference_annotation.gff3 \
    --alt-genome-fasta alternative_genome.fa \
    --ref-genome-fasta reference_genome.fa \
    --output-dir ./transcriptome_output
```

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
    --crosses SWB SF \
    --quant-results-files QUANT_RESULTS_767_vs_SWB QUANT_RESULTS_767_vs_SF \
    --cpm-min 5.0 \
    --output-dir ./processed_output
```

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
    --parent1-samples IM767 \
    --parent2-samples SWB \
    --parent1-name IM767 \
    --parent2-name SWB \
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
    --parent1-samples IM767 \
    --parent2-samples SWB \
    --parent1-name IM767 \
    --parent2-name SWB \
    --output-dir ./DE_results \
    --test-method ttest \
    --fdr-threshold 0.05
```

Key parameters:
- `--expression-file`: Tab-delimited expression matrix (genes as rows, samples as columns)
- `--parent1-samples`: Sample IDs or patterns for parent 1 (space-separated)
- `--parent2-samples`: Sample IDs or patterns for parent 2 (space-separated)
- `--parent1-name`: Name for parent 1 (e.g., IM767)
- `--parent2-name`: Name for parent 2 (e.g., SWB)
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
Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
```

### BibTeX Format
```bibtex
@software{madrigal_roca_2025_salmon_streamer,
  author       = {Madrigal-Roca, Luis Javier and
                  Veltsos, Paris and
                  Kelly, John K.},
  title        = {Salmon Streamer: A Comprehensive Pipeline for RNA-seq to QTL Analysis},
  version      = {1.0.0},
  publisher    = {Zenodo},
  year         = {2025},
  doi          = {10.5281/zenodo.XXXXXXX},
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