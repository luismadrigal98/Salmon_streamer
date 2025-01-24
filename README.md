# Salmon Pipeline

## Overview

This project provides a pipeline for RNA-seq data analysis using the Salmon tool. It includes two main components:

1. `salmon_runner.py`: A script to run the Salmon pipeline for quantifying RNA-seq reads.
2. `salmon_output_processor.py`: A script to combine the quantification results from the Salmon pipeline into a single table.

## Requirements

- Python 3.x
- Salmon
- tqdm
- argparse
- subprocess
- logging

## Installation

1. Clone the repository:

    ```sh
    git clone https://github.com/luismadrigal98/Salmon_streamer
    cd Salmon_streamer
    ```

2. Create and activate a conda environment:

    ```sh
    conda create --name salmon_pipeline python=3.8
    conda activate salmon_pipeline
    ```

3. Install the required Python packages:

    ```sh
    pip install tqdm argparse logging
    ```

4. Install Salmon:

    ```sh
    conda install -c bioconda salmon
    ```

## Usage

### Running the Salmon Pipeline

The `salmon_runner.py` script takes several arguments to configure the pipeline. Below is an example of how to run it:

```sh
python [salmon_runner.py](http://_vscodecontentref_/0) -i <input_directory> --transcriptome <transcriptome.fasta> --reference <reference_genome.fasta> --alternative <alternative_genome.fasta> --working_directory <working_directory> --output <output_directory>

Arguments
-i, --input: Directory with the fastq files (required)
--transcriptome, --tf: Transcriptome in fasta format (required)
--reference, -r: Reference genome in fasta format (required)
--alternative, -a: Alternative genome in fasta format (required)
--working_directory, -w: Working directory (required)
--output, -o: Output directory (required)
--temporal_directory, --temp: Temporal directory (optional, default: ./TEMP)
--reference_name, -rn: Name of the reference genome (optional)
--alternative_name, -an: Name of the alternative genome (optional)
--salmon_index_options, -sio: Salmon index options as a single string (optional, default: --keepDuplicates -k 31)
--quant_options, -qo: Salmon quant options as a single string (optional, default: --noLengthCorrection -l U -p 1)
--threads, -t: Number of threads to use (optional, default: 1)
--chrom_level, -c: Whether to remove scaffolds and contigs from the genomes (optional, default: True)
--memory, -m: Memory to use in the cluster for individual quantification jobs (optional, default: 2)
--clean, -cl: Clean the temporal directory after the run (optional, default: False)
```

Combining Quantification Results
The salmon_output_processor.py script combines the quantification results into a single table. Below is an example of how to run it:

```sh
python [salmon_output_processor.py](http://_vscodecontentref_/1) -o <output_directory> --result_name <result_file_name>

Arguments
-o, --output: Directory with the output of the Salmon pipeline (required)
--result_name, --rn: Name of the output file to store the combined results (optional, default: table.txt)
--mode, -m: Mode to run the commands. Options are "cmd" for command line and "python" for Python code (optional, default: cmd)
```

# Authors
Luis Javier Madrigal-Roca, 
Paris Veltsos, 
John K. Kelly

# License
This project is licensed under the MIT License.