#!/usr/bin/env python3

"""
Main program to run the Salmon pipeline.

This program will take as input the directory with the fastq files that are going to
be mapped to the transcriptome. It will deploy the genome-wide decoys and it going to
produce the table with the quantification of the transcripts derived from Salmon.

It can work using two different approaches:
1) It can quantity the gene expression using a reference genome
2) It can quantify the gene expression using both, a reference genome and an alternative genome, for which case
the table will contain the quantification of the transcripts from the reference genome and the alternative genome.
It could be useful if you care about expression bias in hybrids and other similar cases.

"""

import os
import subprocess
import sys
import logging
from tqdm import tqdm

# Adding src directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

# Add the current directory to the Python path
sys.path.append(os.path.dirname(__file__))

from preprocessing_utilities import *

def main(args):
    # Setting up basic logging configuration
    logging.basicConfig(level=logging.INFO)

    # Define the tqdm format
    tqdm_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]'

    # Parse the arguments
    input_dir = args.input
    transcriptome = args.transcriptome
    reference = args.reference
    alternative = args.alternative
    working_directory = args.working_directory
    output_dir = args.output
    temporal_directory = args.temporal_directory
    reference_name = args.reference_name
    alternative_name = args.alternative_name
    salmon_index_options = args.salmon_index_options.split()
    quant_options = args.quant_options.split()
    threads = args.threads
    chrom_level = args.chrom_level
    memory = args.memory
    clean = args.clean

    # Check if the input directory exists
    if not os.path.isdir(input_dir):
        logging.error('The input directory does not exist')
        sys.exit(1)
    
    # Check if the transcriptome exists and if it is compressed
    if not os.path.isfile(transcriptome):
        logging.error('The transcriptome does not exist')
        sys.exit(1)
    if transcriptome.endswith('.gz'):
        logging.info('Decompressing the transcriptome')
        cmd = f"gunzip {transcriptome}"
        subprocess.run(cmd, shell=True, executable='/bin/bash')

    # Check if the reference genome exists
    if not os.path.isfile(reference):
        logging.error('The reference genome does not exist')
        sys.exit(1)

    # Check if the alternative genome exists
    if alternative is not None and not os.path.isfile(alternative):
        logging.error('The alternative genome does not exist')
        sys.exit(1)
    
    # Check if the working directory exists and make it a full path if '.' is given
    if not os.path.isdir(working_directory):
        logging.error('The working directory does not exist')
        logging.info('Creating the working directory')
        os.mkdir(working_directory)
    if working_directory == '.':
        working_directory = os.getcwd()
        ## DEBUGGING
        print(working_directory)
    
    # Check if the output directory exists
    if not os.path.isdir(output_dir):
        logging.info('Creating the output directory')
        os.mkdir(output_dir)

    # Check if the temporal directory exists
    if not os.path.isdir(temporal_directory):
        logging.info('Creating the temporal directory')
        os.mkdir(temporal_directory)

    # Check if the reference name is given
    if alternative is not None and reference_name is None:
        logging.error('The reference name is not given')
        sys.exit(1)

    # Check if the alternative name is given
    if alternative is not None and alternative_name is None:
        logging.error('The alternative name is not given')
        sys.exit(1)

    # Homogenize the names of the chromosomes and add the specific names
    if alternative is None:
        reference_homogenized = os.path.join(output_dir, os.path.basename(reference).rsplit('.', 1)[0] + '_homogenized.fasta')
        homogenize_headers(reference, reference_homogenized, chrom_level)
        alternative_homogenized = None
    else:
        reference_homogenized, alternative_homogenized = fasta_header_homogenizer(reference, alternative, working_directory, chrom_level)

    # Add the specific names to the chromosomes (Modification in place)
    # This will help to differentiate the chromosomes from the reference and the alternative genomes for further analysis

    if alternative is not None:
        add_name_to_fasta_headers(reference_homogenized, reference_name, remove_bak = True)
        add_name_to_fasta_headers(alternative_homogenized, alternative_name, remove_bak = True)

    # Extract the contig names for the decoys

    decoy_file = os.path.join(temporal_directory, 'decoys.txt')

    get_contig_names_for_decoys(reference_homogenized, alternative_homogenized, decoy_file)

    # Concatenate transcriptome and genomes into a single file

    transcriptome_genome = os.path.join(temporal_directory, 'transcriptome_genome.fasta')

    concatenate_fasta_files(transcriptome, reference_homogenized, alternative_homogenized,
                            transcriptome_genome)
    
    # Create the index for the transcriptome

    run_salmon_index(transcriptome_genome, decoy_file, threads = threads,
                    output_dir = temporal_directory,
                    output_name='salmon_index',
                    salmon_options = salmon_index_options)

    # Run Salmon quantification

    logging.info('Quantifying the files ...')

    files = os.listdir(input_dir)

    # Create the job directory
    job_dir = os.path.join(temporal_directory, 'Jobs')
    if not os.path.isdir(job_dir):
        os.mkdir(job_dir)

    for file in tqdm(files, desc='Quantifying files', bar_format=tqdm_format):
        file_path = os.path.join(input_dir, file)
        
        if file_path.endswith('.gz'):
            # Decompress the file
            cmd_gunzip = f"gunzip {file_path}"
            subprocess.run(cmd_gunzip, shell=True, executable='/bin/bash')
            file_path = file_path.replace('.gz', '')
        
        # Checking that the file includes the fastq extension
        # If not, rename the file on the fly to include the .fastq extension and modify the file_path
        # NOTE: Only fastq files are accepted

        if not file_path.endswith('.fastq'):
            new_file_path = file_path + '.fastq'
            os.rename(file_path, new_file_path)
            file_path = new_file_path

        job = master_script_generator(file_path, working_directory, job_dir, output_dir, threads, memory,
                                    os.path.join(temporal_directory, 'salmon_index'),
                                    quant_options)
        # Run the job
        cmd_launcher = f"sbatch {job}"
        subprocess.run(cmd_launcher, shell=True, executable='/bin/bash')

    # Cleaning all the files and temporal directories
    if clean:
        logging.info('Cleaning the temporal directory')
        cmd_clean = f"rm -r {temporal_directory}"
        subprocess.run(cmd_clean, shell=True, executable='/bin/bash')

    logging.info('Done! Jobs are running in the background')
    logging.info('Check the output directory for the results')
    logging.info('Bye!')

if __name__ == '__main__':
    main()