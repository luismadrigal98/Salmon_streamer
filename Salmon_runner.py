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
import argparse
import subprocess
import sys
import logging
from tqdm import tqdm

# Adding src directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

# Add the current directory to the Python path
sys.path.append(os.path.dirname(__file__))

from preprocessing_utilities import *

def main():

    parser = argparse.ArgumentParser(description='Salmon pipeline. This pipeline allows the quantification of RNA-seq reads using a reference genome, the raw reads in fastq format, and the fast selective alignment algorithm (SA) implemented through salmon software')

    global_group = parser.add_argument_group(
        'Global arguments', 
        description = 'Directory and path settings'
    )

    global_group.add_argument(
        '-i', '--input', 
        required=True, 
        help='Directory with the fastq files'
    )

    global_group.add_argument(
        '--transcriptome', '--tf',
        required=True,
        help='Transcriptome in fasta format'
    )

    global_group.add_argument(
        '--reference', '-r',
        required=True,
        help='Reference genome in fasta format'
    )

    global_group.add_argument(
        '--alternative', '-a',
        required=False,
        default=None,
        help='Alternative genome in fasta format'
    )

    global_group.add_argument(
        '--working_directory', '-w',
        required=True,
        help='Working directory',
        default='.'
    )

    global_group.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory'
    )

    global_group.add_argument(
        '--temporal_directory', '--temp',
        required=False,
        help='Temporal directory',
        default='./TEMP'
    )

    naming_group = parser.add_argument_group(
        'Naming arguments',
        description='Naming settings the chromosomes'
    )

    naming_group.add_argument(
        '--reference_name', '-rn',
        required=False,
        help='Name of the reference genome'
    )

    naming_group.add_argument(
        '--alternative_name', '-an',
        required=False,
        help='Name of the alternative genome. It is advisable to use a different name from the reference genome when both provided'
    )

    salmon_index_group = parser.add_argument_group(
        'Salmon index arguments',
        description='Salmon index settings'
    )

    salmon_index_group.add_argument(
        '--salmon_index_options', '-sio',
        required=False,
        help='Salmon index options as a single string',
        default='--keepDuplicates -k 31'
    )

    salmon_quant_group = parser.add_argument_group(
        'Salmon quant arguments',
        description='Salmon quant settings'
    )

    salmon_quant_group.add_argument(
        '--quant_options', '-qo',
        required=False,
        help='Salmon quant options as a single string',
        default='--noLengthCorrection -l U -p 1'
    )

    miscellanous_group = parser.add_argument_group(
        'Miscellanous arguments',
        description='Miscellanous settings'
    )

    miscellanous_group.add_argument(
        '--threads', '-t',
        type=int,
        required=False,
        help='Number of threads to use',
        default=1
    )

    miscellanous_group.add_argument(
        '--chrom_level', '-c',
        required=False,
        help='Whether you want to remove scaffolds and contigs from the genomes',
        default=True
    )

    miscellanous_group.add_argument(
        '--memory', '-m',
        type=int,
        required=False,
        help='Memory to use in the computing cluster for individual quantification jobs',
        default=2
    )

    miscellanous_group.add_argument(
        '--clean', '-cl',
        action='store_true',
        help='Clean the temporal directory after the run',
        default=False
    )

    # Setting up basic logging configuration
    logging.basicConfig(level=logging.INFO)

    # Define the tqdm format
    tqdm_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]'

    # Parse the arguments
    args = parser.parse_args()

    # Retrieve the arguments from the parser
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
        reference_homogenized = os.path.join(output_dir, os.path.basename(ref_genome).rsplit('.', 1)[0] + '_homogenized.fasta')
        reference_homogenized = homogenize_headers(reference, reference_homogenized, chrom_level)
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