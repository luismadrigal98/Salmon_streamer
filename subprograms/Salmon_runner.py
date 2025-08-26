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
    email = args.email
    partition = args.partition
    conda_env = args.conda_env
    time_limit = args.time_limit
    module_load_cmd = args.module_load_cmd

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

    # Verify the index was created successfully and wait for file system sync
    import time
    salmon_index_path = os.path.join(temporal_directory, 'salmon_index')
    version_info_path = os.path.join(salmon_index_path, 'versionInfo.json')
    
    # Wait up to 30 seconds for the index to be fully written
    max_wait = 30
    wait_time = 0
    while not os.path.exists(version_info_path) and wait_time < max_wait:
        logging.info(f"Waiting for salmon index to be fully created... ({wait_time}s)")
        time.sleep(2)
        wait_time += 2
    
    if not os.path.exists(version_info_path):
        logging.error(f"Salmon index creation failed: {version_info_path} does not exist after {max_wait} seconds")
        sys.exit(1)
    
    logging.info(f"Salmon index verified successfully at {salmon_index_path}")

    # Run Salmon quantification

    logging.info('Quantifying the files ...')

    files = os.listdir(input_dir)

    # Create the job directory
    job_dir = os.path.join(temporal_directory, 'Jobs')
    if not os.path.isdir(job_dir):
        os.mkdir(job_dir)

    for file in tqdm(files, desc='Quantifying files', bar_format=tqdm_format):
        original_file_path = os.path.join(input_dir, file) # Store original path for logging
        file_path = original_file_path # Start with the original path
        decompressed_successfully = True # Flag to track decompression status
        output_dir_path = os.path.dirname(file_path) # Directory where output file will be written

        if file_path.endswith('.gz'):
            decompressed_path = file_path.replace('.gz', '')

            # --- Permission Check ---
            can_read_gz = os.access(file_path, os.R_OK)
            can_write_here = os.access(output_dir_path, os.W_OK)

            if not can_read_gz:
                logging.error(f"Permission denied: Cannot read input file {file_path}. Skipping.")
                decompressed_successfully = False
                continue
            if not can_write_here:
                logging.error(f"Permission denied: Cannot write to directory {output_dir_path} for decompression. Skipping {file_path}.")
                decompressed_successfully = False
                continue
            # --- End Permission Check ---

            # Use gunzip -k to keep the original file, which is safer
            # Use Python 3.6 compatible arguments for subprocess.run
            cmd_gunzip = f"gunzip -kf {file_path}" # -k: keep, -f: force overwrite (if .fq exists)
            logging.info(f"Attempting to decompress: {file_path}")
            try:
                result = subprocess.run(
                    cmd_gunzip,
                    shell=True,
                    executable='/bin/bash',
                    stdout=subprocess.PIPE, # Python 3.6 equivalent of capture_output=True
                    stderr=subprocess.PIPE, # Python 3.6 equivalent of capture_output=True
                    universal_newlines=True # Python 3.6 equivalent of text=True
                )

                if result.returncode == 0:
                    # Check if the decompressed file actually exists now
                    if os.path.exists(decompressed_path):
                        logging.info(f"Successfully decompressed {file_path} to {decompressed_path}")
                        file_path = decompressed_path # Update file_path ONLY on success
                    else:
                        # This might happen if gunzip -f overwrites but then fails, or other edge cases
                        logging.error(f"gunzip command succeeded but output file {decompressed_path} not found. Skipping {original_file_path}.")
                        decompressed_successfully = False
                        continue
                else:
                    logging.error(f"Failed to decompress {file_path}. Return code: {result.returncode}")
                    if result.stderr:
                        logging.error(f"Stderr: {result.stderr.strip()}")
                    if result.stdout: # Log stdout too, might contain info
                        logging.error(f"Stdout: {result.stdout.strip()}")
                    logging.warning(f"Skipping file {original_file_path} due to decompression error.")
                    decompressed_successfully = False
                    continue # Skip to the next file in the loop
            except Exception as e:
                 logging.error(f"An unexpected error occurred during gunzip subprocess for {file_path}: {e}")
                 decompressed_successfully = False
                 continue


        # Proceed only if decompression was successful (or not needed)
        if decompressed_successfully:
            # Checking that the file includes the fastq extension
            # If not, rename the file on the fly to include the .fastq extension and modify the file_path

            # There could be instances where the file extension is fq. Given that upstream fastq is expected,
            # we will rename the file to fastq

            if file_path.endswith('.fq'):
                # Check if the target file exists before attempting rename
                if not os.path.exists(file_path):
                    logging.error(f"File not found before renaming: {file_path}. This might be due to a prior error. Skipping.")
                    continue
                new_file_path = file_path.replace('.fq', '.fastq')
                try:
                    # Check write permissions in the directory before renaming
                    if not os.access(os.path.dirname(file_path), os.W_OK):
                         logging.error(f"Permission denied: Cannot rename file in directory {os.path.dirname(file_path)}. Skipping {original_file_path}.")
                         continue
                    logging.info(f"Renaming {file_path} to {new_file_path}")
                    os.rename(file_path, new_file_path)
                    file_path = new_file_path # Update file_path after successful rename
                except OSError as e:
                    logging.error(f"Failed to rename {file_path} to {new_file_path}: {e}")
                    logging.warning(f"Skipping file {original_file_path} due to renaming error.")
                    continue

            if not file_path.endswith('.fastq'):
                # Check if the target file exists before attempting rename
                if not os.path.exists(file_path):
                    logging.error(f"File not found before renaming: {file_path}. This might be due to a prior error. Skipping.")
                    continue # Skip to the next file

                new_file_path = file_path + '.fastq'
                try:
                    # Check write permissions in the directory before renaming
                    if not os.access(os.path.dirname(file_path), os.W_OK):
                         logging.error(f"Permission denied: Cannot rename file in directory {os.path.dirname(file_path)}. Skipping {original_file_path}.")
                         continue
                    logging.info(f"Renaming {file_path} to {new_file_path}")
                    os.rename(file_path, new_file_path)
                    file_path = new_file_path # Update file_path after successful rename
                except OSError as e:
                    logging.error(f"Failed to rename {file_path} to {new_file_path}: {e}")
                    logging.warning(f"Skipping file {original_file_path} due to renaming error.")
                    continue # Skip to the next file

            # Ensure the final file_path exists before generating the script
            if not os.path.exists(file_path):
                logging.error(f"Final file path {file_path} does not exist before generating job script. Skipping.")
                continue

            job = master_script_generator(file_path, working_directory, job_dir, output_dir, threads, memory,
                                        os.path.join(temporal_directory, 'salmon_index'),
                                        quant_options, email, partition, conda_env, time_limit, module_load_cmd)
            # Run the job
            cmd_launcher = f"sbatch {job}"
            logging.info(f"Submitting job: {cmd_launcher}")
            subprocess.run(cmd_launcher, shell=True, executable='/bin/bash')
        # else: Decompression failed, already logged and skipped

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