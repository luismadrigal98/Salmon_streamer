#!/usr/bin/env python3

"""
This program will handle the annotation transfer from the reference genome to the alternative genome using
liftoff.

@author: Luis Javier Madrigal-Roca & John K. Kelly
@date: 2025-04-30
@version: 1.0.0

"""

import logging
import sys
import subprocess
import os

## Setting up logging (fancy style)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def main(args):
    """
    Run Liftoff to transfer genome annotations from a reference to a target genome.
    This function extracts arguments, validates input/output paths, and executes
    the Liftoff tool to transfer annotations. It creates necessary directories
    if they don't exist and checks for successful execution.
    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments containing:
        - target: Path to the target genome file
        - reference: Path to the reference genome file
        - annotation_gff3: Path to the annotation file in GFF3 format
        - output: Path to the output file
        - intermediate_dir: Directory for intermediate files
        - liftoff_path: Path to the Liftoff executable
        - minimap_path: Path to the minimap executable
        - mm2_options: Options for minimap2
    Returns:
    --------
    None
        The function exits with status code 1 if any input files are missing or
        if the Liftoff command fails.
    Side Effects:
    -------------
    - Creates output and intermediate directories if they don't exist
    - Writes transferred annotations to the output file
    - Logs information, warnings, and errors
    """

    # Decompress the arguments
    
    # Basic input outputs
    target = args.target
    reference = args.reference
    annotation_gff3 = args.annotation_gff3
    output = args.output
    intermediate_dir = args.intermediate_dir
    
    # Executables
    liftoff_path = args.liftoff_path
    minimap_path = args.minimap_path
    mm2_options = args.mm2_options

    output_dir = os.path.dirname(output)

    # Check input, ouptput and intermediate directories
    if not os.path.isfile(target):
        logging.error(f'The target file does not exist: {target}')
        sys.exit(1)
    if not os.path.isfile(reference):
        logging.error(f'The reference file does not exist: {reference}')
        sys.exit(1)
    if not os.path.isfile(annotation_gff3):
        logging.error(f'The annotation file does not exist: {annotation_gff3}')
        sys.exit(1)
    if not os.path.isdir(output_dir):
        logging.warning('The output directory does not exist')
        logging.info('Creating the output directory')
        os.makedirs(output_dir, exist_ok=True)       
    if not os.path.isdir(intermediate_dir):
        logging.warning('The intermediate directory does not exist')
        logging.info('Creating the intermediate directory')
        os.makedirs(intermediate_dir, exist_ok=True)

    # Building the liftoff command
    cmd = f"{liftoff_path} -g {annotation_gff3} -o {output} -dir {intermediate_dir} " + \
            f"-mm2_options {mm2_options} -m {minimap_path} {target} {reference}"
    
    # Run the command
    logging.info(f"Running liftoff command: {cmd}")
    subprocess.run(cmd, shell=True, executable='/bin/bash')

    # Check if the command was successful
    if subprocess.run(cmd, shell=True, executable='/bin/bash').returncode != 0:
        logging.error('Liftoff failed')
        sys.exit(1)
    else:
        logging.info("Liftoff finished successfully")