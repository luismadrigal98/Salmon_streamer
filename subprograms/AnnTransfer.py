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

# Import necessary utilities
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.gff3_transference import *

## Setting up logging
logging.basicConfig(level=logging.INFO)

def main(args):

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

    # Check input, ouptput and intermediate directories
    if not os.path.isdir(target):
        logging.error('The target directory does not exist')
        sys.exit(1)
    if not os.path.isdir(reference):
        logging.error('The reference directory does not exist')
        sys.exit(1)
    if not os.path.isdir(annotation_gff3):
        logging.error('The annotation directory does not exist')
        sys.exit(1)
    if not os.path.isdir(output):
        logging.warning('The output directory does not exist')
        logging.info('Creating the output directory')
        os.makedirs(output, exist_ok=True)       
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