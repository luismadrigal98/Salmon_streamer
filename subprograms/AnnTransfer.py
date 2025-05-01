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
    
    # Executables - Expand user paths here
    liftoff_path = os.path.expanduser(args.liftoff_path)
    minimap_path = os.path.expanduser(args.minimap_path)
    # mm2_options = args.mm2_options # This should contain the string like "-a --eqx ..."

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
    # Only create the output directory if it's specified and doesn't exist
    if output_dir and not os.path.isdir(output_dir):
        logging.warning(f'The output directory does not exist: {output_dir}')
        logging.info(f'Creating the output directory: {output_dir}')
        os.makedirs(output_dir, exist_ok=True)
    if not os.path.isdir(intermediate_dir):
        logging.warning('The intermediate directory does not exist')
        logging.info('Creating the intermediate directory')
        os.makedirs(intermediate_dir, exist_ok=True)

    # Building the liftoff command as a list for shell=False
    # Pass the mm2_options string as the argument value
    cmd_list = [
        liftoff_path,
        '-g', annotation_gff3,
        '-o', output,
        '-dir', intermediate_dir,
        #'-mm2_options', mm2_options, # Pass the string "-a --eqx ..." directly
        '-m', minimap_path,
        target,
        reference
    ]

    # Run the command and capture the result
    # Use shell=False (default) and pass the command as a list
    logging.info(f"Running liftoff command: {' '.join(cmd_list)}")
    try:
        # Use check=True to automatically raise an error for non-zero exit codes
        # Use stdout=subprocess.PIPE, stderr=subprocess.PIPE for Python < 3.7
        # Use universal_newlines=True instead of text=True for Python < 3.7
        result = subprocess.run(cmd_list, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, executable=None)
        logging.info("Liftoff finished successfully")
        # Log stdout and stderr even on success (might contain warnings)
        if result.stdout:
            logging.info("Liftoff stdout:\n" + result.stdout)
        if result.stderr:
            logging.warning("Liftoff stderr:\n" + result.stderr)

    except subprocess.CalledProcessError as e:
        # Log detailed error information if liftoff fails
        logging.error('Liftoff failed')
        logging.error(f"Command: {' '.join(e.cmd)}")
        logging.error(f"Return code: {e.returncode}")
        # Access stdout/stderr from the exception object
        if e.stdout:
            logging.error("stdout:\n" + e.stdout)
        if e.stderr:
            logging.error("stderr:\n" + e.stderr)
        sys.exit(1)
    except FileNotFoundError:
        # Handle case where liftoff executable itself is not found
        # Log the *expanded* path here for clarity
        logging.error(f"Liftoff executable not found at the expanded path: {liftoff_path}")
        sys.exit(1)
    except Exception as e:
        # Catch any other unexpected errors during subprocess execution
        logging.error(f"An unexpected error occurred while running liftoff: {e}")
        sys.exit(1)