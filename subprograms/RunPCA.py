#!/usr/bin/env python3
"""
Wrapper script to execute the PCA_based_QC.R script.
"""
import argparse
import subprocess
import sys
import os
import logging
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def main(args):
    # Construct the command to run the R script
    # Ensure the path to PCA_based_QC.R is correct relative to SalmonStreamer.py
    # For simplicity, assuming SalmonStreamer.py is run from the project root.
    r_script_path = os.path.join(os.path.dirname(__file__), 'subprograms/PCA_based_QC.R')
    if not os.path.exists(r_script_path):
        logging.error(f"R script PCA_based_QC.R not found at expected locations.")
        sys.exit(1)

# Prepare arguments for the R script
    r_args = [
        args.rscript_executable,
        r_script_path,
        args.input_data,
        args.label_rules_file,
        str(args.iqr_multiplier),
        args.output_filtered_data_name,
        str(args.genes_as_rows).upper(), # Pass TRUE/FALSE as string
        args.output_dir,
        str(args.pc_to_retain)
    ]

    try:
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        
        result = subprocess.run(r_args, cwd=project_root, check=True, capture_output=True, text=True)
        logging.info("R script executed successfully.")
        if result.stdout:
            logging.info("R script stdout:\n" + result.stdout)
        if result.stderr:
            # R often prints informational messages to stderr
            logging.info("R script stderr:\n" + result.stderr)
            
    except subprocess.CalledProcessError as e:
        logging.error(f"R script failed with exit code {e.returncode}")
        if e.stdout:
            logging.error("R script stdout:\n" + e.stdout)
        if e.stderr:
            logging.error("R script stderr:\n" + e.stderr)
        sys.exit(1)
    except FileNotFoundError:
        logging.error(f"Rscript executable not found at '{args.rscript_executable}'. Please ensure Rscript is in your PATH or provide the full path.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        sys.exit(1)

