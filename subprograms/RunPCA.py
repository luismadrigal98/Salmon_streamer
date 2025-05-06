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
    # PCA_based_QC.R is in the same directory as this RunPCA.py script
    r_script_path = os.path.join(os.path.dirname(__file__), 'PCA_based_QC.R')
    if not os.path.exists(r_script_path):
        # Fallback if the script is not found, though the above should be correct
        logging.error(f"R script PCA_based_QC.R not found at {r_script_path}.")
        # Attempt path relative to project root if SalmonStreamer.py is in root
        # This part might be redundant if the first path is correct.
        alt_r_script_path = os.path.join("subprograms", "PCA_based_QC.R")
        if os.path.exists(alt_r_script_path):
            r_script_path = alt_r_script_path
            logging.info(f"Found R script at alternative path: {r_script_path}")
        else:
            logging.error(f"R script PCA_based_QC.R not found at expected locations.")
            sys.exit(1)

    r_script_executable_path = os.path.expanduser(args.rscript_executable)

# Prepare arguments for the R script
    r_args = [
        r_script_executable_path,
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

