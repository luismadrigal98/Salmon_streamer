#!/usr/bin/env python3

"""
Calculate raw reads per plant from Salmon outputs.

This program takes Salmon outputs and determines total reads per plant
across all crosses.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2025-07-01

"""

import sys
import os
import argparse

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.postprocessing_utilities import calculate_raw_reads_per_plant

def main(args):
    """
    Main function to calculate raw reads per plant.
    
    Args:
        args: Parsed command line arguments containing:
            - salmon_files: List of Salmon output files
            - output_file: Output file path
    """
    
    calculate_raw_reads_per_plant(
        salmon_files=args.salmon_files,
        output_file=args.output_file
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate raw reads per plant from Salmon outputs')
    parser.add_argument('--salmon-files', nargs='+', required=True,
                       help='List of Salmon output files (e.g., Salmon_outputs.IMlines.updated767.SWB.txt)')
    parser.add_argument('--output-file', default='raw.reads.per.plant.txt',
                       help='Output file path (default: raw.reads.per.plant.txt)')
    
    args = parser.parse_args()
    main(args)
