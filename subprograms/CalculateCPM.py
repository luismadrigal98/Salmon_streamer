#!/usr/bin/env python3

"""
Standardize raw reads to CPM (Counts Per Million) for PCA analysis.

This program takes raw read counts and standardizes them to CPM,
filtering genes based on minimum CPM threshold.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2025-07-01

"""

import sys
import os
import argparse

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.postprocessing_utilities import calculate_cpm_for_pca

def main(args):
    """
    Main function to calculate CPM for PCA analysis.
    
    Args:
        args: Parsed command line arguments containing:
            - raw_reads_file: Path to raw reads per plant file
            - salmon_files: List of Salmon output files
            - cpm_min: Minimum CPM threshold
            - output_file: Output file path
    """
    
    calculate_cpm_for_pca(
        raw_reads_file=args.raw_reads_file,
        salmon_files=args.salmon_files,
        cpm_min=args.cpm_min,
        output_file=args.output_file
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate CPM for PCA analysis')
    parser.add_argument('--raw-reads-file', required=True,
                       help='Path to raw reads per plant file (e.g., raw.reads.per.plant.txt)')
    parser.add_argument('--salmon-files', nargs='+', required=True,
                       help='List of Salmon output files (e.g., Salmon_outputs.IMlines.updated767.SWB.txt)')
    parser.add_argument('--cpm-min', type=float, default=5.0,
                       help='Minimum CPM threshold (default: 5.0)')
    parser.add_argument('--output-file', default='RawSamples_forPCA',
                       help='Output file path (default: RawSamples_forPCA)')
    
    args = parser.parse_args()
    main(args)
