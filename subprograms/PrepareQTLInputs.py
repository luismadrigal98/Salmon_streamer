#!/usr/bin/env python3

"""
Prepare all inputs for QTL analysis including genotype and phenotype files.

This program creates properly formatted R/qtl input files from processed
genotype and phenotype data.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2025-07-01

"""

import sys
import os
import argparse

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.postprocessing_utilities import prepare_qtl_inputs_pipeline

def main(args):
    """
    Main function to prepare QTL analysis inputs.
    
    Args:
        args: Parsed command line arguments containing QTL input preparation parameters
    """
    
    prepare_qtl_inputs_pipeline(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare inputs for QTL analysis')
    parser.add_argument('--cross', required=True,
                       help='Cross identifier (e.g., SWB, SF, 1034)')
    parser.add_argument('--genotype-pp-file', required=True,
                       help='Path to genotype posterior probabilities file (e.g., CROSS.F2_geno_PP.txt)')
    parser.add_argument('--estimates-files', nargs='+', required=True,
                       help='List of estimates files for each chromosome')
    parser.add_argument('--genes-by-cross-file', required=True,
                       help='Path to genes by cross file')
    parser.add_argument('--phenotype-group', required=True,
                       help='Phenotype group identifier (e.g., gene.group1)')
    parser.add_argument('--phenotype-files-dir', default='pfiles',
                       help='Directory containing phenotype files (default: pfiles)')
    parser.add_argument('--output-dir', default='rQTL_files',
                       help='Output directory for R/qtl files (default: rQTL_files)')
    parser.add_argument('--chromosomes', nargs='+', type=int, default=list(range(1, 15)),
                       help='List of chromosome numbers (default: 1-14)')
    
    args = parser.parse_args()
    main(args)
