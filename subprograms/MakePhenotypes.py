#!/usr/bin/env python3

"""
Generate phenotype files from expression data for QTL analysis.

This program processes expression data and generates phenotype files
suitable for QTL analysis, including Box-Cox transformation.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2025-07-01

"""

import sys
import os
import argparse

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.postprocessing_utilities import make_phenotype_files_pipeline

def main(args):
    """
    Main function to generate phenotype files from expression data.
    
    Args:
        args: Parsed command line arguments containing phenotype file generation parameters
    """
    
    make_phenotype_files_pipeline(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate phenotype files from expression data')
    parser.add_argument('--genes-by-cross-file', required=True,
                       help='Path to genes by cross file (e.g., Genes_by_cross.txt)')
    parser.add_argument('--total-reads-files', nargs='+', required=True,
                       help='List of total reads files for each cross (e.g., Total_reads_byplant.SF)')
    parser.add_argument('--readcounts-files', nargs='+', required=True,
                       help='List of readcounts files for each cross')
    parser.add_argument('--f2-lists', nargs='+', required=True,
                       help='List of F2 inclusion files (e.g., SF.included.f2s.txt)')
    parser.add_argument('--crosses', nargs='+', required=True,
                       help='List of cross identifiers')
    parser.add_argument('--im-lines', nargs='+', 
                       default=["62","155","444","502","541","664","909","1034","1192"],
                       help='List of inbred line identifiers')
    parser.add_argument('--output-dir', default='pfiles',
                       help='Output directory for phenotype files (default: pfiles)')
    parser.add_argument('--summary-file', default='f2summaries_by_gene.txt',
                       help='Summary output file (default: f2summaries_by_gene.txt)')
    
    args = parser.parse_args()
    main(args)
