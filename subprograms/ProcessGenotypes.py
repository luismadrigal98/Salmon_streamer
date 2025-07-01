#!/usr/bin/env python3

"""
Process genotypes from transcript mapping data.

This program handles the genotype calling pipeline including transcript mapping,
error rate estimation, and genotype posterior probability calculation.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2025-07-01

"""

import sys
import os
import argparse

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.postprocessing_utilities import process_genotypes_pipeline

def main(args):
    """
    Main function to process genotypes from transcript mapping data.
    
    Args:
        args: Parsed command line arguments containing genotype processing parameters
    """
    
    process_genotypes_pipeline(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process genotypes from transcript mapping data')
    parser.add_argument('--cross', required=True,
                       help='Cross identifier (e.g., SWB, SF, 1034)')
    parser.add_argument('--allele-counts-file', required=True,
                       help='Path to allele counts file (e.g., allele_counts.combined.updated767.CROSS.txt)')
    parser.add_argument('--samples-file', required=True,
                       help='Path to samples file (e.g., CROSS.combined.samples.txt)')
    parser.add_argument('--genes-file', required=True,
                       help='Path to genes mapping file (e.g., Genes_to_updated_767_assembly.txt)')
    parser.add_argument('--min-parental-lines', type=int, default=5,
                       help='Minimum parental lines called (default: 5)')
    parser.add_argument('--mapping-threshold', type=float, default=0.95,
                       help='Minimum frequency of mapping to correct allele (default: 0.95)')
    parser.add_argument('--min-reads-per-plant', type=int, default=100000,
                       help='Minimum reads per plant (default: 100000)')
    parser.add_argument('--min-reads-per-call', type=int, default=6,
                       help='Minimum reads to count an F2 call (default: 6)')
    parser.add_argument('--f2-fraction-threshold', type=float, default=0.5,
                       help='Fraction of F2s required (default: 0.5)')
    parser.add_argument('--homozygous-threshold', type=float, default=0.95,
                       help='Homozygous calling threshold (default: 0.95)')
    parser.add_argument('--het-maf', type=float, default=0.25,
                       help='Heterozygous calling MAF threshold (default: 0.25)')
    parser.add_argument('--output-dir', default='.',
                       help='Output directory (default: current directory)')
    
    args = parser.parse_args()
    main(args)
