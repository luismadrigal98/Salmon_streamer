#!/usr/bin/env python3

"""
Translate Salmon outputs and organize read counts to each allele of each gene per sample.

This program takes salmon outputs (QUANT) and organizes read counts to each allele 
of each gene per sample for a specific cross.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2025-07-01

"""

import sys
import os
import argparse

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.postprocessing_utilities import translate_salmon_outputs

def main(args):
    """
    Main function to translate salmon outputs for a specific cross.
    
    Args:
        args: Parsed command line arguments containing:
            - cross: Cross identifier (e.g., SWB, SF)
            - genes_file: Path to genes mapping file
            - quant_results_file: Path to QUANT results file
            - output_file: Output file path (optional)
    """
    
    # Set default output file if not provided
    if not args.output_file:
        args.output_file = f"Salmon_outputs.IMlines.updated767.{args.cross}.txt"
    
    translate_salmon_outputs(
        cross=args.cross,
        genes_file=args.genes_file,
        quant_results_file=args.quant_results_file,
        output_file=args.output_file
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Translate Salmon outputs for specific cross')
    parser.add_argument('cross', help='Cross identifier (e.g., SWB, SF)')
    parser.add_argument('--genes-file', required=True, 
                       help='Path to genes mapping file (e.g., Genes_to_updated_767_assembly.txt)')
    parser.add_argument('--quant-results-file', required=True,
                       help='Path to QUANT results file (e.g., QUANT_RESULTS_767_vs_SWB)')
    parser.add_argument('--output-file', 
                       help='Output file path (default: Salmon_outputs.IMlines.updated767.[cross].txt)')
    
    args = parser.parse_args()
    main(args)
