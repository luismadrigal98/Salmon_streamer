#!/usr/bin/env python3

"""
This program will handle the transcriptome generation from the reference genome and the alternative genome.

@author: Luis Javier Madrigal-Roca & John K. Kelly
@date: 2025-05-01
@version: 1.0.0

"""

import argparse
import logging
import os
import sys

## Setting up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def main(args):
    
    # Decompress the arguments
    alt_genome_id = args.alt_genome_id
    ref_genome_id = args.ref_genome_id
    liftover_gff = args.liftover_gff
    ref_gene_pos = args.ref_gene_pos
    original_ref_gff = args.original_ref_gff
    alt_genome_fasta = args.alt_genome_fasta
    ref_genome_fasta = args.ref_genome_fasta
    output_dir = args.output_dir
    cov_threshold = args.cov_threshold
    seqid_threshold = args.seqid_threshold
    output_fasta = args.output_fasta

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define intermediate and final file paths
    pos_comparison_file = os.path.join(output_dir, f"{ref_genome_id}_vs_{alt_genome_id}_positions.txt")
    gene_list_file = os.path.join(output_dir, f"genes_in_{ref_genome_id}x{alt_genome_id}_cross.txt")
    filtered_ref_gff = os.path.join(output_dir, f"REF_{ref_genome_id}x{alt_genome_id}_cross.gff3")
    filtered_alt_gff = os.path.join(output_dir, f"ALT_{alt_genome_id}x{ref_genome_id}_cross.gff3")
    final_output_fasta_path = os.path.join(output_dir, output_fasta)

# --- Step 1: Generate Position Comparison (Optional but good for QC) ---
    try:
        generate_position_comparison(alt_genome_id, liftover_gff, ref_gene_pos, pos_comparison_file)
    except Exception as e:
        logging.error(f"Failed during position comparison generation: {e}")
        sys.exit(1)

    # --- Step 2: Filter GFFs by Quality ---
    try:
        filter_gff_by_quality(alt_genome_id, liftover_gff, ref_gene_pos, original_ref_gff,
                                gene_list_file, filtered_ref_gff, filtered_alt_gff,
                                cov_threshold, seqid_threshold)
    except Exception as e:
        logging.error(f"Failed during GFF quality filtering: {e}")
        sys.exit(1)

    # --- Step 3: Extract Sequences ---
    try:
        extract_gene_sequences(alt_genome_id, alt_genome_fasta, ref_genome_fasta,
                                filtered_alt_gff, filtered_ref_gff, final_output_fasta_path)
    except Exception as e:
        logging.error(f"Failed during sequence extraction: {e}")
        sys.exit(1)

    logging.info("Transcriptome generation process completed successfully.")
    logging.info(f"Final combined transcriptome FASTA: {final_output_fasta_path}")

if __name__ == '__main__':
    main()