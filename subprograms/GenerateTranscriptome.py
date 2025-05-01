#!/usr/bin/env python3

"""
This program will handle the transcriptome generation from the reference genome and the alternative genome.

@author: Luis Javier Madrigal-Roca & John K. Kelly
@date: 2025-05-01
@version: 1.0.0

"""

def main(args):
    
    # Decompress the arguments
    alt_genome_id = args.alt_genome_id
    liftover_gff = args.liftover_gff
    ref_gene_pos = args.ref_gene_pos
    original_ref_gff = args.original_ref_gff
    alt_genome_fasta = args.alt_genome_fasta
    ref_genome_fasta = args.ref_genome_fasta
    output_dir = args.output_dir
    cov_threshold = args.cov_threshold
    seqid_threshold = args.seqid_threshold
    output_fasta = args.output_fasta

if __name__ == '__main__':
    main()