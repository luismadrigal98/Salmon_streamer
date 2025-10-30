#!/usr/bin/env python3

"""
Extract Transcriptome from Single Reference Genome

This program extracts transcript sequences from a reference genome using GFF3 annotations,
creating a transcriptome FASTA file that is fully compatible with Salmon quantification.

Unlike the dual-genome GenerateTranscriptome, this tool works with a single reference genome
and extracts mature mRNA sequences by concatenating exons for each transcript.

Features:
- Extracts transcript sequences (mRNA) by concatenating exons
- Handles both gene and transcript-level features
- Supports strand-specific sequence extraction with reverse complement
- Compatible with Salmon quantification requirements
- Flexible feature type selection (mRNA, transcript, etc.)
- Robust GFF3 parsing with error handling

@author: Luis Javier Madrigal-Roca & John K. Kelly  
@date: 2025-10-30
@version: 1.0.0
"""

import argparse
import logging
import os
import sys
from collections import defaultdict

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.transcriptome_utilities import read_fasta

# Setting up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def reverse_complement(sequence):
    """
    Return the reverse complement of a DNA sequence.
    
    Parameters
    ----------
    sequence : str
        DNA sequence string
        
    Returns
    -------
    str
        Reverse complement of the input sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 
                  'g': 'c', 'n': 'n'}
    
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

def parse_gff3_attributes(attribute_string):
    """
    Parse GFF3 attributes string into a dictionary.
    
    Parameters
    ----------
    attribute_string : str
        GFF3 attributes column (column 9)
        
    Returns
    -------
    dict
        Dictionary of attribute key-value pairs
    """
    attributes = {}
    for pair in attribute_string.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attributes[key.strip()] = value.strip()
    return attributes

def read_gff3_transcripts(gff_path, transcript_types=['mRNA'], gene_types=['gene']):
    """
    Read GFF3 file and extract transcript and exon information.
    
    Parameters
    ----------
    gff_path : str
        Path to the GFF3 file
    transcript_types : list
        List of feature types to consider as transcripts (default: ['mRNA'])
    gene_types : list
        List of feature types to consider as genes (default: ['gene'])
        
    Returns
    -------
    tuple
        (genes_dict, transcripts_dict, exons_dict)
        - genes_dict: gene_id -> {chr, start, end, strand, attributes}
        - transcripts_dict: transcript_id -> {chr, start, end, strand, gene_id, attributes}
        - exons_dict: transcript_id -> [list of exon dicts with chr, start, end, strand]
    """
    genes = {}
    transcripts = {}
    exons = defaultdict(list)
    
    logging.info(f"Reading GFF3 file: {gff_path}")
    
    try:
        with open(gff_path, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                line = line.strip()
                
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                
                cols = line.split('\t')
                if len(cols) != 9:
                    logging.warning(f"Line {line_count}: Invalid GFF3 format (expected 9 columns, got {len(cols)})")
                    continue
                
                seqid, source, feature_type, start, end, score, strand, phase, attributes_str = cols
                
                try:
                    start = int(start)
                    end = int(end)
                except ValueError:
                    logging.warning(f"Line {line_count}: Invalid coordinates for {feature_type}")
                    continue
                
                attributes = parse_gff3_attributes(attributes_str)
                
                # Process genes
                if feature_type in gene_types:
                    gene_id = attributes.get('ID')
                    if gene_id:
                        genes[gene_id] = {
                            'chr': seqid,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'attributes': attributes
                        }
                
                # Process transcripts/mRNAs
                elif feature_type in transcript_types:
                    transcript_id = attributes.get('ID')
                    parent_id = attributes.get('Parent')
                    
                    if transcript_id:
                        transcripts[transcript_id] = {
                            'chr': seqid,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'gene_id': parent_id,
                            'attributes': attributes
                        }
                
                # Process exons
                elif feature_type == 'exon':
                    parent_id = attributes.get('Parent')
                    if parent_id:
                        # Handle cases where Parent might reference multiple transcripts
                        parent_ids = parent_id.split(',')
                        for pid in parent_ids:
                            pid = pid.strip()
                            exons[pid].append({
                                'chr': seqid,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'attributes': attributes
                            })
                
    except FileNotFoundError:
        logging.error(f"GFF3 file not found: {gff_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading GFF3 file {gff_path}: {e}")
        sys.exit(1)
    
    logging.info(f"Parsed GFF3: {len(genes)} genes, {len(transcripts)} transcripts, {len(exons)} transcript-exon sets")
    
    return genes, transcripts, exons

def extract_transcript_sequence(genome_sequences, transcript_info, exon_list):
    """
    Extract transcript sequence by concatenating exons.
    
    Parameters
    ----------
    genome_sequences : dict
        Dictionary of chromosome sequences from FASTA
    transcript_info : dict
        Transcript information dictionary
    exon_list : list
        List of exon dictionaries for this transcript
        
    Returns
    -------
    str or None
        Transcript sequence, or None if extraction failed
    """
    if not exon_list:
        return None
    
    chr_name = transcript_info['chr']
    strand = transcript_info['strand']
    
    if chr_name not in genome_sequences:
        logging.warning(f"Chromosome '{chr_name}' not found in genome sequences")
        return None
    
    chromosome_seq = genome_sequences[chr_name]
    
    # Sort exons by start position
    sorted_exons = sorted(exon_list, key=lambda x: x['start'])
    
    # Extract and concatenate exon sequences
    transcript_sequence = ""
    for exon in sorted_exons:
        start = exon['start'] - 1  # Convert to 0-based indexing
        end = exon['end']          # GFF3 end is inclusive, Python slicing is exclusive
        
        if start < 0 or end > len(chromosome_seq):
            logging.warning(f"Exon coordinates ({start+1}-{end}) out of bounds for chromosome {chr_name} (length: {len(chromosome_seq)})")
            continue
        
        exon_seq = chromosome_seq[start:end]
        transcript_sequence += exon_seq
    
    # Apply reverse complement if on negative strand
    if strand == '-':
        transcript_sequence = reverse_complement(transcript_sequence)
    
    return transcript_sequence

def main(args):
    """
    Main function to extract transcriptome from reference genome.
    """
    # Parse arguments
    genome_fasta = args.genome_fasta
    gff_file = args.gff_file
    output_fasta = args.output_fasta
    transcript_types = args.transcript_types.split(',') if args.transcript_types else ['mRNA']
    gene_types = args.gene_types.split(',') if args.gene_types else ['gene']
    id_prefix = args.id_prefix
    include_gene_id = args.include_gene_id
    min_length = args.min_length
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(os.path.abspath(output_fasta)), exist_ok=True)
    
    # Read genome sequences
    logging.info(f"Reading genome FASTA: {genome_fasta}")
    try:
        genome_sequences = read_fasta(genome_fasta)
    except Exception as e:
        logging.error(f"Failed to read genome FASTA file: {e}")
        sys.exit(1)
    
    # Read GFF3 annotations
    genes, transcripts, exons = read_gff3_transcripts(gff_file, transcript_types, gene_types)
    
    if not transcripts:
        logging.error("No transcripts found in GFF3 file. Check your transcript types.")
        sys.exit(1)
    
    # Extract transcript sequences
    logging.info("Extracting transcript sequences...")
    extracted_count = 0
    skipped_count = 0
    short_count = 0
    
    try:
        with open(output_fasta, 'w') as outfile:
            for transcript_id, transcript_info in transcripts.items():
                # Get exons for this transcript
                transcript_exons = exons.get(transcript_id, [])
                
                if not transcript_exons:
                    logging.warning(f"No exons found for transcript {transcript_id}")
                    skipped_count += 1
                    continue
                
                # Extract sequence
                sequence = extract_transcript_sequence(genome_sequences, transcript_info, transcript_exons)
                
                if sequence is None:
                    logging.warning(f"Failed to extract sequence for transcript {transcript_id}")
                    skipped_count += 1
                    continue
                
                # Check minimum length
                if len(sequence) < min_length:
                    logging.debug(f"Transcript {transcript_id} too short ({len(sequence)} bp), skipping")
                    short_count += 1
                    continue
                
                # Create FASTA header
                header_parts = []
                if id_prefix:
                    header_parts.append(id_prefix)
                
                header_parts.append(transcript_id)
                
                if include_gene_id and transcript_info.get('gene_id'):
                    header_parts.append(f"gene:{transcript_info['gene_id']}")
                
                header = '_'.join(header_parts) if len(header_parts) > 1 else header_parts[0]
                
                # Write to FASTA
                outfile.write(f">{header}\n")
                
                # Write sequence in 80-character lines
                for i in range(0, len(sequence), 80):
                    outfile.write(sequence[i:i+80] + '\n')
                
                extracted_count += 1
                
                if extracted_count % 1000 == 0:
                    logging.info(f"Extracted {extracted_count} transcripts...")
    
    except IOError as e:
        logging.error(f"Error writing output FASTA file {output_fasta}: {e}")
        sys.exit(1)
    
    # Summary
    logging.info(f"Transcriptome extraction completed!")
    logging.info(f"Output file: {output_fasta}")
    logging.info(f"Successfully extracted: {extracted_count} transcripts")
    if skipped_count > 0:
        logging.info(f"Skipped (no exons/failed): {skipped_count} transcripts")
    if short_count > 0:
        logging.info(f"Skipped (too short): {short_count} transcripts")
    
    total_transcripts = len(transcripts)
    success_rate = (extracted_count / total_transcripts * 100) if total_transcripts > 0 else 0
    logging.info(f"Success rate: {success_rate:.1f}% ({extracted_count}/{total_transcripts})")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract transcriptome sequences from reference genome using GFF3 annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python ExtractTranscriptome.py --genome genome.fasta --gff annotation.gff3 --output transcriptome.fasta
  
  # Custom transcript types
  python ExtractTranscriptome.py --genome genome.fasta --gff annotation.gff3 --output transcriptome.fasta --transcript-types mRNA,transcript
  
  # Include gene IDs in headers
  python ExtractTranscriptome.py --genome genome.fasta --gff annotation.gff3 --output transcriptome.fasta --include-gene-id
        """
    )
    
    # Required arguments
    parser.add_argument('--genome-fasta', '--genome', required=True,
                        help='Reference genome FASTA file')
    parser.add_argument('--gff-file', '--gff', required=True,
                        help='GFF3 annotation file')
    parser.add_argument('--output-fasta', '--output', required=True,
                        help='Output transcriptome FASTA file')
    
    # Optional arguments
    parser.add_argument('--transcript-types', default='mRNA',
                        help='Comma-separated list of GFF3 feature types to treat as transcripts (default: mRNA)')
    parser.add_argument('--gene-types', default='gene',
                        help='Comma-separated list of GFF3 feature types to treat as genes (default: gene)')
    parser.add_argument('--id-prefix', default=None,
                        help='Prefix to add to transcript IDs in output FASTA headers')
    parser.add_argument('--include-gene-id', action='store_true',
                        help='Include gene ID in FASTA headers')
    parser.add_argument('--min-length', type=int, default=1,
                        help='Minimum transcript length to include (default: 1)')
    
    args = parser.parse_args()
    main(args)