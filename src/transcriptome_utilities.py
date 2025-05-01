"""
Transcriptome Utilities for Salmon Streamer

This module provides utility functions for processing transcriptome data,
particularly focused on handling FASTA files, GFF annotations, and generating
combined transcriptomes from reference and alternative genomes.

The module supports the following key operations:
- Reading FASTA sequence files
- Comparing gene positions between reference and lifted-over annotations
- Filtering GFF annotations based on quality metrics
- Extracting gene sequences from genomes using GFF coordinates

Functions
---------
read_fasta : Reads a FASTA file and returns a dictionary of sequences.
generate_position_comparison : Compares gene positions between reference and liftover GFF.
filter_gff_by_quality : Filters liftover GFF based on quality metrics.
extract_gene_sequences : Extracts gene sequences based on filtered GFF files.

Usage
-----
This module is typically used as part of the Salmon Streamer pipeline for
preparing transcriptomes for RNA-seq analysis.
"""

import logging
import sys
import pandas as pd
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)

def read_fasta(fasta_path):
    """
    Reads a FASTA file and returns a dictionary of sequences.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    dict
        Dictionary where keys are sequence headers (without '>')
        and values are the sequences.
    """
    sequences = {}
    current_seq_id = None
    current_seq = []
    try:
        with open(fasta_path, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq_id:
                        sequences[current_seq_id] = ''.join(current_seq)
                    # Handle headers like >Chr_01 pacid=...
                    current_seq_id = line.split()[0][1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            # Add the last sequence
            if current_seq_id:
                sequences[current_seq_id] = ''.join(current_seq)
    except FileNotFoundError:
        logging.error(f"FASTA file not found: {fasta_path}")
        raise
    except Exception as e:
        logging.error(f"Error reading FASTA file {fasta_path}: {e}")
        raise
    logging.info(f"Read {len(sequences)} sequences from {fasta_path}")
    return sequences

def read_gff3_as_dataframe(gff_path):
    """
    Reads a GFF3 file and returns a DataFrame.

    Parameters
    ----------
    gff_path : str
        Path to the GFF3 file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the GFF3 data.
    """
    try:
        gff_df = pd.read_csv(gff_path, sep='\t', header=None, comment='#')
        gff_df.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        # Filter for gene features
        gff_df = gff_df[gff_df['type'] == 'gene']
        # Extract gene_id from attributes
        gff_df['gene_id'] = gff_df['attributes'].str.extract(r'ID=([^;]+)')
        # Preserve only chromosome, start, and end columns, strand, and gene_id
        gff_df = gff_df[['chr', 'start', 'end', 'strand', 'gene_id']]

        return gff_df
    except FileNotFoundError:
        logging.error(f"GFF3 file not found: {gff_path}")
        raise
    except Exception as e:
        logging.error(f"Error reading GFF3 file {gff_path}: {e}")
        raise   

def generate_position_comparison(ref_gene_id, ref_gff_path,
                                alt_genome_id, liftover_gff_path,
                                output_path, preserve_interm = True):
    """
    Compares gene positions between reference and liftover GFF.
    (Based on original script: liftoff_genes_to_767positions.py from John)
    """
    logging.info(f"Generating position comparison file: {output_path}")
    
    try:
        # Read the annotation file as a DataFrame
        # Chr_01  phytozomev13    gene    13982   16715   .       +       .       ID=MgIM767.01G000100.v2.1;Name=MgIM767.01G000100;
        genes_metadata = read_gff3_as_dataframe(ref_gff_path)
        if preserve_interm:
            genes_metadata.to_csv(os.path.join(output_path, f"Gene_pos_{ref_gene_id}"), 
                                sep='\t', index=False, header=True)
        
    except FileNotFoundError:
        logging.error(f"Reference gene position file not found: {ref_gff_path}")
        sys.exit(1)

    try:
        # Read the liftover GFF file
        # Chr_01  liftoff  gene    13982   16715   .       +       .       ID=MgIM767.01G000100.v2.1;Name=MgIM767.01G000100;
        liftover_gff = read_gff3_as_dataframe(liftover_gff_path)
        if preserve_interm:
            liftover_gff.to_csv(os.path.join(output_path, f"Liftover_{alt_genome_id}"), 
                                sep='\t', index=False, header=True)

    except FileNotFoundError:
        logging.error(f"Liftover GFF file not found: {liftover_gff_path}")
        sys.exit(1)

    try:
        # Let's take advantage of pandas to do the comparison
        # Merge the two DataFrames on gene_id

        merged_df = pd.merge(liftover_gff, genes_metadata, on='gene_id', suffixes=((f'_liftoff_{alt_genome_id}', 
                                                                                    f'_ref_{ref_gene_id}')),
                                                                                    how='left')

        # Save the comparison dataframe to the output file
        merged_df.to_csv(output_path, sep='\t', index=False, header=True)

    except IOError as e:
        logging.error(f"Error writing position comparison file {output_path}: {e}")
        sys.exit(1)
    logging.info("Position comparison file generated.")

def filter_gff_by_quality(alt_genome_id, liftover_gff_path, ref_gene_pos_path, original_ref_gff_path,
                        output_gene_list_path, output_ref_gff_path, output_alt_gff_path,
                        cov_threshold=0.9, seqid_threshold=0.9):
    """
    Filters liftover GFF based on quality metrics and chromosome matching.
    Creates filtered GFFs for both reference and alternative genomes.
    (Based on original script 2)
    """
    logging.info(f"Filtering GFF by quality (Coverage>={cov_threshold}, SeqID>={seqid_threshold})")
    chrom_of_gene = {}
    pos_of_gene = {}
    try:
        with open(ref_gene_pos_path, "r") as infile:
            for line in infile:
                cols = line.strip().split('\t')
                if len(cols) >= 5:
                    # Chr_01	13982	16715	+	MgIM767.01G000100.v2.1
                    geneid = "ID=" + cols[4]
                    try:
                        # Attempt to extract chromosome number, handle potential errors
                        chrom_num = int(cols[0].split('_')[1])
                        chrom_of_gene[geneid] = chrom_num
                    except (IndexError, ValueError):
                        logging.warning(f"Could not parse chromosome number for gene {cols[4]} from {cols[0]}. Skipping chromosome check for this gene.")
                        chrom_of_gene[geneid] = None # Indicate chromosome check is not possible
                    pos_of_gene[geneid] = f"{cols[0]}\t{cols[1]}\t{cols[2]}"
    except FileNotFoundError:
        logging.error(f"Reference gene position file not found: {ref_gene_pos_path}")
        sys.exit(1)

    keepers = {}
    outcomes = [0, 0, 0] # [kept, wrong_chrom, low_qual]
    on_alt = False # Flag to track if current gene block should be written to alt GFF

    try:
        with open(liftover_gff_path, "r") as infile, \
            open(output_alt_gff_path, "w") as out_alt, \
            open(output_gene_list_path, "w") as out_genes:

            for line_idx, line in enumerate(infile):
                if line.startswith('##'): # Write header lines to alt GFF
                    out_alt.write(line)
                    continue
                if line.startswith('#'): # Skip other comment lines
                    continue

                cols = line.strip().split('\t')
                if len(cols) < 9: continue

                if cols[2] == "gene":
                    attributes = {k: v for k, v in (pair.split('=', 1) for pair in cols[8].split(';') if '=' in pair)}
                    geneid = attributes.get('ID')
                    if not geneid: continue
                    geneid_key = "ID=" + geneid # Match key format used above

                    try:
                        cover = float(attributes.get('coverage', 0))
                        sid = float(attributes.get('sequence_ID', 0))
                    except (ValueError, TypeError):
                        logging.warning(f"Could not parse coverage/sequence_ID for gene {geneid}. Skipping.")
                        cover, sid = 0, 0

                    ref_chrom_num = chrom_of_gene.get(geneid_key)
                    alt_chrom_str = cols[0] # e.g., "chr1"

                    # Check if chromosome number matches (if possible)
                    chrom_match = False
                    if ref_chrom_num is not None:
                        try:
                            # Assuming alt chrom format like "chr1", "chr01", "chromosome_1" etc.
                            alt_chrom_num_str = ''.join(filter(str.isdigit, alt_chrom_str))
                            if alt_chrom_num_str:
                                alt_chrom_num = int(alt_chrom_num_str)
                                chrom_match = (alt_chrom_num == ref_chrom_num)
                            else:
                                logging.warning(f"Could not extract number from alt chrom {alt_chrom_str} for gene {geneid}")
                        except ValueError:
                            logging.warning(f"Could not parse alt chrom {alt_chrom_str} for gene {geneid}")
                    else:
                        chrom_match = True # Cannot perform check, assume match for quality filter

                    if cover >= cov_threshold and sid >= seqid_threshold and chrom_match:
                        keepers[geneid_key] = 1
                        out_genes.write(f"{geneid}\t{pos_of_gene.get(geneid_key, 'NA\tNA\tNA')}\n")
                        out_alt.write(line)
                        on_alt = True
                        outcomes[0] += 1
                    elif cover >= cov_threshold and sid >= seqid_threshold: # Good quality, wrong chrom
                        on_alt = False
                        outcomes[1] += 1
                    else: # Low quality
                        on_alt = False
                        outcomes[2] += 1
                elif on_alt: # Write subsequent lines (mRNA, exon, CDS) for kept genes
                    out_alt.write(line)

    except FileNotFoundError:
        logging.error(f"Liftover GFF file not found: {liftover_gff_path}")
        sys.exit(1)
    except IOError as e:
        logging.error(f"Error writing filtered alternative GFF {output_alt_gff_path} or gene list {output_gene_list_path}: {e}")
        sys.exit(1)

    logging.info(f"Liftover outcomes (kept/wrong_chrom/low_qual): {outcomes}")
    logging.info(f"Genes kept based on quality and chromosome: {len(keepers)}")

    # Now filter the original reference GFF based on the keepers list
    logging.info(f"Filtering original reference GFF: {original_ref_gff_path}")
    kept_ref_count = 0
    on_ref = False # Flag to track if current gene block should be written to ref GFF
    try:
        with open(original_ref_gff_path, "r") as infile, open(output_ref_gff_path, "w") as out_ref:
            for line_idx, line in enumerate(infile):
                if line.startswith('##'): # Write header lines
                    out_ref.write(line)
                    continue
                if line.startswith('#'): # Skip other comment lines
                    continue

                cols = line.strip().split('\t')
                if len(cols) < 9: continue

                if cols[2] == "gene":
                    attributes = {k: v for k, v in (pair.split('=', 1) for pair in cols[8].split(';') if '=' in pair)}
                    geneid = attributes.get('ID')
                    if not geneid: continue
                    geneid_key = "ID=" + geneid # Match key format

                    if geneid_key in keepers:
                        out_ref.write(line)
                        on_ref = True
                        kept_ref_count += 1
                    else:
                        on_ref = False
                elif on_ref: # Write subsequent lines for kept genes
                    out_ref.write(line)
    except FileNotFoundError:
        logging.error(f"Original reference GFF file not found: {original_ref_gff_path}")
        sys.exit(1)
    except IOError as e:
        logging.error(f"Error writing filtered reference GFF {output_ref_gff_path}: {e}")
        sys.exit(1)

    logging.info(f"Genes written to filtered reference GFF: {kept_ref_count}")
    if kept_ref_count != len(keepers):
        logging.warning(f"Mismatch between keepers ({len(keepers)}) and genes written to ref GFF ({kept_ref_count})")

def extract_gene_sequences(alt_genome_id, ref_genome_id, alt_genome_fasta_path, ref_genome_fasta_path,
                            filtered_alt_gff_path, filtered_ref_gff_path, output_fasta_path):
    """
    Extracts gene sequences based on filtered GFF files and combines them.
    (Based on original script 3)
    """
    logging.info(f"Extracting sequences for combined transcriptome: {output_fasta_path}")

    try:
        logging.info(f"Reading alternative genome: {alt_genome_fasta_path}")
        alt_sequences = read_fasta(alt_genome_fasta_path)
        logging.info(f"Reading reference genome: {ref_genome_fasta_path}")
        ref_sequences = read_fasta(ref_genome_fasta_path)
    except Exception as e:
        logging.error(f"Failed to read genome FASTA files: {e}")
        sys.exit(1)

    extracted_count = 0
    try:
        with open(output_fasta_path, "w") as outfile:
            # Process filtered alternative GFF
            logging.info(f"Processing filtered alternative GFF: {filtered_alt_gff_path}")
            try:
                with open(filtered_alt_gff_path, "r") as infile:
                    for line in infile:
                        if line.startswith('#'): continue
                        cols = line.strip().split('\t')
                        if len(cols) == 9 and cols[2] == "gene":
                            attributes = {k: v for k, v in (pair.split('=', 1) for pair in cols[8].split(';') if '=' in pair)}
                            gene_id = attributes.get('ID')
                            if not gene_id: continue

                            chrom = cols[0]
                            try:
                                start = int(cols[3]) - 1 # GFF is 1-based, Python is 0-based
                                end = int(cols[4])       # End is inclusive in GFF, exclusive in Python slicing
                                strand = cols[6]
                            except ValueError:
                                logging.warning(f"Could not parse coordinates for gene {gene_id} in {filtered_alt_gff_path}. Skipping.")
                                continue

                            if chrom in alt_sequences:
                                seq = alt_sequences[chrom][start:end]
                                if strand == '-':
                                    # Simple reverse complement
                                    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
                                    seq = "".join(complement.get(base, 'N') for base in reversed(seq.upper()))

                                outfile.write(f">{alt_genome_id}_{gene_id}\n")
                                # Write sequence in lines of 60 characters
                                for i in range(0, len(seq), 60):
                                    outfile.write(seq[i:i+60] + '\n')
                                extracted_count += 1
                            else:
                                logging.warning(f"Chromosome/Scaffold '{chrom}' for gene {gene_id} not found in alternative genome FASTA.")
            except FileNotFoundError:
                logging.error(f"Filtered alternative GFF not found: {filtered_alt_gff_path}")
                sys.exit(1)

            # Process filtered reference GFF
            logging.info(f"Processing filtered reference GFF: {filtered_ref_gff_path}")
            try:
                with open(filtered_ref_gff_path, "r") as infile:
                    for line in infile:
                        if line.startswith('#'): continue
                        cols = line.strip().split('\t')
                        if len(cols) == 9 and cols[2] == "gene":
                            attributes = {k: v for k, v in (pair.split('=', 1) for pair in cols[8].split(';') if '=' in pair)}
                            gene_id = attributes.get('ID')
                            if not gene_id: continue

                            chrom = cols[0]
                            try:
                                start = int(cols[3]) - 1 # GFF is 1-based, Python is 0-based
                                end = int(cols[4])       # End is inclusive in GFF, exclusive in Python slicing
                                strand = cols[6]
                            except ValueError:
                                logging.warning(f"Could not parse coordinates for gene {gene_id} in {filtered_ref_gff_path}. Skipping.")
                                continue

                            if chrom in ref_sequences:
                                seq = ref_sequences[chrom][start:end]
                                if strand == '-':
                                    # Simple reverse complement
                                    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
                                    seq = "".join(complement.get(base, 'N') for base in reversed(seq.upper()))

                                outfile.write(f">{ref_genome_id}_{gene_id}\n") # Use ref_genome_id here
                                # Write sequence in lines of 60 characters
                                for i in range(0, len(seq), 60):
                                    outfile.write(seq[i:i+60] + '\n')
                                # Not incrementing extracted_count here as it was counted from alt GFF
                            else:
                                logging.warning(f"Chromosome/Scaffold '{chrom}' for gene {gene_id} not found in reference genome FASTA.")
            except FileNotFoundError:
                logging.error(f"Filtered reference GFF not found: {filtered_ref_gff_path}")
                sys.exit(1)

    except IOError as e:
        logging.error(f"Error writing output FASTA file {output_fasta_path}: {e}")
        sys.exit(1)

    logging.info(f"Finished extracting sequences. Total unique genes processed (from alt GFF): {extracted_count}")