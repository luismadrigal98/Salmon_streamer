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
    Reads a GFF3 file, filters for 'gene' type, and extracts key columns.
    """
    try:
        # Define column names for GFF
        col_names = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        # Read GFF, skipping comments
        gff_df = pd.read_csv(gff_path, sep='\t', header=None, comment='#', names=col_names, low_memory=False)

        # Filter for gene features
        gff_df = gff_df[gff_df['type'] == 'gene'].copy() # Use .copy() to avoid SettingWithCopyWarning

        # Extract gene_id from attributes more robustly
        gff_df['gene_id'] = gff_df['attributes'].str.extract(r'ID=([^;]+)', expand=False)

        # Drop rows where gene_id could not be extracted
        gff_df.dropna(subset=['gene_id'], inplace=True)

        # Keep only necessary columns
        gff_df = gff_df[['chr', 'start', 'end', 'strand', 'gene_id']].reset_index(drop=True)

        return gff_df
    except FileNotFoundError:
        logging.error(f"GFF3 file not found: {gff_path}")
        raise
    except Exception as e:
        logging.error(f"Error reading or processing GFF3 file {gff_path}: {e}")
        raise

def generate_position_comparison(ref_gene_id, ref_gff_path,
                                alt_genome_id, liftover_gff_path,
                                output_path, preserve_interm = False): # Default preserve_interm to False
    """
    Compares gene positions between reference and liftover GFF, keeping all reference genes.
    """
    logging.info(f"Generating position comparison file: {output_path}")

    try:
        # Read the reference annotation file as a DataFrame
        genes_metadata = read_gff3_as_dataframe(ref_gff_path)
        if genes_metadata.empty:
            logging.error(f"No 'gene' features found or extracted from reference GFF: {ref_gff_path}")
            sys.exit(1)
        # Optional: Save intermediate reference gene data
        if preserve_interm:
            output_dir = os.path.dirname(output_path)
            ref_interm_path = os.path.join(output_dir, f"Intermediate_RefGenes_{ref_gene_id}.tsv")
            try:
                genes_metadata.to_csv(ref_interm_path, sep='\t', index=False, header=True)
                logging.info(f"Saved intermediate reference gene data to {ref_interm_path}")
            except IOError as e:
                logging.warning(f"Could not save intermediate reference gene data: {e}")

    except FileNotFoundError:
        logging.error(f"Reference GFF file not found: {ref_gff_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading reference GFF {ref_gff_path}: {e}")
        sys.exit(1)

    try:
        # Read the liftover GFF file
        liftover_gff = read_gff3_as_dataframe(liftover_gff_path)
        if liftover_gff.empty:
            # This might be okay, just means no genes were lifted over
             logging.warning(f"No 'gene' features found or extracted from liftover GFF: {liftover_gff_path}")
             # Create an empty DataFrame with expected columns if liftover is empty to allow merge
             liftover_gff = pd.DataFrame(columns=['chr', 'start', 'end', 'strand', 'gene_id'])

        # Optional: Save intermediate liftover gene data
        if preserve_interm:
            output_dir = os.path.dirname(output_path)
            lift_interm_path = os.path.join(output_dir, f"Intermediate_LiftoverGenes_{alt_genome_id}.tsv")
            try:
                liftover_gff.to_csv(lift_interm_path, sep='\t', index=False, header=True)
                logging.info(f"Saved intermediate liftover gene data to {lift_interm_path}")
            except IOError as e:
                logging.warning(f"Could not save intermediate liftover gene data: {e}")


    except FileNotFoundError:
        logging.error(f"Liftover GFF file not found: {liftover_gff_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading liftover GFF {liftover_gff_path}: {e}")
        sys.exit(1)

    try:
        # Merge genes_metadata (left) with liftover_gff (right)
        # Keep all rows from genes_metadata (reference)
        merged_df = pd.merge(genes_metadata, liftover_gff, on='gene_id',
                            suffixes=(f'_ref_{ref_gene_id}', f'_liftoff_{alt_genome_id}'),
                            how='left') # Keep all reference genes, add liftover info if match

        # Save the comparison dataframe to the output file
        merged_df.to_csv(output_path, sep='\t', index=False, header=True)

    except IOError as e:
        logging.error(f"Error writing position comparison file {output_path}: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred during merging or writing: {e}")
        sys.exit(1)

    logging.info(f"Position comparison file generated: {output_path}")

def extract_chrom_num(chrom_str): # Renamed from _extract_chrom_num
    """Helper function to extract number from chromosome string."""
    if not isinstance(chrom_str, str):
        return None
    # Handle potential float representation if read from CSV with missing values
    try:
        # Remove '.0' if it's like '1.0'
        chrom_str = str(chrom_str).replace('.0', '')
    except Exception:
        return None # Cannot process non-string types reliably

    num_str = ''.join(filter(str.isdigit, chrom_str))
    return int(num_str) if num_str else None

def filter_gff_by_quality(alt_genome_id, ref_genome_id, liftover_gff_path, pos_comparison_path, original_ref_gff_path,
                        output_gene_list_path, output_ref_gff_path, output_alt_gff_path,
                        cov_threshold=0.9, seqid_threshold=0.9):
    """
    Filters liftover GFF based on quality metrics and chromosome matching using the position comparison file.
    Creates filtered GFFs for both reference and alternative genomes.
    Relies on the comparison file keeping all reference genes.
    """
    logging.info(f"Filtering GFF by quality (Coverage>={cov_threshold}, SeqID>={seqid_threshold}) using {pos_comparison_path}")

    # --- Read Position Comparison File ---
    # This file now contains all reference genes, with liftover info added where available.
    liftover_gene_data = {} # Map gene_id -> {alt_chrom: 'chrXX', alt_start: N, alt_end: N, ref_chrom: 'Chr_YY', ...} or None if no liftover
    try:
        comparison_df = pd.read_csv(pos_comparison_path, sep='\t', low_memory=False)
        if 'gene_id' not in comparison_df.columns:
            logging.error(f"'gene_id' column not found in {pos_comparison_path}")
            sys.exit(1)

        # Define expected columns based on the merge strategy (ref is left, liftoff is right)
        ref_chrom_col = f'chr_ref_{ref_genome_id}'
        ref_start_col = f'start_ref_{ref_genome_id}'
        ref_end_col = f'end_ref_{ref_genome_id}'
        alt_chrom_col = f'chr_liftoff_{alt_genome_id}'
        alt_start_col = f'start_liftoff_{alt_genome_id}'
        alt_end_col = f'end_liftoff_{alt_genome_id}'

        required_cols = ['gene_id', ref_chrom_col, ref_start_col, ref_end_col] # Liftover cols might be NaN
        missing_cols = [col for col in required_cols if col not in comparison_df.columns]
        if missing_cols:
            logging.error(f"Expected reference columns missing in {pos_comparison_path}: {', '.join(missing_cols)}")
            sys.exit(1)

        # Check if liftover columns exist, even if they contain NaNs
        if alt_chrom_col not in comparison_df.columns:
             logging.warning(f"Liftover chromosome column '{alt_chrom_col}' not found in {pos_comparison_path}. Chromosome check might be affected.")
             # Add dummy column to prevent key errors later, but checks will likely fail/be skipped
             comparison_df[alt_chrom_col] = pd.NA


        for _, row in comparison_df.iterrows():
            gene_id = row['gene_id']
            # Store reference info always
            ref_data = {
                'ref_chrom': str(row[ref_chrom_col]),
                'ref_start': int(row[ref_start_col]),
                'ref_end': int(row[ref_end_col])
            }
            # Check if liftover data exists (not NaN) for this reference gene
            if pd.notna(row[alt_chrom_col]):
                liftover_gene_data[gene_id] = {
                    **ref_data, # Include reference data
                    'alt_chrom': str(row[alt_chrom_col]),
                    # Add start/end if needed later, but primarily need alt_chrom for check
                }
            else:
                # Reference gene exists, but no liftover match found
                liftover_gene_data[gene_id] = {**ref_data, 'alt_chrom': None} # Mark alt_chrom as None

    except FileNotFoundError:
        logging.error(f"Position comparison file not found: {pos_comparison_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading or processing position comparison file {pos_comparison_path}: {e}")
        sys.exit(1)

    keepers = {} # Map gene_id -> 1 (Genes that passed filters)
    outcomes = [0, 0, 0, 0] # [kept, wrong_chrom, low_qual, no_liftover]
    on_alt = False # Flag for writing alternative GFF features

    # --- Process Liftover GFF (Original File) ---
    # We still need to read the original liftover GFF to get quality attributes (coverage, seqID)
    # which are not in the comparison file.
    liftover_quality_scores = {} # gene_id -> {'coverage': float, 'sequence_ID': float}
    try:
        with open(liftover_gff_path, "r") as infile:
             for line in infile:
                  if line.startswith('#'): continue
                  cols = line.strip().split('\t')
                  if len(cols) < 9 or cols[2] != 'gene': continue

                  attributes = {k: v for k, v in (pair.split('=', 1) for pair in cols[8].split(';') if '=' in pair)}
                  geneid = attributes.get('ID')
                  if not geneid: continue

                  try:
                      cover = float(attributes.get('coverage', 0))
                      sid = float(attributes.get('sequence_ID', 0))
                      liftover_quality_scores[geneid] = {'coverage': cover, 'sequence_ID': sid}
                  except (ValueError, TypeError):
                      logging.warning(f"Could not parse quality scores for gene {geneid} in {liftover_gff_path}. Skipping quality check for this gene.")
                      liftover_quality_scores[geneid] = {'coverage': 0, 'sequence_ID': 0} # Assign default low scores

    except FileNotFoundError:
         # If liftover GFF doesn't exist, no genes can pass quality checks
         logging.warning(f"Original liftover GFF file not found: {liftover_gff_path}. No genes will pass quality filters.")
         # Proceed, but keepers will remain empty.
    except Exception as e:
         logging.error(f"Error reading liftover GFF for quality scores {liftover_gff_path}: {e}")
         # Proceed cautiously, quality checks might be incomplete
         pass


    # --- Determine Keepers based on Comparison Data and Quality Scores ---
    with open(output_gene_list_path, "w") as out_genes:
        # Write header for gene list file
        out_genes.write("GeneID\tRef_Chrom\tRef_Start\tRef_End\tStatus\n")

        for geneid, data in liftover_gene_data.items():
            ref_chrom_str = data['ref_chrom']
            alt_chrom_str = data.get('alt_chrom') # Will be None if no liftover

            status = "Unknown"

            if alt_chrom_str is None:
                # No liftover entry for this reference gene
                outcomes[3] += 1
                status = "NoLiftover"
                # Write to gene list, but don't add to keepers
                out_genes.write(f"{geneid}\t{data['ref_chrom']}\t{data['ref_start']}\t{data['ref_end']}\t{status}\n")
                continue # Cannot pass filters without liftover

            # Get quality scores obtained from the original liftover GFF
            quality = liftover_quality_scores.get(geneid, {'coverage': 0, 'sequence_ID': 0}) # Default if missing
            cover = quality['coverage']
            sid = quality['sequence_ID']

            # Check chromosome match
            ref_chrom_num = extract_chrom_num(ref_chrom_str)
            alt_chrom_num = extract_chrom_num(alt_chrom_str)
            chrom_match = False
            if ref_chrom_num is not None and alt_chrom_num is not None:
                chrom_match = (alt_chrom_num == ref_chrom_num)
            else:
                # Fallback string comparison
                chrom_match = (alt_chrom_str == ref_chrom_str)
                if not chrom_match:
                     logging.debug(f"Chromosome number extraction failed or mismatch for {geneid}: Ref='{ref_chrom_str}', Alt='{alt_chrom_str}'")


            # Apply filters
            if cover >= cov_threshold and sid >= seqid_threshold and chrom_match:
                keepers[geneid] = 1 # Add to keepers list
                outcomes[0] += 1
                status = "Kept"
            elif cover >= cov_threshold and sid >= seqid_threshold: # Good quality, wrong chrom
                outcomes[1] += 1
                status = "WrongChromosome"
            else: # Low quality
                outcomes[2] += 1
                status = "LowQuality"

            # Write all processed genes (with liftover attempt) to the gene list
            out_genes.write(f"{geneid}\t{data['ref_chrom']}\t{data['ref_start']}\t{data['ref_end']}\t{status}\n")


    logging.info(f"Filter outcomes (kept/wrong_chrom/low_qual/no_liftover): {outcomes}")
    logging.info(f"Genes kept based on quality and chromosome: {len(keepers)}")


    # --- Write Filtered Alternative GFF ---
    # Iterate through the original liftover GFF again, writing only features belonging to 'keepers'
    logging.info(f"Writing filtered alternative GFF: {output_alt_gff_path}")
    try:
        with open(liftover_gff_path, "r") as infile, open(output_alt_gff_path, "w") as out_alt:
            current_gene_id = None
            write_current_gene_block = False
            for line in infile:
                if line.startswith('##'):
                    out_alt.write(line)
                    continue
                if line.startswith('#'):
                    continue

                cols = line.strip().split('\t')
                if len(cols) < 9: continue

                feature_type = cols[2]
                attributes = {k: v for k, v in (pair.split('=', 1) for pair in cols[8].split(';') if '=' in pair)}

                if feature_type == "gene":
                    current_gene_id = attributes.get('ID')
                    if current_gene_id in keepers:
                        write_current_gene_block = True
                        out_alt.write(line) # Write the gene line
                    else:
                        write_current_gene_block = False # Don't write this gene or its children
                elif write_current_gene_block and feature_type in ["mRNA", "exon", "CDS"]:
                    # Basic check: write if the parent gene was a keeper
                    # More robust: check Parent attribute if necessary
                    out_alt.write(line)
                else:
                    # Feature type is not gene or child of a kept gene
                    pass
    except FileNotFoundError:
         logging.warning(f"Original liftover GFF {liftover_gff_path} not found. Filtered alternative GFF will be empty.")
         # Create empty file
         open(output_alt_gff_path, 'w').close()
    except Exception as e:
         logging.error(f"Error writing filtered alternative GFF {output_alt_gff_path}: {e}")
         sys.exit(1)


    # --- Filter Original Reference GFF ---
    # This part remains the same: write features from the original reference GFF if the gene ID is in keepers
    logging.info(f"Filtering original reference GFF: {original_ref_gff_path}")
    kept_ref_count = 0
    try:
        with open(original_ref_gff_path, "r") as infile, open(output_ref_gff_path, "w") as out_ref:
            current_gene_id = None
            write_current_gene_block = False
            for line in infile:
                if line.startswith('##'):
                    out_ref.write(line)
                    continue
                if line.startswith('#'):
                    continue

                cols = line.strip().split('\t')
                if len(cols) < 9: continue

                feature_type = cols[2]
                attributes = {k: v for k, v in (pair.split('=', 1) for pair in cols[8].split(';') if '=' in pair)}

                if feature_type == "gene":
                    current_gene_id = attributes.get('ID')
                    if current_gene_id in keepers:
                        write_current_gene_block = True
                        out_ref.write(line) # Write the gene line
                        kept_ref_count += 1
                    else:
                        write_current_gene_block = False
                elif write_current_gene_block and feature_type in ["mRNA", "exon", "CDS"]:
                    out_ref.write(line)
                else:
                    pass # Skip features not belonging to a kept gene

    except FileNotFoundError:
        logging.error(f"Original reference GFF file not found: {original_ref_gff_path}")
        sys.exit(1)
    except IOError as e:
        logging.error(f"Error writing filtered reference GFF {output_ref_gff_path}: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred while processing {original_ref_gff_path}: {e}")
        sys.exit(1)

    logging.info(f"Genes written to filtered reference GFF: {kept_ref_count}")
    if kept_ref_count != len(keepers):
        logging.warning(f"Mismatch between keepers ({len(keepers)}) and genes written to ref GFF ({kept_ref_count}). This might happen if a kept gene ID was missing in the original reference GFF.")

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