#!/usr/bin/env python3
"""
Salmon streamer is a compilation of a bioinformatic pipeline for RNA-seq data analysis.

@author: Luis Javier Madrigal Roca, Paris Veltsos and John K. Kelly.

@date: 2025-04-30

@version: 1.0.0

"""

import argparse
import sys
import os

# Add the src directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

# Import functions from the existing modules
from subprograms.Salmon_runner import main as run_salmon_main
from subprograms.Salmon_output_processor import main as process_salmon_out_main
from subprograms.voom_from_salmon import main as voom_main
from subprograms.AnnTransfer import main as ann_transfer_main
from subprograms.GenerateTranscriptome import main as generate_transcriptome
from subprograms.RunPCA import main as run_pca_qc_main

def main():
    # Create the main parser
    parser = argparse.ArgumentParser(description="Salmon streamer: RNA-seq data analysis pipeline.")
    subparsers = parser.add_subparsers(dest='command', help='Subcommands')
    
    # Create the parser for annotation transfer
    ann_transfer_parser = subparsers.add_parser('AnnTransfer', help='Transfer annotations from reference to alternative genome')

    # Add arguments for annotation transfer
    ann_transfer_parser.add_argument('--target', required=True, help='Target genome in fasta format')
    ann_transfer_parser.add_argument('--reference', required=True, help='Reference genome in fasta format')
    ann_transfer_parser.add_argument('--annotation_gff3', required=True, help='Annotation file in gff3 format')
    ann_transfer_parser.add_argument('--output', required=True, help='Output file for the annotation transfer. This will be a gff3 file')
    ann_transfer_parser.add_argument('--intermediate_dir', required=True, help='Intermediate directory for the annotation transfer')
    ann_transfer_parser.add_argument('--liftoff_path', required=False, default='~/.conda/envs/salmon/bin/liftoff',
                                        help='Path to the liftoff executable')
    ann_transfer_parser.add_argument('--minimap_path', required=False, default='~/.conda/envs/salmon/bin/minimap2',
                                        help='Path to the minimap2 executable')
    #ann_transfer_parser.add_argument('--mm2_options', default='="-a --eqx -N 50 -p 0.5"', help='Options for minimap2')  # THIS IS BROKEN IN LIFTOFF. THERE IS NO WAY OF PARSING THE OPTIONS

    # Transcriptome builder parser
    # --- Transcriptome Generation Subcommand ---
    trans_gen_parser = subparsers.add_parser('GenerateTranscriptome', help='Generate combined transcriptome from liftover results')
    trans_gen_parser.add_argument('--alt-genome-id', required=True, help="Identifier for the alternative genome (e.g., SWB)")
    trans_gen_parser.add_argument('--ref-genome-id', required=True, help="Identifier for the reference genome (e.g., IM767)")
    trans_gen_parser.add_argument('--liftover-gff', required=True, help="Input GFF from AnnTransfer subcommand")
    trans_gen_parser.add_argument('--original-ref-gff', required=True, help="Original reference GFF file used for AnnTransfer")
    trans_gen_parser.add_argument('--alt-genome-fasta', required=True, help="Alternative (target) genome FASTA file")
    trans_gen_parser.add_argument('--ref-genome-fasta', required=True, help="Reference genome FASTA file")
    trans_gen_parser.add_argument('--output-dir', required=True, help="Directory for intermediate and final output files")
    trans_gen_parser.add_argument('--cov-threshold', type=float, default=0.9, help="Minimum coverage threshold for liftover filtering (default: 0.9)")
    trans_gen_parser.add_argument('--seqid-threshold', type=float, default=0.9, help="Minimum sequence ID threshold for liftover filtering (default: 0.9)")
    trans_gen_parser.add_argument('--output-fasta', default='combined_transcriptome.fasta', help="Name for the final output FASTA file (default: combined_transcriptome.fasta)")
    trans_gen_parser.add_argument('--keep-intermediate', action='store_true', help="Keep intermediate files (default: False)")

    # Create the run subcommand parser (from Salmon_runner.py)
    run_parser = subparsers.add_parser('RunSalmonQuant', help='Run Salmon pipeline for RNA-seq quantification')
    
    # Global arguments for runNote
    global_group = run_parser.add_argument_group('Global arguments', 'Directory and path settings')
    global_group.add_argument('-i', '--input', required=True, help='Directory with the fastq files')
    global_group.add_argument('--transcriptome', '--tf', required=True, help='Transcriptome in fasta format')
    global_group.add_argument('--reference', '-r', required=True, help='Reference genome in fasta format')
    global_group.add_argument('--alternative', '-a', required=False, default=None, help='Alternative genome in fasta format')
    global_group.add_argument('--working_directory', '-w', required=True, help='Working directory', default='.')
    global_group.add_argument('--output', '-o', required=True, help='Output directory')
    global_group.add_argument('--temporal_directory', '--temp', required=False, help='Temporal directory', default='./TEMP')
    
    # Naming arguments
    naming_group = run_parser.add_argument_group('Naming arguments', 'Naming settings the chromosomes')
    naming_group.add_argument('--reference_name', '-rn', required=False, help='Name of the reference genome')
    naming_group.add_argument('--alternative_name', '-an', required=False, help='Name of the alternative genome')
    
    # Salmon index arguments
    salmon_index_group = run_parser.add_argument_group('Salmon index arguments', 'Salmon index settings')
    salmon_index_group.add_argument('--salmon_index_options', '-sio', required=False, 
                            help='Salmon index options as a single string', 
                            default='--keepDuplicates -k 31')
    
    # Salmon quant arguments
    salmon_quant_group = run_parser.add_argument_group('Salmon quant arguments', 'Salmon quant settings')
    salmon_quant_group.add_argument('--quant_options', '-qo', required=False,
                            help='Salmon quant options as a single string',
                            default='--noLengthCorrection -l U -p 1')
    
    # Miscellaneous arguments
    misc_group = run_parser.add_argument_group('Miscellaneous arguments', 'Miscellaneous settings')
    misc_group.add_argument('--threads', '-t', type=int, required=False, help='Number of threads to use', default=1)
    misc_group.add_argument('--chrom_level', '-c', required=False, 
                        help='Whether to remove scaffolds and contigs from the genomes', default=True)
    misc_group.add_argument('--memory', '-m', type=int, required=False, 
                        help='Memory to use in the cluster for individual quantification jobs', default=2)
    misc_group.add_argument('--clean', '-cl', action='store_true', 
                        help='Clean the temporal directory after the run', default=False)
    
    # Create the process subcommand parser (from Salmon_output_processor.py)
    process_parser = subparsers.add_parser('ProcessSalmonOut', help='Process Salmon output into a combined table')
    process_parser.add_argument('-o', '--output', required=True, help='Directory with the output of the Salmon pipeline')
    process_parser.add_argument('--result_name', '--rn', default='table.txt', 
                            help='Name of the output file to store the combined results')
    process_parser.add_argument('--mode', '-m', default='cmd', 
                            help='Mode to run the commands. Options are "cmd" for command line and "python" for Python code.')
    process_parser.add_argument('--includes_alternative_genome', action='store_true', 
                            help='If True, the alternative genome was included in the analysis.')
    
    # --- PCA Based QC Subcommand ---
    pca_qc_parser = subparsers.add_parser('PCA_QC', help='Perform PCA based quality control on expression data using an R script.')
    pca_qc_parser.add_argument('--input-data', required=True, help="Path to the input data file (e.g., combined counts table).")
    pca_qc_parser.add_argument('--label-rules-file', required=True, help="Path to a JSON file defining label grouping rules for 'Source' and 'Group' categorization.")
    pca_qc_parser.add_argument('--iqr-multiplier', type=float, default=1.5, help="Multiplier for IQR outlier detection (default: 1.5).")
    pca_qc_parser.add_argument('--output-filtered-data-name', required=True, help="Filename for the output filtered data table (e.g., 'filtered_counts.tsv').")
    pca_qc_parser.add_argument('--genes-as-rows', action='store_true', help="If set, transpose final filtered data to have genes as rows.")
    pca_qc_parser.add_argument('--output-dir', required=True, help="Directory to save output plots and the filtered data file.")
    pca_qc_parser.add_argument('--pc-to-retain', type=int, default=5, help="Number of principal components to retain for analysis (default: 5).")
    pca_qc_parser.add_argument('--rscript-executable', default='~/.conda/envs/PyR/bin/Rscript', help="Path to the Rscript executable (default: Rscript).")

    # Create the voom subcommand parser (from voom_from_salmon.py)
    voom_parser = subparsers.add_parser('Voom', help='Preprocess Salmon output for voom analysis')
    # Add arguments for voom analysis as needed

    voom_parser.add_argument('-o', '--output', required=True, help='Directory with the output of the Salmon pipeline') # PLACEHOLDER
    
    # Parse arguments
    args = parser.parse_args()
    
    # Execute the appropriate subcommand
    if args.command == 'AnnTransfer':
        ann_transfer_main(args)
    elif args.command == 'GenerateTranscriptome':
        # Call the function to generate the transcriptome
        generate_transcriptome(args)
    elif args.command == 'RunSalmonQuant':
        run_salmon_main(args)
    elif args.command == 'ProcessSalmonOut':
        process_salmon_out_main(args)
    elif args.command == 'PCA_QC':
        run_pca_qc_main(args)
    elif args.command == 'Voom':
        voom_main()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()