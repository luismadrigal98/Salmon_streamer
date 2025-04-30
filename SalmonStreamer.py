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

def main():
    # Create the main parser
    parser = argparse.ArgumentParser(description="Salmon streamer: RNA-seq data analysis pipeline.")
    subparsers = parser.add_subparsers(dest='command', help='Subcommands')
    
    # Create the run subcommand parser (from Salmon_runner.py)
    run_parser = subparsers.add_parser('RunSalmonQuant', help='Run Salmon pipeline for RNA-seq quantification')
    
    # Global arguments for run
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
    process_parser = subparsers.add_parser('process', help='Process Salmon output into a combined table')
    process_parser.add_argument('-o', '--output', required=True, help='Directory with the output of the Salmon pipeline')
    process_parser.add_argument('--result_name', '--rn', default='table.txt', 
                            help='Name of the output file to store the combined results')
    process_parser.add_argument('--mode', '-m', default='cmd', 
                            help='Mode to run the commands. Options are "cmd" for command line and "python" for Python code.')
    process_parser.add_argument('--includes_alternative_genome', action='store_true', 
                            help='If True, the alternative genome was included in the analysis.')
    
    # Create the voom subcommand parser (from voom_from_salmon.py)
    voom_parser = subparsers.add_parser('voom', help='Preprocess Salmon output for voom analysis')
    # Add arguments for voom analysis as needed
    
    # Parse arguments
    args = parser.parse_args()
    
    # Execute the appropriate subcommand
    if args.command == 'RunSalmonQuant':
        run_salmon_main(args)
    elif args.command == 'process':
        sys.argv = [sys.argv[0]] + sys.argv[2:]  # Adjust argv to simulate direct call to the module
        process_main()
    elif args.command == 'voom':
        sys.argv = [sys.argv[0]] + sys.argv[2:]  # Adjust argv to simulate direct call to the module
        voom_main()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()