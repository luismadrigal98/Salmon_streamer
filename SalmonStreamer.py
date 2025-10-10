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
from subprograms.QTL_runner import main as qtl_runner_main
from subprograms.TranslateSalmon import main as translate_salmon_main
from subprograms.CalculateRawReads import main as calculate_raw_reads_main
from subprograms.CalculateCPM import main as calculate_cpm_main
from subprograms.ProcessGenotypes import main as process_genotypes_main
from subprograms.MakePhenotypes import main as make_phenotypes_main
from subprograms.PrepareQTLInputs import main as prepare_qtl_inputs_main
from subprograms.ParentalDE import main as parental_de_main
from src.postprocessing_utilities import process_post_pipeline

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
    misc_group.add_argument('--email', required=False, 
                        help='Email address for SLURM job notifications (optional)', default=None)
    misc_group.add_argument('--partition', required=False, 
                        help='SLURM partition to use', default='sixhour,eeb,kucg,kelly')
    misc_group.add_argument('--conda_env', required=False, 
                        help='Conda environment name to activate', default='salmon')
    misc_group.add_argument('--time_limit', required=False, 
                        help='Time limit for SLURM jobs', default='05:59:00')
    misc_group.add_argument('--module_load_cmd', required=False, 
                        help='Module load command for conda (use "none" to skip)', default='module load conda')
    
    # Create the process subcommand parser (from Salmon_output_processor.py)
    process_parser = subparsers.add_parser('ProcessSalmonOut', help='Process Salmon output into a combined table')
    process_parser.add_argument('-o', '--output', required=True, help='Directory with the output of the Salmon pipeline')
    process_parser.add_argument('--result_name', '--rn', default='table.txt', 
                            help='Name of the output file to store the combined results')
    process_parser.add_argument('--mode', '-m', default='python', 
                            help='Mode to run the commands. Options are "cmd" for command line and "python" (default) for Python code.')
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
    
    # --- Add QTL Runner Subcommand ---
    qtl_parser = subparsers.add_parser('RunQTL', help='Run QTL analysis on expression data')
    
    # Required file paths
    qtl_parser.add_argument('--phenofile-path', nargs="+", required=True, 
                            help="Path to phenotypes file or files (expression data)")
    qtl_parser.add_argument('--genfile-path', nargs="+", required=True, 
                            help="Path to genotypes file or files. Note that if multiple files are provided, is responsibility of the user to pass them in the same order as the previous entries")
    qtl_parser.add_argument('--outdir-base', required=True, 
                            help="Base output directory for QTL results")
    
    # Optional analysis parameters
    qtl_parser.add_argument('--covfile-path', nargs='+', required=False, 
                            help="Path to covariates file (optional) or files. Notice that if multiple files are provided, is responsibility of the user to pass them in the same order as the previous entries. Also, if a given input does to have covariates associates, you still need to pass `None`")
    qtl_parser.add_argument('--qtlmethod', default='mr',
                            help="QTL mapping method to use in R/qtl (default: mr)")
    qtl_parser.add_argument('--modeltype', default='normal',
                            help="Model type for scanone (default: normal)")
    qtl_parser.add_argument('--permnum', type=int, default=1000,
                            help="Number of permutations for significance testing (default: 1000)")
    qtl_parser.add_argument('--crosstype', default='f2',
                            help="Cross type for R/qtl (default: f2)")
    
    # SLURM job settings
    qtl_parser.add_argument('--rscript-executable', default='~/.conda/envs/PyR/bin/Rscript',
                            help="Path to the Rscript executable (default: ~/.conda/envs/PyR/bin/Rscript)")
    qtl_parser.add_argument('--partition', default='sixhour,eeb,kucg,kelly',
                            help="SLURM partition to use (default: sixhour,eeb,kucg,kelly)")
    qtl_parser.add_argument('--nodes', type=int, default=1,
                            help="Number of nodes per job (default: 1)")
    qtl_parser.add_argument('--ntasks-per-node', type=int, default=1,
                            help="Number of tasks per node (default: 1)")
    qtl_parser.add_argument('--cpus-per-task', type=int, default=1,
                            help="Number of CPUs per task (default: 1)")
    qtl_parser.add_argument('--mem-per-cpu', default='4G',
                            help="Memory per CPU (default: 4G)")
    qtl_parser.add_argument('--time-limit', default='05:59:00',
                            help="Time limit for jobs (default: 05:59:00)")
    
    # Job control parameters
    qtl_parser.add_argument('--max-jobs', type=int, default=100,
                            help="Maximum number of jobs to submit at once (default: 100)")
    qtl_parser.add_argument('--job-submission-delay', type=float, default=0.5,
                            help="Delay between job submissions in seconds (default: 0.5)")
    
    # Advanced SLURM job management parameters
    qtl_parser.add_argument('--max-concurrent-jobs', type=int, default=4900,
                            help="Maximum active SLURM jobs allowed (default: 4900)")
    qtl_parser.add_argument('--submission-batch-size', type=int, default=50,
                            help="Number of jobs to submit in each batch (default: 50)")
    qtl_parser.add_argument('--submission-wait-time', type=int, default=10,
                            help="Seconds to wait between job submission batches (default: 10)")
    qtl_parser.add_argument('--poll-interval', type=int, default=60,
                            help="Seconds between SLURM status checks (default: 60)")
    qtl_parser.add_argument('--wait-for-completion', action='store_true',
                            help="Wait for all jobs to complete before exiting (default: False)")
    qtl_parser.add_argument('--max-wait-time', type=int, default=86400,
                            help="Maximum seconds to wait for job completion (default: 86400 = 24h)")
    qtl_parser.add_argument('--job-completion-stringency', type=float, default=0.75,
                            help="Required fraction of successful jobs (default: 0.75)")
    
    # --- Add TranslateSalmon Subcommand ---
    translate_parser = subparsers.add_parser('TranslateSalmon', help='Translate Salmon outputs and organize read counts by allele')
    translate_parser.add_argument('cross', help='Cross identifier (e.g., SWB, SF)')
    translate_parser.add_argument('--genes-file', required=True, 
                                 help='Path to genes mapping file (e.g., Genes_to_updated_767_assembly.txt)')
    translate_parser.add_argument('--quant-results-file', required=True,
                                 help='Path to QUANT results file (e.g., QUANT_RESULTS_767_vs_SWB)')
    translate_parser.add_argument('--output-file', 
                                 help='Output file path (default: Salmon_outputs.IMlines.updated767.[cross].txt)')
    
    # --- Add CalculateRawReads Subcommand ---
    raw_reads_parser = subparsers.add_parser('CalculateRawReads', help='Calculate raw reads per plant from Salmon outputs')
    raw_reads_parser.add_argument('--salmon-files', nargs='+', required=True,
                                 help='List of Salmon output files (e.g., Salmon_outputs.IMlines.updated767.SWB.txt)')
    raw_reads_parser.add_argument('--output-file', default='raw.reads.per.plant.txt',
                                 help='Output file path (default: raw.reads.per.plant.txt)')
    
    # --- Add CalculateCPM Subcommand ---
    cpm_parser = subparsers.add_parser('CalculateCPM', help='Calculate CPM for PCA analysis')
    cpm_parser.add_argument('--raw-reads-file', required=True,
                           help='Path to raw reads per plant file (e.g., raw.reads.per.plant.txt)')
    cpm_parser.add_argument('--salmon-files', nargs='+', required=True,
                           help='List of Salmon output files (e.g., Salmon_outputs.IMlines.updated767.SWB.txt)')
    cpm_parser.add_argument('--cpm-min', type=float, default=5.0,
                           help='Minimum CPM threshold (default: 5.0)')
    cpm_parser.add_argument('--output-file', default='RawSamples_forPCA',
                           help='Output file path (default: RawSamples_forPCA)')
    
    # --- Add ProcessPost Subcommand (Combined pipeline) ---
    process_post_parser = subparsers.add_parser('ProcessPost', help='Run complete post-processing pipeline (TranslateSalmon + CalculateRawReads + CalculateCPM)')
    process_post_parser.add_argument('--crosses', nargs='+', required=True,
                                    help='List of cross identifiers (e.g., SWB SF)')
    process_post_parser.add_argument('--genes-file', required=True,
                                    help='Path to genes mapping file (e.g., Genes_to_updated_767_assembly.txt)')
    process_post_parser.add_argument('--quant-results-files', nargs='+', required=True,
                                    help='List of QUANT results files (in same order as crosses)')
    process_post_parser.add_argument('--cpm-min', type=float, default=5.0,
                                    help='Minimum CPM threshold (default: 5.0)')
    process_post_parser.add_argument('--output-dir', default='.',
                                    help='Output directory (default: current directory)')
    
    # --- Add ProcessGenotypes Subcommand ---
    genotypes_parser = subparsers.add_parser('ProcessGenotypes', help='Process genotypes from transcript mapping data')
    genotypes_parser.add_argument('--cross', required=True,
                                 help='Cross identifier (e.g., SWB, SF, 1034)')
    genotypes_parser.add_argument('--allele-counts-file', required=True,
                                 help='Path to allele counts file (e.g., allele_counts.combined.updated767.CROSS.txt)')
    genotypes_parser.add_argument('--samples-file', required=True,
                                 help='Path to samples file (e.g., CROSS.combined.samples.txt)')
    genotypes_parser.add_argument('--genes-file', required=True,
                                 help='Path to genes mapping file (e.g., Genes_to_updated_767_assembly.txt)')
    genotypes_parser.add_argument('--min-parental-lines', type=int, default=5,
                                 help='Minimum parental lines called (default: 5)')
    genotypes_parser.add_argument('--mapping-threshold', type=float, default=0.95,
                                 help='Minimum frequency of mapping to correct allele (default: 0.95)')
    genotypes_parser.add_argument('--min-reads-per-plant', type=int, default=100000,
                                 help='Minimum reads per plant (default: 100000)')
    genotypes_parser.add_argument('--min-reads-per-call', type=int, default=6,
                                 help='Minimum reads to count an F2 call (default: 6)')
    genotypes_parser.add_argument('--f2-fraction-threshold', type=float, default=0.5,
                                 help='Fraction of F2s required (default: 0.5)')
    genotypes_parser.add_argument('--homozygous-threshold', type=float, default=0.95,
                                 help='Homozygous calling threshold (default: 0.95)')
    genotypes_parser.add_argument('--het-maf', type=float, default=0.25,
                                 help='Heterozygous calling MAF threshold (default: 0.25)')
    genotypes_parser.add_argument('--output-dir', default='.',
                                 help='Output directory (default: current directory)')
    
    # --- Add MakePhenotypes Subcommand ---
    phenotypes_parser = subparsers.add_parser('MakePhenotypes', help='Generate phenotype files from expression data')
    phenotypes_parser.add_argument('--genes-by-cross-file', required=True,
                                  help='Path to genes by cross file (e.g., Genes_by_cross.txt)')
    phenotypes_parser.add_argument('--total-reads-files', nargs='+', required=True,
                                  help='List of total reads files for each cross (e.g., Total_reads_byplant.SF)')
    phenotypes_parser.add_argument('--readcounts-files', nargs='+', required=True,
                                  help='List of readcounts files for each cross')
    phenotypes_parser.add_argument('--f2-lists', nargs='+', required=True,
                                  help='List of F2 inclusion files (e.g., SF.included.f2s.txt)')
    phenotypes_parser.add_argument('--crosses', nargs='+', required=True,
                                  help='List of cross identifiers')
    phenotypes_parser.add_argument('--im-lines', nargs='+', 
                                  default=["62","155","444","502","541","664","909","1034","1192"],
                                  help='List of inbred line identifiers')
    phenotypes_parser.add_argument('--output-dir', default='pfiles',
                                  help='Output directory for phenotype files (default: pfiles)')
    phenotypes_parser.add_argument('--summary-file', default='f2summaries_by_gene.txt',
                                  help='Summary output file (default: f2summaries_by_gene.txt)')
    
    # --- Add PrepareQTLInputs Subcommand ---
    qtl_inputs_parser = subparsers.add_parser('PrepareQTLInputs', help='Prepare inputs for QTL analysis')
    qtl_inputs_parser.add_argument('--cross', required=True,
                                  help='Cross identifier (e.g., SWB, SF, 1034)')
    qtl_inputs_parser.add_argument('--genotype-pp-file', required=True,
                                  help='Path to genotype posterior probabilities file (e.g., CROSS.F2_geno_PP.txt)')
    qtl_inputs_parser.add_argument('--estimates-files', nargs='+', required=True,
                                  help='List of estimates files for each chromosome')
    qtl_inputs_parser.add_argument('--genes-by-cross-file', required=True,
                                  help='Path to genes by cross file')
    qtl_inputs_parser.add_argument('--phenotype-group', required=True,
                                  help='Phenotype group identifier (e.g., gene.group1)')
    qtl_inputs_parser.add_argument('--phenotype-files-dir', default='pfiles',
                                  help='Directory containing phenotype files (default: pfiles)')
    qtl_inputs_parser.add_argument('--output-dir', default='rQTL_files',
                                  help='Output directory for R/qtl files (default: rQTL_files)')
    qtl_inputs_parser.add_argument('--chromosomes', nargs='+', type=int, default=list(range(1, 15)),
                                  help='List of chromosome numbers (default: 1-14)')
    
    # --- Add ParentalDE Subcommand ---
    parental_de_parser = subparsers.add_parser('ParentalDE', 
                                               help='Test for differential expression between parental lines')
    parental_de_parser.add_argument('--expression-file', required=True,
                                   help='Path to expression data file (tab-delimited, genes in rows)')
    parental_de_parser.add_argument('--parent1-samples', nargs='+', required=True,
                                   help='Sample IDs or patterns for parent 1 (space-separated)')
    parental_de_parser.add_argument('--parent2-samples', nargs='+', required=True,
                                   help='Sample IDs or patterns for parent 2 (space-separated)')
    parental_de_parser.add_argument('--parent1-name', required=True,
                                   help='Name for parent 1 (e.g., IM767)')
    parental_de_parser.add_argument('--parent2-name', required=True,
                                   help='Name for parent 2 (e.g., SF, SWB)')
    parental_de_parser.add_argument('--test-method', choices=['ttest', 'mannwhitney'], default='ttest',
                                   help='Statistical test method (default: ttest)')
    parental_de_parser.add_argument('--fdr-threshold', type=float, default=0.05,
                                   help='FDR threshold for significance (default: 0.05)')
    parental_de_parser.add_argument('--output-dir', default='./DE_results',
                                   help='Output directory for results (default: ./DE_results)')
    
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
        voom_main(args)
    elif args.command == 'RunQTL':
        qtl_runner_main(args)  # Call the QTL runner main function
    elif args.command == 'TranslateSalmon':
        translate_salmon_main(args)
    elif args.command == 'CalculateRawReads':
        calculate_raw_reads_main(args)
    elif args.command == 'CalculateCPM':
        calculate_cpm_main(args)
    elif args.command == 'ProcessPost':
        # Run the complete post-processing pipeline
        process_post_pipeline(args)
    elif args.command == 'ProcessGenotypes':
        process_genotypes_main(args)
    elif args.command == 'MakePhenotypes':
        make_phenotypes_main(args)
    elif args.command == 'PrepareQTLInputs':
        prepare_qtl_inputs_main(args)
    elif args.command == 'ParentalDE':
        parental_de_main(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()