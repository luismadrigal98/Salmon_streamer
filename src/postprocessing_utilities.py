"""
Set of functions to postprocess the output files after using Salmon

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2024-12-15

"""

import subprocess
import os
import pandas as pd
import numpy as np

def combine_results(output_dir, result_name = 'table.txt', mode = 'cmd', includes_alternative_genome = False):
    """
    Combine all Salmon quantification results into one table.

    Parameters
    ----------
    output_dir : str
        Path to the output directory containing the quantification results
    result_name : str
        Name of the output file to store the combined results
    mode : str
        Mode to run the commands. Options are 'cmd' for command line and 'python' for Python code.
    
    Returns
    -------
    None
    """
    # Change to the output directory
    os.chdir(output_dir)

    if mode == 'cmd': # Currently broken!!!
        
        # Output column 5 (counts, not TPM column 4) as text file with ID as name
        cmd_cut = "for i in $(ls | grep quant); do cut -f5 $i/quant.sf > $i.txt; done"
        subprocess.run(cmd_cut, shell=True, executable='/bin/bash')

        # Replace first row in file with textfile name ie ID
        cmd_replace = r"for i in $(ls | grep quant.txt); do cat <(echo $i | perl -pe 's/.txt//') <(tail -n +2 $i) > ${i}_count.txt; done"
        subprocess.run(cmd_replace, shell=True, executable='/bin/bash')

        # Get gene names from one of the results, to be used as first column in compiled results file
        # Dynamically get gene names from the first quantification result
        first_quant_dir = subprocess.check_output("ls | grep quant | head -n 1", shell=True, executable='/bin/bash').decode().strip()
        cmd_cut_gene_names = f"cut -f1 {first_quant_dir}/quant.sf > 1count.txt"
        subprocess.run(cmd_cut_gene_names, shell=True, executable='/bin/bash')

        # Put everything together
        cmd_paste = f"paste *count* > {result_name}"
        subprocess.run(cmd_paste, shell=True, executable='/bin/bash')

        if includes_alternative_genome:
            # Manually sort the table so that 2 species are on alternative rows, and remove "quant" from header row
            cmd_sort = r"perl -pe 's/.fastq//g ; s/_quant//g ; s/=/\t/ ; s/aName/aName\taName/' " + result_name + r" | sort -t$'\t' -k2,2 -k1,1 | perl -pe 's/ID\t/ID=/ ; s/\taName//' > clean_" + result_name
            subprocess.run(cmd_sort, shell=True, executable='/bin/bash')

    elif mode == 'python':
        # List to store DataFrames
        df_list = []
        sample_ids = []

        # Corrected directory listing
        quant_dirs = [
            d for d in os.listdir('.')
            if d.endswith('_quant') and os.path.isdir(d)
        ]

        for quant_dir in quant_dirs:
            quant_path = os.path.join(quant_dir, 'quant.sf')
            if os.path.exists(quant_path):
                # Read the quant.sf file
                df = pd.read_csv(
                    quant_path,
                    sep='\t',
                    usecols=['Name', 'NumReads']
                )
                df.rename(
                    columns={'NumReads': quant_dir.replace('_quant', '')},
                    inplace=True
                )
                df_list.append(df.set_index('Name'))
                sample_ids.append(quant_dir.replace('_quant', ''))

        # Combine all DataFrames on the 'Name' column
        combined_df = pd.concat(df_list, axis=1)
        combined_df.reset_index(inplace=True)

        # Write the combined table to a file
        combined_df.to_csv(result_name, sep='\t', index=False)

def translate_salmon_outputs(cross, genes_file, quant_results_file, output_file):
    """
    Take salmon outputs (QUANT) and organize read counts to each allele of each gene per sample.
    
    Parameters
    ----------
    cross : str
        Cross identifier (e.g., SWB, SF)
    genes_file : str
        Path to genes mapping file (e.g., Genes_to_updated_767_assembly.txt)
    quant_results_file : str
        Path to QUANT results file
    output_file : str
        Output file path
        
    Returns
    -------
    None
    """
    
    print(f"Processing cross: {cross}")
    
    # Load included genes
    includedgenes = {}
    with open(genes_file, "r") as inx:
        for line_id, liner in enumerate(inx):
            colx = liner.replace('\n', '').split('\t')
            # Chrom	stpos	endpos	old_name	new_name	62	155	444	502	541	664	909	1034	1192	scored_pops	new_chrom	new_start	new_end
            includedgenes[colx[4] + ".v2.1"] = 1
    
    # Process QUANT results
    data = {}
    plants_in_sequence = []
    fn = [0, 0]  # [found, not_found]
    
    with open(quant_results_file, "r") as src:
        for line_idx, line in enumerate(src):
            cols = line.replace('\n', '').split('\t')
            if line_idx == 0:
                # Header row - extract sample names
                for j in range(1, len(cols)):
                    plants_in_sequence.append(cols[j])
            else:
                # Data rows
                gname = cols[0].split("=")[1]
                try:
                    nn = includedgenes[gname]
                    fn[0] += 1
                    try:
                        ux = data[gname]
                    except KeyError:
                        data[gname] = {"IM767": [], cross: []}
                    
                    allele = cols[0].split("_")[0]
                    
                    for j in range(1, len(cols)):
                        data[gname][allele].append(cols[j])
                        
                except KeyError:
                    fn[1] += 1
    
    # Write output
    with open(output_file, "w") as out1:
        out1.write("GeneID")
        for j in range(len(plants_in_sequence)):
            out1.write(f"\t767_{plants_in_sequence[j]}")
            out1.write(f"\t{cross}_{plants_in_sequence[j]}")
        out1.write("\n")
        
        for gname in data:
            out1.write(gname)
            for j in range(len(plants_in_sequence)):
                out1.write(f"\t{data[gname]['IM767'][j]}")
                out1.write(f"\t{data[gname][cross][j]}")
            out1.write("\n")
    
    print(f"Processed {cross}: Found {fn[0]}, Not found {fn[1]}")

def calculate_raw_reads_per_plant(salmon_files, output_file):
    """
    Calculate total reads per plant from multiple Salmon output files.
    
    Parameters
    ----------
    salmon_files : list
        List of Salmon output files
    output_file : str
        Output file path
        
    Returns
    -------
    None
    """
    
    data = {}
    
    for salmon_file in salmon_files:
        # Determine cross type from filename
        if "SWB" in salmon_file:
            cross_prefix = "SWBcross_"
        elif "SF" in salmon_file:
            cross_prefix = "SFcross_"
        else:
            # Extract cross from filename pattern
            cross_name = salmon_file.split(".")[-2]  # Assumes pattern like ...updated767.CROSS.txt
            cross_prefix = f"{cross_name}cross_"
        
        relevant_numbers = {}
        
        with open(salmon_file, "r") as src:
            for line_idx, line in enumerate(src):
                cols = line.replace('\n', '').split('\t')
                if line_idx == 0:
                    # Header row - extract sample positions
                    for j in range(1, len(cols), 2):
                        plt = cross_prefix + cols[j][4:]  # Remove "767_" prefix
                        relevant_numbers[plt] = [j, j+1]
                        data[plt] = 0.0  # Initialize sum of all reads
                else:
                    # Data rows - sum reads for each plant
                    for plt in relevant_numbers:
                        data[plt] += (float(cols[relevant_numbers[plt][0]]) + 
                                    float(cols[relevant_numbers[plt][1]]))
    
    # Write output
    with open(output_file, "w") as out1:
        for plt in data:
            out1.write(f"{plt}\t{data[plt]}\n")

def calculate_cpm_for_pca(raw_reads_file, salmon_files, cpm_min, output_file):
    """
    Standardize raw reads to CPM (Counts Per Million) for PCA analysis.
    
    Parameters
    ----------
    raw_reads_file : str
        Path to raw reads per plant file
    salmon_files : list
        List of Salmon output files
    cpm_min : float
        Minimum CPM threshold
    output_file : str
        Output file path
        
    Returns
    -------
    None
    """
    
    # Load read counts per plant
    reads = {}
    with open(raw_reads_file, "r") as in1:
        for line_idxx, liner in enumerate(in1):
            colx = liner.replace('\n', '').split('\t')
            sd = colx[0]
            reads[sd] = float(colx[1])
    
    data = {}
    nx = 0
    lastg = None
    
    # Process each salmon file
    for salmon_file in salmon_files:
        # Determine cross type from filename
        if "SWB" in salmon_file:
            cross_prefix = "SWBcross_"
        elif "SF" in salmon_file:
            cross_prefix = "SFcross_"
        else:
            # Extract cross from filename pattern
            cross_name = salmon_file.split(".")[-2]
            cross_prefix = f"{cross_name}cross_"
        
        relevant_numbers = {}
        
        with open(salmon_file, "r") as src:
            for line_idx, line in enumerate(src):
                cols = line.replace('\n', '').split('\t')
                if line_idx == 0:
                    # Header row
                    for j in range(1, len(cols), 2):
                        nx += 1
                        plt = cross_prefix + cols[j][4:]
                        relevant_numbers[plt] = [j, j+1]
                else:
                    # Data rows - calculate CPM
                    geneid = cols[0]
                    if geneid not in data:
                        data[geneid] = {}
                    
                    lastg = geneid
                    for plt in relevant_numbers:
                        raw_count = (float(cols[relevant_numbers[plt][0]]) + 
                                   float(cols[relevant_numbers[plt][1]]))
                        data[geneid][plt] = raw_count * 1000000.0 / reads[plt]
    
    # Write output
    with open(output_file, "w") as out1:
        # Write header
        out1.write("geneid")
        if lastg:
            for plt in data[lastg]:
                out1.write(f'\t{plt}')
        out1.write('\n')
        
        # Write data for genes that meet criteria
        for geneid in data:
            if len(data[geneid]) == nx:  # Gene present in both crosses
                sreads = []
                for plt in data[geneid]:
                    sreads.append(data[geneid][plt])
                
                if np.average(sreads) >= cpm_min:  # Gene meets minimum CPM threshold
                    out1.write(geneid)
                    for plt in data[geneid]:
                        out1.write(f'\t{data[geneid][plt]}')
                    out1.write('\n')

def process_post_pipeline(args):
    """
    Run the complete post-processing pipeline: TranslateSalmon -> CalculateRawReads -> CalculateCPM
    """
    import os
    import sys
    from argparse import Namespace
    
    # Add parent directory to path for imports
    parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)
    
    # Import the functions locally to avoid circular imports
    from subprograms.TranslateSalmon import main as translate_salmon_main
    from subprograms.CalculateRawReads import main as calculate_raw_reads_main
    from subprograms.CalculateCPM import main as calculate_cpm_main
    
    print("Running complete post-processing pipeline...")
    
    # Change to output directory
    original_dir = os.getcwd()
    if args.output_dir != '.':
        os.makedirs(args.output_dir, exist_ok=True)
        os.chdir(args.output_dir)
    
    try:
        # Step 1: Run TranslateSalmon for each cross
        salmon_output_files = []
        for i, cross in enumerate(args.crosses):
            translate_args = Namespace(
                cross=cross,
                genes_file=args.genes_file,
                quant_results_file=args.quant_results_files[i],
                output_file=f"Salmon_outputs.IMlines.updated767.{cross}.txt"
            )
            print(f"Step 1.{i+1}: Translating Salmon outputs for cross {cross}")
            translate_salmon_main(translate_args)
            salmon_output_files.append(translate_args.output_file)
        
        # Step 2: Calculate raw reads per plant
        raw_reads_args = Namespace(
            salmon_files=salmon_output_files,
            output_file="raw.reads.per.plant.txt"
        )
        print("Step 2: Calculating raw reads per plant")
        calculate_raw_reads_main(raw_reads_args)
        
        # Step 3: Calculate CPM for PCA
        cpm_args = Namespace(
            raw_reads_file="raw.reads.per.plant.txt",
            salmon_files=salmon_output_files,
            cpm_min=args.cpm_min,
            output_file="RawSamples_forPCA"
        )
        print("Step 3: Calculating CPM for PCA analysis")
        calculate_cpm_main(cpm_args)
        
        print("Post-processing pipeline completed successfully!")
        print("Output files:")
        for file in salmon_output_files:
            print(f"  - {file}")
        print("  - raw.reads.per.plant.txt")
        print("  - RawSamples_forPCA")
        
    finally:
        # Return to original directory
        os.chdir(original_dir)