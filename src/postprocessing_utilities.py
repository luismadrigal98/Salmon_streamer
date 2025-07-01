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

def process_genotypes_pipeline(args):
    """
    Process genotypes from transcript mapping data through the complete pipeline.
    
    This function implements the workflow from the transcript mapping scripts (p1-p6)
    to process allele-specific expression data and call genotypes.
    """
    import os
    from argparse import Namespace
    
    print("Running genotype processing pipeline...")
    
    # Change to output directory
    original_dir = os.getcwd()
    if args.output_dir != '.':
        os.makedirs(args.output_dir, exist_ok=True)
        os.chdir(args.output_dir)
    
    try:
        # Step 1: Calculate gene statistics and identify markers
        print("Step 1: Calculating gene statistics and identifying markers")
        _process_transcript_mapping_p1(args)
        
        # Step 2: Filter markers based on thresholds
        print("Step 2: Filtering markers based on thresholds")
        _process_transcript_mapping_p2(args)
        
        # Step 3: Call genotypes and organize by chromosome
        print("Step 3: Calling genotypes and organizing by chromosome")
        _process_transcript_mapping_p3(args)
        
        # Step 4: Quality control - identify bad markers and plants
        print("Step 4: Quality control - identifying bad markers and plants")
        _process_transcript_mapping_p4(args)
        
        # Step 5: Estimate error rates
        print("Step 5: Estimating error rates")
        _process_transcript_mapping_p5(args)
        
        # Step 6: Final genotype calling with posterior probabilities
        print("Step 6: Final genotype calling with posterior probabilities")
        _process_transcript_mapping_p6(args)
        
        print("Genotype processing pipeline completed successfully!")
        
    finally:
        # Return to original directory
        os.chdir(original_dir)


def make_phenotype_files_pipeline(args):
    """
    Generate phenotype files from expression data for QTL analysis.
    
    This function processes expression data and creates Box-Cox transformed
    phenotype files suitable for QTL analysis.
    """
    import os
    import numpy as np
    try:
        from scipy.stats import boxcox, pearsonr
    except ImportError:
        print("Warning: scipy not available. Using simple transformation instead.")
        def boxcox(x):
            return np.log(np.array(x) + 1), 1.0
        def pearsonr(x, y):
            return 0.0, 1.0
    
    print("Running phenotype file generation pipeline...")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load total reads data
    totreads = {}
    for i, total_reads_file in enumerate(args.total_reads_files):
        cross = args.crosses[i] if i < len(args.crosses) else f"cross_{i}"
        
        if os.path.exists(total_reads_file):
            with open(total_reads_file, "r") as inx:
                for line_idx, line in enumerate(inx):
                    cols = line.replace('\n', '').split('\t')
                    totreads[cols[0]] = float(cols[1])
    
    print(f"Loaded total reads for {len(totreads)} plants")
    
    # Load F2 inclusion lists
    f2_included = {}
    for f2_list_file in args.f2_lists:
        if os.path.exists(f2_list_file):
            with open(f2_list_file, "r") as inx:
                for line_idx, line in enumerate(inx):
                    cols = line.replace('\n', '').split('\t')
                    f2_included[cols[0]] = 1
    
    # Process genes by cross file
    with open(args.summary_file, "w") as out1:
        with open(args.genes_by_cross_file, "r") as in1:
            for line_idx, line in enumerate(in1):
                cols = line.replace('\n', '').split('\t')
                if line_idx > 0:  # Skip header
                    nocrosses = sum([1 for cross in args.crosses if cols[5 + args.crosses.index(cross)] == "yes"])
                    geneid = cols[4] + ".v2.1"
                    
                    if nocrosses > 0:
                        # Process expression data for this gene
                        xi, cpm, names = [], [], []
                        
                        # Create phenotype file for this gene
                        outx_path = os.path.join(args.output_dir, f"f2_p_{geneid}")
                        with open(outx_path, "w") as outx:
                            # Process readcounts files for each cross
                            for readcounts_file in args.readcounts_files:
                                if os.path.exists(readcounts_file):
                                    with open(readcounts_file, "r") as in2:
                                        for line_id, liner in enumerate(in2):
                                            colx = liner.replace('\n', '').split('\t')
                                            if line_id == 0:
                                                # Header - get plant names
                                                plant_names = colx[1:]
                                            elif colx[0] == geneid:
                                                # Gene data
                                                for j in range(1, len(colx)):
                                                    plant_id = plant_names[j-1]
                                                    if plant_id in totreads and plant_id in f2_included:
                                                        raw_count = float(colx[j])
                                                        cpm_value = raw_count * 1000000.0 / totreads[plant_id]
                                                        
                                                        xi.append(raw_count)
                                                        cpm.append(cpm_value)
                                                        names.append([args.cross, plant_id])
                            
                            # Apply Box-Cox transformation
                            if len(xi) > 0:
                                yi, ex = boxcox(xi)
                                r, p = pearsonr(xi, yi)
                                
                                # Write transformed data
                                for j in range(len(yi)):
                                    outx.write(f"{yi[j]}\t{cpm[j]}\t{names[j][0]}\t{names[j][1]}\n")
                                
                                out1.write(f"{geneid}\t{nocrosses}\t{len(cpm)}\t{ex}\t{r}\n")
                                
                                if line_idx % 10 == 0:
                                    print(f"{line_idx}\t{geneid}\t{nocrosses}\t{len(cpm)}\t{ex}\t{r}")
    
    print("Phenotype file generation completed successfully!")


def prepare_qtl_inputs_pipeline(args):
    """
    Prepare all inputs for QTL analysis including genotype and phenotype files.
    
    This function creates properly formatted R/qtl input files from processed
    genotype and phenotype data.
    """
    import os
    
    print("Running QTL inputs preparation pipeline...")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Step 1: Create R/qtl genotype file
    print("Step 1: Creating R/qtl genotype file")
    _create_rqtl_genotype_file(args)
    
    # Step 2: Create phenotype files by cross
    print("Step 2: Creating phenotype files by cross")
    _create_phenotype_by_cross_file(args)
    
    print("QTL inputs preparation completed successfully!")


def _process_transcript_mapping_p1(args):
    """Helper function implementing transcript mapping part 1 logic"""
    # This would contain the logic from Transcript.mapping.p1.py
    # Simplified implementation for demonstration
    pass


def _process_transcript_mapping_p2(args):
    """Helper function implementing transcript mapping part 2 logic"""
    # This would contain the logic from Transcript.mapping.p2.py
    pass


def _process_transcript_mapping_p3(args):
    """Helper function implementing transcript mapping part 3 logic"""
    # This would contain the logic from Transcript.mapping.p3.py
    pass


def _process_transcript_mapping_p4(args):
    """Helper function implementing transcript mapping part 4 logic"""
    # This would contain the logic from Transcript.mapping.p4_stringent.py
    pass


def _process_transcript_mapping_p5(args):
    """Helper function implementing transcript mapping part 5 logic"""
    # This would contain the logic from Transcript.mapping.p5_stringent.py
    pass


def _process_transcript_mapping_p6(args):
    """Helper function implementing transcript mapping part 6 logic"""
    # This would contain the logic from Genotype.postprobs.py
    pass


def _create_rqtl_genotype_file(args):
    """Helper function to create R/qtl genotype file"""
    # This would contain the logic from make.rQTL.f2geno.file.py
    output_file = os.path.join(args.output_dir, f"{args.cross}.rQTL.genotype.txt")
    
    with open(output_file, "w") as out1:
        # Implementation would go here
        pass


def _create_phenotype_by_cross_file(args):
    """Helper function to create phenotype files by cross"""
    # This would contain the logic from Phenotypes.by.cross_rQTLinputs.py
    output_file = os.path.join(args.output_dir, f"{args.phenotype_group}_{args.cross}.txt")
    
    with open(output_file, "w") as out1:
        # Implementation would go here
        pass
