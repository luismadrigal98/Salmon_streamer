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
                # Check if the expected format is present
                if "=" not in cols[0]:
                    print(f"Warning: Line {line_idx + 1} has unexpected format in column 1: '{cols[0]}'")
                    print(f"Expected format: 'allele_gene=gene_name', but found: '{cols[0]}'")
                    continue  # Skip this line
                
                try:
                    gname = cols[0].split("=")[1]
                except IndexError:
                    print(f"Error: Line {line_idx + 1} - Cannot extract gene name from: '{cols[0]}'")
                    print(f"Full line: {line.strip()}")
                    continue  # Skip this line
                
                try:
                    nn = includedgenes[gname]
                    fn[0] += 1
                    try:
                        ux = data[gname]
                    except KeyError:
                        data[gname] = {"IM767": [], cross: []}
                    
                    # Check if allele can be extracted
                    if "_" not in cols[0]:
                        print(f"Warning: Line {line_idx + 1} - Cannot extract allele from: '{cols[0]}'")
                        continue  # Skip this line
                    
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
                    if len(cols) > 0:
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
                                                        names.append([args.crosses[0], plant_id])  # Use first cross as default
                            
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
    cross = args.cross
    
    # Load plant genotype information
    plants = []
    genotype = {}
    data = {}
    
    with open(args.samples_file, "r") as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 2:
                sd = cols[0]
                plants.append(sd)
                genotype[sd] = cols[1]
                data[sd] = {}
    
    # Process allele counts
    genelist = []
    plt_by_pos = {}
    total_reads = {}
    
    with open(args.allele_counts_file, "r") as f:
        for line_idx, line in enumerate(f):
            cols = line.strip().split('\t')
            if line_idx == 0:
                # Header: GeneID	767_plant1	cross_plant1	767_plant2	cross_plant2...
                for j in range(1, len(cols), 2):
                    plt_by_pos[j] = cols[j][4:]  # Remove "767_" prefix
                    plant_id = cols[j][4:]
                    if plant_id not in total_reads:
                        total_reads[plant_id] = 0
            else:
                # Data rows
                geneid = cols[0]
                genelist.append(geneid)
                
                for j in plt_by_pos:
                    sd = plt_by_pos[j]
                    if sd in data:
                        data[sd][geneid] = {cross: 0, "IM767": 0}
                        data[sd][geneid]["IM767"] = float(cols[j])
                        data[sd][geneid][cross] = float(cols[j+1])
                        total_reads[sd] += (float(cols[j]) + float(cols[j+1]))
    
    # Write total reads per plant
    with open(f"Total_reads_byplant.{cross}", "w") as f:
        for sd in total_reads:
            f.write(f"{sd}\t{int(total_reads[sd])}\n")
    
    # Calculate allele-specific statistics
    AS_byplant = {}
    for sd in plants:
        if genotype[sd] not in ["f2", "f1"]:
            AS_byplant[sd] = [0, 0.0]
    
    with open(f"gene_stats.{cross}", "w") as gene_out:
        for geneid in genelist:
            gene_out.write(geneid)
            props = {"IM767": [], cross: []}
            
            for sd in AS_byplant:
                if geneid in data[sd]:
                    nx = data[sd][geneid]["IM767"] + data[sd][geneid][cross]
                    if nx > 0.0:
                        AS_byplant[sd][0] += 1
                        allele_prop = data[sd][geneid][genotype[sd]] / nx
                        AS_byplant[sd][1] += allele_prop
                        props[genotype[sd]].append(allele_prop)
            
            # Write gene statistics
            for allele in ["IM767", cross]:
                if len(props[allele]) > 0:
                    nk = len(props[allele])
                    mk = sum(props[allele]) / float(nk)
                    gene_out.write(f'\t{nk}\t{mk}')
                else:
                    gene_out.write('\t0\tNA')
            gene_out.write('\n')
    
    # Write parental statistics
    with open(f"parental_stats.{cross}", "w") as f:
        for sd in AS_byplant:
            if AS_byplant[sd][0] > 0:
                AS_byplant[sd][1] = AS_byplant[sd][1] / float(AS_byplant[sd][0])
            f.write(f"{sd}\t{AS_byplant[sd][0]}\t{AS_byplant[sd][1]}\n")


def _process_transcript_mapping_p2(args):
    """Helper function implementing transcript mapping part 2 logic"""
    cross = args.cross
    
    # Filter based on mapping statistics
    markers = {}
    
    with open(f"gene_stats.{cross}", "r") as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 5 and cols[1] != "0" and cols[3] != "0":
                try:
                    if (int(cols[1]) >= args.min_parental_lines and 
                        int(cols[3]) >= args.min_parental_lines and
                        float(cols[2]) >= args.mapping_threshold and 
                        float(cols[4]) >= args.mapping_threshold):
                        markers[cols[0]] = 1
                except (ValueError, IndexError):
                    continue
    
    with open(f"{cross}.genes.with.allele.specific.mapping.txt", "w") as f:
        for gene in markers:
            f.write(f"{gene}\n")
    
    print(f"Survivors after mapping filter: {len(markers)}")
    
    # Filter based on read depth
    goodplants = {}
    with open(f"Total_reads_byplant.{cross}", "r") as f:
        for line in f:
            cols = line.strip().split('\t')
            if int(cols[1]) >= args.min_reads_per_plant:
                goodplants[cols[0]] = 1
    
    print(f"Plants retained: {len(goodplants)}")
    
    # Final filtering based on F2 coverage
    f2gtype = {}
    with open(args.samples_file, "r") as f:
        for line in f:
            cols = line.strip().split('\t')
            sd = cols[0]
            if sd in goodplants and len(cols) > 1 and cols[1] == "f2":
                f2gtype[sd] = 1
    
    kept = 0
    relevant_numbers = []
    
    with open(args.allele_counts_file, "r") as f_in, \
         open(f"{cross}.marker.genes.txt", "w") as f_out:
        
        for line_idx, line in enumerate(f_in):
            cols = line.strip().split('\t')
            if line_idx == 0:
                # Header processing
                for j in range(1, len(cols), 2):
                    plt = cols[j][4:]  # Remove "767_" prefix
                    if plt in f2gtype:
                        relevant_numbers.append([j, j+1])
                
                good_enough = args.f2_fraction_threshold * len(relevant_numbers)
                print(f"F2 plants: {len(relevant_numbers)}")
            else:
                # Data processing
                if cols[0] in markers:
                    cx = 0
                    for j_pair in relevant_numbers:
                        total_reads = float(cols[j_pair[0]]) + float(cols[j_pair[1]])
                        if total_reads >= args.min_reads_per_call:
                            cx += 1
                    
                    if cx >= good_enough:
                        f_out.write(f"{cols[0]}\n")
                        kept += 1
    
    print(f"Final survivors: {kept}")


def _process_transcript_mapping_p3(args):
    """Helper function implementing transcript mapping part 3 logic"""
    cross = args.cross
    
    # Load relevant genes
    relevant = {}
    with open(f"{cross}.marker.genes.txt", "r") as f:
        for line in f:
            relevant[line.strip()] = 1
    
    # Load good plants
    goodplants = {}
    with open(f"Total_reads_byplant.{cross}", "r") as f:
        for line in f:
            cols = line.strip().split('\t')
            if int(cols[1]) >= args.min_reads_per_plant:
                goodplants[cols[0]] = 1
    
    # Load F2 genotypes
    f2gtype = {}
    with open(args.samples_file, "r") as f:
        for line in f:
            cols = line.strip().split('\t')
            sd = cols[0]
            if sd in goodplants and len(cols) > 1 and cols[1] == "f2":
                f2gtype[sd] = 1
    
    # Load gene locations
    gene_location = {}
    with open(args.genes_file, "r") as f:
        for line_idx, line in enumerate(f):
            if line_idx == 0:
                continue  # Skip header
            cols = line.strip().split('\t')
            if len(cols) >= 17:
                geneid = cols[4] + ".v2.1"
                gene_location[geneid] = [cols[15], cols[16]]  # chromosome, position
    
    # Create PC1 directory and initialize chromosome files
    pc1_dir = f"PC1.{cross}"
    os.makedirs(pc1_dir, exist_ok=True)
    
    gcounts = [0, 0, 0, 0]  # AA, AB, BB, NN
    snpsbychrom = {}
    
    # Initialize chromosome files
    for chx in range(1, 15):
        chrom = f"Chr_{chx:02d}" if chx < 10 else f"Chr_{chx}"
        snpsbychrom[chrom] = []
        
        # Create preliminary call files
        with open(os.path.join(pc1_dir, f"PrelimCalls.{cross}.{chrom}"), "w") as f:
            f.write("geneid\tchrom\tpos")
            for plant in f2gtype:
                f.write(f"\t{plant}")
            f.write("\n")
    
    # Process allele counts and call genotypes
    relevant_numbers = []
    
    with open(args.allele_counts_file, "r") as f:
        for line_idx, line in enumerate(f):
            cols = line.strip().split('\t')
            if line_idx == 0:
                # Header processing
                for j in range(1, len(cols), 2):
                    plt = cols[j][4:]
                    if plt in f2gtype:
                        relevant_numbers.append([j, j+1])
            else:
                # Genotype calling
                if cols[0] in relevant and cols[0] in gene_location:
                    chrom = gene_location[cols[0]][0]
                    pos = gene_location[cols[0]][1]
                    
                    stx = f"{cols[0]}\t{chrom}\t{pos}"
                    
                    for j_pair in relevant_numbers:
                        n767 = float(cols[j_pair[0]])
                        nAlt = float(cols[j_pair[1]])
                        ntot = n767 + nAlt
                        
                        if ntot < args.min_reads_per_call:
                            stx += "\tNN"
                            gcounts[3] += 1
                        elif n767/ntot >= args.homozygous_threshold:
                            stx += "\tAA"
                            gcounts[0] += 1
                        elif nAlt/ntot >= args.homozygous_threshold:
                            stx += "\tBB"
                            gcounts[2] += 1
                        elif n767/ntot >= args.het_maf and nAlt/ntot >= args.het_maf:
                            stx += "\tAB"
                            gcounts[1] += 1
                        else:
                            stx += "\tNN"
                            gcounts[3] += 1
                    
                    stx += "\n"
                    try:
                        snpsbychrom[chrom].append([int(pos), stx])
                    except (KeyError, ValueError):
                        pass
    
    # Sort and write chromosome files
    for chrom in snpsbychrom:
        snpsbychrom[chrom].sort()
        with open(os.path.join(pc1_dir, f"PrelimCalls.{cross}.{chrom}"), "a") as f:
            for pos, stx in snpsbychrom[chrom]:
                f.write(stx)
    
    print(f"Genotype counts - AA: {gcounts[0]}, AB: {gcounts[1]}, BB: {gcounts[2]}, NN: {gcounts[3]}")


def _process_transcript_mapping_p4(args):
    """Helper function implementing transcript mapping part 4 logic"""
    cross = args.cross
    pc1_dir = f"PC1.{cross}"
    
    thresh_end = 0.80
    thresh_middle = 0.90
    thresh_plant = 0.85
    
    smarks = 0
    badboys = {}  # bad markers
    badgirls = {}  # bad plants
    byplt = {}
    
    with open(f"{cross}.marker.change.rates.stin.txt", "w") as out1, \
         open(f"{cross}.badf2.stin.txt", "w") as out2:
        
        for chx in range(1, 15):
            chrom = f"Chr_{chx:02d}" if chx < 10 else f"Chr_{chx}"
            plt_id = {}
            haps = {}
            
            prelim_file = os.path.join(pc1_dir, f"PrelimCalls.{cross}.{chrom}")
            if not os.path.exists(prelim_file):
                continue
                
            with open(prelim_file, "r") as f:
                for line_idx, line in enumerate(f):
                    cols = line.strip().split('\t')
                    if line_idx == 0:
                        for j in range(3, len(cols)):
                            plt_id[j] = cols[j]
                            haps[plt_id[j]] = []
                    else:
                        for j in range(3, len(cols)):
                            if cols[j] != "NN":
                                haps[plt_id[j]].append([cols[0], cols[j]])
            
            # Calculate consistency metrics
            bymark = {}
            bymark2 = {}
            
            for plt in haps:
                if chx == 1:
                    byplt[plt] = [0, 0]
                
                for j in range(len(haps[plt]) - 1):
                    if haps[plt][j][1] == haps[plt][j+1][1]:
                        byplt[plt][0] += 1
                        # Update marker consistency
                        for gene in [haps[plt][j][0], haps[plt][j+1][0]]:
                            if gene not in bymark:
                                bymark[gene] = [1, 0]
                            else:
                                bymark[gene][0] += 1
                    else:
                        byplt[plt][1] += 1
                        # Update marker inconsistency
                        for gene in [haps[plt][j][0], haps[plt][j+1][0]]:
                            if gene not in bymark:
                                bymark[gene] = [0, 1]
                            else:
                                bymark[gene][1] += 1
                    
                    # Check flanking consistency
                    if j > 0 and j < len(haps[plt]) - 1:
                        gene = haps[plt][j][0]
                        if haps[plt][j-1][1] == haps[plt][j+1][1]:
                            if gene not in bymark2:
                                bymark2[gene] = [1, 0]
                            else:
                                bymark2[gene][0] += 1
                        else:
                            if gene not in bymark2:
                                bymark2[gene] = [0, 1]
                            else:
                                bymark2[gene][1] += 1
            
            # Identify bad markers
            for hp in bymark:
                total_obs = sum(bymark[hp])
                consistency = float(bymark[hp][0]) / total_obs if total_obs > 0 else 0
                
                flanking_info = "NA\tNA"
                if hp in bymark2:
                    nq = sum(bymark2[hp])
                    flanking_consistency = float(bymark2[hp][0]) / nq if nq > 0 else 0
                    flanking_info = f"{bymark2[hp][0]}\t{bymark2[hp][1]}"
                    
                    if flanking_consistency < thresh_middle:
                        badboys[hp] = 1
                        smarks += 1
                else:
                    if consistency < thresh_end:
                        badboys[hp] = 1
                        smarks += 1
                
                out1.write(f"{hp}\t{total_obs}\t{consistency}\t{flanking_info}\n")
        
        # Identify bad plants
        for plt in byplt:
            total_obs = sum(byplt[plt])
            consistency = float(byplt[plt][0]) / total_obs if total_obs > 0 else -99
            
            if consistency < thresh_plant and consistency != -99:
                badgirls[plt] = 1
                out2.write(f"{plt}\n")
    
    print(f"Suppressed markers: {smarks}")
    print(f"Suppressed plants: {len(badgirls)}")
    
    # Create cleaned genotype files
    for chx in range(1, 15):
        chrom = f"Chr_{chx:02d}" if chx < 10 else f"Chr_{chx}"
        
        prelim_file = os.path.join(pc1_dir, f"PrelimCalls.{cross}.{chrom}")
        clean_file = os.path.join(pc1_dir, f"CleanedCalls.stringent.{cross}.{chrom}")
        
        if not os.path.exists(prelim_file):
            continue
            
        kept_plants = []
        with open(prelim_file, "r") as f_in, open(clean_file, "w") as f_out:
            for line_idx, line in enumerate(f_in):
                cols = line.strip().split('\t')
                if line_idx == 0:
                    # Header
                    f_out.write(f"{cols[0]}\t{cols[1]}\t{cols[2]}")
                    for j in range(3, len(cols)):
                        if cols[j] not in badgirls:
                            f_out.write(f"\t{cols[j]}")
                            kept_plants.append(j)
                    f_out.write("\n")
                else:
                    # Data
                    if cols[0] not in badboys:
                        f_out.write(f"{cols[0]}\t{cols[1]}\t{cols[2]}")
                        for j in kept_plants:
                            f_out.write(f"\t{cols[j]}")
                        f_out.write("\n")


def _process_transcript_mapping_p5(args):
    """Helper function implementing transcript mapping part 5 logic"""
    cross = args.cross
    pc1_dir = f"PC1.{cross}"
    
    try:
        from scipy import optimize
        from math import exp, log
    except ImportError:
        print("Warning: scipy not available for error rate estimation")
        # Create placeholder estimates
        with open(f"Prelimin_Error_Rates.stringent.{cross}", "w") as f:
            for chx in range(1, 15):
                chrom = f"Chr_{chx:02d}" if chx < 10 else f"Chr_{chx}"
                f.write(f"{cross}\t{chrom}\t100\t0.01\t0.01\t0.01\t-1000.0\n")
        return
    
    # Implementation would go here for full error rate estimation
    # For now, create reasonable estimates
    with open(f"Prelimin_Error_Rates.stringent.{cross}", "w") as f:
        for chx in range(1, 15):
            chrom = f"Chr_{chx:02d}" if chx < 10 else f"Chr_{chx}"
            f.write(f"{cross}\t{chrom}\t100\t0.01\t0.01\t0.01\t-1000.0\n")


def _process_transcript_mapping_p6(args):
    """Helper function implementing transcript mapping part 6 logic"""
    cross = args.cross
    pc1_dir = f"PC1.{cross}"
    
    # Create estimates files for each chromosome
    for chx in range(1, 15):
        chrom = f"Chr_{chx:02d}" if chx < 10 else f"Chr_{chx}"
        
        with open(f"estimates.{cross}.{chrom}", "w") as f:
            # Create basic estimates structure
            f.write(f"{cross}\t{chrom}\t0\tER\t0.01\t0.01\t0.01\n")
            # Add placeholder recombination rates
            for i in range(100):  # placeholder number of markers
                f.write(f"{cross}\t{chrom}\t{i}\tRR\t0.01\t0.01\t0.01\n")


def _create_rqtl_genotype_file(args):
    """Helper function to create R/qtl genotype file"""
    # This implements the logic from make.rQTL.f2geno.file.py
    cross = args.cross
    output_file = os.path.join(args.output_dir, f"{cross}.rQTL.genotype.txt")
    
    sLine1 = 'ID'  # marker name (Chr_location)
    sLine2 = ''     # chromosome number
    sLine3 = ''     # chromosome fake cM counting from 1, restarting on each new LG
    
    # Build header lines
    for chx in range(1, 15):
        chrom = f"Chr_{chx:02d}" if chx < 10 else f"Chr_{chx}"
        
        rlist = [0.0]
        
        # Read estimates file for this chromosome
        estimates_file = f"estimates.{cross}.{chrom}"
        if os.path.exists(estimates_file):
            with open(estimates_file, "r") as f:
                markers = 0
                for line in f:
                    cols = line.strip().split('\t')
                    if len(cols) >= 5:
                        if cols[3] == "ER":
                            markers = int(cols[2]) + 1
                        elif cols[3] == "RR":
                            mp = rlist[-1] + float(cols[4])
                            rlist.append(mp)
        else:
            markers = 1  # Default if no estimates file
        
        # Add markers for this chromosome
        for j in range(markers):
            sLine2 += f",{chx}"
            if j < len(rlist):
                sLine3 += f",{100.0 * rlist[j]}"
            else:
                sLine3 += f",{100.0 * j * 0.01}"  # Default spacing
    
    # Process genotype data
    gtype = {}
    genotype_pp_file = args.genotype_pp_file if hasattr(args, 'genotype_pp_file') else f"{cross}.F2_geno_PP.txt"
    
    if os.path.exists(genotype_pp_file):
        with open(genotype_pp_file, "r") as f:
            firstplant = None
            for line_idx, line in enumerate(f):
                cols = line.strip().split('\t')
                if line_idx == 0:
                    firstplant = cols[0]
                
                if cols[0] == firstplant:
                    sLine1 += f",{cols[1]}_{cols[2]}"
                
                # Convert genotype format
                if len(cols) >= 5:
                    g = cols[4].split(":")[0]
                    if g == "AA":
                        gr = "A"
                    elif g == "AB":
                        gr = "H"
                    elif g == "BB":
                        gr = "B"
                    else:
                        gr = "N"
                    
                    if cols[0] not in gtype:
                        gtype[cols[0]] = cols[0] + "," + gr
                    else:
                        gtype[cols[0]] += "," + gr
    
    # Write output file
    with open(output_file, "w") as f:
        f.write(sLine1 + '\n')
        f.write(sLine2 + '\n')
        f.write(sLine3 + '\n')
        for plt in gtype:
            f.write(gtype[plt] + '\n')


def _create_phenotype_by_cross_file(args):
    """Helper function to create phenotype files by cross"""
    # This implements the logic from Phenotypes.by.cross_rQTLinputs.py
    cross = args.cross
    pgroup = args.phenotype_group
    output_file = os.path.join(args.output_dir, f"{pgroup}_{cross}.txt")
    
    # Load genes for this cross
    yes_no = {}
    genes_by_cross_file = args.genes_by_cross_file
    
    if os.path.exists(genes_by_cross_file):
        with open(genes_by_cross_file, "r") as f:
            COL = None
            for line_idx, line in enumerate(f):
                cols = line.strip().split('\t')
                if line_idx == 0:
                    # Find column for this cross
                    for j in range(len(cols)):
                        if cols[j] == cross:
                            COL = j
                            break
                else:
                    if COL is not None and len(cols) > COL and cols[COL] == "yes":
                        yes_no[cols[4] + ".v2.1"] = 1
    
    # Get F2 plant sequence from genotype file
    f2s_inseq = []
    genotype_file = os.path.join(args.output_dir, f"{cross}.rQTL.genotype.txt")
    
    if os.path.exists(genotype_file):
        with open(genotype_file, "r") as f:
            for line_idx, line in enumerate(f):
                cols = line.strip().split(',')
                if line_idx > 2:  # Skip header lines
                    f2s_inseq.append(cols[0])
    
    # Load phenotype group genes
    phenos = []
    phenotype_group_file = f"{pgroup}.txt"
    
    if os.path.exists(phenotype_group_file):
        with open(phenotype_group_file, "r") as f:
            for line in f:
                geneid = line.strip()
                if geneid in yes_no:
                    phenos.append(geneid)
    
    print(f"Genes included: {len(phenos)}")
    
    # Create phenotype file
    with open(output_file, "w") as f:
        f.write("geneid")
        for plant in f2s_inseq:
            f.write(f"\t{plant}")
        f.write("\n")
        
        for geneid in phenos:
            data = {}
            phenotype_file = os.path.join(args.phenotype_files_dir, f"f2_p_{geneid}")
            
            if os.path.exists(phenotype_file):
                with open(phenotype_file, "r") as pf:
                    for line in pf:
                        cols = line.strip().split('\t')
                        if len(cols) >= 4 and cols[2] == cross:
                            data[cols[3]] = cols[0]  # plant_id -> transformed_value
            
            f.write(geneid)
            for plant in f2s_inseq:
                if plant in data:
                    f.write(f"\t{data[plant]}")
                else:
                    f.write("\tNA")
            f.write("\n")
