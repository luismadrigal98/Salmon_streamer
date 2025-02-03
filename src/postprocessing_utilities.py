"""
Set of functions to postprocess the output files after using Salmon

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2024-12-15

"""

import subprocess
import os
import pandas as pd

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