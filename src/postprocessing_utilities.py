def combine_results(output_dir):
    """
    Combine all Salmon quantification results into one table.

    Parameters
    ----------
    output_dir : str
        Path to the output directory containing the quantification results

    Returns
    -------
    None
    """
    # Change to the output directory
    os.chdir(output_dir)

    # Output column 5 (counts, not TPM column 4) as text file with ID as name
    cmd_cut = "for i in $(ls | grep quant); do cut -f5 $i/quant.sf > $i.txt; done"
    subprocess.run(cmd_cut, shell=True, executable='/bin/bash')

    # Replace first row in file with textfile name ie ID
    cmd_replace = r"for i in $(ls | grep quant.txt); do cat <(echo $i | perl -pe 's/.txt//') <(tail -n +2 $i) > ${i}_count.txt; done"
    subprocess.run(cmd_replace, shell=True, executable='/bin/bash')

    # Get gene names from one of the results, to be used as first column in compiled results file
    cmd_cut_gene_names = "cut -f1 s7_767-075_quant/quant.sf > 1count.txt"
    subprocess.run(cmd_cut_gene_names, shell=True, executable='/bin/bash')

    # Put everything together
    cmd_paste = "paste *count* > table.txt"
    subprocess.run(cmd_paste, shell=True, executable='/bin/bash')

    # Manually sort the table so that 2 species are on alternative rows, and remove "quant" from header row
    cmd_sort = r"perl -pe 's/.fastq//g ; s/_quant//g ; s/=/\t/ ; s/aName/aName\taName/' table.txt | sort -t$'\t' -k2,2 -k1,1 | perl -pe 's/ID\t/ID=/ ; s/\taName//' > table2.txt"
    subprocess.run(cmd_sort, shell=True, executable='/bin/bash')