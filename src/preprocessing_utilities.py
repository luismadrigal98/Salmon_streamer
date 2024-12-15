"""
Set of functions to preprocess the input files before using Salmon

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2024-12-14

"""

import subprocess
import logging
import os
import fileinput
from Bio import SeqIO

def homogenize_headers(input_file, output_file, only_chrom):
    """
    This function homogenizes the headers of a FASTA file by removing the characters after the first space in the header.
    This is done to avoid problems with the headers when using Salmon.

    Parameters
    ----------
    input_file : str
        Path to the input FASTA file
    output_file : str
        Path to the output FASTA file with homogenized headers

    Returns
    -------
    None

    """
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fasta'):
            record.id = record.id.split()[0]
            record.description = ''
            if only_chrom and ('scaffold' in record.id or 'contig' in record.id):
                continue
            SeqIO.write(record, out_handle, 'fasta')

def fasta_header_homogenizer(ref_genome, alt_genome, out_dir, only_chrom = True):
    """
    This function homogenizes the headers of the reference genome and the alternative genome.
    The headers are homogenized by removing the characters after the first space in the header.
    This is done to avoid problems with the headers when using Salmon.
    
    Parameters
    ----------
    ref_genome : str
        Path to the reference genome in fasta format
    alt_genome : str
        Path to the alternative genome in fasta format
    out_dir : str
        Path to the directory where the new files are stored
    only_chrom : bool
        If True, scaffolds and contigs are discarded from the fasta file.
        
    Returns
    -------
    str
        Path to the reference genome with homogenized headers
    str
        Path to the alternative genome with homogenized headers
    """
    ref_genome_homogenized = os.path.join(out_dir, os.path.basename(ref_genome).split('.')[0] + '_homogenized.fasta')
    alt_genome_homogenized = os.path.join(out_dir, os.path.basename(alt_genome).split('.')[0] + '_homogenized.fasta')

    # Process reference genome
    homogenize_headers(ref_genome, ref_genome_homogenized, only_chrom)

    # Process alternative genome
    homogenize_headers(alt_genome, alt_genome_homogenized, only_chrom)

    return ref_genome_homogenized, alt_genome_homogenized

def add_name_to_fasta_headers(fasta_file, name):
    """
    This function modifies the headers of a FASTA file in place by adding a specific name after each chromosome.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file to be modified
    name : str
        The name to be added after each chromosome in the header

    Returns
    -------
    None
    """
    with fileinput.FileInput(fasta_file, inplace=True, backup='.bak') as file:
        for line in file:
            if line.startswith('>'):
                line = line.strip() + f"_{name}\n"
            print(line, end='')

def get_contig_names_for_decoys(mg_genome, sf_genome, output_file):
    """
    Extract contig names from the given FASTA files and prepare a list of decoys.

    Parameters
    ----------
    mg_genome : str
        Path to the mg_genome FASTA file
    sf_genome : str
        Path to the sf_genome FASTA file
    output_file : str
        Path to the output file for decoy names

    Returns
    -------
    None
    """
    # Extract contig names
    cmd_grep = f"grep '^>' <(cat {mg_genome} {sf_genome}) | cut -d ' ' -f 1 > {output_file}"
    subprocess.run(cmd_grep, shell=True, executable='/bin/bash')

    # Clean up the decoy list
    cmd_sed = f"sed -i.bak -e 's/>//g' {output_file}"
    subprocess.run(cmd_sed, shell=True)

def concatenate_fasta_files(transcriptome, ref_genome, alt_genome, output_file):
    """
    Concatenate the transcriptome and the genomes into a single FASTA file.

    Parameters
    ----------
    transcriptome : str
        Path to the transcriptome FASTA file
    ref_genome : str
        Path to the reference genome FASTA file
    alt_genome : str
        Path to the alternative genome FASTA file
    output_file : str
        Path to the output file for the concatenated FASTA file

    Returns
    -------
    None
    """
    
    # Concatenate the files
    cmd_cat = f"cat {transcriptome} {ref_genome} {alt_genome} > {output_file}"
    subprocess.run(cmd_cat, shell=True, executable='/bin/bash')

def run_salmon_index(transcriptome_genome, decoy_file, threads, output_dir, output_name, salmon_options):
    """
    Create the Salmon index for the given transcriptome and decoy file.

    Parameters
    ----------
    transcriptome_genome : str
        Path to the concatenated transcriptome and genome file
    decoy_file : str
        Path to the decoy file
    threads : int
        Number of threads to use
    output_dir : str
        Path to the output directory for the index
    output_name : str
        Name of the output index
    salmon_options : list
        List of additional options for Salmon index creation

    Returns
    -------
    None
    """
    index_dir = os.path.join(output_dir, output_name)
    os.makedirs(index_dir, exist_ok=True)

    cmd = [
        "salmon", "index",
        "-t", transcriptome_genome,
        "-d", decoy_file,
        "-i", index_dir,
        "-p", str(threads)
    ] + salmon_options

    logging.info(f"Running Salmon index command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def master_script_generator(file, w_dir, job_dir, output_dir,
                            threads, memory,
                            salmon_index, 
                            quant_options = ['--noLengthCorrection -l A -p 1']):
    """
    Generate a master script to run Salmon quantification on the given FASTQ file.

    Parameters
    ----------
    file : str
        Path to the FASTQ file to be quantified
    w_dir : str
        Path to the working directory
    job_dir : str
        Path to the directory for the job scripts
    output_dir : str
        Path to the output directory for the quantification results
    threads : int
        Number of threads to use
    memory : str
        Memory to use
    salmon_index : str
        Path to the Salmon index
    quant_options : list, optional
        List of additional options for Salmon quantification (default is ['--noLengthCorrection'])

    Returns
    -------
    str
        Path to the job script

    """

    name = os.path.basename(file)

    # Remove the extension from the name

    name = name.split('.')[0]

    out_name = name + '_quant'

    job_script_path = os.path.join(job_dir, f"{name}.sh")
    with open(job_script_path, 'w') as out:
        out.write("#!/bin/bash\n")
        out.write(f"#SBATCH --job-name=salmon_{name}\n")
        out.write(f"#SBATCH --output={os.path.join(job_dir, f'salmon_{name}.out')}\n")
        out.write("#SBATCH --partition=sixhour,eeb,kucg,kelly\n")
        out.write(f"#SBATCH --ntasks={threads}\n")
        out.write(f"#SBATCH --cpus-per-task=1\n")
        out.write("#SBATCH --time=06:00:00\n")
        out.write("#SBATCH --mail-user=l338m483@ku.edu\n")
        out.write("#SBATCH --mail-type=END,FAIL\n")
        out.write(f"#SBATCH --mem-per-cpu={memory}g\n")
        out.write("\n")
        out.write("module load conda\n")
        out.write(r'eval "$(conda shell.bash hook)"')
        out.write("\n")
        out.write("conda activate salmon\n")
        out.write("\n")
        out.write(f"cd {w_dir}\n")
        out.write(f"salmon quant -i {salmon_index} ")
        out.write(f"{' '.join(quant_options)}\n")
        out.write(f"-r {file} -o {os.path.join(output_dir, out_name)} \n")
        

    return job_script_path

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