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
    This is done to avoid problems with the header formats when using Salmon.

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
    This is done to avoid problems with different header formats when using Salmon.
    
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
    ref_genome_homogenized = os.path.join(out_dir, os.path.basename(ref_genome).rsplit('.', 1)[0] + '_homogenized.fasta') # rsplit removes the extension of the input file
    alt_genome_homogenized = os.path.join(out_dir, os.path.basename(alt_genome).rsplit('.', 1)[0] + '_homogenized.fasta') # rsplit removes the extension of the input file

    # Process reference genome
    homogenize_headers(ref_genome, ref_genome_homogenized, only_chrom)

    # Process alternative genome
    homogenize_headers(alt_genome, alt_genome_homogenized, only_chrom)

    return ref_genome_homogenized, alt_genome_homogenized

def add_name_to_fasta_headers(fasta_file, name, remove_bak = True):
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
    
    # If no error occurs, remove the backup file if flag set to True (default)
    if remove_bak:
        os.remove(fasta_file + '.bak')

def get_contig_names_for_decoys(ref_genome, alt_genome = None, output_file = 'decoys.txt'):
    """
    Extract contig names from the given FASTA files and prepare a list of decoys.

    Parameters
    ----------
    ref_genome : str
        Path to the ref_genome FASTA file
    alt_genome : str
        Path to the alt_genome FASTA file. It is none by default.
    output_file : str
        Path to the output file for decoy names

    Returns
    -------
    None
    """
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(ref_genome)

    # Extract contig names
    try:
        if alt_genome is None:
            cmd_grep = f"grep '^>' {ref_genome} | cut -d ' ' -f 1 > {output_file}"
        else:
            cmd_grep = f"grep '^>' <(cat {ref_genome} {alt_genome}) | cut -d ' ' -f 1 > {output_file}"

        print(f"Running command: {cmd_grep}")
        subprocess.run(cmd_grep, shell=True, executable='/bin/bash', check=True)

        # Check if the output file is created and not empty
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            raise FileNotFoundError(f"The decoy file {output_file} was not created or is empty.")
        
        # Clean up the decoy list
        cmd_sed = f"sed -i.bak -e 's/>//g' {output_file}"
        print(f"Running command: {cmd_sed}")
        subprocess.run(cmd_sed, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}")
        raise
    except Exception as e:
        print(f"An error occurred: {e}")
        raise

def concatenate_fasta_files(transcriptome, ref_genome, alt_genome = None, output_file = 'transcriptome_genome.fasta'):
    """
    Concatenate the transcriptome and the genomes into a single FASTA file.

    Parameters
    ----------
    transcriptome : str
        Path to the transcriptome FASTA file
    ref_genome : str
        Path to the reference genome FASTA file
    alt_genome : str
        Path to the alternative genome FASTA file. It is None by default.
    output_file : str
        Path to the output file for the concatenated FASTA file

    Returns
    -------
    None
    """
    
    if alt_genome is None:
        # Concatenate the files
        cmd_cat = f"cat {transcriptome} {ref_genome} > {output_file}"
    else:
        cmd_cat = f"cat {transcriptome} {ref_genome} {alt_genome} > {output_file}"

    # Run the command
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
    
    # Verify that the index was created successfully
    version_info_path = os.path.join(index_dir, "versionInfo.json")
    if not os.path.exists(version_info_path):
        raise FileNotFoundError(f"Salmon index creation failed: {version_info_path} does not exist")
    
    logging.info(f"Salmon index created successfully at {index_dir}")

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
        out.write("#SBATCH --time=05:59:00\n")
        out.write("#SBATCH --mail-user=l338m483@ku.edu\n")
        out.write("#SBATCH --mail-type=FAIL\n")
        out.write(f"#SBATCH --mem-per-cpu={memory}g\n")
        out.write("\n")
        out.write("module load conda\n")
        out.write(r'eval "$(conda shell.bash hook)"')
        out.write("\n")
        out.write("conda activate salmon\n")
        out.write("\n")
        out.write(f"cd {w_dir}\n")
        out.write("\n")
        out.write("# Verify salmon index exists before running quantification\n")
        out.write(f"if [ ! -f {salmon_index}/versionInfo.json ]; then\n")
        out.write(f"    echo \"ERROR: Salmon index not found at {salmon_index}/versionInfo.json\"\n")
        out.write("    exit 1\n")
        out.write("fi\n")
        out.write("\n")
        out.write(f"salmon quant -i {salmon_index} ")
        out.write(f"{' '.join(quant_options)}")
        out.write(f" -r {file} -o {os.path.join(output_dir, out_name)} \n")

    return job_script_path