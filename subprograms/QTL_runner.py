"""
Main subprogram for running QTL analysis.

This script is design to detect individual genes, adn set up individual QTL analysis per genes.

In the future it will be extended to run these analysis also per different crosses defined in a cross file.

@author: Luis Javier Madrigal-Roca
@date: 2025-05-12

"""

import os
import logging
import subprocess
import sys
import time

# Add the path to the subprograms directory to the system path
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

from qtl_mapping_utilities import get_gene_ids

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def main(args):
    # Ensure R script path is absolute
    r_script_full_path = os.path.abspath(args.r_script_path)
    if not os.path.exists(r_script_full_path):
        logger.error(f"R script not found at: {r_script_full_path}")
        return

    gene_ids = get_gene_ids(args.phenofile_path)
    if not gene_ids:
        logger.error("No gene IDs found. Exiting.")
        return

    logger.info(f"Found {len(gene_ids)} genes to process.")

    # Create base output directory if it doesn't exist
    if not os.path.exists(args.outdir_base):
        os.makedirs(args.outdir_base, exist_ok=True)
        logger.info(f"Created base output directory: {args.outdir_base}")

    job_count = 0
    for gene_id in gene_ids:
        # Sanitize gene_id for use in filenames and job names
        safe_gene_id = "".join(c if c.isalnum() or c in ('.', '_', '-') else '_' for c in gene_id)

        # Specific output directory for this gene
        gene_outdir = os.path.join(args.outdir_base, f"gene_{safe_gene_id}")
        if not os.path.exists(gene_outdir):
            os.makedirs(gene_outdir, exist_ok=True)

        job_name = f"eQTL_{safe_gene_id}"
        log_file = os.path.join(gene_outdir, f"{job_name}.log")
        error_file = os.path.join(gene_outdir, f"{job_name}.err")

        r_command_parts = [
            args.rscript_executable,
            r_script_full_path,
            f"--genfile_path={os.path.abspath(args.genfile_path)}",
            f"--phenofile_path={os.path.abspath(args.phenofile_path)}",
            f"--outdir={os.path.abspath(gene_outdir)}",
            f"--gene_id={gene_id}",
            f"--qtlmethod={args.qtlmethod}",
            f"--modeltype={args.modeltype}",
            f"--permnum={args.permnum}",
            f"--crosstype={args.crosstype}"
        ]
        if args.covfile_path:
            r_command_parts.append(f"--covfile_path={os.path.abspath(args.covfile_path)}")

        r_command = " ".join(r_command_parts)

        slurm_script_content = f"""#!/bin/bash
        #SBATCH --job-name={job_name}
        #SBATCH --partition={args.partition}
        #SBATCH --nodes={args.nodes}
        #SBATCH --ntasks-per-node={args.ntasks_per_node}
        #SBATCH --cpus-per-task={args.cpus_per_task}
        #SBATCH --mem-per-cpu={args.mem_per_cpu}
        #SBATCH --time={args.time_limit}
        #SBATCH --output={log_file}
        #SBATCH --error={error_file}

        echo "Starting eQTL analysis for gene: {gene_id}"
        echo "Job ID: $SLURM_JOB_ID"
        echo "Running on host: $(hostname)"
        echo "R command: {r_command}"

        # Load R module if necessary (common on HPC clusters)
        # module load R/your_R_version 

        {r_command}

        echo "Finished eQTL analysis for gene: {gene_id}"
        """
        
        slurm_script_path = os.path.join(gene_outdir, f"{job_name}.slurm")
        with open(slurm_script_path, "w") as f:
            f.write(slurm_script_content)

        # Submit the job
        try:
            # Check current number of running/pending jobs for the user (optional, advanced)
            # For now, just submit up to max_jobs
            if job_count < args.max_jobs :
                submission_command = ["sbatch", slurm_script_path]
                result = subprocess.run(submission_command, capture_output=True, text=True, check=True)
                logger.info(f"Submitted job for gene {gene_id}: {result.stdout.strip()} (Script: {slurm_script_path})")
                job_count += 1
                time.sleep(args.job_submission_delay) # Small delay to avoid overwhelming the scheduler
            else:
                logger.warning(f"Reached max job submission limit ({args.max_jobs}). Waiting for jobs to complete or increase limit.")
                # Basic wait strategy: check job queue or simply pause.
                # A more robust solution would involve monitoring SLURM queue.
                # For now, we'll just stop submitting more.
                # You might want to implement a loop that waits and checks `squeue`
                break # Or implement a waiting loop

        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to submit job for gene {gene_id}: {e}")
            logger.error(f"sbatch stdout: {e.stdout}")
            logger.error(f"sbatch stderr: {e.stderr}")
        except FileNotFoundError:
            logger.error("sbatch command not found. Ensure SLURM tools are in your PATH.")
            return

    logger.info(f"Finished submitting jobs. Total submitted: {job_count}")