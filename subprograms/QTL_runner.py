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

from qtl_mapping_utilities import get_gene_ids, submit_individual_jobs, wait_for_jobs_completion

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def main(args):
    # Ensure R script path is absolute
    r_script_full_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'eQTL_runner.R'))

    # Retrieve the sourcing directory for the R script
    source_dir = os.path.join(os.path.dirname(__file__), '../src/R_src')

    gene_ids = get_gene_ids(args.phenofile_path)
    if not gene_ids:
        logger.error("No gene IDs found. Exiting.")
        return

    logger.info(f"Found {len(gene_ids)} genes to process.")

    # Create base output directory if it doesn't exist
    if not os.path.exists(args.outdir_base):
        os.makedirs(args.outdir_base, exist_ok=True)
        logger.info(f"Created base output directory: {args.outdir_base}")

    # Generate all job scripts before submission
    job_scripts = []
    
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
            f"--crosstype={args.crosstype}",
            f"--source_dir={source_dir}"
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
        
        job_scripts.append(slurm_script_path)

    # Submit the jobs using the advanced job management system
    try:
        logger.info(f"Prepared {len(job_scripts)} job scripts. Beginning submission...")
        submitted_job_ids = submit_individual_jobs(
            job_scripts,
            max_concurrent=args.submission_batch_size, 
            wait_time=args.submission_wait_time,
            max_active_jobs=args.max_concurrent_jobs,
            poll_interval=args.poll_interval
        )
        
        logger.info(f"Successfully submitted {len(submitted_job_ids)} out of {len(job_scripts)} jobs")
        
        # Optionally wait for jobs to complete
        if args.wait_for_completion:
            logger.info(f"Waiting for all submitted jobs to complete (max wait: {args.max_wait_time} seconds)")
            success = wait_for_jobs_completion(
                submitted_job_ids, 
                max_wait_time=args.max_wait_time,
                poll_interval=args.poll_interval,
                stringency=args.job_completion_stringency
            )
            
            if success:
                logger.info("All eQTL jobs completed successfully.")
            else:
                logger.warning("Some eQTL jobs failed or timed out. Check SLURM logs for details.")

    except Exception as e:
        logger.error(f"Error in job submission process: {e}")