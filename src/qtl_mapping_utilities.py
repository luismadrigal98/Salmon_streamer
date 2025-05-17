"""
Utility functions for QTL (Quantitative Trait Loci) mapping analysis.

@author: Luis Javier Madrigal-Roca

"""

import logging
import pandas as pd
import subprocess
import logging
import time
from tqdm import tqdm
import getpass

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def get_gene_ids(phenofile_path, gene_id_column="geneid"):
    """
    Reads gene IDs from the phenotype file.

    Args:
        phenofile_path (str): Path to the phenotype file.
        gene_id_column (str): Column name for gene IDs in the phenotype file.
    Returns:
        list: A list of unique gene IDs.
    
    """
    try:
        # Assuming the phenotype file is tab-separated and has a header
        pheno_df = pd.read_csv(phenofile_path, sep='\t', header=0)
        if gene_id_column not in pheno_df.columns:
            logger.error(f"Gene ID column '{gene_id_column}' not found in {phenofile_path}.")
            logger.info(f"Available columns: {pheno_df.columns.tolist()}")
            return []
        return pheno_df[gene_id_column].unique().tolist()
    except Exception as e:
        logger.error(f"Error reading gene IDs from {phenofile_path}: {e}")
        return []
    
def get_active_job_count(user=None):
    """
    Get the number of active jobs for the specified user.
    
    Parameters:
    user (str): Username to check jobs for. If None, uses the current user.
    
    Returns:
    int: Number of active jobs (PENDING, RUNNING, CONFIGURING states)
    """
    try:
        # If no user specified, get current username
        if user is None:
            user = getpass.getuser()
        
        # Run squeue command to check jobs for this user
        cmd = f"squeue -u {user} -h -o '%T' | grep -E 'PENDING|RUNNING|CONFIGURING' | wc -l"
        output = subprocess.check_output(cmd, shell=True).decode().strip()
        
        # Parse the count
        count = int(output)
        logging.info(f"Current active job count for user {user}: {count}")
        return count
    except Exception as e:
        logging.error(f"Error checking job count: {e}")
        # Return a conservative estimate if we can't check
        return 0

def check_slurm_job_status(job_id):
    """
    Check the status of a SLURM job.
    """
    cmd = ["sacct", "-j", str(job_id), "--format=State", "--noheader", "--parsable2"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logging.error(f"Failed to check job status: {result.stderr}")
        return "UNKNOWN"
    
    output_lines = result.stdout.strip().split('\n')
    if not output_lines:
        return "UNKNOWN"
    
    # Get the first state (primary job state)
    state = output_lines[0].split('|')[0]
    
    # Map common states
    if "PENDING" in state:
        return "PENDING"
    elif "RUNNING" in state:
        return "RUNNING"
    elif "COMPLETED" in state:
        return "COMPLETED"
    elif "FAILED" in state or "TIMEOUT" in state:
        return "FAILED"
    elif "CANCELLED" in state:
        return "CANCELLED"
    else:
        return state

# Implementing a job submission function that respects concurrent job limits
def submit_individual_jobs(job_scripts, max_concurrent=50, wait_time=10, max_active_jobs=4900, poll_interval=60):
    """
    Submit individual job scripts while respecting a maximum active job limit.
    
    Parameters:
    job_scripts (list): List of paths to job script files to submit
    max_concurrent (int): Maximum number of jobs to submit in one batch before checking status
    wait_time (int): Seconds to wait between batches of submissions
    max_active_jobs (int): Maximum number of concurrent active jobs allowed in the cluster
    poll_interval (int): How often (in seconds) to check current job count when waiting
    
    Returns:
    list: List of submitted job IDs
    """
    submitted_job_ids = []
    remaining_scripts = job_scripts.copy()
    
    logging.info(f"Preparing to submit {len(job_scripts)} jobs with maximum {max_active_jobs} active jobs")
    
    with tqdm(total=len(job_scripts), desc="Submitting jobs") as progress_bar:
        while remaining_scripts:
            # Check current active job count
            active_jobs = get_active_job_count()
            
            # Calculate how many more jobs we can submit
            available_slots = max(0, max_active_jobs - active_jobs)
            
            if available_slots <= 0:
                logging.info(f"Active job limit reached ({active_jobs}/{max_active_jobs}). Waiting {poll_interval} seconds...")
                time.sleep(poll_interval)
                continue
            
            # Determine batch size - either max_concurrent or available slots, whichever is smaller
            batch_size = min(max_concurrent, available_slots, len(remaining_scripts))
            
            # Get the next batch of scripts
            batch = remaining_scripts[:batch_size]
            
            # Submit this batch
            for script_path in batch:
                try:
                    output = subprocess.check_output(f"sbatch {script_path}", shell=True).decode().strip()
                    
                    # Extract job ID
                    if "Submitted batch job" in output:
                        job_id = output.split()[-1]
                        # Verify it's a number
                        if job_id.isdigit():
                            submitted_job_ids.append(job_id)
                            logging.info(f"Submitted job {job_id} from {script_path}")
                        else:
                            logging.warning(f"Could not extract job ID from output: {output}")
                    else:
                        logging.warning(f"Unexpected sbatch output: {output}")
                    
                    # Small delay to avoid overwhelming the scheduler
                    time.sleep(0.5)
                    
                except subprocess.CalledProcessError as e:
                    logging.error(f"Error submitting job {script_path}: {e}")
            
            # Update remaining scripts and progress
            remaining_scripts = remaining_scripts[batch_size:]
            progress_bar.update(batch_size)
            
            # If we've submitted a batch, wait before checking again
            if batch:
                if remaining_scripts:
                    logging.info(f"Submitted {batch_size} jobs, waiting {wait_time}s before continuing...")
                    time.sleep(wait_time)
    
    logging.info(f"Total jobs submitted: {len(submitted_job_ids)}")
    return submitted_job_ids

def submit_individual_slurm_jobs(args, num_jobs=None):
    """
    Generate and submit individual SLURM jobs for hyperparameter search instead of array jobs.
    
    Parameters:
    args: Command line arguments object with job parameters
    num_jobs: Optional number of jobs to submit (defaults to args.max_trials)
    
    Returns:
    list: List of submitted job IDs
    """
    # If num_jobs is not provided, use max_trials
    if num_jobs is None:
        num_jobs = args.max_trials
    
    logging.info(f"Preparing {num_jobs} individual job scripts (one per trial)")
    
    # Generate job scripts
    job_scripts = generate_job_scripts(args, num_jobs)
    
    # Debug script contents if requested
    if hasattr(args, 'debug_slurm_submission') and args.debug_slurm_submission:
        logging.info(f"Generated {len(job_scripts)} job scripts")
        logging.info(f"First job script content sample:")
        with open(job_scripts[0], 'r') as f:
            logging.info(f.read())
    
    # Submit the job scripts and get job IDs with concurrent job limit
    logging.info(f"Submitting {len(job_scripts)} individual jobs with max active limit of {args.max_concurrent_jobs}...")
    job_ids = submit_individual_jobs(
        job_scripts, 
        max_concurrent=50,  # Submit in smaller batches
        wait_time=10, 
        max_active_jobs=args.max_concurrent_jobs,
        poll_interval=60  # Check active job count every minute
    )
    
    if job_ids:
        logging.info(f"Successfully submitted {len(job_ids)} jobs")
        return job_ids
    else:
        logging.error("Failed to submit any jobs")
        return []

def wait_for_jobs_completion(job_ids, max_wait_time=86400, poll_interval=60,
                            stringency= 0.75):
    """
    Wait for multiple SLURM jobs to complete.
    
    Parameters:
    job_ids (list): List of job IDs to wait for
    max_wait_time (int): Maximum wait time in seconds
    poll_interval (int): Time between status checks in seconds
    stringency (float): Percentage of jobs that must complete successfully to return True
    
    Returns:
    bool: True if all jobs completed successfully, False otherwise
    """
    if not job_ids:
        logging.warning("No job IDs provided to wait for")
        return True
    
    logging.info(f"Waiting for {len(job_ids)} jobs to complete...")
    start_time = time.time()
    
    # Add initial delay to allow jobs to register in SLURM
    time.sleep(20) ## This is to allow jobs to register in SLURM

    while time.time() - start_time < max_wait_time:
        incomplete_jobs = []
        failed_jobs = []
        unknown_jobs = []
        
        for job_id in job_ids:
            status = check_slurm_job_status(job_id)
            if status == "RUNNING" or status == "PENDING":
                incomplete_jobs.append(job_id)
            elif status == "FAILED" or status == "ERROR":
                failed_jobs.append(job_id)
            elif status == "COMPLETED":
                # Job is complete, don't add to incomplete_jobs
                logging.info(f"Job {job_id} completed successfully")
            else:
                # Only for truly unknown statuses
                unknown_jobs.append(job_id)
                incomplete_jobs.append(job_id)
        
        # Only consider jobs complete when we find none that are incomplete
        if not incomplete_jobs:
            if failed_jobs:
                logging.warning(f"{len(failed_jobs)} jobs failed: {failed_jobs}")
                return len(failed_jobs) / len(job_ids) < stringency
            else:
                logging.info("All jobs completed successfully")
                return True
        
        elapsed = time.time() - start_time
        logging.info(f"Waiting for {len(incomplete_jobs)} jobs to complete... ({elapsed:.1f}s elapsed)")
        time.sleep(poll_interval)
    
    logging.error(f"Exceeded maximum wait time of {max_wait_time}s")
    return False