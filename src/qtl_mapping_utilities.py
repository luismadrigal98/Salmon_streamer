"""
Utility functions for QTL (Quantitative Trait Loci) mapping analysis.

@author: Luis Javier Madrigal-Roca

"""

import logging
import pandas as pd

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