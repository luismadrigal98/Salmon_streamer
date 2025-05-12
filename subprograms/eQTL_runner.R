#' @description This script runs the eQTL using the pipeline developed by Paris.
#' It also implements new functionalities and the possibility of running the analysis
#' using MatrixEQTL.
#' 
#' @author Luis Javier Madrigal-Roca & Paris Veltsos
#' 
#' @date 2025-05-12
#' 
#' #' @version 1.0.0
#' ___________________________________________________________________________________

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

# ARGUMENTS GO HERE

# ******************************************************************************
# 1) Set up environment ----
# ______________________________________________________________________________

# 1.1) Source required scripts ----

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".") # Fallback for non-sourced execution
R_src_dir <- file.path(script_dir, "..", "src", "R_src") # Relative to PCA_based_QC.R

if (!dir.exists(R_src_dir)) {
    # Fallback if script_dir is not determined correctly (e.g. when run from RStudio directly)
    # This assumes the script is run from the project root if the above fails.
    R_src_dir <- "./src/R_src"
}

message(paste("Sourcing R scripts from:", R_src_dir))
for(file in list.files(path = R_src_dir, full.names = TRUE, pattern = "\\.R$"))
{
  source(file)
  message(paste("Sourced:", file))
  rm(file)
}

required_packages <- c('qtl', 'qtl2', 'dplyr')

set_environment(parallel_backend = T, personal_seed = 1998, 
                automatic_download = T, required_pckgs = required_packages)

# ******************************************************************************
# 2) Load the data ----
# ______________________________________________________________________________

# 2.1) Import the genotype data ----
geno_data <- read.table(file = '/home/l338m483/scratch/MG_SF_SWB_project/R_directory/eQTL_Data/62.rQTL.genotype.txt',
                        header = T, sep = ',', stringsAsFactors = F)

# 2.2) Import the phenotype data ----
# Expression profiles per gene per sample

pheno_data <- read.table(file = '/home/l338m483/scratch/MG_SF_SWB_project/R_directory/eQTL_Data/gene.group1_62.txt',
                          header = T, sep = '\t', stringsAsFactors = F)

# ******************************************************************************
# 2) Load the data ----
# ______________________________________________________________________________

