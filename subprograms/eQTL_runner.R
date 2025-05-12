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

R_src_dir <- '.' # PLACEHOLDER: Set this to the actual path where your R scripts are located

for(file in list.files(path = R_src_dir, full.names = TRUE, pattern = "\\.R$"))
{
  source(file)
  message(paste("Sourced:", file))
  rm(file)
}

required_packages <- c('assertthat', 'doFuture', 'ggplot2', 'dplyr',
                       'GGally', 'tibble', 'jsonlite')