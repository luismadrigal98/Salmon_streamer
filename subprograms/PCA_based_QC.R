#' @description This script will take a file with teh sample based expression quantification
#' for each gene and it will check the degree of agreement between samples.
#' 
#' It will use PCA to condense the expression profiles of the genes and it will
#' use an arbitrary (defined by the user or selected interactively) number of PC
#' components to cluster the observations.
#' 
#' QC then will be executed over the cluster to eliminate discordant samples
#' 
#' @author Luis Javier Madrigal-Roca
#' 
#' @date 2025-05-04
#' 
#' @version 1.0.0
#' _____________________________________________________________________________

# THINGS I NEED TO PASS TO THIS SCRIPT

Data_file <- '/home/l338m483/scratch/MG_SF_SWB_project/RawSamples_forPCA'

# ******************************************************************************
# 1) Set up environment ----
# ______________________________________________________________________________

# 1.1) Source required scripts ----

R_src <- "./src/R_src"

for(file in list.files(path = R_src, full.names = T))
{
  source(file)
  rm(file)
}

# 1.2) Define the required packages ----

required_packages <- c('assertthat', 'doFuture', 'ggplot2', 'dplyr')

set_environment(parallel_backend = T, personal_seed = 1998, 
                automatic_download = T, required_pckgs = required_packages)

# ******************************************************************************
# 2) Load the data ----
# ______________________________________________________________________________

# 2.1) Read the data ----
Data <- read.table(file = Data_file, header = T, sep = '\t', 
                   stringsAsFactors = F)

traits <- Data$geneid

# 2.2) Reshape the data ----
Data <- Data |> dplyr::select(-geneid) |> t() |> as.data.frame()
names(Data) <- traits

# ******************************************************************************
# 3) Perform the PCA analysis ----
# ______________________________________________________________________________

PCA_res <- prcomp(x = Data)