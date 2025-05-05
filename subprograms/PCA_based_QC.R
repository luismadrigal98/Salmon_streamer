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