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

required_packages <- c('assertthat', 'doFuture', 'ggplot2', 'dplyr',
                       'GGally')

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

PCA_res <- prcomp(x = Data, center = T, scale = T)

PCA_data <- process_PCA_results(PCA_res, 5) # This should also return the composite variance

# Categorize the samples

color_grouping_rules <- list(
    "SF_Parental" = "SFcross.*_.*P",
    "SWB_Parental" = "SWBcross.*_.*P",
    "SF_Cross" = "^SFcross", # Use ^ to anchor at the beginning
    "SWB_Cross" = "^SWBcross" # Use ^ to anchor at the beginning
  )

group_grouping_rules <- list("D" = "_D.",
                    "J" = "_J.",
                    "767" = "_767.")

PCA_data$Source <- assign_groups_from_labels(default_group = 'Unknown',
                                            rules = color_grouping_rules,
                                            labels = rownames(PCA_data)) |>
  as.factor()

PCA_data$Group <- assign_groups_from_labels(default_group = 'Unknown',
                                            rules = group_grouping_rules,
                                            labels = rownames(PCA_data)) |>
  as.factor()

# Visual exploration before QC

PCA_plotter(PCA_data, x = 'PC1', y = 'PC2', color = 'Source',
            group = 'Group', filename = "./Results/Plots/PCA_1_vs_2.pdf")

PCA_plotter(PCA_data, x = 'PC2', y = 'PC3', color = 'Source',
            group = 'Group', filename = "./Results/Plots/PCA_2_vs_3.pdf")

PCA_plotter(PCA_data, x = 'PC1', y = 'PC3', color = 'Source',
            group = 'Group', filename = "./Results/Plots/PCA_1_vs_3.pdf")

# ******************************************************************************
# 4) QC based on clustering ----
# ______________________________________________________________________________

# 4.1) Preserve observations as a function of absolute deviation from the mean ----

# Merge the two grouping variables

PCA_data$Group_merged <- paste0(PCA_data$Source, "_", PCA_data$Group)

PCA_data <- PCA_data |> dplyr::group_by(Group_merged) |>
  dplyr::filter(!abs(value - mean(value) > 2 * sd(value)))