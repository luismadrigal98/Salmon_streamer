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

PCA_res <- prcomp(x = Data, center = T, scale = F)

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
                                            labels = rownames(PCA_data))

PCA_data$Group <- assign_groups_from_labels(default_group = 'Unknown',
                                            rules = group_grouping_rules,
                                            labels = rownames(PCA_data))

# Visual exploration (multipanel scatterplot, all vs all)

pairs_plot <- ggpairs(
  PCA_data,
  columns = 1:5, # Specify that only PC1-PC5 should form the axes of the pairs
  title = "Scatterplot Matrix of Principal Components",
  upper = list(continuous = wrap("points", mapping = aes(color = Group), 
                                 alpha = 0.5, size = 0.8)), # Map 'Group' column to color
  lower = list(continuous = wrap("points", mapping = aes(color = Source), 
                                 alpha = 0.5, size = 0.8)), # Map 'Source' column to color
  diag = list(continuous = wrap("densityDiag", mapping = aes(fill = Group), 
                                alpha = 0.5)) # Optionally color diagonal by Group as well
) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = 'black'),
    panel.grid = element_blank()
  )

# Ensure the directory exists before saving
plots_dir <- "./Results/Plots/"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

ggsave(paste0(plots_dir, 'PCA_screenen_before_QC.pdf'), plot = pairs_plot, width = 12, height = 12)