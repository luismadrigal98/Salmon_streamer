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

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript PCA_based_QC.R <Data_file_path> <label_rules_file_path> <iqr_multiplier> <output_filtered_data_name> <genes_as_rows_TRUE_FALSE> <output_dir_path> <pc_to_retain>", call. = FALSE)
}

Data_file <- args[1]
label_rules_file_path <- args[2]
iqr_multiplier <- as.numeric(args[3])
output_filtered_data_name <- args[4]
genes_as_rows <- as.logical(args[5]) # Expect "TRUE" or "FALSE"
output_dir <- args[6]
pc_to_retain <- as.integer(args[7]) # Number of PCs to retain

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste("Created output directory:", output_dir))
}
# Ensure Results/Plots subdirectory exists for plots
plots_dir <- file.path(output_dir, "Plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
  message(paste("Created plots directory:", plots_dir))
}

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

# 1.2) Define the required packages ----

required_packages <- c('assertthat', 'doFuture', 'ggplot2', 'dplyr',
                       'GGally', 'tibble', 'jsonlite')

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
# 4) Outlier removal ----
# ______________________________________________________________________________

# 4.1) Preserve observations as a function of absolute deviation from the mean ----

# Merge the two grouping variables

PCA_data$Group_merged <- paste0(PCA_data$Source, "_", PCA_data$Group) |>
  as.factor()

PCA_data <- PCA_data |> dplyr::group_by(Group_merged) |>
  dplyr::filter(!(abs(value - median(value)) > 2 * sd(value)))

# 4.2) Removal of the outliers ---

# Convert rownames to an explicit column BEFORE the dplyr pipeline
PCA_data <- PCA_data |>
  tibble::rownames_to_column(var = "SampleID")

PCA_data_filtered_iqr <- PCA_data |>
  dplyr::group_by(Group_merged) |>
  dplyr::filter(
    dplyr::if_all(
      .cols = dplyr::starts_with("PC"), # Selects all PCs
      .fns = ~ mark_outliers_IQR(., multiplier = iqr_multiplier)
    )
  ) |>
  dplyr::ungroup()

# 4.3) Replot the PCRs ----

PCA_plotter(PCA_data_filtered_iqr, x = 'PC1', y = 'PC2', color = 'Group_merged',
            filename = "./Results/Plots/PCA_1_vs_2_filtered.pdf")

PCA_plotter(PCA_data_filtered_iqr, x = 'PC2', y = 'PC3', color = 'Group_merged', 
            filename = "./Results/Plots/PCA_2_vs_3_filtered.pdf")

PCA_plotter(PCA_data_filtered_iqr, x = 'PC1', y = 'PC3', color = 'Group_merged', 
            filename = "./Results/Plots/PCA_1_vs_3_filtered.pdf")

# 5) Save the filtered set ---

# 5.1) Add a sample ID column ----

if (nrow(PCA_data_filtered_iqr) > 0) {
    subset_original_data <- Data[PCA_data_filtered_iqr$SampleID, ] # Data has samples as rows here

    output_file_path <- file.path(output_dir, output_filtered_data_name)

    if(genes_as_rows) {
      # Transpose back so genes are rows, samples are columns
      write.table(x = t(subset_original_data), file = output_file_path, sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA) # NA for top-left cell
      message(paste("Filtered data (genes as rows) saved to:", output_file_path))
    } else {
      # Samples as rows, genes as columns
      write.table(x = subset_original_data, file = output_file_path, sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
      message(paste("Filtered data (samples as rows) saved to:", output_file_path))
    }
} else {
    message("No data remaining after filtering. Filtered data file will not be created.")
    # Optionally, create an empty file or a file with headers
    output_file_path <- file.path(output_dir, output_filtered_data_name)
    file.create(output_file_path) # Create an empty file
    message(paste("Empty filtered data file created at:", output_file_path))
}

message("PCA QC script finished.")