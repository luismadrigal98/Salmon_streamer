#' @description This script runs the eQTL using the pipeline developed by Paris.
#' It also implements new functionalities and the possibility of running the analysis
#' using MatrixEQTL.
#' 
#' @author Luis Javier Madrigal-Roca & Paris Veltsos
#' 
#' @date 2025-05-12
#' 
#' @version 1.0.0
#' ___________________________________________________________________________________

# ******************************************************************************
# 1) Define and Parse Command Line Arguments (using hardcoded for now) ----
# ______________________________________________________________________________

# --- Argument Parser Setup ---
parser <- argparse::ArgumentParser(description = "Run eQTL analysis for specified genes using R/qtl")
parser$add_argument("--genfile_path", type = "character", required = TRUE,
                    help = "Path to the genotype file (e.g., CSV format for R/qtl)")
parser$add_argument("--phenofile_path", type = "character", required = TRUE,
                    help = "Path to the phenotype (expression) file (e.g., CSV format for R/qtl)")
parser$add_argument("--outdir", type = "character", required = TRUE,
                    help = "Directory to save output files")
parser$add_argument("--gene_id", type = "character", required = FALSE, default = NULL,
                    help = "Specific gene ID to process. If NULL, all genes are processed (for local testing).")
parser$add_argument("--qtlmethod", type = "character", default = "mr",
                    help = "QTL mapping method (e.g., 'mr', 'hk', 'em') [default: %(default)s]")
parser$add_argument("--modeltype", type = "character", default = "normal",
                    help = "Model type for scanone (e.g., 'normal', 'binary') [default: %(default)s]")
parser$add_argument("--permnum", type = "integer", default = 1000,
                    help = "Number of permutations for significance testing [default: %(default)s]")
parser$add_argument("--crosstype", type = "character", default = "f2",
                    help = "Cross type for read.cross (e.g., 'f2', 'bc', 'risib') [default: %(default)s]")
parser$add_argument("--covfile_path", type = "character", required = FALSE, default = NULL,
                    help = "Path to the covariate file (optional, e.g., CSV)")
parser$add_argument("--source_dir", type = "character", required = FALSE, default = NULL,
                    help = "Directory containing the source scripts")

# 2.1) Retrieve the arguments ----
args <- parser$parse_args()

# ******************************************************************************
# 2) Set up environment ----
# ______________________________________________________________________________

# 2.1) Load required packages ----

R_src_dir <- args$source_dir

message(paste("Sourcing R scripts from:", R_src_dir))
for(file in list.files(path = R_src_dir, full.names = TRUE, pattern = "\\.R$"))
{
  source(file)
  message(paste("Sourced:", file))
  rm(file)
}

required_packages <- c('doParallel', 'qtl', 'dplyr', 'doFuture',
                       'argparse')

set_environment(parallel_backend = F, personal_seed = 1998, 
                automatic_download = T, required_pckgs = required_packages)

if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

# ******************************************************************************
# 3) Load and Prepare Raw Data ----
# ______________________________________________________________________________

# 3.1) Import the genotype data ----

message("Loading genotype data from: ", args$genfile_path)
geno_data_raw <- read.table(file = args$genfile_path,
                        header = T, sep = ',', stringsAsFactors = F, check.names = FALSE)
message("Genotype data loaded. Dimensions: ", paste(dim(geno_data_raw), collapse = " x "))

# 3.2) Import the phenotype data ----
message("Loading phenotype data from: ", args$phenofile_path)
pheno_data_raw <- read.table(file = args$phenofile_path,
                          header = T, sep = '\t', stringsAsFactors = F, check.names = FALSE)
message("Phenotype data loaded. Dimensions: ", paste(dim(pheno_data_raw), collapse = " x "))

# 3.3) Prepare genotype data for R/qtl ----
# The geno_data_raw seems to be in the correct format for direct CSV writing
# First column is ID, first data row has chr, second data row has pos.
# Ensure the first column name is suitable or explicitly handle it if needed.
# For read.cross with format="csvs", the first column is implicitly the ID.
# Let's rename the first column to 'id' for clarity, though R/qtl might not strictly require it.

geno_for_qtl <- geno_data_raw

# 3.4) Prepare phenotype data for R/qtl ----

message("Transforming phenotype data...")

# Gene IDs are in the first column, samples are in subsequent column headers
# Transpose so that samples are rows and genes are columns

pheno_transformed <- as.data.frame(t(pheno_data_raw[, -which(colnames(pheno_data_raw) == "geneid")]))
colnames(pheno_transformed) <- pheno_data_raw$geneid

# If a specific gene_id is provided, subset the phenotype data here
if (!is.null(args$gene_id)) {
  if (args$gene_id %in% colnames(pheno_transformed)) {
    message(paste("Processing single specified gene:", args$gene_id))
    pheno_transformed <- pheno_transformed[, args$gene_id, drop = FALSE] # Keep as data.frame
  } else {
    stop(paste("Specified gene_id '", args$gene_id, "' not found in phenotype data.", sep=""))
  }
}

# Add sample IDs as the first column, harmonizing them with genotype IDs (e.g. s3_62.168 -> s3_62-168)
pheno_transformed$id <- rownames(pheno_transformed)
pheno_transformed$id <- gsub("\\.", "-", pheno_transformed$id) # Replace dot with dash

# Reorder to make 'id' the first column
pheno_for_qtl <- pheno_transformed[, c("id", setdiff(names(pheno_transformed), "id"))] |>
  dplyr::rename(ID = id)
message("Phenotype data transformed.")

# 3.5) Write data to temporary CSV files for read.cross ----
message("Writing temporary CSV files for R/qtl...")
temp_geno_file <- tempfile(pattern = "geno_", fileext = ".csv", tmpdir = args$outdir)
temp_pheno_file <- tempfile(pattern = "pheno_", fileext = ".csv", tmpdir = args$outdir)

write.csv(geno_for_qtl, temp_geno_file, row.names = FALSE, quote = FALSE)
write.csv(pheno_for_qtl, temp_pheno_file, row.names = FALSE, quote = FALSE)
message(paste("Temporary genotype CSV:", temp_geno_file))
message(paste("Temporary phenotype CSV:", temp_pheno_file))

# ******************************************************************************
# 4) Create R/qtl cross object ----
# ______________________________________________________________________________
message("Reading data into R/qtl cross object...")

# Ensure your genotype codes are "A", "H", "B" for an F2 cross.
# "A" and "B" are homozygous, "H" is heterozygous.
# estimate.map=TRUE will use the positions in your file and try to re-estimate the map.
# If your map is final, you might set estimate.map=FALSE.
# The original scripts used estimate.map=T.

cross_obj <- read.cross(format = "csvs",
                        dir = dirname(temp_geno_file), # or args$outdir
                        genfile = basename(temp_geno_file),
                        phefile = basename(temp_pheno_file),
                        genotypes = c("A", "H", "B"), # Specify genotype codes
                        alleles = c("A", "B"),       # Specify allele codes
                        estimate.map = TRUE,
                        crosstype = args$crosstype) # IMPORTANT: Verify "f2" or set appropriately

message("Cross object created successfully.")
summary(cross_obj)

# Clean up temporary files
unlink(temp_geno_file)
unlink(temp_pheno_file)
message("Temporary CSV files removed.")

# ******************************************************************************
# 5) QTL Data Processing and Analysis ----
# ______________________________________________________________________________

# 5.1) Initial map processing ----

message("Processing genetic map...")
mapthis <- jittermap(cross_obj) # Add small amount of noise to markers at same position

# plotMap(mapthis) # Optional: visualize map

# If estimate.map=TRUE was used in read.cross, the map is already estimated.
# If you need further refinement or used estimate.map=FALSE:
# mapthis <- est.rf(mapthis)
# newmap <- est.map(mapthis, error.prob=0.01) # As in qtl_01a_setup_short.r
# mapthis <- replace.map(mapthis, newmap)

# 5.2) Calculate genotype probabilities ----

message("Calculating genotype probabilities...")

mapthis <- calc.genoprob(mapthis,
                        step = 1, # Adjust step size as needed
                        error.prob = 0.001,
                        map.function = "kosambi") # Or other map function

message("Genotype probabilities calculated.")

# 5.3) Simulate genotypes (for imputation or fine-mapping) ----
# This step is REQUIRED for effectscan() to work
message("Simulating genotypes...")
mapthis <- sim.geno(mapthis,
                   n.draws = 32, # Number of simulation draws
                   step = 0, # Use 0 for positions at markers
                   off.end = 0.0,
                   error.prob = 1.0e-4,
                   map.function = "kosambi",
                   stepwidth = "fixed")
message("Genotypes simulated.")

# 5.4) Perform eQTL scan for each phenotype ----
message("Starting eQTL analysis...")
pheno_names_in_mapthis <- phenames(mapthis) # These are the gene names available in the cross object

# Determine which phenotypes to process
if (!is.null(args$gene_id)) {
    if (args$gene_id %in% pheno_names_in_mapthis) {
        pheno_names_to_process <- args$gene_id
        message(paste("Targeting single gene from command line:", args$gene_id))
    } else {
        stop(paste("Specified gene_id '", args$gene_id, "' not found in the processed cross object's phenotypes. Available:", paste(pheno_names_in_mapthis, collapse=", ")))
    }
} else {
    pheno_names_to_process <- pheno_names_in_mapthis
    message("No specific gene_id provided, processing all genes found in the phenotype data.")
}

all_results_list <- list()

# Load covariates if you have them
covariates_matrix <- NULL  # Default to NULL
if (!is.null(args$covfile_path)) {
  tryCatch({
    if (file.exists(args$covfile_path)) {
      message(paste("Loading covariates from:", args$covfile_path))
      cov_df <- read.csv(args$covfile_path, row.names = 1, check.names = FALSE)
      covariates_matrix <- as.matrix(cov_df)
      message("Covariates loaded and formatted.")
    } else {
      message(paste("Covariate file not found:", args$covfile_path))
    }
  }, error = function(e) {
    message(paste("Error loading covariates:", e$message))
    covariates_matrix <<- NULL  # Ensure it's set to NULL on error
  })
} else {
  message("No covariate file path provided.")
}

# Ensure covariates_matrix exists before entering the loop
if (!exists("covariates_matrix")) {
  covariates_matrix <- NULL
}

for (current_pheno_name in pheno_names_to_process) { # Modified loop
  message(paste("Analyzing phenotype:", current_pheno_name))
  current_pheno_col <- which(phenames(mapthis) == current_pheno_name)

  # Perform QTL scan
  scanone_result <- scanone(mapthis,
                            pheno.col = current_pheno_col,
                            addcovar = covariates_matrix, # Add if you have covariates
                            method = args$qtlmethod,
                            model = args$modeltype)

  # Perform permutations for significance thresholds
  message(paste("Running permutations for", current_pheno_name, "(n =", args$permnum, ")"))
  perm_result <- scanone(mapthis,
                         pheno.col = current_pheno_col,
                         addcovar = covariates_matrix, # Add if you have covariates
                         method = args$qtlmethod,
                         model = args$modeltype,
                         n.perm = args$permnum)

  # Get LOD threshold (e.g., 5% significance)
  lod_threshold_val <- summary(perm_result, alpha = 0.05)[1] # Standard LOD threshold

  # Get summary of significant QTLs
  significant_qtls <- summary(scanone_result, perms = perm_result, alpha = 0.05, pvalues = TRUE)

  if (nrow(significant_qtls) > 0) {
    significant_qtls$pheno <- current_pheno_name
    significant_qtls$lod_threshold.05 <- lod_threshold_val
    # You might want to add marker names if they are not row names
    # significant_qtls$marker <- rownames(significant_qtls) # If summary output doesn't have it as a column

    # Store or write results
    all_results_list[[current_pheno_name]] <- significant_qtls
    output_file_pheno <- file.path(args$outdir, paste0(gsub("[^a-zA-Z0-9_.-]", "_", current_pheno_name), "_", args$qtlmethod, "_eQTLs.tsv"))
    write.table(significant_qtls, file = output_file_pheno, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("Significant QTL results for", current_pheno_name, "saved to", output_file_pheno))

    # Optional: Effect scan for significant QTLs or all markers
    # effect_scan_result <- effectscan(mapthis, pheno.col = current_pheno_col, get.se = TRUE, draw = FALSE)
    # ... further processing of effect_scan_result ...

  } else {
    message(paste("No significant QTLs found for", current_pheno_name, "at the 0.05 LOD threshold."))
    # Optionally save all scanone results even if not significant
    all_scan_df <- as.data.frame(scanone_result)
    all_scan_df$marker <- rownames(all_scan_df)
    all_scan_df$pheno <- current_pheno_name
    all_scan_df$lod_threshold.05 <- lod_threshold_val
    all_results_list[[paste0(current_pheno_name, "_all")]] <- all_scan_df # Store with a different key
    output_file_pheno_all <- file.path(args$outdir, paste0(gsub("[^a-zA-Z0-9_.-]", "_", current_pheno_name), "_", args$qtlmethod, "_all_markers.tsv"))
    write.table(all_scan_df, file = output_file_pheno_all, sep = "\t", quote = FALSE, row.names = FALSE)

  }

  # Final calculations for all markers

  # 1. Calculate effects for all markers (for lodsAll file)
  message("Calculating QTL effects for all markers...")
  # Create a separate effectscan call that doesn't pass method to plotting functions
  effects_all <- effectscan(mapthis, 
                           pheno.col = current_pheno_col,
                           addcovar = covariates_matrix)  # Remove method parameter here

  # Add this debug information for effects_all
  message("Effects structure:")
  message("Length of effects_all$a: ", length(effects_all$a))
  message("Some sample marker names from effects_all$a: ", 
          paste(head(names(effects_all$a)), collapse=", "))
  message("Some sample marker names from scanone_result: ", 
          paste(head(rownames(scanone_result)), collapse=", "))

  # Add this after calculating effects_all and before creating lods_all_with_effects
  # Calculate p-values for all markers (computationally intensive)
  message("Calculating p-values for all markers...")
  perm_matrix <- if(is.vector(perm_result)) matrix(perm_result, ncol=1) else as.matrix(perm_result)

  # Initialize p-values vector for all markers
  all_p_values <- numeric(nrow(scanone_result))
  names(all_p_values) <- rownames(scanone_result)

  # Calculate empirical p-values for each marker
  for(i in 1:nrow(scanone_result)) {
    marker_lod <- scanone_result[i, 3]  # LOD score column
    # Count proportion of permutations with LOD >= marker_lod
    all_p_values[i] <- mean(perm_matrix[,1] >= marker_lod)
    # If p-value is 0, set to lowest possible value based on permutation count
    if(all_p_values[i] == 0) all_p_values[i] <- 1/nrow(perm_matrix)
  }

  # Define lods_all_with_effects by matching chr and pos
  lods_all_with_effects <- data.frame(
    chr = scanone_result[,"chr"],
    pos = scanone_result[,"pos"],
    lod = scanone_result[,3],
    a = NA,  # Initialize with NA
    d = NA,  # Initialize with NA
    pvalue = all_p_values
  )

  # Match effects for all markers (simple direct match)
  for (i in 1:nrow(lods_all_with_effects)) {
    current_chr <- lods_all_with_effects$chr[i]
    current_pos <- lods_all_with_effects$pos[i]
    
    # Direct exact matching - no tolerance needed
    idx <- which(effects_all$chr == current_chr & effects_all$pos == current_pos)
    
    if (length(idx) > 0) {
      lods_all_with_effects$a[i] <- effects_all$a[idx[1]]
      lods_all_with_effects$d[i] <- effects_all$d[idx[1]]
    }
  }

  rownames(lods_all_with_effects) <- rownames(scanone_result)

  # Add the row names as a new column called "marker"
  lods_all_with_effects_modified <- lods_all_with_effects
  lods_all_with_effects_modified$marker <- rownames(lods_all_with_effects)

  # Reorder columns to put marker first
  lods_all_with_effects_modified <- lods_all_with_effects_modified[, c("marker", 
                                    setdiff(colnames(lods_all_with_effects_modified), "marker"))]

  # Define output file path for lodsAll before using it
  safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", current_pheno_name)
  output_file_lodsAll <- file.path(args$outdir, paste0(safe_name, "_", args$qtlmethod, "_", args$permnum, "_lodsAll.txt"))

  # Write to the output file
  write.table(lods_all_with_effects_modified, file=output_file_lodsAll, sep="\t", quote=FALSE, 
              row.names=FALSE)

  # 2. Get top marker per chromosome (for simpleLods) - improved version with effects
  simple_lods <- c()
  for(chr in unique(scanone_result$chr)) {
    # Get all markers for this chromosome
    chr_markers <- scanone_result[scanone_result$chr == chr,]
    
    if(nrow(chr_markers) > 0) {
      # Find max LOD score for this chromosome
      max_idx <- which.max(chr_markers[,3])
      top_marker <- chr_markers[max_idx,]
      
      # Get the marker's chromosome and position
      current_chr <- top_marker$chr
      current_pos <- top_marker$pos
      
      # Direct matching by chromosome and position (exact match)
      idx <- which(effects_all$chr == current_chr & effects_all$pos == current_pos)
      
      if(length(idx) > 0) {
        # Perfect match found
        top_marker$a <- effects_all$a[idx[1]]
        top_marker$d <- effects_all$d[idx[1]]
      } else {
        top_marker$a <- NA
        top_marker$d <- NA
        message("Warning: No exact match found for marker at chr ", current_chr, " pos ", current_pos)
      }
      
      # Add threshold values
      top_marker$threshold <- lod_threshold_val
      
      simple_lods <- rbind(simple_lods, top_marker)
    }
  }

  # Add the row names as a new column called "marker"
  simple_lods_modified <- simple_lods
  simple_lods_modified$marker <- rownames(simple_lods)

  # Reorder columns to put marker first
  simple_lods_modified <- simple_lods_modified[, c("marker", 
                          setdiff(colnames(simple_lods_modified), "marker"))]

  # 3. Calculate confidence intervals (for ci file)
  ci_results <- NULL
  message("Attempting to calculate confidence intervals...")
  if(nrow(significant_qtls) > 0) {
    message(paste("Found", nrow(significant_qtls), "significant QTLs for CI calculation"))
    ci_results <- c()
    for(i in 1:nrow(significant_qtls)) {
      tryCatch({
        chr_id <- significant_qtls[i, "chr"]
        message(paste("Calculating CI for QTL on chromosome", chr_id))
        # The findqtl function may be needed to identify the QTL first
        qtl_pos <- significant_qtls[i, "pos"]
        lodint_result <- lodint(scanone_result, chr=chr_id, drop=1.5, expandtomarkers=TRUE)
        message(paste("CI calculation complete for chr", chr_id))
        # Add a column for gene and threshold
        lodint_result$gene <- current_pheno_name
        lodint_result$threshold <- lod_threshold_val
        ci_results <- rbind(ci_results, lodint_result)
      }, error = function(e) {
        message(paste("Error calculating CI for QTL on chromosome", chr_id, ":", e$message))
      })
    }
  } else {
    # Optional: Calculate CI for top markers even if not significant
    message("No significant QTLs found. Calculating CI for top markers per chromosome instead.")
    ci_results <- c()
    for(chr in unique(scanone_result$chr)) {
      tryCatch({
        # Find the peak marker for this chromosome from simple_lods
        top_marker_for_chr <- simple_lods[simple_lods$chr == chr,]
        if(nrow(top_marker_for_chr) > 0) {
          message(paste("Calculating CI for top marker on chromosome", chr))
          lodint_result <- lodint(scanone_result, chr=chr, drop=1.5, expandtomarkers=TRUE)
          lodint_result$gene <- current_pheno_name
          lodint_result$threshold <- lod_threshold_val
          lodint_result$significant <- "no"
          ci_results <- rbind(ci_results, lodint_result)
        }
      }, error = function(e) {
        message(paste("Error calculating CI for chromosome", chr, ":", e$message))
      })
    }
  }

  # ci file (confidence intervals) - with better error handling
  if(!is.null(ci_results) && length(ci_results) > 0) {
    if(is.data.frame(ci_results) && nrow(ci_results) > 0) {
      # Add marker column if needed
      if(!"marker" %in% colnames(ci_results) && !is.null(rownames(ci_results))) {
        ci_results$marker <- rownames(ci_results)
      }
      
      output_file_ci <- file.path(args$outdir, paste0(safe_name, "_", args$qtlmethod, "_ci.txt"))
      write.table(ci_results, file=output_file_ci, sep="\t", quote=FALSE, row.names=FALSE)
      message(paste("CI file written to", output_file_ci))
    } else {
      message("CI results structure is incomplete - no file written")
    }
  } else {
    message("No CI results generated to write to file")
  }

  # Save files with expected naming convention
  safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", current_pheno_name)

  # simpleLods file (top marker per chromosome)
  output_file_simpleLods <- file.path(args$outdir, paste0(safe_name, "_", args$qtlmethod, "_simpleLods.txt"))
  write.table(simple_lods_modified, file=output_file_simpleLods, sep="\t", quote=FALSE, 
              row.names=FALSE)

  # ci file (confidence intervals)
  if(exists("ci_results") && nrow(ci_results) > 0) {
    output_file_ci <- file.path(args$outdir, paste0(safe_name, "_", args$qtlmethod, "_ci.txt"))
    write.table(ci_results, file=output_file_ci, sep="\t", quote=FALSE, row.names=FALSE)
  }

  # Regular lods file
  output_file_lods <- file.path(args$outdir, paste0(safe_name, "_", args$qtlmethod, "_", args$permnum, "_lods.txt"))
  # Add the row names as a new column called "marker"
  scanone_result_modified <- as.data.frame(scanone_result)
  scanone_result_modified$marker <- rownames(scanone_result)

  # Reorder columns to put marker first
  scanone_result_modified <- scanone_result_modified[, c("marker", 
                             setdiff(colnames(scanone_result_modified), "marker"))]

  # Write without row names (since they're now in the marker column)
  write.table(scanone_result_modified, file=output_file_lods, 
              sep="\t", quote=FALSE, row.names=FALSE)
}

# Combine all significant results (if any)
if (length(all_results_list) > 0) {
  # Filter for data frames that are not the "_all" marker scans if you only want significant ones combined
  significant_results_df_list <- all_results_list[!grepl("_all$", names(all_results_list))]
  if(length(significant_results_df_list) > 0) {
    final_combined_results <- dplyr::bind_rows(significant_results_df_list)
    output_file_all_significant <- file.path(args$outdir, paste0("all_significant_phenotypes_", args$qtlmethod, "_eQTL_summary.tsv"))
    write.table(final_combined_results, file = output_file_all_significant, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("Combined significant results for all phenotypes saved to", output_file_all_significant))
  } else {
    message("No significant QTLs found across any phenotypes to combine.")
  }
} else {
  message("No results generated for any phenotype.")
}

message("eQTL analysis script finished.")

# Debug information about structures
message("Number of rows in scanone_result: ", nrow(scanone_result))
message("Number of rows in effects_all: ", nrow(effects_all))

# Print first few entries of each to compare
message("First few entries of scanone_result:")
print(head(data.frame(chr=scanone_result$chr, pos=scanone_result$pos)))

message("First few entries of effects_all:")
print(head(data.frame(chr=effects_all$chr, pos=effects_all$pos)))