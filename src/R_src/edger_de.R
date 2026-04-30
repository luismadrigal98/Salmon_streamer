#!/usr/bin/env Rscript

################################################################################
# Differential Expression Analysis using edgeR
# SalmonStreamer Pipeline Module
#
# Author: Luis Javier Madrigal-Roca
# Date: 2026-04-30
#
# Usage (called by SalmonStreamer EdgeRDE):
#   Rscript edger_de.R <input_file> <output_dir> <metadata_file> \
#                      [fdr_threshold] [logfc_threshold] [sample_suffix]
#
# Arguments:
#   input_file       Tab-separated count matrix (genes x samples).
#                    First column: gene/transcript IDs.
#                    Remaining columns: one per sample (column header = sample name).
#
#   output_dir       Directory where all output files are written.
#
#   metadata_file    Tab-separated file describing samples.
#
#                    REQUIRED COLUMNS:
#                      sample_name  Identifier matching the count matrix column
#                                   header (after sample_suffix stripping, if any).
#                      group        Experimental group label used to define
#                                   contrasts (e.g. "SpeciesA_Leaf").
#
#                    ALTERNATIVE: provide 'species' + 'tissue' instead of
#                    'group'; the script will derive group as
#                    paste(species, tissue, sep = "_").
#
#                    OPTIONAL COLUMNS (ignored by this script, may be present):
#                      platform, library, flowcell, lane, replicate, ...
#
#   fdr_threshold    FDR significance cutoff (default: 0.05).
#   logfc_threshold  |log2FC| threshold for significance categories (default: 1).
#   sample_suffix    Regex stripped from count matrix column names before
#                    matching to metadata sample_name values.
#                    Example: "_R1_filtered$". Pass "NULL" to skip.
################################################################################

options(bitmapType = "cairo")
if (!interactive()) options(device = "pdf")

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(readr)
})

# ---------------------------------------------------------------------------
# 0. Argument parsing
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: edger_de.R <input_file> <output_dir> <metadata_file> ",
    "[fdr_threshold] [logfc_threshold] [sample_suffix]"
  )
}

input_file      <- args[1]
output_dir      <- args[2]
metadata_arg    <- args[3]
fdr_threshold   <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
logfc_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 1.0
sample_suffix   <- if (length(args) >= 6 && args[6] != "NULL") args[6] else NULL

metadata_file <- if (metadata_arg == "NULL") NULL else metadata_arg

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== SalmonStreamer EdgeRDE Analysis ===\n")
cat("Input file     :", input_file, "\n")
cat("Output dir     :", output_dir, "\n")
cat("Metadata file  :", ifelse(is.null(metadata_file), "(none)", metadata_file), "\n")
cat("FDR threshold  :", fdr_threshold, "\n")
cat("logFC threshold:", logfc_threshold, "\n\n")

# ---------------------------------------------------------------------------
# 1. Load count matrix
# ---------------------------------------------------------------------------

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

counts_data    <- read_tsv(input_file, show_col_types = FALSE)
transcript_ids <- counts_data[[1]]
count_matrix   <- as.matrix(counts_data[, -1])
rownames(count_matrix) <- transcript_ids

cat("Data dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

if (!is.null(sample_suffix)) {
  colnames(count_matrix) <- gsub(sample_suffix, "", colnames(count_matrix))
  cat("Stripped suffix '", sample_suffix, "' from sample names\n", sep = "")
}
sample_names <- colnames(count_matrix)

# ---------------------------------------------------------------------------
# 2. Load and validate metadata
# ---------------------------------------------------------------------------

if (is.null(metadata_file)) {
  stop(
    "A metadata file is required. Provide it via --metadata-file or ",
    "use --group-samples to let the Python wrapper build one automatically."
  )
}
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

metadata <- read_tsv(metadata_file, show_col_types = FALSE)

# Normalise the sample identifier column name
col_lower <- tolower(colnames(metadata))
if ("sample_name" %in% col_lower) {
  colnames(metadata)[col_lower == "sample_name"] <- "Sample"
} else if (!"sample" %in% col_lower) {
  stop("Metadata must contain a 'sample_name' or 'Sample' column.")
} else {
  colnames(metadata)[col_lower == "sample"] <- "Sample"
}

# Normalise / derive the group column
col_lower <- tolower(colnames(metadata))
if ("group" %in% col_lower) {
  colnames(metadata)[col_lower == "group"] <- "Group"
} else if (all(c("species", "tissue") %in% col_lower)) {
  colnames(metadata)[col_lower == "species"] <- "Species"
  colnames(metadata)[col_lower == "tissue"]  <- "Tissue"
  metadata$Group <- paste(metadata$Species, metadata$Tissue, sep = "_")
  cat("Derived Group from species + tissue columns\n")
} else {
  stop(
    "Metadata must contain either a 'group' column or both ",
    "'species' and 'tissue' columns."
  )
}

# Align metadata to count matrix
in_both   <- intersect(metadata$Sample, sample_names)
only_meta <- setdiff(metadata$Sample, sample_names)
only_mat  <- setdiff(sample_names, metadata$Sample)

if (length(only_meta) > 0)
  warning("Samples in metadata but not in count matrix: ",
          paste(only_meta, collapse = ", "))
if (length(only_mat) > 0)
  warning("Samples in count matrix but not in metadata: ",
          paste(only_mat, collapse = ", "))
if (length(in_both) == 0)
  stop("No samples shared between metadata and count matrix. ",
       "Check sample_name values and any --sample-suffix stripping.")

metadata     <- metadata[metadata$Sample %in% in_both, ]
count_matrix <- count_matrix[, metadata$Sample, drop = FALSE]

cat("Samples retained:", nrow(metadata), "\n")
cat("Groups          :", paste(unique(metadata$Group), collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# 3. DGEList, filtering, TMM normalisation
# ---------------------------------------------------------------------------

metadata$group_factor <- factor(metadata$Group)
dge <- DGEList(counts = count_matrix, samples = metadata,
               group = metadata$group_factor)

cat("Genes before filtering:", nrow(dge), "\n")
keep <- filterByExpr(dge, group = metadata$group_factor)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes after filterByExpr:", nrow(dge), "\n")

dge <- calcNormFactors(dge, method = "TMM")
cat("TMM factors:", paste(round(dge$samples$norm.factors, 3), collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# 4. Exploratory plots
# ---------------------------------------------------------------------------

cat("Generating exploratory plots...\n")
lcpm <- cpm(dge, log = TRUE)

# PCA
pca_result <- prcomp(t(lcpm), scale. = FALSE)
pct_var    <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

pca_df <- data.frame(
  PC1    = pca_result$x[, 1],
  PC2    = pca_result$x[, 2],
  Sample = rownames(pca_result$x)
)
pca_df <- merge(pca_df, metadata[, c("Sample", "Group")], by = "Sample")
hulls  <- pca_df %>% group_by(Group) %>% slice(chull(PC1, PC2))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_polygon(data = hulls, alpha = 0.15, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "PCA of Gene Expression",
    x     = paste0("PC1 (", pct_var[1], "%)"),
    y     = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank())

ggsave(file.path(output_dir, "PCA_plot.pdf"), p_pca,
       width = 8, height = 6, device = "pdf")
ggsave(file.path(output_dir, "PCA_plot.png"), p_pca,
       width = 8, height = 6, dpi = 300, device = "png")

# Sample correlation heatmap
cor_matrix <- cor(lcpm, method = "pearson")
ann_df     <- data.frame(Group = metadata$Group, row.names = metadata$Sample)

pdf(file.path(output_dir, "sample_correlation_heatmap.pdf"), width = 10, height = 10)
pheatmap(cor_matrix, annotation_col = ann_df, annotation_row = ann_df,
         main = "Sample Correlation Heatmap")
dev.off()

# Library sizes
lib_df <- data.frame(
  Sample  = metadata$Sample,
  LibSize = dge$samples$lib.size / 1e6,
  Group   = metadata$Group
)
p_lib <- ggplot(lib_df, aes(x = Sample, y = LibSize, fill = Group)) +
  geom_col(alpha = 0.8) +
  labs(title = "Library Sizes", x = "Sample", y = "Library Size (Millions)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank())

ggsave(file.path(output_dir, "library_sizes.pdf"), p_lib,
       width = 10, height = 6, device = "pdf")
ggsave(file.path(output_dir, "library_sizes.png"), p_lib,
       width = 10, height = 6, dpi = 300, device = "png")

# ---------------------------------------------------------------------------
# 5. Statistical analysis — edgeR quasi-likelihood pipeline
# ---------------------------------------------------------------------------

cat("Fitting edgeR quasi-likelihood model...\n")

# ~ 0 + group_factor gives one explicit column per group (no reference level),
# which makes all pairwise contrasts symmetric and equally easy to specify.
design <- model.matrix(~ 0 + group_factor, data = metadata)
colnames(design) <- make.names(levels(metadata$group_factor))

dge <- estimateDisp(dge, design, robust = TRUE)
cat("Common BCV:", round(sqrt(dge$common.dispersion), 4), "\n")

pdf(file.path(output_dir, "dispersion_plot.pdf"), width = 8, height = 6)
plotBCV(dge, main = "Biological Coefficient of Variation")
dev.off()

fit <- glmQLFit(dge, design, robust = TRUE)

pdf(file.path(output_dir, "QL_dispersion_plot.pdf"), width = 8, height = 6)
plotQLDisp(fit, main = "Quasi-likelihood Dispersion")
dev.off()

# All pairwise contrasts
groups      <- colnames(design)
orig_groups <- levels(metadata$group_factor)  # keep original labels for naming
cat("Groups:", paste(orig_groups, collapse = ", "), "\n\n")

results_list <- list()
for (i in seq_len(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1     <- groups[i]
    g2     <- groups[j]
    g1_lab <- orig_groups[i]
    g2_lab <- orig_groups[j]
    cmp    <- paste0(g1_lab, "_vs_", g2_lab)

    contrast_vec       <- setNames(numeric(length(groups)), groups)
    contrast_vec[g1]   <- 1
    contrast_vec[g2]   <- -1

    qlf <- glmQLFTest(fit, contrast = contrast_vec)
    results_list[[cmp]] <- topTags(qlf, n = Inf)$table
    cat("Completed:", cmp, "\n")
  }
}

# ---------------------------------------------------------------------------
# 6. Results export and visualisations
# ---------------------------------------------------------------------------

cat("\nExporting results...\n")
all_sig_genes <- character(0)

for (cmp in names(results_list)) {
  tbl <- results_list[[cmp]]
  tbl$TranscriptID  <- rownames(tbl)
  tbl <- tbl[, c("TranscriptID", setdiff(colnames(tbl), "TranscriptID"))]

  tbl$Significance <- "Not Significant"
  tbl$Significance[tbl$FDR < fdr_threshold &
                     abs(tbl$logFC) > logfc_threshold] <- "Significant"
  tbl$Significance[tbl$FDR < 0.01 &
                     abs(tbl$logFC) > 2 * logfc_threshold] <- "Highly Significant"
  tbl <- tbl[order(tbl$FDR), ]

  write_tsv(tbl, file.path(output_dir, paste0(cmp, "_DE_results.tsv")))

  sig <- tbl[tbl$FDR < fdr_threshold, ]
  write_tsv(sig, file.path(output_dir, paste0(cmp, "_significant_genes.tsv")))
  all_sig_genes <- union(all_sig_genes, sig$TranscriptID)

  cat(sprintf("  %s -> tested: %d | sig: %d | up: %d | down: %d\n",
              cmp, nrow(tbl), nrow(sig),
              sum(sig$logFC > 0), sum(sig$logFC < 0)))

  # Volcano plot
  vdata <- tbl
  vdata$neg_log10_FDR <- -log10(pmax(vdata$FDR, .Machine$double.xmin))

  p_vol <- ggplot(vdata, aes(x = logFC, y = neg_log10_FDR, color = Significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c(
      "Not Significant"    = "grey70",
      "Significant"        = "#377EB8",
      "Highly Significant" = "#E41A1C"
    )) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(fdr_threshold),
               linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("Volcano:", cmp),
      x     = "log2 Fold Change",
      y     = paste0("-log10(FDR)")
    ) +
    theme_bw() +
    theme(legend.title = element_blank(), panel.grid = element_blank())

  ggsave(file.path(output_dir, paste0(cmp, "_volcano.pdf")), p_vol,
         width = 10, height = 8, device = "pdf")
  ggsave(file.path(output_dir, paste0(cmp, "_volcano.png")), p_vol,
         width = 10, height = 8, dpi = 300, device = "png")
}

# Expression heatmap of top DE genes
if (length(all_sig_genes) > 0) {
  genes_to_plot <- if (length(all_sig_genes) > 100) {
    unique(unlist(lapply(results_list, function(x)
      head(rownames(x)[x$FDR < fdr_threshold], 50))))
  } else {
    all_sig_genes
  }
  genes_to_plot <- intersect(genes_to_plot, rownames(lcpm))

  if (length(genes_to_plot) > 0) {
    hz <- t(scale(t(lcpm[genes_to_plot, , drop = FALSE])))
    pdf(file.path(output_dir, "DE_genes_heatmap.pdf"),
        width = 12, height = max(8, length(genes_to_plot) * 0.15))
    pheatmap(hz, annotation_col = ann_df,
             show_rownames = length(genes_to_plot) <= 50,
             main = "Top Differentially Expressed Genes",
             fontsize_row = 6)
    dev.off()
  }
}

# ---------------------------------------------------------------------------
# 7. Summary report
# ---------------------------------------------------------------------------

summary_lines <- c(
  "=== SalmonStreamer EdgeRDE Analysis Summary ===",
  paste("Date         :", Sys.Date()),
  paste("Input file   :", input_file),
  paste("Metadata     :", metadata_file),
  paste("Genes (post-filter):", nrow(dge)),
  paste("Samples      :", ncol(dge)),
  paste("Groups       :", paste(orig_groups, collapse = ", ")),
  paste("FDR threshold:", fdr_threshold),
  paste("logFC threshold:", logfc_threshold),
  "",
  "Comparisons:"
)
for (cmp in names(results_list)) {
  tbl <- results_list[[cmp]]
  sig <- sum(tbl$FDR < fdr_threshold, na.rm = TRUE)
  summary_lines <- c(summary_lines,
    sprintf("  %-60s %d significant genes", cmp, sig))
}

writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

cat("\nDone. Results written to:", output_dir, "\n")
