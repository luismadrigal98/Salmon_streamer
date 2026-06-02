#!/usr/bin/env Rscript

################################################################################
# Differential Expression from a pre-normalized expression matrix
# SalmonStreamer Pipeline Module (limma-trend pipeline)
#
# Author: Luis Javier Madrigal-Roca
# Date:   2026-06-02
#
# Use this when the input matrix has already been TMM-normalized (or otherwise
# normalized to a sample-comparable scale) and you need to run DE on the
# normalized values — for example, after collapsing paralogs by summing
# TMM-normalized CPMs.
#
# The classic edgeR quasi-likelihood pipeline cannot be reused here because it
# requires raw integer counts. Instead we apply the limma-trend workflow:
#
#   logCPM -> lmFit(design) -> contrasts.fit(contrast) -> eBayes(trend=TRUE)
#
# This is the workflow limma's user guide recommends for already-normalized
# expression data (Section "When to use limma-trend").
#
# Usage (called by SalmonStreamer EdgeRDEFromNormalized):
#   Rscript edger_de_from_normalized.R <input_file> <output_dir> <metadata_file> \
#                                      [fdr_threshold] [logfc_threshold] \
#                                      [sample_suffix] [input_scale] [prior_count]
#
# Arguments:
#   input_file       Tab-separated normalized expression matrix
#                    (genes/transcripts x samples). First column = TranscriptID.
#
#   output_dir       Directory where all output files are written.
#
#   metadata_file    Same format as edger_de.R: requires `sample_name` and
#                    `group` columns (or `species` + `tissue`).
#
#   fdr_threshold    FDR significance cutoff (default: 0.05).
#   logfc_threshold  |log2FC| threshold for significance categories (default: 1).
#   sample_suffix    Regex stripped from column names before matching to
#                    metadata sample_name values. Pass "NULL" to skip.
#   input_scale      "logcpm" (default) or "cpm".
#                    - logcpm: values are used as-is for lmFit.
#                    - cpm   : values are converted via log2(x + prior_count)
#                              before fitting.
#   prior_count      Pseudo-count added before log2 when input_scale == "cpm"
#                    (default: 0.25).
################################################################################

options(bitmapType = "cairo")
if (!interactive()) options(device = "pdf")

suppressPackageStartupMessages({
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
    "Usage: edger_de_from_normalized.R <input_file> <output_dir> <metadata_file> ",
    "[fdr_threshold] [logfc_threshold] [sample_suffix] [input_scale] [prior_count]"
  )
}

input_file      <- args[1]
output_dir      <- args[2]
metadata_arg    <- args[3]
fdr_threshold   <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
logfc_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 1.0
sample_suffix   <- if (length(args) >= 6 && args[6] != "NULL") args[6] else NULL
input_scale     <- if (length(args) >= 7) tolower(args[7]) else "logcpm"
prior_count     <- if (length(args) >= 8) as.numeric(args[8]) else 0.25

if (!input_scale %in% c("logcpm", "cpm")) {
  stop("input_scale must be 'logcpm' or 'cpm' (got '", input_scale, "')")
}

metadata_file <- if (metadata_arg == "NULL") NULL else metadata_arg

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== SalmonStreamer EdgeRDEFromNormalized Analysis ===\n")
cat("Input file     :", input_file, "\n")
cat("Output dir     :", output_dir, "\n")
cat("Metadata file  :", ifelse(is.null(metadata_file), "(none)", metadata_file), "\n")
cat("Input scale    :", input_scale, "\n")
if (input_scale == "cpm") {
  cat("Prior count    :", prior_count, "\n")
}
cat("FDR threshold  :", fdr_threshold, "\n")
cat("logFC threshold:", logfc_threshold, "\n\n")

# ---------------------------------------------------------------------------
# 1. Load expression matrix
# ---------------------------------------------------------------------------

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

expr_data      <- read_tsv(input_file, show_col_types = FALSE)
transcript_ids <- expr_data[[1]]
expr_matrix    <- as.matrix(expr_data[, -1])
rownames(expr_matrix) <- transcript_ids
storage.mode(expr_matrix) <- "double"

cat("Data dimensions:", nrow(expr_matrix), "genes x",
    ncol(expr_matrix), "samples\n")

if (!is.null(sample_suffix)) {
  colnames(expr_matrix) <- gsub(sample_suffix, "", colnames(expr_matrix))
  cat("Stripped suffix '", sample_suffix, "' from sample names\n", sep = "")
}
sample_names <- colnames(expr_matrix)

# Convert raw CPM to log2-CPM with prior
if (input_scale == "cpm") {
  if (any(expr_matrix < 0, na.rm = TRUE)) {
    stop(
      "Input was declared as raw CPM (input_scale='cpm') but contains ",
      "negative values. Did you mean input_scale='logcpm'?"
    )
  }
  cat("Converting raw CPM -> log2(CPM + ", prior_count, ")\n", sep = "")
  expr_matrix <- log2(expr_matrix + prior_count)
}

logcpm <- expr_matrix  # from here on, everything assumes log-scale

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

col_lower <- tolower(colnames(metadata))
if ("sample_name" %in% col_lower) {
  colnames(metadata)[col_lower == "sample_name"] <- "Sample"
} else if (!"sample" %in% col_lower) {
  stop("Metadata must contain a 'sample_name' or 'Sample' column.")
} else {
  colnames(metadata)[col_lower == "sample"] <- "Sample"
}

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

in_both   <- intersect(metadata$Sample, sample_names)
only_meta <- setdiff(metadata$Sample, sample_names)
only_mat  <- setdiff(sample_names, metadata$Sample)

if (length(only_meta) > 0)
  warning("Samples in metadata but not in expression matrix: ",
          paste(only_meta, collapse = ", "))
if (length(only_mat) > 0)
  warning("Samples in expression matrix but not in metadata: ",
          paste(only_mat, collapse = ", "))
if (length(in_both) == 0)
  stop("No samples shared between metadata and expression matrix. ",
       "Check sample_name values and any --sample-suffix stripping.")

metadata <- metadata[metadata$Sample %in% in_both, ]
logcpm   <- logcpm[, metadata$Sample, drop = FALSE]

cat("Samples retained:", nrow(metadata), "\n")
cat("Groups          :", paste(unique(metadata$Group), collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# 3. Design matrix
# ---------------------------------------------------------------------------

metadata$group_factor <- factor(metadata$Group)
design <- model.matrix(~ 0 + group_factor, data = metadata)
colnames(design) <- make.names(levels(metadata$group_factor))

# ---------------------------------------------------------------------------
# 4. Exploratory plots (PCA + sample correlation)
# ---------------------------------------------------------------------------

cat("Generating exploratory plots...\n")

# Drop rows that are constant or all-NA before PCA / correlation
row_sd <- apply(logcpm, 1, sd, na.rm = TRUE)
keep_var <- is.finite(row_sd) & row_sd > 0
if (sum(keep_var) < nrow(logcpm)) {
  cat("Dropping", nrow(logcpm) - sum(keep_var),
      "constant / all-NA rows from PCA & correlation plots\n")
}
logcpm_var <- logcpm[keep_var, , drop = FALSE]

pca_result <- prcomp(t(logcpm_var), scale. = FALSE)
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
    title = "PCA of Normalized Expression",
    x     = paste0("PC1 (", pct_var[1], "%)"),
    y     = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank())

ggsave(file.path(output_dir, "PCA_plot.pdf"), p_pca,
       width = 8, height = 6, device = "pdf")
ggsave(file.path(output_dir, "PCA_plot.png"), p_pca,
       width = 8, height = 6, dpi = 300, device = "png")

cor_matrix <- cor(logcpm_var, method = "pearson")
ann_df     <- data.frame(Group = metadata$Group, row.names = metadata$Sample)

pdf(file.path(output_dir, "sample_correlation_heatmap.pdf"),
    width = 10, height = 10)
pheatmap(cor_matrix, annotation_col = ann_df, annotation_row = ann_df,
         main = "Sample Correlation Heatmap")
dev.off()

# ---------------------------------------------------------------------------
# 5. limma-trend fit
# ---------------------------------------------------------------------------

cat("Fitting limma-trend model on log-CPM...\n")
fit_full <- lmFit(logcpm, design)

# ---------------------------------------------------------------------------
# 6. All pairwise contrasts + per-contrast topTable
# ---------------------------------------------------------------------------

groups      <- colnames(design)
orig_groups <- levels(metadata$group_factor)
cat("Groups:", paste(orig_groups, collapse = ", "), "\n\n")

results_list  <- list()
all_sig_genes <- character(0)

for (i in seq_len(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1     <- groups[i]
    g2     <- groups[j]
    g1_lab <- orig_groups[i]
    g2_lab <- orig_groups[j]
    cmp    <- paste0(g1_lab, "_vs_", g2_lab)

    contrast_vec     <- setNames(numeric(length(groups)), groups)
    contrast_vec[g1] <- 1
    contrast_vec[g2] <- -1

    fit_c <- contrasts.fit(fit_full, contrasts = matrix(contrast_vec, ncol = 1))
    fit_c <- eBayes(fit_c, trend = TRUE)
    tt    <- topTable(fit_c, number = Inf, sort.by = "P")

    # Normalise output column names to match the EdgeRDE schema where possible
    if ("adj.P.Val" %in% colnames(tt)) {
      colnames(tt)[colnames(tt) == "adj.P.Val"] <- "FDR"
    }
    if ("P.Value" %in% colnames(tt)) {
      colnames(tt)[colnames(tt) == "P.Value"] <- "PValue"
    }
    if ("AveExpr" %in% colnames(tt)) {
      colnames(tt)[colnames(tt) == "AveExpr"] <- "logCPM"
    }

    tt$TranscriptID <- rownames(tt)
    tt <- tt[, c("TranscriptID", setdiff(colnames(tt), "TranscriptID"))]
    tt$Significance <- "Not Significant"
    tt$Significance[tt$FDR < fdr_threshold &
                      abs(tt$logFC) > logfc_threshold] <- "Significant"
    tt$Significance[tt$FDR < 0.01 &
                      abs(tt$logFC) > 2 * logfc_threshold] <- "Highly Significant"
    tt <- tt[order(tt$FDR), ]

    results_list[[cmp]] <- tt

    write_tsv(tt, file.path(output_dir, paste0(cmp, "_DE_results.tsv")))

    sig <- tt[tt$FDR < fdr_threshold, ]
    write_tsv(sig, file.path(output_dir, paste0(cmp, "_significant_genes.tsv")))
    all_sig_genes <- union(all_sig_genes, sig$TranscriptID)

    cat(sprintf("  %s -> tested: %d | sig: %d | up: %d | down: %d\n",
                cmp, nrow(tt), nrow(sig),
                sum(sig$logFC > 0), sum(sig$logFC < 0)))

    # Volcano plot
    vdata <- tt
    vdata$neg_log10_FDR <- -log10(pmax(vdata$FDR, .Machine$double.xmin))

    p_vol <- ggplot(vdata, aes(x = logFC, y = neg_log10_FDR,
                               color = Significance)) +
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
        y     = "-log10(FDR)"
      ) +
      theme_bw() +
      theme(legend.title = element_blank(), panel.grid = element_blank())

    ggsave(file.path(output_dir, paste0(cmp, "_volcano.pdf")), p_vol,
           width = 10, height = 8, device = "pdf")
    ggsave(file.path(output_dir, paste0(cmp, "_volcano.png")), p_vol,
           width = 10, height = 8, dpi = 300, device = "png")
  }
}

# Heatmap of top DE genes
if (length(all_sig_genes) > 0) {
  genes_to_plot <- if (length(all_sig_genes) > 100) {
    unique(unlist(lapply(results_list, function(x)
      head(x$TranscriptID[x$FDR < fdr_threshold], 50))))
  } else {
    all_sig_genes
  }
  genes_to_plot <- intersect(genes_to_plot, rownames(logcpm))

  if (length(genes_to_plot) > 0) {
    hz <- t(scale(t(logcpm[genes_to_plot, , drop = FALSE])))
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
  "=== SalmonStreamer EdgeRDEFromNormalized Analysis Summary ===",
  paste("Date            :", Sys.Date()),
  paste("Input file      :", input_file),
  paste("Metadata        :", metadata_file),
  paste("Input scale     :", input_scale),
  if (input_scale == "cpm") paste("Prior count     :", prior_count) else NULL,
  paste("Genes (total)   :", nrow(logcpm)),
  paste("Samples         :", ncol(logcpm)),
  paste("Groups          :", paste(orig_groups, collapse = ", ")),
  paste("FDR threshold   :", fdr_threshold),
  paste("logFC threshold :", logfc_threshold),
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
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))

cat("\nDone. Results written to:", output_dir, "\n")
