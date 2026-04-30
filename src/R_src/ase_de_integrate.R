#!/usr/bin/env Rscript

################################################################################
# ASE + Differential Expression Integration
# SalmonStreamer Pipeline Module
#
# Author: Luis Javier Madrigal-Roca
# Date: 2026-04-30
#
# Usage (called by SalmonStreamer ASEIntegrate):
#   Rscript ase_de_integrate.R \
#       <de_results_dir> <ase_counts_file> <output_dir> \
#       [fdr_threshold] \
#       [expression_file] [metadata_file] \
#       [hybrid_group] [parent1_group] [parent2_group] \
#       [parent1_label] [parent2_label]
#
# Arguments:
#   de_results_dir   Directory containing *_DE_results.tsv files produced by
#                    the EdgeRDE subcommand.  All matching files are processed.
#
#   ase_counts_file  Tab-separated file with per-gene ASE read counts.
#
#                    REQUIRED COLUMNS:
#                      gene_id          Gene identifier (matched against
#                                       TranscriptID in the DE results after
#                                       stripping any ".T[digit]" suffix).
#                      snp_count        Number of diagnostic SNPs in the gene.
#                      bias_ratio       Fraction of SNPs showing allelic bias.
#                      predominant_bias Label of the allele with majority reads
#                                       (matches parent1_label or parent2_label).
#                      <parent1_label>_ratio  Proportion of reads from parent 1.
#                      <parent2_label>_ratio  Proportion of reads from parent 2.
#
#                    This format is produced by ase_read_counter.py from
#                    the EDGE_Penstemon pipeline.
#
#   output_dir       Directory where all output files are written.
#
#   fdr_threshold    FDR significance cutoff (default: 0.05).
#
#   expression_file  (Optional) Count matrix TSV used to compute the Expression
#                    Deviation (Ed) metric.  Required for the Ad/Ed regulatory
#                    classification plot.  Pass "NULL" to skip.
#
#   metadata_file    (Optional) Metadata TSV with sample_name and group columns
#                    (same format as required by EdgeRDE).  Required when
#                    expression_file is provided.  Pass "NULL" to skip.
#
#   hybrid_group     Group label in metadata identifying hybrid/F1 samples
#                    (e.g. "F1_hybrid").  Required when expression_file is
#                    provided.  Pass "NULL" to skip.
#
#   parent1_group    Group label for parent 1 samples.  Pass "NULL" to skip.
#   parent2_group    Group label for parent 2 samples.  Pass "NULL" to skip.
#
#   parent1_label    Value in the predominant_bias column representing parent 1
#                    (default: "parent1").
#   parent2_label    Value in the predominant_bias column representing parent 2
#                    (default: "parent2").
################################################################################

options(bitmapType = "cairo")
if (!interactive()) options(device = "pdf")

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# ---------------------------------------------------------------------------
# 0. Argument parsing
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: ase_de_integrate.R <de_results_dir> <ase_counts_file> <output_dir> ",
    "[fdr_threshold] [expression_file] [metadata_file] ",
    "[hybrid_group] [parent1_group] [parent2_group] ",
    "[parent1_label] [parent2_label]"
  )
}

parse_arg <- function(a, default = NULL) {
  if (is.na(a) || a == "NULL") default else a
}

de_results_dir  <- args[1]
ase_counts_file <- args[2]
output_dir      <- args[3]
fdr_threshold   <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
expression_file <- parse_arg(args[5])
metadata_file   <- parse_arg(args[6])
hybrid_group    <- parse_arg(args[7])
parent1_group   <- parse_arg(args[8])
parent2_group   <- parse_arg(args[9])
parent1_label   <- if (length(args) >= 10 && !is.na(args[10]) && args[10] != "NULL") args[10] else "parent1"
parent2_label   <- if (length(args) >= 11 && !is.na(args[11]) && args[11] != "NULL") args[11] else "parent2"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== SalmonStreamer ASEIntegrate Analysis ===\n")
cat("DE results dir  :", de_results_dir, "\n")
cat("ASE counts file :", ase_counts_file, "\n")
cat("Output dir      :", output_dir, "\n")
cat("FDR threshold   :", fdr_threshold, "\n")
cat("Expression file :", ifelse(is.null(expression_file), "(skipped)", expression_file), "\n")
cat("Parent 1 label  :", parent1_label, "\n")
cat("Parent 2 label  :", parent2_label, "\n\n")

# ---------------------------------------------------------------------------
# 1. Load DE results
# ---------------------------------------------------------------------------

de_files <- list.files(de_results_dir, pattern = "_DE_results\\.tsv$",
                       full.names = TRUE)
if (length(de_files) == 0)
  stop("No *_DE_results.tsv files found in: ", de_results_dir)

cat("Found", length(de_files), "DE results file(s):\n")
for (f in de_files) cat("  ", basename(f), "\n")
cat("\n")

# ---------------------------------------------------------------------------
# 2. Load ASE counts
# ---------------------------------------------------------------------------

if (!file.exists(ase_counts_file))
  stop("ASE counts file not found: ", ase_counts_file)

ase_data <- read_tsv(ase_counts_file, show_col_types = FALSE)

required_ase_cols <- c("gene_id", "snp_count", "bias_ratio", "predominant_bias")
missing_cols <- setdiff(required_ase_cols, colnames(ase_data))
if (length(missing_cols) > 0)
  stop("ASE counts file is missing columns: ", paste(missing_cols, collapse = ", "))

p1_ratio_col <- paste0(parent1_label, "_ratio")
p2_ratio_col <- paste0(parent2_label, "_ratio")
if (!p1_ratio_col %in% colnames(ase_data))
  warning("Column '", p1_ratio_col, "' not found in ASE file. ",
          "Check --parent1-label.")

cat("Loaded ASE data for", nrow(ase_data), "genes\n\n")

# ---------------------------------------------------------------------------
# Helper: classify gene regulatory category from DE + ASE
# ---------------------------------------------------------------------------

categorise_gene <- function(is_de, is_biased, bias_toward_p1, logfc) {
  # positive logFC = higher expression in parent 1 (by EdgeRDE convention)
  if (is_de && is_biased) {
    if ((bias_toward_p1 && logfc > 0) || (!bias_toward_p1 && logfc < 0)) {
      return("cis")             # ASE bias and DE direction agree
    } else {
      return("trans_compensatory")  # ASE bias and DE direction oppose
    }
  } else if (is_de && !is_biased) {
    return("trans_only")
  } else if (!is_de && is_biased) {
    return("compensatory")
  } else {
    return("no_divergence")
  }
}

# ---------------------------------------------------------------------------
# 3. Simple DE + ASE integration for each comparison
# ---------------------------------------------------------------------------

all_cis_candidates <- list()

for (de_file in de_files) {
  cmp <- sub("_DE_results\\.tsv$", "", basename(de_file))
  cat("Processing:", cmp, "\n")

  de_data <- read_tsv(de_file, show_col_types = FALSE)

  # Strip .T[digit] transcript suffix so gene IDs match ASE gene_id
  de_data$gene_id_base <- sub("\\.T[0-9]+$", "", de_data$TranscriptID)

  merged <- de_data %>%
    left_join(ase_data, by = c("gene_id_base" = "gene_id")) %>%
    filter(!is.na(snp_count))

  if (nrow(merged) == 0) {
    cat("  No genes with ASE data found for this comparison.\n\n")
    next
  }

  cat("  Genes with both DE and ASE data:", nrow(merged), "\n")

  merged <- merged %>%
    mutate(
      is_DE      = FDR < fdr_threshold,
      is_biased  = bias_ratio > 0.5,
      bias_p1    = predominant_bias == parent1_label,
      ase_de_category = mapply(
        categorise_gene, is_DE, is_biased, bias_p1, logFC
      )
    )

  # Export combined results
  write_tsv(merged,
            file.path(output_dir, paste0(cmp, "_with_ASE.tsv")))

  # Category summary
  summary_tbl <- merged %>%
    group_by(ase_de_category) %>%
    summarise(
      n_genes           = n(),
      avg_logFC         = mean(logFC, na.rm = TRUE),
      avg_bias_ratio    = mean(bias_ratio, na.rm = TRUE),
      .groups = "drop"
    )
  write_tsv(summary_tbl,
            file.path(output_dir, paste0(cmp, "_ASE_category_summary.tsv")))

  cat("  Category breakdown:\n")
  print(as.data.frame(summary_tbl))

  # Scatter plot: logFC vs ASE bias ratio, coloured by category
  if (p1_ratio_col %in% colnames(merged)) {
    p_scatter <- ggplot(
      merged,
      aes(x = logFC, y = .data[[p1_ratio_col]], color = ase_de_category)
    ) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
      geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey50") +
      scale_color_brewer(palette = "Set2", name = "Category") +
      labs(
        title   = paste("ASE bias vs Differential Expression:", cmp),
        x       = "log2 Fold Change (DE)",
        y       = paste0(parent1_label, " allele ratio (ASE)"),
        caption = "Upper-right and lower-left quadrants indicate cis-regulatory divergence"
      ) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "bottom")

    ggsave(file.path(output_dir, paste0(cmp, "_ASE_vs_DE.pdf")),
           p_scatter, width = 10, height = 7, device = "pdf")
    ggsave(file.path(output_dir, paste0(cmp, "_ASE_vs_DE.png")),
           p_scatter, width = 10, height = 7, dpi = 300, device = "png")
  }

  cis <- merged %>%
    filter(ase_de_category == "cis") %>%
    arrange(FDR)
  if (nrow(cis) > 0)
    all_cis_candidates[[cmp]] <- cis

  cat("  Cis-regulatory candidates:", nrow(cis), "\n\n")
}

# Collated cis candidates across all comparisons
if (length(all_cis_candidates) > 0) {
  cis_all <- bind_rows(all_cis_candidates, .id = "comparison") %>%
    arrange(comparison, FDR)
  write_tsv(cis_all,
            file.path(output_dir, "cis_regulatory_candidates_all.tsv"))
  cat("Saved", nrow(cis_all), "total cis-regulatory candidates across all comparisons\n\n")
}

# ---------------------------------------------------------------------------
# 4. Advanced Ad / Ed analysis (requires expression matrix + metadata)
# ---------------------------------------------------------------------------

run_advanced <- !is.null(expression_file) &&
                !is.null(metadata_file) &&
                !is.null(hybrid_group) &&
                !is.null(parent1_group) &&
                !is.null(parent2_group)

if (run_advanced) {
  cat("--- Advanced Ad/Ed Regulatory Classification ---\n")

  if (!file.exists(expression_file))
    stop("Expression file not found: ", expression_file)
  if (!file.exists(metadata_file))
    stop("Metadata file not found: ", metadata_file)

  # Load and normalise expression
  expr_raw    <- read_tsv(expression_file, show_col_types = FALSE)
  gene_ids    <- expr_raw[[1]]
  expr_matrix <- as.matrix(expr_raw[, -1])
  rownames(expr_matrix) <- gene_ids

  meta <- read_tsv(metadata_file, show_col_types = FALSE)
  col_lower <- tolower(colnames(meta))
  if ("sample_name" %in% col_lower)
    colnames(meta)[col_lower == "sample_name"] <- "Sample"
  if ("group" %in% col_lower)
    colnames(meta)[col_lower == "group"] <- "Group"

  # Retain samples present in the expression matrix
  meta <- meta[meta$Sample %in% colnames(expr_matrix), ]
  expr_matrix <- expr_matrix[, meta$Sample, drop = FALSE]

  # TMM normalisation
  dge_ase <- DGEList(counts = expr_matrix)
  dge_ase <- calcNormFactors(dge_ase, method = "TMM")
  norm_cpm <- cpm(dge_ase, normalized.lib.sizes = TRUE)

  # Per-group averages
  p1_samples <- meta$Sample[meta$Group == parent1_group]
  p2_samples <- meta$Sample[meta$Group == parent2_group]
  hy_samples <- meta$Sample[meta$Group == hybrid_group]

  if (length(p1_samples) == 0) stop("No samples found for parent1_group: ", parent1_group)
  if (length(p2_samples) == 0) stop("No samples found for parent2_group: ", parent2_group)
  if (length(hy_samples) == 0) stop("No samples found for hybrid_group: ",  hybrid_group)

  expr_summary <- data.frame(
    gene_id    = rownames(norm_cpm),
    parent1_avg = rowMeans(norm_cpm[, p1_samples, drop = FALSE]),
    parent2_avg = rowMeans(norm_cpm[, p2_samples, drop = FALSE]),
    hybrid_avg  = rowMeans(norm_cpm[, hy_samples, drop = FALSE]),
    stringsAsFactors = FALSE
  )

  # Merge with ASE data
  combined_ad_ed <- inner_join(ase_data, expr_summary, by = "gene_id") %>%
    mutate(
      # Allele Deviation: deviation from 0.5 expected allele ratio
      # Negative = parent1 bias, positive = parent2 bias
      Ad = 0.5 - .data[[p1_ratio_col]],
      # Parental expression difference
      parental_diff = parent1_avg - parent2_avg,
      # F1 position on the parental spectrum (1 = parent1-like, 0 = parent2-like)
      position_ratio = ifelse(
        parental_diff != 0,
        (hybrid_avg - parent2_avg) / parental_diff,
        NA_real_
      ),
      # Expression Deviation: deviation from midparent expectation
      Ed = 0.5 - position_ratio
    ) %>%
    filter(!is.na(Ad), !is.na(Ed))

  # Classify regulatory mechanism (Wittkopp et al. framework)
  classify_regulatory <- function(ad, ed, cis_thresh = 0.1, trans_thresh = 0.1) {
    has_cis   <- abs(ad) > cis_thresh
    has_trans <- abs(ed) > trans_thresh
    if (has_cis && has_trans) {
      if (sign(ad) == sign(ed)) "Cis + Trans (Reinforcing)" else "Cis x Trans (Opposing)"
    } else if (!has_cis && has_trans) {
      "Trans-only"
    } else {
      "Ambiguous"
    }
  }

  combined_ad_ed$Regulation <- mapply(classify_regulatory,
                                      combined_ad_ed$Ad,
                                      combined_ad_ed$Ed)

  cat("Regulatory classification summary:\n")
  print(table(combined_ad_ed$Regulation))
  cat("\n")

  write_tsv(combined_ad_ed,
            file.path(output_dir, "ase_vs_de_ad_ed_data.tsv"))

  # Regression plot per regulatory class
  plot_data <- combined_ad_ed %>%
    filter(Ed >= -1.5, Ed <= 1.5, Regulation != "Ambiguous")

  if (nrow(plot_data) > 0) {
    group_stats <- plot_data %>%
      group_by(Regulation) %>%
      summarise(
        n      = n(),
        r      = tryCatch(cor.test(Ad, Ed)$estimate, error = function(e) NA_real_),
        p_val  = tryCatch(cor.test(Ad, Ed)$p.value,  error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      mutate(label = sprintf("%s\n(n=%d, r=%.2f, p=%.3g)", Regulation, n, r, p_val))

    plot_data <- left_join(plot_data, group_stats, by = "Regulation")

    p_reg <- ggplot(plot_data,
                    aes(x = Ad, y = Ed, color = label, fill = label)) +
      geom_point(alpha = 0.55, size = 1.8) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "black") +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "black") +
      scale_color_manual(
        values = setNames(
          c("#E41A1C", "#377EB8", "grey50"),
          group_stats$label[match(
            c("Cis + Trans (Reinforcing)", "Cis x Trans (Opposing)", "Trans-only"),
            group_stats$Regulation
          )]
        )
      ) +
      scale_fill_manual(
        values = setNames(
          c("#E41A1C", "#377EB8", "grey50"),
          group_stats$label[match(
            c("Cis + Trans (Reinforcing)", "Cis x Trans (Opposing)", "Trans-only"),
            group_stats$Regulation
          )]
        )
      ) +
      labs(
        title    = "Regulatory Mechanisms: Ad vs Ed",
        subtitle = "Allelic Deviation (cis effect) vs Expression Deviation (total effect)",
        x        = paste0("Allelic Deviation (Ad)\n",
                          "(negative = ", parent1_label, " bias)"),
        y        = "Expression Deviation (Ed)",
        color    = "Regulatory Class",
        fill     = "Regulatory Class"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.text     = element_text(size = 8),
        plot.title      = element_text(face = "bold"),
        axis.title      = element_text(face = "bold")
      )

    ggsave(file.path(output_dir, "regulatory_divergence_regressions.pdf"),
           p_reg, width = 10, height = 7, device = "pdf")
    ggsave(file.path(output_dir, "regulatory_divergence_regressions.png"),
           p_reg, width = 10, height = 7, dpi = 300, device = "png")

    cat("Regulatory divergence plot saved.\n")
  }

} else {
  if (!is.null(expression_file) || !is.null(hybrid_group)) {
    cat("Skipping Ad/Ed analysis: expression_file, metadata_file, hybrid_group,\n")
    cat("parent1_group, and parent2_group must all be provided together.\n")
  }
}

# ---------------------------------------------------------------------------
# 5. Session info
# ---------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
cat("\nDone. Results written to:", output_dir, "\n")
