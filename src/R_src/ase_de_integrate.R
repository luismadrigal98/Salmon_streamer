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
#   ase_counts_file  Tab-separated per-gene ASE table, in either of two formats
#                    (detected from the columns; '#' comment lines are skipped):
#
#                    (a) An ASErtain <prefix>.gene_ase.tsv. A gene is "biased"
#                        when ASErtain called it (ase_call), the direction is
#                        ASErtain's, and the allele fraction is its
#                        mean_<label>_ratio. This is the statistically grounded
#                        input and the one to prefer.
#
#                    (b) The legacy summary table below. A gene is "biased" when
#                        bias_ratio > 0.5 — a majority vote over SNPs, NOT a
#                        test, with no account of depth or replication.
#
#                    LEGACY COLUMNS:
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
#                    These two also say WHICH DE comparison is the parental one
#                    and which way round it is. Cis/trans classification is only
#                    defined for parent 1 versus parent 2, and it needs the sign
#                    of the fold change to mean "higher in parent 1". EdgeRDE
#                    names a comparison <A>_vs_<B> with logFC = A - B, A and B
#                    in alphabetical order, so the sign depends on the group
#                    names. Without the two groups the parental comparison is
#                    inferred from parent1_label / parent2_label appearing in
#                    the group names; every other comparison is skipped.
#
#   parent1_label    Value in the predominant_bias column representing parent 1
#                    (default: "parent1").
#   parent2_label    Value in the predominant_bias column representing parent 2
#                    (default: "parent2").
################################################################################

options(bitmapType = "cairo")
if (!interactive()) options(device = "pdf")

suppressPackageStartupMessages({
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

ase_raw <- read_tsv(ase_counts_file, show_col_types = FALSE, comment = "#")

p1_ratio_col <- paste0(parent1_label, "_ratio")
p2_ratio_col <- paste0(parent2_label, "_ratio")

as_flag <- function(x) as.character(x) %in% c("TRUE", "True", "true", "T", "1")

asertain_cols <- c("gene_id", "ase_call", "direction", "n_snps")
if (all(asertain_cols %in% colnames(ase_raw))) {
  # ASErtain reports the fraction for its "variable" lineage only; work out
  # which of the two parents that is from the column that exists.
  m1 <- paste0("mean_", parent1_label, "_ratio")
  m2 <- paste0("mean_", parent2_label, "_ratio")
  if (m1 %in% colnames(ase_raw)) {
    p1_ratio <- as.numeric(ase_raw[[m1]])
  } else if (m2 %in% colnames(ase_raw)) {
    p1_ratio <- 1 - as.numeric(ase_raw[[m2]])
  } else {
    stop("ASErtain table has neither '", m1, "' nor '", m2, "'. ",
         "--parent1-label / --parent2-label must be the lineage labels used in ",
         "the ASErtain config. Columns found: ",
         paste(colnames(ase_raw), collapse = ", "))
  }
  ase_data <- tibble(
    gene_id          = ase_raw$gene_id,
    snp_count        = ase_raw$n_snps,
    bias_ratio       = NA_real_,
    predominant_bias = as.character(ase_raw$direction),
    is_biased_input  = as_flag(ase_raw$ase_call),
    ase_source       = "asertain"
  )
  ase_data[[p1_ratio_col]] <- p1_ratio
  ase_data[[p2_ratio_col]] <- 1 - p1_ratio
  for (extra in c("q_value", "log2_ratio", "phase_discordance_cause"))
    if (extra %in% colnames(ase_raw))
      ase_data[[paste0("ase_", extra)]] <- ase_raw[[extra]]
  cat("ASE input: ASErtain gene table (bias = ASErtain's ase_call)\n")
} else {
  required_ase_cols <- c("gene_id", "snp_count", "bias_ratio", "predominant_bias")
  missing_cols <- setdiff(required_ase_cols, colnames(ase_raw))
  if (length(missing_cols) > 0)
    stop("ASE counts file is missing columns: ", paste(missing_cols, collapse = ", "))
  if (!p1_ratio_col %in% colnames(ase_raw))
    stop("Column '", p1_ratio_col, "' not found in ASE file. Check --parent1-label.")
  ase_data <- ase_raw
  if (!p2_ratio_col %in% colnames(ase_data))
    ase_data[[p2_ratio_col]] <- 1 - ase_data[[p1_ratio_col]]
  ase_data$is_biased_input <- ase_data$bias_ratio > 0.5
  ase_data$ase_source      <- "legacy"
  cat("ASE input: legacy summary table. NOTE: 'biased' here means bias_ratio > 0.5,\n",
      "  a majority vote over SNPs and not a statistical test.\n", sep = "")
}

# The direction of every bias is read by comparing this column with
# parent1_label, so a label that does not occur in it would silently turn every
# biased gene into "biased toward parent 2". Refuse instead.
bias_values <- unique(as.character(ase_data$predominant_bias))
bias_values <- bias_values[!is.na(bias_values) &
                           !tolower(bias_values) %in% c("", "na", "balanced", "none")]
unknown_values <- setdiff(bias_values, c(parent1_label, parent2_label))
if (length(unknown_values) > 0)
  stop("The ASE table's bias direction column holds value(s) that are neither ",
       "--parent1-label ('", parent1_label, "') nor --parent2-label ('",
       parent2_label, "'): ", paste(head(unknown_values, 5), collapse = ", "),
       ". The two labels must be spelled exactly as in that column.")

cat("Loaded ASE data for", nrow(ase_data), "genes\n\n")

# ---------------------------------------------------------------------------
# Helper: which way round is a DE comparison, relative to parent 1 / parent 2?
# ---------------------------------------------------------------------------
# Returns +1 when the file's logFC already means "parent 1 minus parent 2",
# -1 when it means the reverse, and NA when the comparison is not parent 1
# versus parent 2 at all (parent versus F1, say), for which cis/trans
# classification is not defined.

orient_comparison <- function(cmp) {
  if (!is.null(parent1_group) && !is.null(parent2_group)) {
    if (cmp == paste0(parent1_group, "_vs_", parent2_group)) return(1)
    if (cmp == paste0(parent2_group, "_vs_", parent1_group)) return(-1)
    return(NA_real_)
  }
  parts <- strsplit(cmp, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NA_real_)
  has <- function(x, label) grepl(tolower(label), tolower(x), fixed = TRUE)
  only1 <- function(x) has(x, parent1_label) && !has(x, parent2_label)
  only2 <- function(x) has(x, parent2_label) && !has(x, parent1_label)
  # The two group names must differ ONLY by the parent label. Otherwise
  # "lineA_leaf_vs_lineB_root" would qualify and tissue would be read as a
  # parental difference.
  rest <- function(x, label) sub(tolower(label), "", tolower(x), fixed = TRUE)
  same_condition <- function(a, b)
    rest(a, parent1_label) == rest(b, parent2_label)
  if (only1(parts[1]) && only2(parts[2]) && same_condition(parts[1], parts[2]))
    return(1)
  if (only2(parts[1]) && only1(parts[2]) && same_condition(parts[2], parts[1]))
    return(-1)
  NA_real_
}

# ---------------------------------------------------------------------------
# Helper: classify gene regulatory category from DE + ASE
# ---------------------------------------------------------------------------

categorise_gene <- function(is_de, is_biased, bias_toward_p1, logfc) {
  # `logfc` must already be oriented parent 1 minus parent 2 (see
  # orient_comparison); a positive value then means higher in parent 1.
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
parental_de        <- NULL   # oriented DE of the parental comparison, for step 4
n_classified       <- 0
parental_found     <- character(0)

for (de_file in de_files) {
  cmp <- sub("_DE_results\\.tsv$", "", basename(de_file))
  orientation <- orient_comparison(cmp)
  if (is.na(orientation)) {
    cat("Skipping:", cmp,
        "- not a parent 1 versus parent 2 comparison, so cis/trans",
        "classification does not apply.\n")
    next
  }
  cat("Processing:", cmp,
      if (orientation > 0) "(logFC used as named)" else
        "(logFC sign flipped so that positive = higher in parent 1)", "\n")

  parental_found <- c(parental_found, cmp)
  de_data <- read_tsv(de_file, show_col_types = FALSE)

  # Strip .T[digit] transcript suffix so gene IDs match ASE gene_id
  de_data$gene_id_base <- sub("\\.T[0-9]+$", "", de_data$TranscriptID)
  de_data$logFC_parent1_vs_parent2 <- orientation * de_data$logFC

  if (is.null(parental_de))
    parental_de <- de_data[, c("gene_id_base", "logFC_parent1_vs_parent2", "FDR")]

  merged <- de_data %>%
    left_join(ase_data, by = c("gene_id_base" = "gene_id")) %>%
    filter(!is.na(snp_count))

  if (nrow(merged) == 0) {
    cat("  No genes with ASE data found for this comparison.\n\n")
    next
  }

  cat("  Genes with both DE and ASE data:", nrow(merged), "\n")

  clip <- function(x) pmin(pmax(x, 1e-3), 1 - 1e-3)
  merged <- merged %>%
    mutate(
      is_DE      = FDR < fdr_threshold,
      is_biased  = is_biased_input,
      bias_p1    = predominant_bias == parent1_label,
      ase_de_category = mapply(
        categorise_gene, is_DE, is_biased, bias_p1, logFC_parent1_vs_parent2
      ),
      # The standard decomposition (cis = F1 allelic log-ratio, total = parental
      # log-ratio, trans = total - cis), reported as descriptive quantities.
      cis_log2   = log2(clip(.data[[p1_ratio_col]]) / (1 - clip(.data[[p1_ratio_col]]))),
      trans_log2 = logFC_parent1_vs_parent2 - cis_log2
    ) %>%
    select(-is_biased_input)
  n_classified <- n_classified + 1

  # Export combined results
  write_tsv(merged,
            file.path(output_dir, paste0(cmp, "_with_ASE.tsv")))

  # Category summary
  summary_tbl <- merged %>%
    group_by(ase_de_category) %>%
    summarise(
      n_genes           = n(),
      avg_logFC_parent1_vs_parent2 = mean(logFC_parent1_vs_parent2, na.rm = TRUE),
      avg_parent1_ratio = mean(.data[[p1_ratio_col]], na.rm = TRUE),
      .groups = "drop"
    )
  write_tsv(summary_tbl,
            file.path(output_dir, paste0(cmp, "_ASE_category_summary.tsv")))

  cat("  Category breakdown:\n")
  print(as.data.frame(summary_tbl))

  # Scatter plot: oriented logFC vs parent-1 allele fraction, coloured by category
  p_scatter <- ggplot(
    merged,
    aes(x = logFC_parent1_vs_parent2, y = .data[[p1_ratio_col]],
        color = ase_de_category)
  ) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey50") +
    scale_color_brewer(palette = "Set2", name = "Category") +
    labs(
      title   = paste("ASE bias vs Differential Expression:", cmp),
      x       = paste0("log2 Fold Change (", parent1_label, " / ", parent2_label, ")"),
      y       = paste0(parent1_label, " allele ratio (ASE)"),
      caption = "Upper-right and lower-left quadrants indicate cis-regulatory divergence"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "bottom")

  ggsave(file.path(output_dir, paste0(cmp, "_ASE_vs_DE.pdf")),
         p_scatter, width = 10, height = 7, device = "pdf")
  ggsave(file.path(output_dir, paste0(cmp, "_ASE_vs_DE.png")),
         p_scatter, width = 10, height = 7, dpi = 300, device = "png")

  cis <- merged %>%
    filter(ase_de_category == "cis") %>%
    arrange(FDR)
  if (nrow(cis) > 0)
    all_cis_candidates[[cmp]] <- cis

  cat("  Cis-regulatory candidates:", nrow(cis), "\n\n")
}

if (n_classified == 0 && length(parental_found) > 0) {
  stop(
    "The parental comparison was found (", paste(parental_found, collapse = ", "),
    ") but none of its gene IDs occur in the ASE table, so nothing was ",
    "classified. Compare the TranscriptID column of the DE results (after any ",
    ".T[digit] suffix is removed) with the gene_id column of the ASE table."
  )
}
if (n_classified == 0) {
  stop(
    "No DE comparison could be identified as parent 1 versus parent 2, so nothing ",
    "was classified. Pass --parent1-group and --parent2-group with the two group ",
    "names exactly as they appear in the *_DE_results.tsv file names, or use ",
    "--parent1-label / --parent2-label values that occur in those group names. ",
    "Comparisons found: ",
    paste(sub("_DE_results\\.tsv$", "", basename(de_files)), collapse = ", ")
  )
}

# Collated cis candidates across the parental comparisons
if (length(all_cis_candidates) > 0) {
  cis_all <- bind_rows(all_cis_candidates, .id = "comparison") %>%
    arrange(comparison, FDR)
  write_tsv(cis_all,
            file.path(output_dir, "cis_regulatory_candidates_all.tsv"))
  cat("Saved", nrow(cis_all), "total cis-regulatory candidates across all comparisons\n\n")
}

# ---------------------------------------------------------------------------
# 4. Cis/trans classification with the parents' and the F1's expression
#    (requires expression matrix + metadata)
# ---------------------------------------------------------------------------
#
# The classification is the standard decomposition (Wittkopp et al. 2004;
# McManus et al. 2010):
#
#   cis_log2      = log2 of the F1 allelic ratio, parent 1 / parent 2
#   parental_log2 = log2 of the parental expression ratio, parent 1 / parent 2
#   trans_log2    = parental_log2 - cis_log2
#
#   Cis-only                    allelic imbalance, no trans component
#   Trans-only                  parents differ, alleles balanced in the F1
#   Cis + Trans (Reinforcing)   cis and trans favour the same parent
#   Cis x Trans (Opposing)      cis and trans favour opposite parents
#   Compensatory                allelic imbalance, parents do not differ
#   Conserved                   neither
#   Ambiguous                   parents differ but neither component clears
#                               its threshold
#   Not tested (no DE result)   the gene has no row in the parental DE results
#
# It is descriptive: "cis" means the parent-1 allele fraction is outside
# 0.5 +/- 0.1, "trans" means |trans_log2| >= 1, "parents differ" comes from the
# DE test (see below). There is no test of the trans component itself.
#
# Two further numbers are reported and plotted, and they are NOT the cis and
# trans effects:
#
#   Ad = 0.5 - (parent-1 allele fraction in the F1)
#        the F1's allelic imbalance on the fraction scale. Negative = parent-1
#        bias.
#
#   Ed = 0.5 - (F1 - P2) / (P1 - P2)
#        where the F1's TOTAL expression sits between the two parents. Zero is
#        mid-parent (additive); negative is parent-1-like. Ed measures dominance
#        of expression, which is a different question from cis versus trans: a
#        purely additive trans difference leaves Ed at zero, and an F1 can
#        resemble the parent whose allele it expresses LESS. Earlier versions
#        classified on the signs of Ad and Ed, which labelled a pure-cis,
#        additive gene "Ambiguous" and could label opposing cis and trans
#        effects "Reinforcing".
#
# Ed divides by the parental difference, so it is undefined where the parents do
# not differ and explodes where they barely do. It is computed only for genes
# whose parents are differentially expressed (FDR below the threshold in the
# parent 1 versus parent 2 comparison).

run_advanced <- !is.null(expression_file) &&
                !is.null(metadata_file) &&
                !is.null(hybrid_group) &&
                !is.null(parent1_group) &&
                !is.null(parent2_group)

classify_regulatory <- function(ad, cis_log2, trans_log2, parents_differ,
                                cis_thresh = 0.1, trans_thresh = 1.0) {
  if (is.na(parents_differ)) return("Not tested (no DE result)")
  has_cis   <- abs(ad) > cis_thresh
  has_trans <- abs(trans_log2) >= trans_thresh
  if (!parents_differ)
    return(if (has_cis) "Compensatory" else "Conserved")
  if (has_cis && has_trans) {
    if (sign(cis_log2) == sign(trans_log2)) "Cis + Trans (Reinforcing)"
    else "Cis x Trans (Opposing)"
  } else if (has_cis) {
    "Cis-only"
  } else if (has_trans) {
    "Trans-only"
  } else {
    "Ambiguous"
  }
}

if (run_advanced) {
  cat("--- Cis/trans classification with parental and F1 expression ---\n")

  if (!requireNamespace("edgeR", quietly = TRUE))
    stop("The Ad/Ed analysis needs the Bioconductor package 'edgeR' for TMM ",
         "normalisation. Install it, or drop --expression-file to run the ",
         "DE + ASE integration alone.")
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
  if ("sample" %in% col_lower && !"Sample" %in% colnames(meta))
    colnames(meta)[col_lower == "sample"] <- "Sample"
  if ("group" %in% col_lower) {
    colnames(meta)[col_lower == "group"] <- "Group"
  } else if (all(c("species", "tissue") %in% col_lower)) {
    # Same derivation as EdgeRDE, so the group names match the DE file names.
    meta$Group <- paste(meta[[which(col_lower == "species")]],
                        meta[[which(col_lower == "tissue")]], sep = "_")
  } else {
    stop("Metadata must contain a 'group' column, or both 'species' and 'tissue'.")
  }
  if (!"Sample" %in% colnames(meta))
    stop("Metadata must contain a 'sample_name' or 'Sample' column.")

  # Retain samples present in the expression matrix
  meta <- meta[meta$Sample %in% colnames(expr_matrix), ]
  expr_matrix <- expr_matrix[, meta$Sample, drop = FALSE]

  # TMM normalisation
  dge_ase <- edgeR::DGEList(counts = expr_matrix)
  dge_ase <- edgeR::calcNormFactors(dge_ase, method = "TMM")
  norm_cpm <- edgeR::cpm(dge_ase, normalized.lib.sizes = TRUE)

  # Per-group averages
  p1_samples <- meta$Sample[meta$Group == parent1_group]
  p2_samples <- meta$Sample[meta$Group == parent2_group]
  hy_samples <- meta$Sample[meta$Group == hybrid_group]

  if (length(p1_samples) == 0) stop("No samples found for parent1_group: ", parent1_group)
  if (length(p2_samples) == 0) stop("No samples found for parent2_group: ", parent2_group)
  if (length(hy_samples) == 0) stop("No samples found for hybrid_group: ",  hybrid_group)

  expr_summary <- data.frame(
    # Same .T[digit] stripping as the DE join, so a transcript-level matrix
    # still matches gene-level ASE IDs instead of silently matching nothing.
    gene_id    = sub("\\.T[0-9]+$", "", rownames(norm_cpm)),
    parent1_avg = rowMeans(norm_cpm[, p1_samples, drop = FALSE]),
    parent2_avg = rowMeans(norm_cpm[, p2_samples, drop = FALSE]),
    hybrid_avg  = rowMeans(norm_cpm[, hy_samples, drop = FALSE]),
    stringsAsFactors = FALSE
  )

  combined_ad_ed <- inner_join(rename(ase_data, is_biased = is_biased_input),
                               expr_summary, by = "gene_id")
  if (nrow(combined_ad_ed) == 0)
    stop("No gene IDs are shared between the ASE table and the expression ",
         "matrix (after stripping any .T[digit] suffix). Check the ID format.")

  # Which genes have parents that really differ? Taken from the parental DE
  # test. A gene with no row there (filtered out before testing, typically for
  # low expression) is UNTESTED: parents_differ is NA, not FALSE.
  combined_ad_ed$parental_log2 <- log2((combined_ad_ed$parent1_avg + 0.5) /
                                       (combined_ad_ed$parent2_avg + 0.5))
  de_lookup <- parental_de[!duplicated(parental_de$gene_id_base), ]
  hit <- match(combined_ad_ed$gene_id, de_lookup$gene_id_base)
  combined_ad_ed$parental_FDR   <- de_lookup$FDR[hit]
  combined_ad_ed$parents_differ <- combined_ad_ed$parental_FDR < fdr_threshold
  cat("Parental difference taken from the DE results (FDR <", fdr_threshold, ");",
      sum(is.na(combined_ad_ed$parents_differ)), "gene(s) have no DE result\n")

  clip <- function(x) pmin(pmax(x, 1e-3), 1 - 1e-3)
  combined_ad_ed <- combined_ad_ed %>%
    mutate(
      # Allelic deviation: departure from the 0.5 expected allele fraction.
      # Negative = parent1 bias, positive = parent2 bias
      Ad = 0.5 - .data[[p1_ratio_col]],
      # Parental expression difference
      parental_diff = parent1_avg - parent2_avg,
      # F1 position on the parental spectrum (1 = parent1-like, 0 = parent2-like)
      position_ratio = ifelse(
        !is.na(parents_differ) & parents_differ & parental_diff != 0,
        (hybrid_avg - parent2_avg) / parental_diff,
        NA_real_
      ),
      # Expression-inheritance deviation: departure from mid-parent expression
      Ed = 0.5 - position_ratio,
      # Standard decomposition, for reference
      cis_log2   = log2(clip(.data[[p1_ratio_col]]) / (1 - clip(.data[[p1_ratio_col]]))),
      trans_log2 = parental_log2 - cis_log2
    ) %>%
    filter(!is.na(Ad))

  combined_ad_ed$Regulation <- mapply(classify_regulatory,
                                      combined_ad_ed$Ad,
                                      combined_ad_ed$cis_log2,
                                      combined_ad_ed$trans_log2,
                                      combined_ad_ed$parents_differ)

  cat("Classification summary:\n")
  print(table(combined_ad_ed$Regulation))
  cat("\n")

  write_tsv(combined_ad_ed,
            file.path(output_dir, "ase_vs_de_ad_ed_data.tsv"))

  # Regression plot per class
  plot_data <- combined_ad_ed %>%
    filter(!is.na(Ed), Ed >= -1.5, Ed <= 1.5,
           !Regulation %in% c("Ambiguous", "Not tested (no DE result)"))

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

    class_colours <- c("Cis + Trans (Reinforcing)" = "#E41A1C",
                       "Cis x Trans (Opposing)"    = "#377EB8",
                       "Trans-only"                = "grey50",
                       "Cis-only"                  = "#4DAF4A")
    label_colours <- setNames(class_colours[group_stats$Regulation],
                              group_stats$label)

    p_reg <- ggplot(plot_data,
                    aes(x = Ad, y = Ed, color = label, fill = label)) +
      geom_point(alpha = 0.55, size = 1.8) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "black") +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "black") +
      scale_color_manual(values = label_colours) +
      scale_fill_manual(values = label_colours) +
      labs(
        title    = "Allelic deviation vs expression inheritance",
        subtitle = paste("Ad: F1 allelic imbalance. Ed: F1 departure from mid-parent",
                         "expression (dominance). Parentally DE genes only;",
                         "colour = cis/trans class."),
        x        = paste0("Allelic Deviation (Ad)\n",
                          "(negative = ", parent1_label, " bias)"),
        y        = paste0("Expression-inheritance Deviation (Ed)\n",
                          "(negative = F1 is ", parent1_label, "-like)"),
        color    = "Class",
        fill     = "Class"
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
