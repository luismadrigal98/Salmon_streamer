#!/usr/bin/env Rscript

# Aggregate Salmon transcript quantifications to paralog-aware features with
# tximport.
#
# Wrapped by subprograms/ParalogTximport.py -- see there for the CLI.
#
# Why this exists rather than a hand-summed matrix: summing counts across the
# transcripts of a merged feature loses the average transcript length, and that
# length is what edgeR/DESeq2 need to correct for composition differences between
# samples. tximport recomputes it. The two supported ways of using the result:
#
#   countsFromAbundance = "no"
#       Raw counts plus a length matrix. The most faithful route, but the
#       downstream tool must accept the length matrix as an offset.
#   countsFromAbundance = "lengthScaledTPM"   (the pipeline default)
#       The length correction is folded into the counts themselves, so the matrix
#       can be handed to any count-based DE tool directly, with no offset. This is
#       what makes EdgeRDE work without modification.
#
# Args (positional):
#   1 manifest TSV: sample_name <tab> quant_path
#   2 tx2feature TSV: transcript <tab> feature  (ParalogMerge --out-tx2gene)
#   3 output directory
#   4 countsFromAbundance: no | scaledTPM | lengthScaledTPM
#   5 ignoreTxVersion: TRUE | FALSE
#
# @Author: Luis Javier Madrigal-Roca & John K. Kelly

suppressPackageStartupMessages({
    library(tximport)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("usage: tximport_features.R <manifest> <tx2feature> <outdir> <countsFromAbundance> <ignoreTxVersion>")
}

manifest_path <- args[1]
t2f_path      <- args[2]
out_dir       <- args[3]
cfa           <- args[4]
ignore_ver    <- toupper(args[5]) == "TRUE"

valid_cfa <- c("no", "scaledTPM", "lengthScaledTPM")
if (!cfa %in% valid_cfa) {
    stop(sprintf("countsFromAbundance must be one of %s (got '%s')",
                 paste(valid_cfa, collapse = ", "), cfa))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- inputs ---------------------------------------------------------------
manifest <- read.delim(manifest_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)
if (!all(c("sample_name", "quant_path") %in% colnames(manifest))) {
    stop("manifest must have columns: sample_name, quant_path")
}

files <- manifest$quant_path
names(files) <- manifest$sample_name

missing <- files[!file.exists(files)]
if (length(missing) > 0) {
    stop(sprintf("%d quant.sf file(s) not found, first: %s",
                 length(missing), missing[1]))
}

t2f <- read.delim(t2f_path, header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE)
if (ncol(t2f) < 2) {
    stop("tx2feature must have at least two columns: transcript, feature")
}
t2f <- t2f[, 1:2]
colnames(t2f) <- c("transcript", "feature")

cat(sprintf("samples          : %d\n", length(files)))
cat(sprintf("tx2feature rows  : %d  (%d distinct features)\n",
            nrow(t2f), length(unique(t2f$feature))))
cat(sprintf("countsFromAbundance: %s\n", cfa))

# ---- aggregate ------------------------------------------------------------
txi <- tximport(files,
                type                = "salmon",
                tx2gene             = t2f,
                countsFromAbundance = cfa,
                ignoreTxVersion     = ignore_ver)

write_matrix <- function(mat, path) {
    df <- data.frame(feature = rownames(mat), mat,
                     check.names = FALSE, stringsAsFactors = FALSE)
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

counts_path <- file.path(out_dir, "paralog_feature_counts.tsv")
tpm_path    <- file.path(out_dir, "paralog_feature_abundance_tpm.tsv")
len_path    <- file.path(out_dir, "paralog_feature_length.tsv")

write_matrix(txi$counts,    counts_path)
write_matrix(txi$abundance, tpm_path)
write_matrix(txi$length,    len_path)

# ---- report ---------------------------------------------------------------
n_feat <- nrow(txi$counts)
n_para <- sum(grepl("^PARA_", rownames(txi$counts)))
cat("\n")
cat(sprintf("features written : %d  (%d merged PARA_ features, %d single-gene)\n",
            n_feat, n_para, n_feat - n_para))
cat(sprintf("total counts     : %.0f\n", sum(txi$counts)))
cat(sprintf("counts           -> %s\n", counts_path))
cat(sprintf("abundance (TPM)  -> %s\n", tpm_path))
cat(sprintf("avg tx length    -> %s\n", len_path))

if (cfa == "no") {
    cat("\nNOTE: countsFromAbundance='no', so these counts are NOT length-corrected.\n")
    cat("      Supply paralog_feature_length.tsv to your DE tool as an offset, or\n")
    cat("      re-run with --counts-from-abundance lengthScaledTPM to feed EdgeRDE\n")
    cat("      directly.\n")
} else {
    cat("\nThe length correction is folded into the counts, so this matrix can go\n")
    cat("straight to EdgeRDE (or any count-based DE tool) with no offset.\n")
}
