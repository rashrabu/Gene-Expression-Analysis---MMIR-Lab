#!/usr/bin/env Rscript
# QC, normalization, and batch correction for the annotated/renamed/unlogged
# expression matrix produced by preprocessing.py.
#
# R equivalent of qc_normalize.py. Same steps:
#   1. QC diagnostics on the raw matrix (missingness, distributions,
#      sample correlation, hierarchical clustering, PCA by phenotype/batch)
#   2. Quantile normalization (preprocessCore::normalize.quantiles)
#   3. ComBat batch correction (sva::ComBat), protecting the phenotype
#      column as a covariate
#   4. QC diagnostics again on the corrected matrix
#   5. Writes qc_normalized_expression.csv + qc_report.txt
#
# Required packages (install once):
#   install.packages("optparse")
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("preprocessCore", "sva"))
#
# Usage:
#   Rscript qc_normalize.R \
#     --expression_file=processed_data/annotated_renamed_unlogged.csv \
#     --metadata_file=processed_data/sample_metadata.csv \
#     --output_dir=qc_normalized \
#     --metadata_sample_column=analysis_sample \
#     --phenotype_column=depression_status \
#     --batch_column=batch

suppressPackageStartupMessages({
  library(optparse)
  library(preprocessCore)
  library(limma)
  library(matrixStats)
})

option_list <- list(
  make_option("--expression_file", type = "character"),
  make_option("--metadata_file", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--metadata_sample_column", type = "character", default = "analysis_sample"),
  make_option("--phenotype_column", type = "character", default = "depression_status"),
  make_option("--batch_column", type = "character", default = "batch")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

choose_gene_identifier_column <- function(df) {
  candidates <- c("Gene Symbol", "Gene_Symbol", "gene_symbol", "SYMBOL", "Symbol")
  for (col in candidates) if (col %in% names(df)) return(col)
  fallbacks <- c("ID_REF", "ID", "Probe Set ID")
  for (col in fallbacks) if (col %in% names(df)) return(col)
  stop("Could not identify a gene identifier column.")
}

fill_row_means <- function(mat) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(idx) > 0) mat[idx] <- row_means[idx[, 1]]
  mat
}

plot_distributions <- function(mat, suffix) {
  png(file.path(opt$output_dir, sprintf("qc_sample_distributions_%s.png", suffix)),
      width = max(800, ncol(mat) * 15), height = 600)
  boxplot(mat, outline = FALSE, las = 2, cex.axis = 0.5,
          main = sprintf("Sample Expression Distributions (%s)", suffix),
          ylab = "Expression")
  dev.off()
}

plot_correlation_and_clustering <- function(mat, suffix) {
  filled <- fill_row_means(mat)
  corr <- cor(filled)

  png(file.path(opt$output_dir, sprintf("qc_sample_correlation_%s.png", suffix)),
      width = 900, height = 900)
  heatmap(corr, symm = TRUE, labRow = NA, labCol = NA,
          main = sprintf("Sample-Sample Correlation (%s)", suffix))
  dev.off()

  hc <- hclust(as.dist(1 - corr), method = "average")
  png(file.path(opt$output_dir, sprintf("qc_hierarchical_clustering_%s.png", suffix)),
      width = max(800, ncol(mat) * 15), height = 600)
  plot(hc, main = sprintf("Hierarchical Clustering of Samples (%s)", suffix),
       xlab = "", sub = "", cex = 0.5)
  dev.off()
}

plot_pca <- function(mat, metadata_df, suffix) {
  filled <- fill_row_means(mat)
  filled[is.na(filled)] <- 0
  row_var <- rowVars(filled)
  filled <- filled[row_var > 0, , drop = FALSE]
  scaled <- scale(t(filled))
  pca <- prcomp(scaled, center = FALSE, scale. = FALSE)
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100

  pc_df <- data.frame(sample = colnames(mat), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  available_meta_cols <- intersect(
    c(opt$metadata_sample_column, opt$phenotype_column, opt$batch_column),
    names(metadata_df)
  )
  pc_df <- merge(
    pc_df,
    metadata_df[, available_meta_cols, drop = FALSE],
    by.x = "sample", by.y = opt$metadata_sample_column, all.x = TRUE
  )

  for (color_col in c(opt$phenotype_column, opt$batch_column)) {
    if (!color_col %in% names(pc_df)) next
    groups <- as.factor(pc_df[[color_col]])
    if (length(levels(groups)) == 0) {
      cat(sprintf("Skipping PCA plot colored by '%s': no non-missing values found.\n", color_col))
      next
    }
    png(file.path(opt$output_dir, sprintf("qc_pca_%s_by_%s.png", suffix, color_col)),
        width = 800, height = 600)
    plot(pc_df$PC1, pc_df$PC2, col = as.numeric(groups), pch = 19,
         xlab = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
         ylab = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
         main = sprintf("PCA (%s) colored by %s", suffix, color_col))
    legend("topright", legend = levels(groups), col = seq_along(levels(groups)), pch = 19, cex = 0.7)
    dev.off()
  }
}

# --- Load data ---
expression_df <- read.csv(opt$expression_file, check.names = FALSE, stringsAsFactors = FALSE)
metadata_df   <- read.csv(opt$metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

expression_df$ROW_ID__ <- as.character(seq_len(nrow(expression_df)))

gene_col <- choose_gene_identifier_column(expression_df)

sample_cols <- as.character(metadata_df[[opt$metadata_sample_column]])
sample_cols <- sample_cols[sample_cols %in% names(expression_df)]
if (length(sample_cols) == 0) stop("No metadata sample IDs matched expression columns.")
other_cols <- setdiff(names(expression_df), sample_cols)

matrix_mat <- as.matrix(expression_df[, sample_cols, drop = FALSE])
# rm(expression_df)
# gc()
# result_df <- merge(annotation_df,
#                    corrected_df,
#                    by = gene_col,
#                    all.x = TRUE)
storage.mode(matrix_mat) <- "double"
# rownames(matrix_mat) <- as.character(expression_df[[gene_col]])  # matrices allow duplicate rownames
rownames(matrix_mat) <- expression_df$ROW_ID__

# Drop the bottom half of genes/probes by variance -- same rationale as
# the Python version: they carry little signal and add memory/compute cost.
# gene_var <- apply(matrix_mat, 1, var, na.rm = TRUE)
# keep <- gene_var >= median(gene_var, na.rm = TRUE)
# matrix_mat <- matrix_mat[keep, , drop = FALSE]
gene_var <- rowVars(matrix_mat, na.rm = TRUE)
keep <- gene_var >= median(gene_var, na.rm = TRUE)
matrix_mat <- matrix_mat[keep, , drop = FALSE]

cat(sprintf("1. Matrix built: %d genes x %d samples\n", nrow(matrix_mat), ncol(matrix_mat)))

# --- QC on raw matrix ---
plot_distributions(matrix_mat, "before")
cat("2. Boxplot complete\n")
plot_correlation_and_clustering(matrix_mat, "before")
cat("3. Correlation complete\n")
plot_pca(matrix_mat, metadata_df, "before")
cat("4. PCA complete\n")

# --- Quantile normalization ---
filled_for_norm <- fill_row_means(matrix_mat)
normalized <- normalize.quantiles(filled_for_norm)
dimnames(normalized) <- dimnames(filled_for_norm)
rm(filled_for_norm)
gc()
plot_distributions(normalized, "after_normalization")
cat("5. Quantile normalization complete\n")

# --- Batch correction ---
aligned_meta <- metadata_df[match(colnames(normalized), metadata_df[[opt$metadata_sample_column]]), ]
batch_labels <- as.character(aligned_meta[[opt$batch_column]])

corrected <- normalized
batch_status <- sprintf("Skipped: no usable '%s' column/levels in metadata.", opt$batch_column)

has_batch_col <- opt$batch_column %in% names(metadata_df)
if (has_batch_col && !any(is.na(batch_labels)) && !any(batch_labels == "") &&
    length(unique(batch_labels)) >= 2) {

  # phenotype <- as.character(aligned_meta[[opt$phenotype_column]])
  # mod <- if (length(unique(phenotype)) > 1) model.matrix(~as.factor(phenotype)) else NULL

  # corrected <- tryCatch({
  #   ComBat(dat = normalized, batch = batch_labels, mod = mod, par.prior = TRUE)
  # }, error = function(e) {
  #   cat(sprintf("ComBat failed (%s); falling back to normalized data.\n", conditionMessage(e)))
  #   normalized
  # })
  # batch_status <- sprintf("Applied ComBat using '%s' with %d levels.",
  #                          opt$batch_column, length(unique(batch_labels)))
  phenotype <- as.character(aligned_meta[[opt$phenotype_column]])
  design <- if (length(unique(phenotype)) > 1) model.matrix(~as.factor(phenotype)) else matrix(1, ncol(normalized), 1)

  corrected <- tryCatch({
    removeBatchEffect(normalized, batch = batch_labels, design = design)
  }, error = function(e) {
    cat(sprintf("removeBatchEffect failed (%s); falling back to normalized data.\n", conditionMessage(e)))
    normalized
  })
  rm(normalized)
  gc()
  batch_status <- sprintf("Applied limma::removeBatchEffect using '%s' with %d levels.",
                           opt$batch_column, length(unique(batch_labels)))
}
cat(sprintf("6. Batch correction: %s\n", batch_status))

# --- QC on corrected matrix ---
plot_correlation_and_clustering(corrected, "after")
plot_pca(corrected, metadata_df, "after")
cat("7. QC on corrected matrix complete\n")

# --- Write corrected expression table, preserving non-sample columns ---
# corrected_df <- data.frame(GENE_ID_COL = rownames(corrected),
#                             as.data.frame(corrected, check.names = FALSE),
#                             row.names = NULL, check.names = FALSE)
# names(corrected_df)[1] <- gene_col

# result_df <- merge(expression_df[, other_cols, drop = FALSE], corrected_df,
#                     by = gene_col, all.x = TRUE)

corrected_df <- data.frame(ROW_ID__ = rownames(corrected),
                            as.data.frame(corrected, check.names = FALSE),
                            row.names = NULL, check.names = FALSE)

result_df <- merge(expression_df[, other_cols, drop = FALSE], corrected_df,
                    by = "ROW_ID__", all.x = FALSE)
result_df$ROW_ID__ <- NULL


output_path <- file.path(opt$output_dir, "qc_normalized_expression.csv")
write.csv(result_df, output_path, row.names = FALSE)

report_lines <- c(
  "QC / normalization / batch-correction report",
  "",
  sprintf("Samples: %d", length(sample_cols)),
  sprintf("Genes/probes: %d", nrow(matrix_mat)),
  "",
  "Normalization: quantile normalization applied across all samples.",
  "",
  sprintf("Batch correction: %s", batch_status),
  "",
  "Review before/after PCA and correlation plots to confirm major",
  "variation aligns with phenotype rather than batch after correction."
)
writeLines(report_lines, file.path(opt$output_dir, "qc_report.txt"))

cat(sprintf("Saved QC/normalization outputs to %s\n", opt$output_dir))
cat(output_path, "\n")  # last line: the output CSV path, for the Python wrapper to capture
