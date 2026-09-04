#!/usr/bin/env Rscript
# Estimates cell-type proportions from bulk expression data using two
# reference-based deconvolution methods:
#
#   1) CIBERSORTx - support-vector regression per sample, mirroring the
#      linear nu-SVR core of the CIBERSORT/CIBERSORTx algorithm: for each
#      sample, several candidate nu values are fit, the fit with the lowest
#      residual sum of squares is kept, negative coefficients are clipped to
#      zero, and the result is renormalized to sum to 1 (e1071 package,
#      nu-regression with a linear kernel). Note: this reproduces
#      CIBERSORTx's core regression step only; the official CIBERSORTx web
#      tool additionally offers batch-correction (B-mode/S-mode) and
#      single-cell reference building, which require its licensed Docker
#      image/token and are out of scope for a local, redistributable script.
#   2) Marker-based proxy scoring - a simpler, regression-free method: bulk
#      expression is z-scored gene-wise across samples, and each cell type's
#      score per sample is the mean z-score of that cell type's marker genes
#      (the genes marking it in the signature matrix). Scores are clipped at
#      zero and renormalized to sum to 1 per sample so they are comparable,
#      as rough proportions, to the CIBERSORTx output.
#
# Inputs:
# - Bulk expression CSV: a gene identifier column (auto-detected, e.g.
#   "Gene Symbol" or the first column) plus one numeric column per sample.
#   Non-numeric columns (titles, probe IDs, etc.) are dropped automatically.
# - A reference signature, provided as either:
#   1) --signature_file: a full signature matrix CSV (gene identifier
#      column + one numeric column per cell type), or
#   2) --markers_csv: a marker-gene CSV (long format: cell_type, gene) as
#      used by analyze_cell_population_proxies.py, turned into a binary
#      presence/absence signature matrix (1 if the gene marks that cell
#      type, 0 otherwise) restricted to genes seen in the bulk data.
#
# Both matrices are aligned on their shared genes before deconvolution.
#
# Required packages (install once):
#   install.packages(c("optparse", "e1071"))
#
# Usage:
#   Rscript deconvolution_methods.R \
#     --bulk_file=PBMC/GSE234669/expression_analysis_results/deduplicated_expression_dataset.csv \
#     --markers_csv=markers/pbmc_cell_markers.csv \
#     --output_dir=PBMC/GSE234669/deconvolution_results_r

suppressPackageStartupMessages({
  library(optparse)
  library(e1071)
})

option_list <- list(
  make_option("--bulk_file", type = "character"),
  make_option("--signature_file", type = "character", default = NULL),
  make_option("--markers_csv", type = "character", default = NULL),
  make_option("--output_dir", type = "character"),
  make_option("--bulk_gene_col", type = "character", default = NULL),
  make_option("--signature_gene_col", type = "character", default = NULL),
  make_option("--celltype_column", type = "character", default = "auto"),
  make_option("--gene_column", type = "character", default = "auto")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$signature_file) && is.null(opt$markers_csv)) {
  stop("Provide either --signature_file or --markers_csv")
}
if (!is.null(opt$signature_file) && !is.null(opt$markers_csv)) {
  stop("Provide only one of --signature_file or --markers_csv")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helpers ---------------------------------------------------------------

choose_gene_identifier_column <- function(df, gene_col = NULL) {
  if (!is.null(gene_col)) {
    if (!gene_col %in% names(df)) stop(sprintf("Gene column not found: %s", gene_col))
    return(gene_col)
  }
  candidates <- c("Gene Symbol", "Gene_Symbol", "gene_symbol", "gene")
  for (col in candidates) if (col %in% names(df)) return(col)
  names(df)[1]
}

# Loads a gene x sample (or gene x cell-type) numeric matrix from a CSV,
# collapsing duplicate gene rows by mean, mirroring load_matrix() in
# deconvolution_methods.py.
load_matrix <- function(path, gene_col = NULL) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- choose_gene_identifier_column(df, gene_col)

  genes <- as.character(df[[gene_col]])
  value_cols <- setdiff(names(df), gene_col)
  values <- df[, value_cols, drop = FALSE]
  values <- as.data.frame(lapply(values, function(col) suppressWarnings(as.numeric(col))))
  names(values) <- value_cols

  # Drop columns that are mostly non-numeric (e.g. a free-text "Gene Title"
  # column where a stray literal "0" happens to survive as.numeric()): a
  # real sample column should be numeric throughout, so require at least
  # half of its values to have parsed successfully.
  mostly_non_numeric_cols <- vapply(values, function(col) mean(is.na(col)) > 0.5, logical(1))
  values <- values[, !mostly_non_numeric_cols, drop = FALSE]
  value_cols <- names(values)

  # Collapse duplicate gene rows by mean.
  unique_genes <- unique(genes)
  collapsed <- matrix(NA_real_, nrow = length(unique_genes), ncol = ncol(values),
                       dimnames = list(unique_genes, value_cols))
  for (g in unique_genes) {
    rows <- which(genes == g)
    if (length(rows) == 1) {
      collapsed[g, ] <- as.numeric(values[rows, ])
    } else {
      collapsed[g, ] <- colMeans(values[rows, , drop = FALSE], na.rm = TRUE)
    }
  }
  collapsed
}

# Builds a genes x cell-types binary signature matrix from a long-format
# marker CSV (cell_type, gene), restricted to genes present in the bulk data
# (case-insensitive match), mirroring build_binary_signature_from_markers().
build_binary_signature_from_markers <- function(markers_csv, genes_universe,
                                                  celltype_column = "auto",
                                                  gene_column = "auto") {
  marker_df <- read.csv(markers_csv, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(marker_df) == 0) stop(sprintf("Marker CSV is empty: %s", markers_csv))

  lower_map <- setNames(names(marker_df), tolower(trimws(names(marker_df))))

  if (celltype_column == "auto") {
    for (candidate in c("cell_type", "celltype", "cell", "population", "immune_cell")) {
      if (candidate %in% names(lower_map)) { celltype_column <- lower_map[[candidate]]; break }
    }
  }
  if (gene_column == "auto") {
    for (candidate in c("gene", "gene_symbol", "symbol", "marker", "marker_gene")) {
      if (candidate %in% names(lower_map)) { gene_column <- lower_map[[candidate]]; break }
    }
  }
  if (!celltype_column %in% names(marker_df) || !gene_column %in% names(marker_df)) {
    stop("Could not auto-detect cell_type/gene columns in marker CSV; pass --celltype_column/--gene_column.")
  }

  cell_types <- sort(unique(as.character(marker_df[[celltype_column]])))
  upper_to_gene <- setNames(genes_universe, toupper(genes_universe))

  signature <- matrix(0, nrow = length(genes_universe), ncol = length(cell_types),
                       dimnames = list(genes_universe, cell_types))
  for (i in seq_len(nrow(marker_df))) {
    ct <- as.character(marker_df[[celltype_column]][i])
    marker <- toupper(as.character(marker_df[[gene_column]][i]))
    if (marker %in% names(upper_to_gene)) {
      gene <- upper_to_gene[[marker]]
      signature[gene, ct] <- 1
    }
  }

  keep_cols <- colSums(signature) > 0
  signature <- signature[, keep_cols, drop = FALSE]
  if (ncol(signature) == 0) {
    stop(sprintf("No marker genes from %s matched genes in the bulk dataset.", markers_csv))
  }
  signature
}

# Aligns bulk and signature matrices on their shared, complete-case genes,
# mirroring align_matrices(). Also drops genes that are all-zero across
# every cell type in the signature (i.e. non-marker genes): they are inert
# for the marker-based methods (a zero row carries no marker information)
# but, left in, they act as thousands of near-identical noise points for
# the nu-SVR fit and make it needlessly slow, so it's both faster and more
# principled to restrict deconvolution to genes that actually carry
# signature information.
align_matrices <- function(bulk, signature) {
  shared_genes <- intersect(rownames(bulk), rownames(signature))
  if (length(shared_genes) < 10) {
    stop(sprintf("Only %d genes overlap between bulk and signature matrices; check that gene identifiers match.",
                  length(shared_genes)))
  }
  bulk_aligned <- bulk[shared_genes, , drop = FALSE]
  bulk_aligned <- bulk_aligned[rowSums(is.na(bulk_aligned)) == 0, , drop = FALSE]
  signature_aligned <- signature[rownames(bulk_aligned), , drop = FALSE]
  signature_aligned <- signature_aligned[rowSums(is.na(signature_aligned)) == 0, , drop = FALSE]
  shared_genes <- intersect(rownames(bulk_aligned), rownames(signature_aligned))
  bulk_aligned <- bulk_aligned[shared_genes, , drop = FALSE]
  signature_aligned <- signature_aligned[shared_genes, , drop = FALSE]

  informative_genes <- rowSums(abs(signature_aligned)) > 0
  list(bulk = bulk_aligned[informative_genes, , drop = FALSE],
       signature = signature_aligned[informative_genes, , drop = FALSE])
}

# Method 1: CIBERSORTx (linear nu-SVR core) ----------------------------------
# Standardizes the signature matrix, fits nu-SVR for each candidate nu on
# each sample, keeps the fit with the lowest residual sum of squares, clips
# negative coefficients to zero, and renormalizes to sum to 1.
deconvolve_cibersortx <- function(bulk, signature, nus = c(0.25, 0.5, 0.75)) {
  cell_types <- colnames(signature)
  s_matrix <- signature
  s_mean <- colMeans(s_matrix)
  s_sd <- apply(s_matrix, 2, sd)
  s_sd[s_sd == 0] <- 1.0
  s_scaled <- scale(s_matrix, center = s_mean, scale = s_sd)

  proportions <- matrix(NA_real_, nrow = ncol(bulk), ncol = length(cell_types),
                         dimnames = list(colnames(bulk), cell_types))

  for (sample in colnames(bulk)) {
    b_vector <- bulk[, sample]
    b_mean <- mean(b_vector)
    b_sd <- sd(b_vector); if (b_sd == 0 || is.na(b_sd)) b_sd <- 1.0
    b_scaled <- (b_vector - b_mean) / b_sd

    best_coefs <- NULL
    best_rss <- Inf
    for (nu in nus) {
      model <- e1071::svm(x = s_scaled, y = b_scaled, type = "nu-regression",
                           kernel = "linear", nu = nu, cost = 1.0, scale = FALSE)
      # Linear-kernel primal coefficients: w = t(SV) %*% coefs
      coefs <- as.numeric(t(model$SV) %*% model$coefs)
      predicted <- s_scaled %*% coefs
      rss <- sum((predicted - b_scaled)^2)
      if (rss < best_rss) {
        best_rss <- rss
        best_coefs <- coefs
      }
    }

    best_coefs[best_coefs < 0] <- 0
    total <- sum(best_coefs)
    proportions[sample, ] <- if (total > 0) best_coefs / total else best_coefs
  }
  as.data.frame(proportions, check.names = FALSE)
}

# Method 2: Marker-based proxy scoring ---------------------------------------
# Regression-free scoring: bulk expression is z-scored per gene across
# samples, and each cell type's proxy score per sample is the mean z-score
# of that cell type's marker genes (genes with a positive signature value
# for that cell type). Scores are clipped at zero and renormalized to sum
# to 1 per sample so they read as rough proportions, comparable to the
# CIBERSORTx output, mirroring the marker-averaging approach used in
# analyze_cell_population_proxies.py.
deconvolve_marker_proxy <- function(bulk, signature) {
  cell_types <- colnames(signature)

  gene_means <- rowMeans(bulk)
  gene_sds <- apply(bulk, 1, sd)
  gene_sds[gene_sds == 0 | is.na(gene_sds)] <- 1.0
  z_bulk <- sweep(bulk, 1, gene_means, "-")
  z_bulk <- sweep(z_bulk, 1, gene_sds, "/")

  proportions <- matrix(NA_real_, nrow = ncol(bulk), ncol = length(cell_types),
                         dimnames = list(colnames(bulk), cell_types))

  for (ct in cell_types) {
    marker_genes <- rownames(signature)[signature[, ct] > 0]
    proportions[, ct] <- if (length(marker_genes) > 0) {
      colMeans(z_bulk[marker_genes, , drop = FALSE])
    } else {
      0
    }
  }

  proportions[proportions < 0] <- 0
  row_totals <- rowSums(proportions)
  for (i in seq_len(nrow(proportions))) {
    if (row_totals[i] > 0) proportions[i, ] <- proportions[i, ] / row_totals[i]
  }
  as.data.frame(proportions, check.names = FALSE)
}

save_comparison_plots <- function(cibersortx_props, proxy_props, output_dir) {
  cell_types <- colnames(cibersortx_props)

  cibersortx_means <- colMeans(cibersortx_props)
  proxy_means <- colMeans(proxy_props[, cell_types, drop = FALSE])

  png(file.path(output_dir, "deconvolution_method_comparison_bar.png"),
      width = max(800, length(cell_types) * 90), height = 600)
  bar_matrix <- rbind(CIBERSORTx = cibersortx_means, `Marker Proxy` = proxy_means)
  bp <- barplot(bar_matrix, beside = TRUE, col = c("#66c2a5", "#fc8d62"),
                las = 2, cex.names = 0.8,
                main = "Mean Estimated Cell-Type Proportions by Method",
                ylab = "Mean estimated proportion")
  legend("topright", legend = rownames(bar_matrix), fill = c("#66c2a5", "#fc8d62"))
  dev.off()

  shared_cell_types <- intersect(colnames(cibersortx_props), colnames(proxy_props))
  flat_cibersortx <- as.vector(as.matrix(cibersortx_props[, shared_cell_types, drop = FALSE]))
  flat_proxy <- as.vector(as.matrix(proxy_props[, shared_cell_types, drop = FALSE]))
  corr <- if (length(flat_cibersortx) > 1) cor(flat_cibersortx, flat_proxy) else NA_real_

  png(file.path(output_dir, "deconvolution_method_agreement_scatter.png"),
      width = 700, height = 700)
  lim <- max(c(flat_cibersortx, flat_proxy, 0.01), na.rm = TRUE)
  plot(flat_cibersortx, flat_proxy, xlim = c(0, lim), ylim = c(0, lim), pch = 19,
       col = rgb(0, 0, 0, 0.4),
       xlab = "CIBERSORTx proportion", ylab = "Marker proxy proportion",
       main = sprintf("CIBERSORTx vs. Marker Proxy Estimates (Pearson r = %.2f)", corr))
  abline(a = 0, b = 1, lty = 2, col = "gray")
  dev.off()
}

# --- Run ---------------------------------------------------------------

bulk <- load_matrix(opt$bulk_file, gene_col = opt$bulk_gene_col)
if (!is.null(opt$signature_file)) {
  signature <- load_matrix(opt$signature_file, gene_col = opt$signature_gene_col)
} else {
  signature <- build_binary_signature_from_markers(
    opt$markers_csv, rownames(bulk),
    celltype_column = opt$celltype_column, gene_column = opt$gene_column
  )
}

aligned <- align_matrices(bulk, signature)
bulk_aligned <- aligned$bulk
signature_aligned <- aligned$signature

cibersortx_props <- deconvolve_cibersortx(bulk_aligned, signature_aligned)
proxy_props <- deconvolve_marker_proxy(bulk_aligned, signature_aligned)

write.csv(cbind(sample = rownames(cibersortx_props), cibersortx_props),
          file.path(opt$output_dir, "cibersortx_cell_type_proportions.csv"), row.names = FALSE)
write.csv(cbind(sample = rownames(proxy_props), proxy_props),
          file.path(opt$output_dir, "marker_proxy_cell_type_proportions.csv"), row.names = FALSE)
save_comparison_plots(cibersortx_props, proxy_props, opt$output_dir)

cat(sprintf("Genes used for deconvolution: %d\n", nrow(bulk_aligned)))
cat(sprintf("Saved deconvolution outputs to %s\n", opt$output_dir))
cat("\nCIBERSORTx proportions (head):\n")
print(head(cibersortx_props))
cat("\nMarker proxy proportions (head):\n")
print(head(proxy_props))
