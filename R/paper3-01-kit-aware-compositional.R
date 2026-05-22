#!/usr/bin/env Rscript
# Kit-aware compositional within-sample scoring methods.
#
# Reproducible command:
#   cd /umed-projekty/JAJNIKI/OmicSelector_paper
#   Rscript code/methods/kit_aware_compositional.R
#
# Outputs:
#   results/reviewer_package_v5/supplement/table_S45_kit_aware_compositional_benchmarks.tsv
#   results/reviewer_package_v5/figures/figure_S45_kit_aware_methods.png
#   results/reviewer_package_v5/figures/figure_S45_kit_aware_methods.pdf
#   results/reviewer_package_v5/figures/figure_S45_caption.md

if (!exists(".ka_validate_inputs", mode = "function")) {
  stop("Kit-aware helper functions are unavailable.", call. = FALSE)
}

#' Fit kit-stratified rCLR centering sets.
#'
#' @param expr_matrix Numeric matrix, samples x features, counts or non-negative
#'   abundance values.
#' @param sample_meta data.frame with one row per sample.
#' @param kit_label_col Column in `sample_meta` containing library-kit family.
#'   Missing/unknown kit values fall back to `biofluid_col`, then global.
#' @param biofluid_col Optional fallback column. Default "biofluid".
#' @param trim_upper,trim_lower Feature trimming fractions used to choose
#'   centering features from the training stratum.
#' @param min_centering_size Minimum centering set size before global fallback.
#' @param min_samples_per_stratum Minimum training samples needed for a
#'   stratum-specific centering set.
#' @param exclude_features Features excluded from denominator selection.
#' @return A fit object consumed by `score_kit_stratified_rclr()`.
#' @export
fit_kit_stratified_rclr <- function(expr_matrix,
                                    sample_meta,
                                    kit_label_col,
                                    biofluid_col = "biofluid",
                                    trim_upper = 0.10,
                                    trim_lower = 0.05,
                                    min_centering_size = 8L,
                                    min_samples_per_stratum = 5L,
                                    exclude_features = c("hsa-miR-451a",
                                                         "hsa-miR-16-5p",
                                                         "hsa-miR-486-5p",
                                                         "hsa-miR-144-3p",
                                                         "hsa-miR-223-3p")) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col)
  X <- inp$X
  strata <- .ka_resolve_strata(inp$sample_meta, kit_label_col, biofluid_col)
  global_center <- .ka_trimmed_centering_features(
    X, trim_upper = trim_upper, trim_lower = trim_lower,
    exclude_features = exclude_features, min_centering_size = min_centering_size)

  centers <- list()
  for (s in sort(unique(strata$label))) {
    idx <- which(strata$label == s)
    if (length(idx) < min_samples_per_stratum) next
    cfeat <- .ka_trimmed_centering_features(
      X[idx, , drop = FALSE], trim_upper = trim_upper, trim_lower = trim_lower,
      exclude_features = exclude_features, min_centering_size = min_centering_size)
    if (length(cfeat) >= min_centering_size) centers[[s]] <- cfeat
  }

  list(method = "kit_stratified_rclr",
       features = colnames(X),
       kit_label_col = kit_label_col,
       biofluid_col = biofluid_col,
       global_center = global_center,
       stratum_centers = centers,
       trim_upper = trim_upper,
       trim_lower = trim_lower,
       min_centering_size = min_centering_size,
       min_samples_per_stratum = min_samples_per_stratum,
       exclude_features = exclude_features)
}

#' Score kit-stratified rCLR.
#'
#' @param expr_matrix Numeric matrix, samples x features.
#' @param sample_meta data.frame with one row per sample and `kit_label_col`.
#' @param panel_features Fold-specific panel features, typically top 20.
#' @param kit_label_col Column containing kit family. NA, "unknown", "NR", and
#'   non-NGS sentinels fall back to biofluid, then global centering.
#' @param fit Optional object from `fit_kit_stratified_rclr()`. If omitted, the
#'   fit is estimated from `expr_matrix`, which is convenient for deployment but
#'   should not be used for held-out benchmarking.
#' @param feature_weights Optional named numeric signs/weights for panel
#'   features; default all +1. This lets training folds orient features so
#'   higher scores are more case-like.
#' @param biofluid_col Optional fallback stratum column.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param ... Additional arguments passed to the fitting helper when `fit` is
#'   not supplied.
#' @return Numeric score vector of length nrow(expr_matrix).
#' @export
score_kit_stratified_rclr <- function(expr_matrix,
                                      sample_meta,
                                      panel_features,
                                      kit_label_col,
                                      biofluid_col = "biofluid",
                                      fit = NULL,
                                      feature_weights = NULL,
                                      pseudocount = NULL,
                                      ...) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col, panel_features)
  if (is.null(fit)) {
    fit <- fit_kit_stratified_rclr(inp$X, inp$sample_meta, kit_label_col,
                                   biofluid_col = biofluid_col, ...)
  }
  pw <- .ka_panel_and_weights(inp$panel_features, feature_weights, colnames(inp$X))
  all_centers <- unique(c(fit$global_center, unlist(fit$stratum_centers, use.names = FALSE)))
  needed <- unique(c(pw$panel, all_centers))
  logp <- .ka_log_proportions(inp$X, pseudocount = pseudocount, features = needed)
  strata <- .ka_resolve_strata(inp$sample_meta, kit_label_col, fit$biofluid_col %||% biofluid_col)
  score <- rep(NA_real_, nrow(inp$X))
  for (s in unique(strata$label)) {
    idx <- which(strata$label == s)
    center <- fit$stratum_centers[[s]] %||% fit$global_center
    center <- intersect(center, colnames(logp))
    if (length(center) == 0L) center <- setdiff(colnames(logp), pw$panel)
    if (length(center) == 0L) center <- colnames(logp)
    ctr <- rowMeans(logp[idx, center, drop = FALSE], na.rm = TRUE)
    z_panel <- sweep(logp[idx, pw$panel, drop = FALSE], 1L, ctr, `-`)
    score[idx] <- .ka_weighted_panel_sum(z_panel, pw$panel, pw$weights)
  }
  score
}

#' Fit kit fixed-effect adjusted ALR.
#'
#' @description Per feature, log-CPM values are adjusted by subtracting the
#'   training stratum mean and adding the global training mean. ALR pivots are
#'   selected from features with the lowest between-kit variance. When kit is
#'   unknown, biofluid and then global strata are used.
#' @inheritParams fit_kit_stratified_rclr
#' @param n_pivots Number of low-variance pivot features to select.
#' @param pseudocount Additive pseudocount before log-CPM transformation.
#' @return Fit object consumed by `score_kit_fe_adjusted_alr()`.
#' @export
fit_kit_fe_adjusted_alr <- function(expr_matrix,
                                    sample_meta,
                                    kit_label_col,
                                    biofluid_col = "biofluid",
                                    n_pivots = 6L,
                                    min_samples_per_stratum = 5L,
                                    pseudocount = 0.5,
                                    exclude_features = c("hsa-miR-451a",
                                                         "hsa-miR-16-5p",
                                                         "hsa-miR-486-5p",
                                                         "hsa-miR-144-3p",
                                                         "hsa-miR-223-3p")) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col)
  X <- inp$X
  strata <- .ka_resolve_strata(inp$sample_meta, kit_label_col, biofluid_col)
  logcpm <- .ka_log_cpm(X, pseudocount = pseudocount)
  global_mean <- colMeans(logcpm, na.rm = TRUE)
  overall_var <- apply(logcpm, 2L, stats::var, na.rm = TRUE)

  means <- list()
  for (s in sort(unique(strata$label))) {
    idx <- which(strata$label == s)
    if (length(idx) < min_samples_per_stratum) next
    means[[s]] <- colMeans(logcpm[idx, , drop = FALSE], na.rm = TRUE)
  }
  if (length(means) >= 2L) {
    mean_mat <- do.call(rbind, means)
    kit_var <- apply(mean_mat, 2L, stats::var, na.rm = TRUE)
  } else {
    kit_var <- overall_var
  }
  kit_var[!is.finite(kit_var)] <- Inf
  pivot_pool <- setdiff(names(sort(kit_var, decreasing = FALSE)), exclude_features)
  pivots <- head(pivot_pool, min(n_pivots, length(pivot_pool)))
  if (length(pivots) == 0L) pivots <- head(names(sort(overall_var)), min(n_pivots, ncol(X)))

  list(method = "kit_fe_adjusted_alr",
       features = colnames(X),
       kit_label_col = kit_label_col,
       biofluid_col = biofluid_col,
       global_mean = global_mean,
       stratum_means = means,
       pivot_features = pivots,
       kit_variance = kit_var,
       n_pivots = length(pivots),
       min_samples_per_stratum = min_samples_per_stratum,
       pseudocount = pseudocount,
       exclude_features = exclude_features)
}

#' Score kit fixed-effect adjusted ALR.
#'
#' @inheritParams score_kit_stratified_rclr
#' @param fit Optional object from `fit_kit_fe_adjusted_alr()`.
#' @return Numeric score vector of length nrow(expr_matrix), higher by default
#'   indicating a case-like signed panel sum.
#' @export
score_kit_fe_adjusted_alr <- function(expr_matrix,
                                      sample_meta,
                                      panel_features,
                                      kit_label_col,
                                      biofluid_col = "biofluid",
                                      fit = NULL,
                                      feature_weights = NULL,
                                      pseudocount = NULL,
                                      ...) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col, panel_features)
  if (is.null(fit)) {
    fit <- fit_kit_fe_adjusted_alr(inp$X, inp$sample_meta, kit_label_col,
                                   biofluid_col = biofluid_col, ...)
  }
  pc <- pseudocount %||% fit$pseudocount %||% 0.5
  pw <- .ka_panel_and_weights(inp$panel_features, feature_weights, colnames(inp$X))
  pivots <- intersect(fit$pivot_features, colnames(inp$X))
  if (length(pivots) == 0L) pivots <- setdiff(colnames(inp$X), pw$panel)
  if (length(pivots) == 0L) pivots <- pw$panel
  needed <- unique(c(pw$panel, pivots))
  logcpm <- .ka_log_cpm(inp$X, pseudocount = pc, features = needed)
  strata <- .ka_resolve_strata(inp$sample_meta, kit_label_col, fit$biofluid_col %||% biofluid_col)
  score <- rep(NA_real_, nrow(inp$X))
  for (s in unique(strata$label)) {
    idx <- which(strata$label == s)
    sm <- fit$stratum_means[[s]]
    adj <- logcpm[idx, needed, drop = FALSE]
    common <- intersect(colnames(adj), names(fit$global_mean))
    if (!is.null(sm)) {
      common <- intersect(common, names(sm))
      adj[, common] <- sweep(adj[, common, drop = FALSE], 2L, sm[common], `-`)
      adj[, common] <- sweep(adj[, common, drop = FALSE], 2L, fit$global_mean[common], `+`)
    }
    pivot_now <- intersect(pivots, colnames(adj))
    pivot_ctr <- rowMeans(adj[, pivot_now, drop = FALSE], na.rm = TRUE)
    z_panel <- sweep(adj[, pw$panel, drop = FALSE], 1L, pivot_ctr, `-`)
    score[idx] <- .ka_weighted_panel_sum(z_panel, pw$panel, pw$weights)
  }
  score
}

#' Fit kit-orthogonal ILR-like residualization.
#'
#' @description The training CLR feature space is used to estimate the dominant
#'   kit axis between the two most populated known strata. Scores are computed
#'   after removing the projection onto that axis. If fewer than two strata are
#'   available, scoring falls back to the global CLR panel sum and the fit status
#'   records the limitation.
#' @inheritParams fit_kit_stratified_rclr
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Fit object consumed by `score_kit_orthogonal_ilr()`.
#' @export
fit_kit_orthogonal_ilr <- function(expr_matrix,
                                   sample_meta,
                                   kit_label_col,
                                   biofluid_col = "biofluid",
                                   min_samples_per_stratum = 5L,
                                   pseudocount = NULL) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col)
  X <- inp$X
  strata <- .ka_resolve_strata(inp$sample_meta, kit_label_col, biofluid_col)
  logp <- .ka_log_proportions(X, pseudocount = pseudocount)
  clr <- sweep(logp, 1L, rowMeans(logp, na.rm = TRUE), `-`)
  counts <- sort(table(strata$label), decreasing = TRUE)
  counts <- counts[counts >= min_samples_per_stratum]
  axis <- NULL
  axis_strata <- character(0)
  status <- "fallback_global_no_kit_axis"
  if (length(counts) >= 2L) {
    axis_strata <- names(counts)[seq_len(2L)]
    m1 <- colMeans(clr[strata$label == axis_strata[1L], , drop = FALSE], na.rm = TRUE)
    m2 <- colMeans(clr[strata$label == axis_strata[2L], , drop = FALSE], na.rm = TRUE)
    raw_axis <- m1 - m2
    raw_axis <- raw_axis - mean(raw_axis, na.rm = TRUE)
    norm <- sqrt(sum(raw_axis^2, na.rm = TRUE))
    if (is.finite(norm) && norm > 1e-10) {
      axis <- raw_axis / norm
      status <- "evaluable_kit_axis_removed"
    }
  }
  list(method = "kit_orthogonal_ilr",
       features = colnames(X),
       kit_label_col = kit_label_col,
       biofluid_col = biofluid_col,
       kit_axis = axis,
       axis_strata = axis_strata,
       status = status,
       min_samples_per_stratum = min_samples_per_stratum,
       pseudocount = pseudocount)
}

#' Score kit-orthogonal ILR-like residualized panel.
#'
#' @inheritParams score_kit_stratified_rclr
#' @param fit Optional object from `fit_kit_orthogonal_ilr()`.
#' @return Numeric score vector of length nrow(expr_matrix).
#' @export
score_kit_orthogonal_ilr <- function(expr_matrix,
                                     sample_meta,
                                     panel_features,
                                     kit_label_col,
                                     biofluid_col = "biofluid",
                                     fit = NULL,
                                     feature_weights = NULL,
                                     pseudocount = NULL,
                                     ...) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col, panel_features)
  if (is.null(fit)) {
    fit <- fit_kit_orthogonal_ilr(inp$X, inp$sample_meta, kit_label_col,
                                  biofluid_col = biofluid_col, ...)
  }
  pc <- pseudocount %||% fit$pseudocount
  pw <- .ka_panel_and_weights(inp$panel_features, feature_weights, colnames(inp$X))
  needed <- if (is.null(fit$kit_axis)) pw$panel else unique(c(pw$panel, names(fit$kit_axis)))
  logp <- .ka_log_proportions(inp$X, pseudocount = pc, features = needed)
  clr <- sweep(logp, 1L, rowMeans(logp, na.rm = TRUE), `-`)
  raw_score <- .ka_weighted_panel_sum(clr, pw$panel, pw$weights)
  if (is.null(fit$kit_axis)) return(raw_score)
  common_axis <- intersect(names(fit$kit_axis), colnames(clr))
  if (length(common_axis) == 0L) return(raw_score)
  axis_vec <- fit$kit_axis[common_axis]
  denom <- sum(axis_vec^2, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(raw_score)
  projection <- as.numeric(clr[, common_axis, drop = FALSE] %*% axis_vec) / denom
  panel_axis <- intersect(pw$panel, common_axis)
  correction_weight <- if (length(panel_axis) == 0L) 0 else {
    sum(fit$kit_axis[panel_axis] * pw$weights[panel_axis], na.rm = TRUE)
  }
  raw_score - projection * correction_weight
}

#' Fit kit-residual MAD scoring.
#'
#' @description Per feature, log-CPM values are robustly centered and scaled
#'   within kit/biofluid strata using median and MAD. Unknown kit falls back to
#'   biofluid, then global. This optional method is included as a negative or
#'   robustness comparator.
#' @inheritParams fit_kit_stratified_rclr
#' @param pseudocount Additive pseudocount before log-CPM transformation.
#' @param mad_floor Minimum robust scale used to avoid division by zero.
#' @return Fit object consumed by `score_kit_residual_mad()`.
#' @export
fit_kit_residual_mad <- function(expr_matrix,
                                 sample_meta,
                                 kit_label_col,
                                 biofluid_col = "biofluid",
                                 min_samples_per_stratum = 5L,
                                 pseudocount = 0.5,
                                 mad_floor = 1e-3) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col)
  X <- inp$X
  strata <- .ka_resolve_strata(inp$sample_meta, kit_label_col, biofluid_col)
  logcpm <- .ka_log_cpm(X, pseudocount = pseudocount)
  global_median <- apply(logcpm, 2L, stats::median, na.rm = TRUE)
  global_mad <- apply(logcpm, 2L, stats::mad, na.rm = TRUE)
  global_mad[!is.finite(global_mad) | global_mad < mad_floor] <- mad_floor
  medians <- list()
  mads <- list()
  for (s in sort(unique(strata$label))) {
    idx <- which(strata$label == s)
    if (length(idx) < min_samples_per_stratum) next
    medians[[s]] <- apply(logcpm[idx, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
    m <- apply(logcpm[idx, , drop = FALSE], 2L, stats::mad, na.rm = TRUE)
    m[!is.finite(m) | m < mad_floor] <- mad_floor
    mads[[s]] <- m
  }
  list(method = "kit_residual_mad",
       features = colnames(X),
       kit_label_col = kit_label_col,
       biofluid_col = biofluid_col,
       global_median = global_median,
       global_mad = global_mad,
       stratum_medians = medians,
       stratum_mads = mads,
       min_samples_per_stratum = min_samples_per_stratum,
       pseudocount = pseudocount,
       mad_floor = mad_floor)
}

#' Score kit-residual MAD.
#'
#' @inheritParams score_kit_stratified_rclr
#' @param fit Optional object from `fit_kit_residual_mad()`.
#' @return Numeric score vector of length nrow(expr_matrix).
#' @export
score_kit_residual_mad <- function(expr_matrix,
                                   sample_meta,
                                   panel_features,
                                   kit_label_col,
                                   biofluid_col = "biofluid",
                                   fit = NULL,
                                   feature_weights = NULL,
                                   pseudocount = NULL,
                                   ...) {
  inp <- .ka_validate_inputs(expr_matrix, sample_meta, kit_label_col, panel_features)
  if (is.null(fit)) {
    fit <- fit_kit_residual_mad(inp$X, inp$sample_meta, kit_label_col,
                                biofluid_col = biofluid_col, ...)
  }
  pc <- pseudocount %||% fit$pseudocount %||% 0.5
  pw <- .ka_panel_and_weights(inp$panel_features, feature_weights, colnames(inp$X))
  logcpm <- .ka_log_cpm(inp$X, pseudocount = pc, features = pw$panel)
  strata <- .ka_resolve_strata(inp$sample_meta, kit_label_col, fit$biofluid_col %||% biofluid_col)
  score <- rep(NA_real_, nrow(inp$X))
  for (s in unique(strata$label)) {
    idx <- which(strata$label == s)
    med <- fit$stratum_medians[[s]] %||% fit$global_median
    mad <- fit$stratum_mads[[s]] %||% fit$global_mad
    common <- intersect(pw$panel, colnames(logcpm))
    z <- sweep(logcpm[idx, common, drop = FALSE], 2L, med[common], `-`)
    z <- sweep(z, 2L, mad[common], `/`)
    score[idx] <- .ka_weighted_panel_sum(z, common, pw$weights)
  }
  score
}

# ---- S45 benchmark runner ----------------------------------------------------

.s45_paths <- list(
  table = "results/reviewer_package_v5/supplement/table_S45_kit_aware_compositional_benchmarks.tsv",
  png = "results/reviewer_package_v5/figures/figure_S45_kit_aware_methods.png",
  pdf = "results/reviewer_package_v5/figures/figure_S45_kit_aware_methods.pdf",
  caption = "results/reviewer_package_v5/figures/figure_S45_caption.md",
  log = file.path("results/reviewer_package_v5/delegated_tasks",
                  paste0("worker_E1_kit_aware_compositional_",
                         format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
                         ".log"))
)

.s45_config <- list(seed = 42L, outer_k = 5L, panel_size = 20L,
                    min_valid_folds = 4L, min_pos_per_fold = 5L,
                    min_neg_per_fold = 5L, bootstrap_B = 1000L)

.s45_map_modality <- function(x) {
  x <- trimws(as.character(x))
  data.table::fifelse(grepl("^NGS$|miRNA-seq", x, ignore.case = TRUE), "NGS",
    data.table::fifelse(grepl("Toray|Microarray", x, ignore.case = TRUE), "Toray",
      data.table::fifelse(grepl("NanoString|array-other", x, ignore.case = TRUE), "NanoString",
        data.table::fifelse(grepl("qPCR|FirePlex", x, ignore.case = TRUE), "qPCR", x))))
}

.s45_map_biofluid <- function(anticoagulant) {
  raw <- tolower(.ka_clean_label(anticoagulant))
  out <- rep("unknown", length(raw))
  out[grepl("serum|clot", raw)] <- "serum-clot"
  out[grepl("edta", raw)] <- "plasma-EDTA"
  out[grepl("exosome|\\bev\\b|sev|extracellular vesicle", raw)] <- "exosome"
  out[grepl("plasma", raw) & out == "unknown"] <- "plasma-other"
  out[is.na(raw) | raw %in% c("nr", "paper_not_found", "not reported")] <- "unknown"
  out
}

.s45_map_extraction <- function(extraction_kit) {
  raw <- tolower(.ka_clean_label(extraction_kit))
  out <- rep("NR", length(raw))
  out[grepl("mirneasy|exorneasy|norgen|column|silica", raw)] <- "silica_spin"
  out[grepl("exoquick|precip", raw)] <- "precipitation"
  out[grepl("trizol|mirvana|phenol|chloroform", raw)] <- "phenol_chloroform"
  out[is.na(raw) | raw %in% c("nr", "paper_not_found", "not reported")] <- "NR"
  out
}

.s45_map_library_kit <- function(protocol, modality) {
  p <- tolower(.ka_clean_label(protocol))
  out <- rep("not_NGS", length(p))
  is_ngs <- modality == "NGS"
  out[is_ngs] <- "unknown"
  ok <- is_ngs & p %in% c("custom", "illumina", "nextflex", "qiaseq")
  out[ok] <- p[ok]
  out
}

.s45_primary_cohorts <- function() {
  s5 <- data.table::fread("results/reviewer_package_v5/supplement/table_S5_null_calibration.tsv")
  s22 <- data.table::fread("results/reviewer_package_v5/supplement/table_S22_delong_paired_auc.tsv")
  c5 <- sort(unique(s5$Accession))
  c22 <- sort(unique(s22[status == "evaluable", cohort]))
  if (!identical(c5, c22)) {
    stop("Primary cohort set mismatch between S5 and S22.")
  }
  c5
}

.s45_make_metadata <- function(primary) {
  fig4 <- readRDS("analysis/figures/figure_4_data.rds")
  pc <- data.table::as.data.table(fig4$per_cell)
  fig_meta <- unique(pc[accession %in% primary,
                        .(cohort = accession,
                          cancer_type = cancer_type,
                          modality = .s45_map_modality(modality_family),
                          provenance_block)])
  inv <- data.table::fread("knowledge/union_inventory.tsv", sep = "\t", quote = "")
  inv_meta <- inv[accession %in% primary,
                  .(cohort = accession,
                    present_in, track,
                    modality_raw = modality,
                    body_fluid_inventory = body_fluid)]
  pre <- data.table::fread("knowledge/cohort_preanalytics_metadata.tsv", sep = "\t", quote = "")
  pre_meta <- pre[accession %in% primary,
                  .(cohort = accession,
                    anticoagulant_raw = anticoagulant,
                    extraction_kit_raw = extraction_kit)]
  kit <- data.table::fread("code/smrnaseq/kit_assignments.tsv", sep = "\t", quote = "")
  kit_meta <- kit[, .(cohort = accession, kit_product, smrnaseq_protocol)]
  meta <- Reduce(function(a, b) merge(a, b, by = "cohort", all.x = TRUE, sort = FALSE),
                 list(fig_meta, inv_meta, pre_meta, kit_meta))
  meta[, biofluid := .s45_map_biofluid(anticoagulant_raw)]
  meta[, extraction_kit := .s45_map_extraction(extraction_kit_raw)]
  meta[, library_kit := .s45_map_library_kit(smrnaseq_protocol, modality)]
  meta[, kit_family_for_scoring := data.table::fifelse(
    modality == "NGS" & !(library_kit %in% c("unknown", "not_NGS")),
    library_kit, NA_character_)]
  meta[]
}

.s45_load_data <- function(meta, log_fun) {
  stop("The reviewer-package S45 runner is not available from the installed package; ",
       "use the manuscript repository analysis script instead.", call. = FALSE)
  prep_env <- new.env(parent = globalenv())
  inv <- data.table::fread("knowledge/union_inventory.tsv", sep = "\t", quote = "")
  rows <- inv[accession %in% meta$cohort]
  rows <- rows[match(meta$cohort, accession)]
  sample_filters <- data.table::fread("knowledge/v0_6_sample_filters.tsv", sep = "\t",
                                      quote = "", colClasses = "character")
  get_fun <- function(name) {
    if (exists(name, envir = prep_env, inherits = FALSE)) {
      get(name, envir = prep_env, inherits = FALSE)
    } else {
      get(name, envir = globalenv(), inherits = TRUE)
    }
  }
  outcome_dict <- get_fun("load_outcome_dictionary_v0_4")(strict = TRUE)
  matched_set_audit <- get_fun("load_matched_set_audit")(
    "knowledge/matched_set_audit.tsv", n_tier_a = 35L)
  out <- list()
  failures <- character()
  for (i in seq_len(nrow(rows))) {
    acc <- rows$accession[i]
    log_fun(sprintf("load %02d/%02d %s", i, nrow(rows), acc))
    dat <- tryCatch(
      prep_env$.prepare_cohort_data(rows[i], outcome_dict, sample_filters,
                                    matched_set_audit),
      error = function(e) e)
    if (inherits(dat, "error")) {
      failures <- c(failures, sprintf("%s: %s", acc, dat$message))
      next
    }
    out[[acc]] <- dat
  }
  if (length(failures) > 0L) {
    stop("Cohort load failures: ", paste(failures, collapse = " | "))
  }
  out
}

.s45_sample_meta <- function(dat, meta_row) {
  data.frame(
    sample_id = rownames(dat$X) %||% paste0("sample_", seq_len(nrow(dat$X))),
    cohort = meta_row$cohort,
    biofluid = meta_row$biofluid,
    library_kit = meta_row$library_kit,
    extraction_kit = meta_row$extraction_kit,
    kit_family = meta_row$kit_family_for_scoring,
    global_label = "global",
    stringsAsFactors = FALSE
  )
}

.s45_fit_for_method <- function(method, X_train, meta_train) {
  switch(method,
    kit_stratified_rclr = fit_kit_stratified_rclr(X_train, meta_train, "kit_family"),
    kit_fe_adjusted_alr = fit_kit_fe_adjusted_alr(X_train, meta_train, "kit_family"),
    kit_orthogonal_ilr = fit_kit_orthogonal_ilr(X_train, meta_train, "kit_family"),
    kit_residual_mad = fit_kit_residual_mad(X_train, meta_train, "kit_family"),
    stop("Unknown method: ", method)
  )
}

.s45_score_for_method <- function(method, X_test, meta_test, panel, weights, fit) {
  switch(method,
    kit_stratified_rclr = score_kit_stratified_rclr(
      X_test, meta_test, panel, "kit_family", fit = fit, feature_weights = weights),
    kit_fe_adjusted_alr = score_kit_fe_adjusted_alr(
      X_test, meta_test, panel, "kit_family", fit = fit, feature_weights = weights),
    kit_orthogonal_ilr = score_kit_orthogonal_ilr(
      X_test, meta_test, panel, "kit_family", fit = fit, feature_weights = weights),
    kit_residual_mad = score_kit_residual_mad(
      X_test, meta_test, panel, "kit_family", fit = fit, feature_weights = weights),
    stop("Unknown method: ", method)
  )
}

.s45_run_cohort <- function(acc, dat, meta_row, methods, log_fun) {
  y <- as.integer(dat$y)
  X <- as.matrix(dat$X)
  smeta <- .s45_sample_meta(dat, meta_row)
  folds <- if (!is.null(dat$group_id) && length(dat$group_id) == nrow(X)) {
    .grouped_folds(y, dat$group_id, k = .s45_config$outer_k, seed = .s45_config$seed)
  } else {
    .stratified_folds(y, k = .s45_config$outer_k, seed = .s45_config$seed)
  }
  scores <- as.data.frame(matrix(NA_real_, nrow = nrow(X), ncol = length(methods) + 1L))
  names(scores) <- c("rCLR_baseline", methods)
  valid_fold <- logical(length(folds))
  fold_notes <- character(length(folds))
  for (fi in seq_along(folds)) {
    test_idx <- folds[[fi]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    y_train <- y[train_idx]
    y_test <- y[test_idx]
    if (sum(y_test == 1L) < .s45_config$min_pos_per_fold ||
        sum(y_test == 0L) < .s45_config$min_neg_per_fold ||
        length(unique(y_train)) < 2L) {
      fold_notes[fi] <- "class_floor_failed"
      next
    }
    X_train <- X[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]
    meta_train <- smeta[train_idx, , drop = FALSE]
    meta_test <- smeta[test_idx, , drop = FALSE]
    fit_base <- fit_kit_stratified_rclr(X_train, meta_train, "global_label")
    z_train <- .ka_global_rclr_matrix(X_train, fit_base$global_center)
    weights <- .ka_select_panel_from_matrix(z_train, y_train, panel_size = .s45_config$panel_size)
    panel <- names(weights)
    scores$rCLR_baseline[test_idx] <- score_kit_stratified_rclr(
      X_test, meta_test, panel, "global_label", fit = fit_base, feature_weights = weights)
    for (method in methods) {
      fit_m <- .s45_fit_for_method(method, X_train, meta_train)
      scores[[method]][test_idx] <- .s45_score_for_method(method, X_test, meta_test,
                                                          panel, weights, fit_m)
    }
    valid_fold[fi] <- TRUE
    fold_notes[fi] <- paste0("ok_panel=", paste(panel, collapse = ","))
  }
  if (sum(valid_fold) < .s45_config$min_valid_folds) {
    log_fun(sprintf("%s not evaluable: n_valid_folds=%d", acc, sum(valid_fold)))
  }
  rows <- list()
  for (method in methods) {
    auc_m <- .ka_fast_auc(y, scores[[method]])
    auc_b <- .ka_fast_auc(y, scores$rCLR_baseline)
    n_scored <- sum(is.finite(scores[[method]]) & is.finite(scores$rCLR_baseline))
    n_pos <- sum(y == 1L & is.finite(scores[[method]]) & is.finite(scores$rCLR_baseline))
    n_neg <- sum(y == 0L & is.finite(scores[[method]]) & is.finite(scores$rCLR_baseline))
    status <- if (sum(valid_fold) >= .s45_config$min_valid_folds &&
                  is.finite(auc_m) && is.finite(auc_b)) "evaluable" else "not_evaluable"
    rows[[method]] <- data.table::data.table(
      cohort = acc,
      method = method,
      auc = auc_m,
      var_auc = .ka_hm_var(auc_m, n_pos, n_neg),
      rclr_baseline_auc = auc_b,
      var_rclr_baseline = .ka_hm_var(auc_b, n_pos, n_neg),
      lift = auc_m - auc_b,
      var_lift = .ka_hm_var(auc_m, n_pos, n_neg) + .ka_hm_var(auc_b, n_pos, n_neg),
      n_samples = n_scored,
      n_pos = n_pos,
      n_neg = n_neg,
      n_valid_folds = sum(valid_fold),
      status = status,
      fold_notes = paste(fold_notes, collapse = " | "))
  }
  data.table::rbindlist(rows, fill = TRUE)
}

.s45_pool_metric <- function(dt, metric_col, var_col) {
  dt <- dt[is.finite(get(metric_col)) & is.finite(get(var_col)) & get(var_col) > 0]
  if (nrow(dt) == 0L) return(list(estimate = NA_real_, tau2 = NA_real_, I2 = NA_real_))
  block <- dt[, {
    w <- 1 / pmax(get(var_col), 1e-8)
    .(yi = sum(w * get(metric_col)) / sum(w),
      vi = 1 / sum(w),
      n_cohorts = .N,
      n_samples = sum(n_samples, na.rm = TRUE))
  }, by = provenance_block]
  if (nrow(block) == 1L) {
    return(list(estimate = block$yi[1L], tau2 = NA_real_, I2 = NA_real_))
  }
  fit <- tryCatch(
    metafor::rma.uni(yi = block$yi, vi = block$vi, method = "REML"),
    error = function(e) e)
  if (inherits(fit, "error")) {
    w <- 1 / block$vi
    return(list(estimate = sum(w * block$yi) / sum(w), tau2 = NA_real_, I2 = NA_real_))
  }
  list(estimate = as.numeric(fit$b[1L]), tau2 = as.numeric(fit$tau2), I2 = as.numeric(fit$I2))
}

.s45_boot_ci <- function(dt, metric_col, var_col, B = .s45_config$bootstrap_B) {
  dt <- dt[is.finite(get(metric_col)) & is.finite(get(var_col)) & get(var_col) > 0]
  if (nrow(dt) == 0L) return(c(NA_real_, NA_real_))
  if (nrow(dt) == 1L) return(rep(dt[[metric_col]][1L], 2L))
  vals <- replicate(B, {
    idx <- sample.int(nrow(dt), replace = TRUE)
    .s45_pool_metric(dt[idx], metric_col, var_col)$estimate
  })
  stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

.s45_summary_rows <- function(cohort_dt) {
  axes <- list(
    pooled = cohort_dt[, .(axis = "pooled", stratum = "overall",
                           method, cohort, provenance_block, n_samples,
                           auc, var_auc, rclr_baseline_auc,
                           var_rclr_baseline, lift, var_lift, status)],
    biofluid = cohort_dt[, .(axis = "biofluid", stratum = biofluid,
                             method, cohort, provenance_block, n_samples,
                             auc, var_auc, rclr_baseline_auc,
                             var_rclr_baseline, lift, var_lift, status)],
    library_kit = cohort_dt[modality == "NGS",
                            .(axis = "library_kit", stratum = library_kit,
                              method, cohort, provenance_block, n_samples,
                              auc, var_auc, rclr_baseline_auc,
                              var_rclr_baseline, lift, var_lift, status)],
    extraction_kit = cohort_dt[, .(axis = "extraction_kit", stratum = extraction_kit,
                                   method, cohort, provenance_block, n_samples,
                                   auc, var_auc, rclr_baseline_auc,
                                   var_rclr_baseline, lift, var_lift, status)]
  )
  long <- data.table::rbindlist(axes, fill = TRUE)
  long <- long[status == "evaluable" & !is.na(stratum)]
  out <- long[, {
    auc_pool <- .s45_pool_metric(.SD, "auc", "var_auc")
    lift_pool <- .s45_pool_metric(.SD, "lift", "var_lift")
    auc_ci <- .s45_boot_ci(.SD, "auc", "var_auc")
    lift_ci <- .s45_boot_ci(.SD, "lift", "var_lift")
    n_blocks <- data.table::uniqueN(provenance_block)
    .(n_cohorts = data.table::uniqueN(cohort),
      n_samples = sum(unique(.SD[, .(cohort, n_samples)])$n_samples, na.rm = TRUE),
      AUC_re_pooled = auc_pool$estimate,
      CI_lo = max(0, auc_ci[1L]),
      CI_hi = min(1, auc_ci[2L]),
      tau2 = auc_pool$tau2,
      I2 = auc_pool$I2,
      lift_vs_rCLR_baseline = lift_pool$estimate,
      lift_ci_lo = lift_ci[1L],
      lift_ci_hi = lift_ci[2L],
      status = if (n_blocks < 2L) "descriptive_single_provenance_block" else "evaluable_REML_block_collapsed_cluster_bootstrap")
  }, by = .(method, axis, stratum)]

  flagged_methods <- out[axis != "pooled" & AUC_re_pooled > 0.95,
                         .N, by = method][N >= 2L, method]
  out[method %in% flagged_methods,
      status := paste(status, "plausibility_flag_auc_gt_0.95_multiple_strata", sep = ";")]
  data.table::setorder(out, method, axis, stratum)
  out[]
}

.s45_make_figure <- function(summary_dt) {
  plot_dt <- summary_dt[axis %in% c("pooled", "biofluid", "library_kit", "extraction_kit")]
  plot_dt[, method_label := factor(method,
    levels = c("kit_stratified_rclr", "kit_fe_adjusted_alr",
               "kit_orthogonal_ilr", "kit_residual_mad"),
    labels = c("kit-stratified rCLR", "kit-FE ALR",
               "kit-orthogonal ILR", "kit-residual MAD"))]
  plot_dt[, axis_label := factor(axis,
    levels = c("pooled", "biofluid", "library_kit", "extraction_kit"),
    labels = c("Overall", "Biofluid", "NGS library kit", "Extraction kit"))]
  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = lift_vs_rCLR_baseline, y = stratum,
                 xmin = lift_ci_lo, xmax = lift_ci_hi,
                 colour = method_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.58),
                             linewidth = 0.35, na.rm = TRUE) +
    ggplot2::facet_grid(axis_label ~ method_label, scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "Lift vs same-protocol rCLR baseline (AUC)",
                  y = NULL, colour = "Method") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   strip.text.x = ggplot2::element_text(face = "bold"),
                   strip.text.y = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(.s45_paths$png, p, width = 11, height = 8.5, dpi = 300, bg = "white")
  ggplot2::ggsave(.s45_paths$pdf, p, width = 11, height = 8.5, bg = "white")
}

.s45_write_caption <- function(summary_dt) {
  robust <- summary_dt[!grepl("descriptive_single_provenance_block", status)]
  positive <- robust[lift_ci_lo > 0 & is.finite(lift_ci_lo)]
  negative <- robust[lift_ci_hi < 0 & is.finite(lift_ci_hi)]
  descriptive_positive <- summary_dt[grepl("descriptive_single_provenance_block", status) &
                                       lift_ci_lo > 0 & is.finite(lift_ci_lo)]
  top <- positive[order(-lift_vs_rCLR_baseline)][1L]
  pos_text <- if (nrow(positive) == 0L) {
    "No multi-block kit-aware method-stratum row had a bootstrap 95% lift interval entirely above zero in this provisional run."
  } else {
    paste0(nrow(positive), " multi-block method-stratum row(s) had lift intervals above zero; the largest was ",
           top$method, " in ", top$axis, "=", top$stratum,
           sprintf(" (lift %.3f, 95%% CI %.3f to %.3f).",
                   top$lift_vs_rCLR_baseline, top$lift_ci_lo, top$lift_ci_hi))
  }
  neg_text <- if (nrow(negative) == 0L) {
    "No method-stratum row had a bootstrap 95% lift interval entirely below zero."
  } else {
    paste0(nrow(negative), " method-stratum row(s) had lift intervals below zero and should be reported as negative findings.")
  }
  desc_text <- if (nrow(descriptive_positive) == 0L) {
    "Single-provenance-block descriptive rows did not add positive lift claims."
  } else {
    paste0(nrow(descriptive_positive),
           " single-provenance-block row(s) were positive but are descriptive only and should not be treated as robust lift evidence.")
  }
  leak <- summary_dt[grepl("plausibility_flag", status)]
  leak_text <- if (nrow(leak) == 0L) {
    "Plausibility gate: no method exceeded pooled AUC 0.95 in multiple strata."
  } else {
    paste0("Plausibility gate: ", paste(unique(leak$method), collapse = ", "),
           " exceeded pooled AUC 0.95 in multiple strata; these rows need targeted leakage review before interpretation.")
  }
  cap <- c(
    "# Figure S45. Kit-aware compositional scoring methods",
    "",
    "Forest-style intervals show lift over a same-protocol rCLR baseline for the 22 primary matched-null evaluable cohorts. Each cohort was scored by nested folds with training-fold top-20 rCLR panel selection; the new methods were fitted on training folds only and evaluated on held-out samples. Stratum summaries use REML random-effects pooling after provenance-block collapse, and intervals are cohort-resampling cluster bootstraps with block collapse repeated in each bootstrap.",
    "",
    pos_text,
    neg_text,
    desc_text,
    leak_text,
    "",
    "This Worker E1 artifact is reproducible from `Rscript code/methods/kit_aware_compositional.R` but is not manuscript-accepted. It requires /triple-consensus and /plausibility-check before being used as evidence."
  )
  writeLines(cap, .s45_paths$caption)
}

.s45_main <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor is required.")
  set.seed(.s45_config$seed)
  dir.create(dirname(.s45_paths$table), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(.s45_paths$png), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(.s45_paths$log), recursive = TRUE, showWarnings = FALSE)
  log_lines <- character()
  log_fun <- function(msg) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), msg)
    log_lines <<- c(log_lines, line)
    message(line)
  }

  log_fun("Worker E1 kit-aware compositional benchmark start")
  log_fun("Seed=42; no nested delegation; outputs are provisional pending gates")
  primary <- .s45_primary_cohorts()
  log_fun(sprintf("Primary matched-null universe: %d cohorts", length(primary)))
  meta <- .s45_make_metadata(primary)
  data_list <- .s45_load_data(meta, log_fun)
  methods <- c("kit_stratified_rclr", "kit_fe_adjusted_alr",
               "kit_orthogonal_ilr", "kit_residual_mad")
  cohort_rows <- list()
  for (acc in names(data_list)) {
    log_fun(sprintf("benchmark %s", acc))
    rows <- .s45_run_cohort(acc, data_list[[acc]], meta[cohort == acc][1L],
                            methods, log_fun)
    cohort_rows[[acc]] <- rows
  }
  cohort_dt <- data.table::rbindlist(cohort_rows, fill = TRUE)
  cohort_dt <- merge(cohort_dt,
                     meta[, .(cohort, cancer_type, modality, provenance_block,
                              biofluid, library_kit, extraction_kit)],
                     by = "cohort", all.x = TRUE, sort = FALSE)
  summary_dt <- .s45_summary_rows(cohort_dt)
  data.table::fwrite(summary_dt, .s45_paths$table, sep = "\t", quote = FALSE, na = "")
  .s45_make_figure(summary_dt)
  .s45_write_caption(summary_dt)

  log_fun(sprintf("Wrote table: %s (%d rows)", .s45_paths$table, nrow(summary_dt)))
  log_fun(sprintf("Wrote figure: %s and %s", .s45_paths$png, .s45_paths$pdf))
  log_fun(sprintf("Wrote caption: %s", .s45_paths$caption))
  log_fun("Package export note: this paper repo has no R/ or man/ package layout; downstream OmicSelector package export files were not created in this bounded worker scope.")
  flags <- summary_dt[grepl("plausibility_flag", status)]
  if (nrow(flags) > 0L) {
    log_fun(sprintf("Plausibility flags: %s",
                    paste(unique(paste(flags$method, flags$axis, flags$stratum, sep = "/")),
                          collapse = "; ")))
  } else {
    log_fun("Plausibility flags: none")
  }
  log_fun("Session info follows")
  log_lines <- c(log_lines, capture.output(sessionInfo()))
  writeLines(log_lines, .s45_paths$log)
  message("Run log: ", .s45_paths$log)
  invisible(summary_dt)
}

if (sys.nframe() == 0L) {
  .s45_main()
}
