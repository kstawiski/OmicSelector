#' Provenance-aware within-sample scoring methods
#'
#' These functions implement scoring methods that explicitly account for cohort
#' and specimen-provenance structure. All functions expect an expression matrix
#' with samples in rows and features in columns. `meta` must have one row per
#' sample in the same order as `expr`.
#'
#' The mixed-effects scorer is split into a training step and a prediction step:
#' disease labels are used only by `fit_mixed_effects_scorer()` on the training
#' pool. `predict_mixed_effects_scorer()` and `score_mixed_effects_scorer()`
#' with a supplied `fit` object do not read disease labels from test metadata.
#'
#' @name provenance_aware_scoring
NULL

.pa_check_expr_meta <- function(expr, meta) {
  expr <- as.matrix(expr)
  if (!is.data.frame(meta)) meta <- as.data.frame(meta)
  if (nrow(expr) != nrow(meta)) {
    stop("meta must have one row per expression sample")
  }
  if (is.null(colnames(expr))) {
    stop("expr must have feature names in colnames(expr)")
  }
  storage.mode(expr) <- "double"
  if (any(expr < 0, na.rm = TRUE)) {
    stop("provenance-aware scorers require non-negative abundance-like input")
  }
  list(expr = expr, meta = meta)
}

.pa_panel <- function(expr, panel) {
  panel <- intersect(as.character(panel), colnames(expr))
  if (length(panel) == 0L) {
    stop("none of the panel features are present in expr")
  }
  panel
}

.pa_row_pseudocount <- function(expr, pseudocount = NULL) {
  if (!is.null(pseudocount)) {
    return(rep(as.numeric(pseudocount), nrow(expr)))
  }
  pc <- 1e-6 * rowSums(expr, na.rm = TRUE)
  pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  pc
}

.pa_log_expr <- function(expr, pseudocount = NULL) {
  pc <- .pa_row_pseudocount(expr, pseudocount)
  log(sweep(expr, 1L, pc, FUN = "+"))
}

.pa_mad <- function(x) {
  out <- stats::mad(x, center = stats::median(x, na.rm = TRUE),
                    constant = 1.4826, na.rm = TRUE)
  if (!is.finite(out) || out <= 0) out <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(out) || out <= 0) out <- 1
  out
}

.pa_group_moments <- function(log_expr, groups, features = colnames(log_expr)) {
  groups <- as.character(groups)
  split_idx <- split(seq_len(nrow(log_expr)), groups)
  lapply(split_idx, function(idx) {
    sub <- log_expr[idx, features, drop = FALSE]
    center <- apply(sub, 2L, stats::median, na.rm = TRUE)
    scale <- apply(sub, 2L, .pa_mad)
    scale[!is.finite(scale) | scale <= 0] <- 1
    list(center = center, scale = scale, n = length(idx))
  })
}

.pa_center_only_moments <- function(log_expr, groups, features = colnames(log_expr)) {
  groups <- as.character(groups)
  split_idx <- split(seq_len(nrow(log_expr)), groups)
  lapply(split_idx, function(idx) {
    sub <- log_expr[idx, features, drop = FALSE]
    center <- apply(sub, 2L, stats::median, na.rm = TRUE)
    list(center = center, n = length(idx))
  })
}

.pa_apply_group_z <- function(log_expr, groups, panel, moments,
                              fallback_moments = NULL,
                              allow_new_group_moments = TRUE) {
  groups <- as.character(groups)
  out <- rep(NA_real_, nrow(log_expr))
  fallback_count <- 0L
  for (g in unique(groups)) {
    idx <- which(groups == g)
    mm <- moments[[g]]
    if (is.null(mm) && allow_new_group_moments) {
      sub <- log_expr[idx, panel, drop = FALSE]
      center <- apply(sub, 2L, stats::median, na.rm = TRUE)
      scale <- apply(sub, 2L, .pa_mad)
      scale[!is.finite(scale) | scale <= 0] <- 1
      mm <- list(center = center, scale = scale, n = length(idx))
    }
    if (is.null(mm)) {
      mm <- fallback_moments
      fallback_count <- fallback_count + length(idx)
    }
    if (is.null(mm)) next
    keep <- intersect(panel, names(mm$center))
    if (length(keep) == 0L) next
    z <- sweep(log_expr[idx, keep, drop = FALSE], 2L, mm$center[keep], "-")
    z <- sweep(z, 2L, mm$scale[keep], "/")
    out[idx] <- rowSums(z, na.rm = TRUE)
  }
  attr(out, "fallback_samples") <- fallback_count
  out
}

.pa_apply_block_center <- function(log_expr, groups, panel, moments,
                                   fallback_center = NULL) {
  groups <- as.character(groups)
  out <- rep(NA_real_, nrow(log_expr))
  fallback_count <- 0L
  for (g in unique(groups)) {
    idx <- which(groups == g)
    mm <- moments[[g]]
    center <- if (!is.null(mm)) mm$center else fallback_center
    if (is.null(mm)) fallback_count <- fallback_count + length(idx)
    if (is.null(center)) next
    keep <- intersect(panel, names(center))
    if (length(keep) == 0L) next
    z <- sweep(log_expr[idx, keep, drop = FALSE], 2L, center[keep], "-")
    out[idx] <- rowSums(z, na.rm = TRUE)
  }
  attr(out, "fallback_samples") <- fallback_count
  out
}

#' Fit cohort-wise robust z-score moments
#'
#' @param expr Numeric matrix, samples x features, training or scoring pool.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param cohort_col Column in `meta` identifying cohorts.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return A fit object consumed by `predict_cohort_z_then_pool_score()`.
#' @export
fit_cohort_z_then_pool_score <- function(expr, meta, panel,
                                         cohort_col = "cohort",
                                         pseudocount = NULL) {
  checked <- .pa_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  if (!(cohort_col %in% names(meta))) stop("missing cohort_col in meta")
  panel <- .pa_panel(expr, panel)
  log_expr <- .pa_log_expr(expr, pseudocount)
  moments <- .pa_group_moments(log_expr, meta[[cohort_col]], panel)
  global_center <- apply(log_expr[, panel, drop = FALSE], 2L, stats::median,
                         na.rm = TRUE)
  global_scale <- apply(log_expr[, panel, drop = FALSE], 2L, .pa_mad)
  global_scale[!is.finite(global_scale) | global_scale <= 0] <- 1
  structure(
    list(panel = panel, cohort_col = cohort_col, pseudocount = pseudocount,
         moments = moments,
         fallback = list(center = global_center, scale = global_scale,
                         n = nrow(expr))),
    class = "cohort_z_then_pool_score_fit"
  )
}

#' Score a panel after cohort-wise robust z-standardization
#'
#' Within each cohort, all panel features are log-transformed, centered by the
#' cohort median, scaled by cohort MAD, and summed. When `fit` lacks a new test
#' cohort, the default is to compute unsupervised moments from that test cohort;
#' no disease labels are used for this fallback.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param cohort_col Column in `meta` identifying cohorts.
#' @param block_col Column in `meta` identifying provenance blocks; accepted for
#'   a common scorer signature but not used by this method.
#' @param fit Optional object from `fit_cohort_z_then_pool_score()`.
#' @param allow_new_group_moments Logical; if TRUE, new cohorts are centered on
#'   their own unlabeled samples.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Numeric score vector, one value per sample.
#' @export
score_cohort_z_then_pool_score <- function(expr, meta, panel,
                                           cohort_col = "cohort",
                                           block_col = "provenance_block",
                                           fit = NULL,
                                           allow_new_group_moments = TRUE,
                                           pseudocount = NULL) {
  checked <- .pa_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  if (!(cohort_col %in% names(meta))) stop("missing cohort_col in meta")
  if (is.null(fit)) {
    fit <- fit_cohort_z_then_pool_score(expr, meta, panel, cohort_col,
                                        pseudocount)
  }
  panel <- .pa_panel(expr, fit$panel)
  log_expr <- .pa_log_expr(expr, fit$pseudocount)
  .pa_apply_group_z(log_expr, meta[[cohort_col]], panel, fit$moments,
                    fallback_moments = fit$fallback,
                    allow_new_group_moments = allow_new_group_moments)
}

#' Predict with a fitted cohort-z scorer
#'
#' @inheritParams score_cohort_z_then_pool_score
#' @return Numeric score vector.
#' @export
predict_cohort_z_then_pool_score <- function(fit, expr, meta,
                                             cohort_col = fit$cohort_col,
                                             block_col = "provenance_block",
                                             allow_new_group_moments = TRUE) {
  score_cohort_z_then_pool_score(
    expr = expr, meta = meta, panel = fit$panel, cohort_col = cohort_col,
    block_col = block_col, fit = fit,
    allow_new_group_moments = allow_new_group_moments)
}

#' Fit provenance-block rCLR centering moments
#'
#' Computes one robust log-geometric reference vector per provenance block from
#' a training pool. Prediction subtracts the block-specific reference from each
#' sample's log abundance and sums the panel residuals.
#'
#' @param expr Numeric matrix, samples x features, training pool.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param block_col Column in `meta` identifying provenance blocks.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return A fit object consumed by `predict_block_fe_rclr()`.
#' @export
fit_block_fe_rclr <- function(expr, meta, panel,
                              block_col = "provenance_block",
                              pseudocount = NULL) {
  checked <- .pa_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  if (!(block_col %in% names(meta))) stop("missing block_col in meta")
  panel <- .pa_panel(expr, panel)
  log_expr <- .pa_log_expr(expr, pseudocount)
  moments <- .pa_center_only_moments(log_expr, meta[[block_col]], panel)
  global_center <- apply(log_expr[, panel, drop = FALSE], 2L, stats::median,
                         na.rm = TRUE)
  structure(
    list(panel = panel, block_col = block_col, pseudocount = pseudocount,
         moments = moments, fallback_center = global_center),
    class = "block_fe_rclr_fit"
  )
}

#' Score with provenance-block fixed-effect rCLR centering
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param cohort_col Column in `meta` identifying cohorts; accepted for a common
#'   scorer signature but not used by this method.
#' @param block_col Column in `meta` identifying provenance blocks.
#' @param fit Optional object from `fit_block_fe_rclr()`.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Numeric score vector, one value per sample. Attribute
#'   `fallback_samples` counts samples assigned to the training global center
#'   because their block was absent from the fit.
#' @export
score_block_fe_rclr <- function(expr, meta, panel,
                                cohort_col = "cohort",
                                block_col = "provenance_block",
                                fit = NULL,
                                pseudocount = NULL) {
  checked <- .pa_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  if (!(block_col %in% names(meta))) stop("missing block_col in meta")
  if (is.null(fit)) {
    fit <- fit_block_fe_rclr(expr, meta, panel, block_col, pseudocount)
  }
  panel <- .pa_panel(expr, fit$panel)
  log_expr <- .pa_log_expr(expr, fit$pseudocount)
  .pa_apply_block_center(log_expr, meta[[block_col]], panel, fit$moments,
                         fallback_center = fit$fallback_center)
}

#' Predict with a fitted block rCLR scorer
#'
#' @inheritParams score_block_fe_rclr
#' @return Numeric score vector.
#' @export
predict_block_fe_rclr <- function(fit, expr, meta,
                                  cohort_col = "cohort",
                                  block_col = fit$block_col) {
  score_block_fe_rclr(expr = expr, meta = meta, panel = fit$panel,
                      cohort_col = cohort_col, block_col = block_col,
                      fit = fit)
}

.pa_mixed_formula <- function(has_cohort, has_block) {
  rhs <- "disease"
  if (has_cohort) rhs <- paste(rhs, "+ (1 | cohort)")
  if (has_block) rhs <- paste(rhs, "+ (1 | block)")
  stats::as.formula(paste("value ~", rhs))
}

.pa_lm_formula <- function(has_cohort, has_block) {
  rhs <- "disease"
  if (has_cohort) rhs <- paste(rhs, "+ factor(cohort)")
  if (has_block) rhs <- paste(rhs, "+ factor(block)")
  stats::as.formula(paste("value ~", rhs))
}

.pa_feature_weight <- function(value, disease, cohort, block) {
  df <- data.frame(
    value = as.numeric(value),
    disease = as.numeric(disease),
    cohort = factor(cohort),
    block = factor(block)
  )
  ok <- stats::complete.cases(df)
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 10L || length(unique(df$disease)) < 2L) {
    return(c(weight = NA_real_, model = "not_fit"))
  }
  has_cohort <- nlevels(df$cohort) > 1L
  has_block <- nlevels(df$block) > 1L
  if (requireNamespace("lme4", quietly = TRUE) && (has_cohort || has_block)) {
    fit <- tryCatch(
      lme4::lmer(.pa_mixed_formula(has_cohort, has_block), data = df,
                 REML = FALSE,
                 control = lme4::lmerControl(check.rankX = "silent.drop.cols")),
      error = function(e) NULL)
    if (!is.null(fit)) {
      co <- lme4::fixef(fit)
      if ("disease" %in% names(co) && is.finite(co[["disease"]])) {
        return(c(weight = unname(co[["disease"]]), model = "lmer"))
      }
    }
  }
  fit_lm <- tryCatch(stats::lm(.pa_lm_formula(has_cohort, has_block), data = df),
                     error = function(e) NULL)
  if (is.null(fit_lm)) return(c(weight = NA_real_, model = "not_fit"))
  co <- stats::coef(fit_lm)
  if (!("disease" %in% names(co)) || !is.finite(co[["disease"]])) {
    return(c(weight = NA_real_, model = "not_fit"))
  }
  c(weight = unname(co[["disease"]]), model = "lm_fallback")
}

#' Fit disease weights from a cohort/provenance mixed-effects model
#'
#' Fits one per-feature model on the training pool:
#' `log(feature_count) ~ disease + (1|cohort) + (1|provenance_block)`.
#' Terms with only one level are omitted. If `lme4` cannot fit a feature, a
#' fixed-effect `lm()` fallback is attempted and recorded in the fit audit.
#'
#' @param train_expr Numeric matrix, training samples x features.
#' @param train_meta Data frame aligned to `train_expr`.
#' @param panel Character vector of panel features.
#' @param outcome_col Binary disease label column in `train_meta`.
#' @param cohort_col Training cohort column.
#' @param block_col Training provenance-block column.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return A fit object consumed by `predict_mixed_effects_scorer()`.
#' @export
fit_mixed_effects_scorer <- function(train_expr, train_meta, panel,
                                     outcome_col = "disease",
                                     cohort_col = "cohort",
                                     block_col = "provenance_block",
                                     pseudocount = NULL) {
  checked <- .pa_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  for (cc in c(outcome_col, cohort_col, block_col)) {
    if (!(cc %in% names(train_meta))) stop("missing column in train_meta: ", cc)
  }
  panel <- .pa_panel(train_expr, panel)
  y <- train_meta[[outcome_col]]
  if (is.factor(y)) y <- as.integer(y) - 1L
  y <- as.numeric(y)
  if (length(unique(y[is.finite(y)])) != 2L) {
    stop("outcome_col must contain two disease classes in the training pool")
  }
  log_expr <- .pa_log_expr(train_expr, pseudocount)
  weights <- rep(NA_real_, length(panel))
  names(weights) <- panel
  model_type <- rep(NA_character_, length(panel))
  names(model_type) <- panel
  for (feat in panel) {
    ww <- .pa_feature_weight(log_expr[, feat], y,
                             train_meta[[cohort_col]],
                             train_meta[[block_col]])
    weights[[feat]] <- as.numeric(ww[["weight"]])
    model_type[[feat]] <- as.character(ww[["model"]])
  }
  ok <- is.finite(weights)
  structure(
    list(panel = panel, weights = weights, model_type = model_type,
         pseudocount = pseudocount, outcome_col = outcome_col,
         cohort_col = cohort_col, block_col = block_col,
         n_features_fit = sum(ok),
         n_features_requested = length(panel)),
    class = "mixed_effects_scorer_fit"
  )
}

#' Predict a mixed-effects panel score without reading test labels
#'
#' @param fit Object from `fit_mixed_effects_scorer()`.
#' @param expr Numeric matrix, test samples x features.
#' @param meta Data frame aligned to `expr`; labels are not used.
#' @param cohort_col Test cohort column.
#' @param block_col Test provenance-block column.
#' @param min_weighted_features Minimum finite weights required for scoring.
#' @return Numeric score vector, one value per sample.
#' @export
predict_mixed_effects_scorer <- function(fit, expr, meta,
                                         cohort_col = fit$cohort_col,
                                         block_col = fit$block_col,
                                         min_weighted_features = 3L) {
  checked <- .pa_check_expr_meta(expr, meta)
  expr <- checked$expr
  panel <- intersect(fit$panel, colnames(expr))
  weights <- fit$weights[panel]
  ok <- is.finite(weights)
  min_weighted_features <- max(1L, as.integer(min_weighted_features))
  if (sum(ok) < min_weighted_features) {
    out <- rep(NA_real_, nrow(expr))
    attr(out, "n_weighted_features") <- sum(ok)
    return(out)
  }
  log_expr <- .pa_log_expr(expr, fit$pseudocount)
  out <- as.numeric(log_expr[, panel[ok], drop = FALSE] %*% weights[ok])
  attr(out, "n_weighted_features") <- sum(ok)
  out
}

#' Score with the mixed-effects scorer
#'
#' If `fit` is supplied, this function predicts from the existing fit and does
#' not read disease labels from `meta`. If `fit` is NULL, `train_expr`,
#' `train_meta`, and `outcome_col` are required and the model is fit on that
#' training pool before prediction.
#'
#' @param expr Numeric matrix, test samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param cohort_col Cohort column.
#' @param block_col Provenance-block column.
#' @param fit Optional fit from `fit_mixed_effects_scorer()`.
#' @param train_expr Optional training expression matrix when `fit` is NULL.
#' @param train_meta Optional training metadata when `fit` is NULL.
#' @param outcome_col Binary disease label column in `train_meta`.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Numeric score vector.
#' @export
score_mixed_effects_scorer <- function(expr, meta, panel,
                                       cohort_col = "cohort",
                                       block_col = "provenance_block",
                                       fit = NULL,
                                       train_expr = NULL,
                                       train_meta = NULL,
                                       outcome_col = "disease",
                                       pseudocount = NULL) {
  if (is.null(fit)) {
    if (is.null(train_expr) || is.null(train_meta)) {
      stop("train_expr and train_meta are required when fit is NULL")
    }
    fit <- fit_mixed_effects_scorer(
      train_expr = train_expr, train_meta = train_meta, panel = panel,
      outcome_col = outcome_col, cohort_col = cohort_col,
      block_col = block_col, pseudocount = pseudocount)
  }
  predict_mixed_effects_scorer(fit, expr, meta, cohort_col, block_col)
}

.pa_rclr_matrix <- function(expr, pseudocount = NULL, trim_upper = 0.10,
                            trim_lower = 0.05,
                            exclude_features = c("hsa-miR-451a",
                                                 "hsa-miR-16-5p",
                                                 "hsa-miR-486-5p",
                                                 "hsa-miR-144-3p",
                                                 "hsa-miR-223-3p")) {
  log_expr <- .pa_log_expr(expr, pseudocount)
  feat <- colnames(expr)
  excl <- feat %in% exclude_features
  out <- log_expr
  for (i in seq_len(nrow(expr))) {
    lx <- log_expr[i, ]
    keep <- is.finite(lx) & !excl
    if (sum(keep) >= 8L) {
      vals <- lx[keep]
      lo <- stats::quantile(vals, trim_lower, na.rm = TRUE, names = FALSE)
      hi <- stats::quantile(vals, 1 - trim_upper, na.rm = TRUE, names = FALSE)
      keep <- keep & lx >= lo & lx <= hi
    }
    if (sum(keep) < 2L) keep <- is.finite(lx)
    out[i, ] <- lx - mean(lx[keep], na.rm = TRUE)
  }
  out
}

.pa_alr_matrix <- function(expr, pivot_features = c("hsa-miR-103a-3p",
                                                    "hsa-miR-191-5p",
                                                    "hsa-miR-26a-5p",
                                                    "hsa-miR-30c-5p",
                                                    "hsa-let-7g-5p",
                                                    "hsa-miR-93-5p"),
                           pseudocount = NULL, min_pivot_present = 3L) {
  log_expr <- .pa_log_expr(expr, pseudocount)
  piv <- intersect(pivot_features, colnames(expr))
  if (length(piv) < min_pivot_present) {
    pivot_center <- rowMeans(log_expr, na.rm = TRUE)
  } else {
    pivot_center <- rowMeans(log_expr[, piv, drop = FALSE], na.rm = TRUE)
  }
  sweep(log_expr, 1L, pivot_center, "-")
}

.pa_simple_ilr_score <- function(expr, panel, pseudocount = NULL) {
  log_expr <- .pa_log_expr(expr, pseudocount)
  panel <- intersect(panel, colnames(expr))
  if (length(panel) < 2L) return(rep(NA_real_, nrow(expr)))
  hemo <- intersect(c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                     "hsa-miR-144-3p", "hsa-miR-223-3p"), panel)
  non_hemo <- setdiff(panel, hemo)
  if (length(hemo) > 0L && length(non_hemo) > 0L) {
    coef <- sqrt(length(hemo) * length(non_hemo) /
                   (length(hemo) + length(non_hemo)))
    return(coef * (rowMeans(log_expr[, hemo, drop = FALSE], na.rm = TRUE) -
                     rowMeans(log_expr[, non_hemo, drop = FALSE], na.rm = TRUE)))
  }
  k <- max(1L, floor(length(panel) / 2L))
  row_rank_score <- apply(log_expr[, panel, drop = FALSE], 1L, function(v) {
    ord <- order(v, decreasing = TRUE)
    top <- ord[seq_len(k)]
    tail <- setdiff(seq_along(v), top)
    if (length(tail) == 0L) return(NA_real_)
    coef <- sqrt(length(top) * length(tail) / (length(top) + length(tail)))
    coef * (mean(v[top], na.rm = TRUE) - mean(v[tail], na.rm = TRUE))
  })
  as.numeric(row_rank_score)
}

.pa_robust_scale_vec <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  sc <- .pa_mad(x)
  (x - med) / sc
}

#' Score a cluster-robust ensemble of rCLR, ALR, and ILR-style components
#'
#' The sample-level score is the row mean of robustly scaled rCLR, ALR, and
#' simple ILR component scores. Cluster-robust confidence intervals are not a
#' property of a single score vector; the benchmark script estimates them by
#' resampling cohort/provenance clusters across held-out folds.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param cohort_col Column in `meta` identifying cohorts.
#' @param block_col Column in `meta` identifying provenance blocks.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Numeric score vector. Attribute `component_scores` contains the
#'   finite component matrix.
#' @export
score_cluster_robust_ensemble <- function(expr, meta, panel,
                                          cohort_col = "cohort",
                                          block_col = "provenance_block",
                                          pseudocount = NULL) {
  checked <- .pa_check_expr_meta(expr, meta)
  expr <- checked$expr
  panel <- .pa_panel(expr, panel)
  rclr <- rowSums(.pa_rclr_matrix(expr, pseudocount)[, panel, drop = FALSE],
                  na.rm = TRUE)
  alr <- rowSums(.pa_alr_matrix(expr, pseudocount = pseudocount)[, panel,
                                                                 drop = FALSE],
                 na.rm = TRUE)
  ilr <- .pa_simple_ilr_score(expr, panel, pseudocount)
  comp <- cbind(rCLR = .pa_robust_scale_vec(rclr),
                ALR = .pa_robust_scale_vec(alr),
                ILR = .pa_robust_scale_vec(ilr))
  comp <- comp[, colSums(is.finite(comp)) > 0L, drop = FALSE]
  out <- rowMeans(comp, na.rm = TRUE)
  attr(out, "component_scores") <- comp
  attr(out, "n_components") <- ncol(comp)
  out
}

#' Score baseline robust CLR for provenance-aware benchmark comparisons
#'
#' @inheritParams score_cluster_robust_ensemble
#' @return Numeric score vector.
#' @export
score_baseline_rclr <- function(expr, meta, panel,
                                cohort_col = "cohort",
                                block_col = "provenance_block",
                                pseudocount = NULL) {
  checked <- .pa_check_expr_meta(expr, meta)
  expr <- checked$expr
  panel <- .pa_panel(expr, panel)
  rowSums(.pa_rclr_matrix(expr, pseudocount)[, panel, drop = FALSE],
          na.rm = TRUE)
}
