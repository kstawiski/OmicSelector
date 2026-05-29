#' Technology-aware within-sample scoring methods
#'
#' These scorers keep the platform technology axis separate from kit labels.
#' `technology_col` must identify one of `NGS`, `Toray`, `qPCR`,
#' `NanoString`, or `unknown` for every sample. Training helpers fit only on
#' the supplied training pool; prediction helpers do not read test labels.
#'
#' @name technology_aware_scoring
NULL

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}

.ta_allowed_technologies <- c("NGS", "Toray", "qPCR", "NanoString", "unknown")

.ta_check_expr_meta <- function(expr, meta) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  if (!is.data.frame(meta)) meta <- as.data.frame(meta)
  if (nrow(expr) != nrow(meta)) {
    stop("meta must have one row per expression sample")
  }
  if (is.null(colnames(expr))) {
    stop("expr must have feature names in colnames(expr)")
  }
  if (any(expr < 0, na.rm = TRUE)) {
    stop("technology-aware scorers require non-negative abundance-like input")
  }
  list(expr = expr, meta = meta)
}

.ta_panel <- function(expr, panel) {
  panel <- intersect(as.character(panel), colnames(expr))
  if (length(panel) == 0L) stop("none of the panel features are present in expr")
  panel
}

.ta_clean_technology <- function(x) {
  out <- trimws(as.character(x))
  out[is.na(out) | !nzchar(out)] <- "unknown"
  map <- c(
    "ngs" = "NGS", "mirna-seq" = "NGS", "mirna_seq" = "NGS",
    "small-rna-seq" = "NGS", "small_rna_seq" = "NGS",
    "toray" = "Toray", "array-toray3d" = "Toray",
    "array_toray3d" = "Toray", "toray-3dgene" = "Toray",
    "qpcr" = "qPCR", "q-pcr" = "qPCR", "rt-pcr" = "qPCR",
    "nanostring" = "NanoString", "nano-string" = "NanoString",
    "array-other" = "NanoString"
  )
  key <- tolower(out)
  out <- ifelse(key %in% names(map), unname(map[key]), out)
  out[!(out %in% .ta_allowed_technologies)] <- NA_character_
  if (anyNA(out)) {
    bad <- unique(trimws(as.character(x))[is.na(out)])
    stop("technology_col contains unsupported values: ",
         paste(bad, collapse = ", "),
         ". Allowed values are: ",
         paste(.ta_allowed_technologies, collapse = ", "))
  }
  out
}

.ta_get_technology <- function(meta, technology_col) {
  if (!(technology_col %in% names(meta))) {
    stop("technology_col not found in meta: ", technology_col)
  }
  .ta_clean_technology(meta[[technology_col]])
}

.ta_row_pseudocount <- function(expr, pseudocount = NULL) {
  if (!is.null(pseudocount)) return(rep(as.numeric(pseudocount), nrow(expr)))
  pc <- 1e-6 * rowSums(expr, na.rm = TRUE)
  pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  pc
}

.ta_log_expr <- function(expr, pseudocount = NULL) {
  pc <- .ta_row_pseudocount(expr, pseudocount)
  log(sweep(expr, 1L, pc, `+`))
}

.ta_median_or_na <- function(x) {
  out <- stats::median(x, na.rm = TRUE)
  if (!is.finite(out)) NA_real_ else out
}

.ta_tech_centers <- function(log_expr, technology, features) {
  split_idx <- split(seq_len(nrow(log_expr)), technology)
  lapply(split_idx, function(idx) {
    center <- apply(log_expr[idx, features, drop = FALSE], 2L,
                    .ta_median_or_na)
    list(center = center, n = length(idx))
  })
}

#' Fit technology-stratified rCLR centering moments
#'
#' @param train_expr Numeric matrix, samples x features, training pool.
#' @param train_meta Data frame aligned to `train_expr`.
#' @param panel Character vector of panel features.
#' @param technology_col Column identifying technology class.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Fit object consumed by `predict_tech_stratified_rclr()`.
#' @export
fit_tech_stratified_rclr <- function(train_expr, train_meta, panel,
                                     technology_col = "technology",
                                     pseudocount = NULL) {
  checked <- .ta_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  technology <- .ta_get_technology(train_meta, technology_col)
  panel <- .ta_panel(train_expr, panel)
  log_expr <- .ta_log_expr(train_expr, pseudocount)
  centers <- .ta_tech_centers(log_expr, technology, panel)
  global_center <- apply(log_expr[, panel, drop = FALSE], 2L,
                         .ta_median_or_na)
  structure(
    list(panel = panel, technology_col = technology_col,
         pseudocount = pseudocount, centers = centers,
         global_center = global_center,
         technology_levels = names(centers)),
    class = "tech_stratified_rclr_fit"
  )
}

#' Predict technology-stratified rCLR scores
#'
#' @param fit Fit from `fit_tech_stratified_rclr()`.
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param technology_col Column identifying technology class.
#' @return Numeric score vector. Attribute `fallback_samples` counts samples
#'   scored with the training global center because their technology was absent
#'   from the training fit.
#' @export
predict_tech_stratified_rclr <- function(fit, expr, meta,
                                         technology_col = fit$technology_col) {
  checked <- .ta_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  technology <- .ta_get_technology(meta, technology_col)
  panel <- .ta_panel(expr, fit$panel)
  log_expr <- .ta_log_expr(expr, fit$pseudocount)
  out <- rep(NA_real_, nrow(expr))
  fallback <- 0L
  for (tech in unique(technology)) {
    idx <- which(technology == tech)
    mm <- fit$centers[[tech]]
    center <- if (!is.null(mm)) mm$center else fit$global_center
    if (is.null(mm)) fallback <- fallback + length(idx)
    keep <- intersect(panel, names(center))
    if (length(keep) == 0L) next
    z <- sweep(log_expr[idx, keep, drop = FALSE], 2L, center[keep], "-")
    out[idx] <- rowSums(z, na.rm = TRUE)
  }
  attr(out, "fallback_samples") <- fallback
  out
}

#' Score technology-stratified rCLR
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param technology_col Column identifying technology class.
#' @param fit Optional fit from `fit_tech_stratified_rclr()`.
#' @param train_expr Optional training matrix when `fit` is NULL.
#' @param train_meta Optional training metadata when `fit` is NULL.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Numeric score vector.
#' @export
score_tech_stratified_rclr <- function(expr, meta, panel,
                                       technology_col = "technology",
                                       fit = NULL,
                                       train_expr = NULL,
                                       train_meta = NULL,
                                       pseudocount = NULL) {
  if (is.null(fit)) {
    if (is.null(train_expr) || is.null(train_meta)) {
      train_expr <- expr
      train_meta <- meta
    }
    fit <- fit_tech_stratified_rclr(train_expr, train_meta, panel,
                                    technology_col = technology_col,
                                    pseudocount = pseudocount)
  }
  predict_tech_stratified_rclr(fit, expr, meta, technology_col)
}

#' Select cross-technology anchor features on a training pool
#'
#' A selected anchor must be detected in every required technology class present
#' in `required_technologies`. Selection is based only on `train_expr` and
#' `train_meta`; held-out samples should never be passed here.
#'
#' @param train_expr Numeric matrix, samples x features, training pool.
#' @param train_meta Data frame aligned to `train_expr`.
#' @param technology_col Column identifying technology class.
#' @param required_technologies Technology classes required for selection.
#' @param min_detection_rate Minimum per-technology detection rate.
#' @param n_anchor Maximum number of anchors to select.
#' @return Data frame with detection rates and `selected_yes`.
#' @export
select_cross_tech_anchor_features <- function(
    train_expr, train_meta, technology_col = "technology",
    required_technologies = c("NGS", "Toray", "qPCR", "NanoString"),
    min_detection_rate = 0.05, n_anchor = 20L) {
  checked <- .ta_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  technology <- .ta_get_technology(train_meta, technology_col)
  req <- intersect(required_technologies, .ta_allowed_technologies)
  if (length(req) == 0L) stop("required_technologies has no valid technology labels")
  det <- lapply(req, function(tech) {
    idx <- which(technology == tech)
    if (length(idx) == 0L) {
      rep(NA_real_, ncol(train_expr))
    } else {
      colMeans(is.finite(train_expr[idx, , drop = FALSE]) &
                 train_expr[idx, , drop = FALSE] > 0, na.rm = TRUE)
    }
  })
  det_mat <- do.call(cbind, det)
  colnames(det_mat) <- paste0("detection_rate_", req)
  rownames(det_mat) <- colnames(train_expr)
  n_detected <- rowSums(det_mat >= min_detection_rate, na.rm = TRUE)
  min_det <- apply(det_mat, 1L, function(v) {
    if (any(!is.finite(v))) return(NA_real_)
    min(v)
  })
  mean_det <- rowMeans(det_mat, na.rm = TRUE)
  eligible <- n_detected == length(req) & is.finite(min_det)
  ord <- order(eligible, min_det, mean_det, decreasing = TRUE)
  selected <- rep(FALSE, ncol(train_expr))
  if (any(eligible)) {
    pick <- which(eligible)[order(min_det[eligible], mean_det[eligible],
                                  decreasing = TRUE)]
    selected[pick[seq_len(min(length(pick), as.integer(n_anchor)))]] <- TRUE
  }
  out <- data.frame(
    feature = colnames(train_expr),
    n_techs_detected = as.integer(n_detected),
    det_mat,
    min_detection_rate_across_required = min_det,
    mean_detection_rate_across_required = mean_det,
    selected_yes = ifelse(selected, "yes", "no"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out[ord, , drop = FALSE]
}

#' Fit cross-technology harmonized rCLR anchors
#'
#' @inheritParams select_cross_tech_anchor_features
#' @param panel Character vector of panel features.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Fit object consumed by `predict_cross_tech_harmonized_rclr()`.
#' @export
fit_cross_tech_harmonized_rclr <- function(
    train_expr, train_meta, panel, technology_col = "technology",
    required_technologies = c("NGS", "Toray", "qPCR", "NanoString"),
    min_detection_rate = 0.05, n_anchor = 20L, pseudocount = NULL) {
  checked <- .ta_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  panel <- .ta_panel(train_expr, panel)
  anchors <- select_cross_tech_anchor_features(
    train_expr, train_meta, technology_col = technology_col,
    required_technologies = required_technologies,
    min_detection_rate = min_detection_rate, n_anchor = n_anchor)
  anchor_features <- anchors$feature[anchors$selected_yes == "yes"]
  structure(
    list(panel = panel, technology_col = technology_col,
         required_technologies = required_technologies,
         min_detection_rate = min_detection_rate,
         anchor_features = anchor_features,
         anchor_table = anchors,
         pseudocount = pseudocount),
    class = "cross_tech_harmonized_rclr_fit"
  )
}

#' Predict cross-technology harmonized rCLR scores
#'
#' @param fit Fit from `fit_cross_tech_harmonized_rclr()`.
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`; technology labels are validated but
#'   not used to select anchors.
#' @param technology_col Column identifying technology class.
#' @param min_anchors_present Minimum anchors required in `expr`.
#' @return Numeric score vector with attribute `n_anchors_present`.
#' @export
predict_cross_tech_harmonized_rclr <- function(
    fit, expr, meta, technology_col = fit$technology_col,
    min_anchors_present = 3L) {
  checked <- .ta_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  invisible(.ta_get_technology(meta, technology_col))
  panel <- .ta_panel(expr, fit$panel)
  anchors <- intersect(fit$anchor_features, colnames(expr))
  min_anchors_present <- as.integer(min_anchors_present)
  if (length(anchors) < min_anchors_present) {
    out <- rep(NA_real_, nrow(expr))
    attr(out, "n_anchors_present") <- length(anchors)
    attr(out, "failure_reason") <- "insufficient_cross_technology_anchors"
    return(out)
  }
  log_expr <- .ta_log_expr(expr, fit$pseudocount)
  center <- rowMeans(log_expr[, anchors, drop = FALSE], na.rm = TRUE)
  z <- sweep(log_expr[, panel, drop = FALSE], 1L, center, "-")
  out <- rowSums(z, na.rm = TRUE)
  attr(out, "n_anchors_present") <- length(anchors)
  out
}

#' Score cross-technology harmonized rCLR
#'
#' @inheritParams score_tech_stratified_rclr
#' @param required_technologies Technology classes required for anchor selection.
#' @param min_detection_rate Minimum per-technology detection rate.
#' @param n_anchor Maximum number of anchors to select.
#' @return Numeric score vector.
#' @export
score_cross_tech_harmonized_rclr <- function(
    expr, meta, panel, technology_col = "technology", fit = NULL,
    train_expr = NULL, train_meta = NULL,
    required_technologies = c("NGS", "Toray", "qPCR", "NanoString"),
    min_detection_rate = 0.05, n_anchor = 20L, pseudocount = NULL) {
  if (is.null(fit)) {
    if (is.null(train_expr) || is.null(train_meta)) {
      train_expr <- expr
      train_meta <- meta
    }
    fit <- fit_cross_tech_harmonized_rclr(
      train_expr, train_meta, panel, technology_col = technology_col,
      required_technologies = required_technologies,
      min_detection_rate = min_detection_rate, n_anchor = n_anchor,
      pseudocount = pseudocount)
  }
  predict_cross_tech_harmonized_rclr(fit, expr, meta, technology_col)
}

.ta_default_pivot_pool <- function() {
  c("hsa-miR-103a-3p", "hsa-miR-191-5p", "hsa-miR-26a-5p",
    "hsa-miR-30c-5p", "hsa-let-7g-5p", "hsa-miR-93-5p")
}

.ta_lm_tech_effect <- function(value, disease, technology, cohort, tech_levels) {
  df <- data.frame(value = as.numeric(value),
                   disease = as.numeric(disease),
                   technology = factor(technology, levels = tech_levels),
                   cohort = factor(cohort))
  ok <- stats::complete.cases(df)
  df <- df[ok, , drop = FALSE]
  empty <- stats::setNames(rep(NA_real_, length(tech_levels)), tech_levels)
  if (nrow(df) < 8L || nlevels(droplevels(df$technology)) < 1L) {
    return(list(effect = empty, model = "not_fit"))
  }
  has_disease <- length(unique(df$disease[is.finite(df$disease)])) == 2L
  has_cohort <- nlevels(droplevels(df$cohort)) > 1L
  rhs <- c(if (has_disease) "disease" else NULL, "technology")
  if (requireNamespace("lme4", quietly = TRUE) && has_cohort) {
    form <- stats::as.formula(paste("value ~", paste(rhs, collapse = " + "),
                                    "+ (1 | cohort)"))
    fit <- tryCatch(
      lme4::lmer(form, data = df, REML = TRUE,
                 control = lme4::lmerControl(check.rankX = "silent.drop.cols")),
      error = function(e) NULL)
    if (!is.null(fit)) {
      co <- lme4::fixef(fit)
      eff <- stats::setNames(rep(0, length(tech_levels)), tech_levels)
      for (tech in tech_levels[-1L]) {
        nm <- paste0("technology", tech)
        if (nm %in% names(co) && is.finite(co[[nm]])) eff[[tech]] <- co[[nm]]
      }
      return(list(effect = eff, model = "lmer_REML"))
    }
  }
  form <- stats::as.formula(paste("value ~", paste(c(rhs, if (has_cohort) "cohort"),
                                                   collapse = " + ")))
  fit <- tryCatch(stats::lm(form, data = df), error = function(e) NULL)
  if (is.null(fit)) return(list(effect = empty, model = "not_fit"))
  co <- stats::coef(fit)
  eff <- stats::setNames(rep(0, length(tech_levels)), tech_levels)
  for (tech in tech_levels[-1L]) {
    nm <- paste0("technology", tech)
    if (nm %in% names(co) && is.finite(co[[nm]])) eff[[tech]] <- co[[nm]]
  }
  list(effect = eff, model = "lm_fallback")
}

#' Fit technology-residualized ALR
#'
#' Per feature, the training model is
#' `log(feature_count) ~ disease + technology + (1|cohort)` when `lme4` is
#' available and more than one cohort exists. Prediction subtracts the learned
#' technology fixed effect and then computes an ALR-style panel score.
#'
#' @param train_expr Numeric matrix, training samples x features.
#' @param train_meta Data frame aligned to `train_expr`.
#' @param panel Character vector of panel features.
#' @param technology_col Column identifying technology class.
#' @param outcome_col Binary disease label column in `train_meta`.
#' @param cohort_col Cohort column in `train_meta`.
#' @param pivot_features Candidate ALR pivot features.
#' @param anchor_features Optional train-selected fallback pivots.
#' @param min_pivot_present Minimum pivots required for prediction.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @return Fit object consumed by `predict_tech_residualized_alr()`.
#' @export
fit_tech_residualized_alr <- function(
    train_expr, train_meta, panel, technology_col = "technology",
    outcome_col = "disease", cohort_col = "cohort",
    pivot_features = .ta_default_pivot_pool(), anchor_features = character(),
    min_pivot_present = 3L, pseudocount = NULL) {
  checked <- .ta_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  technology <- .ta_get_technology(train_meta, technology_col)
  if (!(outcome_col %in% names(train_meta))) stop("missing outcome_col in train_meta")
  if (!(cohort_col %in% names(train_meta))) stop("missing cohort_col in train_meta")
  panel <- .ta_panel(train_expr, panel)
  pivots <- intersect(pivot_features, colnames(train_expr))
  pivot_source <- "default_pivot_pool"
  if (length(pivots) < min_pivot_present) {
    pivots <- intersect(anchor_features, colnames(train_expr))
    pivot_source <- "train_anchor_fallback"
  }
  features <- unique(c(panel, pivots))
  log_expr <- .ta_log_expr(train_expr, pseudocount)
  tech_levels <- .ta_allowed_technologies[.ta_allowed_technologies %in% unique(technology)]
  if (length(tech_levels) == 0L) tech_levels <- "unknown"
  effects <- vector("list", length(features))
  names(effects) <- features
  model_type <- stats::setNames(rep("not_fit", length(features)), features)
  for (feat in features) {
    fit_feat <- .ta_lm_tech_effect(log_expr[, feat],
                                   train_meta[[outcome_col]],
                                   technology,
                                   train_meta[[cohort_col]],
                                   tech_levels)
    effects[[feat]] <- fit_feat$effect
    model_type[[feat]] <- fit_feat$model
  }
  structure(
    list(panel = panel, technology_col = technology_col,
         outcome_col = outcome_col, cohort_col = cohort_col,
         pivot_features = pivots, pivot_source = pivot_source,
         min_pivot_present = as.integer(min_pivot_present),
         technology_levels = tech_levels, feature_tech_effects = effects,
         model_type = model_type, pseudocount = pseudocount),
    class = "tech_residualized_alr_fit"
  )
}

#' Predict technology-residualized ALR scores
#'
#' @param fit Fit from `fit_tech_residualized_alr()`.
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`; disease labels are not used.
#' @param technology_col Column identifying technology class.
#' @return Numeric score vector with pivot and model audit attributes.
#' @export
predict_tech_residualized_alr <- function(
    fit, expr, meta, technology_col = fit$technology_col) {
  checked <- .ta_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  technology <- .ta_get_technology(meta, technology_col)
  panel <- .ta_panel(expr, fit$panel)
  pivots <- intersect(fit$pivot_features, colnames(expr))
  if (length(pivots) < fit$min_pivot_present) {
    out <- rep(NA_real_, nrow(expr))
    attr(out, "n_pivots_present") <- length(pivots)
    attr(out, "pivot_source") <- fit$pivot_source
    attr(out, "failure_reason") <- "insufficient_alr_pivots"
    return(out)
  }
  features <- unique(c(panel, pivots))
  log_expr <- .ta_log_expr(expr, fit$pseudocount)
  resid <- log_expr[, features, drop = FALSE]
  fallback <- 0L
  for (feat in features) {
    eff <- fit$feature_tech_effects[[feat]]
    if (is.null(eff)) next
    off <- eff[technology]
    missing_eff <- !is.finite(off)
    if (any(missing_eff)) {
      off[missing_eff] <- 0
      fallback <- fallback + sum(missing_eff)
    }
    resid[, feat] <- resid[, feat] - as.numeric(off)
  }
  pivot_center <- rowMeans(resid[, pivots, drop = FALSE], na.rm = TRUE)
  z <- sweep(resid[, panel, drop = FALSE], 1L, pivot_center, "-")
  out <- rowSums(z, na.rm = TRUE)
  attr(out, "n_pivots_present") <- length(pivots)
  attr(out, "pivot_source") <- fit$pivot_source
  attr(out, "technology_effect_fallback_cells") <- fallback
  attr(out, "model_types") <- sort(unique(fit$model_type))
  out
}

#' Score technology-residualized ALR
#'
#' @inheritParams score_tech_stratified_rclr
#' @param outcome_col Binary disease label column in `train_meta`.
#' @param cohort_col Cohort column in `train_meta`.
#' @param pivot_features Candidate ALR pivot features.
#' @param anchor_features Optional train-selected fallback pivots.
#' @param min_pivot_present Minimum pivots required for prediction.
#' @return Numeric score vector.
#' @export
score_tech_residualized_alr <- function(
    expr, meta, panel, technology_col = "technology", fit = NULL,
    train_expr = NULL, train_meta = NULL, outcome_col = "disease",
    cohort_col = "cohort", pivot_features = .ta_default_pivot_pool(),
    anchor_features = character(), min_pivot_present = 3L,
    pseudocount = NULL) {
  if (is.null(fit)) {
    if (is.null(train_expr) || is.null(train_meta)) {
      stop("train_expr and train_meta are required when fit is NULL")
    }
    fit <- fit_tech_residualized_alr(
      train_expr, train_meta, panel, technology_col = technology_col,
      outcome_col = outcome_col, cohort_col = cohort_col,
      pivot_features = pivot_features, anchor_features = anchor_features,
      min_pivot_present = min_pivot_present, pseudocount = pseudocount)
  }
  predict_tech_residualized_alr(fit, expr, meta, technology_col)
}

#' Paired DeLong SE for technology-aware transfer lift
#'
#' Computes the paired DeLong AUC-difference SE for a technology-aware transfer
#' method against its paired rCLR baseline on the same held-out samples. This
#' uses `orient = "auc"` to match the technology-aware transfer convention
#' `max(AUC, 1 - AUC)`.
#'
#' @param y Binary held-out outcome vector.
#' @param method_scores Numeric method scores on the held-out samples.
#' @param baseline_scores Numeric paired baseline scores on the same samples.
#' @param fold Optional fold identifier for mean-of-folds aggregation.
#'
#' @return A list from `paper3_paired_auc_diff_se()`.
#' @export
paper3_technology_lift_delong <- function(y, method_scores, baseline_scores,
                                          fold = NULL) {
  paper3_paired_auc_diff_se(
    y = y,
    score_method = method_scores,
    score_baseline = baseline_scores,
    fold = fold,
    orient = "auc"
  )
}
