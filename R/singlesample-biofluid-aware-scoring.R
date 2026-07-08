#' Biofluid-aware within-sample scoring methods
#'
#' These scorers treat serum, plasma anticoagulant class, exosome-enriched
#' fractions, and unknown biofluid status as first-class pre-analytical strata.
#' All functions expect expression matrices with samples in rows and miRNAs in
#' columns. Fitting functions use training data only; prediction functions do
#' not read disease labels from test metadata.
#'
#' @name biofluid_aware_scoring
NULL

.bfa_default_biofluid_levels <- c("serum-clot", "plasma-EDTA",
                                  "plasma-other", "exosome", "unknown")

.bfa_default_pivot_pool <- function() {
  c("hsa-miR-103a-3p", "hsa-miR-191-5p", "hsa-miR-26a-5p",
    "hsa-miR-30c-5p", "hsa-let-7g-5p", "hsa-miR-93-5p")
}

.bfa_default_blood_cell_mirnas <- function() {
  c("hsa-miR-451a", "hsa-miR-486-5p", "hsa-miR-21-5p",
    "hsa-miR-16-5p", "hsa-miR-92a-3p", "hsa-miR-191-5p",
    "miR-451a", "miR-486-5p", "miR-21-5p", "miR-16-5p",
    "miR-92a", "miR-92a-3p", "miR-191", "miR-191-5p")
}

.bfa_check_expr_meta <- function(expr, meta) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  if (!is.data.frame(meta)) meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (nrow(expr) != nrow(meta)) stop("meta must have one row per expression sample")
  if (is.null(colnames(expr))) stop("expr must have feature names in colnames(expr)")
  if (nrow(expr) == 0L || ncol(expr) == 0L) stop("expr must be non-empty")
  if (any(expr < 0, na.rm = TRUE)) {
    stop("biofluid-aware scorers require non-negative abundance-like input")
  }
  list(expr = expr, meta = meta)
}

.bfa_panel <- function(expr, panel, min_panel_present = 1L) {
  panel <- intersect(as.character(panel), colnames(expr))
  if (length(panel) < min_panel_present) {
    stop("fewer than min_panel_present panel features are present in expr")
  }
  panel
}

.bfa_canonical_biofluid <- function(x) {
  z <- trimws(tolower(as.character(x)))
  z[is.na(z) | !nzchar(z)] <- "unknown"
  out <- rep("unknown", length(z))
  out[grepl("exosome|extracellular vesicle|\\bev\\b|vesicle", z)] <- "exosome"
  out[grepl("serum|clot", z) & out == "unknown"] <- "serum-clot"
  out[grepl("edta", z) & out == "unknown"] <- "plasma-EDTA"
  out[grepl("citrate|citrat|heparin|plasma|nr-anticoag|plasma-nr", z) &
        out == "unknown"] <- "plasma-other"
  factor(out, levels = .bfa_default_biofluid_levels)
}

.bfa_row_pseudocount <- function(expr, pseudocount = NULL) {
  if (!is.null(pseudocount)) return(rep(as.numeric(pseudocount), nrow(expr)))
  pc <- 1e-6 * rowSums(expr, na.rm = TRUE)
  pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  pc
}

.bfa_log_abundance <- function(expr, pseudocount = NULL) {
  pc <- .bfa_row_pseudocount(expr, pseudocount)
  log(sweep(expr, 1L, pc, "+"))
}

.bfa_log_proportions <- function(expr, pseudocount = NULL) {
  pc <- .bfa_row_pseudocount(expr, pseudocount)
  xpos <- sweep(expr, 1L, pc, "+")
  denom <- rowSums(xpos, na.rm = TRUE)
  denom[!is.finite(denom) | denom <= 0] <- NA_real_
  sweep(log(xpos), 1L, log(denom), "-")
}

.bfa_col_median <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colMedians(as.matrix(x), na.rm = TRUE)
  } else {
    apply(x, 2L, stats::median, na.rm = TRUE)
  }
}

.bfa_trimmed_centering_features <- function(logp,
                                            exclude_features = .bfa_default_blood_cell_mirnas(),
                                            trim_upper = 0.10,
                                            trim_lower = 0.05,
                                            min_centering_size = 8L) {
  if (ncol(logp) == 0L) return(character(0))
  med <- .bfa_col_median(logp)
  names(med) <- colnames(logp)
  candidates <- setdiff(names(med)[is.finite(med)], exclude_features)
  if (length(candidates) < min_centering_size) candidates <- names(med)[is.finite(med)]
  if (length(candidates) <= min_centering_size) return(candidates)
  ord <- candidates[order(med[candidates], decreasing = TRUE)]
  n <- length(ord)
  lo <- min(n, floor(trim_upper * n) + 1L)
  hi <- max(lo, n - floor(trim_lower * n))
  keep <- ord[seq.int(lo, hi)]
  if (length(keep) < min_centering_size) candidates else keep
}

.bfa_weighted_panel_sum <- function(z, panel, feature_weights = NULL) {
  panel <- intersect(panel, colnames(z))
  if (length(panel) == 0L) return(rep(NA_real_, nrow(z)))
  if (is.null(feature_weights)) {
    w <- stats::setNames(rep(1, length(panel)), panel)
  } else {
    w <- as.numeric(feature_weights[panel])
    w[!is.finite(w)] <- 1
    w <- stats::setNames(w, panel)
  }
  as.numeric(z[, panel, drop = FALSE] %*% w)
}

.bfa_global_rclr_matrix <- function(expr, pseudocount = NULL,
                                    exclude_features = .bfa_default_blood_cell_mirnas(),
                                    trim_upper = 0.10,
                                    trim_lower = 0.05,
                                    min_centering_size = 8L) {
  logp <- .bfa_log_proportions(expr, pseudocount)
  out <- logp
  for (i in seq_len(nrow(logp))) {
    row <- logp[i, , drop = FALSE]
    center <- .bfa_trimmed_centering_features(
      row, exclude_features = exclude_features, trim_upper = trim_upper,
      trim_lower = trim_lower, min_centering_size = min_centering_size)
    if (length(center) < 2L) center <- colnames(logp)[is.finite(logp[i, ])]
    ctr <- mean(logp[i, center], na.rm = TRUE)
    out[i, ] <- logp[i, ] - ctr
  }
  out
}

#' Fit biofluid-stratified rCLR centering
#'
#' Selects rCLR centering features separately within each training biofluid
#' class. The held-out score is the signed sum of panel log-ratios centered by
#' the denominator selected for the sample's biofluid class.
#'
#' @param train_expr Numeric matrix, samples x features.
#' @param train_meta Data frame aligned to `train_expr`.
#' @param panel Character vector of panel features.
#' @param biofluid_col Column in `train_meta` with canonical or raw biofluid.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param exclude_features Features excluded from rCLR denominator selection.
#' @param trim_upper,trim_lower Feature trimming fractions for rCLR denominator
#'   selection.
#' @param min_centering_size Minimum denominator size before fallback.
#' @return A `biofluid_stratified_rclr_fit`.
#' @export
fit_biofluid_stratified_rclr <- function(train_expr, train_meta, panel,
                                         biofluid_col = "biofluid",
                                         pseudocount = NULL,
                                         exclude_features = .bfa_default_blood_cell_mirnas(),
                                         trim_upper = 0.10,
                                         trim_lower = 0.05,
                                         min_centering_size = 8L) {
  checked <- .bfa_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  if (!(biofluid_col %in% names(train_meta))) stop("missing biofluid_col in train_meta")
  panel <- .bfa_panel(train_expr, panel)
  bio <- .bfa_canonical_biofluid(train_meta[[biofluid_col]])
  logp <- .bfa_log_proportions(train_expr, pseudocount)
  centering <- list()
  n_by <- table(bio, useNA = "no")
  for (lev in names(n_by)) {
    idx <- which(as.character(bio) == lev)
    centering[[lev]] <- .bfa_trimmed_centering_features(
      logp[idx, , drop = FALSE], exclude_features = exclude_features,
      trim_upper = trim_upper, trim_lower = trim_lower,
      min_centering_size = min_centering_size)
  }
  global <- .bfa_trimmed_centering_features(
    logp, exclude_features = exclude_features, trim_upper = trim_upper,
    trim_lower = trim_lower, min_centering_size = min_centering_size)
  structure(
    list(panel = panel, biofluid_col = biofluid_col, pseudocount = pseudocount,
         centering_features = centering, global_centering_features = global,
         n_by_biofluid = as.integer(n_by), biofluid_levels = names(n_by),
         exclude_features = exclude_features, trim_upper = trim_upper,
         trim_lower = trim_lower, min_centering_size = min_centering_size),
    class = "biofluid_stratified_rclr_fit")
}

#' Predict biofluid-stratified rCLR scores
#'
#' @param fit Object from `fit_biofluid_stratified_rclr()`.
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param biofluid_col Column in `meta` with canonical or raw biofluid.
#' @param feature_weights Optional named numeric panel weights.
#' @return Numeric score vector.
#' @export
predict_biofluid_stratified_rclr <- function(fit, expr, meta,
                                             biofluid_col = fit$biofluid_col,
                                             feature_weights = NULL) {
  checked <- .bfa_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  if (!(biofluid_col %in% names(meta))) stop("missing biofluid_col in meta")
  panel <- .bfa_panel(expr, fit$panel)
  logp <- .bfa_log_proportions(expr, fit$pseudocount)
  bio <- as.character(.bfa_canonical_biofluid(meta[[biofluid_col]]))
  out_z <- matrix(NA_real_, nrow = nrow(expr), ncol = length(panel),
                  dimnames = list(rownames(expr), panel))
  fallback <- 0L
  for (lev in unique(bio)) {
    idx <- which(bio == lev)
    center <- fit$centering_features[[lev]]
    if (is.null(center) || length(center) < 2L) {
      center <- fit$global_centering_features
      fallback <- fallback + length(idx)
    }
    center <- intersect(center, colnames(logp))
    if (length(center) < 2L) center <- colnames(logp)[is.finite(colSums(logp))]
    ctr <- rowMeans(logp[idx, center, drop = FALSE], na.rm = TRUE)
    out_z[idx, ] <- sweep(logp[idx, panel, drop = FALSE], 1L, ctr, "-")
  }
  score <- .bfa_weighted_panel_sum(out_z, panel, feature_weights)
  attr(score, "fallback_samples") <- fallback
  score
}

#' Score biofluid-stratified rCLR
#'
#' @inheritParams predict_biofluid_stratified_rclr
#' @param panel Character vector of panel features.
#' @param fit Optional fit. If NULL, the fit is built on `expr` and `meta`.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param ... Additional arguments passed to `fit_biofluid_stratified_rclr()`.
#' @export
score_biofluid_stratified_rclr <- function(expr, meta, panel,
                                           biofluid_col = "biofluid",
                                           fit = NULL,
                                           feature_weights = NULL,
                                           pseudocount = NULL, ...) {
  if (is.null(fit)) {
    fit <- fit_biofluid_stratified_rclr(expr, meta, panel, biofluid_col,
                                        pseudocount = pseudocount, ...)
  }
  predict_biofluid_stratified_rclr(fit, expr, meta, biofluid_col, feature_weights)
}

.bfa_anchor_table <- function(expr, biofluid, panel = character(0),
                              exclude_features = .bfa_default_blood_cell_mirnas(),
                              pseudocount = NULL,
                              min_detection_rate_per_class = 0.10) {
  logp <- .bfa_log_proportions(expr, pseudocount)
  feats <- colnames(logp)
  rows <- lapply(feats, function(feat) {
    vals <- logp[, feat]
    cls <- split(seq_along(vals), as.character(biofluid))
    med <- vapply(cls, function(idx) {
      det <- mean(is.finite(expr[idx, feat]) & expr[idx, feat] > 0, na.rm = TRUE)
      if (!is.finite(det) || det < min_detection_rate_per_class) return(NA_real_)
      stats::median(vals[idx], na.rm = TRUE)
    }, numeric(1L))
    n_cls <- sum(is.finite(med))
    total_var <- stats::var(vals, na.rm = TRUE)
    between_var <- if (n_cls >= 2L) stats::var(med[is.finite(med)]) else NA_real_
    ratio <- between_var / max(total_var, .Machine$double.eps)
    data.frame(feature = feat, biofluid_variance_ratio = ratio,
               n_biofluid_classes_detected = n_cls,
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[is.finite(tab$biofluid_variance_ratio), , drop = FALSE]
  tab <- tab[!(tab$feature %in% exclude_features), , drop = FALSE]
  tab <- tab[order(tab$biofluid_variance_ratio, -tab$n_biofluid_classes_detected,
                   tab$feature), , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

#' Fit biofluid-stable anchor rCLR
#'
#' Selects denominator anchors with the lowest ratio of between-biofluid
#' variance to total variance in the training pool.
#'
#' @inheritParams fit_biofluid_stratified_rclr
#' @param anchor_n Number of anchor features to select.
#' @param min_anchor_n Minimum anchors required for an evaluable fit.
#' @param min_biofluid_classes Minimum biofluid classes in which an anchor must
#'   be detected.
#' @export
fit_biofluid_anchor_rclr <- function(train_expr, train_meta, panel,
                                     biofluid_col = "biofluid",
                                     pseudocount = NULL,
                                     anchor_n = 12L,
                                     min_anchor_n = 3L,
                                     min_biofluid_classes = 2L,
                                     exclude_features = .bfa_default_blood_cell_mirnas()) {
  checked <- .bfa_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  if (!(biofluid_col %in% names(train_meta))) stop("missing biofluid_col in train_meta")
  panel <- .bfa_panel(train_expr, panel)
  bio <- .bfa_canonical_biofluid(train_meta[[biofluid_col]])
  tab <- .bfa_anchor_table(train_expr, bio, panel = panel,
                           exclude_features = exclude_features,
                           pseudocount = pseudocount)
  candidates <- tab[tab$n_biofluid_classes_detected >= min_biofluid_classes, , drop = FALSE]
  candidates <- candidates[!(candidates$feature %in% panel), , drop = FALSE]
  if (nrow(candidates) < min_anchor_n) {
    candidates <- tab[tab$n_biofluid_classes_detected >= min_biofluid_classes, , drop = FALSE]
  }
  selected <- head(candidates$feature, anchor_n)
  if (length(selected) < min_anchor_n) selected <- head(tab$feature, anchor_n)
  tab$selected <- tab$feature %in% selected
  structure(
    list(panel = panel, biofluid_col = biofluid_col, pseudocount = pseudocount,
         anchor_features = selected, anchor_table = tab,
         min_anchor_n = min_anchor_n, anchor_n = anchor_n,
         exclude_features = exclude_features),
    class = "biofluid_anchor_rclr_fit")
}

#' Predict biofluid-stable anchor rCLR scores
#'
#' @param fit Object from `fit_biofluid_anchor_rclr()`.
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param biofluid_col Column in `meta` with canonical or raw biofluid.
#' @param feature_weights Optional named numeric panel weights.
#' @export
predict_biofluid_anchor_rclr <- function(fit, expr, meta,
                                         biofluid_col = fit$biofluid_col,
                                         feature_weights = NULL) {
  checked <- .bfa_check_expr_meta(expr, meta)
  expr <- checked$expr
  panel <- .bfa_panel(expr, fit$panel)
  logp <- .bfa_log_proportions(expr, fit$pseudocount)
  anchors <- intersect(fit$anchor_features, colnames(logp))
  if (length(anchors) < fit$min_anchor_n) {
    out <- rep(NA_real_, nrow(expr))
    attr(out, "n_anchor_present") <- length(anchors)
    return(out)
  }
  center <- rowMeans(logp[, anchors, drop = FALSE], na.rm = TRUE)
  z <- sweep(logp[, panel, drop = FALSE], 1L, center, "-")
  out <- .bfa_weighted_panel_sum(z, panel, feature_weights)
  attr(out, "n_anchor_present") <- length(anchors)
  out
}

#' Score biofluid-stable anchor rCLR
#'
#' @inheritParams score_biofluid_stratified_rclr
#' @export
score_biofluid_anchor_rclr <- function(expr, meta, panel,
                                       biofluid_col = "biofluid",
                                       fit = NULL,
                                       feature_weights = NULL,
                                       pseudocount = NULL, ...) {
  if (is.null(fit)) {
    fit <- fit_biofluid_anchor_rclr(expr, meta, panel, biofluid_col,
                                    pseudocount = pseudocount, ...)
  }
  predict_biofluid_anchor_rclr(fit, expr, meta, biofluid_col, feature_weights)
}

.bfa_biofluid_effects_one_feature <- function(value, disease, biofluid, cohort) {
  biofluid <- as.character(biofluid)
  ord <- names(sort(table(biofluid), decreasing = TRUE))
  bio_levels <- unique(c(ord, setdiff(.bfa_default_biofluid_levels, ord)))
  df <- data.frame(value = as.numeric(value),
                   disease = as.numeric(disease),
                   biofluid = factor(biofluid, levels = bio_levels),
                   cohort = factor(as.character(cohort)))
  df <- df[stats::complete.cases(df), , drop = FALSE]
  effects <- stats::setNames(rep(0, length(bio_levels)), bio_levels)
  if (nrow(df) < 12L || nlevels(droplevels(df$biofluid)) < 2L) {
    return(list(effect = effects, model = "not_fit"))
  }
  has_disease <- length(unique(df$disease[is.finite(df$disease)])) == 2L
  rhs <- if (has_disease) "disease + biofluid" else "biofluid"
  if (nlevels(droplevels(df$cohort)) > 1L) {
    if (requireNamespace("lme4", quietly = TRUE)) {
      f <- stats::as.formula(paste("value ~", rhs, "+ (1 | cohort)"))
      fit <- tryCatch(lme4::lmer(
        f, data = df, REML = FALSE,
        control = lme4::lmerControl(check.rankX = "silent.drop.cols",
                                    check.conv.singular = "ignore")),
        error = function(e) NULL)
      if (!is.null(fit)) {
        co <- lme4::fixef(fit)
        bio_terms <- grep("^biofluid", names(co), value = TRUE)
        levels_seen <- sub("^biofluid", "", bio_terms)
        effects[levels_seen] <- co[bio_terms]
        effects[!is.finite(effects)] <- 0
        return(list(effect = effects, model = "lmer"))
      }
    }
    f <- stats::as.formula(paste("value ~", rhs, "+ factor(cohort)"))
  } else {
    f <- stats::as.formula(paste("value ~", rhs))
  }
  fit <- tryCatch(stats::lm(f, data = df), error = function(e) NULL)
  if (is.null(fit)) return(list(effect = effects, model = "not_fit"))
  co <- stats::coef(fit)
  bio_terms <- grep("^biofluid", names(co), value = TRUE)
  levels_seen <- sub("^biofluid", "", bio_terms)
  effects[levels_seen] <- co[bio_terms]
  effects[!is.finite(effects)] <- 0
  list(effect = effects, model = "lm")
}

.bfa_low_variance_pivots <- function(log_expr, panel, exclude_features,
                                     n = 6L) {
  pool <- setdiff(colnames(log_expr), unique(c(panel, exclude_features)))
  if (length(pool) < n) pool <- setdiff(colnames(log_expr), exclude_features)
  if (length(pool) == 0L) return(character(0))
  v <- apply(log_expr[, pool, drop = FALSE], 2L, stats::mad, na.rm = TRUE)
  names(sort(v))[seq_len(min(n, length(v)))]
}

#' Fit biofluid-residualized ALR
#'
#' Fits one per-feature training model,
#' `log(feature_count) ~ disease + biofluid + (1 | cohort)`, subtracts the
#' estimated biofluid fixed effect at prediction time, and scores the panel as
#' an ALR against a frozen pivot pool. If `lme4` is unavailable or a feature
#' model fails, an `lm()` fallback is used and recorded.
#'
#' @param train_expr Numeric matrix, training samples x features.
#' @param train_meta Data frame aligned to `train_expr`.
#' @param panel Character vector of panel features.
#' @param biofluid_col Column identifying biofluid class.
#' @param outcome_col Binary disease label column in `train_meta`.
#' @param cohort_col Cohort column in `train_meta`.
#' @param pivot_features Optional ALR pivot features.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param min_pivot_present Minimum pivot features required for prediction.
#' @export
fit_biofluid_residualized_alr <- function(train_expr, train_meta, panel,
                                          biofluid_col = "biofluid",
                                          outcome_col = "disease",
                                          cohort_col = "cohort",
                                          pivot_features = NULL,
                                          pseudocount = NULL,
                                          min_pivot_present = 3L) {
  checked <- .bfa_check_expr_meta(train_expr, train_meta)
  train_expr <- checked$expr
  train_meta <- checked$meta
  for (cc in c(biofluid_col, outcome_col, cohort_col)) {
    if (!(cc %in% names(train_meta))) stop("missing column in train_meta: ", cc)
  }
  panel <- .bfa_panel(train_expr, panel)
  if (is.null(pivot_features)) pivot_features <- .bfa_default_pivot_pool()
  log_expr <- .bfa_log_abundance(train_expr, pseudocount)
  pivot <- intersect(pivot_features, colnames(log_expr))
  pivot_source <- "default_pivot_pool"
  if (length(pivot) < min_pivot_present) {
    pivot <- .bfa_low_variance_pivots(log_expr, panel,
                                      exclude_features = .bfa_default_blood_cell_mirnas(),
                                      n = max(6L, min_pivot_present))
    pivot_source <- "train_low_variance_fallback"
  }
  features <- unique(c(panel, pivot))
  bio <- .bfa_canonical_biofluid(train_meta[[biofluid_col]])
  y <- train_meta[[outcome_col]]
  if (is.factor(y)) y <- as.integer(y) - 1L
  y <- as.numeric(y)
  bio_levels <- .bfa_default_biofluid_levels
  effects <- matrix(0, nrow = length(features), ncol = length(bio_levels),
                    dimnames = list(features, bio_levels))
  model_type <- stats::setNames(rep("not_fit", length(features)), features)
  for (feat in features) {
    if (!(feat %in% colnames(log_expr))) next
    res <- .bfa_biofluid_effects_one_feature(
      log_expr[, feat], y, bio, train_meta[[cohort_col]])
    effects[feat, names(res$effect)] <- res$effect
    model_type[[feat]] <- res$model
  }
  structure(
    list(panel = panel, biofluid_col = biofluid_col, outcome_col = outcome_col,
         cohort_col = cohort_col, pivot_features = pivot,
         pivot_source = pivot_source, min_pivot_present = min_pivot_present,
         pseudocount = pseudocount, effects = effects,
         model_type = model_type,
         n_features_fit = sum(model_type != "not_fit"),
         n_features_requested = length(features)),
    class = "biofluid_residualized_alr_fit")
}

#' Predict biofluid-residualized ALR scores
#'
#' @param fit Object from `fit_biofluid_residualized_alr()`.
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param biofluid_col Column identifying biofluid class.
#' @param feature_weights Optional named numeric panel weights.
#' @export
predict_biofluid_residualized_alr <- function(fit, expr, meta,
                                              biofluid_col = fit$biofluid_col,
                                              feature_weights = NULL) {
  checked <- .bfa_check_expr_meta(expr, meta)
  expr <- checked$expr
  meta <- checked$meta
  if (!(biofluid_col %in% names(meta))) stop("missing biofluid_col in meta")
  panel <- .bfa_panel(expr, fit$panel)
  piv <- intersect(fit$pivot_features, colnames(expr))
  if (length(piv) < fit$min_pivot_present) {
    out <- rep(NA_real_, nrow(expr))
    attr(out, "n_pivot_present") <- length(piv)
    return(out)
  }
  features <- unique(c(panel, piv))
  log_expr <- .bfa_log_abundance(expr[, features, drop = FALSE], fit$pseudocount)
  bio <- as.character(.bfa_canonical_biofluid(meta[[biofluid_col]]))
  for (lev in unique(bio)) {
    idx <- which(bio == lev)
    if (!(lev %in% colnames(fit$effects))) next
    common <- intersect(colnames(log_expr), rownames(fit$effects))
    log_expr[idx, common] <- sweep(log_expr[idx, common, drop = FALSE],
                                   2L, fit$effects[common, lev], "-")
  }
  center <- rowMeans(log_expr[, piv, drop = FALSE], na.rm = TRUE)
  z <- sweep(log_expr[, panel, drop = FALSE], 1L, center, "-")
  out <- .bfa_weighted_panel_sum(z, panel, feature_weights)
  attr(out, "n_pivot_present") <- length(piv)
  attr(out, "pivot_source") <- fit$pivot_source
  out
}

#' Score biofluid-residualized ALR
#'
#' @inheritParams score_biofluid_stratified_rclr
#' @param train_expr Optional training matrix when `fit` is NULL.
#' @param train_meta Optional training metadata when `fit` is NULL.
#' @param outcome_col Binary disease label column in `train_meta`.
#' @param cohort_col Cohort column in `train_meta`.
#' @export
score_biofluid_residualized_alr <- function(expr, meta, panel,
                                            biofluid_col = "biofluid",
                                            fit = NULL,
                                            feature_weights = NULL,
                                            train_expr = NULL,
                                            train_meta = NULL,
                                            outcome_col = "disease",
                                            cohort_col = "cohort",
                                            pseudocount = NULL, ...) {
  if (is.null(fit)) {
    if (is.null(train_expr)) train_expr <- expr
    if (is.null(train_meta)) train_meta <- meta
    fit <- fit_biofluid_residualized_alr(
      train_expr = train_expr, train_meta = train_meta, panel = panel,
      biofluid_col = biofluid_col, outcome_col = outcome_col,
      cohort_col = cohort_col, pseudocount = pseudocount, ...)
  }
  predict_biofluid_residualized_alr(fit, expr, meta, biofluid_col, feature_weights)
}

#' Score rCLR after excluding blood-cell-derived miRNAs
#'
#' Removes the nominated platelet/erythrocyte miRNAs from the expression matrix
#' before computing the rCLR denominator and before scoring the panel.
#'
#' @param expr Numeric matrix, samples x features.
#' @param meta Data frame aligned to `expr`.
#' @param panel Character vector of panel features.
#' @param biofluid_col Biofluid column retained for a common scorer signature.
#' @param excluded_features Features removed before centering and scoring.
#' @param feature_weights Optional named numeric panel weights.
#' @param pseudocount Optional additive pseudocount before log transform.
#' @param trim_upper,trim_lower Feature trimming fractions for rCLR denominator
#'   selection.
#' @param min_centering_size Minimum denominator size before fallback.
#' @param ... Reserved for scorer-interface compatibility.
#' @export
score_platelet_exclusion_rclr <- function(expr, meta, panel,
                                          biofluid_col = "biofluid",
                                          excluded_features = .bfa_default_blood_cell_mirnas(),
                                          feature_weights = NULL,
                                          pseudocount = NULL,
                                          trim_upper = 0.10,
                                          trim_lower = 0.05,
                                          min_centering_size = 8L, ...) {
  checked <- .bfa_check_expr_meta(expr, meta)
  expr <- checked$expr
  keep_features <- setdiff(colnames(expr), excluded_features)
  removed <- intersect(colnames(expr), excluded_features)
  if (length(keep_features) < min_centering_size) {
    out <- rep(NA_real_, nrow(expr))
    attr(out, "features_excluded") <- removed
    return(out)
  }
  expr2 <- expr[, keep_features, drop = FALSE]
  panel2 <- intersect(setdiff(panel, excluded_features), colnames(expr2))
  if (length(panel2) == 0L) {
    out <- rep(NA_real_, nrow(expr))
    attr(out, "features_excluded") <- removed
    return(out)
  }
  z <- .bfa_global_rclr_matrix(expr2, pseudocount = pseudocount,
                               exclude_features = character(0),
                               trim_upper = trim_upper,
                               trim_lower = trim_lower,
                               min_centering_size = min_centering_size)
  out <- .bfa_weighted_panel_sum(z, panel2, feature_weights)
  attr(out, "features_excluded") <- removed
  out
}
