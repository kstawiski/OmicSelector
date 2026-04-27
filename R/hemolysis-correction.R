#' @title Hemolysis-Aware Corrections for Biomarker Panels
#'
#' @description
#' Utilities for simple deployment-facing hemolysis mitigations used in
#' biomarker-panel robustness benchmarks:
#' \itemize{
#'   \item a Blondal-style proxy pre-filter tuned on training controls only;
#'   \item a weighted centered log-ratio (CLR) transform that down-weights
#'     hemolysis-sensitive features while preserving within-sample centering.
#' }
#'
#' The default proxy uses the eLife-14 panel overlap markers
#' `miR-92a-3p / miR-320d` because canonical `miR-23a / miR-451a` are not
#' available in the frozen 14-miRNA panel.
#'
#' @name hemolysis-correction
NULL


#' @keywords internal
.paper1_hemolysis_coefficients <- function() {
  stats::setNames(
    c(0.180, 0.121, 0.085, 0.039),
    c("MIMAT0000092", "MIMAT0006764", "MIMAT0005793", "MIMAT0000451")
  )
}


#' @keywords internal
.paper1_resolve_feature_names <- function(x, feature_names = NULL, original_names = NULL) {
  if (!is.null(feature_names)) {
    if (length(feature_names) != ncol(x)) {
      stop("feature_names must have length ncol(x).", call. = FALSE)
    }
    return(as.character(feature_names))
  }

  if (!is.null(original_names)) {
    if (length(original_names) != ncol(x)) {
      stop("original_names must have length ncol(x).", call. = FALSE)
    }
    return(as.character(original_names))
  }

  if (!is.null(colnames(x))) {
    return(as.character(colnames(x)))
  }

  paste0("feature_", seq_len(ncol(x)))
}


#' @keywords internal
.paper1_match_feature_weights <- function(feature_names, weights = NULL, sensitivity = NULL) {
  feature_names <- as.character(feature_names)

  if (!is.null(weights)) {
    if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
      stop("weights must be a positive numeric vector.", call. = FALSE)
    }
    if (is.null(names(weights))) {
      if (length(weights) != length(feature_names)) {
        stop("Unnamed weights must have length equal to ncol(x).", call. = FALSE)
      }
      return(as.numeric(weights))
    }
    out <- rep(1, length(feature_names))
    idx <- match(feature_names, names(weights))
    keep <- !is.na(idx)
    out[keep] <- as.numeric(weights[idx[keep]])
    if (!any(keep)) {
      stop("Named weights do not overlap feature_names.", call. = FALSE)
    }
    return(out)
  }

  if (is.null(sensitivity)) {
    sensitivity <- .paper1_hemolysis_coefficients()
  }
  if (!is.numeric(sensitivity) || is.null(names(sensitivity))) {
    stop("sensitivity must be a named numeric vector.", call. = FALSE)
  }
  if (any(!is.finite(sensitivity)) || any(sensitivity <= 0)) {
    stop("sensitivity coefficients must be positive and finite.", call. = FALSE)
  }

  idx <- match(feature_names, names(sensitivity))
  if (!any(!is.na(idx))) {
    stop(
      "No overlap between feature_names and sensitivity coefficients; ",
      "provide matching colnames() or explicit weights.",
      call. = FALSE
    )
  }

  matched <- sensitivity[idx[!is.na(idx)]]
  min_coeff <- min(matched)
  out <- rep(1, length(feature_names))
  out[!is.na(idx)] <- pmin(1, min_coeff / sensitivity[idx[!is.na(idx)]])
  as.numeric(out)
}


#' @keywords internal
.paper1_lexicographic_gt <- function(lhs, rhs) {
  for (i in seq_along(lhs)) {
    if (lhs[[i]] > rhs[[i]]) {
      return(TRUE)
    }
    if (lhs[[i]] < rhs[[i]]) {
      return(FALSE)
    }
  }
  FALSE
}


#' @keywords internal
.paper1_proxy_indices <- function(feature_names, numerator, denominator) {
  numerator <- as.character(numerator)
  denominator <- as.character(denominator)
  num_idx <- match(numerator, feature_names)
  den_idx <- match(denominator, feature_names)
  if (is.na(num_idx)) {
    stop("numerator feature '", numerator, "' is not present.", call. = FALSE)
  }
  if (is.na(den_idx)) {
    stop("denominator feature '", denominator, "' is not present.", call. = FALSE)
  }
  c(num_idx, den_idx)
}


#' @keywords internal
.paper1_pseudo_linear_scale_logged <- function(x_log, factor, floor_log2 = -8) {
  x64 <- matrix(as.numeric(x_log), nrow = nrow(x_log), ncol = ncol(x_log))
  factor64 <- matrix(as.numeric(factor), nrow = nrow(x64), ncol = ncol(x64))
  out <- x64
  mask <- abs(factor64 - 1) > 1e-12
  if (!any(mask)) {
    return(matrix(as.numeric(x_log), nrow = nrow(x_log), ncol = ncol(x_log)))
  }
  shifted <- 2^(x64[mask] - floor_log2) - 1
  scaled <- shifted * factor64[mask]
  out[mask] <- log2(scaled + 1) + floor_log2
  out
}


#' @keywords internal
.paper1_simulate_hemolysis_proxy <- function(x, feature_names, strength, sensitivity, seed, floor_log2 = -8) {
  set.seed(as.integer(seed))
  x <- .paper1_require_matrix(x, arg = "x")
  weights <- rep(0, length(feature_names))
  idx <- match(feature_names, names(sensitivity))
  if (any(!is.na(idx))) {
    matched <- sensitivity[idx[!is.na(idx)]]
    weights[!is.na(idx)] <- matched / max(matched)
  }
  h <- stats::runif(nrow(x), min = 0, max = 1)
  factor <- matrix(1, nrow = nrow(x), ncol = ncol(x))
  if (any(weights > 0)) {
    factor <- factor + strength * tcrossprod(h, weights)
  }
  out <- .paper1_pseudo_linear_scale_logged(x, factor, floor_log2 = floor_log2)
  colnames(out) <- feature_names
  out
}


#' @title Hemolysis Proxy Score
#'
#' @description
#' Computes a log-scale hemolysis proxy score for panel robustness checks. By
#' default this is `miR-92a-3p - miR-320d`.
#'
#' @param x Numeric vector or matrix with samples in rows and features in
#'   columns.
#' @param numerator Name of the numerator feature. Defaults to `"MIMAT0000092"`.
#' @param denominator Name of the denominator feature. Defaults to
#'   `"MIMAT0006764"`.
#' @param feature_names Optional feature names. If omitted, `names(x)` or
#'   `colnames(x)` are used.
#'
#' @return Numeric vector of per-sample proxy scores, or a scalar for vector
#'   input.
#'
#' @export
hemolysis_proxy_score <- function(x,
                                  numerator = "MIMAT0000092",
                                  denominator = "MIMAT0006764",
                                  feature_names = NULL) {
  was_vector <- is.vector(x)
  original_names <- if (is.vector(x)) names(x) else colnames(x)
  x <- .paper1_require_matrix(x, arg = "x")
  feature_names <- .paper1_resolve_feature_names(x, feature_names, original_names)
  idx <- .paper1_proxy_indices(feature_names, numerator, denominator)
  scores <- x[, idx[[1L]]] - x[, idx[[2L]]]
  if (was_vector) {
    return(as.numeric(scores[[1L]]))
  }
  as.numeric(scores)
}


#' @title Fit a Hemolysis Pre-Filter
#'
#' @description
#' Tunes a Blondal-style sample rejection threshold using training controls
#' only. The threshold is selected by contrasting clean controls against
#' synthetic proxy-hemolyzed controls and maximizing balanced accuracy.
#'
#' @param X_train_controls Numeric matrix of healthy-control samples only.
#' @param numerator Name of the numerator feature. Defaults to `"MIMAT0000092"`.
#' @param denominator Name of the denominator feature. Defaults to
#'   `"MIMAT0006764"`.
#' @param feature_names Optional feature names. If omitted, `colnames()` are
#'   used.
#' @param sensitivity Named positive numeric vector of hemolysis-sensitivity
#'   coefficients. Defaults to the locked eLife-14 proxy coefficients.
#' @param synthetic_strengths Proxy perturbation strengths used during tuning.
#' @param synthetic_repeats Number of synthetic repeats per strength.
#' @param threshold_grid Quantiles used to generate threshold candidates.
#' @param seed Integer seed for synthetic perturbations.
#' @param floor_log2 Pseudo-linear floor used by the synthetic perturbation.
#'
#' @return An object of class `"omicselector_hemolysis_prefilter"`.
#'
#' @export
fit_hemolysis_prefilter <- function(X_train_controls,
                                    numerator = "MIMAT0000092",
                                    denominator = "MIMAT0006764",
                                    feature_names = NULL,
                                    sensitivity = .paper1_hemolysis_coefficients(),
                                    synthetic_strengths = c(0.5, 1, 2),
                                    synthetic_repeats = 3L,
                                    threshold_grid = seq(0.5, 0.995, length.out = 100),
                                    seed = 42L,
                                    floor_log2 = -8) {
  original_names <- colnames(X_train_controls)
  x <- .paper1_require_matrix(X_train_controls, arg = "X_train_controls")
  feature_names <- .paper1_resolve_feature_names(x, feature_names, original_names)
  colnames(x) <- feature_names

  clean_scores <- hemolysis_proxy_score(
    x,
    numerator = numerator,
    denominator = denominator,
    feature_names = feature_names
  )

  synth_scores <- vector("list", length = as.integer(synthetic_repeats) * length(synthetic_strengths))
  k <- 1L
  for (rep_idx in seq_len(as.integer(synthetic_repeats))) {
    for (strength in synthetic_strengths) {
      synth <- .paper1_simulate_hemolysis_proxy(
        x,
        feature_names = feature_names,
        strength = strength,
        sensitivity = sensitivity,
        seed = as.integer(seed) + rep_idx + round(1000 * strength),
        floor_log2 = floor_log2
      )
      synth_scores[[k]] <- hemolysis_proxy_score(
        synth,
        numerator = numerator,
        denominator = denominator,
        feature_names = feature_names
      )
      k <- k + 1L
    }
  }
  positive_scores <- unlist(synth_scores, use.names = FALSE)
  all_scores <- c(clean_scores, positive_scores)
  labels <- c(rep(0L, length(clean_scores)), rep(1L, length(positive_scores)))
  thresholds <- unique(as.numeric(stats::quantile(all_scores, probs = threshold_grid, na.rm = TRUE, names = FALSE)))

  best <- NULL
  for (threshold in thresholds) {
    pred <- as.integer(all_scores > threshold)
    tp <- sum(pred == 1L & labels == 1L)
    tn <- sum(pred == 0L & labels == 0L)
    fp <- sum(pred == 1L & labels == 0L)
    fn <- sum(pred == 0L & labels == 1L)
    sens <- tp / max(tp + fn, 1L)
    spec <- tn / max(tn + fp, 1L)
    clean_reject <- fp / max(fp + tn, 1L)
    bal_acc <- 0.5 * (sens + spec)
    candidate <- c(bal_acc, spec, sens, -clean_reject, threshold)
    if (is.null(best) || .paper1_lexicographic_gt(candidate, best$candidate)) {
      best <- list(
        candidate = candidate,
        threshold = threshold,
        balanced_accuracy = bal_acc,
        sensitivity = sens,
        specificity = spec,
        clean_reject_rate = clean_reject
      )
    }
  }

  structure(
    list(
      threshold = best$threshold,
      numerator = numerator,
      denominator = denominator,
      feature_names = feature_names,
      sensitivity = sensitivity,
      floor_log2 = floor_log2,
      synthetic_strengths = synthetic_strengths,
      synthetic_repeats = as.integer(synthetic_repeats),
      seed = as.integer(seed),
      tuning = best
    ),
    class = "omicselector_hemolysis_prefilter"
  )
}


#' @title Apply a Hemolysis Pre-Filter
#'
#' @description
#' Scores new samples with a fit created by [fit_hemolysis_prefilter()] and
#' returns a keep/reject decision for each sample.
#'
#' @param fit A fitted `"omicselector_hemolysis_prefilter"` object.
#' @param X_new Numeric vector or matrix of new samples.
#' @param feature_names Optional feature names for `X_new`.
#'
#' @return A `data.frame` with `hemolysis_score`, `keep`, `rejected`, and
#'   `threshold`.
#'
#' @export
apply_hemolysis_prefilter <- function(fit, X_new, feature_names = NULL) {
  if (!inherits(fit, "omicselector_hemolysis_prefilter")) {
    stop("fit must be an object returned by fit_hemolysis_prefilter().", call. = FALSE)
  }

  original_names <- if (is.vector(X_new)) names(X_new) else colnames(X_new)
  x_new <- .paper1_require_matrix(X_new, arg = "X_new")
  resolved_names <- .paper1_resolve_feature_names(x_new, feature_names, original_names)
  colnames(x_new) <- resolved_names
  scores <- hemolysis_proxy_score(
    x_new,
    numerator = fit$numerator,
    denominator = fit$denominator,
    feature_names = resolved_names
  )
  keep <- as.logical(scores <= fit$threshold)
  data.frame(
    hemolysis_score = as.numeric(scores),
    keep = keep,
    rejected = !keep,
    threshold = rep(fit$threshold, length(scores))
  )
}


#' @title Weighted Centered Log-Ratio Transformation
#'
#' @description
#' Computes a weighted within-sample CLR transform. The weighted center uses
#' down-weighted hemolysis-sensitive features, and the transformed coordinates
#' can optionally be shrunk by the same weights before classification.
#'
#' @param x Numeric vector or matrix with samples in rows and features in
#'   columns.
#' @param weights Optional positive numeric vector of feature weights. If named,
#'   names are matched to `feature_names`.
#' @param sensitivity Optional named positive numeric vector of hemolysis
#'   coefficients. Used only when `weights = NULL`.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param pseudocount Small positive constant added before logging raw inputs.
#' @param shrink_features Logical; if `TRUE`, multiply the centered coordinates
#'   by their feature weights after centering.
#' @param feature_names Optional feature names.
#'
#' @return Numeric vector or matrix of weighted CLR features.
#'
#' @export
weighted_clr_transform <- function(x,
                                   weights = NULL,
                                   sensitivity = .paper1_hemolysis_coefficients(),
                                   input_scale = c("log", "raw"),
                                   pseudocount = 1e-8,
                                   shrink_features = TRUE,
                                   feature_names = NULL) {
  input_scale <- match.arg(input_scale)
  was_vector <- is.vector(x)
  original_names <- if (is.vector(x)) names(x) else colnames(x)
  x <- .paper1_require_matrix(x, arg = "x")
  feature_names <- .paper1_resolve_feature_names(x, feature_names, original_names)
  colnames(x) <- feature_names

  if (identical(input_scale, "raw")) {
    if (any(x <= 0, na.rm = TRUE)) {
      stop("Raw weighted CLR inputs must be strictly positive.", call. = FALSE)
    }
    x <- log(x + pseudocount)
  }

  weight_vec <- .paper1_match_feature_weights(
    feature_names = feature_names,
    weights = weights,
    sensitivity = sensitivity
  )
  weighted_center <- rowSums(sweep(x, 2, weight_vec, `*`)) / sum(weight_vec)
  centered <- sweep(x, 1, weighted_center, `-`)
  if (isTRUE(shrink_features)) {
    centered <- sweep(centered, 2, weight_vec, `*`)
  }

  if (was_vector) {
    return(as.numeric(centered[1, ]))
  }
  centered
}


#' @title Train a Weighted CLR + MLP Classifier
#'
#' @description
#' Fits a hemolysis-aware weighted-CLR analogue of [train_clr_mlp()] using the
#' same backends and normalization scheme as the dense CoDA baseline.
#'
#' @inheritParams train_clr_mlp
#' @param weights Optional positive numeric vector of feature weights.
#' @param sensitivity Optional named positive numeric vector of hemolysis
#'   coefficients used when `weights = NULL`.
#' @param shrink_features Logical; whether to shrink centered coordinates by
#'   feature weight after weighted centering.
#'
#' @return An object of class `"omicselector_weighted_clr_mlp"`.
#'
#' @export
train_weighted_clr_mlp <- function(X_train, y_train,
                                   weights = NULL,
                                   sensitivity = .paper1_hemolysis_coefficients(),
                                   shrink_features = TRUE,
                                   epochs = 30L,
                                   lr = 0.001,
                                   batch_size = 256L,
                                   class_weight = TRUE,
                                   device = NULL,
                                   verbose = TRUE,
                                   seed = 42L,
                                   positive = NULL,
                                   input_scale = c("log", "raw"),
                                   backend = c("auto", "torch", "nnet", "glm"),
                                   hidden_units = c(152L, 136L),
                                   dropout = 0.3) {
  input_scale <- match.arg(input_scale)
  backend <- match.arg(backend)
  train <- .paper1_prepare_train_inputs(X_train, y_train, positive = positive)
  device <- .paper1_resolve_device(device)

  if (identical(backend, "auto")) {
    if (requireNamespace("torch", quietly = TRUE)) {
      backend <- "torch"
    } else if (requireNamespace("nnet", quietly = TRUE)) {
      backend <- "nnet"
    } else {
      backend <- "glm"
    }
  }

  wclr_train <- weighted_clr_transform(
    train$x,
    weights = weights,
    sensitivity = sensitivity,
    input_scale = input_scale,
    shrink_features = shrink_features,
    feature_names = colnames(train$x)
  )
  normalized <- .paper1_normalize_train_test(wclr_train)

  fit <- switch(
    backend,
    torch = .fit_clr_mlp_torch(
      x_train = normalized$train,
      y_train = train$y,
      epochs = epochs,
      lr = lr,
      batch_size = batch_size,
      class_weight = class_weight,
      device = device,
      verbose = verbose,
      seed = seed,
      hidden_units = hidden_units,
      dropout = dropout
    ),
    nnet = .fit_clr_mlp_nnet(
      x_train = normalized$train,
      y_train = train$y,
      class_names = train$class_names,
      positive = train$positive,
      epochs = epochs,
      class_weight = class_weight,
      verbose = verbose
    ),
    glm = .fit_clr_mlp_glm(
      x_train = normalized$train,
      y_train = train$y,
      class_weight = class_weight
    )
  )

  structure(
    c(
      fit,
      list(
        class_names = train$class_names,
        positive = train$positive,
        feature_names = colnames(train$x),
        weighted_clr_stats = list(mean = normalized$mean, sd = normalized$sd),
        input_scale = input_scale,
        seed = as.integer(seed),
        hidden_units = as.integer(hidden_units),
        dropout = dropout,
        weights = weights,
        sensitivity = sensitivity,
        shrink_features = shrink_features
      )
    ),
    class = "omicselector_weighted_clr_mlp"
  )
}


#' @title Predict with a Weighted CLR + MLP Classifier
#'
#' @description
#' Generates positive-class probabilities for new samples from a fit created by
#' [train_weighted_clr_mlp()].
#'
#' @param fit A fitted `"omicselector_weighted_clr_mlp"` object.
#' @param X_new Numeric matrix of new samples.
#'
#' @return Numeric vector of positive-class probabilities.
#'
#' @export
predict_weighted_clr_mlp <- function(fit, X_new) {
  if (!inherits(fit, "omicselector_weighted_clr_mlp")) {
    stop("fit must be an object returned by train_weighted_clr_mlp().", call. = FALSE)
  }

  X_new <- .paper1_prepare_new_data(
    X_new,
    n_features = length(fit$feature_names),
    feature_names = fit$feature_names
  )
  wclr_new <- weighted_clr_transform(
    X_new,
    weights = fit$weights,
    sensitivity = fit$sensitivity,
    input_scale = fit$input_scale,
    shrink_features = fit$shrink_features,
    feature_names = fit$feature_names
  )
  if (!is.na(fit$weighted_clr_stats$sd) && fit$weighted_clr_stats$sd > 0) {
    wclr_new <- (wclr_new - fit$weighted_clr_stats$mean) / fit$weighted_clr_stats$sd
  }

  switch(
    fit$backend,
    torch = .predict_clr_mlp_torch(fit, wclr_new),
    nnet = .predict_clr_mlp_nnet(fit, wclr_new),
    glm = .predict_clr_mlp_glm(fit, wclr_new)
  )
}
