#' @title Compute Binary Clinical Operating Points
#'
#' @description
#' Summarizes sensitivity and specificity at prespecified targets. Larger scores
#' are treated as more positive.
#'
#' @param y Binary outcome vector. The second factor level is treated as positive
#'   unless \code{positive} is supplied.
#' @param score Numeric prediction score.
#' @param specificity_targets Specificity thresholds to evaluate.
#' @param sensitivity_targets Sensitivity thresholds to evaluate.
#' @param positive Optional positive-class label.
#'
#' @return An \code{os_operating_points} object.
#' @export
os_operating_points <- function(y, score,
                                specificity_targets = c(0.90, 0.95, 0.98),
                                sensitivity_targets = c(0.90, 0.95),
                                positive = NULL) {
  y_num <- .os_binary_with_positive(y, positive)
  if (length(score) != length(y_num)) stop("`score` and `y` must have the same length.", call. = FALSE)
  score <- as.numeric(score)
  thresholds <- sort(unique(c(-Inf, score, Inf)))
  curve <- do.call(rbind, lapply(thresholds, function(thr) {
    pred <- as.integer(score >= thr)
    tp <- sum(pred == 1L & y_num == 1L)
    fp <- sum(pred == 1L & y_num == 0L)
    tn <- sum(pred == 0L & y_num == 0L)
    fn <- sum(pred == 0L & y_num == 1L)
    sens <- if ((tp + fn) == 0L) NA_real_ else tp / (tp + fn)
    spec <- if ((tn + fp) == 0L) NA_real_ else tn / (tn + fp)
    ppv <- if ((tp + fp) == 0L) NA_real_ else tp / (tp + fp)
    npv <- if ((tn + fn) == 0L) NA_real_ else tn / (tn + fn)
    data.frame(
      threshold = thr,
      sensitivity = sens,
      specificity = spec,
      ppv = ppv,
      npv = npv,
      stringsAsFactors = FALSE
    )
  }))

  by_specificity <- do.call(rbind, lapply(specificity_targets, function(target) {
    candidates <- curve[is.finite(curve$specificity) & curve$specificity >= target, , drop = FALSE]
    .os_best_operating_row(candidates, target, "specificity")
  }))
  by_sensitivity <- do.call(rbind, lapply(sensitivity_targets, function(target) {
    candidates <- curve[is.finite(curve$sensitivity) & curve$sensitivity >= target, , drop = FALSE]
    .os_best_operating_row(candidates, target, "sensitivity")
  }))
  out <- list(
    curve = curve,
    by_specificity = by_specificity,
    by_sensitivity = by_sensitivity,
    prevalence = mean(y_num == 1L),
    n_positive = sum(y_num == 1L),
    n_negative = sum(y_num == 0L)
  )
  class(out) <- c("os_operating_points", "list")
  out
}

#' @title Compute Apparent or Cross-Fitted Calibrated Brier Score
#'
#' @param y Binary outcome vector.
#' @param score Numeric prediction score.
#' @param folds Optional list of held-out row indices for cross-fitting Platt calibration.
#' @param positive Optional positive-class label.
#'
#' @return An \code{os_calibrated_brier} object.
#' @export
os_calibrated_brier <- function(y, score, folds = NULL, positive = NULL) {
  y_num <- .os_binary_with_positive(y, positive)
  if (length(score) != length(y_num)) stop("`score` and `y` must have the same length.", call. = FALSE)
  score <- as.numeric(score)
  prob <- rep(NA_real_, length(score))
  if (is.null(folds)) {
    prob <- .os_platt_predict(y_num, score, seq_along(score), seq_along(score))
    mode <- "apparent"
  } else {
    for (idx in folds) {
      train <- setdiff(seq_along(score), idx)
      prob[idx] <- .os_platt_predict(y_num, score, train, idx)
    }
    mode <- "cross_fitted"
  }
  prob <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  out <- list(
    brier = mean((prob - y_num)^2),
    probabilities = prob,
    mode = mode
  )
  class(out) <- c("os_calibrated_brier", "list")
  out
}

#' @title Fail-Closed Operating-Point Gate
#'
#' @param operating_points Object from \code{os_operating_points}.
#' @param identifiability_status Terminal status from
#'   \code{os_identifiability_gate()} or an equivalent provenance review.
#' @param min_specificity Minimum specificity required.
#' @param min_sensitivity Minimum sensitivity required at \code{min_specificity}.
#'
#' @return A one-row data frame with gate status.
#' @export
os_operating_point_gate <- function(operating_points, identifiability_status,
                                    min_specificity = 0.95,
                                    min_sensitivity = 0.50) {
  if (!inherits(operating_points, "os_operating_points")) {
    stop("`operating_points` must be an `os_operating_points` object.", call. = FALSE)
  }
  op <- operating_points$by_specificity
  op <- op[op$target >= min_specificity, , drop = FALSE]
  op <- op[order(op$target), , drop = FALSE]
  row <- if (nrow(op) == 0L) NULL else op[1L, , drop = FALSE]
  if (!identical(identifiability_status, "PASS")) {
    status <- "FAIL-CLOSED"
    label <- "not_threshold_eligible_after_identifiability_failure"
  } else if (!is.null(row) && isTRUE(row$sensitivity >= min_sensitivity) && isTRUE(row$specificity >= min_specificity)) {
    status <- "PASS"
    label <- "threshold_eligible"
  } else {
    status <- "FAIL"
    label <- "operating_point_below_requirement"
  }
  data.frame(
    status = status,
    claim_label = label,
    min_specificity = min_specificity,
    min_sensitivity = min_sensitivity,
    observed_specificity = if (is.null(row)) NA_real_ else row$specificity,
    observed_sensitivity = if (is.null(row)) NA_real_ else row$sensitivity,
    stringsAsFactors = FALSE
  )
}

.os_best_operating_row <- function(candidates, target, target_type) {
  if (nrow(candidates) == 0L) {
    return(data.frame(
      target_type = target_type,
      target = target,
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      ppv = NA_real_,
      npv = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  if (target_type == "specificity") {
    candidates <- candidates[order(-candidates$sensitivity, -candidates$specificity), , drop = FALSE]
  } else {
    candidates <- candidates[order(-candidates$specificity, -candidates$sensitivity), , drop = FALSE]
  }
  out <- candidates[1L, c("threshold", "sensitivity", "specificity", "ppv", "npv"), drop = FALSE]
  data.frame(target_type = target_type, target = target, out, row.names = NULL, stringsAsFactors = FALSE)
}

.os_platt_predict <- function(y_num, score, train, test) {
  tryCatch({
    if (length(unique(y_num[train])) < 2L) stop("single-class calibration fold", call. = FALSE)
    cal_df <- data.frame(y = y_num[train], score = score[train])
    fit <- suppressWarnings(stats::glm(y ~ score, data = cal_df, family = stats::binomial()))
    as.numeric(stats::predict(fit, newdata = data.frame(score = score[test]), type = "response"))
  }, error = function(e) {
    rep(mean(y_num[train]), length(test))
  })
}
