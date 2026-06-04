#' @title Fit a Top-k Oriented-Pair Panel Scorer
#'
#' @description
#' Fits a deterministic top-k oriented-pair scorer on a locked panel. Pair
#' orientation is selected from the training data so that larger vote fractions
#' indicate the positive class. This is a lightweight, auditable k-TSP-style
#' scorer for small biomarker panels; it is not a full switchBox replacement.
#'
#' @param X Numeric matrix or data frame with samples in rows and features in columns.
#' @param y Binary outcome vector. The second factor level is treated as positive.
#' @param k Number of oriented pairs to retain.
#'
#' @return An object of class \code{"os_ktsp_model"}.
#' @export
os_ktsp_fit <- function(X, y, k = 11L) {
  X <- .os_as_numeric_matrix(X)
  y_num <- .os_binary_y(y)
  if (ncol(X) < 2L) stop("`X` must contain at least two features.", call. = FALSE)
  k <- as.integer(k)
  if (is.na(k) || k < 1L) stop("`k` must be a positive integer.", call. = FALSE)

  pair_rows <- list()
  idx <- 1L
  for (a in seq_len(ncol(X) - 1L)) {
    for (b in (a + 1L):ncol(X)) {
      ctrl_rate <- mean(X[y_num == 0L, a] < X[y_num == 0L, b], na.rm = TRUE)
      case_rate <- mean(X[y_num == 1L, a] < X[y_num == 1L, b], na.rm = TRUE)
      if (!is.finite(ctrl_rate)) ctrl_rate <- 0
      if (!is.finite(case_rate)) case_rate <- 0
      if (case_rate >= ctrl_rate) {
        pair_rows[[idx]] <- data.frame(
          feature_a = colnames(X)[a],
          feature_b = colnames(X)[b],
          orientation = "<",
          score = case_rate - ctrl_rate,
          stringsAsFactors = FALSE
        )
      } else {
        pair_rows[[idx]] <- data.frame(
          feature_a = colnames(X)[b],
          feature_b = colnames(X)[a],
          orientation = "<",
          score = ctrl_rate - case_rate,
          stringsAsFactors = FALSE
        )
      }
      idx <- idx + 1L
    }
  }
  pairs <- do.call(rbind, pair_rows)
  pairs <- pairs[order(-pairs$score, pairs$feature_a, pairs$feature_b), , drop = FALSE]
  pairs <- utils::head(pairs, min(k, nrow(pairs)))
  out <- list(pairs = pairs, features = colnames(X), k = nrow(pairs))
  class(out) <- c("os_ktsp_model", "list")
  out
}

#' @title Predict from a Top-k Oriented-Pair Panel Scorer
#'
#' @param object An \code{os_ktsp_model} from \code{os_ktsp_fit}.
#' @param newdata Numeric matrix or data frame.
#' @param ... Unused.
#'
#' @return Numeric vote fraction per sample.
#' @export
predict.os_ktsp_model <- function(object, newdata, ...) {
  X <- .os_as_numeric_matrix(newdata)
  missing <- setdiff(unique(c(object$pairs$feature_a, object$pairs$feature_b)), colnames(X))
  if (length(missing) > 0L) {
    stop("`newdata` is missing fitted pair features: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  # cbind (not vapply) keeps an nrow x npairs matrix even for a single-row
  # newdata, so the score is always the per-row vote fraction; vapply collapses
  # a 1-row newdata to a length-npairs vector and would return the raw votes.
  votes <- do.call(cbind, lapply(seq_len(nrow(object$pairs)), function(i) {
    a <- object$pairs$feature_a[i]
    b <- object$pairs$feature_b[i]
    as.numeric(X[, a] < X[, b])
  }))
  as.numeric(rowMeans(votes))
}

#' @title Score a Direction-Split Panel by Within-Sample Ranks
#'
#' @description
#' Computes a singscore-style score from within-sample fractional ranks:
#' mean rank of up features minus mean rank of down features. The split
#' cardinality is carried in the return attributes for null-QC checks.
#'
#' @param X Numeric matrix or data frame.
#' @param up Character vector of up-in-positive feature names.
#' @param down Character vector of down-in-positive feature names.
#'
#' @return Numeric score vector.
#' @export
os_singscore <- function(X, up, down) {
  X <- .os_as_numeric_matrix(X)
  missing <- setdiff(c(up, down), colnames(X))
  if (length(missing) > 0L) {
    stop("`X` is missing split features: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(up) < 1L || length(down) < 1L) {
    stop("Both `up` and `down` must contain at least one feature.", call. = FALSE)
  }
  ranks <- t(apply(X, 1L, function(row) rank(row, ties.method = "average") / length(row)))
  score <- rowMeans(ranks[, up, drop = FALSE]) - rowMeans(ranks[, down, drop = FALSE])
  attr(score, "up_n") <- length(up)
  attr(score, "down_n") <- length(down)
  score
}

#' @title Mahalanobis Anomaly Score from a Reference Set
#'
#' @param X Numeric matrix or data frame.
#' @param reference_rows Integer, logical, or character row selector for the reference set.
#' @param ridge Non-negative diagonal ridge added to the covariance matrix.
#'
#' @return Numeric squared Mahalanobis distance per row.
#' @export
os_mahalanobis_score <- function(X, reference_rows, ridge = 1e-6) {
  X <- .os_as_numeric_matrix(X)
  ref_idx <- .os_row_index(reference_rows, X)
  if (length(ref_idx) < 2L) stop("`reference_rows` must select at least two samples.", call. = FALSE)
  X_imp <- .os_impute_by_reference(X, ref_idx)
  center <- colMeans(X_imp[ref_idx, , drop = FALSE])
  cov_mat <- stats::cov(X_imp[ref_idx, , drop = FALSE])
  cov_mat <- cov_mat + diag(ridge, ncol(cov_mat))
  stats::mahalanobis(X_imp, center = center, cov = cov_mat)
}

#' @title Conformal Healthy-Reference Anomaly Score
#'
#' @description
#' Computes a split-conformal anomaly score using Mahalanobis distances from
#' reference-training controls and calibration controls. The returned score is
#' \code{1 - conformal p-value}; larger values are more anomalous.
#'
#' @param X Numeric matrix or data frame.
#' @param train_rows Reference-training row selector.
#' @param calibration_rows Calibration row selector.
#' @param ridge Non-negative diagonal ridge for covariance.
#'
#' @return Numeric score vector.
#' @export
os_conformal_anomaly_score <- function(X, train_rows, calibration_rows, ridge = 1e-6) {
  X <- .os_as_numeric_matrix(X)
  tr <- .os_row_index(train_rows, X)
  cal <- .os_row_index(calibration_rows, X)
  if (length(tr) < 2L || length(cal) < 1L) {
    stop("Conformal scoring requires >=2 training and >=1 calibration rows.", call. = FALSE)
  }
  X_imp <- .os_impute_by_reference(X, tr)
  center <- colMeans(X_imp[tr, , drop = FALSE])
  cov_mat <- stats::cov(X_imp[tr, , drop = FALSE]) + diag(ridge, ncol(X_imp))
  dist_all <- stats::mahalanobis(X_imp, center = center, cov = cov_mat)
  dist_cal <- dist_all[cal]
  p <- vapply(dist_all, function(s) {
    (1 + sum(dist_cal >= s)) / (length(dist_cal) + 1)
  }, FUN.VALUE = numeric(1L))
  1 - p
}

#' @title Matched Random-Panel Null Benchmark
#'
#' @description
#' Benchmarks a locked panel statistic against random panels of the same size.
#' The scorer receives \code{X} subset to a feature set and must return one
#' score per sample. Null panels are sampled without replacement and are not
#' allowed to overlap the locked panel.
#'
#' @param X Numeric matrix or data frame with named columns.
#' @param y Binary outcome vector.
#' @param panel_features Character vector of locked-panel feature names.
#' @param scorer Function with signature \code{scorer(X_sub, features)}.
#' @param n_draws Number of random panels.
#' @param feature_pool Candidate null features. Defaults to all non-panel columns.
#' @param seed Random seed.
#' @param metric Currently only \code{"auc"}.
#' @param ... Additional arguments passed to \code{scorer}.
#'
#' @return Object of class \code{"os_panel_null"}.
#' @export
os_panel_null_benchmark <- function(X, y, panel_features, scorer,
                                    n_draws = 1000L, feature_pool = NULL,
                                    seed = 1L, metric = "auc", ...) {
  X <- .os_as_numeric_matrix(X)
  y_num <- .os_binary_y(y)
  panel_features <- as.character(panel_features)
  missing <- setdiff(panel_features, colnames(X))
  if (length(missing) > 0L) stop("Missing panel features: ", paste(missing, collapse = ", "), call. = FALSE)
  if (metric != "auc") stop("Only metric='auc' is currently supported.", call. = FALSE)
  if (is.null(feature_pool)) feature_pool <- setdiff(colnames(X), panel_features)
  feature_pool <- setdiff(as.character(feature_pool), panel_features)
  if (length(feature_pool) < length(panel_features)) {
    stop("`feature_pool` must contain at least as many non-panel features as the locked panel.", call. = FALSE)
  }
  n_draws <- as.integer(n_draws)
  if (is.na(n_draws) || n_draws < 1L) stop("`n_draws` must be a positive integer.", call. = FALSE)

  observed_score <- scorer(X[, panel_features, drop = FALSE], panel_features, ...)
  observed_stat <- .auc_fast(y_num, observed_score)
  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(rng_state)) .Random.seed <<- rng_state
  }, add = TRUE)
  set.seed(seed)
  null_rows <- vector("list", n_draws)
  null_stats <- numeric(n_draws)
  for (i in seq_len(n_draws)) {
    feats <- sample(feature_pool, length(panel_features), replace = FALSE)
    score <- scorer(X[, feats, drop = FALSE], feats, ...)
    null_stats[i] <- .auc_fast(y_num, score)
    null_rows[[i]] <- data.frame(
      draw = i,
      features = paste(feats, collapse = ","),
      stat = null_stats[i],
      overlaps_locked_panel = any(feats %in% panel_features),
      stringsAsFactors = FALSE
    )
  }
  out <- list(
    metric = metric,
    panel_features = panel_features,
    observed_stat = observed_stat,
    null_stats = null_stats,
    null_draws = do.call(rbind, null_rows),
    null_p95 = unname(stats::quantile(null_stats, 0.95, na.rm = TRUE)),
    null_p99 = unname(stats::quantile(null_stats, 0.99, na.rm = TRUE)),
    permutation_p = (1 + sum(null_stats >= observed_stat, na.rm = TRUE)) /
      (sum(is.finite(null_stats)) + 1),
    n_draws = n_draws
  )
  class(out) <- c("os_panel_null", "list")
  out
}

#' @title Validate a Random-Panel Null Benchmark
#'
#' @param null_result Object from \code{os_panel_null_benchmark}.
#' @param min_draws Minimum accepted number of null draws.
#' @param expected_panel_size Expected locked-panel size.
#'
#' @return A data frame with QC checks and pass/fail values.
#' @export
os_null_qc <- function(null_result, min_draws = 1000L, expected_panel_size = NULL) {
  if (!inherits(null_result, "os_panel_null")) {
    stop("`null_result` must be an `os_panel_null` object.", call. = FALSE)
  }
  panel_size <- length(null_result$panel_features)
  checks <- data.frame(
    check = c("minimum_draws", "finite_nulls", "no_panel_overlap", "panel_size"),
    pass = c(
      null_result$n_draws >= min_draws,
      all(is.finite(null_result$null_stats)),
      !any(null_result$null_draws$overlaps_locked_panel),
      is.null(expected_panel_size) || panel_size == expected_panel_size
    ),
    value = c(
      as.character(null_result$n_draws),
      paste0(sum(is.finite(null_result$null_stats)), "/", length(null_result$null_stats)),
      as.character(any(null_result$null_draws$overlaps_locked_panel)),
      as.character(panel_size)
    ),
    stringsAsFactors = FALSE
  )
  checks
}

#' @title Claim Gate for Panel-Null Benchmarks
#'
#' @description
#' Applies a fail-closed panel-claim policy: if the disease contrast has not
#' passed the provenance/identifiability review, no method can receive a
#' positive primary label. If identifiability passed, p-values are Holm adjusted
#' and a method passes only when adjusted p-value, null95 margin, and null QC
#' all satisfy the prespecified thresholds.
#'
#' @param results Data frame with columns \code{method_id}, \code{p_value},
#'   \code{observed_stat}, \code{null_p95}, and \code{null_qc_pass}.
#' @param identifiability_status One of \code{"PASS"}, \code{"FAIL"}, or
#'   \code{"INCONCLUSIVE"}.
#' @param required_margin Numeric margin required over null95.
#' @param alpha Family-wise alpha for Holm-adjusted p-values.
#'
#' @return Data frame with adjusted p-values, margins, and labels.
#' @export
os_claim_gate <- function(results, identifiability_status, required_margin = 0.02, alpha = 0.05) {
  required <- c("method_id", "p_value", "observed_stat", "null_p95", "null_qc_pass")
  missing <- setdiff(required, names(results))
  if (length(missing) > 0L) stop("`results` missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  out <- results
  out$margin_over_null95 <- out$observed_stat - out$null_p95
  out$holm_p <- stats::p.adjust(out$p_value, method = "holm")
  if (!identical(identifiability_status, "PASS")) {
    out$claim_label <- "not_identifiable_descriptive_only"
    return(out)
  }
  out$claim_label <- ifelse(
    out$null_qc_pass & out$holm_p <= alpha & out$margin_over_null95 >= required_margin,
    "primary_pass",
    "primary_fail"
  )
  out
}

#' @title Required Margin for Matched-Null Claims
#'
#' @param mde80 Minimum detectable effect at 80 percent power.
#' @param w90 Width of the relevant 90 percent confidence interval.
#' @param floor Minimum floor margin.
#'
#' @return Numeric required margin.
#' @export
os_required_margin <- function(mde80, w90, floor = 0.02) {
  vals <- c(floor, mde80, 0.5 * w90)
  if (any(!is.finite(vals))) stop("All margin inputs must be finite.", call. = FALSE)
  max(vals)
}

#' @title Create an MDE/Margin Record
#'
#' @param n Effective sample count.
#' @param n_clusters Effective cluster count.
#' @param method_ids Character vector defining the fixed family.
#' @param null_draws Number of null draws.
#' @param mde80 Minimum detectable effect at 80 percent power.
#' @param w90 Width of the relevant 90 percent confidence interval.
#' @param floor Minimum floor margin.
#'
#' @return List with class \code{"os_mde_margin_record"}.
#' @export
os_mde_margin_record <- function(n, n_clusters, method_ids, null_draws,
                                 mde80, w90, floor = 0.02) {
  required_margin <- os_required_margin(mde80, w90, floor)
  payload <- list(
    n = as.integer(n),
    n_clusters = as.integer(n_clusters),
    method_ids = as.character(method_ids),
    holm_family_size = length(method_ids),
    null_draws = as.integer(null_draws),
    mde80 = as.numeric(mde80),
    w90 = as.numeric(w90),
    required_margin = required_margin
  )
  payload$hash <- digest::digest(payload, algo = "sha256")
  class(payload) <- c("os_mde_margin_record", "list")
  payload
}

.os_as_numeric_matrix <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (is.null(colnames(X))) colnames(X) <- paste0("F", seq_len(ncol(X)))
  if (is.null(rownames(X))) rownames(X) <- paste0("S", seq_len(nrow(X)))
  X
}

.os_binary_y <- function(y) {
  y_fac <- as.factor(y)
  if (nlevels(y_fac) != 2L) stop("`y` must have exactly two levels.", call. = FALSE)
  as.integer(y_fac) - 1L
}

.os_row_index <- function(rows, X) {
  if (is.logical(rows)) {
    if (length(rows) != nrow(X)) stop("Logical row selector has wrong length.", call. = FALSE)
    return(which(rows))
  }
  if (is.character(rows)) {
    idx <- match(rows, rownames(X))
    if (any(is.na(idx))) stop("Unknown row names in selector.", call. = FALSE)
    return(idx)
  }
  rows <- as.integer(rows)
  if (any(is.na(rows)) || any(rows < 1L) || any(rows > nrow(X))) {
    stop("Row selector contains invalid indices.", call. = FALSE)
  }
  rows
}

.os_impute_by_reference <- function(X, ref_idx) {
  out <- X
  med <- apply(out[ref_idx, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  idx <- which(!is.finite(out), arr.ind = TRUE)
  if (nrow(idx) > 0L) out[idx] <- med[idx[, 2L]]
  out
}
