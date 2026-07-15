#' Paired DeLong SE for single-sample AUC lift
#'
#' Computes the empirical DeLong variance of a paired AUC difference
#' (method score minus baseline score). When `fold` is supplied, the per-fold
#' DeLong variances are aggregated as the variance of the mean fold-level AUC
#' difference, matching the `auc_obs_cv = mean(auc_obs_folds)` mean-of-folds estimand.
#'
#' @param y Binary outcome vector, coded 0/1 or a two-level factor.
#' @param score_method Numeric scores from the candidate method.
#' @param score_baseline Numeric scores from the paired baseline method.
#' @param fold Optional fold identifier. If `NULL`, all observations are treated
#'   as one fold.
#' @param orient Orientation rule. `"median"` mirrors pROC's within-cohort
#'   auto-orientation convention; `"auc"` reports `max(AUC, 1 - AUC)` and is
#'   used by the transfer-family technology-aware engines; `"fixed"` preserves
#'   the supplied score direction and therefore does not inspect held-out
#'   outcomes to choose orientation.
#'
#' @return A list with `lift`, `se`, `variance`, `n_folds`, `auc_method`,
#'   `auc_baseline`, and a per-fold `fold_results` data frame.
#'
#' @references DeLong ER, DeLong DM, Clarke-Pearson DL. (1988) Comparing the
#' areas under two or more correlated receiver operating characteristic curves:
#' a nonparametric approach. \emph{Biometrics} 44:837-845.
#'
#' @export
singlesample_paired_auc_diff_se <- function(y, score_method, score_baseline,
                                      fold = NULL,
                                      orient = c("median", "auc", "fixed")) {
  orient <- match.arg(orient)
  if (is.factor(y)) y <- as.integer(y) - 1L
  y <- as.integer(y)
  score_method <- as.numeric(score_method)
  score_baseline <- as.numeric(score_baseline)
  n <- length(y)
  if (length(score_method) != n || length(score_baseline) != n) {
    stop("y, score_method, and score_baseline must have the same length.",
         call. = FALSE)
  }
  if (is.null(fold)) fold <- rep(1L, n)
  if (length(fold) != n) stop("fold must have length equal to y.", call. = FALSE)

  fold_levels <- sort(unique(fold))
  rows <- lapply(fold_levels, function(f) {
    idx <- which(fold == f)
    if (length(unique(y[idx][!is.na(y[idx])])) < 2L) {
      return(data.frame(
        fold = as.character(f), auc_method = NA_real_,
        auc_baseline = NA_real_, lift = NA_real_, variance = NA_real_,
        n_pos = sum(y[idx] == 1L, na.rm = TRUE),
        n_neg = sum(y[idx] == 0L, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
    res <- .singlesample_delong_var_diff_fold(
      y[idx], score_method[idx], score_baseline[idx], orient = orient)
    data.frame(
      fold = as.character(f), auc_method = res$auc_method,
      auc_baseline = res$auc_baseline,
      lift = res$auc_method - res$auc_baseline,
      variance = res$variance,
      n_pos = res$n_pos, n_neg = res$n_neg,
      stringsAsFactors = FALSE
    )
  })
  fold_results <- do.call(rbind, rows)
  ok <- is.finite(fold_results$variance)
  k <- sum(ok)
  if (k == 0L) {
    return(list(
      lift = NA_real_, se = NA_real_, variance = NA_real_, n_folds = 0L,
      auc_method = NA_real_, auc_baseline = NA_real_,
      fold_results = fold_results
    ))
  }
  variance <- sum(fold_results$variance[ok]) / k^2
  list(
    lift = mean(fold_results$lift[ok]),
    se = sqrt(max(0, variance)),
    variance = max(0, variance),
    n_folds = k,
    auc_method = mean(fold_results$auc_method[ok]),
    auc_baseline = mean(fold_results$auc_baseline[ok]),
    fold_results = fold_results
  )
}

.singlesample_delong_placements <- function(y, score,
                                            orient = c("median", "auc", "fixed")) {
  orient <- match.arg(orient)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  pos <- score[y == 1L]
  neg <- score[y == 0L]
  n_pos <- length(pos)
  n_neg <- length(neg)
  if (n_pos == 0L || n_neg == 0L) {
    return(list(v10 = numeric(0), v01 = numeric(0),
                n_pos = n_pos, n_neg = n_neg, auc = NA_real_))
  }
  cmp <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
  v10 <- rowMeans(cmp)
  v01 <- colMeans(cmp)
  auc <- mean(cmp)
  reflect <- if (orient == "fixed") {
    FALSE
  } else if (orient == "auc") {
    is.finite(auc) && auc < 0.5
  } else {
    stats::median(neg) > stats::median(pos)
  }
  if (isTRUE(reflect)) {
    v10 <- 1 - v10
    v01 <- 1 - v01
    auc <- 1 - auc
  }
  list(v10 = v10, v01 = v01, n_pos = n_pos, n_neg = n_neg, auc = auc)
}

.singlesample_delong_var_diff_fold <- function(y, score_method, score_baseline,
                                         orient = c("median", "auc", "fixed")) {
  orient <- match.arg(orient)
  ok <- is.finite(score_method) & is.finite(score_baseline) & !is.na(y)
  y <- y[ok]
  score_method <- score_method[ok]
  score_baseline <- score_baseline[ok]
  pm <- .singlesample_delong_placements(y, score_method, orient = orient)
  pb <- .singlesample_delong_placements(y, score_baseline, orient = orient)
  if (pm$n_pos < 2L || pm$n_neg < 2L) {
    return(list(
      variance = NA_real_, auc_method = pm$auc, auc_baseline = pb$auc,
      covariance = NA_real_, n_pos = pm$n_pos, n_neg = pm$n_neg
    ))
  }
  s10 <- stats::cov(cbind(pm$v10, pb$v10))
  s01 <- stats::cov(cbind(pm$v01, pb$v01))
  var_method <- s10[1L, 1L] / pm$n_pos + s01[1L, 1L] / pm$n_neg
  var_baseline <- s10[2L, 2L] / pm$n_pos + s01[2L, 2L] / pm$n_neg
  covariance <- s10[1L, 2L] / pm$n_pos + s01[1L, 2L] / pm$n_neg
  list(
    variance = max(0, var_method + var_baseline - 2 * covariance),
    auc_method = pm$auc,
    auc_baseline = pb$auc,
    covariance = covariance,
    n_pos = pm$n_pos,
    n_neg = pm$n_neg
  )
}
