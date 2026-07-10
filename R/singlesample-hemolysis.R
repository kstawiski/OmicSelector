#' @title Single-sample robust-regression hemolysis correction (Module B)
#'
#' @description
#' Frozen robust-regression hemolysis correction introduced in the
#' single-sample scoring bank (Module B). Performs per-feature robust
#' regression of log-expression on a hemolysis proxy score (and optionally a
#' platelet score), fitting only on training controls. At deployment, the
#' predicted nuisance contribution is subtracted from each feature.
#'
#' Note: the existing package file \code{hemolysis-correction.R} provides
#' \code{fit_hemolysis_prefilter()} (a different approach based on
#' marker-ratio gating). This file provides the Module-B fit-and-freeze
#' approach and is named \code{singlesample-hemolysis.R} to avoid collision.
#'
#' Methods provided:
#' \itemize{
#'   \item \code{\link{fit_hemolysis_rr}}: fit the per-feature robust-regression
#'     model.
#'   \item \code{\link{apply_hemolysis_rr}}: apply the frozen model to new
#'     samples.
#' }
#'
#' @references
#' Huber PJ. (1964) Robust Estimation of a Location Parameter.
#' \emph{The Annals of Mathematical Statistics} 35(1): 73–101.
#'
#' @name singlesample-hemolysis
NULL


# ----------------------------------------------------------------------------
# hemolysis_index_blondal -- canonical Blondal index helper
# ----------------------------------------------------------------------------

#' @title Blondal hemolysis index (log miR-451a - log miR-23a-3p)
#'
#' @description
#' Computes the eLife-14 Blondal-style hemolysis index on a log-abundance
#' matrix as \eqn{\log(\mathrm{miR\text{-}451a}) - \log(\mathrm{miR\text{-}23a\text{-}3p})}.
#' Both anchor features must be present (canonical names
#' \code{hsa-miR-451a} / \code{hsa-miR-23a-3p}, or the corresponding MIMAT
#' accessions \code{MIMAT0001631} / \code{MIMAT0000078}); otherwise
#' \code{NULL} is returned and the caller should treat the cell as not
#' hemolysis-correctable. This is the convention the single-sample pipeline
#' uses to feed \code{\link{fit_hemolysis_rr}} / \code{\link{apply_hemolysis_rr}}.
#'
#' @param X_log Numeric matrix (samples x features) on log-abundance scale.
#'   Column names must be canonical mature-miRNA names (post-resolver) or
#'   MIMAT accessions.
#'
#' @return Numeric vector of length \code{nrow(X_log)}, or \code{NULL} if
#'   either anchor feature is missing from \code{colnames(X_log)}.
#'
#' @references
#' Blondal T, Jensby Nielsen S, Baker A, Andreasen D, Mouritzen P,
#' Wrang Teilum M, Dahlsveen IK. (2013) Assessing sample and miRNA profile
#' quality in serum and plasma or other biofluids. \emph{Methods} 59(1): S1-S6.
#'
#' @examples
#' X_log <- matrix(c(8, 6, 4, 12, 5, 7, 3, 11), nrow = 2,
#'                 dimnames = list(NULL,
#'                                 c("hsa-miR-451a", "hsa-miR-23a-3p",
#'                                   "hsa-miR-1", "hsa-miR-2")))
#' hemolysis_index_blondal(X_log)
#'
#' @export
hemolysis_index_blondal <- function(X_log) {
  if (!is.matrix(X_log)) X_log <- as.matrix(X_log)
  if (is.null(colnames(X_log))) return(NULL)
  hemo_col <- which(colnames(X_log) %in% c("hsa-miR-451a", "MIMAT0001631"))
  ref_col  <- which(colnames(X_log) %in% c("hsa-miR-23a-3p", "MIMAT0000078"))
  if (length(hemo_col) == 0L || length(ref_col) == 0L) return(NULL)
  as.numeric(X_log[, hemo_col[1L]]) - as.numeric(X_log[, ref_col[1L]])
}


# ----------------------------------------------------------------------------
# fit_hemolysis_rr
# ----------------------------------------------------------------------------

#' @title Fit robust-regression hemolysis nuisance model from training controls
#'
#' @description
#' Fits per-feature Huber M-estimation robust regressions (\code{MASS::rlm})
#' of log-expression on a hemolysis proxy score (and optionally a platelet
#' proxy). The resulting intercept and slope coefficients are frozen at fit
#' time. At deployment, the predicted nuisance contribution is subtracted from
#' test-sample log-expressions via \code{\link{apply_hemolysis_rr}}.
#'
#' Failed fits (convergence failures, NA coefficients) fall back to a zero-slope
#' model (intercept = feature mean) and are flagged in the \code{stable} vector.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features). Ideally a
#'   mix of clean and mildly contaminated training controls.
#' @param hemo_score Numeric vector of length \code{nrow(X_train)}. Hemolysis
#'   proxy score per training sample.
#' @param platelet_score Optional numeric vector of length \code{nrow(X_train)}.
#'   Platelet contamination proxy. Default \code{NULL} (no platelet correction).
#' @param method \code{"huber"} (default, \code{MASS::rlm} with
#'   \code{psi.huber}) or \code{"rlm"} (alias).
#' @param clip_max Maximum absolute correction per feature in log-units.
#'   Default 2. Prevents over-correction from unstable slope estimates.
#'
#' @return Object of class \code{hemolysis_rr_fit} with components:
#'   \code{betas} (features \eqn{\times} coefficients matrix), \code{stable},
#'   \code{has_platelet}, \code{clip_max}, \code{train_hemo_range},
#'   \code{train_platelet_range}, \code{method}, \code{n_train}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rlnorm(100 * 20, 5, 0.8), nrow = 100,
#'             dimnames = list(NULL, paste0("mir-", seq_len(20))))
#' hemo <- runif(100)
#' fit <- fit_hemolysis_rr(X, hemo)
#' fit$n_train
#' }
#'
#' @references
#' Huber PJ. (1964) Robust Estimation of a Location Parameter.
#' \emph{The Annals of Mathematical Statistics} 35(1): 73–101.
#'
#' @export
fit_hemolysis_rr <- function(X_train,
                              hemo_score,
                              platelet_score = NULL,
                              method = c("huber", "rlm"),
                              clip_max = 2.0) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("fit_hemolysis_rr requires the MASS package (MASS::rlm).")
  }
  method <- match.arg(method)
  stopifnot(is.matrix(X_train) || is.data.frame(X_train))
  X <- as.matrix(X_train)
  if (length(hemo_score) != nrow(X)) {
    stop("hemo_score length must match nrow(X_train).")
  }

  log_X <- log(X + 1e-6 * rowSums(X))
  has_pl <- !is.null(platelet_score)
  if (has_pl && length(platelet_score) != nrow(X)) {
    stop("platelet_score length must match nrow(X_train).")
  }

  n_coef <- if (has_pl) 3L else 2L
  betas <- matrix(NA_real_, nrow = ncol(X), ncol = n_coef)
  colnames(betas) <- if (has_pl) c("intercept", "beta_hemo", "beta_plt") else
                                  c("intercept", "beta_hemo")
  rownames(betas) <- colnames(X)

  stable <- logical(ncol(X))

  for (j in seq_len(ncol(X))) {
    df <- data.frame(y = log_X[, j], h = hemo_score)
    if (has_pl) df$p <- platelet_score
    fit_j <- tryCatch(
      MASS::rlm(if (has_pl) y ~ h + p else y ~ h, data = df, method = "M",
                psi = MASS::psi.huber),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    if (is.null(fit_j)) {
      betas[j, ] <- c(mean(df$y, na.rm = TRUE), rep(0, n_coef - 1L))
      stable[j] <- FALSE
    } else {
      coefs <- coef(fit_j)
      betas[j, ] <- coefs
      stable[j] <- !any(is.na(coefs))
    }
  }

  fit <- list(
    betas = betas,
    stable = stable,
    has_platelet = has_pl,
    clip_max = clip_max,
    train_hemo_range = range(hemo_score, na.rm = TRUE),
    train_platelet_range = if (has_pl) range(platelet_score, na.rm = TRUE) else NULL,
    method = method,
    n_train = nrow(X)
  )
  class(fit) <- "hemolysis_rr_fit"
  fit
}


# ----------------------------------------------------------------------------
# apply_hemolysis_rr
# ----------------------------------------------------------------------------

#' @title Apply frozen hemolysis-RR correction to new samples
#'
#' @description
#' Subtracts the predicted hemolysis (and optional platelet) nuisance
#' contribution from one or more test samples using the frozen per-feature
#' regression coefficients from \code{\link{fit_hemolysis_rr}}.
#'
#' Corrections are clipped to \eqn{\pm}\code{clip_max} log-units to prevent
#' over-correction. If \code{report_oob = TRUE} (default), an \code{oob}
#' attribute flags samples whose hemolysis / platelet scores fall outside the
#' training range.
#'
#' @param x Numeric vector (single sample) or matrix (samples \eqn{\times}
#'   features).
#' @param fit Object of class \code{hemolysis_rr_fit} from
#'   \code{\link{fit_hemolysis_rr}}.
#' @param hemo_score Numeric vector (length = number of samples).
#' @param platelet_score Optional numeric vector (same length as
#'   \code{hemo_score}).
#' @param report_oob Logical. If \code{TRUE}, attaches an \code{oob} attribute
#'   flagging out-of-training-range samples. Default \code{TRUE}.
#'
#' @return Same shape as \code{x}. Hemolysis-corrected log-expression values,
#'   with attributes \code{oob} (if \code{report_oob}) and \code{fit_summary}.
#'
#' @export
apply_hemolysis_rr <- function(x, fit, hemo_score, platelet_score = NULL,
                                report_oob = TRUE) {
  stopifnot(inherits(fit, "hemolysis_rr_fit"))

  if (is.vector(x)) {
    return(apply_hemolysis_rr(
      matrix(x, nrow = 1L, dimnames = list(NULL, names(x))),
      fit, hemo_score, platelet_score, report_oob)[1L, ])
  }

  X <- as.matrix(x)
  log_X <- log(X + 1e-6 * rowSums(X))
  betas <- fit$betas
  out <- log_X

  for (i in seq_len(nrow(X))) {
    h <- hemo_score[i]
    pl <- if (fit$has_platelet) platelet_score[i] else 0
    pred <- if (fit$has_platelet)
      betas[, "beta_hemo"] * h + betas[, "beta_plt"] * pl
    else
      betas[, "beta_hemo"] * h
    pred_clipped <- pmin(pmax(pred, -fit$clip_max), fit$clip_max)
    out[i, ] <- log_X[i, ] - pred_clipped
  }

  if (report_oob) {
    oob_h <- hemo_score < fit$train_hemo_range[1] |
              hemo_score > fit$train_hemo_range[2]
    oob_p <- if (fit$has_platelet) {
      platelet_score < fit$train_platelet_range[1] |
        platelet_score > fit$train_platelet_range[2]
    } else {
      rep(FALSE, nrow(X))
    }
    attr(out, "oob") <- oob_h | oob_p
  }
  attr(out, "fit_summary") <- sprintf(
    "hemolysis_rr (n_train=%d, has_platelet=%s, n_stable=%d/%d)",
    fit$n_train, fit$has_platelet, sum(fit$stable), nrow(fit$betas))
  out
}
