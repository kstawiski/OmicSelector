#' @title Fit REO-UCell within-sample rank discriminator
#'
#' @description
#' Learns a frozen up/down feature signature from training data using
#' within-sample ranks, then stores only that signature and fixed
#' hyperparameters for single-sample UCell-style scoring. Training ranks are
#' computed per sample with rank 1 assigned to the most abundant feature. A
#' feature elevated in cases therefore has a lower mean case rank than control
#' rank, and is selected by
#' \code{mean_rank_controls - mean_rank_cases}. No test data or test-time batch
#' statistics are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{k} signature size per
#'   direction (default \code{min(25L, floor(ncol(X_train) / 4))}, at least 1),
#'   \code{maxRank} UCell rank cap (default
#'   \code{min(1500L, ncol(X_train))}), and \code{use_down} logical
#'   (default \code{TRUE}). These are resolved once during fitting and are not
#'   tuned inside \code{fit_reo_ucell()}.
#'
#' @return A plain list of class \code{reo_ucell_model} containing
#'   \code{up_features}, \code{down_features}, \code{maxRank},
#'   \code{feature_universe}, \code{use_down}, \code{k}, and the resolved
#'   \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 20)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
#' model <- fit_reo_ucell(X, y)
#' score <- score_reo_ucell(model, X)
#' }
#'
#' @references
#' Andreatta M, Carmona SJ. (2021) UCell: Robust and scalable single-cell gene
#' signature scoring. \emph{Computational and Structural Biotechnology Journal}
#' 19: 3796-3798.
#'
#' @export
fit_reo_ucell <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_reo_ucell", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_reo_ucell", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train))
  hp <- .reo_resolve_hp(hp, ncol(X_train))

  # rbind (not t(vapply)) preserves an n x p shape even when ncol == 1, where
  # vapply would collapse to a length-n vector and break the dimnames assignment.
  train_ranks <- do.call(rbind, lapply(seq_len(nrow(X_train)), function(i) {
    rank(-X_train[i, ], ties.method = "average")
  }))
  dimnames(train_ranks) <- dimnames(X_train)

  mean_rank_cases <- colMeans(train_ranks[y == 1, , drop = FALSE])
  mean_rank_controls <- colMeans(train_ranks[y == 0, , drop = FALSE])
  elevation_in_cases <- mean_rank_controls - mean_rank_cases

  # method = "radix" pins the feature-name tie-break to locale-independent byte
  # order, so the frozen signature is reproducible across machines/locales.
  up_order <- order(-elevation_in_cases, names(elevation_in_cases),
                    method = "radix")
  up_features <- names(elevation_in_cases)[up_order][seq_len(hp$k)]

  down_order <- order(elevation_in_cases, names(elevation_in_cases),
                      method = "radix")
  down_candidates <- setdiff(names(elevation_in_cases)[down_order], up_features)
  down_features <- utils::head(down_candidates, hp$k)

  model <- list(
    up_features = unname(up_features),
    down_features = unname(down_features),
    maxRank = hp$maxRank,
    feature_universe = colnames(X_train),
    use_down = hp$use_down,
    k = hp$k,
    hp = hp
  )
  class(model) <- "reo_ucell_model"
  model
}


#' @title Score REO-UCell within-sample rank signature
#'
#' @description
#' Scores each row of \code{X} independently with the frozen signature learned
#' by \code{\link{fit_reo_ucell}}. For one sample, features in the frozen
#' training universe that are present in \code{X} are ranked in decreasing
#' abundance using \code{rank(-x, ties.method = "average")}. For a signature
#' \eqn{S} of size \eqn{n}, ranks are capped at \code{maxRank_eff + 1}, where
#' \code{maxRank_eff = min(model$maxRank, number of present universe features)}.
#' The UCell-style score is
#' \deqn{1 - \{ \sum_{f \in S} \min(r_f, maxRank_eff + 1) - n(n+1)/2 \} /
#'       (n \times maxRank_eff).}
#' If no up or down signature features are present, that direction contributes
#' 0. With \code{use_down = TRUE}, the returned score is
#' \code{UCell(up) - UCell(down)}; otherwise it is \code{UCell(up)}.
#'
#' Because both the ranking universe and \code{maxRank_eff} depend on which
#' model-universe features are present in \code{X}, scores are only comparable
#' across specimens evaluated on the same feature panel; score each cohort in a
#' single comparison with one consistent set of features.
#'
#' @param model A \code{reo_ucell_model} object returned by
#'   \code{\link{fit_reo_ucell}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Numeric vector of one finite score per row of \code{X}. Scores use
#'   only each row's own values plus the frozen model, with no scored-batch
#'   centering, quantiles, or renormalization.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 20)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
#' model <- fit_reo_ucell(X, y)
#' score_reo_ucell(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Andreatta M, Carmona SJ. (2021) UCell: Robust and scalable single-cell gene
#' signature scoring. \emph{Computational and Structural Biotechnology Journal}
#' 19: 3796-3798.
#'
#' @export
score_reo_ucell <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_ucell_model")) {
    stop("score_reo_ucell: model must have class reo_ucell_model")
  }
  X <- .reo_check_matrix(X, "score_reo_ucell", "X")
  .reo_check_meta(meta, nrow(X), "score_reo_ucell", "meta")

  universe_present <- intersect(model$feature_universe, colnames(X))
  if (length(universe_present) == 0L) {
    stop("score_reo_ucell: no shared features between model feature_universe and X")
  }

  X_use <- X[, universe_present, drop = FALSE]
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    # Restore feature names: single-row slicing of a 1-column matrix drops them,
    # which would silently empty the signature intersection and return 0.
    x_i <- X_use[i, ]
    names(x_i) <- universe_present
    .reo_ucell_score_row(
      x = x_i,
      up = model$up_features,
      down = model$down_features,
      maxRank = model$maxRank,
      use_down = model$use_down
    )
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_reo_ucell: scorer produced non-finite or wrong-length output")
  }
  out
}


.reo_check_matrix <- function(X, fn, arg) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(fn, ": ", arg, " must be a numeric matrix of samples x features")
  }
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop(fn, ": ", arg, " must have at least one row and one feature")
  }
  feat_names <- colnames(X)
  if (is.null(feat_names) || any(!nzchar(feat_names))) {
    stop(fn, ": ", arg, " must have feature names (colnames)")
  }
  if (any(duplicated(feat_names))) {
    stop(fn, ": ", arg, " must have unique feature names (colnames)")
  }
  if (any(!is.finite(X))) {
    stop(fn, ": ", arg, " must contain only finite numeric values")
  }
  if (any(X < 0)) {
    stop(fn, ": negative values not allowed (compositional input)")
  }
  storage.mode(X) <- "double"
  X
}

.reo_check_labels <- function(y_train, n) {
  if (!is.numeric(y_train)) {
    stop("fit_reo_ucell: y_train must be numeric or integer 0/1 labels")
  }
  if (length(y_train) != n) {
    stop("fit_reo_ucell: length(y_train) must match nrow(X_train)")
  }
  if (any(!is.finite(y_train)) || any(!y_train %in% c(0, 1))) {
    stop("fit_reo_ucell: y_train must contain only 0/1 labels")
  }
  y <- as.integer(y_train)
  if (!any(y == 1L)) stop("fit_reo_ucell: y_train must contain at least one case")
  if (!any(y == 0L)) stop("fit_reo_ucell: y_train must contain at least one control")
  y
}

.reo_check_meta <- function(meta, n, fn, arg) {
  if (is.null(meta)) return(invisible(TRUE))
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (nrow(meta) != n) {
    stop(fn, ": ", arg, " must have one row per row of the expression matrix")
  }
  invisible(TRUE)
}

.reo_resolve_hp <- function(hp, p) {
  if (!is.list(hp)) stop("fit_reo_ucell: hp must be a list")
  allowed <- c("k", "maxRank", "use_down")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_reo_ucell: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  k_default <- max(1L, min(25L, floor(p / 4)))
  k <- hp$k
  if (is.null(k)) k <- k_default
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 1L || k != as.integer(k)) {
    stop("fit_reo_ucell: hp$k must be a positive integer")
  }
  k <- as.integer(k)
  if (k > p) stop("fit_reo_ucell: hp$k cannot exceed ncol(X_train)")

  maxRank <- hp$maxRank
  if (is.null(maxRank)) maxRank <- min(1500L, p)
  if (!is.numeric(maxRank) || length(maxRank) != 1L ||
      !is.finite(maxRank) || maxRank < 1L ||
      maxRank != as.integer(maxRank)) {
    stop("fit_reo_ucell: hp$maxRank must be a positive integer")
  }
  maxRank <- as.integer(maxRank)
  if (maxRank < k) {
    stop("fit_reo_ucell: hp$maxRank must be >= hp$k for valid UCell normalization")
  }

  use_down <- hp$use_down
  if (is.null(use_down)) use_down <- TRUE
  if (!is.logical(use_down) || length(use_down) != 1L || is.na(use_down)) {
    stop("fit_reo_ucell: hp$use_down must be TRUE or FALSE")
  }

  list(k = k, maxRank = maxRank, use_down = use_down)
}

.reo_ucell_score_row <- function(x, up, down, maxRank, use_down) {
  ranks <- rank(-x, ties.method = "average")
  names(ranks) <- names(x)
  maxRank_eff <- min(maxRank, length(ranks))
  up_score <- .reo_ucell_direction(ranks, up, maxRank_eff)
  if (!use_down) return(up_score)
  up_score - .reo_ucell_direction(ranks, down, maxRank_eff)
}

.reo_ucell_direction <- function(ranks, features, maxRank_eff) {
  present <- intersect(features, names(ranks))
  if (length(present) == 0L) return(0)
  r_cap <- pmin(ranks[present], maxRank_eff + 1)
  n <- length(r_cap)
  u <- sum(r_cap) - n * (n + 1) / 2
  1 - u / (n * maxRank_eff)
}
