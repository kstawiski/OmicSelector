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


# allow_na: when FALSE (default; every existing caller) NA/non-finite values are
# rejected -- behaviour is byte-identical to before. When TRUE (the qPCR-Ct entry
# points, where NA encodes a censored/undetected measurement) NA is permitted to pass
# through, and the negative-value check uses na.rm so NA does not break it; NA is mapped
# to a positive floor downstream before any numeric use.
.reo_check_matrix <- function(X, fn, arg, allow_na = FALSE) {
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
  if (!allow_na && any(!is.finite(X))) {
    stop(fn, ": ", arg, " must contain only finite numeric values")
  }
  if (any(X < 0, na.rm = TRUE)) {
    stop(fn, ": negative values not allowed (compositional input)")
  }
  storage.mode(X) <- "double"
  X
}

.reo_check_labels <- function(y_train, n, fn = "fit_reo_ucell") {
  if (!is.numeric(y_train)) {
    stop(fn, ": y_train must be numeric or integer 0/1 labels")
  }
  if (length(y_train) != n) {
    stop(fn, ": length(y_train) must match nrow(X_train)")
  }
  if (any(!is.finite(y_train)) || any(!y_train %in% c(0, 1))) {
    stop(fn, ": y_train must contain only 0/1 labels")
  }
  y <- as.integer(y_train)
  if (!any(y == 1L)) stop(fn, ": y_train must contain at least one case")
  if (!any(y == 0L)) stop(fn, ": y_train must contain at least one control")
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


#' @title Fit REO-singscore within-sample rank discriminator
#'
#' @description
#' Learns a frozen up/down feature signature from training data using
#' within-sample fractional ranks, then stores only that signature and fixed
#' hyperparameters for single-sample singscore-style scoring. Training ranks are
#' computed per sample with ascending ranks divided by the training feature
#' count; higher fractional rank means higher abundance in that sample. A
#' feature elevated in cases therefore has a higher mean case fractional rank
#' than control fractional rank, and is selected by
#' \code{mean_fracrank_cases - mean_fracrank_controls}. No test data or
#' test-time batch statistics are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{k} signature size per
#'   direction (default \code{min(25L, floor(ncol(X_train) / 4))}, at least 1)
#'   and \code{use_down} logical (default \code{TRUE}). These are resolved once
#'   during fitting and are not tuned inside \code{fit_reo_singscore()}.
#'
#' @return A plain list of class \code{reo_singscore_model} containing
#'   \code{up_features}, \code{down_features}, \code{feature_universe},
#'   \code{use_down}, \code{k}, and the resolved \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 20)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
#' model <- fit_reo_singscore(X, y)
#' score <- score_reo_singscore(model, X)
#' }
#'
#' @references
#' Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ. (2018)
#' Single sample scoring of molecular phenotypes. \emph{BMC Bioinformatics}
#' 19: 404. PMID 30400809.
#'
#' @export
fit_reo_singscore <- function(X_train, y_train, meta_train = NULL,
                              hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_reo_singscore", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_reo_singscore",
                  "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_reo_singscore")
  hp <- .reo_resolve_hp_singscore(hp, ncol(X_train))

  # rbind (not t(vapply)) preserves an n x p shape even when ncol == 1, where
  # vapply would collapse to a length-n vector and break the dimnames assignment.
  train_fracranks <- do.call(rbind, lapply(seq_len(nrow(X_train)), function(i) {
    rank(X_train[i, ], ties.method = "average") / ncol(X_train)
  }))
  dimnames(train_fracranks) <- dimnames(X_train)

  mean_fracrank_cases <- colMeans(train_fracranks[y == 1, , drop = FALSE])
  mean_fracrank_controls <- colMeans(train_fracranks[y == 0, , drop = FALSE])
  elevation_in_cases <- mean_fracrank_cases - mean_fracrank_controls

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
    feature_universe = colnames(X_train),
    use_down = hp$use_down,
    k = hp$k,
    hp = hp
  )
  class(model) <- "reo_singscore_model"
  model
}


#' @title Score REO-singscore within-sample rank signature
#'
#' @description
#' Scores each row of \code{X} independently with the frozen signature learned
#' by \code{\link{fit_reo_singscore}}. For one sample, features in the frozen
#' training universe that are present in \code{X} are ranked in increasing
#' abundance using
#' \code{rank(x, ties.method = "average") / length(x)}. The up term is the mean
#' fractional rank of present up features; the down term is the mean fractional
#' rank of present down features. If a signature direction has no present
#' features, it contributes the neutral expected fractional rank \code{0.5}.
#' With \code{use_down = TRUE}, the returned score is
#' \code{mean_rank(up) - mean_rank(down)}; otherwise it is \code{mean_rank(up)}.
#'
#' @param model A \code{reo_singscore_model} object returned by
#'   \code{\link{fit_reo_singscore}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Numeric vector of one finite score per row of \code{X}. Scores use
#'   only each row's own values plus the frozen model, with no scored-batch
#'   centering, quantiles, or renormalization. If no model-universe features are
#'   shared with \code{X}, scoring stops explicitly.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 20)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
#' model <- fit_reo_singscore(X, y)
#' score_reo_singscore(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ. (2018)
#' Single sample scoring of molecular phenotypes. \emph{BMC Bioinformatics}
#' 19: 404. PMID 30400809.
#'
#' @export
score_reo_singscore <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_singscore_model")) {
    stop("score_reo_singscore: model must have class reo_singscore_model")
  }
  X <- .reo_check_matrix(X, "score_reo_singscore", "X")
  .reo_check_meta(meta, nrow(X), "score_reo_singscore", "meta")

  universe_present <- intersect(model$feature_universe, colnames(X))
  if (length(universe_present) == 0L) {
    stop("score_reo_singscore: no shared features between model feature_universe and X")
  }

  X_use <- X[, universe_present, drop = FALSE]
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    # Restore feature names: single-row slicing of a 1-column matrix drops them,
    # which would silently empty the signature intersection and return neutral.
    x_i <- X_use[i, ]
    names(x_i) <- universe_present
    .reo_singscore_score_row(
      x = x_i,
      up = model$up_features,
      down = model$down_features,
      use_down = model$use_down
    )
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_reo_singscore: scorer produced non-finite or wrong-length output")
  }
  out
}

.reo_resolve_hp_singscore <- function(hp, p) {
  if (!is.list(hp)) stop("fit_reo_singscore: hp must be a list")
  allowed <- c("k", "use_down")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_reo_singscore: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  k_default <- max(1L, min(25L, floor(p / 4)))
  k <- hp$k
  if (is.null(k)) k <- k_default
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 1L || k != as.integer(k)) {
    stop("fit_reo_singscore: hp$k must be a positive integer")
  }
  k <- as.integer(k)
  if (k > p) stop("fit_reo_singscore: hp$k cannot exceed ncol(X_train)")

  use_down <- hp$use_down
  if (is.null(use_down)) use_down <- TRUE
  if (!is.logical(use_down) || length(use_down) != 1L || is.na(use_down)) {
    stop("fit_reo_singscore: hp$use_down must be TRUE or FALSE")
  }

  list(k = k, use_down = use_down)
}

.reo_singscore_score_row <- function(x, up, down, use_down) {
  fracranks <- rank(x, ties.method = "average") / length(x)
  names(fracranks) <- names(x)
  up_score <- .reo_singscore_direction(fracranks, up)
  if (!use_down) return(up_score)
  up_score - .reo_singscore_direction(fracranks, down)
}

.reo_singscore_direction <- function(fracranks, features) {
  present <- intersect(features, names(fracranks))
  if (length(present) == 0L) return(0.5)
  mean(fracranks[present])
}


#' @title Fit REO-kTSP within-sample pair-order discriminator
#'
#' @description
#' Learns a frozen k-Top-Scoring-Pairs (k-TSP) discriminator from training
#' samples by calling \code{\link{os_ktsp_fit}}. The fitted primitive stores
#' oriented feature pairs so that larger vote fractions indicate the case class:
#' each retained pair votes 1 for a specimen when
#' \code{feature_a < feature_b} within that specimen. No test data or
#' test-time batch statistics are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{k}, the number of oriented
#'   pairs requested from \code{\link{os_ktsp_fit}} (default \code{11L}).
#'   \code{k} is resolved once during fitting and is not tuned inside
#'   \code{fit_reo_ktsp()}.
#'
#' @return A plain list of class \code{reo_ktsp_model} containing the fitted
#'   \code{os_ktsp_model} in \code{ktsp}, \code{feature_universe}, retained
#'   pair count \code{k}, and the resolved \code{hp}.
#'
#' @examples
#' \dontrun{
#' X <- matrix(stats::rgamma(40 * 12, shape = 2), nrow = 40,
#'             dimnames = list(NULL, paste0("miR-", seq_len(12))))
#' y <- rep(c(0, 1), each = 20)
#' X[y == 0, "miR-1"] <- X[y == 0, "miR-1"] + 5
#' X[y == 1, "miR-2"] <- X[y == 1, "miR-2"] + 5
#' model <- fit_reo_ktsp(X, y, hp = list(k = 3L))
#' score <- score_reo_ktsp(model, X)
#' }
#'
#' @references
#' Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
#' expression profiles from pairwise mRNA comparisons. \emph{Statistical
#' Applications in Genetics and Molecular Biology} 3: Article19.
#'
#' Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. (2005) Simple decision rules
#' for classifying human cancers from gene expression profiles. \emph{Bioinformatics}
#' 21: 3896-3904.
#'
#' @export
fit_reo_ktsp <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_reo_ktsp", "X_train")
  if (ncol(X_train) < 2L) {
    stop("fit_reo_ktsp: X_train must contain at least two features for k-TSP")
  }
  .reo_check_meta(meta_train, nrow(X_train), "fit_reo_ktsp", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_reo_ktsp")
  hp <- .reo_resolve_hp_ktsp(hp, ncol(X_train))

  ktsp <- os_ktsp_fit(X_train, y, k = hp$k)
  model <- list(
    ktsp = ktsp,
    feature_universe = colnames(X_train),
    k = ktsp$k,
    hp = hp
  )
  class(model) <- "reo_ktsp_model"
  model
}


#' @title Score REO-kTSP within-sample pair-order votes
#'
#' @description
#' Scores each row of \code{X} independently with the frozen oriented pairs
#' learned by \code{\link{fit_reo_ktsp}}. For the pairs whose two features are
#' present in \code{X}, the score is the k-TSP vote fraction: the mean of
#' \code{feature_a < feature_b} votes over usable pairs. This is the standard
#' k-TSP decision statistic and supplies a continuous margin score for DeLong
#' comparisons. If no retained pair has both features present, the function
#' returns the neutral vote fraction \code{0.5} for every row.
#'
#' @param model A \code{reo_ktsp_model} object returned by
#'   \code{\link{fit_reo_ktsp}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Numeric vector of one finite vote-fraction score per row of
#'   \code{X}. Scores use only each row's own values plus the frozen model, with
#'   no scored-batch centering, quantiles, or renormalization. On a fully
#'   present feature set the result is identical to
#'   \code{predict(model$ktsp, X)}.
#'
#' @examples
#' \dontrun{
#' X <- matrix(stats::rgamma(40 * 12, shape = 2), nrow = 40,
#'             dimnames = list(NULL, paste0("miR-", seq_len(12))))
#' y <- rep(c(0, 1), each = 20)
#' X[y == 0, "miR-1"] <- X[y == 0, "miR-1"] + 5
#' X[y == 1, "miR-2"] <- X[y == 1, "miR-2"] + 5
#' model <- fit_reo_ktsp(X, y, hp = list(k = 3L))
#' score_reo_ktsp(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
#' expression profiles from pairwise mRNA comparisons. \emph{Statistical
#' Applications in Genetics and Molecular Biology} 3: Article19.
#'
#' Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. (2005) Simple decision rules
#' for classifying human cancers from gene expression profiles. \emph{Bioinformatics}
#' 21: 3896-3904.
#'
#' @export
score_reo_ktsp <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_ktsp_model")) {
    stop("score_reo_ktsp: model must have class reo_ktsp_model")
  }
  X <- .reo_check_matrix(X, "score_reo_ktsp", "X")
  .reo_check_meta(meta, nrow(X), "score_reo_ktsp", "meta")

  usable_pairs <- .reo_ktsp_usable_pairs(model$ktsp$pairs, colnames(X))
  if (nrow(usable_pairs) == 0L) {
    return(rep(0.5, nrow(X)))
  }

  out <- .reo_ktsp_vote_fraction(X, usable_pairs)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_reo_ktsp: scorer produced non-finite or wrong-length output")
  }
  out
}

.reo_resolve_hp_ktsp <- function(hp, p) {
  if (!is.list(hp)) stop("fit_reo_ktsp: hp must be a list")
  allowed <- "k"
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_reo_ktsp: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  k <- hp$k
  if (is.null(k)) k <- 11L
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 1L || k != as.integer(k)) {
    stop("fit_reo_ktsp: hp$k must be a positive integer")
  }
  list(k = as.integer(k))
}

.reo_ktsp_usable_pairs <- function(pairs, feature_names) {
  keep <- pairs$feature_a %in% feature_names & pairs$feature_b %in% feature_names
  pairs[keep, , drop = FALSE]
}

.reo_ktsp_vote_fraction <- function(X, pairs) {
  votes <- do.call(cbind, lapply(seq_len(nrow(pairs)), function(i) {
    a <- pairs$feature_a[i]
    b <- pairs$feature_b[i]
    as.numeric(X[, a] < X[, b])
  }))
  as.numeric(rowMeans(votes))
}


#' @title Fit REO-pairratio within-sample log-ratio discriminator
#'
#' @description
#' Learns a frozen penalized-logistic discriminator over selected pairwise
#' log-ratios. Features are pre-screened using only the training data by the
#' absolute two-sample t-statistic on \code{log(X_train + pseudocount)}. All
#' pairwise log-ratios among the screened features are then generated with
#' \code{\link{ws_logratio}}, and \code{glmnet::cv.glmnet()} fits a binomial
#' elastic-net model with deterministic folds. No test data or test-time batch
#' statistics are used.
#'
#' A retained ratio \eqn{\log(x_a + c) - \log(x_b + c)} is computed from one
#' specimen's two features, so a frozen linear predictor over these ratios is
#' single-sample deployable. Larger returned linear predictors are more
#' case-like.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{pseudocount} (default
#'   \code{0.5}, positive), \code{m_features} (default
#'   \code{min(40L, ncol(X_train))}, at least 2), \code{alpha} (default
#'   \code{1}, in \code{[0, 1]}), \code{nfolds} (default \code{10L}, at least
#'   3), \code{lambda_rule} (default \code{"1se"}, or \code{"min"}), and
#'   \code{seed} (default \code{1L}, used only to make fold assignment
#'   reproducible while restoring the global RNG state).
#'
#' @return A plain list of class \code{reo_pairratio_model} containing
#'   \code{selected_pairs}, nonzero \code{coefficients}, \code{intercept},
#'   \code{pseudocount}, \code{screened_features}, \code{feature_universe},
#'   \code{lambda}, \code{lambda_rule}, and resolved \code{hp}. If the chosen
#'   lambda yields an intercept-only glmnet model, \code{selected_pairs} has
#'   zero rows and scoring returns the intercept for every specimen.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 20, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(20))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, "miR-1"] <- X[y == 1, "miR-1"] + 10
#' X[y == 0, "miR-2"] <- X[y == 0, "miR-2"] + 10
#' model <- fit_reo_pairratio(X, y, hp = list(m_features = 8L, nfolds = 5L))
#' score <- score_reo_pairratio(model, X)
#' }
#'
#' @references
#' Aitchison J. (1986) \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Friedman J, Hastie T, Tibshirani R. (2010) Regularization paths for
#' generalized linear models via coordinate descent. \emph{Journal of
#' Statistical Software} 33: 1-22.
#'
#' @export
fit_reo_pairratio <- function(X_train, y_train, meta_train = NULL,
                              hp = list()) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("fit_reo_pairratio: package 'glmnet' is required for reo-pairratio",
         call. = FALSE)
  }
  X_train <- .reo_check_matrix(X_train, "fit_reo_pairratio", "X_train")
  if (ncol(X_train) < 2L) {
    stop("fit_reo_pairratio: X_train must contain at least two features for log-ratios")
  }
  .reo_check_meta(meta_train, nrow(X_train), "fit_reo_pairratio",
                  "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_reo_pairratio")
  if (min(table(y)) < 2L) {
    stop("fit_reo_pairratio: each class must have at least 2 samples for ",
         "cross-validated glmnet (got ", sum(y == 1L), " case(s), ",
         sum(y == 0L), " control(s))", call. = FALSE)
  }
  hp <- .reo_resolve_hp_pairratio(hp, ncol(X_train), nrow(X_train))

  L <- log(X_train + hp$pseudocount)
  signal <- .reo_pairratio_univariate_signal(L, y)
  feature_order <- order(-abs(signal), names(signal), method = "radix")
  screened_features <- names(signal)[feature_order][seq_len(hp$m_features)]

  L_screened <- L[, screened_features, drop = FALSE]
  R <- ws_logratio(L_screened, feature_names = screened_features)
  if (ncol(R) < 2L) {
    stop("fit_reo_pairratio: need at least 3 screened features to form >= 2 ",
         "log-ratios for glmnet; increase hp$m_features or provide more features",
         call. = FALSE)
  }
  pair_map <- .reo_pairratio_pair_map(screened_features)
  pair_map$ratio_id <- paste0("ratio_", seq_len(nrow(pair_map)))
  colnames(R) <- pair_map$ratio_id

  # RNG hygiene: snapshot/restore the global .Random.seed around ALL RNG-using
  # fitting (seeded folds + cv.glmnet) so fit leaves no global RNG side-effect,
  # even when no seed existed beforehand (cv.glmnet can create one).
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)

  foldid <- .reo_pairratio_foldid(y, hp$nfolds, hp$seed)
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(
      x = R,
      y = y,
      family = "binomial",
      alpha = hp$alpha,
      foldid = foldid,
      nfolds = hp$nfolds
    ),
    error = function(e) {
      stop("fit_reo_pairratio: cv.glmnet failed (training partition likely too ",
           "small or imbalanced for nfolds=", hp$nfolds, "): ",
           conditionMessage(e), call. = FALSE)
    }
  )

  lambda_name <- if (identical(hp$lambda_rule, "1se")) "lambda.1se" else "lambda.min"
  chosen_lambda <- unname(cv_fit[[lambda_name]])
  coef_mat <- as.matrix(stats::coef(cv_fit, s = lambda_name))
  coef_vec <- as.numeric(coef_mat[, 1L])
  names(coef_vec) <- rownames(coef_mat)

  intercept <- unname(coef_vec["(Intercept)"])
  beta <- coef_vec[setdiff(names(coef_vec), "(Intercept)")]
  selected_idx <- which(beta != 0)
  selected_pairs <- pair_map[match(names(beta)[selected_idx], pair_map$ratio_id),
                             c("feature_a", "feature_b", "ratio_name",
                               "ratio_id"),
                             drop = FALSE]
  rownames(selected_pairs) <- NULL
  coefficients <- unname(beta[selected_idx])
  names(coefficients) <- selected_pairs$ratio_id

  model <- list(
    selected_pairs = selected_pairs,
    coefficients = coefficients,
    intercept = intercept,
    pseudocount = hp$pseudocount,
    screened_features = unname(screened_features),
    feature_universe = colnames(X_train),
    lambda = chosen_lambda,
    lambda_rule = hp$lambda_rule,
    hp = hp
  )
  class(model) <- "reo_pairratio_model"
  model
}


#' @title Score REO-pairratio within-sample log-ratio discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen pairwise
#' log-ratio coefficients learned by \code{\link{fit_reo_pairratio}}. For one
#' specimen and each stored pair \code{(a, b)}, if both features are present in
#' \code{X}, the term is
#' \code{log(X[i, a] + pseudocount) - log(X[i, b] + pseudocount)}. If either
#' feature is absent, that pair is dropped for that specimen. The returned
#' score is the frozen linear predictor, not \code{plogis()} of it.
#'
#' If the model is intercept-only, or if none of the stored pairs is usable for
#' a scored specimen, the score is the intercept. Scoring uses no random numbers
#' and no scored-batch statistics.
#'
#' @param model A \code{reo_pairratio_model} object returned by
#'   \code{\link{fit_reo_pairratio}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Numeric vector of one finite linear-predictor score per row of
#'   \code{X}. Larger scores are more case-like.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 20, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(20))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, "miR-1"] <- X[y == 1, "miR-1"] + 10
#' X[y == 0, "miR-2"] <- X[y == 0, "miR-2"] + 10
#' model <- fit_reo_pairratio(X, y, hp = list(m_features = 8L, nfolds = 5L))
#' score_reo_pairratio(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Aitchison J. (1986) \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Friedman J, Hastie T, Tibshirani R. (2010) Regularization paths for
#' generalized linear models via coordinate descent. \emph{Journal of
#' Statistical Software} 33: 1-22.
#'
#' @export
score_reo_pairratio <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_pairratio_model")) {
    stop("score_reo_pairratio: model must have class reo_pairratio_model")
  }
  X <- .reo_check_matrix(X, "score_reo_pairratio", "X")
  .reo_check_meta(meta, nrow(X), "score_reo_pairratio", "meta")

  out <- vapply(seq_len(nrow(X)), function(i) {
    # Restore feature names: single-row slicing of a 1-column matrix drops them,
    # which would otherwise make partial-overlap handling depend on batch shape.
    x_i <- X[i, ]
    names(x_i) <- colnames(X)
    .reo_pairratio_score_row(
      x = x_i,
      pairs = model$selected_pairs,
      coefficients = model$coefficients,
      intercept = model$intercept,
      pseudocount = model$pseudocount
    )
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_reo_pairratio: scorer produced non-finite or wrong-length output")
  }
  out
}

.reo_resolve_hp_pairratio <- function(hp, p, n) {
  if (!is.list(hp)) stop("fit_reo_pairratio: hp must be a list")
  allowed <- c("pseudocount", "m_features", "alpha", "nfolds",
               "lambda_rule", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_reo_pairratio: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  pseudocount <- hp$pseudocount
  if (is.null(pseudocount)) pseudocount <- 0.5
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      !is.finite(pseudocount) || pseudocount <= 0) {
    stop("fit_reo_pairratio: hp$pseudocount must be a positive finite number")
  }

  m_features <- hp$m_features
  if (is.null(m_features)) m_features <- min(40L, p)
  if (!is.numeric(m_features) || length(m_features) != 1L ||
      !is.finite(m_features) || m_features < 2L ||
      m_features != as.integer(m_features)) {
    stop("fit_reo_pairratio: hp$m_features must be an integer >= 2")
  }
  m_features <- min(as.integer(m_features), p)

  alpha <- hp$alpha
  if (is.null(alpha)) alpha <- 1
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("fit_reo_pairratio: hp$alpha must be a finite number in [0, 1]")
  }

  nfolds <- hp$nfolds
  if (is.null(nfolds)) nfolds <- 10L
  if (!is.numeric(nfolds) || length(nfolds) != 1L || !is.finite(nfolds) ||
      nfolds < 3L || nfolds != as.integer(nfolds)) {
    stop("fit_reo_pairratio: hp$nfolds must be an integer >= 3")
  }
  nfolds <- as.integer(nfolds)
  if (nfolds > n) {
    stop("fit_reo_pairratio: hp$nfolds cannot exceed nrow(X_train)")
  }

  lambda_rule <- hp$lambda_rule
  if (is.null(lambda_rule)) lambda_rule <- "1se"
  if (!is.character(lambda_rule) || length(lambda_rule) != 1L ||
      is.na(lambda_rule) || !lambda_rule %in% c("1se", "min")) {
    stop("fit_reo_pairratio: hp$lambda_rule must be '1se' or 'min'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_reo_pairratio: hp$seed must be a non-negative integer")
  }

  list(
    pseudocount = as.numeric(pseudocount),
    m_features = m_features,
    alpha = as.numeric(alpha),
    nfolds = nfolds,
    lambda_rule = lambda_rule,
    seed = as.integer(seed)
  )
}

.reo_pairratio_univariate_signal <- function(L, y) {
  case <- y == 1L
  control <- y == 0L
  mean_case <- colMeans(L[case, , drop = FALSE])
  mean_control <- colMeans(L[control, , drop = FALSE])
  var_case <- apply(L[case, , drop = FALSE], 2L, stats::var)
  var_control <- apply(L[control, , drop = FALSE], 2L, stats::var)
  var_case[is.na(var_case)] <- 0
  var_control[is.na(var_control)] <- 0
  delta <- mean_case - mean_control
  se <- sqrt(var_case / sum(case) + var_control / sum(control))
  signal <- ifelse(se > 0, delta / se, ifelse(delta == 0, 0, sign(delta) * Inf))
  names(signal) <- colnames(L)
  signal
}

.reo_pairratio_pair_map <- function(features) {
  out <- vector("list", length(features) * (length(features) - 1L) / 2L)
  k <- 1L
  for (i in seq_len(length(features) - 1L)) {
    for (j in (i + 1L):length(features)) {
      out[[k]] <- data.frame(
        feature_a = features[i],
        feature_b = features[j],
        ratio_name = paste0(features[i], "_vs_", features[j]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

.reo_pairratio_foldid <- function(y, nfolds, seed) {
  # Seeded, class-stratified fold assignment. The CALLER (fit_reo_pairratio)
  # owns the global-.Random.seed snapshot/restore that spans this call and the
  # subsequent cv.glmnet, so this helper only sets the seed and assigns folds.
  set.seed(seed)
  foldid <- integer(length(y))
  for (cls in c(0L, 1L)) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx), replace = FALSE)
    foldid[idx] <- rep(seq_len(nfolds), length.out = length(idx))
  }
  foldid
}

.reo_pairratio_score_row <- function(x, pairs, coefficients, intercept,
                                     pseudocount) {
  if (nrow(pairs) == 0L) return(intercept)
  usable <- pairs$feature_a %in% names(x) & pairs$feature_b %in% names(x)
  if (!any(usable)) return(intercept)
  idx <- which(usable)
  terms <- vapply(idx, function(k) {
    log(x[[pairs$feature_a[k]]] + pseudocount) -
      log(x[[pairs$feature_b[k]]] + pseudocount)
  }, numeric(1))
  intercept + sum(coefficients[idx] * terms)
}


#' @title Fit REO-MetaKTSP meta-analytic pair-order discriminator
#'
#' @description
#' Learns a meta-analytic k-Top-Scoring-Pairs discriminator across training
#' cohorts named by \code{meta_train[[cohort_col]]}. Each unordered feature
#' pair receives an avgTSP score: the simple equal-cohort-weight mean of
#' per-cohort case-control differences in \code{feature_i < feature_j}. This is
#' a transfer-estimand method at training because pair selection uses
#' cross-cohort information, but it remains single-sample at inference because
#' \code{\link{score_reo_metaktsp}} uses only one specimen's own pair order.
#' When no usable cohort metadata is available, or only one non-missing cohort
#' is present, fitting gracefully degenerates to ordinary single-cohort
#' \code{\link{fit_reo_ktsp}} behavior over the training rows. Rows with
#' missing cohort labels are dropped from meta-analytic selection, cohorts with
#' only cases or only controls are skipped, and a zero-usable-cohort meta split
#' falls back to the single-cohort path. At score time, no usable retained pairs
#' returns the neutral vote fraction \code{0.5}.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. When it contains
#'   \code{cohort_col}, non-missing values define the training cohorts used for
#'   avgTSP pair selection.
#' @param hp List of frozen hyperparameters: \code{k}, the number of oriented
#'   pairs requested (default \code{11L}), and \code{cohort_col}, a single
#'   non-empty character string naming the cohort/study column in
#'   \code{meta_train} (default \code{"accession"}). These are resolved once
#'   during fitting and are not tuned inside \code{fit_reo_metaktsp()}.
#'
#' @return A plain list of class \code{reo_metaktsp_model} containing retained
#'   oriented \code{pairs}, \code{feature_universe}, retained pair count
#'   \code{k}, \code{cohort_col}, number of cohorts used in selection
#'   \code{n_cohorts}, and the resolved \code{hp}. \code{n_cohorts} is the
#'   auditable meta-engagement signal: \code{n_cohorts > 1} means avgTSP
#'   selection ran across multiple cohorts, while \code{n_cohorts == 1} means
#'   fitting degenerated to single-cohort k-TSP (NULL/missing/single-level
#'   \code{cohort_col}, or a zero-usable-cohort fallback). A caller that
#'   intends cross-cohort meta-analysis should check \code{n_cohorts > 1} so a
#'   mis-wired \code{meta_train} cannot silently demote this transfer method to
#'   pooled within-cohort k-TSP.
#'
#' @examples
#' \dontrun{
#' X <- matrix(stats::rgamma(60 * 12, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(12))))
#' y <- rep(c(0, 1), times = 30)
#' meta <- data.frame(accession = rep(paste0("GSE", 1:3), each = 20))
#' X[y == 0, "miR-1"] <- X[y == 0, "miR-1"] + 5
#' X[y == 1, "miR-2"] <- X[y == 1, "miR-2"] + 5
#' model <- fit_reo_metaktsp(X, y, meta_train = meta, hp = list(k = 3L))
#' score <- score_reo_metaktsp(model, X)
#' }
#'
#' @references
#' Kim S, Lin C-W, Tseng GC. (2016) MetaKTSP: a meta-analytically derived
#' k-top-scoring-pair classifier for robust prediction of human cancer
#' subtypes. \emph{Bioinformatics} 32: 1966-1973. PMID: 27153719.
#'
#' Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
#' expression profiles from pairwise mRNA comparisons. \emph{Statistical
#' Applications in Genetics and Molecular Biology} 3: Article19.
#'
#' Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. (2005) Simple decision rules
#' for classifying human cancers from gene expression profiles. \emph{Bioinformatics}
#' 21: 3896-3904.
#'
#' @export
fit_reo_metaktsp <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_reo_metaktsp", "X_train")
  if (ncol(X_train) < 2L) {
    stop("fit_reo_metaktsp: X_train must contain at least two features for k-TSP")
  }
  .reo_check_meta(meta_train, nrow(X_train), "fit_reo_metaktsp", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_reo_metaktsp")
  hp <- .reo_resolve_hp_metaktsp(hp, ncol(X_train))

  if (is.null(meta_train)) {
    return(.reo_metaktsp_single_cohort_model(X_train, y, hp))
  }
  meta_train <- as.data.frame(meta_train, stringsAsFactors = FALSE)
  if (!hp$cohort_col %in% names(meta_train)) {
    return(.reo_metaktsp_single_cohort_model(X_train, y, hp))
  }

  cohort <- as.character(meta_train[[hp$cohort_col]])
  non_missing <- !is.na(cohort)
  cohort_levels <- unique(cohort[non_missing])
  if (length(cohort_levels) <= 1L) {
    return(.reo_metaktsp_single_cohort_model(X_train, y, hp))
  }

  usable_cohorts <- cohort_levels[vapply(cohort_levels, function(cohort_i) {
    rows <- non_missing & cohort == cohort_i
    any(y[rows] == 1L) && any(y[rows] == 0L)
  }, logical(1))]
  if (length(usable_cohorts) == 0L) {
    return(.reo_metaktsp_single_cohort_model(X_train, y, hp,
                                             rows = non_missing))
  }

  pairs <- .reo_metaktsp_pairs(X_train, y, cohort, usable_cohorts)
  pairs <- pairs[order(-pairs$score, pairs$feature_a, pairs$feature_b,
                       method = "radix"), , drop = FALSE]
  pairs <- utils::head(pairs, min(hp$k, nrow(pairs)))
  model <- list(
    pairs = pairs,
    feature_universe = colnames(X_train),
    k = nrow(pairs),
    cohort_col = hp$cohort_col,
    n_cohorts = length(usable_cohorts),
    hp = hp
  )
  class(model) <- "reo_metaktsp_model"
  model
}


#' @title Score REO-MetaKTSP within-sample pair-order votes
#'
#' @description
#' Scores each row of \code{X} independently with the frozen oriented pairs
#' learned by \code{\link{fit_reo_metaktsp}}. The model is trained
#' meta-analytically across cohorts via \code{meta_train[[cohort_col]]} using
#' avgTSP with equal cohort weight, but score-time inference is single-sample:
#' each retained pair votes 1 only when \code{feature_a < feature_b} within
#' that specimen. The graceful single-cohort degeneration, missing-cohort-row
#' drop, and single-class-cohort skip are fit-time behaviors only. If no
#' retained pair has both features present in \code{X}, scoring returns the
#' neutral vote fraction \code{0.5} for every row.
#'
#' @param model A \code{reo_metaktsp_model} object returned by
#'   \code{\link{fit_reo_metaktsp}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Numeric vector of one finite vote-fraction score per row of
#'   \code{X}. Scores use only each row's own values plus the frozen model, with
#'   no scored-batch centering, quantiles, or renormalization.
#'
#' @examples
#' \dontrun{
#' X <- matrix(stats::rgamma(60 * 12, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(12))))
#' y <- rep(c(0, 1), times = 30)
#' meta <- data.frame(accession = rep(paste0("GSE", 1:3), each = 20))
#' X[y == 0, "miR-1"] <- X[y == 0, "miR-1"] + 5
#' X[y == 1, "miR-2"] <- X[y == 1, "miR-2"] + 5
#' model <- fit_reo_metaktsp(X, y, meta_train = meta, hp = list(k = 3L))
#' score_reo_metaktsp(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Kim S, Lin C-W, Tseng GC. (2016) MetaKTSP: a meta-analytically derived
#' k-top-scoring-pair classifier for robust prediction of human cancer
#' subtypes. \emph{Bioinformatics} 32: 1966-1973. PMID: 27153719.
#'
#' Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
#' expression profiles from pairwise mRNA comparisons. \emph{Statistical
#' Applications in Genetics and Molecular Biology} 3: Article19.
#'
#' Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. (2005) Simple decision rules
#' for classifying human cancers from gene expression profiles. \emph{Bioinformatics}
#' 21: 3896-3904.
#'
#' @export
score_reo_metaktsp <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_metaktsp_model")) {
    stop("score_reo_metaktsp: model must have class reo_metaktsp_model")
  }
  X <- .reo_check_matrix(X, "score_reo_metaktsp", "X")
  .reo_check_meta(meta, nrow(X), "score_reo_metaktsp", "meta")

  usable_pairs <- .reo_ktsp_usable_pairs(model$pairs, colnames(X))
  if (nrow(usable_pairs) == 0L) {
    return(rep(0.5, nrow(X)))
  }

  out <- .reo_ktsp_vote_fraction(X, usable_pairs)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_reo_metaktsp: scorer produced non-finite or wrong-length output")
  }
  out
}

.reo_resolve_hp_metaktsp <- function(hp, p) {
  if (!is.list(hp)) stop("fit_reo_metaktsp: hp must be a list")
  allowed <- c("k", "cohort_col")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_reo_metaktsp: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  k <- hp$k
  if (is.null(k)) k <- 11L
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k < 1L || k != as.integer(k)) {
    stop("fit_reo_metaktsp: hp$k must be a positive integer")
  }

  cohort_col <- hp$cohort_col
  if (is.null(cohort_col)) cohort_col <- "accession"
  if (!is.character(cohort_col) || length(cohort_col) != 1L ||
      is.na(cohort_col) || !nzchar(cohort_col)) {
    stop("fit_reo_metaktsp: hp$cohort_col must be a single non-empty character string")
  }

  list(k = as.integer(k), cohort_col = cohort_col)
}

.reo_metaktsp_single_cohort_model <- function(X_train, y, hp, rows = NULL) {
  if (is.null(rows)) {
    rows <- seq_len(nrow(X_train))
  } else {
    rows <- which(rows)
  }
  if (length(rows) == 0L || !any(y[rows] == 1L) || !any(y[rows] == 0L)) {
    rows <- seq_len(nrow(X_train))
  }

  ktsp <- os_ktsp_fit(X_train[rows, , drop = FALSE], y[rows], k = hp$k)
  model <- list(
    pairs = ktsp$pairs,
    feature_universe = colnames(X_train),
    k = ktsp$k,
    cohort_col = hp$cohort_col,
    n_cohorts = 1L,
    hp = hp
  )
  class(model) <- "reo_metaktsp_model"
  model
}

.reo_metaktsp_pairs <- function(X_train, y, cohort, usable_cohorts) {
  features <- colnames(X_train)
  pair_rows <- list()
  idx <- 1L
  for (a in seq_len(ncol(X_train) - 1L)) {
    for (b in (a + 1L):ncol(X_train)) {
      deltas <- vapply(usable_cohorts, function(cohort_i) {
        rows <- !is.na(cohort) & cohort == cohort_i
        case_rate <- mean(X_train[rows & y == 1L, a] <
                            X_train[rows & y == 1L, b], na.rm = TRUE)
        ctrl_rate <- mean(X_train[rows & y == 0L, a] <
                            X_train[rows & y == 0L, b], na.rm = TRUE)
        if (!is.finite(case_rate) || !is.finite(ctrl_rate)) return(NA_real_)
        case_rate - ctrl_rate
      }, numeric(1))
      contributing <- is.finite(deltas)
      n_cohorts_used <- sum(contributing)
      meta_score <- if (n_cohorts_used > 0L) mean(deltas[contributing]) else 0
      if (meta_score >= 0) {
        pair_rows[[idx]] <- data.frame(
          feature_a = features[a],
          feature_b = features[b],
          orientation = "<",
          score = meta_score,
          stringsAsFactors = FALSE
        )
      } else {
        pair_rows[[idx]] <- data.frame(
          feature_a = features[b],
          feature_b = features[a],
          orientation = "<",
          score = -meta_score,
          stringsAsFactors = FALSE
        )
      }
      idx <- idx + 1L
    }
  }
  do.call(rbind, pair_rows)
}
