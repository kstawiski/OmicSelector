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
