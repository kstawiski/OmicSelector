#' @title Fit within-cohort ILR balance discriminator
#'
#' @description
#' Wraps \code{\link{ws_balance_ilr}} as a frozen single-sample discriminator.
#' The raw score is the row sum of valid sequential-binary-partition ILR
#' balances computed from each specimen's own feature values. Fitting learns
#' only a one-bit orientation from the training labels: if the mean raw score is
#' at least as high in cases as controls, the stored sign is \code{+1};
#' otherwise it is \code{-1}. No test data or test-time batch statistic is used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{sbp} sequential binary
#'   partition (default \code{\link{ws_default_sbp}()}), \code{pseudocount}
#'   (default \code{NULL}, using the primitive's per-sample default),
#'   \code{aggregate} (\code{"gmean"} or \code{"trimmed_gmean"}, default
#'   \code{"gmean"}), \code{min_balance_coverage} (default \code{0.8}), and
#'   \code{min_sbp_features} canonical SBP feature floor (default \code{3L}).
#'
#' @return A plain list of class \code{ws_balance_ilr_model} containing the
#'   frozen \code{sbp}, \code{pseudocount}, \code{aggregate},
#'   \code{min_balance_coverage}, \code{min_sbp_features}, orientation
#'   \code{sign}, \code{sbp_canonical_features}, \code{feature_universe}, and
#'   resolved \code{hp}.
#'
#' @examples
#' \dontrun{
#' sbp <- ws_default_sbp()
#' features <- unique(unlist(lapply(sbp, function(b) {
#'   c(b$numerator, b$denominator)
#' }), use.names = FALSE))
#' features <- setdiff(features, c("ALL_NON_RBC", "ALL_NON_PLT",
#'                                 "TOP_K_BY_ABUNDANCE", "TAIL_BY_ABUNDANCE"))
#' X <- matrix(stats::rgamma(40 * length(features), shape = 2), nrow = 40,
#'             dimnames = list(NULL, features))
#' y <- rep(c(0, 1), each = 20)
#' model <- fit_ws_balance_ilr(X, y)
#' score <- score_ws_balance_ilr(model, X)
#' }
#'
#' @references
#' Egozcue J. J., Pawlowsky-Glahn V. (2005) Groups of parts and their
#' balances in compositional data analysis. \emph{Mathematical Geology}
#' 37(7): 795-828.
#'
#' @export
fit_ws_balance_ilr <- function(X_train, y_train, meta_train = NULL,
                               hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ws_balance_ilr", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ws_balance_ilr",
                  "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ws_balance_ilr")
  hp <- .ws_balance_ilr_resolve_hp(hp)
  sbp_canonical_features <- .ws_balance_ilr_canonical_features(hp$sbp)

  raw_train <- .ws_balance_ilr_raw_score(
    X = X_train,
    sbp = hp$sbp,
    pseudocount = hp$pseudocount,
    aggregate = hp$aggregate,
    min_balance_coverage = hp$min_balance_coverage,
    min_sbp_features = hp$min_sbp_features,
    sbp_canonical_features = sbp_canonical_features,
    feature_universe = colnames(X_train)
  )
  if (length(raw_train) != nrow(X_train) || any(!is.finite(raw_train))) {
    stop("fit_ws_balance_ilr: raw scorer produced non-finite or wrong-length output")
  }

  sign <- if (mean(raw_train[y == 1L]) >= mean(raw_train[y == 0L])) 1 else -1

  model <- list(
    sbp = hp$sbp,
    pseudocount = hp$pseudocount,
    aggregate = hp$aggregate,
    min_balance_coverage = hp$min_balance_coverage,
    min_sbp_features = hp$min_sbp_features,
    sign = sign,
    sbp_canonical_features = sbp_canonical_features,
    feature_universe = colnames(X_train),
    hp = hp
  )
  class(model) <- "ws_balance_ilr_model"
  model
}


#' @title Score within-cohort ILR balance discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_ws_balance_ilr}}. The scorer first restricts the specimen
#' panel to features from the frozen training universe that are present in
#' \code{X}. If fewer than \code{model$min_sbp_features} canonical SBP features
#' are present, the specimen/panel receives the documented neutral score
#' \code{0}. Otherwise, \code{\link{ws_balance_ilr}} computes the ILR balance
#' matrix and columns with any non-\code{NA} value are summed per row with
#' \code{na.rm = TRUE}; rows with no valid balance therefore also receive
#' \code{0}. The stored orientation sign is then applied so higher scores are
#' more case-like.
#'
#' A specimen with fewer than the minimum canonical SBP features has a
#' degenerate low-information score. This is the single-sample analogue of the
#' benchmark's panel-coverage ineligibility, surfaced as finite \code{0}
#' instead of \code{NA} so foundation validators can require finite outputs.
#'
#' @param model A \code{ws_balance_ilr_model} object returned by
#'   \code{\link{fit_ws_balance_ilr}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Scores use
#'   only each row's own values plus the frozen model, with no scored-batch
#'   centering, quantiles, or renormalization.
#'
#' @examples
#' \dontrun{
#' sbp <- ws_default_sbp()
#' features <- unique(unlist(lapply(sbp, function(b) {
#'   c(b$numerator, b$denominator)
#' }), use.names = FALSE))
#' features <- setdiff(features, c("ALL_NON_RBC", "ALL_NON_PLT",
#'                                 "TOP_K_BY_ABUNDANCE", "TAIL_BY_ABUNDANCE"))
#' X <- matrix(stats::rgamma(40 * length(features), shape = 2), nrow = 40,
#'             dimnames = list(NULL, features))
#' y <- rep(c(0, 1), each = 20)
#' model <- fit_ws_balance_ilr(X, y)
#' score_ws_balance_ilr(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Egozcue J. J., Pawlowsky-Glahn V. (2005) Groups of parts and their
#' balances in compositional data analysis. \emph{Mathematical Geology}
#' 37(7): 795-828.
#'
#' @export
score_ws_balance_ilr <- function(model, X, meta = NULL) {
  if (!inherits(model, "ws_balance_ilr_model")) {
    stop("score_ws_balance_ilr: model must have class ws_balance_ilr_model")
  }
  X <- .reo_check_matrix(X, "score_ws_balance_ilr", "X")
  .reo_check_meta(meta, nrow(X), "score_ws_balance_ilr", "meta")

  raw <- .ws_balance_ilr_raw_score(
    X = X,
    sbp = model$sbp,
    pseudocount = model$pseudocount,
    aggregate = model$aggregate,
    min_balance_coverage = model$min_balance_coverage,
    min_sbp_features = model$min_sbp_features,
    sbp_canonical_features = model$sbp_canonical_features,
    feature_universe = model$feature_universe
  )
  out <- as.numeric(model$sign * raw)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ws_balance_ilr: scorer produced non-finite or wrong-length output")
  }
  out
}

.ws_balance_ilr_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ws_balance_ilr: hp must be a list")
  allowed <- c("sbp", "pseudocount", "aggregate", "min_balance_coverage",
               "min_sbp_features")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ws_balance_ilr: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  sbp <- hp$sbp
  if (is.null(sbp)) sbp <- ws_default_sbp()
  .ws_balance_ilr_validate_sbp(sbp)

  pseudocount <- hp$pseudocount
  if (!is.null(pseudocount) &&
      (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
       !is.finite(pseudocount) || pseudocount <= 0)) {
    stop("fit_ws_balance_ilr: hp$pseudocount must be NULL or a positive finite number")
  }

  aggregate <- hp$aggregate
  if (is.null(aggregate)) aggregate <- "gmean"
  aggregate <- match.arg(aggregate, c("gmean", "trimmed_gmean"))

  min_balance_coverage <- hp$min_balance_coverage
  if (is.null(min_balance_coverage)) min_balance_coverage <- 0.8
  if (!is.numeric(min_balance_coverage) ||
      length(min_balance_coverage) != 1L ||
      !is.finite(min_balance_coverage) ||
      min_balance_coverage < 0 || min_balance_coverage > 1) {
    stop("fit_ws_balance_ilr: hp$min_balance_coverage must be in [0, 1]")
  }

  min_sbp_features <- hp$min_sbp_features
  if (is.null(min_sbp_features)) min_sbp_features <- 3L
  if (!is.numeric(min_sbp_features) || length(min_sbp_features) != 1L ||
      !is.finite(min_sbp_features) || min_sbp_features < 1L ||
      min_sbp_features != as.integer(min_sbp_features)) {
    stop("fit_ws_balance_ilr: hp$min_sbp_features must be a positive integer")
  }

  list(
    sbp = sbp,
    pseudocount = pseudocount,
    aggregate = aggregate,
    min_balance_coverage = min_balance_coverage,
    min_sbp_features = as.integer(min_sbp_features)
  )
}

.ws_balance_ilr_validate_sbp <- function(sbp) {
  if (!is.list(sbp) || length(sbp) < 1L || is.null(names(sbp)) ||
      any(!nzchar(names(sbp)))) {
    stop("fit_ws_balance_ilr: hp$sbp must be a named list of balances")
  }
  sentinels <- c("ALL_NON_RBC", "ALL_NON_PLT", "TOP_K_BY_ABUNDANCE",
                 "TAIL_BY_ABUNDANCE")
  for (nm in names(sbp)) {
    b <- sbp[[nm]]
    if (!is.list(b) || is.null(b$numerator) || is.null(b$denominator) ||
        !is.character(b$numerator) || !is.character(b$denominator) ||
        length(b$numerator) < 1L || length(b$denominator) < 1L ||
        any(is.na(b$numerator)) || any(is.na(b$denominator)) ||
        any(!nzchar(b$numerator)) || any(!nzchar(b$denominator))) {
      stop("fit_ws_balance_ilr: every hp$sbp balance must have a non-empty ",
           "character numerator and denominator")
    }
    # Reject real-feature overlap between the two sides (sentinels excluded:
    # they denote computed feature sets, not literal features).
    overlap <- intersect(setdiff(b$numerator, sentinels),
                         setdiff(b$denominator, sentinels))
    if (length(overlap) > 0L) {
      stop("fit_ws_balance_ilr: hp$sbp balance '", nm,
           "' has overlapping numerator/denominator feature(s): ",
           paste(overlap, collapse = ", "))
    }
  }
  invisible(TRUE)
}

.ws_balance_ilr_canonical_features <- function(sbp) {
  sentinels <- c("ALL_NON_RBC", "ALL_NON_PLT", "TOP_K_BY_ABUNDANCE",
                 "TAIL_BY_ABUNDANCE")
  features <- unlist(lapply(sbp, function(b) c(b$numerator, b$denominator)),
                     use.names = FALSE)
  unique(features[!features %in% sentinels])
}

.ws_balance_ilr_raw_score <- function(X, sbp, pseudocount, aggregate,
                                      min_balance_coverage,
                                      min_sbp_features,
                                      sbp_canonical_features,
                                      feature_universe = NULL) {
  if (!is.null(feature_universe)) {
    universe_present <- intersect(feature_universe, colnames(X))
    if (length(universe_present) == 0L) return(rep(0, nrow(X)))
    X <- X[, universe_present, drop = FALSE]
  }

  n_sbp <- sum(colnames(X) %in% sbp_canonical_features)
  if (n_sbp < min_sbp_features) return(rep(0, nrow(X)))

  B <- ws_balance_ilr(
    X,
    balances = sbp,
    pseudocount = pseudocount,
    aggregate = aggregate,
    min_balance_coverage = min_balance_coverage
  )

  if (!is.matrix(B)) {
    if (all(is.na(B))) return(rep(0, nrow(X)))
    raw <- as.numeric(B)
  } else {
    valid_cols <- which(colSums(!is.na(B)) > 0L)
    if (length(valid_cols) == 0L) return(rep(0, nrow(X)))
    raw <- rowSums(B[, valid_cols, drop = FALSE], na.rm = TRUE)
  }
  raw <- as.numeric(raw)
  if (length(raw) != nrow(X) || any(!is.finite(raw))) {
    stop(".ws_balance_ilr_raw_score: produced non-finite or wrong-length output")
  }
  raw
}
