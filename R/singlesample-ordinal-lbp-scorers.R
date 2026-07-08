#' @title Fit ordinal Local Binary Pattern (inv-olbp) within-sample discriminator
#'
#' @description
#' Learns a frozen, single-sample discriminator from the \emph{ordinal texture}
#' of a specimen's abundance profile. The construction has three frozen stages,
#' all estimated from training data only:
#' \enumerate{
#'   \item \strong{Canonical feature order.} Features are arranged once, in
#'     \emph{decreasing} training column mean of \code{X_train} (ties broken by
#'     feature index, so the order is deterministic and locale-independent). This
#'     fixed permutation, \code{order_perm}, defines a circular sequence of
#'     feature positions shared by every specimen.
#'   \item \strong{Ordinal Local Binary Pattern (OLBP) histogram.} For one
#'     specimen its values are read off in the canonical order and converted to
#'     within-sample ranks (\code{rank(.,"average")}). At each position the
#'     centre rank is compared to its \code{P} circular neighbours (offsets
#'     \code{c(-(P/2):-1, 1:(P/2))}); each comparison contributes one bit
#'     (\code{1} when the centre rank is \eqn{\ge} the neighbour rank), and the
#'     \code{P} bits are packed -- in the listed offset order, \code{offsets[k]}
#'     carrying weight \eqn{2^{k-1}} -- into an integer code in
#'     \eqn{[0, 2^P-1]}. The length-\eqn{2^P} normalised histogram of these
#'     codes over all positions is the specimen's texture descriptor. Because it
#'     is built from within-sample ranks it is invariant to any per-sample
#'     monotone-increasing transform (library size, input amount, log/qPCR vs
#'     microarray vs NGS scaling).
#'   \item \strong{Frozen linear head.} Per-sample histograms of the training
#'     matrix are standardised by a frozen \code{center}/\code{scale} (each bin's
#'     training mean and standard deviation, the latter floored at \code{eps};
#'     zero-variance bins are deactivated and contribute \code{0}). A logistic
#'     regression (\code{glm.fit}, binomial) of the labels on the standardised
#'     active bins gives the head coefficients; if that fit errors, fails to
#'     converge, returns non-finite coefficients, or exhibits complete or
#'     quasi-complete separation (fitted probabilities pinned at 0 or 1, which
#'     would otherwise freeze a numerically unstable divergent fit), a
#'     diagonal-LDA mean-difference direction (case-minus-control mean on the
#'     standardised active bins, thresholded at the class midpoint) is used
#'     instead. Both orient \emph{larger = more case-like}.
#' }
#' Disease reshapes \emph{which} miRNAs dominate relative to their canonical
#' neighbours; the OLBP histogram is a compact, batch-robust summary of that
#' within-sample shape. No test data and no test-time batch statistic are used:
#' the feature order, neighbour offsets, histogram binning, standardisation, and
#' head are all frozen here from training only.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{neighbors} the number of
#'   circular OLBP neighbours \eqn{P} (default \code{4L}; must be an even integer
#'   in \code{[2, 8]}, giving \eqn{2^P} histogram bins), \code{min_features} the
#'   feature-overlap floor at scoring (default \code{max(neighbors + 1L, 4L)},
#'   positive integer), and \code{eps} the positive standard-deviation floor
#'   applied to active histogram bins in the head (only genuinely zero-variance
#'   bins are deactivated and contribute \code{0}; default \code{1e-8}).
#'   These are resolved once during fitting and are not tuned inside
#'   \code{fit_inv_olbp()}. No randomness is used, so fitting is deterministic.
#'
#' @return A plain list of class \code{inv_olbp_model} containing
#'   \code{feature_universe}, the frozen \code{order_perm} (feature names in
#'   canonical order), \code{neighbor_offsets}, \code{n_codes} (\eqn{= 2^P}), the
#'   frozen \code{head} (\code{type} -- \code{"glm"} or \code{"diag_lda"} --
#'   \code{center}, \code{scale}, \code{active}, \code{coef}, \code{intercept}),
#'   and the resolved \code{hp}.
#'
#' @details
#' The per-specimen score (see \code{\link{score_inv_olbp}}) is the frozen linear
#' predictor on the standardised OLBP histogram,
#' \deqn{S(x) = \beta_0 + \sum_{b} \beta_b \, \frac{h_b(x) - c_b}{s_b},}
#' where \eqn{h(x)} is the specimen's normalised code histogram, \eqn{c, s} are
#' the frozen bin centres/scales, and \eqn{\beta} the frozen head. Larger
#' \eqn{S} means more case-like. Every quantity except \eqn{h(x)} is frozen at
#' fit time, and \eqn{h(x)} depends only on \eqn{x}'s own ranks, so the score is
#' computable from a single specimen and is unchanged by per-sample scaling.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 40
#' X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
#'             dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
#' y <- rep(c(0, 1), each = 60)
#' # cases elevate a mid block, reordering the within-sample profile
#' X[y == 1, 11:16] <- X[y == 1, 11:16] + 150
#' model <- fit_inv_olbp(X, y)
#' score <- score_inv_olbp(model, X)
#' }
#'
#' @references
#' Ojala T, Pietikainen M, Maenpaa T. (2002) Multiresolution gray-scale and
#' rotation invariant texture classification with local binary patterns.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}
#' 24(7): 971-987.
#'
#' Geman D, d'Avignon C, Naiman DQ, Winslow RL. (2004) Classifying gene
#' expression profiles from pairwise mRNA comparisons. \emph{Statistical
#' Applications in Genetics and Molecular Biology} 3: Article19.
#'
#' @export
fit_inv_olbp <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_inv_olbp", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_inv_olbp", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_inv_olbp")
  hp <- .inv_olbp_resolve_hp(hp)

  D <- ncol(X_train)
  if (D < hp$min_features) {
    stop("fit_inv_olbp: X_train must contain at least hp$min_features features")
  }
  if (D <= hp$neighbors) {
    stop("fit_inv_olbp: X_train must contain more than hp$neighbors features")
  }

  # Frozen canonical order: decreasing training column mean, ties broken by
  # feature index (locale-independent, deterministic). order_perm holds the
  # feature NAMES in that fixed order and is reused unchanged at scoring.
  col_mean <- colMeans(X_train)
  order_idx <- order(-col_mean, seq_len(D))
  order_perm <- colnames(X_train)[order_idx]

  offsets <- .inv_olbp_offsets(hp$neighbors)
  n_codes <- 2L^hp$neighbors

  # Per-sample OLBP histograms of the training matrix, built through the SAME
  # helper used at scoring so the descriptor is defined once (no fit/score drift).
  Xo <- X_train[, order_perm, drop = FALSE]
  H <- t(vapply(seq_len(nrow(Xo)), function(i) {
    .inv_olbp_hist(Xo[i, ], offsets, n_codes)
  }, numeric(n_codes)))

  head <- .inv_olbp_fit_head(H, y, hp$eps)

  model <- list(
    feature_universe = colnames(X_train),
    order_perm = order_perm,
    neighbor_offsets = offsets,
    n_codes = n_codes,
    head = head,
    hp = hp
  )
  class(model) <- "inv_olbp_model"
  model
}


#' @title Score ordinal Local Binary Pattern (inv-olbp) within-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_inv_olbp}}. The model-universe features present in \code{X}
#' are taken in the frozen canonical order, and for one specimen its within-sample
#' ranks over exactly those present features are turned into the ordinal Local
#' Binary Pattern histogram (same neighbour offsets, same bit packing, same
#' \eqn{2^P} bins as at fitting; the circular neighbour wrap is recomputed on the
#' present length \eqn{M}). The histogram is standardised by the frozen
#' \code{center}/\code{scale} and passed through the frozen linear head. Larger
#' values are more case-like.
#'
#' Scoring uses only each row's own values plus the frozen model -- no test-batch
#' renormalisation, no cross-row coupling, no statistic estimated from \code{X}.
#' If the number of present model-universe features is below
#' \code{model$hp$min_features} or not greater than the neighbour count \eqn{P}
#' (so a circular OLBP code is undefined), the documented neutral score \code{0}
#' is returned for every row.
#'
#' @param model An \code{inv_olbp_model} object returned by
#'   \code{\link{fit_inv_olbp}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like. The score of a row depends only on that row's values and
#'   the frozen model, and is exactly invariant to per-sample positive scaling.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 40
#' X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
#'             dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
#' y <- rep(c(0, 1), each = 60)
#' X[y == 1, 11:16] <- X[y == 1, 11:16] + 150
#' model <- fit_inv_olbp(X, y)
#' score_inv_olbp(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Ojala T, Pietikainen M, Maenpaa T. (2002) Multiresolution gray-scale and
#' rotation invariant texture classification with local binary patterns.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}
#' 24(7): 971-987.
#'
#' @export
score_inv_olbp <- function(model, X, meta = NULL) {
  if (!inherits(model, "inv_olbp_model")) {
    stop("score_inv_olbp: model must have class inv_olbp_model")
  }
  X <- .reo_check_matrix(X, "score_inv_olbp", "X")
  .reo_check_meta(meta, nrow(X), "score_inv_olbp", "meta")

  P <- length(model$neighbor_offsets)
  # Present model-universe features in the FROZEN canonical order (intersect
  # keeps the order of its first argument, order_perm).
  ordered_present <- intersect(model$order_perm, colnames(X))
  M <- length(ordered_present)
  if (M < model$hp$min_features || M <= P) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, ordered_present, drop = FALSE]
  head <- model$head
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    h <- .inv_olbp_hist(X_use[i, ], model$neighbor_offsets, model$n_codes)
    .inv_olbp_head_predict(h, head)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_inv_olbp: scorer produced non-finite or wrong-length output")
  }
  out
}


.inv_olbp_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_inv_olbp: hp must be a list")
  allowed <- c("neighbors", "min_features", "eps")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_inv_olbp: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  neighbors <- hp$neighbors
  if (is.null(neighbors)) neighbors <- 4L
  if (!is.numeric(neighbors) || length(neighbors) != 1L ||
      !is.finite(neighbors) || neighbors > .Machine$integer.max ||
      neighbors != as.integer(neighbors)) {
    stop("fit_inv_olbp: hp$neighbors must be an even integer in [2, 8]")
  }
  neighbors <- as.integer(neighbors)
  if (neighbors < 2L || neighbors > 8L || neighbors %% 2L != 0L) {
    stop("fit_inv_olbp: hp$neighbors must be an even integer in [2, 8]")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- max(neighbors + 1L, 4L)
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_inv_olbp: hp$min_features must be a positive integer")
  }
  min_features <- as.integer(min_features)

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-8
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_inv_olbp: hp$eps must be a positive finite number")
  }

  list(
    neighbors = neighbors,
    min_features = min_features,
    eps = as.numeric(eps)
  )
}

# Circular neighbour offsets for a P-bit OLBP code: the P/2 positions before and
# the P/2 positions after each centre, in the order (-(P/2):-1, 1:(P/2)). The
# k-th offset carries bit weight 2^(k-1) at packing time.
.inv_olbp_offsets <- function(P) {
  half <- P %/% 2L
  as.integer(c(-(half:1L), 1L:half))
}

# Per-position OLBP codes for one specimen's values ALREADY arranged in the
# canonical order. Uses only this specimen's own values: within-sample ranks
# (average ties), then for each circular offset a bit = 1 when the centre rank is
# >= the neighbour rank. Bits are packed in the listed offset order (offsets[k]
# -> weight 2^(k-1)) into an integer code in [0, 2^P-1]. The circular neighbour
# index of position t under offset o is ((t + o - 1) mod M) + 1; because the
# caller guarantees M > P >= 2*max|o|, no position is ever its own neighbour.
.inv_olbp_codes <- function(a, offsets) {
  M <- length(a)
  r <- rank(a, ties.method = "average")
  pos <- seq_len(M)
  codes <- numeric(M)
  for (k in seq_along(offsets)) {
    o <- offsets[k]
    nb <- ((pos + o - 1L) %% M) + 1L
    codes <- codes + as.integer(r >= r[nb]) * 2^(k - 1L)
  }
  as.integer(codes)
}

# Normalised OLBP histogram (length n_codes = 2^P) for one specimen's ordered
# values. Mass 1/M per position, summing to 1; built from within-sample ranks
# only, so it is invariant to any per-sample monotone-increasing transform.
.inv_olbp_hist <- function(a, offsets, n_codes) {
  codes <- .inv_olbp_codes(a, offsets)
  tabulate(codes + 1L, nbins = n_codes) / length(a)
}

# Frozen linear head over the standardised OLBP histograms. center is the
# training per-bin mean; scale is the training per-bin sd FLOORED at eps so a
# near-constant bin cannot blow up its standardised value. Only genuinely
# zero-variance bins (sd == 0, e.g. an OLBP code that is constant across the
# whole training set) carry no information and are "inactive" -> forced to
# contribute 0. The active standardised bins feed a binomial glm.fit; on failure
# (error, non-convergence, non-finite coefficients, or complete/quasi-complete
# SEPARATION -- detected so an unstable, iteration-dependent separating fit is
# not frozen) the head falls back to a stable diagonal-LDA mean-difference
# direction (consistent with the sibling gsp-gft head). head$type records which
# branch was taken. Both orient larger = case-like. No randomness is used.
.inv_olbp_fit_head <- function(H, y, eps) {
  n_codes <- ncol(H)
  center <- colMeans(H)
  raw_sd <- apply(H, 2L, stats::sd)
  raw_sd[!is.finite(raw_sd)] <- 0
  active <- raw_sd > 0                 # only exactly-zero-variance bins inactive
  scale <- pmax(raw_sd, eps)           # floor the sd of every (incl. active) bin

  Z <- sweep(H, 2L, center, "-")
  Z <- sweep(Z, 2L, scale, "/")
  Z[, !active] <- 0

  coef <- numeric(n_codes)  # aligned to all bins; inactive bins keep 0
  intercept <- 0
  type <- "diag_lda"        # default when no bin is active (degenerate head)
  if (any(active)) {
    Za <- Z[, active, drop = FALSE]
    fit <- .inv_olbp_try_glm(Za, y)
    if (is.null(fit)) {
      mu1 <- colMeans(Za[y == 1L, , drop = FALSE])
      mu0 <- colMeans(Za[y == 0L, , drop = FALSE])
      beta <- mu1 - mu0
      intercept <- -sum(beta * (mu1 + mu0) / 2)
      coef[active] <- beta
      type <- "diag_lda"
    } else {
      intercept <- fit$intercept
      coef[active] <- fit$beta
      type <- "glm"
    }
  }

  list(
    type = type,
    center = center,
    scale = scale,
    active = active,
    coef = coef,
    intercept = as.numeric(intercept)
  )
}

# Attempt a binomial logistic fit of y on cbind(1, Za). Returns NULL (signalling
# the diagonal-LDA fallback) on any error, non-convergence, non-finite
# coefficient, or complete/quasi-complete separation. Separation makes glm.fit
# report convergence with finite but DIVERGENT, iteration-dependent coefficients
# that would freeze an unstable head; it is detected locale-independently by
# fitted probabilities pinned at 0/1, plus the canonical English separation
# warning as a fast path, so the stable diagonal-LDA direction is used instead
# (matching the sibling gsp-gft head).
.inv_olbp_try_glm <- function(Za, y) {
  x <- cbind(1, Za)
  warn <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm.fit(x = x, y = as.numeric(y), family = stats::binomial()),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  co <- fit$coefficients
  if (length(co) != ncol(x) || any(!is.finite(co))) return(NULL)
  if (any(grepl("fitted probabilities numerically 0 or 1|did not converge",
                warn))) {
    return(NULL)
  }
  mu <- fit$fitted.values
  tol <- sqrt(.Machine$double.eps)
  if (length(mu) && (any(mu < tol) || any(mu > 1 - tol))) return(NULL)
  list(intercept = unname(co[1L]), beta = unname(co[-1L]))
}

# Apply the frozen head to one specimen's OLBP histogram. Standardise with the
# frozen center/scale, zero out inactive bins, and return the linear predictor.
.inv_olbp_head_predict <- function(h, head) {
  z <- (h - head$center) / head$scale
  z[!head$active] <- 0
  head$intercept + sum(head$coef * z)
}
