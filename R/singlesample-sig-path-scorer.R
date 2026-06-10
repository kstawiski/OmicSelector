#' @title Truncated path-signature single-sample discriminator (pure base-R, frozen head)
#'
#' @description
#' Single-sample within-cohort discriminator (family M, image/invariance) that
#' maps each specimen's frozen-order per-sample robust-CLR (rCLR) profile to a
#' TIME-AUGMENTED 2D path and scores its \strong{truncated path signature}
#' (Chen 1957; Lyons 1998; the iterated-integral feature map) through a FROZEN
#' ridge-logistic head. A specimen of \eqn{p} universe features is lifted to the
#' 2D path \eqn{P_i = (t_i, v_i)} with the FROZEN time grid \eqn{t_i =
#' (i-1)/(p-1)} along the frozen feature axis and \eqn{v_i} the rCLR value; the
#' increments \eqn{\Delta_i = P_{i+1} - P_i} drive the CLOSED-FORM truncated
#' signature over the truncated tensor algebra \eqn{T^{(L)}(\mathbb{R}^2)},
#' computed in PURE BASE R via Chen's identity, flattened to a fixed-length
#' coefficient vector (the level-0 constant \eqn{1} is dropped), standardized by
#' FROZEN training statistics, and fed to a FROZEN ridge-logistic linear head.
#' Larger = more case-like.
#'
#' \strong{Pure base-R contract (no python at fit OR score).} The PRESPECIFIED
#' \code{dep_route} was \code{reticulate-iisignature}, but iisignature has no
#' wheels and FAILS to build on python 3.14 (build-isolation
#' \code{ModuleNotFoundError: numpy}), and esig's only py3.14 path needs
#' \code{pyrecombine} which also has no wheel. Rather than ship a fragile C++
#' python dependency for a clinic-deployable single-sample method, the PRODUCTION
#' scorer computes the truncated signature in CLOSED FORM in PURE BASE R (Chen's
#' identity for a piecewise-linear path). The fitted model holds ONLY plain R-side
#' state: the frozen feature universe, the signature config (\code{L}, signature
#' dimension \code{D}), the FROZEN per-coordinate coefficient standardization
#' (\code{coef_mean}, \code{coef_sd}), the FROZEN linear head (\code{weights},
#' \code{intercept}), and the resolved hyperparameters. There is NO python
#' dependency, NO GPU, NO external pointer at fit or at score. The signature is a
#' deterministic pure-R function of the path alone, so the default
#' \code{model_digest} is a stable deterministic snapshot suitable for
#' \code{\link{singlesample_assert_row_equivariant}}, and a row's signature ---
#' hence its score --- is bit-identical whether scored alone or in any batch
#' (single-row score == batch score EXACTLY 0).
#'
#' \strong{Why this is inherently single-sample.} The path signature \eqn{S(P)}
#' is a fixed deterministic function of a single specimen's own path with NO
#' learned parameters and NO cross-sample statistics --- each specimen's
#' coefficients depend ONLY on that specimen's own rCLR profile and the frozen
#' time grid. There is no batch coupling and no float32/accumulation artifact
#' (the closed form is exact base-R double arithmetic). The coefficient
#' standardization uses FROZEN training mean/sd (never the query's own), and the
#' head is a frozen linear map, so the whole pipeline is single-sample-safe.
#'
#' \strong{Signature vs log-signature.} This scorer uses the FULL truncated
#' SIGNATURE, not the log-signature. The log-signature is a non-redundant
#' compression (a Lie-algebra element in a Lyndon/Hall basis); the full signature
#' is a complete superset of the same information that the ridge head handles via
#' L2 regularization, and the closed-form full signature avoids the subtle
#' Lyndon/Hall log map. The manifest notes record "path/log-signature".
#'
#' \strong{Signal construction (deterministic, frozen).} For a specimen's
#' universe-aligned per-sample rCLR vector \eqn{v} (length \eqn{p =} size of the
#' FROZEN feature universe, in the FROZEN feature ORDER), the 2D path is the
#' time-augmentation \eqn{P_i = (t_i, v_i)} with \eqn{t_i = (i-1)/(p-1)}. The
#' signature is ORDER-dependent along the feature axis: this is an intended,
#' frozen design choice (the feature order is fixed at fit; there is NO claim of
#' feature-permutation invariance). The row-equivariance harness permutes SAMPLES
#' (rows), never features (columns).
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly
#' invariant to per-specimen positive scaling and uses no cross-row statistic.
#' Universe features absent from a specimen carry the neutral rCLR value \code{0}.
#' The SAME rCLR transform feeds the signature at fit and at score.
#'
#' \strong{Degenerate-neutral.} A query returns the neutral score \code{0} if
#' fewer than \code{hp$min_features} universe features overlap the query columns (a
#' column-overlap floor, batch-independent), or it has empty positive support over
#' the frozen universe on the (pre-rCLR) aligned ORIGINAL abundances
#' (\code{!any(X_use[i, ] > 0)}). A FLAT all-equal-positive composition maps to the
#' rCLR origin (an all-zero \eqn{v}) but is a VALID specimen scored normally (NOT
#' floored): its 2D path is non-trivial because the time coordinate \eqn{t} still
#' moves, so its signature is non-degenerate. The empty-support test is on the
#' ORIGINAL abundances, NOT on \code{all(v == 0)}.
#'
#' \strong{Reference oracle (test only).} roughpy 0.3.0 (py3.14 wheel) is used
#' ONLY as a skip-guarded \eqn{\S}7 reference engine to validate the pure-R
#' signature math + index order to \eqn{< 10^{-9}}; the PRODUCTION scorer has NO
#' python dependency.
#'
#' @references
#' Chen K-T. (1957) Integration of paths, geometric invariants and a generalized
#' Baker--Hausdorff formula. \emph{Annals of Mathematics} 65(1):163-178.
#'
#' Lyons T. (1998) Differential equations driven by rough signals. \emph{Revista
#' Matematica Iberoamericana} 14(2):215-310.
#'
#' Chevyrev I, Kormilitzin A. (2016) A primer on the signature method in machine
#' learning. \emph{arXiv:1603.03788}.
#'
#' @name singlesample-sig-path
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .inv_scatter_rclr.
.sig_path_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Row-wise rCLR of a matrix. rbind (not t(vapply)) keeps the n x p shape even
# when n == 1 or p == 1.
.sig_path_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .sig_path_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.sig_path_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Closed-form truncated path signature (pure base R, d = 2)
# ----------------------------------------------------------------------------
# A truncated tensor-algebra element over R^2 up to level L is stored as a list
# of L+1 numeric vectors: entry [[k+1]] is level k, a length-2^k vector holding
# the coefficients of the words (i_1,...,i_k) flattened ROW-MAJOR (last index
# fastest, i.e. C-order over the k tensor axes). Level 0 is the scalar slot.
# This is exactly roughpy's word order (verified in §7.1 to < 1e-9).

# Signature dimension D = sum_{k=1}^{L} d^k with the level-0 constant dropped.
.sig_path_dim <- function(L, d = 2L) {
  sum(d ^ seq_len(L))
}

# Graded tensor exponential of a single segment increment delta in R^2:
#   exp_otimes(delta) = sum_{k=0}^{L} delta^{otimes k} / k!
# level 0 = 1; level k = the k-fold outer product of delta with itself / k!.
# The k-fold row-major-flattened outer power equals the k-fold Kronecker power
# (as.vector(kronecker(...))); verified against the explicit tensor in §7.1.
.sig_path_seg_exp <- function(delta, L) {
  lev <- vector("list", L + 1L)
  lev[[1L]] <- 1.0
  cur <- 1.0
  for (k in seq_len(L)) {
    cur <- as.vector(kronecker(cur, delta))     # delta^{otimes k}, row-major
    lev[[k + 1L]] <- cur / factorial(k)
  }
  lev
}

# Truncated concatenation (Chen) product of two truncated tensors a, b up to L:
#   (a otimes b)_m = sum_{j+k=m} a_j otimes b_k,  truncated at level L.
# Associative and unital (level-0 scalar acts as the algebra unit). The flattened
# outer product a_j otimes b_k is as.vector(kronecker(a_j, b_k)) (row-major).
.sig_path_chen <- function(a, b, L, d = 2L) {
  out <- vector("list", L + 1L)
  for (m in 0:L) {
    acc <- if (m == 0L) 0.0 else numeric(d ^ m)
    for (j in 0:m) {
      k <- m - j
      acc <- acc + as.vector(kronecker(a[[j + 1L]], b[[k + 1L]]))
    }
    out[[m + 1L]] <- acc
  }
  out
}

# Full truncated signature of a piecewise-linear path given its (s x d) increment
# matrix: the ORDERED Chen product S_1 otimes S_2 otimes ... otimes S_s of the
# per-segment tensor exponentials. Returns the list-of-levels truncated tensor.
.sig_path_signature_levels <- function(incr, L, d = 2L) {
  s <- .sig_path_seg_exp(incr[1L, ], L)
  if (nrow(incr) >= 2L) {
    for (i in 2:nrow(incr)) {
      s <- .sig_path_chen(s, .sig_path_seg_exp(incr[i, ], L), L, d)
    }
  }
  s
}

# Flatten a truncated-signature levels list to the fixed-length D feature vector,
# DROPPING the level-0 constant (always 1). Levels 1..L concatenated in order.
.sig_path_flatten <- function(levels, L) {
  if (L < 1L) return(numeric(0))
  unlist(levels[2:(L + 1L)], use.names = FALSE)
}

# Time-augmented 2D path increments for one rCLR vector v (length p): points
# P_i = (t_i, v_i) with the FROZEN grid t_i = (i-1)/(p-1); Delta_i = P_{i+1}-P_i.
# Requires p >= 2 (a path needs at least one segment); fit enforces this.
.sig_path_increments <- function(v) {
  p <- length(v)
  t <- (seq_len(p) - 1) / (p - 1)
  P <- cbind(t, v)
  P[-1L, , drop = FALSE] - P[-p, , drop = FALSE]
}

# Truncated signature feature vector (length D) of one rCLR vector v.
.sig_path_features_one <- function(v, L) {
  incr <- .sig_path_increments(v)
  lev <- .sig_path_signature_levels(incr, L)
  .sig_path_flatten(lev, L)
}

# n x D signature feature matrix of an n x p rCLR matrix Z, computed ROW-BY-ROW
# (each row independent; no cross-row state). single-row == batch EXACTLY 0.
.sig_path_features <- function(Z, L) {
  D <- .sig_path_dim(L)
  rows <- lapply(seq_len(nrow(Z)), function(i) .sig_path_features_one(Z[i, ], L))
  out <- matrix(unlist(rows, use.names = FALSE), nrow = nrow(Z), ncol = D,
                byrow = TRUE)
  storage.mode(out) <- "double"
  out
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.sig_path_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_sig_path: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_sig_path: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_sig_path: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("L", "lambda", "min_features", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_sig_path: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  L <- hp$L
  if (is.null(L)) L <- 3L
  if (!is.numeric(L) || length(L) != 1L || !is.finite(L) ||
      !L %in% c(1L, 2L, 3L, 4L)) {
    stop("fit_sig_path: hp$L (signature truncation level) must be 1, 2, 3, or 4")
  }

  # NULL lambda => train-only deterministic CV over a fixed grid selects it. A
  # supplied lambda is a single non-negative finite number (ridge penalty).
  lambda <- hp$lambda
  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
        lambda < 0) {
      stop("fit_sig_path: hp$lambda must be NULL or a non-negative finite number")
    }
    lambda <- as.numeric(lambda)
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_sig_path: hp$min_features must be a positive integer")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_sig_path: hp$seed must be a single integer")
  }

  list(
    L = as.integer(L),
    lambda = lambda,                    # NULL => CV-select at fit
    min_features = as.integer(min_features),
    seed = as.integer(seed)
  )
}


# ----------------------------------------------------------------------------
# FROZEN coefficient standardization + FROZEN ridge-logistic head (pure base R)
# ----------------------------------------------------------------------------

# Per-coordinate FROZEN standardization over the TRAINING signatures. A
# zero-(or near-zero-)sd coordinate (constant across training) is guarded to
# sd = 1 so its standardized value is a finite constant 0 for every query.
.sig_path_standardize_fit <- function(C) {
  mu <- colMeans(C)
  sdv <- apply(C, 2L, stats::sd)
  sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  list(coef_mean = as.numeric(mu), coef_sd = as.numeric(sdv))
}

# Apply the FROZEN standardization to an n x D coefficient matrix.
.sig_path_standardize_apply <- function(C, coef_mean, coef_sd) {
  sweep(sweep(C, 2L, coef_mean, "-"), 2L, coef_sd, "/")
}

# FROZEN ridge-logistic head over the standardized signatures, fit by IRLS with
# an L2 (ridge) penalty `lambda` on the WEIGHTS (the intercept is unpenalized).
# D can exceed n, so the ridge is REQUIRED (the un-penalized Hessian would be
# singular). Returns list(intercept, weights) -- a plain pure-R linear map; no RNG.
.sig_path_ridge_logit <- function(Xs, y, lambda, max_iter = 200L, tol = 1e-10) {
  Xa <- cbind(1, Xs)                                 # n x (D+1), col 1 = intercept
  D1 <- ncol(Xa)
  pen <- rep(lambda, D1); pen[1L] <- 0               # do not penalize the intercept
  b <- rep(0, D1)
  for (it in seq_len(max_iter)) {
    eta <- as.vector(Xa %*% b)
    mu <- 1 / (1 + exp(-eta))
    w <- mu * (1 - mu)
    w[w < 1e-9] <- 1e-9                              # IRLS weight floor (stability)
    zwork <- eta + (y - mu) / w
    XtW <- t(Xa * w)
    H <- XtW %*% Xa
    diag(H) <- diag(H) + pen
    g <- XtW %*% zwork
    bn <- tryCatch(solve(H, g), error = function(e) {
      solve(H + diag(1e-8, D1), g)                   # tiny jitter if ill-conditioned
    })
    bn <- as.vector(bn)
    if (max(abs(bn - b)) < tol) { b <- bn; break }
    b <- bn
  }
  list(intercept = b[1L], weights = b[-1L])
}

# Deterministic class-stratified k-fold assignment (NO RNG): within each class,
# rows are dealt to folds round-robin in their existing order. Guarantees each
# fold sees both classes when both are present and k <= min class size.
.sig_path_cv_folds <- function(y, k) {
  fold <- integer(length(y))
  for (cl in c(0L, 1L)) {
    idx <- which(y == cl)
    if (length(idx) > 0L) {
      fold[idx] <- ((seq_along(idx) - 1L) %% k) + 1L
    }
  }
  fold
}

# Select the ridge `lambda` by TRAINING-ONLY deterministic k-fold CV over a fixed
# log-spaced grid, scoring each lambda by mean held-out binomial deviance (lower =
# better). Fully deterministic (no RNG): folds are class-stratified round-robin.
# Returns the chosen lambda (a scalar from the grid).
.sig_path_cv_lambda <- function(Xs, y, seed) {
  grid <- 10^seq(-3, 3, by = 0.5)                    # 13 candidate ridge penalties
  k <- max(2L, min(5L, sum(y == 1L), sum(y == 0L)))  # <= smaller class size, >= 2
  fold <- .sig_path_cv_folds(y, k)
  dev_of <- function(lambda) {
    tot <- 0
    cnt <- 0L
    for (f in seq_len(k)) {
      te <- which(fold == f)
      tr <- which(fold != f)
      if (length(te) == 0L) next
      if (!any(y[tr] == 1L) || !any(y[tr] == 0L)) next  # fold lost a class: skip
      head <- .sig_path_ridge_logit(Xs[tr, , drop = FALSE], y[tr], lambda)
      eta <- head$intercept + as.vector(Xs[te, , drop = FALSE] %*% head$weights)
      mu <- 1 / (1 + exp(-eta))
      eps <- 1e-12
      mu <- pmin(pmax(mu, eps), 1 - eps)
      tot <- tot - 2 * sum(y[te] * log(mu) + (1 - y[te]) * log(1 - mu))
      cnt <- cnt + length(te)
    }
    if (cnt == 0L) return(Inf)
    tot / cnt
  }
  devs <- vapply(grid, dev_of, numeric(1L))
  if (all(!is.finite(devs))) return(1)               # degenerate fallback
  grid[which.min(devs)]                              # smallest-lambda tie-break
}


# ----------------------------------------------------------------------------
# fit_sig_path
# ----------------------------------------------------------------------------

#' @title Fit the truncated path-signature single-sample discriminator
#'
#' @description
#' Maps the TRAINING matrix to the per-sample robust CLR over the frozen feature
#' universe (\code{colnames(X_train)}), computes each training specimen's
#' closed-form truncated path signature (time-augmented 2D path, level \code{L},
#' pure base R, ROW-BY-ROW), freezes a per-coordinate signature-coefficient
#' standardization, and fits a FROZEN ridge-logistic head over the standardized
#' coefficients. If \code{hp$lambda} is \code{NULL} the ridge penalty is selected
#' by TRAINING-ONLY deterministic class-stratified k-fold CV (no RNG). The fitted
#' model holds NO python dependency and NO external pointer.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames). At least 2
#'   features are required (a path needs at least one segment).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{L}
#'   (signature truncation level, one of \code{1L}, \code{2L}, \code{3L},
#'   \code{4L}; default \code{3L}; signature dimension \eqn{D = \sum_{k=1}^{L}
#'   2^k}, so \code{14} at \code{L = 3}), \code{lambda} (ridge penalty for the
#'   logistic head; \code{NULL} (default) selects it by training-only
#'   deterministic k-fold CV, otherwise a non-negative finite number),
#'   \code{min_features} (feature-overlap floor at scoring, positive integer;
#'   default \code{3L}), and \code{seed} (integer; default \code{42L}). NO
#'   \code{device} hp --- the method is pure base R.
#'
#' @return Object of class \code{sig_path_model}: a list with
#'   \code{feature_universe}, \code{L}, \code{D} (signature dimension),
#'   \code{coef_mean}, \code{coef_sd} (frozen standardization), \code{intercept},
#'   \code{weights} (frozen linear head), \code{lambda} (resolved ridge penalty),
#'   \code{seed}, and \code{hp}. No python dependency, no external pointer.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 80; p <- 16; k <- 6
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' model <- fit_sig_path(X, y)
#' score_sig_path(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Chevyrev I, Kormilitzin A. (2016) A primer on the signature method in machine
#' learning. \emph{arXiv:1603.03788}.
#'
#' @export
fit_sig_path <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_sig_path", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_sig_path", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_sig_path")
  hp <- .sig_path_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_sig_path: X_train must contain at least hp$min_features features")
  }
  if (ncol(X_train) < 2L) {
    stop("fit_sig_path: X_train must contain at least 2 features ",
         "(a time-augmented path needs at least one segment)")
  }

  # The signature is a deterministic pure-R function: the head's CV uses NO RNG
  # (class-stratified round-robin folds) and IRLS is deterministic, so fit draws
  # no R randomness. We still save/restore the global R .Random.seed defensively
  # so fit is provably RNG-neutral (§7.9); score likewise consumes no RNG.
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)

  feat_order <- colnames(X_train)
  Z <- .sig_path_rclr_matrix(X_train)                # n x p rCLR over full universe

  C <- .sig_path_features(Z, hp$L)                   # n x D training signatures
  D <- .sig_path_dim(hp$L)

  std <- .sig_path_standardize_fit(C)
  Cs <- .sig_path_standardize_apply(C, std$coef_mean, std$coef_sd)

  lambda <- hp$lambda
  if (is.null(lambda)) lambda <- .sig_path_cv_lambda(Cs, y, hp$seed)
  head <- .sig_path_ridge_logit(Cs, y, lambda)

  model <- list(
    feature_universe = feat_order,
    L = hp$L,
    D = D,
    coef_mean = std$coef_mean,
    coef_sd = std$coef_sd,
    intercept = head$intercept,
    weights = head$weights,
    lambda = lambda,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "sig_path_model"
  model
}


# ----------------------------------------------------------------------------
# score_sig_path
# ----------------------------------------------------------------------------

#' @title Score the truncated path-signature discriminator (pure base R, row-by-row)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_sig_path}}. Each query is mapped to the per-sample robust CLR
#' over the frozen feature universe (absent universe features carry the neutral
#' rCLR value \code{0}), lifted to the time-augmented 2D path, its closed-form
#' truncated path signature computed in PURE BASE R ONE ROW AT A TIME,
#' standardized with the FROZEN training statistics, and passed through the FROZEN
#' ridge-logistic linear head. The score is \eqn{\mathrm{intercept} + \sum_d w_d
#' \cdot \mathrm{std\_sig}_d(x)}; larger = more case-like. The signature is a
#' deterministic pure-R function of the row alone, so a row's score is
#' bit-identical whether scored alone or in any batch.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present
#' (a column-overlap floor), or empty positive support over the frozen universe on
#' the (pre-rCLR) aligned ORIGINAL abundances (\code{!any(X_use[i, ] > 0)}), return
#' the neutral score \code{0}. A FLAT all-equal-positive composition maps to the
#' rCLR origin (all-zero \eqn{v}) but is a VALID specimen scored normally (NOT
#' floored): the time coordinate still moves, so its path signature is
#' non-degenerate. The score of a row depends only on that row and the frozen
#' model, and is exactly invariant to per-specimen positive scaling.
#'
#' @param model A \code{sig_path_model} object from \code{\link{fit_sig_path}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity and
#'   ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values are
#'   more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_sig_path(X, y)
#' score_sig_path(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Chevyrev I, Kormilitzin A. (2016) A primer on the signature method in machine
#' learning. \emph{arXiv:1603.03788}.
#'
#' @export
score_sig_path <- function(model, X, meta = NULL) {
  if (!inherits(model, "sig_path_model")) {
    stop("score_sig_path: model must have class sig_path_model")
  }
  X <- .reo_check_matrix(X, "score_sig_path", "X")
  .reo_check_meta(meta, nrow(X), "score_sig_path", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .sig_path_align(X, feat_order)

  # Rows with empty positive support on the ORIGINAL aligned abundances are floored
  # to neutral 0 (NOT scored). A FLAT composition has an all-zero rCLR but a valid
  # non-trivial path (t still moves): it is scored normally. The test is on X_use,
  # NOT on all(z == 0).
  active <- which(apply(X_use, 1L, function(r) any(r > 0)))

  if (length(active) > 0L) {
    Z <- .sig_path_rclr_matrix(X_use[active, , drop = FALSE])
    C <- .sig_path_features(Z, model$L)              # row-by-row pure-R signatures
    Cs <- .sig_path_standardize_apply(C, model$coef_mean, model$coef_sd)
    eta <- model$intercept + as.vector(Cs %*% model$weights)
    out[active] <- eta
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_sig_path: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen, score-determining R-side state. The model holds NO live
# external-pointer object (the signature is pure base R), so every field is a
# plain R object and this is a deterministic snapshot suitable as the
# `model_digest` for singlesample_assert_row_equivariant() (clause (d): the bytes
# must be identical before and after scoring). Two seed-matched fits produce
# identical digests.
.sig_path_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    L = model$L,
    D = model$D,
    coef_mean = model$coef_mean,
    coef_sd = model$coef_sd,
    intercept = model$intercept,
    weights = model$weights,
    lambda = model$lambda,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
