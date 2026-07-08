#' @title 1D wavelet-scattering single-sample discriminator (python-at-score, frozen head)
#'
#' @description
#' Single-sample within-cohort discriminator (family M, image/invariance) that
#' maps each specimen's frozen-order per-sample robust-CLR (rCLR) profile to a
#' \strong{1D wavelet scattering transform} (Mallat 2012; Anden & Mallat 2014)
#' and scores the scattering coefficients through a FROZEN regularized-logistic
#' head. A specimen of \eqn{p} universe features is laid out as a 1D signal along
#' the frozen feature axis, zero-padded to \eqn{\mathrm{shape} =
#' \mathrm{next\_pow2}(p)}, scattered by a fixed cascade of wavelet convolutions +
#' modulus + low-pass averaging (\code{ScatteringTorch1D}, kymatio), flattened to a
#' fixed-length coefficient vector, standardized by FROZEN training statistics, and
#' fed to a FROZEN ridge-logistic linear head. Larger = more case-like.
#'
#' \strong{Python-at-score contract (no live pointer).} Like the proven
#' python-at-score siblings \code{\link{fit_img_gasfcnn}} / \code{\link{fit_tabicl}}
#' (and unlike the export-to-pure-R trained-torch siblings \code{proto-net} /
#' \code{lrt-deepmaha}), the scattering transform is computed in python at fit AND
#' score, but the fitted model holds ONLY plain R-side state: the scattering CONFIG
#' (\code{J}, \code{Q}, \code{shape}, \code{max_order}, the frozen feature
#' universe), the FROZEN per-coordinate coefficient standardization
#' (\code{coef_mean}, \code{coef_sd}), the FROZEN linear head (\code{weights},
#' \code{intercept}), and the seed/hp. There is NO live python pointer on the
#' model. At SCORE the identical \code{ScatteringTorch1D} is REBUILT in python from
#' the config each call, the coefficients are computed FORCED ROW-BY-ROW in
#' float64, standardized with the frozen stats, and the frozen R-side linear head
#' is applied in pure base R. Because every model field is a plain R object, the
#' default \code{model_digest} is a stable deterministic snapshot suitable for
#' \code{\link{singlesample_assert_row_equivariant}}.
#'
#' \strong{Why this is inherently single-sample (no BatchNorm).} The scattering
#' transform \eqn{S(x)} is a FIXED cascade of wavelet convolutions, modulus, and
#' low-pass averaging with NO learned parameters and NO cross-sample statistics --
#' each specimen's coefficients depend ONLY on that specimen's own signal. There is
#' no BatchNorm and no batch coupling. The only numerical wrinkle is float32
#' convolution accumulation across a batch dimension; this scorer computes the
#' coefficients in \strong{float64} ONE ROW AT A TIME, so the \eqn{n=1} path is used
#' uniformly and a row's coefficients (hence its score) are bit-identical whether
#' scored alone or in any batch -- single-row score == batch score EXACTLY 0. The
#' coefficient standardization uses FROZEN training mean/sd (never the query's own),
#' and the head is a frozen linear map, so the whole pipeline is single-sample-safe.
#'
#' \strong{Signal construction (deterministic, frozen).} For a specimen's
#' universe-aligned per-sample rCLR vector \eqn{v} (length \eqn{p = } size of the
#' FROZEN feature universe, in the FROZEN feature ORDER), the 1D signal is \eqn{v}
#' laid out along the frozen feature axis, zero-padded to \eqn{\mathrm{shape} =
#' \mathrm{next\_pow2}(p)} (FROZEN at fit from the training universe size). The
#' scattering transform is ORDER-dependent along the feature axis: this is an
#' intended, frozen design choice (the feature order is fixed at fit; there is NO
#' claim of feature-permutation invariance). The row-equivariance harness permutes
#' SAMPLES (rows), never features (columns) -- exactly as for img-gasfcnn's GASF.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly
#' invariant to per-specimen positive scaling and uses no cross-row statistic.
#' Universe features absent from a specimen carry the neutral rCLR value \code{0}.
#' The SAME rCLR transform feeds the scattering at fit and at score.
#'
#' \strong{Degenerate-neutral.} A query returns the neutral score \code{0} if fewer
#' than \code{hp$min_features} universe features overlap the query columns (a
#' column-overlap floor, batch-independent), or it has empty positive support over
#' the frozen universe on the (pre-rCLR) aligned ORIGINAL abundances
#' (\code{!any(X_use[i, ] > 0)}). A FLAT all-equal-positive composition maps to the
#' rCLR origin (an all-zero signal) but is a VALID specimen: it yields valid (zero)
#' scattering coefficients and is scored normally (NOT floored) -- the empty-support
#' test is on the ORIGINAL abundances, NOT on \code{all(z == 0)}.
#'
#' \strong{kymatio import workaround.} kymatio 0.3.0's package init
#' (\code{import kymatio.torch}) eagerly imports the 3D harmonic-scattering module,
#' which calls \code{scipy.special.sph_harm} (REMOVED in scipy >= 1.15), so the
#' top-level import RAISES \code{ImportError}. The 1D frontend is therefore imported
#' DIRECTLY from \code{kymatio.scattering1d.frontend.torch_frontend} (verified
#' working: float64, deterministic). The reticulate skip-guard probes \code{torch}
#' AND that direct \code{ScatteringTorch1D} import (NOT \code{import kymatio.torch}).
#'
#' @references
#' Mallat S. (2012) Group Invariant Scattering. \emph{Communications on Pure and
#' Applied Mathematics} 65(10):1331-1398.
#'
#' Anden J, Mallat S. (2014) Deep Scattering Spectrum. \emph{IEEE Transactions on
#' Signal Processing} 62(16):4114-4128.
#'
#' Andreux M, Angles T, Exarchakis G, et al. (2020) Kymatio: Scattering Transforms
#' in Python. \emph{Journal of Machine Learning Research} 21(60):1-6.
#'
#' @name singlesample-inv-scatter
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .img_gasfcnn_rclr / .tabicl_rclr.
.inv_scatter_rclr <- function(v) {
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
.inv_scatter_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .inv_scatter_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.inv_scatter_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}

# Smallest power of two >= p (the FROZEN scattering signal length). p >= 1.
.inv_scatter_next_pow2 <- function(p) {
  s <- 1L
  while (s < p) s <- s * 2L
  as.integer(s)
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.inv_scatter_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_inv_scatter: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_inv_scatter: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_inv_scatter: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("J", "Q", "max_order", "lambda", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_inv_scatter: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  J <- hp$J
  if (is.null(J)) J <- 3L
  if (!is.numeric(J) || length(J) != 1L || !is.finite(J) || J < 1L ||
      J > .Machine$integer.max || J != as.integer(J)) {
    stop("fit_inv_scatter: hp$J must be a positive integer")
  }

  Q <- hp$Q
  if (is.null(Q)) Q <- 4L
  if (!is.numeric(Q) || length(Q) != 1L || !is.finite(Q) || Q < 1L ||
      Q > .Machine$integer.max || Q != as.integer(Q)) {
    stop("fit_inv_scatter: hp$Q must be a positive integer")
  }

  max_order <- hp$max_order
  if (is.null(max_order)) max_order <- 2L
  if (!is.numeric(max_order) || length(max_order) != 1L || !is.finite(max_order) ||
      !max_order %in% c(1L, 2L)) {
    stop("fit_inv_scatter: hp$max_order must be 1 or 2")
  }

  # NULL lambda => train-only deterministic CV over a fixed grid selects it. A
  # supplied lambda is a single non-negative finite number (ridge penalty).
  lambda <- hp$lambda
  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
        lambda < 0) {
      stop("fit_inv_scatter: hp$lambda must be NULL or a non-negative finite number")
    }
    lambda <- as.numeric(lambda)
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_inv_scatter: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_inv_scatter: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_inv_scatter: hp$seed must be a single integer")
  }

  list(
    J = as.integer(J),
    Q = as.integer(Q),
    max_order = as.integer(max_order),
    lambda = lambda,                    # NULL => CV-select at fit
    min_features = as.integer(min_features),
    device = device,
    seed = as.integer(seed)
  )
}


# ----------------------------------------------------------------------------
# Python / torch + kymatio availability + device resolution + scattering build
# ----------------------------------------------------------------------------

# Clear, actionable error if reticulate, torch, or the kymatio 1D frontend is
# unavailable. kymatio's top-level `import kymatio.torch` is BROKEN under
# scipy >= 1.15 (it eagerly imports the 3D module which calls the removed
# scipy.special.sph_harm); we therefore probe the DIRECT 1D-frontend import that
# the scorer actually uses, NOT `import kymatio.torch`.
.inv_scatter_require_modules <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_inv_scatter: package 'reticulate' is required to compute the ",
         "scattering transform. Install with install.packages('reticulate').",
         call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_inv_scatter: Python module 'torch' is not available in the active ",
         "reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  ok <- tryCatch({
    reticulate::import("kymatio.scattering1d.frontend.torch_frontend",
                       delay_load = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    stop("fit_inv_scatter: kymatio's 1D torch frontend ",
         "(kymatio.scattering1d.frontend.torch_frontend.ScatteringTorch1D) is not ",
         "importable in the active reticulate environment. NOTE: the top-level ",
         "'import kymatio.torch' is broken under scipy >= 1.15; install kymatio ",
         "and import the 1D frontend directly.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.inv_scatter_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}

# Build the FROZEN float64 ScatteringTorch1D operator from the config on `device`.
# One place builds the identical operator at fit and at score (rebuilt from the
# stored config). The 1D frontend is imported DIRECTLY (see
# .inv_scatter_require_modules). The operator has NO learned parameters.
.inv_scatter_build_scattering <- function(J, Q, shape, max_order, device) {
  kf <- reticulate::import("kymatio.scattering1d.frontend.torch_frontend",
                           delay_load = FALSE)
  S <- kf$ScatteringTorch1D(
    J = as.integer(J),
    Q = as.integer(Q),
    shape = as.integer(shape),
    max_order = as.integer(max_order)
  )
  S$double()$to(device)
}

# Compute the flattened float64 scattering coefficients of the n x p rCLR matrix
# `Z`, FORCED ROW-BY-ROW. Each row is zero-padded to length `shape` and scattered
# with a batch dimension of 1, so the n=1 path is used uniformly and a row's
# coefficients are bit-identical alone or in any batch (single == batch EXACTLY 0).
# Returns an n x D numeric matrix; D is fixed by the config. torch RNG is NOT
# touched here (the scattering forward consumes no RNG); the caller guards it.
.inv_scatter_coeffs <- function(S, Z, shape, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  n <- nrow(Z)
  p <- ncol(Z)
  no_grad <- torch$no_grad()
  no_grad$`__enter__`()
  on.exit(tryCatch(no_grad$`__exit__`(NULL, NULL, NULL), error = function(e) NULL),
          add = TRUE)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    sig <- numeric(shape)
    sig[seq_len(p)] <- Z[i, ]                       # zero-pad to shape (frozen)
    x <- torch$as_tensor(np$asarray(matrix(sig, nrow = 1L), dtype = "float64"))
    x <- x$to(torch$float64)$to(device)
    co <- S(x)$reshape(list(1L, -1L))               # 1 x D float64
    rows[[i]] <- as.numeric(reticulate::py_to_r(co$cpu()$numpy()))
  }
  out <- do.call(rbind, rows)
  storage.mode(out) <- "double"
  out
}


# ----------------------------------------------------------------------------
# FROZEN coefficient standardization + FROZEN ridge-logistic head (pure base R)
# ----------------------------------------------------------------------------

# Per-coordinate FROZEN standardization over the TRAINING coefficients. A
# zero-(or near-zero-)sd coordinate (constant across training) is guarded to
# sd = 1 so its standardized value is a finite constant 0 for every query.
.inv_scatter_standardize_fit <- function(C) {
  mu <- colMeans(C)
  sdv <- apply(C, 2L, stats::sd)
  sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  list(coef_mean = as.numeric(mu), coef_sd = as.numeric(sdv))
}

# Apply the FROZEN standardization to an n x D coefficient matrix.
.inv_scatter_standardize_apply <- function(C, coef_mean, coef_sd) {
  sweep(sweep(C, 2L, coef_mean, "-"), 2L, coef_sd, "/")
}

# FROZEN ridge-logistic head over the standardized coefficients, fit by IRLS with
# an L2 (ridge) penalty `lambda` on the WEIGHTS (the intercept is unpenalized).
# D can exceed n, so the ridge is REQUIRED (the un-penalized Hessian would be
# singular). Returns list(intercept, weights) -- a plain pure-R linear map; no RNG.
.inv_scatter_ridge_logit <- function(Xs, y, lambda, max_iter = 200L, tol = 1e-10) {
  D <- ncol(Xs)
  if (D > 0L) {
    sv <- svd(Xs, nu = 0L, nv = min(dim(Xs)))
    # Standard numerical-rank cutoff. Exact for exact-null directions (the primal
    # ridge optimum has a strictly-zero component there); for a near-null direction
    # 0 < d_k <= rank_tol that is dropped, the residual vs the primal is O(d_k/lambda)
    # -> negligible at the CV grid floor lambda >= 1e-3 and ~0 for in-distribution
    # queries. So this preserves the ridge optimum up to the numerical-rank tolerance,
    # not literally for every finite matrix (see the lambda=0 note at `force_jitter`).
    rank_tol <- if (length(sv$d) > 0L) {
      sv$d[1L] * max(dim(Xs)) * .Machine$double.eps
    } else {
      0
    }
    keep <- sv$d > rank_tol
    r <- sum(keep)
    V <- if (r > 0L) sv$v[, keep, drop = FALSE] else matrix(numeric(0L), D, 0L)
  } else {
    r <- 0L
    V <- matrix(numeric(0L), 0L, 0L)
  }
  Psi <- if (r > 0L) Xs %*% V else matrix(numeric(0L), nrow(Xs), 0L)
  Xa <- cbind(1, Psi)                                # n x (r+1), col 1 = intercept
  D1 <- ncol(Xa)
  pen <- rep(lambda, D1); pen[1L] <- 0               # do not penalize the intercept
  b <- rep(0, D1)
  # lambda=0 is NOT reachable in the benchmark (CV grid floor is 1e-3; default
  # lambda=NULL -> CV); it is only reachable via explicit user hp$lambda=0. With
  # dropped null directions (D>r) the primal H is singular and falls back to a
  # 1e-8-jittered solve; we force the same jitter here so the EXACT-null case
  # matches. NOTE: for lambda=0 on a RANK-DEFICIENT design the reduced (r+1) jitter
  # and the primal (D+1) jitter condition the singular system differently, so the
  # two are NOT bit-equivalent there (~1e-7); both are ~1e-8-ridge solutions of a
  # singular problem. This corner is documented/tested, not used by W7.
  force_jitter <- isTRUE(lambda == 0) && D > r        # omitted null directions made primal H singular
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
    bn <- if (force_jitter) {
      solve(H + diag(1e-8, D1), g)
    } else tryCatch(solve(H, g), error = function(e) {
      solve(H + diag(1e-8, D1), g)                   # tiny jitter if ill-conditioned
    })
    bn <- as.vector(bn)
    delta <- if (r > 0L) {
      max(abs(c(bn[1L] - b[1L], as.vector(V %*% (bn[-1L] - b[-1L])))))
    } else {
      abs(bn[1L] - b[1L])
    }
    if (delta < tol) { b <- bn; break }
    b <- bn
  }
  weights <- if (r > 0L) as.vector(V %*% b[-1L]) else numeric(D)
  list(intercept = b[1L], weights = weights)
}

# Deterministic class-stratified k-fold assignment (NO RNG): within each class,
# rows are dealt to folds round-robin in their existing order. Guarantees each
# fold sees both classes when both are present and k <= min class size.
.inv_scatter_cv_folds <- function(y, k) {
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
.inv_scatter_cv_lambda <- function(Xs, y, seed) {
  grid <- 10^seq(-3, 3, by = 0.5)                    # 13 candidate ridge penalties
  n <- nrow(Xs)
  k <- max(2L, min(5L, sum(y == 1L), sum(y == 0L)))  # <= smaller class size, >= 2
  fold <- .inv_scatter_cv_folds(y, k)
  dev_of <- function(lambda) {
    tot <- 0
    cnt <- 0L
    for (f in seq_len(k)) {
      te <- which(fold == f)
      tr <- which(fold != f)
      if (length(te) == 0L) next
      if (!any(y[tr] == 1L) || !any(y[tr] == 0L)) next  # fold lost a class: skip
      head <- .inv_scatter_ridge_logit(Xs[tr, , drop = FALSE], y[tr], lambda)
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
# fit_inv_scatter
# ----------------------------------------------------------------------------

#' @title Fit the 1D wavelet-scattering single-sample discriminator
#'
#' @description
#' Maps the TRAINING matrix to the per-sample robust CLR over the frozen feature
#' universe (\code{colnames(X_train)}), freezes the scattering signal length to
#' \code{shape = next_pow2(ncol(X_train))}, computes each training specimen's
#' float64 1D wavelet scattering coefficients (FORCED ROW-BY-ROW) via
#' reticulate-python torch + kymatio, freezes a per-coordinate coefficient
#' standardization, and fits a FROZEN ridge-logistic head over the standardized
#' coefficients. If \code{hp$lambda} is \code{NULL} the ridge penalty is selected
#' by TRAINING-ONLY deterministic class-stratified k-fold CV (no RNG). The fitted
#' model holds NO external pointer -- score rebuilds \code{ScatteringTorch1D} from
#' the stored config.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{J}
#'   (scattering scale, positive integer; default \code{3L}), \code{Q} (wavelets
#'   per octave, positive integer; default \code{4L}), \code{max_order}
#'   (scattering order, \code{1L} or \code{2L}; default \code{2L}), \code{lambda}
#'   (ridge penalty for the logistic head; \code{NULL} (default) selects it by
#'   training-only deterministic k-fold CV, otherwise a non-negative finite
#'   number), \code{min_features} (feature-overlap floor at scoring, positive
#'   integer; default \code{3L}), \code{device} (\code{"cpu"} (default),
#'   \code{"cuda"}, or \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU
#'   with no GPU), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{inv_scatter_model}: a list with
#'   \code{feature_universe}, \code{J}, \code{Q}, \code{shape}, \code{max_order},
#'   \code{coef_mean}, \code{coef_sd} (frozen standardization), \code{intercept},
#'   \code{weights} (frozen linear head), \code{lambda} (resolved ridge penalty),
#'   \code{device} (resolved), \code{seed}, and \code{hp}. No python pointer.
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
#' model <- fit_inv_scatter(X, y)
#' score_inv_scatter(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Anden J, Mallat S. (2014) Deep Scattering Spectrum. \emph{IEEE TSP}
#' 62(16):4114-4128.
#'
#' @export
fit_inv_scatter <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_inv_scatter", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_inv_scatter", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_inv_scatter")
  hp <- .inv_scatter_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_inv_scatter: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's first
  # python init seeds R's RNG. Restoring on exit (incl. the no-prior-seed case)
  # keeps fit globally RNG-neutral. The head's CV uses NO RNG (deterministic folds);
  # this guard makes fit robust even if a future reticulate init draws R randomness.
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

  .inv_scatter_require_modules()

  feat_order <- colnames(X_train)
  p <- ncol(X_train)
  shape <- .inv_scatter_next_pow2(p)
  Z <- .inv_scatter_rclr_matrix(X_train)             # n x p rCLR over full universe

  device <- .inv_scatter_resolve_device(hp$device)

  # Save/restore the torch RNG (CPU + guarded CUDA) around the python work so fit
  # leaves the torch generators byte-unchanged. The scattering forward consumes no
  # RNG; the guard makes RNG-safety robust regardless of build.
  torch <- reticulate::import("torch", delay_load = FALSE)
  torch_state <- tryCatch(torch$random$get_rng_state(), error = function(e) NULL)
  cuda_ok <- tryCatch(isTRUE(torch$cuda$is_available()), error = function(e) FALSE)
  cuda_state <- if (cuda_ok) {
    tryCatch(torch$cuda$get_rng_state_all(), error = function(e) NULL)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(torch_state)) {
      tryCatch(torch$random$set_rng_state(torch_state), error = function(e) NULL)
    }
    if (!is.null(cuda_state)) {
      tryCatch(torch$cuda$set_rng_state_all(cuda_state), error = function(e) NULL)
    }
  }, add = TRUE)

  S <- .inv_scatter_build_scattering(hp$J, hp$Q, shape, hp$max_order, device)
  C <- .inv_scatter_coeffs(S, Z, shape, device)      # n x D training coefficients

  std <- .inv_scatter_standardize_fit(C)
  Cs <- .inv_scatter_standardize_apply(C, std$coef_mean, std$coef_sd)

  lambda <- hp$lambda
  if (is.null(lambda)) lambda <- .inv_scatter_cv_lambda(Cs, y, hp$seed)
  head <- .inv_scatter_ridge_logit(Cs, y, lambda)

  model <- list(
    feature_universe = feat_order,
    J = hp$J,
    Q = hp$Q,
    shape = shape,
    max_order = hp$max_order,
    coef_mean = std$coef_mean,
    coef_sd = std$coef_sd,
    intercept = head$intercept,
    weights = head$weights,
    lambda = lambda,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "inv_scatter_model"
  model
}


# ----------------------------------------------------------------------------
# score_inv_scatter
# ----------------------------------------------------------------------------

#' @title Score the 1D wavelet-scattering discriminator (forced row-by-row)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_inv_scatter}}. The \code{ScatteringTorch1D} operator is REBUILT
#' in python from the stored config (no python pointer survives a fit), then each
#' query is mapped to the per-sample robust CLR over the frozen feature universe
#' (absent universe features carry the neutral rCLR value \code{0}), zero-padded to
#' the frozen \code{shape}, scattered in float64 ONE ROW AT A TIME, standardized
#' with the FROZEN training statistics, and passed through the FROZEN ridge-logistic
#' linear head. The score is \eqn{\mathrm{intercept} + \sum_d w_d \cdot
#' \mathrm{std\_coeff}_d(x)}; larger = more case-like. Forcing the \eqn{n=1}
#' scattering path uniformly (float64) makes a row's coefficients -- hence its
#' score -- bit-identical whether scored alone or in any batch.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present
#' (a column-overlap floor), or empty positive support over the frozen universe on
#' the (pre-rCLR) aligned ORIGINAL abundances (\code{!any(X_use[i, ] > 0)}), return
#' the neutral score \code{0}. A FLAT all-equal-positive composition maps to the
#' rCLR origin (all-zero signal) but is a VALID specimen and is scored normally
#' (NOT floored). The score of a row depends only on that row and the frozen model,
#' and is exactly invariant to per-specimen positive scaling.
#'
#' @param model An \code{inv_scatter_model} object from \code{\link{fit_inv_scatter}}.
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
#' model <- fit_inv_scatter(X, y)
#' score_inv_scatter(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Anden J, Mallat S. (2014) Deep Scattering Spectrum. \emph{IEEE TSP}
#' 62(16):4114-4128.
#'
#' @export
score_inv_scatter <- function(model, X, meta = NULL) {
  if (!inherits(model, "inv_scatter_model")) {
    stop("score_inv_scatter: model must have class inv_scatter_model")
  }
  X <- .reo_check_matrix(X, "score_inv_scatter", "X")
  .reo_check_meta(meta, nrow(X), "score_inv_scatter", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  .inv_scatter_require_modules()
  torch <- reticulate::import("torch", delay_load = FALSE)
  device <- .inv_scatter_resolve_device(model$hp$device)

  X_use <- .inv_scatter_align(X, feat_order)
  shape <- model$shape

  # Rows with empty positive support on the ORIGINAL aligned abundances are floored
  # to neutral 0 (NOT scored). A FLAT composition has an all-zero rCLR but a valid
  # all-zero signal: it is scored normally. The test is on X_use, NOT on all(z==0).
  active <- which(apply(X_use, 1L, function(r) any(r > 0)))

  # Save/restore the torch RNG (CPU + guarded CUDA) around the whole scoring loop so
  # scoring leaves the torch generators byte-unchanged (the scattering forward
  # consumes no RNG, but the guard makes RNG-safety robust). No python pointer
  # survives the call.
  torch_state <- tryCatch(torch$random$get_rng_state(), error = function(e) NULL)
  cuda_ok <- tryCatch(isTRUE(torch$cuda$is_available()), error = function(e) FALSE)
  cuda_state <- if (cuda_ok) {
    tryCatch(torch$cuda$get_rng_state_all(), error = function(e) NULL)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(torch_state)) {
      tryCatch(torch$random$set_rng_state(torch_state), error = function(e) NULL)
    }
    if (!is.null(cuda_state)) {
      tryCatch(torch$cuda$set_rng_state_all(cuda_state), error = function(e) NULL)
    }
  }, add = TRUE)

  if (length(active) > 0L) {
    S <- .inv_scatter_build_scattering(model$J, model$Q, shape, model$max_order,
                                       device)
    Z <- .inv_scatter_rclr_matrix(X_use[active, , drop = FALSE])
    C <- .inv_scatter_coeffs(S, Z, shape, device)    # row-by-row float64 coeffs
    Cs <- .inv_scatter_standardize_apply(C, model$coef_mean, model$coef_sd)
    eta <- model$intercept + as.vector(Cs %*% model$weights)
    out[active] <- eta
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_inv_scatter: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen, score-determining R-side state. The model holds NO live
# external-pointer object (the python ScatteringTorch1D is rebuilt from config at
# score, never stored), so every field is a plain R object and this is a
# deterministic snapshot suitable as the `model_digest` for
# singlesample_assert_row_equivariant() (clause (d): the bytes must be identical
# before and after scoring). Two seed-matched fits produce identical digests.
.inv_scatter_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    J = model$J,
    Q = model$Q,
    shape = model$shape,
    max_order = model$max_order,
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
