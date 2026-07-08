#' @title Fit curvature-scale-space (inv-css) within-sample discriminator
#'
#' @description
#' Learns a frozen, single-sample discriminator from the \emph{multi-scale
#' differential curvature} of a specimen's compositional shape. Each specimen is
#' reduced to a single 1-D curve -- the sorted robust-CLR (rCLR) abundance
#' profile -- and that curve is summarised by a Mokhtarian-Mackworth
#' curvature-scale-space (CSS) descriptor, which is then fed to a frozen linear
#' (LDA) head. The construction has three stages, all estimated from training
#' data only (or, like the scale ladder and kernels, fixed data-free constants):
#' \enumerate{
#'   \item \strong{Single-specimen 1-D curve.} A specimen's non-negative
#'     abundances are mapped to a within-sample robust-CLR vector: over the
#'     specimen's own positive support \eqn{P=\{j: v_j>0\}},
#'     \eqn{z_j = \log v_j - \mathrm{mean}_{k\in P}\log v_k}; structural zeros
#'     (\eqn{v_j=0}) are excluded because they encode detection-rate / sequencing
#'     depth, a batch artifact rather than shape (only structural zeros, keyed on
#'     \eqn{v_j=0}, are dropped; a present coordinate with \eqn{z_j=0} is kept).
#'     The present rCLR coordinates
#'     \eqn{z_{pos}} are resampled to a fixed length \eqn{L} by evaluating their
#'     empirical quantile function at \eqn{u_\ell=(\ell-0.5)/L},
#'     \eqn{\ell=1,\dots,L} (\code{stats::quantile}, \code{type = 7}). The result
#'     \eqn{g\in\mathbb{R}^L} is a monotone non-decreasing curve, approximately
#'     mean-zero (the midpoint quantile resample of the exactly-mean-zero
#'     \eqn{z_{pos}} need not preserve the mean exactly; the curvature descriptors
#'     are in any case invariant to a vertical shift of \eqn{g}), and
#'     \emph{exactly} invariant to per-sample positive scaling (rCLR removes
#'     the per-sample geometric mean, so multiplying \eqn{v} by any \eqn{c>0}
#'     leaves \eqn{z_{pos}}, hence \eqn{g}, unchanged).
#'   \item \strong{Curvature-scale-space descriptor.} A frozen geometric ladder
#'     of \eqn{S} scales \eqn{\sigma_s} runs from \eqn{\sigma_{min}} to
#'     \eqn{\sigma_{max}}. At each scale \eqn{g} is convolved with a discrete,
#'     sum-1 Gaussian kernel (half-width \eqn{h=\lceil 3\sigma_s\rceil},
#'     symmetric/"reflect" boundary padding) to give \eqn{g_s}; unit-grid central
#'     and second differences \eqn{g'_s, g''_s} give the planar curve curvature
#'     \eqn{\kappa_s = g''_s / (1 + g'^2_s)^{3/2}} at the interior grid points,
#'     summarised over the \strong{padding-free sub-interval} only (indices whose
#'     smoothing window lies entirely within the real signal, excluding the
#'     reflect-padding boundary artefact that a monotone curve's mirrored slope
#'     would otherwise inject).
#'     Each scale contributes four descriptors -- total absolute curvature
#'     \eqn{\sum|\kappa_s|}, curvature energy \eqn{\sum\kappa_s^2}, peak curvature
#'     \eqn{\max|\kappa_s|}, and the CSS inflection count (number of sign changes
#'     of \eqn{\kappa_s}) -- concatenated over scales into
#'     \eqn{\varphi\in\mathbb{R}^{4S}}. The scale ladder, kernels and difference
#'     stencils are data-free, so the descriptor is defined identically at fit and
#'     at scoring.
#'   \item \strong{Frozen standardize + LDA head.} Training descriptors are
#'     standardised by a frozen \code{center}/\code{scale} (each descriptor's
#'     training mean and standard deviation, the latter floored at \code{eps};
#'     exactly-zero-variance descriptors are deactivated and contribute \code{0}).
#'     A ridge-regularised linear-discriminant head is fit on the standardised
#'     descriptors: pooled within-class covariance \eqn{S_w}, ridged to positive
#'     definiteness (\eqn{S_w + (\mathrm{shrink}\cdot\overline{\mathrm{diag}} +
#'     \mathrm{eps})I}, the ridge raised until the Cholesky succeeds),
#'     \eqn{w = S_w^{-1}(\mu_{case}-\mu_{ctrl})} and
#'     \eqn{b = -\tfrac12 (\mu_{case}+\mu_{ctrl})^\top w}. The head orients
#'     \emph{larger = more case-like}.
#' }
#' Disease reshapes the \emph{shape} of the sorted compositional profile (a sharp
#' "shoulder" or kink versus a smooth ramp); the CSS descriptor is a compact,
#' batch-robust, scale-free summary of that multi-scale shape. No test data and
#' no test-time batch statistic are used: the rCLR is strictly within-sample, and
#' the scale ladder, kernels, standardisation and head are all frozen here from
#' training only.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters (exact-name reads, strict
#'   allowed-list): \code{n_scales} number of scales \eqn{S} (default \code{6L};
#'   integer \eqn{\ge 2}), \code{resample_len} fixed curve length \eqn{L}
#'   (default \code{64L}; integer \eqn{\ge 16}; a very small \eqn{L} relative to
#'   \code{sigma_max} thins the padding-free interior \eqn{[h+1, L-h-2]} and can
#'   drive the largest scales into the neutral \code{(0,0,0,0)} descriptor -- the
#'   default \eqn{L=64} keeps every default scale well-populated),
#'   \code{sigma_min} smallest scale
#'   (default \code{1.0}; \eqn{>0}), \code{sigma_max} largest scale (default
#'   \code{NULL} \eqn{\to} \code{resample_len/8}; must exceed \code{sigma_min}),
#'   \code{shrink} ridge fraction of the mean covariance diagonal (default
#'   \code{0.1}; \eqn{\ge 0}), \code{min_features} per-specimen floor on the
#'   number of present rCLR coordinates (default \code{8L}; integer \eqn{\ge 3}),
#'   \code{eps} positive standard-deviation / ridge floor (default \code{1e-6}),
#'   and \code{seed} (default \code{1L}; accepted for interface uniformity only --
#'   the method consumes no randomness).
#'
#' @return A plain list of class \code{inv_css_model} containing
#'   \code{feature_universe}, the frozen scale ladder \code{scales} (per-scale
#'   \code{sigma}, half-width \code{h}, normalised \code{kern}),
#'   \code{descriptor_dim} (\eqn{= 4S}), the frozen \code{head}
#'   (\code{center}, \code{scale}, \code{active}, \code{w}, \code{b},
#'   \code{ridge}), and the resolved \code{hp}.
#'
#' @details
#' The per-specimen score (see \code{\link{score_inv_css}}) is the frozen linear
#' predictor on the standardised CSS descriptor,
#' \deqn{S(x) = b + \sum_{c} w_c \, \frac{\varphi_c(x) - \mathrm{center}_c}
#'                                        {\mathrm{scale}_c},}
#' where \eqn{\varphi(x)} is the specimen's CSS descriptor. Every quantity except
#' \eqn{\varphi(x)} is frozen at fit time, and \eqn{\varphi(x)} depends only on
#' \eqn{x}'s own rCLR profile, so the score is computable from a single specimen
#' and is exactly invariant to per-sample positive scaling. Although it reads the
#' same 1-D substrate as the optimal-transport profile methods (the sorted rCLR
#' curve), the descriptor family here is a distinct multi-scale
#' differential-geometric (CSS) representation, not an optimal-transport
#' embedding.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 60
#' features <- paste0("miR-", sprintf("%03d", seq_len(p)))
#' n <- 160
#' y <- rep(c(0, 1), each = n / 2)
#' X <- matrix(0, n, p, dimnames = list(NULL, features))
#' ctrl_mu <- seq(2, 5, length.out = p)                       # smooth ramp shape
#' case_mu <- ifelse(seq_len(p) <= p / 2, 2.75, 4.25)         # sharp step shape
#' for (i in seq_len(n)) {
#'   mu <- if (y[i] == 1) case_mu else ctrl_mu
#'   X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
#' }
#' model <- fit_inv_css(X, y)
#' score <- score_inv_css(model, X)
#' }
#'
#' @references
#' Mokhtarian F, Mackworth AK. (1992) A theory of multiscale, curvature-based
#' shape representation for planar curves. \emph{IEEE Transactions on Pattern
#' Analysis and Machine Intelligence} 14(8): 789-805.
#'
#' Witkin AP. (1983) Scale-space filtering. \emph{Proceedings of the 8th
#' International Joint Conference on Artificial Intelligence} 2: 1019-1022.
#'
#' Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse compositional
#' technique reveals microbial perturbations. \emph{mSystems} 4(1): e00016-19.
#'
#' @export
fit_inv_css <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_inv_css", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_inv_css", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_inv_css")
  hp <- .inv_css_resolve_hp(hp)

  if (ncol(X_train) < hp$min_features) {
    stop("fit_inv_css: X_train must contain at least hp$min_features features")
  }

  # Frozen, data-free scale ladder + Gaussian kernels (identical at fit/score).
  scales <- .inv_css_kernels(hp$n_scales, hp$sigma_min, hp$sigma_max)

  # Per-specimen CSS descriptor over the FULL training universe, built through
  # the SAME helper used at scoring so the descriptor is defined once (no
  # fit/score drift). A specimen with fewer than min_features present rCLR
  # coordinates yields NULL and is dropped from the head fit.
  desc_list <- lapply(seq_len(nrow(X_train)), function(i) {
    .inv_css_descriptor(X_train[i, ], scales, hp$resample_len, hp$min_features)
  })
  keep <- !vapply(desc_list, is.null, logical(1))
  if (sum(keep) < 2L) {
    stop("fit_inv_css: too few training specimens have >= hp$min_features ",
         "present rCLR coordinates")
  }
  Phi <- do.call(rbind, desc_list[keep])
  yk <- y[keep]
  if (sum(yk == 1L) < 1L || sum(yk == 0L) < 1L) {
    stop("fit_inv_css: after the min_features filter both classes must retain ",
         ">= 1 specimen")
  }

  head <- .inv_css_fit_head(Phi, yk, hp$shrink, hp$eps)

  model <- list(
    feature_universe = colnames(X_train),
    scales = scales,
    descriptor_dim = ncol(Phi),
    head = head,
    hp = hp
  )
  class(model) <- "inv_css_model"
  model
}


#' @title Score curvature-scale-space (inv-css) within-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_inv_css}}. The model-universe features present in \code{X}
#' (\code{intersect(model$feature_universe, colnames(X))}) define the present
#' substrate; for one specimen its within-sample rCLR profile over exactly those
#' present features is resampled to the frozen length, summarised by the frozen
#' curvature-scale-space descriptor (same scale ladder, kernels and difference
#' stencils as at fitting), standardised by the frozen \code{center}/\code{scale}
#' and passed through the frozen linear head. Larger values are more case-like.
#'
#' Scoring uses only each row's own values plus the frozen model -- no test-batch
#' renormalisation, no cross-row coupling, no statistic estimated from \code{X}.
#' If fewer than \code{model$hp$min_features} model-universe features are present,
#' the documented neutral score \code{0} is returned for every row; a single
#' specimen whose own number of present rCLR coordinates falls below
#' \code{model$hp$min_features} likewise scores the neutral \code{0}.
#'
#' @param model An \code{inv_css_model} object returned by
#'   \code{\link{fit_inv_css}}.
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
#' p <- 60
#' features <- paste0("miR-", sprintf("%03d", seq_len(p)))
#' n <- 160
#' y <- rep(c(0, 1), each = n / 2)
#' X <- matrix(0, n, p, dimnames = list(NULL, features))
#' ctrl_mu <- seq(2, 5, length.out = p)
#' case_mu <- ifelse(seq_len(p) <= p / 2, 2.75, 4.25)
#' for (i in seq_len(n)) {
#'   mu <- if (y[i] == 1) case_mu else ctrl_mu
#'   X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
#' }
#' model <- fit_inv_css(X, y)
#' score_inv_css(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Mokhtarian F, Mackworth AK. (1992) A theory of multiscale, curvature-based
#' shape representation for planar curves. \emph{IEEE Transactions on Pattern
#' Analysis and Machine Intelligence} 14(8): 789-805.
#'
#' @export
score_inv_css <- function(model, X, meta = NULL) {
  if (!inherits(model, "inv_css_model")) {
    stop("score_inv_css: model must have class inv_css_model")
  }
  X <- .reo_check_matrix(X, "score_inv_css", "X")
  .reo_check_meta(meta, nrow(X), "score_inv_css", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, present, drop = FALSE]
  head <- model$head
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    phi <- .inv_css_descriptor(X_use[i, ], model$scales, model$hp$resample_len,
                               model$hp$min_features)
    if (is.null(phi)) return(0)
    .inv_css_head_predict(phi, head)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_inv_css: scorer produced non-finite or wrong-length output")
  }
  out
}


# ---------------------------------------------------------------------------
# Hyperparameter resolver (strict allowed-list, exact `[[ ]]` reads, per-field
# validation with integer-overflow guards). No randomness is consumed.
# ---------------------------------------------------------------------------
.inv_css_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_inv_css: hp must be a list")
  # Reject fully-unnamed / partially-unnamed / duplicated hp before the
  # allowed-list check: setdiff(names(hp), allowed) silently passes a NULL-names
  # list (e.g. list(8L)) and keeps only the last of a duplicated name.
  if (length(hp) > 0L) {
    nms <- names(hp)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("fit_inv_css: all hp entries must be named")
    }
    if (anyDuplicated(nms)) {
      stop("fit_inv_css: duplicated hp field(s): ",
           paste(unique(nms[duplicated(nms)]), collapse = ", "))
    }
  }
  allowed <- c("n_scales", "resample_len", "sigma_min", "sigma_max", "shrink",
               "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_inv_css: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  n_scales <- hp[["n_scales"]]
  if (is.null(n_scales)) n_scales <- 6L
  if (!is.numeric(n_scales) || length(n_scales) != 1L || !is.finite(n_scales) ||
      n_scales > .Machine$integer.max || n_scales != as.integer(n_scales)) {
    stop("fit_inv_css: hp$n_scales must be an integer >= 2")
  }
  n_scales <- as.integer(n_scales)
  if (n_scales < 2L) stop("fit_inv_css: hp$n_scales must be an integer >= 2")

  resample_len <- hp[["resample_len"]]
  if (is.null(resample_len)) resample_len <- 64L
  if (!is.numeric(resample_len) || length(resample_len) != 1L ||
      !is.finite(resample_len) || resample_len > .Machine$integer.max ||
      resample_len != as.integer(resample_len)) {
    stop("fit_inv_css: hp$resample_len must be an integer >= 16")
  }
  resample_len <- as.integer(resample_len)
  if (resample_len < 16L) {
    stop("fit_inv_css: hp$resample_len must be an integer >= 16")
  }

  sigma_min <- hp[["sigma_min"]]
  if (is.null(sigma_min)) sigma_min <- 1.0
  if (!is.numeric(sigma_min) || length(sigma_min) != 1L ||
      !is.finite(sigma_min) || sigma_min <= 0) {
    stop("fit_inv_css: hp$sigma_min must be a positive finite number")
  }
  sigma_min <- as.numeric(sigma_min)

  sigma_max <- hp[["sigma_max"]]
  if (is.null(sigma_max)) sigma_max <- resample_len / 8
  if (!is.numeric(sigma_max) || length(sigma_max) != 1L ||
      !is.finite(sigma_max) || sigma_max <= 0) {
    stop("fit_inv_css: hp$sigma_max must be a positive finite number")
  }
  sigma_max <- as.numeric(sigma_max)
  if (sigma_max <= sigma_min) {
    stop("fit_inv_css: hp$sigma_max must be greater than hp$sigma_min")
  }

  shrink <- hp[["shrink"]]
  if (is.null(shrink)) shrink <- 0.1
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0) {
    stop("fit_inv_css: hp$shrink must be a non-negative finite number")
  }
  shrink <- as.numeric(shrink)

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 8L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_inv_css: hp$min_features must be an integer >= 3")
  }
  min_features <- as.integer(min_features)
  if (min_features < 3L) {
    stop("fit_inv_css: hp$min_features must be an integer >= 3")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_inv_css: hp$eps must be a positive finite number")
  }
  eps <- as.numeric(eps)

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_inv_css: hp$seed must be a non-negative integer")
  }
  seed <- as.integer(seed)

  list(
    n_scales = n_scales,
    resample_len = resample_len,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    shrink = shrink,
    min_features = min_features,
    eps = eps,
    seed = seed
  )
}


# ---------------------------------------------------------------------------
# Frozen, data-free scale ladder + Gaussian kernels. sigma_s is a geometric
# ladder from sigma_min to sigma_max; each scale carries a discrete, sum-1
# Gaussian kernel of half-width h = ceil(3*sigma). Depends only on hp, so it is
# identical at fit and at score.
# ---------------------------------------------------------------------------
.inv_css_kernels <- function(n_scales, sigma_min, sigma_max) {
  if (n_scales < 2L) {
    sig <- sigma_min
  } else {
    sig <- sigma_min *
      (sigma_max / sigma_min)^((seq_len(n_scales) - 1L) / (n_scales - 1L))
  }
  lapply(sig, function(s) {
    h <- as.integer(ceiling(3 * s))
    if (h < 1L) h <- 1L
    idx <- seq.int(-h, h)
    k <- exp(-(idx^2) / (2 * s^2))
    k <- k / sum(k)
    list(sigma = s, h = h, kern = k)
  })
}


# ---------------------------------------------------------------------------
# Single-specimen within-sample robust-CLR over the specimen's own positive
# support, returning ALL present (v > 0) coordinates. Only STRUCTURAL zeros
# (absent parts, v == 0) are excluded -- NOT coordinates whose centred-log value
# happens to equal the geometric mean (z == 0; e.g. a tied / constant count
# profile), which are kept. Keying on `z != 0` here would wrongly drop those and
# would empty a constant positive profile, distorting the resampled curve.
# Exactly scale-invariant: log(c*v) - mean(log(c*v)) = log(v) - mean(log(v)).
# ---------------------------------------------------------------------------
.inv_css_rclr_pos <- function(v) {
  pos <- which(v > 0)
  if (length(pos) == 0L) return(numeric(0))
  lv <- log(v[pos])
  lv - mean(lv)
}


# Symmetric ("reflect", edge-duplicated) boundary index map with period 2L:
# indices 1..L,L..1 repeat, so any integer index (left/right pad of any width,
# even h > L) maps deterministically into [1, L]. R's `%%` returns a value with
# the sign of the (positive) divisor, so i2 is always in [0, 2L-1].
.inv_css_reflect_index <- function(i, L) {
  if (L == 1L) return(rep.int(1L, length(i)))
  m <- 2L * L
  i2 <- (i - 1L) %% m
  ifelse(i2 < L, i2 + 1L, m - i2)
}

# Reflect-pad g by h on each side using the symmetric index map above.
.inv_css_reflect <- function(g, h) {
  L <- length(g)
  j <- seq_len(L + 2L * h) - h          # original indices 1-h, ..., L+h
  g[.inv_css_reflect_index(j, L)]
}

# Convolve g with a symmetric sum-1 kernel under reflect padding, returning a
# length-L smoothed signal. Because the kernel sums to 1, a constant signal is
# returned unchanged (so a flat profile has ~0 curvature).
.inv_css_smooth <- function(g, kern, h) {
  L <- length(g)
  padded <- .inv_css_reflect(g, h)
  out <- numeric(L)
  K <- length(kern)                      # = 2h + 1
  for (k in seq_len(K)) {
    out <- out + kern[k] * padded[k:(k + L - 1L)]
  }
  out
}

# Unit-grid curvature of the planar curve (u, g_s(u)) at the interior points:
# kappa = g'' / (1 + g'^2)^(3/2), with central first and second differences.
# Returns a length-(L-2) vector.
.inv_css_curvature <- function(g_s) {
  L <- length(g_s)
  if (L < 3L) return(numeric(0))
  t <- 2:(L - 1L)
  g1 <- (g_s[t + 1L] - g_s[t - 1L]) / 2
  g2 <- g_s[t + 1L] - 2 * g_s[t] + g_s[t - 1L]
  g2 / (1 + g1^2)^(3 / 2)
}

# CSS inflection count: number of sign changes of the curvature, ignoring values
# at or below a near-zero noise floor `tol`. A flat / straight region has true
# curvature 0 but its finite-difference value is floating-point noise (~eps*|g|)
# that oscillates in sign; without a tolerance that noise would be counted as
# dozens of spurious inflections. Genuine zero-crossings are preserved (curvature
# is small only AT the crossing and grows away from it with a definite sign).
# `tol` is a machine-eps-relative floor on the curve amplitude, so the count stays
# bit-stable under per-sample positive scaling (the rCLR curve, hence kappa and
# tol, are scale-invariant). Default tol = 0 reproduces the exact-zeros-only form.
.inv_css_sign_changes <- function(kappa, tol = 0) {
  s <- sign(kappa)
  s[abs(kappa) <= tol] <- 0
  s <- s[s != 0]
  if (length(s) < 2L) return(0L)
  as.integer(sum(s[-1L] != s[-length(s)]))
}

# Four per-scale descriptors for a smoothed profile g_s: total absolute
# curvature, curvature energy, peak |curvature|, and the CSS inflection count --
# summarised over the PADDING-FREE INTERIOR only. The smoothed value g_s[i] is
# free of reflect-padding contamination only for i in [h+1, L-h] (its kernel
# window then lies entirely within the real signal), so curvature (which uses
# g_s[i-1,i,i+1]) is padding-free for i in [h+2, L-h-1], i.e. kappa[j] with
# j = i-1 in [h+1, L-h-2]. Restricting to this range drops the boundary artefact
# -- the sorted (monotone) curve's mirrored boundary slope injects a smoothed kink
# at each edge that otherwise DOMINATES sum/max curvature at medium-large scales
# (up to ~85% of total |kappa|, with max|kappa| always at the boundary) -- while
# keeping every genuine real-data curvature point. If the padding-free interior is
# empty (a very large scale relative to L), the scale contributes the neutral
# (0,0,0,0).
.inv_css_scale_descriptors <- function(g_s, h) {
  kappa <- .inv_css_curvature(g_s)
  L <- length(g_s)
  lo <- h + 1L
  hi <- L - h - 2L
  if (length(kappa) == 0L || hi < lo) return(c(0, 0, 0, 0))
  kappa <- kappa[lo:hi]
  ak <- abs(kappa)
  # Noise floor for the inflection count: the second-difference cancellation error
  # is ~ eps * max|g_s|, so curvature below ~1e3*eps*max|g_s| is fp noise, not a
  # genuine zero-crossing. Scale-invariant (g_s is scale-invariant).
  tol <- 1e3 * .Machine$double.eps * max(abs(g_s))
  c(sum(ak), sum(kappa^2), max(ak), .inv_css_sign_changes(kappa, tol))
}


# ---------------------------------------------------------------------------
# Full single-specimen CSS descriptor. Returns NULL when the specimen has fewer
# than min_features present rCLR coordinates (the documented neutral case).
# Used at fit (full universe) and at score (present subset) -- one definition.
# ---------------------------------------------------------------------------
.inv_css_descriptor <- function(v, scales, L, min_features) {
  zp <- .inv_css_rclr_pos(v)
  if (length(zp) < min_features) return(NULL)
  u <- (seq_len(L) - 0.5) / L
  g <- stats::quantile(zp, probs = u, type = 7, names = FALSE)
  d <- length(scales)
  phi <- numeric(d * 4L)
  for (s in seq_len(d)) {
    sc <- scales[[s]]
    g_s <- .inv_css_smooth(g, sc$kern, sc$h)
    phi[((s - 1L) * 4L + 1L):(s * 4L)] <- .inv_css_scale_descriptors(g_s, sc$h)
  }
  phi
}


# ---------------------------------------------------------------------------
# Frozen standardize + ridge-LDA head over the training CSS descriptors.
# center/scale are the per-descriptor training mean / sd (sd floored at eps);
# exactly-zero-variance descriptors are deactivated (forced to contribute 0).
# The pooled within-class covariance is ridged until its Cholesky succeeds, then
# w = Sinv (mu_case - mu_ctrl), b = -0.5 (mu_case + mu_ctrl) . w. Larger = case.
# ---------------------------------------------------------------------------
.inv_css_fit_head <- function(Phi, y, shrink, eps) {
  d <- ncol(Phi)
  center <- colMeans(Phi)
  raw_sd <- apply(Phi, 2L, stats::sd)
  raw_sd[!is.finite(raw_sd)] <- 0
  active <- raw_sd > 0
  scale <- pmax(raw_sd, eps)

  Z <- sweep(Phi, 2L, center, "-")
  Z <- sweep(Z, 2L, scale, "/")
  Z[, !active] <- 0

  Zc <- Z[y == 1L, , drop = FALSE]
  Z0 <- Z[y == 0L, , drop = FALSE]
  mu_case <- colMeans(Zc)
  mu_ctrl <- colMeans(Z0)
  df <- max(nrow(Z) - 2L, 1L)
  Wc <- crossprod(sweep(Zc, 2L, mu_case, "-"))
  W0 <- crossprod(sweep(Z0, 2L, mu_ctrl, "-"))
  Sw <- (Wc + W0) / df
  Sw <- (Sw + t(Sw)) / 2                  # symmetrize against fp asymmetry

  diag_mean <- mean(diag(Sw))
  if (!is.finite(diag_mean) || diag_mean < 0) diag_mean <- 0
  ridge <- shrink * diag_mean + eps
  if (!is.finite(ridge) || ridge <= 0) ridge <- eps

  ch <- NULL
  for (it in seq_len(64L)) {
    Swr <- Sw
    diag(Swr) <- diag(Swr) + ridge
    ch <- tryCatch(chol(Swr), error = function(e) NULL)
    if (!is.null(ch)) break
    ridge <- ridge * 10
  }
  if (is.null(ch)) {
    stop("fit_inv_css: within-class covariance could not be made positive ",
         "definite")
  }
  Sinv <- chol2inv(ch)

  delta <- mu_case - mu_ctrl
  w <- as.numeric(Sinv %*% delta)
  w[!active] <- 0                         # inactive descriptors carry no weight
  b <- -0.5 * sum((mu_case + mu_ctrl) * w)

  list(
    center = center,
    scale = scale,
    active = active,
    w = w,
    b = as.numeric(b),
    ridge = ridge
  )
}

# Apply the frozen head to one specimen's CSS descriptor: standardise with the
# frozen center/scale, zero inactive descriptors, return the linear predictor.
.inv_css_head_predict <- function(phi, head) {
  z <- (phi - head$center) / head$scale
  z[!head$active] <- 0
  head$b + sum(head$w * z)
}
