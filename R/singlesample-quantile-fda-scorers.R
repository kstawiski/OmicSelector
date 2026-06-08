#' @title Fit quantile-function FDA (inv-fdaqf) within-sample discriminator
#'
#' @description
#' Learns a frozen, single-sample discriminator from a functional-PCA (FPCA)
#' decomposition of a specimen's compositional \emph{quantile function} (QF).
#' Each specimen is reduced to a single monotone curve -- the sorted robust-CLR
#' (rCLR) abundance profile, resampled to a fixed length -- and that curve is
#' projected onto a frozen orthonormal FPCA eigenbasis, whose coordinates feed a
#' frozen ridge-regularised linear-discriminant (LDA) head. The construction has
#' four stages, all estimated from training data only (the QF grid is a fixed
#' data-free constant):
#' \enumerate{
#'   \item \strong{Single-specimen quantile function.} A specimen's non-negative
#'     abundances are mapped to a within-sample robust-CLR vector: over the
#'     specimen's own positive support \eqn{P=\{j: v_j>0\}},
#'     \eqn{z_j = \log v_j - \mathrm{mean}_{k\in P}\log v_k}. Structural zeros
#'     (\eqn{v_j=0}) are excluded because they encode detection-rate / sequencing
#'     depth, a batch artifact rather than shape (only structural zeros, keyed on
#'     \eqn{v_j=0}, are dropped; a present coordinate with \eqn{z_j=0}, e.g. a
#'     tied / constant count profile, is kept). The present rCLR coordinates
#'     \eqn{z_{pos}} are resampled to a fixed length \eqn{L} by evaluating their
#'     empirical quantile function at \eqn{u_\ell=(\ell-0.5)/L},
#'     \eqn{\ell=1,\dots,L} (\code{stats::quantile}, \code{type = 7}). The result
#'     \eqn{q\in\mathbb{R}^L} is a monotone non-decreasing curve of fixed length
#'     \eqn{L} regardless of the specimen's support size, and \emph{exactly}
#'     invariant to per-sample positive scaling (rCLR removes the per-sample
#'     geometric mean, so multiplying \eqn{v} by any \eqn{c>0} leaves
#'     \eqn{z_{pos}}, hence \eqn{q}, unchanged). This QF substrate is shared with
#'     the curvature-scale-space (inv-css) and optimal-transport profile
#'     methods, but the mechanism here is a distinct functional-PCA projection,
#'     not a curvature descriptor and not an optimal-transport embedding.
#'   \item \strong{Frozen functional-PCA basis (training only).} The training QF
#'     matrix \eqn{G} (\eqn{n_{kept}\times L}, one row per retained specimen QF)
#'     is centred by its functional mean \eqn{\bar q=\mathrm{colMeans}(G)}:
#'     \eqn{G_c = G - \mathbf{1}\bar q^\top}. The top-\eqn{K} right-singular
#'     vectors of \eqn{G_c} (base-R \code{svd}, no \pkg{fda} dependency) form the
#'     frozen orthonormal loading basis \eqn{V} (\eqn{L\times K}, with
#'     \eqn{V^\top V = I_K} to machine precision). \eqn{K=}\code{n_components},
#'     capped at the numerical rank of \eqn{G_c}: if fewer than \eqn{K} positive
#'     singular values exist only those are kept and the realised \eqn{K} is
#'     recorded (\code{model$fpca$K}). SVD column signs (and, for repeated
#'     singular values, the in-subspace rotation) are arbitrary but \emph{frozen}
#'     in the model and absorbed by the downstream linear head, so the final
#'     scores are sign-robust; the frozen \eqn{V} orientation is nonetheless
#'     LAPACK-build dependent (a reproducibility caveat shared with other frozen
#'     eigenbasis methods such as gsp-gft).
#'   \item \strong{Projection scores.} A specimen's FPCA score vector is the
#'     projection of its centred QF onto the frozen basis,
#'     \eqn{s = V^\top (q - \bar q)\in\mathbb{R}^K}. It uses only the specimen's
#'     own QF and the frozen \eqn{\bar q}, \eqn{V}.
#'   \item \strong{Frozen standardize + ridge-LDA head.} The \eqn{K}-dimensional
#'     score vectors are standardised by a frozen \code{center}/\code{scale}
#'     (each component's training mean and standard deviation, the latter floored
#'     at \code{eps}; exactly-zero-variance components are deactivated and
#'     contribute \code{0}). A ridge-regularised LDA head is then fit: pooled
#'     within-class covariance \eqn{S_w}, ridged to positive definiteness
#'     (\eqn{S_w + (\mathrm{shrink}\cdot\overline{\mathrm{diag}} +
#'     \mathrm{eps})I}, the ridge raised until the Cholesky succeeds),
#'     \eqn{w = S_w^{-1}(\mu_{case}-\mu_{ctrl})} and
#'     \eqn{b = -\tfrac12 (\mu_{case}+\mu_{ctrl})^\top w}. The head orients
#'     \emph{larger = more case-like}.
#' }
#' Disease reshapes the \emph{distributional shape} of the compositional profile
#' (e.g. a heavier upper tail in the rCLR values); the FPCA score vector is a
#' compact, batch-robust, scale-free summary of that shape along the dominant
#' training modes of variation. No test data and no test-time batch statistic are
#' used: the rCLR is strictly within-sample, and the QF grid, functional mean,
#' eigenbasis, standardisation and head are all frozen here from training only.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters (exact-name reads, strict
#'   allowed-list): \code{resample_len} fixed QF length \eqn{L} (default
#'   \code{64L}; integer \eqn{\ge 16}), \code{n_components} number of frozen FPCA
#'   components \eqn{K} (default \code{6L}; integer \eqn{\ge 1}; the realised
#'   \eqn{K} is reduced at fit time when the training centred QFs span a subspace
#'   of rank \eqn{< K}), \code{shrink} ridge fraction of the mean covariance
#'   diagonal (default \code{0.1}; \eqn{\ge 0}), \code{min_features} per-specimen
#'   floor on the number of present rCLR coordinates (default \code{8L}; integer
#'   \eqn{\ge 3}), \code{eps} positive standard-deviation / ridge floor (default
#'   \code{1e-6}), and \code{seed} (default \code{1L}; accepted for interface
#'   uniformity only -- the method consumes no randomness).
#'
#' @return A plain list of class \code{inv_fdaqf_model} containing
#'   \code{feature_universe}, the frozen \code{fpca} basis (functional mean
#'   \code{qbar}, orthonormal loadings \code{V}, realised \code{K}, numerical
#'   \code{rank}, training \code{singular_values}), \code{descriptor_dim}
#'   (\eqn{= K}), the frozen \code{head} (\code{center}, \code{scale},
#'   \code{active}, \code{w}, \code{b}, \code{ridge}), and the resolved \code{hp}.
#'
#' @details
#' The per-specimen score (see \code{\link{score_inv_fdaqf}}) is the frozen
#' linear predictor on the standardised FPCA score vector,
#' \deqn{S(x) = b + \sum_{c} w_c \, \frac{s_c(x) - \mathrm{center}_c}
#'                                        {\mathrm{scale}_c},
#'       \qquad s(x) = V^\top (q(x) - \bar q),}
#' where \eqn{q(x)} is the specimen's quantile function. Every quantity except
#' \eqn{q(x)} is frozen at fit time, and \eqn{q(x)} depends only on \eqn{x}'s own
#' rCLR profile, so the score is computable from a single specimen and is exactly
#' invariant to per-sample positive scaling. Although it reads the same sorted
#' rCLR substrate as the curvature-scale-space and optimal-transport profile
#' methods, the representation here is a distinct functional-PCA projection.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 60
#' features <- paste0("miR-", sprintf("%03d", seq_len(p)))
#' n <- 160
#' y <- rep(c(0, 1), each = n / 2)
#' X <- matrix(0, n, p, dimnames = list(NULL, features))
#' rank_frac <- (seq_len(p) - 1) / (p - 1)
#' ctrl_mu <- 2 + 3 * rank_frac                               # near-linear QF
#' case_mu <- 2 + 3 * rank_frac^3                             # heavier upper tail
#' for (i in seq_len(n)) {
#'   mu <- if (y[i] == 1) case_mu else ctrl_mu
#'   X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
#' }
#' model <- fit_inv_fdaqf(X, y)
#' score <- score_inv_fdaqf(model, X)
#' }
#'
#' @references
#' Ramsay JO, Silverman BW. (2005) Functional Data Analysis, 2nd ed.
#' Springer. (Functional principal-component analysis.)
#'
#' Bolstad BM, Irizarry RA, Astrand M, Speed TP. (2003) A comparison of
#' normalization methods for high density oligonucleotide array data based on
#' variance and bias. \emph{Bioinformatics} 19(2): 185-193. (Quantile
#' representation of a sample distribution.)
#'
#' Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse compositional
#' technique reveals microbial perturbations. \emph{mSystems} 4(1): e00016-19.
#' (Robust CLR.)
#'
#' @export
fit_inv_fdaqf <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_inv_fdaqf", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_inv_fdaqf", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_inv_fdaqf")
  hp <- .inv_fdaqf_resolve_hp(hp)

  if (ncol(X_train) < hp$min_features) {
    stop("fit_inv_fdaqf: X_train must contain at least hp$min_features features")
  }

  # Per-specimen quantile function over the FULL training universe, built through
  # the SAME helper used at scoring so the QF is defined once (no fit/score
  # drift). A specimen with fewer than min_features present rCLR coordinates
  # yields NULL and is dropped from the basis/head fit.
  qf_list <- lapply(seq_len(nrow(X_train)), function(i) {
    .inv_fdaqf_qf(X_train[i, ], hp$resample_len, hp$min_features)
  })
  keep <- !vapply(qf_list, is.null, logical(1))
  if (sum(keep) < 2L) {
    stop("fit_inv_fdaqf: too few training specimens have >= hp$min_features ",
         "present rCLR coordinates")
  }
  G <- do.call(rbind, qf_list[keep])
  yk <- y[keep]
  if (sum(yk == 1L) < 1L || sum(yk == 0L) < 1L) {
    stop("fit_inv_fdaqf: after the min_features filter both classes must retain ",
         ">= 1 specimen")
  }

  # Frozen functional-PCA basis (training only): functional mean + top-K right
  # singular vectors of the centred QF matrix, capped at the numerical rank.
  fpca <- .inv_fdaqf_fit_basis(G, hp$n_components, hp$eps)

  # Per-specimen FPCA score vectors via the SAME projection used at scoring (one
  # definition -> no fit/score drift). Equivalent to the SVD scores Gc %*% V.
  if (fpca$K == 0L) {
    Phi <- matrix(numeric(0), nrow = nrow(G), ncol = 0L)
  } else {
    Phi <- do.call(rbind, lapply(seq_len(nrow(G)), function(i) {
      .inv_fdaqf_project(G[i, ], fpca$qbar, fpca$V)
    }))
  }

  head <- .inv_fdaqf_fit_head(Phi, yk, hp$shrink, hp$eps)

  model <- list(
    feature_universe = colnames(X_train),
    fpca = fpca,
    descriptor_dim = fpca$K,
    head = head,
    hp = hp
  )
  class(model) <- "inv_fdaqf_model"
  model
}


#' @title Score quantile-function FDA (inv-fdaqf) within-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_inv_fdaqf}}. The model-universe features present in \code{X}
#' (\code{intersect(model$feature_universe, colnames(X))}) define the present
#' substrate; for one specimen its within-sample rCLR profile over exactly those
#' present features is resampled to the frozen quantile-function length, projected
#' onto the frozen FPCA eigenbasis (same functional mean and loadings as at
#' fitting), standardised by the frozen \code{center}/\code{scale} and passed
#' through the frozen linear head. Larger values are more case-like.
#'
#' Scoring uses only each row's own values plus the frozen model -- no test-batch
#' renormalisation, no cross-row coupling, no statistic estimated from \code{X}.
#' If fewer than \code{model$hp$min_features} model-universe features are present,
#' the documented neutral score \code{0} is returned for every row; a single
#' specimen whose own number of present rCLR coordinates falls below
#' \code{model$hp$min_features} likewise scores the neutral \code{0}. A degenerate
#' fit whose frozen basis has realised \eqn{K=0} (every kept training QF
#' identical) scores every row the constant \code{0}.
#'
#' @param model An \code{inv_fdaqf_model} object returned by
#'   \code{\link{fit_inv_fdaqf}}.
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
#' rank_frac <- (seq_len(p) - 1) / (p - 1)
#' ctrl_mu <- 2 + 3 * rank_frac
#' case_mu <- 2 + 3 * rank_frac^3
#' for (i in seq_len(n)) {
#'   mu <- if (y[i] == 1) case_mu else ctrl_mu
#'   X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
#' }
#' model <- fit_inv_fdaqf(X, y)
#' score_inv_fdaqf(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Ramsay JO, Silverman BW. (2005) Functional Data Analysis, 2nd ed. Springer.
#'
#' @export
score_inv_fdaqf <- function(model, X, meta = NULL) {
  if (!inherits(model, "inv_fdaqf_model")) {
    stop("score_inv_fdaqf: model must have class inv_fdaqf_model")
  }
  X <- .reo_check_matrix(X, "score_inv_fdaqf", "X")
  .reo_check_meta(meta, nrow(X), "score_inv_fdaqf", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, present, drop = FALSE]
  head <- model$head
  qbar <- model$fpca$qbar
  V <- model$fpca$V
  L <- model$hp$resample_len
  mf <- model$hp$min_features
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    q <- .inv_fdaqf_qf(X_use[i, ], L, mf)
    if (is.null(q)) return(0)
    .inv_fdaqf_head_predict(.inv_fdaqf_project(q, qbar, V), head)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_inv_fdaqf: scorer produced non-finite or wrong-length output")
  }
  out
}


# ---------------------------------------------------------------------------
# Hyperparameter resolver (strict allowed-list, exact `[[ ]]` reads, per-field
# validation with integer-overflow guards). No randomness is consumed. The
# realised number of FPCA components may be reduced below n_components at fit
# time when the training centred QFs span a lower-rank subspace; the resolver
# only validates the requested n_components >= 1.
# ---------------------------------------------------------------------------
.inv_fdaqf_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_inv_fdaqf: hp must be a list")
  # Reject fully-unnamed / partially-unnamed / duplicated hp before the
  # allowed-list check: setdiff(names(hp), allowed) silently passes a NULL-names
  # list (e.g. list(64L)) and keeps only the last of a duplicated name.
  if (length(hp) > 0L) {
    nms <- names(hp)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("fit_inv_fdaqf: all hp entries must be named")
    }
    if (anyDuplicated(nms)) {
      stop("fit_inv_fdaqf: duplicated hp field(s): ",
           paste(unique(nms[duplicated(nms)]), collapse = ", "))
    }
  }
  allowed <- c("resample_len", "n_components", "shrink", "min_features", "eps",
               "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_inv_fdaqf: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  resample_len <- hp[["resample_len"]]
  if (is.null(resample_len)) resample_len <- 64L
  if (!is.numeric(resample_len) || length(resample_len) != 1L ||
      !is.finite(resample_len) || resample_len > .Machine$integer.max ||
      resample_len != as.integer(resample_len)) {
    stop("fit_inv_fdaqf: hp$resample_len must be an integer >= 16")
  }
  resample_len <- as.integer(resample_len)
  if (resample_len < 16L) {
    stop("fit_inv_fdaqf: hp$resample_len must be an integer >= 16")
  }

  n_components <- hp[["n_components"]]
  if (is.null(n_components)) n_components <- 6L
  if (!is.numeric(n_components) || length(n_components) != 1L ||
      !is.finite(n_components) || n_components > .Machine$integer.max ||
      n_components != as.integer(n_components)) {
    stop("fit_inv_fdaqf: hp$n_components must be an integer >= 1")
  }
  n_components <- as.integer(n_components)
  if (n_components < 1L) {
    stop("fit_inv_fdaqf: hp$n_components must be an integer >= 1")
  }

  shrink <- hp[["shrink"]]
  if (is.null(shrink)) shrink <- 0.1
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0) {
    stop("fit_inv_fdaqf: hp$shrink must be a non-negative finite number")
  }
  shrink <- as.numeric(shrink)

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 8L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_inv_fdaqf: hp$min_features must be an integer >= 3")
  }
  min_features <- as.integer(min_features)
  if (min_features < 3L) {
    stop("fit_inv_fdaqf: hp$min_features must be an integer >= 3")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_inv_fdaqf: hp$eps must be a positive finite number")
  }
  eps <- as.numeric(eps)

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_inv_fdaqf: hp$seed must be a non-negative integer")
  }
  seed <- as.integer(seed)

  list(
    resample_len = resample_len,
    n_components = n_components,
    shrink = shrink,
    min_features = min_features,
    eps = eps,
    seed = seed
  )
}


# ---------------------------------------------------------------------------
# Single-specimen within-sample robust-CLR over the specimen's own positive
# support, returning ALL present (v > 0) coordinates. SELF-CONTAINED (does not
# reuse any package *_rclr_matrix / *_rclr helper, several of which centre using
# cross-sample reference statistics and would break single-sample
# deployability). Only STRUCTURAL zeros (absent parts, v == 0) are excluded --
# NOT coordinates whose centred-log value happens to equal the geometric mean
# (z == 0; e.g. a tied / constant count profile), which are kept. Keying on
# `z != 0` here would wrongly drop those and would empty a constant positive
# profile, distorting the resampled curve. Exactly scale-invariant:
# log(c*v) - mean(log(c*v)) = log(v) - mean(log(v)).
# ---------------------------------------------------------------------------
.inv_fdaqf_rclr_pos <- function(v) {
  pos <- which(v > 0)
  if (length(pos) == 0L) return(numeric(0))
  lv <- log(v[pos])
  lv - mean(lv)
}


# ---------------------------------------------------------------------------
# Single-specimen quantile function: the present rCLR coordinates resampled to a
# fixed length L by evaluating their empirical quantile function at the midpoint
# grid u_l = (l - 0.5)/L (type-7). stats::quantile sorts internally, so the QF
# is independent of feature column order and is a monotone non-decreasing curve.
# Returns NULL when the specimen has fewer than min_features present rCLR
# coordinates (the documented neutral case). Used at fit (full universe) and at
# score (present subset) -- one definition, no fit/score drift.
# ---------------------------------------------------------------------------
.inv_fdaqf_qf <- function(v, L, min_features) {
  zp <- .inv_fdaqf_rclr_pos(v)
  if (length(zp) < min_features) return(NULL)
  u <- (seq_len(L) - 0.5) / L
  stats::quantile(zp, probs = u, type = 7, names = FALSE)
}


# ---------------------------------------------------------------------------
# Frozen functional-PCA basis over the training QF matrix G (n_kept x L). The
# functional mean qbar centres G; the top-K right singular vectors of the
# centred matrix Gc are the orthonormal FPCA loadings V (L x K). K is capped at
# the numerical rank of Gc -- the count of singular values above a relative
# tolerance (a robust reading of "positive singular values"; a bare d > 0 test
# would count floating-point noise as genuine components). When Gc is exactly
# rank 0 (all kept training QFs identical) the realised K is 0 and V is an
# L x 0 matrix; the downstream head then degenerates to a constant-0 predictor.
# svd() is deterministic and consumes no RNG. The frozen V column signs /
# in-subspace rotation are LAPACK-build dependent but absorbed by the linear
# head, so final scores are sign-robust.
# ---------------------------------------------------------------------------
.inv_fdaqf_fit_basis <- function(G, n_components, eps) {
  qbar <- colMeans(G)
  Gc <- sweep(G, 2L, qbar, "-")
  sv <- svd(Gc)
  d <- sv$d
  sv_tol <- max(nrow(Gc), ncol(Gc)) * .Machine$double.eps * max(d, 0)
  rank <- sum(d > sv_tol)
  K <- min(as.integer(n_components), rank)
  V <- sv$v[, seq_len(K), drop = FALSE]   # L x K, orthonormal columns (frozen)
  list(
    qbar = qbar,
    V = V,
    K = as.integer(K),
    rank = as.integer(rank),
    singular_values = d
  )
}


# Project a single specimen's quantile function onto the frozen FPCA basis:
# s = V^T (q - qbar), length ncol(V) = K (numeric(0) when K == 0). Uses only the
# specimen's own QF and the frozen qbar, V -- no scored-batch statistic.
.inv_fdaqf_project <- function(q, qbar, V) {
  as.numeric(crossprod(V, q - qbar))
}


# ---------------------------------------------------------------------------
# Frozen standardize + ridge-LDA head over the training FPCA score vectors.
# center/scale are the per-component training mean / sd (sd floored at eps);
# exactly-zero-variance components are deactivated (forced to contribute 0). The
# pooled within-class covariance is ridged until its Cholesky succeeds, then
# w = Sinv (mu_case - mu_ctrl), b = -0.5 (mu_case + mu_ctrl) . w. Larger = case.
# A zero-dimensional descriptor (K == 0, degenerate rank-0 basis) yields a
# well-formed constant-0 head.
# ---------------------------------------------------------------------------
.inv_fdaqf_fit_head <- function(Phi, y, shrink, eps) {
  d <- ncol(Phi)
  if (d == 0L) {
    # No projection coordinates: the linear predictor is the constant 0 for
    # every specimen (.inv_fdaqf_head_predict on a length-0 descriptor -> b = 0).
    return(list(center = numeric(0), scale = numeric(0), active = logical(0),
                w = numeric(0), b = 0, ridge = eps))
  }
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
    stop("fit_inv_fdaqf: within-class covariance could not be made positive ",
         "definite")
  }
  Sinv <- chol2inv(ch)

  delta <- mu_case - mu_ctrl
  w <- as.numeric(Sinv %*% delta)
  w[!active] <- 0                         # inactive components carry no weight
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

# Apply the frozen head to one specimen's FPCA score vector: standardise with
# the frozen center/scale, zero inactive components, return the linear
# predictor. A length-0 descriptor (degenerate K == 0 head) returns b = 0.
.inv_fdaqf_head_predict <- function(phi, head) {
  z <- (phi - head$center) / head$scale
  z[!head$active] <- 0
  head$b + sum(head$w * z)
}
