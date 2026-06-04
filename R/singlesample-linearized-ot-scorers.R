#' @title Fit a linearized-optimal-transport (CDT tangent) LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data by embedding each
#' specimen into \strong{linearized-optimal-transport (LOT) tangent space} -- in
#' one dimension the exact \strong{Cumulative Distribution Transform (CDT)} of
#' Park-Kolouri-Rohde -- and then separating the two classes with a linear
#' (LDA-style, pooled-covariance) Gaussian discriminant in that tangent space.
#' Each specimen is first mapped to a self-contained, per-sample robust
#' centred-log-ratio (rCLR) representation \eqn{z}, where
#' \eqn{z_j = \log v_j - \mathrm{mean}_{k:\,v_k>0}\log v_k} on its nonzero parts
#' and \eqn{z_j = 0} on its zero parts. Because the centring uses only that
#' specimen's own nonzero values, \eqn{z(c\,v) = z(v)} for any \eqn{c>0}: the
#' representation -- and therefore every downstream score -- is exactly (to
#' numerical tolerance) invariant to per-sample scaling (library size / input
#' amount).
#'
#' @section The single-specimen measure (why a 1-D CDT and not a literal W2):
#' A literal squared 2-Wasserstein distance from a \emph{single} specimen (a
#' Dirac point measure) to a class distribution is
#' \eqn{\lVert\cdot\rVert^2 + \text{const}} -- a covariance-blind nearest-centroid
#' that discards the class shape, the same collapse that forces a Gaussian LRT in
#' \code{\link{fit_lrt_bw}}. We avoid it by treating each specimen as a
#' \strong{one-dimensional empirical measure over its own present rCLR coordinate
#' values}. Take the present rCLR coordinates \eqn{z_{\mathrm{pos}} = \{z_j: v_j >
#' 0\}} (rCLR sets absent parts to exactly \eqn{0}; including the structural-zero
#' point mass would leak detection rate / sequencing depth -- a batch artefact --
#' so only structural zeros, keyed on the raw abundance \eqn{v_j = 0}, are
#' excluded; a present coordinate that happens to equal the geometric mean
#' (\eqn{z_j = 0}) is kept). \eqn{z_{\mathrm{pos}}} is
#' mean-zero by rCLR construction, so its \emph{location} is fixed at \eqn{0} and
#' the discriminative content is the \strong{distributional shape} (spread, skew,
#' tails) of the centred-log-ratio profile. The specimen's measure is
#' \eqn{\nu = \tfrac{1}{m}\sum_{j}\delta_{z_{\mathrm{pos},j}}}
#' (\eqn{m = |z_{\mathrm{pos}}|}); its quantile function \eqn{Q_\nu} evaluated at
#' fixed interior levels by the empirical inverse-CDF (\code{type = 1}) is the
#' exact 1-D optimal-transport (monotone rearrangement) representation. This
#' representation is
#' \strong{feature-exchangeable by design} -- it scores the shape of the rCLR
#' value-distribution -- a mechanism distinct from and complementary to the
#' feature-identity-preserving sliced / Gaussian LRTs.
#'
#' @section Fit (training only; no RNG consumed):
#' \enumerate{
#'   \item A fixed interior quantile grid \eqn{t_k = (k-\tfrac12)/K},
#'     \eqn{k = 1,\dots,K} (\eqn{K = } \code{hp$n_quantiles}).
#'   \item For each training specimen the per-sample rCLR is formed over the
#'     feature universe, its present coordinates \eqn{z_{\mathrm{pos}}} are taken,
#'     and -- if at least \code{hp$min_features} of them exist -- its CDT quantile
#'     function \eqn{Q = \mathrm{quantile}(z_{\mathrm{pos}}, t_k, \mathrm{type}=1)}
#'     (the empirical inverse-CDF, i.e. the exact 1-D OT map)
#'     \eqn{\in \mathbb{R}^K} is recorded (a specimen with fewer usable
#'     coordinates is skipped; at least two usable specimens per class are
#'     required).
#'   \item The \strong{frozen reference quantile function} \eqn{Q_{\mathrm{ref}}}
#'     (the LOT tangent base point) is the 1-D Wasserstein barycenter
#'     \eqn{\mathrm{colMeans}(Q_i)} over all usable specimens
#'     (\code{hp$reference = "barycenter"}; the 1-D barycenter quantile function is
#'     the mean of the quantile functions) or the pooled quantile function
#'     \eqn{\mathrm{quantile}(\bigcup_i z_{\mathrm{pos},i}, t_k)}
#'     (\code{hp$reference = "pooled"}).
#'   \item The \strong{CDT tangent embedding} is \eqn{v_i = Q_i - Q_{\mathrm{ref}}}
#'     (the Wasserstein log-map at the reference).
#'   \item A \strong{linear discriminant head} (LDA-style, which keeps the method
#'     faithful to \emph{linear} OT): per-class mean \eqn{\mu_c}, pooled
#'     within-class covariance
#'     \eqn{S = ((n_1-1)\,\mathrm{Cov}_1 + (n_0-1)\,\mathrm{Cov}_0)/(n_1+n_0-2)}.
#'     Because the quantile vector is monotone, adjacent components are highly
#'     correlated and \eqn{S} is near-singular, so regularization is required: a
#'     fixed diagonal load \code{hp$shrink}\eqn{\times\mathrm{mean}(\mathrm{diag}
#'     S)} is added first, then a ridge is raised from \code{hp$eps} (doubling)
#'     until \code{chol} succeeds. With \eqn{S_{\mathrm{inv}} =
#'     \mathrm{chol2inv}(\mathrm{chol}(S+\lambda I))} the head is
#'     \eqn{w = S_{\mathrm{inv}}(\mu_1-\mu_0)},
#'     \eqn{b = -\tfrac12(\mu_1+\mu_0)^\top w}.
#' }
#'
#' @section Reference-origin invariance:
#' The decision value \eqn{w^\top(Q-Q_{\mathrm{ref}})+b} is algebraically
#' independent of the choice of \eqn{Q_{\mathrm{ref}}}: \eqn{w} depends only on the
#' class-mean difference (a constant shift cancels) and the
#' \eqn{Q_{\mathrm{ref}}} terms in \eqn{w^\top(Q-Q_{\mathrm{ref}})} and \eqn{b}
#' cancel, leaving the textbook LDA decision
#' \eqn{w^\top Q - \tfrac12 w^\top(\bar Q_1+\bar Q_0)}. The reference fixes the
#' tangent origin (and is reported for inspection) but does not change any score.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{n_quantiles} CDT grid size
#'   \eqn{K} (positive integer \eqn{\ge 4}, default \code{16L}); \code{reference}
#'   the tangent base point, one of \code{"barycenter"} or \code{"pooled"}
#'   (default \code{"barycenter"}); \code{shrink} fixed diagonal-load fraction
#'   added as \code{shrink * mean(diag(S))} before the ridge loop (non-negative
#'   finite, default \code{0.1}); \code{min_features} minimum present rCLR
#'   coordinates needed to form a quantile function (positive integer \eqn{\ge 2},
#'   default \code{5L}); \code{eps} positive ridge increment / seed (default
#'   \code{1e-6}); and \code{seed} non-negative integer kept for API consistency
#'   only -- no RNG is consumed, and a fit leaves \code{.Random.seed} untouched
#'   (default \code{1L}).
#'
#' @return A plain list of class \code{ot_lot_model} containing \code{t_k} (the
#'   frozen quantile grid), \code{Q_ref} (the frozen reference quantile function),
#'   \code{head} (the frozen linear discriminant: \code{w}, \code{b}, per-class
#'   tangent means \code{mu_case}/\code{mu_control}, \code{Sinv}, the resolved
#'   ridge, and the per-class usable counts), \code{feature_universe},
#'   \code{reference}, the per-class usable counts, and the resolved \code{hp}.
#'
#' @details
#' The single-specimen score (see \code{\link{score_ot_lot}}) is the frozen linear
#' functional \eqn{S(Q) = w^\top(Q - Q_{\mathrm{ref}}) + b} of the specimen's own
#' CDT embedding, larger meaning more case-like. The reference \eqn{Q_{\mathrm{ref}}},
#' the grid \eqn{t_k} and the head \eqn{(w,b)} are all frozen from training; the
#' specimen's quantile function is computed from its own present coordinates, so
#' no scored-batch statistic and no cross-row coupling enter. Because the whole
#' model lives in the fixed-dimension CDT tangent space \eqn{\mathbb{R}^K} (no
#' feature-indexed component beyond which features are admitted), partial feature
#' overlap introduces no systematic batch-statistic drift -- there is no
#' cross-sample renormalisation, and the same specimen scored alone always yields
#' the same score (the inclusion-gate property). Fewer present features do,
#' however, lower the resolution of (and add sampling variance to) the empirical
#' quantile function, which is still evaluated at the frozen \eqn{t_k}; this added
#' variance is intrinsic to the value-distribution representation, not a leak of
#' the present-feature count.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200; p <- 40
#' base <- 6
#' y <- rep(c(0, 1), each = n / 2)
#' E <- matrix(0, n, p, dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' E[y == 0, ] <- stats::rnorm(sum(y == 0) * p, 0, 0.30)   # tight controls
#' E[y == 1, ] <- stats::rnorm(sum(y == 1) * p, 0, 0.60)   # dispersed cases
#' X <- exp(base + E)
#' model <- fit_ot_lot(X, y)
#' score <- score_ot_lot(model, X)
#' }
#'
#' @references
#' Park S.R., Kolouri S., Kundu S., Rohde G.K. (2018) The cumulative distribution
#' transform and linear pattern classification. \emph{Applied and Computational
#' Harmonic Analysis} 45(3): 616-641.
#'
#' Kolouri S., Park S.R., Thorpe M., Slepcev D., Rohde G.K. (2017) Optimal mass
#' transport: signal processing and machine-learning applications. \emph{IEEE
#' Signal Processing Magazine} 34(4): 43-59.
#'
#' Wang W., Slepcev D., Basu S., Ozolek J.A., Rohde G.K. (2013) A linear optimal
#' transportation framework for quantifying and visualizing variation in sets of
#' images. \emph{International Journal of Computer Vision} 101(2): 254-269.
#'
#' Fisher R.A. (1936) The use of multiple measurements in taxonomic problems.
#' \emph{Annals of Eugenics} 7(2): 179-188.
#'
#' @export
fit_ot_lot <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ot_lot", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ot_lot", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ot_lot")
  hp <- .ot_lot_resolve_hp(hp)

  t_k <- .ot_lot_grid(hp$n_quantiles)

  # Per-specimen present rCLR coordinates (the 1-D empirical measure). The same
  # transform is reused at scoring, so the CDT geometry is defined once and there
  # is no fit/score drift.
  zp_list <- lapply(seq_len(nrow(X_train)), function(i) {
    .ot_lot_rclr_pos(X_train[i, ])
  })
  usable <- vapply(zp_list, function(zp) length(zp) >= hp$min_features,
                   logical(1))

  n_usable_case <- sum(usable & y == 1L)
  n_usable_control <- sum(usable & y == 0L)
  if (n_usable_case < 2L || n_usable_control < 2L) {
    stop(paste0("fit_ot_lot: needs at least 2 usable specimens per class (with ",
                ">= hp$min_features present rCLR coordinates) to estimate a ",
                "tangent-space covariance"))
  }

  use_idx <- which(usable)
  Q_all <- do.call(rbind, lapply(use_idx, function(i) {
    .ot_lot_cdt(zp_list[[i]], t_k)
  }))
  y_use <- y[use_idx]

  # Frozen reference quantile function (the LOT tangent base point).
  if (identical(hp$reference, "barycenter")) {
    Q_ref <- colMeans(Q_all)                       # 1-D Wasserstein barycenter
  } else {
    pooled <- unlist(zp_list[use_idx], use.names = FALSE)
    Q_ref <- stats::quantile(pooled, probs = t_k, type = 7, names = FALSE)
  }
  Q_ref <- as.numeric(Q_ref)

  # CDT tangent embedding and the LDA-style linear head.
  V <- sweep(Q_all, 2L, Q_ref, "-")
  V_case <- V[y_use == 1L, , drop = FALSE]
  V_control <- V[y_use == 0L, , drop = FALSE]
  mu_case <- colMeans(V_case)
  mu_control <- colMeans(V_control)
  n_case <- nrow(V_case)
  n_control <- nrow(V_control)
  cov_case <- stats::cov(V_case)
  cov_control <- stats::cov(V_control)
  S <- ((n_case - 1) * cov_case + (n_control - 1) * cov_control) /
    (n_case + n_control - 2)
  S <- (S + t(S)) / 2

  # Monotone-quantile correlation => S near-singular: a fixed diagonal load is
  # added before the eps ridge loop (both documented in hp).
  diag_load <- hp$shrink * mean(diag(S))
  if (is.finite(diag_load) && diag_load > 0) {
    S <- S + diag_load * diag(hp$n_quantiles)
  }
  pd <- .ot_lot_chol_pd(S, hp$eps)

  delta <- mu_case - mu_control
  w <- as.numeric(pd$Sinv %*% delta)
  b <- as.numeric(-0.5 * sum((mu_case + mu_control) * w))

  head <- list(
    w = w,
    b = b,
    mu_case = mu_case,
    mu_control = mu_control,
    Sinv = pd$Sinv,
    ridge = pd$ridge,
    diag_load = diag_load,
    n_case = n_case,
    n_control = n_control
  )

  model <- list(
    t_k = t_k,
    Q_ref = Q_ref,
    head = head,
    feature_universe = colnames(X_train),
    reference = hp$reference,
    n_usable_case = n_usable_case,
    n_usable_control = n_usable_control,
    hp = hp
  )
  class(model) <- "ot_lot_model"
  model
}


#' @title Score a linearized-optimal-transport (CDT tangent) LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_ot_lot}}. For one specimen the abundances are mapped to the
#' self-contained per-sample rCLR representation, its present coordinates
#' \eqn{z_{\mathrm{pos}}} are taken, and its CDT quantile function
#' \eqn{Q = \mathrm{quantile}(z_{\mathrm{pos}}, t_k, \mathrm{type}=1)} (the
#' empirical inverse-CDF) is formed on
#' the frozen grid. The score is the frozen linear functional of the CDT tangent
#' embedding
#' \deqn{S(Q) = w^\top (Q - Q_{\mathrm{ref}}) + b,}
#' larger meaning more case-like.
#'
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}; the per-sample rCLR is
#' formed over exactly that set. The reference \eqn{Q_{\mathrm{ref}}}, the grid
#' \eqn{t_k} and the head \eqn{(w,b)} are frozen from training and live in the
#' fixed-dimension CDT tangent space, so no present-subset re-derivation is needed
#' and no scored-batch statistic is used. The score of a row therefore depends
#' only on that row and the frozen model, and is exactly (to numerical tolerance)
#' invariant to per-sample scaling. If fewer than \code{model$hp$min_features}
#' feature-universe columns are present, the documented neutral score \code{0} is
#' returned for every row; a specimen whose own present rCLR coordinates number
#' fewer than \code{model$hp$min_features} (e.g. a near-empty profile) also scores
#' the neutral \code{0}.
#'
#' @param model An \code{ot_lot_model} object returned by \code{\link{fit_ot_lot}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity and
#'   ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values are
#'   more case-like. Scoring uses only each row's own values plus the frozen
#'   reference, grid and head.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200; p <- 40
#' y <- rep(c(0, 1), each = n / 2)
#' E <- matrix(0, n, p, dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' E[y == 0, ] <- stats::rnorm(sum(y == 0) * p, 0, 0.30)
#' E[y == 1, ] <- stats::rnorm(sum(y == 1) * p, 0, 0.60)
#' X <- exp(6 + E)
#' model <- fit_ot_lot(X, y)
#' score_ot_lot(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Park S.R., Kolouri S., Kundu S., Rohde G.K. (2018) The cumulative distribution
#' transform and linear pattern classification. \emph{Applied and Computational
#' Harmonic Analysis} 45(3): 616-641.
#'
#' @export
score_ot_lot <- function(model, X, meta = NULL) {
  if (!inherits(model, "ot_lot_model")) {
    stop("score_ot_lot: model must have class ot_lot_model")
  }
  X <- .reo_check_matrix(X, "score_ot_lot", "X")
  .reo_check_meta(meta, nrow(X), "score_ot_lot", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, present, drop = FALSE]
  t_k <- model$t_k
  Q_ref <- model$Q_ref
  w <- model$head$w
  b <- model$head$b
  min_features <- model$hp$min_features

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    zp <- .ot_lot_rclr_pos(X_use[i, ])
    if (length(zp) < min_features) {
      return(0)
    }
    Q <- .ot_lot_cdt(zp, t_k)
    sum(w * (Q - Q_ref)) + b
  }, numeric(1))

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ot_lot: scorer produced non-finite or wrong-length output")
  }
  out
}

.ot_lot_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ot_lot: hp must be a list")
  # Reject fully-unnamed / partially-unnamed / duplicated hp before the
  # allowed-list check: setdiff(names(hp), allowed) silently passes a NULL-names
  # list (e.g. list(24L)) and keeps only the last of a duplicated name.
  if (length(hp) > 0L) {
    nms <- names(hp)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("fit_ot_lot: all hp entries must be named")
    }
    if (anyDuplicated(nms)) {
      stop("fit_ot_lot: duplicated hp field(s): ",
           paste(unique(nms[duplicated(nms)]), collapse = ", "))
    }
  }
  allowed <- c("n_quantiles", "reference", "shrink", "min_features", "eps",
               "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ot_lot: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  # EXACT [[ ]] reads only: list `$` does partial matching, which would let a
  # prefix (e.g. a stray field) silently bind to a longer allowed name.
  n_quantiles <- hp[["n_quantiles"]]
  if (is.null(n_quantiles)) n_quantiles <- 16L
  if (!is.numeric(n_quantiles) || length(n_quantiles) != 1L ||
      !is.finite(n_quantiles) || n_quantiles < 4L ||
      n_quantiles > .Machine$integer.max ||
      n_quantiles != as.integer(n_quantiles)) {
    stop("fit_ot_lot: hp$n_quantiles must be a positive integer >= 4")
  }

  reference <- hp[["reference"]]
  if (is.null(reference)) reference <- "barycenter"
  if (!is.character(reference) || length(reference) != 1L || is.na(reference) ||
      !reference %in% c("barycenter", "pooled")) {
    stop("fit_ot_lot: hp$reference must be one of 'barycenter' or 'pooled'")
  }

  shrink <- hp[["shrink"]]
  if (is.null(shrink)) shrink <- 0.1
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0) {
    stop("fit_ot_lot: hp$shrink must be a single non-negative finite number")
  }

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 5L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 2L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_ot_lot: hp$min_features must be a positive integer >= 2")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_ot_lot: hp$eps must be a positive finite number")
  }

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_ot_lot: hp$seed must be a non-negative integer")
  }

  list(
    n_quantiles = as.integer(n_quantiles),
    reference = as.character(reference),
    shrink = as.numeric(shrink),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

# Fixed interior quantile grid t_k = (k - 0.5) / K, k = 1..K.
.ot_lot_grid <- function(K) {
  (seq_len(K) - 0.5) / K
}

# Self-contained per-sample robust CLR of one specimen. The geometric-mean
# centring uses ONLY this specimen's own nonzero parts, so it is a within-sample
# transform: z(c*v) = z(v) for any c > 0 (exact per-sample scale-invariance), and
# zero parts map to 0. This deliberately does NOT reuse any package rclr helper
# that centres using cross-sample / reference statistics.
.ot_lot_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Present (positive-abundance) rCLR coordinates of one specimen -- the support of
# its 1-D empirical measure. Only STRUCTURAL zeros (absent parts, v == 0) are
# excluded (keying on `v > 0`, NOT `z != 0`): dropping the structural-zero point
# mass avoids leaking detection rate / sequencing depth, while a genuinely present
# coordinate whose centred-log value happens to equal the geometric mean (z == 0;
# e.g. a tied / constant count profile) is correctly KEPT. A fully constant
# positive profile therefore yields a degenerate point mass at 0 of size =
# #positive features (scored, not silently dropped), rather than an empty vector.
.ot_lot_rclr_pos <- function(v) {
  z <- .ot_lot_rclr(v)
  z[v > 0]
}

# Exact 1-D CDT quantile function of a present-rCLR support at the frozen grid.
# In one dimension the optimal-transport map to any reference is the monotone
# rearrangement, so the quantile function (which sorts internally) IS the LOT
# embedding coordinate. type = 1 (the empirical inverse-CDF / step quantile) is
# used deliberately: it is the EXACT quantile function of the equal-weight
# empirical measure. Evaluated at K = m interior levels, the squared L2 between
# two equal-m CDT embeddings equals the exact 1-D squared 2-Wasserstein distance
# (= mean of the squared differences of the sorted supports) -- the LOT isometry.
# At the operational K != m, the K-vector L2 is a K-point quadrature of that W2^2,
# NOT exact equality; the head is a linear discriminant on the tangent embedding
# and never computes an exact Wasserstein distance (so the deployed K=16 score is
# not an exact-W2 quantity). type = 7 (piecewise-linear interpolation) would
# instead converge to a DIFFERENT, low-biased limit, breaking the OT identity.
.ot_lot_cdt <- function(zp, t_k) {
  stats::quantile(zp, probs = t_k, type = 1, names = FALSE)
}

# Cholesky of the (loaded) pooled tangent covariance with a ridge guard: try
# chol with no extra ridge first and, on failure, seed the ridge at `eps` and
# double it until chol succeeds. Returns the precision (chol2inv) and the ridge
# actually used. The monotone quantile vector makes the raw S near-singular, so
# the diagonal load + this ridge are what keep S + lambda I positive definite.
.ot_lot_chol_pd <- function(S, eps) {
  K <- nrow(S)
  Imat <- diag(K)
  S <- (S + t(S)) / 2
  ridge <- 0
  ch <- NULL
  for (k in 0:60) {
    Stry <- if (ridge == 0) S else S + ridge * Imat
    ch <- tryCatch(chol(Stry), error = function(e) NULL)
    if (!is.null(ch)) break
    ridge <- if (ridge == 0) eps else ridge * 2
  }
  if (is.null(ch)) {
    stop("fit_ot_lot: failed to regularize the pooled tangent covariance to positive definite")
  }
  list(Sinv = chol2inv(ch), ridge = ridge)
}
