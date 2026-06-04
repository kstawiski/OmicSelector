#' @title Fit a Bures-Wasserstein Gaussian-OT class-conditional LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data using a Gaussian
#' model per class whose covariance is regularized by optimal-transport
#' (Bures-Wasserstein, BW) geometry. Each specimen is first mapped to a
#' self-contained, per-sample robust centred-log-ratio (rCLR) representation
#' \eqn{z}, where \eqn{z_j = \log v_j - \mathrm{mean}_{k:\,v_k>0}\log v_k} on the
#' nonzero parts of the specimen and \eqn{z_j = 0} on its zero parts. Because the
#' centring uses only that specimen's own nonzero values, \eqn{z(c\,v) = z(v)} for
#' any \eqn{c>0}: the representation -- and therefore every downstream score -- is
#' exactly (to numerical tolerance) invariant to per-sample scaling (library size
#' / input amount).
#'
#' @section Why a Gaussian likelihood and not a literal Wasserstein distance:
#' The label \emph{bw} here means \strong{Bures-Wasserstein} -- an
#' optimal-transport metric between Gaussians / SPD covariance matrices -- and not
#' a bandwidth. A literal squared 2-Wasserstein distance from a \emph{single}
#' specimen (a Dirac point, or a fixed-covariance Gaussian) to a class Gaussian
#' \eqn{N(\mu_c,\Sigma_c)} is
#' \eqn{\lVert z - \mu_c\rVert^2 + (\text{class-constant Bures term})}: the class
#' covariance enters \emph{only} through a per-class additive constant, so the
#' W2-distance log-ratio collapses to a Euclidean nearest-centroid (linear)
#' classifier that does \emph{not} use \eqn{\Sigma_c} discriminatively. To build a
#' single-sample LRT that genuinely uses the frozen class covariances we score
#' with the Gaussian log-likelihood ratio (Neyman-Pearson optimal for Gaussian
#' class models); the \strong{Bures-Wasserstein element is carried by how the
#' covariances are estimated and regularized} -- by OT-geometric (BW-geodesic)
#' shrinkage of each class covariance toward a shared barycenter anchor. That is
#' what distinguishes this method from a plain Euclidean-shrinkage Mahalanobis
#' LRT.
#'
#' @section Fit (training only):
#' On the per-sample rCLR \eqn{z} of the frozen class anchors:
#' \enumerate{
#'   \item \eqn{\mu_c} is the class mean (the W2 / Frechet barycenter of the class
#'     points equals the ordinary Euclidean mean in rCLR space) and \eqn{\Sigma_c}
#'     is the class sample covariance (\code{stats::cov}), symmetrized.
#'   \item A shared anchor covariance \eqn{\Sigma_{\mathrm{pool}}} is either the
#'     equal-weight \strong{Bures-Wasserstein barycenter} of
#'     \eqn{\{\Sigma_{\mathrm{case}},\Sigma_{\mathrm{control}}\}} (the midpoint of
#'     the BW geodesic; \code{hp$anchor = "bw_barycenter"}) or the pooled
#'     within-class covariance (\code{hp$anchor = "pooled"}).
#'   \item Each class covariance is BW-geodesically shrunk toward the anchor by
#'     \eqn{t = } \code{hp$shrink} \eqn{\in [0,1)}:
#'     \eqn{\Sigma_c^{\mathrm{reg}} = \mathrm{BWgeo}(\Sigma_c,
#'     \Sigma_{\mathrm{pool}}, t)} (\eqn{t=0} keeps \eqn{\Sigma_c}; \eqn{t=1}
#'     gives the anchor), then a diagonal ridge is added until \code{chol}
#'     succeeds. The Cholesky factor yields \eqn{\Sigma_c^{\mathrm{reg}\,-1}} via
#'     \code{chol2inv} and \eqn{\log\det\Sigma_c^{\mathrm{reg}} = 2\sum\log
#'     \mathrm{diag}}.
#' }
#'
#' @section Bures-Wasserstein geodesic:
#' For SPD matrices \eqn{A,B} and \eqn{t\in[0,1]} the McCann / OT interpolation is
#' \deqn{T = A^{-1/2}\,(A^{1/2} B A^{1/2})^{1/2}\,A^{-1/2}, \quad
#'   M_t = (1-t) I + t\,T, \quad \mathrm{BWgeo}(A,B,t) = M_t\,A\,M_t,}
#' symmetrized. Symmetric matrix square roots use the eigendecomposition
#' \eqn{A^{1/2} = V\,\mathrm{diag}(\sqrt{\max(\lambda,0)})\,V^\top} and
#' \eqn{A^{-1/2} = V\,\mathrm{diag}(1/\sqrt{\max(\lambda,\varepsilon)})\,V^\top}
#' (tiny eigenvalues floored at \code{hp$eps} for the inverse square root). Both
#' \eqn{A} and \eqn{B} are floored to positive definite before forming the
#' geodesic. The equal-weight BW barycenter is \eqn{\mathrm{BWgeo}(A,B,0.5)}.
#'
#' @section High-dimension (p > n) behaviour:
#' When the number of features \eqn{p} exceeds a class's anchor count \eqn{n_c},
#' \eqn{\Sigma_c} is rank-deficient. The eigenvalue floor (null-space directions
#' receive variance \code{hp$eps}), the BW-geodesic shrinkage toward the
#' better-conditioned anchor (with a substantial default \code{hp$shrink = 0.5}),
#' and the final ridge-to-positive-definite loop together guarantee a positive
#' definite \eqn{\Sigma_c^{\mathrm{reg}}}, so the fit and every score remain
#' finite.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{shrink} BW-geodesic shrink of
#'   each class covariance toward the shared anchor, in \eqn{[0,1)} (default
#'   \code{0.5}); \code{anchor} the shared anchor covariance, one of
#'   \code{"bw_barycenter"} or \code{"pooled"} (default \code{"bw_barycenter"});
#'   \code{max_anchors_per_class} positive-integer anchor cap per class (default
#'   \code{200L}); \code{min_features} positive-integer feature-overlap floor at
#'   scoring (default \code{3L}); \code{eps} positive eigenvalue floor / ridge
#'   seed (default \code{1e-6}); and \code{seed} non-negative integer used only for
#'   deterministic anchor subsampling while restoring the global RNG state
#'   (default \code{1L}).
#'
#' @return A plain list of class \code{lrt_bw_model} containing
#'   \code{case_anchors}, \code{control_anchors}, \code{feature_universe},
#'   \code{repr} (the full-universe frozen class representation: per-class mean,
#'   regularized precision and log determinant), and the resolved \code{hp}.
#'
#' @details
#' The single-specimen score (see \code{\link{score_lrt_bw}}) is the Gaussian
#' log-likelihood ratio
#' \deqn{S(z) = \ell_{\mathrm{case}}(z) - \ell_{\mathrm{control}}(z), \qquad
#'   \ell_c(z) = -\tfrac{1}{2}(z-\mu_c)^\top \Sigma_c^{\mathrm{reg}\,-1}(z-\mu_c)
#'              - \tfrac{1}{2}\log\det\Sigma_c^{\mathrm{reg}},}
#' larger meaning more case-like. The shared \eqn{-\tfrac{p}{2}\log(2\pi)} term
#' cancels in the difference (both classes are scored over the same present
#' feature set) and is omitted.
#'
#' Fitting stores the raw case/control anchor rows (so the class representation
#' can be re-derived consistently over any present-feature subset at scoring), the
#' frozen training feature universe, the resolved hyperparameters, and -- for
#' inspection and tests -- the full-universe class representation. No test data and
#' no test-time batch statistic are used.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 160; p <- 30
#' L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' shift <- rep(c(1.2, -1.2), length.out = 12)
#' L[y == 1, 1:12] <- L[y == 1, 1:12] + matrix(shift, sum(y == 1), 12, byrow = TRUE)
#' X <- exp(L)
#' model <- fit_lrt_bw(X, y)
#' score <- score_lrt_bw(model, X)
#' }
#'
#' @references
#' Bhatia R., Jain T., Lim Y. (2019) On the Bures-Wasserstein distance between
#' positive definite matrices. \emph{Expositiones Mathematicae} 37(2): 165-191.
#'
#' Agueh M., Carlier G. (2011) Barycenters in the Wasserstein space.
#' \emph{SIAM Journal on Mathematical Analysis} 43(2): 904-924.
#'
#' Takatsu A. (2011) Wasserstein geometry of Gaussian measures.
#' \emph{Osaka Journal of Mathematics} 48(4): 1005-1026.
#'
#' Neyman J., Pearson E.S. (1933) On the problem of the most efficient tests of
#' statistical hypotheses. \emph{Philosophical Transactions of the Royal Society A}
#' 231: 289-337.
#'
#' @export
fit_lrt_bw <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_lrt_bw", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_lrt_bw", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_lrt_bw")
  hp <- .lrt_bw_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_lrt_bw: X_train must contain at least hp$min_features features")
  }
  if (ncol(X_train) < 2L) {
    stop("fit_lrt_bw: requires at least 2 features to estimate a class covariance")
  }

  anchor_idx <- .lrt_bw_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  if (length(anchor_idx$case) < 2L || length(anchor_idx$control) < 2L) {
    stop("fit_lrt_bw: needs at least 2 case and 2 control anchors to estimate a class covariance")
  }
  case_anchors <- X_train[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_train[anchor_idx$control, , drop = FALSE]

  # Full-universe class representation: stored for inspection / tests. The score
  # path recomputes this over the present-feature subset through the same helper,
  # so there is one geometry definition and no fit/score drift.
  repr <- .lrt_bw_class_repr(case_anchors, control_anchors, hp)

  model <- list(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    feature_universe = colnames(X_train),
    repr = repr,
    hp = hp
  )
  class(model) <- "lrt_bw_model"
  model
}


#' @title Score a Bures-Wasserstein Gaussian-OT class-conditional LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_lrt_bw}}. For one specimen the abundances are mapped to the
#' self-contained per-sample rCLR representation \eqn{z}, and the score is the
#' Gaussian log-likelihood ratio of case vs control under the frozen
#' BW-regularized class Gaussians
#' \deqn{S(z) = -\tfrac{1}{2}(z-\mu_{\mathrm{case}})^\top
#'              \Sigma_{\mathrm{case}}^{\mathrm{reg}\,-1}(z-\mu_{\mathrm{case}})
#'              - \tfrac{1}{2}\log\det\Sigma_{\mathrm{case}}^{\mathrm{reg}}
#'            + \tfrac{1}{2}(z-\mu_{\mathrm{control}})^\top
#'              \Sigma_{\mathrm{control}}^{\mathrm{reg}\,-1}
#'              (z-\mu_{\mathrm{control}})
#'            + \tfrac{1}{2}\log\det\Sigma_{\mathrm{control}}^{\mathrm{reg}}.}
#'
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. The class
#' representation (per-class mean, BW-regularized precision and log determinant)
#' is re-derived over exactly that present set from the frozen raw anchors --
#' using only frozen training data, never a scored-batch statistic -- which keeps
#' partial-overlap scores consistent and equal to the fit-time representation at
#' full overlap. The score of a row therefore depends only on that row and the
#' frozen model, and is exactly (to numerical tolerance) invariant to per-sample
#' scaling. If fewer than \code{model$hp$min_features} features are present, the
#' documented neutral score \code{0} is returned for every row.
#'
#' @param model An \code{lrt_bw_model} object returned by \code{\link{fit_lrt_bw}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values are
#'   more case-like. Scoring uses only each row's own values plus the frozen
#'   anchors and hyperparameters.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 160; p <- 30
#' L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' shift <- rep(c(1.2, -1.2), length.out = 12)
#' L[y == 1, 1:12] <- L[y == 1, 1:12] + matrix(shift, sum(y == 1), 12, byrow = TRUE)
#' X <- exp(L)
#' model <- fit_lrt_bw(X, y)
#' score_lrt_bw(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Bhatia R., Jain T., Lim Y. (2019) On the Bures-Wasserstein distance between
#' positive definite matrices. \emph{Expositiones Mathematicae} 37(2): 165-191.
#'
#' Takatsu A. (2011) Wasserstein geometry of Gaussian measures.
#' \emph{Osaka Journal of Mathematics} 48(4): 1005-1026.
#'
#' @export
score_lrt_bw <- function(model, X, meta = NULL) {
  if (!inherits(model, "lrt_bw_model")) {
    stop("score_lrt_bw: model must have class lrt_bw_model")
  }
  X <- .reo_check_matrix(X, "score_lrt_bw", "X")
  .reo_check_meta(meta, nrow(X), "score_lrt_bw", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  repr <- .lrt_bw_class_repr(
    model$case_anchors[, present, drop = FALSE],
    model$control_anchors[, present, drop = FALSE],
    model$hp
  )
  X_use <- X[, present, drop = FALSE]

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    z <- .lrt_bw_rclr(X_use[i, ])
    d_case <- z - repr$mu_case
    d_control <- z - repr$mu_control
    l_case <- -0.5 * sum(d_case * (repr$Rinv_case %*% d_case)) -
      0.5 * repr$logdet_case
    l_control <- -0.5 * sum(d_control * (repr$Rinv_control %*% d_control)) -
      0.5 * repr$logdet_control
    l_case - l_control
  }, numeric(1))

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_lrt_bw: scorer produced non-finite or wrong-length output")
  }
  out
}

.lrt_bw_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_lrt_bw: hp must be a list")
  allowed <- c("shrink", "anchor", "max_anchors_per_class", "min_features",
               "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_lrt_bw: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  # EXACT [[ ]] reads only: list `$` does partial matching, which would let a
  # prefix (e.g. a stray field) silently bind to a longer allowed name.
  shrink <- hp[["shrink"]]
  if (is.null(shrink)) shrink <- 0.5
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0 || shrink >= 1) {
    stop("fit_lrt_bw: hp$shrink must be a single number in [0, 1)")
  }

  anchor <- hp[["anchor"]]
  if (is.null(anchor)) anchor <- "bw_barycenter"
  if (!is.character(anchor) || length(anchor) != 1L || is.na(anchor) ||
      !anchor %in% c("bw_barycenter", "pooled")) {
    stop("fit_lrt_bw: hp$anchor must be one of 'bw_barycenter' or 'pooled'")
  }

  max_anchors_per_class <- hp[["max_anchors_per_class"]]
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_lrt_bw: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_lrt_bw: hp$min_features must be a positive integer")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_lrt_bw: hp$eps must be a positive finite number")
  }

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_lrt_bw: hp$seed must be a non-negative integer")
  }

  list(
    shrink = as.numeric(shrink),
    anchor = as.character(anchor),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

# Deterministic per-class anchor subsampling with global-RNG save/restore
# (mirrors the lrt-copula pattern). When neither class exceeds the cap the data
# are used as-is and the RNG is never touched.
.lrt_bw_anchor_indices <- function(y, max_anchors_per_class, seed) {
  case_idx <- which(y == 1L)
  control_idx <- which(y == 0L)
  needs_subsample <- length(case_idx) > max_anchors_per_class ||
    length(control_idx) > max_anchors_per_class
  if (!needs_subsample) {
    return(list(case = case_idx, control = control_idx))
  }

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

  set.seed(seed)
  if (length(case_idx) > max_anchors_per_class) {
    case_idx <- sort(sample(case_idx, max_anchors_per_class, replace = FALSE))
  }
  if (length(control_idx) > max_anchors_per_class) {
    control_idx <- sort(sample(control_idx, max_anchors_per_class,
                               replace = FALSE))
  }
  list(case = case_idx, control = control_idx)
}

# Self-contained per-sample robust CLR of one specimen. The geometric-mean
# centring uses ONLY this specimen's own nonzero parts, so it is a within-sample
# transform: z(c*v) = z(v) for any c > 0 (exact per-sample scale-invariance), and
# zero parts map to 0. An all-zero (or single-nonzero) specimen maps to the
# all-zero vector. This deliberately does NOT reuse any package rclr helper that
# centres using cross-sample / reference statistics.
.lrt_bw_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Row-wise rCLR of an anchor matrix. rbind (not t(vapply)) keeps an n x p shape
# even when n == 1 or p == 1, where the collapsed forms would break alignment.
.lrt_bw_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .lrt_bw_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Numerically symmetric copy: cleans the float asymmetry of products that are
# symmetric in exact arithmetic.
.lrt_bw_symm <- function(S) {
  (S + t(S)) / 2
}

# Symmetric PSD matrix square root via eigendecomposition:
# A^{1/2} = V diag(sqrt(max(lambda, 0))) V^T. (d * t(V) scales row i of V^T by
# d[i], i.e. equals diag(d) %*% V^T.)
.lrt_bw_sqrt <- function(S) {
  e <- eigen(.lrt_bw_symm(S), symmetric = TRUE)
  d <- sqrt(pmax(e$values, 0))
  .lrt_bw_symm(e$vectors %*% (d * t(e$vectors)))
}

# Symmetric inverse square root via eigendecomposition with tiny eigenvalues
# floored at `floor` so the inverse stays finite:
# A^{-1/2} = V diag(1/sqrt(max(lambda, floor))) V^T.
.lrt_bw_invsqrt <- function(S, floor) {
  e <- eigen(.lrt_bw_symm(S), symmetric = TRUE)
  d <- 1 / sqrt(pmax(e$values, floor))
  .lrt_bw_symm(e$vectors %*% (d * t(e$vectors)))
}

# Floor a symmetric matrix to positive definite by clamping its eigenvalues at
# `eps` and reconstructing. Used to make the rank-deficient (p > n) sample
# covariances admissible bases for the BW geodesic operations.
.lrt_bw_make_pd <- function(S, eps) {
  e <- eigen(.lrt_bw_symm(S), symmetric = TRUE)
  lam <- pmax(e$values, eps)
  .lrt_bw_symm(e$vectors %*% (lam * t(e$vectors)))
}

# Bures-Wasserstein (McCann / OT) geodesic between SPD A and B at parameter t:
#   T  = A^{-1/2} (A^{1/2} B A^{1/2})^{1/2} A^{-1/2}   (OT map A -> B)
#   Mt = (1 - t) I + t T
#   BWgeo(A, B, t) = Mt A Mt   (symmetrized)
# Endpoints are returned exactly (t<=0 -> A, t>=1 -> B). A and B must be PD; the
# caller floors them first. The equal-weight BW barycenter is BWgeo(A, B, 0.5).
.lrt_bw_geo <- function(A, B, t, floor) {
  A <- .lrt_bw_symm(A)
  B <- .lrt_bw_symm(B)
  if (t <= 0) return(A)
  if (t >= 1) return(B)
  A_sqrt <- .lrt_bw_sqrt(A)
  A_isqrt <- .lrt_bw_invsqrt(A, floor)
  inner <- .lrt_bw_sqrt(.lrt_bw_symm(A_sqrt %*% B %*% A_sqrt))
  Tmap <- .lrt_bw_symm(A_isqrt %*% inner %*% A_isqrt)
  p <- nrow(A)
  Mt <- (1 - t) * diag(p) + t * Tmap
  .lrt_bw_symm(Mt %*% A %*% Mt)
}

# Cholesky of a (regularized, already-PD) covariance with a ridge guard: starts
# at no ridge and, on failure, seeds the ridge at `eps` and doubles it until chol
# succeeds. Geometric growth (rather than the additive step used for the bounded
# correlation matrices of lrt-copula) keeps the loop short regardless of the
# covariance scale. Returns the precision (chol2inv) and log determinant
# (2 * sum(log(diag(L)))).
.lrt_bw_chol_pd <- function(S, eps) {
  p <- nrow(S)
  Imat <- diag(p)
  S <- .lrt_bw_symm(S)
  ridge <- 0
  ch <- NULL
  for (k in 0:60) {
    Stry <- if (ridge == 0) S else S + ridge * Imat
    ch <- tryCatch(chol(Stry), error = function(e) NULL)
    if (!is.null(ch)) break
    ridge <- if (ridge == 0) eps else ridge * 2
  }
  if (is.null(ch)) {
    stop("fit_lrt_bw: failed to regularize a class covariance to positive definite")
  }
  list(Rinv = chol2inv(ch), logdet = 2 * sum(log(diag(ch))))
}

# Derive the two-class Gaussian representation (per-class mean + BW-regularized
# precision and log determinant) from raw anchor matrices over their current
# (full or present) feature set. Used by both fit (full universe) and score
# (present subset), so the geometry is defined once and there is no fit/score
# drift.
.lrt_bw_class_repr <- function(case_anchors, control_anchors, hp) {
  Z_case <- .lrt_bw_rclr_matrix(case_anchors)
  Z_control <- .lrt_bw_rclr_matrix(control_anchors)
  mu_case <- colMeans(Z_case)
  mu_control <- colMeans(Z_control)
  Sigma_case <- .lrt_bw_symm(stats::cov(Z_case))
  Sigma_control <- .lrt_bw_symm(stats::cov(Z_control))

  # PD bases for the geodesic (null-space directions floored at hp$eps).
  A_case <- .lrt_bw_make_pd(Sigma_case, hp$eps)
  A_control <- .lrt_bw_make_pd(Sigma_control, hp$eps)

  # Shared anchor covariance.
  if (identical(hp$anchor, "bw_barycenter")) {
    Sigma_pool <- .lrt_bw_geo(A_case, A_control, 0.5, hp$eps)
  } else {
    n_case <- nrow(Z_case)
    n_control <- nrow(Z_control)
    Sigma_pool <- .lrt_bw_symm(
      ((n_case - 1) * Sigma_case + (n_control - 1) * Sigma_control) /
        (n_case + n_control - 2)
    )
  }
  B_pool <- .lrt_bw_make_pd(Sigma_pool, hp$eps)

  # BW-geodesic shrink of each class covariance toward the shared anchor.
  Sigma_case_reg <- .lrt_bw_geo(A_case, B_pool, hp$shrink, hp$eps)
  Sigma_control_reg <- .lrt_bw_geo(A_control, B_pool, hp$shrink, hp$eps)

  pd_case <- .lrt_bw_chol_pd(Sigma_case_reg, hp$eps)
  pd_control <- .lrt_bw_chol_pd(Sigma_control_reg, hp$eps)

  list(
    mu_case = mu_case,
    mu_control = mu_control,
    Rinv_case = pd_case$Rinv,
    logdet_case = pd_case$logdet,
    Rinv_control = pd_control$Rinv,
    logdet_control = pd_control$logdet
  )
}
