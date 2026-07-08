#' @title Fit class-conditional Student-t-copula LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data using a
#' \emph{Student-t copula} for each class. It is the heavy-tailed generalization
#' of the Gaussian-copula discriminator (\code{\link{fit_lrt_copula}}): the
#' per-sample representation, the empirical marginals, and the shrink/ridge
#' correlation are identical, but the copula scores are Student-t quantiles
#' \eqn{q = t_\nu^{-1}(u)} (shared degrees of freedom \eqn{\nu = } \code{hp$df})
#' and the dependence is scored by the Student-t-copula log-density rather than
#' the Gaussian one. Because \eqn{t_\nu^{-1} \to \Phi^{-1}} as \eqn{\nu \to
#' \infty}, this method reduces to \code{lrt-copula} in the Gaussian limit.
#'
#' Each specimen is first mapped to a self-contained, per-sample robust
#' centred-log-ratio (rCLR) representation \eqn{z}, where \eqn{z_j = \log v_j -
#' \mathrm{mean}_{k:\,v_k>0}\log v_k} on the nonzero parts of the specimen and
#' \eqn{z_j = 0} on its zero parts. Because the centring uses only that
#' specimen's own nonzero values, \eqn{z(c\,v) = z(v)} for any \eqn{c>0}: the
#' representation -- and therefore every downstream score -- is exactly invariant
#' to per-sample scaling (library size / input amount).
#'
#' Each class \eqn{c \in \{\mathrm{case}, \mathrm{control}\}} is modelled by a
#' Student-t copula built from frozen ingredients estimated on \emph{training
#' only}:
#' \enumerate{
#'   \item a per-feature empirical marginal used as a probability-integral
#'     transform (PIT). The marginal is the \emph{linearly interpolated}
#'     plotting-position ECDF through the training order statistics with Weibull
#'     positions \eqn{u_{(i)} = i / (n_c + 1)} (constant extrapolation past the
#'     training min/max), clamped to \eqn{[\delta, 1-\delta]},
#'     \eqn{\delta = } \code{hp$pit_clamp}. It equals the step plotting-position
#'     ECDF \eqn{\#\{z^{\mathrm{train}} \le x\} / (n_c + 1)} at the data points but
#'     is continuous in \eqn{x}, so the score is exactly (to numerical tolerance)
#'     invariant to per-sample scaling instead of jumping when the irreducible
#'     floating-point noise of rescaling pushes \eqn{z} across a breakpoint; and
#'   \item a per-class correlation matrix \eqn{R_c} of the Student-t scores
#'     \eqn{q = t_\nu^{-1}(u)} of the training anchors, regularized by shrinkage
#'     toward the identity (\eqn{R_c \leftarrow (1-s)R_c + sI}, \eqn{s = }
#'     \code{hp$shrink}) and a ridge that is increased until \eqn{R_c} is positive
#'     definite.
#' }
#' The marginals are matched out by the PIT, so the copula captures only the
#' class-specific \emph{dependence} structure between features. The single
#' specimen score (see \code{\link{score_lrt_tcopula}}) is the Student-t-copula
#' log-density ratio of case vs control; larger means more case-like.
#'
#' Fitting stores the raw case/control anchor rows (so the class copula can be
#' re-derived consistently over any present-feature subset at scoring), the
#' frozen training feature universe, the resolved hyperparameters, and -- for
#' inspection and tests -- the full-universe class representation. No test data
#' and no test-time batch statistic are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{shrink} correlation shrinkage
#'   toward the identity in \eqn{[0, 1)} (default \code{0.05}); \code{pit_clamp}
#'   PIT clamp \eqn{\delta} in \eqn{(0, 0.5)} keeping \eqn{t_\nu^{-1}(u)} finite
#'   (default \code{1e-4}); \code{df} Student-t degrees of freedom \eqn{\nu}, a
#'   single finite number \eqn{\ge 1} (default \code{5}; \eqn{\nu=1} is the
#'   Cauchy copula, \eqn{\nu \to \infty} the Gaussian copula);
#'   \code{max_anchors_per_class} positive-integer anchor cap per class (default
#'   \code{200L}); \code{min_features} positive-integer feature-overlap floor at
#'   scoring (default \code{3L}); \code{eps} positive ridge increment / dispersion
#'   floor (default \code{1e-6}); and \code{seed} non-negative integer used only
#'   for deterministic anchor subsampling while restoring the global RNG state
#'   (default \code{1L}).
#'
#' @return A plain list of class \code{lrt_tcopula_model} containing
#'   \code{case_anchors}, \code{control_anchors}, \code{feature_universe},
#'   \code{repr} (the full-universe frozen class representation), and the
#'   resolved \code{hp}.
#'
#' @details
#' The score is the Student-t-copula log-density ratio
#' \deqn{S(z) = \ell_{\mathrm{case}}(q_{\mathrm{case}}) -
#'   \ell_{\mathrm{control}}(q_{\mathrm{control}}),}
#' where \eqn{q_c = t_\nu^{-1}(\mathrm{PIT}_c(z))} are the class-\eqn{c} Student-t
#' scores of the specimen over its \eqn{d} present features and
#' \deqn{\ell_c(q) = \log\Gamma\!\Big(\tfrac{\nu+d}{2}\Big)
#'   + (d-1)\log\Gamma\!\Big(\tfrac{\nu}{2}\Big)
#'   - d\,\log\Gamma\!\Big(\tfrac{\nu+1}{2}\Big)
#'   - \tfrac{1}{2}\log\det R_c
#'   - \tfrac{\nu+d}{2}\,\log\!\Big(1 + \tfrac{q^\top R_c^{-1} q}{\nu}\Big)
#'   + \tfrac{\nu+1}{2}\sum_j \log\!\Big(1 + \tfrac{q_j^2}{\nu}\Big).}
#' This is the multivariate-t log-density minus the sum of the univariate-t
#' marginal log-densities, i.e. the copula (dependence-only) log-density
#' \eqn{\log f_{\mathrm{MVt}}(q; R_c, \nu) - \sum_j \log f_t(q_j; \nu)}; the
#' \eqn{\tfrac{d}{2}\log(\nu\pi)} terms cancel between the joint and the
#' marginals. The \eqn{\nu/d}-only \eqn{\log\Gamma} constant is identical for the
#' two classes (shared \eqn{\nu}, same present \eqn{d}) and cancels in the
#' difference; it is kept in each \eqn{\ell_c} for clarity. As \eqn{\nu \to
#' \infty}, \eqn{-\tfrac{\nu+d}{2}\log(1+\tfrac{q^\top R_c^{-1} q}{\nu}) \to
#' -\tfrac{1}{2} q^\top R_c^{-1} q} and \eqn{\tfrac{\nu+1}{2}\sum_j
#' \log(1+\tfrac{q_j^2}{\nu}) \to \tfrac{1}{2} q^\top q}, so \eqn{\ell_c \to
#' -\tfrac{1}{2} q^\top (R_c^{-1} - I) q - \tfrac{1}{2}\log\det R_c}, the
#' Gaussian-copula log-density.
#'
#' This method is a deliberately \emph{tail-aware} dependence discriminator and is
#' the intended complement to \code{\link{fit_lrt_copula}} (the Gaussian-copula
#' limit). At finite \eqn{\nu} the Student-t copula with \eqn{R_c = I} is \emph{not}
#' the independence copula: the marginal-correction term
#' \eqn{\tfrac{\nu+1}{2}\sum_j \log(1+q_j^2/\nu)} does not cancel the joint term
#' exactly (it does only as \eqn{\nu \to \infty}), so \eqn{\ell_c} assigns high
#' density to \emph{jointly tail-extreme} score vectors. The density is
#' sign-invariant in each coordinate and permutation-invariant (it depends on the
#' \eqn{q_j} only through \eqn{q_j^2}), but it is \emph{not} purely radial: only the
#' joint-t term depends on \eqn{q^\top R_c^{-1} q} alone, whereas the
#' coordinatewise marginal term is concave in each \eqn{q_j^2}, so at a fixed radius
#' it \emph{rewards multi-coordinate tail extremeness} (mass spread across several
#' features) over a single dominant coordinate. For example, at \eqn{R_c = I},
#' \eqn{\nu = 5}, the vectors \eqn{(3, 0)} and \eqn{(\sqrt{4.5}, \sqrt{4.5})} share
#' the same quadratic form \eqn{q^\top R_c^{-1} q = 9} (Euclidean radius
#' \eqn{\lVert q \rVert = 3}) but score \eqn{-0.415} and \eqn{+0.347}
#' respectively. The LRT
#' therefore responds to class-specific \emph{joint-tail} behaviour as well as
#' class-specific correlation, which is the entire point of the heavy-tailed
#' generalization. A practical consequence is that, unlike the Gaussian copula
#' (whose exact \eqn{-I} term removes the marginal contribution), this score is
#' \emph{not} designed to capture a pure per-feature mean/location shift between
#' classes: such a shift is matched out by the per-class PIT, and a specimen that is
#' marginally extreme under the opposite class's marginals can be rewarded by that
#' class's tail. The score is best read as a tail-aware dependence contrast; its
#' empirical discrimination on any given cohort is established by the benchmark, not
#' assumed. The orientation is the prespecified likelihood-ratio contrast
#' (\eqn{\ell_{\mathrm{case}} - \ell_{\mathrm{control}}}, larger = more
#' case-like) and is never flipped post hoc.
#'
#' Degenerate features are handled without any caller-side special-casing. A
#' feature that is constant across a class's training anchors (e.g. zero in all of
#' them) collapses to a single-knot marginal, so its PIT maps directly to that
#' knot's plotting position (Student-t score \eqn{\approx 0}) with no
#' interpolation; and its zero-variance score column -- which makes
#' \code{stats::cor} return \code{NA} -- is forced to the independent (identity)
#' row before the correlation is shrunk and ridge-regularized to positive
#' definite. The fit and every score therefore remain finite. \code{hp$pit_clamp}
#' is validated to the open interval \eqn{(0, 0.5)} so \eqn{t_\nu^{-1}} stays
#' finite, and \code{hp$df} is validated to \eqn{\ge 1} for tail stability.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30
#' L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' f <- stats::rnorm(sum(y == 1))
#' L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(c(1.4, -1.4), 4))  # case block
#' X <- exp(L)
#' model <- fit_lrt_tcopula(X, y, hp = list(df = 5))
#' score <- score_lrt_tcopula(model, X)
#' }
#'
#' @references
#' Sklar A. (1959) Fonctions de repartition a n dimensions et leurs marges.
#' \emph{Publications de l'Institut de Statistique de l'Universite de Paris}
#' 8: 229-231.
#'
#' Demarta S., McNeil A.J. (2005) The t copula and related copulas.
#' \emph{International Statistical Review} 73(1): 111-129.
#'
#' Joe H. (2014) \emph{Dependence Modeling with Copulas}. Chapman and Hall/CRC.
#'
#' @export
fit_lrt_tcopula <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_lrt_tcopula", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_lrt_tcopula", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_lrt_tcopula")
  hp <- .lrt_tcopula_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_lrt_tcopula: X_train must contain at least hp$min_features features")
  }
  if (ncol(X_train) < 2L) {
    stop("fit_lrt_tcopula: requires at least 2 features to estimate a copula correlation")
  }

  anchor_idx <- .lrt_tcopula_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  if (length(anchor_idx$case) < 2L || length(anchor_idx$control) < 2L) {
    stop("fit_lrt_tcopula: needs at least 2 case and 2 control anchors to estimate a copula correlation")
  }
  case_anchors <- X_train[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_train[anchor_idx$control, , drop = FALSE]

  # Full-universe class representation: stored for inspection / tests. The score
  # path recomputes this over the present-feature subset through the same helper,
  # so there is one copula definition and no fit/score drift.
  repr <- .lrt_tcopula_class_repr(case_anchors, control_anchors, hp)

  model <- list(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    feature_universe = colnames(X_train),
    repr = repr,
    hp = hp
  )
  class(model) <- "lrt_tcopula_model"
  model
}


#' @title Score class-conditional Student-t-copula LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_lrt_tcopula}}. For one specimen the abundances are mapped to
#' the self-contained per-sample rCLR representation \eqn{z}, the present features
#' are PIT-transformed under each class's frozen marginal and mapped to Student-t
#' scores \eqn{q_c = t_\nu^{-1}(\mathrm{PIT}_c(z))}, and the score is the
#' Student-t-copula log-density ratio \eqn{S(z) = \ell_{\mathrm{case}}(q_{\mathrm{
#' case}}) - \ell_{\mathrm{control}}(q_{\mathrm{control}})} (see
#' \code{\link{fit_lrt_tcopula}} for \eqn{\ell_c}).
#'
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. The class copula
#' (per-feature marginals and correlation) is re-derived over exactly that present
#' set from the frozen raw anchors -- using only frozen training data, never a
#' scored-batch statistic -- which keeps partial-overlap scores consistent and
#' equal to the fit-time representation at full overlap. The score of a row
#' therefore depends only on that row and the frozen model, and is exactly
#' invariant to per-sample scaling. If fewer than \code{model$hp$min_features}
#' features are present, the documented neutral score \code{0} is returned for
#' every row.
#'
#' @param model An \code{lrt_tcopula_model} object returned by
#'   \code{\link{fit_lrt_tcopula}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like (the specimen's feature dependence matches the case copula
#'   better than the control copula). Scoring uses only each row's own values plus
#'   the frozen anchors and hyperparameters.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30
#' L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' f <- stats::rnorm(sum(y == 1))
#' L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(c(1.4, -1.4), 4))
#' X <- exp(L)
#' model <- fit_lrt_tcopula(X, y, hp = list(df = 5))
#' score_lrt_tcopula(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Demarta S., McNeil A.J. (2005) The t copula and related copulas.
#' \emph{International Statistical Review} 73(1): 111-129.
#'
#' Joe H. (2014) \emph{Dependence Modeling with Copulas}. Chapman and Hall/CRC.
#'
#' @export
score_lrt_tcopula <- function(model, X, meta = NULL) {
  if (!inherits(model, "lrt_tcopula_model")) {
    stop("score_lrt_tcopula: model must have class lrt_tcopula_model")
  }
  X <- .reo_check_matrix(X, "score_lrt_tcopula", "X")
  .reo_check_meta(meta, nrow(X), "score_lrt_tcopula", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  nu <- model$hp$df
  repr <- .lrt_tcopula_class_repr(
    model$case_anchors[, present, drop = FALSE],
    model$control_anchors[, present, drop = FALSE],
    model$hp
  )
  X_use <- X[, present, drop = FALSE]

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    z <- .lrt_tcopula_rclr(X_use[i, ])
    q_case <- .lrt_tcopula_pit_q(z, repr$marg_case, model$hp$pit_clamp, nu)
    q_control <- .lrt_tcopula_pit_q(z, repr$marg_control, model$hp$pit_clamp, nu)
    l_case <- .lrt_tcopula_logdens(q_case, repr$Rinv_case, repr$logdet_case, nu)
    l_control <- .lrt_tcopula_logdens(q_control, repr$Rinv_control,
                                      repr$logdet_control, nu)
    l_case - l_control
  }, numeric(1))

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_lrt_tcopula: scorer produced non-finite or wrong-length output")
  }
  out
}

# Student-t-copula log-density of one specimen over its d present features.
# q = t-scores (qt of the PIT), Rinv = class precision, logdet = log det of the
# class correlation, nu = degrees of freedom. This is the multivariate-t log-
# density minus the sum of univariate-t marginal log-densities (the copula =
# joint / product-of-marginals); the (d/2) log(nu*pi) terms cancel between joint
# and marginals. log1p is used for numerical safety. The nu/d-only lgamma
# constant is identical across classes (shared nu, same d) and cancels in the
# case-control score difference; it is retained here for clarity.
.lrt_tcopula_logdens <- function(q, Rinv, logdet, nu) {
  d <- length(q)
  quad <- as.numeric(crossprod(q, Rinv %*% q))   # q^T Rinv q
  lgamma((nu + d) / 2) + (d - 1) * lgamma(nu / 2) - d * lgamma((nu + 1) / 2) -
    0.5 * logdet -
    ((nu + d) / 2) * log1p(quad / nu) +
    ((nu + 1) / 2) * sum(log1p(q^2 / nu))
}

.lrt_tcopula_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_lrt_tcopula: hp must be a list")
  allowed <- c("shrink", "pit_clamp", "df", "max_anchors_per_class",
               "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_lrt_tcopula: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  shrink <- hp$shrink
  if (is.null(shrink)) shrink <- 0.05
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0 || shrink >= 1) {
    stop("fit_lrt_tcopula: hp$shrink must be a single number in [0, 1)")
  }

  pit_clamp <- hp$pit_clamp
  if (is.null(pit_clamp)) pit_clamp <- 1e-4
  if (!is.numeric(pit_clamp) || length(pit_clamp) != 1L ||
      !is.finite(pit_clamp) || pit_clamp <= 0 || pit_clamp >= 0.5) {
    stop("fit_lrt_tcopula: hp$pit_clamp must be a single number in (0, 0.5)")
  }

  df <- hp$df
  if (is.null(df)) df <- 5
  # Real-valued df is allowed (qt/dt accept non-integer df); df >= 1 is required
  # for tail stability (df = 1 is the Cauchy copula; df < 1 gives heavier-than-
  # Cauchy marginals that destabilize qt and the quadratic-form weighting).
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1) {
    stop("fit_lrt_tcopula: hp$df must be a single finite number >= 1")
  }

  max_anchors_per_class <- hp$max_anchors_per_class
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_lrt_tcopula: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_lrt_tcopula: hp$min_features must be a positive integer")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_lrt_tcopula: hp$eps must be a positive finite number")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_lrt_tcopula: hp$seed must be a non-negative integer")
  }

  list(
    shrink = as.numeric(shrink),
    pit_clamp = as.numeric(pit_clamp),
    df = as.numeric(df),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

# Deterministic per-class anchor subsampling with global-RNG save/restore
# (mirrors the ig-fisherrao / lrt-copula pattern). When neither class exceeds the
# cap the data are used as-is and the RNG is never touched.
.lrt_tcopula_anchor_indices <- function(y, max_anchors_per_class, seed) {
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
.lrt_tcopula_rclr <- function(v) {
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
.lrt_tcopula_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .lrt_tcopula_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Per-feature continuous plotting-position empirical marginal. For feature j the
# training rCLR values are sorted and assigned Weibull plotting positions
# u_(i) = i / (n + 1); exact ties (e.g. several anchors with rCLR 0 for a feature
# that is zero in all of them) are collapsed to a single knot carrying the mean of
# their plotting positions, so the knots are strictly increasing in x and the
# interpolant is monotone. The knots define the linearly interpolated empirical
# CDF used as the PIT (see .lrt_tcopula_pit_q): it equals the step plotting-
# position ECDF u = #{train <= x}/(n+1) at the data points but is CONTINUOUS in x,
# which is what makes the score exactly (to 1e-8) invariant to per-sample scaling.
# A hard step ECDF would jump by O(1/(n+1)) whenever the irreducible float noise of
# rescaling (fl(c*v) != c*fl(v) bit-for-bit, so rCLR(c*v) differs from rCLR(v) at
# ~1e-16) pushes z across a training breakpoint.
.lrt_tcopula_marginal <- function(zj) {
  n <- length(zj)
  s <- sort(zj)
  pp <- seq_len(n) / (n + 1)
  ux <- unique(s)
  if (length(ux) == n) {
    list(x = s, u = pp, n = n)
  } else {
    uu <- vapply(ux, function(val) mean(pp[s == val]), numeric(1))
    list(x = ux, u = uu, n = n)
  }
}

# Probability-integral transform of a single specimen's rCLR vector under the
# frozen per-class marginals, mapped to Student-t scores q = qt(u, df). u is the
# linearly interpolated plotting-position ECDF (constant extrapolation past the
# training min/max via approx rule = 2), then clamped to [delta, 1-delta] so qt
# stays finite. qt is a smooth monotone transform of u, so per-sample scale-
# invariance is preserved exactly as in lrt-copula. A degenerate feature whose
# training column is constant has a single knot and PITs to that knot's plotting
# position (q = 0 when it is the midpoint).
.lrt_tcopula_pit_q <- function(z, marg, pit_clamp, df) {
  hi <- 1 - pit_clamp
  knots <- marg$knots
  vapply(seq_along(z), function(j) {
    k <- knots[[j]]
    if (length(k$x) == 1L) {
      u <- k$u[1L]
    } else {
      u <- stats::approx(k$x, k$u, xout = z[[j]], rule = 2)$y
    }
    if (u < pit_clamp) {
      u <- pit_clamp
    } else if (u > hi) {
      u <- hi
    }
    stats::qt(u, df = df)
  }, numeric(1))
}

# Sanitized Pearson correlation of the Student-t-score matrix. A constant column
# (e.g. a feature that is zero across all anchors -> constant rCLR) yields NA in
# cor(); those entries are forced to the independent value (0 off-diagonal, 1 on
# the diagonal). The result is positive semidefinite by construction.
.lrt_tcopula_corr <- function(Q) {
  R <- suppressWarnings(stats::cor(Q))
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  R
}

# Shrink a correlation matrix toward the identity and return its inverse and log
# determinant, increasing a diagonal ridge by hp$eps until the matrix is positive
# definite (chol succeeds). For a PSD correlation matrix and shrink > 0 the ridge
# is never needed; the loop is the documented guard for shrink == 0 or numerical
# rank deficiency. logdet = 2 * sum(log(diag(L))) for the Cholesky factor L.
.lrt_tcopula_regularize <- function(R, shrink, eps) {
  p <- nrow(R)
  Imat <- diag(p)
  Rreg <- (1 - shrink) * R + shrink * Imat
  ridge <- 0
  ch <- NULL
  for (k in 0:50) {
    Rtry <- Rreg + ridge * Imat
    ch <- tryCatch(chol(Rtry), error = function(e) NULL)
    if (!is.null(ch)) break
    ridge <- ridge + eps
  }
  if (is.null(ch)) {
    stop("fit_lrt_tcopula: failed to regularize a class correlation to positive definite")
  }
  list(Rinv = chol2inv(ch), logdet = 2 * sum(log(diag(ch))))
}

# Derive the two-class Student-t-copula representation (per-feature frozen
# marginals + regularized correlation) from raw anchor matrices over their
# current (full or present) feature set. Used by both fit (full universe) and
# score (present subset), so the copula is defined once and there is no fit/score
# drift. The correlation is estimated as the Pearson correlation of the Student-t
# scores q = qt(u, df) of the training anchors -- a normal-scores-style Pearson
# plug-in, chosen for sibling consistency with lrt-copula (same R estimator, only
# the copula family differs) and because it is what makes the score converge to
# lrt-copula as df -> Inf; it is NOT the rank-based Kendall-tau elliptical estimator
# R = sin(pi*tau/2). Returns the per-feature frozen plotting-position marginals
# needed to PIT a new specimen, plus each class's regularized precision and logdet.
.lrt_tcopula_class_repr <- function(case_anchors, control_anchors, hp) {
  nu <- hp$df
  build <- function(A) {
    Z <- .lrt_tcopula_rclr_matrix(A)
    p <- ncol(Z)
    n <- nrow(Z)
    knots <- lapply(seq_len(p), function(j) .lrt_tcopula_marginal(Z[, j]))
    marg <- list(knots = knots, n = n)
    Q <- do.call(rbind, lapply(seq_len(n), function(i) {
      .lrt_tcopula_pit_q(Z[i, ], marg, hp$pit_clamp, nu)
    }))
    R <- .lrt_tcopula_corr(Q)
    pd <- .lrt_tcopula_regularize(R, hp$shrink, hp$eps)
    list(marg = marg, Rinv = pd$Rinv, logdet = pd$logdet)
  }
  rc <- build(case_anchors)
  r0 <- build(control_anchors)
  list(
    marg_case = rc$marg,
    marg_control = r0$marg,
    Rinv_case = rc$Rinv,
    Rinv_control = r0$Rinv,
    logdet_case = rc$logdet,
    logdet_control = r0$logdet
  )
}
