#' @title Fit class-conditional vine-copula LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data using a
#' \emph{vine copula} (pair-copula construction) for each class. It is the
#' flexible generalization of the Gaussian-copula (\code{\link{fit_lrt_copula}})
#' and Student-t-copula (\code{\link{fit_lrt_tcopula}}) discriminators: instead of
#' a single elliptical copula whose dependence is one correlation matrix, each
#' class is modelled by a vine copula -- a cascade of bivariate pair-copulas, each
#' free to take its own parametric family (Gaussian, Clayton, Gumbel, Frank) or
#' collapse to independence -- so the class dependence structure can be
#' asymmetric, tail-dependent, and pair-specific rather than purely elliptical.
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
#' vine copula built from two frozen ingredients estimated on \emph{training
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
#'   \item a per-class vine copula fitted to the pseudo-observations
#'     \eqn{u = \mathrm{PIT}_c(z)} of the training anchors via
#'     \code{rvinecopulib::vinecop}, with a truncation level
#'     (\code{hp$trunc_lvl}) that keeps only the strongest dependence trees and a
#'     restricted, parametric \code{hp$family_set} (including \code{"indep"}) so
#'     weak pairs collapse to independence. The marginals are matched out by the
#'     PIT, so the vine captures only the class-specific \emph{dependence}
#'     structure between features.
#' }
#' The single specimen score (see \code{\link{score_lrt_vinecopula}}) is the
#' vine-copula log-density ratio of case vs control; larger means more case-like.
#'
#' \strong{Design decisions (flagged for the gate).}
#' \itemize{
#'   \item \emph{Frozen top-\eqn{d} vine feature universe.} A vine copula has
#'     \eqn{O(d^2)} pair-copulas and is intractable / heavily overfit on the full
#'     hundreds-of-feature panel with only \eqn{\sim 200} anchors. The model is
#'     therefore restricted to the \eqn{d = } \code{hp$n_vine_features} features
#'     with the largest \emph{training column mean} (ties broken by column index).
#'     These highest-abundance miRNAs are also the most platform-portable, which
#'     suits the single-sample, batch-robust contract. This frozen feature set is
#'     the model's \code{feature_universe}; the vine is never evaluated outside it.
#'   \item \emph{Present-subset re-derivation.} As in the sibling copula scorers,
#'     the per-class marginals and vine are re-derived from the frozen raw anchors
#'     over exactly the feature subset present at scoring time
#'     (\code{intersect(feature_universe, colnames(X))}) -- using only frozen
#'     training data, never a scored-batch statistic -- so partial-overlap scores
#'     stay consistent and equal the fit-time representation at full overlap.
#'   \item \emph{Truncated, restricted, parametric vine.} Truncation
#'     (\code{trunc_lvl}, default \code{3L}) and a parametric \code{family_set}
#'     with an explicit independence family control overfitting and keep the
#'     log-density closed-form and deterministic.
#' }
#'
#' Fitting stores the raw case/control anchor rows restricted to the frozen
#' \eqn{d}-feature universe (so the class copula can be re-derived consistently
#' over any present-feature subset at scoring), the frozen feature universe, the
#' resolved hyperparameters, and -- for inspection and tests -- the full-universe
#' class representation. No test data and no test-time batch statistic are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{n_vine_features} positive
#'   integer \eqn{\ge 2} size of the frozen top-mean feature universe the vine is
#'   built on (default \code{10L}); \code{trunc_lvl} positive-integer vine
#'   truncation level (default \code{3L}); \code{family_set} non-empty character
#'   vector of \code{rvinecopulib} pair-copula family names (default
#'   \code{c("gaussian", "clayton", "gumbel", "frank", "indep")}); \code{pit_clamp}
#'   PIT clamp \eqn{\delta} in \eqn{(0, 0.5)} keeping the pseudo-observations in
#'   the open unit interval (default \code{1e-4}); \code{max_anchors_per_class}
#'   positive-integer anchor cap per class (default \code{200L});
#'   \code{min_features} positive-integer feature-overlap floor at scoring
#'   (default \code{3L}); \code{eps} positive numerical floor reserved for
#'   interface parity with the sibling scorers (default \code{1e-6}); and
#'   \code{seed} non-negative integer used only for deterministic anchor
#'   subsampling and as a guard around the (deterministic) vine fit while
#'   restoring the global RNG state (default \code{1L}).
#'
#' @return A plain list of class \code{lrt_vinecopula_model} containing
#'   \code{case_anchors}, \code{control_anchors} (raw anchor rows restricted to the
#'   frozen \eqn{d}-feature universe), \code{feature_universe} (the frozen top-\eqn{d}
#'   feature ids, in their training-mean order), \code{repr} (the full-universe
#'   frozen class representation), and the resolved \code{hp}.
#'
#' @details
#' The score is the vine-copula log-density ratio
#' \deqn{S(z) = \log c_{\mathrm{case}}(u_{\mathrm{case}})
#'            - \log c_{\mathrm{control}}(u_{\mathrm{control}}),}
#' where \eqn{u_c = \mathrm{PIT}_c(z)} are the class-\eqn{c} pseudo-observations of
#' the specimen over its present features and \eqn{c_c} is the class-\eqn{c} vine
#' copula density (\code{rvinecopulib::dvinecop}). The vine density is by
#' construction a density on uniform margins, i.e. it is the \emph{copula}
#' (dependence-only) density: the per-feature marginals are matched out by the PIT
#' and contribute no term to the score. This matches the sibling copula scorers,
#' which are likewise copula-only -- \code{lrt-copula} via the \eqn{-I} term that
#' removes the standard-normal marginal contribution, \code{lrt-tcopula} by
#' explicitly subtracting the univariate-t marginal log-densities. A specimen is
#' therefore scored purely by how well its feature \emph{dependence} matches each
#' class. The orientation is the prespecified likelihood-ratio contrast
#' (\eqn{\log c_{\mathrm{case}} - \log c_{\mathrm{control}}}, larger = more
#' case-like) and is never flipped post hoc.
#'
#' Degenerate features are handled without any caller-side special-casing. A
#' feature that is constant across a class's training anchors (e.g. zero in all of
#' them) collapses to a single-knot marginal, so its PIT maps directly to that
#' knot's plotting position with no interpolation; its near-constant
#' pseudo-observation column carries no usable dependence and is selected to the
#' independence family by \code{vinecop}. The fit and every score therefore remain
#' finite. \code{hp$pit_clamp} is validated to the open interval \eqn{(0, 0.5)} so
#' the pseudo-observations stay strictly inside the unit interval.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30
#' L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' f <- stats::rnorm(sum(y == 1))
#' L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(1.4, 8))  # case-only block
#' X <- exp(L)
#' model <- fit_lrt_vinecopula(X, y)
#' score <- score_lrt_vinecopula(model, X)
#' }
#'
#' @references
#' Sklar A. (1959) Fonctions de repartition a n dimensions et leurs marges.
#' \emph{Publications de l'Institut de Statistique de l'Universite de Paris}
#' 8: 229-231.
#'
#' Aas K., Czado C., Frigessi A., Bakken H. (2009) Pair-copula constructions of
#' multiple dependence. \emph{Insurance: Mathematics and Economics} 44(2):
#' 182-198.
#'
#' Joe H. (2014) \emph{Dependence Modeling with Copulas}. Chapman and Hall/CRC.
#'
#' @export
fit_lrt_vinecopula <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
    stop("fit_lrt_vinecopula: package 'rvinecopulib' is required for lrt-vinecopula")
  }
  X_train <- .reo_check_matrix(X_train, "fit_lrt_vinecopula", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_lrt_vinecopula", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_lrt_vinecopula")
  hp <- .lrt_vinecopula_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_lrt_vinecopula: X_train must contain at least hp$min_features features")
  }
  if (ncol(X_train) < 2L) {
    stop("fit_lrt_vinecopula: requires at least 2 features to estimate a vine copula")
  }

  # Frozen top-d vine feature universe: the d features with the largest training
  # column mean, ties broken by column index (a stable sort on -mean keeps the
  # lower index first). The vine is built and evaluated only on these features.
  universe <- .lrt_vinecopula_feature_universe(X_train, hp$n_vine_features)
  if (length(universe) < 2L) {
    stop("fit_lrt_vinecopula: requires at least 2 features in the vine universe")
  }
  X_uni <- X_train[, universe, drop = FALSE]

  anchor_idx <- .lrt_vinecopula_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  if (length(anchor_idx$case) < 2L || length(anchor_idx$control) < 2L) {
    stop("fit_lrt_vinecopula: needs at least 2 case and 2 control anchors to estimate a vine copula")
  }
  case_anchors <- X_uni[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_uni[anchor_idx$control, , drop = FALSE]

  # Full-universe class representation: stored for inspection / tests. The score
  # path recomputes this over the present-feature subset through the same helper,
  # so there is one vine definition and no fit/score drift. Wrapped in the RNG
  # guard because vinecop() initialises the global RNG when none exists.
  repr <- .lrt_vinecopula_rng_guard(
    .lrt_vinecopula_class_repr(case_anchors, control_anchors, hp)
  )

  model <- list(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    feature_universe = universe,
    repr = repr,
    hp = hp
  )
  class(model) <- "lrt_vinecopula_model"
  model
}


#' @title Score class-conditional vine-copula LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_lrt_vinecopula}}. For one specimen the abundances are mapped to
#' the self-contained per-sample rCLR representation \eqn{z}, the present features
#' are PIT-transformed under each class's frozen marginal to pseudo-observations
#' \eqn{u_c = \mathrm{PIT}_c(z)}, and the score is the vine-copula log-density
#' ratio
#' \deqn{S(z) = \log c_{\mathrm{case}}(u_{\mathrm{case}})
#'            - \log c_{\mathrm{control}}(u_{\mathrm{control}}).}
#'
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. The class vine
#' (per-feature marginals and pair-copulas) is re-derived over exactly that present
#' set from the frozen raw anchors -- using only frozen training data, never a
#' scored-batch statistic -- which keeps partial-overlap scores consistent and
#' equal to the fit-time representation at full overlap. The score of a row
#' therefore depends only on that row and the frozen model, and is exactly
#' invariant to per-sample scaling. If fewer than \code{model$hp$min_features}
#' features are present, the documented neutral score \code{0} is returned for
#' every row.
#'
#' @param model An \code{lrt_vinecopula_model} object returned by
#'   \code{\link{fit_lrt_vinecopula}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like (the specimen's feature dependence matches the case vine
#'   better than the control vine). Scoring uses only each row's own values plus
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
#' L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(1.4, 8))
#' X <- exp(L)
#' model <- fit_lrt_vinecopula(X, y)
#' score_lrt_vinecopula(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Aas K., Czado C., Frigessi A., Bakken H. (2009) Pair-copula constructions of
#' multiple dependence. \emph{Insurance: Mathematics and Economics} 44(2):
#' 182-198.
#'
#' Joe H. (2014) \emph{Dependence Modeling with Copulas}. Chapman and Hall/CRC.
#'
#' @export
score_lrt_vinecopula <- function(model, X, meta = NULL) {
  if (!inherits(model, "lrt_vinecopula_model")) {
    stop("score_lrt_vinecopula: model must have class lrt_vinecopula_model")
  }
  if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
    stop("score_lrt_vinecopula: package 'rvinecopulib' is required for lrt-vinecopula")
  }
  X <- .reo_check_matrix(X, "score_lrt_vinecopula", "X")
  .reo_check_meta(meta, nrow(X), "score_lrt_vinecopula", "meta")

  # Preserve the frozen training-mean order of the universe; intersect() with the
  # model universe as the first argument keeps that order, never colnames(X)'s.
  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, present, drop = FALSE]

  # rvinecopulib's vinecop() (in the present-subset re-derivation) AND dvinecop()
  # (in the per-row loop) both initialise the global RNG if none exists, so the
  # whole RNG-touching body is wrapped in a save/restore guard: the score is a
  # deterministic function of the inputs and must leave the caller's RNG state
  # exactly as it was (or absent, if it was absent).
  out <- .lrt_vinecopula_rng_guard({
    repr <- .lrt_vinecopula_class_repr(
      model$case_anchors[, present, drop = FALSE],
      model$control_anchors[, present, drop = FALSE],
      model$hp
    )
    vapply(seq_len(nrow(X_use)), function(i) {
      z <- .lrt_vinecopula_rclr(X_use[i, ])
      u_case <- .lrt_vinecopula_pit_u(z, repr$marg_case, model$hp$pit_clamp)
      u_control <- .lrt_vinecopula_pit_u(z, repr$marg_control, model$hp$pit_clamp)
      l_case <- .lrt_vinecopula_logdens(u_case, repr$vine_case)
      l_control <- .lrt_vinecopula_logdens(u_control, repr$vine_control)
      l_case - l_control
    }, numeric(1))
  })

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_lrt_vinecopula: scorer produced non-finite or wrong-length output")
  }
  out
}

# Evaluate an expression while guaranteeing the global RNG state is left exactly
# as found. rvinecopulib's vinecop()/dvinecop() initialise .Random.seed when none
# exists and advance it otherwise, so any fit/score that calls them must restore
# the caller's state (or remove .Random.seed if it was absent) to stay a pure,
# side-effect-free function of its inputs. Mirrors the anchor-subsampling
# save/restore idiom used by the sibling copula scorers.
.lrt_vinecopula_rng_guard <- function(expr) {
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
  force(expr)
}

# Vine-copula log-density of one specimen over its d present features. u =
# pseudo-observations (PIT of the rCLR), vine = a frozen rvinecopulib vinecop
# object. dvinecop returns the copula density on uniform margins; log of it is the
# dependence-only log-density (the per-feature marginals are matched out by the
# PIT and contribute nothing), matching the copula-only convention of lrt-copula
# and lrt-tcopula. A pmax floor at .Machine$double.xmin keeps log finite if a
# clamped tail pushes the density to a denormalized value.
.lrt_vinecopula_logdens <- function(u, vine) {
  dens <- rvinecopulib::dvinecop(matrix(u, nrow = 1L), vine)
  log(max(dens, .Machine$double.xmin))
}

.lrt_vinecopula_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_lrt_vinecopula: hp must be a list")
  if (length(hp) > 0L) {
    if (is.null(names(hp)) || any(!nzchar(names(hp)))) {
      stop("fit_lrt_vinecopula: hp fields must be named")
    }
    if (any(duplicated(names(hp)))) {
      stop("fit_lrt_vinecopula: hp fields must be unique")
    }
  }
  allowed <- c("n_vine_features", "trunc_lvl", "family_set", "pit_clamp",
               "max_anchors_per_class", "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_lrt_vinecopula: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  n_vine_features <- hp$n_vine_features
  if (is.null(n_vine_features)) n_vine_features <- 10L
  if (!is.numeric(n_vine_features) || length(n_vine_features) != 1L ||
      !is.finite(n_vine_features) || n_vine_features < 2L ||
      n_vine_features > .Machine$integer.max ||
      n_vine_features != as.integer(n_vine_features)) {
    stop("fit_lrt_vinecopula: hp$n_vine_features must be an integer >= 2")
  }

  trunc_lvl <- hp$trunc_lvl
  if (is.null(trunc_lvl)) trunc_lvl <- 3L
  if (!is.numeric(trunc_lvl) || length(trunc_lvl) != 1L ||
      !is.finite(trunc_lvl) || trunc_lvl < 1L ||
      trunc_lvl > .Machine$integer.max ||
      trunc_lvl != as.integer(trunc_lvl)) {
    stop("fit_lrt_vinecopula: hp$trunc_lvl must be a positive integer")
  }

  family_set <- hp$family_set
  if (is.null(family_set)) {
    family_set <- c("gaussian", "clayton", "gumbel", "frank", "indep")
  }
  if (!is.character(family_set) || length(family_set) < 1L ||
      any(is.na(family_set)) || any(!nzchar(family_set))) {
    stop("fit_lrt_vinecopula: hp$family_set must be a non-empty character vector of family names")
  }
  bad <- setdiff(family_set, .lrt_vinecopula_known_families())
  if (length(bad) > 0L) {
    stop("fit_lrt_vinecopula: hp$family_set has unknown family name(s): ",
         paste(bad, collapse = ", "))
  }

  pit_clamp <- hp$pit_clamp
  if (is.null(pit_clamp)) pit_clamp <- 1e-4
  if (!is.numeric(pit_clamp) || length(pit_clamp) != 1L ||
      !is.finite(pit_clamp) || pit_clamp <= 0 || pit_clamp >= 0.5) {
    stop("fit_lrt_vinecopula: hp$pit_clamp must be a single number in (0, 0.5)")
  }

  max_anchors_per_class <- hp$max_anchors_per_class
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_lrt_vinecopula: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_lrt_vinecopula: hp$min_features must be a positive integer")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_lrt_vinecopula: hp$eps must be a positive finite number")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_lrt_vinecopula: hp$seed must be a non-negative integer")
  }

  list(
    n_vine_features = as.integer(n_vine_features),
    trunc_lvl = as.integer(trunc_lvl),
    family_set = as.character(family_set),
    pit_clamp = as.numeric(pit_clamp),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

# The rvinecopulib pair-copula family names accepted by vinecop()'s family_set,
# including the group aliases. Validated here so a bad family is rejected with a
# package-style error before vinecop() is ever called.
.lrt_vinecopula_known_families <- function() {
  c("indep", "gaussian", "t", "clayton", "gumbel", "frank", "joe",
    "bb1", "bb6", "bb7", "bb8", "tll",
    "all", "parametric", "nonparametric", "onepar", "twopar", "elliptical",
    "archimedean", "itau", "ellipt")
}

# Frozen top-d feature universe: the d features with the largest TRAINING COLUMN
# MEAN, ties broken by column index. order(-mean) on a stable sort keeps the lower
# index first among ties; the returned ids are in decreasing-mean order.
.lrt_vinecopula_feature_universe <- function(X_train, d) {
  cm <- colMeans(X_train)
  ord <- order(cm, seq_along(cm), decreasing = c(TRUE, FALSE),
               method = "radix")
  keep <- ord[seq_len(min(d, length(ord)))]
  colnames(X_train)[keep]
}

# Deterministic per-class anchor subsampling with global-RNG save/restore
# (mirrors the lrt-copula / lrt-tcopula pattern). When neither class exceeds the
# cap the data are used as-is and the RNG is never touched.
.lrt_vinecopula_anchor_indices <- function(y, max_anchors_per_class, seed) {
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
.lrt_vinecopula_rclr <- function(v) {
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
.lrt_vinecopula_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .lrt_vinecopula_rclr(A[i, ])
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
# CDF used as the PIT (see .lrt_vinecopula_pit_u): it equals the step plotting-
# position ECDF u = #{train <= x}/(n+1) at the data points but is CONTINUOUS in x,
# which is what makes the score exactly (to 1e-8) invariant to per-sample scaling.
# A hard step ECDF would jump by O(1/(n+1)) whenever the irreducible float noise of
# rescaling (fl(c*v) != c*fl(v) bit-for-bit, so rCLR(c*v) differs from rCLR(v) at
# ~1e-16) pushes z across a training breakpoint.
.lrt_vinecopula_marginal <- function(zj) {
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
# frozen per-class marginals, returning pseudo-observations u in (0, 1) for the
# vine. u is the linearly interpolated plotting-position ECDF (constant
# extrapolation past the training min/max via approx rule = 2), then clamped to
# [delta, 1-delta] so dvinecop stays in the open unit cube. Unlike the elliptical
# siblings (which map u onto Gaussian / Student-t scores), the vine consumes the
# uniform u directly. A degenerate feature whose training column is constant has a
# single knot and PITs to that knot's plotting position.
.lrt_vinecopula_pit_u <- function(z, marg, pit_clamp) {
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
    u
  }, numeric(1))
}

# Fit one class vine on the anchors' pseudo-observations. The PITs are clamped to
# (delta, 1-delta) so the pseudo-observations are strictly inside the unit cube.
# Selection (selcrit / MST structure) is deterministic given inputs; the seed is
# saved/restored anyway as a guard. trunc_lvl truncates to the strongest trees and
# the restricted family_set (with "indep") lets weak pairs collapse to
# independence -- both control overfitting on the d-dimensional, ~200-anchor vine.
.lrt_vinecopula_fit_vine <- function(U, hp) {
  d <- ncol(U)
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
  set.seed(hp$seed)
  rvinecopulib::vinecop(
    data = U,
    var_types = rep("c", d),
    family_set = hp$family_set,
    trunc_lvl = hp$trunc_lvl,
    keep_data = FALSE
  )
}

# Derive the two-class vine-copula representation (per-feature frozen marginals +
# fitted vine) from raw anchor matrices over their current (full or present)
# feature set. Used by both fit (full universe) and score (present subset), so the
# vine is defined once and there is no fit/score drift. Returns the per-feature
# frozen plotting-position marginals needed to PIT a new specimen, plus each
# class's fitted vinecop object.
.lrt_vinecopula_class_repr <- function(case_anchors, control_anchors, hp) {
  build <- function(A) {
    Z <- .lrt_vinecopula_rclr_matrix(A)
    p <- ncol(Z)
    n <- nrow(Z)
    knots <- lapply(seq_len(p), function(j) .lrt_vinecopula_marginal(Z[, j]))
    marg <- list(knots = knots, n = n)
    U <- do.call(rbind, lapply(seq_len(n), function(i) {
      .lrt_vinecopula_pit_u(Z[i, ], marg, hp$pit_clamp)
    }))
    vine <- .lrt_vinecopula_fit_vine(U, hp)
    list(marg = marg, vine = vine)
  }
  rc <- build(case_anchors)
  r0 <- build(control_anchors)
  list(
    marg_case = rc$marg,
    marg_control = r0$marg,
    vine_case = rc$vine,
    vine_control = r0$vine
  )
}
