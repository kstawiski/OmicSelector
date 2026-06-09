#' @title Forced-single-sample entropic-OT (Sinkhorn) negative control
#'
#' @description
#' Forced \eqn{n = 1}, row-by-row entropic optimal-transport (Sinkhorn) scorer
#' (family NC), a deliberately weak \strong{negative control} on the
#' single-sample discriminator roster. It scores each specimen INDEPENDENTLY (as
#' its own one-atom / Dirac distribution) against two FROZEN training reference
#' point clouds -- the case cloud and the control cloud -- in per-sample robust
#' CLR space, and reports the difference of the two entropic-kernel transport
#' energies oriented so that larger = more case-like.
#'
#' \strong{Why a separate method, and what out-of-N coupling it removes.} The
#' package already ships a MULTI-ROW Sinkhorn-OT scorer
#' (\code{\link{fit_sinkhorn_ot_scorer}} / \code{\link{score_sinkhorn_ot_scorer}}).
#' That primitive is transductive / out-of-N: its batch projection
#' (\code{.sot_project_scaled}, the \code{nrow(z_scaled) > 1} branch) builds ONE
#' joint Sinkhorn transport plan whose source marginal is
#' \code{a <- rep(1 / nrow(z_scaled), nrow(z_scaled))} -- it couples ALL scored
#' rows to the barycenter atoms and then renormalises each row by its plan mass
#' \code{row_mass <- rowSums(plan)}, so a specimen's projection (hence its score)
#' depends on which other specimens are scored alongside it. Its cost
#' normalisation (\code{.sot_cost}) also divides by the MEDIAN pairwise cost of
#' the scored batch -- a second cross-row statistic. Both make the primitive fail
#' the single-sample deployability gate
#' (\code{\link{singlesample_assert_row_equivariant}}). This method removes BOTH
#' couplings: every specimen is transported on its own as a one-atom measure
#' (the row-independent entropic-kernel softmin the primitive itself uses only in
#' its \code{nrow == 1} branch), and the cost scale is FROZEN at fit from the
#' training clouds, never recomputed from the scored set.
#'
#' \strong{Entropic-kernel transport energy.} For one specimen with rCLR vector
#' \eqn{z} (a Dirac \eqn{\delta_z}) and a frozen reference cloud
#' \eqn{\{a_j\}_{j=1}^m} with uniform weights, the entropic-kernel (Gibbs /
#' Sinkhorn) transport energy is the softmin of squared-Euclidean transport
#' costs
#' \deqn{E_\varepsilon(z, A) = -\varepsilon \log\!\Big( \tfrac{1}{m}
#'   \sum_{j=1}^m \exp(-\,C(z, a_j) / \varepsilon) \Big),\qquad
#'   C(z, a_j) = \lVert z - a_j \rVert_2^2 / \tau,}
#' where \eqn{\varepsilon} is the entropic regularisation and \eqn{\tau} the
#' frozen cost scale (median training pairwise cost). This is exactly the
#' single-source entropic-OT / Sinkhorn-potential value (the log-sum-exp the
#' multi-row primitive evaluates per atom in its single-sample branch); as
#' \eqn{\varepsilon \to 0} it tends to the nearest-atom transport cost. The score
#' is \eqn{\mathrm{orientation} \cdot (E_\varepsilon(z, \mathrm{control}) -
#'   E_\varepsilon(z, \mathrm{case}))}: a specimen closer (in entropic-OT energy)
#' to the case cloud than to the control cloud scores higher.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the
#' self-contained per-sample robust CLR over its OWN strictly-positive support
#' (geometric-mean centring on \code{v > 0}); this is exactly invariant to
#' per-specimen scaling and uses no cross-row statistic, exactly as the other
#' roster scorers do. Universe features absent from a specimen are dropped from
#' BOTH the specimen and the frozen reference atoms for that one transport (the
#' cost is computed on the present-feature subspace), so a missing feature
#' contributes no transport cost.
#'
#' \strong{Negative-control honesty.} This is EXPECTED to discriminate weakly on
#' real data and is not engineered to look strong; it nonetheless produces a
#' finite, deterministic, exactly row-equivariant score and separates a strongly
#' planted case-vs-control shift (so the wiring is testable). Degenerate inputs
#' (no present features, fewer than \code{hp$min_features} present, an all-zero
#' specimen, or an empty reference cloud) return the neutral score \code{0}.
#'
#' @references
#' Cuturi M. (2013) Sinkhorn Distances: Lightspeed Computation of Optimal
#' Transport. \emph{Advances in Neural Information Processing Systems (NeurIPS)}
#' 26:2292-2300.
#'
#' Peyre G, Cuturi M. (2019) Computational Optimal Transport.
#' \emph{Foundations and Trends in Machine Learning} 11(5-6):355-607.
#'
#' @name singlesample-sinkhorn-single
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Deliberately does NOT reuse a package
# rclr helper that centres on a fixed pseudocount or cross-sample reference
# statistics (either would break exact scale-invariance / single-sample
# deployability). Mirrors .ecod_copod_rclr.
.sinkhorn_single_rclr <- function(v) {
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
.sinkhorn_single_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .sinkhorn_single_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}


# ----------------------------------------------------------------------------
# Entropic-kernel single-source transport energy
# ----------------------------------------------------------------------------

# Entropic-OT (Sinkhorn) transport energy of a single specimen (Dirac at z) to a
# frozen reference cloud `atoms` (m x d, columns = the present-feature subspace),
# with uniform reference weights. Squared-Euclidean cost scaled by the frozen
# `tau`, then the softmin (log-sum-exp) at temperature epsilon:
#   E = -epsilon * log( mean_j exp( -||z - a_j||^2 / (tau * epsilon) ) ).
# Implemented via a numerically stable log-sum-exp (subtract the min cost), so it
# is finite for any inputs. This is the single-source entropic-OT / Sinkhorn-
# potential value -- the same Gibbs-kernel softmin the multi-row primitive uses
# in its nrow == 1 branch -- and depends only on z, the frozen atoms, and the
# frozen (tau, epsilon); no statistic of any other scored specimen enters.
.sinkhorn_single_energy <- function(z, atoms, tau, epsilon) {
  d2 <- colSums((t(atoms) - z)^2)            # ||z - a_j||^2, length m
  cost <- d2 / tau                            # frozen cost scale
  t_arg <- cost / epsilon                     # -log Gibbs kernel per atom
  m <- length(t_arg)
  t_min <- min(t_arg)
  # log( mean_j exp(-t_j) ) = -t_min + log( mean_j exp(-(t_j - t_min)) )
  lse <- -t_min + log(mean(exp(-(t_arg - t_min))))
  -epsilon * lse
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.sinkhorn_single_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_sinkhorn_single: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_sinkhorn_single: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_sinkhorn_single: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("epsilon", "min_features", "max_atoms", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_sinkhorn_single: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  epsilon <- hp$epsilon
  if (is.null(epsilon)) epsilon <- 0.1
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0) {
    stop("fit_sinkhorn_single: hp$epsilon must be a positive finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_sinkhorn_single: hp$min_features must be a positive integer")
  }

  max_atoms <- hp$max_atoms
  if (is.null(max_atoms)) max_atoms <- 500L
  if (!is.numeric(max_atoms) || length(max_atoms) != 1L ||
      !is.finite(max_atoms) || max_atoms < 1L ||
      max_atoms > .Machine$integer.max ||
      max_atoms != as.integer(max_atoms)) {
    stop("fit_sinkhorn_single: hp$max_atoms must be a positive integer")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-8
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_sinkhorn_single: hp$eps must be a positive finite number")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_sinkhorn_single: hp$seed must be a single integer")
  }

  list(
    epsilon = as.numeric(epsilon),
    min_features = as.integer(min_features),
    max_atoms = as.integer(max_atoms),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}


# Median pairwise squared-Euclidean cost scale of the pooled reference atoms,
# frozen at fit. Mirrors the spirit of .sot_cost's median normalisation but is a
# TRAINING-ONLY constant (never recomputed from a scored batch). Falls back to 1
# for a degenerate (empty / single-atom / non-positive-median) pool.
.sinkhorn_single_cost_scale <- function(atoms) {
  if (is.null(atoms) || nrow(atoms) < 2L) return(1)
  a2 <- rowSums(atoms * atoms)
  d2 <- outer(a2, a2, `+`) - 2 * tcrossprod(atoms)
  d2 <- d2[upper.tri(d2)]
  d2 <- d2[is.finite(d2) & d2 > 0]
  med <- if (length(d2)) stats::median(d2) else 1
  if (!is.finite(med) || med <= 0) 1 else med
}


# Mann-Whitney AUC of `score` against 0/1 labels `y` (P(score_case >
# score_ctrl), ties at 0.5). Self-contained so the fit does not depend on pROC.
.sinkhorn_single_auc <- function(y, score) {
  case <- score[y == 1L]
  ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}


# ----------------------------------------------------------------------------
# fit_sinkhorn_single
# ----------------------------------------------------------------------------

#' @title Fit the forced-single-sample entropic-OT (Sinkhorn) negative control
#'
#' @description
#' Fits the forced \eqn{n = 1} entropic-OT negative control. The training matrix
#' is mapped to the per-sample robust CLR (over each row's own positive support);
#' the rows are split into a frozen CASE reference cloud (\code{y_train == 1}) and
#' a frozen CONTROL reference cloud (\code{y_train == 0}). The frozen model stores
#' only these two reference clouds (capped at \code{hp$max_atoms} atoms each by a
#' deterministic seeded subsample), the frozen squared-cost scale, the entropic
#' regularisation \code{hp$epsilon}, and a single orientation scalar chosen so
#' larger = more case-like. Scoring is pure R and transports each specimen
#' independently against these frozen clouds.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{epsilon}
#'   (positive entropic-OT regularisation; default \code{0.1}),
#'   \code{min_features} (positive integer feature-overlap floor at scoring;
#'   default \code{3L}), \code{max_atoms} (positive integer cap on the atoms kept
#'   per reference cloud, seeded subsample; default \code{500L}), \code{eps}
#'   (positive finite numeric guard / scale floor; default \code{1e-8}), and
#'   \code{seed} (integer seed for the deterministic atom subsample; default
#'   \code{42L}).
#'
#' @return Object of class \code{sinkhorn_single_model}: a list with
#'   \code{feature_universe}, \code{case_atoms} / \code{control_atoms} (frozen
#'   rCLR reference clouds aligned to \code{feature_universe}), \code{cost_scale},
#'   \code{orientation}, \code{n_case} / \code{n_control}, and \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30; k <- 8
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' model <- fit_sinkhorn_single(X, y)
#' score_sinkhorn_single(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Cuturi M. (2013) Sinkhorn Distances. \emph{NeurIPS} 26:2292-2300.
#'
#' @export
fit_sinkhorn_single <- function(X_train, y_train, meta_train = NULL,
                                hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_sinkhorn_single", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_sinkhorn_single", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_sinkhorn_single")
  hp <- .sinkhorn_single_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_sinkhorn_single: X_train must contain at least hp$min_features ",
         "features")
  }

  feat_order <- colnames(X_train)
  Z <- .sinkhorn_single_rclr_matrix(X_train)        # n x p, full universe

  # Deterministic seeded subsample of each class cloud to <= max_atoms atoms.
  # RNG is local: the global .Random.seed is saved and restored around set.seed.
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
  set.seed(hp$seed)

  subsample <- function(idx) {
    if (length(idx) > hp$max_atoms) idx <- sort(sample(idx, hp$max_atoms))
    idx
  }
  case_idx <- subsample(which(y == 1L))
  control_idx <- subsample(which(y == 0L))

  case_atoms <- Z[case_idx, , drop = FALSE]
  control_atoms <- Z[control_idx, , drop = FALSE]

  # Frozen squared-cost scale from the pooled reference atoms (training only).
  cost_scale <- .sinkhorn_single_cost_scale(rbind(case_atoms, control_atoms))
  cost_scale <- max(cost_scale, hp$eps)

  # Orientation: difference of entropic energies (control - case) on the TRAINING
  # rows, oriented so the score is positively associated with case. Computed with
  # +1 orientation, then flipped if anti-correlated (training AUC < 0.5).
  full_idx <- seq_len(length(feat_order))
  train_raw <- vapply(seq_len(nrow(Z)), function(i) {
    # Same empty-support floor as the score path, so fit and score use the
    # IDENTICAL convention (empty support -> 0; flat / all other rows -> scored).
    if (!any(X_train[i, ] > 0)) return(0)            # empty positive support
    .sinkhorn_single_score_one(Z[i, ], full_idx, case_atoms, control_atoms,
                               cost_scale, hp$epsilon)
  }, numeric(1L))
  auc_pos <- .sinkhorn_single_auc(y, train_raw)
  orientation <- if (is.finite(auc_pos) && auc_pos < 0.5) -1 else 1

  model <- list(
    feature_universe = feat_order,
    case_atoms = case_atoms,
    control_atoms = control_atoms,
    cost_scale = cost_scale,
    orientation = orientation,
    n_case = nrow(case_atoms),
    n_control = nrow(control_atoms),
    hp = hp
  )
  class(model) <- "sinkhorn_single_model"
  model
}


# Raw (un-oriented) single-specimen score from the frozen clouds, on the
# present-feature subspace `present_idx` (indices into feature_universe). The
# specimen rCLR `z` and both reference clouds are restricted to those columns, so
# a feature absent from the specimen contributes no transport cost. Returns the
# difference of entropic transport energies (control - case): smaller energy to a
# cloud means "closer to it", so a larger control-minus-case difference means
# closer to the case cloud. Returns 0 (neutral) if either cloud is empty.
.sinkhorn_single_score_one <- function(z, present_idx, case_atoms, control_atoms,
                                       cost_scale, epsilon) {
  if (nrow(case_atoms) < 1L || nrow(control_atoms) < 1L) return(0)
  z_p <- z[present_idx]
  ca <- case_atoms[, present_idx, drop = FALSE]
  co <- control_atoms[, present_idx, drop = FALSE]
  e_case <- .sinkhorn_single_energy(z_p, ca, cost_scale, epsilon)
  e_ctrl <- .sinkhorn_single_energy(z_p, co, cost_scale, epsilon)
  e_ctrl - e_case
}


# ----------------------------------------------------------------------------
# score_sinkhorn_single
# ----------------------------------------------------------------------------

#' @title Score the forced-single-sample entropic-OT (Sinkhorn) negative control
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_sinkhorn_single}}. For one specimen the abundances are mapped
#' to the per-sample robust CLR over its own positive support; its entropic-OT
#' (Sinkhorn) transport energy to the frozen case cloud and to the frozen control
#' cloud is computed in pure R on the present-feature subspace, and the frozen
#' orientation is applied to \eqn{(E_{\mathrm{control}} - E_{\mathrm{case}})} so
#' larger = more case-like. No statistic of any other scored row enters, so the
#' score of a row depends only on that row and the frozen model and is exactly
#' invariant to per-specimen scaling.
#'
#' The present feature set is \code{intersect(model$feature_universe,
#' colnames(X))}; features missing from a specimen contribute no transport cost.
#' If fewer than \code{model$hp$min_features} universe features are present, or a
#' reference cloud is empty, or a specimen has empty positive support, the neutral
#' score \code{0} is returned (for that row, or for every row in the
#' below-floor / empty-cloud cases).
#'
#' @param model A \code{sinkhorn_single_model} object from
#'   \code{\link{fit_sinkhorn_single}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values
#'   are more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_sinkhorn_single(X, y)
#' score_sinkhorn_single(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Cuturi M. (2013) Sinkhorn Distances. \emph{NeurIPS} 26:2292-2300.
#'
#' @export
score_sinkhorn_single <- function(model, X, meta = NULL) {
  if (!inherits(model, "sinkhorn_single_model")) {
    stop("score_sinkhorn_single: model must have class sinkhorn_single_model")
  }
  X <- .reo_check_matrix(X, "score_sinkhorn_single", "X")
  .reo_check_meta(meta, nrow(X), "score_sinkhorn_single", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))
  if (length(present) < model$hp$min_features ||
      model$n_case < 1L || model$n_control < 1L) {
    return(rep(0, nrow(X)))
  }
  present_idx <- match(present, feat_order)

  # Align X to the frozen universe order; absent features stay 0 (no support).
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  X_use[, present] <- X[, present, drop = FALSE]

  case_atoms <- model$case_atoms
  control_atoms <- model$control_atoms
  cost_scale <- model$cost_scale
  epsilon <- model$hp$epsilon
  orientation <- model$orientation

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (it maps
    # to the origin) but is a valid specimen that MUST be transported to the frozen
    # clouds, not floored. Using all(z == 0) here would conflate the two and break
    # fit/score consistency (the fit orientation path scores flat rows).
    if (!any(X_use[i, ] > 0)) return(0)              # empty positive support
    z <- .sinkhorn_single_rclr(X_use[i, ])
    orientation * .sinkhorn_single_score_one(z, present_idx, case_atoms,
                                              control_atoms, cost_scale, epsilon)
  }, numeric(1L))

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_sinkhorn_single: scorer produced non-finite or wrong-length ",
         "output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen, score-determining R-side state of the model (feature
# order, the two reference clouds, the frozen cost scale, orientation, the cloud
# sizes, and hp). The model holds no live external-pointer object, so this is a
# deterministic snapshot suitable as the `model_digest` for
# singlesample_assert_row_equivariant() (clause (d): the bytes must be identical
# before and after scoring).
.sinkhorn_single_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    case_atoms = model$case_atoms,
    control_atoms = model$control_atoms,
    cost_scale = model$cost_scale,
    orientation = model$orientation,
    n_case = model$n_case,
    n_control = model$n_control,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
