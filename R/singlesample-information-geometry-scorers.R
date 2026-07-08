#' @title Fit Fisher-Rao geodesic class-conditional LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data using the
#' information geometry of compositions. Each specimen is closed to proportions
#' \eqn{p = v / \sum v} (optionally additively smoothed toward the uniform
#' composition by the pseudocount, \eqn{(p + pc) / \sum (p + pc)}) and mapped to
#' the square-root embedding \eqn{u = \sqrt{p}}. Because \eqn{\sum_i p_i = 1} we
#' have \eqn{\|u\| = 1}, so every composition lands on the positive orthant of
#' the unit sphere \eqn{S^{D-1}}, the manifold on which the Fisher-Rao
#' (multinomial information) metric is the round metric. The Fisher-Rao geodesic
#' distance between two compositions is the great-circle distance of their
#' embeddings,
#' \deqn{d(p, q) = 2 \, \arccos\!\Big(\sum_i \sqrt{p_i \, q_i}\Big),}
#' the term inside being the Bhattacharyya coefficient \eqn{= \langle u_p, u_q
#' \rangle}. Closing to proportions first -- and applying the pseudocount to the
#' already-closed proportions -- makes the embedding (and hence every distance
#' and the final score) \emph{exactly} invariant to per-sample scaling (library
#' size / input amount) for any \code{pc}, because \eqn{p(c v) = p(v)} for
#' \eqn{c > 0} and the smoothing is a fixed function of the scale-invariant
#' \eqn{p}.
#'
#' Each class \eqn{c \in \{case, control\}} is modelled as a Riemannian-Gaussian
#' on the sphere: a class mean \eqn{\mu_c} and a scalar Fisher-Rao dispersion
#' \eqn{\sigma_c^2}. The mean is the extrinsic (projected) Frechet mean
#' \eqn{\mu_c = \bar u_c / \|\bar u_c\|} with \eqn{\bar u_c} the average of the
#' class embeddings; compositional data are confined to one orthant of the
#' sphere (any two embeddings have non-negative inner product, so the
#' round-sphere angle between them, \eqn{\arccos\langle u, w \rangle}, is
#' \eqn{\le \pi/2} -- half the injectivity radius \eqn{\pi}; the factor-2
#' Fisher-Rao distance \eqn{d = 2\arccos} is then \eqn{\le \pi}), where the
#' extrinsic mean is a consistent, deterministic, closed-form proxy for the
#' iterative Karcher mean. The dispersion is the mean squared geodesic distance
#' of the class embeddings to \eqn{\mu_c} (the Frechet variance), floored at
#' \code{hp$eps} so a near-degenerate class cannot produce an infinite weight.
#'
#' Fitting stores the raw case/control anchor rows (so the class representation
#' can be re-derived consistently over any present-feature subset at scoring),
#' the frozen training feature universe, the resolved hyperparameters, and -- for
#' inspection and tests -- the full-universe class means and dispersions. No test
#' data and no test-time batch statistic are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{pc} non-negative pseudocount
#'   applied as additive smoothing of the closed proportions toward the uniform
#'   composition (default \code{0}, the pure Fisher-Rao geometry; a positive
#'   value moves boundary/zero compositions off the simplex boundary while, by
#'   acting on the already-closed proportions, preserving exact scale-invariance
#'   for any \code{pc}),
#'   \code{shared_dispersion} logical (default \code{FALSE}; \code{TRUE} pools a
#'   single within-class dispersion, reducing the score to a Fisher-Rao
#'   nearest-centroid contrast), \code{eps} positive dispersion floor (default
#'   \code{1e-6}), \code{max_anchors_per_class} positive integer anchor cap per
#'   class (default \code{200L}), \code{min_features} positive integer
#'   feature-overlap floor at scoring (default \code{2L}), and \code{seed}
#'   non-negative integer used only for deterministic anchor subsampling while
#'   restoring the global RNG state (default \code{1L}).
#'
#' @return A plain list of class \code{ig_fisherrao_model} containing
#'   \code{case_anchors}, \code{control_anchors}, \code{feature_universe},
#'   \code{mu_case}, \code{mu_control}, \code{sigma2_case},
#'   \code{sigma2_control}, \code{sigma2_pooled} (all over the full training
#'   universe), and the resolved \code{hp}.
#'
#' @details
#' The single-specimen score (see \code{\link{score_ig_fisherrao}}) is the
#' class-conditional Riemannian-Gaussian log-likelihood ratio
#' \deqn{S(x) = \frac{d(x, \mu_{ctrl})^2}{2 \sigma_{ctrl}^2}
#'            - \frac{d(x, \mu_{case})^2}{2 \sigma_{case}^2},}
#' i.e. \eqn{\log p_{case}(x) - \log p_{ctrl}(x)} up to an additive constant (the
#' class log-normalizers and the log-dispersion-ratio term) that does not depend
#' on \eqn{x} and so leaves ranking, AUC, and any monotone threshold calibration
#' unchanged. Because each dispersion is estimated in the same squared-distance
#' units as the numerator, the arbitrary Fisher-Rao distance constant (the factor
#' 2 above) cancels in every \eqn{d^2/\sigma^2} ratio whenever the dispersion is
#' the estimated Frechet variance; the fixed \code{hp$eps} floor is a numerical
#' regularizer that does not rescale with the convention, so for a degenerate
#' (floor-bound) class that convention-independence is only approximate. Larger
#' \eqn{S} means geometrically closer (dispersion-weighted) to the case class.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
#' X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
#' model <- fit_ig_fisherrao(X, y)
#' score <- score_ig_fisherrao(model, X)
#' }
#'
#' @references
#' Rao CR. (1945) Information and accuracy attainable in the estimation of
#' statistical parameters. \emph{Bulletin of the Calcutta Mathematical Society}
#' 37: 81-91.
#'
#' Amari S, Nagaoka H. (2000) \emph{Methods of Information Geometry}. American
#' Mathematical Society / Oxford University Press.
#'
#' @export
fit_ig_fisherrao <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ig_fisherrao", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ig_fisherrao", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ig_fisherrao")
  hp <- .ig_fisherrao_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_ig_fisherrao: X_train must contain at least hp$min_features features")
  }

  anchor_idx <- .ig_fisherrao_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  case_anchors <- X_train[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_train[anchor_idx$control, , drop = FALSE]

  # Full-universe class representation: stored for inspection / tests. The score
  # path recomputes this over the present-feature subset through the same helper,
  # so there is one geometry definition and no fit/score drift.
  repr <- .ig_fisherrao_class_repr(case_anchors, control_anchors, hp)

  model <- list(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    feature_universe = colnames(X_train),
    mu_case = repr$mu_case,
    mu_control = repr$mu_control,
    sigma2_case = repr$sigma2_case,
    sigma2_control = repr$sigma2_control,
    sigma2_pooled = repr$sigma2_pooled,
    hp = hp
  )
  class(model) <- "ig_fisherrao_model"
  model
}


#' @title Score Fisher-Rao geodesic class-conditional LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_ig_fisherrao}}. For one specimen \eqn{x} the score is the
#' class-conditional Riemannian-Gaussian log-likelihood ratio
#' \deqn{S(x) = \frac{d(x, \mu_{ctrl})^2}{2 \sigma_{ctrl}^2}
#'            - \frac{d(x, \mu_{case})^2}{2 \sigma_{case}^2},}
#' where \eqn{d} is the Fisher-Rao geodesic distance and
#' \eqn{\mu_c, \sigma_c^2} are the class mean and dispersion. With
#' \code{model$hp$shared_dispersion = TRUE} a single pooled dispersion replaces
#' both \eqn{\sigma_c^2}, reducing the score to
#' \eqn{(d_{ctrl}^2 - d_{case}^2) / (2 \sigma_{pooled}^2)}, a Fisher-Rao
#' nearest-centroid contrast.
#'
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. Both the specimen and
#' the frozen raw anchors are closed to proportions over exactly that present set
#' and square-root embedded (applying the same pseudocount smoothing as at
#' fitting), and the class means and dispersions are re-derived
#' over the present set from the frozen anchors -- using only frozen training
#' data, never a scored-batch statistic -- which keeps partial-overlap distances
#' consistent and the score exactly invariant to per-sample scaling. The score of
#' a row therefore depends only on that row and the frozen model. If fewer than
#' \code{model$hp$min_features} features are present, the documented neutral score
#' \code{0} is returned for every row.
#'
#' @param model An \code{ig_fisherrao_model} object returned by
#'   \code{\link{fit_ig_fisherrao}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like (geometrically closer, dispersion-weighted, to the case
#'   class). Scoring uses only each row's own values plus the frozen anchors and
#'   hyperparameters.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
#' X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
#' model <- fit_ig_fisherrao(X, y)
#' score_ig_fisherrao(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Rao CR. (1945) Information and accuracy attainable in the estimation of
#' statistical parameters. \emph{Bulletin of the Calcutta Mathematical Society}
#' 37: 81-91.
#'
#' Amari S, Nagaoka H. (2000) \emph{Methods of Information Geometry}. American
#' Mathematical Society / Oxford University Press.
#'
#' @export
score_ig_fisherrao <- function(model, X, meta = NULL) {
  if (!inherits(model, "ig_fisherrao_model")) {
    stop("score_ig_fisherrao: model must have class ig_fisherrao_model")
  }
  X <- .reo_check_matrix(X, "score_ig_fisherrao", "X")
  .reo_check_meta(meta, nrow(X), "score_ig_fisherrao", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  repr <- .ig_fisherrao_class_repr(
    model$case_anchors[, present, drop = FALSE],
    model$control_anchors[, present, drop = FALSE],
    model$hp
  )
  X_use <- X[, present, drop = FALSE]

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    # Single-row slicing drops names; restore for downstream name-based helpers.
    x_i <- X_use[i, ]
    names(x_i) <- present
    u <- .ig_fisherrao_sqrt_embed(x_i, pc = model$hp$pc)
    d2_case <- .ig_fisherrao_geodesic(u, repr$mu_case)^2
    d2_ctrl <- .ig_fisherrao_geodesic(u, repr$mu_control)^2
    if (isTRUE(model$hp$shared_dispersion)) {
      (d2_ctrl - d2_case) / (2 * repr$sigma2_pooled)
    } else {
      d2_ctrl / (2 * repr$sigma2_control) - d2_case / (2 * repr$sigma2_case)
    }
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ig_fisherrao: scorer produced non-finite or wrong-length output")
  }
  out
}

.ig_fisherrao_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ig_fisherrao: hp must be a list")
  allowed <- c("pc", "shared_dispersion", "eps", "max_anchors_per_class",
               "min_features", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ig_fisherrao: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  pc <- hp$pc
  if (is.null(pc)) pc <- 0
  if (!is.numeric(pc) || length(pc) != 1L || !is.finite(pc) || pc < 0) {
    stop("fit_ig_fisherrao: hp$pc must be a non-negative finite number")
  }

  shared_dispersion <- hp$shared_dispersion
  if (is.null(shared_dispersion)) shared_dispersion <- FALSE
  if (!is.logical(shared_dispersion) || length(shared_dispersion) != 1L ||
      is.na(shared_dispersion)) {
    stop("fit_ig_fisherrao: hp$shared_dispersion must be a single TRUE/FALSE")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_ig_fisherrao: hp$eps must be a positive finite number")
  }

  max_anchors_per_class <- hp$max_anchors_per_class
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_ig_fisherrao: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 2L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_ig_fisherrao: hp$min_features must be a positive integer")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_ig_fisherrao: hp$seed must be a non-negative integer")
  }

  list(
    pc = as.numeric(pc),
    shared_dispersion = shared_dispersion,
    eps = as.numeric(eps),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    seed = as.integer(seed)
  )
}

.ig_fisherrao_anchor_indices <- function(y, max_anchors_per_class, seed) {
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

# Square-root embedding of one specimen. The specimen is FIRST closed to
# proportions p = v / sum(v); the pseudocount (when pc > 0) is then applied to
# those ALREADY-CLOSED proportions as additive simplex smoothing toward the
# uniform composition, (p + pc) / sum(p + pc). Closing first -- and smoothing a
# quantity that is already a fixed function of the scale-invariant p -- makes the
# embedding (and every distance and the score) EXACTLY invariant to per-sample
# scaling for ANY pc, because p(c*v) = p(v) for c > 0. (Adding pc to the raw
# abundances BEFORE closure would instead break exact invariance for pc > 0.) An
# all-zero specimen maps to the uniform composition (embedding sqrt(1/D) in every
# coordinate, for any pc). The result is a unit vector on the positive orthant of
# the sphere (sum p = 1 => ||u|| = 1). The pc = 0 default path is unchanged.
.ig_fisherrao_sqrt_embed <- function(x, pc) {
  s <- sum(x)
  p <- if (s > 0) x / s else rep.int(1 / length(x), length(x))
  if (pc > 0) p <- (p + pc) / sum(p + pc)
  sqrt(p)
}

.ig_fisherrao_sqrt_embed_matrix <- function(X, pc) {
  rs <- rowSums(X)
  safe <- ifelse(rs > 0, rs, 1)
  P <- X / safe
  if (any(rs <= 0)) P[rs <= 0, ] <- 1 / ncol(X)
  # Smoothing the closed proportions row-wise via sweep() over each row's own sum
  # (NOT a recycled `(P+pc)/rowSums(P+pc)`, which would divide column-wise).
  if (pc > 0) P <- sweep(P + pc, 1L, rowSums(P + pc), "/")
  sqrt(P)
}

# Fisher-Rao geodesic distance between two square-root embeddings (unit vectors
# on the sphere): 2 * arccos(<u, w>). The inner product is the Bhattacharyya
# coefficient and is clamped to [-1, 1] to absorb floating-point overshoot before
# acos (an overshoot to 1 + eps would otherwise yield NaN). For two compositions
# the inner product is non-negative, so the distance is in [0, pi].
.ig_fisherrao_geodesic <- function(u, w) {
  ip <- sum(u * w)
  if (ip > 1) ip <- 1 else if (ip < -1) ip <- -1
  2 * acos(ip)
}

# Extrinsic (projected) class mean on the sphere: normalize the average of the
# class square-root embeddings. The average of unit vectors confined to one
# orthant is non-zero in general; if it is numerically zero (e.g. all anchors
# zero over the present features) fall back to the uniform unit embedding.
.ig_fisherrao_class_mean <- function(U) {
  ubar <- colMeans(U)
  nrm <- sqrt(sum(ubar^2))
  if (!is.finite(nrm) || nrm <= 0) {
    return(rep.int(1 / sqrt(ncol(U)), ncol(U)))
  }
  ubar / nrm
}

# Squared Fisher-Rao geodesic distances of the class embeddings (rows of U) to
# the class mean mu. The raw vector (not just its mean) is returned so the pooled
# dispersion can be formed over both classes' distances.
.ig_fisherrao_sq_dists <- function(U, mu) {
  vapply(seq_len(nrow(U)), function(i) {
    .ig_fisherrao_geodesic(U[i, ], mu)^2
  }, numeric(1))
}

# Derive the two-class Fisher-Rao representation (means, per-class and pooled
# dispersions) from raw anchor matrices over their current (full or present)
# feature set. Used by both fit (full universe) and score (present subset), so
# the geometry is defined once. Dispersions are the mean squared geodesic
# distance to the class mean (the Frechet variance), floored at eps so a
# near-degenerate class cannot produce an infinite 1/sigma^2 weight.
.ig_fisherrao_class_repr <- function(case_anchors, control_anchors, hp) {
  U_case <- .ig_fisherrao_sqrt_embed_matrix(case_anchors, pc = hp$pc)
  U_control <- .ig_fisherrao_sqrt_embed_matrix(control_anchors, pc = hp$pc)
  mu_case <- .ig_fisherrao_class_mean(U_case)
  mu_control <- .ig_fisherrao_class_mean(U_control)

  d2_case <- .ig_fisherrao_sq_dists(U_case, mu_case)
  d2_control <- .ig_fisherrao_sq_dists(U_control, mu_control)

  sigma2_case <- max(mean(d2_case), hp$eps)
  sigma2_control <- max(mean(d2_control), hp$eps)
  sigma2_pooled <- max(mean(c(d2_case, d2_control)), hp$eps)

  list(
    mu_case = mu_case,
    mu_control = mu_control,
    sigma2_case = sigma2_case,
    sigma2_control = sigma2_control,
    sigma2_pooled = sigma2_pooled
  )
}
