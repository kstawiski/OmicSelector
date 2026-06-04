#' @title Fit per-feature KDE naive-Bayes class-conditional LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data using a
#' \emph{per-feature kernel density estimate} (KDE) of each class-conditional
#' marginal under a naive-Bayes (feature-independence) assumption. Each specimen
#' is first mapped to a self-contained, per-sample robust centred-log-ratio
#' (rCLR) representation \eqn{z}, where \eqn{z_j = \log v_j -
#' \mathrm{mean}_{k:\,v_k>0}\log v_k} on the nonzero parts of the specimen and
#' \eqn{z_j = 0} on its zero parts. Because the centring uses only that
#' specimen's own nonzero values, \eqn{z(c\,v) = z(v)} for any \eqn{c>0}: the
#' representation -- and therefore every downstream score -- is exactly invariant
#' to per-sample scaling (library size / input amount).
#'
#' This is the \emph{marginal} complement to the Gaussian-copula LRT. The copula
#' matches out the per-feature marginals via a probability-integral transform and
#' scores only the residual feature \emph{dependence}; \code{lrt-nbkde} does the
#' opposite -- it scores the per-feature class-conditional \emph{marginal} density
#' ratio and ignores dependence (the naive-Bayes assumption). For each feature
#' \eqn{j} and class \eqn{c \in \{\mathrm{case}, \mathrm{control}\}} the frozen
#' marginal is a Gaussian-kernel KDE of the training rCLR values
#' \eqn{z^{\mathrm{train}}_{c,j}}:
#' \deqn{\log f_{c,j}(t) = \mathrm{logSumExp}_{i}\!\left(
#'   -\tfrac{1}{2}\Big(\tfrac{t - z^{\mathrm{train}}_{c,j,i}}{h_j}\Big)^2\right)
#'   - \log n_{c} - \log h_j - \tfrac{1}{2}\log(2\pi),}
#' evaluated directly at the specimen's own coordinate \eqn{t = z_j} (not on a
#' fixed grid), so it is exact and computable from a single specimen.
#'
#' The per-feature bandwidth \eqn{h_j} is frozen at fit. By default it is the
#' Silverman rule-of-thumb \code{stats::bw.nrd0} of the \emph{pooled} rCLR values
#' of feature \eqn{j} across \emph{both} classes (so the case and control
#' marginals of a feature are smoothed comparably), scaled by \code{hp$bw_adjust}
#' and floored at \code{hp$eps}: \eqn{h_j = \max(\mathrm{adjust}\cdot
#' \mathrm{bw.nrd0}(\mathrm{pool}_j),\, \epsilon)}. If \code{hp$bw} is supplied,
#' that positive scalar replaces \code{bw.nrd0} for every feature.
#'
#' Fitting stores the raw case/control anchor rows (so the per-feature KDEs can be
#' re-derived consistently over any present-feature subset at scoring), the frozen
#' training feature universe, the resolved hyperparameters, and -- for inspection
#' and tests -- the full-universe class representation (per-class rCLR anchor
#' values and the per-feature bandwidths). No test data and no test-time batch
#' statistic are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{bw} either \code{NULL}
#'   (default; the per-feature pooled \code{stats::bw.nrd0} bandwidth is used) or
#'   a single positive finite scalar bandwidth applied to every feature;
#'   \code{bw_adjust} positive finite multiplier on the bandwidth analogous to
#'   \code{density(adjust=)} (default \code{1.0}); \code{max_anchors_per_class}
#'   positive-integer anchor cap per class (default \code{200L});
#'   \code{min_features} positive-integer feature-overlap floor at scoring
#'   (default \code{3L}); \code{eps} positive bandwidth/density floor (default
#'   \code{1e-6}); and \code{seed} non-negative integer used only for
#'   deterministic anchor subsampling while restoring the global RNG state
#'   (default \code{1L}).
#'
#' @return A plain list of class \code{lrt_nbkde_model} containing
#'   \code{case_anchors}, \code{control_anchors}, \code{feature_universe},
#'   \code{repr} (the full-universe frozen class representation: \code{z_case},
#'   \code{z_control}, and the per-feature \code{bw}), and the resolved \code{hp}.
#'
#' @details
#' The single-specimen score (see \code{\link{score_lrt_nbkde}}) is the naive-Bayes
#' class-conditional marginal log-likelihood ratio
#' \deqn{S(z) = \sum_{j} \big[\log f_{\mathrm{case},j}(z_j)
#'   - \log f_{\mathrm{control},j}(z_j)\big],}
#' a sum of per-feature class-conditional marginal log-density ratios; larger
#' means more case-like. The bandwidth and the \eqn{-\tfrac12\log(2\pi)} constant
#' cancel within each feature's ratio, so only the kernel sums and the per-class
#' log-counts \eqn{\log n_c} survive.
#'
#' Degenerate features are handled without any caller-side special-casing. A
#' feature that is constant across a class's training anchors (e.g. zero in all of
#' them) still yields a finite KDE: the pooled \code{bw.nrd0} stays positive (it
#' falls back to a unit scale for a constant input) and is in any case floored at
#' \code{hp$eps}, and the Gaussian kernel sum is evaluated through a stable
#' log-sum-exp so each \eqn{\log f_{c,j}} is finite. The fit and every score
#' therefore remain finite.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30
#' L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, 1:8] <- L[y == 1, 1:8] + 1.2   # case-only marginal location shift
#' X <- exp(L)
#' model <- fit_lrt_nbkde(X, y)
#' score <- score_lrt_nbkde(model, X)
#' }
#'
#' @references
#' Silverman B.W. (1986) \emph{Density Estimation for Statistics and Data
#' Analysis}. Chapman and Hall, London.
#'
#' Scott D.W. (1992) \emph{Multivariate Density Estimation: Theory, Practice, and
#' Visualization}. Wiley, New York.
#'
#' Hand D.J., Yu K. (2001) Idiot's Bayes -- not so stupid after all?
#' \emph{International Statistical Review} 69(3): 385-398.
#'
#' @export
fit_lrt_nbkde <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_lrt_nbkde", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_lrt_nbkde", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_lrt_nbkde")
  hp <- .lrt_nbkde_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_lrt_nbkde: X_train must contain at least hp$min_features features")
  }

  anchor_idx <- .lrt_nbkde_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  if (length(anchor_idx$case) < 1L || length(anchor_idx$control) < 1L) {
    stop("fit_lrt_nbkde: needs at least 1 case and 1 control anchor")
  }
  case_anchors <- X_train[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_train[anchor_idx$control, , drop = FALSE]

  # Full-universe class representation: stored for inspection / tests. The score
  # path recomputes this over the present-feature subset through the same helper,
  # so there is one KDE definition and no fit/score drift.
  repr <- .lrt_nbkde_class_repr(case_anchors, control_anchors, hp)

  model <- list(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    feature_universe = colnames(X_train),
    repr = repr,
    hp = hp
  )
  class(model) <- "lrt_nbkde_model"
  model
}


#' @title Score per-feature KDE naive-Bayes class-conditional LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_lrt_nbkde}}. For one specimen the abundances are mapped to the
#' self-contained per-sample rCLR representation \eqn{z}, and the score is the
#' naive-Bayes class-conditional marginal log-likelihood ratio
#' \deqn{S(z) = \sum_{j \in \mathrm{present}} \big[\log f_{\mathrm{case},j}(z_j)
#'   - \log f_{\mathrm{control},j}(z_j)\big],}
#' where each \eqn{\log f_{c,j}} is the frozen Gaussian-kernel KDE log-density of
#' class \eqn{c} for feature \eqn{j} evaluated at the specimen's own coordinate
#' \eqn{z_j} via a numerically stable log-sum-exp.
#'
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. The per-feature KDEs
#' (per-class rCLR anchor values and the pooled bandwidths) are re-derived over
#' exactly that present set from the frozen raw anchors -- using only frozen
#' training data, never a scored-batch statistic -- which keeps partial-overlap
#' scores consistent and equal to the fit-time representation at full overlap. The
#' score of a row therefore depends only on that row and the frozen model, and is
#' exactly invariant to per-sample scaling. If fewer than
#' \code{model$hp$min_features} features are present, the documented neutral score
#' \code{0} is returned for every row.
#'
#' @param model An \code{lrt_nbkde_model} object returned by \code{\link{fit_lrt_nbkde}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity and
#'   ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values are
#'   more case-like (the specimen's per-feature marginals match the case KDEs
#'   better than the control KDEs). Scoring uses only each row's own values plus
#'   the frozen anchors, bandwidths and hyperparameters.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30
#' L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, 1:8] <- L[y == 1, 1:8] + 1.2
#' X <- exp(L)
#' model <- fit_lrt_nbkde(X, y)
#' score_lrt_nbkde(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Silverman B.W. (1986) \emph{Density Estimation for Statistics and Data
#' Analysis}. Chapman and Hall, London.
#'
#' Hand D.J., Yu K. (2001) Idiot's Bayes -- not so stupid after all?
#' \emph{International Statistical Review} 69(3): 385-398.
#'
#' @export
score_lrt_nbkde <- function(model, X, meta = NULL) {
  if (!inherits(model, "lrt_nbkde_model")) {
    stop("score_lrt_nbkde: model must have class lrt_nbkde_model")
  }
  X <- .reo_check_matrix(X, "score_lrt_nbkde", "X")
  .reo_check_meta(meta, nrow(X), "score_lrt_nbkde", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  repr <- .lrt_nbkde_class_repr(
    model$case_anchors[, present, drop = FALSE],
    model$control_anchors[, present, drop = FALSE],
    model$hp
  )
  X_use <- X[, present, drop = FALSE]

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    z <- .lrt_nbkde_rclr(X_use[i, ])
    lf_case <- .lrt_nbkde_class_logdens(z, repr$z_case, repr$bw)
    lf_control <- .lrt_nbkde_class_logdens(z, repr$z_control, repr$bw)
    sum(lf_case - lf_control)
  }, numeric(1))

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_lrt_nbkde: scorer produced non-finite or wrong-length output")
  }
  out
}

.lrt_nbkde_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_lrt_nbkde: hp must be a list")
  allowed <- c("bw", "bw_adjust", "max_anchors_per_class", "min_features",
               "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_lrt_nbkde: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  # Use exact [[ ]] extraction (NOT $) throughout: $ does partial matching on
  # lists, so hp$bw would silently resolve to hp$bw_adjust when a caller supplies
  # bw_adjust but not bw (bw is a unique prefix of bw_adjust). [[ "bw" ]] returns
  # NULL unless an exact "bw" element exists, which is the intended behaviour.
  bw <- hp[["bw"]]
  if (!is.null(bw)) {
    if (!is.numeric(bw) || length(bw) != 1L || !is.finite(bw) || bw <= 0) {
      stop("fit_lrt_nbkde: hp$bw must be NULL or a single positive finite number")
    }
    bw <- as.numeric(bw)
  }

  bw_adjust <- hp[["bw_adjust"]]
  if (is.null(bw_adjust)) bw_adjust <- 1.0
  if (!is.numeric(bw_adjust) || length(bw_adjust) != 1L ||
      !is.finite(bw_adjust) || bw_adjust <= 0) {
    stop("fit_lrt_nbkde: hp$bw_adjust must be a single positive finite number")
  }

  max_anchors_per_class <- hp[["max_anchors_per_class"]]
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_lrt_nbkde: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_lrt_nbkde: hp$min_features must be a positive integer")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_lrt_nbkde: hp$eps must be a positive finite number")
  }

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_lrt_nbkde: hp$seed must be a non-negative integer")
  }

  list(
    bw = bw,
    bw_adjust = as.numeric(bw_adjust),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

# Deterministic per-class anchor subsampling with global-RNG save/restore
# (mirrors the lrt-copula pattern). When neither class exceeds the cap the data
# are used as-is and the RNG is never touched.
.lrt_nbkde_anchor_indices <- function(y, max_anchors_per_class, seed) {
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
.lrt_nbkde_rclr <- function(v) {
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
.lrt_nbkde_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .lrt_nbkde_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Per-feature frozen bandwidths. By default each feature's bandwidth is the
# Silverman rule-of-thumb bw.nrd0 of the POOLED rCLR values of that feature across
# BOTH classes (so the two class marginals of a feature are smoothed comparably),
# scaled by hp$bw_adjust and floored at hp$eps. If hp$bw is a positive scalar it
# replaces bw.nrd0 for every feature. bw.nrd0 needs >= 2 points and returns a
# positive value even for a constant input (it falls back to a unit scale); the
# eps floor is the backstop that keeps every bandwidth strictly positive.
.lrt_nbkde_bandwidths <- function(Zc, Z0, hp) {
  p <- ncol(Zc)
  bw_fixed <- hp[["bw"]]
  if (!is.null(bw_fixed)) {
    base <- rep(bw_fixed, p)
  } else {
    base <- vapply(seq_len(p), function(j) {
      pooled <- c(Zc[, j], Z0[, j])
      if (length(pooled) < 2L) hp[["eps"]] else stats::bw.nrd0(pooled)
    }, numeric(1))
  }
  pmax(hp[["bw_adjust"]] * base, hp[["eps"]])
}

# Column-wise numerically stable log-sum-exp of an n x P matrix: returns a length-P
# vector with entry j = log(sum_i exp(E[i, j])). The per-column max is subtracted
# before exponentiating so the largest term is exp(0) = 1 and the sum is >= 1.
.lrt_nbkde_logsumexp_cols <- function(E) {
  m <- apply(E, 2L, max)
  unname(m + log(colSums(exp(sweep(E, 2L, m, FUN = "-")))))
}

# Per-feature Gaussian-kernel KDE log-density of one query rCLR vector z (length P)
# against a class rCLR anchor matrix Zmat (n x P) with per-feature bandwidths bw
# (length P). Returns a length-P vector whose j-th entry is
#   log f_j(z_j) = logSumExp_i(-0.5*((z_j - Zmat[i,j])/bw_j)^2)
#                  - log(n) - log(bw_j) - 0.5*log(2*pi).
# The squared argument makes the sign of (z_j - Zmat[i,j]) irrelevant. Every term
# is finite (bw_j > 0, n >= 1, inputs finite), so each log-density is finite.
.lrt_nbkde_class_logdens <- function(z, Zmat, bw) {
  n <- nrow(Zmat)
  diff <- sweep(Zmat, 2L, z, FUN = "-")
  diff <- sweep(diff, 2L, bw, FUN = "/")
  E <- -0.5 * diff * diff
  lse <- .lrt_nbkde_logsumexp_cols(E)
  lse - log(n) - log(bw) - 0.5 * log(2 * pi)
}

# Derive the two-class per-feature KDE representation (per-class rCLR anchor values
# + per-feature pooled bandwidths) from raw anchor matrices over their current
# (full or present) feature set. Used by both fit (full universe) and score
# (present subset), so the KDEs are defined once and there is no fit/score drift.
.lrt_nbkde_class_repr <- function(case_anchors, control_anchors, hp) {
  Zc <- .lrt_nbkde_rclr_matrix(case_anchors)
  Z0 <- .lrt_nbkde_rclr_matrix(control_anchors)
  bw <- .lrt_nbkde_bandwidths(Zc, Z0, hp)
  list(z_case = Zc, z_control = Z0, bw = bw)
}
