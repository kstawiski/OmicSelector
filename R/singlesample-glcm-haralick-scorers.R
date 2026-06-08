#' @title Fit GLCM Haralick texture (inv-glcm) within-sample discriminator
#'
#' @description
#' Learns a frozen, single-sample discriminator from the 1-D Haralick texture of
#' a specimen's canonical-order robust-CLR (rCLR) abundance landscape. The
#' construction has five stages:
#' \enumerate{
#'   \item \strong{Frozen canonical feature order.} Features are sorted once
#'     by decreasing training column mean (\code{colMeans(X_train)}), ties broken
#'     by feature index (locale-independent, deterministic). The resulting
#'     permutation \code{order_perm} is reused unchanged at scoring.
#'   \item \strong{Per-specimen rCLR landscape.} Over each specimen's own positive
#'     support \eqn{\{j : v_j > 0\}}, the within-sample robust-CLR is
#'     \eqn{g_j = \log v_j - \mathrm{mean}_{k \in P}\log v_k}. Structural zeros
#'     (\eqn{v_j = 0}) are excluded; present coordinates with \eqn{g_j = 0} are
#'     kept (keyed on \code{v > 0}, NOT \code{g != 0}). Values are placed in the
#'     frozen canonical order, yielding a sequence \eqn{g_1, \dots, g_m}. The
#'     rCLR landscape \eqn{g} is exactly scale-invariant because
#'     \eqn{\log(c v_j) - \mathrm{mean}\log(c v) = \log v_j - \mathrm{mean}\log v}
#'     for any \eqn{c > 0}; the frozen-edge quantised levels, hence the score, are
#'     therefore unchanged (bit-identical in practice -- the only exception is the
#'     measure-zero set where an rCLR value lies within floating-point rounding of
#'     a frozen edge and a positive rescaling could flip its bin).
#'   \item \strong{Frozen gray-level quantisation.} \code{K = hp$n_levels} gray
#'     levels are defined by \eqn{K-1} interior quantile edges learned from the
#'     pooled training landscape values. Each specimen's landscape is quantised
#'     via left-closed binning: \code{findInterval(g, edges) + 1} (intervals
#'     \eqn{[e_i, e_{i+1})}, so a value exactly equal to an edge falls into the
#'     higher-numbered level), yielding levels in \eqn{\{1, \dots, K\}}. Because
#'     both \eqn{g} and the frozen \code{edges} live on the scale-invariant rCLR
#'     scale, the quantised level sequence is unchanged by positive per-sample
#'     scaling.
#'   \item \strong{Per-specimen 1-D GLCM + Haralick descriptors.} From the
#'     quantised sequence \eqn{q_1, \dots, q_m} a distance-1 co-occurrence matrix
#'     \eqn{C[i,j]} is built (\eqn{C[q_t, q_{t+1}]} incremented for
#'     \eqn{t = 1, \dots, m-1}), symmetrised (\eqn{C \leftarrow C + C^\top}), and
#'     normalised to a probability matrix \eqn{P = C / \mathrm{sum}(C)}. Five
#'     Haralick features are extracted: contrast, homogeneity (IDM), energy (ASM),
#'     entropy, and correlation (with degeneracy guard). This helper is defined
#'     once and used identically at fit and at score.
#'   \item \strong{Frozen standardize + ridge-LDA head.} Per-descriptor frozen
#'     \code{center}/\code{scale} (sd floored at \code{eps}; zero-variance
#'     descriptors deactivated), pooled within-class covariance ridged to PD,
#'     \eqn{w = S_w^{-1}(\mu_\text{case} - \mu_\text{ctrl})},
#'     \eqn{b = -\frac12 (\mu_\text{case} + \mu_\text{ctrl})^\top w}.
#'     Larger score = more case-like.
#' }
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease, 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity; ignored.
#' @param hp List of frozen hyperparameters (strict allowed-list, exact
#'   \code{[[ ]]} reads): \code{n_levels} K (default \code{8L}; integer in
#'   \code{[2, 32]}) -- GLCM gray levels; \code{shrink} (default \code{0.1};
#'   non-negative finite) -- ridge fraction of mean covariance diagonal;
#'   \code{min_features} (default \code{8L}; integer \eqn{\ge 3}) --
#'   per-specimen positive-support floor; \code{eps} (default \code{1e-6};
#'   positive finite) -- sd / ridge / correlation-degeneracy floor; \code{seed}
#'   (default \code{1L}; non-negative integer) -- accepted for interface
#'   uniformity, no randomness is used.
#'
#' @return A plain list of class \code{inv_glcm_model} containing
#'   \code{feature_universe}, \code{order_perm} (feature names in frozen
#'   canonical order), \code{edges} (frozen \eqn{K-1} quantile breakpoints),
#'   \code{K} (number of gray levels), \code{descriptor_dim} (5), the frozen
#'   \code{head} (\code{center}, \code{scale}, \code{active}, \code{w}, \code{b},
#'   \code{ridge}), and the resolved \code{hp}.
#'
#' @details
#' The per-specimen score is the frozen linear predictor on the standardised
#' Haralick descriptor:
#' \deqn{S(x) = b + \sum_c w_c \frac{\phi_c(x) - \mathrm{center}_c}
#'                                   {\mathrm{scale}_c},}
#' where \eqn{\phi(x)} is the 5-dim descriptor. Every quantity except
#' \eqn{\phi(x)} is frozen at fit time; \eqn{\phi(x)} depends only on
#' \eqn{x}'s own rCLR landscape, so the score is computable from a single
#' specimen and is unchanged by per-sample positive scaling (bit-identical in
#' practice; see \code{\link{score_inv_glcm}} for the measure-zero quantiser-edge
#' qualification).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 50
#' X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
#'             dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
#' y <- rep(c(0, 1), each = 60)
#' X[y == 1, 11:20] <- X[y == 1, 11:20] + 150
#' model <- fit_inv_glcm(X, y)
#' score <- score_inv_glcm(model, X)
#' }
#'
#' @references
#' Haralick RM, Shanmugam K, Dinstein I. (1973) Textural features for image
#' classification. \emph{IEEE Transactions on Systems, Man, and Cybernetics}
#' SMC-3(6): 610-621.
#'
#' Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse compositional
#' technique reveals microbial perturbations. \emph{mSystems} 4(1): e00016-19.
#'
#' @export
fit_inv_glcm <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_inv_glcm", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_inv_glcm", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_inv_glcm")
  hp <- .inv_glcm_resolve_hp(hp)

  if (ncol(X_train) < hp$min_features) {
    stop("fit_inv_glcm: X_train must contain at least hp$min_features features")
  }

  # Stage 1: frozen canonical order (decreasing training column mean, ties by
  # feature index so the result is locale-independent and deterministic).
  col_mean <- colMeans(X_train)
  order_idx <- order(-col_mean, seq_len(ncol(X_train)))
  order_perm <- colnames(X_train)[order_idx]

  # Stage 2+3: per-specimen rCLR landscapes, then pool all values to learn edges.
  # Only training specimens are used; present features are taken in frozen order.
  Xo <- X_train[, order_perm, drop = FALSE]
  pooled <- unlist(lapply(seq_len(nrow(Xo)), function(i) {
    .inv_glcm_rclr_landscape(Xo[i, ])
  }), use.names = FALSE)

  if (length(pooled) == 0L) {
    stop("fit_inv_glcm: no positive-support values in training data")
  }
  K <- hp$n_levels
  # K-1 interior quantile edges from pooled training landscape values (type 7).
  # Consecutive ties produce equal edges; left-closed binning
  # findInterval(value, edges) + 1 (intervals [e_i, e_{i+1}); a value equal to an
  # edge falls into the higher level) still maps to a valid level in {1..K}, so
  # no additional guard is needed beyond what findInterval provides.
  edges <- if (K > 1L) {
    as.numeric(stats::quantile(pooled, probs = seq_len(K - 1L) / K,
                               type = 7, names = FALSE))
  } else {
    numeric(0)
  }

  # Stage 4: Haralick descriptors for every training specimen.
  desc_list <- lapply(seq_len(nrow(Xo)), function(i) {
    .inv_glcm_descriptor(Xo[i, ], edges, K, hp$min_features, hp$eps)
  })
  keep <- !vapply(desc_list, is.null, logical(1))
  if (sum(keep) < 2L) {
    stop("fit_inv_glcm: too few training specimens have >= hp$min_features ",
         "positive-support coordinates")
  }
  Phi <- do.call(rbind, desc_list[keep])
  yk <- y[keep]
  if (sum(yk == 1L) < 1L || sum(yk == 0L) < 1L) {
    stop("fit_inv_glcm: after the min_features filter both classes must ",
         "retain >= 1 specimen")
  }

  # Stage 5: frozen standardize + ridge-LDA head.
  head <- .inv_glcm_fit_head(Phi, yk, hp$shrink, hp$eps)

  model <- list(
    feature_universe = colnames(X_train),
    order_perm       = order_perm,
    edges            = edges,
    K                = K,
    descriptor_dim   = 5L,
    head             = head,
    hp               = hp
  )
  class(model) <- "inv_glcm_model"
  model
}


#' @title Score GLCM Haralick texture (inv-glcm) within-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_inv_glcm}}. The model-universe features present in \code{X}
#' are taken in the frozen canonical order; for one specimen its within-sample
#' rCLR landscape is quantised with the frozen edges, a 1-D GLCM is built, the
#' Haralick descriptor is computed, standardised by the frozen head, and the
#' frozen LDA score is returned. Larger = more case-like.
#'
#' Scoring uses only each row's own values plus the frozen model; no test-batch
#' renormalisation, no cross-row coupling, no statistic estimated from \code{X}.
#' If fewer than \code{model$hp$min_features} model-universe features are present
#' in \code{X}, or a single specimen has fewer than \code{model$hp$min_features}
#' positive values, the documented neutral score \code{0} is returned.
#'
#' @param model An \code{inv_glcm_model} object returned by
#'   \code{\link{fit_inv_glcm}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity;
#'   ignored.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like. Each score depends only on that row's values and the
#'   frozen model. It is invariant to per-sample positive scaling: the rCLR
#'   landscape is exactly scale-invariant, so the frozen-edge quantised levels and
#'   hence the score are bit-identical, except on the measure-zero set where an
#'   rCLR value lies within floating-point rounding of a frozen edge.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 50
#' X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
#'             dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
#' y <- rep(c(0, 1), each = 60)
#' X[y == 1, 11:20] <- X[y == 1, 11:20] + 150
#' model <- fit_inv_glcm(X, y)
#' score_inv_glcm(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Haralick RM, Shanmugam K, Dinstein I. (1973) Textural features for image
#' classification. \emph{IEEE Transactions on Systems, Man, and Cybernetics}
#' SMC-3(6): 610-621.
#'
#' @export
score_inv_glcm <- function(model, X, meta = NULL) {
  if (!inherits(model, "inv_glcm_model")) {
    stop("score_inv_glcm: model must have class inv_glcm_model")
  }
  X <- .reo_check_matrix(X, "score_inv_glcm", "X")
  .reo_check_meta(meta, nrow(X), "score_inv_glcm", "meta")

  # Present model-universe features in the FROZEN canonical order.
  ordered_present <- intersect(model$order_perm, colnames(X))
  M <- length(ordered_present)
  if (M < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, ordered_present, drop = FALSE]
  head  <- model$head
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    phi <- .inv_glcm_descriptor(X_use[i, ], model$edges, model$K,
                                model$hp$min_features, model$hp$eps)
    if (is.null(phi)) return(0)
    .inv_glcm_head_predict(phi, head)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_inv_glcm: scorer produced non-finite or wrong-length output")
  }
  out
}


# ---------------------------------------------------------------------------
# Hyperparameter resolver (strict allowed-list, exact [[ ]] reads, integer-
# overflow guards; unnamed/duplicate names rejected BEFORE the allowed-list
# check -- mirroring .inv_css_resolve_hp).
# ---------------------------------------------------------------------------
.inv_glcm_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_inv_glcm: hp must be a list")
  if (length(hp) > 0L) {
    nms <- names(hp)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("fit_inv_glcm: all hp entries must be named")
    }
    if (anyDuplicated(nms)) {
      stop("fit_inv_glcm: duplicated hp field(s): ",
           paste(unique(nms[duplicated(nms)]), collapse = ", "))
    }
  }
  allowed <- c("n_levels", "shrink", "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_inv_glcm: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  n_levels <- hp[["n_levels"]]
  if (is.null(n_levels)) n_levels <- 8L
  if (!is.numeric(n_levels) || length(n_levels) != 1L ||
      !is.finite(n_levels) || n_levels > .Machine$integer.max ||
      n_levels != as.integer(n_levels)) {
    stop("fit_inv_glcm: hp$n_levels must be an integer in [2, 32]")
  }
  n_levels <- as.integer(n_levels)
  if (n_levels < 2L || n_levels > 32L) {
    stop("fit_inv_glcm: hp$n_levels must be an integer in [2, 32]")
  }

  shrink <- hp[["shrink"]]
  if (is.null(shrink)) shrink <- 0.1
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0) {
    stop("fit_inv_glcm: hp$shrink must be a non-negative finite number")
  }
  shrink <- as.numeric(shrink)

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 8L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_inv_glcm: hp$min_features must be an integer >= 3")
  }
  min_features <- as.integer(min_features)
  if (min_features < 3L) {
    stop("fit_inv_glcm: hp$min_features must be an integer >= 3")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_inv_glcm: hp$eps must be a positive finite number")
  }
  eps <- as.numeric(eps)

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_inv_glcm: hp$seed must be a non-negative integer")
  }
  seed <- as.integer(seed)

  list(
    n_levels     = n_levels,
    shrink       = shrink,
    min_features = min_features,
    eps          = eps,
    seed         = seed
  )
}


# ---------------------------------------------------------------------------
# Within-sample rCLR landscape for one specimen: ordered values in the frozen
# canonical order (caller passes the sub-matrix already column-ordered). Returns
# a named numeric vector of rCLR values over the POSITIVE support (v > 0).
# Keyed on v > 0, NOT g != 0, to preserve coordinates whose rCLR equals 0 (the
# C1 lesson from inv-css). Scale-invariant: log(c*v) - mean(log(c*v)) = rCLR(v).
# ---------------------------------------------------------------------------
.inv_glcm_rclr_landscape <- function(v_ordered) {
  pos <- which(v_ordered > 0)
  if (length(pos) == 0L) return(numeric(0))
  lv <- log(v_ordered[pos])
  lv - mean(lv)
}


# ---------------------------------------------------------------------------
# Quantise a rCLR landscape vector g to levels {1..K} using frozen edges.
# Right-closed binning: findInterval(g, edges) + 1. With K=1 (edges = numeric(0))
# every value maps to level 1. Levels are always in {1..K} by construction.
# ---------------------------------------------------------------------------
.inv_glcm_quantise <- function(g, edges, K) {
  as.integer(findInterval(g, edges) + 1L)
}


# ---------------------------------------------------------------------------
# 1-D distance-1 GLCM from a quantised sequence q (integer vector in {1..K}).
# For t = 1..(m-1): C[q[t], q[t+1]] += 1. Symmetrised (C + t(C)) and
# normalised to P = C / sum(C). Returns a K x K probability matrix.
# If the sequence has length < 2, no adjacencies exist; returns the K x K zero
# matrix (the descriptor helper handles this via the degeneracy guard).
# ---------------------------------------------------------------------------
.inv_glcm_build_P <- function(q, K) {
  C <- matrix(0, K, K)
  m <- length(q)
  if (m >= 2L) {
    for (t in seq_len(m - 1L)) {
      C[q[t], q[t + 1L]] <- C[q[t], q[t + 1L]] + 1
    }
  }
  C <- C + t(C)
  s <- sum(C)
  if (s > 0) C / s else C
}


# ---------------------------------------------------------------------------
# Five Haralick descriptors from a K x K normalised GLCM probability matrix P.
# Let i, j index levels 1..K. Definitions:
#   contrast    = sum_{ij} (i-j)^2 * P_{ij}
#   homogeneity = sum_{ij} P_{ij} / (1 + (i-j)^2)     [IDM]
#   energy      = sum_{ij} P_{ij}^2                    [ASM]
#   entropy     = - sum_{ij} P_{ij} * log(P_{ij})  (0*log(0) := 0)
#   correlation = (sum_{ij} i*j*P_{ij} - mu_i*mu_j) / (sigma_i * sigma_j)
#                 where mu_i  = sum_i  i * rowSums(P)
#                       mu_j  = sum_j  j * colSums(P)
#                       sig_i = sqrt(sum_i (i-mu_i)^2 * rowSums(P))
#                       sig_j = sqrt(sum_j (j-mu_j)^2 * colSums(P))
#                 If sig_i < eps OR sig_j < eps (degenerate: all mass on one
#                 level, i.e. constant landscape) then correlation := 0.
# Returns a length-5 numeric vector (contrast, homogeneity, energy, entropy,
# correlation). A constant level sequence (q all equal) yields a P with all
# mass on the diagonal block for that level: contrast 0, energy 1 (if one cell
# = 1), entropy 0, correlation 0 by the guard -- all well-defined, no errors.
# ---------------------------------------------------------------------------
.inv_glcm_haralick <- function(P, K, eps) {
  lv <- seq_len(K)
  pi_i <- rowSums(P)   # marginal over rows (levels 1..K)
  pi_j <- colSums(P)   # marginal over cols (levels 1..K)

  # Contrast and homogeneity: depend only on (i-j)^2 weight
  contrast    <- 0
  homogeneity <- 0
  energy      <- sum(P^2)
  entropy     <- 0
  cross_mean  <- 0

  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      p_ij  <- P[i, j]
      d2    <- (i - j)^2
      contrast    <- contrast    + d2 * p_ij
      homogeneity <- homogeneity + p_ij / (1 + d2)
      if (p_ij > 0) entropy <- entropy - p_ij * log(p_ij)
      cross_mean <- cross_mean + i * j * p_ij
    }
  }

  mu_i <- sum(lv * pi_i)
  mu_j <- sum(lv * pi_j)
  sig_i <- sqrt(sum((lv - mu_i)^2 * pi_i))
  sig_j <- sqrt(sum((lv - mu_j)^2 * pi_j))

  if (sig_i < eps || sig_j < eps) {
    correlation <- 0
  } else {
    correlation <- (cross_mean - mu_i * mu_j) / (sig_i * sig_j)
  }

  c(contrast, homogeneity, energy, entropy, correlation)
}


# ---------------------------------------------------------------------------
# Full single-specimen GLCM+Haralick descriptor. Returns NULL when the specimen
# has fewer than min_features positive-support values. Called identically at fit
# (full ordered universe passed) and at score (present ordered subset passed),
# so there is no fit/score drift.
# v_ordered: the specimen's values ALREADY placed in the frozen canonical order
#            (caller selects ordered_present columns before calling).
# ---------------------------------------------------------------------------
.inv_glcm_descriptor <- function(v_ordered, edges, K, min_features, eps) {
  g <- .inv_glcm_rclr_landscape(v_ordered)
  if (length(g) < min_features) return(NULL)
  q <- .inv_glcm_quantise(g, edges, K)
  P <- .inv_glcm_build_P(q, K)
  .inv_glcm_haralick(P, K, eps)
}


# ---------------------------------------------------------------------------
# Frozen standardize + ridge-LDA head -- identical construction to inv-css.
# center/scale from training Phi; sd floored at eps; zero-variance descriptors
# deactivated. Pooled within-class covariance ridged until Cholesky succeeds.
# w = Sinv (mu_case - mu_ctrl), b = -0.5 (mu_case + mu_ctrl).w. Larger = case.
# ---------------------------------------------------------------------------
.inv_glcm_fit_head <- function(Phi, y, shrink, eps) {
  d <- ncol(Phi)
  center  <- colMeans(Phi)
  raw_sd  <- apply(Phi, 2L, stats::sd)
  raw_sd[!is.finite(raw_sd)] <- 0
  active  <- raw_sd > 0
  scale   <- pmax(raw_sd, eps)

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
  Sw <- (Sw + t(Sw)) / 2

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
    stop("fit_inv_glcm: within-class covariance could not be made PD")
  }
  Sinv <- chol2inv(ch)

  delta <- mu_case - mu_ctrl
  w <- as.numeric(Sinv %*% delta)
  w[!active] <- 0
  b <- -0.5 * sum((mu_case + mu_ctrl) * w)

  list(
    center = center,
    scale  = scale,
    active = active,
    w      = w,
    b      = as.numeric(b),
    ridge  = ridge
  )
}

# Apply frozen head: standardise, zero inactive descriptors, return LDA score.
.inv_glcm_head_predict <- function(phi, head) {
  z <- (phi - head$center) / head$scale
  z[!head$active] <- 0
  head$b + sum(head$w * z)
}
