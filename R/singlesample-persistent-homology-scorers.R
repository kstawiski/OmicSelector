#' @title Fit topological persistence-image (tda-ph) within-sample discriminator
#'
#' @description
#' Learns a frozen, single-sample discriminator from the 0-dimensional
#' sublevel-set persistent homology of a specimen's within-sample rCLR landscape.
#' The construction has five stages, all estimated from training data only (or,
#' like the persistence-image Gaussian integral stencil, fixed data-free
#' constants):
#' \enumerate{
#'   \item \strong{Canonical feature order.} Features are arranged once, in
#'     \emph{decreasing} training column mean of \code{X_train}; ties are broken
#'     by feature index, making \code{order_perm} deterministic and
#'     locale-independent. This fixed order defines the path graph on which every
#'     specimen's landscape is read.
#'   \item \strong{Single-specimen rCLR landscape.} For one specimen, the
#'     model-universe features present in \code{X} are taken in the frozen order
#'     (\code{intersect(order_perm, colnames(X))}). Over that specimen's own
#'     positive support \eqn{P=\{j: v_j>0\}}, the landscape is
#'     \eqn{g_j=\log v_j-\mathrm{mean}_{k\in P}\log v_k}, keeping the frozen
#'     order and dropping only structural zeros. A coordinate with
#'     \eqn{g_j=0} is kept; exclusion keys on \eqn{v_j=0}, never on
#'     \eqn{g_j=0}. This rCLR landscape is exactly invariant to per-sample
#'     positive scaling.
#'   \item \strong{Path lower-star persistence.} The 0-dimensional finite
#'     sublevel-set classes of \eqn{g} on the path graph are computed by a
#'     deterministic union-find merge tree. Vertices enter in increasing
#'     \eqn{g}; ties are broken by ascending position. When two existing
#'     components merge, the younger component (higher birth, ties by later
#'     birth position) dies at the saddle value. The essential global-minimum
#'     class is excluded by design.
#'   \item \strong{Frozen persistence image.} Finite classes are mapped to
#'     birth-persistence coordinates and rasterised onto a frozen
#'     \eqn{R\times R} grid by the exact Gaussian pixel integral using
#'     \code{stats::pnorm}. The weight is \eqn{w(p)=p}. Vectorisation is
#'     row-major: pixel \eqn{(a,c)} (birth pixel \eqn{a}, persistence pixel
#'     \eqn{c}) is stored at \code{(a - 1) * R + c}.
#'   \item \strong{Frozen standardize + LDA head.} Training images are
#'     standardised by frozen per-pixel \code{center}/\code{scale}; exactly
#'     zero-variance pixels are inactive and contribute \code{0}. A
#'     ridge-regularised LDA head is fit on active standardised pixels, oriented
#'     so larger scores are more case-like.
#' }
#' The key topological design decision is that only \emph{finite} 0-dimensional
#' classes are used. The essential global-minimum class is omitted because it
#' primarily re-encodes global rCLR spread; finite classes isolate secondary
#' valley structure in the canonical landscape. A strictly monotone landscape has
#' no finite classes and gives an all-zero persistence image. That all-zero image
#' is \emph{not} short-circuited: it is standardised and passed through the frozen
#' head exactly as every other descriptor in this family, so a monotone specimen
#' receives the head's learned, constant linear readout for the \dQuote{no
#' secondary valleys} descriptor -- a single deterministic baseline shared by
#' every empty-diagram specimen. This is the intended behaviour (the absence of
#' finite valleys is itself a discriminative cue, not an abstention); it is
#' \emph{not} the bare intercept \eqn{b}, and it is distinct from the hard neutral
#' \code{0} that \code{\link{score_tda_ph}} returns only when a specimen falls
#' below the \code{min_features} floor (a genuine abstention on insufficient
#' support). In the degenerate case where \emph{every} training specimen is
#' monotone, all training images are zero, every descriptor is deactivated, and
#' the head collapses to the constant-\code{0} scorer (\code{w = 0}, \code{b =
#' 0}).
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters (exact-name reads, strict
#'   allowed-list): \code{grid_res} persistence-image side length \eqn{R}
#'   (default \code{8L}; integer in \code{[3, 16]}), \code{sigma_frac} Gaussian
#'   bandwidth as a multiple of the smaller pixel width (default \code{1.0};
#'   positive finite), \code{grid_pad} fractional padding of frozen grid ranges
#'   (default \code{0.05}; non-negative finite), \code{shrink} ridge fraction of
#'   the mean covariance diagonal (default \code{0.1}; non-negative finite),
#'   \code{min_features} per-specimen positive-support floor (default \code{8L};
#'   integer \eqn{\ge 3}), \code{eps} positive standard-deviation / ridge floor
#'   (default \code{1e-6}), and \code{seed} (default \code{1L}; accepted for
#'   interface uniformity only -- no randomness is used).
#'
#' @return A plain list of class \code{tda_ph_model} containing
#'   \code{feature_universe}, the frozen canonical \code{order_perm}, frozen
#'   persistence-image \code{grid}, \code{descriptor_dim} (\eqn{= R^2}), frozen
#'   ridge-LDA \code{head} (\code{center}, \code{scale}, \code{active},
#'   \code{w}, \code{b}, \code{ridge}), and the resolved \code{hp}.
#'
#' @details
#' The per-specimen score (see \code{\link{score_tda_ph}}) is the frozen linear
#' predictor on the standardised persistence image,
#' \deqn{S(x) = b + \sum_c w_c \frac{\phi_c(x)-\mathrm{center}_c}
#'                                      {\mathrm{scale}_c},}
#' where \eqn{\phi(x)} is the specimen's finite-class persistence image. Every
#' quantity except \eqn{\phi(x)} is frozen at fit time, and \eqn{\phi(x)} depends
#' only on \eqn{x}'s own rCLR landscape, its own 0-dimensional finite diagram,
#' and the frozen image grid. The score is therefore computable from one
#' specimen alone, with no test-batch normalisation, no cross-row coupling, and
#' exact invariance to per-sample positive scaling. A specimen whose finite
#' (birth, persistence) coordinates fall outside the frozen training grid has part
#' of its persistence-image mass truncated by the grid edges; the score is still
#' deterministic, finite, and rank/AUC-valid, but its raw magnitude is not
#' calibrated across that grid boundary.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 72
#' n <- 160
#' features <- paste0("miR-", sprintf("%03d", seq_len(p)))
#' y <- rep(c(0, 1), each = n / 2)
#' X <- matrix(0, n, p, dimnames = list(NULL, features))
#' base_mu <- seq(10, 2, length.out = p)
#' block <- 24:48
#' pattern <- rep(c(0.45, -0.75, 0.35, -0.65), length.out = length(block))
#' for (i in seq_len(n)) {
#'   mu <- base_mu
#'   if (y[i] == 1) mu[block] <- mu[block] + pattern
#'   X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.02))
#' }
#' model <- fit_tda_ph(X, y)
#' score <- score_tda_ph(model, X)
#' }
#'
#' @references
#' Edelsbrunner H, Letscher D, Zomorodian A. (2002) Topological persistence and
#' simplification. \emph{Discrete & Computational Geometry} 28(4): 511-533.
#'
#' Adams H, Emerson T, Kirby M, et al. (2017) Persistence images: a stable
#' vector representation of persistent homology. \emph{Journal of Machine
#' Learning Research} 18(8): 1-35.
#'
#' Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse compositional
#' technique reveals microbial perturbations. \emph{mSystems} 4(1): e00016-19.
#'
#' @export
fit_tda_ph <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_tda_ph", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_tda_ph", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_tda_ph")
  hp <- .tda_ph_resolve_hp(hp)

  if (ncol(X_train) < hp$min_features) {
    stop("fit_tda_ph: X_train must contain at least hp$min_features features")
  }

  col_mean <- colMeans(X_train)
  order_idx <- order(-col_mean, seq_len(ncol(X_train)))
  order_perm <- colnames(X_train)[order_idx]
  Xo <- X_train[, order_perm, drop = FALSE]

  landscape_list <- lapply(seq_len(nrow(Xo)), function(i) {
    .tda_ph_rclr_landscape(Xo[i, ])
  })
  keep <- vapply(landscape_list, length, integer(1)) >= hp$min_features
  if (sum(keep) < 2L) {
    stop("fit_tda_ph: too few training specimens have >= hp$min_features ",
         "present rCLR coordinates")
  }
  yk <- y[keep]
  if (sum(yk == 1L) < 1L || sum(yk == 0L) < 1L) {
    stop("fit_tda_ph: after the min_features filter both classes must retain ",
         ">= 1 specimen")
  }

  diagram_list <- lapply(landscape_list[keep], .tda_ph_diagram_0d)
  grid <- .tda_ph_fit_grid(diagram_list, hp)
  descriptor_dim <- hp$grid_res * hp$grid_res
  Phi <- t(vapply(diagram_list, .tda_ph_persistence_image,
                  numeric(descriptor_dim), grid = grid))

  head <- .tda_ph_fit_head(Phi, yk, hp$shrink, hp$eps)

  model <- list(
    feature_universe = colnames(X_train),
    order_perm = order_perm,
    grid = grid,
    descriptor_dim = descriptor_dim,
    head = head,
    hp = hp
  )
  class(model) <- "tda_ph_model"
  model
}


#' @title Score topological persistence-image (tda-ph) within-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_tda_ph}}. The model-universe features present in \code{X} are
#' taken in the frozen canonical order, a per-specimen rCLR landscape is built
#' over that row's own positive support, finite 0-dimensional lower-star classes
#' are rasterised through the frozen exact-integral persistence-image grid, and
#' the resulting vector is passed through the frozen standardisation and
#' ridge-LDA head. Larger values are more case-like.
#'
#' Scoring uses only each row's own values plus the frozen model -- no test-batch
#' renormalisation, no cross-row coupling, no statistic estimated from \code{X}.
#' If fewer than \code{model$hp$min_features} model-universe features are
#' present, the documented neutral score \code{0} is returned for every row; a
#' specimen whose own positive-support count is below \code{min_features}
#' likewise scores the neutral \code{0}. A specimen with \emph{sufficient} support
#' but an empty finite diagram (e.g. a strictly monotone landscape) is \emph{not}
#' abstained: its all-zero persistence image is standardised and scored by the
#' frozen head, yielding the constant empty-diagram baseline described in
#' \code{\link{fit_tda_ph}} (which is generally not the bare intercept \eqn{b}).
#'
#' @param model A \code{tda_ph_model} object returned by
#'   \code{\link{fit_tda_ph}}.
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
#' X <- matrix(stats::rgamma(80 * 40, shape = 10, rate = 1), nrow = 80,
#'             dimnames = list(NULL, paste0("miR-", seq_len(40))))
#' y <- rep(c(0, 1), each = 40)
#' X[y == 1, 15:24] <- X[y == 1, 15:24] * rep(c(1.8, 0.55), 5)
#' model <- fit_tda_ph(X, y)
#' score_tda_ph(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Edelsbrunner H, Letscher D, Zomorodian A. (2002) Topological persistence and
#' simplification. \emph{Discrete & Computational Geometry} 28(4): 511-533.
#'
#' Adams H, Emerson T, Kirby M, et al. (2017) Persistence images: a stable
#' vector representation of persistent homology. \emph{Journal of Machine
#' Learning Research} 18(8): 1-35.
#'
#' @export
score_tda_ph <- function(model, X, meta = NULL) {
  if (!inherits(model, "tda_ph_model")) {
    stop("score_tda_ph: model must have class tda_ph_model")
  }
  X <- .reo_check_matrix(X, "score_tda_ph", "X")
  .reo_check_meta(meta, nrow(X), "score_tda_ph", "meta")

  ordered_present <- intersect(model$order_perm, colnames(X))
  if (length(ordered_present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, ordered_present, drop = FALSE]
  head <- model$head
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    phi <- .tda_ph_descriptor(X_use[i, ], model$grid, model$hp$min_features)
    if (is.null(phi)) return(0)
    .tda_ph_head_predict(phi, head)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_tda_ph: scorer produced non-finite or wrong-length output")
  }
  out
}


# ---------------------------------------------------------------------------
# Hyperparameter resolver (strict allowed-list, exact `[[ ]]` reads, per-field
# validation with integer-overflow guards). No randomness is consumed.
# ---------------------------------------------------------------------------
.tda_ph_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_tda_ph: hp must be a list")
  if (length(hp) > 0L) {
    nms <- names(hp)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("fit_tda_ph: all hp entries must be named")
    }
    if (anyDuplicated(nms)) {
      stop("fit_tda_ph: duplicated hp field(s): ",
           paste(unique(nms[duplicated(nms)]), collapse = ", "))
    }
  }
  allowed <- c("grid_res", "sigma_frac", "grid_pad", "shrink",
               "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_tda_ph: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  grid_res <- hp[["grid_res"]]
  if (is.null(grid_res)) grid_res <- 8L
  if (!is.numeric(grid_res) || length(grid_res) != 1L ||
      !is.finite(grid_res) || grid_res > .Machine$integer.max ||
      grid_res != as.integer(grid_res)) {
    stop("fit_tda_ph: hp$grid_res must be an integer in [3, 16]")
  }
  grid_res <- as.integer(grid_res)
  if (grid_res < 3L || grid_res > 16L) {
    stop("fit_tda_ph: hp$grid_res must be an integer in [3, 16]")
  }

  sigma_frac <- hp[["sigma_frac"]]
  if (is.null(sigma_frac)) sigma_frac <- 1.0
  if (!is.numeric(sigma_frac) || length(sigma_frac) != 1L ||
      !is.finite(sigma_frac) || sigma_frac <= 0) {
    stop("fit_tda_ph: hp$sigma_frac must be a positive finite number")
  }
  sigma_frac <- as.numeric(sigma_frac)

  grid_pad <- hp[["grid_pad"]]
  if (is.null(grid_pad)) grid_pad <- 0.05
  if (!is.numeric(grid_pad) || length(grid_pad) != 1L ||
      !is.finite(grid_pad) || grid_pad < 0) {
    stop("fit_tda_ph: hp$grid_pad must be a non-negative finite number")
  }
  grid_pad <- as.numeric(grid_pad)

  shrink <- hp[["shrink"]]
  if (is.null(shrink)) shrink <- 0.1
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0) {
    stop("fit_tda_ph: hp$shrink must be a non-negative finite number")
  }
  shrink <- as.numeric(shrink)

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 8L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_tda_ph: hp$min_features must be an integer >= 3")
  }
  min_features <- as.integer(min_features)
  if (min_features < 3L) {
    stop("fit_tda_ph: hp$min_features must be an integer >= 3")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_tda_ph: hp$eps must be a positive finite number")
  }
  eps <- as.numeric(eps)

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_tda_ph: hp$seed must be a non-negative integer")
  }
  seed <- as.integer(seed)

  list(
    grid_res = grid_res,
    sigma_frac = sigma_frac,
    grid_pad = grid_pad,
    shrink = shrink,
    min_features = min_features,
    eps = eps,
    seed = seed
  )
}


# ---------------------------------------------------------------------------
# Single-specimen within-sample robust-CLR over the specimen's own positive
# support, retaining the existing vector order. Only STRUCTURAL zeros (v == 0)
# are excluded; coordinates with centred log value exactly 0 are kept.
# ---------------------------------------------------------------------------
.tda_ph_rclr_landscape <- function(v) {
  pos <- which(v > 0)
  if (length(pos) == 0L) return(numeric(0))
  lv <- log(v[pos])
  lv - mean(lv)
}


.tda_ph_empty_diagram <- function() {
  matrix(numeric(0), ncol = 2L,
         dimnames = list(NULL, c("birth", "death")))
}


# ---------------------------------------------------------------------------
# 0-dimensional lower-star persistence on a path. Vertices are processed by
# increasing landscape value, ties by ascending position. The essential class
# (global minimum) is deliberately excluded; only finite birth/death pairs with
# positive persistence are returned.
# ---------------------------------------------------------------------------
.tda_ph_diagram_0d <- function(g) {
  m <- length(g)
  if (m < 2L) return(.tda_ph_empty_diagram())

  ord <- order(g, seq_len(m))
  parent <- seq_len(m)
  birth <- rep(NA_real_, m)
  birth_pos <- rep(NA_integer_, m)
  added <- rep(FALSE, m)
  births <- numeric(0)
  deaths <- numeric(0)

  find_root <- function(x) {
    while (parent[x] != x) x <- parent[x]
    x
  }

  is_older <- function(a, b) {
    birth[a] < birth[b] ||
      (birth[a] == birth[b] && birth_pos[a] < birth_pos[b])
  }

  for (i in ord) {
    added[i] <- TRUE
    parent[i] <- i
    birth[i] <- g[i]
    birth_pos[i] <- i

    nbrs <- integer(0)
    if (i > 1L && added[i - 1L]) nbrs <- c(nbrs, i - 1L)
    if (i < m && added[i + 1L]) nbrs <- c(nbrs, i + 1L)
    if (length(nbrs) == 0L) next

    root_current <- find_root(nbrs[1L])
    parent[i] <- root_current
    if (length(nbrs) == 1L) next

    for (nb in nbrs[-1L]) {
      root_other <- find_root(nb)
      root_current <- find_root(root_current)
      if (root_other == root_current) next

      if (is_older(root_current, root_other)) {
        elder <- root_current
        younger <- root_other
      } else {
        elder <- root_other
        younger <- root_current
      }
      if (g[i] > birth[younger]) {
        births <- c(births, birth[younger])
        deaths <- c(deaths, g[i])
      }
      parent[younger] <- elder
      root_current <- elder
    }
  }

  if (length(births) == 0L) return(.tda_ph_empty_diagram())
  out <- cbind(birth = births, death = deaths)
  keep <- is.finite(out[, "birth"]) & is.finite(out[, "death"]) &
    out[, "death"] > out[, "birth"]
  if (!any(keep)) return(.tda_ph_empty_diagram())
  out[keep, , drop = FALSE]
}


.tda_ph_birth_persistence <- function(diagram) {
  if (nrow(diagram) == 0L) {
    return(matrix(numeric(0), ncol = 2L,
                  dimnames = list(NULL, c("birth", "persistence"))))
  }
  p <- diagram[, "death"] - diagram[, "birth"]
  out <- cbind(birth = diagram[, "birth"], persistence = p)
  keep <- is.finite(out[, "birth"]) & is.finite(out[, "persistence"]) &
    out[, "persistence"] > 0
  out[keep, , drop = FALSE]
}


# ---------------------------------------------------------------------------
# Frozen persistence-image grid learned from TRAINING finite diagrams only.
# If training contains no finite classes, a deterministic non-degenerate fallback
# grid is used; all training descriptors are then all-zero and the ridge-LDA head
# becomes the neutral intercept-only direction.
# ---------------------------------------------------------------------------
.tda_ph_fit_grid <- function(diagram_list, hp) {
  bp_list <- lapply(diagram_list, .tda_ph_birth_persistence)
  points <- if (length(bp_list) > 0L) do.call(rbind, bp_list) else NULL
  if (is.null(points) || nrow(points) == 0L) {
    b_center <- 0
    span <- 1
    p_max <- 0
  } else {
    b_lo_raw <- min(points[, "birth"])
    b_hi_raw <- max(points[, "birth"])
    if (b_hi_raw == b_lo_raw) {
      b_center <- b_lo_raw
      span <- 1
    } else {
      b_center <- (b_lo_raw + b_hi_raw) / 2
      span <- b_hi_raw - b_lo_raw
    }
    p_max <- max(points[, "persistence"])
  }

  b_lo <- b_center - span / 2 - hp$grid_pad * span
  b_hi <- b_center + span / 2 + hp$grid_pad * span
  p_hi <- max(hp$eps, (1 + hp$grid_pad) * p_max)
  R <- hp$grid_res
  birth_edges <- seq(b_lo, b_hi, length.out = R + 1L)
  pers_edges <- seq(0, p_hi, length.out = R + 1L)
  wb <- (b_hi - b_lo) / R
  wp <- p_hi / R
  sigma <- hp$sigma_frac * min(wb, wp)
  if (!is.finite(sigma) || sigma <= 0) {
    stop("fit_tda_ph: persistence-image bandwidth is non-positive")
  }

  list(
    R = R,
    b_lo = b_lo,
    b_hi = b_hi,
    p_hi = p_hi,
    wb = wb,
    wp = wp,
    sigma = sigma,
    birth_edges = birth_edges,
    pers_edges = pers_edges,
    weight = "linear_persistence",
    vector_order = "row-major birth pixel, then persistence pixel"
  )
}


# Exact pixel-integral persistence image for one finite diagram, vectorised
# row-major: index (a, c) -> (a - 1) * R + c.
.tda_ph_persistence_image <- function(diagram, grid) {
  R <- grid$R
  phi <- numeric(R * R)
  bp <- .tda_ph_birth_persistence(diagram)
  if (nrow(bp) == 0L) return(phi)

  for (k in seq_len(nrow(bp))) {
    b <- bp[k, "birth"]
    p <- bp[k, "persistence"]
    db <- diff(stats::pnorm((grid$birth_edges - b) / grid$sigma))
    dp <- diff(stats::pnorm((grid$pers_edges - p) / grid$sigma))
    for (a in seq_len(R)) {
      idx <- ((a - 1L) * R + 1L):(a * R)
      phi[idx] <- phi[idx] + p * db[a] * dp
    }
  }
  phi
}


# Full single-specimen descriptor. Returns NULL below the positive-support floor.
# Used at fit (full ordered universe) and at score (ordered present subset).
.tda_ph_descriptor <- function(v, grid, min_features) {
  g <- .tda_ph_rclr_landscape(v)
  if (length(g) < min_features) return(NULL)
  diagram <- .tda_ph_diagram_0d(g)
  .tda_ph_persistence_image(diagram, grid)
}


# ---------------------------------------------------------------------------
# Frozen standardize + ridge-LDA head over training persistence images.
# ---------------------------------------------------------------------------
.tda_ph_fit_head <- function(Phi, y, shrink, eps) {
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
    stop("fit_tda_ph: within-class covariance could not be made positive ",
         "definite")
  }
  Sinv <- chol2inv(ch)

  delta <- mu_case - mu_ctrl
  w <- as.numeric(Sinv %*% delta)
  w[!active] <- 0
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


.tda_ph_head_predict <- function(phi, head) {
  z <- (phi - head$center) / head$scale
  z[!head$active] <- 0
  head$b + sum(head$w * z)
}
