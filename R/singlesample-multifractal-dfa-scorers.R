#' @title Fit MFDFA spectrum (frac-mfdfa) within-sample discriminator
#'
#' @description
#' Learns a frozen, single-sample discriminator from the multifractal detrended
#' fluctuation analysis (MFDFA) spectrum of a specimen's canonical-order
#' robust-CLR (rCLR) abundance landscape. The construction has six stages:
#' \enumerate{
#'   \item \strong{Frozen canonical feature order.} Features are sorted once by
#'     decreasing training column mean, with ties broken by feature index. The
#'     resulting \code{order_perm} is reused unchanged at scoring.
#'   \item \strong{Per-specimen rCLR landscape.} For one specimen, model-universe
#'     features present in \code{X} are placed in the frozen order. Over that
#'     specimen's own positive support \eqn{P=\{j:v_j>0\}}, the within-sample rCLR
#'     is \eqn{g_j=\log v_j-\mathrm{mean}_{k\in P}\log v_k}. Structural zeros are
#'     excluded, but coordinates whose centred-log value is exactly zero are kept.
#'     This landscape is exactly invariant to multiplying the specimen by any
#'     positive scalar.
#'   \item \strong{Fixed-length resample.} MFDFA needs one frozen absolute scale
#'     ladder. Because the positive support size can vary by specimen, the rCLR
#'     landscape is linearly interpolated to \code{hp$resample_len} equally spaced
#'     positions over the native index range. This imposes a fixed low-pass at
#'     the smallest scales, so the default ladder starts at \code{s_min >= 8} and
#'     avoids the smallest windows.
#'   \item \strong{Base-R MFDFA.} The length-\eqn{L} series is profiled as
#'     \eqn{Y(i)=\sum_{k\le i}(x_k-\bar x)}. For each frozen scale, non-overlap
#'     windows are taken from both the start and the end of the profile. A
#'     polynomial trend of order \code{hp$m_poly} is removed by least squares, the
#'     window mean residual square is floored at \code{hp$eps}, and fluctuation
#'     functions \eqn{F_q(s)} are formed for the frozen \code{hp$q} grid. The
#'     generalized Hurst exponent \eqn{h(q)} is the slope of
#'     \eqn{\log F_q(s)} against \eqn{\log s}.
#'   \item \strong{Multifractal descriptors.} The fixed descriptor contains all
#'     \eqn{h(q)} values, spectrum width \code{dAlpha}, \eqn{h(2)}, and a spectrum
#'     asymmetry term
#'     \eqn{\alpha(q_{max})+\alpha(q_{min})-2\alpha((q_{min}+q_{max})/2)}, with
#'     the midpoint alpha obtained by linear interpolation on the frozen q grid.
#'   \item \strong{Frozen standardize + ridge-LDA head.} Descriptor center,
#'     scale, active flags, ridge, weight vector, and intercept are learned from
#'     training descriptors only. The head is oriented so larger scores are more
#'     case-like.
#' }
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp Optional list of frozen MFDFA hyperparameters. Supported names are
#'   validated by the internal resolver: resample_len, n_scales, s_min, s_max,
#'   m_poly, q, shrink, min_features, eps, and seed.
#'
#' @return A plain list of class \code{frac_mfdfa_model} containing
#'   \code{feature_universe}, \code{order_perm}, frozen \code{scales},
#'   \code{q}, \code{descriptor_dim}, \code{descriptor_names}, frozen
#'   \code{head}, and resolved \code{hp}.
#'
#' @details
#' The per-specimen score is the frozen linear predictor
#' \deqn{S(x)=b+\sum_c w_c\{\phi_c(x)-center_c\}/scale_c,}
#' where \eqn{\phi(x)} is the fixed MFDFA descriptor. Every quantity except
#' \eqn{\phi(x)} is frozen at fit time. At scoring, \eqn{\phi(x)} is computed
#' from only that row's present positive values in the frozen feature order; no
#' scored-batch statistic, cross-row coupling, or test-batch normalization is
#' used. The rCLR landscape, fixed-length interpolation, MFDFA fluctuation
#' ratios, descriptors, and final score are invariant to per-sample positive
#' scaling.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 180
#' n <- 160
#' features <- paste0("miR-", sprintf("%03d", seq_len(p)))
#' y <- rep(c(0, 1), each = n / 2)
#' t <- seq(0, 1, length.out = p)
#' X <- matrix(0, n, p, dimnames = list(NULL, features))
#' for (i in seq_len(n)) {
#'   rough <- if (y[i] == 1) 0.8 * sin(48 * pi * t) else 0
#'   mu <- seq(5, 2, length.out = p) + rough
#'   X[i, ] <- exp(stats::rnorm(p, mu, 0.08))
#' }
#' model <- fit_frac_mfdfa(X, y)
#' score <- score_frac_mfdfa(model, X)
#' }
#'
#' @references
#' Kantelhardt JW, Zschiegner SA, Koscielny-Bunde E, et al. (2002)
#' Multifractal detrended fluctuation analysis of nonstationary time series.
#' \emph{Physica A} 316(1-4): 87-114.
#'
#' Peng CK, Buldyrev SV, Havlin S, et al. (1994) Mosaic organization of DNA
#' nucleotides. \emph{Physical Review E} 49(2): 1685-1689.
#'
#' Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse compositional
#' technique reveals microbial perturbations. \emph{mSystems} 4(1): e00016-19.
#'
#' @export
fit_frac_mfdfa <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_frac_mfdfa", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_frac_mfdfa", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_frac_mfdfa")
  hp <- .frac_mfdfa_resolve_hp(hp)

  if (ncol(X_train) < hp$min_features) {
    stop("fit_frac_mfdfa: X_train must contain at least hp$min_features features")
  }

  col_mean <- colMeans(X_train)
  order_idx <- order(-col_mean, seq_len(ncol(X_train)))
  order_perm <- colnames(X_train)[order_idx]

  Xo <- X_train[, order_perm, drop = FALSE]
  desc_list <- lapply(seq_len(nrow(Xo)), function(i) {
    .frac_mfdfa_descriptor(Xo[i, ], hp)
  })
  keep <- !vapply(desc_list, is.null, logical(1))
  if (sum(keep) < 2L) {
    stop("fit_frac_mfdfa: too few training specimens have >= hp$min_features ",
         "positive-support coordinates")
  }
  Phi <- do.call(rbind, desc_list[keep])
  yk <- y[keep]
  if (sum(yk == 1L) < 1L || sum(yk == 0L) < 1L) {
    stop("fit_frac_mfdfa: after the min_features filter both classes must ",
         "retain >= 1 specimen")
  }

  head <- .frac_mfdfa_fit_head(Phi, yk, hp$shrink, hp$eps)

  model <- list(
    feature_universe = colnames(X_train),
    order_perm = order_perm,
    scales = hp$scales,
    q = hp$q,
    descriptor_dim = ncol(Phi),
    descriptor_names = colnames(Phi),
    head = head,
    hp = hp
  )
  class(model) <- "frac_mfdfa_model"
  model
}


#' @title Score MFDFA spectrum (frac-mfdfa) within-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_frac_mfdfa}}. Model-universe features present in \code{X} are
#' placed in the frozen canonical order; for one specimen its within-sample rCLR
#' landscape is resampled to the frozen length, transformed into a fixed MFDFA
#' spectrum descriptor using the frozen scales and q grid, standardized by the
#' frozen training center/scale, and passed through the frozen ridge-LDA head.
#' Larger scores are more case-like.
#'
#' Scoring uses only each row's own values plus the frozen model. If fewer than
#' \code{model$hp$min_features} model-universe features are present in \code{X},
#' or a specimen's own positive support is below that floor, the documented
#' neutral score \code{0} is returned.
#'
#' @param model A \code{frac_mfdfa_model} returned by
#'   \code{\link{fit_frac_mfdfa}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like.
#'
#' @examples
#' \dontrun{
#' score_frac_mfdfa(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Kantelhardt JW, Zschiegner SA, Koscielny-Bunde E, et al. (2002)
#' Multifractal detrended fluctuation analysis of nonstationary time series.
#' \emph{Physica A} 316(1-4): 87-114.
#'
#' Peng CK, Buldyrev SV, Havlin S, et al. (1994) Mosaic organization of DNA
#' nucleotides. \emph{Physical Review E} 49(2): 1685-1689.
#'
#' @export
score_frac_mfdfa <- function(model, X, meta = NULL) {
  if (!inherits(model, "frac_mfdfa_model")) {
    stop("score_frac_mfdfa: model must have class frac_mfdfa_model")
  }
  X <- .reo_check_matrix(X, "score_frac_mfdfa", "X")
  .reo_check_meta(meta, nrow(X), "score_frac_mfdfa", "meta")

  ordered_present <- intersect(model$order_perm, colnames(X))
  if (length(ordered_present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, ordered_present, drop = FALSE]
  out <- vapply(seq_len(nrow(X_use)), function(i) {
    phi <- .frac_mfdfa_descriptor(X_use[i, ], model$hp)
    if (is.null(phi)) return(0)
    .frac_mfdfa_head_predict(phi, model$head)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_frac_mfdfa: scorer produced non-finite or wrong-length output")
  }
  out
}


.frac_mfdfa_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_frac_mfdfa: hp must be a list")
  if (length(hp) > 0L) {
    nms <- names(hp)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("fit_frac_mfdfa: all hp entries must be named")
    }
    if (anyDuplicated(nms)) {
      stop("fit_frac_mfdfa: duplicated hp field(s): ",
           paste(unique(nms[duplicated(nms)]), collapse = ", "))
    }
  }
  allowed <- c("resample_len", "n_scales", "s_min", "s_max", "m_poly", "q",
               "shrink", "min_features", "eps", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_frac_mfdfa: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  resample_len <- hp[["resample_len"]]
  if (is.null(resample_len)) resample_len <- 256L
  if (!is.numeric(resample_len) || length(resample_len) != 1L ||
      !is.finite(resample_len) || resample_len > .Machine$integer.max ||
      resample_len != as.integer(resample_len)) {
    stop("fit_frac_mfdfa: hp$resample_len must be an integer >= 32")
  }
  resample_len <- as.integer(resample_len)
  if (resample_len < 32L) {
    stop("fit_frac_mfdfa: hp$resample_len must be an integer >= 32")
  }

  n_scales <- hp[["n_scales"]]
  if (is.null(n_scales)) n_scales <- 12L
  if (!is.numeric(n_scales) || length(n_scales) != 1L ||
      !is.finite(n_scales) || n_scales > .Machine$integer.max ||
      n_scales != as.integer(n_scales)) {
    stop("fit_frac_mfdfa: hp$n_scales must be an integer >= 4")
  }
  n_scales <- as.integer(n_scales)
  if (n_scales < 4L) {
    stop("fit_frac_mfdfa: hp$n_scales must be an integer >= 4")
  }

  s_min <- hp[["s_min"]]
  if (is.null(s_min)) s_min <- 8L
  if (!is.numeric(s_min) || length(s_min) != 1L || !is.finite(s_min) ||
      s_min > .Machine$integer.max || s_min != as.integer(s_min)) {
    stop("fit_frac_mfdfa: hp$s_min must be an integer >= 4")
  }
  s_min <- as.integer(s_min)
  if (s_min < 4L) stop("fit_frac_mfdfa: hp$s_min must be an integer >= 4")

  s_max <- hp[["s_max"]]
  if (is.null(s_max)) s_max <- resample_len %/% 4L
  if (!is.numeric(s_max) || length(s_max) != 1L || !is.finite(s_max) ||
      s_max > .Machine$integer.max || s_max != as.integer(s_max)) {
    stop("fit_frac_mfdfa: hp$s_max must be an integer > hp$s_min and <= L/2")
  }
  s_max <- as.integer(s_max)
  if (s_max <= s_min || s_max > resample_len %/% 2L) {
    stop("fit_frac_mfdfa: hp$s_max must be > hp$s_min and <= resample_len %/% 2")
  }

  m_poly <- hp[["m_poly"]]
  if (is.null(m_poly)) m_poly <- 1L
  if (!is.numeric(m_poly) || length(m_poly) != 1L || !is.finite(m_poly) ||
      m_poly > .Machine$integer.max || m_poly != as.integer(m_poly)) {
    stop("fit_frac_mfdfa: hp$m_poly must be an integer in [1, 3]")
  }
  m_poly <- as.integer(m_poly)
  if (m_poly < 1L || m_poly > 3L) {
    stop("fit_frac_mfdfa: hp$m_poly must be an integer in [1, 3]")
  }
  # The smallest MFDFA window has s_min points and is detrended by an order-m_poly
  # polynomial (m_poly + 1 coefficients). Require s_min >= m_poly + 2 so that
  # window is OVER-determined: at s_min == m_poly + 1 the detrend is an exact fit,
  # forcing F2 = 0 -> the eps floor -> a constant-floored smallest scale that
  # biases the log F_q(s) vs log s slope. The defaults (s_min = 8, m_poly = 1) are
  # far from this corner; this guard only rejects the degenerate combination.
  if (s_min < m_poly + 2L) {
    stop("fit_frac_mfdfa: hp$s_min must be >= hp$m_poly + 2 so the smallest ",
         "MFDFA window is not exactly fit by the detrending polynomial")
  }

  q <- hp[["q"]]
  if (is.null(q)) q <- c(-5, -3, -1, 2, 3, 5)
  if (!is.numeric(q) || length(q) < 2L || any(!is.finite(q))) {
    stop("fit_frac_mfdfa: hp$q must be a finite numeric vector of length >= 2")
  }
  q <- sort(as.numeric(q))
  if (anyDuplicated(q)) stop("fit_frac_mfdfa: hp$q must not contain duplicates")
  if (!any(q == 2)) stop("fit_frac_mfdfa: hp$q must include 2")

  shrink <- hp[["shrink"]]
  if (is.null(shrink)) shrink <- 0.1
  if (!is.numeric(shrink) || length(shrink) != 1L || !is.finite(shrink) ||
      shrink < 0) {
    stop("fit_frac_mfdfa: hp$shrink must be a non-negative finite number")
  }
  shrink <- as.numeric(shrink)

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- max(32L, 2L * s_max)
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_frac_mfdfa: hp$min_features must be an integer >= 2 * hp$s_max")
  }
  min_features <- as.integer(min_features)
  if (min_features < 2L * s_max) {
    stop("fit_frac_mfdfa: hp$min_features must be >= 2 * hp$s_max")
  }

  eps <- hp[["eps"]]
  if (is.null(eps)) eps <- 1e-8
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_frac_mfdfa: hp$eps must be a positive finite number")
  }
  eps <- as.numeric(eps)

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_frac_mfdfa: hp$seed must be a non-negative integer")
  }
  seed <- as.integer(seed)

  scales <- .frac_mfdfa_scale_ladder(resample_len, s_min, s_max, n_scales)
  if (length(scales) < 4L) {
    stop("fit_frac_mfdfa: frozen scale ladder must contain at least 4 scales")
  }

  list(
    resample_len = resample_len,
    n_scales = n_scales,
    s_min = s_min,
    s_max = s_max,
    m_poly = m_poly,
    q = q,
    scales = scales,
    shrink = shrink,
    min_features = min_features,
    eps = eps,
    seed = seed
  )
}


.frac_mfdfa_scale_ladder <- function(L, s_min, s_max, n_scales) {
  raw <- exp(seq(log(s_min), log(s_max), length.out = n_scales))
  sort(unique(as.integer(round(raw))))
}


.frac_mfdfa_rclr_landscape <- function(v_ordered) {
  pos <- which(v_ordered > 0)
  if (length(pos) == 0L) return(numeric(0))
  lv <- log(v_ordered[pos])
  lv - mean(lv)
}


.frac_mfdfa_resample <- function(g, L) {
  m <- length(g)
  stats::approx(
    x = seq_len(m),
    y = as.numeric(g),
    xout = seq(1, m, length.out = L),
    ties = "ordered",
    rule = 2
  )$y
}


.frac_mfdfa_descriptor <- function(v_ordered, hp) {
  g <- .frac_mfdfa_rclr_landscape(v_ordered)
  if (length(g) < hp$min_features) return(NULL)
  x <- .frac_mfdfa_resample(g, hp$resample_len)
  h <- .frac_mfdfa_hq(x, hp$scales, hp$q, hp$m_poly, hp$eps)
  .frac_mfdfa_spectrum_descriptor(h, hp$q)
}


.frac_mfdfa_hq <- function(x, scales, q, m_poly, eps) {
  Y <- cumsum(x - mean(x))
  log_s <- log(scales)
  h <- numeric(length(q))
  for (qi in seq_along(q)) {
    log_fq <- vapply(scales, function(s) {
      F2 <- .frac_mfdfa_segment_f2(Y, s, m_poly)
      F2 <- pmax(F2, eps)
      if (q[qi] == 0) {
        0.5 * mean(log(F2))
      } else {
        a <- (q[qi] / 2) * log(F2)
        amax <- max(a)
        (log(mean(exp(a - amax))) + amax) / q[qi]
      }
    }, numeric(1))
    fit <- stats::lm.fit(cbind(1, log_s), log_fq)
    h[qi] <- fit$coefficients[2L]
  }
  names(h) <- paste0("h_q", .frac_mfdfa_q_label(q))
  h
}


.frac_mfdfa_segment_f2 <- function(Y, s, m_poly) {
  L <- length(Y)
  Ns <- floor(L / s)
  if (Ns < 1L) return(numeric(0))
  design <- .frac_mfdfa_design(s, m_poly)
  out <- numeric(2L * Ns)
  k <- 1L
  for (nu in seq_len(Ns)) {
    idx <- ((nu - 1L) * s + 1L):(nu * s)
    out[k] <- .frac_mfdfa_window_f2(Y[idx], design)
    k <- k + 1L
  }
  for (nu in seq_len(Ns)) {
    idx <- (L - nu * s + 1L):(L - (nu - 1L) * s)
    out[k] <- .frac_mfdfa_window_f2(Y[idx], design)
    k <- k + 1L
  }
  out
}


.frac_mfdfa_design <- function(s, m_poly) {
  u <- seq(-1, 1, length.out = s)
  outer(u, seq.int(0L, m_poly), "^")
}


.frac_mfdfa_window_f2 <- function(y, design) {
  fit <- stats::lm.fit(design, y)
  mean(fit$residuals^2)
}


.frac_mfdfa_spectrum_descriptor <- function(h, q) {
  tau <- q * h - 1
  alpha <- .frac_mfdfa_alpha(q, tau)
  d_alpha <- max(alpha) - min(alpha)
  h2 <- h[which(q == 2)[1L]]
  q_mid <- (min(q) + max(q)) / 2
  alpha_mid <- stats::approx(q, alpha, xout = q_mid, ties = "ordered")$y
  asym <- alpha[length(alpha)] + alpha[1L] - 2 * alpha_mid
  out <- c(as.numeric(h), dAlpha = d_alpha, h2 = as.numeric(h2),
           asymmetry = as.numeric(asym))
  names(out)[seq_along(q)] <- paste0("h_q", .frac_mfdfa_q_label(q))
  out
}


.frac_mfdfa_alpha <- function(q, tau) {
  n <- length(q)
  alpha <- numeric(n)
  if (n == 2L) {
    alpha[] <- (tau[2L] - tau[1L]) / (q[2L] - q[1L])
    return(alpha)
  }
  alpha[1L] <- (tau[2L] - tau[1L]) / (q[2L] - q[1L])
  alpha[n] <- (tau[n] - tau[n - 1L]) / (q[n] - q[n - 1L])
  for (i in 2L:(n - 1L)) {
    alpha[i] <- (tau[i + 1L] - tau[i - 1L]) / (q[i + 1L] - q[i - 1L])
  }
  alpha
}


.frac_mfdfa_q_label <- function(q) {
  lab <- format(q, scientific = FALSE, trim = TRUE)
  lab <- gsub("-", "m", lab, fixed = TRUE)
  gsub("\\.", "p", lab)
}


.frac_mfdfa_fit_head <- function(Phi, y, shrink, eps) {
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
    stop("fit_frac_mfdfa: within-class covariance could not be made positive ",
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


.frac_mfdfa_head_predict <- function(phi, head) {
  z <- (phi - head$center) / head$scale
  z[!head$active] <- 0
  head$b + sum(head$w * z)
}
