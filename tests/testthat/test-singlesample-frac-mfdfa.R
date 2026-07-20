library(testthat)

.frac_hp <- function(...) {
  hp <- list(
    resample_len = 128L,
    n_scales = 8L,
    s_min = 8L,
    s_max = 32L,
    min_features = 64L
  )
  utils::modifyList(hp, list(...))
}

.make_frac_mfdfa_data <- function(n = 160L, p = 180L, seed = 71L) {
  set.seed(seed)
  y <- rep(c(0L, 1L), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  t <- seq(0, 1, length.out = p)
  base <- seq(5.2, 2.0, length.out = p)
  smooth <- 0.35 * sin(2 * pi * t) + 0.18 * cos(4 * pi * t)
  pulse <- rep(c(rep(1.0, 4L), rep(-0.25, 12L)), length.out = p)
  rough <- smooth + 0.95 * pulse + 0.30 * sin(44 * pi * t)
  X <- matrix(0, nrow = n, ncol = p,
              dimnames = list(paste0("S", seq_len(n)), features))
  for (i in seq_len(n)) {
    shape <- if (y[i] == 1L) rough else smooth
    mu <- base + shape + stats::rnorm(1L, 0, 0.25)
    X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.06))
  }
  list(X = X, y = y)
}

.auc_wilcox_frac <- function(y, score) {
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2L) / (n1 * n0)
}

.ind_mfdfa_hq <- function(x, scales, q, m_poly, eps) {
  Y <- cumsum(x - mean(x))
  log_s <- log(scales)
  h <- numeric(length(q))
  for (iq in seq_along(q)) {
    log_fq <- numeric(length(scales))
    for (is in seq_along(scales)) {
      s <- scales[is]
      Ns <- floor(length(Y) / s)
      f2 <- numeric(2L * Ns)
      u <- seq(-1, 1, length.out = s)
      design <- outer(u, seq.int(0L, m_poly), "^")
      k <- 1L
      for (nu in seq_len(Ns)) {
        idx <- ((nu - 1L) * s + 1L):(nu * s)
        fit <- stats::lm.fit(design, Y[idx])
        f2[k] <- mean(fit$residuals^2)
        k <- k + 1L
      }
      for (nu in seq_len(Ns)) {
        idx <- (length(Y) - nu * s + 1L):(length(Y) - (nu - 1L) * s)
        fit <- stats::lm.fit(design, Y[idx])
        f2[k] <- mean(fit$residuals^2)
        k <- k + 1L
      }
      f2 <- pmax(f2, eps)
      if (q[iq] == 0) {
        log_fq[is] <- 0.5 * mean(log(f2))
      } else {
        a <- (q[iq] / 2) * log(f2)
        amax <- max(a)
        log_fq[is] <- (log(mean(exp(a - amax))) + amax) / q[iq]
      }
    }
    h[iq] <- stats::lm.fit(cbind(1, log_s), log_fq)$coefficients[2L]
  }
  h
}


test_that("frac-mfdfa fit/score roundtrip has the right shape and class", {
  dat <- .make_frac_mfdfa_data(n = 120L, seed = 72L)
  model <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())
  score <- score_frac_mfdfa(model, dat$X)

  expect_s3_class(model, "frac_mfdfa_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_setequal(model$order_perm, colnames(dat$X))
  expect_identical(model$descriptor_dim, length(model$hp$q) + 3L)
  expect_identical(model$scales, model$hp$scales)
  expect_length(model$head$center, model$descriptor_dim)
  expect_length(model$head$scale, model$descriptor_dim)
  expect_length(model$head$w, model$descriptor_dim)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
  expect_null(names(score))
})


test_that("frac-mfdfa separates a planted intermittent landscape contrast", {
  dat <- .make_frac_mfdfa_data(seed = 73L)
  tr <- c(1:60, 81:140)
  te <- c(61:80, 141:160)
  model <- fit_frac_mfdfa(dat$X[tr, ], dat$y[tr], hp = .frac_hp())
  score <- score_frac_mfdfa(model, dat$X[te, ])

  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[te] == 1L]), mean(score[dat$y[te] == 0L]))
  expect_gt(.auc_wilcox_frac(dat$y[te], score), 0.7)
})


test_that("frac-mfdfa passes the canonical row-equivariance gate", {
  dat <- .make_frac_mfdfa_data(seed = 74L)
  tr <- c(1:60, 81:140)
  te <- c(61:80, 141:160)
  model <- fit_frac_mfdfa(dat$X[tr, ], dat$y[tr], hp = .frac_hp())
  X_test <- dat$X[te, ]
  score_fun <- function(model, X, meta) score_frac_mfdfa(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})


test_that("frac-mfdfa scoring is exactly invariant to per-sample positive scaling", {
  dat <- .make_frac_mfdfa_data(seed = 75L)
  model <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_frac_mfdfa(model, row * 7),
               score_frac_mfdfa(model, row), tolerance = 1e-8)
  expect_equal(score_frac_mfdfa(model, row * 1e6),
               score_frac_mfdfa(model, row), tolerance = 1e-8)

  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e5
  expect_equal(score_frac_mfdfa(model, scaled),
               score_frac_mfdfa(model, batch), tolerance = 1e-8)

  mixed <- batch * c(2, 0.5, 10, 1, 1e3, 7)
  expect_equal(score_frac_mfdfa(model, mixed),
               score_frac_mfdfa(model, batch), tolerance = 1e-8)
})


test_that("frac-mfdfa handles partial feature overlap and column reorder", {
  dat <- .make_frac_mfdfa_data(seed = 76L)
  model <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())

  keep <- rev(sample(colnames(dat$X), 90L))
  X_part <- dat$X[1:12, keep, drop = FALSE]
  s_part <- score_frac_mfdfa(model, X_part)
  expect_length(s_part, nrow(X_part))
  expect_true(all(is.finite(s_part)))
  expect_equal(score_frac_mfdfa(model, X_part[5, , drop = FALSE]),
               s_part[5], tolerance = 1e-12)

  X_reperm <- X_part[, sample(ncol(X_part)), drop = FALSE]
  expect_equal(score_frac_mfdfa(model, X_reperm), s_part, tolerance = 1e-12)
})


test_that("frac-mfdfa is deployable through singlesample_score_call", {
  dat <- .make_frac_mfdfa_data(seed = 77L)
  tr <- c(1:60, 81:140)
  te <- c(61:80, 141:160)
  model <- fit_frac_mfdfa(dat$X[tr, ], dat$y[tr], hp = .frac_hp())
  X_test <- dat$X[te, ]
  direct <- score_frac_mfdfa(model, X_test)

  roster <- singlesample_method_roster()
  if (!"frac-mfdfa" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "frac-mfdfa"
    tmpl$fit_fn <- "fit_frac_mfdfa"
    tmpl$score_fn <- "score_frac_mfdfa"
    tmpl$pkg_status <- "new-secondary"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_frac_mfdfa", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "frac-mfdfa",
      function(model, X, meta) score_frac_mfdfa(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("frac-mfdfa", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})


test_that("frac-mfdfa returns neutral 0 below the feature floor", {
  dat <- .make_frac_mfdfa_data(seed = 78L)
  model <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())

  X_few <- dat$X[1:6, colnames(dat$X)[1:20], drop = FALSE]
  expect_equal(score_frac_mfdfa(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  X_none <- matrix(1, 4L, 3L,
                   dimnames = list(paste0("N", seq_len(4L)),
                                   c("other-a", "other-b", "other-c")))
  expect_equal(score_frac_mfdfa(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)

  Xs <- dat$X[1:5, , drop = FALSE]
  Xs[2, ] <- 0
  Xs[2, 1:20] <- dat$X[2, 1:20]
  s <- score_frac_mfdfa(model, Xs)
  expect_equal(s[2], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s)))
  expect_true(any(s[-2] != 0))
})


test_that("frac-mfdfa fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_frac_mfdfa_data(seed = 79L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())
  expect_identical(model_1, model_2)

  s1 <- score_frac_mfdfa(model_1, dat$X[1:20, ])
  s2 <- score_frac_mfdfa(model_2, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
})


test_that("frac-mfdfa hp validation errors are explicit", {
  dat <- .make_frac_mfdfa_data(seed = 80L)

  expect_error(fit_frac_mfdfa(dat$X, dat$y, hp = list(foo = 1)),
               "unknown hp")
  expect_error(fit_frac_mfdfa(dat$X, dat$y, hp = list(128L)),
               "must be named")
  expect_error(
    fit_frac_mfdfa(dat$X, dat$y,
                   hp = list(n_scales = 6L, n_scales = 8L)),
    "duplicated hp"
  )
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(q = c(-3, -1, 1, 3))),
               "include 2")
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(q = c(-3, -1, 2, 2))),
               "duplicates")
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(s_max = 8L)),
               "s_max")
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(m_poly = 0L)),
               "m_poly")
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(m_poly = 4L)),
               "m_poly")
  # s_min must exceed the detrend DOF so the smallest window is over-determined
  # (s_min = 4 with m_poly = 3 would be an exact fit -> F2 = 0 -> eps-floor bias).
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(s_min = 4L, m_poly = 3L)),
               "s_min")
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(resample_len = 3e9)),
               "resample_len")
  expect_error(fit_frac_mfdfa(dat$X, dat$y,
                              hp = .frac_hp(min_features = 3e9)),
               "min_features")
  expect_error(score_frac_mfdfa(list(), dat$X), "frac_mfdfa_model")
})


test_that("frac-mfdfa scores single rows identically standalone and in batch", {
  dat <- .make_frac_mfdfa_data(seed = 81L)
  model <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())
  s_all <- score_frac_mfdfa(model, dat$X[1:8, ])
  for (i in c(1L, 4L, 8L)) {
    expect_equal(score_frac_mfdfa(model, dat$X[i, , drop = FALSE]),
                 s_all[i], tolerance = 1e-12)
  }
})


test_that("frac-mfdfa rCLR keeps coordinates at the geometric mean", {
  g <- .frac_mfdfa_rclr_landscape(c(a = 0, b = 1, c = 2, d = 4, e = 0))
  expect_length(g, 3L)
  expect_equal(as.numeric(sort(g)), c(-log(2), 0, log(2)), tolerance = 1e-12)
  expect_length(.frac_mfdfa_rclr_landscape(rep(5, 9)), 9L)
  expect_true(all(.frac_mfdfa_rclr_landscape(rep(5, 9)) == 0))
})


test_that("frac-mfdfa MFDFA descriptors are finite on monofractal and degenerate series", {
  hp <- .frac_mfdfa_resolve_hp(.frac_hp(q = c(-4, -2, 2, 4)))
  set.seed(123)
  x <- stats::rnorm(hp$resample_len)
  h <- .frac_mfdfa_hq(x, hp$scales, hp$q, hp$m_poly, hp$eps)
  desc <- .frac_mfdfa_spectrum_descriptor(h, hp$q)

  expect_true(all(is.finite(h)))
  expect_lt(max(h) - min(h), 0.45)
  expect_lt(desc["dAlpha"], 0.9)

  f2 <- .frac_mfdfa_segment_f2(seq_len(64L), s = 16L, m_poly = 1L)
  expect_lt(max(f2), 1e-24)
  h_flat <- .frac_mfdfa_hq(rep(3, hp$resample_len), hp$scales, hp$q,
                           hp$m_poly, hp$eps)
  expect_true(all(is.finite(h_flat)))
  expect_equal(as.numeric(h_flat), rep(0, length(h_flat)), tolerance = 1e-12)
})


test_that("frac-mfdfa h(2) orders smooth and rough deterministic series as documented", {
  hp <- .frac_mfdfa_resolve_hp(.frac_hp(q = c(-4, -2, 2, 4)))
  L <- hp$resample_len
  t <- seq(0, 1, length.out = L)
  smooth <- sin(2 * pi * t)
  rough <- sin(2 * pi * t) + 0.8 * sin(42 * pi * t)
  h_smooth <- .frac_mfdfa_hq(smooth, hp$scales, hp$q, hp$m_poly, hp$eps)
  h_rough <- .frac_mfdfa_hq(rough, hp$scales, hp$q, hp$m_poly, hp$eps)

  # With this fixed ladder, added high-frequency intermittency lowers h(2).
  expect_gt(h_smooth[hp$q == 2], h_rough[hp$q == 2])
})


test_that("frac-mfdfa internal h(q) matches an independent MFDFA recompute", {
  hp <- .frac_mfdfa_resolve_hp(.frac_hp(q = c(-5, -3, -1, 2, 3, 5)))
  t <- seq(0, 1, length.out = hp$resample_len)
  x <- sin(6 * pi * t) + 0.35 * cos(18 * pi * t) + 0.2 * t

  h_pkg <- .frac_mfdfa_hq(x, hp$scales, hp$q, hp$m_poly, hp$eps)
  h_ind <- .ind_mfdfa_hq(x, hp$scales, hp$q, hp$m_poly, hp$eps)
  expect_equal(as.numeric(h_pkg), h_ind, tolerance = 1e-8)
})


test_that("frac-mfdfa cached scale windows are bit-identical to the legacy q loop", {
  hp <- .frac_mfdfa_resolve_hp(.frac_hp(q = c(-5, -3, -1, 0, 2, 3, 5)))
  inputs <- list(
    sin(seq(0, 14 * pi, length.out = hp$resample_len)),
    rep(3, hp$resample_len),
    seq(-2, 2, length.out = hp$resample_len)^3
  )
  for (x in inputs) {
    cached <- .frac_mfdfa_hq(x, hp$scales, hp$q, hp$m_poly, hp$eps)
    legacy <- .ind_mfdfa_hq(x, hp$scales, hp$q, hp$m_poly, hp$eps)
    expect_identical(unname(cached), legacy)
  }
})


test_that("frac-mfdfa ridge-LDA head solves the normal equations", {
  set.seed(5)
  Phi <- matrix(stats::rnorm(40 * 6), 40, 6)
  y <- rep(c(0L, 1L), each = 20L)
  Phi[y == 1L, ] <- Phi[y == 1L, ] + 0.8
  head <- .frac_mfdfa_fit_head(Phi, y, shrink = 0.1, eps = 1e-8)

  center <- colMeans(Phi)
  raw_sd <- apply(Phi, 2L, stats::sd)
  active <- raw_sd > 0
  scale <- pmax(raw_sd, 1e-8)
  Z <- sweep(sweep(Phi, 2L, center, "-"), 2L, scale, "/")
  Z[, !active] <- 0
  mu1 <- colMeans(Z[y == 1L, , drop = FALSE])
  mu0 <- colMeans(Z[y == 0L, , drop = FALSE])
  Sw <- (crossprod(sweep(Z[y == 1L, , drop = FALSE], 2L, mu1, "-")) +
         crossprod(sweep(Z[y == 0L, , drop = FALSE], 2L, mu0, "-"))) /
    (nrow(Z) - 2L)
  Sw <- (Sw + t(Sw)) / 2
  Swr <- Sw
  diag(Swr) <- diag(Swr) + head$ridge
  w_check <- as.numeric(chol2inv(chol(Swr)) %*% (mu1 - mu0))
  b_check <- -0.5 * sum((mu1 + mu0) * w_check)

  expect_equal(head$w, w_check, tolerance = 1e-8)
  expect_equal(head$b, b_check, tolerance = 1e-8)

  dat <- .make_frac_mfdfa_data(seed = 82L)
  model <- fit_frac_mfdfa(dat$X, dat$y, hp = .frac_hp())
  phi <- .frac_mfdfa_descriptor(dat$X[7, model$order_perm], model$hp)
  z <- (phi - model$head$center) / model$head$scale
  z[!model$head$active] <- 0
  manual <- model$head$b + sum(model$head$w * z)
  expect_equal(score_frac_mfdfa(model, dat$X[7, , drop = FALSE]), manual,
               tolerance = 1e-9)
})


test_that("frac-mfdfa ridge makes a collinear descriptor covariance positive definite", {
  set.seed(9)
  Phi <- matrix(stats::rnorm(8 * 6, sd = 1e-7), 8, 6)
  Phi[, 1] <- Phi[, 2]
  y <- rep(c(0L, 1L), each = 4L)
  head <- .frac_mfdfa_fit_head(Phi, y, shrink = 0.1, eps = 1e-8)
  expect_true(all(is.finite(head$w)))
  expect_true(is.finite(head$b))
  expect_true(head$ridge > 0)

  center <- colMeans(Phi)
  scale <- pmax(apply(Phi, 2L, stats::sd), 1e-8)
  Z <- sweep(sweep(Phi, 2L, center, "-"), 2L, scale, "/")
  mu1 <- colMeans(Z[y == 1L, , drop = FALSE])
  mu0 <- colMeans(Z[y == 0L, , drop = FALSE])
  Sw <- (crossprod(sweep(Z[y == 1L, , drop = FALSE], 2L, mu1, "-")) +
         crossprod(sweep(Z[y == 0L, , drop = FALSE], 2L, mu0, "-"))) /
    max(nrow(Z) - 2L, 1L)
  Sw <- (Sw + t(Sw)) / 2
  Swr <- Sw
  diag(Swr) <- diag(Swr) + head$ridge
  expect_error(chol(Swr), NA)
})
