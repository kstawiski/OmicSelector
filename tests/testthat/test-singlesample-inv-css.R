library(testthat)

# ---------------------------------------------------------------------------
# Synthetic data with a planted MULTI-SCALE SHAPE difference in the sorted rCLR
# profile (the substrate inv-css reads). Both classes share a comparable level
# and range; only the *shape* of the sorted profile differs:
#   - Controls: per-feature log-means follow a smooth LINEAR ramp -> the sorted
#     rCLR profile is nearly linear -> low curvature at every scale.
#   - Cases: per-feature log-means follow a sharp STEP (two plateaus) -> the
#     sorted rCLR profile has a steep shoulder/kink -> high curvature (large
#     energy, large peak, an inflection pair) across scales.
# The signal is a pure rank-shape difference (scale-free), exactly what the CSS
# descriptor is designed to capture, and the class means are matched (3.5 each)
# so the contrast is shape, not level.
# ---------------------------------------------------------------------------
.make_inv_css_data <- function(n = 200L, p = 60L, seed = 71L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(0, nrow = n, ncol = p,
              dimnames = list(paste0("S", seq_len(n)), features))
  ctrl_mu <- seq(2, 5, length.out = p)                       # smooth ramp shape
  case_mu <- ifelse(seq_len(p) <= p / 2L, 2.75, 4.25)        # sharp step shape
  for (i in seq_len(n)) {
    mu <- if (y[i] == 1L) case_mu else ctrl_mu
    X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
  }
  list(X = X, y = y)
}

.auc_or_wilcoxon <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


test_that("inv-css fit/score roundtrip has the right shape and class", {
  dat <- .make_inv_css_data()
  model <- fit_inv_css(dat$X, dat$y)
  score <- score_inv_css(model, dat$X)

  expect_s3_class(model, "inv_css_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_identical(model$descriptor_dim, 6L * 4L)
  expect_length(model$scales, 6L)
  expect_length(model$head$center, model$descriptor_dim)
  expect_length(model$head$scale, model$descriptor_dim)
  expect_length(model$head$w, model$descriptor_dim)

  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})


test_that("inv-css separates a planted multi-scale shape difference", {
  dat <- .make_inv_css_data(seed = 72L)
  tr <- c(1:70, 101:170)
  te <- c(71:100, 171:200)
  model <- fit_inv_css(dat$X[tr, ], dat$y[tr])
  score <- score_inv_css(model, dat$X[te, ])

  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[te] == 1]), mean(score[dat$y[te] == 0]))
  expect_gt(.auc_or_wilcoxon(dat$y[te], score), 0.7)
})


test_that("inv-css passes the canonical row-equivariance gate", {
  dat <- .make_inv_css_data(seed = 73L)
  tr <- c(1:70, 101:170)
  te <- c(71:100, 171:200)
  model <- fit_inv_css(dat$X[tr, ], dat$y[tr])
  X_test <- dat$X[te, ]
  score_fun <- function(model, X, meta) score_inv_css(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})


test_that("inv-css scoring is exactly invariant to per-sample positive scaling", {
  dat <- .make_inv_css_data(seed = 74L)
  model <- fit_inv_css(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_inv_css(model, row * 7),
               score_inv_css(model, row), tolerance = 1e-8)
  expect_equal(score_inv_css(model, row * 1e6),
               score_inv_css(model, row), tolerance = 1e-8)

  # Scaling some rows of a batch leaves the other rows' scores untouched.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e5
  expect_equal(score_inv_css(model, scaled),
               score_inv_css(model, batch), tolerance = 1e-8)

  # Mixed: each row scaled by its own positive factor (column-major recycling of
  # a length-nrow vector scales row i by factor[i]).
  mixed <- batch * c(2, 0.5, 10, 1, 1e3, 7)
  expect_equal(score_inv_css(model, mixed),
               score_inv_css(model, batch), tolerance = 1e-8)
})


test_that("inv-css handles partial feature overlap (subset + reorder) consistently", {
  dat <- .make_inv_css_data(seed = 75L)
  model <- fit_inv_css(dat$X, dat$y)

  keep <- rev(sample(colnames(dat$X), 30L))     # subset and reorder columns
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_inv_css(model, X_partial)

  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # A single row scored alone equals its score within the batch.
  expect_equal(score_inv_css(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # The descriptor sorts values, so re-permuting present columns cannot change it.
  X_reperm <- X_partial[, sample(ncol(X_partial)), drop = FALSE]
  expect_equal(score_inv_css(model, X_reperm), s_partial, tolerance = 1e-12)
})


test_that("inv-css is deployable through the canonical dispatch", {
  skip_if_not(exists("singlesample_score_call"))
  dat <- .make_inv_css_data(seed = 76L)
  tr <- c(1:70, 101:170)
  te <- c(71:100, 171:200)
  model <- fit_inv_css(dat$X[tr, ], dat$y[tr])
  X_test <- dat$X[te, ]
  direct <- score_inv_css(model, X_test)

  # Post-integration the manifest carries the inv-css row and score_inv_css lives
  # in the OmicSelector namespace, so the canonical adapter resolves it directly.
  # In the pre-integration staged self-test score_inv_css is not yet in the
  # namespace, so we clone the roster row if absent and register a thin canonical
  # adapter to exercise the same dispatch wiring. Either way the dispatched score
  # must equal the direct call.
  roster <- singlesample_method_roster()
  if (!"inv-css" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "inv-css"
    tmpl$fit_fn <- "fit_inv_css"
    tmpl$score_fn <- "score_inv_css"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_inv_css", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "inv-css",
      function(model, X, meta) score_inv_css(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("inv-css", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})


test_that("inv-css returns neutral 0 below the feature-overlap floor", {
  dat <- .make_inv_css_data(seed = 77L)
  model <- fit_inv_css(dat$X, dat$y)            # default min_features = 8

  # Fewer than min_features model-universe features present -> neutral 0.
  X_few <- dat$X[1:6, colnames(dat$X)[1:5], drop = FALSE]
  expect_equal(score_inv_css(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", seq_len(4L)),
                                   c("other-a", "other-b", "other-c")))
  expect_equal(score_inv_css(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)

  # Enough features present, but ONE specimen has too few nonzero coordinates ->
  # that row alone scores the neutral 0; the others score normally.
  Xs <- dat$X[1:5, , drop = FALSE]
  Xs[2, ] <- 0
  Xs[2, 1:4] <- dat$X[2, 1:4]                   # only 4 nonzero < min_features 8
  s <- score_inv_css(model, Xs)
  expect_equal(s[2], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s)))
  expect_true(any(s[-2] != 0))
})


test_that("inv-css fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_inv_css_data(seed = 78L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_inv_css(dat$X, dat$y)
  # The method consumes no randomness, so the global RNG is untouched by a fit.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_inv_css(dat$X, dat$y)
  expect_identical(model_1, model_2)

  s1 <- score_inv_css(model_1, dat$X[1:20, ])
  s2 <- score_inv_css(model_2, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
})


test_that("inv-css hp validation errors are explicit", {
  dat <- .make_inv_css_data(seed = 79L)

  expect_error(fit_inv_css(dat$X, dat$y, hp = list(foo = 1)), "unknown hp")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(n_scales = 6.5)), "n_scales")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(n_scales = 1L)), "n_scales")
  # Out-of-integer-range must hit the explicit error, not an as.integer overflow.
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(n_scales = 3e9)), "n_scales")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(resample_len = 4L)),
               "resample_len")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(sigma_max = 0.5,
                                                   sigma_min = 1)), "sigma_max")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(sigma_min = 0)), "sigma_min")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(shrink = -0.1)), "shrink")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(seed = -1)), "hp\\$seed")
  # `$` partial matching would let "min" bind to "min_features"; the exact [[ ]]
  # reads + allowed-list reject the unknown field instead.
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(min = 3L)), "unknown hp")
  # Fully-unnamed hp (names(hp) == NULL) must be rejected, not silently ignored:
  # setdiff(NULL, allowed) is empty, so the allowed-list check alone would pass it.
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(6L)), "must be named")
  # A partially-unnamed entry is rejected too.
  expect_error(fit_inv_css(dat$X, dat$y, hp = list(n_scales = 6L, 99)),
               "must be named")
  # Duplicated field names must be rejected, not silently last-wins.
  expect_error(
    fit_inv_css(dat$X, dat$y, hp = list(n_scales = 4L, n_scales = 6L)),
    "duplicated hp")
})

test_that("inv-css keeps present coordinates at the geometric mean (no z!=0 drop)", {
  # Structural-zero exclusion must key on v > 0, NOT z != 0: a present coordinate
  # whose centred-log equals the geometric mean (z == 0) must be kept, and a
  # constant positive profile must yield a length-(#positive) point mass at 0
  # (so the resampled curve is well-defined), not an empty vector.
  rclr_pos <- get(".inv_css_rclr_pos", envir = asNamespace("OmicSelector"))
  # v = (1, 2, 4): geomean = 2, so the middle present coordinate has rCLR exactly 0.
  zp <- rclr_pos(c(1, 2, 4))
  expect_length(zp, 3L)
  expect_equal(sort(zp), c(-log(2), 0, log(2)), tolerance = 1e-12)
  # Structural zeros (v == 0) are excluded; the present middle-0 is kept.
  expect_length(rclr_pos(c(0, 1, 2, 4, 0)), 3L)
  # A constant positive profile -> point mass at 0 of size = #positive features.
  expect_length(rclr_pos(rep(2, 7)), 7L)
  expect_true(all(rclr_pos(rep(2, 7)) == 0))
  # End to end: a specimen constant across all present features scores finite
  # (not dropped to neutral 0 because z != 0 emptied it).
  dat <- .make_inv_css_data(seed = 80L)
  model <- fit_inv_css(dat$X, dat$y)
  Xc <- dat$X[1:3, , drop = FALSE]
  Xc[1, ] <- 5
  expect_true(all(is.finite(score_inv_css(model, Xc))))
})


test_that("inv-css scores a single row identically standalone and in batch", {
  dat <- .make_inv_css_data(seed = 80L)
  model <- fit_inv_css(dat$X, dat$y)
  s_all <- score_inv_css(model, dat$X[1:8, ])
  for (i in c(1L, 4L, 8L)) {
    s_one <- score_inv_css(model, dat$X[i, , drop = FALSE])
    expect_length(s_one, 1L)
    expect_equal(s_one, s_all[i], tolerance = 1e-12)
  }
})


test_that("inv-css CSS curvature is correct on hand-built profiles", {
  L <- 64L
  scales <- .inv_css_kernels(6L, 1.0, 8.0)

  # (1) Convex parabola -> constant-sign curvature -> 0 inflections. For
  #     g = a (t - c)^2 the second difference is exactly 2a > 0 everywhere.
  g_par <- 0.01 * (seq_len(L) - (L + 1) / 2)^2
  kap_par <- .inv_css_curvature(g_par)
  expect_true(all(kap_par > 0))
  expect_identical(.inv_css_sign_changes(kap_par), 0L)

  # (2) S-shaped (logistic) profile -> exactly one curvature sign change (the
  #     single inflection), preserved after coarse-scale smoothing.
  g_sig <- 1 / (1 + exp(-seq(-6, 6, length.out = L)))
  expect_identical(.inv_css_sign_changes(.inv_css_curvature(g_sig)), 1L)
  for (s in c(2L, 3L, 4L)) {                    # several coarse scales
    g_sig_s <- .inv_css_smooth(g_sig, scales[[s]]$kern, scales[[s]]$h)
    expect_identical(.inv_css_sign_changes(.inv_css_curvature(g_sig_s)), 1L)
  }

  # (3) Flat-line profile -> Gaussian smoothing returns it unchanged (sum-1
  #     kernel) and curvature is ~0 at every scale.
  g_flat <- rep(3.3, L)
  for (s in seq_along(scales)) {
    g_flat_s <- .inv_css_smooth(g_flat, scales[[s]]$kern, scales[[s]]$h)
    expect_equal(g_flat_s, g_flat, tolerance = 1e-12)
    expect_lt(max(abs(.inv_css_curvature(g_flat_s))), 1e-9)
  }

  # Reflect padding is exactly scale/shift-equivariant on a constant and matches
  # the simple rev() form for h <= L.
  expect_identical(.inv_css_reflect(1:5, 2L),
                   c(2L, 1L, 1L, 2L, 3L, 4L, 5L, 5L, 4L))
})


test_that("inv-css descriptors exclude the reflect-padding boundary artifact", {
  # The sorted rCLR curve is monotone, so reflect padding mirrors a non-zero
  # boundary slope into a smoothed kink at each edge -- a curvature artifact that
  # would otherwise dominate the descriptors at medium-large scales. The per-scale
  # descriptors must summarise only the PADDING-FREE interior [h+1, L-h-2], so a
  # perfectly straight (monotone) curve -- which has zero true curvature -- yields
  # near-zero descriptors even though its boundary-included curvature is large.
  L <- 64L
  scales <- .inv_css_kernels(6L, 1.0, 8.0)
  g_line <- 3.5 * seq_len(L) + 1                  # straight monotone curve
  for (s in scales) {
    g_s <- .inv_css_smooth(g_line, s$kern, s$h)
    full <- sum(abs(.inv_css_curvature(g_s)))     # boundary-included (artifact)
    d <- .inv_css_scale_descriptors(g_s, s$h)     # padding-free interior only
    expect_lt(d[1], 1e-6)                          # total |kappa| over interior ~ 0
    expect_lt(d[3], 1e-6)                          # peak |kappa| ~ 0
    expect_identical(as.integer(d[4]), 0L)         # no spurious inflections
    expect_gt(full, 1e-2)                          # the excluded boundary held the artifact
  }
  # Guard branch: a scale whose half-width leaves no padding-free interior
  # contributes the neutral (0,0,0,0) rather than erroring.
  big <- .inv_css_kernels(1L, 30.0, 31.0)[[1]]     # h = ceil(3*30) = 90 > L
  g_s_big <- .inv_css_smooth(g_line, big$kern, big$h)
  expect_equal(.inv_css_scale_descriptors(g_s_big, big$h), c(0, 0, 0, 0))
})


test_that("inv-css head solves the ridge-LDA normal equations", {
  set.seed(5)
  Phi <- matrix(stats::rnorm(40 * 6), 40, 6)
  y <- rep(c(0, 1), each = 20)
  Phi[y == 1, ] <- Phi[y == 1, ] + 0.8         # genuine class separation
  head <- .inv_css_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)

  # Independent recompute of the standardize + ridge-LDA head.
  center <- colMeans(Phi)
  raw_sd <- apply(Phi, 2L, stats::sd)
  active <- raw_sd > 0
  scale <- pmax(raw_sd, 1e-6)
  Z <- sweep(sweep(Phi, 2L, center, "-"), 2L, scale, "/")
  Z[, !active] <- 0
  mu1 <- colMeans(Z[y == 1, , drop = FALSE])
  mu0 <- colMeans(Z[y == 0, , drop = FALSE])
  df <- nrow(Z) - 2L
  Sw <- (crossprod(sweep(Z[y == 1, , drop = FALSE], 2L, mu1, "-")) +
         crossprod(sweep(Z[y == 0, , drop = FALSE], 2L, mu0, "-"))) / df
  Sw <- (Sw + t(Sw)) / 2
  Swr <- Sw
  diag(Swr) <- diag(Swr) + head$ridge
  w_check <- as.numeric(chol2inv(chol(Swr)) %*% (mu1 - mu0))
  b_check <- -0.5 * sum((mu1 + mu0) * w_check)

  expect_equal(head$w, w_check, tolerance = 1e-8)
  expect_equal(head$b, b_check, tolerance = 1e-8)

  # And the scorer's linear predictor equals the hand-computed LDA value.
  dat <- .make_inv_css_data(seed = 81L)
  model <- fit_inv_css(dat$X, dat$y)
  phi <- .inv_css_descriptor(dat$X[7, ], model$scales, model$hp$resample_len,
                             model$hp$min_features)
  z <- (phi - model$head$center) / model$head$scale
  z[!model$head$active] <- 0
  manual <- model$head$b + sum(model$head$w * z)
  expect_equal(score_inv_css(model, dat$X[7, , drop = FALSE]), manual,
               tolerance = 1e-9)
})


test_that("inv-css head stays positive definite after ridge on a near-degenerate covariance", {
  set.seed(9)
  # 8 samples, 6 descriptors, tiny variance + an exactly collinear pair: the
  # pooled within-class covariance is (near-)singular, so the ridge-to-PD loop
  # must still produce a finite Cholesky-based head.
  Phi <- matrix(stats::rnorm(8 * 6, sd = 1e-7), 8, 6)
  Phi[, 1] <- Phi[, 2]
  y <- rep(c(0, 1), each = 4)
  head <- .inv_css_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)
  expect_true(all(is.finite(head$w)))
  expect_true(is.finite(head$b))
  expect_true(head$ridge > 0)

  # Integration: a fit on a p > n regime (few rows after the min_features filter)
  # also stays finite end to end.
  dat <- .make_inv_css_data(n = 12L, p = 40L, seed = 82L)
  model <- fit_inv_css(dat$X, dat$y)
  sc <- score_inv_css(model, dat$X)
  expect_length(sc, nrow(dat$X))
  expect_true(all(is.finite(sc)))
})
