library(testthat)

# Runs inside the package test harness (devtools::test / R CMD check), where
# fit_ecod_copod / score_ecod_copod and all package helpers are already loaded.
# Fitting needs the Python pyod module via reticulate; scoring is pure R. Tests
# that fit are skipped when pyod is unavailable.
skip_if_no_pyod <- function() {
  venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
  if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = venv)
  }
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(reticulate::py_module_available("pyod"),
                        "Python module 'pyod' not available")
}

# Synthetic generator with a CASE-ONLY tail/location shift in a feature block.
# The first `k` of `p` log-abundance features are elevated by `shift` for case
# samples, so after exp() and the per-sample rCLR centring those features sit in
# the upper tail relative to the control reference -- exactly the per-feature
# extreme-tail signal that ECOD/COPOD aggregate. Abundances are exp() so the
# input is strictly positive (compositional).
.make_ecod_copod_data <- function(n = 160L, p = 30L, k = 8L, shift = 1.3,
                                  sd = 0.5, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case <- which(y == 1L)
  L[case, seq_len(k)] <- L[case, seq_len(k)] + shift
  X <- exp(L)
  list(X = X, y = y, block = features[seq_len(k)])
}

.auc_mw_ecod <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Independent (loop-based) recompute of pyod's ECOD / COPOD per-specimen
# aggregates against a frozen reference, used only by the §7 anchor test to
# cross-check the package helper without sharing its code path.
.manual_ecod_copod <- function(z, ref_cols, skew_sign) {
  ecod <- 0; copod <- 0
  for (j in seq_along(z)) {
    ref <- ref_cols[[j]]
    n <- length(ref)
    u_l <- -log(max(sum(ref <= z[j]), 0.5) / n)
    u_r <- -log(max(sum(ref >= z[j]), 0.5) / n)
    s <- skew_sign[[j]]
    u_skew <- if (s > 0) u_r else if (s < 0) u_l else (u_l + u_r)
    ecod <- ecod + max(u_l, u_r, u_skew)
    copod <- copod + max(u_skew, (u_l + u_r) / 2)
  }
  c(ecod = ecod, copod = copod)
}


test_that("ecod-copod fit/score roundtrip has the right shape and types", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data()
  model <- fit_ecod_copod(dat$X, dat$y)
  expect_s3_class(model, "ecod_copod_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_length(model$ref_cols, ncol(dat$X))
  expect_length(model$skew_sign, ncol(dat$X))
  expect_equal(model$n_control, sum(dat$y == 0L))
  expect_true(model$orientation %in% c(-1, 1))
  expect_named(model$score_mean, c("ecod", "copod"))
  expect_named(model$score_sd, c("ecod", "copod"))

  score <- score_ecod_copod(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("ecod-copod frozen reimplementation matches pyod decision_scores_ (<1e-8)", {
  skip_if_no_pyod()
  # §7 EXTERNAL GROUND TRUTH. On the control TRAINING rows the frozen per-feature
  # reference ECDF equals the fit-time ECDF, so the package's pure-R ECOD/COPOD
  # aggregates must reproduce pyod's decision_scores_ to machine precision. This
  # is the arbiter of exact ECDF/skew-convention correctness and is locked here
  # as a regression.
  dat <- .make_ecod_copod_data(seed = 202L)
  model <- fit_ecod_copod(dat$X, dat$y)

  Zc <- .ecod_copod_rclr_matrix(dat$X[dat$y == 0L, , drop = FALSE])
  ecod_mod <- reticulate::import("pyod.models.ecod", delay_load = FALSE)
  copod_mod <- reticulate::import("pyod.models.copod", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  Zc_py <- np$asarray(Zc, dtype = "float64")
  ef <- ecod_mod$ECOD(); ef$fit(Zc_py)
  cf <- copod_mod$COPOD(); cf$fit(Zc_py)
  ds_ecod <- as.numeric(ef$decision_scores_)
  ds_copod <- as.numeric(cf$decision_scores_)

  full_idx <- seq_len(ncol(Zc))
  agg <- t(vapply(seq_len(nrow(Zc)), function(i) {
    .ecod_copod_aggregates(Zc[i, ], model$feature_universe, full_idx,
                           model$ref_cols, model$skew_sign)
  }, numeric(2L)))

  expect_lt(max(abs(agg[, "ecod"] - ds_ecod)), 1e-8)
  expect_lt(max(abs(agg[, "copod"] - ds_copod)), 1e-8)

  # And the independent loop-based recompute agrees with the package helper.
  agg_manual <- t(vapply(seq_len(nrow(Zc)), function(i) {
    .manual_ecod_copod(Zc[i, ], model$ref_cols, model$skew_sign)
  }, numeric(2L)))
  expect_equal(agg, agg_manual, tolerance = 1e-12)
})

test_that("ecod-copod separates a planted case-vs-control tail shift", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data(n = 300L, seed = 303L)
  ctrl <- which(dat$y == 0L); case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_ecod_copod(dat$X[train, ], dat$y[train])
  score <- score_ecod_copod(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_mw_ecod(dat$y[test], score), 0.7)
})

test_that("ecod-copod passes the canonical row-equivariance gate", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data(seed = 43L)
  model <- fit_ecod_copod(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_ecod_copod(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(score_fun, model, X_test,
                                        model_digest = .ecod_copod_model_digest))
  expect_true(
    singlesample_is_row_equivariant(score_fun, model, X_test,
                                    model_digest = .ecod_copod_model_digest))
})

test_that("ecod-copod scoring is exactly invariant to per-sample scaling", {
  skip_if_no_pyod()
  # The frozen ECDF is a STEP function of the rCLR value, so a control row scored
  # against a reference that contains its own value sits exactly on step edges.
  # The tie tolerance keeps the inclusive count -- and hence the score -- exactly
  # invariant to the round-off introduced by per-sample scaling.
  dat <- .make_ecod_copod_data(seed = 44L)
  model <- fit_ecod_copod(dat$X, dat$y)
  for (i in c(1L, 3L, 81L, 140L)) {            # controls and cases
    row <- dat$X[i, , drop = FALSE]
    expect_equal(score_ecod_copod(model, row * 7),
                 score_ecod_copod(model, row), tolerance = 1e-8)
    expect_equal(score_ecod_copod(model, row * 1e6),
                 score_ecod_copod(model, row), tolerance = 1e-8)
    expect_equal(score_ecod_copod(model, row * 1e-3),
                 score_ecod_copod(model, row), tolerance = 1e-8)
  }
  # Per-row random positive scaling within a batch leaves every score fixed.
  batch <- dat$X[1:8, , drop = FALSE]
  set.seed(5)
  scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  expect_equal(score_ecod_copod(model, scaled),
               score_ecod_copod(model, batch), tolerance = 1e-8)
})

test_that("ecod-copod tie tolerance preserves genuinely distinct fine-grained ranks", {
  skip_if_no_pyod()
  # Regression for an over-snapping tie tolerance. An earlier fixed 1e-9 tolerance
  # merged GENUINELY DISTINCT rCLR values: controls whose raw log-abundances differ
  # by only 1e-10 have real nonzero rCLR gaps (~6.7e-11, ~3.3e-11) -- five orders of
  # magnitude above the measured ~4e-15 per-sample-scaling round-off, yet below the
  # 1e-9 blanket, so the snap made the frozen pure-R aggregates diverge from pyod's
  # decision_scores_ and the verify_pyod fit hard-failed (maxdiff ~4.16). The
  # ulp-relative tolerance (~4.5e-13) must leave these ranks distinct, fit cleanly,
  # and keep the anchor exact.
  y <- rep(c(0, 1), each = 4L)
  d <- c(0, 1e-10, 2e-10, 3e-10, 1, 1.1, 1.2, 1.3)
  X <- cbind(exp(d), 1, 1)
  colnames(X) <- paste0("f", 1:3)
  rownames(X) <- paste0("S", seq_len(8L))

  ctrl_z <- .ecod_copod_rclr_matrix(X[y == 0L, , drop = FALSE])
  min_gap <- min(vapply(seq_len(ncol(ctrl_z)), function(j) {
    g <- diff(sort(unique(ctrl_z[, j]))); if (length(g)) min(g) else Inf
  }, numeric(1L)))
  expect_gt(min_gap, 0)      # the ranks are genuinely distinct
  expect_lt(min_gap, 1e-9)   # ... yet finer than the old blanket tolerance

  # Default fit enforces the pyod anchor internally; over-snapping would error here.
  model <- fit_ecod_copod(X, y)
  expect_s3_class(model, "ecod_copod_model")

  # The frozen aggregates still reproduce pyod's decision_scores_ on the controls.
  ecod_mod <- reticulate::import("pyod.models.ecod", delay_load = FALSE)
  copod_mod <- reticulate::import("pyod.models.copod", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  Zc_py <- np$asarray(ctrl_z, dtype = "float64")
  ef <- ecod_mod$ECOD(); ef$fit(Zc_py)
  cf <- copod_mod$COPOD(); cf$fit(Zc_py)
  full_idx <- seq_len(ncol(ctrl_z))
  agg <- t(vapply(seq_len(nrow(ctrl_z)), function(i) {
    .ecod_copod_aggregates(ctrl_z[i, ], model$feature_universe, full_idx,
                           model$ref_cols, model$skew_sign)
  }, numeric(2L)))
  expect_lt(max(abs(agg[, "ecod"] - as.numeric(ef$decision_scores_))), 1e-8)
  expect_lt(max(abs(agg[, "copod"] - as.numeric(cf$decision_scores_))), 1e-8)
})

test_that("ecod-copod handles partial feature overlap (consistent + single-row)", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data(seed = 45L)
  model <- fit_ecod_copod(dat$X, dat$y)
  keep <- rev(c(dat$block[1:5], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_ecod_copod(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_ecod_copod(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_ecod_copod(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("ecod-copod matches the canonical dispatch", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data(seed = 46L)
  model <- fit_ecod_copod(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  direct <- score_ecod_copod(model, X_test)

  # Post-integration the manifest carries the nc-ecod-copod row and
  # score_ecod_copod lives in the OmicSelector namespace, so the canonical
  # adapter resolves it directly. In the pre-integration staged self-test the
  # namespace entry does not yet exist; clone a one-row roster from a committed
  # within-discriminator row and register a thin adapter so the same dispatch
  # wiring is exercised. Either way the dispatched score must equal the direct
  # call.
  roster <- singlesample_method_roster()
  if (!"nc-ecod-copod" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "nc-ecod-copod"
    tmpl$fit_fn <- "fit_ecod_copod"
    tmpl$score_fn <- "score_ecod_copod"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_ecod_copod", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "nc-ecod-copod",
      function(model, X, meta) score_ecod_copod(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("nc-ecod-copod", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("ecod-copod returns the neutral 0 on degenerate input", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data(seed = 47L)
  model <- fit_ecod_copod(dat$X, dat$y, hp = list(min_features = 3L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, colnames(dat$X)[1:2], drop = FALSE]
  expect_equal(score_ecod_copod(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_ecod_copod(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
  # An all-zero specimen (empty positive support) scores the neutral 0; the
  # remaining non-degenerate rows still score finite.
  X_zero <- dat$X[1:3, , drop = FALSE]
  X_zero[1, ] <- 0
  s_zero <- score_ecod_copod(model, X_zero)
  expect_equal(s_zero[1], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s_zero)))
})

test_that("ecod-copod orientation flips when anomaly anti-correlates with case", {
  skip_if_no_pyod()
  # Controls are DISPERSED (tail-heavy reference) and cases sit tightly at the
  # control centre, so cases look in-distribution and the controls are the
  # anomalous ones. The +1-oriented anomaly score then anti-correlates with case
  # (training AUC < 0.5), so fit must store orientation = -1; after orientation
  # the oriented score is again positively associated with case.
  set.seed(7)
  n <- 160L; p <- 30L
  y <- rep(c(0, 1), each = n / 2L)
  feat <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(0, n, p, dimnames = list(paste0("S", seq_len(n)), feat))
  L[y == 0L, ] <- stats::rnorm(sum(y == 0L) * p, 4, 1.4)   # wide controls
  L[y == 1L, ] <- stats::rnorm(sum(y == 1L) * p, 4, 0.35)  # tight cases
  X <- exp(L)
  model <- fit_ecod_copod(X, y)
  expect_equal(model$orientation, -1)
  expect_gt(.auc_mw_ecod(y, score_ecod_copod(model, X)), 0.5)
})

test_that("ecod-copod fitting/scoring is deterministic and RNG-safe", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data(seed = 48L)
  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_ecod_copod(dat$X, dat$y)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_ecod_copod(dat$X, dat$y)
  expect_identical(model_1, model_2)

  s1 <- score_ecod_copod(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_ecod_copod(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  # Scoring never touches the global RNG.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
  # The model holds no external pointer, so the digest is a stable snapshot.
  expect_identical(.ecod_copod_model_digest(model_1),
                   .ecod_copod_model_digest(model_2))
})

test_that("ecod-copod hyperparameter validation is strict", {
  skip_if_no_pyod()
  dat <- .make_ecod_copod_data(seed = 50L)
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(min_features = 1.5)),
               "min_features")
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(eps = -1)), "hp\\$eps")
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(verify_pyod = NA)),
               "verify_pyod")
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(verify_pyod = "yes")),
               "verify_pyod")
  # Unnamed / duplicate hp fields are rejected before the allow-list.
  expect_error(fit_ecod_copod(dat$X, dat$y, hp = list(3L)),
               "must be named")
  expect_error(
    fit_ecod_copod(dat$X, dat$y,
                   hp = structure(list(1e-8, 1e-8), names = c("eps", "eps"))),
    "duplicate hp")
})
