library(testthat)

# Runs inside the package test harness (devtools::test / R CMD check), where
# fit_sinkhorn_single / score_sinkhorn_single and all package helpers are already
# loaded. This negative-control method is pure R (no Python); the only optional
# dependency below is the multi-row Sinkhorn-OT primitive used by the NEGATIVE
# test (that the out-of-N scorer is NOT row-equivariant), which is part of the
# package.

# Synthetic generator with a CASE-ONLY location shift in a feature block. The
# first `k` of `p` log-abundance features are elevated by `shift` for case
# samples, so after exp() and the per-sample rCLR centring the case specimens sit
# closer (in entropic-OT energy) to the frozen case reference cloud than to the
# control cloud. Abundances are exp() so the input is strictly positive
# (compositional).
.make_sinkhorn_single_data <- function(n = 160L, p = 30L, k = 8L, shift = 1.3,
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

.auc_mw_sks <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Independent (loop-based) recompute of the single-source entropic transport
# energy, used by the §7 anchor test to cross-check the package helper without
# sharing its code path.
.manual_sinkhorn_energy <- function(z, atoms, tau, epsilon) {
  m <- nrow(atoms)
  cost <- numeric(m)
  for (j in seq_len(m)) cost[j] <- sum((z - atoms[j, ])^2) / tau
  -epsilon * log(mean(exp(-cost / epsilon)))
}


test_that("sinkhorn-single fit/score roundtrip has the right shape and types", {
  dat <- .make_sinkhorn_single_data()
  model <- fit_sinkhorn_single(dat$X, dat$y)
  expect_s3_class(model, "sinkhorn_single_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_equal(ncol(model$case_atoms), ncol(dat$X))
  expect_equal(ncol(model$control_atoms), ncol(dat$X))
  expect_equal(model$n_case, sum(dat$y == 1L))
  expect_equal(model$n_control, sum(dat$y == 0L))
  expect_true(model$orientation %in% c(-1, 1))
  expect_true(is.finite(model$cost_scale) && model$cost_scale > 0)

  score <- score_sinkhorn_single(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("sinkhorn-single energy matches an independent recompute (§7 anchor)", {
  # §7 GROUND TRUTH. The single-source entropic-OT energy is a log-sum-exp of
  # squared transport costs; the package's stable (min-subtracted) implementation
  # must equal a direct, independent loop recompute to machine precision. This is
  # the arbiter of the entropic-kernel formula and is locked here as a regression.
  dat <- .make_sinkhorn_single_data(seed = 202L)
  model <- fit_sinkhorn_single(dat$X, dat$y)
  Z <- .sinkhorn_single_rclr_matrix(dat$X)
  full_idx <- seq_len(ncol(Z))
  tau <- model$cost_scale; eps <- model$hp$epsilon
  for (i in c(1L, 17L, 90L, 160L)) {
    z <- Z[i, ]
    pkg_case <- .sinkhorn_single_energy(z, model$case_atoms, tau, eps)
    man_case <- .manual_sinkhorn_energy(z, model$case_atoms, tau, eps)
    expect_equal(pkg_case, man_case, tolerance = 1e-12)
    pkg_ctrl <- .sinkhorn_single_energy(z, model$control_atoms, tau, eps)
    man_ctrl <- .manual_sinkhorn_energy(z, model$control_atoms, tau, eps)
    expect_equal(pkg_ctrl, man_ctrl, tolerance = 1e-12)
    # And the oriented per-specimen score equals orientation * (ctrl - case).
    pkg_score <- score_sinkhorn_single(model, dat$X[i, , drop = FALSE])
    expect_equal(pkg_score, model$orientation * (man_ctrl - man_case),
                 tolerance = 1e-10)
  }
})

test_that("sinkhorn-single separates a strongly planted case-vs-control shift", {
  dat <- .make_sinkhorn_single_data(n = 300L, seed = 303L)
  ctrl <- which(dat$y == 0L); case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_sinkhorn_single(dat$X[train, ], dat$y[train])
  score <- score_sinkhorn_single(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  # Negative control, but the wiring must separate a STRONG plant clearly.
  expect_gt(.auc_mw_sks(dat$y[test], score), 0.7)
})

test_that("sinkhorn-single single-row == batch (decisive single-sample property)", {
  # The defining property of this n=1 negative control: a specimen's score is
  # bit-identical whether scored alone or inside any batch. maxdiff must be 0.
  dat <- .make_sinkhorn_single_data(seed = 71L)
  model <- fit_sinkhorn_single(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  batch <- score_sinkhorn_single(model, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_sinkhorn_single(model, X_test[i, , drop = FALSE])
  }, numeric(1L))
  expect_equal(max(abs(batch - singles)), 0, tolerance = 1e-6)
})

test_that("sinkhorn-single passes the canonical row-equivariance gate", {
  dat <- .make_sinkhorn_single_data(seed = 43L)
  model <- fit_sinkhorn_single(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_sinkhorn_single(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, model, X_test,
      model_digest = .sinkhorn_single_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, model, X_test,
      model_digest = .sinkhorn_single_model_digest))
})

test_that("sinkhorn-single is row-permutation invariant", {
  dat <- .make_sinkhorn_single_data(seed = 91L)
  model <- fit_sinkhorn_single(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  base <- score_sinkhorn_single(model, X_test)
  perm <- sample(nrow(X_test))
  permuted <- score_sinkhorn_single(model, X_test[perm, , drop = FALSE])
  expect_equal(permuted, base[perm], tolerance = 1e-12)
})

test_that("sinkhorn-single scoring is exactly invariant to per-sample scaling", {
  # The per-sample rCLR cancels any positive per-specimen scalar exactly (up to
  # log round-off), and the frozen cost scale / clouds / orientation never see a
  # scored-batch statistic, so the score is invariant to per-sample scaling.
  dat <- .make_sinkhorn_single_data(seed = 44L)
  model <- fit_sinkhorn_single(dat$X, dat$y)
  for (i in c(1L, 3L, 81L, 140L)) {            # controls and cases
    row <- dat$X[i, , drop = FALSE]
    expect_equal(score_sinkhorn_single(model, row * 7),
                 score_sinkhorn_single(model, row), tolerance = 1e-8)
    expect_equal(score_sinkhorn_single(model, row * 1e6),
                 score_sinkhorn_single(model, row), tolerance = 1e-8)
    expect_equal(score_sinkhorn_single(model, row * 1e-3),
                 score_sinkhorn_single(model, row), tolerance = 1e-8)
  }
  # Per-row random positive scaling within a batch leaves every score fixed.
  batch <- dat$X[1:8, , drop = FALSE]
  set.seed(5)
  scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  expect_equal(score_sinkhorn_single(model, scaled),
               score_sinkhorn_single(model, batch), tolerance = 1e-8)
})

test_that("sinkhorn-single handles partial feature overlap (consistent + single-row)", {
  dat <- .make_sinkhorn_single_data(seed = 45L)
  model <- fit_sinkhorn_single(dat$X, dat$y)
  keep <- rev(c(dat$block[1:5], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_sinkhorn_single(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_sinkhorn_single(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_sinkhorn_single(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("sinkhorn-single matches the canonical dispatch", {
  dat <- .make_sinkhorn_single_data(seed = 46L)
  model <- fit_sinkhorn_single(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  direct <- score_sinkhorn_single(model, X_test)

  # Post-integration the manifest carries the nc-sinkhorn-single row and
  # score_sinkhorn_single lives in the OmicSelector namespace, so the canonical
  # adapter resolves it directly. In the pre-integration staged self-test the
  # namespace entry does not yet exist; clone a one-row roster from a committed
  # within-discriminator row and register a thin adapter so the same dispatch
  # wiring is exercised. Either way the dispatched score must equal the direct
  # call.
  roster <- singlesample_method_roster()
  if (!"nc-sinkhorn-single" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "nc-sinkhorn-single"
    tmpl$fit_fn <- "fit_sinkhorn_single"
    tmpl$score_fn <- "score_sinkhorn_single"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_sinkhorn_single", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "nc-sinkhorn-single",
      function(model, X, meta) score_sinkhorn_single(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("nc-sinkhorn-single", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("sinkhorn-single returns the neutral 0 on degenerate input", {
  dat <- .make_sinkhorn_single_data(seed = 47L)
  model <- fit_sinkhorn_single(dat$X, dat$y, hp = list(min_features = 3L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, colnames(dat$X)[1:2], drop = FALSE]
  expect_equal(score_sinkhorn_single(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_sinkhorn_single(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
  # An all-zero specimen (empty positive support) scores the neutral 0; the
  # remaining non-degenerate rows still score finite.
  X_zero <- dat$X[1:3, , drop = FALSE]
  X_zero[1, ] <- 0
  s_zero <- score_sinkhorn_single(model, X_zero)
  expect_equal(s_zero[1], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s_zero)))
})

test_that("sinkhorn-single scores a FLAT composition (not floored as degenerate)", {
  # Regression: an earlier empty-support check used all(z == 0), which ALSO floored
  # a valid FLAT composition (full positive support, EQUAL abundances) -- its
  # per-sample rCLR maps to the ORIGIN (all zeros) but it is a legitimate specimen
  # that must be transported to the frozen case/control clouds, not floored to
  # neutral 0. Only a genuinely empty positive support (no positive abundance) is
  # degenerate. The fix tests the ORIGINAL abundances (any(X_use[i, ] > 0)).
  dat <- .make_sinkhorn_single_data(seed = 51L)
  model <- fit_sinkhorn_single(dat$X, dat$y, hp = list(min_features = 3L))
  p <- length(model$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", model$feature_universe))
  z_flat <- .sinkhorn_single_rclr(flat[1, ])
  expect_true(all(flat > 0))         # full positive support
  expect_true(all(z_flat == 0))      # ... yet rCLR is the origin (all zeros)

  # It must be SCORED (entropic energy difference from the origin to the frozen
  # clouds), reproducing the raw per-row scorer -- NOT floored to 0.
  s_flat <- score_sinkhorn_single(model, flat)
  manual <- model$orientation *
    .sinkhorn_single_score_one(z_flat, seq_len(p), model$case_atoms,
                               model$control_atoms, model$cost_scale,
                               model$hp$epsilon)
  expect_equal(s_flat, manual, tolerance = 1e-12)
  expect_true(is.finite(s_flat) && abs(s_flat) > 1e-9)   # genuinely scored, not 0

  # Contrast: a genuinely empty positive support IS floored to neutral 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", model$feature_universe))
  expect_equal(score_sinkhorn_single(model, empty), 0, tolerance = 1e-12)
})

test_that("sinkhorn-single applies the frozen orientation scalar", {
  # The label-defined case / control clouds make the raw discriminant
  # (control_energy - case_energy) point toward case by construction, so the
  # fitted orientation is +1 here. The orientation is nonetheless a genuine,
  # applied frozen scalar: negating it negates every score. (A flip to -1 would
  # require the raw discriminant to anti-correlate with case on training, which
  # this self-referential cloud construction structurally prevents -- documented
  # honestly rather than fabricated.)
  dat <- .make_sinkhorn_single_data(seed = 52L)
  model <- fit_sinkhorn_single(dat$X, dat$y)
  expect_equal(model$orientation, 1)
  s_pos <- score_sinkhorn_single(model, dat$X[1:10, ])
  model_neg <- model; model_neg$orientation <- -1L
  s_neg <- score_sinkhorn_single(model_neg, dat$X[1:10, ])
  expect_equal(s_neg, -s_pos, tolerance = 1e-12)
})

test_that("sinkhorn-single fitting/scoring is deterministic and RNG-safe", {
  dat <- .make_sinkhorn_single_data(seed = 48L)
  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_sinkhorn_single(dat$X, dat$y)
  # Fit uses a LOCAL seed for the atom subsample and restores the global RNG.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_sinkhorn_single(dat$X, dat$y)
  expect_identical(model_1, model_2)

  s1 <- score_sinkhorn_single(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_sinkhorn_single(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  # Scoring never touches the global RNG.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
  # The model holds no external pointer, so the digest is a stable snapshot.
  expect_identical(.sinkhorn_single_model_digest(model_1),
                   .sinkhorn_single_model_digest(model_2))
})

test_that("sinkhorn-single fit restores RNG when no global seed pre-exists", {
  # If the session has never drawn a random number, .Random.seed does not exist;
  # the local seeding must not leak a .Random.seed into the global environment.
  dat <- .make_sinkhorn_single_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  model <- fit_sinkhorn_single(dat$X, dat$y)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "sinkhorn_single_model")
})

test_that("sinkhorn-single hyperparameter validation is strict", {
  dat <- .make_sinkhorn_single_data(seed = 50L)
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(epsilon = 0)),
               "hp\\$epsilon")
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(epsilon = -1)),
               "hp\\$epsilon")
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(min_features = 1.5)),
               "min_features")
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(max_atoms = 0)),
               "max_atoms")
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(seed = 1.5)), "seed")
  # Unnamed / duplicate hp fields are rejected before the allow-list.
  expect_error(fit_sinkhorn_single(dat$X, dat$y, hp = list(3L)),
               "must be named")
  expect_error(
    fit_sinkhorn_single(dat$X, dat$y,
                        hp = structure(list(0.1, 0.1),
                                       names = c("epsilon", "epsilon"))),
    "duplicate hp")
})

test_that("the multi-row Sinkhorn-OT primitive is NOT row-equivariant (negative test)", {
  # The whole reason this method exists: the multi-row primitive
  # score_sinkhorn_ot_scorer couples the scored rows (joint transport plan +
  # batch-median cost scale), so its score for a specimen depends on the batch.
  # This negative test pins that contrast -- the out-of-N scorer must FAIL the
  # single-sample gate that nc-sinkhorn-single passes. Skipped only if the
  # primitive cannot be fit in this environment (it needs >=2 train cohorts).
  skip_if_not(exists("fit_sinkhorn_ot_scorer", mode = "function") &&
                exists("score_sinkhorn_ot_scorer", mode = "function"),
              "multi-row Sinkhorn-OT primitive not available")
  dat <- .make_sinkhorn_single_data(n = 200L, seed = 61L)
  meta <- data.frame(
    cohort = rep(c("A", "B", "C", "D"), length.out = nrow(dat$X)),
    stringsAsFactors = FALSE)
  prim <- tryCatch(
    fit_sinkhorn_ot_scorer(dat$X, dat$y, meta, cohort_col = "cohort",
                           tune = FALSE, n_barycenter_atoms = 60L,
                           max_cohort_atoms = 30L, barycenter_iter = 2L,
                           max_sinkhorn_iter = 80L),
    error = function(e) e)
  if (inherits(prim, "error")) {
    skip(paste("primitive fit unavailable:", conditionMessage(prim)))
  }
  X_test <- dat$X[1:40, ]
  prim_score_fun <- function(model, X, meta) {
    as.numeric(score_sinkhorn_ot_scorer(X, model))
  }
  expect_false(
    singlesample_is_row_equivariant(prim_score_fun, prim, X_test))
})
