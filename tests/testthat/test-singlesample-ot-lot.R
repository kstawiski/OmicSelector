library(testthat)

# Synthetic generator with a CASE-specific DISTRIBUTIONAL-SHAPE difference (the
# design domain of a CDT-LDA): case specimens' log-abundances are drawn from a
# MORE DISPERSED (and, optionally, heavier-tailed) generator than controls. A
# constant per-feature baseline keeps abundances strictly positive (compositional)
# while leaving the per-sample rCLR geometry driven by the dispersion of the noise
# -- exactly the spread/tail signal the 1-D CDT embedding reads. A pure location
# shift would be removed by the within-sample rCLR centring, so it is the SHAPE
# (not the mean) of the centred-log-ratio profile that separates the classes.
.make_ot_lot_data <- function(n = 200L, p = 40L, base = 6, sd_ctrl = 0.30,
                              sd_case = 0.60, df = NULL, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  E <- matrix(0, nrow = n, ncol = p,
              dimnames = list(paste0("S", seq_len(n)), features))
  ctrl <- which(y == 0L)
  case <- which(y == 1L)
  E[ctrl, ] <- stats::rnorm(length(ctrl) * p, mean = 0, sd = sd_ctrl)
  if (is.null(df)) {
    E[case, ] <- stats::rnorm(length(case) * p, mean = 0, sd = sd_case)
  } else {
    # Heavier-tailed (Student-t) case profiles, scaled to sd_case.
    E[case, ] <- stats::rt(length(case) * p, df = df) *
      (sd_case / sqrt(df / (df - 2)))
  }
  X <- exp(base + E)
  list(X = X, y = y, features = features)
}

.auc_or_wilcoxon_lot <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Standalone per-sample rCLR + its present support (matches the scorer's frozen
# transform) for the manual CDT / LOT-isometry cross-checks.
.rclr0_lot <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}
# Present (v > 0) coordinates -- only structural zeros excluded, mirroring the
# scorer's .ot_lot_rclr_pos (a present coordinate that equals the geometric mean,
# z == 0, is kept).
.rclr_pos_lot <- function(v) {
  z <- .rclr0_lot(v)
  z[v > 0]
}
# type = 1 (empirical inverse-CDF) = the exact 1-D OT quantile map, mirroring
# .ot_lot_cdt; the squared L2 between two such embeddings is the exact 1-D W2^2.
.cdt_q_lot <- function(zp, t) {
  stats::quantile(zp, probs = t, type = 1, names = FALSE)
}
.grid_lot <- function(K) (seq_len(K) - 0.5) / K

test_that("ot-lot fit/score roundtrip has the right shape and types", {
  dat <- .make_ot_lot_data()
  model <- fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 16L))
  expect_s3_class(model, "ot_lot_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  K <- 16L
  expect_length(model$t_k, K)
  expect_length(model$Q_ref, K)
  expect_length(model$head$w, K)
  expect_length(model$head$b, 1L)
  expect_equal(model$t_k, (seq_len(K) - 0.5) / K)
  # Frozen tangent precision is symmetric, the right size and finite.
  expect_equal(dim(model$head$Sinv), c(K, K))
  expect_equal(model$head$Sinv, t(model$head$Sinv), tolerance = 1e-10)
  expect_true(all(is.finite(model$head$w)) && is.finite(model$head$b))
  expect_equal(model$n_usable_case, sum(dat$y == 1L))
  expect_equal(model$n_usable_control, sum(dat$y == 0L))

  score <- score_ot_lot(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("ot-lot separates a planted distributional-shape difference", {
  dat <- .make_ot_lot_data(n = 300L, seed = 202L)
  ctrl <- which(dat$y == 0L)
  case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_ot_lot(dat$X[train, ], dat$y[train])
  score <- score_ot_lot(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_lot(dat$y[test], score), 0.7)

  # A heavier-tailed (Student-t) case generator is also separated.
  dat_t <- .make_ot_lot_data(n = 300L, sd_case = 0.45, df = 3, seed = 303L)
  train_t <- c(utils::head(which(dat_t$y == 0L), 112L),
               utils::head(which(dat_t$y == 1L), 112L))
  test_t <- c(utils::tail(which(dat_t$y == 0L), 38L),
              utils::tail(which(dat_t$y == 1L), 38L))
  model_t <- fit_ot_lot(dat_t$X[train_t, ], dat_t$y[train_t])
  score_t <- score_ot_lot(model_t, dat_t$X[test_t, ])
  expect_gt(.auc_or_wilcoxon_lot(dat_t$y[test_t], score_t), 0.7)
})

test_that("ot-lot passes the canonical row-equivariance gate", {
  dat <- .make_ot_lot_data(seed = 43L)
  model <- fit_ot_lot(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:200, ]
  score_fun <- function(model, X, meta) score_ot_lot(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("ot-lot scoring is exactly invariant to per-sample scaling", {
  dat <- .make_ot_lot_data(seed = 44L)
  model <- fit_ot_lot(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]
  expect_equal(score_ot_lot(model, row * 7),
               score_ot_lot(model, row), tolerance = 1e-8)
  expect_equal(score_ot_lot(model, row * 1e6),
               score_ot_lot(model, row), tolerance = 1e-8)
  # Mixed per-row scaling: each row scaled by its own factor leaves scores fixed.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch * c(1, 7, 100, 1e-3, 1e6, 0.5)
  expect_equal(score_ot_lot(model, scaled),
               score_ot_lot(model, batch), tolerance = 1e-8)
  # Row-equivariant batch: scaling some rows leaves the others' scores fixed.
  scaled2 <- batch
  scaled2[c(2L, 5L), ] <- scaled2[c(2L, 5L), ] * 1e6
  expect_equal(score_ot_lot(model, scaled2),
               score_ot_lot(model, batch), tolerance = 1e-8)
  # A row with mostly zeros (rCLR over its own nonzero parts) still scores finite.
  X_sparse <- dat$X[1:3, , drop = FALSE]
  X_sparse[1, 8:ncol(X_sparse)] <- 0
  expect_true(all(is.finite(score_ot_lot(model, X_sparse))))
})

test_that("ot-lot handles partial feature overlap (consistent + single-row)", {
  dat <- .make_ot_lot_data(seed = 45L)
  model <- fit_ot_lot(dat$X, dat$y)
  # Subset and REORDER the columns: present = intersect(universe, colnames) keeps
  # feature-universe order, and the CDT sorts values, so the score is invariant to
  # the X column order (feature-exchangeable by design).
  keep <- rev(c(dat$features[1:12], dat$features[20:30]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_ot_lot(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_ot_lot(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_ot_lot(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("ot-lot matches the canonical dispatch", {
  dat <- .make_ot_lot_data(seed = 46L)
  model <- fit_ot_lot(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:200, ]
  direct <- score_ot_lot(model, X_test)

  # Post-integration the manifest carries the ot-lot row and score_ot_lot lives in
  # the OmicSelector namespace, so the canonical adapter resolves it directly (the
  # real roster -> package dispatch path). In the pre-integration staged self-test
  # score_ot_lot is not yet in the namespace, so we clone the roster row if absent
  # and register a thin canonical adapter to exercise the same dispatch wiring.
  # Either way the dispatched score must equal the direct call.
  roster <- singlesample_method_roster()
  if (!"ot-lot" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "ot-lot"
    tmpl$fit_fn <- "fit_ot_lot"
    tmpl$score_fn <- "score_ot_lot"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_ot_lot", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "ot-lot",
      function(model, X, meta) score_ot_lot(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("ot-lot", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("ot-lot returns the neutral 0 below the feature-overlap floor", {
  dat <- .make_ot_lot_data(seed = 47L)
  model <- fit_ot_lot(dat$X, dat$y, hp = list(min_features = 5L))
  # Only three universe features present -> below min_features=5 -> neutral 0.
  X_three <- dat$X[1:6, colnames(dat$X)[1:3], drop = FALSE]
  expect_equal(score_ot_lot(model, X_three), rep(0, nrow(X_three)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_ot_lot(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
  # All universe columns present but a specimen with too FEW nonzero coordinates
  # (a near-empty profile) -> its own row scores the neutral 0.
  X_sparse <- dat$X[1:4, , drop = FALSE]
  X_sparse[1, ] <- 0
  X_sparse[1, 1:3] <- dat$X[1, 1:3]      # only 3 nonzero coords < min_features 5
  s_sparse <- score_ot_lot(model, X_sparse)
  expect_equal(s_sparse[1], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s_sparse)))
})

test_that("ot-lot fitting/scoring is deterministic and consumes no RNG", {
  dat <- .make_ot_lot_data(seed = 48L)
  hp <- list(n_quantiles = 20L, reference = "pooled")

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_ot_lot(dat$X, dat$y, hp = hp)
  # A fit must not touch the global RNG (no sample()/rnorm() anywhere).
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_ot_lot(dat$X, dat$y, hp = hp)
  expect_identical(model_1, model_2)

  s1 <- score_ot_lot(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_ot_lot(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("ot-lot fitting leaves no RNG state in a fresh session", {
  # With no prior global seed, a fit must NOT create one (it consumes no RNG).
  dat <- .make_ot_lot_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  model <- fit_ot_lot(dat$X, dat$y, hp = list(seed = 29L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "ot_lot_model")
})

test_that("ot-lot tangent precision stays PD when K exceeds the class count", {
  # Few usable specimens per class (8) but K = 16 quantile components -> the raw
  # pooled tangent covariance is rank-deficient, and the monotone quantile vector
  # makes adjacent components highly correlated. The fixed diagonal load + the
  # ridge-to-PD loop must keep the precision positive definite and every score
  # finite.
  dat <- .make_ot_lot_data(n = 16L, p = 40L, seed = 71L)
  model <- fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 16L))
  expect_equal(model$n_usable_case, 8L)
  expect_equal(model$n_usable_control, 8L)
  expect_true(all(is.finite(model$head$Sinv)))
  ev <- eigen(model$head$Sinv, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))                       # precision PD <=> covariance PD
  score <- score_ot_lot(model, dat$X)
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("ot-lot pooled reference and barycenter give identical scores", {
  # Reference-origin invariance: the LDA decision value is algebraically
  # independent of the tangent base point, so the barycenter and pooled references
  # (different Q_ref vectors) must yield the SAME scores from the same training
  # data.
  dat <- .make_ot_lot_data(seed = 52L)
  m_bary <- fit_ot_lot(dat$X, dat$y, hp = list(reference = "barycenter"))
  m_pool <- fit_ot_lot(dat$X, dat$y, hp = list(reference = "pooled"))
  expect_false(isTRUE(all.equal(m_bary$Q_ref, m_pool$Q_ref)))   # bases differ
  expect_equal(score_ot_lot(m_bary, dat$X),
               score_ot_lot(m_pool, dat$X), tolerance = 1e-8)
})

test_that("ot-lot hyperparameter validation is strict", {
  dat <- .make_ot_lot_data(seed = 50L)
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(reference = "banana")),
               "hp\\$reference")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 1.5)),
               "n_quantiles")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 3L)),
               "n_quantiles")
  # Out-of-integer-range must hit the explicit error, not an as.integer overflow.
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 3e9)),
               "n_quantiles")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(min_features = 1)),
               "min_features")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(shrink = -0.1)),
               "hp\\$shrink")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(seed = -1)), "hp\\$seed")
  # `$` partial matching would let "min" bind to "min_features"; the exact [[ ]]
  # reads + allowed-list reject the unknown field instead.
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(min = 5L)), "unknown hp")
  # Fully-unnamed hp (names(hp) == NULL) must be rejected, not silently ignored:
  # setdiff(NULL, allowed) is empty, so the allowed-list check alone would pass it.
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(24L)), "must be named")
  # A partially-unnamed entry is rejected too.
  expect_error(fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 8L, 99)),
               "must be named")
  # Duplicated field names must be rejected, not silently last-wins.
  expect_error(
    fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 8L, n_quantiles = 16L)),
    "duplicated hp")
})

test_that("ot-lot keeps present coordinates at the geometric mean (no z!=0 drop)", {
  # Structural-zero exclusion must key on v > 0 (raw structural zeros), NOT z != 0:
  # a genuinely present coordinate whose centred-log value lands exactly at the
  # geometric mean (z == 0) must be kept, and a fully-constant positive profile
  # must yield a point mass at 0 of size = #positive features (not an empty vector
  # that silently scores neutral 0).
  # v = (1, 2, 4): geomean = 2, so the middle present coordinate has rCLR exactly 0.
  zp <- .rclr_pos_lot(c(1, 2, 4))
  expect_length(zp, 3L)                            # all three present coords kept
  expect_equal(sort(zp), c(-log(2), 0, log(2)), tolerance = 1e-12)
  # Structural zeros (v == 0) ARE excluded; the present middle-0 is still kept.
  zp2 <- .rclr_pos_lot(c(0, 1, 2, 4, 0))
  expect_length(zp2, 3L)
  expect_equal(sort(zp2), c(-log(2), 0, log(2)), tolerance = 1e-12)
  # A constant positive profile -> point mass at 0 of size = #positive features.
  expect_length(.rclr_pos_lot(rep(2, 5)), 5L)
  expect_true(all(.rclr_pos_lot(rep(2, 5)) == 0))
  # End to end: a specimen that is constant across all (>= min_features) present
  # features scores finite (not dropped to neutral 0 because z != 0 emptied it).
  dat <- .make_ot_lot_data(seed = 51L)
  model <- fit_ot_lot(dat$X, dat$y, hp = list(min_features = 5L))
  Xc <- dat$X[1:3, , drop = FALSE]
  Xc[1, ] <- 5                                     # constant positive profile (m = p)
  s <- score_ot_lot(model, Xc)
  expect_true(all(is.finite(s)))
})

test_that("ot-lot CDT is the exact 1-D OT map and is LOT-isometric to W2", {
  # (a) The 1-D optimal-transport map is the sorted (monotone) rearrangement, and
  # the CDT uses the empirical inverse-CDF (type = 1), so the quantile function is
  # invariant to input order and -- at K = m interior levels -- returns exactly the
  # sorted support. For zp = {0, 1, 2, 3} on grid t = (k - 0.5)/4 the type-1
  # quantiles are the order statistics themselves: 0, 1, 2, 3.
  t4 <- .grid_lot(4L)
  hand <- c(0, 1, 2, 3)
  expect_equal(.cdt_q_lot(c(3, 1, 0, 2), t4), hand, tolerance = 1e-12)
  expect_equal(.cdt_q_lot(c(3, 1, 0, 2), t4), .cdt_q_lot(sort(c(3, 1, 0, 2)), t4),
               tolerance = 1e-12)

  # (b) GENUINE (non-circular) LOT isometry. For two equal-size empirical measures
  # the exact 1-D squared 2-Wasserstein distance is mean((sort(z_i) - sort(z_j))^2)
  # (sorted-support matching). With the empirical inverse-CDF (type = 1), sampling
  # at K = m interior levels returns the sorted support exactly, so the squared L2
  # between the two CDT embeddings EQUALS that exact empirical W2^2 to machine
  # precision -- the defining isometry. The earlier draft compared type-7 at K=64
  # vs type-7 at K=8192 (same interpolating helper, so it only measured type-7
  # quadrature stability, NOT the W2 isometry); type-7 converges to a different,
  # low-biased limit and would fail the assertions below.
  dat <- .make_ot_lot_data(seed = 61L)
  zp_i <- .rclr_pos_lot(dat$X[1, ])
  zp_j <- .rclr_pos_lot(dat$X[151, ])
  m <- ncol(dat$X)
  expect_length(zp_i, m)                          # full (all-present) profiles, equal m
  expect_length(zp_j, m)
  # type=1 at K=m IS the sorted support.
  expect_equal(.cdt_q_lot(zp_i, .grid_lot(m)), sort(zp_i), tolerance = 1e-12)
  emp_w2 <- mean((sort(zp_i) - sort(zp_j))^2)      # exact 1-D empirical W2^2
  lot_w2 <- mean((.cdt_q_lot(zp_i, .grid_lot(m)) -
                    .cdt_q_lot(zp_j, .grid_lot(m)))^2)
  expect_equal(lot_w2, emp_w2, tolerance = 1e-10)  # exact isometry (not approx)
  # A fine grid (K >> m) still converges to the same exact empirical W2.
  fine <- .grid_lot(20000L)
  lot_fine <- mean((.cdt_q_lot(zp_i, fine) - .cdt_q_lot(zp_j, fine))^2)
  expect_equal(lot_fine, emp_w2, tolerance = 0.02)
  # Regression guard: a type-7 (piecewise-linear) embedding is materially different
  # from the exact empirical W2 (here ~12% low), so reverting .ot_lot_cdt to type=7
  # would break the machine-precision isometry above.
  t7 <- function(z, t) stats::quantile(z, t, type = 7, names = FALSE)
  lot7_fine <- mean((t7(zp_i, fine) - t7(zp_j, fine))^2)
  expect_false(isTRUE(all.equal(lot7_fine, emp_w2, tolerance = 1e-3)))

  # (c) The score equals the hand-computed linear decision value on the frozen
  # head, score = w . (Q - Q_ref) + b, AND equals the textbook LDA decision value
  # w . Q - 0.5 w . (Qbar_case + Qbar_control) (the Q_ref origin cancels).
  model <- fit_ot_lot(dat$X, dat$y, hp = list(n_quantiles = 16L))
  w <- model$head$w
  b <- model$head$b
  Q_ref <- model$Q_ref
  t_k <- model$t_k
  Qbar_case <- model$head$mu_case + Q_ref        # back to Q-space class means
  Qbar_control <- model$head$mu_control + Q_ref
  for (i in c(1L, 7L, 60L, 151L, 199L)) {
    Q <- .cdt_q_lot(.rclr_pos_lot(dat$X[i, ]), t_k)
    manual <- sum(w * (Q - Q_ref)) + b
    lda <- sum(w * Q) - 0.5 * sum(w * (Qbar_case + Qbar_control))
    got <- score_ot_lot(model, dat$X[i, , drop = FALSE])
    expect_equal(got, manual, tolerance = 1e-9)
    expect_equal(got, lda, tolerance = 1e-9)
  }
})

test_that("ot-lot rejects degenerate training data explicitly", {
  # Too few usable specimens per class (all-but-one specimen sparse below the
  # min-features floor) -> the fit must stop, not return a rank-0 head.
  dat <- .make_ot_lot_data(n = 20L, p = 40L, seed = 55L)
  X <- dat$X
  # Knock all but one specimen per class down to a single nonzero coordinate.
  keep_full <- c(which(dat$y == 0L)[1], which(dat$y == 1L)[1])
  for (i in setdiff(seq_len(nrow(X)), keep_full)) {
    X[i, ] <- 0
    X[i, 1] <- dat$X[i, 1]
  }
  expect_error(fit_ot_lot(X, dat$y, hp = list(min_features = 5L)),
               "usable specimens per class")
})
