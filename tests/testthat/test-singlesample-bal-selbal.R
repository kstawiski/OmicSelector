library(testthat)

# Runs inside the package test harness (devtools::test / R CMD check), where all
# fit_bal_selbal/score_bal_selbal and package helpers are already loaded.

skip_if_no_selbal <- function() {
  testthat::skip_if_not_installed("selbal")
}

# Planted single-balance data: a known numerator group is elevated in cases and
# a known denominator group depressed (otherwise identical baseline). Moderate
# multipliers so selbal's accuracy measure does not instantly saturate.
.make_bal_selbal_data <- function(n = 140L, p = 16L, seed = 7L, mult = 1.8) {
  set.seed(seed)
  feat <- paste0("f", sprintf("%02d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 30, rate = 2), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), feat))
  y <- rep(c(0, 1), each = n / 2L)
  num_pl <- c("f01", "f02", "f03")
  den_pl <- c("f14", "f15", "f16")
  X[y == 1L, num_pl] <- X[y == 1L, num_pl] * mult
  X[y == 1L, den_pl] <- X[y == 1L, den_pl] / mult
  list(X = X, y = y, num_pl = num_pl, den_pl = den_pl, feat = feat)
}

.auc_mw <- function(y, s) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, s, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(s, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("bal-selbal fit/score returns finite numeric vector of right shape", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data()
  model <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                          hp = list(logit_acc = "Dev"))
  score <- score_bal_selbal(model, dat$X[101:140, ])

  expect_s3_class(model, "bal_selbal_model")
  expect_type(model$numerator, "character")
  expect_type(model$denominator, "character")
  expect_true(length(model$numerator) >= 1L)
  expect_true(length(model$denominator) >= 1L)
  expect_length(intersect(model$numerator, model$denominator), 0L)
  expect_type(score, "double")
  expect_length(score, 40L)
  expect_true(all(is.finite(score)))
})

test_that("bal-selbal selects planted groups and separates classes", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data()
  tr <- c(1:50, 71:120)
  te <- setdiff(seq_len(nrow(dat$X)), tr)
  model <- fit_bal_selbal(dat$X[tr, ], dat$y[tr], hp = list(logit_acc = "Dev"))

  # selbal may orient the planted groups onto either side; count overlap with
  # both orientations (a swapped balance is the same balance up to sign, which
  # the logistic calibration absorbs).
  same <- length(intersect(model$numerator, dat$num_pl)) +
    length(intersect(model$denominator, dat$den_pl))
  swap <- length(intersect(model$numerator, dat$den_pl)) +
    length(intersect(model$denominator, dat$num_pl))
  expect_gte(max(same, swap), 3L)

  score <- score_bal_selbal(model, dat$X[te, ])
  expect_gt(mean(score[dat$y[te] == 1]), mean(score[dat$y[te] == 0]))
  expect_gt(.auc_mw(dat$y[te], score), 0.9)
})

test_that("bal-selbal balance is exactly invariant to per-sample positive scaling", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 11L)
  model <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                          hp = list(logit_acc = "Dev"))
  X_test <- dat$X[101:140, ]
  s_base <- score_bal_selbal(model, X_test)

  # constant scale
  expect_equal(score_bal_selbal(model, X_test * 13.7), s_base,
               tolerance = 1e-6)
  # per-row random positive scale
  set.seed(5)
  scal <- stats::runif(nrow(X_test), 0.05, 40)
  expect_equal(score_bal_selbal(model, X_test * scal), s_base,
               tolerance = 1e-6)
})

test_that("bal-selbal passes the canonical row-equivariance gate", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 12L)
  model <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                          hp = list(logit_acc = "Dev"))
  X_test <- dat$X[101:140, ]
  score_fun <- function(model, X, meta) score_bal_selbal(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("bal-selbal is column-permutation invariant and single-row == batch", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 13L)
  model <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                          hp = list(logit_acc = "Dev"))
  X_test <- dat$X[101:140, ]
  s_base <- score_bal_selbal(model, X_test)

  set.seed(9)
  perm <- sample(ncol(X_test))
  expect_equal(score_bal_selbal(model, X_test[, perm, drop = FALSE]), s_base,
               tolerance = 1e-10)
  # single-row == its position in the batch
  expect_equal(score_bal_selbal(model, X_test[7, , drop = FALSE]),
               s_base[7], tolerance = 1e-12)
})

test_that("bal-selbal handles partial frozen-group overlap consistently", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 14L)
  model <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                          hp = list(logit_acc = "Dev"))
  # drop ONE numerator and ONE denominator feature (keep >= 1 each side so the
  # balance is over the present subset, not neutral)
  drop_one <- c(model$numerator[1], model$denominator[1])
  keep <- setdiff(colnames(dat$X), drop_one)
  X_partial <- dat$X[101:112, keep, drop = FALSE]

  s_partial <- score_bal_selbal(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # single-row consistency under partial overlap (no batch coupling)
  expect_equal(score_bal_selbal(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
})

test_that("bal-selbal returns neutral 0 below coverage and on degenerate selection", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 15L)
  model <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                          hp = list(logit_acc = "Dev"))

  # no shared features at all
  X_none <- matrix(1, nrow = 4L, ncol = 1L,
                   dimnames = list(NULL, "unrelated"))
  expect_equal(score_bal_selbal(model, X_none), rep(0, 4L))

  # only the numerator side present (denominator absent) -> neutral 0
  X_numonly <- dat$X[101:104, model$numerator, drop = FALSE]
  expect_equal(score_bal_selbal(model, X_numonly), rep(0, 4L))

  # coverage floor too high to ever be met -> all 0
  model_hi <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                             hp = list(logit_acc = "Dev",
                                       min_group_coverage = 999L))
  expect_equal(score_bal_selbal(model_hi, dat$X[101:110, ]), rep(0, 10L))

  # degenerate (empty) selection -> all 0
  model_degen <- model
  model_degen$degenerate <- TRUE
  expect_equal(score_bal_selbal(model_degen, dat$X[101:110, ]), rep(0, 10L))

  # an all-zero specimen row is finite neutral 0, not NaN/NA
  X_zero <- dat$X[101:103, , drop = FALSE]
  X_zero[1, ] <- 0
  s_zero <- score_bal_selbal(model, X_zero)
  expect_true(all(is.finite(s_zero)))
  expect_equal(s_zero[1], 0)
})

test_that("bal-selbal refit is deterministic and leaves the global RNG untouched", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 16L)

  set.seed(2024)
  before <- get(".Random.seed", envir = globalenv())
  m1 <- fit_bal_selbal(dat$X, dat$y, hp = list(logit_acc = "Dev"))
  after_fit <- get(".Random.seed", envir = globalenv())
  s1 <- score_bal_selbal(m1, dat$X[1:20, ])
  after_score <- get(".Random.seed", envir = globalenv())
  m2 <- fit_bal_selbal(dat$X, dat$y, hp = list(logit_acc = "Dev"))
  s2 <- score_bal_selbal(m2, dat$X[1:20, ])

  expect_identical(m1$numerator, m2$numerator)
  expect_identical(m1$denominator, m2$denominator)
  expect_identical(m1$intercept, m2$intercept)
  expect_identical(m1$slope, m2$slope)
  expect_identical(s1, s2)
  # RNG untouched by fit AND score
  expect_identical(before, after_fit)
  expect_identical(after_fit, after_score)

  # fit with NO pre-existing seed leaves none
  rm(list = ".Random.seed", envir = globalenv())
  m3 <- fit_bal_selbal(dat$X, dat$y, hp = list(logit_acc = "Dev"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m3, "bal_selbal_model")
})

test_that("bal-selbal is deployable through the canonical dispatch", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 17L)
  model <- fit_bal_selbal(dat$X[1:100, ], dat$y[1:100],
                          hp = list(logit_acc = "Dev"))
  X_test <- dat$X[101:140, ]

  # Build a minimal one-row roster so singlesample_score_call resolves
  # bal-selbal without the (not-yet-integrated) package manifest, and register
  # a canonical adapter pointing at score_bal_selbal. The dispatched score must
  # equal the direct call.
  roster <- data.frame(
    method_id = "bal-selbal", family = "K", estimand = "within",
    role = "discriminator", tier = "R1", dep_route = "selbal-github",
    fit_fn = "fit_bal_selbal", score_fn = "score_bal_selbal",
    pkg_status = "new", notes = "", row_source = "per_cohort_rows",
    lopbo_mechanism = "within_cohort", stringsAsFactors = FALSE)
  singlesample_register_score_adapter(
    "bal-selbal", function(model, X, meta) score_bal_selbal(model, X, meta))

  expect_equal(
    singlesample_score_call("bal-selbal", model, X_test, roster = roster),
    score_bal_selbal(model, X_test),
    tolerance = 1e-12)
})

test_that("bal-selbal strict hp resolver rejects malformed hp", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 18L)

  expect_error(fit_bal_selbal(dat$X, dat$y, hp = "notalist"),
               "hp must be a list")
  # unknown field
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(maxV = 5L)),
               "unknown hp")
  # unnamed field
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(10L)),
               "must be named")
  # duplicate field
  expect_error(
    fit_bal_selbal(dat$X, dat$y,
                   hp = stats::setNames(list(5L, 6L), c("max_vars", "max_vars"))),
    "duplicate hp")
  # out-of-range max_vars
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(max_vars = 1L)),
               "max_vars")
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(max_vars = 2.5)),
               "max_vars")
  # bad logit_acc
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(logit_acc = "bad")),
               "logit_acc")
  # bad zero_rep
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(zero_rep = "nope")),
               "zero_rep")
  # bad min_group_coverage
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(min_group_coverage = 0L)),
               "min_group_coverage")
  # bad seed
  expect_error(fit_bal_selbal(dat$X, dat$y, hp = list(seed = -1L)),
               "seed")
})

test_that("bal-selbal input validation errors are explicit", {
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 19L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_bal_selbal(X_unnamed, dat$y), "feature names")
  expect_error(fit_bal_selbal(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_bal_selbal(dat$X, rep(0, nrow(dat$X))), "at least one case")
  expect_error(fit_bal_selbal(dat$X[, 1, drop = FALSE], dat$y),
               "at least two features")
  # score guards on model class
  expect_error(score_bal_selbal(list(), dat$X), "class bal_selbal_model")
})

test_that("bal-selbal frozen default logit_acc is 'Dev' and the default config discriminates", {
  # The shipped default must equal the validated/benchmarked config: deviance keeps
  # adding informative features after class separation (AUC saturates early at a
  # tiny balance). This exercises the DEFAULT path (no logit_acc passed).
  skip_if_no_selbal()
  dat <- .make_bal_selbal_data(seed = 23L)
  tr <- c(1:50, 71:120)          # balanced: 50 controls + 50 cases (rows are y-ordered)
  te <- setdiff(seq_len(nrow(dat$X)), tr)
  model <- fit_bal_selbal(dat$X[tr, ], dat$y[tr])          # default config
  expect_identical(model$hp$logit_acc, "Dev")
  s <- score_bal_selbal(model, dat$X[te, ])
  expect_true(all(is.finite(s)))
  # held-out Mann-Whitney AUC (independent of scorer internals)
  yte <- dat$y[te]; r <- rank(s); n1 <- sum(yte == 1); n0 <- sum(yte == 0)
  auc <- (sum(r[yte == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  expect_gt(auc, 0.7)
})
