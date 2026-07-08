library(testthat)

# Runs inside the package test harness (devtools::test / R CMD check), where
# fit_sel_stablemate/score_sel_stablemate and all package helpers are already
# loaded. StableMate is a GitHub dependency, so every test skips when it is
# not installed.

skip_if_no_stablemate <- function() {
  testthat::skip_if_not_installed("StableMate")
}

# Planted multi-environment data. Stable features (miR-1, miR-2, miR-3) are
# elevated in cases CONSISTENTLY across all environments (multiplicative shift,
# so the signal lives in the per-sample rCLR). Spurious features (miR-4, miR-5,
# miR-6) are elevated in cases in ONE environment each only -- predictive within
# that cohort but unstable across cohorts. The remaining features are noise.
.make_sel_stablemate_data <- function(seed = 3L, n_per = 70L,
                                      envs = c("e1", "e2", "e3"),
                                      stable_mult = 2.2, spur_mult = 3.0) {
  set.seed(seed)
  p <- 12L
  n <- n_per * length(envs)
  env <- rep(envs, each = n_per)
  y <- rep(rep(c(0, 1), each = n_per / 2L), length(envs))
  X <- matrix(stats::rgamma(n * p, shape = 8, rate = 1), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("miR-", seq_len(p))))
  for (f in 1:3) X[y == 1L, f] <- X[y == 1L, f] * stable_mult
  spur_idx <- c(4L, 5L, 6L)
  for (j in seq_along(envs)) {
    if (j > length(spur_idx)) break
    rows <- env == envs[j] & y == 1L
    X[rows, spur_idx[j]] <- X[rows, spur_idx[j]] * spur_mult
  }
  list(X = X, y = y, env = env,
       meta = data.frame(accession = env, stringsAsFactors = FALSE),
       stable = paste0("miR-", 1:3), spurious = paste0("miR-", 4:6))
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

test_that("sel-stablemate fit/score returns finite numeric vector of right shape", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data()
  model <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                              hp = list(K = 60L, sigthresh = 0.8))
  score <- score_sel_stablemate(model, dat$X)

  expect_s3_class(model, "sel_stablemate_model")
  expect_type(model$selected_features, "character")
  expect_type(model$coefficients, "double")
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("sel-stablemate selects stable predictors over spurious ones and generalizes to a held-out environment", {
  skip_if_no_stablemate()
  # Train on e1 + e2, hold out e3 entirely (never seen at fit). The stable set
  # must dominate the selection, and the frozen logistic must discriminate the
  # held-out environment.
  tr <- .make_sel_stablemate_data(seed = 3L, envs = c("e1", "e2"))
  te <- .make_sel_stablemate_data(seed = 99L, envs = c("e3"))
  model <- fit_sel_stablemate(tr$X, tr$y, meta_train = tr$meta,
                              hp = list(K = 80L, sigthresh = 0.8))

  n_stable <- length(intersect(model$selected_features, tr$stable))
  n_spur <- length(intersect(model$selected_features, tr$spurious))
  expect_equal(model$n_environments, 2L)
  expect_gte(n_stable, 2L)            # recovers the stable set
  expect_gt(n_stable, n_spur)         # stable overlap dominates spurious

  score_te <- score_sel_stablemate(model, te$X)
  expect_gt(mean(score_te[te$y == 1L]), mean(score_te[te$y == 0L]))
  expect_gt(.auc_mw(te$y, score_te), 0.8)   # generalizes to held-out env
})

test_that("sel-stablemate selects NOTHING when no predictor is stable across environments", {
  # A TRANSFER/stability method must NOT fall back to predictivity-only,
  # environment-specific features when the stability ensemble finds no stable
  # predictor: that would freeze spurious batch-driven signal while reporting
  # n_environments > 1. The honest outcome is an empty selection -> degenerate
  # (neutral-0) model. (Regression for the R1 "stab-empty fallback to pred" defect.)
  skip_if_no_stablemate()
  set.seed(7L)
  p <- 10L; feat <- paste0("miR-", seq_len(p))
  mk <- function(env) {
    n <- 80L
    X <- matrix(stats::rlnorm(n * p, log(40), 0.5), n, p, dimnames = list(NULL, feat))
    y <- rep(c(0L, 1L), each = n / 2L)
    if (env == "e1") X[y == 1L, 1] <- X[y == 1L, 1] * 2.5  # miR-1 case-shifted in e1 only
    if (env == "e2") X[y == 1L, 2] <- X[y == 1L, 2] * 2.5  # miR-2 case-shifted in e2 only
    list(X = X, y = y)
  }
  e1 <- mk("e1"); e2 <- mk("e2")
  X <- rbind(e1$X, e2$X); y <- c(e1$y, e2$y)
  meta <- data.frame(accession = c(rep("e1", 80L), rep("e2", 80L)),
                     stringsAsFactors = FALSE)
  model <- fit_sel_stablemate(X, y, meta_train = meta,
                              hp = list(K = 100L, max_features = 10L,
                                        sigthresh = 0.9, seed = 1L))
  expect_equal(model$n_environments, 2L)        # the env-aware path DID engage
  expect_length(model$selected_features, 0L)    # but found NO stable predictor
  expect_true(model$degenerate)
  expect_true(all(score_sel_stablemate(model, X) == 0))  # neutral-0, not spurious
})

test_that("sel-stablemate score equals the frozen logistic linear predictor over selected rCLR features", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 5L)
  model <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                              hp = list(K = 60L, sigthresh = 0.8))
  s <- score_sel_stablemate(model, dat$X)

  # Independent recompute: per-sample rCLR over the frozen universe, restricted
  # to the selected features, frozen intercept + sum(beta * rclr).
  rclr <- function(v) {
    z <- rep(0, length(v)); pos <- which(v > 0)
    if (length(pos)) { lv <- log(v[pos]); z[pos] <- lv - mean(lv) }
    z
  }
  manual <- vapply(seq_len(nrow(dat$X)), function(i) {
    xi <- dat$X[i, ]; up <- intersect(model$feature_universe, names(xi))
    v <- xi[up]; z <- rclr(v); names(z) <- up
    use <- intersect(model$selected_features, up)
    if (length(use) < model$min_selected) return(0)
    model$intercept + sum(model$coefficients[use] * z[use])
  }, numeric(1))
  expect_equal(s, manual, tolerance = 1e-8)
})

test_that("sel-stablemate is exactly invariant to per-sample positive scaling", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 11L)
  model <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                              hp = list(K = 60L, sigthresh = 0.8))
  X_test <- dat$X[1:60, , drop = FALSE]
  s_base <- score_sel_stablemate(model, X_test)

  expect_equal(score_sel_stablemate(model, X_test * 13.7), s_base,
               tolerance = 1e-6)
  set.seed(5)
  scal <- stats::runif(nrow(X_test), 0.05, 40)
  expect_equal(score_sel_stablemate(model, X_test * scal), s_base,
               tolerance = 1e-6)
})

test_that("sel-stablemate passes the canonical row-equivariance gate", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 12L)
  model <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                              hp = list(K = 60L, sigthresh = 0.8))
  X_test <- dat$X[1:50, , drop = FALSE]
  score_fun <- function(model, X, meta) score_sel_stablemate(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("sel-stablemate is column-permutation invariant and single-row == batch", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 13L)
  model <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                              hp = list(K = 60L, sigthresh = 0.8))
  X_test <- dat$X[1:50, , drop = FALSE]
  s_base <- score_sel_stablemate(model, X_test)

  set.seed(9)
  perm <- sample(ncol(X_test))
  expect_equal(score_sel_stablemate(model, X_test[, perm, drop = FALSE]), s_base,
               tolerance = 1e-10)
  expect_equal(score_sel_stablemate(model, X_test[7, , drop = FALSE]),
               s_base[7], tolerance = 1e-12)
})

test_that("sel-stablemate reports n_environments and degenerates gracefully without usable cohorts", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 7L)

  # multi-environment -> engaged
  m_multi <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                                hp = list(K = 60L, sigthresh = 0.8))
  expect_equal(m_multi$n_environments, 3L)

  # NULL meta -> degenerate predictivity-only selection, n_environments == 1
  m_null <- fit_sel_stablemate(dat$X, dat$y, meta_train = NULL,
                               hp = list(K = 60L, sigthresh = 0.8))
  expect_equal(m_null$n_environments, 1L)
  expect_gte(length(intersect(m_null$selected_features, dat$stable)), 2L)
  expect_true(all(is.finite(score_sel_stablemate(m_null, dat$X))))

  # single cohort -> degenerate, n_environments == 1
  one_meta <- data.frame(accession = rep("one", nrow(dat$X)),
                         stringsAsFactors = FALSE)
  m_one <- fit_sel_stablemate(dat$X, dat$y, meta_train = one_meta,
                              hp = list(K = 60L, sigthresh = 0.8))
  expect_equal(m_one$n_environments, 1L)

  # missing cohort column -> degenerate, n_environments == 1
  m_miss <- fit_sel_stablemate(dat$X, dat$y,
                               meta_train = data.frame(other = dat$env),
                               hp = list(K = 60L, sigthresh = 0.8))
  expect_equal(m_miss$n_environments, 1L)
})

test_that("sel-stablemate returns neutral 0 below min_selected and on empty/degenerate selection", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 15L)
  model <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                              hp = list(K = 60L, sigthresh = 0.8))

  # no selected features present at all -> neutral 0
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(NULL, c("zzz-a", "zzz-b")))
  expect_equal(score_sel_stablemate(model, X_none), rep(0, 4L))

  # below min_selected (require 2, present only 1) -> neutral 0
  keep1 <- model$selected_features[1]
  X1 <- dat$X[1:5, keep1, drop = FALSE]
  model_hi <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                                 hp = list(K = 60L, sigthresh = 0.8,
                                           min_selected = 2L))
  expect_equal(score_sel_stablemate(model_hi, X1), rep(0, 5L))

  # degenerate (empty) selection -> all 0
  model_degen <- model
  model_degen$degenerate <- TRUE
  expect_equal(score_sel_stablemate(model_degen, dat$X[1:10, ]), rep(0, 10L))

  # an all-zero specimen row is finite neutral, not NaN/NA
  X_zero <- dat$X[1:3, , drop = FALSE]
  X_zero[1, ] <- 0
  s_zero <- score_sel_stablemate(model, X_zero)
  expect_true(all(is.finite(s_zero)))
})

test_that("sel-stablemate refit is deterministic and leaves the global RNG untouched", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 16L)

  set.seed(2024)
  before <- get(".Random.seed", envir = globalenv())
  m1 <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                           hp = list(K = 60L, sigthresh = 0.8))
  after_fit <- get(".Random.seed", envir = globalenv())
  s1 <- score_sel_stablemate(m1, dat$X[1:20, ])
  after_score <- get(".Random.seed", envir = globalenv())
  m2 <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                           hp = list(K = 60L, sigthresh = 0.8))
  s2 <- score_sel_stablemate(m2, dat$X[1:20, ])

  expect_identical(m1$selected_features, m2$selected_features)
  expect_identical(m1$intercept, m2$intercept)
  expect_identical(m1$coefficients, m2$coefficients)
  expect_identical(s1, s2)
  expect_identical(before, after_fit)      # RNG untouched by fit
  expect_identical(after_fit, after_score) # RNG untouched by score

  # fit with NO pre-existing seed leaves none
  rm(list = ".Random.seed", envir = globalenv())
  m3 <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                           hp = list(K = 60L, sigthresh = 0.8))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m3, "sel_stablemate_model")
})

test_that("sel-stablemate is deployable through the canonical dispatch", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 17L)
  model <- fit_sel_stablemate(dat$X, dat$y, meta_train = dat$meta,
                              hp = list(K = 60L, sigthresh = 0.8))
  X_test <- dat$X[1:40, , drop = FALSE]

  # Self-contained one-row roster so singlesample_score_call resolves
  # sel-stablemate without depending on manifest integration order. The
  # dispatched score must equal the direct call.
  roster <- data.frame(
    method_id = "sel-stablemate", family = "N", estimand = "transfer",
    role = "discriminator", tier = "R1", dep_route = "StableMate-github",
    fit_fn = "fit_sel_stablemate", score_fn = "score_sel_stablemate",
    pkg_status = "new", notes = "", row_source = "transfer_fold_rows",
    lopbo_mechanism = "m4_decision", stringsAsFactors = FALSE)

  expect_equal(
    singlesample_score_call("sel-stablemate", model, X_test, roster = roster),
    score_sel_stablemate(model, X_test),
    tolerance = 1e-12)
})

test_that("sel-stablemate strict hp resolver rejects malformed hp", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 18L)

  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = "notalist"),
               "hp must be a list")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(foo = 5L)),
               "unknown hp")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(10L)),
               "must be named")
  expect_error(
    fit_sel_stablemate(dat$X, dat$y,
                       hp = stats::setNames(list(5L, 6L), c("K", "K"))),
    "duplicate hp")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(K = 1L)), "K")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(K = 2.5)), "K")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(max_features = 1L)),
               "max_features")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(min_selected = 0L)),
               "min_selected")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(sigthresh = 0)),
               "sigthresh")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(sigthresh = 1)),
               "sigthresh")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(cohort_col = 1)),
               "cohort_col")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(cohort_col = "")),
               "cohort_col")
  expect_error(fit_sel_stablemate(dat$X, dat$y, hp = list(seed = -1L)),
               "seed")
})

test_that("sel-stablemate input validation errors are explicit", {
  skip_if_no_stablemate()
  dat <- .make_sel_stablemate_data(seed = 19L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_sel_stablemate(X_unnamed, dat$y, dat$meta), "feature names")
  expect_error(fit_sel_stablemate(dat$X[, 1, drop = FALSE], dat$y, dat$meta),
               "at least two features")
  expect_error(fit_sel_stablemate(dat$X, dat$y[-1], dat$meta),
               "length\\(y_train\\)")
  expect_error(fit_sel_stablemate(dat$X, dat$y, dat$meta[-1, , drop = FALSE]),
               "one row per row")
  expect_error(fit_sel_stablemate(dat$X, rep(0, nrow(dat$X)), dat$meta),
               "at least one case")
  expect_error(fit_sel_stablemate(dat$X, rep(1, nrow(dat$X)), dat$meta),
               "at least one control")
  expect_error(score_sel_stablemate(list(), dat$X),
               "class sel_stablemate_model")
})
