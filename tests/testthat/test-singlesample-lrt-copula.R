library(testthat)

# Synthetic generator with a CASE-ONLY feature-correlation block. A shared latent
# factor loads onto the first `k` features for case samples only, so those
# features are correlated across cases but independent across controls. The
# loadings ALTERNATE in sign (sum to zero), i.e. the factor is a contrast among
# the block features. A contrast is orthogonal to the within-sample rCLR centring
# direction, so unlike a uniform shift it is fully preserved by the rCLR and
# produces a strong class-specific dependence (first-half block features move up
# while second-half move down together, within cases only). Because the Gaussian-
# copula LRT matches out marginals via the PIT and scores a specimen purely by how
# well its feature DEPENDENCE matches each class, this is exactly the signal the
# method is built to detect. Abundances are exp() of the log-scale construction so
# the input is strictly positive (compositional).
.make_lrt_copula_data <- function(n = 160L, p = 30L, k = 16L, load = 2.0,
                                  sd = 0.5, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case <- which(y == 1L)
  f <- stats::rnorm(length(case))                 # shared latent factor (cases)
  loadings <- rep(c(load, -load), length.out = k) # mean-zero contrast loadings
  L[case, seq_len(k)] <- L[case, seq_len(k)] + outer(f, loadings)
  X <- exp(L)
  list(X = X, y = y, block = features[seq_len(k)])
}

.auc_or_wilcoxon_cop <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("lrt-copula fit/score roundtrip has the right shape and types", {
  dat <- .make_lrt_copula_data()
  model <- fit_lrt_copula(dat$X, dat$y,
                          hp = list(max_anchors_per_class = 40L, seed = 7L))
  expect_s3_class(model, "lrt_copula_model")
  expect_equal(nrow(model$case_anchors), 40L)
  expect_equal(nrow(model$control_anchors), 40L)
  expect_identical(model$feature_universe, colnames(dat$X))
  # Frozen full-universe precision matrices are symmetric and the right size.
  p <- ncol(dat$X)
  expect_equal(dim(model$repr$Rinv_case), c(p, p))
  expect_equal(dim(model$repr$Rinv_control), c(p, p))

  score <- score_lrt_copula(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("lrt-copula separates the planted class-specific dependence block", {
  # Larger cohort so the held-out AUC is a stable estimate of the (genuine, copula-
  # only) dependence signal rather than a small-sample fluctuation.
  dat <- .make_lrt_copula_data(n = 300L, seed = 202L)
  ctrl <- which(dat$y == 0L)
  case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_lrt_copula(dat$X[train, ], dat$y[train],
                          hp = list(max_anchors_per_class = 112L, seed = 11L))
  score <- score_lrt_copula(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_cop(dat$y[test], score), 0.7)
})

test_that("lrt-copula passes the canonical row-equivariance gate", {
  dat <- .make_lrt_copula_data(seed = 43L)
  model <- fit_lrt_copula(dat$X[1:120, ], dat$y[1:120],
                          hp = list(max_anchors_per_class = 30L, seed = 13L))
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_lrt_copula(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("lrt-copula scoring is exactly invariant to per-sample scaling", {
  dat <- .make_lrt_copula_data(seed = 44L)
  model <- fit_lrt_copula(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]
  expect_equal(score_lrt_copula(model, row * 7),
               score_lrt_copula(model, row), tolerance = 1e-8)
  expect_equal(score_lrt_copula(model, row * 1e6),
               score_lrt_copula(model, row), tolerance = 1e-8)
  # Row-equivariant batch: scaling some rows leaves the others' scores fixed.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e6
  expect_equal(score_lrt_copula(model, scaled),
               score_lrt_copula(model, batch), tolerance = 1e-8)
  # A row with mostly zeros (rCLR over its own nonzero parts) still scores finite.
  X_sparse <- dat$X[1:3, , drop = FALSE]
  X_sparse[1, 3:ncol(X_sparse)] <- 0
  expect_true(all(is.finite(score_lrt_copula(model, X_sparse))))
})

test_that("lrt-copula handles partial feature overlap (consistent + single-row)", {
  dat <- .make_lrt_copula_data(seed = 45L)
  model <- fit_lrt_copula(dat$X, dat$y)
  # Subset and REORDER the columns: present = intersect(universe, colnames) keeps
  # feature-universe order, so the score is invariant to the X column order.
  keep <- rev(c(dat$block[1:5], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_lrt_copula(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_lrt_copula(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_lrt_copula(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("lrt-copula matches the canonical dispatch", {
  dat <- .make_lrt_copula_data(seed = 46L)
  model <- fit_lrt_copula(dat$X[1:120, ], dat$y[1:120],
                          hp = list(max_anchors_per_class = 30L, seed = 17L))
  X_test <- dat$X[121:160, ]
  # In-package, singlesample_score_call() resolves the manifest score_fn
  # (score_lrt_copula) directly from the OmicSelector namespace -- no adapter
  # registration needed; this exercises the real roster->package dispatch path.
  expect_equal(
    singlesample_score_call("lrt-copula", model, X_test),
    score_lrt_copula(model, X_test),
    tolerance = 1e-12
  )
})

test_that("lrt-copula returns the neutral 0 below the feature-overlap floor", {
  dat <- .make_lrt_copula_data(seed = 47L)
  model <- fit_lrt_copula(dat$X, dat$y, hp = list(min_features = 3L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, colnames(dat$X)[1:2], drop = FALSE]
  expect_equal(score_lrt_copula(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_lrt_copula(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
})

test_that("lrt-copula fitting/scoring is deterministic and RNG-safe", {
  dat <- .make_lrt_copula_data(seed = 48L)
  hp <- list(max_anchors_per_class = 25L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_lrt_copula(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_lrt_copula(dat$X, dat$y, hp = hp)
  expect_identical(model_1, model_2)

  s1 <- score_lrt_copula(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_lrt_copula(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("lrt-copula fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_lrt_copula_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  # max_anchors_per_class > class size -> no subsampling -> RNG never touched.
  model <- fit_lrt_copula(dat$X, dat$y,
                          hp = list(max_anchors_per_class = 200L, seed = 29L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "lrt_copula_model")
})

test_that("lrt-copula leaves no RNG state when subsampling fits in a fresh session", {
  # Companion to the no-subsampling case above, exercising the OTHER branch of the
  # anchor-subsampling RNG save/restore: here a class EXCEEDS the anchor cap, so
  # set.seed() IS called inside the fit. With no prior global seed, the on.exit
  # restore must rm() the .Random.seed that set.seed() created -- not leave it
  # behind. (80 anchors/class > cap 20 -> subsampling path taken.)
  dat <- .make_lrt_copula_data(seed = 53L)        # this set.seed() seeds the RNG
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())       # wipe to a genuine no-seed state
  }
  model <- fit_lrt_copula(dat$X, dat$y,
                          hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_lte(nrow(model$case_anchors), 20L)        # subsampling actually happened
  expect_lte(nrow(model$control_anchors), 20L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # And subsampling is deterministic across fits (frozen seed).
  model2 <- fit_lrt_copula(dat$X, dat$y,
                           hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_identical(rownames(model$case_anchors), rownames(model2$case_anchors))
})

test_that("lrt-copula handles a constant/degenerate training feature without NA", {
  # A feature that is zero across ALL anchors of a class has a constant (all-zero)
  # rCLR column. Its empirical marginal collapses to a SINGLE knot (no approx call,
  # PIT -> the knot plotting position, q ~ 0) and its Gaussian-score column has zero
  # variance, so cor() yields NA for that row/column. The marginal single-knot guard
  # and the correlation NA->identity guard must keep the fit and every score finite.
  dat <- .make_lrt_copula_data(seed = 54L)
  X <- dat$X
  X[, ncol(X)] <- 0                       # feature constant (zero) in BOTH classes
  X[dat$y == 1L, ncol(X) - 1L] <- 0       # feature constant (zero) in cases only
  model <- fit_lrt_copula(X, dat$y, hp = list(max_anchors_per_class = 40L))
  expect_true(all(is.finite(model$repr$Rinv_case)))
  expect_true(all(is.finite(model$repr$Rinv_control)))
  expect_true(is.finite(model$repr$logdet_case))
  expect_true(is.finite(model$repr$logdet_control))
  score <- score_lrt_copula(model, X)
  expect_length(score, nrow(X))
  expect_true(all(is.finite(score)))
  # Scoring a specimen that is itself zero on the degenerate features is still finite.
  expect_true(all(is.finite(score_lrt_copula(model, X[1:4, , drop = FALSE]))))
})

test_that("lrt-copula hyperparameter validation is strict", {
  dat <- .make_lrt_copula_data(seed = 50L)
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(shrink = 1)),
               "hp\\$shrink")
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(pit_clamp = 0.7)),
               "hp\\$pit_clamp")
  # Both open-interval endpoints are rejected: pit_clamp = 0 would make
  # qnorm(0) = -Inf and pit_clamp = 0.5 collapses the clamp window.
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(pit_clamp = 0)),
               "hp\\$pit_clamp")
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(pit_clamp = 0.5)),
               "hp\\$pit_clamp")
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
               "max_anchors_per_class")
  # Out-of-integer-range must hit the explicit error, not an as.integer overflow.
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
               "max_anchors_per_class")
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_lrt_copula(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
})

test_that("lrt-copula copula log-density ratio is finite with the correct sign", {
  # Tiny 3-feature example. The latent factor loads on a SUBSET (f1, f2) of the
  # three features in cases: a uniform shift across all features would be removed
  # by the within-sample rCLR centring, so the case-specific dependence survives
  # only as a positive z1-z2 correlation relative to controls. A case-consistent
  # specimen (f1, f2 jointly elevated above f3) must score above a case-
  # inconsistent specimen (f1 up while f2 down).
  set.seed(123)
  n <- 80L
  y <- rep(c(0, 1), each = n / 2L)
  feats <- c("f1", "f2", "f3")
  L <- matrix(stats::rnorm(n * 3L, 3, 0.5), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), feats))
  case <- which(y == 1L)
  fac <- stats::rnorm(length(case))
  L[case, 1:2] <- L[case, 1:2] + outer(fac, c(1.7, 1.7))  # block on f1, f2 only
  X <- exp(L)
  model <- fit_lrt_copula(X, y)

  x_caselike <- matrix(exp(c(5.0, 5.0, 3.0)), nrow = 1,
                       dimnames = list(NULL, feats))   # z1,z2 jointly high
  x_ctrllike <- matrix(exp(c(5.0, 3.0, 5.0)), nrow = 1,
                       dimnames = list(NULL, feats))   # z1 high, z2 low
  s_case <- score_lrt_copula(model, x_caselike)
  s_ctrl <- score_lrt_copula(model, x_ctrllike)
  expect_true(is.finite(s_case) && is.finite(s_ctrl))
  expect_length(s_case, 1L)
  expect_gt(s_case, s_ctrl)

  # On the training batch the case group scores higher on average.
  s_all <- score_lrt_copula(model, X)
  expect_true(all(is.finite(s_all)))
  expect_gt(mean(s_all[y == 1]), mean(s_all[y == 0]))
})
