library(testthat)

.make_reo_singscore_data <- function(n = 96L, p = 80L, seed = 11L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  up <- features[1:8]
  down <- features[9:16]
  X[y == 1, up] <- X[y == 1, up] + 8
  X[y == 0, down] <- X[y == 0, down] + 8
  list(X = X, y = y, up = up, down = down)
}

.auc_or_wilcoxon_singscore <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("reo-singscore fit/score separates planted up/down signal", {
  dat <- .make_reo_singscore_data()
  train <- c(seq_len(36L), 49L:84L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_reo_singscore(dat$X[train, ], dat$y[train], hp = list(k = 8L))
  score <- score_reo_singscore(model, dat$X[test, ])

  expect_s3_class(model, "reo_singscore_model")
  expect_equal(model$k, 8L)
  expect_true(all(dat$up %in% model$up_features))
  expect_true(any(dat$down %in% model$down_features))
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_singscore(dat$y[test], score), 0.7)
})

test_that("reo-singscore passes the canonical row-equivariance gate", {
  dat <- .make_reo_singscore_data(seed = 12L)
  model <- fit_reo_singscore(dat$X[1:72, ], dat$y[1:72], hp = list(k = 8L))
  X_test <- dat$X[73:96, ]
  score_fun <- function(model, X, meta) score_reo_singscore(model, X, meta)

  expect_invisible(paper3_assert_row_equivariant(score_fun, model, X_test))
  expect_true(paper3_is_row_equivariant(score_fun, model, X_test))
})

test_that("reo-singscore equals os_singscore on a fully-present signature", {
  dat <- .make_reo_singscore_data(seed = 13L)
  model <- fit_reo_singscore(dat$X, dat$y, hp = list(k = 8L))
  X_full <- dat$X[1:20, model$feature_universe, drop = FALSE]

  score <- score_reo_singscore(model, X_full)
  primitive <- OmicSelector:::os_singscore(
    X_full,
    model$up_features,
    model$down_features
  )

  expect_equal(score, as.numeric(primitive), tolerance = 1e-10)
})

test_that("reo-singscore scores are deterministic finite numeric vectors", {
  dat <- .make_reo_singscore_data(seed = 14L)
  model <- fit_reo_singscore(dat$X, dat$y, hp = list(k = 8L))

  s1 <- score_reo_singscore(model, dat$X[1:20, ])
  s2 <- score_reo_singscore(model, dat$X[1:20, ])

  expect_identical(s1, s2)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("reo-singscore handles single-column / single-shared-feature inputs", {
  # (a) fit + score must not crash on a 1-column matrix, and the empty down
  # direction must contribute the documented neutral fractional rank.
  X1 <- matrix(c(10, 5, 2, 8), nrow = 4,
               dimnames = list(c("S1", "S2", "S3", "S4"), "miR-001"))
  y1 <- c(1, 0, 1, 0)
  m1 <- fit_reo_singscore(X1, y1, hp = list(k = 1L))
  s1 <- score_reo_singscore(m1, X1)
  expect_length(s1, 4L)
  expect_true(all(is.finite(s1)))
  expect_equal(s1, rep(0.5, 4L))

  # (b) reachable in single-sample deployment: a specimen sharing only ONE
  # feature with a larger frozen model universe. With a single present feature
  # the fractional rank is forced to 1 for every specimen, so the score is the
  # constant 0.5 (up = 1, empty down = neutral 0.5), value-invariant by
  # construction. This guards no-crash / finiteness / singleton == batch, and
  # that the names-restore fix returns 0.5 rather than the pre-fix silent 0
  # (a dropped-name empty intersect would have given 0.5 - 0.5 = 0).
  dat <- .make_reo_singscore_data(seed = 19L)
  model <- fit_reo_singscore(dat$X, dat$y, hp = list(k = 8L))
  shared <- model$up_features[1]
  X_one <- dat$X[1:5, shared, drop = FALSE]
  s_one <- score_reo_singscore(model, X_one)
  expect_length(s_one, 5L)
  expect_true(all(is.finite(s_one)))
  expect_equal(s_one, rep(0.5, 5L))
  expect_equal(score_reo_singscore(model, X_one[3, , drop = FALSE]), s_one[3],
               tolerance = 1e-12)
})

test_that("reo-singscore input validation errors are explicit", {
  dat <- .make_reo_singscore_data(seed = 15L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_reo_singscore(X_unnamed, dat$y), "feature names")

  expect_error(fit_reo_singscore(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_reo_singscore(dat$X, rep(0, nrow(dat$X))),
               "at least one case")
  expect_error(fit_reo_singscore(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")

  model <- fit_reo_singscore(dat$X, dat$y, hp = list(k = 8L))
  X_no_shared <- dat$X[1:4, , drop = FALSE]
  colnames(X_no_shared) <- paste0("other-", seq_len(ncol(X_no_shared)))
  expect_error(score_reo_singscore(model, X_no_shared), "no shared features")
})
