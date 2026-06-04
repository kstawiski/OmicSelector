library(testthat)

.make_kme_witness_data <- function(n = 120L, p = 40L, seed = 61L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 30, rate = 1), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case_high <- features[1:8]
  control_high <- features[9:16]
  n0 <- sum(y == 0L)
  n1 <- sum(y == 1L)

  X[y == 1L, case_high] <- 140 +
    matrix(stats::runif(n1 * length(case_high), 0, 3), nrow = n1)
  X[y == 1L, control_high] <- 8 +
    matrix(stats::runif(n1 * length(control_high), 0, 3), nrow = n1)
  X[y == 0L, case_high] <- 8 +
    matrix(stats::runif(n0 * length(case_high), 0, 3), nrow = n0)
  X[y == 0L, control_high] <- 140 +
    matrix(stats::runif(n0 * length(control_high), 0, 3), nrow = n0)

  list(X = X, y = y, case_high = case_high, control_high = control_high)
}

.auc_or_wilcoxon_kme <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("kme-witness fit/score separates planted CLR signal case-high", {
  dat <- .make_kme_witness_data()
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_kme_witness(
    dat$X[train, ], dat$y[train],
    hp = list(max_anchors_per_class = 20L, seed = 17L)
  )
  score <- score_kme_witness(model, dat$X[test, ])

  expect_s3_class(model, "kme_witness_model")
  expect_equal(nrow(model$case_anchors), 20L)
  expect_equal(nrow(model$control_anchors), 20L)
  expect_true(is.finite(model$gamma))
  expect_gt(model$gamma, 0)
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_kme(dat$y[test], score), 0.7)
})

test_that("kme-witness passes the canonical row-equivariance gate", {
  dat <- .make_kme_witness_data(seed = 62L)
  model <- fit_kme_witness(
    dat$X[1:90, ], dat$y[1:90],
    hp = list(max_anchors_per_class = 18L, seed = 19L)
  )
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_kme_witness(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("kme-witness is deployable through the canonical dispatch", {
  dat <- .make_kme_witness_data(seed = 63L)
  model <- fit_kme_witness(
    dat$X[1:90, ], dat$y[1:90],
    hp = list(max_anchors_per_class = 18L, seed = 19L)
  )
  X_test <- dat$X[91:120, ]

  expect_equal(
    singlesample_score_call("kme-witness", model, X_test),
    score_kme_witness(model, X_test),
    tolerance = 1e-12
  )
})

test_that("kme-witness fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_kme_witness_data(seed = 64L)
  hp <- list(max_anchors_per_class = 10L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_kme_witness(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_kme_witness(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  expect_identical(model_1, model_2)

  s1 <- score_kme_witness(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_kme_witness(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("kme-witness fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_kme_witness_data(seed = 65L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  model <- fit_kme_witness(
    dat$X, dat$y,
    hp = list(max_anchors_per_class = 10L, seed = 29L)
  )

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "kme_witness_model")
})

test_that("kme-witness scoring is exactly scale-invariant at the default pc", {
  # Closing to proportions before the pseudocount CLR makes the score exactly
  # invariant to per-sample scaling for ANY pc -- at the DEFAULT pc and ordinary
  # magnitudes, not only in a contrived tiny-pc / pre-scaled regime.
  dat <- .make_kme_witness_data(seed = 66L)
  model <- fit_kme_witness(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_kme_witness(model, row * 7),
               score_kme_witness(model, row), tolerance = 1e-8)
  expect_equal(score_kme_witness(model, row * 1e6),
               score_kme_witness(model, row), tolerance = 1e-8)
  # The batch is scored row-equivariantly: scaling some rows leaves others fixed.
  batch <- dat$X[1:5, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e6
  expect_equal(score_kme_witness(model, scaled),
               score_kme_witness(model, batch), tolerance = 1e-8)

  # An all-zero specimen (closure -> uniform composition) scores finite.
  X_zero <- dat$X[1:3, , drop = FALSE]
  X_zero[1, ] <- 0
  expect_true(all(is.finite(score_kme_witness(model, X_zero))))
})

test_that("kme-witness handles partial feature overlap as finite neutral-aware scores", {
  dat <- .make_kme_witness_data(seed = 67L)
  model <- fit_kme_witness(dat$X, dat$y)

  keep <- c(dat$case_high[1:4], dat$control_high[1:4],
            colnames(dat$X)[25:30])
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_kme_witness(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  expect_equal(score_kme_witness(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)

  X_one <- dat$X[1:6, dat$case_high[1], drop = FALSE]
  expect_equal(score_kme_witness(model, X_one),
               rep(0, nrow(X_one)), tolerance = 1e-12)

  X_no_shared <- matrix(1, nrow = 4L, ncol = 1L,
                        dimnames = list(paste0("N", seq_len(4L)), "other"))
  expect_equal(score_kme_witness(model, X_no_shared),
               rep(0, nrow(X_no_shared)), tolerance = 1e-12)
})

test_that("kme-witness input validation errors are explicit", {
  dat <- .make_kme_witness_data(seed = 68L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_kme_witness(X_unnamed, dat$y), "feature names")

  expect_error(fit_kme_witness(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_kme_witness(dat$X, rep(0, nrow(dat$X))),
               "at least one case")
  expect_error(fit_kme_witness(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")

  expect_error(fit_kme_witness(dat$X, dat$y, hp = list(pc = 0)),
               "hp\\$pc")
  expect_error(fit_kme_witness(dat$X, dat$y, hp = list(gamma = -1)),
               "hp\\$gamma")
  expect_error(
    fit_kme_witness(dat$X, dat$y,
                    hp = list(max_anchors_per_class = 1.5)),
    "max_anchors_per_class"
  )
  # Out-of-integer-range values must hit the explicit validation error, not a
  # generic as.integer() overflow ("missing value where TRUE/FALSE needed").
  expect_error(
    fit_kme_witness(dat$X, dat$y,
                    hp = list(max_anchors_per_class = 3e9)),
    "max_anchors_per_class"
  )
  expect_error(
    fit_kme_witness(dat$X, dat$y, hp = list(min_features = 3e9)),
    "min_features"
  )
})
