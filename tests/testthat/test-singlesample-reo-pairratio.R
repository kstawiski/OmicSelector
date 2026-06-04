library(testthat)

.make_reo_pairratio_data <- function(n = 120L, p_noise = 10L,
                                     seed = 31L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- c("A", "B", paste0("N", sprintf("%02d", seq_len(p_noise))))
  X <- matrix(stats::rgamma(n * length(features), shape = 20, rate = 1),
              nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))

  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  X[y == 0, "A"] <- 100 + stats::runif(n0, 0, 20)
  X[y == 0, "B"] <- 500 + stats::runif(n0, 0, 20)
  X[y == 1, "A"] <- 500 + stats::runif(n1, 0, 20)
  X[y == 1, "B"] <- 100 + stats::runif(n1, 0, 20)

  list(X = X, y = y)
}

.auc_or_wilcoxon_pairratio <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("reo-pairratio fit/score separates planted log-ratio signal", {
  skip_if_not_installed("glmnet")
  dat <- .make_reo_pairratio_data()
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_reo_pairratio(
    dat$X[train, ], dat$y[train],
    hp = list(m_features = 8L, nfolds = 5L, seed = 17L)
  )
  score <- score_reo_pairratio(model, dat$X[test, ])

  expect_s3_class(model, "reo_pairratio_model")
  expect_gt(nrow(model$selected_pairs), 0L)
  expect_true(any(model$selected_pairs$feature_a %in% c("A", "B") &
                    model$selected_pairs$feature_b %in% c("A", "B")))
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_pairratio(dat$y[test], score), 0.7)
})

test_that("reo-pairratio passes the canonical row-equivariance gate", {
  skip_if_not_installed("glmnet")
  dat <- .make_reo_pairratio_data(seed = 32L)
  model <- fit_reo_pairratio(
    dat$X[1:90, ], dat$y[1:90],
    hp = list(m_features = 8L, nfolds = 5L, seed = 17L)
  )
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_reo_pairratio(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("reo-pairratio fitting and scoring are deterministic", {
  skip_if_not_installed("glmnet")
  dat <- .make_reo_pairratio_data(seed = 33L)
  hp <- list(m_features = 8L, nfolds = 5L, seed = 19L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_reo_pairratio(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_reo_pairratio(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)

  expect_identical(model_1$selected_pairs, model_2$selected_pairs)
  expect_identical(model_1$coefficients, model_2$coefficients)
  expect_identical(model_1$intercept, model_2$intercept)
  expect_identical(model_1$lambda, model_2$lambda)

  score_1 <- score_reo_pairratio(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  score_2 <- score_reo_pairratio(model_1, dat$X[1:20, ])
  expect_identical(score_1, score_2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("reo-pairratio fitting leaves the global RNG untouched (no prior seed)", {
  skip_if_not_installed("glmnet")
  dat <- .make_reo_pairratio_data(seed = 37L)
  hp <- list(m_features = 8L, nfolds = 5L, seed = 23L)

  # No .Random.seed exists before fit: cv.glmnet can create one, so the fit
  # must restore the absent state (global-RNG hygiene, not just value-identity).
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  model <- fit_reo_pairratio(dat$X, dat$y, hp = hp)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "reo_pairratio_model")
})

test_that("reo-pairratio scores are finite and nearly scale-invariant", {
  skip_if_not_installed("glmnet")
  dat <- .make_reo_pairratio_data(seed = 34L)
  model <- fit_reo_pairratio(
    dat$X, dat$y,
    hp = list(m_features = 8L, nfolds = 5L, seed = 17L)
  )

  score <- score_reo_pairratio(model, dat$X[1:20, ])
  expect_type(score, "double")
  expect_length(score, 20L)
  expect_true(all(is.finite(score)))

  row_large <- dat$X[3, , drop = FALSE] * 1e6
  expect_equal(score_reo_pairratio(model, row_large * 7),
               score_reo_pairratio(model, row_large),
               tolerance = 1e-8)
})

test_that("reo-pairratio handles partial feature overlap", {
  skip_if_not_installed("glmnet")
  dat <- .make_reo_pairratio_data(seed = 35L)
  model <- fit_reo_pairratio(
    dat$X, dat$y,
    hp = list(m_features = 8L, nfolds = 5L, seed = 17L)
  )
  expect_gt(nrow(model$selected_pairs), 0L)
  p1 <- model$selected_pairs[1, ]

  X_partial <- dat$X[1:12, setdiff(colnames(dat$X), p1$feature_a),
                     drop = FALSE]
  score_partial <- score_reo_pairratio(model, X_partial)
  expect_length(score_partial, nrow(X_partial))
  expect_true(all(is.finite(score_partial)))
  expect_equal(score_reo_pairratio(model, X_partial[5, , drop = FALSE]),
               score_partial[5], tolerance = 1e-12)

  X_one <- dat$X[1:6, p1$feature_a, drop = FALSE]
  expect_equal(score_reo_pairratio(model, X_one),
               rep(model$intercept, nrow(X_one)), tolerance = 1e-12)

  X_no_shared <- dat$X[1:4, 1:2, drop = FALSE]
  colnames(X_no_shared) <- c("other-a", "other-b")
  expect_equal(score_reo_pairratio(model, X_no_shared),
               rep(model$intercept, nrow(X_no_shared)), tolerance = 1e-12)
})

test_that("reo-pairratio input validation errors are explicit", {
  skip_if_not_installed("glmnet")
  dat <- .make_reo_pairratio_data(seed = 36L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_reo_pairratio(X_unnamed, dat$y), "feature names")

  expect_error(fit_reo_pairratio(dat$X[, 1, drop = FALSE], dat$y),
               "at least two features")
  expect_error(fit_reo_pairratio(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_reo_pairratio(dat$X, rep(0, nrow(dat$X))),
               "at least one case")
  expect_error(fit_reo_pairratio(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")

  # Too few log-ratios for glmnet: 2 features -> 1 ratio -> clear method error
  # (not an opaque glmnet "x should be a matrix with 2 or more columns" crash).
  expect_error(
    fit_reo_pairratio(dat$X[, 1:2, drop = FALSE], dat$y),
    "log-ratios"
  )

  # Single-member minority class -> clear method error before glmnet.
  one_case <- c(which(dat$y == 0L), which(dat$y == 1L)[1])
  expect_error(
    fit_reo_pairratio(dat$X[one_case, , drop = FALSE], dat$y[one_case]),
    "at least 2 samples"
  )
})
