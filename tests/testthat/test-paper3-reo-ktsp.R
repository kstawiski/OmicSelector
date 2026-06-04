library(testthat)

.make_reo_ktsp_data <- function(n = 96L, p_noise = 12L, seed = 21L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  pair_a <- paste0("A", seq_len(3L))
  pair_b <- paste0("B", seq_len(3L))
  noise <- paste0("N", sprintf("%02d", seq_len(p_noise)))
  features <- c(pair_a, pair_b, noise)
  X <- matrix(stats::rgamma(n * length(features), shape = 2, rate = 1),
              nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))

  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  X[y == 0, pair_a] <- 12 + matrix(stats::runif(n0 * length(pair_a), 0, 0.1),
                                   nrow = n0)
  X[y == 0, pair_b] <- 1 + matrix(stats::runif(n0 * length(pair_b), 0, 0.1),
                                  nrow = n0)
  X[y == 1, pair_a] <- 1 + matrix(stats::runif(n1 * length(pair_a), 0, 0.1),
                                  nrow = n1)
  X[y == 1, pair_b] <- 12 + matrix(stats::runif(n1 * length(pair_b), 0, 0.1),
                                   nrow = n1)

  list(X = X, y = y, pair_a = pair_a, pair_b = pair_b)
}

.auc_or_wilcoxon_ktsp <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("reo-ktsp fit/score separates planted pair-order signal", {
  dat <- .make_reo_ktsp_data()
  train <- c(seq_len(36L), 49L:84L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_reo_ktsp(dat$X[train, ], dat$y[train], hp = list(k = 5L))
  score <- score_reo_ktsp(model, dat$X[test, ])

  expect_s3_class(model, "reo_ktsp_model")
  expect_s3_class(model$ktsp, "os_ktsp_model")
  expect_equal(model$hp$k, 5L)
  expect_equal(model$k, 5L)
  expect_true(any(model$ktsp$pairs$feature_a == "A1" &
                    model$ktsp$pairs$feature_b == "B1"))
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_ktsp(dat$y[test], score), 0.7)
})

test_that("reo-ktsp passes the canonical row-equivariance gate", {
  dat <- .make_reo_ktsp_data(seed = 22L)
  model <- fit_reo_ktsp(dat$X[1:72, ], dat$y[1:72], hp = list(k = 5L))
  X_test <- dat$X[73:96, ]
  score_fun <- function(model, X, meta) score_reo_ktsp(model, X, meta)

  expect_invisible(paper3_assert_row_equivariant(score_fun, model, X_test))
  expect_true(paper3_is_row_equivariant(score_fun, model, X_test))
})

test_that("reo-ktsp equals predict.os_ktsp_model on fully-present features", {
  dat <- .make_reo_ktsp_data(seed = 23L)
  model <- fit_reo_ktsp(dat$X, dat$y, hp = list(k = 5L))
  X_full <- dat$X[1:20, model$feature_universe, drop = FALSE]

  score <- score_reo_ktsp(model, X_full)
  primitive <- predict(model$ktsp, X_full)

  expect_equal(score, as.numeric(primitive), tolerance = 1e-12)

  # single-row equivalence: predict.os_ktsp_model is now single-sample safe, so
  # a 1-row specimen scores identically through the wrapper and the primitive.
  one <- score_reo_ktsp(model, X_full[1, , drop = FALSE])
  expect_length(one, 1L)
  expect_equal(one, as.numeric(predict(model$ktsp, X_full[1, , drop = FALSE])),
               tolerance = 1e-12)
})

test_that("reo-ktsp scores correctly with exactly one usable pair", {
  dat <- .make_reo_ktsp_data(seed = 27L)
  model <- fit_reo_ktsp(dat$X, dat$y, hp = list(k = 5L))
  p1 <- model$ktsp$pairs[1, ]
  # Restrict X to a single retained pair's two features -> exactly one usable
  # pair, exercising the do.call(cbind, ...) single-column path directly.
  keep <- unique(c(p1$feature_a, p1$feature_b))
  X_1pair <- dat$X[1:6, keep, drop = FALSE]
  s <- score_reo_ktsp(model, X_1pair)
  expect_length(s, 6L)
  expect_true(all(is.finite(s)))
  # vote fraction over one pair == that pair's raw within-row vote
  expect_equal(s, as.numeric(X_1pair[, p1$feature_a] < X_1pair[, p1$feature_b]))
  # and the single-row case stays scalar and consistent
  expect_equal(score_reo_ktsp(model, X_1pair[3, , drop = FALSE]), s[3])
})

test_that("reo-ktsp scores are deterministic finite numeric vectors", {
  dat <- .make_reo_ktsp_data(seed = 24L)
  model <- fit_reo_ktsp(dat$X, dat$y, hp = list(k = 5L))

  s1 <- score_reo_ktsp(model, dat$X[1:20, ])
  s2 <- score_reo_ktsp(model, dat$X[1:20, ])

  expect_identical(s1, s2)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("reo-ktsp handles partial overlap and neutral no-pair cases", {
  dat <- .make_reo_ktsp_data(seed = 25L)
  model <- fit_reo_ktsp(dat$X, dat$y, hp = list(k = 5L))
  dropped <- model$ktsp$pairs$feature_a[1]
  X_partial <- dat$X[1:12, setdiff(colnames(dat$X), dropped), drop = FALSE]

  score_partial <- score_reo_ktsp(model, X_partial)
  expect_length(score_partial, nrow(X_partial))
  expect_true(all(is.finite(score_partial)))
  usable <- model$ktsp$pairs$feature_a %in% colnames(X_partial) &
    model$ktsp$pairs$feature_b %in% colnames(X_partial)
  expect_true(any(usable))

  row_id <- 5L
  expect_equal(score_reo_ktsp(model, X_partial[row_id, , drop = FALSE]),
               score_partial[row_id], tolerance = 1e-12)

  X_one <- dat$X[1:6, model$ktsp$pairs$feature_a[1], drop = FALSE]
  expect_equal(score_reo_ktsp(model, X_one), rep(0.5, nrow(X_one)))

  X_no_shared <- dat$X[1:4, 1:2, drop = FALSE]
  colnames(X_no_shared) <- c("other-a", "other-b")
  expect_equal(score_reo_ktsp(model, X_no_shared),
               rep(0.5, nrow(X_no_shared)))
})

test_that("reo-ktsp input validation errors are explicit", {
  dat <- .make_reo_ktsp_data(seed = 26L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_reo_ktsp(X_unnamed, dat$y), "feature names")

  expect_error(fit_reo_ktsp(dat$X[, 1, drop = FALSE], dat$y),
               "at least two features")
  expect_error(fit_reo_ktsp(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_reo_ktsp(dat$X, rep(0, nrow(dat$X))),
               "at least one case")
  expect_error(fit_reo_ktsp(dat$X, dat$y, hp = list(k = 1.5)),
               "positive integer")
  expect_error(fit_reo_ktsp(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")
})
