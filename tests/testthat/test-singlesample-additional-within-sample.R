library(testthat)

set.seed(42L)

# ============================================================================
# Shared data builder
# ============================================================================

.make_ws_data <- function(n = 50L, p = 20L, seed = 42L) {
  set.seed(seed)
  feats <- paste0("miR-", seq_len(p))
  x <- matrix(abs(stats::rnorm(n * p, 100, 30)), nrow = n,
              dimnames = list(NULL, feats))
  x
}

# ============================================================================
# fit_logistic_normal_eb / apply_logistic_normal_eb
# ============================================================================

test_that("fit_logistic_normal_eb: low-depth samples shrink harder toward prior (headline)", {
  set.seed(42L)
  x_train <- .make_ws_data(n = 50L, p = 20L)
  fit <- fit_logistic_normal_eb(x_train)
  expect_s3_class(fit, "logistic_normal_eb_fit")
  expect_length(fit$prior_var, 20L)
  expect_true(all(fit$prior_var > 0))

  set.seed(10L)
  x_lo <- matrix(abs(stats::rnorm(10 * 20, 5,  2)),  nrow = 10, dimnames = list(NULL, paste0("miR-", seq_len(20))))
  x_hi <- matrix(abs(stats::rnorm(10 * 20, 200, 40)), nrow = 10, dimnames = list(NULL, paste0("miR-", seq_len(20))))

  out_lo <- apply_logistic_normal_eb(fit, x_lo)
  out_hi <- apply_logistic_normal_eb(fit, x_hi)

  prior_mat <- matrix(fit$prior_mean, nrow = 10, ncol = 20, byrow = TRUE)
  d_lo <- mean((out_lo - prior_mat)^2)
  d_hi <- mean((out_hi - prior_mat)^2)
  expect_lt(d_lo, d_hi)   # low-count samples shrink harder (closer to prior)
})

test_that("fit_logistic_normal_eb: negative input rejected", {
  x_train <- .make_ws_data()
  x_bad <- x_train; x_bad[1, 1] <- -1
  expect_error(fit_logistic_normal_eb(x_bad), "negative values")
})

test_that("fit_logistic_normal_eb: missing features error by default; fill_with_prior pads", {
  x_train <- .make_ws_data()
  fit <- fit_logistic_normal_eb(x_train)

  x_short <- x_train[1:5, 1:15, drop = FALSE]
  expect_error(apply_logistic_normal_eb(fit, x_short), "missing")

  out_pad <- apply_logistic_normal_eb(fit, x_short, missing_features = "fill_with_prior")
  expect_equal(ncol(out_pad), 20L)
})

test_that("fit_logistic_normal_eb: matrix without colnames raises error", {
  x_anon <- matrix(abs(stats::rnorm(50 * 20, 100, 30)), nrow = 50, ncol = 20)
  expect_error(fit_logistic_normal_eb(x_anon), "feature names")
})


# ============================================================================
# fit_frozen_quantile / apply_frozen_quantile
# ============================================================================

test_that("fit_frozen_quantile: apply aligns test distribution to training scale (headline)", {
  set.seed(42L)
  feats <- paste0("miR-", seq_len(20))
  x_train <- matrix(abs(stats::rnorm(50 * 20, 100, 30)), nrow = 50,
                    dimnames = list(NULL, feats))
  x_test  <- matrix(abs(stats::rnorm(10 * 20, 150, 50)), nrow = 10,
                    dimnames = list(NULL, feats))
  fit <- fit_frozen_quantile(x_train, n_quantiles = 100L)
  out <- apply_frozen_quantile(fit, x_test)

  expect_equal(dim(out), c(10L, 20L))
  expect_true(all(is.finite(out)))
  # Post-mapping median should be closer to training median than raw test median.
  expect_lt(abs(stats::median(out) - stats::median(x_train)),
            abs(stats::median(x_test) - stats::median(x_train)))
})

test_that("fit_frozen_quantile: negative values rejected at fit and apply", {
  x_train <- .make_ws_data()
  x_bad <- x_train; x_bad[1, 1] <- -1
  expect_error(fit_frozen_quantile(x_bad), "negative values")

  fit <- fit_frozen_quantile(x_train, n_quantiles = 50L)
  x_test_bad <- .make_ws_data(n = 5L, seed = 99L); x_test_bad[1, 1] <- -1
  expect_error(apply_frozen_quantile(fit, x_test_bad), "negative values")
})

test_that("fit_frozen_quantile: missing features error by default; skip returns subset", {
  feats <- paste0("miR-", seq_len(20))
  x_train <- matrix(abs(stats::rnorm(50 * 20, 100, 30)), nrow = 50,
                    dimnames = list(NULL, feats))
  fit <- fit_frozen_quantile(x_train, n_quantiles = 50L)

  # Only 10 of 20 features present.
  x_test_short <- matrix(abs(stats::rnorm(5 * 10, 100, 30)), nrow = 5,
                          dimnames = list(NULL, feats[1:10]))
  expect_error(apply_frozen_quantile(fit, x_test_short, missing_features = "error"),
               "training features missing")

  out_skip <- apply_frozen_quantile(fit, x_test_short, missing_features = "skip")
  expect_equal(ncol(out_skip), 10L)   # only shared features returned
})

test_that("fit_frozen_quantile: matrix without colnames raises error", {
  x_anon <- matrix(abs(stats::rnorm(50 * 20)), nrow = 50, ncol = 20)
  expect_error(fit_frozen_quantile(x_anon), "feature names")
})
