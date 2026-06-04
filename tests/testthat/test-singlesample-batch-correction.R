library(testthat)

set.seed(42L)

skip_if_not_installed("MASS")

# ============================================================================
# Shared synthetic data builder
# ============================================================================

.make_batch_data <- function(n_train = 150L, n_test = 40L, p = 50L, seed = 42L) {
  set.seed(seed)
  batch <- stats::rnorm(n_train + n_test, sd = 1.5)
  feat_loadings <- stats::rnorm(p, sd = 0.5)
  X_clean <- matrix(stats::rnorm((n_train + n_test) * p, sd = 0.5),
                    nrow = n_train + n_test, ncol = p)
  X <- X_clean + batch %o% feat_loadings +
    matrix(stats::rnorm((n_train + n_test) * p, sd = 0.1),
            nrow = n_train + n_test)
  ctrl_names <- c(
    "hsa-miR-93-5p", "hsa-miR-26a-5p", "hsa-miR-191-5p", "hsa-miR-92a-3p",
    "hsa-miR-30c-5p", "hsa-let-7g-5p", "hsa-miR-103a-3p", "hsa-miR-23a-3p",
    "hsa-miR-21-5p"
  )
  colnames(X) <- c(ctrl_names, paste0("mir-", sprintf("%03d", 10:p)))
  list(
    X_train = X[seq_len(n_train), , drop = FALSE],
    X_test  = X[(n_train + 1):(n_train + n_test), , drop = FALSE],
    batch_test = batch[(n_train + 1):(n_train + n_test)]
  )
}

# ============================================================================
# fit_frozen_ruv / apply_frozen_ruv
# ============================================================================

test_that("fit_frozen_ruv: batch correlation reduced after correction (headline behaviour)", {
  dat <- .make_batch_data()
  fit <- fit_frozen_ruv(dat$X_train, k = 2:4)
  expect_s3_class(fit, "frozen_ruv_fit")
  expect_true(fit$k_selected %in% 2:4)

  X_corr <- apply_frozen_ruv(dat$X_test, fit)
  cor_pre  <- mean(abs(stats::cor(dat$X_test, dat$batch_test,
                                   use = "complete.obs")))
  cor_post <- mean(abs(stats::cor(X_corr, dat$batch_test,
                                   use = "complete.obs")))
  expect_lt(cor_post, cor_pre)
})

test_that("fit_frozen_ruv: falls back to low-variance proxy when controls not found", {
  dat <- .make_batch_data()
  # Strip control feature names so none are matched.
  X_anon <- dat$X_train
  colnames(X_anon) <- paste0("anon-", seq_len(ncol(X_anon)))
  expect_warning(
    fit <- fit_frozen_ruv(X_anon, k = 2:3),
    "control features found"
  )
  expect_s3_class(fit, "frozen_ruv_fit")
})

test_that("apply_frozen_ruv: returns same dimensions as input", {
  dat <- .make_batch_data()
  fit <- fit_frozen_ruv(dat$X_train, k = 2:3)
  out <- apply_frozen_ruv(dat$X_test, fit)
  expect_equal(dim(out), dim(dat$X_test))
})

test_that("apply_frozen_ruv: vector input (single sample) is accepted", {
  dat <- .make_batch_data()
  fit <- fit_frozen_ruv(dat$X_train, k = 2L, select = "fixed")
  x_vec <- dat$X_test[1, ]
  out <- apply_frozen_ruv(x_vec, fit)
  expect_length(out, ncol(dat$X_test))
})


# ============================================================================
# fit_robust_pca_residual / apply_robust_pca_residual
# ============================================================================

.make_pca_data <- function(n = 80L, p = 40L, seed = 42L) {
  set.seed(seed)
  feats <- paste0("miR-", seq_len(p))
  bio <- matrix(stats::rnorm(n * p, sd = 0.3), nrow = n)
  batch <- matrix(rep(stats::rnorm(p, sd = 1.5), n), nrow = n, byrow = TRUE)
  x <- exp(bio + batch) * 100
  colnames(x) <- feats
  x
}

test_that("fit_robust_pca_residual: returns fit object with correct dimensions (headline)", {
  x_train <- .make_pca_data()
  fit <- fit_robust_pca_residual(x_train, n_components = 3L)
  expect_s3_class(fit, "robust_pca_residual_fit")
  expect_equal(ncol(fit$v), 3L)
  expect_equal(nrow(fit$v), ncol(x_train))
})

test_that("fit_robust_pca_residual: rejects negative values", {
  x_train <- .make_pca_data()
  x_bad <- x_train; x_bad[1, 1] <- -1
  expect_error(fit_robust_pca_residual(x_bad), "negative values")
})

test_that("fit_robust_pca_residual: residual + reconstruction = original z (round-trip)", {
  x_train <- .make_pca_data()
  x_test  <- .make_pca_data(n = 10L, seed = 99L)
  fit <- fit_robust_pca_residual(x_train, n_components = 3L)
  res <- apply_robust_pca_residual(fit, x_test, return_reconstructed = TRUE)
  expect_equal(dim(res$scores),       c(10L, 3L))
  expect_equal(dim(res$residual),     c(10L, ncol(x_test)))
  expect_equal(dim(res$reconstructed), c(10L, ncol(x_test)))
  # residual + reconstructed must sum to z_test (within floating-point error)
  z_sum <- res$residual + res$reconstructed
  z_ref <- log(x_test + 0.5)
  z_ref <- sweep(z_ref, 2L, apply(log(x_train + 0.5), 2L, stats::median), "-")
  z_ref <- sweep(z_ref, 2L, apply(log(x_train + 0.5), 2L, stats::mad), "/")
  expect_lt(max(abs(z_sum - z_ref)), 1e-10)
})

test_that("apply_robust_pca_residual: missing training features raise error", {
  x_train <- .make_pca_data()
  fit <- fit_robust_pca_residual(x_train, n_components = 2L)
  x_short <- x_train[1:5, 1:20, drop = FALSE]
  expect_error(apply_robust_pca_residual(fit, x_short), "training features missing")
})
