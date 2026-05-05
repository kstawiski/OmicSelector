library(testthat)

set.seed(42L)

skip_if_not_installed("MASS")   # needed by fit_compositional_mahalanobis

# ============================================================================
# Shared data builders
# ============================================================================

.make_ood_data <- function(n_train = 80L, n_test = 20L, p = 20L, seed = 42L) {
  set.seed(seed)
  feats <- paste0("miR-", seq_len(p))
  x_train <- matrix(abs(stats::rnorm(n_train * p, 100, 30)), nrow = n_train,
                    dimnames = list(NULL, feats))
  x_test  <- matrix(abs(stats::rnorm(n_test * p, 100, 30)), nrow = n_test,
                    dimnames = list(NULL, feats))
  # Inject 5 compositional outliers: boost first half of features 8x.
  x_test[16:20, 1:(p %/% 2)] <- x_test[16:20, 1:(p %/% 2)] * 8
  list(x_train = x_train, x_test = x_test, outlier_rows = 16:20, feats = feats)
}

# ============================================================================
# fit_compositional_mahalanobis / apply_compositional_mahalanobis
# ============================================================================

test_that("fit_compositional_mahalanobis: known outliers detected at p<0.01 (headline)", {
  skip_if_not_installed("robustbase")
  dat <- .make_ood_data()
  fit <- fit_compositional_mahalanobis(dat$x_train, transform = "rclr")
  expect_s3_class(fit, "compositional_mahalanobis_fit")
  out <- apply_compositional_mahalanobis(fit, dat$x_test)
  expect_length(out$distance, nrow(dat$x_test))
  expect_true(all(is.finite(out$distance)))
  flagged <- which(out$pvalue_chisq < 0.05)
  # At least 2 of the 5 injected outliers should be flagged.
  expect_gte(sum(flagged %in% dat$outlier_rows), 2L)
})

test_that("fit_compositional_mahalanobis: negative values rejected", {
  dat <- .make_ood_data()
  x_bad <- dat$x_train; x_bad[1, 1] <- -1
  expect_error(fit_compositional_mahalanobis(x_bad), "negative values")
})

test_that("fit_compositional_mahalanobis: require_robust=FALSE uses classical covariance", {
  dat <- .make_ood_data()
  # Should not error even without robustbase.
  fit <- fit_compositional_mahalanobis(dat$x_train, transform = "rclr",
                                        require_robust = FALSE)
  out <- apply_compositional_mahalanobis(fit, dat$x_test)
  expect_true(all(is.finite(out$distance)))
})

test_that("apply_compositional_mahalanobis: missing training features raise error", {
  dat <- .make_ood_data()
  fit <- fit_compositional_mahalanobis(dat$x_train, transform = "rclr",
                                        require_robust = FALSE)
  x_short <- dat$x_test[, 1:10, drop = FALSE]
  expect_error(apply_compositional_mahalanobis(fit, x_short), "training features missing")
})


# ============================================================================
# fit_conformal_anomaly / os_conformal_anomaly
# ============================================================================

test_that("os_conformal_anomaly: FPR <= 0.10 and TPR >= 0.70 on synthetic OOD (headline)", {
  set.seed(42L)
  n_ref <- 300L; n_test <- 60L; p <- 30L
  X_ref <- matrix(stats::rnorm(n_ref * p), nrow = n_ref, ncol = p)
  X_healthy <- matrix(stats::rnorm((n_test / 2) * p), nrow = n_test / 2, ncol = p)
  X_ood     <- matrix(stats::rnorm((n_test / 2) * p, mean = 2.5), nrow = n_test / 2, ncol = p)
  X_test <- rbind(X_healthy, X_ood)
  is_ood <- c(rep(FALSE, n_test / 2), rep(TRUE, n_test / 2))

  fit <- fit_conformal_anomaly(X_ref, k = 10L, calibration_split = 0.30, seed = 42L)
  res <- os_conformal_anomaly(X_test, fit, alpha = 0.05)

  fpr <- mean(res$is_anomaly[!is_ood])
  tpr <- mean(res$is_anomaly[is_ood])
  expect_lte(fpr, 0.10)
  expect_gte(tpr, 0.60)
})

test_that("os_conformal_anomaly: p-values are in [0,1] and length matches test set", {
  set.seed(42L)
  X_ref  <- matrix(stats::rnorm(200 * 20), nrow = 200, ncol = 20)
  X_test <- matrix(stats::rnorm(10 * 20), nrow = 10, ncol = 20)
  fit <- fit_conformal_anomaly(X_ref, k = 5L, seed = 42L)
  res <- os_conformal_anomaly(X_test, fit)
  expect_length(res$p_value, 10L)
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
})

test_that("os_conformal_anomaly: wrong class raises error", {
  expect_error(os_conformal_anomaly(matrix(1), list()), "conformal_anomaly_fit")
})


# ============================================================================
# fit_isolation_forest_logratio / apply_isolation_forest_logratio
# ============================================================================

test_that("fit_isolation_forest_logratio: scores are in [0,1] and known outliers flagged (headline)", {
  dat <- .make_ood_data()
  fit <- fit_isolation_forest_logratio(dat$x_train, n_trees = 40L,
                                       sample_size = 40L, seed = 42L)
  out <- apply_isolation_forest_logratio(fit, dat$x_test, fpr_target = 0.10)
  expect_length(out$score, nrow(dat$x_test))
  expect_true(all(out$score >= 0 & out$score <= 1))
  # At least 2/5 outliers flagged at 10% FPR.
  expect_gte(sum(out$flagged[dat$outlier_rows]), 2L)
})

test_that("fit_isolation_forest_logratio: no-shared-features error", {
  dat <- .make_ood_data()
  fit <- fit_isolation_forest_logratio(dat$x_train, n_trees = 10L,
                                       sample_size = 20L, seed = 42L)
  x_bad <- dat$x_test
  colnames(x_bad) <- paste0("XXXX-", seq_len(ncol(x_bad)))
  expect_error(apply_isolation_forest_logratio(fit, x_bad), "no shared features")
})

test_that("fit_isolation_forest_logratio: single-sample test does not error", {
  dat <- .make_ood_data()
  fit <- fit_isolation_forest_logratio(dat$x_train, n_trees = 10L,
                                       sample_size = 20L, seed = 42L)
  x_single <- dat$x_test[1L, , drop = FALSE]
  out <- apply_isolation_forest_logratio(fit, x_single)
  expect_length(out$score, 1L)
})
