library(testthat)

set.seed(42L)

skip_if_not_installed("MASS")

# ============================================================================
# Shared synthetic dataset
# ============================================================================

.make_hemo_data <- function(n = 100L, p = 50L, seed = 42L) {
  set.seed(seed)
  X_clean <- matrix(stats::rlnorm(n * p, meanlog = 5, sdlog = 0.8),
                    nrow = n, ncol = p)
  feat_names <- c(paste0("mir-", sprintf("%03d", seq_len(p - 5L))),
                  "hsa-miR-451a", "hsa-miR-16-5p",
                  "hsa-miR-486-5p", "hsa-miR-144-3p", "hsa-miR-223-3p")
  colnames(X_clean) <- feat_names
  hemo <- stats::runif(n, 0, 1)
  shift_idx <- which(colnames(X_clean) %in%
                       c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                         "hsa-miR-144-3p", "hsa-miR-223-3p"))
  X_obs <- X_clean
  X_obs[, shift_idx] <- X_obs[, shift_idx] * exp(2 * hemo)
  list(X_clean = X_clean, X_obs = X_obs, hemo = hemo,
       shift_idx = shift_idx, feat_names = feat_names)
}

# ============================================================================
# fit_hemolysis_rr
# ============================================================================

test_that("fit_hemolysis_rr: fitted object has correct structure (headline behaviour)", {
  dat <- .make_hemo_data()
  fit <- fit_hemolysis_rr(dat$X_obs, dat$hemo)
  expect_s3_class(fit, "hemolysis_rr_fit")
  expect_equal(nrow(fit$betas), ncol(dat$X_obs))
  expect_true("beta_hemo" %in% colnames(fit$betas))
  expect_false(fit$has_platelet)
  expect_equal(fit$n_train, nrow(dat$X_obs))
})

test_that("fit_hemolysis_rr: hemo_score length mismatch raises error", {
  dat <- .make_hemo_data()
  expect_error(fit_hemolysis_rr(dat$X_obs, dat$hemo[1:5]),
               "hemo_score length")
})

test_that("fit_hemolysis_rr: platelet proxy is accepted and stored", {
  dat <- .make_hemo_data(n = 60L, p = 20L)
  plt_score <- stats::runif(60L)
  fit <- fit_hemolysis_rr(dat$X_obs, dat$hemo, platelet_score = plt_score)
  expect_true(fit$has_platelet)
  expect_true("beta_plt" %in% colnames(fit$betas))
  expect_equal(ncol(fit$betas), 3L)
})

# ============================================================================
# apply_hemolysis_rr
# ============================================================================

test_that("apply_hemolysis_rr: post-correction error < 50% of pre-correction (headline behaviour)", {
  dat <- .make_hemo_data()
  fit <- fit_hemolysis_rr(dat$X_obs, dat$hemo)
  X_corr <- apply_hemolysis_rr(dat$X_obs, fit, dat$hemo)
  err_pre  <- mean(abs(log(dat$X_obs[, dat$shift_idx] + 1) -
                         log(dat$X_clean[, dat$shift_idx] + 1)))
  err_post <- mean(abs(X_corr[, dat$shift_idx] -
                         log(dat$X_clean[, dat$shift_idx] + 1)))
  expect_lt(err_post, err_pre * 0.75)
})

test_that("apply_hemolysis_rr: out-of-bound hemo scores flagged when report_oob=TRUE", {
  dat <- .make_hemo_data(n = 60L, p = 20L)
  fit <- fit_hemolysis_rr(dat$X_obs, dat$hemo)
  hemo_test <- c(rep(0.5, 10L), rep(999, 5L))   # 5 OOB samples
  X_test <- dat$X_obs[1:15, , drop = FALSE]
  out <- apply_hemolysis_rr(X_test, fit, hemo_test)
  oob <- attr(out, "oob")
  expect_length(oob, 15L)
  expect_true(all(oob[11:15]))    # OOB samples flagged
  expect_false(any(oob[1:10]))    # in-range samples not flagged
})

test_that("apply_hemolysis_rr: wrong class raises error", {
  expect_error(apply_hemolysis_rr(matrix(1), list(), hemo_score = 0),
               "hemolysis_rr_fit")
})


# ============================================================================
# hemolysis_index_blondal -- Blondal canonical index helper
# ============================================================================

test_that("hemolysis_index_blondal computes log(miR-451a) - log(miR-23a-3p)", {
  X_log <- matrix(c(8, 6, 4, 12, 5, 7, 3, 11), nrow = 2, byrow = FALSE,
                  dimnames = list(NULL, c("hsa-miR-451a", "hsa-miR-23a-3p",
                                          "hsa-miR-1", "hsa-miR-2")))
  idx <- hemolysis_index_blondal(X_log)
  expect_equal(idx, X_log[, "hsa-miR-451a"] - X_log[, "hsa-miR-23a-3p"])
})

test_that("hemolysis_index_blondal returns NULL when anchors missing", {
  X_log <- matrix(1:6, nrow = 2,
                  dimnames = list(NULL, c("hsa-miR-1", "hsa-miR-2", "hsa-miR-3")))
  expect_null(hemolysis_index_blondal(X_log))
})

test_that("hemolysis_index_blondal returns NULL when only one anchor present", {
  X_log <- matrix(1:6, nrow = 2,
                  dimnames = list(NULL, c("hsa-miR-451a", "hsa-miR-1", "hsa-miR-2")))
  expect_null(hemolysis_index_blondal(X_log))
})

test_that("hemolysis_index_blondal accepts MIMAT aliases", {
  X_log <- matrix(c(8, 6, 4, 12), nrow = 2,
                  dimnames = list(NULL, c("MIMAT0001631", "MIMAT0000078")))
  idx <- hemolysis_index_blondal(X_log)
  expect_equal(idx, X_log[, "MIMAT0001631"] - X_log[, "MIMAT0000078"])
})

test_that("hemolysis_index_blondal returns NULL when input has no colnames", {
  expect_null(hemolysis_index_blondal(matrix(1:6, nrow = 2)))
})
