library(testthat)

set.seed(42L)

# ============================================================================
# preprocess_inverse_log
# ============================================================================

test_that("preprocess_inverse_log: output is strictly positive for finite log2 input", {
  set.seed(42L)
  log_mat <- matrix(rnorm(50 * 20, mean = 0, sd = 1.5), nrow = 50, ncol = 20)
  colnames(log_mat) <- paste0("S", seq_len(20))
  rownames(log_mat) <- paste0("miR-", seq_len(50))
  out <- preprocess_inverse_log(log_mat, base = 2)
  expect_true(all(out > 0))
  expect_equal(dim(out), dim(log_mat))
  expect_equal(attr(out, "preprocessing"), "inverse_log_base_2")
})

test_that("preprocess_inverse_log: round-trip x -> 2^x -> log2 recovers x to machine precision", {
  set.seed(42L)
  log_mat <- matrix(rnorm(30 * 15, mean = 0, sd = 1.5), nrow = 30, ncol = 15,
                     dimnames = list(paste0("f", seq_len(30)), paste0("S", seq_len(15))))
  out <- preprocess_inverse_log(log_mat, base = 2)
  rt <- log2(out)
  expect_lt(max(abs(rt - log_mat)), 1e-10)
})

test_that("preprocess_inverse_log: non-finite values trigger warning and are imputed", {
  set.seed(42L)
  mat <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
  mat[1, 3] <- Inf
  mat[5, 7] <- NA
  expect_warning(
    out <- preprocess_inverse_log(mat, base = 2),
    "non-finite values"
  )
  expect_equal(dim(out), dim(mat))
  expect_true(all(is.finite(out[1, ])))
})

test_that("preprocess_inverse_log: base 10 produces correct output", {
  set.seed(42L)
  x <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  out <- preprocess_inverse_log(x, base = 10)
  expect_equal(out[1, 1], 10^1)
  expect_equal(out[2, 2], 10^4)
})
