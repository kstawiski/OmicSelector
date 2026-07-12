library(testthat)

skip_if_not_installed("MASS")

.make_rank_test_compositions <- function(n = 80L, p = 13L, seed = 20260712L) {
  set.seed(seed)
  matrix(
    stats::rgamma(n * p, shape = 3, rate = 0.05),
    nrow = n,
    dimnames = list(NULL, paste0("miR-", seq_len(p)))
  )
}

test_that("rCLR Mahalanobis fits in an explicit D-1 coordinate basis", {
  skip_if_not_installed("robustbase")
  x <- .make_rank_test_compositions()

  expect_silent(
    fit <- fit_compositional_mahalanobis(x, transform = "rclr")
  )

  expect_equal(dim(fit$projection_basis), c(ncol(x), ncol(x) - 1L))
  expect_equal(dim(fit$cov), c(ncol(x) - 1L, ncol(x) - 1L))
  expect_equal(qr(fit$cov)$rank, ncol(x) - 1L)
  expect_equal(fit$df, qr(fit$cov)$rank)

  z_rclr <- .singlesample_compose_transform(
    x, transform = "rclr", pseudocount = fit$pseudocount
  ) %*% fit$projection_basis
  z_ilr <- .singlesample_compose_transform(
    x, transform = "ilr", pseudocount = fit$pseudocount
  )
  expect_equal(z_rclr, z_ilr, tolerance = 1e-12)
})

test_that("robust MCD fails clearly when n is not greater than dimension", {
  skip_if_not_installed("robustbase")
  x <- .make_rank_test_compositions(n = 10L, p = 13L)

  expect_error(
    fit_compositional_mahalanobis(
      x, transform = "rclr", require_robust = TRUE
    ),
    "robust MCD requires n_train > transformed dimension"
  )
})

test_that("classical fallback records the actual covariance rank", {
  skip_if_not_installed("robustbase")
  x <- .make_rank_test_compositions(n = 10L, p = 13L)

  expect_warning(
    fit <- fit_compositional_mahalanobis(
      x, transform = "rclr", require_robust = FALSE
    ),
    "falling back to classical covariance"
  )

  expect_equal(fit$df, qr(fit$cov)$rank)
  expect_lte(fit$df, nrow(x) - 1L)
})

test_that("projected rCLR scoring is single-row equivalent and serializable", {
  skip_if_not_installed("robustbase")
  x <- .make_rank_test_compositions()
  fit <- fit_compositional_mahalanobis(x, transform = "rclr")
  x_test <- .make_rank_test_compositions(n = 6L, seed = 20260713L)

  batch <- apply_compositional_mahalanobis(fit, x_test)
  singleton <- lapply(seq_len(nrow(x_test)), function(i) {
    apply_compositional_mahalanobis(fit, x_test[i, , drop = FALSE])
  })

  expect_equal(
    batch$distance_sq,
    vapply(singleton, function(z) z$distance_sq, numeric(1)),
    tolerance = 1e-12
  )

  restored <- unserialize(serialize(fit, NULL))
  expect_equal(
    apply_compositional_mahalanobis(restored, x_test)$distance_sq,
    batch$distance_sq,
    tolerance = 0
  )
})
