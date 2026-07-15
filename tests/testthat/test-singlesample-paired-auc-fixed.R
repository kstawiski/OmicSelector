test_that("paired AUC fixed orientation never reflects from held-out labels", {
  y <- c(0L, 0L, 1L, 1L)
  inverse_method <- c(4, 3, 2, 1)
  forward_baseline <- c(1, 2, 3, 4)

  fixed <- singlesample_paired_auc_diff_se(
    y, inverse_method, forward_baseline, orient = "fixed"
  )
  automatic <- singlesample_paired_auc_diff_se(
    y, inverse_method, forward_baseline, orient = "median"
  )

  expect_equal(fixed$auc_method, 0)
  expect_equal(fixed$auc_baseline, 1)
  expect_equal(fixed$lift, -1)
  expect_equal(automatic$auc_method, 1)
  expect_equal(automatic$auc_baseline, 1)
  expect_equal(automatic$lift, 0)
})

test_that("paired AUC defaults to the frozen supplied direction", {
  y <- c(0L, 0L, 1L, 1L)
  inverse_method <- c(4, 3, 2, 1)
  forward_baseline <- c(1, 2, 3, 4)

  default <- singlesample_paired_auc_diff_se(
    y, inverse_method, forward_baseline
  )

  expect_equal(default$auc_method, 0)
  expect_equal(default$auc_baseline, 1)
  expect_equal(default$lift, -1)
})

test_that("fixed orientation is respected independently within supplied folds", {
  y <- rep(c(0L, 0L, 1L, 1L), 2L)
  fold <- rep(1:2, each = 4L)
  method <- rep(c(4, 3, 2, 1), 2L)
  baseline <- rep(c(1, 2, 3, 4), 2L)

  out <- singlesample_paired_auc_diff_se(
    y, method, baseline, fold = fold, orient = "fixed"
  )

  expect_equal(out$n_folds, 2L)
  expect_equal(out$fold_results$auc_method, c(0, 0))
  expect_equal(out$fold_results$auc_baseline, c(1, 1))
  expect_equal(out$lift, -1)
})
