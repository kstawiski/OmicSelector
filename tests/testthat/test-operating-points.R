library(OmicSelector)

test_that("operating points report sensitivity at specificity targets", {
  y <- c(rep(0, 100), rep(1, 40))
  score <- c(seq(0, 0.8, length.out = 100), seq(0.4, 1, length.out = 40))
  op <- os_operating_points(y, score, specificity_targets = c(0.90, 0.95))
  expect_s3_class(op, "os_operating_points")
  expect_equal(nrow(op$by_specificity), 2L)
  expect_true(all(op$by_specificity$specificity >= op$by_specificity$target))
})

test_that("calibrated brier supports grouped cross-fitting", {
  y <- rep(c(0, 1), each = 20)
  score <- y + rnorm(length(y), sd = 0.1)
  folds <- os_make_grouped_stratified_folds(y, n_folds = 4, seed = 2)
  brier <- os_calibrated_brier(y, score, folds = folds)
  expect_s3_class(brier, "os_calibrated_brier")
  expect_true(is.finite(brier$brier))
  expect_equal(length(brier$probabilities), length(y))
})

test_that("operating point gate refuses thresholds after identifiability failure", {
  y <- c(rep(0, 100), rep(1, 40))
  score <- c(seq(0, 0.8, length.out = 100), seq(0.4, 1, length.out = 40))
  op <- os_operating_points(y, score)
  gate <- os_operating_point_gate(op, identifiability_status = "FAIL-DEMOTED", min_specificity = 0.95)
  expect_equal(gate$status, "FAIL-CLOSED")
  expect_equal(gate$claim_label, "not_threshold_eligible_after_identifiability_failure")
})
