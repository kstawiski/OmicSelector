test_that("Platt fits logistic regression on the clipped probability at fit and prediction", {
  probs <- rep(c(0, 0.2, 0.5, 0.8, 1), each = 4)
  labels <- c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1)
  clipped <- pmax(pmin(probs, 1 - 1e-10), 1e-10)
  expected <- glm(labels ~ clipped, family = binomial())
  calibrator <- fit_platt_scaling(probs, labels)
  expect_equal(unname(attr(calibrator, "coefficients")), unname(coef(expected)))
  new_probs <- c(0, 1e-12, 0.25, 0.7, 1 - 1e-12, 1)
  expected_probs <- plogis(coef(expected)[1] + coef(expected)[2] *
                            pmax(pmin(new_probs, 1 - 1e-10), 1e-10))
  expect_equal(calibrator(new_probs), as.numeric(expected_probs))
  # predict.glm uses a matrix product; optimized BLAS may round one row and many rows differently.
  expect_equal(calibrator(new_probs), vapply(new_probs, calibrator, numeric(1)), tolerance = 1e-12)
  expect_error(fit_platt_scaling(c(-1, 0, 1), c(0, 0, 1)))
  expect_error(calibrator(c(-1, 2)))
})
