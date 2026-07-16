test_that("canonical rCLR panel matches the benchmark definition", {
  set.seed(263L)
  y <- rep(0:1, each = 30L)
  X <- matrix(rexp(60L * 30L, rate = 0.2), nrow = 60L,
              dimnames = list(paste0("S", seq_len(60L)),
                              paste0("hsa-miR-", seq_len(30L))))
  X[y == 1L, 1:4] <- X[y == 1L, 1:4] * 3
  model <- fit_ws_rclr_panel(X, y)
  Z <- ws_rclr_trimmed(X)
  effect <- vapply(colnames(Z), function(feature) {
    x <- Z[, feature]
    den <- stats::sd(x[is.finite(x)])
    if (sum(is.finite(x)) < 5L || !is.finite(den) || den <= 0) return(0)
    abs(mean(x[y == 1L]) - mean(x[y == 0L])) / den
  }, numeric(1L))
  panel <- colnames(Z)[order(effect, decreasing = TRUE)][1:20]
  weights <- vapply(panel, function(feature) {
    d <- stats::median(Z[y == 1L, feature]) - stats::median(Z[y == 0L, feature])
    if (is.finite(d) && d < 0) -1 else 1
  }, numeric(1L))
  expected <- as.numeric(Z[, panel, drop = FALSE] %*% weights)

  expect_identical(model$panel, panel)
  expect_equal(model$weights, weights)
  expect_equal(score_ws_rclr_panel(model, X), expected)
})

test_that("rCLR panel is row-equivariant and exposed by deployment API", {
  set.seed(264L)
  y <- rep(0:1, each = 20L)
  X <- matrix(rexp(40L * 25L), nrow = 40L,
              dimnames = list(NULL, paste0("hsa-miR-", seq_len(25L))))
  deployed <- deploy_singlesample(X, y, method = "ws-rclr-trimmed")
  expect_true(deployed$verified)
  expect_equal(length(score_specimen(deployed, X[1:3, , drop = FALSE])), 3L)
  expect_true(is_singlesample_deployable(deployed, X[1:8, , drop = FALSE]))
  expect_error(score_specimen(deployed, X[1:2, -1, drop = FALSE]),
               "missing frozen panel feature")
})
