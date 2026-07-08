library(OmicSelector)

test_that("k-TSP scorer learns oriented pair votes", {
  X <- data.frame(
    A = c(5, 4, 3, 1, 1, 0),
    B = c(1, 1, 0, 4, 5, 6),
    C = c(0, 1, 0, 1, 0, 1)
  )
  y <- c(0, 0, 0, 1, 1, 1)
  fit <- os_ktsp_fit(X, y, k = 1)
  pred <- predict(fit, X)
  expect_s3_class(fit, "os_ktsp_model")
  expect_equal(length(pred), nrow(X))
  expect_gt(mean(pred[y == 1]), mean(pred[y == 0]))

  # Single-sample deployability: a 1-row newdata must return one scalar vote
  # fraction equal to its position in the batch (not the raw per-pair votes).
  fit5 <- os_ktsp_fit(X, y, k = 3)
  pred_batch <- predict(fit5, X)
  pred_one <- predict(fit5, X[1, , drop = FALSE])
  expect_length(pred_one, 1L)
  expect_equal(pred_one, pred_batch[1])
})

test_that("singscore preserves direction-split cardinality", {
  X <- matrix(c(
    5, 4, 1, 0,
    4, 3, 2, 1,
    1, 0, 4, 5
  ), nrow = 3, byrow = TRUE)
  colnames(X) <- c("up1", "up2", "down1", "down2")
  score <- os_singscore(X, up = c("up1", "up2"), down = c("down1", "down2"))
  expect_gt(score[1], score[3])
  expect_equal(attr(score, "up_n"), 2L)
  expect_equal(attr(score, "down_n"), 2L)
  expect_error(os_singscore(X, up = "up1", down = character()), "at least one")
})

test_that("Mahalanobis and conformal scorers use reference rows", {
  set.seed(1)
  controls <- matrix(rnorm(80, 0, 1), ncol = 4)
  cases <- matrix(rnorm(40, 3, 1), ncol = 4)
  X <- rbind(controls, cases)
  colnames(X) <- paste0("F", 1:4)
  mahal <- os_mahalanobis_score(X, reference_rows = 1:nrow(controls))
  conformal <- os_conformal_anomaly_score(
    X,
    train_rows = 1:10,
    calibration_rows = 11:20
  )
  expect_gt(mean(mahal[21:30]), mean(mahal[1:20]))
  expect_gt(mean(conformal[21:30]), mean(conformal[1:20]))
})

test_that("matched random-panel null excludes locked panel and emits QC", {
  set.seed(2)
  n <- 80L
  y <- c(rep(0, 40), rep(1, 40))
  X <- matrix(rnorm(n * 20), nrow = n)
  colnames(X) <- paste0("F", 1:20)
  X[y == 1, 1:3] <- X[y == 1, 1:3] + 2
  X[y == 0, 1:3] <- X[y == 0, 1:3] - 2
  scorer <- function(X_sub, features) rowMeans(X_sub)
  null <- os_panel_null_benchmark(
    X, y,
    panel_features = c("F1", "F2", "F3"),
    scorer = scorer,
    n_draws = 25,
    seed = 10
  )
  qc <- os_null_qc(null, min_draws = 25, expected_panel_size = 3)
  expect_s3_class(null, "os_panel_null")
  expect_true(all(qc$pass))
  expect_false(any(null$null_draws$overlaps_locked_panel))
  expect_gt(null$observed_stat, null$null_p95)
  expect_lte(null$permutation_p, 1)
})

test_that("claim gate fails closed under identifiability failure and applies Holm after pass", {
  results <- data.frame(
    method_id = c("m1", "m2"),
    p_value = c(0.01, 0.20),
    observed_stat = c(0.74, 0.62),
    null_p95 = c(0.68, 0.61),
    null_qc_pass = c(TRUE, TRUE)
  )
  closed <- os_claim_gate(results, identifiability_status = "FAIL", required_margin = 0.02)
  expect_true(all(closed$claim_label == "not_identifiable_descriptive_only"))
  open <- os_claim_gate(results, identifiability_status = "PASS", required_margin = 0.02)
  expect_equal(open$claim_label, c("primary_pass", "primary_fail"))
})

test_that("MDE margin record hashes fixed family inputs", {
  rec1 <- os_mde_margin_record(
    n = 100, n_clusters = 90,
    method_ids = c("m1", "m2"),
    null_draws = 1000,
    mde80 = 0.03,
    w90 = 0.08
  )
  rec2 <- os_mde_margin_record(
    n = 100, n_clusters = 90,
    method_ids = c("m1", "m2"),
    null_draws = 1000,
    mde80 = 0.03,
    w90 = 0.08
  )
  expect_s3_class(rec1, "os_mde_margin_record")
  expect_equal(rec1$required_margin, 0.04)
  expect_identical(rec1$hash, rec2$hash)
})
