library(OmicSelector)

# -----------------------------------------------------------------------------
# os_log_transform_adaptive
# -----------------------------------------------------------------------------

test_that("os_log_transform_adaptive: half_min produces per-feature pseudocounts", {
  x <- matrix(c(0.001, 0.05, 0.13,
                0.002, 0.01, 0.20,
                0,     0.03, 1500),
              nrow = 3, byrow = TRUE)
  colnames(x) <- c("ratio_a", "ratio_b", "count_c")
  out <- os_log_transform_adaptive(x, method = "half_min")
  ps <- attr(out, "pseudocount")
  expect_equal(unname(ps[1]), 0.0005)             # half of 0.001 (smallest positive)
  expect_equal(unname(ps[2]), 0.005)              # half of 0.01
  expect_equal(unname(ps[3]), 0.065)              # half of 0.13 (smallest positive in count_c)
  expect_equal(unname(out[1, 1]), log2(0.001 + 0.0005))
  expect_equal(dim(out), dim(x))
})

test_that("os_log_transform_adaptive: fixed mode is backward-compatible with +1", {
  x <- matrix(c(0, 1, 100, 0, 9, 99), nrow = 2, byrow = TRUE)
  out <- os_log_transform_adaptive(x, method = "fixed", value = 1)
  expect_equal(out, log2(x + 1), ignore_attr = TRUE)
  expect_equal(attr(out, "pseudocount"), rep(1, ncol(x)))
})

test_that("os_log_transform_adaptive: epsilon mode uses 1e-6 everywhere", {
  x <- matrix(c(0.001, 1e6), nrow = 1)
  out <- os_log_transform_adaptive(x, method = "epsilon")
  expect_equal(unname(attr(out, "pseudocount")), c(1e-6, 1e-6))
})

test_that("os_log_transform_adaptive: column with no positives falls back to 1e-6", {
  x <- matrix(c(0, 5, 0, 10, 0, 7), nrow = 3, byrow = TRUE)
  colnames(x) <- c("zero_only", "ok")
  out <- os_log_transform_adaptive(x, method = "half_min")
  expect_equal(unname(attr(out, "pseudocount"))[1], 1e-6)
})

test_that("os_log_transform_adaptive: rejects negative inputs", {
  x <- matrix(c(-1, 1, 2, 3), nrow = 2)
  expect_error(os_log_transform_adaptive(x), "negative values")
})


# -----------------------------------------------------------------------------
# os_plate_median_frozen
# -----------------------------------------------------------------------------

test_that("os_plate_median_frozen: removes per-plate location shift", {
  set.seed(7)
  plate <- rep(c("P1", "P2", "P3"), each = 20)
  x <- matrix(stats::rnorm(60 * 4, mean = 0, sd = 0.3), nrow = 60, ncol = 4)
  shift <- c(P1 = 0, P2 = 5, P3 = -3)[plate]
  x_shifted <- x + matrix(shift, nrow = 60, ncol = 4)
  fc <- os_plate_median_frozen$new()
  fc$fit(x_shifted, plate)
  corrected <- fc$transform(x_shifted, plate)
  # After correction, per-plate medians should be approximately equal across plates.
  med_p1 <- apply(corrected[plate == "P1", ], 2L, median)
  med_p2 <- apply(corrected[plate == "P2", ], 2L, median)
  med_p3 <- apply(corrected[plate == "P3", ], 2L, median)
  expect_equal(med_p1, med_p2, tolerance = 1e-8)
  expect_equal(med_p1, med_p3, tolerance = 1e-8)
})

test_that("os_plate_median_frozen: leaves identical-batch data approximately unchanged", {
  set.seed(11)
  plate <- rep("P1", 30)
  x <- matrix(stats::rnorm(30 * 5, mean = 2, sd = 1), nrow = 30, ncol = 5)
  fc <- os_plate_median_frozen$new()
  fc$fit(x, plate)
  out <- fc$transform(x, plate)
  # With one plate, transform is identity (per-plate median == grand median).
  expect_equal(out, x, tolerance = 1e-12)
})

test_that("os_plate_median_frozen: scale='iqr' equalises per-plate IQR", {
  set.seed(13)
  plate <- rep(c("A", "B"), each = 50)
  x <- matrix(stats::rnorm(100 * 3), nrow = 100)
  x[plate == "B", ] <- x[plate == "B", ] * 4   # B has 4x the spread
  fc <- os_plate_median_frozen$new(scale = "iqr")
  fc$fit(x, plate)
  out <- fc$transform(x, plate)
  iqr_a <- apply(out[plate == "A", ], 2, IQR)
  iqr_b <- apply(out[plate == "B", ], 2, IQR)
  # After IQR-correction the two plates should have similar spreads.
  expect_lt(max(abs(iqr_a - iqr_b) / iqr_a), 0.05)
})

test_that("os_plate_median_frozen: unseen plate triggers warning + grand-median fallback", {
  set.seed(17)
  plate_train <- rep(c("P1", "P2"), each = 15)
  x_train <- matrix(stats::rnorm(30 * 2), nrow = 30)
  fc <- os_plate_median_frozen$new()
  fc$fit(x_train, plate_train)
  x_test <- matrix(stats::rnorm(5 * 2), nrow = 5)
  expect_warning(fc$transform(x_test, rep("P_unseen", 5)), "not seen during fit")
})

test_that("os_plate_median_frozen: transform() before fit() errors", {
  fc <- os_plate_median_frozen$new()
  expect_error(fc$transform(matrix(1, 2, 2), c("a", "b")), "before")
})

test_that("os_plate_median_frozen: get_plate_levels() returns sorted training levels", {
  fc <- os_plate_median_frozen$new()
  fc$fit(matrix(stats::rnorm(40), nrow = 10), rep(c("P_b", "P_a"), each = 5))
  expect_equal(fc$get_plate_levels(), c("P_a", "P_b"))
})


# -----------------------------------------------------------------------------
# os_paired_delong
# -----------------------------------------------------------------------------

test_that("os_paired_delong: identical scores yield p ~ 1 / NA, delta == 0", {
  skip_if_not_installed("pROC")
  set.seed(7)
  y <- rep(c(0, 1), each = 50)
  s <- y + stats::rnorm(100, sd = 0.5)
  out <- os_paired_delong(y, s, s)
  expect_equal(out$delta, 0)
  expect_true(is.na(out$p_value) || out$p_value > 0.99)
})

test_that("os_paired_delong: informative beats random with p<0.05 and positive delta", {
  skip_if_not_installed("pROC")
  set.seed(11)
  y <- rep(c(0, 1), each = 100)
  s_good <- y + stats::rnorm(200, sd = 0.5)
  s_rand <- stats::rnorm(200)
  out <- os_paired_delong(y, s_good, s_rand, alternative = "greater")
  expect_gt(out$delta, 0)
  expect_lt(out$p_value, 0.05)
  expect_equal(out$auc_a, as.numeric(pROC::auc(pROC::roc(y, s_good, quiet = TRUE, direction = "<"))))
})

test_that("os_paired_delong: requires both classes", {
  skip_if_not_installed("pROC")
  expect_error(os_paired_delong(rep(1, 10), stats::rnorm(10), stats::rnorm(10)),
               "both classes")
})


# -----------------------------------------------------------------------------
# os_oof_pipeline_compare
# -----------------------------------------------------------------------------

make_synth <- function(N = 120, seed = 7) {
  set.seed(seed)
  groups <- as.character(seq_len(N))
  y <- rep(c(0, 1), length.out = N)
  X_signal <- matrix(stats::rnorm(N * 3, mean = 0, sd = 1), nrow = N)
  X_signal[y == 1, ] <- X_signal[y == 1, ] + 0.8
  list(y = y, groups = groups, X = X_signal,
       batch = sample(c("b1", "b2", "b3"), N, replace = TRUE),
       plate = sample(c("p1", "p2"), N, replace = TRUE))
}

pipe_signal <- function(Xtr, ytr, Xte, batch_train, batch_test, plate_train, plate_test) {
  # Use feature 1 (signal) as the score. No fitting; mimics a pre-trained scorer.
  list(train = Xtr[, 1L], test = Xte[, 1L])
}
pipe_random <- function(Xtr, ytr, Xte, batch_train, batch_test, plate_train, plate_test) {
  list(train = stats::runif(nrow(Xtr)), test = stats::runif(nrow(Xte)))
}

test_that("os_oof_pipeline_compare: signal pipeline beats random; OOF count == N", {
  skip_if_not_installed("pROC")
  d <- make_synth()
  res <- os_oof_pipeline_compare(
    y = d$y, groups = d$groups,
    pipelines = list(signal = pipe_signal, random = pipe_random),
    X = d$X, batch = d$batch, plate = d$plate,
    n_folds = 5L, seeds = c(7L, 42L), reference = "random",
    verbose = FALSE
  )
  # OOF prediction count == N for each pipeline per seed.
  for (seed in c("7", "42")) {
    expect_equal(length(res$oof_predictions[[seed]]$signal), length(d$y))
    expect_equal(length(res$oof_predictions[[seed]]$random), length(d$y))
    # No NA in OOF predictions (every sample is in exactly one test fold).
    expect_true(!anyNA(res$oof_predictions[[seed]]$signal))
  }
  # Signal must beat random in expectation.
  s_signal <- res$summary$mean_auc[res$summary$pipeline == "signal"]
  s_random <- res$summary$mean_auc[res$summary$pipeline == "random"]
  expect_gt(s_signal, s_random)
  # Pairwise DeLong table contains 1 row per (non-reference pipeline, seed).
  expect_equal(nrow(res$pairwise_delong), 2L)
  expect_equal(unique(res$pairwise_delong$pipeline), "signal")
})

test_that("os_oof_pipeline_compare: errors when reference name is unknown", {
  d <- make_synth()
  expect_error(
    os_oof_pipeline_compare(
      y = d$y, groups = d$groups,
      pipelines = list(signal = pipe_signal),
      X = d$X, reference = "no_such_pipeline"
    ),
    "not in names"
  )
})

test_that("os_oof_pipeline_compare: pairwise_summary is NULL when reference is NULL", {
  skip_if_not_installed("pROC")
  d <- make_synth(N = 60)
  res <- os_oof_pipeline_compare(
    y = d$y, groups = d$groups,
    pipelines = list(signal = pipe_signal),
    X = d$X, n_folds = 3L, seeds = 7L, reference = NULL
  )
  expect_null(res$pairwise_delong)
  expect_null(res$pairwise_summary)
})


# -----------------------------------------------------------------------------
# Internal helper: stratified-grouped folds
# -----------------------------------------------------------------------------

test_that(".grouped_stratified_folds: each row in exactly one fold", {
  set.seed(7)
  N <- 100
  y <- rep(c(0, 1), 50)
  groups <- as.character(seq_len(N))
  folds <- OmicSelector:::.grouped_stratified_folds(y, groups, n_folds = 5L)
  expect_equal(unname(sort(unlist(folds))), seq_len(N))
  expect_equal(length(folds), 5L)
})
