# Tests for Signature Selection Module

test_that("select_best_signature with constrained_1se mode works", {
  skip_if_not_installed("mlr3")

  # Create mock NestedCVResult
  mock_result <- create_mock_nested_cv_result()

  # Test 1SE selection
  best <- select_best_signature(mock_result, mode = "constrained_1se")

  expect_s3_class(best, "data.table")
  expect_true("selected" %in% names(best))
  expect_true("mean_metric" %in% names(best))
  expect_true("se_metric" %in% names(best))
  expect_equal(sum(best$selected), 1)  # Exactly one selected
})

test_that("select_best_signature with weighted mode works", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result()

  best <- select_best_signature(
    mock_result,
    mode = "weighted",
    weights = c(performance = 0.5, stability = 0.3, parsimony = 0.2)
  )

  expect_s3_class(best, "data.table")
  expect_true("score" %in% names(best))
  expect_true("selected" %in% names(best))
  expect_equal(sum(best$selected), 1)

  # Score should be between 0 and 1 (normalized)
  expect_true(all(best$score >= 0 & best$score <= 1))
})

test_that("select_best_signature with pareto mode works", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result(n_learners = 5)

  best <- select_best_signature(mock_result, mode = "pareto")

  expect_s3_class(best, "data.table")
  expect_true("pareto" %in% names(best))
  expect_true(any(best$pareto))  # At least one Pareto-optimal point
})

test_that("select_best_signature applies AUC constraint correctly", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result()

  # Get all candidates first
  all_candidates <- select_best_signature(mock_result, mode = "constrained_1se")
  min_auc <- min(all_candidates$mean_metric)
  max_auc <- max(all_candidates$mean_metric)

  # Apply constraint that filters some candidates
  threshold <- (min_auc + max_auc) / 2
  filtered <- select_best_signature(mock_result, mode = "constrained_1se", auc_min = threshold)

  expect_true(nrow(filtered) <= nrow(all_candidates))
  expect_true(all(filtered$mean_metric >= threshold))
})

test_that("select_best_signature applies k_max constraint correctly", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result(with_features = TRUE)

  # Apply strict feature count constraint
  filtered <- select_best_signature(mock_result, mode = "constrained_1se", k_max = 5)

  if (nrow(filtered) > 0) {
    expect_true(all(filtered$mean_k <= 5 | is.na(filtered$mean_k)))
  }
})

test_that("Pareto frontier computation is correct", {
  # Create test data with known Pareto points
  dt <- data.table::data.table(
    id = c("A", "B", "C", "D", "E"),
    metric = c(0.9, 0.85, 0.8, 0.7, 0.6),  # higher better
    stability = c(0.5, 0.7, 0.8, 0.9, 0.85),  # higher better
    k = c(20, 15, 10, 5, 8)  # lower better
  )

  # A dominates nothing, B dominates nothing, C dominates E, D dominates nothing
  # Pareto set should be: A, B, C, D (E is dominated by C)
  idx <- .pareto_frontier_idx(dt, cols = c("metric", "stability", "k"),
                              maximize = c(TRUE, TRUE, FALSE))

  pareto_ids <- dt$id[idx]
  expect_true("A" %in% pareto_ids)
  expect_true("B" %in% pareto_ids)
  expect_true("D" %in% pareto_ids)
  expect_false("E" %in% pareto_ids)  # E is dominated by C (same or better on all)
})

test_that("Pareto frontier handles simple 2D case", {
  dt <- data.table::data.table(
    id = c("P1", "P2", "P3"),
    x = c(1, 2, 3),  # higher better
    y = c(3, 2, 1)   # higher better
  )

  # All points are on Pareto frontier (trade-off between x and y)
  idx <- .pareto_frontier_idx(dt, cols = c("x", "y"), maximize = c(TRUE, TRUE))
  expect_equal(length(idx), 3)
})

test_that("minmax normalization works correctly", {
  x <- c(10, 20, 30, 40, 50)
  normalized <- .minmax_normalize(x)

  expect_equal(min(normalized), 0)
  expect_equal(max(normalized), 1)
  expect_equal(normalized, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("minmax normalization handles constant values", {
  x <- c(5, 5, 5, 5)
  normalized <- .minmax_normalize(x)

  # All values should be 1 when constant
  expect_equal(normalized, c(1, 1, 1, 1))
})

test_that("minmax normalization handles NA values", {
  x <- c(10, NA, 30, 40, 50)
  normalized <- .minmax_normalize(x)

  expect_equal(length(normalized), 5)
  expect_false(any(is.na(normalized[c(1, 3, 4, 5)])))
})

test_that("get_consensus_features works with mock data", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result(with_features = TRUE)
  best <- select_best_signature(mock_result, mode = "constrained_1se")

  if (nrow(best) > 0 && "learner_id" %in% names(best)) {
    learner_id <- best[selected == TRUE]$learner_id[1]

    consensus <- get_consensus_features(mock_result, learner_id, min_frequency = 0.5)

    expect_s3_class(consensus, "data.table")
    expect_true("feature" %in% names(consensus))
    expect_true("frequency" %in% names(consensus))
    expect_true("selected" %in% names(consensus))
  }
})

test_that("select_best_signature handles empty result gracefully", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result()

  # Apply impossible constraint
  result <- select_best_signature(mock_result, mode = "constrained_1se", auc_min = 2.0)

  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 0)
})

test_that("weights must sum to positive value", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result()

  # Weights are normalized internally, so any positive weights work
  result <- select_best_signature(
    mock_result,
    mode = "weighted",
    weights = c(performance = 1, stability = 0.1, parsimony = 0.1)
  )

  expect_s3_class(result, "data.table")
})

test_that("weighted mode requires correct weight names", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result()

  expect_error(
    select_best_signature(
      mock_result,
      mode = "weighted",
      weights = c(auc = 0.5, stab = 0.3)  # Wrong names
    ),
    "performance"
  )
})

test_that("selection attributes are preserved", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result()

  result <- select_best_signature(
    mock_result,
    mode = "weighted",
    auc_min = 0.5,
    k_max = 20
  )

  expect_equal(attr(result, "mode"), "weighted")
  expect_equal(attr(result, "constraints")$auc_min, 0.5)
  expect_equal(attr(result, "constraints")$k_max, 20)
})

test_that("NA stability handling works in 1SE mode", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result(with_na_stability = TRUE)

  # Should not error with NA stability when not constrained
  result <- select_best_signature(mock_result, mode = "constrained_1se", na_stability = "worst")

  expect_s3_class(result, "data.table")
})

test_that("1SE selection keeps candidates within 1 SE of best", {
  skip_if_not_installed("mlr3")

  mock_result <- create_mock_nested_cv_result(n_learners = 10)

  result <- select_best_signature(mock_result, mode = "constrained_1se")

  if (nrow(result) > 1) {
    best <- result[which.max(mean_metric)]
    threshold <- best$mean_metric - best$se_metric

    if (!is.na(threshold)) {
      expect_true(all(result$mean_metric >= threshold - 1e-10))
    }
  }
})


# Helper function to create mock NestedCVResult
create_mock_nested_cv_result <- function(n_learners = 3, n_folds = 5,
                                          with_features = FALSE,
                                          with_na_stability = FALSE) {
  set.seed(42)

  # Create mock benchmark result
  task <- mlr3::tsk("iris")
  learner_ids <- paste0("learner_", seq_len(n_learners))

  # Create mock score data
  score_data <- data.table::data.table(
    learner_id = rep(learner_ids, each = n_folds),
    task_id = "iris",
    resampling_id = "cv",
    iteration = rep(seq_len(n_folds), n_learners),
    classif.auc = runif(n_learners * n_folds, 0.6, 0.95),
    classif.acc = runif(n_learners * n_folds, 0.7, 0.95)
  )

  # Create mock benchmark result object
  mock_bmr <- list(
    score = function() score_data,
    n_resample_results = n_learners,
    resample_result = function(i) {
      list(
        learner = list(id = learner_ids[i]),
        task = task,
        iters = n_folds
      )
    }
  )
  class(mock_bmr) <- c("BenchmarkResult", "R6")

  # Create stability info
  if (with_na_stability) {
    stab_vals <- c(runif(n_learners - 1, 0.5, 0.9), NA)
  } else {
    stab_vals <- runif(n_learners, 0.5, 0.9)
  }

  stability <- list(
    nogueira_index = data.table::data.table(
      learner_id = learner_ids,
      nogueira_index = stab_vals
    ),
    selected_features_per_fold = if (with_features) {
      # Create mock feature selection data
      sf_data <- data.table::data.table(
        learner_id = rep(learner_ids, each = n_folds),
        outer_fold = rep(seq_len(n_folds), n_learners),
        selected_features = lapply(seq_len(n_learners * n_folds), function(i) {
          n_feat <- sample(5:15, 1)
          paste0("gene_", sample(100, n_feat))
        })
      )
      sf_data
    } else {
      list()
    }
  )

  # Build NestedCVResult
  result <- list(
    benchmark_result = mock_bmr,
    performance = score_data[, .(
      mean_auc = mean(classif.auc),
      mean_acc = mean(classif.acc)
    ), by = .(learner_id, task_id, resampling_id)],
    per_fold = list(),
    stability = stability,
    outer_folds = n_folds,
    inner_folds = 3,
    seed = 42,
    timestamp = Sys.time()
  )

  class(result) <- c("NestedCVResult", "list")
  result
}
