# Tests for Gold Standard Synthetic Dataset
# These tests verify the dataset structure and leakage detection capabilities

# Source the generator from data-raw (works in development and R CMD check)
# Try multiple paths to find the source file
source_paths <- c(
  "../../data-raw/make_gold_standard.R",  # From testthat context
  file.path(getwd(), "data-raw/make_gold_standard.R"),  # From package root
  file.path(dirname(getwd()), "data-raw/make_gold_standard.R")  # Fallback
)

for (path in source_paths) {
  if (file.exists(path)) {
    source(path)
    break
  }
}

# Skip tests if generator not available (e.g., in installed package without data-raw)
skip_if_not(
  exists("generate_gold_standard"),
  "Gold standard generator not available (data-raw not installed)"
)

test_that("Gold Standard dataset generates deterministically", {
  gs1 <- generate_gold_standard(seed = 123, n_patients = 20, n_features = 100)
  gs2 <- generate_gold_standard(seed = 123, n_patients = 20, n_features = 100)

  expect_equal(gs1$data, gs2$data)
  expect_equal(gs1$train_indices, gs2$train_indices)
  expect_equal(gs1$test_indices, gs2$test_indices)
})

test_that("Gold Standard dataset has correct structure", {
  gs <- generate_gold_standard(
    n_patients = 30,
    samples_per_patient = 2,
    n_features = 150,
    n_causal = 10,
    n_trap = 15,
    seed = 42
  )

  # Check dimensions
  expect_equal(nrow(gs$data), 30 * 2)
  expect_equal(nrow(gs$train_data) + nrow(gs$test_data), nrow(gs$data))

  # Check ground truth exists

  expect_length(gs$ground_truth$causal_features, 10)
  expect_length(gs$ground_truth$trap_features, 15)

  # Check no overlap between train and test patient IDs
  expect_length(intersect(gs$train_patient_ids, gs$test_patient_ids), 0)

  # Check all samples assigned to either train or test
  all_indices <- c(gs$train_indices, gs$test_indices)
  expect_equal(sort(all_indices), 1:nrow(gs$data))
})

test_that("Causal features correlate with outcome", {
  gs <- generate_gold_standard(
    n_patients = 100,
    samples_per_patient = 2,
    n_features = 100,
    n_causal = 10,
    n_trap = 10,
    seed = 42
  )

  # Get causal feature columns
  causal_cols <- gs$ground_truth$causal_features
  outcome <- gs$data$outcome_numeric

  # Test correlation for each causal feature
  for (feat in causal_cols) {
    cor_value <- cor(gs$data[[feat]], outcome)
    # Causal features should have meaningful correlation
    expect_true(abs(cor_value) > 0.1,
                info = paste("Causal feature", feat, "has weak correlation:", cor_value))
  }
})

test_that("Trap features detect leakage when selected", {
  gs <- generate_gold_standard(
    n_patients = 50,
    samples_per_patient = 2,
    n_features = 100,
    n_causal = 10,
    n_trap = 15,
    seed = 42
  )

  # Simulate a "leaky" feature selection that picks trap features
  leaky_selection <- c(
    gs$ground_truth$causal_features[1:5],
    gs$ground_truth$trap_features[1:3]  # This is leakage!
  )

  eval_result <- evaluate_feature_selection(leaky_selection, gs)

  # Should fail due to trap feature selection
  expect_false(eval_result$pass)
  expect_equal(eval_result$n_trap_selected, 3)
  expect_match(eval_result$verdict, "leakage")
})

test_that("Correct feature selection passes evaluation", {
  gs <- generate_gold_standard(
    n_patients = 50,
    samples_per_patient = 2,
    n_features = 100,
    n_causal = 10,
    n_trap = 15,
    seed = 42
  )

  # Simulate correct feature selection (only causal features)
  correct_selection <- gs$ground_truth$causal_features[1:7]

  eval_result <- evaluate_feature_selection(correct_selection, gs)

  # Should pass
  expect_true(eval_result$pass)
  expect_equal(eval_result$n_trap_selected, 0)
  expect_equal(eval_result$n_causal_selected, 7)
  expect_equal(eval_result$leakage_score, 0)
})

test_that("Patient-level split prevents sample leakage", {
  gs <- generate_gold_standard(
    n_patients = 50,
    samples_per_patient = 3,
    n_features = 50,
    seed = 42
  )

  # Get patient IDs for train and test samples
  train_patients <- unique(gs$train_data$patient_id)
  test_patients <- unique(gs$test_data$patient_id)

  # No patient should appear in both sets
  expect_length(intersect(train_patients, test_patients), 0)

  # All samples from same patient should be in same set
  for (pid in unique(gs$data$patient_id)) {
    sample_indices <- which(gs$data$patient_id == pid)
    in_train <- sample_indices %in% gs$train_indices
    in_test <- sample_indices %in% gs$test_indices

    # Either all in train or all in test
    expect_true(all(in_train) | all(in_test),
                info = paste("Patient", pid, "has samples in both train and test"))
  }
})

test_that("Trap features behave differently in train vs test", {
  gs <- generate_gold_standard(
    n_patients = 100,
    samples_per_patient = 2,
    n_features = 100,
    n_causal = 10,
    n_trap = 15,
    seed = 42
  )

  trap_cols <- gs$ground_truth$trap_features
  outcome_train <- gs$train_data$outcome_numeric
  outcome_test <- gs$test_data$outcome_numeric

  # For Type 1 traps (inverted in test), correlation should flip
  # For Type 2 traps (random in test), correlation should drop
  for (i in seq_along(trap_cols)) {
    feat <- trap_cols[i]
    cor_train <- cor(gs$train_data[[feat]], outcome_train)
    cor_test <- cor(gs$test_data[[feat]], outcome_test)

    # TRAP DESIGN (v2.0): Training data should show NO/WEAK correlation
    # This ensures a correct method using only training data won't select traps
    # Only when test data leaks in should traps be selected
    # Threshold 0.30 allows for random sampling variation
    expect_true(abs(cor_train) < 0.30,
                info = paste("Trap", feat, "should have weak/no train correlation, got:", cor_train))

    # Test correlation should be STRONG (this is where the signal is)
    # When train+test are combined, traps appear significant (leakage indicator)
    expect_true(abs(cor_test) > 0.2,
                info = paste("Trap", feat, "should have strong test correlation, got:", cor_test))
  }
})

test_that("Batch effects are present", {
  gs <- generate_gold_standard(
    n_patients = 50,
    samples_per_patient = 2,
    n_features = 150,
    n_batches = 4,
    seed = 42
  )

  # Check batch column exists and has correct levels
  expect_true("batch" %in% names(gs$data))
  expect_equal(length(unique(gs$data$batch)), 4)

  # Check batch features if they exist
  if (!is.null(gs$ground_truth$batch_features)) {
    batch_cols <- gs$ground_truth$batch_features
    batch <- as.numeric(as.character(gs$data$batch))

    # Batch features should correlate with batch
    for (feat in batch_cols[1:min(5, length(batch_cols))]) {
      if (feat %in% names(gs$data)) {
        # Use ANOVA to test batch effect
        model <- aov(gs$data[[feat]] ~ gs$data$batch)
        p_value <- summary(model)[[1]][["Pr(>F)"]][1]
        expect_true(p_value < 0.1,
                    info = paste("Batch feature", feat, "shows no batch effect"))
      }
    }
  }
})
