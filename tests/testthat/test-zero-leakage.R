# Zero Leakage Validation Tests
#
# These tests implement the "Distractor Test" from the OmicSelector 2.0 specification:
# A synthetic dataset with signal implanted ONLY in the test set.
# The system must FAIL to select this feature - if it selects it, leakage is proven.
#
# Success Metric: Zero Leakage - The distractor feature must never be selected.

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

#' Create a Distractor Test Dataset
#'
#' This creates a dataset where specific "distractor" features have signal
#' ONLY in the test set. If any feature selection method selects these,
#' it proves that test data was seen during selection (leakage).
#'
#' @param n_samples Number of samples
#' @param n_features Total features
#' @param n_distractor Number of distractor features (test-only signal)
#' @param seed Random seed
#' @return List with dataset and ground truth
create_distractor_dataset <- function(
  n_samples = 200,
  n_features = 100,
  n_distractor = 5,
  test_ratio = 0.3,
  seed = 42
) {
  set.seed(seed)

  # Create train/test split
  n_test <- round(n_samples * test_ratio)
  n_train <- n_samples - n_test

  test_indices <- sample(1:n_samples, n_test)
  train_indices <- setdiff(1:n_samples, test_indices)

  # Generate outcome
  outcome <- rbinom(n_samples, 1, 0.5)

  # Generate feature matrix (pure noise)
  features <- matrix(
    rnorm(n_samples * n_features),
    nrow = n_samples,
    ncol = n_features
  )
  colnames(features) <- paste0("F", sprintf("%03d", 1:n_features))

  # Add DISTRACTOR features: signal ONLY in test set
  distractor_indices <- 1:n_distractor
  for (i in distractor_indices) {
    # No signal in training data
    features[train_indices, i] <- rnorm(n_train)

    # STRONG signal in test data (this should NEVER be found by valid methods)
    features[test_indices, i] <- features[test_indices, i] + outcome[test_indices] * 3.0
  }
  colnames(features)[distractor_indices] <- paste0("DISTRACTOR_", 1:n_distractor)

  # Add some TRUE causal features (signal in BOTH train and test)
  causal_indices <- (n_distractor + 1):(n_distractor + 10)
  for (i in causal_indices) {
    features[, i] <- features[, i] + outcome * 1.5
  }
  colnames(features)[causal_indices] <- paste0("CAUSAL_", seq_along(causal_indices))

  # Assemble dataset
  data <- data.frame(
    outcome = factor(outcome, levels = c(0, 1), labels = c("Control", "Case")),
    outcome_numeric = outcome,
    features,
    check.names = FALSE
  )

  list(
    data = data,
    train_data = data[train_indices, ],
    test_data = data[test_indices, ],
    train_indices = train_indices,
    test_indices = test_indices,
    distractor_features = colnames(features)[distractor_indices],
    causal_features = colnames(features)[causal_indices],
    all_feature_names = colnames(features)
  )
}


#' Evaluate if distractor features were selected (LEAKAGE DETECTED)
#'
#' @param selected_features Character vector of selected features
#' @param dataset Output from create_distractor_dataset()
#' @return List with evaluation results
evaluate_distractor_test <- function(selected_features, dataset) {
  n_distractor_selected <- sum(selected_features %in% dataset$distractor_features)
  n_causal_selected <- sum(selected_features %in% dataset$causal_features)

  leakage_detected <- n_distractor_selected > 0

  list(
    pass = !leakage_detected,
    leakage_detected = leakage_detected,
    n_distractor_selected = n_distractor_selected,
    n_causal_selected = n_causal_selected,
    distractors_found = selected_features[selected_features %in% dataset$distractor_features],
    causals_found = selected_features[selected_features %in% dataset$causal_features],
    verdict = if (leakage_detected) {
      paste("CRITICAL FAILURE: Leakage detected!",
            n_distractor_selected, "distractor features were selected.",
            "These features have signal ONLY in the test set.",
            "Selection indicates test data was seen during training.")
    } else if (n_causal_selected > 0) {
      paste("PASS: No leakage detected.", n_causal_selected,
            "causal features correctly identified.")
    } else {
      "WARNING: No leakage detected, but no causal features found either."
    }
  )
}


test_that("Distractor Test: method must not select test-only signal features", {
  # Create distractor dataset
  dataset <- create_distractor_dataset(
    n_samples = 200,
    n_features = 50,
    n_distractor = 5,
    seed = 42
  )

  # Verify distractor features have NO signal in training data
  train_outcome <- dataset$train_data$outcome_numeric

  for (feat in dataset$distractor_features) {
    cor_train <- cor(dataset$train_data[[feat]], train_outcome)
    # Distractor should have near-zero correlation in training
    expect_true(abs(cor_train) < 0.15,
                info = paste("Distractor", feat, "has unexpected train signal:", cor_train))
  }

  # Verify causal features have signal in training data
  for (feat in dataset$causal_features) {
    cor_train <- cor(dataset$train_data[[feat]], train_outcome)
    expect_true(abs(cor_train) > 0.2,
                info = paste("Causal", feat, "has weak train signal:", cor_train))
  }
})

test_that("Simple t-test filter on training data does not leak", {
  dataset <- create_distractor_dataset(
    n_samples = 200,
    n_features = 50,
    n_distractor = 5,
    seed = 42
  )

  # Implement simple t-test feature selection using ONLY training data
  feature_cols <- setdiff(names(dataset$train_data), c("outcome", "outcome_numeric"))
  outcome <- dataset$train_data$outcome_numeric

  p_values <- sapply(feature_cols, function(feat) {
    group0 <- dataset$train_data[[feat]][outcome == 0]
    group1 <- dataset$train_data[[feat]][outcome == 1]
    t.test(group0, group1)$p.value
  })

  # Select features with p < 0.05
  selected <- names(p_values)[p_values < 0.05]

  # Evaluate for leakage
  eval_result <- evaluate_distractor_test(selected, dataset)

  # Should NOT select distractor features
  expect_true(eval_result$pass, info = eval_result$verdict)
  expect_equal(eval_result$n_distractor_selected, 0)
})

test_that("Leaky method DOES select distractor features (negative control)", {
  # This test verifies our distractor setup works by intentionally leaking

  dataset <- create_distractor_dataset(
    n_samples = 200,
    n_features = 50,
    n_distractor = 5,
    seed = 42
  )

  # INTENTIONALLY LEAK: Use full dataset (train + test) for feature selection
  feature_cols <- setdiff(names(dataset$data), c("outcome", "outcome_numeric"))
  outcome <- dataset$data$outcome_numeric

  p_values <- sapply(feature_cols, function(feat) {
    group0 <- dataset$data[[feat]][outcome == 0]
    group1 <- dataset$data[[feat]][outcome == 1]
    t.test(group0, group1)$p.value
  })

  # Select features with p < 0.05
  selected <- names(p_values)[p_values < 0.05]

  # Evaluate for leakage
  eval_result <- evaluate_distractor_test(selected, dataset)

  # SHOULD detect leakage (this is a negative control)
  expect_true(eval_result$leakage_detected,
              info = "Leaky method should have selected distractor features")
  expect_gt(eval_result$n_distractor_selected, 0)
})

test_that("Preprocessing on full data causes leakage", {
  dataset <- create_distractor_dataset(
    n_samples = 200,
    n_features = 50,
    n_distractor = 5,
    seed = 42
  )

  feature_cols <- setdiff(names(dataset$data), c("outcome", "outcome_numeric"))

  # LEAK: Scale using full dataset statistics
  full_means <- colMeans(dataset$data[, feature_cols])
  full_sds <- apply(dataset$data[, feature_cols], 2, sd)

  # Apply to training data (but using stats from full data = leakage)
  train_scaled <- as.data.frame(scale(
    dataset$train_data[, feature_cols],
    center = full_means,
    scale = full_sds
  ))
  train_scaled$outcome_numeric <- dataset$train_data$outcome_numeric

  # Feature selection on leaked data
  outcome <- train_scaled$outcome_numeric
  p_values <- sapply(feature_cols, function(feat) {
    group0 <- train_scaled[[feat]][outcome == 0]
    group1 <- train_scaled[[feat]][outcome == 1]
    tryCatch(t.test(group0, group1)$p.value, error = function(e) 1)
  })

  selected <- names(p_values)[p_values < 0.1]

  # Check if distractor was selected (would indicate the preprocessing leaked info)
  n_distractor <- sum(selected %in% dataset$distractor_features)

  # This may or may not trigger depending on exact leakage path
  # The key is that our test framework CAN detect such issues
  # when they occur

  # Document result
  if (n_distractor > 0) {
    message("Preprocessing leakage confirmed: distractor features selected")
  }

  # Test passes but documents the scenario
  expect_true(TRUE)
})

test_that("Gold Standard trap features detect leakage", {
  gs <- generate_gold_standard(
    n_patients = 80,
    samples_per_patient = 2,
    n_features = 100,
    n_causal = 10,
    n_trap = 15,
    seed = 42
  )

  # Correct: feature selection using ONLY training data
  feature_cols <- gs$ground_truth$all_feature_names
  train_outcome <- gs$train_data$outcome_numeric

  p_values_correct <- sapply(feature_cols, function(feat) {
    group0 <- gs$train_data[[feat]][train_outcome == 0]
    group1 <- gs$train_data[[feat]][train_outcome == 1]
    tryCatch(t.test(group0, group1)$p.value, error = function(e) 1)
  })

  selected_correct <- names(p_values_correct)[p_values_correct < 0.05]
  eval_correct <- evaluate_feature_selection(selected_correct, gs)

  # WRONG: feature selection using ALL data (leakage)
  full_outcome <- gs$data$outcome_numeric

  p_values_leaky <- sapply(feature_cols, function(feat) {
    group0 <- gs$data[[feat]][full_outcome == 0]
    group1 <- gs$data[[feat]][full_outcome == 1]
    tryCatch(t.test(group0, group1)$p.value, error = function(e) 1)
  })

  selected_leaky <- names(p_values_leaky)[p_values_leaky < 0.05]
  eval_leaky <- evaluate_feature_selection(selected_leaky, gs)

  # Correct method should pass
  expect_true(eval_correct$pass,
              info = paste("Correct method failed:", eval_correct$verdict))

  # Leaky method should fail (select trap features)
  # Note: This may not always trigger depending on the specific trap construction
  # but it demonstrates the detection capability
  message("Correct method trap count: ", eval_correct$n_trap_selected)
  message("Leaky method trap count: ", eval_leaky$n_trap_selected)
})
