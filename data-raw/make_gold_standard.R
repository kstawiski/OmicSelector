#' Generate Gold Standard Synthetic Dataset for Leakage Detection
#'
#' This script generates the OmicSelector Gold Standard synthetic dataset designed
#' to detect data leakage in feature selection and cross-validation pipelines.
#'
#' Dataset Structure:
#' - 10,000 features total:
#'   - 20 TRUE causal features (known ground truth)
#'   - 50 correlated pathway blocks (nuisance but realistic)
#'   - 30 TRAP features (leakage detectors - correlate with outcome only if test data leaks)
#'   - ~9,900 pure noise features
#' - Patient hierarchies (multiple samples per patient)
#' - Batch/site effects
#' - Time/longitudinal structure
#'
#' @author OmicSelector 2.0 Team

library(MASS)  # For mvrnorm

#' Generate Gold Standard Synthetic Dataset
#'
#' @param n_patients Number of patients (default: 100)
#' @param samples_per_patient Number of samples per patient (default: 3)
#' @param n_features Total number of features (default: 1000 for CI speed)
#' @param n_causal Number of true causal features (default: 20)
#' @param n_trap Number of trap/leakage detector features (default: 30)
#' @param n_batches Number of batch effects (default: 4)
#' @param seed Random seed for reproducibility (default: 42)
#' @param test_ratio Ratio of patients for test set (default: 0.3)
#'
#' @return A list containing:
#'   - data: Full dataset (data.frame)
#'   - train_data: Training split (data.frame)
#'   - test_data: Test split (data.frame)
#'   - train_indices: Row indices for training
#'   - test_indices: Row indices for testing
#'   - ground_truth: Named list of feature indices
#'   - patient_splits: Patient IDs per split
#'   - metadata: Dataset generation parameters
generate_gold_standard <- function(
  n_patients = 100,
  samples_per_patient = 3,
  n_features = 1000,
  n_causal = 20,
  n_trap = 30,
  n_batches = 4,
  seed = 42,
  test_ratio = 0.3
) {
  set.seed(seed)

  n_samples <- n_patients * samples_per_patient

  # =========================================================================
  # 1. PATIENT STRUCTURE
  # =========================================================================
  patient_id <- rep(1:n_patients, each = samples_per_patient)
  sample_id <- 1:n_samples

  # Random patient-level outcome (binary classification)
  patient_outcome <- rbinom(n_patients, 1, 0.5)
  outcome <- patient_outcome[patient_id]

  # =========================================================================
  # 2. SPLIT BY PATIENT (prevents leakage from same-patient samples)
  # =========================================================================
  n_test_patients <- round(n_patients * test_ratio)
  test_patient_ids <- sample(1:n_patients, n_test_patients)
  train_patient_ids <- setdiff(1:n_patients, test_patient_ids)

  train_indices <- which(patient_id %in% train_patient_ids)
  test_indices <- which(patient_id %in% test_patient_ids)

  # =========================================================================
  # 3. BATCH/SITE EFFECTS
  # =========================================================================
  batch <- sample(1:n_batches, n_samples, replace = TRUE)
  batch_effect_strength <- rnorm(n_batches, mean = 0, sd = 1.5)

  # =========================================================================
  # 4. TIME/LONGITUDINAL STRUCTURE
  # =========================================================================
  # Samples within patient are timepoints
  timepoint <- rep(1:samples_per_patient, n_patients)

  # =========================================================================
  # 5. FEATURE GENERATION
  # =========================================================================

  # Initialize feature matrix
  feature_matrix <- matrix(
    rnorm(n_samples * n_features, mean = 0, sd = 1),
    nrow = n_samples,
    ncol = n_features
  )

  # Feature names
  feature_names <- paste0("feature_", sprintf("%04d", 1:n_features))
  colnames(feature_matrix) <- feature_names

  # --- 5.1 TRUE CAUSAL FEATURES (1:n_causal) ---
  causal_indices <- 1:n_causal
  causal_effect_sizes <- runif(n_causal, min = 0.8, max = 2.0) * sample(c(-1, 1), n_causal, replace = TRUE)

  for (i in seq_along(causal_indices)) {
    idx <- causal_indices[i]
    # Add outcome-dependent signal
    feature_matrix[, idx] <- feature_matrix[, idx] + outcome * causal_effect_sizes[i]
    # Add patient-level random effect
    patient_effects <- rnorm(n_patients, mean = 0, sd = 0.5)
    feature_matrix[, idx] <- feature_matrix[, idx] + patient_effects[patient_id]
  }

  # Rename causal features
  colnames(feature_matrix)[causal_indices] <- paste0("CAUSAL_", sprintf("%02d", 1:n_causal))

  # --- 5.2 TRAP FEATURES (LEAKAGE DETECTORS) ---
  # These correlate with outcome ONLY if you improperly include test data in training
  # CRITICAL: Training data alone should show NO signal for trap features
  trap_start <- n_causal + 1
  trap_indices <- trap_start:(trap_start + n_trap - 1)

  for (i in seq_along(trap_indices)) {
    idx <- trap_indices[i]
    trap_type <- (i %% 3) + 1  # Cycle through 3 trap types

    if (trap_type == 1) {
      # TYPE 1: NO signal in training, STRONG signal in test only
      # A correct method using only training data will NOT select this
      # A leaky method using train+test will select this
      feature_matrix[train_indices, idx] <- rnorm(length(train_indices))  # Pure noise
      feature_matrix[test_indices, idx] <- feature_matrix[test_indices, idx] +
        outcome[test_indices] * 2.5  # Strong test-only signal
    } else if (trap_type == 2) {
      # TYPE 2: Weak NEGATIVE correlation in train, STRONG POSITIVE in test
      # Combined: appears correlated, but train-only shows negative/weak
      feature_matrix[train_indices, idx] <- feature_matrix[train_indices, idx] -
        outcome[train_indices] * 0.3  # Weak negative
      feature_matrix[test_indices, idx] <- feature_matrix[test_indices, idx] +
        outcome[test_indices] * 3.0  # Strong positive - dominates when combined
    } else {
      # TYPE 3: Correlates with split membership (detects data mixing)
      # Pure noise in both sets when analyzed separately
      # But row-index correlation when combined
      feature_matrix[train_indices, idx] <- rnorm(length(train_indices))
      feature_matrix[test_indices, idx] <- rnorm(length(test_indices)) +
        outcome[test_indices] * 2.0  # Test-only signal
    }
  }

  # Rename trap features
  colnames(feature_matrix)[trap_indices] <- paste0("TRAP_", sprintf("%02d", 1:n_trap))

  # --- 5.3 BATCH-CORRELATED NUISANCE FEATURES ---
  batch_start <- trap_indices[length(trap_indices)] + 1
  n_batch_features <- 50
  batch_indices <- batch_start:(batch_start + n_batch_features - 1)

  for (i in seq_along(batch_indices)) {
    idx <- batch_indices[i]
    if (idx <= n_features) {
      # Correlates with batch, weak/no correlation with outcome
      feature_matrix[, idx] <- feature_matrix[, idx] +
        batch_effect_strength[batch] * runif(1, 0.5, 2.0)
    }
  }

  if (max(batch_indices) <= n_features) {
    colnames(feature_matrix)[batch_indices] <- paste0("BATCH_", sprintf("%02d", 1:n_batch_features))
  }

  # --- 5.4 CORRELATED PATHWAY BLOCKS ---
  # Create groups of correlated features (mimics biological pathways)
  n_pathway_blocks <- 10
  features_per_block <- 5
  pathway_start <- ifelse(max(batch_indices) <= n_features, max(batch_indices) + 1, batch_start)

  pathway_indices <- list()
  current_idx <- pathway_start

  for (block in 1:n_pathway_blocks) {
    if (current_idx + features_per_block - 1 <= n_features) {
      block_indices <- current_idx:(current_idx + features_per_block - 1)
      pathway_indices[[block]] <- block_indices

      # Create correlated block using Cholesky decomposition
      correlation_strength <- runif(1, 0.5, 0.9)
      sigma <- matrix(correlation_strength, nrow = features_per_block, ncol = features_per_block)
      diag(sigma) <- 1

      # Generate correlated features
      correlated_block <- mvrnorm(n = n_samples, mu = rep(0, features_per_block), Sigma = sigma)
      feature_matrix[, block_indices] <- correlated_block

      # Name them
      colnames(feature_matrix)[block_indices] <- paste0("PATHWAY", block, "_", 1:features_per_block)

      current_idx <- current_idx + features_per_block
    }
  }

  # =========================================================================
  # 6. ADD BATCH EFFECT TO ALL FEATURES (realistic preprocessing challenge)
  # =========================================================================
  for (j in 1:n_features) {
    feature_matrix[, j] <- feature_matrix[, j] + batch_effect_strength[batch] * 0.3
  }

  # =========================================================================
  # 7. ASSEMBLE FINAL DATASET
  # =========================================================================
  data <- data.frame(
    sample_id = sample_id,
    patient_id = patient_id,
    outcome = factor(outcome, levels = c(0, 1), labels = c("Control", "Case")),
    outcome_numeric = outcome,
    batch = factor(batch),
    timepoint = timepoint,
    feature_matrix,
    check.names = FALSE
  )

  # =========================================================================
  # 8. DEFINE GROUND TRUTH
  # =========================================================================
  ground_truth <- list(
    causal_features = colnames(feature_matrix)[causal_indices],
    causal_indices = causal_indices,
    causal_effect_sizes = causal_effect_sizes,
    trap_features = colnames(feature_matrix)[trap_indices],
    trap_indices = trap_indices,
    batch_features = if (max(batch_indices) <= n_features) colnames(feature_matrix)[batch_indices] else NULL,
    pathway_blocks = pathway_indices,
    all_feature_names = colnames(feature_matrix)
  )

  # =========================================================================
  # 9. CREATE TRAIN/TEST SPLITS
  # =========================================================================
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]

  # =========================================================================
  # 10. METADATA

  # =========================================================================
  metadata <- list(
    n_patients = n_patients,
    samples_per_patient = samples_per_patient,
    n_features = n_features,
    n_causal = n_causal,
    n_trap = n_trap,
    n_batches = n_batches,
    seed = seed,
    test_ratio = test_ratio,
    generation_date = Sys.time(),
    description = paste(
      "OmicSelector Gold Standard Synthetic Dataset v2.0",
      "Designed to detect data leakage in feature selection pipelines.",
      "TRUE CAUSAL features should be selected by correct methods.",
      "TRAP features should NOT be selected - selecting them indicates leakage.",
      sep = "\n"
    )
  )

  # =========================================================================
  # 11. RETURN RESULTS
  # =========================================================================
  result <- list(
    data = data,
    train_data = train_data,
    test_data = test_data,
    train_indices = train_indices,
    test_indices = test_indices,
    train_patient_ids = train_patient_ids,
    test_patient_ids = test_patient_ids,
    ground_truth = ground_truth,
    metadata = metadata
  )

  class(result) <- c("OmicSelectorGoldStandard", "list")

  return(result)
}


#' Evaluate Feature Selection for Leakage
#'
#' @param selected_features Character vector of selected feature names
#' @param gold_standard Output from generate_gold_standard()
#'
#' @return A list with evaluation metrics
evaluate_feature_selection <- function(selected_features, gold_standard) {
  gt <- gold_standard$ground_truth

  # Calculate metrics
  n_causal_selected <- sum(selected_features %in% gt$causal_features)
  n_trap_selected <- sum(selected_features %in% gt$trap_features)
  n_total_selected <- length(selected_features)

  # Sensitivity: proportion of causal features found
  sensitivity <- n_causal_selected / length(gt$causal_features)

  # Leakage score: proportion of selected features that are traps
  leakage_score <- if (n_total_selected > 0) n_trap_selected / n_total_selected else 0

  # Precision: proportion of selected features that are truly causal
  precision <- if (n_total_selected > 0) n_causal_selected / n_total_selected else 0

  # Verdict
  if (n_trap_selected > 0) {
    verdict <- "FAIL: Data leakage detected - TRAP features were selected"
    pass <- FALSE
  } else if (n_causal_selected == 0) {
    verdict <- "FAIL: No true causal features found"
    pass <- FALSE
  } else if (sensitivity < 0.5) {
    verdict <- "WARNING: Low sensitivity - many causal features missed"
    pass <- TRUE
  } else {
    verdict <- "PASS: Good feature selection without leakage"
    pass <- TRUE
  }

  result <- list(
    selected_features = selected_features,
    n_total_selected = n_total_selected,
    n_causal_selected = n_causal_selected,
    n_trap_selected = n_trap_selected,
    sensitivity = sensitivity,
    precision = precision,
    leakage_score = leakage_score,
    verdict = verdict,
    pass = pass,
    causal_found = selected_features[selected_features %in% gt$causal_features],
    traps_found = selected_features[selected_features %in% gt$trap_features]
  )

  class(result) <- c("OmicSelectorEvaluation", "list")
  return(result)
}


#' Print method for OmicSelectorGoldStandard
print.OmicSelectorGoldStandard <- function(x, ...) {
  cat("=== OmicSelector Gold Standard Dataset ===\n")
  cat("\nDimensions:\n")
  cat(sprintf("  Total samples: %d\n", nrow(x$data)))
  cat(sprintf("  Total features: %d\n", x$metadata$n_features))
  cat(sprintf("  Patients: %d\n", x$metadata$n_patients))
  cat(sprintf("  Samples per patient: %d\n", x$metadata$samples_per_patient))
  cat("\nSplit:\n")
  cat(sprintf("  Training: %d samples (%d patients)\n",
              nrow(x$train_data), length(x$train_patient_ids)))
  cat(sprintf("  Test: %d samples (%d patients)\n",
              nrow(x$test_data), length(x$test_patient_ids)))
  cat("\nGround Truth:\n")
  cat(sprintf("  Causal features: %d\n", length(x$ground_truth$causal_features)))
  cat(sprintf("  Trap features: %d (leakage detectors)\n", length(x$ground_truth$trap_features)))
  cat(sprintf("  Batch features: %d\n", length(x$ground_truth$batch_features)))
  cat("\nSeed: ", x$metadata$seed, "\n")
  invisible(x)
}


#' Print method for OmicSelectorEvaluation
print.OmicSelectorEvaluation <- function(x, ...) {
  cat("=== Feature Selection Evaluation ===\n\n")
  cat(sprintf("Features selected: %d\n", x$n_total_selected))
  cat(sprintf("  - Causal (true positives): %d / %d\n",
              x$n_causal_selected, 20))  # Hardcoded for simplicity
  cat(sprintf("  - Trap (leakage indicators): %d\n", x$n_trap_selected))
  cat("\nMetrics:\n")
  cat(sprintf("  Sensitivity: %.2f\n", x$sensitivity))
  cat(sprintf("  Precision: %.2f\n", x$precision))
  cat(sprintf("  Leakage Score: %.2f\n", x$leakage_score))
  cat("\n", x$verdict, "\n")
  invisible(x)
}


# =========================================================================
# GENERATE AND SAVE DATASETS
# =========================================================================

if (interactive() || !exists("SKIP_DATA_GENERATION")) {
  cat("Generating Gold Standard datasets...\n")

  # Small dataset for quick CI tests
  gold_standard_small <- generate_gold_standard(
    n_patients = 50,
    samples_per_patient = 2,
    n_features = 200,
    n_causal = 10,
    n_trap = 15,
    seed = 42
  )

  # Full dataset for thorough testing
  gold_standard_full <- generate_gold_standard(
    n_patients = 100,
    samples_per_patient = 3,
    n_features = 1000,
    n_causal = 20,
    n_trap = 30,
    seed = 42
  )

  cat("Saving datasets to data/...\n")

  save(gold_standard_small, file = "data/gold_standard_small.rda", compress = "xz")
  save(gold_standard_full, file = "data/gold_standard_full.rda", compress = "xz")

  # Also export the generator functions
  save(
    generate_gold_standard,
    evaluate_feature_selection,
    file = "data/gold_standard_functions.rda",
    compress = "xz"
  )

  cat("Done!\n")
  print(gold_standard_small)
}
