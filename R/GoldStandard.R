#' @title Gold Standard Synthetic Dataset for Leakage Detection
#'
#' @description
#' Generates and evaluates synthetic datasets designed to detect data leakage
#' in feature selection and cross-validation pipelines.
#'
#' @details
#' The Gold Standard dataset includes:
#' - **Causal features**: Known ground truth signal correlated with outcome
#' - **Trap features**: Leakage detectors that appear informative only when
#'   test data improperly leaks into training
#' - **Batch features**: Nuisance features confounded with batch effects
#' - **Pathway blocks**: Correlated feature groups mimicking biological pathways
#'
#' A correctly implemented cross-validation pipeline should:
#' 1. Select mostly CAUSAL features
#' 2. Select ZERO TRAP features (selecting any indicates leakage)
#'
#' @name GoldStandard
NULL


#' Generate Gold Standard Synthetic Dataset
#'
#' Creates a synthetic dataset with known ground truth for validating
#' feature selection pipelines against data leakage.
#'
#' @param n_patients Number of patients (default: 100)
#' @param samples_per_patient Number of samples per patient (default: 3)
#' @param n_features Total number of features (default: 1000)
#' @param n_causal Number of true causal features (default: 20)
#' @param n_trap Number of trap/leakage detector features (default: 30)
#' @param n_batches Number of batch effects (default: 4)
#' @param seed Random seed for reproducibility (default: 42)
#' @param test_ratio Ratio of patients for test set (default: 0.3)
#'
#' @return A list of class `OmicSelectorGoldStandard` containing:
#' \describe{
#'   \item{data}{Full dataset (data.frame)}
#'   \item{train_data}{Training split (data.frame)}
#'   \item{test_data}{Test split (data.frame)}
#'   \item{train_indices}{Row indices for training}
#'   \item{test_indices}{Row indices for testing}
#'   \item{ground_truth}{Named list of feature indices and names}
#'   \item{train_patient_ids, test_patient_ids}{Patient IDs per split}
#'   \item{metadata}{Dataset generation parameters}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate a small dataset for testing
#' gs <- generate_gold_standard(
#'   n_patients = 50,
#'   samples_per_patient = 2,
#'   n_features = 200,
#'   seed = 42
#' )
#' print(gs)
#' }
#'
#' @export
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

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for generating correlated pathway blocks.", call. = FALSE)
  }

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
  timepoint <- rep(1:samples_per_patient, n_patients)

  # =========================================================================
  # 5. FEATURE GENERATION
  # =========================================================================
  feature_matrix <- matrix(
    rnorm(n_samples * n_features, mean = 0, sd = 1),
    nrow = n_samples,
    ncol = n_features
  )
  feature_names <- paste0("feature_", sprintf("%04d", 1:n_features))
  colnames(feature_matrix) <- feature_names

  # --- 5.1 TRUE CAUSAL FEATURES (1:n_causal) ---
  causal_indices <- 1:n_causal
  causal_effect_sizes <- runif(n_causal, min = 0.8, max = 2.0) * sample(c(-1, 1), n_causal, replace = TRUE)

  for (i in seq_along(causal_indices)) {
    idx <- causal_indices[i]
    feature_matrix[, idx] <- feature_matrix[, idx] + outcome * causal_effect_sizes[i]
    patient_effects <- rnorm(n_patients, mean = 0, sd = 0.5)
    feature_matrix[, idx] <- feature_matrix[, idx] + patient_effects[patient_id]
  }
  colnames(feature_matrix)[causal_indices] <- paste0("CAUSAL_", sprintf("%02d", 1:n_causal))

  # --- 5.2 TRAP FEATURES (LEAKAGE DETECTORS) ---
  trap_start <- n_causal + 1
  trap_indices <- trap_start:(trap_start + n_trap - 1)

  for (i in seq_along(trap_indices)) {
    idx <- trap_indices[i]
    trap_type <- (i %% 3) + 1

    if (trap_type == 1) {
      feature_matrix[train_indices, idx] <- rnorm(length(train_indices))
      feature_matrix[test_indices, idx] <- feature_matrix[test_indices, idx] + outcome[test_indices] * 2.5
    } else if (trap_type == 2) {
      feature_matrix[train_indices, idx] <- feature_matrix[train_indices, idx] - outcome[train_indices] * 0.3
      feature_matrix[test_indices, idx] <- feature_matrix[test_indices, idx] + outcome[test_indices] * 3.0
    } else {
      feature_matrix[train_indices, idx] <- rnorm(length(train_indices))
      feature_matrix[test_indices, idx] <- rnorm(length(test_indices)) + outcome[test_indices] * 2.0
    }
  }
  colnames(feature_matrix)[trap_indices] <- paste0("TRAP_", sprintf("%02d", 1:n_trap))

  # --- 5.3 BATCH-CORRELATED NUISANCE FEATURES ---
  batch_start <- trap_indices[length(trap_indices)] + 1
  n_batch_features <- 50
  batch_indices <- batch_start:(batch_start + n_batch_features - 1)

  for (i in seq_along(batch_indices)) {
    idx <- batch_indices[i]
    if (idx <= n_features) {
      feature_matrix[, idx] <- feature_matrix[, idx] + batch_effect_strength[batch] * runif(1, 0.5, 2.0)
    }
  }
  if (max(batch_indices) <= n_features) {
    colnames(feature_matrix)[batch_indices] <- paste0("BATCH_", sprintf("%02d", 1:n_batch_features))
  }

  # --- 5.4 CORRELATED PATHWAY BLOCKS ---
  n_pathway_blocks <- 10
  features_per_block <- 5
  pathway_start <- ifelse(max(batch_indices) <= n_features, max(batch_indices) + 1, batch_start)

  pathway_indices <- list()
  current_idx <- pathway_start

  for (block in 1:n_pathway_blocks) {
    if (current_idx + features_per_block - 1 <= n_features) {
      block_indices <- current_idx:(current_idx + features_per_block - 1)
      pathway_indices[[block]] <- block_indices

      correlation_strength <- runif(1, 0.5, 0.9)
      sigma <- matrix(correlation_strength, nrow = features_per_block, ncol = features_per_block)
      diag(sigma) <- 1

      correlated_block <- MASS::mvrnorm(n = n_samples, mu = rep(0, features_per_block), Sigma = sigma)
      feature_matrix[, block_indices] <- correlated_block
      colnames(feature_matrix)[block_indices] <- paste0("PATHWAY", block, "_", 1:features_per_block)

      current_idx <- current_idx + features_per_block
    }
  }

  # =========================================================================
  # 6. ADD BATCH EFFECT TO ALL FEATURES
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
#' Assesses whether a feature selection method has data leakage by checking
#' if TRAP features were selected.
#'
#' @param selected_features Character vector of selected feature names
#' @param gold_standard Output from [generate_gold_standard()]
#'
#' @return A list of class `OmicSelectorEvaluation` containing:
#' \describe{
#'   \item{n_causal_selected}{Number of true causal features found}
#'   \item{n_trap_selected}{Number of trap features found (should be 0)}
#'   \item{sensitivity}{Proportion of causal features recovered}
#'   \item{precision}{Proportion of selected features that are causal}
#'   \item{leakage_score}{Proportion of selected features that are traps}
#'   \item{verdict}{Human-readable assessment}
#'   \item{pass}{Logical indicating if selection passed (no traps)}
#' }
#'
#' @examples
#' \dontrun{
#' gs <- generate_gold_standard(n_features = 200, n_causal = 10, n_trap = 15)
#'
#' # A good selection: only causal features
#' good <- evaluate_feature_selection(gs$ground_truth$causal_features[1:5], gs)
#' print(good)  # PASS
#'
#' # A leaky selection: includes trap features
#' bad <- evaluate_feature_selection(c("CAUSAL_01", "TRAP_01"), gs)
#' print(bad)  # FAIL: Data leakage detected
#' }
#'
#' @export
evaluate_feature_selection <- function(selected_features, gold_standard) {
  gt <- gold_standard$ground_truth

  n_causal_selected <- sum(selected_features %in% gt$causal_features)
  n_trap_selected <- sum(selected_features %in% gt$trap_features)
  n_total_selected <- length(selected_features)
  n_causal_total <- length(gt$causal_features)

  sensitivity <- n_causal_selected / n_causal_total
  leakage_score <- if (n_total_selected > 0) n_trap_selected / n_total_selected else 0
  precision <- if (n_total_selected > 0) n_causal_selected / n_total_selected else 0

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
    n_causal_total = n_causal_total,
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


#' @export
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


#' @export
print.OmicSelectorEvaluation <- function(x, ...) {
  cat("=== Feature Selection Evaluation ===\n\n")
  cat(sprintf("Features selected: %d\n", x$n_total_selected))
  cat(sprintf("  - Causal (true positives): %d / %d\n",
              x$n_causal_selected, x$n_causal_total))
  cat(sprintf("  - Trap (leakage indicators): %d\n", x$n_trap_selected))
  cat("\nMetrics:\n")
  cat(sprintf("  Sensitivity: %.2f\n", x$sensitivity))
  cat(sprintf("  Precision: %.2f\n", x$precision))
  cat(sprintf("  Leakage Score: %.2f\n", x$leakage_score))
  cat("\n", x$verdict, "\n")
  invisible(x)
}
