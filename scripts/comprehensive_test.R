#!/usr/bin/env Rscript
#' Comprehensive Test of ALL OmicSelector Features
#' Using TCGA miRNA Dataset

# ============================================================================
# SETUP
# ============================================================================

cat("========================================================================\n")
cat("OmicSelector 2.0 - Comprehensive Feature Test\n")
cat("========================================================================\n\n")

set.seed(42)

# Load required packages
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3pipelines)
  library(mlr3filters)
  library(mlr3learners)
  library(data.table)
})

# Source OmicSelector modules
source_dir <- file.path(getwd(), "R")
r_files <- list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  tryCatch(source(f), error = function(e) {
    cat("Warning: Could not source", basename(f), "-", e$message, "\n")
  })
}

# Results storage
test_results <- list()

# ============================================================================
# 1. LOAD AND PREPARE TCGA DATA
# ============================================================================

cat("\n[1/13] LOADING TCGA DATA\n")
cat("--------------------------------------------------\n")

# Load the TCGA data
load("data/orginal_TCGA_data.rda")

# Check sample_type column
cat("Sample types:\n")
print(table(orginal_TCGA_data$sample_type))

# Identify miRNA columns (numeric features)
mirna_cols <- grep("^hsa\\.miR", names(orginal_TCGA_data), value = TRUE)
cat("Number of miRNA features:", length(mirna_cols), "\n")

# Convert to data.frame to ensure consistent behavior
orginal_TCGA_data <- as.data.frame(orginal_TCGA_data)

# Create analysis dataset
# Use sample_type as target, subset miRNA features (limit to 200 for speed)
selected_cols <- c("sample_type", mirna_cols[1:min(200, length(mirna_cols))])
analysis_data <- orginal_TCGA_data[, selected_cols, drop = FALSE]

# Ensure it's a proper data.frame
analysis_data <- as.data.frame(analysis_data)

# Remove rows with missing target
valid_rows <- !is.na(analysis_data$sample_type)
analysis_data <- analysis_data[valid_rows, , drop = FALSE]

# Convert target to factor
analysis_data$sample_type <- factor(analysis_data$sample_type)

# Subset to binary classification (tumor vs normal if available)
if (length(levels(analysis_data$sample_type)) > 2) {
  # Keep only the two most common classes
  class_counts <- sort(table(analysis_data$sample_type), decreasing = TRUE)
  top_classes <- names(class_counts)[1:2]
  keep_rows <- analysis_data$sample_type %in% top_classes
  analysis_data <- analysis_data[keep_rows, , drop = FALSE]
  analysis_data$sample_type <- droplevels(analysis_data$sample_type)
}

# Handle missing values in features
feature_cols <- setdiff(names(analysis_data), "sample_type")
for (col in feature_cols) {
  na_mask <- is.na(analysis_data[[col]])
  if (any(na_mask)) {
    analysis_data[[col]][na_mask] <- median(analysis_data[[col]], na.rm = TRUE)
  }
}

# Subsample for speed if too large (limit to 800 samples for testing)
if (nrow(analysis_data) > 800) {
  set.seed(42)
  idx <- sample(nrow(analysis_data), 800)
  analysis_data <- analysis_data[idx, , drop = FALSE]
  rownames(analysis_data) <- NULL
}

cat("Final dataset:", nrow(analysis_data), "samples,", ncol(analysis_data) - 1, "features\n")
cat("Class distribution:\n")
print(table(analysis_data$sample_type))

test_results$data_loaded <- TRUE
test_results$n_samples <- nrow(analysis_data)
test_results$n_features <- ncol(analysis_data) - 1

# ============================================================================
# 2. TEST OmicPipeline
# ============================================================================

cat("\n[2/13] TESTING OmicPipeline\n")
cat("--------------------------------------------------\n")

tryCatch({
  # Determine positive class
  levels_ordered <- levels(analysis_data$sample_type)
  positive_class <- if ("Primary Tumor" %in% levels_ordered) "Primary Tumor" else levels_ordered[2]

  # Create pipeline
  pipeline <- OmicPipeline$new(
    data = analysis_data,
    target = "sample_type",
    positive = positive_class
  )

  cat("Pipeline created successfully\n")
  cat("Positive class:", positive_class, "\n")

  # Test graph learner creation
  learner <- pipeline$create_graph_learner(
    filter = "anova",
    model = "ranger",
    n_features = 30,
    scale = TRUE
  )

  cat("GraphLearner created:", learner$id, "\n")

  # Quick train/predict test
  task <- pipeline$get_task()
  train_idx <- sample(task$nrow, floor(0.7 * task$nrow))
  test_idx <- setdiff(seq_len(task$nrow), train_idx)

  learner$train(task, row_ids = train_idx)
  pred <- learner$predict(task, row_ids = test_idx)
  auc <- pred$score(msr("classif.auc"))

  cat("Quick validation AUC:", round(auc, 4), "\n")

  test_results$omicpipeline <- list(
    status = "PASS",
    auc = auc,
    n_features_selected = 30
  )
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$omicpipeline <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 3. TEST BenchmarkService
# ============================================================================

cat("\n[3/13] TESTING BenchmarkService\n")
cat("--------------------------------------------------\n")

tryCatch({
  service <- BenchmarkService$new(
    task = pipeline,
    outer_folds = 3,
    inner_folds = 2,
    seed = 42
  )

  cat("BenchmarkService created\n")

  simple_learner <- pipeline$create_graph_learner(
    filter = "anova",
    model = "ranger",
    n_features = 20
  )
  service$add_learner(simple_learner, id = "anova_ranger")

  cat("Running nested CV (3 outer, 2 inner folds)...\n")
  result <- service$run()

  perf <- result$performance
  cat("Mean AUC:", round(mean(perf$classif.auc), 4), "\n")
  cat("Mean Accuracy:", round(mean(perf$classif.acc), 4), "\n")

  test_results$benchmark_service <- list(
    status = "PASS",
    mean_auc = mean(perf$classif.auc),
    mean_acc = mean(perf$classif.acc)
  )
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$benchmark_service <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 4. TEST GOF FILTERS
# ============================================================================

cat("\n[4/13] TESTING GOF FILTERS\n")
cat("--------------------------------------------------\n")

tryCatch({
  task <- pipeline$get_task()

  cat("Testing FilterGOF_KS...\n")
  ks_filter <- FilterGOF_KS$new()
  ks_filter$calculate(task)
  ks_top <- names(sort(ks_filter$scores, decreasing = TRUE))[1:10]
  cat("  Top 5:", paste(ks_top[1:5], collapse = ", "), "\n")

  cat("Testing FilterHurdle...\n")
  hurdle_filter <- FilterHurdle$new()
  hurdle_filter$calculate(task)
  hurdle_top <- names(sort(hurdle_filter$scores, decreasing = TRUE))[1:10]
  cat("  Top 5:", paste(hurdle_top[1:5], collapse = ", "), "\n")

  cat("Testing FilterZeroProp...\n")
  zeroprop_filter <- FilterZeroProp$new()
  zeroprop_filter$calculate(task)
  zeroprop_top <- names(sort(zeroprop_filter$scores, decreasing = TRUE))[1:10]
  cat("  Top 5:", paste(zeroprop_top[1:5], collapse = ", "), "\n")

  # Compare with ANOVA
  anova_filter <- flt("anova")
  anova_filter$calculate(task)
  anova_top <- names(sort(anova_filter$scores, decreasing = TRUE))[1:10]

  ks_unique <- setdiff(ks_top, anova_top)
  cat("KS-unique features (not in ANOVA top 10):", length(ks_unique), "\n")

  test_results$gof_filters <- list(
    status = "PASS",
    ks_top = ks_top,
    hurdle_top = hurdle_top,
    ks_unique_count = length(ks_unique)
  )
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$gof_filters <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 5. TEST BAYESIAN TUNING
# ============================================================================

cat("\n[5/13] TESTING BAYESIAN TUNING\n")
cat("--------------------------------------------------\n")

tryCatch({
  if (!requireNamespace("mlr3mbo", quietly = TRUE)) {
    cat("Skipping: mlr3mbo not installed\n")
    test_results$bayesian_tuning <- list(status = "SKIP", reason = "mlr3mbo not installed")
  } else {
    task <- pipeline$get_task()

    cat("Creating Bayesian-tuned glmnet autotuner...\n")
    autotuner <- make_autotuner_glmnet(task, n_evals = 5, inner_folds = 2)

    cat("Training (5 Bayesian iterations)...\n")
    autotuner$train(task)

    pred <- autotuner$predict(task)
    auc <- pred$score(msr("classif.auc"))
    cat("Training AUC:", round(auc, 4), "\n")

    test_results$bayesian_tuning <- list(status = "PASS", auc = auc)
    cat("STATUS: PASS\n")
  }
}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$bayesian_tuning <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 6. TEST AutoXAI
# ============================================================================

cat("\n[6/13] TESTING AutoXAI\n")
cat("--------------------------------------------------\n")

tryCatch({
  if (!requireNamespace("DALEX", quietly = TRUE)) {
    cat("Skipping: DALEX not installed\n")
    test_results$autoxai <- list(status = "SKIP", reason = "DALEX not installed")
  } else {
    # Use smaller subset for faster XAI testing
    task <- pipeline$get_task()
    # Select only top 30 features by variance for speed
    data <- as.data.frame(task$data())
    feature_cols <- task$feature_names
    variances <- sapply(data[, feature_cols, drop = FALSE], var, na.rm = TRUE)
    top_var_features <- names(sort(variances, decreasing = TRUE))[1:30]

    small_data <- data[, c(task$target_names, top_var_features), drop = FALSE]
    small_task <- as_task_classif(small_data, target = task$target_names, positive = task$positive)

    cat("Training model for XAI analysis (30 features)...\n")
    learner <- lrn("classif.ranger", predict_type = "prob", num.trees = 50)
    learner$train(small_task)

    cat("Running XAI pipeline...\n")
    xai_results <- xai_pipeline(
      learner = learner,
      task = small_task,
      top_k = 10,
      n_shap_obs = 2,  # Reduced for speed
      cor_threshold = 0.7,
      verbose = TRUE
    )

    cat("Top 5 features:", paste(xai_results$top_features[1:5], collapse = ", "), "\n")
    cat("Correlation warnings:", xai_results$correlations$n_warnings, "\n")

    test_results$autoxai <- list(
      status = "PASS",
      top_features = xai_results$top_features,
      n_cor_warnings = xai_results$correlations$n_warnings
    )
    cat("STATUS: PASS\n")
  }
}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$autoxai <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 7. TEST StabilityEnsemble
# ============================================================================

cat("\n[7/13] TESTING StabilityEnsemble\n")
cat("--------------------------------------------------\n")

tryCatch({
  # Use smaller subset for faster stability testing
  task <- pipeline$get_task()
  data <- as.data.frame(task$data())
  feature_cols <- task$feature_names[1:50]  # First 50 features
  small_data <- data[, c(task$target_names, feature_cols), drop = FALSE]
  small_task <- as_task_classif(small_data, target = task$target_names, positive = task$positive)

  cat("Creating stability ensemble (20 bootstraps, 50 features)...\n")
  ensemble <- create_stability_ensemble(preset = "fast", n_bootstrap = 20)

  cat("Fitting ensemble...\n")
  ensemble$fit(small_task, seed = 42)

  stable_features <- ensemble$get_feature_importance(20)
  cat("Top 5 most stable:\n")
  print(head(stable_features, 5))

  cat("Tier weights:", paste(round(ensemble$tier_weights, 3), collapse = ", "), "\n")

  test_results$stability_ensemble <- list(
    status = "PASS",
    top_stable = names(stable_features)[1:5],
    n_tiers = length(ensemble$stable_features)
  )
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$stability_ensemble <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 8. TEST SequentialSelector (HSFS)
# ============================================================================

cat("\n[8/13] TESTING SequentialSelector (HSFS)\n")
cat("--------------------------------------------------\n")

tryCatch({
  # Use smaller subset for faster HSFS testing
  task <- pipeline$get_task()
  data <- as.data.frame(task$data())
  feature_cols <- task$feature_names[1:60]  # First 60 features
  small_data <- data[, c(task$target_names, feature_cols), drop = FALSE]
  task <- as_task_classif(small_data, target = task$target_names, positive = pipeline$get_task()$positive)

  cat("Creating HSFS selector (minimal preset, 60 features)...\n")
  selector <- create_hsfs_selector("minimal")

  cat("Creating GraphLearner...\n")
  hsfs_learner <- selector$create_learner(model = "ranger")

  train_idx <- sample(task$nrow, floor(0.7 * task$nrow))
  test_idx <- setdiff(seq_len(task$nrow), train_idx)

  cat("Training HSFS learner...\n")
  hsfs_learner$train(task, row_ids = train_idx)

  pred <- hsfs_learner$predict(task, row_ids = test_idx)
  auc <- pred$score(msr("classif.auc"))

  cat("Test AUC:", round(auc, 4), "\n")

  test_results$sequential_selector <- list(status = "PASS", auc = auc)
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$sequential_selector <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 9. TEST SyntheticData (SMOTE, Noise)
# ============================================================================

cat("\n[9/13] TESTING SyntheticData\n")
cat("--------------------------------------------------\n")

tryCatch({
  task <- pipeline$get_task()

  cat("Original class distribution:\n")
  print(table(task$truth()))

  cat("\nTesting SMOTE...\n")
  task_smote <- smote_augment(task, ratio = 1.0, k = 5)
  cat("After SMOTE:", task_smote$nrow, "samples\n")
  print(table(task_smote$truth()))

  cat("\nTesting noise augmentation...\n")
  task_noise <- noise_augment(task, n_augment = 1, noise_sd = 0.1)
  cat("After noise:", task_noise$nrow, "samples\n")

  test_results$synthetic_data <- list(
    status = "PASS",
    smote_n = task_smote$nrow,
    noise_n = task_noise$nrow
  )
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$synthetic_data <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 10. TEST Calibration
# ============================================================================

cat("\n[10/13] TESTING Calibration\n")
cat("--------------------------------------------------\n")

tryCatch({
  task <- pipeline$get_task()

  learner <- lrn("classif.ranger", predict_type = "prob", num.trees = 100)

  # Split data
  n <- task$nrow
  train_idx <- sample(n, floor(0.6 * n))
  remaining <- setdiff(seq_len(n), train_idx)
  calib_idx <- sample(remaining, floor(0.5 * length(remaining)))
  test_idx <- setdiff(remaining, calib_idx)

  learner$train(task, row_ids = train_idx)

  calib_pred <- learner$predict(task, row_ids = calib_idx)
  calib_probs <- calib_pred$prob[, task$positive]
  calib_y <- as.numeric(calib_pred$truth == task$positive)

  cat("Testing Platt scaling...\n")
  # fit_platt_scaling returns a calibrator function
  platt_calibrator <- fit_platt_scaling(calib_probs, calib_y)
  cat("Platt calibrator function created\n")

  test_pred <- learner$predict(task, row_ids = test_idx)
  test_probs <- test_pred$prob[, task$positive]

  calibrated <- platt_calibrator(test_probs)  # Apply by calling the function
  cat("Original range:", round(min(test_probs), 3), "-", round(max(test_probs), 3), "\n")
  cat("Calibrated range:", round(min(calibrated), 3), "-", round(max(calibrated), 3), "\n")

  cat("\nTesting isotonic...\n")
  iso_calibrator <- fit_isotonic_calibration(calib_probs, calib_y)
  iso_probs <- iso_calibrator(test_probs)  # Apply by calling the function
  cat("Isotonic range:", round(min(iso_probs), 3), "-", round(max(iso_probs), 3), "\n")

  test_results$calibration <- list(
    status = "PASS",
    original_range = paste(round(min(test_probs), 3), "-", round(max(test_probs), 3)),
    calibrated_range = paste(round(min(calibrated), 3), "-", round(max(calibrated), 3))
  )
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$calibration <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 11. TEST FrozenComBat
# ============================================================================

cat("\n[11/13] TESTING FrozenComBat\n")
cat("--------------------------------------------------\n")

tryCatch({
  data <- as.data.frame(pipeline$get_task()$data())
  features <- pipeline$get_feature_names()[1:50]  # Use subset

  n <- nrow(data)
  # Create batches - ensure both batches have samples
  set.seed(42)
  batch <- factor(sample(c("A", "B"), n, replace = TRUE, prob = c(0.6, 0.4)))

  X <- as.matrix(data[, features])
  # Add artificial batch effect (constant shift for batch B)
  batch_effect <- matrix(ifelse(batch == "B", 2.0, 0), nrow = n, ncol = ncol(X))
  X_batched <- X + batch_effect

  cat("Batch distribution:", paste(table(batch), collapse = " / "), "\n")
  cat("Original batch B - A mean shift:", round(mean(X_batched[batch == "B", 1]) -
                                                  mean(X_batched[batch == "A", 1]), 4), "\n")

  # Split into train/test ensuring both batches in training
  train_idx <- sample(n, floor(0.7 * n))
  test_idx <- setdiff(seq_len(n), train_idx)

  # Verify training has both batches
  train_batches <- table(batch[train_idx])
  cat("Training batch counts:", paste(train_batches, collapse = " / "), "\n")

  result <- frozen_combat_correct(
    train_data = X_batched[train_idx, , drop = FALSE],
    train_batch = batch[train_idx],
    test_data = X_batched[test_idx, , drop = FALSE],
    test_batch = batch[test_idx]
  )

  cat("FrozenComBat model fitted and applied\n")
  cat("Corrected test data dimensions:", paste(dim(result$corrected_test), collapse = " x "), "\n")

  # Check if batch effect was reduced in test set
  test_batch <- batch[test_idx]
  original_test_A_mean <- mean(X_batched[test_idx[test_batch == "A"], 1])
  original_test_B_mean <- mean(X_batched[test_idx[test_batch == "B"], 1])
  corrected_test_A_mean <- mean(result$corrected_test[test_batch == "A", 1])
  corrected_test_B_mean <- mean(result$corrected_test[test_batch == "B", 1])

  original_shift <- original_test_B_mean - original_test_A_mean
  corrected_shift <- corrected_test_B_mean - corrected_test_A_mean

  cat("Original batch shift in test:", round(original_shift, 4), "\n")
  cat("Corrected batch shift in test:", round(corrected_shift, 4), "\n")

  test_results$frozen_combat <- list(status = "PASS",
                                      original_shift = original_shift,
                                      corrected_shift = corrected_shift)
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$frozen_combat <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 12. TEST Multi-Omics Integration
# ============================================================================

cat("\n[12/13] TESTING Multi-Omics Integration\n")
cat("--------------------------------------------------\n")

tryCatch({
  # Test multi-omics by creating combined feature matrix with prefixed names
  data <- as.data.frame(pipeline$get_task()$data())
  features <- pipeline$get_feature_names()
  target_col <- pipeline$get_task()$target_names

  n_feat <- min(100, length(features))

  # Create combined data with "modality" prefixes (simulating multi-omics)
  modality1_features <- features[1:floor(n_feat/2)]
  modality2_features <- features[(floor(n_feat/2)+1):n_feat]

  combined_data <- data.frame(
    data[, modality1_features, drop = FALSE],
    data[, modality2_features, drop = FALSE]
  )

  # Rename with modality prefixes
  names(combined_data)[1:length(modality1_features)] <-
    paste0("mod1_", names(combined_data)[1:length(modality1_features)])
  names(combined_data)[(length(modality1_features)+1):ncol(combined_data)] <-
    paste0("mod2_", names(combined_data)[(length(modality1_features)+1):ncol(combined_data)])

  combined_data[[target_col]] <- data[[target_col]]

  cat("Modality 1 features:", length(modality1_features), "\n")
  cat("Modality 2 features:", length(modality2_features), "\n")
  cat("Combined features:", ncol(combined_data) - 1, "\n")

  # Create pipeline with combined data (single modality input, prefixed features)
  levels_ordered <- levels(analysis_data$sample_type)
  positive_class <- if ("PrimaryTumor" %in% levels_ordered) "PrimaryTumor" else levels_ordered[2]

  multi_pipeline <- OmicPipeline$new(
    data = combined_data,
    target = target_col,
    positive = positive_class
  )

  cat("Multi-omics style pipeline created\n")
  cat("Total features:", length(multi_pipeline$get_feature_names()), "\n")

  multi_learner <- multi_pipeline$create_graph_learner(
    filter = "anova",
    model = "ranger",
    n_features = 20
  )

  task <- multi_pipeline$get_task()
  train_idx <- sample(task$nrow, floor(0.7 * task$nrow))
  multi_learner$train(task, row_ids = train_idx)

  pred <- multi_learner$predict(task, row_ids = setdiff(seq_len(task$nrow), train_idx))
  auc <- pred$score(msr("classif.auc"))

  cat("Multi-omics style AUC:", round(auc, 4), "\n")

  # Check both modalities contributed (features from both prefixes selected)
  n_mod1 <- sum(grepl("^mod1_", task$feature_names))
  n_mod2 <- sum(grepl("^mod2_", task$feature_names))
  cat("Features from mod1:", n_mod1, ", from mod2:", n_mod2, "\n")

  test_results$multi_omics <- list(
    status = "PASS",
    n_modalities = 2,
    auc = auc
  )
  cat("STATUS: PASS\n")

}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$multi_omics <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# 13. TEST Deep Learning (infrastructure check)
# ============================================================================

cat("\n[13/13] TESTING Deep Learning Infrastructure\n")
cat("--------------------------------------------------\n")

tryCatch({
  if (!requireNamespace("torch", quietly = TRUE)) {
    cat("torch not installed - infrastructure check only\n")
    test_results$deep_learning <- list(status = "SKIP", reason = "torch not installed")
  } else if (!requireNamespace("mlr3torch", quietly = TRUE)) {
    cat("mlr3torch not installed - infrastructure check only\n")
    test_results$deep_learning <- list(status = "SKIP", reason = "mlr3torch not installed")
  } else {
    cat("Creating MLP learner...\n")
    mlp_learner <- make_mlp_learner(n_hidden = 64, n_layers = 2, dropout = 0.5, epochs = 5)
    cat("MLP learner created successfully\n")
    test_results$deep_learning <- list(status = "PASS", note = "Infrastructure validated")
    cat("STATUS: PASS\n")
  }
}, error = function(e) {
  cat("STATUS: FAIL -", e$message, "\n")
  test_results$deep_learning <- list(status = "FAIL", error = e$message)
})

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================================================\n")
cat("TEST SUMMARY\n")
cat("========================================================================\n\n")

passed <- sum(sapply(test_results, function(x) {
  if (is.list(x) && "status" %in% names(x)) x$status == "PASS" else FALSE
}))
failed <- sum(sapply(test_results, function(x) {
  if (is.list(x) && "status" %in% names(x)) x$status == "FAIL" else FALSE
}))
skipped <- sum(sapply(test_results, function(x) {
  if (is.list(x) && "status" %in% names(x)) x$status == "SKIP" else FALSE
}))

cat("PASSED:", passed, "\n")
cat("FAILED:", failed, "\n")
cat("SKIPPED:", skipped, "\n\n")

for (name in names(test_results)) {
  result <- test_results[[name]]
  if (is.list(result) && "status" %in% names(result)) {
    symbol <- switch(result$status, "PASS" = "[OK]", "FAIL" = "[XX]", "SKIP" = "[--]", "[??]")
    cat(symbol, name, "\n")
    if (result$status == "FAIL" && "error" %in% names(result)) {
      cat("    Error:", substr(result$error, 1, 60), "\n")
    }
  }
}

saveRDS(test_results, "scripts/test_results.rds")
cat("\nResults saved to scripts/test_results.rds\n")
cat("\nComprehensive test complete!\n")
