#!/usr/bin/env Rscript
#' @title Test Bayesian Hyperparameter Tuning
#' @description
#' Tests the BayesianTuner implementation on TCGA kidney cancer data.

cat("=", rep("=", 70), "\n", sep = "")
cat("Bayesian Hyperparameter Tuning Test\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Load packages
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3pipelines)
  library(mlr3learners)
  library(mlr3tuning)
  library(data.table)
})

# Source BayesianTuner
pkg_dir <- "/projekty/OmicSelector2/OmicSelector/R"
source(file.path(pkg_dir, "BayesianTuner.R"))

set.seed(42)

# ============================================================================
# 1. Prepare Data (subset for faster testing)
# ============================================================================

cat("[1/4] Preparing TCGA Kidney Data...\n")
cat("-", rep("-", 50), "\n", sep = "")

load("/projekty/OmicSelector2/OmicSelector/data/orginal_TCGA_data.rda")
df <- as.data.frame(orginal_TCGA_data)

# Get miRNA columns
mirna_cols <- grep("^hsa", names(df), value = TRUE)

# Filter kidney samples
kidney_tumor <- df[df$primary_site == "Kidney" & df$sample_type == "PrimaryTumor", ]
kidney_normal <- df[df$primary_site == "Kidney" & df$sample_type == "SolidTissueNormal", ]

# Subsample for speed (100 tumor, 50 normal)
set.seed(42)
kidney_tumor <- kidney_tumor[sample(nrow(kidney_tumor), 100), ]
kidney_normal <- kidney_normal[sample(nrow(kidney_normal), 50), ]

# Combine
kidney_data <- rbind(
  cbind(kidney_tumor[, mirna_cols, drop = FALSE], outcome = "Tumor"),
  cbind(kidney_normal[, mirna_cols, drop = FALSE], outcome = "Normal")
)
kidney_data$outcome <- factor(kidney_data$outcome, levels = c("Normal", "Tumor"))

# Remove zero variance
var_check <- sapply(kidney_data[, mirna_cols], var, na.rm = TRUE)
zero_var <- names(var_check[var_check == 0 | is.na(var_check)])
if (length(zero_var) > 0) {
  kidney_data <- kidney_data[, !names(kidney_data) %in% zero_var]
  mirna_cols <- setdiff(mirna_cols, zero_var)
}

# Pre-filter to top 100 features for speed
feature_vars <- sapply(kidney_data[, mirna_cols], var, na.rm = TRUE)
top_features <- names(sort(feature_vars, decreasing = TRUE))[1:100]
kidney_data <- kidney_data[, c(top_features, "outcome")]

cat(sprintf("  Samples: %d (Tumor: %d, Normal: %d)\n",
            nrow(kidney_data), sum(kidney_data$outcome == "Tumor"),
            sum(kidney_data$outcome == "Normal")))
cat(sprintf("  Features: %d\n", length(top_features)))

# Create task
task <- as_task_classif(kidney_data, target = "outcome", positive = "Tumor")
task$id <- "TCGA_Kidney_Test"

cat("  Task created successfully\n\n")

# ============================================================================
# 2. Test Package Check Function
# ============================================================================

cat("[2/4] Testing Package Dependencies...\n")
cat("-", rep("-", 50), "\n", sep = "")

tryCatch({
  check_mbo_packages()
  cat("  All required packages available\n\n")
}, error = function(e) {
  cat(sprintf("  Missing packages: %s\n\n", e$message))
})

# ============================================================================
# 3. Test AutoTuner Creation
# ============================================================================

cat("[3/4] Testing AutoTuner Creation...\n")
cat("-", rep("-", 50), "\n", sep = "")

# Test glmnet (fast)
tryCatch({
  at_glmnet <- make_autotuner_glmnet(task, n_evals = 5, inner_folds = 2)
  cat(sprintf("  glmnet AutoTuner: OK (id=%s)\n", at_glmnet$id))
}, error = function(e) {
  cat(sprintf("  glmnet AutoTuner: FAILED - %s\n", e$message))
})

# Test ranger
tryCatch({
  at_ranger <- make_autotuner_ranger(task, n_evals = 5, inner_folds = 2)
  cat(sprintf("  ranger AutoTuner: OK (id=%s)\n", at_ranger$id))
}, error = function(e) {
  cat(sprintf("  ranger AutoTuner: FAILED - %s\n", e$message))
})

# Test xgboost
tryCatch({
  at_xgb <- make_autotuner_xgboost(task, n_evals = 5, inner_folds = 2)
  cat(sprintf("  xgboost AutoTuner: OK (id=%s)\n", at_xgb$id))
}, error = function(e) {
  cat(sprintf("  xgboost AutoTuner: FAILED - %s\n", e$message))
})

cat("\n")

# ============================================================================
# 4. Quick Benchmark Test
# ============================================================================

cat("[4/4] Running Quick Benchmark (glmnet only, 5 evals)...\n")
cat("-", rep("-", 50), "\n", sep = "")

tryCatch({
  # Use only glmnet for speed
  at <- make_autotuner_glmnet(task, n_evals = 5, inner_folds = 2)
  at$id <- "tuned_glmnet"

  # Single train/predict test
  cat("  Training AutoTuner...\n")
  at$train(task)

  cat("  Extracting optimal parameters...\n")
  opt_params <- get_optimal_params(at)
  cat(sprintf("    alpha = %.4f\n", opt_params$alpha))
  cat(sprintf("    s = %.6f\n", opt_params$s))

  # Quick CV
  cat("\n  Running 3-fold CV...\n")
  resampling <- rsmp("cv", folds = 3)
  rr <- resample(task, at, resampling, store_models = TRUE)

  auc <- rr$aggregate(msr("classif.auc"))
  acc <- rr$aggregate(msr("classif.acc"))

  cat(sprintf("\n  Results:\n"))
  cat(sprintf("    AUC: %.4f\n", auc))
  cat(sprintf("    Accuracy: %.4f\n", acc))

}, error = function(e) {
  cat(sprintf("\n  ERROR: %s\n", e$message))
  cat(sprintf("  Traceback: %s\n", paste(capture.output(traceback()), collapse = "\n")))
})

# ============================================================================
# Summary
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("Bayesian Tuner Test Complete\n")
cat(rep("=", 70), "\n", sep = "")
