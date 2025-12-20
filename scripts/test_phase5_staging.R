#!/usr/bin/env Rscript
#' @title Test Phase 5 Features on Clinical Staging (Harder Task)
#' @description
#' Tests Phase 5 on Early (I/II) vs Late (III/IV) stage classification
#' This is a clinically relevant and harder task than tumor vs normal.

cat("=" , rep("=", 70), "\n", sep = "")
cat("Phase 5 Clinical Staging Test: Early vs Late Stage Kidney Cancer\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

# Load packages
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3pipelines)
  library(mlr3learners)
  library(mlr3filters)
  library(data.table)
  library(R6)
})

# Source Phase 5 modules
pkg_dir <- "/projekty/OmicSelector2/OmicSelector/R"
source(file.path(pkg_dir, "SequentialSelector.R"))
source(file.path(pkg_dir, "StabilityEnsemble.R"))

set.seed(42)

# ============================================================================
# 1. Prepare Clinical Staging Data
# ============================================================================

cat("[1/5] Preparing Clinical Staging Data...\n")
cat("-" , rep("-", 50), "\n", sep = "")

load("/projekty/OmicSelector2/OmicSelector/data/orginal_TCGA_data.rda")
df <- as.data.frame(orginal_TCGA_data)

# Get miRNA columns
mirna_cols <- grep("^hsa", names(df), value = TRUE)

# Filter kidney tumors with stage info
kidney <- df[df$primary_site == "Kidney" &
             df$sample_type == "PrimaryTumor" &
             !is.na(df$tumor_stage), ]

cat(sprintf("  Kidney tumor samples with staging: %d\n", nrow(kidney)))

# Create binary staging outcome
kidney$stage_binary <- ifelse(
  kidney$tumor_stage %in% c("stagei", "stageii"),
  "Early",
  "Late"
)
kidney$stage_binary <- factor(kidney$stage_binary, levels = c("Early", "Late"))

cat(sprintf("  Early stage (I/II): %d\n", sum(kidney$stage_binary == "Early")))
cat(sprintf("  Late stage (III/IV): %d\n", sum(kidney$stage_binary == "Late")))

# Create dataset
staging_data <- kidney[, c(mirna_cols, "stage_binary")]
names(staging_data)[ncol(staging_data)] <- "outcome"

# Remove zero-variance features
var_check <- sapply(staging_data[, mirna_cols], var, na.rm = TRUE)
zero_var <- names(var_check[var_check == 0 | is.na(var_check)])
if (length(zero_var) > 0) {
  staging_data <- staging_data[, !names(staging_data) %in% zero_var]
  mirna_cols <- setdiff(mirna_cols, zero_var)
}

cat(sprintf("  Features after filtering: %d\n", length(mirna_cols)))

# Create mlr3 task
task <- as_task_classif(staging_data, target = "outcome", positive = "Late")
task$id <- "TCGA_Kidney_Staging"

cat("  Task created successfully\n\n")

# ============================================================================
# 2. Baseline: Simple ANOVA + Random Forest
# ============================================================================

cat("[2/5] Baseline (ANOVA + Random Forest)...\n")
cat("-" , rep("-", 50), "\n", sep = "")

baseline_graph <- po("filter",
                     filter = flt("anova"),
                     filter.nfeat = 50) %>>%
  po("learner", lrn("classif.ranger", predict_type = "prob", num.trees = 200))

baseline_learner <- as_learner(baseline_graph)
baseline_learner$id <- "baseline_anova_rf"

resampling <- rsmp("cv", folds = 5)
rr_baseline <- resample(task, baseline_learner, resampling, store_models = TRUE)

baseline_auc <- rr_baseline$aggregate(msr("classif.auc"))
baseline_acc <- rr_baseline$aggregate(msr("classif.acc"))

cat(sprintf("  Baseline Results (5-fold CV):\n"))
cat(sprintf("    - AUC: %.4f\n", baseline_auc))
cat(sprintf("    - Accuracy: %.4f\n", baseline_acc))

# ============================================================================
# 3. HSFS (Hybrid Sequential Feature Selection)
# ============================================================================

cat("\n[3/5] HSFS (Variance → ANOVA → RF-importance → glmnet)...\n")
cat("-" , rep("-", 50), "\n", sep = "")

tryCatch({
  hsfs <- SequentialSelector$new(
    variance_threshold = 0.01,
    univariate_n = 200,
    univariate_method = "anova",
    rfe_n = 50,
    lasso_n = NULL,
    verbose = FALSE
  )
  hsfs$stages$lasso$enabled <- FALSE  # Skip LASSO for now

  hsfs_learner <- hsfs$create_learner(model = "ranger", scale = TRUE)

  rr_hsfs <- resample(task, hsfs_learner, resampling, store_models = TRUE)

  hsfs_auc <- rr_hsfs$aggregate(msr("classif.auc"))
  hsfs_acc <- rr_hsfs$aggregate(msr("classif.acc"))

  cat(sprintf("  HSFS Results (5-fold CV):\n"))
  cat(sprintf("    - AUC: %.4f\n", hsfs_auc))
  cat(sprintf("    - Accuracy: %.4f\n", hsfs_acc))

}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
  hsfs_auc <<- NA
  hsfs_acc <<- NA
})

# ============================================================================
# 4. Stability Ensemble
# ============================================================================

cat("\n[4/5] StabilityEnsemble (Bootstrap-based selection)...\n")
cat("-" , rep("-", 50), "\n", sep = "")

tryCatch({
  ensemble <- create_stability_ensemble(
    preset = "default",
    n_bootstrap = 100,
    n_features = 30,
    verbose = FALSE
  )

  cat("  Running 100 bootstrap resamples...\n")
  ensemble$fit(task, seed = 42)

  # Get stability results
  top_features <- ensemble$get_feature_importance(20)

  cat(sprintf("\n  Features with >50%% stability: %d\n",
              sum(ensemble$frequencies >= 0.5)))
  cat(sprintf("  Features with >70%% stability: %d\n",
              sum(ensemble$frequencies >= 0.7)))
  cat(sprintf("  Features with >90%% stability: %d\n",
              sum(ensemble$frequencies >= 0.9)))

  cat("\n  Top 10 Most Stable Features for Staging:\n")
  for (i in 1:min(10, length(top_features))) {
    cat(sprintf("    %2d. %s (%.1f%%)\n",
                i, names(top_features)[i], top_features[i] * 100))
  }

  stability_results <- list(
    n_50 = sum(ensemble$frequencies >= 0.5),
    n_70 = sum(ensemble$frequencies >= 0.7),
    n_90 = sum(ensemble$frequencies >= 0.9),
    top_features = top_features
  )

}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
  stability_results <<- list(error = e$message)
})

# ============================================================================
# 5. Summary
# ============================================================================

cat("\n[5/5] Summary\n")
cat("=" , rep("=", 70), "\n", sep = "")

cat("\nCLINICAL STAGING CLASSIFICATION (Early vs Late)\n")
cat("-" , rep("-", 50), "\n", sep = "")
cat(sprintf("  Samples: %d (Early: %d, Late: %d)\n",
            nrow(staging_data),
            sum(staging_data$outcome == "Early"),
            sum(staging_data$outcome == "Late")))

cat("\nPERFORMANCE COMPARISON (5-fold CV):\n")
cat(sprintf("  Baseline (ANOVA+RF):  AUC = %.4f, Acc = %.4f\n",
            baseline_auc, baseline_acc))
if (exists("hsfs_auc") && !is.na(hsfs_auc)) {
  cat(sprintf("  HSFS (Sequential):    AUC = %.4f, Acc = %.4f\n",
              hsfs_auc, hsfs_acc))
}

cat("\nINTERPRETATION:\n")
if (baseline_auc < 0.75) {
  cat("  - AUC < 0.75: This is a challenging task (expected for staging)\n")
} else if (baseline_auc < 0.85) {
  cat("  - AUC 0.75-0.85: Moderate predictive signal detected\n")
} else {
  cat("  - AUC > 0.85: Strong predictive signal for staging\n")
}

cat("\n" , rep("=", 70), "\n", sep = "")
cat("Clinical staging test complete.\n")
