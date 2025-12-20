#!/usr/bin/env Rscript
#' @title Test Auto XAI Module
#' @description Tests the AutoXAI implementation on TCGA kidney cancer data.

cat("=", rep("=", 70), "\n", sep = "")
cat("Auto XAI Test\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Load packages
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3learners)
  library(data.table)
  library(ranger)
  library(DALEX)
  library(ggplot2)
})

# Source AutoXAI
pkg_dir <- "/projekty/OmicSelector2/OmicSelector/R"
source(file.path(pkg_dir, "AutoXAI.R"))

set.seed(42)

# ============================================================================
# 1. Prepare Data
# ============================================================================

cat("[1/5] Preparing TCGA Kidney Data...\n")
cat("-", rep("-", 50), "\n", sep = "")

load("/projekty/OmicSelector2/OmicSelector/data/orginal_TCGA_data.rda")
df <- as.data.frame(orginal_TCGA_data)
mirna_cols <- grep("^hsa", names(df), value = TRUE)

# Filter kidney
kidney_tumor <- df[df$primary_site == "Kidney" & df$sample_type == "PrimaryTumor", ]
kidney_normal <- df[df$primary_site == "Kidney" & df$sample_type == "SolidTissueNormal", ]

# Subsample for speed
set.seed(42)
kidney_tumor <- kidney_tumor[sample(nrow(kidney_tumor), 60), ]
kidney_normal <- kidney_normal[sample(nrow(kidney_normal), 30), ]

kidney_data <- rbind(
  cbind(kidney_tumor[, mirna_cols, drop = FALSE], outcome = "Tumor"),
  cbind(kidney_normal[, mirna_cols, drop = FALSE], outcome = "Normal")
)
kidney_data$outcome <- factor(kidney_data$outcome, levels = c("Normal", "Tumor"))

# Remove zero variance and limit to top 50 features
var_check <- sapply(kidney_data[, mirna_cols], var, na.rm = TRUE)
zero_var <- names(var_check[var_check == 0 | is.na(var_check)])
kidney_data <- kidney_data[, !names(kidney_data) %in% zero_var]
mirna_cols <- setdiff(mirna_cols, zero_var)

# Top 50 by variance
feature_vars <- sapply(kidney_data[, mirna_cols], var, na.rm = TRUE)
top_features <- names(sort(feature_vars, decreasing = TRUE))[1:50]
kidney_data <- kidney_data[, c(top_features, "outcome")]

cat(sprintf("  Samples: %d, Features: %d\n\n", nrow(kidney_data), length(top_features)))

# Create task
task <- as_task_classif(kidney_data, target = "outcome", positive = "Tumor")
task$id <- "TCGA_Kidney"

# ============================================================================
# 2. Train Model
# ============================================================================

cat("[2/5] Training Random Forest model...\n")
cat("-", rep("-", 50), "\n", sep = "")

learner <- lrn("classif.ranger",
               predict_type = "prob",
               num.trees = 200,
               importance = "impurity")
learner$train(task)

cat("  Model trained successfully\n\n")

# ============================================================================
# 3. Test Explainer Creation
# ============================================================================

cat("[3/5] Testing DALEX explainer creation...\n")
cat("-", rep("-", 50), "\n", sep = "")

tryCatch({
  explainer <- xai_explainer_mlr3(learner, task)
  cat(sprintf("  Explainer created: %s\n", explainer$label))
  cat(sprintf("  Data dimensions: %d x %d\n", nrow(explainer$data), ncol(explainer$data)))
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
})

# ============================================================================
# 4. Test Correlation Diagnostics
# ============================================================================

cat("\n[4/5] Testing correlation diagnostics...\n")
cat("-", rep("-", 50), "\n", sep = "")

tryCatch({
  X <- as.data.frame(task$data())[, top_features]
  cor_results <- xai_correlations(X, threshold = 0.7)

  cat(sprintf("  High correlation pairs (|r| > 0.7): %d\n", cor_results$n_warnings))
  cat(sprintf("  Correlation clusters: %d\n", length(cor_results$clusters)))

  if (nrow(cor_results$high_cor_pairs) > 0) {
    cat("\n  Top 5 correlated pairs:\n")
    for (i in 1:min(5, nrow(cor_results$high_cor_pairs))) {
      cat(sprintf("    %s <-> %s (r = %.3f)\n",
                  cor_results$high_cor_pairs$feature1[i],
                  cor_results$high_cor_pairs$feature2[i],
                  cor_results$high_cor_pairs$correlation[i]))
    }
  }
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
})

# ============================================================================
# 5. Test Full XAI Pipeline
# ============================================================================

cat("\n[5/5] Testing full XAI pipeline...\n")
cat("-", rep("-", 50), "\n", sep = "")

tryCatch({
  xai_results <- xai_pipeline(
    learner = learner,
    task = task,
    top_k = 10,
    n_shap_obs = 3,
    cor_threshold = 0.7,
    verbose = TRUE
  )

  # Print summary
  print_xai_summary(xai_results)

  # Test importance plot
  cat("\nGenerating importance plot...\n")
  p <- plot_xai_importance(xai_results, top_n = 10)
  ggsave("/projekty/OmicSelector2/OmicSelector/scripts/xai_importance_test.png",
         p, width = 8, height = 6)
  cat("  Saved to: scripts/xai_importance_test.png\n")

}, error = function(e) {
  cat(sprintf("\n  ERROR: %s\n", e$message))
  traceback()
})

cat("\n", rep("=", 70), "\n", sep = "")
cat("Auto XAI Test Complete\n")
cat(rep("=", 70), "\n", sep = "")
