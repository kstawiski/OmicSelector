#!/usr/bin/env Rscript
#' @title Test GOF Filters
#' @description Tests GOF filters on TCGA kidney cancer data.

cat("=", rep("=", 70), "\n", sep = "")
cat("GOF Filters Test\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Load packages
suppressPackageStartupMessages({
  library(mlr3)
  library(mlr3filters)
  library(data.table)
  library(R6)
})

# Source GOF filters
pkg_dir <- "/projekty/OmicSelector2/OmicSelector/R"
source(file.path(pkg_dir, "FilterGOF.R"))

set.seed(42)

# ============================================================================
# 1. Prepare Data
# ============================================================================

cat("[1/4] Preparing TCGA Kidney Data...\n")
cat("-", rep("-", 50), "\n", sep = "")

load("/projekty/OmicSelector2/OmicSelector/data/orginal_TCGA_data.rda")
df <- as.data.frame(orginal_TCGA_data)
mirna_cols <- grep("^hsa", names(df), value = TRUE)

# Filter kidney
kidney_tumor <- df[df$primary_site == "Kidney" & df$sample_type == "PrimaryTumor", ]
kidney_normal <- df[df$primary_site == "Kidney" & df$sample_type == "SolidTissueNormal", ]

# Subsample for speed
set.seed(42)
kidney_tumor <- kidney_tumor[sample(nrow(kidney_tumor), 100), ]
kidney_normal <- kidney_normal[sample(nrow(kidney_normal), 50), ]

kidney_data <- rbind(
  cbind(kidney_tumor[, mirna_cols, drop = FALSE], outcome = "Tumor"),
  cbind(kidney_normal[, mirna_cols, drop = FALSE], outcome = "Normal")
)
kidney_data$outcome <- factor(kidney_data$outcome, levels = c("Normal", "Tumor"))

# Keep all features (including zero-heavy ones)
cat(sprintf("  Samples: %d, Features: %d\n\n", nrow(kidney_data), length(mirna_cols)))

# Create task
task <- as_task_classif(kidney_data, target = "outcome", positive = "Tumor")
task$id <- "TCGA_Kidney"

# ============================================================================
# 2. Test Individual Filters
# ============================================================================

cat("[2/4] Testing Individual Filters...\n")
cat("-", rep("-", 50), "\n", sep = "")

# Test KS filter
cat("\n  FilterGOF_KS:\n")
tryCatch({
  ks_filter <- FilterGOF_KS$new()
  ks_filter$calculate(task)
  ks_scores <- ks_filter$scores
  ks_top <- sort(ks_scores, decreasing = TRUE)[1:10]

  for (i in 1:10) {
    cat(sprintf("    %2d. %s (D = %.4f)\n", i, names(ks_top)[i], ks_top[i]))
  }
}, error = function(e) {
  cat(sprintf("    ERROR: %s\n", e$message))
})

# Test Hurdle filter
cat("\n  FilterHurdle:\n")
tryCatch({
  hurdle_filter <- FilterHurdle$new()
  hurdle_filter$calculate(task)
  hurdle_scores <- hurdle_filter$scores
  hurdle_top <- sort(hurdle_scores, decreasing = TRUE)[1:10]

  for (i in 1:10) {
    cat(sprintf("    %2d. %s (score = %.2f)\n", i, names(hurdle_top)[i], hurdle_top[i]))
  }
}, error = function(e) {
  cat(sprintf("    ERROR: %s\n", e$message))
})

# Test ZeroProp filter
cat("\n  FilterZeroProp:\n")
tryCatch({
  zp_filter <- FilterZeroProp$new()
  zp_filter$calculate(task)
  zp_scores <- zp_filter$scores
  zp_top <- sort(zp_scores, decreasing = TRUE)[1:10]

  for (i in 1:10) {
    cat(sprintf("    %2d. %s (diff = %.4f)\n", i, names(zp_top)[i], zp_top[i]))
  }
}, error = function(e) {
  cat(sprintf("    ERROR: %s\n", e$message))
})

# ============================================================================
# 3. Compare with ANOVA
# ============================================================================

cat("\n[3/4] Comparing with ANOVA...\n")
cat("-", rep("-", 50), "\n", sep = "")

tryCatch({
  anova_filter <- flt("anova")
  anova_filter$calculate(task)
  anova_scores <- anova_filter$scores
  anova_top <- sort(anova_scores, decreasing = TRUE)[1:10]

  cat("  ANOVA Top 10:\n")
  for (i in 1:10) {
    cat(sprintf("    %2d. %s (F = %.2f)\n", i, names(anova_top)[i], anova_top[i]))
  }

  # Compare overlap
  ks_top10 <- names(sort(ks_scores, decreasing = TRUE))[1:10]
  hurdle_top10 <- names(sort(hurdle_scores, decreasing = TRUE))[1:10]
  anova_top10 <- names(anova_top)

  cat("\n  Overlap Analysis:\n")
  cat(sprintf("    KS vs ANOVA overlap: %d/10\n", length(intersect(ks_top10, anova_top10))))
  cat(sprintf("    Hurdle vs ANOVA overlap: %d/10\n", length(intersect(hurdle_top10, anova_top10))))
  cat(sprintf("    KS vs Hurdle overlap: %d/10\n", length(intersect(ks_top10, hurdle_top10))))

}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
})

# ============================================================================
# 4. Test with mlr3pipelines
# ============================================================================

cat("\n[4/4] Testing Pipeline Integration...\n")
cat("-", rep("-", 50), "\n", sep = "")

tryCatch({
  library(mlr3pipelines)
  library(mlr3learners)

  # Use make_gof_filter convenience function
  ks_filt <- make_gof_filter("ks")

  # Create pipeline with KS filter
  graph <- po("filter", filter = ks_filt, filter.nfeat = 20) %>>%
    po("learner", lrn("classif.ranger", predict_type = "prob", num.trees = 100))

  graph_learner <- as_learner(graph)
  graph_learner$id <- "gof_ks_ranger"

  # Quick CV
  resampling <- rsmp("cv", folds = 3)
  rr <- resample(task, graph_learner, resampling)

  auc <- rr$aggregate(msr("classif.auc"))
  cat(sprintf("  KS Filter + RF Results (3-fold CV):\n"))
  cat(sprintf("    AUC: %.4f\n", auc))

  # Compare with ANOVA baseline
  graph_anova <- po("filter", filter = flt("anova"), filter.nfeat = 20) %>>%
    po("learner", lrn("classif.ranger", predict_type = "prob", num.trees = 100))

  gl_anova <- as_learner(graph_anova)
  gl_anova$id <- "anova_ranger"

  rr_anova <- resample(task, gl_anova, resampling)
  auc_anova <- rr_anova$aggregate(msr("classif.auc"))

  cat(sprintf("\n  ANOVA + RF Baseline:\n"))
  cat(sprintf("    AUC: %.4f\n", auc_anova))

  cat(sprintf("\n  Difference: %.4f\n", auc - auc_anova))

}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message))
})

# ============================================================================
# Summary
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("GOF Filters Test Complete\n")
cat(rep("=", 70), "\n", sep = "")
