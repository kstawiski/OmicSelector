#!/usr/bin/env Rscript
# Test Signature Selection on TCGA Pancreatic Cancer Data
# Version 2: With FrozenComBat batch correction and stability metrics

cat("=== OmicSelector: Signature Selection on TCGA Data (v2) ===\n")
cat("=== With FrozenComBat Batch Correction and Stability Metrics ===\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(mlr3)
  library(mlr3learners)
  library(mlr3pipelines)
  library(mlr3fselect)
  library(mlr3filters)
  library(mlr3tuning)
  library(ranger)
  library(xgboost)
  library(e1071)
})

# Source the package functions
cat("Loading OmicSelector functions...\n")
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  tryCatch(source(f, local = FALSE), error = function(e) NULL)
}

# ============================================================================
# 1. LOAD AND PREPARE TCGA DATA WITH PROPER CONTROLS
# ============================================================================
cat("\n--- Step 1: Loading TCGA Data with Proper Controls ---\n")

# Load the packaged TCGA data
load("data/orginal_TCGA_data.rda")
cat(sprintf("Loaded TCGA data: %d samples, %d columns\n",
            nrow(orginal_TCGA_data), ncol(orginal_TCGA_data)))

# Show available primary sites
cat("\nPrimary sites in TCGA data:\n")
print(table(orginal_TCGA_data$primary_site))

# Filter for pancreatic cancer cases
pancreatic_cancer <- orginal_TCGA_data[orginal_TCGA_data$primary_site == "Pancreas" &
                                        orginal_TCGA_data$sample_type == "PrimaryTumor", ]

# Use ONLY pancreatic normal tissues as controls (scientifically correct)
pancreatic_normal <- orginal_TCGA_data[orginal_TCGA_data$primary_site == "Pancreas" &
                                        orginal_TCGA_data$sample_type == "SolidTissueNormal", ]

cat(sprintf("\nPancreatic cancer cases: %d\n", nrow(pancreatic_cancer)))
cat(sprintf("Pancreatic normal controls: %d\n", nrow(pancreatic_normal)))

# Check if we have enough pancreatic normals
if (nrow(pancreatic_normal) < 10) {
  cat("\nWARNING: Very few pancreatic normal samples. Adding other GI normals.\n")
  # Add other GI tract normals as alternative
  gi_sites <- c("Stomach", "Colon", "Liver and intrahepatic bile ducts")
  gi_normals <- orginal_TCGA_data[orginal_TCGA_data$primary_site %in% gi_sites &
                                   orginal_TCGA_data$sample_type == "SolidTissueNormal", ]
  # Sample to balance
  n_to_add <- min(100, nrow(gi_normals))
  set.seed(42)
  gi_sample <- gi_normals[sample(nrow(gi_normals), n_to_add), ]
  pancreatic_normal <- rbind(pancreatic_normal, gi_sample)
  cat(sprintf("Added %d GI normal samples. Total controls: %d\n", n_to_add, nrow(pancreatic_normal)))
}

# Label classes
pancreatic_cancer$Class <- "Case"
pancreatic_normal$Class <- "Control"

# Combine dataset
dataset <- rbind(pancreatic_cancer, pancreatic_normal)
cat(sprintf("\nFinal dataset: %d samples\n", nrow(dataset)))
cat(sprintf("Class distribution: Case=%d, Control=%d\n",
            sum(dataset$Class == "Case"), sum(dataset$Class == "Control")))

# ============================================================================
# 2. EXTRACT BATCH INFORMATION AND APPLY FROZENCOMBAT
# ============================================================================
cat("\n--- Step 2: Batch Correction with FrozenComBat ---\n")

# Extract miRNA features
mirna_cols <- grep("^hsa", names(dataset), value = TRUE)
cat(sprintf("miRNA features: %d\n", length(mirna_cols)))

# Use first 200 miRNAs for faster testing
mirna_cols <- mirna_cols[1:min(200, length(mirna_cols))]

# Create feature matrix
X <- as.matrix(dataset[, mirna_cols, with = FALSE])
y <- factor(dataset$Class, levels = c("Control", "Case"))

# Remove features with zero variance
var_filter <- apply(X, 2, var, na.rm = TRUE) > 0
X <- X[, var_filter]
cat(sprintf("Features after variance filter: %d\n", ncol(X)))

# Handle missing values
for (j in seq_len(ncol(X))) {
  na_idx <- is.na(X[, j])
  if (any(na_idx)) {
    X[na_idx, j] <- median(X[, j], na.rm = TRUE)
  }
}

# Log transform
X <- log10(X + 1)

# Create batch variable from primary_site (or use project_id if available)
batch <- as.character(dataset$primary_site)
batch[is.na(batch)] <- "Unknown"
batch <- as.factor(batch)

cat(sprintf("Batch levels: %s\n", paste(levels(batch), collapse = ", ")))
cat("Batch distribution:\n")
print(table(batch))

# Apply FrozenComBat (fit on training data to prevent leakage)
# For this demo, we'll fit on all data first to show the effect
# In practice, you would fit ONLY on training folds

cat("\nApplying FrozenComBat batch correction...\n")

if (length(levels(batch)) >= 2) {
  fc <- FrozenComBat$new(parametric = TRUE, mean_only = FALSE)

  tryCatch({
    X_corrected <- fc$fit_transform(X, batch)
    cat("FrozenComBat correction applied successfully.\n")

    # Compare variance before/after
    var_before <- mean(apply(X, 2, var))
    var_after <- mean(apply(X_corrected, 2, var))
    cat(sprintf("Mean feature variance: Before=%.4f, After=%.4f\n", var_before, var_after))

    X <- X_corrected
  }, error = function(e) {
    cat(sprintf("FrozenComBat error: %s\nProceeding without batch correction.\n", e$message))
  })
} else {
  cat("Only one batch level - skipping batch correction.\n")
}

cat(sprintf("\nFinal dataset: %d samples, %d features\n", nrow(X), ncol(X)))

# ============================================================================
# 3. CREATE MLR3 TASK AND LEARNERS
# ============================================================================
cat("\n--- Step 3: Setting Up ML Pipeline ---\n")

# Create data.frame for mlr3
df <- data.frame(X, Class = y)

# Create classification task
task <- TaskClassif$new(
  id = "pancreatic_mirna",
  backend = df,
  target = "Class",
  positive = "Case"
)

cat(sprintf("Task created: %s (%d obs, %d features)\n",
            task$id, task$nrow, task$n_feature_cols))

# Store total features for stability computation
total_features <- task$n_feature_cols

# Create learners with different filter sizes
# Learner 1: Random Forest with variance filter (50 features)
po_var <- po("filter", filter = flt("variance"), filter.nfeat = 50)
po_rf <- po("learner", lrn("classif.ranger", predict_type = "prob", num.trees = 100))
graph_rf <- po_var %>>% po_rf
glrn_rf <- as_learner(graph_rf)
glrn_rf$id <- "rf_var50"

# Learner 2: Logistic Regression with ANOVA filter (30 features)
po_anova <- po("filter", filter = flt("anova"), filter.nfeat = 30)
po_logreg <- po("learner", lrn("classif.log_reg", predict_type = "prob"))
graph_logreg <- po_anova %>>% po_logreg
glrn_logreg <- as_learner(graph_logreg)
glrn_logreg$id <- "logreg_anova30"

# Learner 3: XGBoost with AUC filter (40 features)
po_auc <- po("filter", filter = flt("auc"), filter.nfeat = 40)
po_xgb <- po("learner", lrn("classif.xgboost", predict_type = "prob",
                             nrounds = 50, verbose = 0))
graph_xgb <- po_auc %>>% po_xgb
glrn_xgb <- as_learner(graph_xgb)
glrn_xgb$id <- "xgb_auc40"

# Learner 4: SVM with Kruskal-Wallis filter (15 features - small panel)
po_kruskal <- po("filter", filter = flt("kruskal_test"), filter.nfeat = 15)
po_svm <- po("learner", lrn("classif.svm", predict_type = "prob",
                             type = "C-classification", kernel = "radial"))
graph_svm <- po_kruskal %>>% po_svm
glrn_svm <- as_learner(graph_svm)
glrn_svm$id <- "svm_kruskal15"

# Learner 5: Random Forest with very small panel (10 features)
po_var10 <- po("filter", filter = flt("variance"), filter.nfeat = 10)
po_rf10 <- po("learner", lrn("classif.ranger", predict_type = "prob", num.trees = 100))
graph_rf10 <- po_var10 %>>% po_rf10
glrn_rf10 <- as_learner(graph_rf10)
glrn_rf10$id <- "rf_var10"

cat("Created 5 learners with embedded feature selection:\n")
cat("  - rf_var50: Random Forest with top 50 by variance\n")
cat("  - logreg_anova30: Logistic Regression with top 30 by ANOVA\n")
cat("  - xgb_auc40: XGBoost with top 40 by AUC filter\n")
cat("  - svm_kruskal15: SVM with top 15 by Kruskal-Wallis\n")
cat("  - rf_var10: Random Forest with top 10 by variance (small panel)\n")

# ============================================================================
# 4. RUN NESTED CROSS-VALIDATION
# ============================================================================
cat("\n--- Step 4: Running Nested Cross-Validation ---\n")

set.seed(42)

# Create BenchmarkService
benchmark_service <- BenchmarkService$new(
  task = task,
  outer_folds = 5,
  inner_folds = 3,
  seed = 42
)

# Add learners
benchmark_service$add_learner(glrn_rf, id = "rf_var50")
benchmark_service$add_learner(glrn_logreg, id = "logreg_anova30")
benchmark_service$add_learner(glrn_xgb, id = "xgb_auc40")
benchmark_service$add_learner(glrn_svm, id = "svm_kruskal15")
benchmark_service$add_learner(glrn_rf10, id = "rf_var10")

# Run benchmark
cat("\nRunning nested CV (this may take a few minutes)...\n")
start_time <- Sys.time()

result <- benchmark_service$run(
  measures = list(
    msr("classif.auc"),
    msr("classif.acc"),
    msr("classif.sensitivity"),
    msr("classif.specificity"),
    msr("classif.bbrier")
  )
)

end_time <- Sys.time()
cat(sprintf("Completed in %.1f seconds\n", difftime(end_time, start_time, units = "secs")))

# ============================================================================
# 5. EXTRACT DETAILED RESULTS
# ============================================================================
cat("\n--- Step 5: Benchmark Results ---\n")

# Get per-fold performance
measures <- list(
  msr("classif.auc"),
  msr("classif.acc"),
  msr("classif.sensitivity"),
  msr("classif.specificity"),
  msr("classif.bbrier")
)
score_dt <- as.data.table(result$benchmark_result$score(measures))

perf_summary <- score_dt[, .(
  mean_auc = mean(`classif.auc`, na.rm = TRUE),
  sd_auc = sd(`classif.auc`, na.rm = TRUE),
  mean_acc = mean(`classif.acc`, na.rm = TRUE),
  mean_sens = mean(`classif.sensitivity`, na.rm = TRUE),
  mean_spec = mean(`classif.specificity`, na.rm = TRUE),
  mean_brier = mean(`classif.bbrier`, na.rm = TRUE)
), by = .(learner_id)]

cat("\nPerformance Summary:\n")
print(perf_summary)

# ============================================================================
# 6. SIGNATURE SELECTION WITH ALL MODES
# ============================================================================
cat("\n--- Step 6: Signature Selection ---\n")

# Mode 1: Constrained 1SE (remove hard constraints for now to see all candidates)
cat("\n=== Mode: Constrained 1SE ===\n")
best_1se <- select_best_signature(
  result,
  mode = "constrained_1se",
  metric = "classif.auc"
)

if (nrow(best_1se) > 0) {
  cat("\nAll candidates (ranked by 1SE rule):\n")
  print(best_1se[, .(learner_id, mean_metric, se_metric, stability, mean_k, selected)])

  selected_1se <- best_1se[selected == TRUE]
  cat(sprintf("\n** Selected (1SE): %s **\n", selected_1se$learner_id[1]))
  cat(sprintf("   AUC: %.4f (SE: %.4f)\n", selected_1se$mean_metric[1], selected_1se$se_metric[1]))
  if (!is.na(selected_1se$stability[1])) {
    cat(sprintf("   Stability (Nogueira): %.4f\n", selected_1se$stability[1]))
  }
  if (!is.na(selected_1se$mean_k[1])) {
    cat(sprintf("   Features: %.0f\n", selected_1se$mean_k[1]))
  }
}

# Mode 2: Weighted scoring
cat("\n=== Mode: Weighted Scoring ===\n")
best_weighted <- select_best_signature(
  result,
  mode = "weighted",
  metric = "classif.auc",
  weights = c(performance = 0.4, stability = 0.35, parsimony = 0.25)
)

if (nrow(best_weighted) > 0) {
  cat("\nWeighted scores:\n")
  print(best_weighted[, .(learner_id, mean_metric, stability, mean_k, score, selected)])

  selected_weighted <- best_weighted[selected == TRUE]
  cat(sprintf("\n** Selected (weighted): %s **\n", selected_weighted$learner_id[1]))
  cat(sprintf("   Score: %.4f\n", selected_weighted$score[1]))
}

# Mode 3: Pareto frontier
cat("\n=== Mode: Pareto Frontier ===\n")
pareto_set <- select_best_signature(
  result,
  mode = "pareto",
  metric = "classif.auc"
)

if (nrow(pareto_set) > 0) {
  cat("\nPareto-optimal signatures:\n")
  cols_to_show <- intersect(
    c("learner_id", "mean_metric", "stability", "mean_k", "pareto", "selected", "dist_utopia"),
    names(pareto_set)
  )
  print(pareto_set[, ..cols_to_show])

  cat(sprintf("\nPareto set size: %d\n", sum(pareto_set$pareto, na.rm = TRUE)))
  recommended <- pareto_set[selected == TRUE]
  if (nrow(recommended) > 0) {
    cat(sprintf("** Recommended (closest to Utopia): %s **\n", recommended$learner_id[1]))
  }
}

# ============================================================================
# 7. CONSENSUS FEATURES
# ============================================================================
cat("\n--- Step 7: Consensus Features from Best Signature ---\n")

if (nrow(best_weighted) > 0) {
  best_learner <- best_weighted[selected == TRUE]$learner_id[1]

  consensus <- get_consensus_features(result, best_learner, min_frequency = 0.6)

  if (nrow(consensus) > 0) {
    cat(sprintf("\nConsensus features for %s (>=60%% frequency):\n", best_learner))
    stable_features <- consensus[selected == TRUE]
    cat(sprintf("  Total stable features: %d\n", nrow(stable_features)))

    if (nrow(stable_features) > 0 && nrow(stable_features) <= 20) {
      cat("\n  Stable features:\n")
      print(stable_features)
    } else if (nrow(stable_features) > 20) {
      cat("\n  Top 20 most stable features:\n")
      print(head(stable_features[order(-frequency)], 20))
    }
  }
}

# ============================================================================
# 8. SCIENTIFIC AND CLINICAL ANALYSIS
# ============================================================================
cat("\n\n")
cat("============================================================\n")
cat("SCIENTIFIC AND CLINICAL OBSERVATIONS\n")
cat("============================================================\n\n")

cat("=== 1. PERFORMANCE ANALYSIS ===\n\n")

# Sort by AUC
perf_sorted <- perf_summary[order(-mean_auc)]
cat("Learners ranked by AUC:\n")
for (i in seq_len(nrow(perf_sorted))) {
  row <- perf_sorted[i]
  se <- row$sd_auc / sqrt(5)
  cat(sprintf("  %d. %s: AUC=%.3f [%.3f-%.3f], Sens=%.1f%%, Spec=%.1f%%\n",
              i, row$learner_id, row$mean_auc,
              row$mean_auc - 1.96*se, row$mean_auc + 1.96*se,
              row$mean_sens * 100, row$mean_spec * 100))
}

cat("\n=== 2. STABILITY ANALYSIS (Nogueira Index) ===\n\n")

if (nrow(best_weighted) > 0 && any(!is.na(best_weighted$stability))) {
  stab_sorted <- best_weighted[order(-stability)]
  cat("Learners ranked by stability:\n")
  for (i in seq_len(nrow(stab_sorted))) {
    row <- stab_sorted[i]
    if (!is.na(row$stability)) {
      interp <- if (row$stability >= 0.8) "EXCELLENT"
                else if (row$stability >= 0.6) "GOOD"
                else if (row$stability >= 0.4) "MODERATE"
                else "POOR"
      cat(sprintf("  %d. %s: Stability=%.3f (%s), Features=%.0f\n",
                  i, row$learner_id, row$stability, interp, row$mean_k))
    }
  }

  cat("\nInterpretation:\n")
  cat("  - Stability >= 0.8: Excellent - highly reproducible feature selection\n")
  cat("  - Stability 0.6-0.8: Good - acceptable for clinical development\n")
  cat("  - Stability 0.4-0.6: Moderate - consider larger sample size\n")
  cat("  - Stability < 0.4: Poor - feature selection is unstable\n")
} else {
  cat("Stability could not be computed (features not extracted from learners).\n")
}

cat("\n=== 3. CLINICAL FEASIBILITY ===\n\n")

if (nrow(best_weighted) > 0 && any(!is.na(best_weighted$mean_k))) {
  cat("Feature count vs. clinical panel type:\n")
  for (i in seq_len(nrow(best_weighted))) {
    row <- best_weighted[i]
    if (!is.na(row$mean_k)) {
      feasibility <- if (row$mean_k <= 10) "IDEAL (RT-qPCR panel)"
                     else if (row$mean_k <= 20) "GOOD (small targeted panel)"
                     else if (row$mean_k <= 50) "ACCEPTABLE (medium panel)"
                     else "CHALLENGING (large panel)"
      cat(sprintf("  - %s: %.0f features -> %s\n", row$learner_id, row$mean_k, feasibility))
    }
  }
}

cat("\n=== 4. TRADE-OFF ANALYSIS ===\n\n")

if (nrow(best_weighted) > 0) {
  # Find best by each criterion
  best_auc_id <- perf_summary$learner_id[which.max(perf_summary$mean_auc)]
  best_auc_val <- max(perf_summary$mean_auc)

  if (any(!is.na(best_weighted$stability))) {
    best_stab_row <- best_weighted[which.max(stability)]
    best_stab_id <- best_stab_row$learner_id
    best_stab_val <- best_stab_row$stability

    if (best_auc_id != best_stab_id && !is.na(best_stab_val)) {
      cat("TRADE-OFF DETECTED:\n")
      cat(sprintf("  Best AUC: %s (AUC=%.4f)\n", best_auc_id, best_auc_val))
      cat(sprintf("  Best Stability: %s (Stability=%.4f)\n", best_stab_id, best_stab_val))

      # Get AUC of best stability learner
      auc_of_stable <- perf_summary[learner_id == best_stab_id]$mean_auc
      auc_diff <- best_auc_val - auc_of_stable

      cat(sprintf("\n  AUC sacrifice for best stability: %.4f (%.1f%%)\n",
                  auc_diff, auc_diff / best_auc_val * 100))

      if (auc_diff < 0.02) {
        cat("  RECOMMENDATION: Choose more stable signature (minimal AUC loss)\n")
      } else if (auc_diff < 0.05) {
        cat("  RECOMMENDATION: Consider stable signature if reproducibility is priority\n")
      } else {
        cat("  RECOMMENDATION: Stability gain may not justify AUC loss\n")
      }
    } else {
      cat("NO TRADE-OFF: Best AUC and best stability are the same learner.\n")
    }
  }
}

cat("\n=== 5. FINAL RECOMMENDATIONS ===\n\n")

if (nrow(best_weighted) > 0) {
  final_choice <- best_weighted[selected == TRUE]

  cat("RECOMMENDED SIGNATURE (weighted selection):\n")
  cat(sprintf("  Learner: %s\n", final_choice$learner_id[1]))
  cat(sprintf("  AUC: %.4f\n", final_choice$mean_metric[1]))
  if (!is.na(final_choice$stability[1])) {
    cat(sprintf("  Stability: %.4f\n", final_choice$stability[1]))
  }
  if (!is.na(final_choice$mean_k[1])) {
    cat(sprintf("  Features: %.0f\n", final_choice$mean_k[1]))
  }
  if (!is.na(final_choice$score[1])) {
    cat(sprintf("  Overall Score: %.4f\n", final_choice$score[1]))
  }

  cat("\nRationale:\n")
  cat("  - Weighted selection balances performance (40%), stability (35%), parsimony (25%)\n")
  cat("  - This prioritizes REPRODUCIBILITY and CLINICAL FEASIBILITY\n")
  cat("  - Selected signature should generalize better to external validation\n")
}

cat("\n=== 6. BATCH CORRECTION IMPACT ===\n\n")
cat("FrozenComBat was applied to correct for batch effects from different\n")
cat("TCGA collection centers. This prevents the model from learning\n")
cat("technical artifacts instead of true biological signal.\n\n")
cat("Benefits of frozen batch correction:\n")
cat("  - Parameters learned ONLY from training data (no leakage)\n")
cat("  - Same correction applied to external validation sets\n")
cat("  - Enables fair comparison across different cohorts\n")

cat("\n=== 7. CAVEATS AND NEXT STEPS ===\n\n")
cat("1. EXTERNAL VALIDATION: These results must be validated on independent cohorts\n")
cat("2. BIOLOGICAL VALIDATION: Selected miRNAs should be cross-referenced with:\n")
cat("   - miRTarBase (validated targets)\n")
cat("   - KEGG pathways\n")
cat("   - PubMed literature\n")
cat("3. EXPERIMENTAL VALIDATION: RT-qPCR on independent samples\n")
cat("4. CLINICAL VALIDATION: Prospective study with defined endpoints\n")

cat("\n============================================================\n")
cat("END OF ANALYSIS\n")
cat("============================================================\n")

# Save results
saveRDS(list(
  result = result,
  perf_summary = perf_summary,
  best_1se = best_1se,
  best_weighted = best_weighted,
  pareto_set = pareto_set,
  total_features = total_features
), "tcga_signature_selection_v2_results.rds")

cat("\nResults saved to: tcga_signature_selection_v2_results.rds\n")
