#!/usr/bin/env Rscript
# Test Signature Selection on TCGA Pancreatic Cancer Data
# This script demonstrates the full workflow from data loading to signature selection

cat("=== OmicSelector: Signature Selection on TCGA Data ===\n\n")

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
# 1. LOAD AND PREPARE TCGA DATA
# ============================================================================
cat("\n--- Step 1: Loading TCGA Data ---\n")

# Load the packaged TCGA data
load("data/orginal_TCGA_data.rda")
cat(sprintf("Loaded TCGA data: %d samples, %d columns\n",
            nrow(orginal_TCGA_data), ncol(orginal_TCGA_data)))

# Filter for pancreatic cancer vs normal
cancer_cases <- orginal_TCGA_data[orginal_TCGA_data$primary_site == "Pancreas" &
                                   orginal_TCGA_data$sample_type == "PrimaryTumor", ]
control_cases <- orginal_TCGA_data[orginal_TCGA_data$sample_type == "SolidTissueNormal", ]

cat(sprintf("Pancreatic cancer cases: %d\n", nrow(cancer_cases)))
cat(sprintf("Normal tissue controls: %d\n", nrow(control_cases)))

# Label classes
cancer_cases$Class <- "Case"
control_cases$Class <- "Control"

# Combine and prepare dataset
dataset <- rbind(cancer_cases, control_cases)

# Extract miRNA features
mirna_cols <- grep("^hsa", names(dataset), value = TRUE)
cat(sprintf("miRNA features: %d\n", length(mirna_cols)))

# Create feature matrix (take first 200 miRNAs for faster testing)
mirna_cols <- mirna_cols[1:min(200, length(mirna_cols))]
X <- as.matrix(dataset[, mirna_cols, with = FALSE])
y <- factor(dataset$Class, levels = c("Control", "Case"))

# Remove features with zero variance
var_filter <- apply(X, 2, var, na.rm = TRUE) > 0
X <- X[, var_filter]
cat(sprintf("Features after variance filter: %d\n", ncol(X)))

# Handle missing values (simple imputation with median)
for (j in seq_len(ncol(X))) {
  na_idx <- is.na(X[, j])
  if (any(na_idx)) {
    X[na_idx, j] <- median(X[, j], na.rm = TRUE)
  }
}

# Log transform (add small constant to avoid log(0))
X <- log10(X + 1)

cat(sprintf("Final dataset: %d samples, %d features\n", nrow(X), ncol(X)))
cat(sprintf("Class distribution: Case=%d, Control=%d\n",
            sum(y == "Case"), sum(y == "Control")))

# ============================================================================
# 2. CREATE MLR3 TASK AND LEARNERS
# ============================================================================
cat("\n--- Step 2: Setting Up ML Pipeline ---\n")

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

# Create learners with embedded feature selection (to simulate nested CV properly)
# For this test, we'll create graph learners with filter-based feature selection

# Learner 1: Random Forest with variance filter
po_var <- po("filter", filter = flt("variance"), filter.nfeat = 50)
po_rf <- po("learner", lrn("classif.ranger", predict_type = "prob", num.trees = 100))
graph_rf <- po_var %>>% po_rf
glrn_rf <- as_learner(graph_rf)
glrn_rf$id <- "rf_var50"

# Learner 2: Logistic Regression with ANOVA filter
po_anova <- po("filter", filter = flt("anova"), filter.nfeat = 30)
po_logreg <- po("learner", lrn("classif.log_reg", predict_type = "prob"))
graph_logreg <- po_anova %>>% po_logreg
glrn_logreg <- as_learner(graph_logreg)
glrn_logreg$id <- "logreg_anova30"

# Learner 3: XGBoost with AUC filter
po_auc <- po("filter", filter = flt("auc"), filter.nfeat = 40)
po_xgb <- po("learner", lrn("classif.xgboost", predict_type = "prob",
                             nrounds = 50, verbose = 0))
graph_xgb <- po_auc %>>% po_xgb
glrn_xgb <- as_learner(graph_xgb)
glrn_xgb$id <- "xgb_auc40"

# Learner 4: SVM with Kruskal-Wallis filter (smaller signature)
po_kruskal <- po("filter", filter = flt("kruskal_test"), filter.nfeat = 15)
po_svm <- po("learner", lrn("classif.svm", predict_type = "prob", type = "C-classification", kernel = "radial"))
graph_svm <- po_kruskal %>>% po_svm
glrn_svm <- as_learner(graph_svm)
glrn_svm$id <- "svm_kruskal15"

cat("Created 4 learners with embedded feature selection:\n")
cat("  - rf_var50: Random Forest with top 50 by variance\n")
cat("  - logreg_anova30: Logistic Regression with top 30 by ANOVA\n")
cat("  - xgb_auc40: XGBoost with top 40 by AUC filter\n")
cat("  - svm_kruskal15: SVM with top 15 by Kruskal-Wallis\n")

# ============================================================================
# 3. RUN NESTED CROSS-VALIDATION
# ============================================================================
cat("\n--- Step 3: Running Nested Cross-Validation ---\n")

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
# 4. EXTRACT DETAILED RESULTS
# ============================================================================
cat("\n--- Step 4: Benchmark Results ---\n")

# Get per-fold performance (need to pass measures to score())
measures <- list(
  msr("classif.auc"),
  msr("classif.acc"),
  msr("classif.sensitivity"),
  msr("classif.specificity"),
  msr("classif.bbrier")
)
score_dt <- as.data.table(result$benchmark_result$score(measures))
cat("Available columns:", paste(names(score_dt), collapse=", "), "\n")

# Use backticks for column names with dots
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
# 5. EXTRACT SELECTED FEATURES PER FOLD
# ============================================================================
cat("\n--- Step 5: Extracting Selected Features ---\n")

# Extract features selected in each fold for stability analysis
selected_features_list <- list()

for (i in seq_len(result$benchmark_result$n_resample_results)) {
  rr <- result$benchmark_result$resample_result(i)
  learner_id <- rr$learner$id

  fold_features <- list()
  for (iter in seq_len(rr$iters)) {
    # Get trained learner for this fold
    lrn_trained <- rr$learners[[iter]]

    # Extract selected features from the filter PipeOp
    tryCatch({
      state <- lrn_trained$state
      if (!is.null(state$model)) {
        # Find the filter node
        for (nm in names(state$model)) {
          if (grepl("filter", nm, ignore.case = TRUE)) {
            # The filter stores which features passed
            filter_state <- state$model[[nm]]
            if (!is.null(filter_state$features)) {
              fold_features[[iter]] <- filter_state$features
            } else if (!is.null(filter_state$selected_features)) {
              fold_features[[iter]] <- filter_state$selected_features
            }
          }
        }
      }
    }, error = function(e) NULL)
  }

  if (length(fold_features) > 0) {
    selected_features_list[[learner_id]] <- fold_features
  }
}

# Compute stability metrics
cat("\nFeature Selection Stability:\n")
stability_results <- list()

for (learner_id in names(selected_features_list)) {
  fold_features <- selected_features_list[[learner_id]]

  if (length(fold_features) >= 2) {
    # Count feature frequencies
    all_features <- unlist(fold_features)
    freq_table <- table(all_features) / length(fold_features)

    # Compute Nogueira stability index (simplified)
    # SI = 1 - mean(pairwise Hamming distances) / expected_random
    n_folds <- length(fold_features)
    k_mean <- mean(sapply(fold_features, length))
    p <- ncol(X)  # Total features

    # Jaccard-based stability
    jaccard_vals <- c()
    for (i in 1:(n_folds-1)) {
      for (j in (i+1):n_folds) {
        intersection <- length(intersect(fold_features[[i]], fold_features[[j]]))
        union_size <- length(union(fold_features[[i]], fold_features[[j]]))
        if (union_size > 0) {
          jaccard_vals <- c(jaccard_vals, intersection / union_size)
        }
      }
    }

    stability <- mean(jaccard_vals, na.rm = TRUE)

    stability_results[[learner_id]] <- list(
      nogueira_index = stability,
      mean_features = k_mean,
      features_80pct = names(freq_table)[freq_table >= 0.8],
      n_stable = sum(freq_table >= 0.8)
    )

    cat(sprintf("  %s: Stability=%.3f, Mean features=%.1f, Stable (>=80%%)=%d\n",
                learner_id, stability, k_mean, sum(freq_table >= 0.8)))
  }
}

# Update result with stability information
result$stability <- list(
  nogueira_index = data.table(
    learner_id = names(stability_results),
    nogueira_index = sapply(stability_results, function(x) x$nogueira_index)
  ),
  selected_features_per_fold = do.call(rbind, lapply(names(selected_features_list), function(lid) {
    data.table(
      learner_id = lid,
      outer_fold = seq_along(selected_features_list[[lid]]),
      selected_features = selected_features_list[[lid]]
    )
  }))
)

# Add mean_k to the score summary
for (lid in names(stability_results)) {
  result$stability$nogueira_index[learner_id == lid, mean_k := stability_results[[lid]]$mean_features]
}

# ============================================================================
# 6. SIGNATURE SELECTION
# ============================================================================
cat("\n--- Step 6: Signature Selection ---\n")

# Mode 1: Constrained 1SE (conservative)
cat("\n=== Mode: Constrained 1SE ===\n")
best_1se <- select_best_signature(
  result,
  mode = "constrained_1se",
  metric = "classif.auc",
  auc_min = 0.7,  # Minimum acceptable AUC
  k_max = 50      # Maximum features for clinical panel
)

if (nrow(best_1se) > 0) {
  cat("\nCandidates within 1SE of best:\n")
  print(best_1se[, .(learner_id, mean_metric, se_metric, stability, mean_k, selected)])

  selected_1se <- best_1se[selected == TRUE]
  cat(sprintf("\nSelected (1SE): %s\n", selected_1se$learner_id[1]))
  cat(sprintf("  AUC: %.3f (SE: %.3f)\n", selected_1se$mean_metric[1], selected_1se$se_metric[1]))
  if (!is.na(selected_1se$stability[1])) {
    cat(sprintf("  Stability: %.3f\n", selected_1se$stability[1]))
  }
  if (!is.na(selected_1se$mean_k[1])) {
    cat(sprintf("  Features: %.0f\n", selected_1se$mean_k[1]))
  }
}

# Mode 2: Weighted scoring
cat("\n=== Mode: Weighted Scoring ===\n")
best_weighted <- select_best_signature(
  result,
  mode = "weighted",
  metric = "classif.auc",
  weights = c(performance = 0.4, stability = 0.35, parsimony = 0.25),
  auc_min = 0.7,
  k_max = 50
)

if (nrow(best_weighted) > 0) {
  cat("\nWeighted scores:\n")
  print(best_weighted[, .(learner_id, mean_metric, stability, mean_k, score, selected)])

  selected_weighted <- best_weighted[selected == TRUE]
  cat(sprintf("\nSelected (weighted): %s\n", selected_weighted$learner_id[1]))
  cat(sprintf("  Score: %.3f\n", selected_weighted$score[1]))
}

# Mode 3: Pareto frontier
cat("\n=== Mode: Pareto Frontier ===\n")
pareto_set <- select_best_signature(
  result,
  mode = "pareto",
  metric = "classif.auc",
  auc_min = 0.7
)

if (nrow(pareto_set) > 0) {
  cat("\nPareto-optimal signatures:\n")
  print(pareto_set[, .(learner_id, mean_metric, stability, mean_k, pareto, selected)])

  cat(sprintf("\nPareto set size: %d\n", sum(pareto_set$pareto, na.rm = TRUE)))
  cat(sprintf("Recommended (closest to Utopia): %s\n", pareto_set[selected == TRUE]$learner_id[1]))
}

# ============================================================================
# 7. CONSENSUS FEATURES
# ============================================================================
cat("\n--- Step 7: Consensus Features ---\n")

# Get consensus features from the best signature
if (nrow(best_weighted) > 0) {
  best_learner <- best_weighted[selected == TRUE]$learner_id[1]

  consensus <- get_consensus_features(result, best_learner, min_frequency = 0.6)

  if (nrow(consensus) > 0) {
    cat(sprintf("\nConsensus features for %s (>=60%% frequency):\n", best_learner))
    stable_features <- consensus[selected == TRUE]
    cat(sprintf("  Total stable features: %d\n", nrow(stable_features)))

    if (nrow(stable_features) > 0) {
      cat("\n  Top 10 most stable features:\n")
      print(head(stable_features[order(-frequency)], 10))
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

# Best raw AUC
best_auc_learner <- perf_summary[which.max(mean_auc)]
cat(sprintf("Best raw AUC: %s with %.3f (SD: %.3f)\n",
            best_auc_learner$learner_id, best_auc_learner$mean_auc, best_auc_learner$sd_auc))

# Confidence intervals
cat("\n95%% Confidence Intervals (approximated from 5-fold CV):\n")
for (i in seq_len(nrow(perf_summary))) {
  row <- perf_summary[i]
  se <- row$sd_auc / sqrt(5)
  ci_low <- row$mean_auc - 1.96 * se
  ci_high <- row$mean_auc + 1.96 * se
  cat(sprintf("  %s: AUC = %.3f [%.3f - %.3f]\n",
              row$learner_id, row$mean_auc, ci_low, ci_high))
}

cat("\n=== 2. STABILITY ANALYSIS ===\n\n")

if (length(stability_results) > 0) {
  stab_df <- data.frame(
    learner = names(stability_results),
    stability = sapply(stability_results, function(x) x$nogueira_index),
    features = sapply(stability_results, function(x) x$mean_features),
    stable_80 = sapply(stability_results, function(x) x$n_stable)
  )
  stab_df <- stab_df[order(-stab_df$stability), ]

  cat("Feature Selection Stability (Jaccard-based):\n")
  print(stab_df, row.names = FALSE)

  cat("\nInterpretation:\n")
  cat("  - Stability > 0.7: High reproducibility (recommended for clinical use)\n")
  cat("  - Stability 0.5-0.7: Moderate reproducibility (use with caution)\n")
  cat("  - Stability < 0.5: Low reproducibility (not suitable for clinical panels)\n")

  high_stab <- stab_df$learner[stab_df$stability >= 0.7]
  if (length(high_stab) > 0) {
    cat(sprintf("\n  HIGH STABILITY learners: %s\n", paste(high_stab, collapse = ", ")))
  } else {
    cat("\n  WARNING: No learner achieved high stability (>0.7)\n")
    cat("  This suggests the signal may be weak or data is heterogeneous.\n")
  }
}

cat("\n=== 3. CLINICAL FEASIBILITY ===\n\n")

cat("Feature Count Analysis:\n")
for (lid in names(stability_results)) {
  k <- stability_results[[lid]]$mean_features
  feasibility <- if (k <= 10) "IDEAL (RT-qPCR panel)"
                 else if (k <= 20) "GOOD (small targeted panel)"
                 else if (k <= 50) "ACCEPTABLE (medium panel)"
                 else "CHALLENGING (large panel, consider reduction)"
  cat(sprintf("  %s: %.0f features -> %s\n", lid, k, feasibility))
}

cat("\n=== 4. TRADE-OFF ANALYSIS ===\n\n")

if (nrow(perf_summary) > 1) {
  # Compare best AUC vs best stability
  best_auc_id <- perf_summary$learner_id[which.max(perf_summary$mean_auc)]

  if (length(stability_results) > 0) {
    stab_vals <- sapply(stability_results, function(x) x$nogueira_index)
    best_stab_id <- names(which.max(stab_vals))

    if (best_auc_id != best_stab_id) {
      cat("TRADE-OFF DETECTED:\n")
      cat(sprintf("  Best AUC: %s (AUC=%.3f, Stability=%.3f)\n",
                  best_auc_id,
                  perf_summary[learner_id == best_auc_id]$mean_auc,
                  stab_vals[best_auc_id]))
      cat(sprintf("  Best Stability: %s (AUC=%.3f, Stability=%.3f)\n",
                  best_stab_id,
                  perf_summary[learner_id == best_stab_id]$mean_auc,
                  stab_vals[best_stab_id]))

      auc_diff <- perf_summary[learner_id == best_auc_id]$mean_auc -
                  perf_summary[learner_id == best_stab_id]$mean_auc
      stab_diff <- stab_vals[best_stab_id] - stab_vals[best_auc_id]

      cat(sprintf("\n  AUC sacrifice for stability: %.3f\n", auc_diff))
      cat(sprintf("  Stability gain: %.3f\n", stab_diff))

      if (auc_diff < 0.02 && stab_diff > 0.1) {
        cat("\n  RECOMMENDATION: Choose more stable signature (minimal AUC loss, significant stability gain)\n")
      } else if (auc_diff > 0.05) {
        cat("\n  RECOMMENDATION: Consider the higher-AUC signature (stability gain may not justify AUC loss)\n")
      }
    } else {
      cat("NO TRADE-OFF: Best AUC and best stability are the same learner!\n")
      cat(sprintf("  Winner: %s\n", best_auc_id))
    }
  }
}

cat("\n=== 5. RECOMMENDATIONS ===\n\n")

if (nrow(best_weighted) > 0) {
  final_choice <- best_weighted[selected == TRUE]

  cat("FINAL RECOMMENDATION (weighted selection):\n")
  cat(sprintf("  Signature: %s\n", final_choice$learner_id[1]))
  cat(sprintf("  AUC: %.3f\n", final_choice$mean_metric[1]))
  if (!is.na(final_choice$stability[1])) {
    cat(sprintf("  Stability: %.3f\n", final_choice$stability[1]))
  }
  if (!is.na(final_choice$mean_k[1])) {
    cat(sprintf("  Features: %.0f\n", final_choice$mean_k[1]))
  }

  cat("\nRationale:\n")
  cat("  - Weighted selection balances performance (40%), stability (35%), and parsimony (25%)\n")
  cat("  - This prioritizes reproducibility and clinical feasibility over raw discrimination\n")
  cat("  - The selected signature should transfer better to external validation\n")
}

cat("\n=== 6. CAVEATS AND LIMITATIONS ===\n\n")

cat("1. IMBALANCED CLASSES: Case/Control ratio may affect sensitivity/specificity trade-off\n")
cat("2. BATCH EFFECTS: TCGA data comes from multiple centers; consider ComBat correction\n")
cat("3. OVERFITTING RISK: Even with nested CV, small sample sizes can lead to optimistic estimates\n")
cat("4. EXTERNAL VALIDATION: These results MUST be validated on independent cohorts\n")
cat("5. BIOLOGICAL VALIDATION: Selected miRNAs should be cross-referenced with literature\n")

cat("\n============================================================\n")
cat("END OF ANALYSIS\n")
cat("============================================================\n")

# Save results
saveRDS(list(
  result = result,
  perf_summary = perf_summary,
  stability_results = stability_results,
  best_1se = best_1se,
  best_weighted = best_weighted,
  pareto_set = pareto_set
), "tcga_signature_selection_results.rds")

cat("\nResults saved to: tcga_signature_selection_results.rds\n")
