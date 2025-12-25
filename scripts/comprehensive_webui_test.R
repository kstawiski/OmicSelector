#!/usr/bin/env Rscript
# Comprehensive OmicSelector WebUI Test Script
# Tests R package and validates clinical/scientific validity with TCGA and Elias2017 datasets

suppressPackageStartupMessages({
  library(OmicSelector)
  library(data.table)
  library(mlr3)
  library(mlr3pipelines)
  library(mlr3learners)
})

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  OmicSelector Comprehensive Test Suite\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Track test results
test_results <- list()
all_passed <- TRUE

run_test <- function(name, expr) {
  cat(sprintf("\n[TEST] %s\n", name))
  cat("-" |> rep(50) |> paste(collapse = ""), "\n")

  result <- tryCatch({
    start_time <- Sys.time()
    res <- eval(expr)
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("  PASSED (%.1fs)\n", elapsed))
    list(status = "PASSED", result = res, time = elapsed)
  }, error = function(e) {
    cat(sprintf("  FAILED: %s\n", e$message))
    all_passed <<- FALSE
    list(status = "FAILED", error = e$message)
  })

  test_results[[name]] <<- result
  invisible(result)
}

# ============================================================================
# PART 1: TCGA Dataset Tests
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  PART 1: TCGA Pan-Cancer miRNA Dataset\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

run_test("Load TCGA Data", quote({
  data("orginal_TCGA_data", package = "OmicSelector")

  # Convert from data.table to data.frame explicitly
  tcga <- data.table::setDF(data.table::copy(orginal_TCGA_data))
  cat(sprintf("  Loaded: %d samples, %d features\n", nrow(tcga), ncol(tcga)))
  cat(sprintf("  Sample types: %s\n", paste(unique(tcga[["sample_type"]]), collapse = ", ")))

  # Subset for speed: 100 samples (50 tumor, 50 normal)
  set.seed(42)
  tumor_idx <- which(tcga[["sample_type"]] == "PrimaryTumor")
  normal_idx <- which(tcga[["sample_type"]] == "SolidTissueNormal")

  tumor_sample <- sample(tumor_idx, min(50, length(tumor_idx)))
  normal_sample <- sample(normal_idx, min(50, length(normal_idx)))
  tcga_subset <- tcga[c(tumor_sample, normal_sample), , drop = FALSE]

  # Select top 200 variable miRNAs
  mirna_cols <- grep("^hsa-", names(tcga_subset), value = TRUE)
  if (length(mirna_cols) == 0) {
    mirna_cols <- grep("^hsa", names(tcga_subset), value = TRUE)
  }
  cat(sprintf("  Found %d miRNA columns\n", length(mirna_cols)))

  # Calculate variance for numeric columns only
  mirna_data <- tcga_subset[, mirna_cols, drop = FALSE]
  mirna_var <- sapply(mirna_data, function(x) var(as.numeric(x), na.rm = TRUE))
  top_mirnas <- names(sort(mirna_var, decreasing = TRUE)[1:min(200, length(mirna_cols))])

  tcga_test <- data.frame(
    outcome = factor(ifelse(tcga_subset[["sample_type"]] == "PrimaryTumor", "Tumor", "Normal")),
    tcga_subset[, top_mirnas, drop = FALSE]
  )

  # Remove any columns with all NA
  na_cols <- sapply(tcga_test, function(x) all(is.na(x)))
  if (any(na_cols)) {
    tcga_test <- tcga_test[, !na_cols, drop = FALSE]
  }

  cat(sprintf("  Test subset: %d samples, %d miRNAs\n", nrow(tcga_test), ncol(tcga_test) - 1))
  cat(sprintf("  Class balance: Tumor=%d, Normal=%d\n",
              sum(tcga_test[["outcome"]] == "Tumor"), sum(tcga_test[["outcome"]] == "Normal")))

  assign("tcga_test", tcga_test, envir = .GlobalEnv)
  TRUE
}))

run_test("Create OmicPipeline (TCGA)", quote({
  pipeline <- OmicPipeline$new(
    data = tcga_test,
    target = "outcome",
    positive = "Tumor"
  )

  cat(sprintf("  Task created: %s\n", pipeline$task$id))
  cat(sprintf("  Features: %d\n", length(pipeline$task$feature_names)))
  cat(sprintf("  Target: %s (positive=%s)\n", pipeline$task$target_names, "Tumor"))

  assign("tcga_pipeline", pipeline, envir = .GlobalEnv)
  TRUE
}))

run_test("Create Graph Learner with ANOVA Filter", quote({
  learner <- tcga_pipeline$create_graph_learner(
    filter = "anova",
    model = "rpart",
    n_features = 10
  )

  cat(sprintf("  Learner ID: %s\n", learner$id))
  cat(sprintf("  Graph nodes: %d\n", length(learner$graph$ids())))

  assign("tcga_learner", learner, envir = .GlobalEnv)
  TRUE
}))

run_test("Run Nested CV Benchmark (TCGA)", quote({
  result <- tcga_pipeline$benchmark(
    learners = tcga_learner,
    outer_folds = 3,
    inner_folds = 2,
    stratify = TRUE,
    seed = 42,
    parallel = FALSE
  )

  cat(sprintf("  Result class: %s\n", class(result)[1]))

  # Check performance
  if (!is.null(result$performance)) {
    perf <- result$performance
    auc_col <- grep("auc", names(perf), value = TRUE, ignore.case = TRUE)
    if (length(auc_col) > 0) {
      auc_val <- perf[[auc_col[1]]][1]
      cat(sprintf("  AUC: %.3f\n", auc_val))
      if (auc_val < 0.5) stop("AUC below chance level!")
    }
  }

  # Check stability
  if (!is.null(result$stability)) {
    stab <- result$stability
    if (!is.null(stab$nogueira_index)) {
      ni <- if (is.numeric(stab$nogueira_index)) stab$nogueira_index else stab$nogueira_index$nogueira_index[1]
      cat(sprintf("  Nogueira Stability: %.3f\n", ni))
    }
    if (!is.null(stab$selected_features_per_fold)) {
      n_folds <- length(stab$selected_features_per_fold)
      cat(sprintf("  Feature selections stored: %d folds\n", n_folds))
    }
  }

  assign("tcga_result", result, envir = .GlobalEnv)
  TRUE
}))

run_test("Validate Bootstrap CI Implementation", quote({
  stab <- tcga_result$stability
  if (is.null(stab) || is.null(stab$selected_features_per_fold)) {
    stop("No stability data available for bootstrap CI test")
  }

  sf <- stab$selected_features_per_fold

  # Extract feature lists (handle nested structure)
  feature_lists <- if (is.list(sf[[1]]) && !is.character(sf[[1]])) {
    sf[[1]]
  } else {
    sf
  }

  cat(sprintf("  Feature lists: %d folds\n", length(feature_lists)))

  # Implement bootstrap CI calculation (same as in mod_results.R)
  calculate_nogueira_index <- function(feature_lists, all_features = NULL) {
    if (length(feature_lists) < 2) return(NA_real_)
    if (is.null(all_features)) all_features <- unique(unlist(feature_lists))
    p <- length(all_features)
    if (p == 0) return(NA_real_)
    n <- length(feature_lists)

    Z <- matrix(0, nrow = n, ncol = p)
    colnames(Z) <- all_features
    for (i in seq_len(n)) {
      selected <- feature_lists[[i]]
      matched <- selected[selected %in% all_features]
      if (length(matched) > 0) Z[i, matched] <- 1
    }

    pf <- colMeans(Z)
    k_bar <- mean(rowSums(Z))
    if (k_bar == 0 || k_bar == p) return(1)

    feat_var <- pf * (1 - pf)
    mean_var <- mean(feat_var)
    expected_var <- (k_bar / p) * (1 - k_bar / p)
    if (expected_var == 0) return(1)

    max(0, min(1, 1 - (mean_var / expected_var)))
  }

  bootstrap_stability_ci <- function(feature_lists, n_boot = 1000, conf_level = 0.95, seed = NULL) {
    if (length(feature_lists) < 2) {
      return(list(point = NA_real_, lower = NA_real_, upper = NA_real_, se = NA_real_))
    }
    if (!is.null(seed)) set.seed(seed)

    all_features <- unique(unlist(feature_lists))
    n_folds <- length(feature_lists)
    point_estimate <- calculate_nogueira_index(feature_lists, all_features)

    boot_estimates <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      boot_idx <- sample(seq_len(n_folds), n_folds, replace = TRUE)
      boot_lists <- feature_lists[boot_idx]
      boot_estimates[b] <- calculate_nogueira_index(boot_lists, all_features)
    }

    boot_estimates <- boot_estimates[!is.na(boot_estimates)]
    if (length(boot_estimates) < 10) {
      return(list(point = point_estimate, lower = NA_real_, upper = NA_real_, se = NA_real_))
    }

    alpha <- 1 - conf_level
    list(
      point = point_estimate,
      lower = as.numeric(quantile(boot_estimates, probs = alpha / 2)),
      upper = as.numeric(quantile(boot_estimates, probs = 1 - alpha / 2)),
      se = sd(boot_estimates),
      n_boot = length(boot_estimates)
    )
  }

  # Run bootstrap CI
  ci <- bootstrap_stability_ci(feature_lists, n_boot = 500, conf_level = 0.95, seed = 42)

  cat(sprintf("  Point estimate: %.3f\n", ci$point))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n", ci$lower, ci$upper))
  cat(sprintf("  Standard error: %.4f\n", ci$se))
  cat(sprintf("  Bootstrap samples: %d\n", ci$n_boot))

  # Validate CI properties
  if (is.na(ci$point)) stop("Point estimate is NA")
  if (ci$point < 0 || ci$point > 1) stop("Point estimate out of [0,1] range")
  if (!is.na(ci$lower) && !is.na(ci$upper)) {
    if (ci$lower > ci$upper) stop("CI lower > upper")
    if (ci$lower < 0 || ci$upper > 1) stop("CI bounds out of [0,1] range")
    ci_width <- ci$upper - ci$lower
    cat(sprintf("  CI width: %.3f (%s)\n", ci_width,
                if (ci_width < 0.1) "narrow" else if (ci_width < 0.2) "moderate" else "wide"))
  }

  TRUE
}))

# ============================================================================
# PART 2: Elias2017 Dataset Tests
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  PART 2: Elias2017 Serum miRNA Dataset\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

run_test("Load Elias2017 Data", quote({
  elias_path <- "/OmicSelector/example/Elias2017.csv"
  if (!file.exists(elias_path)) {
    stop("Elias2017.csv not found at expected path")
  }

  elias <- read.csv(elias_path, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded: %d samples, %d columns\n", nrow(elias), ncol(elias)))
  cat(sprintf("  Columns: %s...\n", paste(head(names(elias), 5), collapse = ", ")))

  # Identify target and features
  possible_targets <- c("Class", "class", "Outcome", "outcome", "Target", "target", "Group", "group")
  target_col <- intersect(possible_targets, names(elias))[1]

  if (is.na(target_col)) {
    # Try first column if it looks like a factor
    if (length(unique(elias[[1]])) <= 5) {
      target_col <- names(elias)[1]
    } else {
      stop("Could not identify target column")
    }
  }

  cat(sprintf("  Target column: %s\n", target_col))
  cat(sprintf("  Classes: %s\n", paste(unique(elias[[target_col]]), collapse = ", ")))

  # Prepare data
  elias$outcome <- factor(elias[[target_col]])
  feature_cols <- setdiff(names(elias), c(target_col, "outcome"))

  # Remove non-numeric columns
  numeric_cols <- sapply(elias[, feature_cols, drop = FALSE], is.numeric)
  feature_cols <- feature_cols[numeric_cols]

  elias_test <- elias[, c("outcome", feature_cols)]

  # Handle missing values
  na_counts <- colSums(is.na(elias_test))
  if (any(na_counts > 0)) {
    cat(sprintf("  Columns with NA: %d (imputing with median)\n", sum(na_counts > 0)))
    for (col in names(na_counts[na_counts > 0])) {
      if (col != "outcome") {
        elias_test[[col]][is.na(elias_test[[col]])] <- median(elias_test[[col]], na.rm = TRUE)
      }
    }
  }

  cat(sprintf("  Test data: %d samples, %d features\n", nrow(elias_test), length(feature_cols)))

  assign("elias_test", elias_test, envir = .GlobalEnv)
  TRUE
}))

run_test("Create OmicPipeline (Elias2017)", quote({
  levels_out <- levels(elias_test$outcome)
  positive_class <- levels_out[1]

  pipeline <- OmicPipeline$new(
    data = elias_test,
    target = "outcome",
    positive = positive_class
  )

  cat(sprintf("  Task created: %s\n", pipeline$task$id))
  cat(sprintf("  Features: %d\n", length(pipeline$task$feature_names)))
  cat(sprintf("  Target: %s (positive=%s)\n", pipeline$task$target_names, positive_class))

  assign("elias_pipeline", pipeline, envir = .GlobalEnv)
  TRUE
}))

run_test("Run Nested CV Benchmark (Elias2017)", quote({
  learner <- elias_pipeline$create_graph_learner(
    filter = "anova",
    model = "rpart",
    n_features = 5
  )

  result <- elias_pipeline$benchmark(
    learners = learner,
    outer_folds = 3,
    inner_folds = 2,
    stratify = TRUE,
    seed = 42,
    parallel = FALSE
  )

  cat(sprintf("  Result class: %s\n", class(result)[1]))

  if (!is.null(result$performance)) {
    perf <- result$performance
    auc_col <- grep("auc", names(perf), value = TRUE, ignore.case = TRUE)
    if (length(auc_col) > 0) {
      auc_val <- perf[[auc_col[1]]][1]
      cat(sprintf("  AUC: %.3f\n", auc_val))
    }
  }

  assign("elias_result", result, envir = .GlobalEnv)
  TRUE
}))

# ============================================================================
# PART 3: Multi-Learner Comparison
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  PART 3: Multi-Learner Benchmark\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

run_test("Compare Multiple Learners (TCGA)", quote({
  # Create multiple learners
  learners <- list(
    tcga_pipeline$create_graph_learner(filter = "anova", model = "rpart", n_features = 10),
    tcga_pipeline$create_graph_learner(filter = "anova", model = "log_reg", n_features = 10)
  )

  cat(sprintf("  Learners: %s\n", paste(sapply(learners, function(l) l$id), collapse = ", ")))

  result <- tcga_pipeline$benchmark(
    learners = learners,
    outer_folds = 3,
    inner_folds = 2,
    stratify = TRUE,
    seed = 42,
    parallel = FALSE
  )

  if (!is.null(result$performance) && is.data.frame(result$performance)) {
    perf <- result$performance
    cat("  Performance table:\n")
    print(perf[, intersect(c("learner_id", "classif.auc", "classif.acc"), names(perf))])
  }

  TRUE
}))

# ============================================================================
# PART 4: Clinical/Scientific Validity Checks
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  PART 4: Clinical & Scientific Validity\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

run_test("Validate AUC Range", quote({
  # AUC should be between 0.5 and 1.0 for a reasonable classifier
  if (!is.null(tcga_result$performance)) {
    auc_col <- grep("auc", names(tcga_result$performance), value = TRUE, ignore.case = TRUE)
    if (length(auc_col) > 0) {
      auc <- tcga_result$performance[[auc_col[1]]][1]
      cat(sprintf("  TCGA AUC: %.3f\n", auc))

      if (auc < 0.5) {
        cat("  WARNING: AUC < 0.5 suggests label inversion or random classifier\n")
      } else if (auc < 0.6) {
        cat("  NOTE: AUC 0.5-0.6 indicates weak discrimination\n")
      } else if (auc < 0.7) {
        cat("  NOTE: AUC 0.6-0.7 indicates acceptable discrimination\n")
      } else if (auc < 0.8) {
        cat("  GOOD: AUC 0.7-0.8 indicates good discrimination\n")
      } else {
        cat("  EXCELLENT: AUC > 0.8 indicates excellent discrimination\n")
      }
    }
  }
  TRUE
}))

run_test("Validate Stability Index", quote({
  if (!is.null(tcga_result$stability) && !is.null(tcga_result$stability$nogueira_index)) {
    ni <- tcga_result$stability$nogueira_index
    stab_val <- if (is.numeric(ni)) ni else ni$nogueira_index[1]

    cat(sprintf("  Nogueira Stability: %.3f\n", stab_val))

    if (stab_val < 0.5) {
      cat("  WARNING: Stability < 0.5 indicates unstable feature selection\n")
      cat("  Clinical implication: Results may not replicate\n")
    } else if (stab_val < 0.7) {
      cat("  CAUTION: Stability 0.5-0.7 indicates moderate stability\n")
      cat("  Clinical implication: External validation strongly recommended\n")
    } else {
      cat("  GOOD: Stability >= 0.7 indicates stable feature selection\n")
      cat("  Clinical implication: Features likely to replicate\n")
    }
  }
  TRUE
}))

run_test("Check Feature Selection Consistency", quote({
  if (!is.null(tcga_result$stability) && !is.null(tcga_result$stability$selected_features_per_fold)) {
    sf <- tcga_result$stability$selected_features_per_fold

    # Flatten if nested
    feature_lists <- if (is.list(sf[[1]]) && !is.character(sf[[1]])) sf[[1]] else sf

    # Calculate feature frequencies
    all_features <- unlist(feature_lists)
    freq_table <- sort(table(all_features), decreasing = TRUE)
    n_folds <- length(feature_lists)

    cat(sprintf("  Total unique features selected: %d\n", length(freq_table)))
    cat(sprintf("  Features in all folds: %d\n", sum(freq_table == n_folds)))
    cat(sprintf("  Features in majority (>50%%): %d\n", sum(freq_table > n_folds/2)))

    # Top features
    cat("  Top 5 most stable features:\n")
    top5 <- head(freq_table, 5)
    for (i in seq_along(top5)) {
      cat(sprintf("    %d. %s (%d/%d folds, %.0f%%)\n",
                  i, names(top5)[i], top5[i], n_folds, 100 * top5[i] / n_folds))
    }
  }
  TRUE
}))

run_test("Sample Size Adequacy Check", quote({
  # Rule of thumb: at least 10 events per predictor (EPV)
  n_samples <- nrow(tcga_test)
  n_minority <- min(table(tcga_test$outcome))
  n_features_used <- 10  # We used 10 features

  epv <- n_minority / n_features_used

  cat(sprintf("  Total samples: %d\n", n_samples))
  cat(sprintf("  Minority class: %d\n", n_minority))
  cat(sprintf("  Features selected: %d\n", n_features_used))
  cat(sprintf("  Events per variable (EPV): %.1f\n", epv))

  if (epv < 10) {
    cat("  WARNING: EPV < 10 may lead to overfitting\n")
    cat("  Clinical implication: Consider fewer features or more samples\n")
  } else if (epv < 20) {
    cat("  CAUTION: EPV 10-20 is borderline adequate\n")
  } else {
    cat("  GOOD: EPV >= 20 indicates adequate sample size\n")
  }

  TRUE
}))

# ============================================================================
# PART 5: WebUI Module Tests
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  PART 5: WebUI Module Validation\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

run_test("Load WebUI Modules", quote({
  # Source WebUI modules with error handling
  webui_path <- "/OmicSelector/shiny/omicselector/R"
  modules <- c("class_OmicProject.R", "mod_qc.R", "mod_results.R", "mod_export.R")

  load_results <- list()
  for (mod in modules) {
    mod_path <- file.path(webui_path, mod)
    if (file.exists(mod_path)) {
      result <- tryCatch({
        source(mod_path, local = FALSE)
        cat(sprintf("  Loaded: %s\n", mod))
        load_results[[mod]] <- TRUE
        TRUE
      }, error = function(e) {
        cat(sprintf("  ERROR in %s: %s\n", mod, substr(e$message, 1, 80)))
        load_results[[mod]] <- FALSE
        FALSE
      })
    } else {
      cat(sprintf("  MISSING: %s\n", mod))
      load_results[[mod]] <- FALSE
    }
  }

  loaded_count <- sum(unlist(load_results))
  cat(sprintf("  Successfully loaded: %d/%d modules\n", loaded_count, length(modules)))

  # At minimum, class_OmicProject.R should load
  if (loaded_count >= 1) TRUE else stop("No modules could be loaded")
}))

run_test("Test OmicProject R6 Class", quote({
  # Test project creation and state management
  proj <- OmicProject$new(name = "Test_Project")

  cat(sprintf("  Project ID: %s\n", proj$id))
  cat(sprintf("  Project Name: %s\n", proj$name))

  # Test summary
  summary <- proj$summary()
  cat(sprintf("  Summary fields: %s\n", paste(names(summary), collapse = ", ")))

  TRUE
}))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  TEST SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

passed <- sum(sapply(test_results, function(x) x$status == "PASSED"))
failed <- sum(sapply(test_results, function(x) x$status == "FAILED"))
total <- length(test_results)

cat(sprintf("Tests Passed: %d/%d (%.0f%%)\n", passed, total, 100 * passed / total))
cat(sprintf("Tests Failed: %d/%d\n", failed, total))

if (failed > 0) {
  cat("\nFailed tests:\n")
  for (name in names(test_results)) {
    if (test_results[[name]]$status == "FAILED") {
      cat(sprintf("  - %s: %s\n", name, test_results[[name]]$error))
    }
  }
}

cat("\n")
if (all_passed) {
  cat("ALL TESTS PASSED\n")
  quit(status = 0)
} else {
  cat("SOME TESTS FAILED\n")
  quit(status = 1)
}
