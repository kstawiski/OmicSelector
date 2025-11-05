# Test suite for feature_clustering.R
# Tests OmicSelector_cluster_features() and related functions

library(testthat)

# ===== Test Data Setup =====

# Create synthetic correlated features
create_correlated_data <- function(n_samples = 100, n_clusters = 5, features_per_cluster = 5, seed = 123) {
  set.seed(seed)

  # Create cluster centers
  cluster_centers <- matrix(rnorm(n_samples * n_clusters), ncol = n_clusters)
  colnames(cluster_centers) <- paste0("cluster_", 1:n_clusters)

  # Create correlated features within each cluster
  data_list <- list()
  for (i in 1:n_clusters) {
    for (j in 1:features_per_cluster) {
      # Add noise to cluster center
      feature <- cluster_centers[, i] + rnorm(n_samples, sd = 0.3)
      feature_name <- paste0("feat_c", i, "_", j)
      data_list[[feature_name]] <- feature
    }
  }

  # Add outcome
  outcome <- ifelse(rowSums(cluster_centers[, 1:2]) > median(rowSums(cluster_centers[, 1:2])),
                    "High", "Low")
  data_list$outcome <- factor(outcome)

  # Convert to data frame
  data <- as.data.frame(data_list)

  return(data)
}

test_data <- create_correlated_data()

# ===== Test 1: Input Validation =====

test_that("OmicSelector_cluster_features validates inputs correctly", {

  # Test missing data
  expect_error(
    OmicSelector_cluster_features(
      data = NULL
    ),
    "data"
  )

  # Test invalid method
  expect_error(
    OmicSelector_cluster_features(
      data = test_data,
      method = "invalid_method"
    ),
    "should be one of"
  )

  # Test invalid correlation threshold
  expect_error(
    OmicSelector_cluster_features(
      data = test_data,
      min_correlation = 1.5
    ),
    "min_correlation"
  )

  # Test invalid distance metric
  expect_error(
    OmicSelector_cluster_features(
      data = test_data,
      distance_metric = "invalid"
    ),
    "distance_metric"
  )
})

# ===== Test 2: Correlation-based Clustering =====

test_that("Correlation-based clustering works correctly", {

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],  # Remove outcome
    method = "correlation",
    min_correlation = 0.6,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_cluster_features")

  # Check required elements
  expect_true("clusters" %in% names(result))
  expect_true("representatives" %in% names(result))
  expect_true("cluster_assignments" %in% names(result))
  expect_true("correlation_matrix" %in% names(result))
  expect_true("method" %in% names(result))

  # Check that representatives exist
  expect_type(result$representatives, "character")
  expect_true(length(result$representatives) > 0)

  # Check that all features are assigned to clusters
  expect_equal(length(result$cluster_assignments),
               ncol(test_data) - 1)  # Excluding outcome

  # Method should be recorded
  expect_equal(result$method, "correlation")
})

# ===== Test 3: Hierarchical Clustering =====

test_that("Hierarchical clustering works correctly", {

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "hierarchical",
    n_clusters = 5,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_cluster_features")

  # Check number of clusters
  n_unique_clusters <- length(unique(result$cluster_assignments))
  expect_equal(n_unique_clusters, 5)

  # Check that each cluster has a representative
  expect_equal(length(result$representatives), 5)
})

# ===== Test 4: Expression Pattern Clustering =====

test_that("Expression pattern clustering works correctly", {

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "expression_pattern",
    n_clusters = 5,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_cluster_features")

  # Check structure
  expect_true(length(result$representatives) > 0)
  expect_true(all(result$representatives %in% colnames(test_data)))
})

# ===== Test 5: K-medoids Clustering =====

test_that("K-medoids clustering works correctly", {
  skip_if_not_installed("cluster")

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "kmedoids",
    n_clusters = 5,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(result, "OmicSelector_cluster_features")

  # Check number of representatives
  expect_equal(length(result$representatives), 5)

  # Representatives should be actual features (medoids)
  expect_true(all(result$representatives %in% colnames(test_data)))
})

# ===== Test 6: Auto-detection of Number of Clusters =====

test_that("Auto-detection of clusters works", {

  # When n_clusters = NULL, should auto-detect
  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    n_clusters = NULL,
    min_correlation = 0.7,
    verbose = FALSE
  )

  # Should determine some number of clusters
  n_clusters <- length(unique(result$cluster_assignments))
  expect_true(n_clusters > 0)
  expect_true(n_clusters <= ncol(test_data) - 1)

  # Number of representatives should match number of clusters
  expect_equal(length(result$representatives), n_clusters)
})

# ===== Test 7: Representative Selection Strategies =====

test_that("Different representative selection strategies work", {

  # Strategy: highest mean correlation
  result_mean <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    n_clusters = 5,
    representative_method = "mean_cor",
    verbose = FALSE
  )

  expect_true(length(result_mean$representatives) == 5)

  # Strategy: highest variance
  result_var <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    n_clusters = 5,
    representative_method = "highest_var",
    verbose = FALSE
  )

  expect_true(length(result_var$representatives) == 5)

  # Results might differ
  # (not necessarily, but possible depending on data)
})

# ===== Test 8: Replacement Map Generation =====

test_that("Replacement maps are generated correctly", {

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    n_clusters = 5,
    verbose = FALSE
  )

  # Check replacement map exists
  expect_true("replacement_map" %in% names(result) ||
              "alternatives" %in% names(result))

  # Each non-representative feature should have alternatives
  non_reps <- setdiff(colnames(test_data)[-ncol(test_data)],
                      result$representatives)

  if ("replacement_map" %in% names(result)) {
    for (feat in non_reps) {
      if (feat %in% names(result$replacement_map)) {
        alternatives <- result$replacement_map[[feat]]
        expect_type(alternatives, "character")
        expect_true(length(alternatives) > 0)
      }
    }
  }
})

# ===== Test 9: Dimensionality Reduction Metrics =====

test_that("Dimensionality reduction is measured", {

  n_original <- ncol(test_data) - 1  # Exclude outcome

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    n_clusters = 5,
    verbose = FALSE
  )

  n_representatives <- length(result$representatives)

  # Check metrics
  expect_true("dimensionality_reduction" %in% names(result) ||
              "metrics" %in% names(result))

  # Representatives should be fewer than original features
  expect_true(n_representatives < n_original)

  # Reduction percentage should be calculable
  reduction_pct <- (1 - n_representatives / n_original) * 100
  expect_true(reduction_pct > 0 && reduction_pct < 100)
})

# ===== Test 10: S3 Methods =====

test_that("S3 methods work correctly", {

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    n_clusters = 5,
    verbose = FALSE
  )

  # Test print method
  expect_output(print(result), "OmicSelector")
  expect_output(print(result), "Feature Clustering")

  # Test summary method (if exists)
  if ("summary.OmicSelector_cluster_features" %in% methods("summary")) {
    expect_output(summary(result), "Cluster")
  }

  # Test plot method (if exists)
  if ("plot.OmicSelector_cluster_features" %in% methods("plot")) {
    expect_silent(p <- plot(result))
    expect_true(inherits(p, "ggplot") ||
                inherits(p, "recordedplot") ||
                is.list(p))
  }
})

# ===== Test 11: Distance Metrics =====

test_that("Different distance metrics work", {

  # Pearson correlation
  result_pearson <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    distance_metric = "pearson",
    n_clusters = 5,
    verbose = FALSE
  )

  expect_s3_class(result_pearson, "OmicSelector_cluster_features")

  # Spearman correlation
  result_spearman <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    distance_metric = "spearman",
    n_clusters = 5,
    verbose = FALSE
  )

  expect_s3_class(result_spearman, "OmicSelector_cluster_features")

  # Euclidean distance
  result_euclidean <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "hierarchical",
    distance_metric = "euclidean",
    n_clusters = 5,
    verbose = FALSE
  )

  expect_s3_class(result_euclidean, "OmicSelector_cluster_features")
})

# ===== Test 12: Edge Cases =====

test_that("Handles edge cases gracefully", {

  # Very few features
  small_data <- test_data[, 1:5]
  result_small <- OmicSelector_cluster_features(
    data = small_data,
    method = "correlation",
    n_clusters = 2,
    verbose = FALSE
  )

  expect_s3_class(result_small, "OmicSelector_cluster_features")

  # All features perfectly correlated (within cluster)
  perfect_data <- data.frame(
    feat1 = rnorm(100),
    feat2 = rnorm(100),
    feat3 = rnorm(100)
  )
  perfect_data$feat1_dup <- perfect_data$feat1 + rnorm(100, sd = 0.01)
  perfect_data$feat2_dup <- perfect_data$feat2 + rnorm(100, sd = 0.01)

  result_perfect <- OmicSelector_cluster_features(
    data = perfect_data,
    method = "correlation",
    min_correlation = 0.9,
    verbose = FALSE
  )

  expect_s3_class(result_perfect, "OmicSelector_cluster_features")

  # Should identify the high correlation
  expect_true(length(result_perfect$representatives) < ncol(perfect_data))
})

# ===== Test 13: Integration with Feature Selection =====

test_that("Results can be used with feature selection", {
  skip_if_not_installed("glmnet")

  # First, cluster features
  cluster_result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    n_clusters = 10,
    verbose = FALSE
  )

  # Then, do feature selection on representatives
  rep_data <- test_data[, c(cluster_result$representatives, "outcome")]

  if (requireNamespace("caret", quietly = TRUE)) {
    # Should be able to build a model
    set.seed(123)
    model <- caret::train(
      outcome ~ .,
      data = rep_data,
      method = "glmnet",
      trControl = caret::trainControl(method = "cv", number = 3),
      tuneLength = 1
    )

    expect_s3_class(model, "train")
  }
})

# ===== Test 14: Biomarker Replaceability =====

test_that("Biomarker replacement suggestions work", {

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    min_correlation = 0.7,
    verbose = FALSE
  )

  # Pick a representative
  if (length(result$representatives) > 0) {
    rep_feature <- result$representatives[1]

    # Get cluster members
    cluster_id <- result$cluster_assignments[rep_feature]
    cluster_members <- names(result$cluster_assignments)[
      result$cluster_assignments == cluster_id
    ]

    # There should be multiple members in at least some clusters
    if (length(cluster_members) > 1) {
      # Alternative biomarkers exist
      alternatives <- setdiff(cluster_members, rep_feature)
      expect_true(length(alternatives) > 0)
    }
  }
})

# ===== Test 15: Correlation Matrix Calculation =====

test_that("Correlation matrix is calculated correctly", {

  result <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "correlation",
    verbose = FALSE
  )

  # Correlation matrix should exist
  expect_true("correlation_matrix" %in% names(result))

  cor_mat <- result$correlation_matrix

  # Should be a matrix
  expect_true(is.matrix(cor_mat))

  # Should be square
  expect_equal(nrow(cor_mat), ncol(cor_mat))

  # Diagonal should be 1 (or very close)
  expect_true(all(abs(diag(cor_mat) - 1) < 1e-10))

  # Should be symmetric
  expect_equal(cor_mat, t(cor_mat))

  # Values should be between -1 and 1
  expect_true(all(cor_mat >= -1 & cor_mat <= 1))
})

# ===== Test 16: Reproducibility =====

test_that("Results are reproducible with seed", {

  result1 <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "kmedoids",
    n_clusters = 5,
    seed = 42,
    verbose = FALSE
  )

  result2 <- OmicSelector_cluster_features(
    data = test_data[, -ncol(test_data)],
    method = "kmedoids",
    n_clusters = 5,
    seed = 42,
    verbose = FALSE
  )

  # Should get identical results
  expect_identical(result1$representatives, result2$representatives)
  expect_identical(result1$cluster_assignments, result2$cluster_assignments)
})

# ===== Test 17: Performance with Larger Data =====

test_that("Function handles moderately large datasets", {
  skip_on_cran()  # Performance tests can be slow

  # Create larger dataset
  large_data <- create_correlated_data(
    n_samples = 300,
    n_clusters = 20,
    features_per_cluster = 10
  )

  # Should complete without errors
  expect_error(
    result <- OmicSelector_cluster_features(
      data = large_data[, -ncol(large_data)],
      method = "correlation",
      n_clusters = 20,
      verbose = FALSE
    ),
    NA  # No error expected
  )
})

# ===== Test 18: Missing Value Handling =====

test_that("Handles missing values appropriately", {

  # Create data with missing values
  data_with_na <- test_data
  data_with_na[sample(1:nrow(data_with_na), 10), sample(1:(ncol(data_with_na)-1), 3)] <- NA

  # Should either handle or give informative error
  result <- tryCatch({
    OmicSelector_cluster_features(
      data = data_with_na[, -ncol(data_with_na)],
      method = "correlation",
      n_clusters = 5,
      verbose = FALSE
    )
  }, error = function(e) {
    # If it errors, message should mention missing values
    expect_true(grepl("NA|missing|complete", e$message, ignore.case = TRUE))
    NULL
  })

  # If it succeeded, check result is valid
  if (!is.null(result)) {
    expect_s3_class(result, "OmicSelector_cluster_features")
  }
})
