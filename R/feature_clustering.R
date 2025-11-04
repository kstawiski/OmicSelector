

#' Feature Clustering for Biomarker Replaceability
#'
#' Groups features based on their similarity to identify alternative biomarkers.
#' This is critical for platform transitions (e.g., RNA-seq to qPCR) where
#' not all features may be available or measurable on the new platform.
#'
#' @param data A data frame containing features (columns should be numeric)
#' @param features Character vector of feature names to cluster. If NULL, uses all
#'   numeric columns except those in `exclude_cols`.
#' @param exclude_cols Character vector of column names to exclude (e.g., outcome, IDs)
#' @param n_clusters Integer, desired number of clusters. If NULL, will use
#'   optimal clustering (gap statistic or silhouette)
#' @param method Character string specifying clustering method:
#'   \itemize{
#'     \item "correlation" - Cluster based on Pearson correlation
#'     \item "hierarchical" - Hierarchical clustering
#'     \item "kmeans" - K-means clustering
#'     \item "pam" - Partitioning Around Medoids
#'   }
#' @param distance_metric For hierarchical clustering: "pearson", "spearman", "euclidean"
#' @param min_correlation Numeric (0-1), minimum correlation to consider features
#'   as similar. Only used for correlation method. Default: 0.7
#' @param linkage For hierarchical clustering: "complete", "average", "single", "ward.D2"
#' @param select_representatives Logical, whether to automatically select one
#'   representative feature from each cluster. Default: TRUE
#' @param representative_method How to select representatives:
#'   \itemize{
#'     \item "highest_variance" - Feature with highest variance
#'     \item "most_central" - Feature most correlated with cluster center
#'     \item "median" - Feature closest to cluster median
#'   }
#' @param plot Logical, whether to generate visualization. Default: TRUE
#' @param seed Integer for reproducibility
#'
#' @return A list object of class "OmicSelector_feature_clusters" containing:
#'   \item{clusters}{Named integer vector of cluster assignments}
#'   \item{n_clusters}{Number of clusters}
#'   \item{representatives}{Character vector of representative features (one per cluster)}
#'   \item{cluster_summary}{Data frame summarizing each cluster}
#'   \item{correlation_matrix}{Correlation matrix (if method = "correlation")}
#'   \item{dendrogram}{Dendrogram object (if method = "hierarchical")}
#'   \item{replacements}{List showing which features can replace which}
#'   \item{method}{Method used}
#'
#' @examples
#' \dontrun{
#' # Cluster miRNA features
#' clusters <- OmicSelector_cluster_features(
#'   data = miR_data,
#'   exclude_cols = c("Class", "patient_id"),
#'   method = "correlation",
#'   min_correlation = 0.75,
#'   n_clusters = 20
#' )
#'
#' # View results
#' print(clusters)
#' plot(clusters)
#'
#' # Get representative features
#' representatives <- clusters$representatives
#'
#' # Find replacements for a specific feature
#' replacements <- clusters$replacements[["hsa-miR-21-5p"]]
#' }
#'
#' @export
OmicSelector_cluster_features <- function(
  data,
  features = NULL,
  exclude_cols = NULL,
  n_clusters = NULL,
  method = c("correlation", "hierarchical", "kmeans", "pam"),
  distance_metric = c("pearson", "spearman", "euclidean"),
  min_correlation = 0.7,
  linkage = c("complete", "average", "single", "ward.D2"),
  select_representatives = TRUE,
  representative_method = c("highest_variance", "most_central", "median"),
  plot = TRUE,
  seed = 123
) {

  # Validate inputs
  method <- match.arg(method)
  distance_metric <- match.arg(distance_metric)
  linkage <- match.arg(linkage)
  representative_method <- match.arg(representative_method)

  set.seed(seed)

  # Determine features to cluster
  if (is.null(features)) {
    # Use all numeric columns except excluded ones
    numeric_cols <- names(data)[sapply(data, is.numeric)]
    features <- setdiff(numeric_cols, exclude_cols)
  }

  if (length(features) < 2) {
    stop("Need at least 2 features to cluster")
  }

  # Extract feature data
  feature_data <- as.matrix(data[, features, drop = FALSE])

  # Remove features with zero variance
  feature_vars <- apply(feature_data, 2, var, na.rm = TRUE)
  valid_features <- features[feature_vars > 0 & !is.na(feature_vars)]

  if (length(valid_features) < length(features)) {
    warning(paste0("Removed ", length(features) - length(valid_features),
                   " features with zero or missing variance"))
    features <- valid_features
    feature_data <- feature_data[, features, drop = FALSE]
  }

  message(paste0("Clustering ", length(features), " features using method: ", method))

  # Perform clustering based on method
  if (method == "correlation") {
    result <- .cluster_by_correlation(
      feature_data = feature_data,
      features = features,
      min_correlation = min_correlation,
      n_clusters = n_clusters,
      distance_metric = distance_metric
    )
  } else if (method == "hierarchical") {
    result <- .cluster_hierarchical(
      feature_data = feature_data,
      features = features,
      n_clusters = n_clusters,
      distance_metric = distance_metric,
      linkage = linkage
    )
  } else if (method == "kmeans") {
    result <- .cluster_kmeans(
      feature_data = feature_data,
      features = features,
      n_clusters = n_clusters
    )
  } else if (method == "pam") {
    result <- .cluster_pam(
      feature_data = feature_data,
      features = features,
      n_clusters = n_clusters,
      distance_metric = distance_metric
    )
  }

  # Select representative features from each cluster
  if (select_representatives) {
    message("Selecting representative features...")
    representatives <- .select_cluster_representatives(
      feature_data = feature_data,
      clusters = result$clusters,
      method = representative_method
    )
    result$representatives <- representatives
  }

  # Create cluster summary
  cluster_summary <- .summarize_clusters(
    feature_data = feature_data,
    clusters = result$clusters,
    representatives = if (select_representatives) result$representatives else NULL
  )
  result$cluster_summary <- cluster_summary

  # Create replacement map (which features can replace which)
  result$replacements <- .create_replacement_map(
    clusters = result$clusters,
    representatives = if (select_representatives) result$representatives else NULL
  )

  # Add metadata
  result$method <- method
  result$n_features <- length(features)
  result$parameters <- list(
    min_correlation = min_correlation,
    distance_metric = distance_metric,
    linkage = linkage,
    representative_method = representative_method
  )

  class(result) <- c("OmicSelector_feature_clusters", "list")

  message(paste0("✓ Created ", result$n_clusters, " feature clusters"))

  # Plot if requested
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    print(plot(result))
  }

  return(result)
}


#' Internal: Cluster by Correlation
#'
#' @keywords internal
#' @noRd
.cluster_by_correlation <- function(
  feature_data, features, min_correlation, n_clusters, distance_metric
) {

  # Calculate correlation matrix
  if (distance_metric == "pearson") {
    cor_matrix <- cor(feature_data, use = "pairwise.complete.obs")
  } else if (distance_metric == "spearman") {
    cor_matrix <- cor(feature_data, use = "pairwise.complete.obs", method = "spearman")
  } else {
    cor_matrix <- cor(feature_data, use = "pairwise.complete.obs")
  }

  # Convert correlation to distance
  dist_matrix <- as.dist(1 - abs(cor_matrix))

  # Hierarchical clustering on correlation distance
  hclust_result <- hclust(dist_matrix, method = "complete")

  # Determine number of clusters
  if (is.null(n_clusters)) {
    # Cut tree at correlation threshold
    n_clusters <- length(unique(cutree(hclust_result, h = 1 - min_correlation)))
    message(paste0("Auto-detected ", n_clusters, " clusters at correlation threshold ",
                   min_correlation))
  }

  # Cut tree
  clusters <- cutree(hclust_result, k = n_clusters)
  names(clusters) <- features

  result <- list(
    clusters = clusters,
    n_clusters = n_clusters,
    correlation_matrix = cor_matrix,
    dendrogram = hclust_result,
    distance_matrix = dist_matrix
  )

  return(result)
}


#' Internal: Hierarchical Clustering
#'
#' @keywords internal
#' @noRd
.cluster_hierarchical <- function(
  feature_data, features, n_clusters, distance_metric, linkage
) {

  # Calculate distance
  if (distance_metric == "pearson") {
    cor_matrix <- cor(feature_data, use = "pairwise.complete.obs")
    dist_matrix <- as.dist(1 - abs(cor_matrix))
  } else if (distance_metric == "spearman") {
    cor_matrix <- cor(feature_data, use = "pairwise.complete.obs", method = "spearman")
    dist_matrix <- as.dist(1 - abs(cor_matrix))
  } else {
    dist_matrix <- dist(t(feature_data), method = distance_metric)
    cor_matrix <- NULL
  }

  # Hierarchical clustering
  hclust_result <- hclust(dist_matrix, method = linkage)

  # Determine number of clusters if not specified
  if (is.null(n_clusters)) {
    # Use silhouette or gap statistic
    if (requireNamespace("cluster", quietly = TRUE)) {
      library(cluster)
      # Try different numbers of clusters
      max_k <- min(15, floor(length(features) / 2))
      sil_scores <- sapply(2:max_k, function(k) {
        clusters_temp <- cutree(hclust_result, k = k)
        sil <- silhouette(clusters_temp, dist_matrix)
        mean(sil[, 3])
      })
      n_clusters <- which.max(sil_scores) + 1
      message(paste0("Auto-detected ", n_clusters, " clusters using silhouette method"))
    } else {
      n_clusters <- max(2, floor(sqrt(length(features) / 2)))
      message(paste0("Using default: ", n_clusters, " clusters"))
    }
  }

  # Cut tree
  clusters <- cutree(hclust_result, k = n_clusters)
  names(clusters) <- features

  result <- list(
    clusters = clusters,
    n_clusters = n_clusters,
    dendrogram = hclust_result,
    distance_matrix = dist_matrix
  )

  if (!is.null(cor_matrix)) {
    result$correlation_matrix <- cor_matrix
  }

  return(result)
}


#' Internal: K-means Clustering
#'
#' @keywords internal
#' @noRd
.cluster_kmeans <- function(feature_data, features, n_clusters) {

  # Transpose so features are rows
  feature_data_t <- t(feature_data)

  # Determine number of clusters if not specified
  if (is.null(n_clusters)) {
    # Use elbow method
    max_k <- min(15, floor(nrow(feature_data_t) / 2))
    wss <- sapply(1:max_k, function(k) {
      kmeans(feature_data_t, centers = k, nstart = 10)$tot.withinss
    })

    # Find elbow (simple method)
    diffs <- diff(wss)
    n_clusters <- which.min(abs(diff(diffs))) + 1
    message(paste0("Auto-detected ", n_clusters, " clusters using elbow method"))
  }

  # K-means clustering
  km_result <- kmeans(feature_data_t, centers = n_clusters, nstart = 25)

  clusters <- km_result$cluster
  names(clusters) <- features

  result <- list(
    clusters = clusters,
    n_clusters = n_clusters,
    centers = km_result$centers,
    kmeans_result = km_result
  )

  return(result)
}


#' Internal: PAM Clustering
#'
#' @keywords internal
#' @noRd
.cluster_pam <- function(feature_data, features, n_clusters, distance_metric) {

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("cluster package required for PAM. Install with: install.packages('cluster')")
  }

  library(cluster)

  # Calculate distance
  if (distance_metric == "pearson") {
    cor_matrix <- cor(feature_data, use = "pairwise.complete.obs")
    dist_matrix <- as.dist(1 - abs(cor_matrix))
  } else if (distance_metric == "spearman") {
    cor_matrix <- cor(feature_data, use = "pairwise.complete.obs", method = "spearman")
    dist_matrix <- as.dist(1 - abs(cor_matrix))
  } else {
    dist_matrix <- dist(t(feature_data), method = distance_metric)
  }

  # Determine number of clusters if not specified
  if (is.null(n_clusters)) {
    max_k <- min(15, floor(length(features) / 2))
    sil_scores <- sapply(2:max_k, function(k) {
      pam_temp <- pam(dist_matrix, k = k)
      pam_temp$silinfo$avg.width
    })
    n_clusters <- which.max(sil_scores) + 1
    message(paste0("Auto-detected ", n_clusters, " clusters using silhouette"))
  }

  # PAM clustering
  pam_result <- pam(dist_matrix, k = n_clusters)

  clusters <- pam_result$clustering
  names(clusters) <- features

  result <- list(
    clusters = clusters,
    n_clusters = n_clusters,
    medoids = features[pam_result$medoids],
    pam_result = pam_result
  )

  return(result)
}


#' Internal: Select Cluster Representatives
#'
#' @keywords internal
#' @noRd
.select_cluster_representatives <- function(feature_data, clusters, method) {

  unique_clusters <- sort(unique(clusters))
  representatives <- character(length(unique_clusters))
  names(representatives) <- paste0("Cluster_", unique_clusters)

  for (i in seq_along(unique_clusters)) {
    cluster_id <- unique_clusters[i]
    cluster_features <- names(clusters)[clusters == cluster_id]

    if (length(cluster_features) == 1) {
      representatives[i] <- cluster_features
    } else {
      cluster_data <- feature_data[, cluster_features, drop = FALSE]

      if (method == "highest_variance") {
        vars <- apply(cluster_data, 2, var, na.rm = TRUE)
        representatives[i] <- names(which.max(vars))
      } else if (method == "most_central") {
        # Feature most correlated with cluster mean
        cluster_mean <- rowMeans(cluster_data, na.rm = TRUE)
        cors <- apply(cluster_data, 2, function(x) cor(x, cluster_mean, use = "complete.obs"))
        representatives[i] <- names(which.max(cors))
      } else if (method == "median") {
        # Feature closest to median profile
        cluster_median <- apply(cluster_data, 1, median, na.rm = TRUE)
        dists <- apply(cluster_data, 2, function(x) {
          sqrt(mean((x - cluster_median)^2, na.rm = TRUE))
        })
        representatives[i] <- names(which.min(dists))
      }
    }
  }

  return(representatives)
}


#' Internal: Summarize Clusters
#'
#' @keywords internal
#' @noRd
.summarize_clusters <- function(feature_data, clusters, representatives) {

  unique_clusters <- sort(unique(clusters))

  summary_df <- data.frame(
    Cluster = unique_clusters,
    N_Features = sapply(unique_clusters, function(c) sum(clusters == c)),
    Representative = if (!is.null(representatives)) representatives else NA,
    stringsAsFactors = FALSE
  )

  # Add cluster members
  summary_df$Members <- sapply(unique_clusters, function(c) {
    paste(names(clusters)[clusters == c], collapse = ", ")
  })

  return(summary_df)
}


#' Internal: Create Replacement Map
#'
#' @keywords internal
#' @noRd
.create_replacement_map <- function(clusters, representatives) {

  replacement_map <- list()

  for (feature in names(clusters)) {
    cluster_id <- clusters[feature]
    # All features in same cluster can replace each other
    cluster_members <- names(clusters)[clusters == cluster_id]
    # Exclude the feature itself
    replacements <- setdiff(cluster_members, feature)

    # Prioritize representative if available
    if (!is.null(representatives)) {
      rep_feature <- representatives[paste0("Cluster_", cluster_id)]
      if (rep_feature %in% replacements) {
        # Move representative to front
        replacements <- c(rep_feature, setdiff(replacements, rep_feature))
      }
    }

    replacement_map[[feature]] <- replacements
  }

  return(replacement_map)
}


#' Print Method for OmicSelector_feature_clusters
#'
#' @param x An OmicSelector_feature_clusters object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_feature_clusters <- function(x, ...) {
  cat("OmicSelector Feature Clustering\n")
  cat("================================\n\n")

  cat("Method:", x$method, "\n")
  cat("Features clustered:", x$n_features, "\n")
  cat("Number of clusters:", x$n_clusters, "\n\n")

  cat("Cluster Summary:\n")
  print(x$cluster_summary[, c("Cluster", "N_Features", "Representative")])

  cat("\n")
  cat("Largest clusters:\n")
  top_5 <- head(x$cluster_summary[order(x$cluster_summary$N_Features, decreasing = TRUE), ], 5)
  for (i in 1:nrow(top_5)) {
    cat(sprintf("  Cluster %d: %d features (Rep: %s)\n",
                top_5$Cluster[i], top_5$N_Features[i], top_5$Representative[i]))
  }

  cat("\nUse summary() for detailed cluster membership.\n")
}


#' Summary Method for OmicSelector_feature_clusters
#'
#' @param object An OmicSelector_feature_clusters object
#' @param ... Additional arguments (not used)
#' @export
summary.OmicSelector_feature_clusters <- function(object, ...) {
  print(object)

  cat("\n")
  cat("Detailed Cluster Membership:\n")
  cat("============================\n\n")

  for (i in 1:min(object$n_clusters, 10)) {
    cluster_features <- names(object$clusters)[object$clusters == i]
    cat(sprintf("Cluster %d (%d features):\n", i, length(cluster_features)))
    cat("  ", paste(cluster_features, collapse=", "), "\n\n")
  }

  if (object$n_clusters > 10) {
    cat(sprintf("... and %d more clusters\n", object$n_clusters - 10))
  }
}


#' Plot Method for OmicSelector_feature_clusters
#'
#' @param x An OmicSelector_feature_clusters object
#' @param ... Additional arguments (not used)
#' @export
plot.OmicSelector_feature_clusters <- function(x, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }

  library(ggplot2)

  # Plot cluster sizes
  p1 <- ggplot(x$cluster_summary, aes(x = reorder(factor(Cluster), -N_Features),
                                       y = N_Features)) +
    geom_col(fill = "steelblue") +
    labs(title = "Feature Cluster Sizes",
         subtitle = paste0("Method: ", x$method, ", Total: ", x$n_clusters, " clusters"),
         x = "Cluster", y = "Number of Features") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p1)

  # If dendrogram available, could plot it
  if (!is.null(x$dendrogram) && requireNamespace("ggdendro", quietly = TRUE)) {
    library(ggdendro)
    p2 <- ggdendrogram(x$dendrogram) +
      labs(title = "Hierarchical Clustering Dendrogram")
    print(p2)
  }
}
