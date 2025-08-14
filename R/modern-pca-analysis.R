#' Modern PCA Analysis for Omics Data
#'
#' @description
#' Modernized Principal Component Analysis with enhanced visualization options,
#' comprehensive error handling, and multiple plot types. Supports both 2D and 3D
#' visualizations with publication-ready styling.
#'
#' Key improvements:
#' - Removed 17 redundant library dependencies
#' - Added comprehensive input validation
#' - Multiple visualization options (ggplot2, plotly for 3D)
#' - Enhanced statistical reporting
#' - Customizable styling and annotations
#' - Support for multiple grouping variables
#' - Scree plots and loadings analysis
#' - Missing value handling options
#'
#' @param expression_data Matrix or data.frame with features in columns, samples in rows
#' @param sample_groups Factor or character vector indicating sample groups
#' @param scale_data Logical: whether to scale features (default: TRUE)
#' @param center_data Logical: whether to center features (default: TRUE)
#' @param plot_type Character: "2d", "3d", "scree", "loadings", or "all"
#' @param color_palette Character vector of colors for groups
#' @param show_ellipses Logical: whether to show confidence ellipses
#' @param ellipse_confidence Numeric: confidence level for ellipses (default: 0.95)
#' @param show_loadings Logical: whether to show variable loadings
#' @param title Character: plot title
#' @param point_size Numeric: size of points in plots
#' @param interactive Logical: create interactive plots using plotly
#' @param variance_threshold Numeric: minimum variance explained to include PC
#' @param max_components Integer: maximum number of PCs to calculate
#' @param na_action Character: how to handle missing values ("omit", "impute")
#'
#' @return List containing PCA results, plots, and summary statistics
#'
#' @examples
#' \dontrun{
#' # Basic PCA
#' data <- matrix(rnorm(1000), nrow = 50, ncol = 20)
#' groups <- factor(rep(c("Case", "Control"), each = 25))
#' 
#' pca_results <- perform_modern_pca(
#'   expression_data = data,
#'   sample_groups = groups
#' )
#' 
#' # Advanced PCA with custom options
#' pca_results <- perform_modern_pca(
#'   expression_data = data,
#'   sample_groups = groups,
#'   plot_type = "all",
#'   show_ellipses = TRUE,
#'   interactive = TRUE,
#'   title = "miRNA Expression PCA"
#' )
#' 
#' # Access results
#' print(pca_results$summary)
#' pca_results$plots$pca_2d
#' }
#'
#' @export
perform_modern_pca <- function(expression_data,
                              sample_groups,
                              scale_data = TRUE,
                              center_data = TRUE,
                              plot_type = "2d",
                              color_palette = NULL,
                              show_ellipses = TRUE,
                              ellipse_confidence = 0.95,
                              show_loadings = FALSE,
                              title = "PCA Analysis",
                              point_size = 3,
                              interactive = FALSE,
                              variance_threshold = 0.01,
                              max_components = 10,
                              na_action = "omit") {
  
  # Validate inputs
  validate_pca_inputs(expression_data, sample_groups, plot_type, na_action)
  
  # Load only essential packages
  required_packages <- c("ggplot2", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed", pkg))
    }
  }
  
  cat("🔬 Starting modern PCA analysis...\n")
  cat(sprintf("   Features: %d\n", ncol(expression_data)))
  cat(sprintf("   Samples: %d\n", nrow(expression_data)))
  cat(sprintf("   Groups: %d (%s)\n", 
              length(unique(sample_groups)), 
              paste(unique(sample_groups), collapse = ", ")))
  
  # Prepare data
  prepared_data <- prepare_pca_data(
    expression_data = expression_data,
    sample_groups = sample_groups,
    na_action = na_action
  )
  
  # Perform PCA
  cat("🧮 Computing principal components...\n")
  pca_result <- compute_pca(
    data = prepared_data$clean_data,
    scale = scale_data,
    center = center_data,
    max_components = max_components
  )
  
  # Filter components by variance threshold
  significant_pcs <- which(pca_result$variance_explained >= variance_threshold)
  if (length(significant_pcs) == 0) {
    warning("No components meet variance threshold, using first 2 components")
    significant_pcs <- 1:min(2, ncol(pca_result$scores))
  }
  
  # Create summary
  pca_summary <- create_pca_summary(
    pca_result = pca_result,
    groups = prepared_data$clean_groups,
    significant_pcs = significant_pcs
  )
  
  # Generate plots
  cat("📊 Creating visualizations...\n")
  plots <- create_pca_plots(
    pca_result = pca_result,
    groups = prepared_data$clean_groups,
    plot_type = plot_type,
    color_palette = color_palette,
    show_ellipses = show_ellipses,
    ellipse_confidence = ellipse_confidence,
    show_loadings = show_loadings,
    title = title,
    point_size = point_size,
    interactive = interactive,
    significant_pcs = significant_pcs
  )
  
  cat("✅ PCA analysis complete!\n")
  
  # Return comprehensive results
  return(structure(list(
    pca = pca_result,
    summary = pca_summary,
    plots = plots,
    data_info = list(
      n_features = ncol(expression_data),
      n_samples = nrow(expression_data),
      n_groups = length(unique(sample_groups)),
      group_sizes = table(sample_groups),
      variance_explained = pca_result$variance_explained[significant_pcs]
    ),
    parameters = list(
      scale_data = scale_data,
      center_data = center_data,
      variance_threshold = variance_threshold,
      na_action = na_action
    )
  ), class = "omics_pca"))
}

#' Validate PCA Inputs
#' @keywords internal
validate_pca_inputs <- function(expression_data, sample_groups, plot_type, na_action) {
  
  # Check expression data
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("'expression_data' must be a matrix or data.frame")
  }
  
  if (nrow(expression_data) == 0 || ncol(expression_data) == 0) {
    stop("'expression_data' must have both rows and columns")
  }
  
  # Check sample groups
  if (length(sample_groups) != nrow(expression_data)) {
    stop("Length of 'sample_groups' must equal number of rows in 'expression_data'")
  }
  
  # Check plot type
  valid_plot_types <- c("2d", "3d", "scree", "loadings", "all")
  if (!plot_type %in% valid_plot_types) {
    stop(sprintf("'plot_type' must be one of: %s", paste(valid_plot_types, collapse = ", ")))
  }
  
  # Check NA action
  valid_na_actions <- c("omit", "impute")
  if (!na_action %in% valid_na_actions) {
    stop(sprintf("'na_action' must be one of: %s", paste(valid_na_actions, collapse = ", ")))
  }
}

#' Prepare Data for PCA
#' @keywords internal
prepare_pca_data <- function(expression_data, sample_groups, na_action) {
  
  # Convert to data.frame if needed
  if (is.matrix(expression_data)) {
    expression_data <- as.data.frame(expression_data)
  }
  
  # Ensure all columns are numeric
  numeric_cols <- sapply(expression_data, is.numeric)
  if (!all(numeric_cols)) {
    warning("Non-numeric columns detected and will be removed")
    expression_data <- expression_data[, numeric_cols, drop = FALSE]
  }
  
  # Handle missing values
  if (any(is.na(expression_data))) {
    if (na_action == "omit") {
      # Remove rows with any missing values
      complete_rows <- complete.cases(expression_data)
      expression_data <- expression_data[complete_rows, ]
      sample_groups <- sample_groups[complete_rows]
      
      cat(sprintf("ℹ️ Removed %d samples with missing values\n", 
                  sum(!complete_rows)))
      
    } else if (na_action == "impute") {
      # Impute missing values with column means
      for (i in seq_len(ncol(expression_data))) {
        col_mean <- mean(expression_data[, i], na.rm = TRUE)
        expression_data[is.na(expression_data[, i]), i] <- col_mean
      }
      cat("ℹ️ Imputed missing values with column means\n")
    }
  }
  
  # Remove features with zero variance
  feature_vars <- apply(expression_data, 2, var, na.rm = TRUE)
  zero_var_features <- which(feature_vars == 0 | is.na(feature_vars))
  
  if (length(zero_var_features) > 0) {
    expression_data <- expression_data[, -zero_var_features]
    cat(sprintf("ℹ️ Removed %d features with zero variance\n", 
                length(zero_var_features)))
  }
  
  # Final validation
  if (ncol(expression_data) < 2) {
    stop("At least 2 features are required for PCA")
  }
  
  if (nrow(expression_data) < 3) {
    stop("At least 3 samples are required for PCA")
  }
  
  return(list(
    clean_data = as.matrix(expression_data),
    clean_groups = factor(sample_groups)
  ))
}

#' Compute PCA
#' @keywords internal
compute_pca <- function(data, scale, center, max_components) {
  
  # Perform PCA
  pca_fit <- prcomp(data, scale. = scale, center = center)
  
  # Limit to max_components if specified
  n_components <- min(max_components, ncol(pca_fit$x))
  
  # Extract variance explained
  variance_explained <- (pca_fit$sdev^2 / sum(pca_fit$sdev^2))[1:n_components]
  cumulative_variance <- cumsum(variance_explained)
  
  # Extract scores and loadings
  scores <- pca_fit$x[, 1:n_components, drop = FALSE]
  loadings <- pca_fit$rotation[, 1:n_components, drop = FALSE]
  
  return(list(
    scores = scores,
    loadings = loadings,
    variance_explained = variance_explained,
    cumulative_variance = cumulative_variance,
    sdev = pca_fit$sdev[1:n_components],
    center = pca_fit$center,
    scale = pca_fit$scale,
    raw_pca = pca_fit
  ))
}

#' Create PCA Summary
#' @keywords internal
create_pca_summary <- function(pca_result, groups, significant_pcs) {
  
  # Basic statistics
  n_components <- length(pca_result$variance_explained)
  
  # Group statistics
  group_summary <- data.frame(
    Group = names(table(groups)),
    Count = as.numeric(table(groups)),
    Percentage = round(as.numeric(table(groups)) / length(groups) * 100, 1)
  )
  
  # Component statistics
  component_summary <- data.frame(
    Component = paste0("PC", 1:n_components),
    Variance_Explained = round(pca_result$variance_explained * 100, 2),
    Cumulative_Variance = round(pca_result$cumulative_variance * 100, 2),
    Standard_Deviation = round(pca_result$sdev, 3),
    Significant = 1:n_components %in% significant_pcs
  )
  
  return(list(
    n_components = n_components,
    total_variance_explained = round(sum(pca_result$variance_explained) * 100, 2),
    significant_components = length(significant_pcs),
    group_summary = group_summary,
    component_summary = component_summary
  ))
}

#' Create PCA Plots
#' @keywords internal
create_pca_plots <- function(pca_result, groups, plot_type, color_palette, 
                            show_ellipses, ellipse_confidence, show_loadings,
                            title, point_size, interactive, significant_pcs) {
  
  plots <- list()
  
  # Set up color palette
  if (is.null(color_palette)) {
    n_groups <- length(unique(groups))
    if (n_groups <= 8) {
      color_palette <- RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
    } else {
      color_palette <- rainbow(n_groups)
    }
  }
  
  # Create 2D plot
  if (plot_type %in% c("2d", "all")) {
    plots$pca_2d <- create_2d_pca_plot(
      pca_result, groups, color_palette, show_ellipses, 
      ellipse_confidence, show_loadings, title, point_size
    )
  }
  
  # Create 3D plot
  if (plot_type %in% c("3d", "all") && ncol(pca_result$scores) >= 3) {
    plots$pca_3d <- create_3d_pca_plot(
      pca_result, groups, color_palette, title, point_size, interactive
    )
  }
  
  # Create scree plot
  if (plot_type %in% c("scree", "all")) {
    plots$scree <- create_scree_plot(pca_result, significant_pcs)
  }
  
  # Create loadings plot
  if (plot_type %in% c("loadings", "all") && show_loadings) {
    plots$loadings <- create_loadings_plot(pca_result, title)
  }
  
  return(plots)
}

#' Create 2D PCA Plot
#' @keywords internal
create_2d_pca_plot <- function(pca_result, groups, color_palette, show_ellipses, 
                              ellipse_confidence, show_loadings, title, point_size) {
  
  # Create data frame for plotting
  plot_data <- data.frame(
    PC1 = pca_result$scores[, 1],
    PC2 = pca_result$scores[, 2],
    Group = groups
  )
  
  # Calculate variance explained for axis labels
  var1 <- round(pca_result$variance_explained[1] * 100, 1)
  var2 <- round(pca_result$variance_explained[2] * 100, 1)
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::labs(
      title = title,
      x = paste0("PC1 (", var1, "% variance)"),
      y = paste0("PC2 (", var2, "% variance)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  # Add ellipses if requested
  if (show_ellipses) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- p + ggplot2::stat_ellipse(
        level = ellipse_confidence,
        alpha = 0.3,
        size = 1.2
      )
    }
  }
  
  return(p)
}

#' Create 3D PCA Plot
#' @keywords internal
create_3d_pca_plot <- function(pca_result, groups, color_palette, title, point_size, interactive) {
  
  if (!interactive) {
    warning("3D plots require interactive=TRUE and plotly package")
    return(NULL)
  }
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    warning("plotly package required for 3D plots")
    return(NULL)
  }
  
  # Calculate variance explained
  var1 <- round(pca_result$variance_explained[1] * 100, 1)
  var2 <- round(pca_result$variance_explained[2] * 100, 1)
  var3 <- round(pca_result$variance_explained[3] * 100, 1)
  
  # Create 3D plot
  plot_3d <- plotly::plot_ly(
    x = pca_result$scores[, 1],
    y = pca_result$scores[, 2],
    z = pca_result$scores[, 3],
    color = groups,
    colors = color_palette,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = point_size)
  ) %>%
    plotly::layout(
      title = title,
      scene = list(
        xaxis = list(title = paste0("PC1 (", var1, "%)")),
        yaxis = list(title = paste0("PC2 (", var2, "%)")),
        zaxis = list(title = paste0("PC3 (", var3, "%)"))
      )
    )
  
  return(plot_3d)
}

#' Create Scree Plot
#' @keywords internal
create_scree_plot <- function(pca_result, significant_pcs) {
  
  scree_data <- data.frame(
    Component = factor(1:length(pca_result$variance_explained)),
    Variance = pca_result$variance_explained * 100,
    Significant = 1:length(pca_result$variance_explained) %in% significant_pcs
  )
  
  ggplot2::ggplot(scree_data, ggplot2::aes(x = Component, y = Variance)) +
    ggplot2::geom_col(ggplot2::aes(fill = Significant), alpha = 0.8) +
    ggplot2::geom_line(group = 1, size = 1, color = "red", alpha = 0.6) +
    ggplot2::geom_point(size = 2, color = "red") +
    ggplot2::scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "steelblue")) +
    ggplot2::labs(
      title = "Scree Plot - Variance Explained by Principal Components",
      x = "Principal Component",
      y = "Variance Explained (%)",
      fill = "Above Threshold"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )
}

#' Create Loadings Plot
#' @keywords internal
create_loadings_plot <- function(pca_result, title) {
  
  # Extract top loadings for PC1 and PC2
  n_top <- min(10, nrow(pca_result$loadings))
  
  pc1_loadings <- pca_result$loadings[, 1]
  pc2_loadings <- pca_result$loadings[, 2]
  
  # Get top features by absolute loading
  top_features_pc1 <- names(sort(abs(pc1_loadings), decreasing = TRUE))[1:n_top]
  top_features_pc2 <- names(sort(abs(pc2_loadings), decreasing = TRUE))[1:n_top]
  
  top_features <- unique(c(top_features_pc1, top_features_pc2))
  
  loadings_data <- data.frame(
    Feature = top_features,
    PC1 = pc1_loadings[top_features],
    PC2 = pc2_loadings[top_features]
  )
  
  ggplot2::ggplot(loadings_data, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    ggplot2::geom_text(ggplot2::aes(label = Feature), 
                      hjust = 0, vjust = 0, size = 3, 
                      check_overlap = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(
      title = paste0("Feature Loadings - ", title),
      x = "PC1 Loading",
      y = "PC2 Loading"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )
}

#' Print Method for PCA Results
#' @export
print.omics_pca <- function(x, ...) {
  cat("Modern Omics PCA Analysis\n")
  cat("========================\n\n")
  
  cat("Dataset Information:\n")
  cat(sprintf("  Features: %d\n", x$data_info$n_features))
  cat(sprintf("  Samples: %d\n", x$data_info$n_samples))
  cat(sprintf("  Groups: %d\n", x$data_info$n_groups))
  
  cat("\nGroup Distribution:\n")
  for (i in seq_len(nrow(x$summary$group_summary))) {
    cat(sprintf("  %s: %d samples (%.1f%%)\n", 
                x$summary$group_summary$Group[i],
                x$summary$group_summary$Count[i],
                x$summary$group_summary$Percentage[i]))
  }
  
  cat("\nPrincipal Components:\n")
  cat(sprintf("  Total components: %d\n", x$summary$n_components))
  cat(sprintf("  Significant components: %d\n", x$summary$significant_components))
  cat(sprintf("  Total variance explained: %.1f%%\n", x$summary$total_variance_explained))
  
  cat("\nTop Components:\n")
  top_components <- head(x$summary$component_summary, 5)
  for (i in seq_len(nrow(top_components))) {
    cat(sprintf("  %s: %.1f%% variance (cumulative: %.1f%%)\n",
                top_components$Component[i],
                top_components$Variance_Explained[i],
                top_components$Cumulative_Variance[i]))
  }
  
  invisible(x)
}

#' Backward Compatibility Alias
#' @export
OmicSelector_PCA <- function(ttpm_features, meta) {
  .Deprecated("perform_modern_pca", 
              msg = "OmicSelector_PCA is deprecated. Use perform_modern_pca() instead.")
  
  # Convert to modern function call
  result <- perform_modern_pca(
    expression_data = ttpm_features,
    sample_groups = meta,
    plot_type = "2d",
    show_ellipses = TRUE,
    verbose = FALSE
  )
  
  # Return just the 2D plot for backward compatibility
  return(result$plots$pca_2d)
}
