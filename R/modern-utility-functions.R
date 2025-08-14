#' Modern Utility Functions for OmicSelector
#'
#' @description
#' Collection of modernized utility functions with improved error handling,
#' performance optimizations, and comprehensive validation.
#'

#' Modern Formula Merging Function
#'
#' @description
#' Intelligently merge multiple formulas while handling duplicates and validation.
#' Improved version of OmicSelector_merge_formulas with better error handling.
#'
#' @param formula_list List of formula objects to merge
#' @param class_column Character: name of the class/response column (default: "Class")
#' @param remove_duplicates Logical: remove duplicate features (default: TRUE)
#' @param max_features Integer: maximum number of features to include (default: NULL)
#' @param priority_method Character: how to prioritize features when max_features is set
#' @param validate_formulas Logical: validate input formulas (default: TRUE)
#'
#' @return Single formula object with merged features
#'
#' @examples
#' \dontrun{
#' # Basic merging
#' f1 <- Class ~ feature1 + feature2
#' f2 <- Class ~ feature2 + feature3
#' merged <- merge_omics_formulas(list(f1, f2))
#' 
#' # With feature limit
#' merged <- merge_omics_formulas(
#'   list(f1, f2), 
#'   max_features = 5,
#'   priority_method = "frequency"
#' )
#' }
#'
#' @export
merge_omics_formulas <- function(formula_list,
                                class_column = "Class",
                                remove_duplicates = TRUE,
                                max_features = NULL,
                                priority_method = "frequency",
                                validate_formulas = TRUE) {
  
  # Validate inputs
  if (!is.list(formula_list) || length(formula_list) == 0) {
    stop("'formula_list' must be a non-empty list")
  }
  
  if (validate_formulas) {
    # Check that all elements are formulas
    is_formula <- sapply(formula_list, function(x) inherits(x, "formula"))
    if (!all(is_formula)) {
      stop("All elements in 'formula_list' must be formula objects")
    }
  }
  
  # Extract features from all formulas
  all_features <- character(0)
  feature_counts <- list()
  
  for (i in seq_along(formula_list)) {
    formula_obj <- formula_list[[i]]
    
    # Parse formula to extract features
    formula_terms <- attr(terms(formula_obj), "term.labels")
    
    # Clean feature names (remove backticks if present)
    clean_features <- gsub("`", "", formula_terms)
    
    all_features <- c(all_features, clean_features)
    
    # Count feature occurrences for prioritization
    for (feature in clean_features) {
      feature_counts[[feature]] <- (feature_counts[[feature]] %||% 0) + 1
    }
  }
  
  # Handle empty feature set
  if (length(all_features) == 0) {
    warning("No features found in any formula")
    return(as.formula(paste(class_column, "~ 1")))
  }
  
  # Remove duplicates if requested
  if (remove_duplicates) {
    unique_features <- unique(all_features)
  } else {
    unique_features <- all_features
  }
  
  # Apply feature limit with prioritization
  if (!is.null(max_features) && length(unique_features) > max_features) {
    
    if (priority_method == "frequency") {
      # Prioritize by frequency of occurrence across formulas
      feature_freq <- sapply(unique_features, function(f) feature_counts[[f]] %||% 1)
      priority_order <- order(feature_freq, decreasing = TRUE)
      
    } else if (priority_method == "alphabetical") {
      # Alphabetical order
      priority_order <- order(unique_features)
      
    } else if (priority_method == "first") {
      # Order of first appearance
      priority_order <- seq_along(unique_features)
      
    } else {
      stop("Invalid priority_method. Must be 'frequency', 'alphabetical', or 'first'")
    }
    
    selected_features <- unique_features[priority_order][1:max_features]
    
  } else {
    selected_features <- unique_features
  }
  
  # Create merged formula
  create_omics_formula(selected_features, class_column = class_column)
}

#' Modern Signature Overlap Analysis
#'
#' @description
#' Analyze overlap between different feature signatures with comprehensive statistics
#' and visualization options. Modernized version of OmicSelector_signiture_overlap.
#'
#' @param signature_list Named list of character vectors representing signatures
#' @param method Character: analysis method ("jaccard", "overlap", "dice", "all")
#' @param min_overlap Integer: minimum overlap to report (default: 1)
#' @param create_plot Logical: create visualization (default: TRUE)
#' @param plot_type Character: "heatmap", "network", "upset", or "all"
#' @param return_details Logical: return detailed overlap information
#'
#' @return Data frame with overlap statistics or list with plots and statistics
#'
#' @examples
#' \dontrun{
#' # Basic overlap analysis
#' sig1 <- c("hsa_miR_1", "hsa_miR_2", "hsa_miR_3")
#' sig2 <- c("hsa_miR_2", "hsa_miR_3", "hsa_miR_4")
#' sig3 <- c("hsa_miR_1", "hsa_miR_4", "hsa_miR_5")
#' 
#' signatures <- list("Method1" = sig1, "Method2" = sig2, "Method3" = sig3)
#' overlap_results <- analyze_signature_overlap(signatures)
#' }
#'
#' @export
analyze_signature_overlap <- function(signature_list,
                                     method = "all",
                                     min_overlap = 1,
                                     create_plot = TRUE,
                                     plot_type = "heatmap",
                                     return_details = TRUE) {
  
  # Validate inputs
  if (!is.list(signature_list) || length(signature_list) < 2) {
    stop("'signature_list' must be a list with at least 2 signatures")
  }
  
  if (is.null(names(signature_list))) {
    names(signature_list) <- paste0("Signature_", seq_along(signature_list))
  }
  
  # Ensure all signatures are character vectors
  signature_list <- lapply(signature_list, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  # Calculate pairwise overlap statistics
  n_signatures <- length(signature_list)
  signature_names <- names(signature_list)
  
  # Initialize results matrices
  overlap_matrix <- matrix(0, nrow = n_signatures, ncol = n_signatures,
                          dimnames = list(signature_names, signature_names))
  jaccard_matrix <- overlap_matrix
  dice_matrix <- overlap_matrix
  
  detailed_results <- list()
  
  for (i in seq_len(n_signatures)) {
    for (j in seq_len(n_signatures)) {
      
      sig_i <- signature_list[[i]]
      sig_j <- signature_list[[j]]
      
      if (i == j) {
        # Self-comparison
        overlap_count <- length(sig_i)
        jaccard_coef <- 1.0
        dice_coef <- 1.0
        
      } else {
        # Calculate overlap
        intersection <- intersect(sig_i, sig_j)
        union_set <- union(sig_i, sig_j)
        
        overlap_count <- length(intersection)
        
        # Jaccard coefficient
        jaccard_coef <- if (length(union_set) > 0) {
          length(intersection) / length(union_set)
        } else {
          0
        }
        
        # Dice coefficient
        dice_coef <- if (length(sig_i) + length(sig_j) > 0) {
          2 * length(intersection) / (length(sig_i) + length(sig_j))
        } else {
          0
        }
        
        # Store detailed results
        if (overlap_count >= min_overlap) {
          detailed_results[[paste(signature_names[i], signature_names[j], sep = "_vs_")]] <- list(
            signature1 = signature_names[i],
            signature2 = signature_names[j],
            size1 = length(sig_i),
            size2 = length(sig_j),
            intersection = intersection,
            overlap_count = overlap_count,
            jaccard = jaccard_coef,
            dice = dice_coef
          )
        }
      }
      
      # Fill matrices
      overlap_matrix[i, j] <- overlap_count
      jaccard_matrix[i, j] <- jaccard_coef
      dice_matrix[i, j] <- dice_coef
    }
  }
  
  # Create summary statistics
  summary_stats <- create_overlap_summary(signature_list, overlap_matrix, jaccard_matrix, dice_matrix)
  
  # Create plots if requested
  plots <- NULL
  if (create_plot) {
    plots <- create_overlap_plots(
      overlap_matrix = overlap_matrix,
      jaccard_matrix = jaccard_matrix,
      dice_matrix = dice_matrix,
      signature_list = signature_list,
      plot_type = plot_type
    )
  }
  
  # Prepare return object
  if (return_details) {
    return(list(
      summary = summary_stats,
      overlap_matrix = overlap_matrix,
      jaccard_matrix = jaccard_matrix,
      dice_matrix = dice_matrix,
      detailed_results = detailed_results,
      plots = plots
    ))
  } else {
    return(summary_stats)
  }
}

#' Create Overlap Summary Statistics
#' @keywords internal
create_overlap_summary <- function(signature_list, overlap_matrix, jaccard_matrix, dice_matrix) {
  
  signature_names <- names(signature_list)
  n_signatures <- length(signature_list)
  
  # Signature size statistics
  signature_sizes <- sapply(signature_list, length)
  
  # Pairwise statistics (excluding diagonal)
  upper_tri_idx <- upper.tri(overlap_matrix)
  
  summary_df <- data.frame(
    Signature = signature_names,
    Size = signature_sizes,
    Mean_Overlap = apply(overlap_matrix, 1, function(x) mean(x[-which.max(x)])),
    Max_Overlap = apply(overlap_matrix, 1, function(x) max(x[-which.max(x)])),
    Mean_Jaccard = apply(jaccard_matrix, 1, function(x) mean(x[-which.max(x)])),
    Max_Jaccard = apply(jaccard_matrix, 1, function(x) max(x[-which.max(x)])),
    stringsAsFactors = FALSE
  )
  
  # Overall statistics
  overall_stats <- list(
    total_signatures = n_signatures,
    mean_signature_size = mean(signature_sizes),
    median_signature_size = median(signature_sizes),
    mean_pairwise_overlap = mean(overlap_matrix[upper_tri_idx]),
    mean_jaccard = mean(jaccard_matrix[upper_tri_idx]),
    mean_dice = mean(dice_matrix[upper_tri_idx]),
    max_pairwise_overlap = max(overlap_matrix[upper_tri_idx]),
    total_unique_features = length(unique(unlist(signature_list)))
  )
  
  return(list(
    signature_summary = summary_df,
    overall_statistics = overall_stats
  ))
}

#' Create Overlap Visualization Plots
#' @keywords internal
create_overlap_plots <- function(overlap_matrix, jaccard_matrix, dice_matrix, 
                                signature_list, plot_type) {
  
  plots <- list()
  
  # Heatmap of overlap
  if (plot_type %in% c("heatmap", "all")) {
    plots$overlap_heatmap <- create_overlap_heatmap(overlap_matrix, "Feature Overlap Count")
    plots$jaccard_heatmap <- create_overlap_heatmap(jaccard_matrix, "Jaccard Coefficient")
    plots$dice_heatmap <- create_overlap_heatmap(dice_matrix, "Dice Coefficient")
  }
  
  # UpSet plot for multi-way intersections
  if (plot_type %in% c("upset", "all")) {
    plots$upset_plot <- create_upset_plot(signature_list)
  }
  
  # Network plot
  if (plot_type %in% c("network", "all")) {
    plots$network_plot <- create_network_plot(jaccard_matrix, threshold = 0.1)
  }
  
  return(plots)
}

#' Create Overlap Heatmap
#' @keywords internal
create_overlap_heatmap <- function(matrix_data, title) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for heatmap visualization")
    return(NULL)
  }
  
  # Convert matrix to long format
  matrix_df <- expand.grid(
    Signature1 = rownames(matrix_data),
    Signature2 = colnames(matrix_data),
    stringsAsFactors = FALSE
  )
  matrix_df$Value <- as.vector(matrix_data)
  
  ggplot2::ggplot(matrix_df, ggplot2::aes(x = Signature1, y = Signature2, fill = Value)) +
    ggplot2::geom_tile(color = "white", size = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(Value, 3)), color = "black", size = 3) +
    ggplot2::scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                                 midpoint = max(matrix_data) / 2) +
    ggplot2::labs(title = title, x = "Signature", y = "Signature", fill = "Value") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

#' Create UpSet Plot for Multi-way Intersections
#' @keywords internal
create_upset_plot <- function(signature_list) {
  
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    warning("UpSetR package required for UpSet plots")
    return(NULL)
  }
  
  # Convert to UpSetR format
  all_features <- unique(unlist(signature_list))
  
  upset_data <- data.frame(
    Feature = all_features,
    stringsAsFactors = FALSE
  )
  
  # Add columns for each signature
  for (sig_name in names(signature_list)) {
    upset_data[[sig_name]] <- as.integer(upset_data$Feature %in% signature_list[[sig_name]])
  }
  
  # Create UpSet plot
  tryCatch({
    UpSetR::upset(upset_data, 
                  sets = names(signature_list),
                  order.by = "freq",
                  mainbar.y.label = "Intersection Size",
                  sets.x.label = "Signature Size")
  }, error = function(e) {
    warning("Could not create UpSet plot: ", e$message)
    return(NULL)
  })
}

#' Create Network Plot
#' @keywords internal
create_network_plot <- function(similarity_matrix, threshold = 0.1) {
  
  if (!requireNamespace("igraph", quietly = TRUE) || 
      !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("igraph and ggplot2 packages required for network plots")
    return(NULL)
  }
  
  # Create adjacency matrix with threshold
  adj_matrix <- similarity_matrix
  adj_matrix[adj_matrix < threshold] <- 0
  diag(adj_matrix) <- 0
  
  # Create igraph object
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, 
                                          mode = "undirected", 
                                          weighted = TRUE)
  
  if (igraph::ecount(g) == 0) {
    warning("No connections above threshold for network plot")
    return(NULL)
  }
  
  # Extract layout
  layout <- igraph::layout_with_fr(g)
  
  # Create edge and vertex data frames
  edges <- igraph::as_data_frame(g, what = "edges")
  vertices <- data.frame(
    name = igraph::V(g)$name,
    x = layout[, 1],
    y = layout[, 2],
    stringsAsFactors = FALSE
  )
  
  # Create ggplot
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, 
                         ggplot2::aes(x = vertices$x[match(from, vertices$name)],
                                     y = vertices$y[match(from, vertices$name)],
                                     xend = vertices$x[match(to, vertices$name)],
                                     yend = vertices$y[match(to, vertices$name)],
                                     alpha = weight),
                         color = "gray50") +
    ggplot2::geom_point(data = vertices, 
                       ggplot2::aes(x = x, y = y), 
                       size = 8, color = "steelblue") +
    ggplot2::geom_text(data = vertices, 
                      ggplot2::aes(x = x, y = y, label = name), 
                      color = "white", size = 3, fontface = "bold") +
    ggplot2::labs(title = "Signature Similarity Network",
                 alpha = "Jaccard\nCoefficient") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# Null-coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Backward Compatibility Aliases
#' @export
OmicSelector_merge_formulas <- function(...) {
  .Deprecated("merge_omics_formulas", 
              msg = "OmicSelector_merge_formulas is deprecated. Use merge_omics_formulas() instead.")
  
  # Convert old call to new function
  formulas <- list(...)
  merge_omics_formulas(formulas, validate_formulas = FALSE)
}

#' @export  
OmicSelector_signiture_overlap <- function(...) {
  .Deprecated("analyze_signature_overlap", 
              msg = "OmicSelector_signiture_overlap is deprecated. Use analyze_signature_overlap() instead.")
  
  # Convert old call to new function
  signatures <- list(...)
  result <- analyze_signature_overlap(signatures, return_details = FALSE)
  
  # Return just the summary for backward compatibility
  return(result$signature_summary)
}
