#' Modern Heatmap Visualization for Omics Data
#'
#' @description
#' Create publication-ready heatmaps for omics data using modern ggplot2 and
#' pheatmap approaches. Supports multiple annotation tracks and customization options.
#'
#' @param expression_data Matrix or data frame with features in columns and samples in rows
#' @param sample_annotations Data frame with sample metadata for annotation tracks
#' @param feature_annotations Data frame with feature metadata (optional)
#' @param scale_data Character indicating how to scale data: "none", "row", "column", "zscore" (default: "zscore")
#' @param cluster_rows Logical indicating whether to cluster rows (default: TRUE)
#' @param cluster_cols Logical indicating whether to cluster columns (default: TRUE)
#' @param show_row_names Logical indicating whether to show row names (default: FALSE)
#' @param show_col_names Logical indicating whether to show column names (default: TRUE)
#' @param color_palette Character vector of colors or function to generate colors
#' @param title Character string for plot title
#' @param trim_range Numeric vector of length 2 for trimming extreme values
#' @param method Character indicating plotting method: "ggplot2" or "pheatmap" (default: "pheatmap")
#' @param ... Additional arguments passed to plotting functions
#'
#' @return ggplot object (if method="ggplot2") or pheatmap object (if method="pheatmap")
#'
#' @examples
#' \dontrun{
#' # Basic heatmap
#' data <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' colnames(data) <- paste0("Feature_", 1:10)
#' rownames(data) <- paste0("Sample_", 1:20)
#' 
#' # Simple heatmap
#' create_modern_heatmap(data)
#' 
#' # With annotations
#' annotations <- data.frame(
#'   Class = factor(rep(c("Case", "Control"), each = 10)),
#'   Batch = factor(rep(c("A", "B"), 10))
#' )
#' create_modern_heatmap(data, sample_annotations = annotations)
#' }
#'
#' @export
create_modern_heatmap <- function(expression_data,
                                sample_annotations = NULL,
                                feature_annotations = NULL,
                                scale_data = "zscore",
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                show_row_names = FALSE,
                                show_col_names = TRUE,
                                color_palette = NULL,
                                title = NULL,
                                trim_range = NULL,
                                method = "pheatmap",
                                ...) {
  
  # Input validation
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("'expression_data' must be a matrix or data frame", call. = FALSE)
  }
  
  if (!is.null(sample_annotations) && nrow(sample_annotations) != nrow(expression_data)) {
    stop("Number of rows in 'sample_annotations' must match rows in 'expression_data'", call. = FALSE)
  }
  
  # Convert to matrix if needed
  if (is.data.frame(expression_data)) {
    numeric_cols <- sapply(expression_data, is.numeric)
    if (!all(numeric_cols)) {
      warning("Non-numeric columns detected and will be removed", call. = FALSE)
      expression_data <- expression_data[, numeric_cols, drop = FALSE]
    }
    expression_data <- as.matrix(expression_data)
  }
  
  # Handle missing values
  if (any(is.na(expression_data))) {
    warning("Missing values detected; they will be replaced with column means", call. = FALSE)
    for (i in seq_len(ncol(expression_data))) {
      col_mean <- mean(expression_data[, i], na.rm = TRUE)
      expression_data[is.na(expression_data[, i]), i] <- col_mean
    }
  }
  
  # Scale data
  scaled_data <- switch(scale_data,
    "none" = expression_data,
    "row" = t(scale(t(expression_data))),
    "column" = scale(expression_data),
    "zscore" = scale(expression_data),
    stop("Invalid scale_data option. Use 'none', 'row', 'column', or 'zscore'", call. = FALSE)
  )
  
  # Trim extreme values if requested
  if (!is.null(trim_range)) {
    if (length(trim_range) != 2 || trim_range[1] >= trim_range[2]) {
      stop("'trim_range' must be a numeric vector of length 2 with min < max", call. = FALSE)
    }
    scaled_data[scaled_data < trim_range[1]] <- trim_range[1]
    scaled_data[scaled_data > trim_range[2]] <- trim_range[2]
  }
  
  # Set default color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  }
  
  # Choose plotting method
  if (method == "pheatmap") {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      stop("Package 'pheatmap' is required for method='pheatmap'", call. = FALSE)
    }
    
    # Prepare annotation data frame
    annotation_row <- if (!is.null(sample_annotations)) {
      sample_annotations
    } else {
      NULL
    }
    
    annotation_col <- if (!is.null(feature_annotations)) {
      feature_annotations
    } else {
      NULL
    }
    
    # Create pheatmap
    p <- pheatmap::pheatmap(
      scaled_data,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      show_rownames = show_row_names,
      show_colnames = show_col_names,
      annotation_row = annotation_row,
      annotation_col = annotation_col,
      color = color_palette,
      main = title,
      ...
    )
    
    return(p)
    
  } else if (method == "ggplot2") {
    if (!requireNamespace("ggplot2", quietly = TRUE) || 
        !requireNamespace("reshape2", quietly = TRUE) ||
        !requireNamespace("rlang", quietly = TRUE)) {
      stop("Packages 'ggplot2', 'reshape2', and 'rlang' are required for method='ggplot2'", call. = FALSE)
    }
    
    # Convert to long format
    scaled_df <- as.data.frame(scaled_data)
    scaled_df$Sample <- rownames(scaled_df)
    melted_data <- reshape2::melt(scaled_df, id.vars = "Sample", 
                                 variable.name = "Feature", value.name = "Expression")
    
    # Add sample annotations if provided
    if (!is.null(sample_annotations)) {
      sample_annotations$Sample <- rownames(sample_annotations)
      melted_data <- merge(melted_data, sample_annotations, by = "Sample")
    }
    
    # Create base plot with proper NSE handling
    p <- ggplot2::ggplot(melted_data, ggplot2::aes(x = !!rlang::sym("Feature"), 
                                                   y = !!rlang::sym("Sample"), 
                                                   fill = !!rlang::sym("Expression"))) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradientn(colors = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = if (show_col_names) ggplot2::element_text(angle = 45, hjust = 1) else ggplot2::element_blank(),
        axis.text.y = if (show_row_names) ggplot2::element_text() else ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(title = title, x = "Features", y = "Samples")
    
    return(p)
    
  } else {
    stop("Invalid method. Use 'pheatmap' or 'ggplot2'", call. = FALSE)
  }
}

#' @rdname create_modern_heatmap
#' @export
OmicSelector_heatmap <- function(x, 
                                rlab = NULL, 
                                zscore = FALSE, 
                                margins = c(10, 10), 
                                expression_name = "log10(TPM)", 
                                trim_min = NULL, 
                                trim_max = NULL, 
                                centered_on = NULL, 
                                legend_pos = "topright", 
                                legend_cex = 0.8, 
                                ...) {
  
  # Backward compatibility wrapper
  .Deprecated("create_modern_heatmap", 
              msg = "OmicSelector_heatmap is deprecated. Use create_modern_heatmap instead.")
  
  # Handle legacy parameter mapping
  scale_method <- if (zscore) "zscore" else "none"
  trim_range <- if (!is.null(trim_min) && !is.null(trim_max)) c(trim_min, trim_max) else NULL
  
  # Create modern heatmap with legacy-compatible settings
  create_modern_heatmap(
    expression_data = x,
    sample_annotations = rlab,
    scale_data = scale_method,
    trim_range = trim_range,
    method = "pheatmap",
    ...
  )
}
