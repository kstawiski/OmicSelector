#' Modern Volcano Plot for Differential Expression Analysis
#'
#' @description
#' Create publication-ready volcano plots for differential expression analysis
#' using modern ggplot2 with customizable highlighting and labeling.
#'
#' @param de_results Data frame with differential expression results containing columns:
#'   log2FoldChange, pvalue (and optionally padj for adjusted p-values)
#' @param feature_column Character string indicating column name with feature identifiers (default: "feature")
#' @param fc_column Character string indicating fold change column name (default: "log2FoldChange")
#' @param pval_column Character string indicating p-value column name (default: "pvalue")
#' @param padj_column Character string indicating adjusted p-value column name (default: "padj")
#' @param selected_features Character vector of features to highlight
#' @param label_features Character vector of features to label (subset of selected_features)
#' @param use_adjusted_p Logical indicating whether to use adjusted p-values (default: TRUE)
#' @param fc_threshold Numeric fold change threshold for significance (default: 1)
#' @param p_threshold Numeric p-value threshold for significance (default: 0.05)
#' @param point_size Numeric size for points (default: 1.5)
#' @param label_size Numeric size for labels (default: 3)
#' @param colors Named character vector with colors for different point types
#' @param title Character string for plot title
#' @param ... Additional arguments passed to ggplot
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#' # Create sample differential expression data
#' de_data <- data.frame(
#'   feature = paste0("miR_", 1:100),
#'   log2FoldChange = rnorm(100, 0, 2),
#'   pvalue = runif(100, 0, 0.1),
#'   padj = p.adjust(runif(100, 0, 0.1), method = "BH")
#' )
#' 
#' # Basic volcano plot
#' create_volcano_plot(de_data)
#' 
#' # Highlight specific features
#' selected <- de_data$feature[1:10]
#' create_volcano_plot(de_data, selected_features = selected)
#' }
#'
#' @export
create_volcano_plot <- function(de_results,
                               feature_column = "feature",
                               fc_column = "log2FoldChange", 
                               pval_column = "pvalue",
                               padj_column = "padj",
                               selected_features = NULL,
                               label_features = NULL,
                               use_adjusted_p = TRUE,
                               fc_threshold = 1,
                               p_threshold = 0.05,
                               point_size = 1.5,
                               label_size = 3,
                               colors = c("significant" = "#E31A1C", 
                                        "selected" = "#1F78B4", 
                                        "ns" = "grey60"),
                               title = "Volcano Plot",
                               ...) {
  
  # Input validation
  if (!is.data.frame(de_results)) {
    stop("'de_results' must be a data frame", call. = FALSE)
  }
  
  required_cols <- c(feature_column, fc_column, pval_column)
  missing_cols <- setdiff(required_cols, colnames(de_results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  if (use_adjusted_p) {
    if (!padj_column %in% colnames(de_results)) {
      stop("Column '", padj_column, "' not found for adjusted p-values", call. = FALSE)
    }
    p_col <- padj_column
  } else {
    p_col <- pval_column
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required", call. = FALSE)
  }
  
  # Prepare data
  plot_data <- de_results
  plot_data$log10p <- -log10(plot_data[[p_col]])
  plot_data$feature_id <- plot_data[[feature_column]]
  plot_data$fold_change <- plot_data[[fc_column]]
  
  # Remove infinite and missing values
  plot_data <- plot_data[is.finite(plot_data$log10p) & 
                        is.finite(plot_data$fold_change) & 
                        !is.na(plot_data$log10p) & 
                        !is.na(plot_data$fold_change), ]
  
  if (nrow(plot_data) == 0) {
    stop("No valid data points after removing missing/infinite values", call. = FALSE)
  }
  
  # Classify points
  plot_data$significance <- "ns"
  significant_idx <- abs(plot_data$fold_change) >= fc_threshold & 
                   plot_data[[p_col]] <= p_threshold
  plot_data$significance[significant_idx] <- "significant"
  
  # Mark selected features
  if (!is.null(selected_features)) {
    selected_idx <- plot_data$feature_id %in% selected_features
    plot_data$significance[selected_idx] <- "selected"
  }
  
  # Prepare labels
  plot_data$label <- NA_character_
  if (!is.null(label_features)) {
    label_idx <- plot_data$feature_id %in% label_features
    plot_data$label[label_idx] <- plot_data$feature_id[label_idx]
  } else if (!is.null(selected_features)) {
    selected_idx <- plot_data$feature_id %in% selected_features
    plot_data$label[selected_idx] <- plot_data$feature_id[selected_idx]
  }
  
  # Create plot - using strings to avoid NSE issues
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "fold_change", y = "log10p")) +
    ggplot2::geom_point(ggplot2::aes_string(color = "significance"), 
                       size = point_size, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors, 
                               name = "Significance",
                               labels = c("ns" = "Not significant",
                                        "significant" = "Significant", 
                                        "selected" = "Selected")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = paste("log2(Fold Change)"),
      y = paste0("-log10(", if (use_adjusted_p) "adjusted " else "", "p-value)")
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  
  # Add significance thresholds
  if (fc_threshold > 0) {
    p <- p + 
      ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
                         linetype = "dashed", alpha = 0.7) 
  }
  
  if (p_threshold < 1) {
    p <- p + 
      ggplot2::geom_hline(yintercept = -log10(p_threshold), 
                         linetype = "dashed", alpha = 0.7)
  }
  
  # Add labels
  if (any(!is.na(plot_data$label))) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      # Fallback to basic text labels
      label_data <- plot_data[!is.na(plot_data$label), ]
      p <- p + ggplot2::geom_text(data = label_data,
                                 ggplot2::aes_string(label = "label"),
                                 size = label_size, 
                                 hjust = 0, vjust = 0,
                                 nudge_x = 0.1, nudge_y = 0.1)
    } else {
      # Use ggrepel for better label positioning
      label_data <- plot_data[!is.na(plot_data$label), ]
      p <- p + ggrepel::geom_text_repel(data = label_data,
                                       ggplot2::aes_string(label = "label"),
                                       size = label_size,
                                       max.overlaps = 20)
    }
  }
  
  return(p)
}

#' @rdname create_volcano_plot
#' @export
OmicSelector_vulcano_plot <- function(selected_miRNAs, 
                                     DE = NULL, 
                                     only_label = NULL, 
                                     take_adjusted_p = FALSE) {
  
  # Backward compatibility wrapper
  .Deprecated("create_volcano_plot", 
              msg = "OmicSelector_vulcano_plot is deprecated. Use create_volcano_plot instead.")
  
  if (is.null(DE)) {
    stop("DE parameter is required for legacy function", call. = FALSE)
  }
  
  # Handle legacy column names - try to detect them
  feature_col <- if ("miR" %in% colnames(DE)) "miR" else 
                if ("feature" %in% colnames(DE)) "feature" else
                colnames(DE)[1]
  
  fc_col <- if ("log2FC" %in% colnames(DE)) "log2FC" else
           if ("logFC" %in% colnames(DE)) "logFC" else
           if ("log2FoldChange" %in% colnames(DE)) "log2FoldChange" else
           stop("Could not identify fold change column")
  
  pval_col <- if ("pvalue" %in% colnames(DE)) "pvalue" else
             if ("P.Value" %in% colnames(DE)) "P.Value" else
             if ("p.value" %in% colnames(DE)) "p.value" else
             stop("Could not identify p-value column")
  
  padj_col <- if ("padj" %in% colnames(DE)) "padj" else
             if ("adj.P.Val" %in% colnames(DE)) "adj.P.Val" else
             if ("p.adjust" %in% colnames(DE)) "p.adjust" else
             "padj"  # Default, may not exist
  
  # Create modern volcano plot
  create_volcano_plot(
    de_results = DE,
    feature_column = feature_col,
    fc_column = fc_col,
    pval_column = pval_col,
    padj_column = padj_col,
    selected_features = selected_miRNAs,
    label_features = only_label,
    use_adjusted_p = take_adjusted_p
  )
}
