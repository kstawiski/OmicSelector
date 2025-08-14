#' S3 Methods for OmicSelector Objects
#'
#' @description
#' S3 methods for OmicSelector result objects
#'
#' @name omicselector-methods
NULL

#' Print Method for OmicSelector Objects
#'
#' @param x OmicSelector object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector <- function(x, ...) {
  cat("OmicSelector Results\n")
  cat("===================\n\n")
  
  cat("Run ID:", x$stamp, "\n")
  cat("Started:", format(x$started_at), "\n")
  
  if (!is.null(x$completed_at)) {
    cat("Completed:", format(x$completed_at), "\n")
    cat("Runtime:", round(x$runtime, 2), "minutes\n")
  } else {
    cat("Status: Running or incomplete\n")
  }
  
  cat("\nData Summary:\n")
  cat("- Training samples:", x$data_info$train_samples, "\n")
  cat("- Features:", x$data_info$train_features, "\n")
  cat("- Test samples:", x$data_info$test_samples, "\n")
  
  cat("\nMethods:\n")
  cat("- Requested:", length(x$methods), "\n")
  cat("- Completed:", length(x$formulas), "\n")
  cat("- Failed:", length(x$errors), "\n")
  
  if (length(x$formulas) > 0) {
    cat("\nTop formulas:\n")
    for (i in seq_len(min(3, length(x$formulas)))) {
      formula_str <- x$formulas[[i]]
      if (nchar(formula_str) > 60) {
        formula_str <- paste0(substr(formula_str, 1, 57), "...")
      }
      cat(" ", i, ":", formula_str, "\n")
    }
    if (length(x$formulas) > 3) {
      cat("  ... and", length(x$formulas) - 3, "more\n")
    }
  }
  
  cat("\nUse summary() for more details\n")
  invisible(x)
}

#' Summary Method for OmicSelector Objects
#'
#' @param object OmicSelector object
#' @param ... Additional arguments (not used)
#' @export
summary.OmicSelector <- function(object, ...) {
  cat("OmicSelector Results Summary\n")
  cat("===========================\n\n")
  
  # Basic info
  print(object)
  
  # Configuration details
  cat("\nConfiguration:\n")
  config_items <- c("max_iterations", "cores", "prefer_no_features", "timeout_sec")
  for (item in config_items) {
    if (item %in% names(object$config)) {
      cat("- ", item, ":", object$config[[item]], "\n")
    }
  }
  
  # Method details
  if (length(object$method_results) > 0) {
    cat("\nMethod Results:\n")
    for (method_id in names(object$method_results)) {
      result <- object$method_results[[method_id]]
      n_features <- length(result$features)
      cat("- Method", method_id, ":", n_features, "features selected\n")
    }
  }
  
  # Errors and warnings
  if (length(object$errors) > 0) {
    cat("\nErrors:\n")
    for (method_id in names(object$errors)) {
      cat("- Method", method_id, ":", object$errors[[method_id]], "\n")
    }
  }
  
  invisible(object)
}

#' Extract Formulas from OmicSelector Object
#'
#' @param x OmicSelector object
#' @param ... Additional arguments (not used)
#' @return Character vector of formulas
#' @export
formulas <- function(x, ...) {
  UseMethod("formulas")
}

#' @export
formulas.OmicSelector <- function(x, ...) {
  return(unlist(x$formulas))
}

#' Extract Features from OmicSelector Object
#'
#' @param x OmicSelector object
#' @param method_id Optional specific method ID to extract features from
#' @param ... Additional arguments (not used)
#' @return List of character vectors containing selected features
#' @export
features <- function(x, method_id = NULL, ...) {
  UseMethod("features")
}

#' @export
features.OmicSelector <- function(x, method_id = NULL, ...) {
  if (is.null(method_id)) {
    # Return all features from all methods
    all_features <- list()
    for (mid in names(x$method_results)) {
      all_features[[mid]] <- x$method_results[[mid]]$features
    }
    return(all_features)
  } else {
    # Return features from specific method
    if (as.character(method_id) %in% names(x$method_results)) {
      return(x$method_results[[as.character(method_id)]]$features)
    } else {
      warning("Method ", method_id, " not found in results")
      return(NULL)
    }
  }
}

#' Plot Method for OmicSelector Objects
#'
#' @param x OmicSelector object
#' @param type Character string specifying plot type
#' @param ... Additional arguments passed to plotting functions
#' @export
plot.OmicSelector <- function(x, type = "overview", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  switch(type,
    "overview" = plot_overview(x, ...),
    "features" = plot_feature_counts(x, ...),
    "methods" = plot_method_success(x, ...),
    stop("Unknown plot type: ", type)
  )
}

#' Plot Overview of Results
#'
#' @param x OmicSelector object
#' @param ... Additional arguments
#' @keywords internal
plot_overview <- function(x, ...) {
  
  # Create summary data
  method_data <- data.frame(
    method = names(x$method_results),
    features = sapply(x$method_results, function(r) length(r$features)),
    stringsAsFactors = FALSE
  )
  
  if (nrow(method_data) == 0) {
    message("No completed methods to plot")
    return(invisible(NULL))
  }
  
  p <- ggplot2::ggplot(method_data, ggplot2::aes(x = reorder(method, features), y = features)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Number of Features Selected by Method",
      x = "Method",
      y = "Number of Features",
      subtitle = paste("Run:", x$stamp)
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  invisible(p)
}

#' Plot Feature Count Distribution
#'
#' @param x OmicSelector object
#' @param ... Additional arguments
#' @keywords internal
plot_feature_counts <- function(x, ...) {
  
  feature_counts <- sapply(x$method_results, function(r) length(r$features))
  
  if (length(feature_counts) == 0) {
    message("No completed methods to plot")
    return(invisible(NULL))
  }
  
  count_data <- data.frame(
    count = feature_counts,
    method = names(feature_counts),
    stringsAsFactors = FALSE
  )
  
  p <- ggplot2::ggplot(count_data, ggplot2::aes(x = count)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "steelblue", alpha = 0.7) +
    ggplot2::labs(
      title = "Distribution of Feature Counts",
      x = "Number of Features Selected",
      y = "Number of Methods",
      subtitle = paste("Run:", x$stamp)
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  invisible(p)
}

#' Plot Method Success Rate
#'
#' @param x OmicSelector object
#' @param ... Additional arguments
#' @keywords internal
plot_method_success <- function(x, ...) {
  
  success_data <- data.frame(
    status = c("Completed", "Failed"),
    count = c(length(x$formulas), length(x$errors)),
    stringsAsFactors = FALSE
  )
  
  p <- ggplot2::ggplot(success_data, ggplot2::aes(x = status, y = count, fill = status)) +
    ggplot2::geom_col(alpha = 0.7) +
    ggplot2::scale_fill_manual(values = c("Completed" = "green", "Failed" = "red")) +
    ggplot2::labs(
      title = "Method Success Rate",
      x = "Status",
      y = "Number of Methods",
      subtitle = paste("Run:", x$stamp)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  
  print(p)
  invisible(p)
}
