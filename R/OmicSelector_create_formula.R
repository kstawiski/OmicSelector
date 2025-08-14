#' Create Formula for Feature Selection
#'
#' @description
#' Modern helper function to create R formulas based on selected features.
#' Includes improved validation and error handling.
#'
#' @param selected_features Character vector of selected feature names
#' @param class_column Character string specifying the class column name (default: "Class")
#' @param validate_features Logical indicating whether to validate feature names (default: TRUE)
#' 
#' @return Formula object of the form "Class ~ feature1 + feature2 + ..."
#' 
#' @examples 
#' \dontrun{
#' # Basic usage
#' features <- c("hsa_miR_1", "hsa_miR_2", "hsa_miR_3")
#' formula <- create_omics_formula(features)
#' 
#' # Custom class column
#' formula <- create_omics_formula(features, class_column = "Response")
#' }
#'
#' @export
create_omics_formula <- function(selected_features, 
                                class_column = "Class", 
                                validate_features = TRUE) {
  
  # Input validation
  if (!is.character(selected_features)) {
    stop("'selected_features' must be a character vector", call. = FALSE)
  }
  
  if (!is.character(class_column) || length(class_column) != 1) {
    stop("'class_column' must be a single character string", call. = FALSE)
  }
  
  # Remove missing values
  selected_features <- selected_features[!is.na(selected_features)]
  
  # Check if any features remain
  if (length(selected_features) == 0) {
    stop("No valid features provided after removing NAs", call. = FALSE)
  }
  
  # Optional feature name validation
  if (validate_features) {
    invalid_features <- selected_features[!grepl("^[a-zA-Z][a-zA-Z0-9._]*$", selected_features)]
    if (length(invalid_features) > 0) {
      warning("Some feature names may be invalid: ", 
              paste(head(invalid_features, 3), collapse = ", "),
              if (length(invalid_features) > 3) "...", 
              call. = FALSE)
    }
  }
  
  # Create formula with proper escaping
  feature_string <- paste0("`", selected_features, "`", collapse = " + ")
  formula_string <- paste(class_column, "~", feature_string)
  
  tryCatch({
    as.formula(formula_string)
  }, error = function(e) {
    stop("Failed to create formula: ", e$message, call. = FALSE)
  })
}

#' @rdname create_omics_formula
#' @export
OmicSelector_create_formula <- function(selected_features) {
  # Backward compatibility wrapper
  .Deprecated("create_omics_formula", 
              msg = "OmicSelector_create_formula is deprecated. Use create_omics_formula instead.")
  
  if (length(selected_features) == 1 && selected_features[1] == ".") {
    stop("Too few features selected.", call. = FALSE)
  }
  
  create_omics_formula(selected_features)
}
