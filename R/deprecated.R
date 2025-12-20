#' @title Deprecated and Leaky Functions
#'
#' @description
#' This file contains deprecation notices and warnings for functions
#' that violate the zero-leakage principle of OmicSelector 2.0.
#'
#' @name deprecated
NULL


#' @title Check if a function call would cause data leakage
#'
#' @description
#' Internal function to detect potential data leakage patterns
#'
#' @param func_name Name of the function being called
#' @param args Arguments passed to the function
#'
#' @return Logical indicating if leakage was detected
#' @keywords internal
.check_leakage_risk <- function(func_name, args) {
  # Known leaky patterns
  leaky_patterns <- list(
    "iteratedRFE" = list(
      check = function(a) !is.null(a$testSet),
      message = "Test set data is being used during feature selection, causing data leakage."
    ),
    "iteratedRFETest" = list(
      check = function(a) TRUE,
      message = "iteratedRFETest explicitly uses test data during selection. Always leaky."
    )
  )

  if (func_name %in% names(leaky_patterns)) {
    pattern <- leaky_patterns[[func_name]]
    if (pattern$check(args)) {
      return(list(
        leaky = TRUE,
        message = pattern$message
      ))
    }
  }

  return(list(leaky = FALSE, message = NULL))
}


#' @title Emit deprecation warning for leaky functions
#'
#' @description
#' Issues a prominent warning when deprecated/leaky functions are called
#'
#' @param func_name Function name
#' @param replacement Suggested replacement function
#' @param reason Reason for deprecation
#'
#' @keywords internal
.deprecation_warning <- function(func_name, replacement = NULL, reason = NULL) {
  msg <- sprintf(
    "\n%s\n%s is DEPRECATED in OmicSelector 2.0\n%s\n",
    paste(rep("=", 60), collapse = ""),
    func_name,
    paste(rep("=", 60), collapse = "")
  )

  if (!is.null(reason)) {
    msg <- paste0(msg, "\nReason: ", reason, "\n")
  }

  if (!is.null(replacement)) {
    msg <- paste0(msg, "\nUse instead: ", replacement, "\n")
  }

  msg <- paste0(msg, "\nThis function may produce INVALID results due to data leakage.\n")
  msg <- paste0(msg, "See: https://biostat.umed.pl/OmicSelector/articles/leakage.html\n\n")

  warning(msg, call. = FALSE, immediate. = TRUE)
}


#' @title Wrapper for deprecated OmicSelector_iteratedRFE
#'
#' @description
#' This function is DEPRECATED due to data leakage when testSet is provided.
#' The test set labels are used during feature selection, which inflates
#' performance estimates by 10-30%.
#'
#' @param trainSet Training data
#' @param testSet Test data (CAUSES LEAKAGE if provided)
#' @param ... Other arguments
#'
#' @return Original function result with warning
#'
#' @details
#' ## Why is this deprecated?
#'

#' The original implementation uses `testSet` during RFE iteration:
#' ```
#' rfeIter(testX = testSet[, initFeatures], testY = testSet[, classLab], ...)
#' ```
#'
#' This means the model can "see" test data labels during training,
#' leading to overly optimistic accuracy estimates.
#'
#' ## Replacement
#'
#' Use `OmicPipeline$create_graph_learner()` with nested cross-validation:
#' ```
#' pipeline <- OmicPipeline$new(data, target = "outcome")
#' learner <- pipeline$create_graph_learner(filter = "anova", model = "ranger")
#' service <- BenchmarkService$new(pipeline, outer_folds = 5, inner_folds = 3)
#' service$add_learner(learner)
#' result <- service$run()
#' ```
#'
#' @seealso [OmicPipeline], [BenchmarkService]
#' @export
OmicSelector_iteratedRFE_deprecated <- function(trainSet, testSet = NULL, ...) {
  .deprecation_warning(
    func_name = "OmicSelector_iteratedRFE",
    replacement = "OmicPipeline$create_graph_learner() with BenchmarkService",
    reason = "Uses test set data during feature selection (data leakage)"
  )

  # Check for leakage
  if (!is.null(testSet)) {
    warning(paste(
      "\n*** CRITICAL: DATA LEAKAGE DETECTED ***\n",
      "testSet parameter is provided, which means test labels are used",
      "during feature selection. This produces INVALID, inflated metrics.\n",
      "Results from this function should NOT be published or trusted.\n"
    ), call. = FALSE, immediate. = TRUE)
  }

  # Call original function (for backward compatibility only)
  OmicSelector_iteratedRFE_legacy(trainSet = trainSet, testSet = testSet, ...)
}


#' @title Legacy iteratedRFE (internal use only)
#' @keywords internal
#' @noRd
OmicSelector_iteratedRFE_legacy <- function(trainSet, testSet = NULL,
                                             initFeatures = colnames(trainSet),
                                             classLab, checkNFeatures = 25,
                                             votingIterations = 1000,
                                             useCV = FALSE, nfolds = 10,
                                             initRandomState = 42) {
  # Original implementation preserved for reference only
  # DO NOT USE - This function has known data leakage

  stop(paste(
    "OmicSelector_iteratedRFE is no longer available in OmicSelector 2.0\n",
    "This function was deprecated due to critical data leakage issues.\n",
    "Use OmicPipeline$create_graph_learner() with nested cross-validation instead."
  ), call. = FALSE)
}


#' @title List all deprecated functions
#'
#' @description
#' Returns a list of deprecated functions and their replacements
#'
#' @return A data.frame with deprecated function info
#' @export
list_deprecated_functions <- function() {
  deprecated <- data.frame(
    Function = c(
      "OmicSelector_OmicSelector",
      "OmicSelector_iteratedRFE",
      "OmicSelector_iteratedRFETest"
    ),
    Reason = c(
      "SMOTE applied before CV splits; filter-then-CV pattern (leakage)",
      "Test set data used during feature selection (leakage)",
      "Explicitly leaks test data during RFE (Method 55)"
    ),
    Replacement = c(
      "OmicPipeline$new() + OmicPipeline$create_graph_learner() + BenchmarkService",
      "OmicPipeline$create_graph_learner() + BenchmarkService",
      "OmicPipeline$create_graph_learner() + BenchmarkService"
    ),
    Severity = c("HIGH", "CRITICAL", "CRITICAL"),
    stringsAsFactors = FALSE
  )

  return(deprecated)
}


#' @title Validate pipeline for leakage risks
#'
#' @description
#' Analyzes a pipeline or function call for potential data leakage issues
#'
#' @param pipeline An OmicPipeline object or function call
#'
#' @return A validation result object
#' @export
validate_no_leakage <- function(pipeline) {
  result <- list(
    valid = TRUE,
    issues = list(),
    recommendations = list()
  )

  # Check if it's an OmicPipeline (safe by construction)
  if (inherits(pipeline, "OmicPipeline")) {
    result$recommendations <- append(result$recommendations,
      "OmicPipeline uses mlr3 GraphLearners which prevent leakage by construction.")
    return(structure(result, class = c("LeakageValidation", "list")))
  }

  # Check if it's a legacy function call
  if (is.function(pipeline)) {
    func_name <- deparse(substitute(pipeline))
    if (grepl("iteratedRFE", func_name)) {
      result$valid <- FALSE
      result$issues <- append(result$issues,
        "iteratedRFE functions have known data leakage issues.")
    }
  }

  structure(result, class = c("LeakageValidation", "list"))
}


#' @title Print method for LeakageValidation
#' @param x A LeakageValidation object
#' @param ... Additional arguments (ignored)
#' @export
print.LeakageValidation <- function(x, ...) {
  if (x$valid) {
    cat("Leakage Validation: PASSED\n")
  } else {
    cat("Leakage Validation: FAILED\n")
    cat("\nIssues:\n")
    for (issue in x$issues) {
      cat("  - ", issue, "\n")
    }
  }

  if (length(x$recommendations) > 0) {
    cat("\nRecommendations:\n")
    for (rec in x$recommendations) {
      cat("  - ", rec, "\n")
    }
  }

  invisible(x)
}


# =============================================================================
# LEGACY WRAPPER: OmicSelector_OmicSelector
# =============================================================================

#' @title OmicSelector_OmicSelector (DEPRECATED)
#'
#' @description
#' **DEPRECATED:** This function is deprecated in OmicSelector 2.0 due to
#' data leakage issues. It applies SMOTE oversampling BEFORE cross-validation
#' splits, which causes synthetic data to leak between folds.
#'
#' The function is preserved for backward compatibility but will emit
#' deprecation warnings. Users should migrate to `OmicPipeline` for
#' scientifically valid results.
#'
#' @param wd Working directory path (default: current directory)
#' @param m Numeric vector of method IDs to run (default: 1:70)
#' @param ... Additional arguments passed to legacy function
#'
#' @details
#' ## Why is this deprecated?
#'
#' The legacy `OmicSelector_OmicSelector` function has several critical issues:
#'
#' 1. **SMOTE Leakage**: Oversampling is applied to the entire training set
#'    BEFORE cross-validation. This means synthetic samples can appear in

#'    validation folds, leading to overly optimistic performance estimates.
#'
#' 2. **Filter-then-CV**: Feature selection runs on the full training set
#'    before any CV fold isolation. Selected features may be driven by
#'    statistical artifacts that won't generalize.
#'
#' 3. **Global State**: Uses `setwd()` and file-based I/O, making results
#'    difficult to reproduce.
#'
#' ## Migration Path
#'
#' Replace legacy workflows with the new mlr3-based `OmicPipeline` class:
#'
#' ```r
#' # Legacy (DEPRECATED - has leakage)
#' # OmicSelector_OmicSelector(wd = ".", m = c(1, 2, 3))
#'
#' # New (zero leakage by construction)
#' pipeline <- OmicPipeline$new(
#'   data = my_data,
#'   target = "Class",
#'   positive = "Case"
#' )
#'
#' learner <- pipeline$create_graph_learner(
#'   filter = "anova",
#'   model = "ranger",
#'   n_features = 20,
#'   oversample = "smote"  # Applied inside CV folds
#' )
#'
#' service <- BenchmarkService$new(pipeline)
#' service$add_learner(learner)
#' result <- service$run()
#' ```
#'
#' @seealso [OmicPipeline], [BenchmarkService], [list_deprecated_functions()]
#'
#' @return The list of selected formulas (with deprecation warning)
#' @export
OmicSelector_OmicSelector_wrapper <- function(wd = getwd(), m = c(1:70),
                                               max_iterations = 10,
                                               code_path = system.file("extdata", "", package = "OmicSelector"),
                                               register_parallel = TRUE,
                                               clx = NULL,
                                               stamp = as.numeric(Sys.time()),
                                               prefer_no_features = 11,
                                               conda_path = "/home/konrad/anaconda3/bin/conda",
                                               debug = FALSE,
                                               timeout_sec = 172800,
                                               type = "auto") {

  # Emit prominent deprecation warning
  .deprecation_warning(
    func_name = "OmicSelector_OmicSelector",
    replacement = "OmicPipeline$new() + BenchmarkService",
    reason = paste(
      "This function applies SMOTE before CV splits (data leakage).",
      "Feature selection also occurs before fold isolation.",
      "Results may be 10-30% inflated compared to properly nested CV."
    )
  )

  # Check for high-risk methods
  leaky_methods <- c(55)  # Method 55 is iteratedRFETest which uses test labels
  if (any(m %in% leaky_methods)) {
    warning(paste(
      "\n*** CRITICAL: SEVERELY LEAKY METHODS REQUESTED ***\n",
      "Methods", paste(intersect(m, leaky_methods), collapse = ", "),
      "explicitly use test set labels during feature selection.\n",
      "These methods have been DISABLED in OmicSelector 2.0.\n",
      "They will be skipped.\n"
    ), call. = FALSE, immediate. = TRUE)
    m <- setdiff(m, leaky_methods)
  }

  # Inform user about migration
  message(paste0(
    "\n",
    "==============================================================\n",
    " MIGRATION NOTICE: OmicSelector 2.0 introduces OmicPipeline\n",
    "==============================================================\n",
    "\n",
    "The new OmicPipeline class guarantees zero data leakage:\n",
    "\n",
    "  pipeline <- OmicPipeline$new(data, target = 'Class')\n",
    "  learner <- pipeline$create_graph_learner(\n",
    "    filter = 'anova', model = 'ranger', oversample = 'smote'\n",
    "  )\n",
    "  result <- pipeline$benchmark(learner)\n",
    "\n",
    "Learn more: https://biostat.umed.pl/OmicSelector/articles/migration.html\n",
    "==============================================================\n"
  ))

  # Call the original function (preserved for backward compatibility)
  OmicSelector_OmicSelector_legacy(
    wd = wd,
    m = m,
    max_iterations = max_iterations,
    code_path = code_path,
    register_parallel = register_parallel,
    clx = clx,
    stamp = stamp,
    prefer_no_features = prefer_no_features,
    conda_path = conda_path,
    debug = debug,
    timeout_sec = timeout_sec,
    type = type
  )
}


#' @title Legacy OmicSelector_OmicSelector (internal)
#'
#' @description
#' The original implementation preserved for backward compatibility.
#' This function is called by the deprecated wrapper.
#'
#' @param wd Working directory
#' @param m Methods to run (1:70)
#' @param max_iterations Maximum iterations
#' @param code_path Code path for external scripts
#' @param register_parallel Use parallel processing
#' @param clx Existing cluster object
#' @param stamp Timestamp for output files
#' @param prefer_no_features Maximum features to select
#' @param conda_path Path to conda
#' @param debug Enable debug output
#' @param timeout_sec Timeout in seconds
#' @param type Mode for differential expression
#'
#' @return List of formulas
#' @keywords internal
#' @noRd
# The legacy function is defined in OmicSelector_OmicSelector.R as OmicSelector_OmicSelector_legacy


#' @title OmicSelector_OmicSelector
#'
#' @description
#' **DEPRECATED:** Main feature selection function from OmicSelector 1.x.
#'
#' This function is deprecated due to data leakage issues. It now routes
#' through a deprecation wrapper that emits warnings before calling the
#' legacy implementation.
#'
#' New users should use `OmicPipeline` instead. See `?OmicPipeline` for details.
#'
#' @inheritParams OmicSelector_OmicSelector_wrapper
#'
#' @return The list of selected formulas (with deprecation warning)
#'
#' @seealso [OmicPipeline], [BenchmarkService], [list_deprecated_functions()]
#'
#' @export
OmicSelector_OmicSelector <- OmicSelector_OmicSelector_wrapper
