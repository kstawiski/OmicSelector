#' @title TRIPOD+AI Report Generation for OmicSelector 2.0
#'
#' @description
#' Functions for generating reproducible, TRIPOD+AI compliant reports.
#' Reports are generated as immutable HTML/PDF artifacts.
#'
#' @details
#' TRIPOD+AI checklist compliance requires:
#' - Transparent reporting of model development and validation
#' - Calibration metrics (not just discrimination)
#' - Decision curve analysis for clinical utility
#' - Feature stability across resamples
#'
#' @name report
NULL


#' @title Create Report Data Schema
#'
#' @description
#' Creates a standardized data structure containing all information
#' required for TRIPOD+AI compliant reporting. This separates
#' computation from rendering.
#'
#' @param benchmark_result Result from BenchmarkService$run()
#' @param stability Optional NogueiraStability object
#' @param pipeline OmicPipeline object used in analysis
#' @param config Optional configuration list for provenance
#'
#' @return A ReportData object (list) with standardized structure
#'
#' @export
create_report_data <- function(benchmark_result,
                                stability = NULL,
                                pipeline = NULL,
                                config = NULL) {

  report_data <- list(
    # Metadata
    meta = list(
      timestamp = Sys.time(),
      omicselector_version = as.character(utils::packageVersion("OmicSelector")),
      r_version = R.version.string,
      session_info = utils::sessionInfo()
    ),

    # Data characteristics (TRIPOD items 4a-4b)
    data = list(
      n_samples = NULL,
      n_features = NULL,
      n_modalities = NULL,
      target_name = NULL,
      class_distribution = NULL,
      missing_data = NULL
    ),

    # Validation strategy (TRIPOD items 5a-5c, 10a-10d)
    validation = list(
      type = "internal_nested_cv",
      outer_folds = NULL,
      inner_folds = NULL,
      stratified = TRUE,
      grouped_by = NULL
    ),

    # Model specification (TRIPOD items 7a-7b, 10a)
    model = list(
      filter_method = NULL,
      learner_type = NULL,
      n_features_selected = NULL,
      hyperparameters = list()
    ),

    # Discrimination metrics (TRIPOD items 10a-10d)
    discrimination = list(
      auc = NULL,
      auc_ci = NULL,
      sensitivity = NULL,
      specificity = NULL,
      ppv = NULL,
      npv = NULL,
      brier_score = NULL
    ),

    # Calibration metrics (TRIPOD items 10a-10d)
    calibration = list(
      slope = NULL,
      intercept = NULL,
      ece = NULL,  # Expected Calibration Error
      calibration_data = NULL  # For plotting
    ),

    # Stability metrics (OmicSelector-specific)
    stability = list(
      nogueira_index = NULL,
      consensus_features = NULL,
      selection_frequency = NULL,
      interpretation = NULL
    ),

    # Reproducibility (TRIPOD+AI addition)
    reproducibility = list(
      seed = NULL,
      renv_hash = NULL,
      docker_tag = NULL,
      runtime_seconds = NULL
    )
  )

  # Populate from benchmark result
  if (!is.null(benchmark_result)) {
    report_data <- .populate_from_benchmark(report_data, benchmark_result)
  }

  # Populate from stability
  if (!is.null(stability)) {
    report_data <- .populate_from_stability(report_data, stability)
  }

  # Populate from pipeline
  if (!is.null(pipeline)) {
    report_data <- .populate_from_pipeline(report_data, pipeline)
  }

  # Add config for provenance
  if (!is.null(config)) {
    report_data$meta$config <- config
  }

  class(report_data) <- c("ReportData", "list")
  return(report_data)
}


#' @title Populate Report Data from Benchmark Result
#' @keywords internal
.populate_from_benchmark <- function(report_data, bmr) {

  # Handle both raw mlr3 BenchmarkResult and wrapped OmicBenchmarkResult
  if (inherits(bmr, "OmicBenchmarkResult")) {
    actual_bmr <- bmr$benchmark_result
    report_data$reproducibility$runtime_seconds <- as.numeric(
      difftime(bmr$timestamp, Sys.time(), units = "secs")
    )
  } else {
    actual_bmr <- bmr
  }

  # Extract performance metrics if mlr3 result
  if (inherits(actual_bmr, "BenchmarkResult")) {
    tryCatch({
      # Aggregate metrics
      agg <- actual_bmr$aggregate()

      # Try to get AUC
      if ("classif.auc" %in% names(agg)) {
        report_data$discrimination$auc <- agg$classif.auc
      }

      # Try to get other metrics
      if ("classif.bbrier" %in% names(agg)) {
        report_data$discrimination$brier_score <- agg$classif.bbrier
      }

    }, error = function(e) {
      warning("Could not extract metrics from BenchmarkResult: ", e$message)
    })
  }

  report_data
}


#' @title Populate Report Data from Stability Analysis
#' @keywords internal
.populate_from_stability <- function(report_data, stability) {

  if (inherits(stability, "NogueiraStability")) {
    report_data$stability$nogueira_index <- stability$nogueira_index
    report_data$stability$consensus_features <- stability$consensus_features
    report_data$stability$selection_frequency <- stability$selection_frequency
    report_data$stability$interpretation <- stability$interpretation
  }

  report_data
}


#' @title Populate Report Data from Pipeline
#' @keywords internal
.populate_from_pipeline <- function(report_data, pipeline) {

  if (inherits(pipeline, "OmicPipeline")) {
    task <- pipeline$get_task()

    report_data$data$n_samples <- task$nrow
    report_data$data$n_features <- length(pipeline$get_feature_names())
    report_data$data$target_name <- task$target_names

    # Class distribution
    if (inherits(task, "TaskClassif")) {
      report_data$data$class_distribution <- table(task$truth())
    }

    # Multi-omics info
    if (pipeline$is_multi_omics()) {
      mod_info <- pipeline$get_modality_info()
      report_data$data$n_modalities <- nrow(mod_info)
    } else {
      report_data$data$n_modalities <- 1
    }
  }

  report_data
}


#' @title Generate TRIPOD+AI Report
#'
#' @description
#' Renders a TRIPOD+AI compliant report from a ReportData object.
#'
#' @param results Benchmark results or ReportData object
#' @param output_file Output file path
#' @param format Output format: "html" (default) or "pdf"
#' @param template Custom Rmd template path (optional)
#'
#' @return Path to generated report (invisibly)
#'
#' @export
generate_tripod_report <- function(results,
                                    output_file = "tripod_report.html",
                                    format = c("html", "pdf"),
                                    template = NULL) {

  format <- match.arg(format)

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' required for report generation. Install with: install.packages('rmarkdown')")
  }

  # Convert to ReportData if needed
  if (!inherits(results, "ReportData")) {
    report_data <- create_report_data(results)
  } else {
    report_data <- results
  }

  # Use built-in template if not provided
  if (is.null(template)) {
    template <- system.file("rmd", "tripod_template.Rmd", package = "OmicSelector")
    if (template == "" || !file.exists(template)) {
      # Fallback: generate simple report
      return(.generate_simple_report(report_data, output_file))
    }
  }

  # Render report
  output_format <- switch(format,
    html = rmarkdown::html_document(
      toc = TRUE,
      toc_float = TRUE,
      code_folding = "hide",
      theme = "cosmo"
    ),
    pdf = rmarkdown::pdf_document(
      toc = TRUE,
      number_sections = TRUE
    )
  )

  output_path <- rmarkdown::render(
    input = template,
    output_format = output_format,
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    params = list(report_data = report_data),
    quiet = TRUE
  )

  message("Report generated: ", output_path)
  invisible(output_path)
}


#' @title Generate Simple Text Report
#' @keywords internal
.generate_simple_report <- function(report_data, output_file) {

  lines <- c(
    "=" %rep% 60,
    "OmicSelector 2.0 - TRIPOD+AI Analysis Report",
    "=" %rep% 60,
    "",
    sprintf("Generated: %s", format(report_data$meta$timestamp, "%Y-%m-%d %H:%M:%S")),
    sprintf("OmicSelector Version: %s", report_data$meta$omicselector_version),
    sprintf("R Version: %s", report_data$meta$r_version),
    "",
    "-" %rep% 40,
    "DATA CHARACTERISTICS",
    "-" %rep% 40,
    sprintf("Samples: %s", report_data$data$n_samples %||% "N/A"),
    sprintf("Features: %s", report_data$data$n_features %||% "N/A"),
    sprintf("Modalities: %s", report_data$data$n_modalities %||% "1"),
    sprintf("Target: %s", report_data$data$target_name %||% "N/A"),
    "",
    "-" %rep% 40,
    "DISCRIMINATION METRICS",
    "-" %rep% 40,
    sprintf("AUC: %s", format_metric(report_data$discrimination$auc)),
    sprintf("Brier Score: %s", format_metric(report_data$discrimination$brier_score)),
    "",
    "-" %rep% 40,
    "STABILITY METRICS",
    "-" %rep% 40,
    sprintf("Nogueira Stability Index: %s", format_metric(report_data$stability$nogueira_index)),
    sprintf("Interpretation: %s", report_data$stability$interpretation %||% "N/A"),
    sprintf("Consensus Features: %d", length(report_data$stability$consensus_features %||% character(0))),
    "",
    "=" %rep% 60,
    "End of Report",
    "=" %rep% 60
  )

  writeLines(lines, output_file)
  message("Simple report generated: ", output_file)
  invisible(output_file)
}


#' @title Format Metric for Display
#' @keywords internal
format_metric <- function(x, digits = 4) {
  if (is.null(x) || is.na(x)) {
    return("N/A")
  }
  sprintf("%.*f", digits, x)
}


#' @title String Repetition Operator
#' @keywords internal
`%rep%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}


#' @title Null Coalescing Operator
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' @title Print method for ReportData
#' @param x A ReportData object
#' @param ... Additional arguments (ignored)
#' @export
print.ReportData <- function(x, ...) {
  cat("=== OmicSelector Report Data ===\n\n")

  cat("Metadata:\n")
  cat(sprintf("  Timestamp: %s\n", format(x$meta$timestamp)))
  cat(sprintf("  Version: %s\n", x$meta$omicselector_version))

  cat("\nData:\n")
  cat(sprintf("  Samples: %s\n", x$data$n_samples %||% "N/A"))
  cat(sprintf("  Features: %s\n", x$data$n_features %||% "N/A"))
  cat(sprintf("  Modalities: %s\n", x$data$n_modalities %||% "1"))

  cat("\nDiscrimination:\n")
  cat(sprintf("  AUC: %s\n", format_metric(x$discrimination$auc)))
  cat(sprintf("  Brier: %s\n", format_metric(x$discrimination$brier_score)))

  cat("\nStability:\n")
  cat(sprintf("  Nogueira Index: %s\n", format_metric(x$stability$nogueira_index)))

  invisible(x)
}
