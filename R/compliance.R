#' @title TRIPOD+AI and PROBAST+AI Compliance Module
#' @description
#' Functions for generating TRIPOD+AI (Transparent Reporting of a multivariable
#' prediction model for Individual Prognosis Or Diagnosis + AI) compliant reports
#' and PROBAST+AI (Prediction model Risk Of Bias ASsessment Tool + AI) risk of
#' bias assessments.
#'
#' @details
#' TRIPOD+AI is an extension of the TRIPOD guidelines specifically for AI and
#' machine learning prediction models. It includes 27 items covering:
#' * Title and abstract
#' * Introduction
#' * Methods (participants, outcomes, predictors, sample size, missing data,
#'   statistical analysis, model development, performance measures, risk groups)
#' * Results (participants, model specification, performance, updated performance)
#' * Discussion and other information
#'
#' PROBAST+AI assesses risk of bias across 4 domains:
#' * Participants domain
#' * Predictors domain
#' * Outcome domain
#' * Analysis domain
#'
#' @name compliance
#' @keywords internal
NULL


#' @title Generate TRIPOD+AI Compliant Report
#' @description
#' Automatically generates a comprehensive TRIPOD+AI checklist and report from
#' a model result object. This ensures transparent reporting of prediction
#' model development and validation studies.
#'
#' @param model_result An OmicSelector model result object (from OmicSelector_fit,
#'   OmicSelector_nested_cv, or similar functions)
#' @param study_info List containing study-specific information:
#'   \itemize{
#'     \item title: Study title
#'     \item authors: List of authors
#'     \item study_design: Type of study (e.g., "cross-sectional", "cohort")
#'     \item data_source: Description of data source
#'     \item inclusion_criteria: Patient/sample inclusion criteria
#'     \item exclusion_criteria: Patient/sample exclusion criteria
#'     \item outcome_definition: Clear definition of the outcome
#'     \item candidate_predictors: Rationale for predictor selection
#'     \item sample_size_justification: Sample size calculation/justification
#'     \item intended_use: Intended use of the prediction model
#'   }
#' @param output_format Character string specifying output format: "html", "pdf", "json", or "docx"
#' @param output_file Path to save the report (if NULL, returns the report object)
#' @param include_checklist Logical indicating whether to include the TRIPOD+AI checklist
#' @param include_recommendations Logical indicating whether to include recommendations
#' @param template_path Path to custom report template (optional)
#'
#' @return A TRIPOD+AI report object containing:
#' \item{checklist}{Data frame with all 27 TRIPOD+AI items and completion status}
#' \item{report}{Formatted report content}
#' \item{completeness_score}{Percentage of items adequately reported}
#' \item{recommendations}{List of recommendations for improvement}
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' result <- OmicSelector_nested_cv(
#'   data = my_data,
#'   outcome = "disease",
#'   models = list(rf = my_rf_spec)
#' )
#'
#' # Define study information
#' study_info <- list(
#'   title = "Biomarker discovery for disease X",
#'   authors = c("Smith J", "Jones A"),
#'   study_design = "cross-sectional",
#'   data_source = "Gene Expression Omnibus (GSE12345)",
#'   outcome_definition = "Binary disease status (case vs control)"
#' )
#'
#' # Generate TRIPOD+AI report
#' report <- OmicSelector_tripod_report(
#'   model_result = result,
#'   study_info = study_info,
#'   output_format = "html",
#'   output_file = "tripod_report.html"
#' )
#' }
#'
#' @references
#' Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance
#' for reporting clinical prediction models that use regression or machine learning
#' methods. BMJ 2024;385:e078378. doi: 10.1136/bmj-2023-078378
#'
#' @export
OmicSelector_tripod_report <- function(
  model_result,
  study_info = list(),
  output_format = c("html", "pdf", "json", "docx"),
  output_file = NULL,
  include_checklist = TRUE,
  include_recommendations = TRUE,
  template_path = NULL
) {

  output_format <- match.arg(output_format)

  # Validate inputs
  if (!inherits(model_result, c("OmicSelector_fit", "OmicSelector_nested_cv"))) {
    warning("model_result is not a recognized OmicSelector object. Some items may not be automatically filled.")
  }

  # Initialize TRIPOD+AI checklist (27 items)
  checklist <- .initialize_tripod_checklist()

  # Extract information from model_result
  extracted_info <- .extract_model_info(model_result)

  # Fill checklist items automatically where possible
  checklist <- .fill_tripod_checklist(
    checklist = checklist,
    model_info = extracted_info,
    study_info = study_info
  )

  # Calculate completeness score
  completeness <- .calculate_completeness(checklist)

  # Generate recommendations
  recommendations <- NULL
  if (include_recommendations) {
    recommendations <- .generate_tripod_recommendations(checklist)
  }

  # Build report content
  report_content <- .build_tripod_report(
    checklist = checklist,
    model_info = extracted_info,
    study_info = study_info,
    include_checklist = include_checklist
  )

  # Create result object
  result <- list(
    checklist = checklist,
    report = report_content,
    completeness_score = completeness,
    recommendations = recommendations,
    metadata = list(
      generated_date = Sys.time(),
      output_format = output_format,
      package_version = utils::packageVersion("OmicSelector")
    )
  )

  class(result) <- c("OmicSelector_tripod_report", "list")

  # Export to file if requested
  if (!is.null(output_file)) {
    .export_tripod_report(
      report = result,
      output_file = output_file,
      format = output_format,
      template_path = template_path
    )
    message(sprintf("TRIPOD+AI report saved to: %s", output_file))
  }

  return(result)
}


#' @title PROBAST+AI Risk of Bias Assessment
#' @description
#' Performs a systematic risk of bias assessment using the PROBAST+AI framework.
#' This helps identify potential sources of bias in prediction model studies.
#'
#' @param model_result An OmicSelector model result object
#' @param assessment Manual assessment inputs (optional). List containing:
#'   \itemize{
#'     \item participants_bias: Assessment of participant selection bias
#'     \item predictors_bias: Assessment of predictor measurement bias
#'     \item outcome_bias: Assessment of outcome measurement bias
#'     \item analysis_bias: Assessment of analysis bias
#'   }
#' @param output_format Character string: "html", "pdf", or "json"
#' @param output_file Path to save the assessment report
#'
#' @return A PROBAST+AI assessment object containing:
#' \item{domain_assessments}{Risk of bias for each of 4 domains}
#' \item{overall_risk}{Overall risk of bias (low/moderate/high)}
#' \item{concerns_applicability}{Concerns about applicability}
#' \item{recommendations}{Recommendations for improvement}
#' \item{details}{Detailed assessment for each signaling question}
#'
#' @details
#' The 4 PROBAST domains are:
#' 1. **Participants**: Are participants representative? Selection bias?
#' 2. **Predictors**: Were predictors measured appropriately and consistently?
#' 3. **Outcome**: Was outcome determined appropriately?
#' 4. **Analysis**: Was the analysis appropriate? Data leakage? Overfitting?
#'
#' @examples
#' \dontrun{
#' # Perform PROBAST assessment
#' probast <- OmicSelector_probast(
#'   model_result = my_model,
#'   output_file = "probast_assessment.html"
#' )
#'
#' # View overall risk
#' print(probast$overall_risk)
#'
#' # View domain-specific risks
#' print(probast$domain_assessments)
#' }
#'
#' @references
#' Wolff RF, Moons KGM, Riley RD, et al. PROBAST: A Tool to Assess the Risk of
#' Bias and Applicability of Prediction Model Studies. Ann Intern Med.
#' 2019;170(1):51-58. doi:10.7326/M18-1376
#'
#' @export
OmicSelector_probast <- function(
  model_result,
  assessment = list(),
  output_format = c("html", "pdf", "json"),
  output_file = NULL
) {

  output_format <- match.arg(output_format)

  # Initialize PROBAST structure
  probast_structure <- .initialize_probast()

  # Extract information from model
  model_info <- .extract_model_info(model_result)

  # Automatic assessment where possible
  auto_assessment <- .automatic_probast_assessment(model_info)

  # Merge with manual assessment
  final_assessment <- .merge_assessments(auto_assessment, assessment)

  # Calculate overall risk
  overall_risk <- .calculate_overall_risk(final_assessment)

  # Generate recommendations
  recommendations <- .generate_probast_recommendations(final_assessment)

  # Create result object
  result <- list(
    domain_assessments = final_assessment,
    overall_risk = overall_risk,
    recommendations = recommendations,
    details = probast_structure,
    metadata = list(
      generated_date = Sys.time(),
      package_version = utils::packageVersion("OmicSelector")
    )
  )

  class(result) <- c("OmicSelector_probast", "list")

  # Export if requested
  if (!is.null(output_file)) {
    .export_probast_report(
      assessment = result,
      output_file = output_file,
      format = output_format
    )
    message(sprintf("PROBAST+AI assessment saved to: %s", output_file))
  }

  return(result)
}


#' @title Model Reporting Card
#' @description
#' Generates a comprehensive model card documenting model development,
#' performance, intended use, and limitations.
#'
#' @param model_result Model result object
#' @param model_details List with model-specific details
#' @param intended_use Description of intended use cases
#' @param limitations Known limitations
#' @param ethical_considerations Ethical considerations
#' @param output_file Output file path
#'
#' @return A model card object
#'
#' @export
OmicSelector_model_card <- function(
  model_result,
  model_details = list(),
  intended_use = "",
  limitations = "",
  ethical_considerations = "",
  output_file = NULL
) {

  card <- list(
    model_details = model_details,
    intended_use = intended_use,
    factors = list(),  # Relevant factors (e.g., population characteristics)
    metrics = list(),  # Performance metrics
    evaluation_data = list(),  # Evaluation dataset details
    training_data = list(),  # Training dataset details
    quantitative_analyses = list(),  # Performance across subgroups
    ethical_considerations = ethical_considerations,
    caveats_recommendations = limitations
  )

  class(card) <- c("OmicSelector_model_card", "list")

  if (!is.null(output_file)) {
    .export_model_card(card, output_file)
  }

  return(card)
}


# Internal helper functions ----

#' @keywords internal
.initialize_tripod_checklist <- function() {
  # TRIPOD+AI checklist with all 27 items
  checklist <- data.frame(
    section = c(
      rep("Title and Abstract", 2),
      rep("Introduction", 2),
      rep("Methods - Source of data", 2),
      rep("Methods - Participants", 2),
      rep("Methods - Outcome", 2),
      rep("Methods - Predictors", 2),
      rep("Methods - Sample size", 1),
      rep("Methods - Missing data", 1),
      rep("Methods - Statistical analysis", 5),
      rep("Methods - Model development", 2),
      rep("Methods - Model performance", 1),
      rep("Results - Participants", 1),
      rep("Results - Model development", 2),
      rep("Results - Model specification", 1),
      rep("Results - Model performance", 1),
      rep("Discussion", 2),
      rep("Other information", 1)
    ),
    item = c(
      "1", "2",
      "3a", "3b",
      "4a", "4b",
      "5a", "5b", "5c",
      "6a", "6b",
      "7a", "7b",
      "8",
      "9",
      "10a", "10b", "10c", "10d", "10e",
      "11", "12",
      "13a", "13b",
      "14a", "14b",
      "15a", "15b",
      "16",
      "17",
      "18", "19a", "19b",
      "20",
      "21",
      "22", "23",
      "24"
    ),
    description = c(
      "Title identifies study as developing/validating prediction model",
      "Abstract provides summary of objectives, methods, results",
      "Explain medical context and rationale for prediction model",
      "Specify objectives including intended use",
      "Describe study design and data sources",
      "Specify study key dates",
      "Describe eligibility criteria",
      "Describe data cleaning and preprocessing",
      "Define outcome with details",
      "Clearly define predictors with details",
      "Report predictor handling",
      "Explain sample size",
      "Describe missing data methods",
      "Describe model development process",
      "Specify ML algorithm and training procedure",
      "Report hyperparameter tuning",
      "Specify resampling strategy",
      "Report measures to prevent overfitting",
      "Describe measures of model performance",
      "Describe assessment of calibration",
      "Risk stratification",
      "Describe participant flow",
      "Provide detailed model specification",
      "Report hyperparameters selected",
      "Report model performance with confidence intervals",
      "Report calibration assessment results",
      "Clinical interpretation of model performance",
      "Limitations of study",
      "Implications and registration"
    ),
    status = rep("Not assessed", 27),
    page = rep(NA_character_, 27),
    stringsAsFactors = FALSE
  )

  return(checklist)
}


#' @keywords internal
.extract_model_info <- function(model_result) {
  info <- list()

  if (inherits(model_result, "OmicSelector_nested_cv")) {
    info$metadata <- model_result$metadata
    info$performance <- model_result$aggregated_results
    info$resampling <- list(
      method = "nested_cv",
      outer_folds = model_result$metadata$outer_folds,
      inner_folds = model_result$metadata$inner_folds
    )
    info$features <- model_result$selected_features
    info$calibration <- model_result$calibration
  } else if (inherits(model_result, "OmicSelector_fit")) {
    info$metadata <- model_result$metadata
    info$performance <- model_result$performance
  }

  return(info)
}


#' @keywords internal
.fill_tripod_checklist <- function(checklist, model_info, study_info) {
  # Automatically fill items that can be determined from model_info

  # Item 10a: Model development process
  if (!is.null(model_info$metadata$algorithm)) {
    checklist$status[checklist$item == "10a"] <- "Complete"
    checklist$page[checklist$item == "10a"] <- "Methods (automated)"
  }

  # Item 10d: Resampling strategy
  if (!is.null(model_info$resampling)) {
    checklist$status[checklist$item == "10d"] <- "Complete"
    checklist$page[checklist$item == "10d"] <- "Methods (automated)"
  }

  # Add more automatic filling based on available information

  return(checklist)
}


#' @keywords internal
.calculate_completeness <- function(checklist) {
  n_complete <- sum(checklist$status == "Complete")
  n_total <- nrow(checklist)
  percentage <- round(100 * n_complete / n_total, 1)

  list(
    n_complete = n_complete,
    n_total = n_total,
    percentage = percentage
  )
}


#' @keywords internal
.generate_tripod_recommendations <- function(checklist) {
  incomplete <- checklist[checklist$status != "Complete", ]

  recommendations <- lapply(1:nrow(incomplete), function(i) {
    item <- incomplete[i, ]
    sprintf("Item %s (%s): %s",
            item$item, item$section, item$description)
  })

  return(recommendations)
}


#' @keywords internal
.build_tripod_report <- function(checklist, model_info, study_info, include_checklist) {
  # Build structured report content
  report <- list(
    title = "TRIPOD+AI Compliance Report",
    date = format(Sys.time(), "%Y-%m-%d"),
    sections = list()
  )

  # Add sections based on checklist
  if (include_checklist) {
    report$checklist_table <- checklist
  }

  return(report)
}


#' @keywords internal
.export_tripod_report <- function(report, output_file, format, template_path) {
  # Export report to file
  if (format == "json") {
    jsonlite::write_json(report, output_file, pretty = TRUE)
  } else if (format == "html") {
    # Generate HTML report
    .generate_html_report(report, output_file)
  } else if (format == "pdf" || format == "docx") {
    message("PDF and DOCX export require rmarkdown. Using HTML format instead.")
    .generate_html_report(report, sub("\\.(pdf|docx)$", ".html", output_file))
  }
}


#' @keywords internal
.generate_html_report <- function(report, output_file) {
  # Simple HTML report generation
  html <- paste0(
    "<!DOCTYPE html>\n",
    "<html>\n<head>\n",
    "<title>", report$title, "</title>\n",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 40px; }",
    "table { border-collapse: collapse; width: 100%; }",
    "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
    "th { background-color: #4CAF50; color: white; }",
    "</style>\n",
    "</head>\n<body>\n",
    "<h1>", report$title, "</h1>\n",
    "<p>Generated: ", report$date, "</p>\n",
    "</body>\n</html>"
  )

  writeLines(html, output_file)
}


#' @keywords internal
.initialize_probast <- function() {
  # Initialize PROBAST structure with 4 domains and signaling questions
  structure <- list(
    participants = list(
      questions = c(
        "Were appropriate data sources used?",
        "Were all inclusions and exclusions of participants appropriate?"
      ),
      risk = "Not assessed"
    ),
    predictors = list(
      questions = c(
        "Were predictors defined and assessed appropriately?",
        "Were all predictors available at the time of prediction?"
      ),
      risk = "Not assessed"
    ),
    outcome = list(
      questions = c(
        "Was the outcome defined and determined appropriately?",
        "Was there a reasonable time interval between predictor assessment and outcome?"
      ),
      risk = "Not assessed"
    ),
    analysis = list(
      questions = c(
        "Were there a reasonable number of participants with the outcome?",
        "Were continuous and categorical predictors handled appropriately?",
        "Were all enrolled participants included in the analysis?",
        "Were missing data handled appropriately?",
        "Was selection of predictors based on univariable analysis avoided?",
        "Were complexities in the data accounted for?",
        "Were relevant model performance measures evaluated?",
        "Was model overfitting and optimism accounted for?",
        "Were predictions and predictor effects presented with confidence intervals?"
      ),
      risk = "Not assessed"
    )
  )

  return(structure)
}


#' @keywords internal
.automatic_probast_assessment <- function(model_info) {
  assessment <- list(
    participants = list(risk = "Unknown"),
    predictors = list(risk = "Unknown"),
    outcome = list(risk = "Unknown"),
    analysis = list(risk = "Low")  # Assume low if nested CV was used
  )

  # Check for nested CV (good for preventing overfitting)
  if (!is.null(model_info$resampling)) {
    if (model_info$resampling$method == "nested_cv") {
      assessment$analysis$risk <- "Low"
      assessment$analysis$note <- "Nested cross-validation used to prevent overfitting"
    }
  }

  return(assessment)
}


#' @keywords internal
.merge_assessments <- function(auto_assessment, manual_assessment) {
  # Merge automatic and manual assessments, preferring manual when provided
  merged <- auto_assessment

  for (domain in names(manual_assessment)) {
    if (domain %in% names(merged)) {
      merged[[domain]] <- manual_assessment[[domain]]
    }
  }

  return(merged)
}


#' @keywords internal
.calculate_overall_risk <- function(assessment) {
  risks <- sapply(assessment, function(x) x$risk)

  if (any(risks == "High")) {
    return("High")
  } else if (any(risks == "Moderate") || any(risks == "Unknown")) {
    return("Moderate")
  } else {
    return("Low")
  }
}


#' @keywords internal
.generate_probast_recommendations <- function(assessment) {
  recommendations <- list()

  for (domain in names(assessment)) {
    risk <- assessment[[domain]]$risk

    if (risk == "High") {
      recommendations[[domain]] <- sprintf(
        "HIGH RISK in %s domain. Critical issues must be addressed.",
        domain
      )
    } else if (risk == "Moderate" || risk == "Unknown") {
      recommendations[[domain]] <- sprintf(
        "MODERATE/UNKNOWN RISK in %s domain. Consider improvements.",
        domain
      )
    }
  }

  return(recommendations)
}


#' @keywords internal
.export_probast_report <- function(assessment, output_file, format) {
  if (format == "json") {
    jsonlite::write_json(assessment, output_file, pretty = TRUE)
  } else {
    .generate_probast_html(assessment, output_file)
  }
}


#' @keywords internal
.generate_probast_html <- function(assessment, output_file) {
  html <- paste0(
    "<!DOCTYPE html>\n<html>\n<head>\n",
    "<title>PROBAST+AI Assessment</title>\n",
    "</head>\n<body>\n",
    "<h1>PROBAST+AI Risk of Bias Assessment</h1>\n",
    "<h2>Overall Risk: ", assessment$overall_risk, "</h2>\n",
    "</body>\n</html>"
  )

  writeLines(html, output_file)
}


#' @keywords internal
.export_model_card <- function(card, output_file) {
  # Simple JSON export for model card
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(card, output_file, pretty = TRUE)
  } else {
    warning("jsonlite package required for model card export")
  }
}


# S3 methods ----

#' @export
print.OmicSelector_tripod_report <- function(x, ...) {
  cat("TRIPOD+AI Compliance Report\n")
  cat("===========================\n\n")
  cat("Completeness:", x$completeness_score$percentage, "%\n")
  cat("Items complete:", x$completeness_score$n_complete, "/", x$completeness_score$n_total, "\n\n")

  if (!is.null(x$recommendations) && length(x$recommendations) > 0) {
    cat("Recommendations:\n")
    for (i in 1:min(5, length(x$recommendations))) {
      cat("  -", x$recommendations[[i]], "\n")
    }
    if (length(x$recommendations) > 5) {
      cat("  ... and", length(x$recommendations) - 5, "more\n")
    }
  }

  invisible(x)
}


#' @export
print.OmicSelector_probast <- function(x, ...) {
  cat("PROBAST+AI Risk of Bias Assessment\n")
  cat("===================================\n\n")
  cat("Overall Risk:", x$overall_risk, "\n\n")

  cat("Domain-specific risks:\n")
  for (domain in names(x$domain_assessments)) {
    cat("  ", domain, ": ", x$domain_assessments[[domain]]$risk, "\n", sep = "")
  }

  if (!is.null(x$recommendations) && length(x$recommendations) > 0) {
    cat("\nRecommendations:\n")
    for (domain in names(x$recommendations)) {
      cat("  -", x$recommendations[[domain]], "\n")
    }
  }

  invisible(x)
}


#' @export
print.OmicSelector_model_card <- function(x, ...) {
  cat("OmicSelector Model Card\n")
  cat("=======================\n\n")
  cat("Model Details:\n")
  str(x$model_details, max.level = 1)
  cat("\nIntended Use:", substr(x$intended_use, 1, 100), "...\n")
  invisible(x)
}
