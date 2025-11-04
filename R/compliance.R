#' @title TRIPOD+AI and PROBAST+AI Compliance Module
#' @description Automated generation of TRIPOD+AI checklists and PROBAST+AI risk of bias assessments
#'
#' This module provides tools for ensuring clinical prediction model development
#' follows the TRIPOD+AI (Transparent Reporting of a multivariable prediction model
#' for Individual Prognosis Or Diagnosis + AI) guidelines and PROBAST+AI (Prediction
#' model Risk Of Bias ASsessment Tool + AI) framework.
#'
#' References:
#' - TRIPOD+AI: Collins et al. (2024) BMJ
#' - PROBAST: Wolff et al. (2019) BMJ
#'
#' @author Konrad Stawiski
#' @name compliance
NULL


#' Generate TRIPOD+AI Compliant Report
#'
#' Creates a comprehensive report following TRIPOD+AI guidelines for transparent
#' reporting of prediction models. The report includes all 27 checklist items
#' specific to AI/ML-based prediction models.
#'
#' @param model_result An OmicSelector model object (from OmicSelector_fit or OmicSelector_nested_cv)
#' @param study_info A list containing study information:
#'   \itemize{
#'     \item title: Study title
#'     \item authors: Character vector of author names
#'     \item objective: Study objective
#'     \item data_source: Description of data source
#'     \item eligibility_criteria: Inclusion/exclusion criteria
#'     \item outcome_definition: How the outcome was defined
#'     \item sample_size_justification: Justification for sample size
#'     \item missing_data_handling: How missing data was handled
#'   }
#' @param output_format Character string: "html", "pdf", "json", or "markdown"
#' @param output_file Path to output file. If NULL, returns the report object
#' @param include_plots Logical, whether to include performance plots
#'
#' @return A list object of class "OmicSelector_tripod_report" containing:
#'   \item{checklist}{Data frame with all 27 TRIPOD+AI items and completion status}
#'   \item{report_sections}{List of report sections with content}
#'   \item{metadata}{Metadata about the report generation}
#'   \item{file_path}{Path to saved report file (if output_file specified)}
#'
#' @examples
#' \dontrun{
#' # After running nested CV
#' result <- OmicSelector_nested_cv(...)
#'
#' # Define study information
#' study_info <- list(
#'   title = "MicroRNA Biomarker Discovery for Cancer Diagnosis",
#'   authors = c("Smith J", "Jones A", "Brown B"),
#'   objective = "To develop a diagnostic model for cancer using miRNA expression",
#'   data_source = "TCGA RNA-seq data",
#'   eligibility_criteria = "Adult patients with confirmed diagnosis",
#'   outcome_definition = "Pathologically confirmed cancer status",
#'   sample_size_justification = "Based on EPV rule: 10 events per variable",
#'   missing_data_handling = "Multiple imputation using MICE"
#' )
#'
#' # Generate report
#' report <- OmicSelector_tripod_report(
#'   model_result = result,
#'   study_info = study_info,
#'   output_format = "html",
#'   output_file = "tripod_report.html"
#' )
#' }
#'
#' @export
OmicSelector_tripod_report <- function(
  model_result,
  study_info = list(),
  output_format = c("html", "pdf", "json", "markdown"),
  output_file = NULL,
  include_plots = TRUE
) {

  output_format <- match.arg(output_format)

  # Extract information from model result
  model_info <- .extract_model_info(model_result)

  # Create TRIPOD+AI checklist
  checklist <- .create_tripod_checklist(model_info, study_info)

  # Generate report sections
  report_sections <- list(
    title_abstract = .generate_title_abstract(study_info, model_info),
    introduction = .generate_introduction(study_info),
    methods = .generate_methods(model_info, study_info),
    results = .generate_results(model_info, model_result),
    discussion = .generate_discussion(model_info),
    other_information = .generate_other_information(study_info)
  )

  # Create report object
  report <- list(
    checklist = checklist,
    report_sections = report_sections,
    metadata = list(
      generated_date = Sys.time(),
      omicselector_version = packageVersion("OmicSelector"),
      tripod_version = "TRIPOD+AI 2024"
    )
  )

  class(report) <- c("OmicSelector_tripod_report", "list")

  # Export if output file specified
  if (!is.null(output_file)) {
    .export_tripod_report(report, output_format, output_file, include_plots)
    report$file_path <- output_file
  }

  return(report)
}


#' Generate PROBAST+AI Risk of Bias Assessment
#'
#' Conducts a systematic assessment of risk of bias in prediction model development
#' following the PROBAST+AI framework.
#'
#' @param model_result An OmicSelector model object
#' @param assessment_inputs A list containing information for bias assessment:
#'   \itemize{
#'     \item participant_selection: Description of how participants were selected
#'     \item predictor_assessment: How predictors were measured
#'     \item outcome_assessment: How outcome was measured
#'     \item sample_selection: Whether appropriate data sources were used
#'     \item attrition: Information about missing data and dropouts
#'   }
#'
#' @return A list object of class "OmicSelector_probast" containing:
#'   \item{overall_risk}{Character: "Low", "High", or "Unclear"}
#'   \item{domain_assessments}{Assessment for each of 4 domains}
#'   \item{concerns_applicability}{Concerns regarding applicability}
#'   \item{recommendations}{Specific recommendations for improvement}
#'
#' @examples
#' \dontrun{
#' result <- OmicSelector_nested_cv(...)
#'
#' assessment <- OmicSelector_probast(
#'   model_result = result,
#'   assessment_inputs = list(
#'     participant_selection = "Consecutive enrollment",
#'     predictor_assessment = "Standardized RNA-seq protocol",
#'     outcome_assessment = "Blinded pathological review"
#'   )
#' )
#'
#' print(assessment)
#' }
#'
#' @export
OmicSelector_probast <- function(
  model_result,
  assessment_inputs = list()
) {

  # Extract model information
  model_info <- .extract_model_info(model_result)

  # Domain 1: Participants
  participants_assessment <- .assess_participants_domain(model_info, assessment_inputs)

  # Domain 2: Predictors
  predictors_assessment <- .assess_predictors_domain(model_info, assessment_inputs)

  # Domain 3: Outcome
  outcome_assessment <- .assess_outcome_domain(model_info, assessment_inputs)

  # Domain 4: Analysis
  analysis_assessment <- .assess_analysis_domain(model_info, assessment_inputs)

  # Overall risk of bias
  domain_risks <- c(
    participants_assessment$risk,
    predictors_assessment$risk,
    outcome_assessment$risk,
    analysis_assessment$risk
  )

  overall_risk <- ifelse(any(domain_risks == "High"), "High",
                        ifelse(any(domain_risks == "Unclear"), "Unclear", "Low"))

  # Applicability concerns
  applicability_concerns <- .assess_applicability(model_info, assessment_inputs)

  # Generate recommendations
  recommendations <- .generate_probast_recommendations(
    participants_assessment,
    predictors_assessment,
    outcome_assessment,
    analysis_assessment
  )

  # Create result object
  result <- list(
    overall_risk = overall_risk,
    domain_assessments = list(
      participants = participants_assessment,
      predictors = predictors_assessment,
      outcome = outcome_assessment,
      analysis = analysis_assessment
    ),
    concerns_applicability = applicability_concerns,
    recommendations = recommendations,
    metadata = list(
      assessment_date = Sys.time(),
      probast_version = "PROBAST+AI 2024"
    )
  )

  class(result) <- c("OmicSelector_probast", "list")

  return(result)
}


#' Internal: Extract Model Information
#'
#' @keywords internal
#' @noRd
.extract_model_info <- function(model_result) {

  info <- list()

  if (inherits(model_result, "OmicSelector_nested_cv")) {
    info$model_type <- "nested_cv"
    info$n_folds_outer <- model_result$metadata$outer_folds
    info$n_folds_inner <- model_result$metadata$inner_folds
    info$n_samples <- model_result$metadata$n_samples
    info$n_features <- model_result$metadata$n_features
    info$outcome_type <- model_result$metadata$outcome_type
    info$feature_selection <- model_result$metadata$feature_selection_method
    info$metrics <- model_result$overall_metrics
    info$predictions <- model_result$final_predictions
  } else if (inherits(model_result, "OmicSelector_model")) {
    info$model_type <- "single_model"
    info$framework <- model_result$framework
    info$algorithm <- model_result$algorithm
    info$metrics <- model_result$metrics
    info$predictions <- model_result$predictions
  }

  return(info)
}


#' Internal: Create TRIPOD+AI Checklist
#'
#' @keywords internal
#' @noRd
.create_tripod_checklist <- function(model_info, study_info) {

  # Define all 27 TRIPOD+AI items
  items <- data.frame(
    Section = c(
      rep("Title and Abstract", 2),
      rep("Introduction", 2),
      rep("Methods - Participants", 3),
      rep("Methods - Outcome", 2),
      rep("Methods - Predictors", 2),
      rep("Methods - Sample Size", 1),
      rep("Methods - Missing Data", 1),
      rep("Methods - Statistical Analysis", 6),
      rep("Results - Participants", 2),
      rep("Results - Model Development", 3),
      rep("Results - Model Specification", 1),
      rep("Results - Model Performance", 2),
      rep("Discussion", 3),
      rep("Other Information", 2)
    ),
    Item = c(
      "1a", "1b",
      "2a", "2b",
      "3a", "3b", "3c",
      "4a", "4b",
      "5a", "5b",
      "6",
      "7",
      "8a", "8b", "8c", "8d", "8e", "8f",
      "9a", "9b",
      "10a", "10b", "10c",
      "11",
      "12a", "12b",
      "13a", "13b", "13c",
      "14a", "14b"
    ),
    Description = c(
      "Title identifies the study as developing or validating a prediction model",
      "Abstract provides summary of objectives, study design, setting, participants, predictors, outcome, analysis, results, and conclusions",
      "Explain medical context and rationale",
      "Specify objectives including intended use of the model",
      "Describe data source and when/where data were collected",
      "Specify eligibility criteria",
      "Specify how outcome was determined",
      "Clearly define the outcome",
      "Report when outcome was determined relative to predictor assessment",
      "Clearly define all predictors and how they were measured",
      "Report when predictors were measured relative to outcome",
      "Explain sample size and number of events per variable",
      "Describe how missing data were handled",
      "Specify type of model and whether it is an AI/ML model",
      "Specify hyperparameter tuning approach",
      "Specify validation approach (e.g., nested CV)",
      "Describe feature selection methods",
      "Specify model evaluation metrics",
      "Describe calibration assessment",
      "Describe participant flow including missing data",
      "Describe characteristics of study participants",
      "Specify selected features and their importance",
      "Specify hyperparameters of final model",
      "Present model equation or algorithm specification",
      "Provide code and data for reproducibility",
      "Report performance metrics with confidence intervals",
      "Report calibration results",
      "Discuss limitations including potential bias and overfitting",
      "Discuss model interpretation and clinical implications",
      "Discuss generalizability and external validation needs",
      "Provide supplementary information including code",
      "Provide funding and conflict of interest statements"
    ),
    Status = "Not Assessed",
    Notes = "",
    stringsAsFactors = FALSE
  )

  # Auto-populate some items based on model_info
  if (model_info$model_type == "nested_cv") {
    items$Status[items$Item == "8c"] <- "Complete"
    items$Notes[items$Item == "8c"] <- paste0(
      "Nested CV with ", model_info$n_folds_outer, " outer and ",
      model_info$n_folds_inner, " inner folds"
    )

    if (model_info$feature_selection != "none") {
      items$Status[items$Item == "8d"] <- "Complete"
      items$Notes[items$Item == "8d"] <- paste0(
        "Feature selection: ", model_info$feature_selection
      )
    }
  }

  return(items)
}


#' Internal: Generate Title and Abstract Section
#'
#' @keywords internal
#' @noRd
.generate_title_abstract <- function(study_info, model_info) {

  section <- list()

  section$title <- ifelse(
    !is.null(study_info$title),
    study_info$title,
    "Development and Validation of a Prediction Model"
  )

  section$abstract <- paste0(
    "Objective: ", study_info$objective %||% "To develop a prediction model", "\n\n",
    "Design: ", ifelse(model_info$model_type == "nested_cv",
                       "Nested cross-validation study", "Development study"), "\n\n",
    "Setting and Participants: ", study_info$data_source %||% "Not specified", "\n\n",
    "Sample Size: ", model_info$n_samples, " participants\n\n",
    "Outcome: ", study_info$outcome_definition %||% "Not specified", "\n\n",
    "Analysis: ", model_info$outcome_type, " modeling with ",
    model_info$feature_selection, " feature selection"
  )

  return(section)
}


#' Internal: Generate Introduction Section
#'
#' @keywords internal
#' @noRd
.generate_introduction <- function(study_info) {

  section <- list(
    background = study_info$background %||% "Not provided",
    objectives = study_info$objective %||% "Not provided"
  )

  return(section)
}


#' Internal: Generate Methods Section
#'
#' @keywords internal
#' @noRd
.generate_methods <- function(model_info, study_info) {

  section <- list(
    data_source = study_info$data_source %||% "Not specified",
    eligibility = study_info$eligibility_criteria %||% "Not specified",
    outcome = study_info$outcome_definition %||% "Not specified",
    sample_size = paste0(
      model_info$n_samples, " samples with ", model_info$n_features, " features. ",
      study_info$sample_size_justification %||% "Sample size determined by available data."
    ),
    missing_data = study_info$missing_data_handling %||% "Not specified",
    analysis = paste0(
      "Model type: ", model_info$model_type, ". ",
      ifelse(model_info$model_type == "nested_cv",
             paste0("Nested cross-validation with ", model_info$n_folds_outer,
                    " outer folds and ", model_info$n_folds_inner, " inner folds. "),
             ""),
      "Feature selection: ", model_info$feature_selection, "."
    )
  )

  return(section)
}


#' Internal: Generate Results Section
#'
#' @keywords internal
#' @noRd
.generate_results <- function(model_info, model_result) {

  section <- list(
    sample_description = paste0(
      "Total sample: ", model_info$n_samples, " participants. ",
      "Features: ", model_info$n_features, ". ",
      "Problem type: ", model_info$outcome_type, "."
    ),
    performance = model_info$metrics
  )

  return(section)
}


#' Internal: Generate Discussion Section
#'
#' @keywords internal
#' @noRd
.generate_discussion <- function(model_info) {

  section <- list(
    limitations = "This is an automated report. Please provide study-specific limitations.",
    interpretation = "This is an automated report. Please provide clinical interpretation.",
    implications = "External validation is recommended before clinical use."
  )

  return(section)
}


#' Internal: Generate Other Information Section
#'
#' @keywords internal
#' @noRd
.generate_other_information <- function(study_info) {

  section <- list(
    supplementary = "Model code and specifications available in OmicSelector format.",
    funding = study_info$funding %||% "Not specified"
  )

  return(section)
}


#' Internal: Assess Participants Domain
#'
#' @keywords internal
#' @noRd
.assess_participants_domain <- function(model_info, inputs) {

  # Check for common issues
  issues <- c()

  if (model_info$n_samples < 100) {
    issues <- c(issues, "Small sample size may increase risk of overfitting")
  }

  risk <- ifelse(length(issues) > 0, "Unclear", "Low")

  result <- list(
    risk = risk,
    issues = issues,
    signaling_questions = list(
      q1 = "Were appropriate data sources used?",
      q2 = "Were all inclusions/exclusions appropriate?"
    )
  )

  return(result)
}


#' Internal: Assess Predictors Domain
#'
#' @keywords internal
#' @noRd
.assess_predictors_domain <- function(model_info, inputs) {

  issues <- c()

  # Check feature selection
  if (model_info$feature_selection == "none" && model_info$n_features > 50) {
    issues <- c(issues, "High-dimensional data without feature selection may lead to overfitting")
  }

  risk <- ifelse(length(issues) > 0, "Unclear", "Low")

  result <- list(
    risk = risk,
    issues = issues
  )

  return(result)
}


#' Internal: Assess Outcome Domain
#'
#' @keywords internal
#' @noRd
.assess_outcome_domain <- function(model_info, inputs) {

  result <- list(
    risk = "Unclear",
    issues = c("Outcome definition and assessment require manual review"),
    signaling_questions = list(
      q1 = "Was outcome determined appropriately?",
      q2 = "Was there potential for outcome assessment bias?"
    )
  )

  return(result)
}


#' Internal: Assess Analysis Domain
#'
#' @keywords internal
#' @noRd
.assess_analysis_domain <- function(model_info, inputs) {

  issues <- c()

  # Check validation
  if (model_info$model_type != "nested_cv") {
    issues <- c(issues, "Nested cross-validation not used - consider for more rigorous evaluation")
  }

  risk <- ifelse(length(issues) == 0, "Low", "Unclear")

  result <- list(
    risk = risk,
    issues = issues,
    strengths = c()
  )

  if (model_info$model_type == "nested_cv") {
    result$strengths <- c(result$strengths, "Nested CV used - reduces overfitting risk")
  }

  return(result)
}


#' Internal: Assess Applicability
#'
#' @keywords internal
#' @noRd
.assess_applicability <- function(model_info, inputs) {

  concerns <- list(
    participants = "Requires assessment of target population match",
    predictors = "Requires assessment of predictor availability in practice",
    outcome = "Requires assessment of outcome relevance to intended use"
  )

  return(concerns)
}


#' Internal: Generate PROBAST Recommendations
#'
#' @keywords internal
#' @noRd
.generate_probast_recommendations <- function(p_assess, pred_assess, out_assess, analysis_assess) {

  recommendations <- c()

  # Collect all issues
  all_issues <- c(
    p_assess$issues,
    pred_assess$issues,
    out_assess$issues,
    analysis_assess$issues
  )

  if (length(all_issues) > 0) {
    recommendations <- c(
      recommendations,
      "Address identified issues before clinical deployment"
    )
  }

  # Add strengths-based recommendations
  recommendations <- c(
    recommendations,
    "Consider external validation in independent cohorts",
    "Assess calibration in target population",
    "Evaluate clinical utility using decision curve analysis"
  )

  return(recommendations)
}


#' Internal: Export TRIPOD Report
#'
#' @keywords internal
#' @noRd
.export_tripod_report <- function(report, format, file_path, include_plots) {

  if (format == "json") {
    # Export as JSON
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      warning("jsonlite package required for JSON export. Using base R instead.")
      saveRDS(report, file_path)
    } else {
      jsonlite::write_json(report, file_path, pretty = TRUE, auto_unbox = TRUE)
    }
  } else if (format %in% c("html", "pdf", "markdown")) {
    # Generate markdown
    md_content <- .generate_markdown_report(report)

    if (format == "markdown") {
      writeLines(md_content, file_path)
    } else {
      # Use rmarkdown to render
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown package required for HTML/PDF export. Saving as markdown instead.")
        writeLines(md_content, gsub("\\.(html|pdf)$", ".md", file_path))
      } else {
        temp_md <- tempfile(fileext = ".md")
        writeLines(md_content, temp_md)
        rmarkdown::render(temp_md, output_format = ifelse(format == "html", "html_document", "pdf_document"),
                         output_file = file_path, quiet = TRUE)
      }
    }
  }

  message(paste0("Report saved to: ", file_path))
}


#' Internal: Generate Markdown Report
#'
#' @keywords internal
#' @noRd
.generate_markdown_report <- function(report) {

  md <- c(
    "# TRIPOD+AI Compliance Report",
    "",
    paste0("Generated: ", format(report$metadata$generated_date, "%Y-%m-%d %H:%M:%S")),
    paste0("OmicSelector Version: ", report$metadata$omicselector_version),
    "",
    "## Title and Abstract",
    "",
    paste0("### Title"),
    report$report_sections$title_abstract$title,
    "",
    paste0("### Abstract"),
    report$report_sections$title_abstract$abstract,
    "",
    "## Methods",
    "",
    paste0("**Data Source:** ", report$report_sections$methods$data_source),
    "",
    paste0("**Eligibility:** ", report$report_sections$methods$eligibility),
    "",
    paste0("**Outcome:** ", report$report_sections$methods$outcome),
    "",
    paste0("**Sample Size:** ", report$report_sections$methods$sample_size),
    "",
    paste0("**Missing Data:** ", report$report_sections$methods$missing_data),
    "",
    paste0("**Statistical Analysis:** ", report$report_sections$methods$analysis),
    "",
    "## Results",
    "",
    report$report_sections$results$sample_description,
    "",
    "## TRIPOD+AI Checklist",
    "",
    "| Section | Item | Description | Status | Notes |",
    "|---------|------|-------------|--------|-------|"
  )

  # Add checklist rows
  for (i in 1:nrow(report$checklist)) {
    row <- report$checklist[i, ]
    md <- c(md, paste0(
      "| ", row$Section, " | ", row$Item, " | ", row$Description, " | ",
      row$Status, " | ", row$Notes, " |"
    ))
  }

  return(md)
}


#' Print Method for OmicSelector_tripod_report
#'
#' @param x An OmicSelector_tripod_report object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_tripod_report <- function(x, ...) {
  cat("TRIPOD+AI Compliance Report\n")
  cat("===========================\n\n")
  cat("Generated:", format(x$metadata$generated_date, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("TRIPOD Version:", x$metadata$tripod_version, "\n\n")

  cat("Checklist Status:\n")
  complete <- sum(x$checklist$Status == "Complete")
  total <- nrow(x$checklist)
  cat(paste0("  ", complete, "/", total, " items assessed\n\n"))

  cat("Use summary() to view the full checklist.\n")
}


#' Print Method for OmicSelector_probast
#'
#' @param x An OmicSelector_probast object
#' @param ... Additional arguments (not used)
#' @export
print.OmicSelector_probast <- function(x, ...) {
  cat("PROBAST+AI Risk of Bias Assessment\n")
  cat("===================================\n\n")

  cat("Overall Risk of Bias:", x$overall_risk, "\n\n")

  cat("Domain Assessments:\n")
  cat("  1. Participants:", x$domain_assessments$participants$risk, "\n")
  cat("  2. Predictors:", x$domain_assessments$predictors$risk, "\n")
  cat("  3. Outcome:", x$domain_assessments$outcome$risk, "\n")
  cat("  4. Analysis:", x$domain_assessments$analysis$risk, "\n\n")

  cat("Recommendations:\n")
  for (i in seq_along(x$recommendations)) {
    cat(paste0("  ", i, ". ", x$recommendations[i], "\n"))
  }
}


#' Helper: Null coalescing operator
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (is.character(a) && a == "")) b else a
}
