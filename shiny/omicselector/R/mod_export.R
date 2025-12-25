#' Export Module
#'
#' @description Model export and deployment preparation

# Module UI
mod_export_ui <- function(id) {
  ns <- NS(id)

  tagList(
    div(
      class = "step-header mb-4",
      h3(bsicons::bs_icon("box-arrow-up-right"), "Step 6: Export & Deploy"),
      p(class = "text-muted", "Export your model for deployment and generate reports")
    ),

    # Export Readiness Card
    card(
      class = "mb-4",
      card_header("Export Readiness"),
      card_body(
        uiOutput(ns("readiness_checks"))
      )
    ),

    # Export Options
    layout_column_wrap(
      width = 1/3,
      fill = FALSE,

      # RDS Bundle
      card(
        card_header(
          bsicons::bs_icon("file-earmark-binary"),
          span("RDS Bundle", class = "ms-2")
        ),
        card_body(
          p("Complete R object including fitted model, feature list, and preprocessing recipe."),
          textInput(ns("rds_name"), "Filename", value = "omic_model"),
          downloadButton(ns("download_rds"), "Download RDS", class = "btn-primary w-100")
        )
      ),

      # API Export
      card(
        card_header(
          bsicons::bs_icon("cloud-arrow-up"),
          span("Plumber API", class = "ms-2")
        ),
        card_body(
          p("Generate Plumber API template for RESTful deployment."),
          textInput(ns("api_name"), "API Name", value = "predict_api"),
          actionButton(ns("generate_api"), "Generate API Template",
                       class = "btn-secondary w-100",
                       icon = icon("server"))
        )
      ),

      # ONNX Export
      card(
        card_header(
          bsicons::bs_icon("diagram-2"),
          span("ONNX (Beta)", class = "ms-2")
        ),
        card_body(
          p("Export to ONNX format for cross-platform deployment."),
          actionButton(ns("export_onnx"), "Export ONNX",
                       class = "btn-outline-secondary w-100",
                       icon = icon("share-nodes"))
        )
      )
    ),

    # Report Generation
    card(
      class = "mt-4",
      card_header(
        bsicons::bs_icon("file-earmark-text"),
        span("Report Generation", class = "ms-2")
      ),
      card_body(
        layout_column_wrap(
          width = 1/2,

          div(
            h6("TRIPOD-Compliant Report"),
            p(class = "text-muted small",
              "Generate a comprehensive report following TRIPOD guidelines for transparent reporting of prediction models."),
            textInput(ns("report_title"), "Report Title", value = "Biomarker Discovery Report"),
            textInput(ns("report_author"), "Author"),
            actionButton(ns("generate_report"), "Generate Report",
                         class = "btn-success",
                         icon = icon("file-pdf"))
          ),

          div(
            h6("Report Sections"),
            checkboxGroupInput(ns("report_sections"), NULL,
                               choices = c(
                                 "Data Summary" = "data",
                                 "QC Results" = "qc",
                                 "Pipeline Configuration" = "pipeline",
                                 "Benchmark Results" = "benchmark",
                                 "Stability Analysis" = "stability",
                                 "Feature Importance" = "importance",
                                 "Model Coefficients" = "coefficients"
                               ),
                               selected = c("data", "qc", "pipeline", "benchmark", "stability"))
          )
        )
      )
    ),

    # Project Summary
    card(
      class = "mt-4",
      card_header("Project Summary"),
      card_body(
        verbatimTextOutput(ns("project_summary")),
        hr(),
        layout_column_wrap(
          width = 1/2,
          downloadButton(ns("download_project"), "Download Project (.rds)",
                         class = "btn-outline-primary"),
          actionButton(ns("save_project"), "Save Project",
                       class = "btn-primary",
                       icon = icon("save"))
        )
      )
    ),

    # Audit Log
    card(
      class = "mt-4",
      card_header(
        bsicons::bs_icon("journal-text"),
        span("Audit Log", class = "ms-2")
      ),
      card_body(
        div(style = "max-height: 200px; overflow-y: auto; font-family: monospace; font-size: 0.85em;",
            verbatimTextOutput(ns("audit_log")))
      )
    )
  )
}

# Module Server
mod_export_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Readiness checks
    output$readiness_checks <- renderUI({
      proj <- project()
      if (is.null(proj)) {
        return(div(class = "text-muted", "No project loaded"))
      }

      summary <- proj$summary()

      checks <- list(
        list(
          name = "Data Loaded",
          pass = !is.na(summary$n_samples),
          detail = if (!is.na(summary$n_samples)) sprintf("%d samples", summary$n_samples) else "Not loaded"
        ),
        list(
          name = "Pipeline Created",
          pass = summary$has_pipeline,
          detail = if (summary$has_pipeline) "Configured" else "Not created"
        ),
        list(
          name = "Benchmark Complete",
          pass = summary$has_benchmark,
          detail = if (summary$has_benchmark) "Complete" else "Not run"
        ),
        list(
          name = "Final Model Fitted",
          pass = summary$has_final_fit,
          detail = if (summary$has_final_fit) "Fitted" else "Not fitted"
        )
      )

      check_items <- lapply(checks, function(check) {
        icon_class <- if (check$pass) "text-success" else "text-danger"
        icon_name <- if (check$pass) "check-circle-fill" else "x-circle-fill"

        div(
          class = "d-flex align-items-center mb-2",
          bsicons::bs_icon(icon_name, class = paste(icon_class, "me-2")),
          strong(check$name),
          span(class = "ms-auto text-muted small", check$detail)
        )
      })

      do.call(tagList, check_items)
    })

    # Project summary
    output$project_summary <- renderPrint({
      proj <- project()
      if (is.null(proj)) {
        cat("No project loaded")
        return(invisible(NULL))
      }

      summary <- proj$summary()
      cat("Project:", summary$name, "\n")
      cat("ID:", summary$id, "\n")
      cat("Path:", summary$path, "\n")
      cat("Created:", format(summary$created_at), "\n")
      cat("Modified:", format(summary$modified_at), "\n")
      cat("\n")
      cat("Samples:", summary$n_samples, "\n")
      cat("Features:", summary$n_features, "\n")
      cat("Multi-omics:", summary$is_multi_omics, "\n")
      cat("Target:", summary$target, "\n")
      cat("\n")
      cat("Pipeline:", if (summary$has_pipeline) "Yes" else "No", "\n")
      cat("Benchmark:", if (summary$has_benchmark) "Yes" else "No", "\n")
      cat("Final Fit:", if (summary$has_final_fit) "Yes" else "No", "\n")
      cat("QC Status:", summary$qc_overall, "\n")
    })

    # Audit log
    output$audit_log <- renderPrint({
      proj <- project()
      if (is.null(proj)) {
        cat("No project loaded")
        return(invisible(NULL))
      }

      log <- proj$get_audit_log()
      if (length(log) == 0) {
        cat("No events logged")
      } else {
        cat(paste(log, collapse = "\n"))
      }
    })

    # Download RDS
    output$download_rds <- downloadHandler(
      filename = function() {
        paste0(input$rds_name, ".rds")
      },
      content = function(file) {
        proj <- project()
        if (is.null(proj)) return()

        final_fit <- proj$get_final_fit()
        if (is.null(final_fit)) {
          showNotification("No fitted model to export", type = "error")
          return()
        }

        # Handle both old (direct fit) and new (bundle) structures
        if (is.list(final_fit) && !is.null(final_fit$fit)) {
          # New bundle structure - include consensus features
          bundle <- list(
            fit = final_fit$fit,
            learner = final_fit$learner,
            consensus_features = final_fit$consensus_features,
            learner_id = final_fit$learner_id,
            stability = final_fit$stability,
            auc = final_fit$auc,
            mapping = proj$get_mapping(),
            summary = proj$summary(),
            created = Sys.time()
          )
        } else {
          # Legacy structure
          bundle <- list(
            fit = final_fit,
            mapping = proj$get_mapping(),
            summary = proj$summary(),
            created = Sys.time()
          )
        }

        saveRDS(bundle, file)
      }
    )

    # Generate API
    observeEvent(input$generate_api, {
      proj <- project()
      req(proj)

      fit <- proj$get_final_fit()
      if (is.null(fit)) {
        showNotification("Train final model first", type = "error")
        return()
      }

      # Generate Plumber API template with robust learner handling
      api_code <- sprintf('
# OmicSelector Plumber API
# Generated: %s

library(plumber)

# Load model bundle
# Expects the RDS file to be named 'model.rds' or passed as an environment variable
model_path <- Sys.getenv("OMICSELECTOR_MODEL_PATH", "model.rds")
if (!file.exists(model_path)) {
  stop(sprintf("Model file '%s' not found. Please place the exported RDS file in this directory and rename it to 'model.rds' or set OMICSELECTOR_MODEL_PATH.", model_path))
}
bundle <- readRDS(model_path)

# Extract the fitted model - handle different structures
get_learner <- function(bundle) {
  if (!is.null(bundle$fit$learner)) {
    return(bundle$fit$learner)
  } else if (!is.null(bundle$fit$model)) {
    return(bundle$fit$model)
  } else if (inherits(bundle$fit, "Learner")) {
    return(bundle$fit)
  } else if (!is.null(bundle$learner)) {
    return(bundle$learner)
  } else {
    stop("Could not extract learner from bundle")
  }
}

learner <- get_learner(bundle)

#* @apiTitle %s
#* @apiDescription Biomarker prediction API generated by OmicSelector

#* Predict endpoint
#* @param data:object Input data as JSON object or array
#* @post /predict
function(req, data) {
  tryCatch({
    newdata <- as.data.frame(data)

    # Handle prediction with error checking
    if (!is.null(learner$predict_newdata)) {
      pred <- learner$predict_newdata(newdata)
      list(
        success = TRUE,
        prediction = as.character(pred$response),
        probability = if (!is.null(pred$prob)) as.list(as.data.frame(pred$prob)) else NULL
      )
    } else if (!is.null(learner$predict)) {
      pred <- learner$predict(newdata)
      list(
        success = TRUE,
        prediction = as.character(pred)
      )
    } else {
      list(success = FALSE, error = "Learner does not support prediction")
    }
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

#* Health check
#* @get /health
function() {
  list(
    status = "healthy",
    timestamp = as.character(Sys.time()),
    model_type = class(learner)[1]
  )
}

#* Model info
#* @get /info
function() {
  list(
    mapping = bundle$mapping,
    created = as.character(bundle$created),
    summary = bundle$summary
  )
}
', Sys.time(), input$api_name)

      # Save to project directory
      # Security: Sanitize filename to prevent path traversal attacks
      clean_name <- gsub("[^a-zA-Z0-9_-]", "", input$api_name)
      if (!nzchar(clean_name)) clean_name <- "predict_api"
      api_path <- file.path(proj$path, paste0(clean_name, ".R"))

      # P0 Fix: Add error handling for file write operations
      tryCatch({
        # Use atomic write pattern: write to temp file first, then rename
        temp_path <- paste0(api_path, ".tmp")
        writeLines(api_code, temp_path)
        if (file.exists(temp_path)) {
          file.rename(temp_path, api_path)
          showNotification(sprintf("API template saved to: %s", api_path), type = "message")
        } else {
          stop("Failed to write temporary file")
        }
      }, error = function(e) {
        showNotification(
          sprintf("Failed to write API template: %s. Check directory permissions.", e$message),
          type = "error",
          duration = 8
        )
      })
    })

    # Export ONNX
    observeEvent(input$export_onnx, {
      proj <- project()
      req(proj)

      final_fit <- proj$get_final_fit()
      if (is.null(final_fit)) {
        showNotification("Train final model first", type = "error")
        return()
      }

      # Extract the actual fit from the bundle (new structure)
      learner <- if (is.list(final_fit) && !is.null(final_fit$fit)) {
        final_fit$fit
      } else {
        final_fit
      }

      # Get the task from the pipeline
      pipeline <- proj$get_pipeline()
      task <- if (!is.null(pipeline) && !is.null(pipeline$task)) {
        pipeline$task
      } else {
        NULL
      }

      if (is.null(task)) {
        showNotification("ONNX export requires a task object. Pipeline may not support this.", type = "warning")
        # Continue anyway - export_onnx may handle NULL task
      }

      # Create output directory if needed
      onnx_dir <- file.path(proj$path, "onnx_export")
      if (!dir.exists(onnx_dir)) dir.create(onnx_dir)

      # Determine basic filename
      base_name <- paste0("omic_model_", format(Sys.time(), "%Y%m%d_%H%M%S"))
      onnx_path <- file.path(onnx_dir, paste0(base_name, ".onnx"))

      withProgress(message = "Exporting to ONNX...", value = 0.5, {
        tryCatch({
          # Check if export_onnx function exists
          if (!exists("export_onnx", mode = "function")) {
            showNotification("ONNX export function not available. Install required packages.", type = "error")
            return()
          }

          # Call backend export function
          result <- export_onnx(
            learner = learner,
            path = onnx_path,
            task = task,
            export_preprocessing = TRUE
          )

          if (result$success) {
            showNotification(sprintf("Export successful! Saved to: %s", onnx_dir), type = "message")
          } else {
            showNotification(paste("Export warning:", result$message), type = "warning")
          }
        }, error = function(e) {
          showNotification(paste("Export failed:", e$message), type = "error")
        })
      })
    })

    # Generate report
    observeEvent(input$generate_report, {
      proj <- project()
      req(proj)

      # Check requirements
      if (!proj$can_proceed("results")) {
        showNotification("Complete benchmark and results steps first", type = "warning")
        return()
      }

      # P0 Fix: Check if report generation functions are available
      if (!exists("create_report_data", mode = "function") ||
          !exists("generate_tripod_report", mode = "function")) {
        showNotification(
          "Report generation is not available in this build. The functions create_report_data() and generate_tripod_report() are not defined. Please ensure the OmicSelector package is properly installed with report generation support.",
          type = "error",
          duration = 10
        )
        return()
      }

      report_dir <- file.path(proj$path, "reports")
      if (!dir.exists(report_dir)) dir.create(report_dir)

      output_file <- file.path(report_dir, paste0("report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"))

      withProgress(message = "Generating TRIPOD+AI Report...", value = 0.2, {
        tryCatch({
          # Gather data
          incProgress(0.2, detail = "Gathering data...")

          # Create ReportData object
          report_data <- create_report_data(
            benchmark_result = proj$get_benchmark_result(),
            stability = if(!is.null(proj$get_benchmark_result())) proj$get_benchmark_result()$stability else NULL,
            pipeline = proj$get_pipeline(),
            config = list(
              project_name = proj$summary()$name,
              project_id = proj$summary()$id
            )
          )

          incProgress(0.4, detail = "Rendering HTML...")

          # Generate report
          report_path <- generate_tripod_report(
            results = report_data,
            output_file = output_file,
            format = "html"
          )

          incProgress(0.2, detail = "Done!")
          showNotification(sprintf("Report generated: %s", basename(report_path)), type = "message")

          # Optionally open (local)
          # utils::browseURL(report_path)

        }, error = function(e) {
          showNotification(paste("Report generation failed:", e$message), type = "error")
        })
      })
    })

    # Save project
    observeEvent(input$save_project, {
      proj <- project()
      req(proj)

      tryCatch({
        proj$save()
        showNotification("Project saved!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })

    # Download project
    output$download_project <- downloadHandler(
      filename = function() {
        proj <- project()
        if (is.null(proj)) return("project.rds")
        paste0(proj$id, "_project.rds")
      },
      content = function(file) {
        proj <- project()
        if (!is.null(proj)) {
          proj$save(file)
        }
      }
    )
  })
}
