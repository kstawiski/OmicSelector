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

        fit <- proj$get_final_fit()
        if (is.null(fit)) {
          showNotification("No fitted model to export", type = "error")
          return()
        }

        bundle <- list(
          fit = fit,
          mapping = proj$get_mapping(),
          summary = proj$summary(),
          created = Sys.time()
        )

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

      # Generate Plumber API template
      api_code <- sprintf('
# OmicSelector Plumber API
# Generated: %s

library(plumber)

# Load model
model <- readRDS("model.rds")

#* @apiTitle %s
#* @apiDescription Biomarker prediction API

#* Predict endpoint
#* @param data:object Input data
#* @post /predict
function(req, data) {
  newdata <- as.data.frame(data)
  pred <- model$fit$learner$predict_newdata(newdata)
  list(
    prediction = pred$response,
    probability = as.list(pred$prob)
  )
}

#* Health check
#* @get /health
function() {
  list(status = "healthy", timestamp = Sys.time())
}
', Sys.time(), input$api_name)

      # Save to project directory
      api_path <- file.path(proj$path, paste0(input$api_name, ".R"))
      writeLines(api_code, api_path)

      showNotification(sprintf("API template saved to: %s", api_path), type = "message")
    })

    # Export ONNX
    observeEvent(input$export_onnx, {
      showNotification("ONNX export is in beta. Please check torch compatibility.", type = "warning")
    })

    # Generate report
    observeEvent(input$generate_report, {
      showNotification("Report generation triggered. This feature is under development.", type = "message")
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
