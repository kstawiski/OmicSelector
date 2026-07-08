#' Benchmark Module
#'
#' @description Run nested cross-validation with progress tracking

# Module UI
mod_benchmark_ui <- function(id) {
  ns <- NS(id)

  tagList(
    div(
      class = "step-header mb-4",
      h3(bsicons::bs_icon("speedometer2"), "Step 4: Benchmark"),
      p(class = "text-muted", "Run nested cross-validation for unbiased performance estimation")
    ),

    layout_columns(
      col_widths = c(4, 8),

      # Configuration Card
      card(
        card_header("Benchmark Configuration"),
        card_body(
          h6("Cross-Validation Settings"),
          numericInput(ns("outer_folds"), "Outer Folds (Evaluation)", value = 5, min = 2, max = 20),
          numericInput(ns("inner_folds"), "Inner Folds (Tuning)", value = 3, min = 2, max = 10),
          checkboxInput(ns("stratify"), "Stratified CV", value = TRUE),
          numericInput(ns("seed"), "Random Seed", value = 42, min = 1),

          hr(),

          h6("Execution Settings"),
          checkboxInput(ns("parallel"), "Parallel Execution", value = TRUE),
          conditionalPanel(
            condition = sprintf("input['%s']", ns("parallel")),
            uiOutput(ns("workers_ui"))
          ),

          hr(),

          actionButton(ns("run_benchmark"), "Start Benchmark",
                       class = "btn-danger btn-lg w-100",
                       icon = icon("play")),

          div(class = "mt-3",
              uiOutput(ns("benchmark_status_badge")))
        )
      ),

      # Progress & Results Card
      card(
        card_header("Benchmark Progress"),
        card_body(
          uiOutput(ns("progress_ui")),

          hr(),

          h5("Performance Summary"),
          tableOutput(ns("performance_table")),

          hr(),

          h5("Stability Analysis"),
          verbatimTextOutput(ns("stability_summary")),

          hr(),

          actionButton(ns("next_step"), "View Results",
                       class = "btn-success btn-lg",
                       icon = icon("arrow-right"))
        )
      )
    ),

    # Information Card
    card(
      class = "mt-4",
      card_header(
        bsicons::bs_icon("info-circle"),
        span("About Nested Cross-Validation", class = "ms-2")
      ),
      card_body(
        div(
          class = "row",
          div(
            class = "col-md-6",
            h6("Why Nested CV?"),
            p("Nested cross-validation provides unbiased performance estimates by separating:"),
            tags$ul(
              tags$li("Outer loop: Model evaluation (test performance)"),
              tags$li("Inner loop: Hyperparameter tuning / feature selection")
            ),
            p("This prevents data leakage and overly optimistic results.")
          ),
          div(
            class = "col-md-6",
            h6("Stability Index"),
            p("The Nogueira Stability Index measures how consistently features are selected across folds:"),
            tags$ul(
              tags$li("1.0 = Perfect stability (same features every fold)"),
              tags$li("0.0 = No stability (random selection)"),
              tags$li("> 0.7 = Good stability (recommended)")
            )
          )
        )
      )
    )
  )
}

# Module Server
mod_benchmark_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Detect available cores for worker validation
    available_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(available_cores) || available_cores < 1) available_cores <- 4
    recommended_cores <- max(1, available_cores - 1)  # Leave one core free

    # Workers UI with system-aware validation
    output$workers_ui <- renderUI({
      tagList(
        numericInput(ns("workers"), "Workers",
                     value = min(2, recommended_cores),
                     min = 1,
                     max = available_cores),
        p(class = "text-muted small",
          sprintf("Detected %d cores. Recommended: %d (leaves 1 free)",
                  available_cores, recommended_cores))
      )
    })

    # Benchmark status
    benchmark_running <- reactiveVal(FALSE)

    # Status badge
    output$benchmark_status_badge <- renderUI({
      proj <- project()
      if (is.null(proj)) return(NULL)

      result <- proj$get_benchmark_result()

      if (benchmark_running()) {
        div(
          class = "alert alert-warning d-flex align-items-center",
          div(class = "spinner-border spinner-border-sm me-2"),
          "Benchmark running..."
        )
      } else if (!is.null(result)) {
        div(
          class = "alert alert-success",
          bsicons::bs_icon("check-circle"),
          " Benchmark complete!"
        )
      } else {
        div(
          class = "alert alert-info",
          "Ready to run benchmark"
        )
      }
    })

    # Progress UI
    output$progress_ui <- renderUI({
      if (benchmark_running()) {
        div(
          class = "progress-container",
          div(class = "progress mb-2",
              div(class = "progress-bar progress-bar-striped progress-bar-animated",
                  role = "progressbar",
                  style = "width: 100%")),
          p(class = "text-muted text-center", "Running nested cross-validation...")
        )
      } else {
        NULL
      }
    })

    # Run benchmark
    observeEvent(input$run_benchmark, {
      proj <- project()
      req(proj)

      pipeline <- proj$get_pipeline()
      learners <- proj$get_learner()

      if (is.null(pipeline) || is.null(learners)) {
        showNotification("Create pipeline first", type = "error")
        return()
      }

      # Handle both single learner and list of learners
      # Check if it's a single Learner object (not a list of learners)
      if (inherits(learners, "Learner")) {
        # Single learner object - wrap in list
        learners <- list(learners)
      } else if (is.list(learners) && length(learners) > 0) {
        # Already a list - check if first element is a Learner to validate
        if (!inherits(learners[[1]], "Learner")) {
          showNotification("Invalid learner configuration", type = "error")
          return()
        }
        # Keep as-is, it's a valid list of learners
      } else if (!is.list(learners)) {
        # Not a list and not a Learner - invalid
        showNotification("No valid learners configured", type = "error")
        return()
      }

      n_learners <- length(learners)
      message(sprintf("[mod_benchmark] Starting benchmark with %d learners", n_learners))

      benchmark_running(TRUE)
      # Ensure running state is reset even on error
      on.exit(benchmark_running(FALSE), add = TRUE)

      # Use withProgress for basic progress indication
      withProgress(message = sprintf("Running nested CV with %d learners...", n_learners), value = 0, {
        tryCatch({
          incProgress(0.1, detail = "Setting up...")

          # Safely get worker count (may not be available on first render)
          n_workers <- if (input$parallel) {
            # Use input$workers if available, otherwise use recommended default
            workers_val <- input$workers
            if (is.null(workers_val) || is.na(workers_val) || workers_val < 1) {
              recommended_cores
            } else {
              min(workers_val, available_cores)
            }
          } else {
            1
          }

          # Run benchmark with all learners
          result <- pipeline$benchmark(
            learners = learners,
            outer_folds = input$outer_folds,
            inner_folds = input$inner_folds,
            stratify = input$stratify,
            seed = input$seed,
            parallel = input$parallel,
            threads = n_workers
          )

          incProgress(0.8, detail = "Processing results...")

          proj$set_benchmark_result(result)

          incProgress(0.1, detail = "Done!")

          showNotification(sprintf("Benchmark completed! Compared %d learners.", n_learners), type = "message")
        }, error = function(e) {
          message("[mod_benchmark] ERROR: ", e$message)
          showNotification(paste("Error:", e$message), type = "error", duration = 10)
        })
      })
    })

    # Performance table
    output$performance_table <- renderTable({
      proj <- project()
      if (is.null(proj)) return(NULL)

      result <- proj$get_benchmark_result()
      if (is.null(result)) {
        return(data.frame(Message = "Run benchmark to see results"))
      }

      if (!is.null(result$performance) && is.data.frame(result$performance)) {
        result$performance
      } else {
        data.frame(Message = "No performance data available")
      }
    })

    # Stability summary
    output$stability_summary <- renderPrint({
      proj <- project()
      if (is.null(proj)) return(invisible(NULL))

      result <- proj$get_benchmark_result()
      if (is.null(result)) {
        cat("Run benchmark to see stability analysis")
        return(invisible(NULL))
      }

      if (!is.null(result$stability)) {
        print(result$stability)
      } else {
        cat("Stability analysis not available")
      }
    })

    # Continue to results
    observeEvent(input$next_step, {
      proj <- project()
      if (is.null(proj) || !proj$can_proceed("results")) {
        showNotification("Run benchmark first", type = "error")
        return()
      }
      switch_tab("results")
    })
  })
}
