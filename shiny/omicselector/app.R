# OmicSelector WebUI - Main Application
# Biomarker Discovery with Zero-Leakage Machine Learning

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(bslib)
  library(bsicons)
  library(shinyWidgets)
  library(OmicSelector)
  library(data.table)
  library(plotly)
  library(R6)
})

options(shiny.maxRequestSize = 100 * 1024^2)

# Source modules
source("R/class_OmicProject.R")
source("R/mod_home.R")
source("R/mod_data.R")
source("R/mod_qc.R")
source("R/mod_pipeline.R")
source("R/mod_benchmark.R")
source("R/mod_results.R")
source("R/mod_export.R")

# Analysis directory management
analysis_root <- function() {
  root <- Sys.getenv("OMICSELECTOR_ANALYSES_DIR", unset = "/OmicSelector/analyses")
  if (!nzchar(root)) {
    root <- "/OmicSelector/analyses"
  }
  if (!dir.exists(root)) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(root)) {
    root <- file.path(getwd(), "analyses")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

# Theme definition
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#2563eb",
  secondary = "#475569",
  success = "#10b981",
  warning = "#f59e0b",
  danger = "#ef4444",
  "font-family-base" = "Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"
)

# Workflow steps definition
workflow_steps <- list(
  list(id = "home", num = "0", title = "Home", subtitle = "Project setup", icon = "house"),
  list(id = "data", num = "1", title = "Data Import", subtitle = "Load & map variables", icon = "database"),
  list(id = "qc", num = "2", title = "QC & EDA", subtitle = "Quality control", icon = "clipboard-check"),
  list(id = "pipeline", num = "3", title = "Pipeline", subtitle = "Configure ML pipeline", icon = "gear"),
  list(id = "benchmark", num = "4", title = "Benchmark", subtitle = "Nested cross-validation", icon = "speedometer2"),
  list(id = "results", num = "5", title = "Results", subtitle = "Select signature", icon = "trophy"),
  list(id = "export", num = "6", title = "Export", subtitle = "Deploy & report", icon = "box-arrow-up-right")
)

# UI
ui <- page_sidebar(
  title = div(
    class = "d-flex align-items-center",
    bsicons::bs_icon("activity", size = "1.5em"),
    span("OmicSelector", class = "ms-2 fw-bold"),
    span(class = "badge bg-primary ms-2", "v2")
  ),
  theme = app_theme,
  fillable = TRUE,

  # Sidebar with workflow stepper
  sidebar = sidebar(
    id = "main_sidebar",
    width = 280,
    bg = "#f8fafc",

    # Project info section
    uiOutput("project_info_sidebar"),

    hr(class = "my-3"),

    # Workflow steps
    div(
      class = "workflow-steps",

      lapply(workflow_steps, function(step) {
        actionLink(
          inputId = paste0("step_", step$id),
          label = NULL,
          class = "workflow-step",
          div(
            class = "step-number",
            if (step$num == "0") bsicons::bs_icon(step$icon, size = "0.9em") else step$num
          ),
          div(
            class = "step-content",
            div(class = "step-title", step$title),
            div(class = "step-subtitle", step$subtitle)
          )
        )
      })
    ),

    hr(class = "my-3"),

    # Dark mode toggle
    div(
      class = "d-flex justify-content-between align-items-center px-2",
      span("Dark Mode", class = "small"),
      input_dark_mode(id = "dark_mode", mode = "light")
    )
  ),

  # Main content area with conditional panels
  div(
    id = "main_content",
    class = "p-3",

    # Home module
    conditionalPanel(
      condition = "output.current_tab == 'home'",
      mod_home_ui("home")
    ),

    # Data module
    conditionalPanel(
      condition = "output.current_tab == 'data'",
      mod_data_ui("data")
    ),

    # QC module
    conditionalPanel(
      condition = "output.current_tab == 'qc'",
      mod_qc_ui("qc")
    ),

    # Pipeline module
    conditionalPanel(
      condition = "output.current_tab == 'pipeline'",
      mod_pipeline_ui("pipeline")
    ),

    # Benchmark module
    conditionalPanel(
      condition = "output.current_tab == 'benchmark'",
      mod_benchmark_ui("benchmark")
    ),

    # Results module
    conditionalPanel(
      condition = "output.current_tab == 'results'",
      mod_results_ui("results")
    ),

    # Export module
    conditionalPanel(
      condition = "output.current_tab == 'export'",
      mod_export_ui("export")
    )
  ),

  # Include custom CSS and JS
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('updateStepClasses', function(data) {
        $('.workflow-step').removeClass('active completed');
        var stepOrder = data.steps;
        var activeIndex = stepOrder.indexOf(data.active);

        stepOrder.forEach(function(stepId, index) {
          var el = $('#step_' + stepId);
          if (stepId === data.active) {
            el.addClass('active');
          } else if (index < activeIndex) {
            el.addClass('completed');
          }
        });
      });
    "))
  ),

  # Enable shinyjs
  useShinyjs()
)

# Server
server <- function(input, output, session) {

  # Current project (shared state)
  project <- reactiveVal(NULL)

  # Current active tab
  current_tab <- reactiveVal("home")

  # Output for conditional panels
  output$current_tab <- renderText({
    current_tab()
  })
  outputOptions(output, "current_tab", suspendWhenHidden = FALSE)

  # Switch tab function (passed to modules)
  switch_tab <- function(tab_id) {
    current_tab(tab_id)

    # Update step classes via JavaScript
    shinyjs::runjs(sprintf("
      $('.workflow-step').removeClass('active');
      $('#step_%s').addClass('active');
    ", tab_id))
  }

  # Handle step clicks
  lapply(workflow_steps, function(step) {
    observeEvent(input[[paste0("step_", step$id)]], {
      proj <- project()

      # Check if navigation is allowed
      if (step$id == "home") {
        current_tab("home")
      } else if (!is.null(proj) && proj$can_proceed(step$id)) {
        current_tab(step$id)
      } else if (is.null(proj) && step$id != "home") {
        showNotification("Please create or load a project first", type = "warning")
        current_tab("home")
      } else {
        showNotification("Complete previous steps first", type = "warning")
      }
    }, ignoreInit = TRUE)
  })

  # Update sidebar step classes based on current tab
  observe({
    tab <- current_tab()

    # Use shinyjs to update classes
    session$sendCustomMessage("updateStepClasses", list(
      active = tab,
      steps = sapply(workflow_steps, function(s) s$id)
    ))
  })

  # Project info in sidebar
  output$project_info_sidebar <- renderUI({
    proj <- project()

    if (is.null(proj)) {
      div(
        class = "project-info",
        div(class = "project-name text-muted", "No project loaded"),
        div(class = "project-meta", "Create or load a project to begin")
      )
    } else {
      summary <- proj$summary()
      div(
        class = "project-info",
        div(class = "project-name", summary$name),
        div(class = "project-meta",
            sprintf("%s samples | %s features",
                    if (!is.na(summary$n_samples)) summary$n_samples else "0",
                    if (!is.na(summary$n_features)) summary$n_features else "0")),
        if (!is.na(summary$qc_overall) && summary$qc_overall != "unknown") {
          # Map green/yellow/red to CSS classes pass/warn/fail
          css_class <- switch(summary$qc_overall,
            "green" = "pass",
            "yellow" = "warn",
            "red" = "fail",
            "unknown"
          )
          icon_name <- switch(summary$qc_overall,
            "green" = "check-circle-fill",
            "yellow" = "exclamation-triangle-fill",
            "red" = "x-circle-fill",
            "question-circle-fill"
          )
          label_text <- switch(summary$qc_overall,
            "green" = "QC Passed",
            "yellow" = "QC Warnings",
            "red" = "QC Failed",
            "QC Unknown"
          )
          div(
            class = paste0("qc-summary-badge mt-2 ", css_class),
            bsicons::bs_icon(icon_name),
            span(class = "ms-1", label_text)
          )
        }
      )
    }
  })

  # Initialize modules
  mod_home_server("home", project, switch_tab)
  mod_data_server("data", project, switch_tab)
  mod_qc_server("qc", project, switch_tab)
  mod_pipeline_server("pipeline", project, switch_tab)
  mod_benchmark_server("benchmark", project, switch_tab)
  mod_results_server("results", project, switch_tab)
  mod_export_server("export", project, switch_tab)
}

# Run application
shinyApp(ui = ui, server = server)
