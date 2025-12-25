#' Home/Project Module
#'
#' @description Landing page for project selection and management

# Module UI
mod_home_ui <- function(id) {
  ns <- NS(id)

  tagList(
    div(
      class = "home-header text-center mb-4",
      tags$img(src = "https://biostat.umed.pl/OmicSelector/logo.png",
               height = "80px", class = "mb-3"),
      h1("OmicSelector", class = "display-4"),
      p(class = "lead text-muted",
        "Rigorous Biomarker Discovery with Zero-Leakage Machine Learning")
    ),

    layout_column_wrap(
      width = 1/2,
      fill = FALSE,

      # New Project Card
      card(
        card_header(
          class = "bg-primary text-white",
          div(class = "d-flex align-items-center",
              bsicons::bs_icon("plus-circle", size = "1.5em"),
              span("New Project", class = "ms-2 fw-bold"))
        ),
        card_body(
          textInput(ns("new_project_name"), "Project Name",
                    placeholder = "e.g., Glioblastoma_Signature_v1"),
          p(class = "text-muted small",
            "Start a new biomarker discovery analysis."),
          actionButton(ns("create_project"), "Create Project",
                       class = "btn-primary w-100",
                       icon = icon("rocket"))
        )
      ),

      # Load Project Card
      card(
        card_header(
          class = "bg-secondary text-white",
          div(class = "d-flex align-items-center",
              bsicons::bs_icon("folder2-open", size = "1.5em"),
              span("Load Project", class = "ms-2 fw-bold"))
        ),
        card_body(
          fileInput(ns("load_project_file"), "Project File (.rds)",
                    accept = ".rds"),
          p(class = "text-muted small",
            "Continue from a saved project state."),
          actionButton(ns("load_project"), "Load Project",
                       class = "btn-secondary w-100",
                       icon = icon("upload"))
        )
      )
    ),

    # Recent Projects
    card(
      class = "mt-4",
      card_header("Recent Projects"),
      card_body(
        uiOutput(ns("recent_projects"))
      )
    ),

    # Quick Start Guide
    card(
      class = "mt-4",
      card_header(
        div(class = "d-flex align-items-center",
            bsicons::bs_icon("lightbulb", size = "1.2em"),
            span("Quick Start Guide", class = "ms-2"))
      ),
      card_body(
        div(
          class = "workflow-steps",
          div(class = "workflow-step",
              span(class = "step-number", "1"),
              span(class = "step-text", "Upload your omics data (CSV, TSV, Excel, or RDS)")),
          div(class = "workflow-step",
              span(class = "step-number", "2"),
              span(class = "step-text", "Map target variable and validate data quality")),
          div(class = "workflow-step",
              span(class = "step-number", "3"),
              span(class = "step-text", "Configure feature selection and ML pipeline")),
          div(class = "workflow-step",
              span(class = "step-number", "4"),
              span(class = "step-text", "Run nested cross-validation benchmark")),
          div(class = "workflow-step",
              span(class = "step-number", "5"),
              span(class = "step-text", "Select optimal signature using stability metrics")),
          div(class = "workflow-step",
              span(class = "step-number", "6"),
              span(class = "step-text", "Export model for deployment"))
        )
      )
    )
  )
}

# Module Server
mod_home_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Helper function to determine appropriate tab based on project state
    get_route_for_project <- function(proj) {
      summary <- proj$summary()
      if (summary$has_final_fit) {
        "export"
      } else if (summary$has_benchmark) {
        "results"
      } else if (summary$has_pipeline) {
        "benchmark"
      } else if (!is.na(summary$n_samples)) {
        "qc"
      } else {
        "data"
      }
    }

    # Create new project
    observeEvent(input$create_project, {
      name <- input$new_project_name
      if (is.null(name) || !nzchar(trimws(name))) {
        name <- paste("Project", format(Sys.time(), "%Y%m%d-%H%M"))
      }

      new_proj <- OmicProject$new(name = name)
      project(new_proj)

      showNotification(
        sprintf("Project '%s' created!", name),
        type = "message",
        duration = 3
      )

      # Switch to Data Import tab
      switch_tab("data")
    })

    # Load existing project
    observeEvent(input$load_project, {
      req(input$load_project_file)

      tryCatch({
        proj <- OmicProject$new()
        proj$load(input$load_project_file$datapath)
        project(proj)

        showNotification(
          sprintf("Project '%s' loaded!", proj$name),
          type = "message",
          duration = 3
        )

        # Switch to appropriate tab based on project state
        switch_tab(get_route_for_project(proj))
      }, error = function(e) {
        showNotification(
          paste("Error loading project:", e$message),
          type = "error",
          duration = 5
        )
      })
    })

    # Recent projects list
    output$recent_projects <- renderUI({
      root <- Sys.getenv("OMICSELECTOR_ANALYSES_DIR", "/OmicSelector/analyses")
      if (!dir.exists(root)) {
        root <- file.path(getwd(), "analyses")
      }

      if (!dir.exists(root)) {
        return(p(class = "text-muted", "No recent projects found."))
      }

      # Find project state files (limit depth to avoid scanning very deep directories)
      # Only scan immediate subdirectories (project folders)
      project_dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
      project_files <- character()

      # Limit search to direct children and one level deep to avoid performance issues
      for (dir in c(root, project_dirs)) {
        files_in_dir <- list.files(dir, pattern = "project_state\\.rds$", full.names = TRUE)
        project_files <- c(project_files, files_in_dir)
        # Stop early if we have enough files
        if (length(project_files) >= 20) break
      }

      if (length(project_files) == 0) {
        return(p(class = "text-muted", "No recent projects found."))
      }

      # Get info for up to 5 most recent
      project_files <- head(project_files[order(file.mtime(project_files),
                                                 decreasing = TRUE)], 5)

      project_items <- lapply(seq_along(project_files), function(i) {
        pf <- project_files[i]
        info <- tryCatch({
          state <- readRDS(pf)
          list(
            name = state$name %||% basename(dirname(pf)),
            modified = state$modified_at,
            path = pf
          )
        }, error = function(e) {
          list(name = basename(dirname(pf)), modified = file.mtime(pf), path = pf)
        })

        tags$button(
          class = "list-group-item list-group-item-action d-flex justify-content-between align-items-center",
          onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})",
                            ns("quick_load"), info$path),
          span(info$name),
          span(class = "badge bg-secondary",
               format(info$modified, "%Y-%m-%d %H:%M"))
        )
      })

      div(class = "list-group", project_items)
    })

    # Quick load from recent projects
    observeEvent(input$quick_load, {
      req(input$quick_load)

      tryCatch({
        proj <- OmicProject$new()
        proj$load(input$quick_load)
        project(proj)

        showNotification(
          sprintf("Project '%s' loaded!", proj$name),
          type = "message",
          duration = 3
        )

        # Switch to appropriate tab based on project state
        switch_tab(get_route_for_project(proj))
      }, error = function(e) {
        showNotification(
          paste("Error loading project:", e$message),
          type = "error",
          duration = 5
        )
      })
    })
  })
}
