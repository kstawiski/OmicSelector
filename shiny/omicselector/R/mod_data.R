#' Data Import Module
#'
#' @description Data upload and variable role mapping with validation

# Module UI
mod_data_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # Step indicator
    div(
      class = "step-header mb-4",
      h3(bsicons::bs_icon("database"), "Step 1: Data Import"),
      p(class = "text-muted", "Upload your omics data and map column roles")
    ),

    layout_column_wrap(
      width = 1/2,
      fill = FALSE,

      # Upload Card
      card(
        card_header("Data Upload"),
        card_body(
          radioButtons(
            ns("data_source"),
            "Data Source",
            choices = c("Upload File" = "upload", "Built-in Dataset" = "builtin"),
            selected = "upload",
            inline = TRUE
          ),

          conditionalPanel(
            condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
            radioButtons(
              ns("upload_mode"),
              "Upload Mode",
              choices = c("Single Dataset" = "single", "Multi-omics (Multiple Files)" = "multi"),
              selected = "single",
              inline = TRUE
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'single'", ns("upload_mode")),
              fileInput(
                ns("data_file"),
                "Upload Data File",
                # Security: .rds/.rda removed - deserialization can execute arbitrary code
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
              ),
              p(class = "text-muted small", "Supported: CSV, TSV, TXT, Excel")
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'multi'", ns("upload_mode")),
              fileInput(
                ns("multi_files"),
                "Upload Multiple Files",
                multiple = TRUE,
                # Security: .rds/.rda removed - deserialization can execute arbitrary code
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
              ),
              p(class = "text-muted small", "Supported: CSV, TSV, TXT, Excel")
            )
          ),

          conditionalPanel(
            condition = sprintf("input['%s'] == 'builtin'", ns("data_source")),
            selectInput(
              ns("builtin_dataset"),
              "Built-in Dataset",
              choices = c(
                "TCGA miRNA subset" = "tcga",
                "Tutorial balanced" = "tutorial_balanced",
                "Tutorial mixed" = "tutorial_mixed"
              ),
              selected = "tcga"
            )
          ),

          actionButton(ns("load_data"), "Load Data",
                       class = "btn-success w-100 mt-3",
                       icon = icon("upload"))
        )
      ),

      # Data Preview Card
      card(
        card_header("Data Preview"),
        card_body(
          uiOutput(ns("data_info")),
          div(style = "overflow-x: auto; max-height: 300px;",
              tableOutput(ns("data_preview")))
        )
      )
    ),

    # Variable Role Mapping Card
    card(
      class = "mt-4",
      card_header(
        div(class = "d-flex align-items-center",
            bsicons::bs_icon("diagram-3"),
            span("Variable Role Mapping", class = "ms-2"),
            span(class = "ms-auto badge bg-info", "Required"))
      ),
      card_body(
        p(class = "text-muted mb-3",
          "Assign roles to your data columns. This ensures proper leakage prevention."),

        layout_column_wrap(
          width = 1/3,
          fill = FALSE,

          div(
            class = "role-box role-required",
            h6(bsicons::bs_icon("bullseye"), " Target Variable"),
            uiOutput(ns("target_ui")),
            uiOutput(ns("positive_class_ui"))
          ),

          div(
            class = "role-box role-recommended",
            h6(bsicons::bs_icon("person-badge"), " Patient ID"),
            uiOutput(ns("patient_id_ui")),
            p(class = "text-muted small",
              "Recommended for grouped CV (prevents leakage)")
          ),

          div(
            class = "role-box role-optional",
            h6(bsicons::bs_icon("collection"), " Batch Column"),
            uiOutput(ns("batch_ui")),
            p(class = "text-muted small",
              "Optional for batch effect correction")
          )
        ),

        hr(),

        layout_column_wrap(
          width = 1/2,
          fill = FALSE,

          div(
            textInput(ns("feature_prefix"), "Feature Prefix Filter",
                      value = "", placeholder = "e.g., hsa"),
            p(class = "text-muted small",
              "Filter features by prefix (leave empty for all numeric columns)")
          ),

          div(
            uiOutput(ns("row_id_ui")),
            p(class = "text-muted small",
              "Optional: use a column as row identifiers")
          )
        )
      ),
      card_footer(
        div(class = "d-flex justify-content-between",
            actionButton(ns("apply_mapping"), "Apply Mapping",
                         class = "btn-primary",
                         icon = icon("check")),
            actionButton(ns("next_step"), "Continue to QC",
                         class = "btn-outline-success",
                         icon = icon("arrow-right")))
      )
    ),

    # File Manifest (for multi-omics)
    uiOutput(ns("file_manifest_ui"))
  )
}

# Module Server
mod_data_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive trigger to force UI updates when R6 object state changes
    data_trigger <- reactiveVal(0)

    # Helper to read table files
    read_table_file <- function(path, clean_names = TRUE) {
      ext <- tolower(tools::file_ext(path))

      # Security: Reject RDS/RDA files - deserialization can execute arbitrary code
      if (ext %in% c("rds", "rda", "rdata")) {
        stop("RDS/RDA files are not supported for security reasons. Please use CSV, TSV, or Excel format.", call. = FALSE)
      }

      if (ext %in% c("csv", "tsv", "txt")) {
        df <- data.table::fread(path, data.table = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        if (requireNamespace("readxl", quietly = TRUE)) {
          df <- as.data.frame(readxl::read_excel(path))
        } else {
          stop("Package 'readxl' required for Excel files", call. = FALSE)
        }
      } else {
        stop("Unsupported file type. Use CSV, TSV, TXT, or Excel.", call. = FALSE)
      }

      df <- as.data.frame(df)
      if (isTRUE(clean_names)) {
        names(df) <- make.names(names(df), unique = TRUE)
      }
      df
    }

    # Load built-in dataset
    load_builtin <- function(name) {
      if (name == "tcga") {
        data("original_TCGA_data", package = "OmicSelector", envir = environment())
        df <- as.data.frame(get("original_TCGA_data", envir = environment()))
        # Subset for demo - limit features for faster exploration
        feature_cols <- grep("^hsa", names(df), value = TRUE)
        total_features <- length(feature_cols)
        demo_limit <- 200
        feature_cols <- head(feature_cols, demo_limit)
        meta_cols <- c("patient", "sample_type")
        meta_cols <- meta_cols[meta_cols %in% names(df)]
        df <- df[, c(meta_cols, feature_cols), drop = FALSE]
        # Notify user about demo limitation
        showNotification(
          sprintf("TCGA demo: Using %d of %d features for faster exploration. Upload full data for complete analysis.",
                  demo_limit, total_features),
          type = "message",
          duration = 8
        )
        return(df)
      } else if (name == "tutorial_balanced") {
        data("miRNAselector_tutorial_balanced_dataset", package = "OmicSelector",
             envir = environment())
        return(as.data.frame(get("miRNAselector_tutorial_balanced_dataset",
                                 envir = environment())))
      } else if (name == "tutorial_mixed") {
        data("miRNAselector_tutorial_balanced_mixed", package = "OmicSelector",
             envir = environment())
        return(as.data.frame(get("miRNAselector_tutorial_balanced_mixed",
                                 envir = environment())))
      }
      stop("Unknown built-in dataset", call. = FALSE)
    }

    # Load data
    observeEvent(input$load_data, {
      message("[mod_data] Load button clicked")
      proj <- project()
      if (is.null(proj)) {
        message("[mod_data] ERROR: No project loaded")
        showNotification("Please create a project first", type = "error")
        return()
      }
      message("[mod_data] Project exists: ", proj$id)

      withProgress(message = "Loading data...", value = 0.2, {
        tryCatch({
          data_loaded <- NULL
          is_multi <- FALSE

          message("[mod_data] Data source: ", input$data_source)
          if (input$data_source == "builtin") {
            message("[mod_data] Loading built-in: ", input$builtin_dataset)
            data_loaded <- load_builtin(input$builtin_dataset)
          } else if (input$upload_mode == "single") {
            req(input$data_file)
            data_loaded <- read_table_file(input$data_file$datapath)
          } else {
            req(input$multi_files)
            files <- input$multi_files
            data_loaded <- list()
            for (i in seq_len(nrow(files))) {
              mod_name <- tools::file_path_sans_ext(files$name[i])
              mod_name <- make.names(mod_name)
              df <- read_table_file(files$datapath[i])
              data_loaded[[mod_name]] <- df
            }
            is_multi <- TRUE

            # P0 Fix: Validate sample alignment across modalities
            if (length(data_loaded) > 1) {
              nrows <- sapply(data_loaded, nrow)
              if (length(unique(nrows)) != 1) {
                # Different row counts - warn user
                row_info <- paste(names(nrows), nrows, sep = "=", collapse = ", ")
                showNotification(
                  sprintf("Warning: Multi-omics files have different row counts (%s). Please ensure samples are aligned across modalities, or use a common Patient ID column for alignment.", row_info),
                  type = "warning",
                  duration = 10
                )
              }
            }
          }

          incProgress(0.5)
          # Log dimensions correctly for both single and multi-omics
          if (is.data.frame(data_loaded)) {
            message("[mod_data] Data loaded, dims: ", nrow(data_loaded), " x ", ncol(data_loaded))
          } else if (is.list(data_loaded)) {
            dims_str <- paste(sapply(names(data_loaded), function(m) {
              sprintf("%s: %dx%d", m, nrow(data_loaded[[m]]), ncol(data_loaded[[m]]))
            }), collapse = ", ")
            message("[mod_data] Multi-omics loaded: ", dims_str)
          }

          project()$set_data(data_loaded, is_multi)
          message("[mod_data] set_data completed")

          incProgress(0.3)

          # Trigger UI refresh (R6 state changes don't trigger reactivity)
          data_trigger(data_trigger() + 1)

          showNotification("Data loaded successfully!", type = "message")
          message("[mod_data] SUCCESS")
        }, error = function(e) {
          message("[mod_data] ERROR: ", e$message)
          showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })

    # Get column names reactively
    columns <- reactive({
      # Depend on trigger to refresh when data changes
      data_trigger()
      proj <- project()
      if (is.null(proj)) return(character(0))
      data <- proj$get_data()
      if (is.null(data)) return(character(0))
      if (is.data.frame(data)) {
        names(data)
      } else if (is.list(data)) {
        unique(unlist(lapply(data, names)))
      } else {
        character(0)
      }
    })

    # Data info display
    output$data_info <- renderUI({
      # Depend on trigger to refresh when data changes
      data_trigger()
      proj <- project()
      if (is.null(proj)) return(NULL)
      data <- proj$get_data()
      if (is.null(data)) {
        return(div(class = "text-muted", "No data loaded"))
      }

      if (is.data.frame(data)) {
        tags$div(
          class = "alert alert-info",
          sprintf("Loaded: %d samples x %d columns", nrow(data), ncol(data))
        )
      } else {
        tags$div(
          class = "alert alert-info",
          sprintf("Multi-omics: %d modalities", length(data)),
          tags$ul(
            lapply(names(data), function(m) {
              tags$li(sprintf("%s: %d x %d", m, nrow(data[[m]]), ncol(data[[m]])))
            })
          )
        )
      }
    })

    # Data preview
    output$data_preview <- renderTable({
      # Depend on trigger to refresh when data changes
      data_trigger()
      proj <- project()
      if (is.null(proj)) return(NULL)
      data <- proj$get_data()
      if (is.null(data)) return(NULL)

      if (is.data.frame(data)) {
        head(data, 6)
      } else if (is.list(data)) {
        # For multi-omics, show combined preview with modality prefix
        # Take first few columns from each modality
        preview_dfs <- lapply(names(data), function(mod_name) {
          mod_data <- data[[mod_name]]
          # Take first 3 columns from each modality
          n_cols <- min(3, ncol(mod_data))
          sub_df <- mod_data[, 1:n_cols, drop = FALSE]
          # Add modality prefix to column names
          names(sub_df) <- paste0(mod_name, "::", names(sub_df))
          sub_df
        })
        combined <- do.call(cbind, preview_dfs)
        head(combined, 6)
      }
    }, rownames = TRUE)

    # Target column UI
    output$target_ui <- renderUI({
      cols <- columns()
      if (length(cols) == 0) {
        return(p(class = "text-muted", "Load data first"))
      }
      selectInput(ns("target_col"), NULL, choices = c("Select..." = "", cols))
    })

    # Positive class UI
    output$positive_class_ui <- renderUI({
      proj <- project()
      req(input$target_col)
      if (is.null(proj) || !nzchar(input$target_col)) return(NULL)

      data <- proj$get_data()
      if (is.null(data)) return(NULL)

      target_vals <- NULL
      if (is.data.frame(data) && input$target_col %in% names(data)) {
        target_vals <- data[[input$target_col]]
      } else if (is.list(data)) {
        for (m in names(data)) {
          if (input$target_col %in% names(data[[m]])) {
            target_vals <- data[[m]][[input$target_col]]
            break
          }
        }
      }

      if (is.null(target_vals)) return(NULL)

      if (is.numeric(target_vals)) {
        return(p(class = "text-muted small", "Numeric target (regression)"))
      }

      classes <- unique(target_vals)
      if (length(classes) > 10) {
        return(p(class = "text-warning small", "Many classes detected"))
      }

      selectInput(ns("positive_class"), "Positive Class", choices = classes)
    })

    # Patient ID UI
    output$patient_id_ui <- renderUI({
      cols <- columns()
      if (length(cols) == 0) return(NULL)
      selectInput(ns("patient_id_col"), NULL,
                  choices = c("None" = "", cols))
    })

    # Batch UI
    output$batch_ui <- renderUI({
      cols <- columns()
      if (length(cols) == 0) return(NULL)
      selectInput(ns("batch_col"), NULL,
                  choices = c("None" = "", cols))
    })

    # Row ID UI
    output$row_id_ui <- renderUI({
      cols <- columns()
      if (length(cols) == 0) return(NULL)
      selectInput(ns("row_id_col"), "Row ID Column",
                  choices = c("None" = "", cols))
    })

    # File manifest for multi-omics
    output$file_manifest_ui <- renderUI({
      proj <- project()
      if (is.null(proj)) return(NULL)
      summary <- proj$summary()
      if (!summary$is_multi_omics) return(NULL)

      data <- proj$get_data()
      if (!is.list(data)) return(NULL)

      card(
        class = "mt-4",
        card_header("Multi-omics File Manifest"),
        card_body(
          tableOutput(ns("file_manifest_table"))
        )
      )
    })

    output$file_manifest_table <- renderTable({
      proj <- project()
      if (is.null(proj)) return(NULL)
      data <- proj$get_data()
      if (!is.list(data)) return(NULL)

      data.frame(
        Modality = names(data),
        Samples = sapply(data, nrow),
        Features = sapply(data, ncol),
        stringsAsFactors = FALSE
      )
    })

    # Apply mapping
    observeEvent(input$apply_mapping, {
      req(project())

      target <- input$target_col
      if (is.null(target) || !nzchar(target)) {
        showNotification("Please select a target column", type = "error")
        return()
      }

      # Check if target is numeric - if so, don't pass positive class
      data <- project()$get_data()
      positive_class <- input$positive_class

      if (!is.null(data)) {
        target_vals <- NULL
        if (is.data.frame(data) && target %in% names(data)) {
          target_vals <- data[[target]]
        } else if (is.list(data)) {
          for (m in names(data)) {
            if (target %in% names(data[[m]])) {
              target_vals <- data[[m]][[target]]
              break
            }
          }
        }

        # For numeric targets, don't set positive class (regression task)
        if (!is.null(target_vals) && is.numeric(target_vals)) {
          positive_class <- NULL
        }
      }

      project()$set_mapping(
        target = target,
        positive = positive_class,
        patient_id = input$patient_id_col,
        batch = input$batch_col,
        feature_prefix = input$feature_prefix
      )

      showNotification("Mapping applied!", type = "message")
    })

    # Continue to QC
    observeEvent(input$next_step, {
      proj <- project()
      if (is.null(proj)) {
        showNotification("No project loaded", type = "error")
        return()
      }

      if (is.null(proj$get_data())) {
        showNotification("Please load data first", type = "error")
        return()
      }

      mapping <- proj$get_mapping()
      if (is.null(mapping$target) || !nzchar(mapping$target)) {
        showNotification("Please select a target column", type = "error")
        return()
      }

      switch_tab("qc")
    })
  })
}
