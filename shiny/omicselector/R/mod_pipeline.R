#' Pipeline Designer Module
#'
#' @description Configure ML pipeline with guided and advanced modes

# Module UI
mod_pipeline_ui <- function(id) {
  ns <- NS(id)

  tagList(
    div(
      class = "step-header mb-4",
      h3(bsicons::bs_icon("gear"), "Step 3: Pipeline Configuration"),
      p(class = "text-muted", "Design your biomarker discovery pipeline")
    ),

    # Pipeline Designer Cards
    layout_column_wrap(
      width = 1/3,
      fill = FALSE,

      # Preprocessing Card
      card(
        card_header(
          class = "bg-light",
          div(class = "d-flex align-items-center",
              bsicons::bs_icon("funnel"),
              span("1. Preprocessing", class = "ms-2 fw-bold"))
        ),
        card_body(
          selectInput(ns("impute_method"), "Imputation",
                      choices = c("Median" = "median", "Mean" = "mean", "Sample" = "sample"),
                      selected = "median"),

          checkboxInput(ns("scale_features"), "Scale Features", value = TRUE),

          selectInput(ns("oversample"), "Class Balancing",
                      choices = c("None" = "none", "SMOTE" = "smote", "ROSE" = "rose"),
                      selected = "none"),

          accordion(
            id = ns("preproc_advanced"),
            open = FALSE,
            accordion_panel(
              "Advanced Options",
              checkboxInput(ns("batch_correct"), "Batch Correction (FrozenComBat)", value = FALSE),
              checkboxInput(ns("screening"), "Variance Screening", value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s']", ns("screening")),
                numericInput(ns("screening_frac"), "Screening Fraction", value = 0.2, min = 0.01, max = 1)
              )
            )
          )
        )
      ),

      # Feature Selection Card
      card(
        card_header(
          class = "bg-light",
          div(class = "d-flex align-items-center",
              bsicons::bs_icon("funnel-fill"),
              span("2. Feature Selection", class = "ms-2 fw-bold"))
        ),
        card_body(
          selectizeInput(ns("filter_methods"), "Selection Methods (multi-select)",
                      choices = c(
                        "Univariate (ANOVA)" = "anova",
                        "Kruskal-Wallis" = "kruskal",
                        "Information Gain" = "information_gain",
                        "mRMR" = "mrmr",
                        "Correlation" = "correlation",
                        "AUC" = "auc",
                        "Model Importance" = "importance",
                        "Permutation" = "permutation"
                      ),
                      selected = c("anova", "information_gain"),
                      multiple = TRUE,
                      options = list(plugins = list('remove_button'))),

          numericInput(ns("n_features"), "Number of Features", value = 20, min = 1, max = 500),

          accordion(
            id = ns("fs_advanced"),
            open = FALSE,
            accordion_panel(
              "Auto-Tuning",
              checkboxInput(ns("tune_nfeat"), "Auto-tune Feature Count", value = FALSE),
              textInput(ns("tune_values"), "Candidate Values", value = "5,10,20,50")
            )
          )
        )
      ),

      # Model Selection Card
      card(
        card_header(
          class = "bg-light",
          div(class = "d-flex align-items-center",
              bsicons::bs_icon("cpu"),
              span("3. Model", class = "ms-2 fw-bold"))
        ),
        card_body(
          selectizeInput(ns("models"), "Algorithms (multi-select)",
                      choices = c(
                        "Random Forest" = "ranger",
                        "Elastic Net" = "glmnet",
                        "Logistic Regression" = "log_reg",
                        "SVM" = "svm",
                        "XGBoost" = "xgboost",
                        "k-NN" = "kknn",
                        "Neural Network (MLP)" = "mlp",
                        "TabTransformer" = "tabtransformer"
                      ),
                      selected = c("ranger", "glmnet"),
                      multiple = TRUE,
                      options = list(plugins = list('remove_button'))),

          accordion(
            id = ns("model_advanced"),
            open = FALSE,
            accordion_panel(
              "Deep Learning Settings",
              # Show DL settings when MLP or TabTransformer is selected (multi-select check)
              conditionalPanel(
                condition = sprintf("input['%s'] && (input['%s'].includes('mlp') || input['%s'].includes('tabtransformer'))",
                                    ns("models"), ns("models"), ns("models")),
                textInput(ns("dl_hidden"), "Hidden Layers", value = "64,32"),
                numericInput(ns("dl_epochs"), "Epochs", value = 100, min = 10),
                numericInput(ns("dl_dropout"), "Dropout", value = 0.2, min = 0, max = 0.9, step = 0.05),
                selectInput(ns("dl_device"), "Device", choices = c("CPU" = "cpu", "GPU" = "cuda"))
              ),
              # Show message when no DL models selected
              conditionalPanel(
                condition = sprintf("!input['%s'] || (!input['%s'].includes('mlp') && !input['%s'].includes('tabtransformer'))",
                                    ns("models"), ns("models"), ns("models")),
                p(class = "text-muted small", "Select MLP or TabTransformer to configure deep learning settings")
              )
            ),
            accordion_panel(
              "Autoencoder",
              checkboxInput(ns("use_autoencoder"), "Enable Autoencoder", value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s']", ns("use_autoencoder")),
                numericInput(ns("ae_latent"), "Latent Dimension", value = 32, min = 2),
                numericInput(ns("ae_epochs"), "AE Epochs", value = 50, min = 10)
              )
            )
          )
        )
      )
    ),

    # Pipeline Summary Card
    card(
      class = "mt-4",
      card_header(
        div(class = "d-flex align-items-center justify-content-between",
            div(
              bsicons::bs_icon("diagram-3"),
              span("Pipeline Summary", class = "ms-2")
            ),
            actionButton(ns("view_code"), "View R Code",
                         class = "btn-sm btn-outline-secondary",
                         icon = icon("code")))
      ),
      card_body(
        verbatimTextOutput(ns("pipeline_summary")),

        hr(),

        layout_column_wrap(
          width = 1/2,
          actionButton(ns("create_pipeline"), "Create Pipeline",
                       class = "btn-primary btn-lg w-100",
                       icon = icon("wrench")),
          actionButton(ns("next_step"), "Continue to Benchmark",
                       class = "btn-success btn-lg w-100",
                       icon = icon("arrow-right"))
        )
      )
    ),

    # Code Modal
    uiOutput(ns("code_modal"))
  )
}

# Module Server
mod_pipeline_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Pipeline summary
    output$pipeline_summary <- renderText({
      proj <- project()
      if (is.null(proj)) return("No project loaded")

      mapping <- proj$get_mapping()
      if (is.null(mapping$target)) return("Target variable not set")

      n_filters <- length(input$filter_methods)
      n_models <- length(input$models)
      n_combos <- n_filters * n_models

      # Get data dimensions for cost estimate
      data <- proj$get_data()
      n_samples <- if (!is.null(data)) {
        if (is.data.frame(data)) nrow(data) else nrow(data[[1]])
      } else 0

      lines <- c(
        sprintf("Target: %s", mapping$target),
        sprintf("Preprocessing: %s imputation, %s",
                input$impute_method,
                if (input$scale_features) "scaled" else "unscaled"),
        sprintf("Class Balancing: %s", input$oversample),
        sprintf("Feature Selection: %s (top %d)", paste(input$filter_methods, collapse = ", "), input$n_features),
        sprintf("Models: %s", paste(input$models, collapse = ", ")),
        sprintf("Total Learner Combinations: %d", n_combos)
      )

      if (input$batch_correct) lines <- c(lines, "Batch correction: FrozenComBat")
      if (input$screening) lines <- c(lines, sprintf("Screening: %.0f%%", input$screening_frac * 100))
      if (input$use_autoencoder) lines <- c(lines, sprintf("Autoencoder: %d latent dims", input$ae_latent))

      # Add computational cost warning
      lines <- c(lines, "")
      if (n_combos > 20) {
        lines <- c(lines, sprintf("WARNING: %d combinations may take significant time!", n_combos))
      }
      # Rough estimate: ~30s per learner per 100 samples for basic models
      has_dl <- any(c("mlp", "tabtransformer") %in% input$models)
      has_xgb <- "xgboost" %in% input$models
      base_time <- if (has_dl) 120 else if (has_xgb) 45 else 30
      est_seconds <- n_combos * (n_samples / 100) * base_time / 5  # Divide by typical CV folds
      if (est_seconds > 60) {
        est_minutes <- round(est_seconds / 60)
        if (est_minutes > 60) {
          lines <- c(lines, sprintf("Estimated time: ~%d hours (single-threaded)", round(est_minutes / 60)))
        } else {
          lines <- c(lines, sprintf("Estimated time: ~%d minutes (single-threaded)", est_minutes))
        }
      }

      paste(lines, collapse = "\n")
    })

    # Create pipeline
    observeEvent(input$create_pipeline, {
      proj <- project()
      req(proj)

      data <- proj$get_data()
      mapping <- proj$get_mapping()

      if (is.null(data) || is.null(mapping$target)) {
        showNotification("Data and target required", type = "error")
        return()
      }

      filter_methods <- input$filter_methods
      models <- input$models

      if (length(filter_methods) == 0 || length(models) == 0) {
        showNotification("Select at least one feature selection method and one model", type = "error")
        return()
      }

      n_combos <- length(filter_methods) * length(models)

      withProgress(message = sprintf("Creating %d learner combinations...", n_combos), value = 0, {
        tryCatch({
          pipeline <- OmicPipeline$new(
            data = data,
            target = mapping$target,
            positive = mapping$positive,
            patient_id = mapping$patient_id,
            batch = mapping$batch
          )

          incProgress(0.2, detail = "Pipeline initialized")

          # Build learners for all combinations
          oversample <- if (input$oversample == "none") NULL else input$oversample
          ae_params <- NULL
          if (input$use_autoencoder) {
            ae_params <- list(
              latent_dim = input$ae_latent,
              epochs = input$ae_epochs
            )
          }

          learners <- list()
          combo_idx <- 0
            for (flt in filter_methods) {
            for (mdl in models) {
              combo_idx <- combo_idx + 1
              incProgress(0.6 / n_combos, detail = sprintf("Creating %s + %s (%d/%d)", flt, mdl, combo_idx, n_combos))

              # Handle Deep Learning learner configuration
              if (mdl == "mlp" || mdl == "tabtransformer") {
                # Parse hidden layers string "64,32" -> c(64, 32)
                hidden_layers <- tryCatch({
                  as.numeric(trimws(strsplit(input$dl_hidden, ",")[[1]]))
                }, error = function(e) c(64, 32))
                if (any(is.na(hidden_layers))) hidden_layers <- c(64, 32)
                
                if (mdl == "mlp") {
                  # Manually create the MLP learner with user settings
                  learner_obj <- create_mlp_learner(
                    hidden_layers = hidden_layers,
                    epochs = input$dl_epochs,
                    dropout = input$dl_dropout,
                    device = input$dl_device
                  )
                  
                  # Use the object instead of the string name
                  learner <- pipeline$create_graph_learner(
                    filter = flt,
                    model = learner_obj, # Pass the object!
                    n_features = input$n_features,
                    impute_method = input$impute_method,
                    scale = input$scale_features,
                    oversample = oversample,
                    screening = input$screening,
                    screening_frac = input$screening_frac,
                    batch_correct = input$batch_correct,
                    autoencoder = ae_params
                  )
                  
                  # Explicitly set learner ID
                  learner$id <- paste(flt, mdl, input$n_features, sep = "_")
                  
                  learners[[length(learners) + 1]] <- learner
                  next # Skip the standard creation call below
                }
                 # Future: TabTransformer specific handling if needed
              }

              # Standard learner creation
              learner <- pipeline$create_graph_learner(
                filter = flt,
                model = mdl,
                n_features = input$n_features,
                impute_method = input$impute_method,
                scale = input$scale_features,
                oversample = oversample,
                screening = input$screening,
                screening_frac = input$screening_frac,
                batch_correct = input$batch_correct,
                autoencoder = ae_params
              )
              
              # Explicitly set learner ID for reliable retrieval in Results module
              learner$id <- paste(flt, mdl, input$n_features, sep = "_")
              
              learners[[length(learners) + 1]] <- learner
            }
          }

          incProgress(0.2, detail = "Finalizing...")

          proj$set_pipeline(pipeline)
          proj$set_learner(learners)  # Now stores a list of learners

          showNotification(sprintf("Created %d learner combinations!", n_combos), type = "message")
        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })

    # View code modal
    observeEvent(input$view_code, {
      proj <- project()
      code <- if (!is.null(proj)) proj$generate_code() else "# No project loaded"

      showModal(modalDialog(
        title = "Generated R Code",
        tags$pre(code, style = "background: #1e293b; color: #f8fafc; padding: 1rem; border-radius: 0.5rem; overflow-x: auto;"),
        footer = tagList(
          downloadButton(ns("download_code"), "Download", class = "btn-primary"),
          modalButton("Close")
        ),
        size = "l"
      ))
    })

    output$download_code <- downloadHandler(
      filename = function() "omicselector_pipeline.R",
      content = function(file) {
        proj <- project()
        code <- if (!is.null(proj)) proj$generate_code() else "# No project"
        writeLines(code, file)
      }
    )

    # Continue to benchmark
    observeEvent(input$next_step, {
      proj <- project()
      if (is.null(proj) || !proj$can_proceed("benchmark")) {
        showNotification("Please create pipeline first", type = "error")
        return()
      }
      switch_tab("benchmark")
    })
  })
}
