#' QC & EDA Module
#'
#' @description Quality control with traffic light system and exploratory analysis

# Module UI
mod_qc_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # Step indicator
    div(
      class = "step-header mb-4",
      h3(bsicons::bs_icon("clipboard-check"), "Step 2: Quality Control & EDA"),
      p(class = "text-muted", "Validate data quality before proceeding to analysis")
    ),

    # Traffic Light QC Panel
    card(
      card_header(
        div(class = "d-flex align-items-center",
            bsicons::bs_icon("stoplights"),
            span("Quality Control Checks", class = "ms-2"))
      ),
      card_body(
        actionButton(ns("run_qc"), "Run QC Checks",
                     class = "btn-warning mb-3",
                     icon = icon("stethoscope")),

        uiOutput(ns("qc_traffic_lights")),

        hr(),

        div(
          class = "d-flex justify-content-between align-items-center",
          uiOutput(ns("qc_overall_status")),
          actionButton(ns("next_step"), "Continue to Pipeline",
                       class = "btn-outline-success",
                       icon = icon("arrow-right"))
        )
      )
    ),

    # EDA Tabs
    navset_card_tab(
      title = "Exploratory Data Analysis",

      # PCA Tab
      nav_panel(
        "PCA",
        icon = bsicons::bs_icon("graph-up"),
        layout_column_wrap(
          width = 1/2,
          fill = FALSE,

          card(
            card_header("PCA 2D"),
            card_body(
              uiOutput(ns("pca_color_ui")),
              actionButton(ns("plot_pca"), "Update PCA", class = "btn-sm btn-outline-primary mb-2"),
              plotOutput(ns("pca_2d"), height = "350px")
            )
          ),

          card(
            card_header("PCA 3D (Interactive)"),
            card_body(
              plotly::plotlyOutput(ns("pca_3d"), height = "400px")
            )
          )
        )
      ),

      # Distribution Tab
      nav_panel(
        "Distributions",
        icon = bsicons::bs_icon("bar-chart"),
        layout_column_wrap(
          width = 1/2,
          fill = FALSE,

          card(
            card_header("Target Distribution"),
            card_body(
              plotOutput(ns("target_dist"), height = "250px")
            )
          ),

          card(
            card_header("Missing Values"),
            card_body(
              plotOutput(ns("missing_plot"), height = "250px")
            )
          )
        )
      ),

      # Heatmap Tab (New Feature)
      nav_panel(
        "Heatmap",
        icon = bsicons::bs_icon("grid-3x3"),
        card(
          card_body(
            layout_column_wrap(
              width = 1/4,
              numericInput(ns("heatmap_top_n"), "Top N Variable Features",
                           value = 50, min = 10, max = 500, step = 10),
              selectInput(ns("heatmap_scale"), "Scale",
                          choices = c("Row (Z-score)" = "row",
                                      "Column" = "column",
                                      "None" = "none")),
              selectInput(ns("heatmap_cluster"), "Clustering",
                          choices = c("Both" = "both",
                                      "Rows only" = "row",
                                      "Columns only" = "column",
                                      "None" = "none")),
              actionButton(ns("plot_heatmap"), "Generate Heatmap",
                           class = "btn-primary mt-4")
            ),
            plotOutput(ns("heatmap_plot"), height = "500px")
          )
        )
      ),

      # Missing Values Analysis Tab (New Feature)
      nav_panel(
        "Missing Values",
        icon = bsicons::bs_icon("dash-circle"),
        layout_columns(
          col_widths = c(4, 8),
          card(
            card_header("Missing Value Analysis"),
            card_body(
              uiOutput(ns("missing_summary")),
              hr(),
              h6("Imputation"),
              selectInput(ns("impute_method"), "Method",
                          choices = c("None" = "none",
                                      "Mean" = "mean",
                                      "Median" = "median",
                                      "KNN" = "knn",
                                      "PMM (Advanced)" = "pmm")),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'knn'", ns("impute_method")),
                numericInput(ns("impute_knn_k"), "K Neighbors", value = 5, min = 1, max = 20)
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'pmm'", ns("impute_method")),
                div(class = "alert alert-info small mb-2", role = "alert",
                    bsicons::bs_icon("info-circle"),
                    " PMM (Predictive Mean Matching) uses single imputation (m=1) to preserve ",
                    "data distribution. Requires 'mice' package. Slower on high-dimensional data. ",
                    "Note: This is NOT multiple imputation with pooling."),
                numericInput(ns("impute_pmm_maxit"), "Iterations", value = 5, min = 1, max = 20)
              ),
              actionButton(ns("apply_imputation"), "Apply Imputation",
                           class = "btn-warning w-100",
                           icon = icon("magic")),
              hr(),
              actionButton(ns("reset_imputation"), "Reset to Original",
                           class = "btn-outline-secondary w-100",
                           icon = icon("undo"))
            )
          ),
          card(
            card_header("Missing Pattern"),
            card_body(
              plotOutput(ns("missing_pattern_plot"), height = "400px")
            )
          )
        ),

        layout_column_wrap(
          width = 1/2,
          fill = FALSE,
          class = "mt-3",
          card(
            card_header("Missing by Feature (Top 30)"),
            card_body(
              plotOutput(ns("missing_by_feature"), height = "300px")
            )
          ),
          card(
            card_header("Missing by Sample"),
            card_body(
              plotOutput(ns("missing_by_sample"), height = "300px")
            )
          )
        )
      ),

      # Correlation Tab
      nav_panel(
        "Correlation",
        icon = bsicons::bs_icon("diagram-2"),
        card(
          card_body(
            layout_column_wrap(
              width = 1/4,
              uiOutput(ns("cor_x_ui")),
              uiOutput(ns("cor_y_ui")),
              actionButton(ns("plot_cor"), "Plot", class = "btn-primary mt-4")
            ),
            plotOutput(ns("correlation_plot"), height = "400px")
          )
        )
      ),

      # DE Tab
      nav_panel(
        "Differential Expression",
        icon = bsicons::bs_icon("sort-numeric-down"),
        layout_columns(
          col_widths = c(3, 9),
          card(
            card_header("DE Settings"),
            card_body(
              selectInput(ns("de_mode"), "Value Type",
                          choices = c("Log TPM" = "logtpm", "Delta Ct" = "deltact")),
              numericInput(ns("de_pvalue"), "P-value Threshold", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
              numericInput(ns("de_log2fc"), "Log2 FC Threshold", value = 1, min = 0.1, max = 5, step = 0.1),
              numericInput(ns("de_top_n"), "Label Top N", value = 10, min = 0, max = 50),
              actionButton(ns("run_de"), "Run DE Analysis",
                           class = "btn-warning w-100",
                           icon = icon("calculator")),
              hr(),
              downloadButton(ns("download_de"), "Download Table", class = "btn-outline-secondary w-100")
            )
          ),

          card(
            card_header("Volcano Plot"),
            card_body(
              plotOutput(ns("volcano_plot"), height = "450px")
            )
          )
        ),

        card(
          class = "mt-3",
          card_header("DE Results Table"),
          card_body(
            div(style = "max-height: 300px; overflow-y: auto;",
                tableOutput(ns("de_table")))
          )
        )
      )
    ),

    # Help/FAQ Card
    card(
      class = "mt-4",
      card_header(
        bsicons::bs_icon("question-circle"),
        span("Help & FAQ", class = "ms-2")
      ),
      card_body(
        accordion(
          id = ns("qc_help_accordion"),
          open = FALSE,

          accordion_panel(
            "What do the QC traffic lights mean?",
            icon = bsicons::bs_icon("stoplights"),
            p("The QC system uses a traffic light metaphor:"),
            tags$ul(
              tags$li(tags$strong("Green:"), " Check passed - no issues detected."),
              tags$li(tags$strong("Yellow:"), " Warning - proceed with caution, results may be affected."),
              tags$li(tags$strong("Red:"), " Issue detected - may block progression if critical.")
            ),
            p("Key checks include: class balance, sample size per class, feature/sample ratio (p/n),
               missing values, zero-variance features, and patient grouping for leakage prevention.")
          ),

          accordion_panel(
            "Why is patient grouping important?",
            icon = bsicons::bs_icon("people"),
            p("If your data contains multiple samples from the same patient (e.g., repeated measures,
               longitudinal samples, or multiple biopsies), you MUST set the Patient ID column."),
            p("Without patient grouping, cross-validation may place samples from the same patient in
               both training and test sets, causing ", tags$strong("data leakage"), " and overly
               optimistic performance estimates."),
            p(class = "text-danger", "Leakage can make a model appear to work well when it actually
               won't generalize to new patients.")
          ),

          accordion_panel(
            "Which imputation method should I use?",
            icon = bsicons::bs_icon("magic"),
            tags$ul(
              tags$li(tags$strong("Median:"), " Fast, robust to outliers. Good default for most cases."),
              tags$li(tags$strong("Mean:"), " Simple, but sensitive to outliers."),
              tags$li(tags$strong("KNN:"), " Uses similar samples to impute. Good when features are correlated."),
              tags$li(tags$strong("PMM (Advanced):"), " Predictive Mean Matching preserves data distribution.
                      Best for maintaining realistic values, but slower. Uses single imputation (m=1),
                      not full multiple imputation with pooling.")
            ),
            p("For high-dimensional data with >20% missing values, consider removing features with
               excessive missingness before imputation.")
          ),

          accordion_panel(
            "What is p/n ratio and why does it matter?",
            icon = bsicons::bs_icon("graph-up-arrow"),
            p("The p/n ratio is the number of features (p) divided by the number of samples (n)."),
            tags$ul(
              tags$li("When p/n > 10, you have 'high-dimensional' data."),
              tags$li("When p/n > 50, overfitting risk is very high."),
              tags$li("Many ML models struggle when p >> n.")
            ),
            p("Solutions: Use feature filtering (in Pipeline step), apply regularization,
               or reduce features using domain knowledge or prefix filtering.")
          ),

          accordion_panel(
            "How do I interpret the PCA plot?",
            icon = bsicons::bs_icon("scatter-chart"),
            p("PCA (Principal Component Analysis) reduces high-dimensional data to 2-3 dimensions
               for visualization."),
            tags$ul(
              tags$li("Samples that cluster together are similar in their feature profiles."),
              tags$li("Clear separation by target class suggests the features contain discriminative information."),
              tags$li("Batch effects may cause samples to cluster by batch rather than biology."),
              tags$li("Outliers may appear as isolated points far from other samples.")
            )
          ),

          accordion_panel(
            "What is Differential Expression (DE) analysis?",
            icon = bsicons::bs_icon("sort-numeric-down"),
            p("DE analysis identifies features that differ significantly between classes
               using t-tests. This is useful for:"),
            tags$ul(
              tags$li("Exploring which features distinguish classes."),
              tags$li("Initial feature filtering before ML."),
              tags$li("Biological interpretation of results.")
            ),
            p(class = "text-warning", "Note: DE is univariate and doesn't account for feature
               interactions. ML-based feature selection (in Pipeline step) often selects different
               features than DE.")
          )
        )
      )
    )
  )
}

# Module Server
mod_qc_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive trigger to force UI updates when R6 object state changes
    qc_trigger <- reactiveVal(0)

    # Reactive for DE results
    de_results <- reactiveVal(NULL)

    # Run QC checks
    observeEvent(input$run_qc, {
      proj <- project()
      req(proj)

      message("[mod_qc] Running QC checks...")
      withProgress(message = "Running QC checks...", value = 0.5, {
        proj$run_qc()
      })

      # Trigger UI refresh
      qc_trigger(qc_trigger() + 1)
      message("[mod_qc] QC complete, trigger updated")
      showNotification("QC checks complete!", type = "message")
    })

    # Traffic light display
    output$qc_traffic_lights <- renderUI({
      # Depend on trigger to refresh when QC runs
      qc_trigger()
      proj <- project()
      if (is.null(proj)) return(NULL)

      qc <- proj$get_qc_status()
      if (is.null(qc)) {
        return(div(class = "text-muted", "Click 'Run QC Checks' to validate your data"))
      }

      # Create traffic light items
      lights <- lapply(names(qc), function(check_name) {
        check <- qc[[check_name]]
        status <- check$status
        message <- check$message
        blocking <- check$blocking

        icon_name <- switch(status,
                            "green" = "check-circle-fill",
                            "yellow" = "exclamation-triangle-fill",
                            "red" = "x-circle-fill",
                            "question-circle")

        icon_class <- switch(status,
                             "green" = "text-success",
                             "yellow" = "text-warning",
                             "red" = "text-danger",
                             "text-secondary")

        blocking_badge <- if (blocking && status == "red") {
          span(class = "badge bg-danger ms-2", "Blocking")
        } else NULL

        div(
          class = "d-flex align-items-center mb-2 p-2 border rounded",
          bsicons::bs_icon(icon_name, class = paste(icon_class, "me-2"), size = "1.5em"),
          div(
            strong(gsub("_", " ", check_name)),
            blocking_badge,
            br(),
            span(class = "text-muted small", message)
          )
        )
      })

      do.call(tagList, lights)
    })

    # Overall QC status
    output$qc_overall_status <- renderUI({
      # Depend on trigger to refresh when QC runs
      qc_trigger()
      proj <- project()
      if (is.null(proj)) return(NULL)

      qc <- proj$get_qc_status()
      if (is.null(qc)) {
        return(span(class = "badge bg-secondary", "QC Not Run"))
      }

      summary <- proj$summary()
      overall <- summary$qc_overall

      badge_class <- switch(overall,
                            "green" = "bg-success",
                            "yellow" = "bg-warning text-dark",
                            "red" = "bg-danger",
                            "bg-secondary")

      badge_text <- switch(overall,
                           "green" = "All Checks Passed",
                           "yellow" = "Warnings Present",
                           "red" = "Issues Found",
                           "Unknown")

      span(class = paste("badge", badge_class, "fs-6"), badge_text)
    })

    # Helper to get feature matrix (supports both single and multi-omics)
    get_feature_matrix <- function() {
      proj <- project()
      if (is.null(proj)) return(NULL)

      data <- proj$get_data()
      mapping <- proj$get_mapping()

      if (is.null(data) || is.null(mapping$target)) return(NULL)

      meta_cols <- c(mapping$target, mapping$patient_id, mapping$batch, "mix")
      meta_cols <- meta_cols[!is.null(meta_cols) & nzchar(meta_cols)]

      # Helper to process a single data frame
      process_df <- function(df, prefix_filter = TRUE) {
        keep <- setdiff(names(df), meta_cols)
        x <- df[, keep, drop = FALSE]

        # Filter by prefix if requested
        if (prefix_filter && !is.null(mapping$feature_prefix) && nzchar(mapping$feature_prefix)) {
          escaped_prefix <- gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", mapping$feature_prefix)
          keep_cols <- grep(paste0("^", escaped_prefix), names(x), value = TRUE)
          if (length(keep_cols) > 0) {
            x <- x[, keep_cols, drop = FALSE]
          } else {
            x <- x[, character(0), drop = FALSE]
          }
        }

        # Keep only numeric columns
        numeric_cols <- sapply(x, is.numeric)
        x[, numeric_cols, drop = FALSE]
      }

      if (is.data.frame(data)) {
        process_df(data)
      } else if (is.list(data)) {
        # Multi-omics: combine feature matrices from all modalities
        # Try to align by patient_id if available
        patient_col <- mapping$patient_id

        if (!is.null(patient_col) && nzchar(patient_col)) {
          # Align modalities by patient ID
          all_ids <- NULL
          for (mod in names(data)) {
            if (patient_col %in% names(data[[mod]])) {
              mod_ids <- data[[mod]][[patient_col]]
              if (is.null(all_ids)) {
                all_ids <- mod_ids
              } else {
                all_ids <- intersect(all_ids, mod_ids)
              }
            }
          }

          if (!is.null(all_ids) && length(all_ids) > 0) {
            # Merge modalities on common patient IDs
            combined_dfs <- lapply(names(data), function(mod_name) {
              mod_data <- data[[mod_name]]
              if (!(patient_col %in% names(mod_data))) return(NULL)

              # Subset and order by patient ID
              mod_data <- mod_data[mod_data[[patient_col]] %in% all_ids, , drop = FALSE]
              mod_data <- mod_data[order(match(mod_data[[patient_col]], all_ids)), , drop = FALSE]

              x <- process_df(mod_data, prefix_filter = TRUE)
              if (ncol(x) > 0) {
                names(x) <- paste(mod_name, names(x), sep = "::")
              }
              x
            })
            combined_dfs <- combined_dfs[!sapply(combined_dfs, is.null)]
            if (length(combined_dfs) > 0) {
              return(do.call(cbind, combined_dfs))
            }
          }
        }

        # Fallback: combine assuming same row order (with warning)
        combined_dfs <- lapply(names(data), function(mod_name) {
          x <- process_df(data[[mod_name]], prefix_filter = TRUE)
          if (ncol(x) > 0) {
            names(x) <- paste(mod_name, names(x), sep = "::")
          }
          x
        })
        if (length(combined_dfs) > 0 && nrow(combined_dfs[[1]]) > 0) {
          return(do.call(cbind, combined_dfs))
        }
        NULL
      } else {
        NULL
      }
    }

    # Helper to get target vector (supports both single and multi-omics)
    get_target_vector <- function() {
      proj <- project()
      if (is.null(proj)) return(NULL)

      data <- proj$get_data()
      mapping <- proj$get_mapping()

      if (is.null(data) || is.null(mapping$target)) return(NULL)

      if (is.data.frame(data)) {
        if (mapping$target %in% names(data)) {
          return(data[[mapping$target]])
        }
      } else if (is.list(data)) {
        # Multi-omics: search for target in all modalities
        for (mod in names(data)) {
          if (mapping$target %in% names(data[[mod]])) {
            return(data[[mod]][[mapping$target]])
          }
        }
      }
      NULL
    }

    # Column choices for plots - P1 Fix: Support multi-omics data
    columns <- reactive({
      proj <- project()
      if (is.null(proj)) return(character(0))
      data <- proj$get_data()
      if (is.null(data)) return(character(0))

      # Handle both single data.frame and multi-omics list
      if (is.data.frame(data)) {
        names(data)
      } else if (is.list(data)) {
        # Multi-omics: collect unique column names from all modalities
        unique(unlist(lapply(data, names)))
      } else {
        character(0)
      }
    })

    # PCA color UI
    output$pca_color_ui <- renderUI({
      cols <- columns()
      proj <- project()
      if (is.null(proj)) return(NULL)

      mapping <- proj$get_mapping()
      default <- mapping$target %||% ""

      selectInput(ns("pca_color"), "Color By",
                  choices = c("None" = "", cols),
                  selected = default)
    })

    # PCA 2D plot
    pca_data <- eventReactive(input$plot_pca, {
      x <- get_feature_matrix()
      if (is.null(x) || ncol(x) < 2 || nrow(x) < 2) return(NULL)

      # Handle missing values with user notification
      original_n <- nrow(x)
      complete_idx <- complete.cases(x)
      x <- x[complete_idx, , drop = FALSE]
      dropped_n <- original_n - nrow(x)

      if (dropped_n > 0) {
        pct_dropped <- round(dropped_n / original_n * 100, 1)
        showNotification(
          sprintf("PCA: Excluded %d samples (%.1f%%) with missing values", dropped_n, pct_dropped),
          type = if (pct_dropped > 20) "warning" else "message",
          duration = 5
        )
      }

      if (nrow(x) < 2) {
        showNotification("Not enough complete samples for PCA", type = "error")
        return(NULL)
      }

      pca <- prcomp(x, scale. = TRUE)

      proj <- project()
      data <- proj$get_data()

      # P1 Fix: Handle PCA grouping for both single and multi-omics data
      groups <- NULL
      if (!is.null(input$pca_color) && nzchar(input$pca_color)) {
        if (is.data.frame(data)) {
          # Single data.frame case
          if (input$pca_color %in% names(data)) {
            groups <- data[[input$pca_color]][complete_idx]
          }
        } else if (is.list(data)) {
          # Multi-omics: find the column in any modality
          for (mod in names(data)) {
            if (input$pca_color %in% names(data[[mod]])) {
              # Get the grouping vector and apply the same complete_idx
              groups <- data[[mod]][[input$pca_color]][complete_idx]
              break
            }
          }
        }
      }

      list(scores = pca$x, groups = groups, var_exp = summary(pca)$importance[2, 1:2],
           n_samples = nrow(x), n_dropped = dropped_n)
    })

    output$pca_2d <- renderPlot({
      pca <- pca_data()
      if (is.null(pca)) {
        plot.new()
        text(0.5, 0.5, "Not enough data for PCA", cex = 1.2)
        return()
      }

      scores <- pca$scores[, 1:2]
      var_exp <- pca$var_exp * 100

      if (is.null(pca$groups)) {
        plot(scores, pch = 19, col = "#2563eb",
             xlab = sprintf("PC1 (%.1f%%)", var_exp[1]),
             ylab = sprintf("PC2 (%.1f%%)", var_exp[2]),
             main = "PCA")
      } else {
        groups <- as.factor(pca$groups)
        cols <- rainbow(length(levels(groups)))
        plot(scores, pch = 19, col = cols[groups],
             xlab = sprintf("PC1 (%.1f%%)", var_exp[1]),
             ylab = sprintf("PC2 (%.1f%%)", var_exp[2]),
             main = "PCA")
        legend("topright", legend = levels(groups), col = cols, pch = 19, cex = 0.8)
      }
    })

    # PCA 3D plot - P1 Fix: Support multi-omics data
    output$pca_3d <- plotly::renderPlotly({
      x <- get_feature_matrix()
      if (is.null(x) || ncol(x) < 3 || nrow(x) < 3) return(NULL)

      # Track complete cases for consistent indexing
      complete_idx <- complete.cases(x)
      x <- x[complete_idx, , drop = FALSE]
      if (nrow(x) < 3) return(NULL)

      pca <- prcomp(x, scale. = TRUE)
      scores <- as.data.frame(pca$x[, 1:3])

      proj <- project()
      data <- proj$get_data()

      # P1 Fix: Handle grouping for both single and multi-omics data
      scores$group <- "All"
      if (!is.null(input$pca_color) && nzchar(input$pca_color)) {
        if (is.data.frame(data)) {
          if (input$pca_color %in% names(data)) {
            scores$group <- as.factor(data[[input$pca_color]][complete_idx])
          }
        } else if (is.list(data)) {
          # Multi-omics: find the column in any modality
          for (mod in names(data)) {
            if (input$pca_color %in% names(data[[mod]])) {
              scores$group <- as.factor(data[[mod]][[input$pca_color]][complete_idx])
              break
            }
          }
        }
      }

      plotly::plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3,
                      color = ~group, colors = "Set2",
                      type = "scatter3d", mode = "markers",
                      marker = list(size = 5)) |>
        plotly::layout(title = "PCA 3D")
    })

    # Target distribution
    output$target_dist <- renderPlot({
      target <- get_target_vector()
      if (is.null(target)) {
        plot.new()
        text(0.5, 0.5, "No target variable", cex = 1.2)
        return()
      }

      tbl <- table(target)
      barplot(tbl, col = c("#10b981", "#ef4444"),
              main = "Target Distribution",
              ylab = "Count")
    })

    # Missing values plot
    output$missing_plot <- renderPlot({
      x <- get_feature_matrix()
      if (is.null(x)) {
        plot.new()
        text(0.5, 0.5, "No data", cex = 1.2)
        return()
      }

      missing_pct <- colMeans(is.na(x)) * 100
      missing_pct <- sort(missing_pct, decreasing = TRUE)
      missing_pct <- head(missing_pct[missing_pct > 0], 20)

      if (length(missing_pct) == 0) {
        plot.new()
        text(0.5, 0.5, "No missing values!", cex = 1.2, col = "#10b981")
        return()
      }

      barplot(missing_pct, las = 2, col = "#f59e0b",
              main = "Top Missing Features (%)",
              ylab = "% Missing", cex.names = 0.7)
    })

    # Correlation UI
    output$cor_x_ui <- renderUI({
      x <- get_feature_matrix()
      if (is.null(x)) return(NULL)
      selectInput(ns("cor_x"), "X Variable", choices = names(x))
    })

    output$cor_y_ui <- renderUI({
      x <- get_feature_matrix()
      if (is.null(x)) return(NULL)
      default <- if (ncol(x) >= 2) names(x)[2] else names(x)[1]
      selectInput(ns("cor_y"), "Y Variable", choices = names(x), selected = default)
    })

    # Correlation plot - use eventReactive pattern for better performance
    cor_plot_trigger <- reactiveVal(0)

    observeEvent(input$plot_cor, {
      cor_plot_trigger(cor_plot_trigger() + 1)
    })

    output$correlation_plot <- renderPlot({
      # Depend on trigger to update only when button clicked
      cor_plot_trigger()

      x <- get_feature_matrix()
      req(x, input$cor_x, input$cor_y)

      x_vals <- x[[input$cor_x]]
      y_vals <- x[[input$cor_y]]

      cor_val <- cor(x_vals, y_vals, use = "complete.obs")

      plot(x_vals, y_vals, pch = 19, col = "#2563eb80",
           xlab = input$cor_x, ylab = input$cor_y,
           main = sprintf("Correlation: %.3f", cor_val))
      abline(lm(y_vals ~ x_vals), col = "#ef4444", lwd = 2)
    })

    # DE Analysis
    observeEvent(input$run_de, {
      x <- get_feature_matrix()
      target <- get_target_vector()

      if (is.null(x) || is.null(target)) {
        showNotification("Data or target not available", type = "error")
        return()
      }

      if (length(unique(target)) != 2) {
        showNotification("DE requires binary target", type = "error")
        return()
      }

      withProgress(message = "Running DE analysis...", value = 0.5, {
        classes <- as.factor(target)
        case_idx <- classes == levels(classes)[2]
        ctrl_idx <- classes == levels(classes)[1]

        # Vectorized computation of log2FC using matrix operations
        x_matrix <- as.matrix(x)
        case_means <- colMeans(x_matrix[case_idx, , drop = FALSE], na.rm = TRUE)
        ctrl_means <- colMeans(x_matrix[ctrl_idx, , drop = FALSE], na.rm = TRUE)
        log2fc_vals <- case_means - ctrl_means

        # Vectorized p-value computation using apply (more efficient than sapply with function)
        # Pre-split the data for better performance
        case_data <- x_matrix[case_idx, , drop = FALSE]
        ctrl_data <- x_matrix[ctrl_idx, , drop = FALSE]

        pvalue_vals <- vapply(seq_len(ncol(x_matrix)), function(i) {
          tryCatch(
            t.test(case_data[, i], ctrl_data[, i])$p.value,
            error = function(e) NA_real_
          )
        }, FUN.VALUE = numeric(1))

        results <- data.frame(
          feature = names(x),
          log2FC = log2fc_vals,
          pvalue = pvalue_vals,
          stringsAsFactors = FALSE
        )

        results$padj <- p.adjust(results$pvalue, method = "BH")
        results$neg_log10_p <- -log10(results$pvalue)

        if (input$de_mode == "deltact") {
          results$log2FC <- -results$log2FC
        }

        de_results(results[order(results$pvalue), ])
      })

      showNotification("DE analysis complete!", type = "message")
    })

    # Volcano plot
    output$volcano_plot <- renderPlot({
      de <- de_results()
      if (is.null(de)) {
        plot.new()
        text(0.5, 0.5, "Run DE analysis first", cex = 1.2)
        return()
      }

      # Use configurable thresholds
      pval_threshold <- input$de_pvalue %||% 0.05
      log2fc_threshold <- input$de_log2fc %||% 1

      plot(de$log2FC, de$neg_log10_p,
           pch = 19, col = "#64748b80",
           xlab = "log2 Fold Change",
           ylab = "-log10(p-value)",
           main = sprintf("Volcano Plot (p < %.3f, |log2FC| > %.1f)", pval_threshold, log2fc_threshold))

      # Highlight significant using user-defined thresholds
      sig <- de$padj < pval_threshold & abs(de$log2FC) > log2fc_threshold
      points(de$log2FC[sig], de$neg_log10_p[sig], pch = 19, col = "#ef4444")

      # Label top N
      top_n <- input$de_top_n
      if (top_n > 0 && nrow(de) >= top_n) {
        top_de <- head(de, top_n)
        text(top_de$log2FC, top_de$neg_log10_p,
             labels = top_de$feature, pos = 3, cex = 0.7)
      }

      # Draw threshold lines using user-defined values
      abline(h = -log10(pval_threshold), col = "gray", lty = 2)
      abline(v = c(-log2fc_threshold, log2fc_threshold), col = "gray", lty = 2)

      # Add legend showing counts
      n_sig <- sum(sig, na.rm = TRUE)
      legend("topright",
             legend = c(sprintf("Significant: %d", n_sig), sprintf("Non-significant: %d", sum(!sig, na.rm = TRUE))),
             col = c("#ef4444", "#64748b80"),
             pch = 19, cex = 0.8)
    })

    # DE table
    output$de_table <- renderTable({
      de <- de_results()
      if (is.null(de)) return(NULL)
      head(de, 50)
    })

    # Download DE
    output$download_de <- downloadHandler(
      filename = function() paste0("de_results_", Sys.Date(), ".csv"),
      content = function(file) {
        de <- de_results()
        if (!is.null(de)) {
          write.csv(de, file, row.names = FALSE)
        }
      }
    )

    # =========================================================================
    # HEATMAP (New Feature)
    # =========================================================================

    heatmap_trigger <- reactiveVal(0)

    observeEvent(input$plot_heatmap, {
      heatmap_trigger(heatmap_trigger() + 1)
    })

    output$heatmap_plot <- renderPlot({
      heatmap_trigger()

      x <- get_feature_matrix()
      if (is.null(x) || ncol(x) < 2) {
        plot.new()
        text(0.5, 0.5, "Not enough features for heatmap", cex = 1.2)
        return()
      }

      # Use complete cases only
      complete_idx <- complete.cases(x)
      x <- x[complete_idx, , drop = FALSE]

      if (nrow(x) < 2) {
        plot.new()
        text(0.5, 0.5, "Not enough complete samples", cex = 1.2)
        return()
      }

      # Select top variable features
      top_n <- min(input$heatmap_top_n %||% 50, ncol(x))
      feature_vars <- apply(x, 2, var, na.rm = TRUE)
      top_features <- names(sort(feature_vars, decreasing = TRUE))[1:top_n]
      x <- x[, top_features, drop = FALSE]

      # Convert to matrix
      mat <- as.matrix(x)

      # Scale if requested
      scale_opt <- input$heatmap_scale %||% "row"
      if (scale_opt == "row") {
        mat <- t(scale(t(mat)))
      } else if (scale_opt == "column") {
        mat <- scale(mat)
      }

      # Replace any NaN/Inf from scaling
      mat[!is.finite(mat)] <- 0

      # Determine clustering
      cluster_opt <- input$heatmap_cluster %||% "both"
      cluster_rows <- cluster_opt %in% c("both", "row")
      cluster_cols <- cluster_opt %in% c("both", "column")

      # Get target for row annotation if available
      target <- get_target_vector()
      row_colors <- NULL
      if (!is.null(target)) {
        target <- target[complete_idx]
        target_factor <- as.factor(target)
        target_cols <- c("#10b981", "#ef4444")[as.numeric(target_factor)]
        row_colors <- target_cols
      }

      # Create heatmap using base R heatmap function
      tryCatch({
        if (!is.null(row_colors)) {
          heatmap(mat,
                  Rowv = if (cluster_rows) as.dendrogram(hclust(dist(mat))) else NA,
                  Colv = if (cluster_cols) as.dendrogram(hclust(dist(t(mat)))) else NA,
                  scale = "none",  # We already scaled
                  col = colorRampPalette(c("#2563eb", "white", "#ef4444"))(100),
                  RowSideColors = row_colors,
                  margins = c(8, 4),
                  labRow = NA,  # Hide row labels for large datasets
                  cexCol = 0.6)
        } else {
          heatmap(mat,
                  Rowv = if (cluster_rows) as.dendrogram(hclust(dist(mat))) else NA,
                  Colv = if (cluster_cols) as.dendrogram(hclust(dist(t(mat)))) else NA,
                  scale = "none",
                  col = colorRampPalette(c("#2563eb", "white", "#ef4444"))(100),
                  margins = c(8, 4),
                  labRow = NA,
                  cexCol = 0.6)
        }
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error generating heatmap:", e$message), cex = 0.9)
      })
    })

    # =========================================================================
    # MISSING VALUES ANALYSIS (New Feature)
    # =========================================================================

    # Missing summary statistics
    output$missing_summary <- renderUI({
      x <- get_feature_matrix()
      if (is.null(x)) {
        return(p(class = "text-muted", "No data available"))
      }

      total_cells <- prod(dim(x))
      missing_cells <- sum(is.na(x))
      missing_pct <- round(missing_cells / total_cells * 100, 2)

      features_with_missing <- sum(colSums(is.na(x)) > 0)
      samples_with_missing <- sum(rowSums(is.na(x)) > 0)

      complete_samples <- sum(complete.cases(x))

      tagList(
        div(class = "mb-2",
            strong("Total Cells: "), format(total_cells, big.mark = ",")),
        div(class = "mb-2",
            strong("Missing Cells: "), format(missing_cells, big.mark = ","),
            sprintf(" (%.2f%%)", missing_pct)),
        div(class = "mb-2",
            strong("Features with Missing: "), features_with_missing,
            sprintf(" / %d", ncol(x))),
        div(class = "mb-2",
            strong("Samples with Missing: "), samples_with_missing,
            sprintf(" / %d", nrow(x))),
        div(class = "mb-2",
            strong("Complete Samples: "), complete_samples,
            sprintf(" (%.1f%%)", complete_samples / nrow(x) * 100)),
        if (missing_pct == 0) {
          div(class = "alert alert-success mt-3 mb-0", role = "alert",
              bsicons::bs_icon("check-circle"), " No missing values detected!")
        } else if (missing_pct > 20) {
          div(class = "alert alert-danger mt-3 mb-0", role = "alert",
              bsicons::bs_icon("exclamation-triangle"),
              " High missing rate. Consider imputation or feature filtering.")
        } else if (missing_pct > 5) {
          div(class = "alert alert-warning mt-3 mb-0", role = "alert",
              bsicons::bs_icon("exclamation-circle"),
              " Moderate missing rate. Imputation recommended.")
        } else {
          div(class = "alert alert-info mt-3 mb-0", role = "alert",
              bsicons::bs_icon("info-circle"),
              " Low missing rate. Imputation optional.")
        }
      )
    })

    # Missing pattern plot (heatmap of missing vs observed)
    output$missing_pattern_plot <- renderPlot({
      x <- get_feature_matrix()
      if (is.null(x) || ncol(x) < 1) {
        plot.new()
        text(0.5, 0.5, "No data", cex = 1.2)
        return()
      }

      # Create binary missing matrix
      missing_mat <- is.na(x)

      if (sum(missing_mat) == 0) {
        plot.new()
        text(0.5, 0.5, "No missing values!", cex = 1.2, col = "#10b981")
        return()
      }

      # Limit to features with missing values for clarity
      cols_with_missing <- colSums(missing_mat) > 0
      if (sum(cols_with_missing) == 0) {
        plot.new()
        text(0.5, 0.5, "No missing values!", cex = 1.2, col = "#10b981")
        return()
      }

      missing_mat <- missing_mat[, cols_with_missing, drop = FALSE]

      # Limit to top 50 features with most missing for readability
      if (ncol(missing_mat) > 50) {
        top_missing <- order(colSums(missing_mat), decreasing = TRUE)[1:50]
        missing_mat <- missing_mat[, top_missing, drop = FALSE]
      }

      # Convert to numeric for image plot
      mat_numeric <- matrix(as.numeric(missing_mat), nrow = nrow(missing_mat))
      colnames(mat_numeric) <- colnames(missing_mat)

      # Order rows by missingness pattern
      row_order <- order(rowSums(missing_mat), decreasing = TRUE)
      mat_numeric <- mat_numeric[row_order, , drop = FALSE]

      # Plot
      image(t(mat_numeric[nrow(mat_numeric):1, , drop = FALSE]),
            col = c("#e5e7eb", "#ef4444"),
            axes = FALSE,
            main = "Missing Pattern (red = missing)")

      # Add column labels if not too many
      if (ncol(mat_numeric) <= 30) {
        axis(1, at = seq(0, 1, length.out = ncol(mat_numeric)),
             labels = colnames(mat_numeric), las = 2, cex.axis = 0.6)
      }
    })

    # Missing by feature plot
    output$missing_by_feature <- renderPlot({
      x <- get_feature_matrix()
      if (is.null(x)) {
        plot.new()
        text(0.5, 0.5, "No data", cex = 1.2)
        return()
      }

      missing_pct <- colMeans(is.na(x)) * 100
      missing_pct <- sort(missing_pct, decreasing = TRUE)
      missing_pct <- head(missing_pct[missing_pct > 0], 30)

      if (length(missing_pct) == 0) {
        plot.new()
        text(0.5, 0.5, "No missing values!", cex = 1.2, col = "#10b981")
        return()
      }

      # Color code by severity
      bar_cols <- ifelse(missing_pct > 20, "#ef4444",
                         ifelse(missing_pct > 5, "#f59e0b", "#10b981"))

      barplot(missing_pct, las = 2, col = bar_cols,
              main = "Missing % by Feature",
              ylab = "% Missing", cex.names = 0.6)
      abline(h = c(5, 20), col = c("#f59e0b", "#ef4444"), lty = 2)
    })

    # Missing by sample plot
    output$missing_by_sample <- renderPlot({
      x <- get_feature_matrix()
      if (is.null(x)) {
        plot.new()
        text(0.5, 0.5, "No data", cex = 1.2)
        return()
      }

      missing_pct <- rowMeans(is.na(x)) * 100

      if (all(missing_pct == 0)) {
        plot.new()
        text(0.5, 0.5, "No missing values!", cex = 1.2, col = "#10b981")
        return()
      }

      # Histogram of sample missingness
      hist(missing_pct, breaks = 30, col = "#2563eb",
           main = "Distribution of Missing % per Sample",
           xlab = "% Missing per Sample",
           ylab = "Number of Samples")
      abline(v = c(5, 20), col = c("#f59e0b", "#ef4444"), lty = 2, lwd = 2)

      # Add legend
      legend("topright",
             legend = c("5% threshold", "20% threshold"),
             col = c("#f59e0b", "#ef4444"),
             lty = 2, lwd = 2, cex = 0.8)
    })

    # Apply imputation
    observeEvent(input$apply_imputation, {
      proj <- project()
      req(proj)

      method <- input$impute_method
      if (is.null(method) || method == "none") {
        showNotification("Please select an imputation method", type = "warning")
        return()
      }

      data <- proj$get_data()
      if (is.null(data)) {
        showNotification("No data available", type = "error")
        return()
      }

      # PMM requires mice package - check availability first
      if (method == "pmm") {
        if (!requireNamespace("mice", quietly = TRUE)) {
          showNotification(
            "PMM imputation requires the 'mice' package. Install with: install.packages('mice')",
            type = "error", duration = 10
          )
          return()
        }
      }

      withProgress(message = sprintf("Applying %s imputation...", method), value = 0.5, {

        # Get mapping to exclude metadata columns from imputation
        mapping <- proj$get_mapping()
        meta_cols <- c(mapping$target, mapping$patient_id, mapping$batch, "mix")
        meta_cols <- meta_cols[!is.null(meta_cols) & nzchar(meta_cols)]

        # PMM imputation function (handles df holistically, not per-column)
        impute_pmm <- function(df) {
          numeric_cols <- sapply(df, is.numeric)
          feature_cols <- setdiff(names(df)[numeric_cols], meta_cols)

          # Check if there are columns with missing values to impute
          cols_with_na <- feature_cols[sapply(df[, feature_cols, drop = FALSE], function(x) any(is.na(x)))]

          if (length(cols_with_na) == 0) {
            return(df)  # Nothing to impute
          }

          # Pre-validation: remove constant columns (mice fails on these)
          non_constant <- sapply(df[, cols_with_na, drop = FALSE], function(x) {
            vals <- x[!is.na(x)]
            length(unique(vals)) > 1
          })
          valid_cols <- cols_with_na[non_constant]

          if (length(valid_cols) == 0) {
            showNotification("No valid columns for PMM (all constant). Using median fallback.",
                             type = "warning")
            # Fallback to median for constant columns
            for (col in cols_with_na) {
              df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
            }
            return(df)
          }

          # Check minimum observations per column (need at least 10 donors for PMM)
          sufficient_obs <- sapply(df[, valid_cols, drop = FALSE], function(x) sum(!is.na(x)) >= 10)
          pmm_cols <- valid_cols[sufficient_obs]
          fallback_cols <- valid_cols[!sufficient_obs]

          # Handle columns with insufficient observations via median
          if (length(fallback_cols) > 0) {
            for (col in fallback_cols) {
              df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
            }
          }

          if (length(pmm_cols) == 0) {
            return(df)  # All handled by fallback
          }

          # Extract feature subset for PMM
          feature_df <- df[, pmm_cols, drop = FALSE]

          # Warn if high-dimensional (>500 features)
          if (ncol(feature_df) > 500) {
            showNotification(
              sprintf("PMM on %d features may be slow. Consider filtering features first.", ncol(feature_df)),
              type = "warning", duration = 5
            )
          }

          # Run PMM with mice
          tryCatch({
            # Use quickpred to limit predictors for high-dimensional data
            pred_matrix <- mice::quickpred(feature_df, mincor = 0.1, minpuc = 0.5)

            # Suppress mice output, single imputation (m=1), limited iterations
            maxit <- input$impute_pmm_maxit %||% 5
            imp <- mice::mice(
              feature_df,
              method = "pmm",
              predictorMatrix = pred_matrix,
              m = 1,
              maxit = maxit,
              printFlag = FALSE
            )

            # Extract completed dataset
            completed <- mice::complete(imp, 1)

            # Merge back into original df
            df[, pmm_cols] <- completed

            df
          }, error = function(e) {
            showNotification(
              sprintf("PMM failed: %s. Falling back to median imputation.", e$message),
              type = "warning", duration = 8
            )
            # Fallback to median
            for (col in pmm_cols) {
              if (any(is.na(df[[col]]))) {
                df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
              }
            }
            df
          })
        }

        # Standard imputation function (per-column for mean/median/knn)
        impute_df <- function(df) {
          numeric_cols <- sapply(df, is.numeric)
          feature_cols <- setdiff(names(df)[numeric_cols], meta_cols)

          for (col in feature_cols) {
            if (any(is.na(df[[col]]))) {
              if (method == "mean") {
                df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
              } else if (method == "median") {
                df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
              } else if (method == "knn") {
                # Simple KNN imputation using column means weighted by correlation
                k <- input$impute_knn_k %||% 5
                col_data <- df[[col]]
                missing_idx <- which(is.na(col_data))

                if (length(missing_idx) > 0 && sum(!is.na(col_data)) >= k) {
                  complete_idx <- which(!is.na(col_data))
                  other_numeric <- df[, feature_cols[feature_cols != col], drop = FALSE]
                  other_numeric <- other_numeric[, sapply(other_numeric, is.numeric), drop = FALSE]

                  for (mi in missing_idx) {
                    if (ncol(other_numeric) > 0) {
                      dists <- apply(other_numeric[complete_idx, , drop = FALSE], 1, function(row) {
                        sqrt(sum((row - other_numeric[mi, ])^2, na.rm = TRUE))
                      })
                      nearest <- complete_idx[order(dists)[1:min(k, length(complete_idx))]]
                    } else {
                      nearest <- complete_idx[1:min(k, length(complete_idx))]
                    }
                    df[[col]][mi] <- mean(col_data[nearest], na.rm = TRUE)
                  }
                } else {
                  # Fallback to mean
                  df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
                }
              }
            }
          }
          df
        }

        # Apply appropriate imputation method
        if (method == "pmm") {
          if (is.data.frame(data)) {
            imputed_data <- impute_pmm(data)
          } else if (is.list(data)) {
            imputed_data <- lapply(data, impute_pmm)
          } else {
            showNotification("Unsupported data format", type = "error")
            return()
          }
        } else {
          if (is.data.frame(data)) {
            imputed_data <- impute_df(data)
          } else if (is.list(data)) {
            imputed_data <- lapply(data, impute_df)
          } else {
            showNotification("Unsupported data format", type = "error")
            return()
          }
        }

        # Update project data
        proj$update_data(imputed_data, paste0("imputation_", method))
      })

      # Trigger UI refresh
      qc_trigger(qc_trigger() + 1)
      showNotification(sprintf("%s imputation applied!", tools::toTitleCase(method)),
                       type = "message")
    })

    # Reset to original data
    observeEvent(input$reset_imputation, {
      proj <- project()
      req(proj)

      proj$reset_data()
      qc_trigger(qc_trigger() + 1)
      showNotification("Data reset to original values", type = "message")
    })

    # Continue to Pipeline
    observeEvent(input$next_step, {
      proj <- project()
      if (is.null(proj) || !proj$can_proceed("pipeline")) {
        showNotification("Please resolve blocking QC issues first", type = "error")
        return()
      }
      switch_tab("pipeline")
    })
  })
}
