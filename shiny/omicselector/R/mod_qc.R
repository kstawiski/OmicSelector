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
    )
  )
}

# Module Server
mod_qc_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive for DE results
    de_results <- reactiveVal(NULL)

    # Run QC checks
    observeEvent(input$run_qc, {
      proj <- project()
      req(proj)

      withProgress(message = "Running QC checks...", value = 0.5, {
        proj$run_qc()
      })
    })

    # Traffic light display
    output$qc_traffic_lights <- renderUI({
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

    # Helper to get feature matrix
    get_feature_matrix <- function() {
      proj <- project()
      if (is.null(proj)) return(NULL)

      data <- proj$get_data()
      mapping <- proj$get_mapping()

      if (is.null(data) || is.null(mapping$target)) return(NULL)

      if (is.data.frame(data)) {
        meta_cols <- c(mapping$target, mapping$patient_id, mapping$batch, "mix")
        meta_cols <- meta_cols[!is.null(meta_cols) & nzchar(meta_cols)]
        keep <- setdiff(names(data), meta_cols)
        x <- data[, keep, drop = FALSE]

        # Filter by prefix
        if (!is.null(mapping$feature_prefix) && nzchar(mapping$feature_prefix)) {
          keep_cols <- grep(paste0("^", mapping$feature_prefix), names(x), value = TRUE)
          if (length(keep_cols) > 0) {
            x <- x[, keep_cols, drop = FALSE]
          }
        }

        # Keep only numeric
        numeric_cols <- sapply(x, is.numeric)
        x[, numeric_cols, drop = FALSE]
      } else {
        NULL
      }
    }

    # Helper to get target vector
    get_target_vector <- function() {
      proj <- project()
      if (is.null(proj)) return(NULL)

      data <- proj$get_data()
      mapping <- proj$get_mapping()

      if (is.null(data) || is.null(mapping$target)) return(NULL)

      if (is.data.frame(data) && mapping$target %in% names(data)) {
        data[[mapping$target]]
      } else {
        NULL
      }
    }

    # Column choices for plots
    columns <- reactive({
      proj <- project()
      if (is.null(proj)) return(character(0))
      data <- proj$get_data()
      if (is.null(data) || !is.data.frame(data)) return(character(0))
      names(data)
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

      # Handle missing values
      x <- x[complete.cases(x), , drop = FALSE]
      if (nrow(x) < 2) return(NULL)

      pca <- prcomp(x, scale. = TRUE)

      proj <- project()
      data <- proj$get_data()

      groups <- NULL
      if (!is.null(input$pca_color) && nzchar(input$pca_color) &&
          input$pca_color %in% names(data)) {
        groups <- data[[input$pca_color]][complete.cases(get_feature_matrix())]
      }

      list(scores = pca$x, groups = groups, var_exp = summary(pca)$importance[2, 1:2])
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

    # PCA 3D plot
    output$pca_3d <- plotly::renderPlotly({
      x <- get_feature_matrix()
      if (is.null(x) || ncol(x) < 3 || nrow(x) < 3) return(NULL)

      x <- x[complete.cases(x), , drop = FALSE]
      if (nrow(x) < 3) return(NULL)

      pca <- prcomp(x, scale. = TRUE)
      scores <- as.data.frame(pca$x[, 1:3])

      proj <- project()
      data <- proj$get_data()

      if (!is.null(input$pca_color) && nzchar(input$pca_color) &&
          input$pca_color %in% names(data)) {
        scores$group <- as.factor(data[[input$pca_color]][complete.cases(get_feature_matrix())])
      } else {
        scores$group <- "All"
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

    # Correlation plot
    observeEvent(input$plot_cor, {
      output$correlation_plot <- renderPlot({
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

        results <- data.frame(
          feature = names(x),
          stringsAsFactors = FALSE
        )

        results$log2FC <- sapply(seq_len(ncol(x)), function(i) {
          mean(x[case_idx, i], na.rm = TRUE) - mean(x[ctrl_idx, i], na.rm = TRUE)
        })

        results$pvalue <- sapply(seq_len(ncol(x)), function(i) {
          tryCatch(
            t.test(x[case_idx, i], x[ctrl_idx, i])$p.value,
            error = function(e) NA_real_
          )
        })

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

      plot(de$log2FC, de$neg_log10_p,
           pch = 19, col = "#64748b80",
           xlab = "log2 Fold Change",
           ylab = "-log10(p-value)",
           main = "Volcano Plot")

      # Highlight significant
      sig <- de$padj < 0.05 & abs(de$log2FC) > 1
      points(de$log2FC[sig], de$neg_log10_p[sig], pch = 19, col = "#ef4444")

      # Label top N
      top_n <- input$de_top_n
      if (top_n > 0 && nrow(de) >= top_n) {
        top_de <- head(de, top_n)
        text(top_de$log2FC, top_de$neg_log10_p,
             labels = top_de$feature, pos = 3, cex = 0.7)
      }

      abline(h = -log10(0.05), col = "gray", lty = 2)
      abline(v = c(-1, 1), col = "gray", lty = 2)
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
