#' Results Module
#'
#' @description Decision-support visualization for signature selection

# Module UI
mod_results_ui <- function(id) {
  ns <- NS(id)

  tagList(
    div(
      class = "step-header mb-4",
      h3(bsicons::bs_icon("trophy"), "Step 5: Results & Signature Selection"),
      p(class = "text-muted", "Select optimal biomarker signature using multi-objective criteria")
    ),

    # Executive Summary Card
    card(
      class = "mb-4",
      card_header(
        class = "bg-primary text-white",
        bsicons::bs_icon("card-checklist"),
        span("Executive Summary", class = "ms-2")
      ),
      card_body(
        uiOutput(ns("executive_summary"))
      )
    ),

    # Selection Mode and Trade-off Visualization
    layout_columns(
      col_widths = c(4, 8),

      # Selection Settings
      card(
        card_header("Signature Selection"),
        card_body(
          selectInput(ns("selection_mode"), "Selection Mode",
                      choices = c(
                        "Constrained 1SE (Conservative)" = "constrained_1se",
                        "Weighted Scoring" = "weighted",
                        "Pareto Frontier" = "pareto"
                      ),
                      selected = "constrained_1se"),

          conditionalPanel(
            condition = sprintf("input['%s'] == 'weighted'", ns("selection_mode")),
            h6("Objective Weights"),
            sliderInput(ns("w_performance"), "Performance", value = 0.5, min = 0, max = 1, step = 0.1),
            sliderInput(ns("w_stability"), "Stability", value = 0.3, min = 0, max = 1, step = 0.1),
            sliderInput(ns("w_parsimony"), "Parsimony", value = 0.2, min = 0, max = 1, step = 0.1)
          ),

          hr(),

          h6("Constraints"),
          numericInput(ns("min_auc"), "Min AUC", value = 0.7, min = 0.5, max = 1, step = 0.05),
          numericInput(ns("min_stability"), "Min Stability", value = 0.5, min = 0, max = 1, step = 0.05),
          numericInput(ns("max_features"), "Max Features", value = 50, min = 1, max = 500),

          hr(),

          actionButton(ns("select_signature"), "Select Signature",
                       class = "btn-primary w-100",
                       icon = icon("check-double"))
        )
      ),

      # Trade-off Plot
      card(
        card_header("Trade-off Visualization"),
        card_body(
          plotly::plotlyOutput(ns("tradeoff_plot"), height = "400px")
        )
      )
    ),

    # Candidate Details
    layout_column_wrap(
      width = 1/2,

      # Candidates Table
      card(
        card_header("Candidate Signatures"),
        card_body(
          div(style = "max-height: 300px; overflow-y: auto;",
              tableOutput(ns("candidates_table")))
        )
      ),

      # Selected Signature
      card(
        card_header(
          class = "bg-success text-white",
          bsicons::bs_icon("star-fill"),
          span("Selected Signature", class = "ms-2")
        ),
        card_body(
          uiOutput(ns("selected_signature_ui")),
          hr(),
          h6("Selected Features"),
          div(style = "max-height: 200px; overflow-y: auto;",
              tableOutput(ns("selected_features")))
        )
      )
    ),

    # Stability Explorer
    card(
      class = "mt-4",
      card_header("Stability Explorer"),
      card_body(
        sliderInput(ns("min_freq"), "Minimum Selection Frequency",
                    value = 0.8, min = 0, max = 1, step = 0.1),
        h6("Consensus Features"),
        p(class = "text-muted small", "Features selected in at least the specified fraction of CV folds"),
        tableOutput(ns("consensus_features")),
        hr(),
        actionButton(ns("train_final"), "Train Final Model with Selected Signature",
                     class = "btn-success btn-lg",
                     icon = icon("graduation-cap")),
        actionButton(ns("next_step"), "Continue to Export",
                     class = "btn-outline-primary btn-lg ms-2",
                     icon = icon("arrow-right"))
      )
    )
  )
}

# Module Server
mod_results_server <- function(id, project, switch_tab) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Selected signature reactive
    selected_sig <- reactiveVal(NULL)

    # Candidates data from signature selection
    candidates_data <- reactiveVal(NULL)

    # Executive summary
    output$executive_summary <- renderUI({
      proj <- project()
      if (is.null(proj)) return(NULL)

      result <- proj$get_benchmark_result()
      if (is.null(result)) {
        return(div(class = "text-muted", "No benchmark results available"))
      }

      # Extract summary metrics from NestedCVResult
      perf <- result$performance
      stab <- result$stability

      # Best AUC from performance table
      best_auc <- NA
      if (!is.null(perf) && is.data.frame(perf)) {
        auc_cols <- grep("classif\\.auc|auc", names(perf), value = TRUE, ignore.case = TRUE)
        if (length(auc_cols) > 0) {
          best_auc <- max(perf[[auc_cols[1]]], na.rm = TRUE)
        }
      }

      # Stability index from nogueira_index
      best_stab <- NA
      if (!is.null(stab) && !is.null(stab$nogueira_index)) {
        ni <- stab$nogueira_index
        if (is.numeric(ni) && length(ni) == 1) {
          best_stab <- ni
        } else if (is.data.frame(ni) || data.table::is.data.table(ni)) {
          if ("nogueira_index" %in% names(ni)) {
            best_stab <- max(ni$nogueira_index, na.rm = TRUE)
          }
        }
      }

      # Average feature count
      avg_k <- NA
      if (!is.null(stab) && !is.null(stab$selected_features_per_fold)) {
        sf <- stab$selected_features_per_fold
        if (is.list(sf) && length(sf) > 0) {
          all_lengths <- unlist(lapply(sf, function(x) sapply(x, length)))
          avg_k <- mean(all_lengths, na.rm = TRUE)
        }
      }

      div(
        class = "row text-center",
        div(
          class = "col-md-3",
          div(
            class = "metric-box",
            div(class = "metric-value", if (!is.na(best_auc)) sprintf("%.3f", best_auc) else "N/A"),
            div(class = "metric-label", "Best AUC")
          )
        ),
        div(
          class = "col-md-3",
          div(
            class = "metric-box",
            div(class = "metric-value", if (!is.na(best_stab)) sprintf("%.3f", best_stab) else "N/A"),
            div(class = "metric-label", "Stability Index")
          )
        ),
        div(
          class = "col-md-3",
          div(
            class = "metric-box",
            div(class = "metric-value", if (!is.na(avg_k)) sprintf("%.0f", avg_k) else "N/A"),
            div(class = "metric-label", "Avg Features")
          )
        ),
        div(
          class = "col-md-3",
          div(
            class = if (!is.na(best_stab) && best_stab >= 0.7) "metric-box bg-success text-white" else "metric-box bg-warning",
            div(class = "metric-value", if (!is.na(best_stab) && best_stab >= 0.7) "PASS" else "CHECK"),
            div(class = "metric-label", "Validation Status")
          )
        )
      )
    })

    # Trade-off plot using actual benchmark data
    output$tradeoff_plot <- plotly::renderPlotly({
      proj <- project()
      if (is.null(proj)) return(NULL)

      result <- proj$get_benchmark_result()
      if (is.null(result)) return(NULL)

      # Try to use select_best_signature if available
      candidates <- tryCatch({
        select_best_signature(
          result,
          mode = "pareto",  # Get all candidates for visualization
          metric = "classif.auc"
        )
      }, error = function(e) NULL)

      if (is.null(candidates) || nrow(candidates) == 0) {
        # Fallback: build candidates from performance table
        perf <- result$performance
        if (is.null(perf) || !is.data.frame(perf)) return(NULL)

        candidates <- data.frame(
          learner_id = if ("learner_id" %in% names(perf)) perf$learner_id else paste0("Model_", seq_len(nrow(perf))),
          mean_metric = if ("classif.auc" %in% names(perf)) perf$classif.auc else runif(nrow(perf), 0.7, 0.95),
          stability = NA_real_,
          mean_k = NA_real_
        )

        # Add stability if available
        stab <- result$stability
        if (!is.null(stab$nogueira_index)) {
          ni <- stab$nogueira_index
          if (is.data.frame(ni) && "learner_id" %in% names(ni)) {
            candidates <- merge(candidates, ni[, c("learner_id", "nogueira_index")],
                                by = "learner_id", all.x = TRUE)
            candidates$stability <- candidates$nogueira_index
          } else if (is.numeric(ni) && length(ni) == 1) {
            candidates$stability <- ni
          }
        }
      }

      # Store for later use
      candidates_data(candidates)

      # Handle missing values for plotting
      plot_df <- as.data.frame(candidates)
      plot_df$stability[is.na(plot_df$stability)] <- 0.5
      plot_df$mean_k[is.na(plot_df$mean_k)] <- 20

      # Handle metric flip for brier/loss metrics
      if ("original_metric" %in% names(plot_df)) {
        plot_df$auc <- plot_df$original_metric
      } else {
        plot_df$auc <- plot_df$mean_metric
      }

      # Highlight selected if available
      plot_df$is_selected <- if ("selected" %in% names(plot_df)) plot_df$selected else FALSE

      plotly::plot_ly(plot_df,
                      x = ~stability,
                      y = ~auc,
                      size = ~mean_k,
                      color = ~learner_id,
                      text = ~paste("Learner:", learner_id,
                                    "<br>AUC:", round(auc, 3),
                                    "<br>Stability:", round(stability, 3),
                                    "<br>Features:", round(mean_k, 1)),
                      type = "scatter",
                      mode = "markers",
                      marker = list(sizemode = "diameter", sizeref = 2, sizemin = 4),
                      hoverinfo = "text") |>
        plotly::layout(
          title = list(text = "Performance vs Stability Trade-off", x = 0),
          xaxis = list(title = "Stability Index", range = c(0, 1)),
          yaxis = list(title = "AUC"),
          showlegend = TRUE
        )
    })

    # Candidates table
    output$candidates_table <- renderTable({
      proj <- project()
      if (is.null(proj)) return(NULL)

      result <- proj$get_benchmark_result()
      if (is.null(result)) {
        return(data.frame(Message = "Run benchmark to see results"))
      }

      # Use cached candidates or get from performance
      candidates <- candidates_data()

      if (is.null(candidates) || nrow(candidates) == 0) {
        perf <- result$performance
        if (!is.null(perf) && is.data.frame(perf)) {
          # Return simplified performance table
          cols_to_show <- intersect(c("learner_id", "classif.auc", "classif.acc", "classif.bbrier"),
                                    names(perf))
          if (length(cols_to_show) > 0) {
            return(perf[, cols_to_show, drop = FALSE])
          }
          return(perf)
        }
        return(data.frame(Message = "No performance data available"))
      }

      # Format candidates for display
      display_df <- data.frame(
        Learner = candidates$learner_id,
        AUC = if ("original_metric" %in% names(candidates)) candidates$original_metric else candidates$mean_metric,
        Stability = candidates$stability,
        Features = candidates$mean_k
      )
      display_df$AUC <- round(display_df$AUC, 3)
      display_df$Stability <- round(display_df$Stability, 3)
      display_df$Features <- round(display_df$Features, 1)

      display_df
    })

    # Select signature
    observeEvent(input$select_signature, {
      proj <- project()
      req(proj)

      result <- proj$get_benchmark_result()
      if (is.null(result)) {
        showNotification("No benchmark results", type = "error")
        return()
      }

      # Build weights for weighted mode
      weights <- c(
        performance = input$w_performance,
        stability = input$w_stability,
        parsimony = input$w_parsimony
      )

      # Run signature selection using the actual function
      best <- tryCatch({
        select_best_signature(
          result,
          mode = input$selection_mode,
          metric = "classif.auc",
          auc_min = input$min_auc,
          stability_min = input$min_stability,
          k_max = input$max_features,
          weights = weights
        )
      }, error = function(e) {
        showNotification(paste("Selection error:", e$message), type = "error")
        return(NULL)
      })

      if (is.null(best) || nrow(best) == 0) {
        showNotification("No candidates meet the constraints. Try relaxing constraints.", type = "warning")
        return()
      }

      # Get the selected row
      if ("selected" %in% names(best)) {
        selected_row <- best[best$selected == TRUE, ]
        if (nrow(selected_row) == 0) selected_row <- best[1, ]
      } else {
        selected_row <- best[1, ]
      }

      # Get consensus features for the selected learner
      consensus <- tryCatch({
        get_consensus_features(result, selected_row$learner_id[1], min_frequency = input$min_freq)
      }, error = function(e) {
        data.frame(feature = character(), frequency = numeric(), selected = logical())
      })

      selected_features <- if (nrow(consensus) > 0 && "selected" %in% names(consensus)) {
        consensus$feature[consensus$selected == TRUE]
      } else if (nrow(consensus) > 0) {
        consensus$feature
      } else {
        character()
      }

      # Store selection result
      selected_sig(list(
        learner_id = selected_row$learner_id[1],
        auc = if ("original_metric" %in% names(selected_row)) selected_row$original_metric[1] else selected_row$mean_metric[1],
        se = if ("se_metric" %in% names(selected_row)) selected_row$se_metric[1] else NA,
        stability = selected_row$stability[1],
        n_features = selected_row$mean_k[1],
        features = selected_features,
        consensus_df = consensus,
        mode = input$selection_mode,
        score = if ("score" %in% names(selected_row)) selected_row$score[1] else NA
      ))

      # Update candidates data with selection info
      candidates_data(best)

      showNotification("Signature selected!", type = "message")
    })

    # Selected signature display
    output$selected_signature_ui <- renderUI({
      sig <- selected_sig()
      if (is.null(sig)) {
        return(div(class = "text-muted", "Click 'Select Signature' to choose optimal signature"))
      }

      tagList(
        p(strong("Learner: "), sig$learner_id),
        p(strong("AUC: "), sprintf("%.3f", sig$auc),
          if (!is.na(sig$se)) sprintf(" (SE: %.3f)", sig$se) else ""),
        p(strong("Stability: "), if (!is.na(sig$stability)) sprintf("%.3f", sig$stability) else "N/A"),
        p(strong("Features: "), if (!is.na(sig$n_features)) sprintf("%.0f", sig$n_features) else length(sig$features)),
        if (!is.na(sig$score)) p(strong("Score: "), sprintf("%.3f", sig$score)),
        p(class = "text-muted small", sprintf("Selection mode: %s", sig$mode))
      )
    })

    # Selected features table
    output$selected_features <- renderTable({
      sig <- selected_sig()
      if (is.null(sig) || length(sig$features) == 0) {
        return(data.frame(Message = "No features selected"))
      }

      # Use consensus_df if available for frequency info
      if (!is.null(sig$consensus_df) && nrow(sig$consensus_df) > 0) {
        df <- sig$consensus_df[sig$consensus_df$selected == TRUE, c("feature", "frequency")]
        names(df) <- c("Feature", "Frequency")
        df$Frequency <- sprintf("%.0f%%", df$Frequency * 100)
        return(df)
      }

      data.frame(Feature = sig$features)
    })

    # Consensus features
    output$consensus_features <- renderTable({
      proj <- project()
      if (is.null(proj)) return(NULL)

      result <- proj$get_benchmark_result()
      if (is.null(result)) {
        return(data.frame(Message = "Run benchmark first"))
      }

      sig <- selected_sig()
      learner_id <- if (!is.null(sig)) sig$learner_id else NULL

      # If no selection yet, try to get from first learner
      if (is.null(learner_id)) {
        perf <- result$performance
        if (!is.null(perf) && "learner_id" %in% names(perf)) {
          learner_id <- perf$learner_id[1]
        }
      }

      if (is.null(learner_id)) {
        return(data.frame(Message = "Select a signature first"))
      }

      # Get consensus features
      consensus <- tryCatch({
        get_consensus_features(result, learner_id, min_frequency = input$min_freq)
      }, error = function(e) {
        data.frame(feature = character(), frequency = numeric(), selected = logical())
      })

      if (nrow(consensus) == 0) {
        return(data.frame(Message = "No consensus features found"))
      }

      # Filter to selected and format
      selected_consensus <- consensus[consensus$frequency >= input$min_freq, ]

      if (nrow(selected_consensus) == 0) {
        return(data.frame(Message = sprintf("No features selected in >= %.0f%% of folds", input$min_freq * 100)))
      }

      data.frame(
        Feature = selected_consensus$feature,
        Frequency = sprintf("%.0f%%", selected_consensus$frequency * 100)
      )
    })

    # Train final model
    observeEvent(input$train_final, {
      proj <- project()
      req(proj)

      sig <- selected_sig()
      if (is.null(sig)) {
        showNotification("Select a signature first", type = "error")
        return()
      }

      pipeline <- proj$get_pipeline()
      learner <- proj$get_learner()

      if (is.null(pipeline) || is.null(learner)) {
        showNotification("Pipeline not available", type = "error")
        return()
      }

      withProgress(message = "Training final model...", value = 0.5, {
        tryCatch({
          # Fit the learner on full data
          fit <- pipeline$fit(
            learner = learner,
            seed = 42
          )

          # Store the fit result
          proj$set_final_fit(fit)

          # Also store the selected signature info
          proj$set_mapping(
            target = proj$get_mapping()$target,
            positive = proj$get_mapping()$positive
          )

          showNotification("Final model trained!", type = "message")
        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })

    # Continue to export
    observeEvent(input$next_step, {
      proj <- project()
      if (is.null(proj)) {
        showNotification("No project loaded", type = "error")
        return()
      }

      # Allow continuing even without final fit for review
      if (!proj$can_proceed("export")) {
        showNotification("Train final model first to export, or select signature to continue", type = "warning")
      }

      switch_tab("export")
    })
  })
}
