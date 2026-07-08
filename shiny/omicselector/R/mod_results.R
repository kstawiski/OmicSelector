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
        # Stability CI Analysis Section
        div(
          class = "stability-ci-section mb-4",
          h6(bsicons::bs_icon("graph-up-arrow"), " Stability Confidence Interval"),
          uiOutput(ns("stability_ci_details")),
          div(
            class = "accordion mt-2",
            id = ns("stability_ci_help_accordion"),
            div(
              class = "accordion-item",
              h2(class = "accordion-header",
                 tags$button(
                   class = "accordion-button collapsed py-2 small",
                   type = "button",
                   `data-bs-toggle` = "collapse",
                   `data-bs-target` = paste0("#", ns("stability_ci_help_body")),
                   "How to interpret stability CIs?"
                 )
              ),
              div(
                id = ns("stability_ci_help_body"),
                class = "accordion-collapse collapse",
                `data-bs-parent` = paste0("#", ns("stability_ci_help_accordion")),
                div(
                  class = "accordion-body small",
                  tags$dl(
                    tags$dt("Narrow CI (e.g., [0.72, 0.78])"),
                    tags$dd("Stable feature selection across folds. Reliable for clinical use."),
                    tags$dt("Wide CI (e.g., [0.45, 0.85])"),
                    tags$dd("High variability. Consider collecting more data or simplifying the model."),
                    tags$dt("CI below 0.7"),
                    tags$dd("Feature selection is unstable. Results may not replicate."),
                    tags$dt("CI spans 0.7 threshold"),
                    tags$dd("Uncertain stability. External validation strongly recommended.")
                  )
                )
              )
            )
          )
        ),

        hr(),

        sliderInput(ns("min_freq"), "Minimum Selection Frequency",
                    value = 0.8, min = 0, max = 1, step = 0.1),
        h6("Consensus Features"),
        p(class = "text-muted small", "Features selected in at least the specified fraction of CV folds"),
        tableOutput(ns("consensus_features")),
        hr(),
        layout_column_wrap(
          width = 1/2,
          numericInput(ns("final_seed"), "Random Seed for Final Model", value = 42, min = 1),
          div()
        ),
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

    # Helper function: Calculate Nogueira Stability Index from binary feature selection matrix
    # Based on Nogueira et al. (2018) "Stability of Feature Selection Algorithms"
    # Formula: S = 1 - (mean variance across features) / (expected variance under random selection)
    calculate_nogueira_index <- function(feature_lists, all_features = NULL) {
      if (length(feature_lists) < 2) return(NA_real_)

      # Get all unique features
      if (is.null(all_features)) {
        all_features <- unique(unlist(feature_lists))
      }
      p <- length(all_features)  # Total number of features
      if (p == 0) return(NA_real_)

      n <- length(feature_lists)  # Number of folds/bootstrap samples

      # Create binary selection matrix (folds x features)
      Z <- matrix(0, nrow = n, ncol = p)
      colnames(Z) <- all_features
      for (i in seq_len(n)) {
        selected <- feature_lists[[i]]
        matched <- selected[selected %in% all_features]
        if (length(matched) > 0) {
          Z[i, matched] <- 1
        }
      }

      # Calculate selection frequencies per feature
      pf <- colMeans(Z)  # proportion of times each feature is selected

      # Mean number of selected features per fold
      k_bar <- mean(rowSums(Z))
      if (k_bar == 0 || k_bar == p) return(1)  # Edge case: all or none selected

      # Calculate Nogueira stability index
      # Variance for each feature: pf * (1 - pf)
      feat_var <- pf * (1 - pf)
      mean_var <- mean(feat_var)

      # Expected variance under random selection
      expected_var <- (k_bar / p) * (1 - k_bar / p)
      if (expected_var == 0) return(1)

      # Stability = 1 - (observed variance / expected variance)
      stability <- 1 - (mean_var / expected_var)

      # Clamp to [0, 1]
      max(0, min(1, stability))
    }

    # Helper function: Bootstrap confidence interval for Nogueira stability
    # Resamples fold results with replacement to estimate uncertainty
    bootstrap_stability_ci <- function(feature_lists, n_boot = 1000, conf_level = 0.95, seed = NULL) {
      if (length(feature_lists) < 2) {
        return(list(point = NA_real_, lower = NA_real_, upper = NA_real_, se = NA_real_))
      }

      if (!is.null(seed)) set.seed(seed)

      # Get all unique features for consistent calculation
      all_features <- unique(unlist(feature_lists))
      n_folds <- length(feature_lists)

      # Calculate point estimate
      point_estimate <- calculate_nogueira_index(feature_lists, all_features)

      # Bootstrap resampling of folds
      boot_estimates <- numeric(n_boot)
      for (b in seq_len(n_boot)) {
        # Resample fold indices with replacement
        boot_idx <- sample(seq_len(n_folds), n_folds, replace = TRUE)
        boot_lists <- feature_lists[boot_idx]
        boot_estimates[b] <- calculate_nogueira_index(boot_lists, all_features)
      }

      # Remove NAs
      boot_estimates <- boot_estimates[!is.na(boot_estimates)]
      if (length(boot_estimates) < 10) {
        return(list(point = point_estimate, lower = NA_real_, upper = NA_real_, se = NA_real_))
      }

      # Calculate confidence interval (percentile method)
      alpha <- 1 - conf_level
      lower <- quantile(boot_estimates, probs = alpha / 2, na.rm = TRUE)
      upper <- quantile(boot_estimates, probs = 1 - alpha / 2, na.rm = TRUE)
      se <- sd(boot_estimates, na.rm = TRUE)

      list(
        point = point_estimate,
        lower = as.numeric(lower),
        upper = as.numeric(upper),
        se = se,
        n_boot = length(boot_estimates)
      )
    }

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

      # Stability index from nogueira_index with bootstrap CI
      best_stab <- NA
      stab_ci <- list(point = NA_real_, lower = NA_real_, upper = NA_real_, se = NA_real_)

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

      # Calculate bootstrap CI from selected_features_per_fold
      avg_k <- NA
      if (!is.null(stab) && !is.null(stab$selected_features_per_fold)) {
        sf <- stab$selected_features_per_fold
        if (is.list(sf) && length(sf) > 0) {
          # Calculate average feature count
          all_lengths <- unlist(lapply(sf, function(x) sapply(x, length)))
          avg_k <- mean(all_lengths, na.rm = TRUE)

          # Flatten feature lists for bootstrap CI calculation
          # sf can be nested (learner -> folds) or flat (folds)
          feature_lists <- if (is.list(sf[[1]]) && !is.character(sf[[1]])) {
            # Nested: take first learner's folds
            sf[[1]]
          } else {
            # Flat list of character vectors
            sf
          }

          # Compute bootstrap CI (use 500 for speed in interactive context)
          stab_ci <- tryCatch({
            bootstrap_stability_ci(feature_lists, n_boot = 500, conf_level = 0.95, seed = 42)
          }, error = function(e) {
            list(point = best_stab, lower = NA_real_, upper = NA_real_, se = NA_real_)
          })

          # Use bootstrap point estimate if original is NA
          if (is.na(best_stab) && !is.na(stab_ci$point)) {
            best_stab <- stab_ci$point
          }
        }
      }

      # Determine if CI contains the 0.7 threshold (uncertain validation)
      ci_uncertain <- FALSE
      if (!is.na(stab_ci$lower) && !is.na(stab_ci$upper)) {
        ci_uncertain <- stab_ci$lower < 0.7 && stab_ci$upper >= 0.7
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
            div(class = "metric-label", "Stability Index"),
            # Show 95% CI if available
            if (!is.na(stab_ci$lower) && !is.na(stab_ci$upper)) {
              div(class = "metric-ci text-muted small",
                  sprintf("95%% CI: [%.3f, %.3f]", stab_ci$lower, stab_ci$upper))
            }
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
            class = if (!is.na(best_stab) && best_stab >= 0.7 && !ci_uncertain) {
              "metric-box bg-success text-white"
            } else if (ci_uncertain) {
              "metric-box bg-warning"
            } else {
              "metric-box bg-warning"
            },
            div(class = "metric-value",
                if (!is.na(best_stab) && best_stab >= 0.7 && !ci_uncertain) {
                  "PASS"
                } else if (ci_uncertain) {
                  "UNCERTAIN"
                } else {
                  "CHECK"
                }),
            div(class = "metric-label", "Validation Status"),
            if (ci_uncertain) {
              div(class = "small", "CI spans threshold")
            }
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

      # P0 Fix: Check if select_best_signature function exists before calling
      candidates <- NULL
      if (exists("select_best_signature", mode = "function")) {
        candidates <- tryCatch({
          select_best_signature(
            result,
            mode = "pareto",  # Get all candidates for visualization
            metric = "classif.auc"
          )
        }, error = function(e) NULL)
      }

      if (is.null(candidates) || nrow(candidates) == 0) {
        # Fallback: build candidates from performance table with proper values
        perf <- result$performance
        if (is.null(perf) || !is.data.frame(perf)) {
          showNotification("No performance data available for visualization", type = "warning")
          return(NULL)
        }

        # Extract AUC column
        auc_col <- intersect(c("classif.auc", "auc", "mean_auc"), names(perf))
        auc_vals <- if (length(auc_col) > 0) perf[[auc_col[1]]] else rep(NA_real_, nrow(perf))

        candidates <- data.frame(
          learner_id = if ("learner_id" %in% names(perf)) perf$learner_id else paste0("Model_", seq_len(nrow(perf))),
          mean_metric = auc_vals,
          stability = NA_real_,
          mean_k = NA_real_,
          stringsAsFactors = FALSE
        )

        # Add stability if available from result
        stab <- result$stability
        if (!is.null(stab)) {
          if (!is.null(stab$nogueira_index)) {
            ni <- stab$nogueira_index
            if (is.data.frame(ni) && "learner_id" %in% names(ni) && "nogueira_index" %in% names(ni)) {
              candidates <- merge(candidates, ni[, c("learner_id", "nogueira_index"), drop = FALSE],
                                  by = "learner_id", all.x = TRUE)
              candidates$stability <- candidates$nogueira_index
              candidates$nogueira_index <- NULL
            } else if (is.numeric(ni) && length(ni) == 1) {
              candidates$stability <- ni
            }
          }
          # Extract mean_k from selected features if available
          if (!is.null(stab$selected_features_per_fold)) {
            sf <- stab$selected_features_per_fold
            if (is.list(sf) && length(sf) > 0) {
              avg_k <- mean(unlist(lapply(sf, function(x) sapply(x, length))), na.rm = TRUE)
              candidates$mean_k <- avg_k
            }
          }
        }

        # Notify user about fallback
        showNotification("Using fallback visualization (some metrics may be unavailable)", type = "message", duration = 3)
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

      # P0 Fix: Check if select_best_signature function exists before calling
      if (!exists("select_best_signature", mode = "function")) {
        showNotification(
          "Signature selection function (select_best_signature) is not available. Please ensure OmicSelector package is properly installed.",
          type = "error",
          duration = 10
        )
        return()
      }

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

    # Update selected signature features when frequency slider changes
    observeEvent(input$min_freq, {
      sig <- selected_sig()
      if (is.null(sig)) return()  # Only update if a signature is already selected

      proj <- project()
      if (is.null(proj)) return()

      result <- proj$get_benchmark_result()
      if (is.null(result)) return()

      # Recalculate consensus features with new threshold
      consensus <- tryCatch({
        get_consensus_features(result, sig$learner_id, min_frequency = input$min_freq)
      }, error = function(e) {
        data.frame(feature = character(), frequency = numeric(), selected = logical())
      })

      selected_features <- if (nrow(consensus) > 0 && "selected" %in% names(consensus)) {
        consensus$feature[consensus$selected == TRUE]
      } else if (nrow(consensus) > 0) {
        consensus$feature[consensus$frequency >= input$min_freq]
      } else {
        character()
      }

      # Update the selected signature with new features
      updated_sig <- sig
      updated_sig$features <- selected_features
      updated_sig$consensus_df <- consensus
      updated_sig$n_features <- length(selected_features)

      selected_sig(updated_sig)
    }, ignoreInit = TRUE)

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

    # Stability CI details for Stability Explorer
    output$stability_ci_details <- renderUI({
      proj <- project()
      if (is.null(proj)) return(NULL)

      result <- proj$get_benchmark_result()
      if (is.null(result)) {
        return(div(class = "text-muted", "Run benchmark to see stability analysis"))
      }

      stab <- result$stability
      if (is.null(stab) || is.null(stab$selected_features_per_fold)) {
        return(div(class = "text-muted", "No feature selection data available"))
      }

      # Extract feature lists and compute bootstrap CI
      sf <- stab$selected_features_per_fold
      feature_lists <- if (is.list(sf[[1]]) && !is.character(sf[[1]])) {
        sf[[1]]  # Nested: first learner's folds
      } else {
        sf  # Flat list
      }

      # Calculate bootstrap CI with 1000 iterations for detailed view
      stab_ci <- tryCatch({
        bootstrap_stability_ci(feature_lists, n_boot = 1000, conf_level = 0.95, seed = 42)
      }, error = function(e) {
        list(point = NA_real_, lower = NA_real_, upper = NA_real_, se = NA_real_, n_boot = 0)
      })

      n_folds <- length(feature_lists)

      if (is.na(stab_ci$point)) {
        return(div(class = "text-muted", "Could not compute stability index"))
      }

      # Determine interpretation
      ci_width <- if (!is.na(stab_ci$upper) && !is.na(stab_ci$lower)) {
        stab_ci$upper - stab_ci$lower
      } else NA

      interpretation <- if (is.na(stab_ci$lower) || is.na(stab_ci$upper)) {
        list(badge = "secondary", text = "CI unavailable")
      } else if (stab_ci$lower >= 0.7) {
        list(badge = "success", text = "Reliably stable")
      } else if (stab_ci$upper < 0.7) {
        list(badge = "danger", text = "Consistently unstable")
      } else {
        list(badge = "warning", text = "Uncertain - spans threshold")
      }

      ci_precision <- if (is.na(ci_width)) {
        ""
      } else if (ci_width < 0.1) {
        " (narrow - high precision)"
      } else if (ci_width < 0.2) {
        " (moderate precision)"
      } else {
        " (wide - consider more data)"
      }

      div(
        class = "stability-ci-result",
        div(
          class = "d-flex align-items-center mb-2",
          span(class = "h4 mb-0 me-2", sprintf("%.3f", stab_ci$point)),
          span(class = sprintf("badge bg-%s", interpretation$badge), interpretation$text)
        ),
        if (!is.na(stab_ci$lower) && !is.na(stab_ci$upper)) {
          tagList(
            p(class = "mb-1",
              strong("95% Confidence Interval: "),
              sprintf("[%.3f, %.3f]", stab_ci$lower, stab_ci$upper),
              span(class = "text-muted small", ci_precision)
            ),
            if (!is.na(stab_ci$se)) {
              p(class = "text-muted small mb-1",
                sprintf("Standard Error: %.4f | Bootstrap samples: %d | CV folds: %d",
                        stab_ci$se, stab_ci$n_boot %||% 1000, n_folds))
            }
          )
        },
        # Visual CI bar
        if (!is.na(stab_ci$lower) && !is.na(stab_ci$upper)) {
          div(
            class = "ci-visual mt-2",
            style = "position: relative; height: 24px; background: linear-gradient(to right, #dc3545 0%, #dc3545 70%, #198754 70%, #198754 100%); border-radius: 4px;",
            # 0.7 threshold marker
            div(
              style = "position: absolute; left: 70%; top: 0; bottom: 0; width: 2px; background: white;",
              title = "0.7 threshold"
            ),
            # CI range bar
            div(
              style = sprintf(
                "position: absolute; left: %.1f%%; width: %.1f%%; top: 4px; height: 16px; background: rgba(255,255,255,0.9); border: 2px solid #333; border-radius: 3px;",
                max(0, stab_ci$lower * 100),
                min(100, (stab_ci$upper - stab_ci$lower) * 100)
              ),
              title = sprintf("95%% CI: [%.3f, %.3f]", stab_ci$lower, stab_ci$upper)
            ),
            # Point estimate marker
            div(
              style = sprintf(
                "position: absolute; left: calc(%.1f%% - 4px); top: 2px; width: 8px; height: 20px; background: #333; border-radius: 2px;",
                stab_ci$point * 100
              ),
              title = sprintf("Point estimate: %.3f", stab_ci$point)
            ),
            # Scale labels
            div(
              class = "d-flex justify-content-between mt-1 small text-muted",
              span("0"), span("0.7"), span("1.0")
            )
          )
        }
      )
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
          incProgress(0.1, detail = "Preparing learner...")

          # Get consensus features from the selected signature
          consensus_features <- sig$features
          if (is.null(consensus_features) || length(consensus_features) == 0) {
            showNotification("No consensus features available. Select a signature first.", type = "error")
            return()
          }

          message(sprintf("[mod_results] Training final model with %d consensus features", length(consensus_features)))

          # Find the specific learner matching the selected signature
          selected_learner <- learner
          if (is.list(learner) && length(learner) > 0) {
            learner_ids <- sapply(learner, function(l) {
              if (!is.null(l$id)) l$id else ""
            })
            match_idx <- which(learner_ids == sig$learner_id)
            if (length(match_idx) > 0) {
              selected_learner <- learner[[match_idx[1]]]
            } else {
              partial_match <- grep(sig$learner_id, learner_ids, fixed = TRUE)
              if (length(partial_match) > 0) {
                selected_learner <- learner[[partial_match[1]]]
              } else {
                selected_learner <- learner[[1]]
                showNotification(sprintf("Learner '%s' not found, using first learner", sig$learner_id), type = "warning")
              }
            }
          }

          incProgress(0.2, detail = "Creating fixed-feature learner...")

          # Create a learner that uses ONLY the consensus features (no re-selection)
          # This prevents the learner from selecting different features than what was validated
          seed_val <- input$final_seed %||% 42
          set.seed(seed_val)

          # Try to use pipeline's fit_with_features if available, otherwise fit normally
          fit <- tryCatch({
            pipeline$fit_with_features(
              learner = selected_learner,
              features = consensus_features,
              seed = seed_val
            )
          }, error = function(e) {
            # Fallback: use standard fit but warn about potential feature mismatch
            message("[mod_results] fit_with_features not available, using standard fit")
            showNotification(
              "Note: Using standard fit. Features may differ from consensus set.",
              type = "warning",
              duration = 5
            )
            pipeline$fit(
              learner = selected_learner,
              seed = seed_val
            )
          })

          incProgress(0.5, detail = "Storing results...")

          # Store the fit result along with the consensus features
          fit_bundle <- list(
            fit = fit,
            learner = selected_learner,
            consensus_features = consensus_features,
            learner_id = sig$learner_id,
            stability = sig$stability,
            auc = sig$auc,
            seed = seed_val
          )
          proj$set_final_fit(fit_bundle)

          showNotification(
            sprintf("Final model trained with %d consensus features!", length(consensus_features)),
            type = "message"
          )
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
