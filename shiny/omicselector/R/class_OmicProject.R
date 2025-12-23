#' @title OmicProject: R6 State Management Class for WebUI
#'
#' @description
#' Centralized state manager for the OmicSelector WebUI. Handles data storage,
#' validation, persistence, and workflow gating. Separates business logic from
#' UI reactivity to prevent "reactive spaghetti".
#'
#' @details
#' This class is designed to work with Shiny modules and allows:
#' - Project persistence (save/load to RDS)
#' - Workflow gating (can_proceed() checks)
#' - Zero-leakage validation
#' - Audit trail logging
#'
#' @importFrom R6 R6Class
#' @export
OmicProject <- R6::R6Class(
  "OmicProject",

  public = list(
    #' @field id Project identifier
    id = NULL,

    #' @field name Human-readable project name
    name = NULL,

    #' @field path Project directory path
    path = NULL,

    #' @field created_at Creation timestamp
    created_at = NULL,

    #' @field modified_at Last modification timestamp
    modified_at = NULL,

    #' @description Create a new OmicProject
    #' @param id Project ID (auto-generated if NULL)
    #' @param name Project name
    #' @param path Project directory path
    initialize = function(id = NULL, name = "Untitled Project", path = NULL) {
      self$id <- id %||% private$.generate_id()
      self$name <- name
      self$created_at <- Sys.time()
      self$modified_at <- Sys.time()

      if (is.null(path)) {
        root <- Sys.getenv("OMICSELECTOR_ANALYSES_DIR", "/OmicSelector/analyses")
        if (!dir.exists(root)) {
          root <- file.path(getwd(), "analyses")
        }
        self$path <- file.path(root, self$id)
      } else {
        self$path <- path
      }

      if (!dir.exists(self$path)) {
        dir.create(self$path, recursive = TRUE, showWarnings = FALSE)
      }

      private$.log_event("Project created")
    },

    # =========================================================================
    # DATA MANAGEMENT
    # =========================================================================

    #' @description Set raw data
    #' @param data Data frame or list of data frames (multi-omics)
    #' @param is_multi Logical, is this multi-omics data?
    set_data = function(data, is_multi = FALSE) {
      private$.data_raw <- data
      private$.data_current <- data
      private$.is_multi_omics <- is_multi
      private$.qc_status <- NULL
      private$.log_event("Data loaded")
      self$modified_at <- Sys.time()
      self$save_checkpoint("data_loaded")
      invisible(self)
    },

    #' @description Get current working data
    get_data = function() {
      private$.data_current
    },

    #' @description Get raw (original) data
    get_raw_data = function() {
      private$.data_raw
    },

    #' @description Update current data (after preprocessing)
    #' @param data Preprocessed data
    #' @param step_name Name of preprocessing step
    update_data = function(data, step_name = "preprocessing") {
      private$.data_current <- data
      private$.log_event(paste0("Data updated: ", step_name))
      self$modified_at <- Sys.time()
      invisible(self)
    },

    #' @description Reset data to raw
    reset_data = function() {
      private$.data_current <- private$.data_raw
      private$.log_event("Data reset to raw")
      self$modified_at <- Sys.time()
      invisible(self)
    },

    # =========================================================================
    # COLUMN MAPPING
    # =========================================================================

    #' @description Set column mappings
    #' @param target Target column name
    #' @param positive Positive class label
    #' @param patient_id Patient ID column (for grouped CV)
    #' @param batch Batch column (for batch correction)
    #' @param feature_prefix Feature prefix filter
    set_mapping = function(target = NULL, positive = NULL, patient_id = NULL,
                           batch = NULL, feature_prefix = NULL) {
      if (!is.null(target)) private$.target_col <- target
      if (!is.null(positive)) private$.positive_class <- positive
      if (!is.null(patient_id)) private$.patient_id_col <- patient_id
      if (!is.null(batch)) private$.batch_col <- batch
      if (!is.null(feature_prefix)) private$.feature_prefix <- feature_prefix

      private$.log_event("Column mapping updated")
      self$modified_at <- Sys.time()
      invisible(self)
    },

    #' @description Get column mappings
    get_mapping = function() {
      list(
        target = private$.target_col,
        positive = private$.positive_class,
        patient_id = private$.patient_id_col,
        batch = private$.batch_col,
        feature_prefix = private$.feature_prefix
      )
    },

    # =========================================================================
    # QC & VALIDATION (Traffic Light System)
    # =========================================================================

    #' @description Run QC checks and update status
    run_qc = function() {
      data <- private$.data_current
      target <- private$.target_col

      checks <- list(
        data_loaded = list(
          status = "red",
          message = "No data loaded",
          blocking = TRUE
        ),
        target_valid = list(
          status = "red",
          message = "Target column not set",
          blocking = TRUE
        ),
        class_balance = list(
          status = "green",
          message = "Classes balanced",
          blocking = FALSE
        ),
        missing_values = list(
          status = "green",
          message = "No missing values",
          blocking = FALSE
        ),
        patient_grouping = list(
          status = "yellow",
          message = "Patient ID not set (may cause leakage)",
          blocking = FALSE
        )
      )

      # Check: Data loaded
      if (!is.null(data)) {
        n_samples <- if (is.data.frame(data)) nrow(data) else nrow(data[[1]])
        n_features <- if (is.data.frame(data)) ncol(data) else sum(sapply(data, ncol))
        checks$data_loaded <- list(
          status = "green",
          message = sprintf("%d samples, %d features", n_samples, n_features),
          blocking = FALSE
        )
      }

      # Check: Target valid
      if (!is.null(target) && nzchar(target)) {
        target_vals <- private$.get_target_vector()
        if (!is.null(target_vals)) {
          n_classes <- length(unique(target_vals))
          checks$target_valid <- list(
            status = "green",
            message = sprintf("Target '%s' (%d classes)", target, n_classes),
            blocking = FALSE
          )

          # Check: Class balance
          if (n_classes == 2) {
            tbl <- table(target_vals)
            ratio <- min(tbl) / max(tbl)
            # Check severe imbalance first (more restrictive threshold)
            if (ratio < 0.1) {
              checks$class_balance <- list(
                status = "red",
                message = sprintf("Severe class imbalance (ratio: %.2f)", ratio),
                blocking = FALSE
              )
            } else if (ratio < 0.3) {
              checks$class_balance <- list(
                status = "yellow",
                message = sprintf("Class imbalance (ratio: %.2f)", ratio),
                blocking = FALSE
              )
            }
          }
        }
      }

      # Check: Missing values
      if (!is.null(data)) {
        x <- private$.get_feature_matrix()
        if (!is.null(x)) {
          missing_pct <- mean(is.na(x)) * 100
          if (missing_pct > 20) {
            checks$missing_values <- list(
              status = "red",
              message = sprintf("%.1f%% missing values", missing_pct),
              blocking = FALSE
            )
          } else if (missing_pct > 5) {
            checks$missing_values <- list(
              status = "yellow",
              message = sprintf("%.1f%% missing values", missing_pct),
              blocking = FALSE
            )
          } else if (missing_pct > 0) {
            checks$missing_values <- list(
              status = "green",
              message = sprintf("%.1f%% missing values", missing_pct),
              blocking = FALSE
            )
          }
        }
      }

      # Check: Patient ID
      if (!is.null(private$.patient_id_col) && nzchar(private$.patient_id_col)) {
        checks$patient_grouping <- list(
          status = "green",
          message = "Patient grouping enabled",
          blocking = FALSE
        )
      }

      private$.qc_status <- checks
      private$.log_event("QC completed")
      self$modified_at <- Sys.time()

      invisible(checks)
    },

    #' @description Get QC status
    get_qc_status = function() {
      private$.qc_status
    },

    #' @description Check if a workflow step can proceed
    #' @param step Step name: "qc", "pipeline", "benchmark", "results", "export"
    can_proceed = function(step) {
      switch(step,
        "qc" = !is.null(private$.data_current) && !is.null(private$.target_col),
        "pipeline" = {
          qc <- private$.qc_status
          !is.null(qc) && !any(sapply(qc, function(x) x$blocking && x$status == "red"))
        },
        "benchmark" = !is.null(private$.pipeline) && !is.null(private$.learner),
        "results" = !is.null(private$.benchmark_result),
        "export" = !is.null(private$.final_fit),
        FALSE
      )
    },

    # =========================================================================
    # PIPELINE & ML
    # =========================================================================

    #' @description Store pipeline object
    #' @param pipeline OmicPipeline object
    set_pipeline = function(pipeline) {
      private$.pipeline <- pipeline
      private$.log_event("Pipeline created")
      self$modified_at <- Sys.time()
      invisible(self)
    },

    #' @description Get pipeline
    get_pipeline = function() {
      private$.pipeline
    },

    #' @description Store learner object
    #' @param learner mlr3 learner
    set_learner = function(learner) {
      private$.learner <- learner
      private$.log_event("Learner configured")
      self$modified_at <- Sys.time()
      invisible(self)
    },

    #' @description Get learner
    get_learner = function() {
      private$.learner
    },

    #' @description Store benchmark results
    #' @param result Benchmark result object
    set_benchmark_result = function(result) {
      private$.benchmark_result <- result
      private$.log_event("Benchmark completed")
      self$modified_at <- Sys.time()
      self$save_checkpoint("benchmark_complete")
      invisible(self)
    },

    #' @description Get benchmark results
    get_benchmark_result = function() {
      private$.benchmark_result
    },

    #' @description Store final fitted model
    #' @param fit Final fit object
    set_final_fit = function(fit) {
      private$.final_fit <- fit
      private$.log_event("Final model fitted")
      self$modified_at <- Sys.time()
      self$save_checkpoint("final_fit")
      invisible(self)
    },

    #' @description Get final fit
    get_final_fit = function() {
      private$.final_fit
    },

    # =========================================================================
    # PERSISTENCE
    # =========================================================================

    #' @description Save project state to RDS
    #' @param filename Optional filename (default: project state)
    save = function(filename = NULL) {
      if (is.null(filename)) {
        filename <- file.path(self$path, "project_state.rds")
      }

      state <- list(
        id = self$id,
        name = self$name,
        path = self$path,
        created_at = self$created_at,
        modified_at = self$modified_at,
        data_raw = private$.data_raw,
        data_current = private$.data_current,
        is_multi_omics = private$.is_multi_omics,
        target_col = private$.target_col,
        positive_class = private$.positive_class,
        patient_id_col = private$.patient_id_col,
        batch_col = private$.batch_col,
        feature_prefix = private$.feature_prefix,
        qc_status = private$.qc_status,
        pipeline = private$.pipeline,
        learner = private$.learner,
        benchmark_result = private$.benchmark_result,
        final_fit = private$.final_fit,
        audit_log = private$.audit_log
      )

      saveRDS(state, filename)
      private$.log_event(paste0("Project saved: ", filename))
      invisible(filename)
    },

    #' @description Load project state from RDS
    #' @param filename Path to saved project
    load = function(filename) {
      if (!file.exists(filename)) {
        stop("Project file not found: ", filename, call. = FALSE)
      }

      state <- readRDS(filename)

      self$id <- state$id
      self$name <- state$name
      self$path <- state$path
      self$created_at <- state$created_at
      self$modified_at <- state$modified_at
      private$.data_raw <- state$data_raw
      private$.data_current <- state$data_current
      private$.is_multi_omics <- state$is_multi_omics
      private$.target_col <- state$target_col
      private$.positive_class <- state$positive_class
      private$.patient_id_col <- state$patient_id_col
      private$.batch_col <- state$batch_col
      private$.feature_prefix <- state$feature_prefix
      private$.qc_status <- state$qc_status
      private$.pipeline <- state$pipeline
      private$.learner <- state$learner
      private$.benchmark_result <- state$benchmark_result
      private$.final_fit <- state$final_fit
      private$.audit_log <- state$audit_log %||% character()

      private$.log_event(paste0("Project loaded: ", filename))
      invisible(self)
    },

    #' @description Save a checkpoint
    #' @param name Checkpoint name
    save_checkpoint = function(name) {
      filename <- file.path(self$path, sprintf("checkpoint_%s_%s.rds",
                                                name, format(Sys.time(), "%Y%m%d_%H%M%S")))
      self$save(filename)
      invisible(filename)
    },

    # =========================================================================
    # AUDIT LOG
    # =========================================================================

    #' @description Get audit log
    get_audit_log = function() {
      private$.audit_log
    },

    # =========================================================================
    # CODE GENERATION (for reproducibility)
    # =========================================================================

    #' @description Generate R code snippet for the current pipeline
    generate_code = function() {
      lines <- character()
      lines <- c(lines, "# OmicSelector Analysis Code")
      lines <- c(lines, "# Generated:", format(Sys.time()))
      lines <- c(lines, "library(OmicSelector)")
      lines <- c(lines, "")

      # Data loading (placeholder)
      lines <- c(lines, "# Load your data")
      lines <- c(lines, "# data <- read.csv('your_data.csv')")
      lines <- c(lines, "")

      # Pipeline creation
      if (!is.null(private$.pipeline)) {
        mapping <- self$get_mapping()
        lines <- c(lines, "# Create pipeline")
        lines <- c(lines, sprintf("pipeline <- OmicPipeline$new("))
        lines <- c(lines, sprintf("  data = data,"))
        lines <- c(lines, sprintf("  target = '%s',", mapping$target %||% "outcome"))
        if (!is.null(mapping$positive)) {
          lines <- c(lines, sprintf("  positive = '%s',", mapping$positive))
        }
        if (!is.null(mapping$patient_id) && nzchar(mapping$patient_id)) {
          lines <- c(lines, sprintf("  patient_id = '%s',", mapping$patient_id))
        }
        if (!is.null(mapping$batch) && nzchar(mapping$batch)) {
          lines <- c(lines, sprintf("  batch = '%s',", mapping$batch))
        }
        lines <- c(lines, ")")
        lines <- c(lines, "")
      }

      paste(lines, collapse = "\n")
    },

    # =========================================================================
    # SUMMARY
    # =========================================================================

    #' @description Get project summary
    summary = function() {
      data <- private$.data_current
      n_samples <- if (!is.null(data)) {
        if (is.data.frame(data)) nrow(data) else nrow(data[[1]])
      } else NA

      n_features <- if (!is.null(data)) {
        if (is.data.frame(data)) ncol(data) else sum(sapply(data, ncol))
      } else NA

      list(
        id = self$id,
        name = self$name,
        path = self$path,
        created_at = self$created_at,
        modified_at = self$modified_at,
        n_samples = n_samples,
        n_features = n_features,
        is_multi_omics = private$.is_multi_omics,
        target = private$.target_col,
        has_pipeline = !is.null(private$.pipeline),
        has_benchmark = !is.null(private$.benchmark_result),
        has_final_fit = !is.null(private$.final_fit),
        qc_overall = private$.get_overall_qc_status()
      )
    }
  ),

  private = list(
    .data_raw = NULL,
    .data_current = NULL,
    .is_multi_omics = FALSE,
    .target_col = NULL,
    .positive_class = NULL,
    .patient_id_col = NULL,
    .batch_col = NULL,
    .feature_prefix = NULL,
    .qc_status = NULL,
    .pipeline = NULL,
    .learner = NULL,
    .benchmark_result = NULL,
    .final_fit = NULL,
    .audit_log = character(),

    .generate_id = function() {
      stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
      rand <- substr(digest::digest(runif(1)), 1, 6)
      paste("project", stamp, rand, sep = "-")
    },

    .log_event = function(message) {
      entry <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message)
      private$.audit_log <- c(private$.audit_log, entry)
    },

    .get_target_vector = function() {
      data <- private$.data_current
      target <- private$.target_col
      if (is.null(data) || is.null(target)) return(NULL)

      if (is.data.frame(data)) {
        if (target %in% names(data)) return(data[[target]])
      } else if (is.list(data)) {
        for (mod in names(data)) {
          if (target %in% names(data[[mod]])) {
            return(data[[mod]][[target]])
          }
        }
      }
      NULL
    },

    .get_feature_matrix = function() {
      data <- private$.data_current
      target <- private$.target_col
      if (is.null(data)) return(NULL)

      if (is.data.frame(data)) {
        meta_cols <- c(target, private$.patient_id_col, private$.batch_col, "mix")
        meta_cols <- meta_cols[!is.null(meta_cols) & nzchar(meta_cols)]
        keep <- setdiff(names(data), meta_cols)
        x <- data[, keep, drop = FALSE]
        numeric_cols <- sapply(x, is.numeric)
        return(x[, numeric_cols, drop = FALSE])
      }
      NULL
    },

    .get_overall_qc_status = function() {
      qc <- private$.qc_status
      if (is.null(qc)) return("unknown")

      statuses <- sapply(qc, function(x) x$status)
      if (any(statuses == "red")) return("red")
      if (any(statuses == "yellow")) return("yellow")
      "green"
    }
  )
)

#' @title Null coalescing operator
#' @param x Left value
#' @param y Default value if x is NULL
#' @return x if not NULL, else y
`%||%` <- function(x, y) if (is.null(x)) y else x
