#' @title Multi-Omics Support for OmicSelector 2.0
#'
#' @description
#' Functions and utilities for handling multi-omics data containers.
#' Designed to future-proof the API for Phase 4 integration.
#'
#' @details
#' This module provides:
#' - Validation and normalization of multi-omics input
#' - Modality namespacing (e.g., rna::TP53, prot::TP53)
#' - Late integration support for ensemble models
#'
#' @name multi-omics
NULL


#' @title Validate Multi-Omics Input
#'
#' @description
#' Validates and normalizes multi-omics input data. Accepts either a single
#' data.frame (backward compatible) or a named list of data matrices.
#'
#' @param data Either a data.frame or a named list of data.frames/matrices
#' @param target Name of the target column
#' @param require_common_samples Logical, whether to require sample alignment
#'
#' @return A normalized list structure with validated omics data
#'
#' @examples
#' \dontrun{
#' # Single data.frame (backward compatible)
#' result <- validate_omics_input(my_data, target = "outcome")
#'
#' # Multi-omics list
#' multi_data <- list(
#'   rna = rna_expression,
#'   mirna = mirna_expression,
#'   prot = protein_levels
#' )
#' result <- validate_omics_input(multi_data, target = "outcome")
#' }
#'
#' @export
validate_omics_input <- function(data, target, require_common_samples = TRUE) {

  result <- list(
    is_multi_omics = FALSE,
    modalities = list(),
    target_data = NULL,
    sample_ids = NULL,
    feature_map = list(),
    warnings = character(0)
  )

  # Case 1: Single data.frame
  if (is.data.frame(data)) {
    # Validate target column exists
    if (!target %in% names(data)) {
      stop(sprintf("Target column '%s' not found in data", target), call. = FALSE)
    }

    result$is_multi_omics <- FALSE
    result$modalities$main <- .prepare_modality(data, target, modality_name = "main")
    result$target_data <- data[[target]]
    result$sample_ids <- rownames(data)
    result$feature_map$main <- setdiff(names(data), target)
    return(structure(result, class = c("OmicsInput", "list")))
  }

  # Case 2: Named list of matrices/data.frames
  if (is.list(data) && !is.data.frame(data)) {

    # Validate list has names
    if (is.null(names(data)) || any(names(data) == "")) {
      stop("Multi-omics data must be a named list (e.g., list(rna = ..., mirna = ...))",
           call. = FALSE)
    }

    result$is_multi_omics <- TRUE

    # Find target column across modalities
    target_found <- FALSE
    target_modality <- NULL

    for (mod_name in names(data)) {
      mod_data <- data[[mod_name]]

      # Convert matrix to data.frame if needed
      if (is.matrix(mod_data)) {
        mod_data <- as.data.frame(mod_data)
      }

      if (!is.data.frame(mod_data)) {
        stop(sprintf("Modality '%s' must be a data.frame or matrix", mod_name),
             call. = FALSE)
      }

      # Check for target in this modality
      if (target %in% names(mod_data)) {
        if (target_found) {
          stop(sprintf("Target '%s' found in multiple modalities. Please include target in only one modality.",
                       target), call. = FALSE)
        }
        target_found <- TRUE
        target_modality <- mod_name
        result$target_data <- mod_data[[target]]
      }

      # Namespace features with modality prefix
      feature_cols <- setdiff(names(mod_data), target)
      namespaced_features <- paste0(mod_name, "::", feature_cols)
      result$feature_map[[mod_name]] <- setNames(namespaced_features, feature_cols)

      # Store modality data
      result$modalities[[mod_name]] <- .prepare_modality(
        mod_data, target, modality_name = mod_name
      )
    }

    if (!target_found) {
      stop(sprintf("Target column '%s' not found in any modality", target),
           call. = FALSE)
    }

    # Sample alignment check
    if (require_common_samples && length(result$modalities) > 1) {
      sample_sets <- lapply(result$modalities, function(m) m$sample_ids)
      common_samples <- Reduce(intersect, sample_sets)

      if (length(common_samples) == 0) {
        stop("No common samples found across modalities", call. = FALSE)
      }

      if (length(common_samples) < min(sapply(sample_sets, length))) {
        result$warnings <- c(result$warnings,
          sprintf("Only %d common samples found across %d modalities",
                  length(common_samples), length(sample_sets)))
      }

      result$sample_ids <- common_samples
    } else {
      # Use samples from target modality
      result$sample_ids <- result$modalities[[target_modality]]$sample_ids
    }

    return(structure(result, class = c("OmicsInput", "list")))
  }

  stop("data must be a data.frame or a named list of data.frames", call. = FALSE)
}


#' @title Prepare Single Modality Data
#' @keywords internal
.prepare_modality <- function(data, target, modality_name) {
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  list(
    name = modality_name,
    data = data,
    features = setdiff(names(data), target),
    n_features = length(setdiff(names(data), target)),
    sample_ids = rownames(data),
    has_target = target %in% names(data)
  )
}


#' @title Merge Multi-Omics Data for Analysis
#'
#' @description
#' Merges multiple omics modalities into a single data.frame with
#' namespaced features for downstream analysis.
#'
#' @param omics_input A validated OmicsInput object
#' @param sample_subset Optional character vector of sample IDs to include
#'
#' @return A merged data.frame with namespaced feature columns
#'
#' @export
merge_omics_data <- function(omics_input, sample_subset = NULL) {

  if (!inherits(omics_input, "OmicsInput")) {
    stop("Input must be a validated OmicsInput object from validate_omics_input()",
         call. = FALSE)
  }

  # Use common samples if not specified
  if (is.null(sample_subset)) {
    sample_subset <- omics_input$sample_ids
  }

  # For single modality, just return the data
  if (!omics_input$is_multi_omics) {
    mod <- omics_input$modalities$main
    return(mod$data[sample_subset, , drop = FALSE])
  }

  # Merge multiple modalities
  merged <- NULL

  for (mod_name in names(omics_input$modalities)) {
    mod <- omics_input$modalities[[mod_name]]
    mod_data <- mod$data

    # Subset to common samples
    common <- intersect(sample_subset, mod$sample_ids)
    mod_data <- mod_data[common, , drop = FALSE]

    # Rename features with namespace
    feature_cols <- mod$features
    new_names <- paste0(mod_name, "::", feature_cols)
    names(mod_data)[match(feature_cols, names(mod_data))] <- new_names

    if (is.null(merged)) {
      merged <- mod_data
    } else {
      # Merge by row names (sample IDs)
      merged <- merge(merged, mod_data, by = "row.names", all = FALSE)
      rownames(merged) <- merged$Row.names
      merged$Row.names <- NULL
    }
  }

  return(merged)
}


#' @title Get Modality Information
#'
#' @description
#' Returns summary information about each modality in a multi-omics dataset.
#'
#' @param omics_input A validated OmicsInput object
#'
#' @return A data.frame summarizing each modality
#'
#' @export
get_modality_info <- function(omics_input) {

  if (!inherits(omics_input, "OmicsInput")) {
    stop("Input must be a validated OmicsInput object", call. = FALSE)
  }

  info <- data.frame(
    modality = character(0),
    n_features = integer(0),
    n_samples = integer(0),
    has_target = logical(0),
    stringsAsFactors = FALSE
  )

  for (mod_name in names(omics_input$modalities)) {
    mod <- omics_input$modalities[[mod_name]]
    info <- rbind(info, data.frame(
      modality = mod_name,
      n_features = mod$n_features,
      n_samples = length(mod$sample_ids),
      has_target = mod$has_target,
      stringsAsFactors = FALSE
    ))
  }

  return(info)
}


#' @title Print method for OmicsInput
#' @param x An OmicsInput object
#' @param ... Additional arguments (ignored)
#' @export
print.OmicsInput <- function(x, ...) {
  cat("=== OmicsInput ===\n")
  cat(sprintf("Type: %s\n", if (x$is_multi_omics) "Multi-Omics" else "Single Modality"))
  cat(sprintf("Total samples: %d\n", length(x$sample_ids)))

  cat("\nModalities:\n")
  for (mod_name in names(x$modalities)) {
    mod <- x$modalities[[mod_name]]
    cat(sprintf("  %s: %d features, %d samples%s\n",
                mod_name, mod$n_features, length(mod$sample_ids),
                if (mod$has_target) " (contains target)" else ""))
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat("  - ", w, "\n")
    }
  }

  invisible(x)
}


# =============================================================================
# Multi-Omics Late Integration (Stacking) - Super Learner Methodology
# =============================================================================

#' @title OmicModalitySpec R6 Class
#'
#' @description
#' Specification for a single omics modality including data, preprocessing, and learner.
#' Encapsulates everything needed to train/predict on one data type.
#'
#' @export
OmicModalitySpec <- R6::R6Class(
  "OmicModalitySpec",

  public = list(
    #' @field id Unique identifier for this modality (e.g., "mRNA", "miRNA")
    id = NULL,

    #' @field x Data matrix (samples x features)
    x = NULL,

    #' @field learner mlr3 Learner object
    learner = NULL,

    #' @field preproc Optional mlr3pipelines Graph or PipeOp for preprocessing
    preproc = NULL,

    #' @description
    #' Create a new modality specification
    #'
    #' @param id Character, unique modality identifier
    #' @param x Matrix or data.frame with samples in rows, features in columns.
    #'        Rownames should be sample IDs.
    #' @param learner mlr3 Learner or character shortcut (e.g., "classif.ranger")
    #' @param preproc Optional preprocessing pipeline (Graph or PipeOp)
    #'
    #' @return An OmicModalitySpec object
    initialize = function(id, x, learner, preproc = NULL) {
      checkmate::assert_string(id)
      checkmate::assert_true(is.matrix(x) || is.data.frame(x))

      self$id <- id
      self$x <- as.matrix(x)

      # Accept learner as string or Learner object
      if (is.character(learner)) {
        self$learner <- mlr3::lrn(learner)
      } else {
        checkmate::assert_class(learner, "Learner")
        self$learner <- learner$clone(deep = TRUE)
      }

      self$preproc <- preproc

      # Validate sample IDs
      if (is.null(rownames(self$x))) {
        warning(sprintf("Modality '%s' has no rownames. Using row indices as sample IDs.", id))
        rownames(self$x) <- as.character(seq_len(nrow(self$x)))
      }

      invisible(self)
    },

    #' @description
    #' Get sample IDs for this modality
    #' @return Character vector of sample IDs
    get_sample_ids = function() {
      rownames(self$x)
    },

    #' @description
    #' Get number of features
    #' @return Integer
    n_features = function() {
      ncol(self$x)
    },

    #' @description
    #' Create a GraphLearner combining preprocessing and learner
    #' @return GraphLearner or Learner
    make_graph_learner = function() {
      if (is.null(self$preproc)) {
        lrn <- self$learner$clone(deep = TRUE)
        lrn$id <- paste0(self$id, "_learner")
        return(lrn)
      }

      # Combine preproc and learner into GraphLearner
      if (requireNamespace("mlr3pipelines", quietly = TRUE)) {
        graph <- self$preproc %>>% self$learner
        gl <- mlr3pipelines::as_learner(graph)
        gl$id <- paste0(self$id, "_graph")
        return(gl)
      } else {
        warning("mlr3pipelines not available. Returning learner without preprocessing.")
        return(self$learner$clone(deep = TRUE))
      }
    },

    #' @description
    #' Create an mlr3 Task from this modality's data
    #'
    #' @param y Target vector (must align with rownames of x)
    #' @param target_name Name for target column
    #' @param task_type "classif" or "regr"
    #' @param sample_ids Optional subset of sample IDs to use
    #'
    #' @return TaskClassif or TaskRegr
    as_task = function(y, target_name = "target", task_type = c("classif", "regr"),
                       sample_ids = NULL) {
      task_type <- match.arg(task_type)

      # Handle sample subsetting
      if (!is.null(sample_ids)) {
        x_subset <- self$x[sample_ids, , drop = FALSE]
        y_subset <- y[sample_ids]
      } else {
        x_subset <- self$x
        y_subset <- y
      }

      # Create data frame
      dat <- data.frame(y_subset, x_subset, check.names = FALSE)
      names(dat)[1] <- target_name

      # Create task
      if (task_type == "classif") {
        mlr3::TaskClassif$new(
          id = paste0("task_", self$id),
          backend = dat,
          target = target_name
        )
      } else {
        mlr3::TaskRegr$new(
          id = paste0("task_", self$id),
          backend = dat,
          target = target_name
        )
      }
    },

    #' @description
    #' Print method
    print = function() {
      cat(sprintf("=== OmicModalitySpec: %s ===\n", self$id))
      cat(sprintf("  Samples: %d\n", nrow(self$x)))
      cat(sprintf("  Features: %d\n", ncol(self$x)))
      cat(sprintf("  Learner: %s\n", self$learner$id))
      cat(sprintf("  Preprocessing: %s\n", if (is.null(self$preproc)) "None" else "Yes"))
      invisible(self)
    }
  )
)


#' @title OmicStackedEnsemble R6 Class
#'
#' @description
#' Orchestrates multi-omics stacking with leakage-free OOF prediction generation.
#' Trains base learners per modality and combines using a meta-learner.
#'
#' Based on Super Learner methodology (van der Laan et al., 2007).
#'
#' @references
#' van der Laan, M.J., Polley, E.C., Hubbard, A.E. (2007). Super Learner.
#' Statistical Applications in Genetics and Molecular Biology, 6(1).
#'
#' @export
OmicStackedEnsemble <- R6::R6Class(
  "OmicStackedEnsemble",

  public = list(
    #' @field modalities List of OmicModalitySpec objects
    modalities = NULL,

    #' @field meta_learner Learner for combining base predictions
    meta_learner = NULL,

    #' @field resampling Resampling strategy for OOF predictions
    resampling = NULL,

    #' @field prediction_type "prob" for probabilities, "response" for class labels
    prediction_type = NULL,

    #' @field task_type "classif" or "regr"
    task_type = NULL,

    #' @description
    #' Create a new stacked ensemble
    #'
    #' @param modalities List of OmicModalitySpec objects
    #' @param meta_learner Learner for stacking (default: logistic regression)
    #' @param resampling Resampling strategy (default: 5-fold CV)
    #' @param prediction_type "prob" or "response"
    #' @param task_type "classif" or "regr"
    #'
    #' @return An OmicStackedEnsemble object
    initialize = function(modalities,
                          meta_learner = NULL,
                          resampling = NULL,
                          prediction_type = "prob",
                          task_type = "classif") {

      # Validate modalities
      checkmate::assert_list(modalities, min.len = 2)
      for (m in modalities) {
        checkmate::assert_class(m, "OmicModalitySpec")
      }
      self$modalities <- modalities

      # Default meta-learner: logistic regression (interpretable)
      if (is.null(meta_learner)) {
        if (task_type == "classif") {
          self$meta_learner <- mlr3::lrn("classif.log_reg", predict_type = "prob")
        } else {
          self$meta_learner <- mlr3::lrn("regr.lm")
        }
      } else if (is.character(meta_learner)) {
        self$meta_learner <- mlr3::lrn(meta_learner)
      } else {
        self$meta_learner <- meta_learner$clone(deep = TRUE)
      }

      # Default resampling: 5-fold CV
      if (is.null(resampling)) {
        self$resampling <- mlr3::rsmp("cv", folds = 5)
      } else if (is.character(resampling)) {
        self$resampling <- mlr3::rsmp(resampling)
      } else {
        self$resampling <- resampling
      }

      self$prediction_type <- prediction_type
      self$task_type <- task_type

      invisible(self)
    },

    #' @description
    #' Train the stacked ensemble
    #'
    #' @param y Target vector (named with sample IDs)
    #' @param target_name Name for target column
    #'
    #' @return Self (invisibly), for method chaining
    train = function(y, target_name = "target") {
      # 1. Align samples across all modalities
      aligned <- private$.align_samples(y)
      common_ids <- aligned$common_ids
      y_aligned <- aligned$y

      private$.common_ids <- common_ids
      private$.target_name <- target_name
      private$.classes <- if (self$task_type == "classif") levels(as.factor(y_aligned)) else NULL

      message(sprintf("Training stacked ensemble on %d samples across %d modalities",
                      length(common_ids), length(self$modalities)))

      # 2. Generate OOF predictions for each modality
      stack_parts <- list()

      for (m in self$modalities) {
        message(sprintf("  Processing modality: %s (%d features)",
                        m$id, m$n_features()))

        # Subset to common samples
        x_m <- m$x[common_ids, , drop = FALSE]

        # Create temporary modality with aligned data
        m_aligned <- OmicModalitySpec$new(
          id = m$id,
          x = x_m,
          learner = m$learner$clone(deep = TRUE),
          preproc = m$preproc
        )

        # Create task
        task_m <- m_aligned$as_task(
          y = y_aligned,
          target_name = target_name,
          task_type = self$task_type
        )

        # Get GraphLearner
        lrn_m <- m_aligned$make_graph_learner()

        # Set prediction type for classification
        if (self$task_type == "classif" && self$prediction_type == "prob") {
          lrn_m$predict_type <- "prob"
        }

        # Run resampling to get OOF predictions
        rr <- mlr3::resample(
          task = task_m,
          learner = lrn_m,
          resampling = self$resampling,
          store_models = FALSE
        )

        # Extract OOF predictions
        pred <- rr$prediction()

        # Build stacking features
        if (self$task_type == "classif" && self$prediction_type == "prob") {
          prob <- as.data.frame(pred$prob)
          colnames(prob) <- paste0(m$id, "_prob_", colnames(prob))
          stack_parts[[m$id]] <- prob
        } else {
          stack_parts[[m$id]] <- data.frame(
            setNames(list(pred$response), paste0(m$id, "_pred"))
          )
        }
      }

      # 3. Combine into stacking matrix
      X_stack <- do.call(cbind, stack_parts)
      private$.stack_features <- colnames(X_stack)

      # Create meta-task
      stack_df <- data.frame(
        target = y_aligned,
        X_stack,
        check.names = FALSE
      )
      names(stack_df)[1] <- target_name

      private$.oof_stack <- stack_df

      message("  Training meta-learner on OOF predictions")

      # 4. Fit meta-learner on OOF stacking features
      if (self$task_type == "classif") {
        meta_task <- mlr3::TaskClassif$new(
          id = "meta_task",
          backend = stack_df,
          target = target_name
        )
      } else {
        meta_task <- mlr3::TaskRegr$new(
          id = "meta_task",
          backend = stack_df,
          target = target_name
        )
      }

      if (self$task_type == "classif" && self$prediction_type == "prob") {
        self$meta_learner$predict_type <- "prob"
      }

      private$.fitted_meta <- self$meta_learner$clone(deep = TRUE)
      private$.fitted_meta$train(meta_task)

      # 5. Fit final base learners on full data (for prediction)
      message("  Training final base models on full data")

      private$.fitted_base <- lapply(self$modalities, function(m) {
        x_m <- m$x[common_ids, , drop = FALSE]

        m_aligned <- OmicModalitySpec$new(
          id = m$id,
          x = x_m,
          learner = m$learner$clone(deep = TRUE),
          preproc = m$preproc
        )

        task_m <- m_aligned$as_task(
          y = y_aligned,
          target_name = target_name,
          task_type = self$task_type
        )

        lrn_m <- m_aligned$make_graph_learner()

        if (self$task_type == "classif" && self$prediction_type == "prob") {
          lrn_m$predict_type <- "prob"
        }

        lrn_m$train(task_m)
        lrn_m
      })
      names(private$.fitted_base) <- sapply(self$modalities, function(m) m$id)

      private$.fitted <- TRUE

      message("Stacked ensemble training complete")
      invisible(self)
    },

    #' @description
    #' Predict using the stacked ensemble
    #'
    #' @param newdata_list Named list of data matrices, one per modality.
    #'        Names must match modality IDs. Missing modalities can be NULL.
    #'
    #' @return Prediction object
    predict = function(newdata_list) {
      if (!private$.fitted) {
        stop("Ensemble not trained. Call $train() first.")
      }

      checkmate::assert_list(newdata_list, names = "named")

      # Check for missing modalities
      modality_ids <- sapply(self$modalities, function(m) m$id)
      provided <- names(newdata_list)
      missing <- setdiff(modality_ids, provided)

      if (length(missing) > 0) {
        warning(sprintf("Missing modalities: %s. Will impute with training prior.",
                        paste(missing, collapse = ", ")))
      }

      # Align samples across provided modalities
      available_data <- newdata_list[!sapply(newdata_list, is.null)]

      if (length(available_data) == 0) {
        stop("No modality data provided")
      }

      ids_list <- lapply(available_data, rownames)
      common_ids <- Reduce(intersect, ids_list)

      if (length(common_ids) == 0) {
        stop("No overlapping sample IDs across provided modalities")
      }

      # Generate base predictions
      stack_parts <- list()

      for (m_id in modality_ids) {
        if (!is.null(newdata_list[[m_id]])) {
          # Get base model predictions
          x_new <- newdata_list[[m_id]][common_ids, , drop = FALSE]

          # Create dummy task for prediction
          dummy_y <- rep(NA, length(common_ids))
          if (self$task_type == "classif") {
            dummy_y <- factor(dummy_y, levels = private$.classes)
          }

          dummy_df <- data.frame(
            target = dummy_y,
            x_new,
            check.names = FALSE
          )
          names(dummy_df)[1] <- private$.target_name

          if (self$task_type == "classif") {
            pred_task <- mlr3::TaskClassif$new(
              id = paste0("pred_", m_id),
              backend = dummy_df,
              target = private$.target_name
            )
          } else {
            pred_task <- mlr3::TaskRegr$new(
              id = paste0("pred_", m_id),
              backend = dummy_df,
              target = private$.target_name
            )
          }

          pred <- private$.fitted_base[[m_id]]$predict(pred_task)

          # Extract predictions
          if (self$task_type == "classif" && self$prediction_type == "prob") {
            prob <- as.data.frame(pred$prob)
            colnames(prob) <- paste0(m_id, "_prob_", colnames(prob))
            stack_parts[[m_id]] <- prob
          } else {
            stack_parts[[m_id]] <- data.frame(
              setNames(list(pred$response), paste0(m_id, "_pred"))
            )
          }
        } else {
          # Impute missing modality with training prior
          n_samples <- length(common_ids)

          if (self$task_type == "classif" && self$prediction_type == "prob") {
            # Use class proportions from training
            prior <- table(private$.oof_stack[[private$.target_name]]) / nrow(private$.oof_stack)
            prob <- matrix(rep(as.numeric(prior), each = n_samples),
                           nrow = n_samples, ncol = length(prior))
            prob <- as.data.frame(prob)
            colnames(prob) <- paste0(m_id, "_prob_", names(prior))
            stack_parts[[m_id]] <- prob
          } else {
            # Use mean response from training
            mean_pred <- mean(as.numeric(private$.oof_stack[[private$.target_name]]))
            stack_parts[[m_id]] <- data.frame(
              setNames(list(rep(mean_pred, n_samples)), paste0(m_id, "_pred"))
            )
          }
        }
      }

      # Combine into stacking matrix
      X_stack <- do.call(cbind, stack_parts)

      # Ensure column order matches training
      X_stack <- X_stack[, private$.stack_features, drop = FALSE]

      # Create meta-task for prediction
      dummy_y <- rep(NA, nrow(X_stack))
      if (self$task_type == "classif") {
        dummy_y <- factor(dummy_y, levels = private$.classes)
      }

      meta_df <- data.frame(
        target = dummy_y,
        X_stack,
        check.names = FALSE
      )
      names(meta_df)[1] <- private$.target_name

      if (self$task_type == "classif") {
        meta_task <- mlr3::TaskClassif$new(
          id = "meta_pred",
          backend = meta_df,
          target = private$.target_name
        )
      } else {
        meta_task <- mlr3::TaskRegr$new(
          id = "meta_pred",
          backend = meta_df,
          target = private$.target_name
        )
      }

      # Meta-learner prediction
      private$.fitted_meta$predict(meta_task)
    },

    #' @description
    #' Check if ensemble is fitted
    #' @return Logical
    is_fitted = function() {
      private$.fitted
    },

    #' @description
    #' Get OOF stacking features (for auditing)
    #' @return Data frame with OOF predictions used for meta-learner training
    get_oof_stack = function() {
      if (!private$.fitted) {
        stop("Ensemble not trained")
      }
      private$.oof_stack
    },

    #' @description
    #' Get modality contribution (meta-learner coefficients if linear)
    #' @return Named vector or NULL if meta-learner is not linear
    get_modality_weights = function() {
      if (!private$.fitted) {
        stop("Ensemble not trained")
      }

      # Try to extract coefficients (works for log_reg, lm)
      tryCatch({
        model <- private$.fitted_meta$model
        if (inherits(model, "glm") || inherits(model, "lm")) {
          coef(model)
        } else {
          NULL
        }
      }, error = function(e) NULL)
    },

    #' @description
    #' Print method
    print = function() {
      cat("=== OmicStackedEnsemble ===\n")
      cat(sprintf("  Modalities: %d\n", length(self$modalities)))
      for (m in self$modalities) {
        cat(sprintf("    - %s: %d features, %s\n",
                    m$id, m$n_features(), m$learner$id))
      }
      cat(sprintf("  Meta-learner: %s\n", self$meta_learner$id))
      cat(sprintf("  Resampling: %s\n", self$resampling$id))
      cat(sprintf("  Fitted: %s\n", private$.fitted))
      if (private$.fitted) {
        cat(sprintf("  Training samples: %d\n", length(private$.common_ids)))
      }
      invisible(self)
    }
  ),

  private = list(
    .fitted = FALSE,
    .fitted_base = NULL,
    .fitted_meta = NULL,
    .oof_stack = NULL,
    .common_ids = NULL,
    .target_name = NULL,
    .stack_features = NULL,
    .classes = NULL,

    .align_samples = function(y) {
      # Get sample IDs from all modalities
      ids_list <- lapply(self$modalities, function(m) m$get_sample_ids())
      modality_names <- sapply(self$modalities, function(m) m$id)

      # Find intersection
      common_ids <- Reduce(intersect, ids_list)

      if (length(common_ids) == 0) {
        stop("No overlapping sample IDs across modalities. Check rownames.")
      }

      # Report alignment
      for (i in seq_along(self$modalities)) {
        n_orig <- length(ids_list[[i]])
        n_common <- length(common_ids)
        if (n_orig != n_common) {
          message(sprintf("  Modality '%s': %d -> %d samples after alignment",
                          modality_names[i], n_orig, n_common))
        }
      }

      # Align y
      if (!is.null(names(y))) {
        # Named vector: subset by common_ids
        y_aligned <- y[common_ids]
        if (any(is.na(y_aligned))) {
          stop("Missing y values for some common samples")
        }
      } else {
        # Assume y is already aligned with first modality
        warning("y is not named. Assuming alignment with first modality's rownames.")
        first_ids <- ids_list[[1]]
        idx <- match(common_ids, first_ids)
        y_aligned <- y[idx]
      }

      list(common_ids = common_ids, y = y_aligned)
    }
  )
)


#' @title Simple Weighted Averaging Ensemble
#'
#' @description
#' Alternative to stacking: weighted average of base model predictions.
#' Simpler and more interpretable, often performs comparably.
#'
#' @export
OmicWeightedEnsemble <- R6::R6Class(
  "OmicWeightedEnsemble",

  public = list(
    #' @field modalities List of OmicModalitySpec objects
    modalities = NULL,

    #' @field weights Named numeric vector of modality weights
    weights = NULL,

    #' @field resampling Resampling strategy for weight estimation
    resampling = NULL,

    #' @field task_type "classif" or "regr"
    task_type = NULL,

    #' @description
    #' Create a new weighted ensemble
    #'
    #' @param modalities List of OmicModalitySpec objects
    #' @param weights Optional named weights. If NULL, estimated from OOF performance.
    #' @param resampling Resampling strategy
    #' @param task_type "classif" or "regr"
    #'
    #' @return An OmicWeightedEnsemble object
    initialize = function(modalities,
                          weights = NULL,
                          resampling = NULL,
                          task_type = "classif") {
      checkmate::assert_list(modalities, min.len = 2)
      self$modalities <- modalities
      self$weights <- weights

      if (is.null(resampling)) {
        self$resampling <- mlr3::rsmp("cv", folds = 5)
      } else {
        self$resampling <- resampling
      }

      self$task_type <- task_type

      invisible(self)
    },

    #' @description
    #' Train the weighted ensemble
    #'
    #' @param y Target vector
    #' @param target_name Name for target column
    #' @param weight_method How to estimate weights: "inverse_logloss", "auc", "uniform"
    #'
    #' @return Self (invisibly)
    train = function(y, target_name = "target", weight_method = c("inverse_logloss", "auc", "uniform")) {
      weight_method <- match.arg(weight_method)

      # Align samples
      ids_list <- lapply(self$modalities, function(m) m$get_sample_ids())
      common_ids <- Reduce(intersect, ids_list)

      if (!is.null(names(y))) {
        y_aligned <- y[common_ids]
      } else {
        idx <- match(common_ids, ids_list[[1]])
        y_aligned <- y[idx]
      }

      private$.common_ids <- common_ids
      private$.target_name <- target_name
      private$.classes <- if (self$task_type == "classif") levels(as.factor(y_aligned)) else NULL

      # If weights not provided, estimate from OOF performance
      if (is.null(self$weights)) {
        message("Estimating weights from OOF performance")

        performances <- c()

        for (m in self$modalities) {
          x_m <- m$x[common_ids, , drop = FALSE]

          m_aligned <- OmicModalitySpec$new(
            id = m$id,
            x = x_m,
            learner = m$learner$clone(deep = TRUE),
            preproc = m$preproc
          )

          task_m <- m_aligned$as_task(
            y = y_aligned,
            target_name = target_name,
            task_type = self$task_type
          )

          lrn_m <- m_aligned$make_graph_learner()
          if (self$task_type == "classif") {
            lrn_m$predict_type <- "prob"
          }

          rr <- mlr3::resample(
            task = task_m,
            learner = lrn_m,
            resampling = self$resampling,
            store_models = FALSE
          )

          # Compute performance metric
          if (weight_method == "auc" && self$task_type == "classif") {
            perf <- rr$aggregate(mlr3::msr("classif.auc"))
          } else if (weight_method == "inverse_logloss" && self$task_type == "classif") {
            perf <- 1 / (rr$aggregate(mlr3::msr("classif.logloss")) + 0.01)
          } else {
            perf <- 1  # Uniform
          }

          performances[m$id] <- perf
        }

        # Normalize to sum to 1
        self$weights <- performances / sum(performances)
        message(sprintf("Estimated weights: %s",
                        paste(sprintf("%s=%.3f", names(self$weights), self$weights), collapse = ", ")))
      }

      # Train final base models
      private$.fitted_base <- lapply(self$modalities, function(m) {
        x_m <- m$x[common_ids, , drop = FALSE]

        m_aligned <- OmicModalitySpec$new(
          id = m$id,
          x = x_m,
          learner = m$learner$clone(deep = TRUE),
          preproc = m$preproc
        )

        task_m <- m_aligned$as_task(
          y = y_aligned,
          target_name = target_name,
          task_type = self$task_type
        )

        lrn_m <- m_aligned$make_graph_learner()
        if (self$task_type == "classif") {
          lrn_m$predict_type <- "prob"
        }

        lrn_m$train(task_m)
        lrn_m
      })
      names(private$.fitted_base) <- sapply(self$modalities, function(m) m$id)

      private$.fitted <- TRUE
      invisible(self)
    },

    #' @description
    #' Predict using weighted averaging
    #'
    #' @param newdata_list Named list of data matrices
    #'
    #' @return Data frame with weighted average predictions
    predict = function(newdata_list) {
      if (!private$.fitted) {
        stop("Ensemble not trained")
      }

      # Align samples
      available_data <- newdata_list[!sapply(newdata_list, is.null)]
      ids_list <- lapply(available_data, rownames)
      common_ids <- Reduce(intersect, ids_list)

      modality_ids <- sapply(self$modalities, function(m) m$id)

      # Collect predictions
      all_probs <- list()

      for (m_id in modality_ids) {
        if (!is.null(newdata_list[[m_id]])) {
          x_new <- newdata_list[[m_id]][common_ids, , drop = FALSE]

          dummy_y <- rep(NA, length(common_ids))
          if (self$task_type == "classif") {
            dummy_y <- factor(dummy_y, levels = private$.classes)
          }

          dummy_df <- data.frame(target = dummy_y, x_new, check.names = FALSE)
          names(dummy_df)[1] <- private$.target_name

          if (self$task_type == "classif") {
            pred_task <- mlr3::TaskClassif$new(
              id = paste0("pred_", m_id),
              backend = dummy_df,
              target = private$.target_name
            )
          } else {
            pred_task <- mlr3::TaskRegr$new(
              id = paste0("pred_", m_id),
              backend = dummy_df,
              target = private$.target_name
            )
          }

          pred <- private$.fitted_base[[m_id]]$predict(pred_task)

          if (self$task_type == "classif") {
            all_probs[[m_id]] <- as.matrix(pred$prob)
          } else {
            all_probs[[m_id]] <- pred$response
          }
        }
      }

      # Weighted average
      if (self$task_type == "classif") {
        # Average probabilities
        weighted_probs <- Reduce(`+`, lapply(names(all_probs), function(m_id) {
          all_probs[[m_id]] * self$weights[m_id]
        }))

        # Normalize (in case some modalities were missing)
        weight_sum <- sum(self$weights[names(all_probs)])
        weighted_probs <- weighted_probs / weight_sum

        data.frame(
          row_id = common_ids,
          prob = weighted_probs,
          response = factor(
            colnames(weighted_probs)[apply(weighted_probs, 1, which.max)],
            levels = private$.classes
          )
        )
      } else {
        # Average responses
        weighted_response <- Reduce(`+`, lapply(names(all_probs), function(m_id) {
          all_probs[[m_id]] * self$weights[m_id]
        }))

        weight_sum <- sum(self$weights[names(all_probs)])
        weighted_response <- weighted_response / weight_sum

        data.frame(row_id = common_ids, response = weighted_response)
      }
    },

    #' @description
    #' Print method
    print = function() {
      cat("=== OmicWeightedEnsemble ===\n")
      cat(sprintf("  Modalities: %d\n", length(self$modalities)))
      for (m in self$modalities) {
        w <- if (!is.null(self$weights)) sprintf(" (weight: %.3f)", self$weights[m$id]) else ""
        cat(sprintf("    - %s: %d features%s\n", m$id, m$n_features(), w))
      }
      cat(sprintf("  Fitted: %s\n", private$.fitted))
      invisible(self)
    }
  ),

  private = list(
    .fitted = FALSE,
    .fitted_base = NULL,
    .common_ids = NULL,
    .target_name = NULL,
    .classes = NULL
  )
)


#' @title Create Multi-Omics Stacked Ensemble (Convenience Function)
#'
#' @description
#' High-level interface for creating and optionally training a stacked ensemble.
#'
#' @param data_list Named list of data matrices (one per modality)
#' @param learner_list Optional named list of learners (or single learner for all)
#' @param y Optional target vector (if provided, trains immediately)
#' @param meta_learner Learner for stacking
#' @param resampling Resampling strategy
#' @param task_type "classif" or "regr"
#'
#' @return OmicStackedEnsemble object
#'
#' @examples
#' \dontrun{
#' # Create and train ensemble
#' ensemble <- stack_omics(
#'   data_list = list(mRNA = rna_matrix, miRNA = mirna_matrix),
#'   learner_list = list(mRNA = "classif.ranger", miRNA = "classif.glmnet"),
#'   y = y_vector,
#'   meta_learner = "classif.log_reg"
#' )
#'
#' # Predict
#' preds <- ensemble$predict(list(mRNA = new_rna, miRNA = new_mirna))
#' }
#'
#' @export
stack_omics <- function(data_list,
                         learner_list = NULL,
                         y = NULL,
                         meta_learner = NULL,
                         resampling = NULL,
                         task_type = "classif") {

  checkmate::assert_list(data_list, names = "named", min.len = 2)

  # Default learners
  if (is.null(learner_list)) {
    default_lrn <- if (task_type == "classif") "classif.ranger" else "regr.ranger"
    learner_list <- setNames(rep(list(default_lrn), length(data_list)), names(data_list))
  } else if (length(learner_list) == 1 && !is.list(learner_list)) {
    # Single learner for all modalities
    learner_list <- setNames(rep(list(learner_list), length(data_list)), names(data_list))
  }

  # Create modality specs
  modalities <- lapply(names(data_list), function(id) {
    OmicModalitySpec$new(
      id = id,
      x = data_list[[id]],
      learner = learner_list[[id]]
    )
  })

  # Create ensemble
  ensemble <- OmicStackedEnsemble$new(
    modalities = modalities,
    meta_learner = meta_learner,
    resampling = resampling,
    task_type = task_type
  )

  # Train if y provided
  if (!is.null(y)) {
    ensemble$train(y)
  }

  ensemble
}
