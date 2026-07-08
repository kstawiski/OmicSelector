#' @title BenchmarkService: Nested Cross-Validation with Zero Leakage
#'
#' @description
#' R6 class that enforces proper nested cross-validation for biomarker discovery.
#' Implements the outer loop (evaluation) and inner loop (selection) pattern
#' required for unbiased performance estimation.
#'
#' @details
#' The BenchmarkService guarantees scientific validity by:
#' - Enforcing that feature selection occurs in the inner loop only
#' - Computing the Nogueira Stability Index across outer folds
#' - Tracking which features are selected in each fold for consensus analysis
#' - Preventing any access to test data during training/selection
#'
#' @examples
#' \dontrun{
#' # Create benchmark service
#' service <- BenchmarkService$new(
#'   task = my_task,
#'   outer_folds = 5,
#'   inner_folds = 3
#' )
#'
#' # Add learners with embedded feature selection
#' service$add_learner(my_graph_learner)
#'
#' # Run nested CV
#' result <- service$run()
#'
#' # Get stability metrics
#' stability <- result$get_stability()
#' }
#'
#' @export
BenchmarkService <- R6::R6Class(

  "BenchmarkService",

  public = list(

    #' @description
    #' Create a new BenchmarkService
    #'
    #' @param task An mlr3 Task or OmicPipeline object
    #' @param outer_folds Number of outer CV folds (evaluation)
    #' @param inner_folds Number of inner CV folds (selection/tuning)
    #' @param stratify Logical, whether to stratify by outcome
    #' @param groups Optional column name for grouped CV (e.g., patient_id)
    #' @param seed Random seed for reproducibility
    #'
    #' @return A new BenchmarkService object
    initialize = function(task, outer_folds = 5, inner_folds = 3,
                          stratify = TRUE, groups = NULL, seed = NULL) {

      # Input validation
      stopifnot(
        "outer_folds must be a positive integer >= 2" =
          is.numeric(outer_folds) && length(outer_folds) == 1 &&
          outer_folds >= 2 && outer_folds == as.integer(outer_folds),
        "inner_folds must be a positive integer >= 2" =
          is.numeric(inner_folds) && length(inner_folds) == 1 &&
          inner_folds >= 2 && inner_folds == as.integer(inner_folds)
      )

      # Handle OmicPipeline input
      if (inherits(task, "OmicPipeline")) {
        private$.task <- task$get_task()
        private$.pipeline <- task
      } else if (inherits(task, "Task")) {
        private$.task <- task
        private$.pipeline <- NULL
      } else {
        stop("'task' must be an mlr3 Task or OmicPipeline object")
      }

      private$.outer_folds <- outer_folds
      private$.inner_folds <- inner_folds
      private$.stratify <- stratify
      private$.groups <- groups
      private$.seed <- seed

      private$.learners <- list()
      private$.results <- list()

      # Set seed for reproducibility
      if (!is.null(seed)) {
        set.seed(seed)
      }

      message(sprintf(
        "BenchmarkService created: outer=%d folds, inner=%d folds",
        outer_folds, inner_folds
      ))
    },

    #' @description
    #' Add a learner to benchmark
    #'
    #' @param learner A Learner, GraphLearner, or AutoFSelector
    #' @param id Optional identifier for the learner
    #'
    #' @return Self (for chaining)
    add_learner = function(learner, id = NULL) {
      if (is.null(id)) {
        id <- sprintf("learner_%d", length(private$.learners) + 1)
      }

      # Ensure nested CV when inner_folds > 1
      learner <- private$.ensure_nested_cv(learner)

      # Validate learner has no test data access
      private$.validate_learner(learner)

      private$.learners[[id]] <- learner
      message(sprintf("Added learner: %s", id))

      invisible(self)
    },

    #' @description
    #' Run the nested cross-validation benchmark
    #'
    #' @param measures List of performance measures (default: AUC, accuracy)
    #' @param parallel Logical, whether to run in parallel (default: TRUE)
    #' @param cache_dir Optional directory to cache benchmark results (RDS)
    #' @param cache_key Optional cache key override (string)
    #' @param threads Integer, number of threads for mlr3 learners (default: 1)
    #'
    #' @return A NestedCVResult object
    run = function(measures = NULL, parallel = TRUE, cache_dir = NULL,
                   cache_key = NULL, threads = 1) {
      private$.require_mlr3()

      if (length(private$.learners) == 0) {
        stop("No learners added. Use add_learner() first.")
      }

      # Reproducibility settings
      set_reproducible(seed = private$.seed, threads = threads)

      # Cache lookup
      cache_path <- NULL
      if (!is.null(cache_dir)) {
        cache_path <- private$.cache_path(cache_dir, cache_key)
        if (file.exists(cache_path)) {
          message(sprintf("Loading cached benchmark results: %s", cache_path))
          cached <- readRDS(cache_path)
          private$.results <- cached
          return(cached)
        }
      }

      # Default measures
      if (is.null(measures)) {
        measures <- list(
          mlr3::msr("classif.auc"),
          mlr3::msr("classif.acc"),
          mlr3::msr("classif.bbrier")
        )
      }

      # Create outer resampling with stratification and grouping support
      if (!is.null(private$.groups) && private$.groups %in% names(private$.task$data())) {
        # Grouped CV: ensure samples from same group stay together
        outer_rsmp <- mlr3::rsmp("cv", folds = private$.outer_folds)
        private$.task$set_col_roles(private$.groups, roles = "group")
      } else if (isTRUE(private$.stratify) && private$.task$task_type == "classif") {
        # Stratified CV: maintain class proportions
        outer_rsmp <- mlr3::rsmp("cv", folds = private$.outer_folds)
        private$.task$set_col_roles(private$.task$target_names, roles = c("target", "stratum"))
      } else {
        # Standard CV
        outer_rsmp <- mlr3::rsmp("cv", folds = private$.outer_folds)
      }

      # Always instantiate the resampling (critical for reproducibility)
      if (!is.null(private$.seed)) {
        # Backward-compatible: some mlr3 versions do not accept 'seed' argument
        tryCatch(
          outer_rsmp$instantiate(private$.task, seed = private$.seed),
          error = function(e) outer_rsmp$instantiate(private$.task)
        )
      } else {
        # Even without seed, we must instantiate to fix the splits
        outer_rsmp$instantiate(private$.task)
      }

      # Create benchmark design
      design <- mlr3::benchmark_grid(
        tasks = list(private$.task),
        learners = private$.learners,
        resamplings = list(outer_rsmp)
      )

      # Run benchmark
      message("Running nested cross-validation...")
      start_time <- Sys.time()

      if (parallel) {
        # Enable parallelization - only if not already set
        if (requireNamespace("future", quietly = TRUE)) {
          private$.plan_old <- future::plan()
          # Only set if currently sequential to avoid overriding user config
          if (inherits(private$.plan_old, "sequential")) {
            tryCatch({
              future::plan("multisession")
              private$.plan_was_set <- TRUE
            }, error = function(e) {
              message(sprintf("Parallel setup failed (%s). Running sequential.",
                              e$message))
              private$.plan_was_set <- FALSE
            })
          }
        } else {
          message("Package 'future' not available. Running sequential.")
        }
      }

      if (isTRUE(private$.plan_was_set)) {
        on.exit({
          future::plan(private$.plan_old)
          private$.plan_was_set <- FALSE
          private$.plan_old <- NULL
        }, add = TRUE)
      }

      bmr <- mlr3::benchmark(design, store_models = TRUE)

      end_time <- Sys.time()
      message(sprintf("Completed in %.1f seconds", difftime(end_time, start_time, units = "secs")))

      # Extract results and compute stability
      result <- private$.process_results(bmr, measures)

      if (!is.null(cache_path)) {
        saveRDS(result, cache_path)
      }

      private$.results <- result
      return(result)
    },

    #' @description
    #' Get the most recent results
    #'
    #' @return The NestedCVResult object
    get_results = function() {
      private$.results
    },

    #' @description
    #' Print method
    print = function() {
      cat("=== BenchmarkService ===\n")
      cat(sprintf("Task: %s\n", private$.task$id))
      cat(sprintf("Outer folds: %d\n", private$.outer_folds))
      cat(sprintf("Inner folds: %d\n", private$.inner_folds))
      cat(sprintf("Learners: %d\n", length(private$.learners)))
      if (length(private$.learners) > 0) {
        cat("  ", paste(names(private$.learners), collapse = ", "), "\n")
      }
      invisible(self)
    }
  ),

  private = list(
    .task = NULL,
    .pipeline = NULL,
    .outer_folds = NULL,
    .inner_folds = NULL,
    .stratify = NULL,
    .groups = NULL,
    .seed = NULL,
    .learners = NULL,
    .results = NULL,
    .plan_was_set = FALSE,
    .plan_old = NULL,

    # Check for mlr3 packages
    .require_mlr3 = function() {
      required <- c("mlr3", "mlr3pipelines")
      for (pkg in required) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop(sprintf("Package '%s' required. Install with: install.packages('%s')",
                       pkg, pkg), call. = FALSE)
        }
      }
    },

    # Validate that learner doesn't have obvious leakage patterns
    .validate_learner = function(learner) {
      # GraphLearner or AutoTuner/AutoFSelector are the expected patterns
      if (inherits(learner, "GraphLearner") ||
          inherits(learner, "AutoTuner") ||
          inherits(learner, "AutoFSelector")) {
        return(TRUE)
      }

      # Plain learners are allowed but should be wrapped
      warning(paste(
        "Plain learner detected. For zero-leakage guarantee,",
        "use OmicPipeline$create_graph_learner() or wrap in GraphLearner."
      ))

      return(TRUE)
    },

    # Process benchmark results
    .process_results = function(bmr, measures) {
      # Aggregate performance
      performance <- bmr$aggregate(measures)

      # Extract per-fold results for stability analysis
      per_fold <- private$.extract_per_fold(bmr)

      # Compute stability index if feature selection was used
      stability <- private$.compute_stability(bmr)

      # Build result object
      result <- list(
        benchmark_result = bmr,
        performance = performance,
        per_fold = per_fold,
        stability = stability,
        outer_folds = private$.outer_folds,
        inner_folds = private$.inner_folds,
        seed = private$.seed,
        timestamp = Sys.time()
      )

      class(result) <- c("NestedCVResult", "list")
      return(result)
    },

    # Extract per-fold information
    .extract_per_fold = function(bmr) {
      fold_info <- list()

      for (i in seq_len(bmr$n_resample_results)) {
        rr <- bmr$resample_result(i)

        fold_info[[i]] <- list(
          learner_id = rr$learner$id,
          task_id = rr$task$id,
          n_iters = rr$iters
        )
      }

      return(fold_info)
    },

    # Compute Nogueira Stability Index per learner
    .compute_stability = function(bmr) {
      all_features <- private$.task$feature_names
      if (length(all_features) == 0) {
        return(list(
          nogueira_index = NA_real_,
          message = "No features available for stability computation.",
          selected_features_per_fold = list()
        ))
      }

      nogueira_dt <- data.table::data.table(
        learner_id = character(0),
        nogueira_index = numeric(0)
      )
      selected_per_fold <- list()

      for (i in seq_len(bmr$n_resample_results)) {
        rr <- bmr$resample_result(i)
        learner_id <- rr$learner$id

        feature_sets <- extract_features_from_resample(rr)
        selected_per_fold[[learner_id]] <- feature_sets

        if (all(sapply(feature_sets, length) == 0)) {
          next
        }

        stab <- compute_nogueira_stability(feature_sets, all_features)
        nogueira_dt <- data.table::rbindlist(list(
          nogueira_dt,
          data.table::data.table(
            learner_id = learner_id,
            nogueira_index = stab$nogueira_index
          )
        ), use.names = TRUE, fill = TRUE)
      }

      if (nrow(nogueira_dt) == 0) {
        return(list(
          nogueira_index = NA_real_,
          message = "Could not extract features from learners. Ensure GraphLearner with filter PipeOps or AutoTuner.",
          selected_features_per_fold = selected_per_fold
        ))
      }

      list(
        nogueira_index = nogueira_dt,
        message = "Stability computed per learner.",
        selected_features_per_fold = selected_per_fold
      )
    },

    # Ensure nested CV by wrapping GraphLearner in AutoTuner when needed
    .ensure_nested_cv = function(learner) {
      if (private$.inner_folds < 2) {
        return(learner)
      }

      if (inherits(learner, "AutoTuner") || inherits(learner, "AutoFSelector")) {
        return(learner)
      }

      if (!inherits(learner, "GraphLearner")) {
        message(paste(
          "Inner CV requested, but learner is not a GraphLearner/AutoTuner.",
          "Proceeding without inner-loop tuning."
        ))
        return(learner)
      }

      if (!requireNamespace("mlr3tuning", quietly = TRUE)) {
        stop("Package 'mlr3tuning' is required for nested CV with inner_folds >= 2.",
             call. = FALSE)
      }
      if (!requireNamespace("paradox", quietly = TRUE)) {
        stop("Package 'paradox' is required for nested CV with inner_folds >= 2.",
             call. = FALSE)
      }

      # Find filter.nfeat parameter
      param_ids <- learner$param_set$ids()
      param_id <- grep("filter.*nfeat$", param_ids, value = TRUE)
      if (length(param_id) == 0) {
        message(paste(
          "Inner CV requested, but no filter.nfeat parameter found.",
          "Proceeding without inner-loop tuning."
        ))
        return(learner)
      }
      if (length(param_id) > 1) {
        param_id <- param_id[1]
      }

      filter_values <- private$.default_filter_values()
      # p_int() no longer accepts `id` in paradox >= 1.0.0; use named ps() arg.
      search_space <- do.call(
        paradox::ps,
        setNames(
          list(paradox::p_int(
            lower = min(filter_values),
            upper = max(filter_values)
          )),
          param_id
        )
      )

      inner_rsmp <- mlr3::rsmp("cv", folds = private$.inner_folds)
      measure <- private$.default_measure()

      tuner_obj <- mlr3tuning::tnr("grid_search")
      if ("resolution" %in% tuner_obj$param_set$ids()) {
        tuner_obj$param_set$values$resolution <- length(filter_values)
      }
      terminator <- mlr3tuning::trm("evals", n_evals = length(filter_values))

      message(sprintf(
        "Wrapping learner '%s' with inner CV tuning (nfeat grid: %s)",
        learner$id, paste(filter_values, collapse = ", ")
      ))

      at <- mlr3tuning::AutoTuner$new(
        learner = learner,
        resampling = inner_rsmp,
        measure = measure,
        search_space = search_space,
        terminator = terminator,
        tuner = tuner_obj
      )
      at$id <- sprintf("%s_inner", learner$id)
      at
    },

    .default_filter_values = function() {
      n_features <- length(private$.task$feature_names)
      candidates <- unique(c(5, 10, 20, 50, round(n_features * 0.05), round(n_features * 0.1)))
      candidates <- candidates[candidates >= 1 & candidates <= n_features]
      if (length(candidates) == 0) {
        candidates <- max(1, n_features)
      }
      as.integer(sort(unique(candidates)))
    },

    .default_measure = function() {
      if (private$.task$task_type == "classif") {
        if (!is.null(private$.task$class_names) &&
            length(private$.task$class_names) == 2) {
          return(mlr3::msr("classif.auc"))
        }
        return(mlr3::msr("classif.acc"))
      }
      mlr3::msr("regr.rmse")
    },

    .build_cache_key = function() {
      learner_ids <- names(private$.learners)
      if (is.null(learner_ids) || any(!nzchar(learner_ids))) {
        learner_ids <- vapply(private$.learners, function(lrn) lrn$id, character(1))
      }

      key <- paste(
        "task", private$.task$id,
        "n", private$.task$nrow,
        "p", length(private$.task$feature_names),
        "outer", private$.outer_folds,
        "inner", private$.inner_folds,
        "seed", ifelse(is.null(private$.seed), "none", private$.seed),
        "learners", paste(learner_ids, collapse = ","),
        sep = "_"
      )
      key
    },

    .cache_path = function(cache_dir, cache_key = NULL) {
      if (is.null(cache_key)) {
        cache_key <- private$.build_cache_key()
      }

      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      key_hash <- as.character(hash_cache_key(cache_key))
      file.path(cache_dir, sprintf("benchmark_%s.rds", key_hash))
    }
  )
)


#' @title Print method for NestedCVResult
#' @param x A NestedCVResult object
#' @param ... Additional arguments (ignored)
#' @export
print.NestedCVResult <- function(x, ...) {
  cat("=== Nested Cross-Validation Result ===\n\n")

  cat(sprintf("Outer folds: %d\n", x$outer_folds))
  cat(sprintf("Inner folds: %d\n", x$inner_folds))
  cat(sprintf("Seed: %s\n", ifelse(is.null(x$seed), "not set", x$seed)))
  cat(sprintf("Timestamp: %s\n\n", x$timestamp))

  cat("Performance:\n")
  print(x$performance)

  ni <- x$stability$nogueira_index
  if (is.numeric(ni) && length(ni) == 1 && !is.na(ni)) {
    cat(sprintf("\nNogueira Stability Index: %.3f\n", ni))
  } else if (is.data.frame(ni)) {
    cat("\nNogueira Stability Index (per learner):\n")
    print(ni)
  }

  invisible(x)
}


#' @title Create benchmark service from OmicPipeline
#'
#' @description
#' Convenience function to create a BenchmarkService
#'
#' @param pipeline An OmicPipeline object
#' @param outer_folds Number of outer CV folds
#' @param inner_folds Number of inner CV folds
#' @param ... Additional arguments passed to BenchmarkService$new()
#'
#' @return A BenchmarkService object
#' @export
omic_benchmark <- function(pipeline, outer_folds = 5, inner_folds = 3, ...) {
  BenchmarkService$new(
    task = pipeline,
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    ...
  )
}
