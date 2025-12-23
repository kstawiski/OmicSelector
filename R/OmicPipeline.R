#' @title OmicPipeline: Zero-Leakage Feature Selection Pipeline
#'
#' @description
#' R6 class that encapsulates the mlr3 pipeline for biomarker discovery.
#' Guarantees zero data leakage by enforcing all preprocessing, feature selection,
#' and model training within proper cross-validation folds.
#'
#' @importFrom mlr3pipelines `%>>%`
#'
#' @details
#' OmicPipeline is the central class for OmicSelector 2.0. It replaces the legacy
#' script-based approach with a rigorous, composable, and reproducible architecture.
#'
#' Key features:
#' - All preprocessing (imputation, scaling) occurs inside CV folds
#' - Feature selection is embedded in the inner loop of nested CV
#' - Oversampling (SMOTE/ROSE) is applied only to training data per fold
#' - Factory methods generate configured GraphLearners
#'
#' @examples
#' \dontrun{
#' # Create pipeline from data
#' pipeline <- OmicPipeline$new(
#'   data = my_data,
#'   target = "outcome",
#'   positive = "Case"
#' )
#'
#' # Create a graph learner with feature selection
#' learner <- pipeline$create_graph_learner(
#'   filter = "anova",
#'   model = "ranger",
#'   n_features = 20
#' )
#'
#' # Run nested cross-validation
#' result <- pipeline$benchmark(learner, outer_folds = 5, inner_folds = 3)
#' }
#'
#' @export
OmicPipeline <- R6::R6Class(

"OmicPipeline",

  public = list(

    #' @description
    #' Create a new OmicPipeline object
    #'
    #' @param data Either a data.frame or a named list of data.frames for multi-omics.
    #'   For multi-omics, use named list: list(rna = rna_data, mirna = mirna_data).
    #'   Features will be namespaced: rna::gene1, mirna::hsa-miR-21.
    #' @param target Name of the target column
    #' @param positive Positive class label (for binary classification)
    #' @param patient_id Optional column name for patient grouping (prevents leakage)
    #' @param batch Optional column name for batch information
    #' @param id Unique identifier for this pipeline
    #'
    #' @return A new OmicPipeline object
    initialize = function(data, target, positive = NULL, patient_id = NULL,
                          batch = NULL, id = "omic_task") {

      # Validate and normalize multi-omics input
      private$.omics_input <- validate_omics_input(data, target)

      # Store metadata
      private$.target <- target
      private$.positive <- positive
      private$.patient_id <- patient_id
      private$.batch <- batch
      private$.id <- id
      private$.is_multi_omics <- private$.omics_input$is_multi_omics

      # Merge modalities into single data.frame for mlr3 Task
      merged_data <- merge_omics_data(private$.omics_input)

      # Re-add target column
      merged_data[[target]] <- private$.omics_input$target_data

      # Identify feature columns (exclude metadata columns)
      # Note: c() automatically drops NULL values, so we just need to filter empty strings
      meta_cols <- c(target, patient_id, batch)
      meta_cols <- meta_cols[!is.na(meta_cols) & nzchar(meta_cols)]
      private$.feature_names <- setdiff(names(merged_data), meta_cols)

      # Create mlr3 Task
      private$.create_task(merged_data, target, positive)

      # Log creation
      n_modalities <- length(private$.omics_input$modalities)
      if (private$.is_multi_omics) {
        message(sprintf(
          "OmicPipeline created: %d samples, %d features (%d modalities), target='%s'",
          nrow(merged_data), length(private$.feature_names), n_modalities, target
        ))
      } else {
        message(sprintf(
          "OmicPipeline created: %d samples, %d features, target='%s'",
          nrow(merged_data), length(private$.feature_names), target
        ))
      }
    },

    #' @description
    #' Create a GraphLearner with proper leakage prevention
    #'
    #' @param filter Filter method name (e.g., "anova", "mrmr", "correlation")
    #' @param model Model type (e.g., "ranger", "glmnet", "svm", "mlp",
    #'   "tabtransformer") or an mlr3 Learner object.
    #' @param n_features Number of features to select (or proportion if < 1)
    #' @param impute_method Imputation method ("median", "mean", "sample")
    #' @param scale Logical, whether to scale features
    #' @param oversample Oversampling method (NULL, "smote", "rose")
    #' @param screening Logical, whether to apply a fast variance-based screening
    #'   filter before the main selection step (default: FALSE)
    #' @param screening_nfeat Integer, number of features to keep during screening.
    #'   If NULL, uses screening_frac.
    #' @param screening_frac Numeric (0,1], fraction of features to keep during
    #'   screening when screening_nfeat is NULL (default: 0.2)
    #' @param batch_correct Logical or character. If TRUE, adds FrozenComBat batch
    #'   correction using the batch column specified in pipeline creation. If a
    #'   character string, uses that as the batch column name. Default: FALSE.
    #' @param autoencoder Logical or list. If TRUE, inserts a torch autoencoder
    #'   PipeOp after scaling. If a list, its contents are passed to
    #'   \code{create_autoencoder_pipeop()} to configure latent dimension,
    #'   training, and transfer options. Default: NULL (disabled).
    #'
    #' @return A mlr3 GraphLearner object
    create_graph_learner = function(filter = "anova", model = "ranger",
                                     n_features = 20, impute_method = "median",
                                     scale = TRUE, oversample = NULL,
                                     screening = FALSE, screening_nfeat = NULL,
                                     screening_frac = 0.2,
                                     batch_correct = FALSE,
                                     autoencoder = NULL) {
      private$.require_mlr3()

      # Build pipeline components
      ops <- list()

      # 1. Imputation (always first)
      ops$impute <- private$.create_impute_op(impute_method)

      # 2. Batch correction (AFTER imputation, BEFORE scaling)
      # Must be early to correct systematic batch effects before other transforms
      if (!isFALSE(batch_correct)) {
        batch_col <- if (isTRUE(batch_correct)) private$.batch else batch_correct
        if (!is.null(batch_col)) {
          ops$batch <- create_frozen_combat_pipeop(batch_col = batch_col)
          message(sprintf("Adding FrozenComBat batch correction (column: %s)", batch_col))
        } else {
          warning("batch_correct=TRUE but no batch column specified. Skipping batch correction.")
        }
      }

      # 3. Screening (variance-based, optional, before scaling/selection)
      if (isTRUE(screening)) {
        ops$screen <- private$.create_screening_op(
          screening_nfeat = screening_nfeat,
          screening_frac = screening_frac
        )
      }

      # 4. Scaling (optional, before selection)
      if (scale) {
        ops$scale <- mlr3pipelines::po("scale")
      }

      # 5. Autoencoder (optional, after scaling)
      if (isTRUE(autoencoder) || is.list(autoencoder)) {
        ae_params <- if (is.list(autoencoder)) autoencoder else list()
        ops$autoencoder <- do.call(create_autoencoder_pipeop, ae_params)
      }

      # 6. Oversampling (AFTER scaling/autoencoder, BEFORE selection)
      # RATIONALE: SMOTE uses k-NN which requires scaled features for proper
      # distance calculations. If autoencoder is enabled, SMOTE operates on
      # latent features to reduce sparsity and noise.
      # Order: Impute → Batch → Screen → Scale → Autoencoder → SMOTE → Filter → Learner
      if (!is.null(oversample)) {
        ops$oversample <- private$.create_oversample_op(oversample)
      }

      # 7. Feature Selection
      ops$filter <- private$.create_filter_op(filter, n_features)

      # 8. Learner
      ops$learner <- private$.create_learner(model)

      # Compose the graph
      graph <- private$.compose_graph(ops)

      # Wrap as GraphLearner
      graph_learner <- mlr3pipelines::GraphLearner$new(graph)
      model_id <- if (inherits(model, "Learner")) model$id else model
      graph_learner$id <- sprintf("omic_%s_%s", filter, model_id)

      return(graph_learner)
    },

    #' @description
    #' Create an AutoTuner that tunes filter.nfeat in the inner CV loop
    #'
    #' @param learner A GraphLearner produced by create_graph_learner()
    #' @param filter_values Integer vector of n_features values to try
    #' @param inner_resampling Inner resampling strategy (default: 3-fold CV)
    #' @param measure Performance measure (default: classif.auc)
    #' @param tuner Tuning strategy (default: "grid_search")
    #'
    #' @return An mlr3tuning::AutoTuner object
    create_auto_fselector = function(learner, filter_values = c(5, 10, 20, 50),
                                      inner_resampling = NULL, measure = NULL,
                                      tuner = "grid_search") {
      private$.require_mlr3()

      if (!requireNamespace("mlr3tuning", quietly = TRUE)) {
        stop("Package 'mlr3tuning' is required for create_auto_fselector().",
             call. = FALSE)
      }
      if (!requireNamespace("paradox", quietly = TRUE)) {
        stop("Package 'paradox' is required for create_auto_fselector().",
             call. = FALSE)
      }

      if (!inherits(learner, "Learner")) {
        stop("learner must be an mlr3 Learner or GraphLearner", call. = FALSE)
      }

      filter_values <- unique(as.integer(filter_values))
      filter_values <- filter_values[filter_values > 0]
      if (length(filter_values) == 0) {
        stop("filter_values must contain positive integers.", call. = FALSE)
      }

      # Default inner resampling: 3-fold CV
      if (is.null(inner_resampling)) {
        inner_resampling <- mlr3::rsmp("cv", folds = 3)
      }

      # Default measure: AUC for classification
      if (is.null(measure)) {
        measure <- mlr3::msr("classif.auc")
      }

      # Find filter.nfeat parameter id
      param_ids <- learner$param_set$ids()
      param_id <- grep("filter.*nfeat$", param_ids, value = TRUE)
      if (length(param_id) == 0) {
        stop("No filter.nfeat parameter found in learner. ",
             "Ensure the GraphLearner includes a filter PipeOp.",
             call. = FALSE)
      }
      if (length(param_id) > 1) {
        param_id <- param_id[1]
      }

      # paradox::p_int() in newer versions does not accept an `id` argument.
      # The parameter id is provided by naming the argument in paradox::ps().
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

      terminator <- mlr3tuning::trm("evals", n_evals = length(filter_values))
      tuner_obj <- mlr3tuning::tnr(tuner)
      if ("resolution" %in% tuner_obj$param_set$ids()) {
        tuner_obj$param_set$values$resolution <- length(filter_values)
      }

      mlr3tuning::AutoTuner$new(
        learner = learner,
        resampling = inner_resampling,
        measure = measure,
        search_space = search_space,
        terminator = terminator,
        tuner = tuner_obj
      )
    },

    #' @description
    #' Run benchmark with proper nested cross-validation
    #'
    #' @param learners List of learners to benchmark
    #' @param outer_folds Number of outer CV folds
    #' @param inner_folds Number of inner CV folds
    #' @param stratify Logical, whether to stratify by outcome
    #' @param seed Random seed for reproducibility
    #' @param parallel Logical, whether to run in parallel (default: TRUE)
    #' @param threads Integer, number of threads for mlr3 learners (default: 1)
    #' @param cache_dir Optional directory to cache benchmark results (RDS)
    #' @param cache_key Optional cache key override (string)
    #' @param screening Logical, apply variance screening before nested CV
    #' @param screening_nfeat Integer, number of features to keep in screening
    #' @param screening_frac Numeric (0,1], fraction of features to keep in screening
    #'
    #' @return A BenchmarkResult object
    benchmark = function(learners, outer_folds = 5, inner_folds = 3,
                         stratify = TRUE, seed = NULL, parallel = TRUE,
                         threads = 1, cache_dir = NULL, cache_key = NULL,
                         screening = FALSE, screening_nfeat = NULL,
                         screening_frac = 0.2) {
      private$.require_mlr3()

      set_reproducible(seed = seed, threads = threads)

      # Ensure learners is a list
      if (!is.list(learners)) {
        learners <- list(learners)
      }

      if (isTRUE(screening)) {
        learners <- lapply(learners, function(learner) {
          private$.apply_screening(
            learner = learner,
            screening_nfeat = screening_nfeat,
            screening_frac = screening_frac
          )
        })
      }

      # Use BenchmarkService for nested CV
      service <- BenchmarkService$new(
        task = self,
        outer_folds = outer_folds,
        inner_folds = inner_folds,
        stratify = stratify,
        seed = seed
      )

      learner_names <- names(learners)
      for (i in seq_along(learners)) {
        if (!is.null(learner_names) && nzchar(learner_names[i])) {
          service$add_learner(learners[[i]], id = learner_names[i])
        } else {
          service$add_learner(learners[[i]])
        }
      }

      service$run(
        parallel = parallel,
        cache_dir = cache_dir,
        cache_key = cache_key,
        threads = threads
      )
    },

    #' @description
    #' Fit a learner on the full dataset and return the trained model
    #' and selected features (for deployment).
    #'
    #' @param learner A Learner or GraphLearner
    #' @param seed Random seed for reproducibility
    #' @param threads Integer, number of threads for mlr3 learners (default: 1)
    #' @param screening Logical, apply variance screening before selection
    #' @param screening_nfeat Integer, number of features to keep in screening
    #' @param screening_frac Numeric (0,1], fraction of features to keep in screening
    #'
    #' @return A list with trained learner and selected features
    fit = function(learner, seed = NULL, threads = 1,
                   screening = FALSE, screening_nfeat = NULL,
                   screening_frac = 0.2) {
      private$.require_mlr3()

      set_reproducible(seed = seed, threads = threads)

      learner <- private$.apply_screening(
        learner = learner,
        screening_nfeat = screening_nfeat,
        screening_frac = screening_frac,
        screening = screening
      )

      learner$train(private$.task)
      selected <- extract_selected_features(learner)

      result <- list(
        learner = learner,
        selected_features = selected,
        task_id = private$.task$id,
        n_features = length(private$.feature_names)
      )
      class(result) <- c("OmicFit", "list")
      result
    },

    #' @description
    #' Get the underlying mlr3 Task
    #'
    #' @return The mlr3 Task object
    get_task = function() {
      private$.task
    },

    #' @description
    #' Get feature names
    #'
    #' @return Character vector of feature names (namespaced for multi-omics)
    get_feature_names = function() {
      private$.feature_names
    },

    #' @description
    #' Check if this is a multi-omics pipeline
    #'
    #' @return Logical
    is_multi_omics = function() {
      private$.is_multi_omics
    },

    #' @description
    #' Get modality information for multi-omics data
    #'
    #' @return A data.frame with modality details, or NULL for single-modality
    get_modality_info = function() {
      get_modality_info(private$.omics_input)
    },

    #' @description
    #' Get features for a specific modality
    #'
    #' @param modality Name of the modality (e.g., "rna", "mirna")
    #' @return Character vector of feature names for that modality
    get_modality_features = function(modality) {
      if (!modality %in% names(private$.omics_input$modalities)) {
        stop(sprintf("Modality '%s' not found. Available: %s",
                     modality,
                     paste(names(private$.omics_input$modalities), collapse = ", ")),
             call. = FALSE)
      }
      private$.omics_input$feature_map[[modality]]
    },

    #' @description
    #' Print method
    print = function() {
      cat("=== OmicPipeline ===\n")
      cat(sprintf("ID: %s\n", private$.id))
      cat(sprintf("Samples: %d\n", private$.task$nrow))
      cat(sprintf("Features: %d\n", length(private$.feature_names)))
      cat(sprintf("Target: %s\n", private$.target))
      if (!is.null(private$.positive)) {
        cat(sprintf("Positive class: %s\n", private$.positive))
      }
      if (!is.null(private$.patient_id)) {
        cat(sprintf("Patient grouping: %s\n", private$.patient_id))
      }

      # Multi-omics info
      if (private$.is_multi_omics) {
        cat("\nModalities:\n")
        info <- self$get_modality_info()
        for (i in seq_len(nrow(info))) {
          cat(sprintf("  %s: %d features%s\n",
                      info$modality[i], info$n_features[i],
                      if (info$has_target[i]) " (target)" else ""))
        }
      }

      invisible(self)
    }
  ),

  private = list(
    .task = NULL,
    .target = NULL,
    .positive = NULL,
    .patient_id = NULL,
    .batch = NULL,
    .id = NULL,
    .feature_names = NULL,
    .omics_input = NULL,
    .is_multi_omics = FALSE,

    # Check if mlr3 packages are available
    .require_mlr3 = function() {
      required <- c("mlr3", "mlr3pipelines", "mlr3learners")
      for (pkg in required) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop(sprintf("Package '%s' required but not installed. Install with: install.packages('%s')",
                       pkg, pkg), call. = FALSE)
        }
      }
    },

    # Create the mlr3 Task
    .create_task = function(data, target, positive) {
      private$.require_mlr3()

      # Determine task type based on target
      target_values <- data[[target]]

      if (is.factor(target_values) || is.character(target_values)) {
        # Classification task
        if (!is.factor(target_values)) {
          data[[target]] <- as.factor(target_values)
        }

        private$.task <- mlr3::TaskClassif$new(
          id = private$.id,
          backend = data[, c(private$.feature_names, target)],
          target = target
        )

        # Set positive class if specified
        if (!is.null(positive)) {
          private$.task$positive <- positive
        }
      } else if (is.numeric(target_values)) {
        # Regression task
        private$.task <- mlr3::TaskRegr$new(
          id = private$.id,
          backend = data[, c(private$.feature_names, target)],
          target = target
        )
      } else {
        stop("Target must be factor/character (classification) or numeric (regression)")
      }
    },

    # Create imputation PipeOp
    .create_impute_op = function(method) {
      if (method == "median") {
        mlr3pipelines::po("imputemedian")
      } else if (method == "mean") {
        mlr3pipelines::po("imputemean")
      } else if (method == "sample") {
        mlr3pipelines::po("imputesample")
      } else {
        stop(sprintf("Unknown imputation method: %s", method))
      }
    },

    # Create oversampling PipeOp
    .create_oversample_op = function(method) {
      if (method == "smote") {
        # Note: Actual implementation requires mlr3pipelines::po("smote")
        # which may need additional packages
        if (requireNamespace("themis", quietly = TRUE)) {
          mlr3pipelines::po("smote")
        } else {
          warning("SMOTE requires 'themis' package. Skipping oversampling.")
          return(NULL)
        }
      } else if (method == "rose") {
        warning("ROSE oversampling not yet implemented. Skipping.")
        return(NULL)
      } else {
        stop(sprintf("Unknown oversampling method: %s", method))
      }
    },

    # Create screening (variance) PipeOp
    .create_screening_op = function(screening_nfeat = NULL, screening_frac = 0.2) {
      n_features <- length(private$.feature_names)

      if (n_features == 0) {
        return(NULL)
      }

      if (is.null(screening_nfeat)) {
        if (!is.numeric(screening_frac) || length(screening_frac) != 1 ||
            screening_frac <= 0 || screening_frac > 1) {
          stop("screening_frac must be in (0, 1].", call. = FALSE)
        }
        screening_nfeat <- as.integer(round(n_features * screening_frac))
      }

      screening_nfeat <- as.integer(screening_nfeat)
      if (is.na(screening_nfeat) || screening_nfeat <= 0) {
        stop("screening_nfeat must be a positive integer.", call. = FALSE)
      }

      if (screening_nfeat >= n_features) {
        return(NULL)
      }

      mlr3pipelines::po(
        "filter",
        id = "screen",
        filter = mlr3filters::flt("variance"),
        filter.nfeat = screening_nfeat
      )
    },

    # Create filter PipeOp
    .create_filter_op = function(filter, n_features) {
      # Map filter names to mlr3filters
      filter_map <- list(
        # Univariate statistical tests
        "anova" = "anova",                       # ANOVA F-test (default)
        "kruskal" = "kruskal_test",              # Non-parametric alternative to ANOVA
        "chi_squared" = "chi_squared",           # Chi-squared test for discrete features

        # Variance-based
        "variance" = "variance",                 # Remove low-variance features

        # Correlation-based
        "correlation" = "correlation",           # Correlation with target

        # Information-theoretic
        "information_gain" = "information_gain", # Entropy-based importance
        "gain_ratio" = "gain_ratio",             # Normalized information gain
        "mrmr" = "mrmr",                         # Minimum Redundancy Maximum Relevance
        "cmim" = "cmim",                         # Conditional Mutual Information Maximization
        "jmim" = "jmim",                         # Joint Mutual Information Maximization
        "jmi" = "jmi",                           # Joint Mutual Information

        # Model-based importance
        "auc" = "auc",                           # AUC of univariate models
        "relief" = "relief",                     # Relief algorithm
        "importance" = "importance",             # Random Forest importance
        "permutation" = "permutation"            # Permutation importance
      )

      if (!filter %in% names(filter_map)) {
        available <- paste(names(filter_map), collapse = ", ")
        stop(sprintf("Unknown filter '%s'. Available: %s", filter, available))
      }

      # Create the filter PipeOp
      mlr3pipelines::po("filter",
                        filter = mlr3filters::flt(filter_map[[filter]]),
                        filter.nfeat = n_features)
    },

    # Apply screening to an existing GraphLearner
    .apply_screening = function(learner, screening_nfeat = NULL,
                                screening_frac = 0.2, screening = TRUE) {
      if (!isTRUE(screening)) {
        return(learner)
      }

      if (!inherits(learner, "GraphLearner")) {
        warning("Screening requested, but learner is not a GraphLearner. Skipping.")
        return(learner)
      }

      screen_op <- private$.create_screening_op(
        screening_nfeat = screening_nfeat,
        screening_frac = screening_frac
      )

      if (is.null(screen_op)) {
        return(learner)
      }

      graph <- screen_op %>>% learner$graph
      screened <- mlr3pipelines::GraphLearner$new(graph)
      screened$id <- sprintf("%s_screen", learner$id)
      screened
    },

    # Create learner
    .create_learner = function(model) {
      if (inherits(model, "Learner")) {
        return(model)
      }

      # Supported learners - map shorthand names to mlr3 learner IDs
      learner_map <- list(
        # Tree-based ensemble methods
        "ranger" = "classif.ranger",           # Random Forest (fast implementation)
        "xgboost" = "classif.xgboost",         # XGBoost gradient boosting
        "lightgbm" = "classif.lightgbm",       # LightGBM gradient boosting
        "rpart" = "classif.rpart",             # Single decision tree

        # Linear models
        "glmnet" = "classif.glmnet",           # Elastic net (L1/L2 regularization)
        "log_reg" = "classif.log_reg",         # Logistic regression

        # Distance-based methods
        "svm" = "classif.svm",                 # Support Vector Machine (RBF kernel)
        "kknn" = "classif.kknn",               # K-nearest neighbors

        # Probabilistic classifiers
        "naive_bayes" = "classif.naive_bayes", # Naive Bayes

        # Discriminant analysis
        "lda" = "classif.lda",                 # Linear Discriminant Analysis
        "qda" = "classif.qda",                 # Quadratic Discriminant Analysis

        # Neural networks
        "nnet" = "classif.nnet"                # Single-layer neural network
      )

      if (!model %in% names(learner_map)) {
        # Deep learning shortcuts (mlr3torch)
        if (identical(model, "mlp")) {
          return(create_mlp_learner())
        }
        if (identical(model, "tabtransformer")) {
          return(make_tabtransformer_learner())
        }

        available <- paste(c(names(learner_map), "mlp", "tabtransformer"), collapse = ", ")
        stop(sprintf("Unknown model '%s'. Available: %s", model, available))
      }

      mlr3::lrn(learner_map[[model]], predict_type = "prob")
    },

    # Compose graph from operations
    .compose_graph = function(ops) {
      # Remove NULL operations
      ops <- ops[!sapply(ops, is.null)]

      if (length(ops) < 2) {
        stop("Need at least imputation + learner")
      }

      # Chain operations using %>>%
      graph <- ops[[1]]
      for (i in 2:length(ops)) {
        graph <- graph %>>% ops[[i]]
      }

      return(graph)
    },

    # Wrap benchmark result with additional analysis
    .wrap_benchmark_result = function(bmr) {
      result <- list(
        benchmark_result = bmr,
        timestamp = Sys.time(),
        pipeline_id = private$.id,
        n_features = length(private$.feature_names),
        class = "OmicBenchmarkResult"
      )

      class(result) <- c("OmicBenchmarkResult", "list")
      return(result)
    }
  )
)


#' @title Quick pipeline creation from data
#'
#' @description
#' Convenience function to create an OmicPipeline from a data.frame
#'
#' @param data A data.frame containing features and target
#' @param target Name of the target column
#' @param positive Positive class label (for binary classification)
#' @param ... Additional arguments passed to OmicPipeline$new()
#'
#' @return An OmicPipeline object
#' @export
omic_pipeline <- function(data, target, positive = NULL, ...) {
  OmicPipeline$new(data = data, target = target, positive = positive, ...)
}


#' @title Print method for OmicBenchmarkResult
#' @param x An OmicBenchmarkResult object
#' @param ... Additional arguments (ignored)
#' @export
print.OmicBenchmarkResult <- function(x, ...) {
  cat("=== OmicSelector Benchmark Result ===\n")
  cat(sprintf("Pipeline: %s\n", x$pipeline_id))
  cat(sprintf("Features: %d\n", x$n_features))
  cat(sprintf("Timestamp: %s\n", x$timestamp))
  cat("\nPerformance:\n")
  if (!is.null(x$benchmark_result)) {
    print(x$benchmark_result$aggregate())
  }
  invisible(x)
}

#' @title Print method for OmicFit
#' @param x An OmicFit object
#' @param ... Additional arguments (ignored)
#' @export
print.OmicFit <- function(x, ...) {
  cat("=== OmicSelector Fit ===\n")
  cat(sprintf("Task: %s\n", x$task_id))
  cat(sprintf("Selected features: %d\n",
              ifelse(is.null(x$selected_features), 0L, length(x$selected_features))))
  if (!is.null(x$learner)) {
    cat(sprintf("Learner: %s\n", x$learner$id))
  }
  invisible(x)
}
