#' @title Bayesian Hyperparameter Optimization for Omics
#'
#' @description
#' Provides Bayesian optimization for hyperparameter tuning using mlr3mbo.
#' Optimized for high-dimensional omics data with sensible default search spaces
#' for XGBoost, LightGBM, Random Forest, and glmnet.
#'
#' @details
#' Key features:
#' - Omics-optimized search spaces (strong regularization, shallow trees)
#' - Automatic budget scaling based on sample size
#' - RF surrogate model (robust for mixed integer/continuous params)
#' - Expected Improvement (EI) acquisition function
#'
#' @references
#' Bischl et al. (2023). mlr3mbo: Bayesian Optimization in mlr3.
#'
#' @name BayesianTuner
NULL


#' @title Create Bayesian-Optimized AutoTuner for XGBoost
#'
#' @description
#' Creates an AutoTuner with XGBoost learner and omics-optimized search space.
#'
#' @param task An mlr3 Task (used to set data-aware defaults)
#' @param n_evals Number of Bayesian optimization iterations (default: auto based on n)
#' @param inner_folds Number of inner CV folds for tuning (default: 3)
#' @param measure Performance measure (default: classif.auc)
#' @param early_stopping Use early stopping with nrounds (default: TRUE)
#' @param nthread Number of threads for XGBoost (default: 4)
#'
#' @return An AutoTuner object ready for nested CV
#'
#' @examples
#' \dontrun{
#' at <- make_autotuner_xgboost(task, n_evals = 30)
#' rr <- resample(task, at, rsmp("cv", folds = 5))
#' }
#'
#' @export
make_autotuner_xgboost <- function(task,
                                    n_evals = NULL,
                                    inner_folds = 3,
                                    measure = "classif.auc",
                                    early_stopping = TRUE,
                                    nthread = 4L) {

  check_mbo_packages()

  # Auto-scale n_evals based on sample size
  if (is.null(n_evals)) {
    n <- task$nrow
    n_evals <- if (n < 100) 25 else if (n < 300) 40 else 60
  }

  # Create learner
  learner <- mlr3::lrn("classif.xgboost",
    predict_type = "prob",
    eval_metric = "logloss",
    nthread = nthread,
    verbosity = 0
  )

  # Omics-optimized search space
  # Strong regularization, shallow trees, feature subsampling
  if (early_stopping) {
    search_space <- paradox::ps(
      nrounds = paradox::p_int(lower = 100, upper = 2000),
      eta = paradox::p_dbl(lower = 0.01, upper = 0.3, logscale = TRUE),
      max_depth = paradox::p_int(lower = 2, upper = 6),
      min_child_weight = paradox::p_dbl(lower = 1, upper = 50, logscale = TRUE),
      subsample = paradox::p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = paradox::p_dbl(lower = 0.05, upper = 0.8),
      lambda = paradox::p_dbl(lower = 1e-3, upper = 100, logscale = TRUE),
      alpha = paradox::p_dbl(lower = 1e-3, upper = 10, logscale = TRUE)
    )
  } else {
    search_space <- paradox::ps(
      nrounds = paradox::p_int(lower = 200, upper = 1000),
      eta = paradox::p_dbl(lower = 0.01, upper = 0.3, logscale = TRUE),
      max_depth = paradox::p_int(lower = 2, upper = 6),
      min_child_weight = paradox::p_dbl(lower = 1, upper = 50, logscale = TRUE),
      subsample = paradox::p_dbl(lower = 0.5, upper = 1),
      colsample_bytree = paradox::p_dbl(lower = 0.05, upper = 0.8),
      lambda = paradox::p_dbl(lower = 1e-3, upper = 100, logscale = TRUE),
      alpha = paradox::p_dbl(lower = 1e-3, upper = 10, logscale = TRUE)
    )
  }

  create_autotuner(learner, search_space, n_evals, inner_folds, measure)
}


#' @title Create Bayesian-Optimized AutoTuner for LightGBM
#'
#' @description
#' Creates an AutoTuner with LightGBM learner and omics-optimized search space.
#'
#' @param task An mlr3 Task
#' @param n_evals Number of Bayesian optimization iterations
#' @param inner_folds Number of inner CV folds
#' @param measure Performance measure
#' @param nthread Number of threads
#'
#' @return An AutoTuner object
#'
#' @export
make_autotuner_lightgbm <- function(task,
                                     n_evals = NULL,
                                     inner_folds = 3,
                                     measure = "classif.auc",
                                     nthread = 4L) {

  check_mbo_packages()

  if (!requireNamespace("lightgbm", quietly = TRUE)) {
    stop("Package 'lightgbm' is required. Install with: install.packages('lightgbm')",
         call. = FALSE)
  }

  # Auto-scale n_evals
  if (is.null(n_evals)) {
    n <- task$nrow
    n_evals <- if (n < 100) 25 else if (n < 300) 40 else 60
  }

 learner <- mlr3::lrn("classif.lightgbm",
    predict_type = "prob",
    verbose = -1,
    num_threads = nthread
  )

  # Omics-optimized: small leaves, strong regularization, feature subsampling
  search_space <- paradox::ps(
    num_iterations = paradox::p_int(lower = 100, upper = 2000),
    learning_rate = paradox::p_dbl(lower = 0.01, upper = 0.2, logscale = TRUE),
    num_leaves = paradox::p_int(lower = 8, upper = 128),
    max_depth = paradox::p_int(lower = 3, upper = 8),
    min_data_in_leaf = paradox::p_int(lower = 5, upper = 100),
    feature_fraction = paradox::p_dbl(lower = 0.05, upper = 0.8),
    bagging_fraction = paradox::p_dbl(lower = 0.5, upper = 1),
    lambda_l1 = paradox::p_dbl(lower = 1e-3, upper = 10, logscale = TRUE),
    lambda_l2 = paradox::p_dbl(lower = 1e-3, upper = 100, logscale = TRUE)
  )

  create_autotuner(learner, search_space, n_evals, inner_folds, measure)
}


#' @title Create Bayesian-Optimized AutoTuner for Random Forest
#'
#' @description
#' Creates an AutoTuner with ranger and omics-optimized search space.
#'
#' @param task An mlr3 Task
#' @param n_evals Number of Bayesian optimization iterations
#' @param inner_folds Number of inner CV folds
#' @param measure Performance measure
#' @param num_trees Number of trees (fixed, usually not worth tuning)
#'
#' @return An AutoTuner object
#'
#' @export
make_autotuner_ranger <- function(task,
                                   n_evals = NULL,
                                   inner_folds = 3,
                                   measure = "classif.auc",
                                   num_trees = 500L) {

  check_mbo_packages()

  # Auto-scale n_evals (RF is faster, can do more)
  if (is.null(n_evals)) {
    n <- task$nrow
    n_evals <- if (n < 100) 30 else if (n < 300) 50 else 80
  }

  # Data-aware mtry bounds
  p <- length(task$feature_names)
  mtry_lower <- max(1, floor(sqrt(p) / 2))
  mtry_upper <- min(p, floor(p * 0.5))

  learner <- mlr3::lrn("classif.ranger",
    predict_type = "prob",
    num.trees = num_trees,
    importance = "impurity"
  )

  search_space <- paradox::ps(
    mtry = paradox::p_int(lower = mtry_lower, upper = mtry_upper),
    min.node.size = paradox::p_int(lower = 1, upper = 20),
    sample.fraction = paradox::p_dbl(lower = 0.5, upper = 1)
  )

  create_autotuner(learner, search_space, n_evals, inner_folds, measure)
}


#' @title Create Bayesian-Optimized AutoTuner for glmnet
#'
#' @description
#' Creates an AutoTuner with glmnet (elastic net) and search space for alpha.
#'
#' @param task An mlr3 Task
#' @param n_evals Number of Bayesian optimization iterations
#' @param inner_folds Number of inner CV folds
#' @param measure Performance measure
#'
#' @return An AutoTuner object
#'
#' @export
make_autotuner_glmnet <- function(task,
                                   n_evals = NULL,
                                   inner_folds = 3,
                                   measure = "classif.auc") {

  check_mbo_packages()

  # glmnet is fast, fewer evals needed
  if (is.null(n_evals)) {
    n_evals <- 20
  }

  learner <- mlr3::lrn("classif.glmnet",
    predict_type = "prob"
  )

  # Tune alpha (mixing parameter) and s (lambda)
  search_space <- paradox::ps(
    alpha = paradox::p_dbl(lower = 0, upper = 1),
    s = paradox::p_dbl(lower = 1e-4, upper = 1, logscale = TRUE)
  )

  create_autotuner(learner, search_space, n_evals, inner_folds, measure)
}


#' @title Create AutoTuner with MBO
#'
#' @description
#' Internal helper to create an AutoTuner with mlr3mbo.
#'
#' @param learner An mlr3 Learner
#' @param search_space A paradox ParamSet
#' @param n_evals Number of evaluations
#' @param inner_folds Inner CV folds
#' @param measure Measure name or object
#'
#' @return An AutoTuner object
#'
#' @keywords internal
create_autotuner <- function(learner, search_space, n_evals, inner_folds, measure) {

  # Get measure object
  if (is.character(measure)) {
    measure <- mlr3::msr(measure)
  }

  # Inner resampling
  inner_rsmp <- mlr3::rsmp("cv", folds = inner_folds)

  # Terminator
  terminator <- mlr3tuning::trm("evals", n_evals = n_evals)

  # MBO tuner (uses RF surrogate by default, which is robust)
  tuner <- mlr3tuning::tnr("mbo")

  # Create AutoTuner
  mlr3tuning::AutoTuner$new(
    learner = learner,
    resampling = inner_rsmp,
    measure = measure,
    search_space = search_space,
    terminator = terminator,
    tuner = tuner
  )
}


#' @title Check MBO Package Dependencies
#'
#' @description
#' Verifies that required packages for Bayesian optimization are installed.
#'
#' @keywords internal
check_mbo_packages <- function() {
  if (!requireNamespace("mlr3mbo", quietly = TRUE)) {
    stop("Package 'mlr3mbo' is required for Bayesian optimization. ",
         "Install with: install.packages('mlr3mbo')", call. = FALSE)
  }
  if (!requireNamespace("mlr3tuning", quietly = TRUE)) {
    stop("Package 'mlr3tuning' is required. Install with: install.packages('mlr3tuning')",
         call. = FALSE)
  }
  if (!requireNamespace("paradox", quietly = TRUE)) {
    stop("Package 'paradox' is required. Install with: install.packages('paradox')",
         call. = FALSE)
  }
}


#' @title Run Bayesian Optimization Benchmark
#'
#' @description
#' Convenience function to run Bayesian-optimized learners in nested CV
#' and compare performance.
#'
#' @param task An mlr3 Task
#' @param learners Character vector of learner names: "xgboost", "lightgbm", "ranger", "glmnet"
#' @param outer_folds Number of outer CV folds (default: 5)
#' @param inner_folds Number of inner CV folds (default: 3)
#' @param n_evals Number of MBO iterations per learner (default: auto)
#' @param measure Performance measure (default: classif.auc)
#' @param seed Random seed
#'
#' @return A benchmark result with tuned learners
#'
#' @examples
#' \dontrun{
#' result <- run_bayesian_benchmark(
#'   task,
#'   learners = c("xgboost", "ranger", "glmnet"),
#'   outer_folds = 5,
#'   n_evals = 30
#' )
#' result$aggregate(msr("classif.auc"))
#' }
#'
#' @export
run_bayesian_benchmark <- function(task,
                                    learners = c("xgboost", "ranger", "glmnet"),
                                    outer_folds = 5,
                                    inner_folds = 3,
                                    n_evals = NULL,
                                    measure = "classif.auc",
                                    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Create AutoTuners for each learner
  autotuners <- list()

  for (lrn_name in learners) {
    at <- switch(lrn_name,
      "xgboost" = make_autotuner_xgboost(task, n_evals, inner_folds, measure),
      "lightgbm" = make_autotuner_lightgbm(task, n_evals, inner_folds, measure),
      "ranger" = make_autotuner_ranger(task, n_evals, inner_folds, measure),
      "glmnet" = make_autotuner_glmnet(task, n_evals, inner_folds, measure),
      stop(sprintf("Unknown learner: %s", lrn_name), call. = FALSE)
    )
    at$id <- sprintf("tuned_%s", lrn_name)
    autotuners[[lrn_name]] <- at
  }

  # Outer resampling
  outer_rsmp <- mlr3::rsmp("cv", folds = outer_folds)

  # Create benchmark design
  design <- mlr3::benchmark_grid(
    tasks = list(task),
    learners = autotuners,
    resamplings = list(outer_rsmp)
  )

  # Run benchmark
  mlr3::benchmark(design, store_models = TRUE)
}


#' @title Get Optimal Hyperparameters from AutoTuner
#'
#' @description
#' Extracts the best hyperparameters found by Bayesian optimization.
#'
#' @param autotuner A trained AutoTuner object
#'
#' @return A list of optimal hyperparameter values
#'
#' @export
get_optimal_params <- function(autotuner) {
  if (!inherits(autotuner, "AutoTuner")) {
    stop("Input must be an AutoTuner object", call. = FALSE)
  }

  if (is.null(autotuner$tuning_result)) {
    stop("AutoTuner has not been trained yet", call. = FALSE)
  }

  autotuner$tuning_result$learner_param_vals[[1]]
}
