#' @title OmicSelector: Zero-Leakage Biomarker Discovery Toolkit
#'
#' @description
#' OmicSelector 2.0 is a PhD-level toolkit for high-dimensional biomarker discovery
#' that guarantees scientific validity through rigorous machine learning methodology.
#'
#' @details
#' ## Key Features
#'
#' - **Zero Data Leakage**: All preprocessing, feature selection, and model training
#'   occurs within proper cross-validation folds via mlr3 GraphLearners.
#'
#' - **Nested Cross-Validation**: Inner loop for feature selection and hyperparameter
#'   tuning, outer loop for unbiased performance estimation.
#'
#' - **Stability Metrics**: Nogueira Stability Index to ensure selected features are
#'   robust across resamples, not just high accuracy.
#'
#' - **Reproducibility**: renv lockfiles, Docker containers, and deterministic pipelines.
#'
#' ## Core Classes
#'
#' - [OmicPipeline]: Central R6 class for creating leakage-free pipelines
#' - [BenchmarkService]: Enforces proper nested cross-validation
#'
#' ## Quick Start
#'
#' ```r
#' library(OmicSelector)
#'
#' # Create pipeline from data
#' pipeline <- OmicPipeline$new(
#'   data = my_data,
#'   target = "outcome",
#'   positive = "Case"
#' )
#'
#' # Create graph learner with feature selection
#' learner <- pipeline$create_graph_learner(
#'   filter = "anova",
#'   model = "ranger",
#'   n_features = 20
#' )
#'
#' # Run nested CV benchmark
#' service <- BenchmarkService$new(pipeline, outer_folds = 5, inner_folds = 3)
#' service$add_learner(learner)
#' result <- service$run()
#' ```
#'
#' ## Philosophy
#'
#' "Optimization without validation is hallucination."
#'
#' OmicSelector 2.0 prioritizes **zero data leakage** and **feature stability**
#' above raw accuracy metrics. High accuracy with unstable feature sets
#' (the "Rashomon Effect") indicates overfitting, not real signal.
#'
#' @docType package
#' @name OmicSelector-package
#' @aliases OmicSelector
#'
#' @import R6
#' @importFrom stats cor rnorm rbinom runif aov setNames
#' @importFrom utils head tail
#'
#' @seealso
#' - [OmicPipeline]: Main pipeline class
#' - [BenchmarkService]: Nested CV service
#'
NULL

#' @title OmicSelector 2.0 Package Startup
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "OmicSelector 2.0 - Zero-Leakage Biomarker Discovery\n",
    "Use OmicPipeline$new() to start. Run ?OmicSelector for help.\n"
  )
}
