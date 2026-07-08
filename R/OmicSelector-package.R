#' @title OmicSelector: Zero-Leakage Biomarker Discovery Toolkit
#'
#' @description
#' OmicSelector 2.6.0 is a toolkit for high-dimensional biomarker discovery with
#' leakage-resistant validation, stable feature selection, within-sample panel
#' methods, frozen-reference denoising add-ons, qPCR non-detect handling,
#' matched-null benchmarks, and provenance-aware audit helpers.
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
#' - **Panel Governance**: Group-aware resampling, provenance-floor diagnostics,
#'   identifiability gates, matched-null benchmarks, and operating-point summaries
#'   for biomarker panels evaluated with public omics data.
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
#' OmicSelector 2.6.0 prioritizes leakage-resistant validation, feature
#' stability, and auditable provenance checks above raw apparent accuracy.
#' High accuracy with unstable feature sets or strong provenance-only prediction
#' is treated as a warning signal rather than as direct biomarker evidence.
#'
#' @name OmicSelector-package
#' @aliases OmicSelector
#'
#' @import R6
#' @importFrom data.table as.data.table
#' @importFrom stats aggregate aov as.stepfun binomial coef cor glm isoreg
#'   .lm.fit mad median optimize pnorm predict qnorm quantile rbinom rnorm
#'   p.adjust runif sd setNames weighted.mean
#' @importFrom methods new
#' @importFrom utils capture.output head modifyList sessionInfo tail
#'
#' @seealso
#' - [OmicPipeline]: Main pipeline class
#' - [BenchmarkService]: Nested CV service
#'
"_PACKAGE"

#' @title OmicSelector 2.6.0 Package Startup
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "OmicSelector 2.6.0 - Leakage-Resistant Biomarker Discovery\n",
    "Use OmicPipeline$new() to start. Run ?OmicSelector for help.\n"
  )
}
