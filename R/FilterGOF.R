#' @title Goodness-of-Fit Filters for Sparse Omics Data
#'
#' @description
#' Implements feature selection filters designed for sparse/zero-inflated omics data.
#' Standard variance-based filters often discard biologically relevant features
#' that are "turned off" (zero) in one condition but expressed in another.
#'
#' @details
#' Two filters are provided:
#' - **FilterGOF_KS**: Kolmogorov-Smirnov test comparing distributions between classes
#' - **FilterHurdle**: Two-part hurdle model testing both zero-frequency and magnitude
#'
#' These filters are designed for binary classification tasks in biomarker discovery.
#'
#' @name FilterGOF
NULL


#' @title Kolmogorov-Smirnov GOF Filter
#'
#' @description
#' Uses the KS test D-statistic as feature importance score.
#' Captures differences in distribution shape, location, and zero-frequency
#' between two classes without assuming a parametric distribution.
#'
#' @details
#' The KS test compares Empirical Cumulative Distribution Functions (ECDF).
#' If one class has 80% zeros and the other has 20% zeros, the ECDFs show
#' a massive difference captured by the D statistic (bounded [0, 1]).
#'
#' Higher D values indicate better class separation.
#'
#' @section Dictionary:
#' This Filter can be instantiated via the dictionary `mlr_filters` or with
#' the associated sugar function `flt()`:
#' ```
#' flt("gof_ks")
#' mlr_filters$get("gof_ks")
#' ```
#'
#' @references
#' Kolmogorov, A. N. (1933). Sulla determinazione empirica di una legge di distribuzione.
#'
#' @export
FilterGOF_KS <- R6::R6Class(
  "FilterGOF_KS",
  inherit = mlr3filters::Filter,

  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "gof_ks",
        task_types = "classif",
        param_set = paradox::ps(
          min_nonzero_prop = paradox::p_dbl(lower = 0, upper = 1, default = 0.01)
        ),
        feature_types = c("numeric", "integer"),
        packages = character(),
        label = "Kolmogorov-Smirnov GOF",
        man = "OmicSelector::FilterGOF_KS"
      )
    }
  ),

  private = list(
    .calculate = function(task, nfeat) {
      # Validate task is binary classification
      if (length(task$class_names) != 2) {
        stop("FilterGOF_KS requires a binary classification task", call. = FALSE)
      }

      # Get data
      data <- task$data()
      target <- task$target_names
      features <- task$feature_names

      # Get class labels
      y <- data[[target]]
      classes <- levels(y)

      # Minimum proportion of non-zero values
      min_nonzero <- self$param_set$values$min_nonzero_prop %||% 0.01

      # Compute KS D-statistic for each feature
      scores <- vapply(features, function(feat) {
        x <- as.numeric(data[[feat]])

        # Split by class
        x1 <- x[y == classes[1]]
        x2 <- x[y == classes[2]]

        # Check minimum non-zero proportion
        nonzero_prop <- mean(x != 0, na.rm = TRUE)
        if (nonzero_prop < min_nonzero) {
          return(0)
        }

        # KS test (suppress warnings for ties)
        ks_result <- tryCatch({
          suppressWarnings(stats::ks.test(x1, x2))
        }, error = function(e) {
          return(list(statistic = 0))
        })

        # Return D statistic (higher = more separation)
        as.numeric(ks_result$statistic)
      }, FUN.VALUE = numeric(1))

      scores
    }
  )
)


#' @title Hurdle Filter for Zero-Inflated Data
#'
#' @description
#' Two-part filter designed for zero-inflated count/expression data.
#' Tests both the difference in zero-frequency AND the difference in
#' magnitude among non-zero values.
#'
#' @details
#' The Hurdle filter combines two tests:
#' 1. **Binary component**: Chi-square test on zero vs non-zero contingency table
#' 2. **Continuous component**: Wilcoxon test on non-zero values
#'
#' The score is computed using Fisher's method to combine p-values:
#' \code{chi2 = -2 * sum(log(p_values))}
#'
#' This captures genes that are good biomarkers because:
#' - They are detectable in one class but not the other (zero-frequency difference)
#' - They are upregulated in one class among expressed samples (magnitude difference)
#'
#' @section Dictionary:
#' This Filter can be instantiated via the dictionary `mlr_filters` or with
#' the associated sugar function `flt()`:
#' ```
#' flt("hurdle")
#' mlr_filters$get("hurdle")
#' ```
#'
#' @export
FilterHurdle <- R6::R6Class(
  "FilterHurdle",
  inherit = mlr3filters::Filter,

  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "hurdle",
        task_types = "classif",
        param_set = paradox::ps(
          min_nonzero_per_class = paradox::p_int(lower = 1, upper = 100, default = 3),
          combine_method = paradox::p_fct(levels = c("fisher", "min"), default = "fisher")
        ),
        feature_types = c("numeric", "integer"),
        packages = character(),
        label = "Hurdle (Zero-Inflated)",
        man = "OmicSelector::FilterHurdle"
      )
    }
  ),

  private = list(
    .calculate = function(task, nfeat) {
      # Validate task is binary classification
      if (length(task$class_names) != 2) {
        stop("FilterHurdle requires a binary classification task", call. = FALSE)
      }

      # Get data
      data <- task$data()
      target <- task$target_names
      features <- task$feature_names

      # Get class labels
      y <- data[[target]]
      classes <- levels(y)

      # Parameters
      min_nonzero <- self$param_set$values$min_nonzero_per_class %||% 3
      combine_method <- self$param_set$values$combine_method %||% "fisher"

      # Compute hurdle score for each feature
      scores <- vapply(features, function(feat) {
        x <- as.numeric(data[[feat]])

        # Split by class
        x1 <- x[y == classes[1]]
        x2 <- x[y == classes[2]]

        # Binary component: zero vs non-zero contingency table
        n1_zero <- sum(x1 == 0, na.rm = TRUE)
        n1_nonzero <- sum(x1 != 0, na.rm = TRUE)
        n2_zero <- sum(x2 == 0, na.rm = TRUE)
        n2_nonzero <- sum(x2 != 0, na.rm = TRUE)

        # Check if we have enough non-zero values
        if (n1_nonzero < min_nonzero || n2_nonzero < min_nonzero) {
          # Not enough non-zero values for continuous test
          # Use only binary component
          p_binary <- tryCatch({
            cont_table <- matrix(c(n1_zero, n1_nonzero, n2_zero, n2_nonzero),
                                 nrow = 2, byrow = TRUE)
            # Use Fisher's exact test for small samples
            if (min(cont_table) < 5) {
              fisher.test(cont_table)$p.value
            } else {
              chisq.test(cont_table, correct = FALSE)$p.value
            }
          }, error = function(e) 1)

          # Return -log10(p) as score (higher = more significant)
          return(-log10(max(p_binary, 1e-300)))
        }

        # Binary component
        p_binary <- tryCatch({
          cont_table <- matrix(c(n1_zero, n1_nonzero, n2_zero, n2_nonzero),
                               nrow = 2, byrow = TRUE)
          if (min(cont_table) < 5) {
            fisher.test(cont_table)$p.value
          } else {
            chisq.test(cont_table, correct = FALSE)$p.value
          }
        }, error = function(e) 1)

        # Continuous component: Wilcoxon test on non-zero values
        x1_nonzero <- x1[x1 != 0]
        x2_nonzero <- x2[x2 != 0]

        p_continuous <- tryCatch({
          suppressWarnings(wilcox.test(x1_nonzero, x2_nonzero)$p.value)
        }, error = function(e) 1)

        # Combine p-values
        if (combine_method == "fisher") {
          # Fisher's method: chi2 = -2 * sum(log(p))
          p_values <- c(p_binary, p_continuous)
          p_values <- pmax(p_values, 1e-300)  # Avoid log(0)
          fisher_stat <- -2 * sum(log(p_values))
          return(fisher_stat)
        } else {
          # Minimum p-value approach
          p_min <- min(p_binary, p_continuous)
          return(-log10(max(p_min, 1e-300)))
        }
      }, FUN.VALUE = numeric(1))

      scores
    }
  )
)


#' @title Zero-Proportion Filter
#'
#' @description
#' Simple filter based on the difference in zero-proportion between classes.
#' Useful for quickly identifying on/off switches in gene expression.
#'
#' @details
#' Computes: |P(X=0|Class1) - P(X=0|Class2)|
#' Higher values indicate features that are expressed in one class but not the other.
#'
#' @export
FilterZeroProp <- R6::R6Class(
  "FilterZeroProp",
  inherit = mlr3filters::Filter,

  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "zero_prop",
        task_types = "classif",
        param_set = paradox::ps(),
        feature_types = c("numeric", "integer"),
        packages = character(),
        label = "Zero-Proportion Difference",
        man = "OmicSelector::FilterZeroProp"
      )
    }
  ),

  private = list(
    .calculate = function(task, nfeat) {
      if (length(task$class_names) != 2) {
        stop("FilterZeroProp requires a binary classification task", call. = FALSE)
      }

      data <- task$data()
      target <- task$target_names
      features <- task$feature_names

      y <- data[[target]]
      classes <- levels(y)

      scores <- vapply(features, function(feat) {
        x <- as.numeric(data[[feat]])

        # Zero proportions per class
        prop1 <- mean(x[y == classes[1]] == 0, na.rm = TRUE)
        prop2 <- mean(x[y == classes[2]] == 0, na.rm = TRUE)

        # Absolute difference
        abs(prop1 - prop2)
      }, FUN.VALUE = numeric(1))

      scores
    }
  )
)


#' @title Register GOF Filters in mlr3
#'
#' @description
#' Registers the GOF filters in the mlr3filters dictionary.
#' Call this function to make filters available via `flt("gof_ks")`, etc.
#'
#' @export
register_gof_filters <- function() {
  if (!requireNamespace("mlr3filters", quietly = TRUE)) {
    stop("Package 'mlr3filters' is required", call. = FALSE)
  }

  # Register filters
  mlr3filters::mlr_filters$add("gof_ks", FilterGOF_KS)
  mlr3filters::mlr_filters$add("hurdle", FilterHurdle)
  mlr3filters::mlr_filters$add("zero_prop", FilterZeroProp)

  invisible(TRUE)
}


#' @title Create GOF Filter
#'
#' @description
#' Convenience function to create GOF filters.
#'
#' @param type Filter type: "ks", "hurdle", or "zero_prop"
#' @param ... Additional parameters passed to the filter
#'
#' @return A Filter object
#'
#' @examples
#' \dontrun{
#' filt <- make_gof_filter("ks")
#' filt <- make_gof_filter("hurdle", min_nonzero_per_class = 5)
#' }
#'
#' @export
make_gof_filter <- function(type = c("ks", "hurdle", "zero_prop"), ...) {
  type <- match.arg(type)

  filter <- switch(type,
    "ks" = FilterGOF_KS$new(),
    "hurdle" = FilterHurdle$new(),
    "zero_prop" = FilterZeroProp$new()
  )

  # Set parameters
  params <- list(...)
  if (length(params) > 0) {
    valid_params <- intersect(names(params), filter$param_set$ids())
    if (length(valid_params) > 0) {
      filter$param_set$set_values(.values = params[valid_params])
    }
  }

  filter
}


#' @title Compare GOF Filters on Task
#'
#' @description
#' Compares different GOF filter rankings on a task.
#'
#' @param task An mlr3 classification task
#' @param top_n Number of top features to compare
#'
#' @return A data frame with filter rankings
#'
#' @export
compare_gof_filters <- function(task, top_n = 20) {

  filters <- list(
    gof_ks = FilterGOF_KS$new(),
    hurdle = FilterHurdle$new(),
    zero_prop = FilterZeroProp$new()
  )

  # Add ANOVA for comparison if available
  if (requireNamespace("mlr3filters", quietly = TRUE)) {
    filters$anova <- mlr3filters::flt("anova")
  }

  results <- lapply(names(filters), function(name) {
    filt <- filters[[name]]
    scores <- filt$calculate(task)
    scores <- sort(scores, decreasing = TRUE)

    data.frame(
      filter = name,
      rank = seq_along(scores),
      feature = names(scores),
      score = as.numeric(scores),
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, results)

  # Pivot to show top features per filter
  top_df <- result_df[result_df$rank <= top_n, ]

  top_df
}
