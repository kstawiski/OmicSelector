#' @title CoDA-Aware Feature Selection for Biomarker Panels
#'
#' @description
#' Feature-selection methods tailored to compositional omics data. The exported
#' functions return per-feature importance scores and selected features, while
#' the accompanying `mlr3filters::Filter` classes expose the methods inside the
#' existing `OmicPipeline` feature-selection path.
#'
#' The implemented selectors cover five complementary CoDA-aware strategies:
#' - `codaFS_plr_variance()`: pairwise log-ratio screening with feature
#'   aggregation.
#' - `codaFS_selbal_wrapper()`: sparse forward-balance selection using the
#'   official `selbal` package when available, or a native fallback.
#' - `codaFS_codacore_wrapper()`: sparse balance selection via CoDaCoRe using
#'   the existing OmicSelector interface.
#' - `codaFS_logcontrast_lasso()`: CLR-space lasso with post-fit zero-sum
#'   projection.
#' - `codaFS_stability_logratio()`: stability selection on pairwise log-ratio
#'   features with elastic-net aggregation back to individual miRNAs.
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Rivera-Pinto J, Egozcue JJ, Pawlowsky-Glahn V, et al. (2018). Balances:
#' a new perspective for microbiome analysis. \emph{mSystems}, 3(4),
#' e00053-18.
#'
#' Gordon-Rodriguez E, Susin A, McEwen JD, et al. (2022). CoDaCoRe: learning
#' sparse log-ratios for high-throughput sequencing data.
#' \emph{Bioinformatics}, 38(12), 3179-3186.
#'
#' Lin W, Shi P, Feng R, Li H. (2014). Variable selection in regression with
#' compositional covariates. \emph{Biometrika}, 101(4), 785-797.
#'
#' Meinshausen N, Buhlmann P. (2010). Stability selection.
#' \emph{Journal of the Royal Statistical Society Series B}, 72(4), 417-473.
#'
#' @name coda-feature-selection
NULL


#' @title Pairwise Log-Ratio Variance Feature Filter
#'
#' @description
#' Scores pairwise log-ratios by combining their across-sample variance with
#' their association to the binary outcome, then aggregates the highest-scoring
#' ratios back to individual features by weighted frequency.
#'
#' @param X Numeric matrix with samples in rows and features in columns.
#' @param y Binary outcome vector.
#' @param n_features Optional target panel size. Used only to derive the default
#'   number of pairwise log-ratios retained during aggregation.
#' @param positive Optional positive-class label for non-numeric outcomes.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param top_pairs Number of top pairwise log-ratios to aggregate. If `NULL`,
#'   this is set to `pair_multiplier * n_features` or a conservative default.
#' @param pair_multiplier Multiplier used when `top_pairs` is `NULL`.
#' @param max_pair_features Maximum number of features retained before computing
#'   all pairwise log-ratios. Prescreening uses training-fold-only CLR effects.
#' @param prescreen_metric Either `"clr_effect"` or `"variance"`.
#'
#' @return A list with per-feature `scores`, `selected_features`, and the
#'   aggregated top pairwise log-ratio table.
#'
#' @export
codaFS_plr_variance <- function(X, y,
                                n_features = NULL,
                                positive = NULL,
                                input_scale = c("log", "raw"),
                                top_pairs = NULL,
                                pair_multiplier = 5L,
                                max_pair_features = 200L,
                                prescreen_metric = c("clr_effect", "variance")) {
  input_scale <- match.arg(input_scale)
  prescreen_metric <- match.arg(prescreen_metric)
  prep <- .coda_fs_prepare_inputs(X, y, positive = positive, input_scale = input_scale)
  candidate <- .coda_fs_prescreen(
    prep$log_x,
    prep$y,
    max_pair_features = max_pair_features,
    metric = prescreen_metric
  )
  local_x <- prep$log_x[, candidate, drop = FALSE]
  pair_stats <- .coda_fs_pair_summary(local_x, prep$y)

  n_features_eff <- .coda_fs_target_size(n_features, ncol(local_x))
  if (is.null(top_pairs)) {
    top_pairs <- max(10L, as.integer(pair_multiplier) * n_features_eff)
  }
  top_pairs <- min(as.integer(top_pairs), nrow(pair_stats))
  if (top_pairs > 0L) {
    pair_stats <- pair_stats[order(pair_stats$score, decreasing = TRUE), , drop = FALSE]
    pair_stats <- pair_stats[seq_len(top_pairs), , drop = FALSE]
  }

  local_scores <- .coda_fs_aggregate_pair_scores(
    pair_table = pair_stats,
    feature_names = colnames(local_x),
    weight_col = "score"
  )
  scores <- setNames(rep(0, ncol(prep$log_x)), prep$feature_names)
  scores[names(local_scores)] <- local_scores

  list(
    method = "coda_plr_variance",
    scores = scores,
    selected_features = .coda_fs_select_top(scores, n_features),
    pair_table = .coda_fs_label_pair_table(pair_stats, colnames(local_x))
  )
}


#' @title Selbal-Style Forward Balance Feature Selection
#'
#' @description
#' Wraps the official `selbal` implementation when it is available. In
#' environments where `selbal` is not installed, or when `backend = "fallback"`
#' is requested, the function uses a native forward balance-selection routine
#' that greedily improves univariate logistic deviance on training data only.
#'
#' @param X Numeric matrix with samples in rows and features in columns.
#' @param y Binary outcome vector.
#' @param n_features Optional target panel size. Used as the default
#'   `max_terms`.
#' @param positive Optional positive-class label for non-numeric outcomes.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param backend One of `"auto"`, `"official"`, or `"fallback"`.
#' @param max_terms Maximum number of features in the selected balance. If
#'   `NULL`, defaults to `n_features`.
#' @param min_improvement Minimum deviance improvement required to add a new
#'   feature during the fallback forward search.
#' @param max_pair_features Maximum number of prescreened features considered by
#'   the forward search.
#' @param prescreen_metric Either `"clr_effect"` or `"variance"`.
#'
#' @return A list with per-feature `scores`, `selected_features`, balance
#'   membership, and backend details.
#'
#' @export
codaFS_selbal_wrapper <- function(X, y,
                                  n_features = NULL,
                                  positive = NULL,
                                  input_scale = c("log", "raw"),
                                  backend = c("auto", "official", "fallback"),
                                  max_terms = NULL,
                                  min_improvement = 1e-4,
                                  max_pair_features = 80L,
                                  prescreen_metric = c("clr_effect", "variance")) {
  input_scale <- match.arg(input_scale)
  backend <- match.arg(backend)
  prescreen_metric <- match.arg(prescreen_metric)
  prep <- .coda_fs_prepare_inputs(X, y, positive = positive, input_scale = input_scale)

  max_terms <- .coda_fs_target_size(max_terms %||% n_features, ncol(prep$log_x))
  result <- NULL
  selected_backend <- "fallback"

  if (backend %in% c("auto", "official") && requireNamespace("selbal", quietly = TRUE)) {
    official <- tryCatch(
      {
        y_factor <- factor(
          prep$y,
          levels = c(0L, 1L),
          labels = c(prep$negative, prep$positive)
        )
        fit <- selbal::selbal(
          x = prep$log_x,
          y = y_factor,
          th.imp = min_improvement,
          logt = FALSE,
          tab = FALSE,
          draw = FALSE
        )
        num <- as.character(fit$numerator)
        den <- as.character(fit$denominator)
        scores <- setNames(rep(0, length(prep$feature_names)), prep$feature_names)
        selected <- unique(c(num, den))
        if (length(selected) > 0L) {
          scores[selected] <- rev(seq_along(selected))
        }
        list(
          method = "coda_selbal",
          backend = "official",
          scores = scores,
          numerator = num,
          denominator = den,
          selected_features = .coda_fs_select_top(scores, n_features),
          history = NULL
        )
      },
      error = function(e) {
        if (identical(backend, "official")) {
          stop("Official selbal backend failed: ", conditionMessage(e), call. = FALSE)
        }
        NULL
      }
    )
    if (!is.null(official)) {
      result <- official
      selected_backend <- "official"
    }
  }

  if (is.null(result)) {
    candidate <- .coda_fs_prescreen(
      prep$log_x,
      prep$y,
      max_pair_features = max_pair_features,
      metric = prescreen_metric
    )
    fallback <- .coda_fs_forward_balance(
      prep$log_x[, candidate, drop = FALSE],
      prep$y,
      max_terms = max_terms,
      min_improvement = min_improvement
    )
    local_names <- colnames(prep$log_x)[candidate]
    numerator <- local_names[fallback$numerator]
    denominator <- local_names[fallback$denominator]
    ordered <- c(
      local_names[fallback$selection_order[fallback$selection_side == "numerator"]],
      local_names[fallback$selection_order[fallback$selection_side == "denominator"]]
    )
    ordered <- unique(ordered)
    scores <- setNames(rep(0, length(prep$feature_names)), prep$feature_names)
    if (length(ordered) > 0L) {
      scores[ordered] <- rev(seq_along(ordered))
    }
    result <- list(
      method = "coda_selbal",
      backend = selected_backend,
      scores = scores,
      numerator = numerator,
      denominator = denominator,
      selected_features = .coda_fs_select_top(scores, n_features),
      history = fallback$history
    )
  }

  result
}


#' @title CoDaCoRe-Derived Feature Selection
#'
#' @description
#' Selects features appearing in sparse CoDaCoRe balances. The official
#' `codacore` backend is used when available and functional; otherwise the
#' existing OmicSelector sparse-balance fallback is reused. Feature scores are
#' aggregated across balances using absolute slopes when available.
#'
#' @param X Numeric matrix with samples in rows and features in columns.
#' @param y Binary outcome vector.
#' @param n_features Optional target panel size.
#' @param positive Optional positive-class label for non-numeric outcomes.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param backend One of `"auto"`, `"official"`, or `"fallback"`.
#' @param max_balances Maximum number of balances to extract.
#' @param top_k Initial sparse balance size for the fallback backend.
#' @param class_weight Logical; passed to the fallback backend.
#' @param verbose Logical; print backend messages.
#' @param seed Integer random seed.
#'
#' @return A list with per-feature `scores`, `selected_features`, extracted
#'   balances, and backend details.
#'
#' @export
codaFS_codacore_wrapper <- function(X, y,
                                    n_features = NULL,
                                    positive = NULL,
                                    input_scale = c("log", "raw"),
                                    backend = c("auto", "official", "fallback"),
                                    max_balances = NULL,
                                    top_k = 3L,
                                    class_weight = FALSE,
                                    verbose = FALSE,
                                    seed = 42L) {
  input_scale <- match.arg(input_scale)
  backend <- match.arg(backend)
  prep <- .coda_fs_prepare_inputs(X, y, positive = positive, input_scale = input_scale)

  max_balances <- .coda_fs_default_max_balances(max_balances, n_features, ncol(prep$log_x))
  train_backend <- switch(
    backend,
    auto = "auto",
    official = "codacore",
    fallback = "fallback"
  )
  fit <- train_codacore(
    X_train = prep$log_x,
    y_train = factor(prep$y, levels = c(0L, 1L), labels = c(prep$negative, prep$positive)),
    class_weight = class_weight,
    verbose = verbose,
    seed = seed,
    input_scale = "log",
    backend = train_backend,
    max_balances = max_balances,
    top_k = top_k
  )
  if (identical(backend, "official") && !identical(fit$backend, "codacore")) {
    stop("Official codacore backend was requested but did not complete successfully.", call. = FALSE)
  }

  scores <- setNames(rep(0, length(prep$feature_names)), prep$feature_names)
  balances <- list()

  if (identical(fit$backend, "codacore")) {
    n_lr <- codacore::getNumLogRatios(fit$model)
    slopes <- abs(codacore::getSlopes(fit$model))
    for (i in seq_len(n_lr)) {
      num <- as.character(codacore::getNumeratorParts(fit$model, baseLearnerIndex = i, boolean = FALSE))
      den <- as.character(codacore::getDenominatorParts(fit$model, baseLearnerIndex = i, boolean = FALSE))
      weight <- slopes[[i]] %||% 1
      scores[num] <- scores[num] + weight
      scores[den] <- scores[den] + weight
      balances[[i]] <- list(numerator = num, denominator = den, weight = weight)
    }
  } else {
    for (i in seq_along(fit$balance_spec)) {
      spec <- fit$balance_spec[[i]]
      num <- prep$feature_names[spec$numerator]
      den <- prep$feature_names[spec$denominator]
      weight <- 1 / i
      scores[num] <- scores[num] + weight
      scores[den] <- scores[den] + weight
      balances[[i]] <- list(numerator = num, denominator = den, weight = weight)
    }
  }

  list(
    method = "coda_codacore",
    backend = fit$backend,
    scores = scores,
    selected_features = .coda_fs_select_top(scores, n_features),
    balances = balances
  )
}


#' @title Log-Contrast Lasso on CLR Features
#'
#' @description
#' Fits a lasso or elastic-net logistic model on centered log-ratio (CLR)
#' features, then projects the coefficient vector onto the zero-sum hyperplane
#' to obtain a valid log-contrast representation. Feature importance is the
#' absolute projected coefficient magnitude.
#'
#' @param X Numeric matrix with samples in rows and features in columns.
#' @param y Binary outcome vector.
#' @param n_features Optional target panel size.
#' @param positive Optional positive-class label for non-numeric outcomes.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param alpha Elastic-net mixing parameter passed to `glmnet`.
#' @param nfolds Number of cross-validation folds used by `cv.glmnet`.
#' @param nlambda Number of lambda values used by `cv.glmnet`.
#' @param lambda_rule Either `"lambda.1se"` (default) or `"lambda.min"`.
#'
#' @return A list with per-feature `scores`, the raw and zero-sum projected
#'   coefficients, and `selected_features`.
#'
#' @export
codaFS_logcontrast_lasso <- function(X, y,
                                     n_features = NULL,
                                     positive = NULL,
                                     input_scale = c("log", "raw"),
                                     alpha = 1,
                                     nfolds = 5L,
                                     nlambda = 100L,
                                     lambda_rule = c("lambda.1se", "lambda.min")) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for codaFS_logcontrast_lasso().", call. = FALSE)
  }
  input_scale <- match.arg(input_scale)
  lambda_rule <- match.arg(lambda_rule)
  prep <- .coda_fs_prepare_inputs(X, y, positive = positive, input_scale = input_scale)
  clr_x <- clr_transform(prep$log_x, input_scale = "log")
  cv_nfolds <- .coda_fs_cv_folds(prep$y, nfolds)

  cv_fit <- glmnet::cv.glmnet(
    x = clr_x,
    y = prep$y,
    family = "binomial",
    alpha = alpha,
    nfolds = cv_nfolds,
    nlambda = as.integer(nlambda),
    standardize = FALSE,
    intercept = TRUE,
    type.measure = "auc"
  )
  coef_raw <- as.numeric(stats::coef(cv_fit, s = cv_fit[[lambda_rule]]))[-1L]
  coef_projected <- coef_raw - mean(coef_raw)
  names(coef_raw) <- names(coef_projected) <- prep$feature_names
  scores <- abs(coef_projected)

  list(
    method = "coda_logcontrast_lasso",
    scores = scores,
    selected_features = .coda_fs_select_top(scores, n_features),
    coefficients_raw = coef_raw,
    coefficients = coef_projected,
    zero_sum_error = sum(coef_projected),
    active_features = names(coef_raw)[abs(coef_raw) > 0]
  )
}


#' @title Stability Selection on Pairwise Log-Ratios
#'
#' @description
#' Applies subsampling-based stability selection to pairwise log-ratio features
#' using elastic-net logistic models. Stable ratios are aggregated back to
#' individual features by their incidence-weighted selection frequency.
#'
#' @param X Numeric matrix with samples in rows and features in columns.
#' @param y Binary outcome vector.
#' @param n_features Optional target panel size.
#' @param positive Optional positive-class label for non-numeric outcomes.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param alpha Elastic-net mixing parameter.
#' @param subsample_frac Fraction of samples used in each stability iteration.
#' @param n_subsamples Number of stability iterations.
#' @param nfolds Number of cross-validation folds within each subsample.
#' @param nlambda Number of lambda values used by `cv.glmnet`.
#' @param lambda_rule Either `"lambda.1se"` (default) or `"lambda.min"`.
#' @param selection_threshold Frequency threshold used to define stable ratios.
#' @param max_pair_features Maximum number of prescreened features retained
#'   before constructing pairwise ratios.
#' @param prescreen_metric Either `"clr_effect"` or `"variance"`.
#'
#' @return A list with per-feature `scores`, stable pairwise log-ratio
#'   frequencies, and `selected_features`.
#'
#' @export
codaFS_stability_logratio <- function(X, y,
                                      n_features = NULL,
                                      positive = NULL,
                                      input_scale = c("log", "raw"),
                                      alpha = 0.9,
                                      subsample_frac = 0.7,
                                      n_subsamples = 20L,
                                      nfolds = 3L,
                                      nlambda = 60L,
                                      lambda_rule = c("lambda.1se", "lambda.min"),
                                      selection_threshold = 0.6,
                                      max_pair_features = 60L,
                                      prescreen_metric = c("clr_effect", "variance")) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for codaFS_stability_logratio().", call. = FALSE)
  }
  input_scale <- match.arg(input_scale)
  lambda_rule <- match.arg(lambda_rule)
  prescreen_metric <- match.arg(prescreen_metric)
  prep <- .coda_fs_prepare_inputs(X, y, positive = positive, input_scale = input_scale)
  candidate <- .coda_fs_prescreen(
    prep$log_x,
    prep$y,
    max_pair_features = max_pair_features,
    metric = prescreen_metric
  )
  local_x <- prep$log_x[, candidate, drop = FALSE]
  pairs <- utils::combn(ncol(local_x), 2L)
  pair_matrix <- .coda_fs_build_pair_matrix(local_x, pairs)

  n <- nrow(pair_matrix)
  sub_n <- max(4L, floor(subsample_frac * n))
  sub_n <- min(sub_n, n)
  pair_freq <- rep(0, ncol(pair_matrix))

  for (b in seq_len(as.integer(n_subsamples))) {
    idx <- sort(sample.int(n, size = sub_n, replace = FALSE))
    y_sub <- prep$y[idx]
    if (length(unique(y_sub)) < 2L) {
      next
    }
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x = pair_matrix[idx, , drop = FALSE],
        y = y_sub,
        family = "binomial",
        alpha = alpha,
        nfolds = .coda_fs_cv_folds(y_sub, nfolds),
        nlambda = as.integer(nlambda),
        standardize = TRUE,
        intercept = TRUE,
        type.measure = "deviance"
      ),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) {
      next
    }
    beta <- as.numeric(stats::coef(cv_fit, s = cv_fit[[lambda_rule]]))[-1L]
    pair_freq[abs(beta) > 0] <- pair_freq[abs(beta) > 0] + 1
  }
  pair_freq <- pair_freq / max(1L, as.integer(n_subsamples))

  pair_table <- data.frame(
    i = pairs[1, ],
    j = pairs[2, ],
    frequency = pair_freq
  )
  stable_pairs <- pair_table[pair_table$frequency >= selection_threshold, , drop = FALSE]
  if (nrow(stable_pairs) == 0L) {
    stable_pairs <- pair_table
  }

  local_scores <- .coda_fs_aggregate_pair_scores(
    pair_table = stable_pairs,
    feature_names = colnames(local_x),
    weight_col = "frequency"
  )
  scores <- setNames(rep(0, length(prep$feature_names)), prep$feature_names)
  scores[names(local_scores)] <- local_scores

  list(
    method = "coda_stability_logratio",
    scores = scores,
    selected_features = .coda_fs_select_top(scores, n_features),
    pair_table = .coda_fs_label_pair_table(stable_pairs, colnames(local_x))
  )
}


#' @title Pairwise Log-Ratio Variance Filter
#'
#' @description
#' `mlr3filters::Filter` wrapper around [codaFS_plr_variance()].
#'
#' @section Dictionary:
#' This Filter can be instantiated via the dictionary `mlr_filters` or with the
#' associated sugar function `flt()`:
#' ```
#' flt("coda_plr_variance")
#' mlr_filters$get("coda_plr_variance")
#' ```
#'
#' @export
FilterCoDA_PLRVariance <- R6::R6Class(
  "FilterCoDA_PLRVariance",
  inherit = mlr3filters::Filter,
  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "coda_plr_variance",
        task_types = "classif",
        param_set = paradox::ps(
          input_scale = paradox::p_fct(levels = c("log", "raw"), default = "log"),
          top_pairs = paradox::p_int(lower = 0L, default = 0L),
          pair_multiplier = paradox::p_int(lower = 1L, default = 5L),
          max_pair_features = paradox::p_int(lower = 2L, default = 200L),
          prescreen_metric = paradox::p_fct(
            levels = c("clr_effect", "variance"),
            default = "clr_effect"
          )
        ),
        feature_types = c("numeric", "integer"),
        packages = character(),
        label = "CoDA PLR Variance",
        man = "OmicSelector::FilterCoDA_PLRVariance"
      )
    }
  ),
  private = list(
    .calculate = function(task, nfeat) {
      pars <- self$param_set$values
      res <- codaFS_plr_variance(
        X = .coda_fs_task_matrix(task),
        y = task$truth(),
        n_features = nfeat,
        positive = task$positive,
        input_scale = pars$input_scale %||% "log",
        top_pairs = if ((pars$top_pairs %||% 0L) > 0L) pars$top_pairs else NULL,
        pair_multiplier = pars$pair_multiplier %||% 5L,
        max_pair_features = pars$max_pair_features %||% 200L,
        prescreen_metric = pars$prescreen_metric %||% "clr_effect"
      )
      res$scores
    }
  )
)


#' @title Selbal-Style Forward Balance Filter
#'
#' @description
#' `mlr3filters::Filter` wrapper around [codaFS_selbal_wrapper()].
#'
#' @section Dictionary:
#' ```
#' flt("coda_selbal")
#' mlr_filters$get("coda_selbal")
#' ```
#'
#' @export
FilterCoDA_Selbal <- R6::R6Class(
  "FilterCoDA_Selbal",
  inherit = mlr3filters::Filter,
  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "coda_selbal",
        task_types = "classif",
        param_set = paradox::ps(
          input_scale = paradox::p_fct(levels = c("log", "raw"), default = "log"),
          backend = paradox::p_fct(
            levels = c("auto", "official", "fallback"),
            default = "auto"
          ),
          max_terms = paradox::p_int(lower = 0L, default = 0L),
          min_improvement = paradox::p_dbl(lower = 0, default = 1e-4),
          max_pair_features = paradox::p_int(lower = 2L, default = 80L),
          prescreen_metric = paradox::p_fct(
            levels = c("clr_effect", "variance"),
            default = "clr_effect"
          )
        ),
        feature_types = c("numeric", "integer"),
        packages = character(),
        label = "CoDA Selbal",
        man = "OmicSelector::FilterCoDA_Selbal"
      )
    }
  ),
  private = list(
    .calculate = function(task, nfeat) {
      pars <- self$param_set$values
      res <- codaFS_selbal_wrapper(
        X = .coda_fs_task_matrix(task),
        y = task$truth(),
        n_features = nfeat,
        positive = task$positive,
        input_scale = pars$input_scale %||% "log",
        backend = pars$backend %||% "auto",
        max_terms = if ((pars$max_terms %||% 0L) > 0L) pars$max_terms else NULL,
        min_improvement = pars$min_improvement %||% 1e-4,
        max_pair_features = pars$max_pair_features %||% 80L,
        prescreen_metric = pars$prescreen_metric %||% "clr_effect"
      )
      res$scores
    }
  )
)


#' @title CoDaCoRe Feature Filter
#'
#' @description
#' `mlr3filters::Filter` wrapper around [codaFS_codacore_wrapper()].
#'
#' @section Dictionary:
#' ```
#' flt("coda_codacore")
#' mlr_filters$get("coda_codacore")
#' ```
#'
#' @export
FilterCoDA_CoDaCoRe <- R6::R6Class(
  "FilterCoDA_CoDaCoRe",
  inherit = mlr3filters::Filter,
  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "coda_codacore",
        task_types = "classif",
        param_set = paradox::ps(
          input_scale = paradox::p_fct(levels = c("log", "raw"), default = "log"),
          backend = paradox::p_fct(
            levels = c("auto", "official", "fallback"),
            default = "auto"
          ),
          max_balances = paradox::p_int(lower = 0L, default = 0L),
          top_k = paradox::p_int(lower = 1L, default = 3L),
          class_weight = paradox::p_lgl(default = FALSE),
          verbose = paradox::p_lgl(default = FALSE),
          seed = paradox::p_int(lower = 1L, default = 42L)
        ),
        feature_types = c("numeric", "integer"),
        packages = character(),
        label = "CoDA CoDaCoRe",
        man = "OmicSelector::FilterCoDA_CoDaCoRe"
      )
    }
  ),
  private = list(
    .calculate = function(task, nfeat) {
      pars <- self$param_set$values
      res <- codaFS_codacore_wrapper(
        X = .coda_fs_task_matrix(task),
        y = task$truth(),
        n_features = nfeat,
        positive = task$positive,
        input_scale = pars$input_scale %||% "log",
        backend = pars$backend %||% "auto",
        max_balances = if ((pars$max_balances %||% 0L) > 0L) pars$max_balances else NULL,
        top_k = pars$top_k %||% 3L,
        class_weight = pars$class_weight %||% FALSE,
        verbose = pars$verbose %||% FALSE,
        seed = pars$seed %||% 42L
      )
      res$scores
    }
  )
)


#' @title Log-Contrast Lasso Filter
#'
#' @description
#' `mlr3filters::Filter` wrapper around [codaFS_logcontrast_lasso()].
#'
#' @section Dictionary:
#' ```
#' flt("coda_logcontrast_lasso")
#' mlr_filters$get("coda_logcontrast_lasso")
#' ```
#'
#' @export
FilterCoDA_LogcontrastLasso <- R6::R6Class(
  "FilterCoDA_LogcontrastLasso",
  inherit = mlr3filters::Filter,
  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "coda_logcontrast_lasso",
        task_types = "classif",
        param_set = paradox::ps(
          input_scale = paradox::p_fct(levels = c("log", "raw"), default = "log"),
          alpha = paradox::p_dbl(lower = 0, upper = 1, default = 1),
          nfolds = paradox::p_int(lower = 2L, default = 5L),
          nlambda = paradox::p_int(lower = 10L, default = 100L),
          lambda_rule = paradox::p_fct(
            levels = c("lambda.1se", "lambda.min"),
            default = "lambda.1se"
          )
        ),
        feature_types = c("numeric", "integer"),
        packages = "glmnet",
        label = "CoDA Log-Contrast Lasso",
        man = "OmicSelector::FilterCoDA_LogcontrastLasso"
      )
    }
  ),
  private = list(
    .calculate = function(task, nfeat) {
      pars <- self$param_set$values
      res <- codaFS_logcontrast_lasso(
        X = .coda_fs_task_matrix(task),
        y = task$truth(),
        n_features = nfeat,
        positive = task$positive,
        input_scale = pars$input_scale %||% "log",
        alpha = pars$alpha %||% 1,
        nfolds = pars$nfolds %||% 5L,
        nlambda = pars$nlambda %||% 100L,
        lambda_rule = pars$lambda_rule %||% "lambda.1se"
      )
      res$scores
    }
  )
)


#' @title Stability Log-Ratio Filter
#'
#' @description
#' `mlr3filters::Filter` wrapper around [codaFS_stability_logratio()].
#'
#' @section Dictionary:
#' ```
#' flt("coda_stability_logratio")
#' mlr_filters$get("coda_stability_logratio")
#' ```
#'
#' @export
FilterCoDA_StabilityLogratio <- R6::R6Class(
  "FilterCoDA_StabilityLogratio",
  inherit = mlr3filters::Filter,
  public = list(
    #' @description
    #' Creates a new instance of this Filter.
    initialize = function() {
      super$initialize(
        id = "coda_stability_logratio",
        task_types = "classif",
        param_set = paradox::ps(
          input_scale = paradox::p_fct(levels = c("log", "raw"), default = "log"),
          alpha = paradox::p_dbl(lower = 0, upper = 1, default = 0.9),
          subsample_frac = paradox::p_dbl(lower = 0.5, upper = 1, default = 0.7),
          n_subsamples = paradox::p_int(lower = 2L, default = 20L),
          nfolds = paradox::p_int(lower = 2L, default = 3L),
          nlambda = paradox::p_int(lower = 10L, default = 60L),
          lambda_rule = paradox::p_fct(
            levels = c("lambda.1se", "lambda.min"),
            default = "lambda.1se"
          ),
          selection_threshold = paradox::p_dbl(lower = 0, upper = 1, default = 0.6),
          max_pair_features = paradox::p_int(lower = 2L, default = 60L),
          prescreen_metric = paradox::p_fct(
            levels = c("clr_effect", "variance"),
            default = "clr_effect"
          )
        ),
        feature_types = c("numeric", "integer"),
        packages = "glmnet",
        label = "CoDA Stability Log-Ratio",
        man = "OmicSelector::FilterCoDA_StabilityLogratio"
      )
    }
  ),
  private = list(
    .calculate = function(task, nfeat) {
      pars <- self$param_set$values
      res <- codaFS_stability_logratio(
        X = .coda_fs_task_matrix(task),
        y = task$truth(),
        n_features = nfeat,
        positive = task$positive,
        input_scale = pars$input_scale %||% "log",
        alpha = pars$alpha %||% 0.9,
        subsample_frac = pars$subsample_frac %||% 0.7,
        n_subsamples = pars$n_subsamples %||% 20L,
        nfolds = pars$nfolds %||% 3L,
        nlambda = pars$nlambda %||% 60L,
        lambda_rule = pars$lambda_rule %||% "lambda.1se",
        selection_threshold = pars$selection_threshold %||% 0.6,
        max_pair_features = pars$max_pair_features %||% 60L,
        prescreen_metric = pars$prescreen_metric %||% "clr_effect"
      )
      res$scores
    }
  )
)


#' @title Register CoDA Feature Selection Filters
#'
#' @description
#' Registers the CoDA-aware `mlr3filters` feature filters provided by this
#' package. The function is idempotent and is also called on package load.
#'
#' @export
register_coda_feature_selection_filters <- function() {
  if (!requireNamespace("mlr3filters", quietly = TRUE)) {
    stop("Package 'mlr3filters' is required.", call. = FALSE)
  }

  .paper1_register_filter("coda_plr_variance", FilterCoDA_PLRVariance)
  .paper1_register_filter("coda_selbal", FilterCoDA_Selbal)
  .paper1_register_filter("coda_codacore", FilterCoDA_CoDaCoRe)
  .paper1_register_filter("coda_logcontrast_lasso", FilterCoDA_LogcontrastLasso)
  .paper1_register_filter("coda_stability_logratio", FilterCoDA_StabilityLogratio)

  invisible(TRUE)
}


#' @keywords internal
.coda_fs_prepare_inputs <- function(X, y, positive, input_scale) {
  train <- .paper1_prepare_train_inputs(X, y, positive = positive)
  list(
    log_x = .codacore_logged_matrix(train$x, input_scale = input_scale),
    y = train$y,
    feature_names = colnames(train$x),
    positive = train$positive,
    negative = train$negative
  )
}


#' @keywords internal
.coda_fs_task_matrix <- function(task) {
  x <- as.matrix(task$data(cols = task$feature_names))
  storage.mode(x) <- "numeric"
  x
}


#' @keywords internal
.coda_fs_target_size <- function(n_features, p) {
  if (is.null(n_features) || is.na(n_features)) {
    return(min(10L, p))
  }
  max(1L, min(as.integer(n_features), p))
}


#' @keywords internal
.coda_fs_default_max_balances <- function(max_balances, n_features, p) {
  if (!is.null(max_balances) && !is.na(max_balances)) {
    return(max(1L, min(as.integer(max_balances), p)))
  }
  target <- max(1L, ceiling(.coda_fs_target_size(n_features, p) / 4))
  min(target, max(1L, floor(p / 2L)))
}


#' @keywords internal
.coda_fs_select_top <- function(scores, n_features) {
  if (is.null(n_features)) {
    return(names(scores)[scores > 0])
  }
  ord <- order(scores, decreasing = TRUE)
  head(names(scores)[ord], .coda_fs_target_size(n_features, length(scores)))
}


#' @keywords internal
.coda_fs_prescreen <- function(log_x, y, max_pair_features, metric = c("clr_effect", "variance")) {
  metric <- match.arg(metric)
  p <- ncol(log_x)
  max_pair_features <- max(2L, min(as.integer(max_pair_features), p))
  if (max_pair_features >= p) {
    return(seq_len(p))
  }

  clr_x <- clr_transform(log_x, input_scale = "log")
  if (identical(metric, "variance")) {
    score <- apply(clr_x, 2, stats::sd)
  } else {
    effect <- abs(
      colMeans(clr_x[y == 1L, , drop = FALSE]) -
        colMeans(clr_x[y == 0L, , drop = FALSE])
    )
    score <- effect * pmax(apply(clr_x, 2, stats::sd), 1e-8)
  }
  score[!is.finite(score)] <- 0
  order(score, decreasing = TRUE)[seq_len(max_pair_features)]
}


#' @keywords internal
.coda_fs_pair_summary <- function(log_x, y) {
  pairs <- utils::combn(ncol(log_x), 2L)
  centered_x <- scale(log_x, center = TRUE, scale = FALSE)
  cov_x <- crossprod(centered_x) / max(1, nrow(log_x) - 1L)
  var_x <- diag(cov_x)
  y_centered <- y - mean(y)
  cov_xy <- as.vector(crossprod(centered_x, y_centered) / max(1, nrow(log_x) - 1L))
  y_var <- stats::var(y)

  i <- pairs[1, ]
  j <- pairs[2, ]
  pair_var <- pmax(var_x[i] + var_x[j] - 2 * cov_x[cbind(i, j)], 0)
  pair_cov <- cov_xy[i] - cov_xy[j]
  pair_cor <- abs(pair_cov) / sqrt(pmax(pair_var, 1e-8) * max(y_var, 1e-8))
  score <- pair_var * pair_cor

  data.frame(
    i = i,
    j = j,
    variance = pair_var,
    effect = abs(pair_cov),
    score = score
  )
}


#' @keywords internal
.coda_fs_aggregate_pair_scores <- function(pair_table, feature_names, weight_col = "score") {
  scores <- setNames(rep(0, length(feature_names)), feature_names)
  if (nrow(pair_table) == 0L) {
    return(scores)
  }
  weights <- pair_table[[weight_col]]
  weights[!is.finite(weights)] <- 0
  for (k in seq_len(nrow(pair_table))) {
    scores[[feature_names[pair_table$i[[k]]]]] <- scores[[feature_names[pair_table$i[[k]]]]] + weights[[k]]
    scores[[feature_names[pair_table$j[[k]]]]] <- scores[[feature_names[pair_table$j[[k]]]]] + weights[[k]]
  }
  scores
}


#' @keywords internal
.coda_fs_label_pair_table <- function(pair_table, feature_names) {
  if (nrow(pair_table) == 0L) {
    return(pair_table)
  }
  out <- pair_table
  out$feature_i <- feature_names[out$i]
  out$feature_j <- feature_names[out$j]
  out
}


#' @keywords internal
.coda_fs_balance <- function(log_x, numerator, denominator) {
  s1 <- length(numerator)
  s2 <- length(denominator)
  sqrt((s1 * s2) / (s1 + s2)) * (
    rowMeans(log_x[, numerator, drop = FALSE]) -
      rowMeans(log_x[, denominator, drop = FALSE])
  )
}


#' @keywords internal
.coda_fs_binomial_objective <- function(balance, y) {
  fit <- tryCatch(
    suppressWarnings(
      stats::glm.fit(
        x = cbind(`(Intercept)` = 1, balance = balance),
        y = y,
        family = stats::binomial()
      )
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || anyNA(fit$coefficients) || !is.finite(fit$deviance)) {
    return(-Inf)
  }
  fit$null.deviance - fit$deviance
}


#' @keywords internal
.coda_fs_forward_balance <- function(log_x, y, max_terms, min_improvement) {
  pair_stats <- .coda_fs_pair_summary(log_x, y)
  best_pair <- pair_stats[which.max(pair_stats$score), , drop = FALSE]
  numerator <- best_pair$i[[1]]
  denominator <- best_pair$j[[1]]
  current_obj <- .coda_fs_binomial_objective(
    .coda_fs_balance(log_x, numerator, denominator),
    y
  )

  history <- data.frame(
    step = c(1L, 2L),
    feature = c(numerator, denominator),
    side = c("numerator", "denominator"),
    objective = c(current_obj, current_obj)
  )

  while (length(c(numerator, denominator)) < max_terms) {
    remaining <- setdiff(seq_len(ncol(log_x)), c(numerator, denominator))
    if (length(remaining) == 0L) {
      break
    }

    best_feature <- NA_integer_
    best_side <- NA_character_
    best_obj <- current_obj

    for (idx in remaining) {
      obj_num <- .coda_fs_binomial_objective(
        .coda_fs_balance(log_x, c(numerator, idx), denominator),
        y
      )
      if (obj_num > best_obj) {
        best_obj <- obj_num
        best_feature <- idx
        best_side <- "numerator"
      }

      obj_den <- .coda_fs_binomial_objective(
        .coda_fs_balance(log_x, numerator, c(denominator, idx)),
        y
      )
      if (obj_den > best_obj) {
        best_obj <- obj_den
        best_feature <- idx
        best_side <- "denominator"
      }
    }

    if (!is.finite(best_obj) || (best_obj - current_obj) <= min_improvement || is.na(best_feature)) {
      break
    }

    if (identical(best_side, "numerator")) {
      numerator <- c(numerator, best_feature)
    } else {
      denominator <- c(denominator, best_feature)
    }
    current_obj <- best_obj
    history <- rbind(
      history,
      data.frame(
        step = nrow(history) + 1L,
        feature = best_feature,
        side = best_side,
        objective = current_obj
      )
    )
  }

  list(
    numerator = numerator,
    denominator = denominator,
    history = history,
    selection_order = history$feature,
    selection_side = history$side
  )
}


#' @keywords internal
.coda_fs_build_pair_matrix <- function(log_x, pairs) {
  out <- log_x[, pairs[1, ], drop = FALSE] - log_x[, pairs[2, ], drop = FALSE]
  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1L)
  }
  colnames(out) <- paste0(
    colnames(log_x)[pairs[1, ]],
    "__lr__",
    colnames(log_x)[pairs[2, ]]
  )
  out
}


#' @keywords internal
.coda_fs_cv_folds <- function(y, nfolds) {
  min_class <- min(table(y))
  max(2L, min(as.integer(nfolds), as.integer(min_class)))
}
