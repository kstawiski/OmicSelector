#' @title Within-Sample Perturbation Benchmark
#'
#' @description
#' Benchmarks a set of classification learners under the Paper 1 v2.2
#' perturbation conditions: baseline, per-sample additive shift, per-sample
#' scaling, and the combined perturbation. Learners are trained on the original
#' training fold and evaluated on a perturbed copy of the held-out fold.
#'
#' @param task An `mlr3::TaskClassif` with numeric features.
#' @param learners A list of `mlr3` learners or learner IDs understood by
#'   [mlr3::lrn()].
#' @param conditions Character vector of perturbation conditions. Defaults to
#'   `c("baseline", "shift", "scale", "both")`.
#' @param resampling An `mlr3` resampling object. Defaults to 5-fold CV.
#' @param shift_sd Standard deviation of the per-sample additive perturbation.
#'   Defaults to `2`.
#' @param scale_range Length-2 numeric vector defining the uniform scaling
#'   range. Defaults to `c(0.5, 2.0)`.
#' @param seed Integer random seed controlling perturbation draws.
#'
#' @return A list with class `"ws_perturbation_benchmark"` containing
#'   per-fold scores, aggregated AUC by learner and condition, and per-fold
#'   predictions.
#'
#' @examples
#' \dontrun{
#' task <- mlr3::tsk("sonar")
#' learners <- list("classif.rpart", "classif.clr_mlp")
#' ws_perturbation_benchmark(task, learners, seed = 1)
#' }
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017).
#' Microbiome datasets are compositional: and this is not optional.
#' \emph{Frontiers in Microbiology}, 8, 2224.
#'
#' Quinn TP, Erb I, Richardson MF, Crowley TM. (2020).
#' Understanding sequencing data as compositions: an outlook and review.
#' \emph{Bioinformatics}, 36(16), 4424-4432.
#'
#' @export
ws_perturbation_benchmark <- function(task,
                                      learners,
                                      conditions = c("baseline", "shift", "scale", "both"),
                                      resampling = mlr3::rsmp("cv", folds = 5L),
                                      shift_sd = 2,
                                      scale_range = c(0.5, 2.0),
                                      seed = 42L) {
  if (!inherits(task, "TaskClassif")) {
    stop("task must inherit from mlr3::TaskClassif.", call. = FALSE)
  }
  if (!all(task$feature_types$type %in% c("integer", "numeric"))) {
    stop("ws_perturbation_benchmark() requires numeric/integer features only.", call. = FALSE)
  }
  if (!length(learners)) {
    stop("learners must contain at least one learner.", call. = FALSE)
  }

  measure <- mlr3::msr("classif.auc")
  feature_names <- task$feature_names
  target_name <- task$target_names[[1L]]
  instantiated <- resampling$clone(deep = TRUE)
  instantiated$instantiate(task)

  learner_objects <- lapply(
    learners,
    function(learner) {
      if (is.character(learner) && length(learner) == 1L) {
        return(mlr3::lrn(learner, predict_type = "prob"))
      }
      if (!inherits(learner, "Learner")) {
        stop("Each learner must be an mlr3 learner or learner ID.", call. = FALSE)
      }
      learner$clone(deep = TRUE)
    }
  )

  score_rows <- list()
  pred_rows <- list()
  row_counter <- 1L

  for (fold_idx in seq_len(instantiated$iters)) {
    train_ids <- instantiated$train_set(fold_idx)
    test_ids <- instantiated$test_set(fold_idx)

    train_task <- task$clone(deep = TRUE)
    train_task$filter(train_ids)

    test_df <- task$data(
      rows = test_ids,
      cols = c(feature_names, target_name)
    )
    baseline_features <- as.matrix(test_df[, ..feature_names])
    storage.mode(baseline_features) <- "numeric"

    for (learner in learner_objects) {
      fitted <- learner$clone(deep = TRUE)
      fitted$predict_type <- "prob"
      fitted$train(train_task)

      for (condition in conditions) {
        perturbed_features <- .paper1_shift_scale(
          baseline_features,
          condition = condition,
          seed = as.integer(seed) + fold_idx,
          shift_sd = shift_sd,
          scale_range = scale_range
        )

        perturbed_df <- test_df
        perturbed_df[, (feature_names) := as.data.table(perturbed_features)]
        test_task <- mlr3::TaskClassif$new(
          id = sprintf("%s_%s_fold%s", task$id, condition, fold_idx),
          backend = perturbed_df,
          target = target_name,
          positive = task$positive %||% tail(task$class_names, 1L)
        )
        prediction <- fitted$predict(test_task)
        auc <- .paper1_auc_summary(prediction, measure)

        score_rows[[row_counter]] <- data.frame(
          learner_id = fitted$id,
          condition = condition,
          fold = fold_idx,
          auc = auc,
          stringsAsFactors = FALSE
        )

        positive <- task$positive %||% tail(task$class_names, 1L)
        prob_positive <- prediction$prob[, positive]
        pred_rows[[row_counter]] <- data.frame(
          learner_id = fitted$id,
          condition = condition,
          fold = fold_idx,
          row_id = prediction$row_ids,
          truth = as.character(prediction$truth),
          prob_positive = as.numeric(prob_positive),
          stringsAsFactors = FALSE
        )
        row_counter <- row_counter + 1L
      }
    }
  }

  score_table <- data.table::rbindlist(score_rows)
  prediction_table <- data.table::rbindlist(pred_rows)
  score_df <- as.data.frame(score_table)
  auc_mean <- stats::aggregate(
    auc ~ learner_id + condition,
    data = score_df,
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  auc_sd <- stats::aggregate(
    auc ~ learner_id + condition,
    data = score_df,
    FUN = function(z) stats::sd(z, na.rm = TRUE)
  )
  auc_folds <- stats::aggregate(
    auc ~ learner_id + condition,
    data = score_df,
    FUN = function(z) sum(!is.na(z))
  )
  names(auc_mean)[[3L]] <- "auc_mean"
  names(auc_sd)[[3L]] <- "auc_sd"
  names(auc_folds)[[3L]] <- "auc_folds"
  aggregate_auc <- Reduce(
    function(x, y) merge(x, y, by = c("learner_id", "condition"), all = TRUE),
    list(auc_mean, auc_sd, auc_folds)
  )

  structure(
    list(
      scores = score_table,
      aggregate_auc = aggregate_auc,
      predictions = prediction_table,
      conditions = conditions,
      resampling = instantiated,
      task_id = task$id
    ),
    class = "ws_perturbation_benchmark"
  )
}


#' @references
#' Quinn TP, Erb I, Richardson MF, Crowley TM. (2020).
#' Understanding sequencing data as compositions: an outlook and review.
#' \emph{Bioinformatics}, 36(16), 4424-4432.
#'
#' @export
print.ws_perturbation_benchmark <- function(x, ...) {
  cat("Within-sample perturbation benchmark\n")
  cat("Task:", x$task_id, "\n")
  print(x$aggregate_auc)
  invisible(x)
}
