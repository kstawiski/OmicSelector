#' @title mlr3torch Checkpoint Utilities
#'
#' @description
#' Export, import, and fine-tune mlr3torch learners using torch checkpoints.
#' These helpers operate on the underlying torch module inside a trained learner.
#'
#' @name torch-checkpoint
NULL


.require_torch_mlr3torch <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("Package 'torch' is required. Install with torch::install_torch().",
         call. = FALSE)
  }
  if (!requireNamespace("mlr3torch", quietly = TRUE)) {
    stop("Package 'mlr3torch' is required. Install with install.packages('mlr3torch').",
         call. = FALSE)
  }
}


.resolve_base_learner <- function(learner) {
  if (inherits(learner, "GraphLearner")) {
    pipeops <- learner$graph$pipeops
    if ("learner" %in% names(pipeops)) {
      return(pipeops$learner$learner)
    }
    if (requireNamespace("mlr3pipelines", quietly = TRUE)) {
      for (po in pipeops) {
        if (inherits(po, "PipeOpLearner")) {
          return(po$learner)
        }
      }
    }
  }
  learner
}


.extract_torch_model <- function(learner) {
  if (!inherits(learner, "Learner")) {
    stop("Expected an mlr3 Learner.", call. = FALSE)
  }

 # Check if learner is trained
  # mlr3torch uses state$model presence, not state$trained
  is_trained <- !is.null(learner$state) &&
    (!is.null(learner$state$model) || isTRUE(learner$state$trained))

  if (!is_trained) {
    stop("Learner must be trained before exporting a checkpoint.", call. = FALSE)
  }

  model <- learner$model
  if (inherits(model, "nn_module")) {
    return(model)
  }
  if (is.list(model)) {
    # mlr3torch 0.3+ stores network in model$network
    for (key in c("network", "model", "module", "net")) {
      if (!is.null(model[[key]]) && inherits(model[[key]], "nn_module")) {
        return(model[[key]])
      }
    }
  }

  stop("Could not locate torch model inside learner$model.", call. = FALSE)
}


.load_state_into_learner <- function(learner, state, strict = TRUE) {
  if (!inherits(learner, "Learner")) {
    stop("Expected an mlr3 Learner.", call. = FALSE)
  }

  # Check if learner is trained (mlr3torch uses state$model, not state$trained)
  is_trained <- !is.null(learner$state) &&
    (!is.null(learner$state$model) || isTRUE(learner$state$trained))

  if (!is_trained) {
    return(FALSE)
  }

  model <- .extract_torch_model(learner)
  model$load_state_dict(state, strict = strict)
  TRUE
}


.attach_pretrained_state <- function(learner, state) {
  ids <- learner$param_set$ids()
  if ("pretrained_state" %in% ids) {
    learner$param_set$values$pretrained_state <- state
    return(TRUE)
  }
  if ("init_state" %in% ids) {
    learner$param_set$values$init_state <- state
    return(TRUE)
  }
  if ("state_dict" %in% ids) {
    learner$param_set$values$state_dict <- state
    return(TRUE)
  }

  attr(learner, "omic_pretrained_state") <- state
  FALSE
}


#' @title Export mlr3torch Checkpoint
#'
#' @description
#' Saves a trained mlr3torch learner's torch module state to disk.
#'
#' @param learner Trained mlr3 learner or GraphLearner.
#' @param path Output file path (recommended: .pt).
#' @param include_optimizer Logical; attempts to save optimizer state if available.
#' @param metadata Optional list of metadata to store alongside weights.
#'
#' @export
export_mlr3torch_checkpoint <- function(learner,
                                        path,
                                        include_optimizer = FALSE,
                                        metadata = NULL) {
  .require_torch_mlr3torch()

  base <- .resolve_base_learner(learner)
  model <- .extract_torch_model(base)

  state <- list(
    model_state = model$state_dict(),
    meta = list(
      learner_id = base$id,
      learner_class = class(base),
      params = base$param_set$values,
      timestamp = Sys.time()
    )
  )

  if (!is.null(metadata)) {
    state$meta$user_metadata <- metadata
  }

  if (isTRUE(include_optimizer) && is.list(base$model)) {
    opt <- base$model$optimizer
    if (!is.null(opt) && is.function(opt$state_dict)) {
      state$optimizer_state <- opt$state_dict()
    }
  }

  torch::torch_save(state, path)
  invisible(path)
}


#' @title Import mlr3torch Checkpoint
#'
#' @description
#' Loads a torch checkpoint into a trained learner for inference, or attaches
#' it for potential fine-tuning if the learner is not yet trained.
#'
#' @param learner mlr3 learner or GraphLearner.
#' @param path Path to a checkpoint created by \code{export_mlr3torch_checkpoint()}.
#' @param strict Logical, enforce exact layer matching.
#'
#' @return The updated learner.
#' @export
import_mlr3torch_checkpoint <- function(learner, path, strict = TRUE) {
  .require_torch_mlr3torch()

  if (!file.exists(path)) {
    stop("Checkpoint file not found.", call. = FALSE)
  }

  state <- torch::torch_load(path)
  base <- .resolve_base_learner(learner)

  if (!is.list(state) || is.null(state$model_state)) {
    stop("Invalid checkpoint format: missing model_state.", call. = FALSE)
  }

  loaded <- FALSE
  # Check if learner is trained (mlr3torch uses state$model, not state$trained)
  is_trained <- !is.null(base$state) &&
    (!is.null(base$state$model) || isTRUE(base$state$trained))

  if (is_trained) {
    loaded <- .load_state_into_learner(base, state$model_state, strict = strict)
  }

  if (!loaded) {
    .attach_pretrained_state(base, state$model_state)
  }

  learner
}


#' @title Fine-tune from mlr3torch Checkpoint
#'
#' @description
#' Loads a checkpoint and trains the learner on a new task, attempting to
#' warm-start when supported by the mlr3torch learner.
#'
#' @param learner mlr3 learner or GraphLearner.
#' @param task mlr3 Task to train on.
#' @param path Checkpoint path.
#' @param epochs Optional number of epochs to set before training.
#' @param seed Optional seed.
#' @param device Optional device ("cpu" or "cuda") if supported by learner.
#' @param strict Logical, enforce exact layer matching.
#'
#' @return Trained learner.
#' @export
finetune_mlr3torch_checkpoint <- function(learner,
                                          task,
                                          path,
                                          epochs = NULL,
                                          seed = NULL,
                                          device = NULL,
                                          strict = TRUE) {
  .require_torch_mlr3torch()

  base <- .resolve_base_learner(learner)
  if (!inherits(task, "Task")) {
    stop("task must be an mlr3 Task.", call. = FALSE)
  }

  learner <- import_mlr3torch_checkpoint(learner, path, strict = strict)

  if (!is.null(epochs) && "epochs" %in% base$param_set$ids()) {
    base$param_set$values$epochs <- epochs
  }
  if (!is.null(device) && "device" %in% base$param_set$ids()) {
    base$param_set$values$device <- device
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  base$train(task)
  learner
}


#' @title Export mlr3torch Checkpoint from OmicFit
#'
#' @description
#' Convenience wrapper that exports a checkpoint from an OmicFit object.
#'
#' @param fit OmicFit object returned by \code{OmicPipeline$fit()}.
#' @param path Output file path (recommended: .pt).
#' @param include_optimizer Logical; attempts to save optimizer state if available.
#' @param metadata Optional list of metadata to store alongside weights.
#'
#' @export
export_omicfit_checkpoint <- function(fit,
                                      path,
                                      include_optimizer = FALSE,
                                      metadata = NULL) {
  if (!inherits(fit, "OmicFit")) {
    stop("fit must be an OmicFit object.", call. = FALSE)
  }
  export_mlr3torch_checkpoint(
    learner = fit$learner,
    path = path,
    include_optimizer = include_optimizer,
    metadata = metadata
  )
}


#' @title Import mlr3torch Checkpoint into OmicFit
#'
#' @description
#' Convenience wrapper that loads a checkpoint into the learner inside an OmicFit object.
#'
#' @param fit OmicFit object returned by \code{OmicPipeline$fit()}.
#' @param path Checkpoint path created by \code{export_mlr3torch_checkpoint()}.
#' @param strict Logical, enforce exact layer matching.
#'
#' @return Updated OmicFit object.
#' @export
import_omicfit_checkpoint <- function(fit, path, strict = TRUE) {
  if (!inherits(fit, "OmicFit")) {
    stop("fit must be an OmicFit object.", call. = FALSE)
  }
  fit$learner <- import_mlr3torch_checkpoint(
    learner = fit$learner,
    path = path,
    strict = strict
  )
  fit
}


#' @title Fine-tune OmicFit from mlr3torch Checkpoint
#'
#' @description
#' Convenience wrapper that fine-tunes an OmicFit learner and returns
#' a new OmicFit object.
#'
#' @param fit OmicFit object returned by \code{OmicPipeline$fit()}.
#' @param task mlr3 Task to train on.
#' @param path Checkpoint path.
#' @param epochs Optional number of epochs to set before training.
#' @param seed Optional seed.
#' @param device Optional device ("cpu" or "cuda") if supported by learner.
#' @param strict Logical, enforce exact layer matching.
#'
#' @return A new OmicFit object.
#' @export
finetune_omicfit_checkpoint <- function(fit,
                                        task,
                                        path,
                                        epochs = NULL,
                                        seed = NULL,
                                        device = NULL,
                                        strict = TRUE) {
  if (!inherits(fit, "OmicFit")) {
    stop("fit must be an OmicFit object.", call. = FALSE)
  }
  if (!inherits(task, "Task")) {
    stop("task must be an mlr3 Task.", call. = FALSE)
  }

  finetuned <- finetune_mlr3torch_checkpoint(
    learner = fit$learner,
    task = task,
    path = path,
    epochs = epochs,
    seed = seed,
    device = device,
    strict = strict
  )

  result <- fit
  result$learner <- finetuned
  result
}
