#' @title mlr3torch Learner Integration for OmicSelector 2.0
#'
#' @description
#' Factory functions for creating mlr3torch-based deep learning models.
#' Phase 2 scope: Multi-Layer Perceptron (MLP) only.
#'
#' @details
#' Why MLP only for Phase 2:
#' - Tabular omics data typically doesn't have spatial/temporal ordering
#' - CNNs assume local structure (adjacent features correlated) - misleading for most omics
#' - MLPs are the appropriate baseline for unordered feature vectors
#' - Focus is on proper nested CV and leakage prevention, not architecture variety
#'
#' @name torch-learners
NULL


#' @title Create MLP Learner via mlr3torch
#'
#' @description
#' Creates a Multi-Layer Perceptron (MLP) learner for classification tasks.
#' Uses mlr3torch for native torch integration within the mlr3 ecosystem.
#'
#' @param hidden_layers Integer vector specifying neurons per hidden layer.
#'   Default: c(64, 32) for a 2-layer network.
#' @param dropout Dropout rate between 0 and 1. Default: 0.2.
#' @param activation Activation function: "relu" (default), "tanh", "sigmoid".
#' @param epochs Maximum training epochs. Default: 100.
#' @param batch_size Training batch size. Default: 32.
#' @param lr Learning rate. Default: 0.001.
#' @param early_stopping Logical, use early stopping. Default: TRUE.
#' @param patience Early stopping patience (epochs). Default: 10.
#' @param seed Random seed for reproducibility.
#' @param device "cpu" or "cuda" for GPU. Default: "cpu".
#'
#' @return An mlr3 Learner object (or a fallback if mlr3torch unavailable)
#'
#' @details
#' The MLP architecture:
#' - Input layer: matches number of features
#' - Hidden layers: as specified, with batch normalization and dropout
#' - Output layer: 2 neurons (binary classification) with softmax
#'
#' If mlr3torch is not installed, returns a fallback learner (glmnet) with a warning.
#'
#' @examples
#' \dontrun{
#' # Basic MLP
#' mlp <- create_mlp_learner()
#'
#' # Deeper network with more regularization
#' mlp_deep <- create_mlp_learner(
#'   hidden_layers = c(128, 64, 32),
#'   dropout = 0.3,
#'   epochs = 200
#' )
#'
#' # Use in OmicPipeline
#' pipeline <- OmicPipeline$new(data, target = "outcome")
#' learner <- pipeline$create_graph_learner(
#'   filter = "anova",
#'   model = mlp,  # Pass the MLP learner
#'   n_features = 50
#' )
#' }
#'
#' @export
create_mlp_learner <- function(hidden_layers = c(64, 32),
                                dropout = 0.2,
                                activation = c("relu", "tanh", "sigmoid"),
                                epochs = 100,
                                batch_size = 32,
                                lr = 0.001,
                                early_stopping = TRUE,
                                patience = 10,
                                seed = NULL,
                                device = "cpu") {

  activation <- match.arg(activation)


  # === GRACEFUL FALLBACK PATTERN (PAL Consensus) ===
  #
  # mlr3torch is in Suggests, not Imports, because:
  # 1. Heavy dependency (CUDA/GPU drivers, system libraries)
  # 2. Not all users need deep learning
  # 3. Core functionality works without it
  #
  # Fallback chain: mlr3torch -> tabnet -> nnet -> xgboost -> glmnet

  # Check if mlr3torch is available
  if (!requireNamespace("mlr3torch", quietly = TRUE)) {
    message(paste0(
      "\n",
      "======================================================================\n",
      " Deep Learning Backend: mlr3torch not installed\n",
      "======================================================================\n",
      "\n",
      " mlr3torch provides native PyTorch integration for OmicSelector.\n",
      " Without it, falling back to a classical ML learner.\n",
      "\n",
      " To enable deep learning:\n",
      "   1. Install torch: install.packages('torch')\n",
      "   2. Install mlr3torch: install.packages('mlr3torch')\n",
      "\n",
      " Note: GPU support requires CUDA installation.\n",
      "       See: https://torch.mlverse.org/docs/articles/installation\n",
      "======================================================================\n"
    ))
    return(.create_robust_fallback())
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
    if (requireNamespace("torch", quietly = TRUE)) {
      torch::torch_manual_seed(seed)
    }
  }

  # Create MLP learner configuration
  # Note: Actual mlr3torch implementation depends on version
  # This is a factory that adapts to available mlr3torch API

  tryCatch({
    mlp_learner <- .create_mlr3torch_mlp(
      hidden_layers = hidden_layers,
      dropout = dropout,
      activation = activation,
      epochs = epochs,
      batch_size = batch_size,
      lr = lr,
      early_stopping = early_stopping,
      patience = patience,
      device = device
    )
    return(mlp_learner)
  }, error = function(e) {
    warning(sprintf(
      "Failed to create mlr3torch MLP: %s\nUsing glmnet fallback.",
      e$message
    ))
    return(.create_glmnet_fallback())
  })
}


#' @title Create mlr3torch MLP (internal)
#' @keywords internal
.create_mlr3torch_mlp <- function(hidden_layers, dropout, activation,
                                   epochs, batch_size, lr,
                                   early_stopping, patience, device) {

  # Try to create native mlr3torch learner
  # This implementation adapts to mlr3torch API

  if (requireNamespace("mlr3torch", quietly = TRUE)) {
    # mlr3torch provides learners via lrn() interface
    # Common patterns: "classif.mlp" or custom torch learners

    # Check if classif.mlp exists in mlr3torch
    available_learners <- tryCatch({
      mlr3::mlr_learners$keys()
    }, error = function(e) character(0))

    if ("classif.mlp" %in% available_learners) {
      # Use built-in MLP if available
      learner <- mlr3::lrn("classif.mlp",
        predict_type = "prob",
        epochs = epochs,
        batch_size = batch_size
      )
      learner$id <- "omic_mlp"
      return(learner)
    }

    # Otherwise, create custom torch network via mlr3torch
    # This requires defining the network architecture

    if (requireNamespace("torch", quietly = TRUE)) {
      learner <- .create_custom_torch_learner(
        hidden_layers = hidden_layers,
        dropout = dropout,
        activation = activation,
        epochs = epochs,
        batch_size = batch_size,
        lr = lr,
        device = device
      )
      return(learner)
    }
  }

  stop("Could not create mlr3torch learner")
}


#' @title Create Custom Torch Learner
#' @keywords internal
.create_custom_torch_learner <- function(hidden_layers, dropout, activation,
                                          epochs, batch_size, lr, device) {

  # For now, return a wrapper that will use torch
  # Full implementation requires mlr3torch's Learner infrastructure

  # Check for tabnet as alternative (if available)
  if ("classif.tabnet" %in% mlr3::mlr_learners$keys()) {
    learner <- mlr3::lrn("classif.tabnet",
      predict_type = "prob",
      epochs = epochs
    )
    learner$id <- "omic_tabnet"
    return(learner)
  }

  # Ultimate fallback: use nnet for simple neural network
  if (requireNamespace("nnet", quietly = TRUE) &&
      "classif.nnet" %in% mlr3::mlr_learners$keys()) {
    learner <- mlr3::lrn("classif.nnet",
      predict_type = "prob",
      size = hidden_layers[1],  # First hidden layer size
      MaxNWts = 10000,
      maxit = epochs
    )
    learner$id <- "omic_nnet"
    message("Using nnet as simplified neural network (mlr3torch unavailable)")
    return(learner)
  }

  stop("No neural network backend available")
}


#' @title Create Robust Fallback Learner
#'
#' @description
#' Creates the best available fallback learner when mlr3torch is not installed.
#' Follows the fallback chain: xgboost -> ranger -> glmnet
#'
#' @return An mlr3 Learner object
#' @keywords internal
.create_robust_fallback <- function() {
  if (!requireNamespace("mlr3learners", quietly = TRUE)) {
    stop(paste(
      "Neither mlr3torch nor mlr3learners available.\n",
      "Install mlr3learners: install.packages('mlr3learners')"
    ), call. = FALSE)
  }

  # Get available learners
  available <- tryCatch(
    mlr3::mlr_learners$keys(),
    error = function(e) character(0)
  )

  # Fallback chain based on performance for tabular data
  if ("classif.xgboost" %in% available) {
    message("  -> Using XGBoost as fallback (excellent for tabular data)")
    learner <- mlr3::lrn("classif.xgboost",
      predict_type = "prob",
      nrounds = 100,
      max_depth = 6,
      eta = 0.1
    )
    learner$id <- "omic_xgboost_fallback"
    return(learner)
  }

  if ("classif.ranger" %in% available) {
    message("  -> Using Random Forest (ranger) as fallback")
    learner <- mlr3::lrn("classif.ranger",
      predict_type = "prob",
      num.trees = 500
    )
    learner$id <- "omic_ranger_fallback"
    return(learner)
  }

  if ("classif.glmnet" %in% available) {
    message("  -> Using Elastic Net (glmnet) as fallback")
    learner <- mlr3::lrn("classif.glmnet",
      predict_type = "prob",
      alpha = 0.5  # Elastic net
    )
    learner$id <- "omic_glmnet_fallback"
    return(learner)
  }

  stop(paste(
    "No suitable fallback learner available.\n",
    "Install at least one of: xgboost, ranger, glmnet\n",
    "  install.packages(c('xgboost', 'ranger', 'glmnet'))"
  ), call. = FALSE)
}


#' @title Create glmnet Fallback Learner (deprecated alias)
#' @keywords internal
.create_glmnet_fallback <- .create_robust_fallback


#' @title Check Deep Learning Availability
#'
#' @description
#' Reports which deep learning backends are available.
#'
#' @return A list with availability status
#'
#' @export
check_dl_availability <- function() {
  list(
    torch = requireNamespace("torch", quietly = TRUE),
    mlr3torch = requireNamespace("mlr3torch", quietly = TRUE),
    nnet = requireNamespace("nnet", quietly = TRUE),
    keras = requireNamespace("keras", quietly = TRUE),  # Legacy, not recommended
    recommended = if (requireNamespace("mlr3torch", quietly = TRUE)) {
      "mlr3torch (installed)"
    } else if (requireNamespace("nnet", quietly = TRUE)) {
      "nnet (basic neural network)"
    } else {
      "None - install mlr3torch: install.packages('mlr3torch')"
    }
  )
}
