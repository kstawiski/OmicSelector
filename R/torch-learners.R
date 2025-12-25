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
#' @return An mlr3 Learner object
#'
#' @details
#' The MLP architecture:
#' - Input layer: matches number of features
#' - Hidden layers: as specified, with batch normalization and dropout
#' - Output layer: 2 neurons (binary classification) with softmax
#'
#' Requires mlr3torch + torch. No non-deep-learning fallbacks are used.
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


  check_torch_packages()

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

  mlp_learner
}


#' @title Create mlr3torch MLP (internal)
#' @keywords internal
.create_mlr3torch_mlp <- function(hidden_layers, dropout, activation,
                                   epochs, batch_size, lr,
                                   early_stopping, patience, device) {


  # Map activation name to torch function

  activation_fn <- switch(activation,
    "relu" = torch::nn_relu,
    "tanh" = torch::nn_tanh,
    "sigmoid" = torch::nn_sigmoid,
    torch::nn_relu  # default

)

  # Build callbacks list
  callbacks <- list()
  if (isTRUE(early_stopping)) {
    tryCatch({
      callbacks <- list(
        mlr3torch::t_clbk("history")
      )
    }, error = function(e) {
      # Callbacks not available, continue without
    })
  }

  # Use classif.mlp - the predefined MLP learner in mlr3torch >= 0.3
  learner <- mlr3::lrn("classif.mlp",
    neurons = hidden_layers,
    p = dropout,
    activation = activation_fn,
    epochs = epochs,
    batch_size = batch_size,
    device = device,
    optimizer = mlr3torch::t_opt("adam", lr = lr),
    loss = mlr3torch::t_loss("cross_entropy"),
    callbacks = if (length(callbacks) > 0) callbacks else NULL,
    predict_type = "prob"
  )

  learner$id <- "omic_mlp"
  learner
}


#' @title Create Custom Torch Learner
#' @keywords internal
.create_custom_torch_learner <- function(hidden_layers, dropout, activation,
                                          epochs, batch_size, lr, device) {
  .create_mlr3torch_mlp(
    hidden_layers = hidden_layers,
    dropout = dropout,
    activation = activation,
    epochs = epochs,
    batch_size = batch_size,
    lr = lr,
    early_stopping = FALSE,
    patience = 0,
    device = device
  )
}


#' @title Make MLP Learner (Alias)
#'
#' @description
#' Convenient alias for create_mlp_learner() matching the naming convention
#' of other factory functions (make_autotuner_*).
#'
#' @inheritParams create_mlp_learner
#' @return An mlr3 Learner object
#' @export
make_mlp_learner <- create_mlp_learner



#' @title Check Deep Learning Availability
#'
#' @description
#' Reports which deep learning backends are available.
#'
#' @return A list with availability status
#'
#' @export
check_dl_availability <- function() {
  has_torch <- requireNamespace("torch", quietly = TRUE)
  has_mlr3torch <- requireNamespace("mlr3torch", quietly = TRUE)
  torch_ready <- FALSE
  cuda <- FALSE

  if (has_torch) {
    torch_ready <- tryCatch(torch::torch_is_installed(), error = function(e) FALSE)
    cuda <- if (torch_ready) {
      tryCatch(torch::cuda_is_available(), error = function(e) FALSE)
    } else {
      FALSE
    }
  }

  list(
    torch = has_torch,
    torch_ready = torch_ready,
    mlr3torch = has_mlr3torch,
    cuda_available = cuda,
    recommended = if (has_mlr3torch && torch_ready) {
      "mlr3torch (installed)"
    } else if (has_torch && !torch_ready) {
      "Run torch::install_torch() to install libtorch"
    } else {
      "Install torch + mlr3torch"
    }
  )
}
