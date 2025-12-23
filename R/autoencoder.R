#' @title Autoencoder Utilities (torch)
#'
#' @description
#' Train and use a torch autoencoder for representation learning on omics data.
#' Provides a PipeOp for leakage-safe integration with mlr3pipelines.
#'
#' @name autoencoder
NULL


check_torch_available <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(
      "Package 'torch' is required for autoencoder support.\n",
      "Install with:\n",
      "  install.packages('torch')\n",
      "  torch::install_torch()",
      call. = FALSE
    )
  }
}


.torch_device <- function(device) {
  if (identical(device, "cuda")) {
    if (!torch::cuda_is_available()) {
      warning("CUDA requested but not available. Falling back to CPU.", call. = FALSE)
      return(torch::torch_device("cpu"))
    }
  }
  torch::torch_device(device)
}


.activation_layer <- function(activation) {
  if (identical(activation, "gelu") && !exists("nn_gelu", asNamespace("torch"), inherits = FALSE)) {
    activation <- "relu"
  }
  switch(activation,
    relu = torch::nn_relu(),
    tanh = torch::nn_tanh(),
    sigmoid = torch::nn_sigmoid(),
    gelu = torch::nn_gelu(),
    torch::nn_relu()
  )
}


.build_autoencoder <- function(input_dim,
                               hidden_layers,
                               latent_dim,
                               activation,
                               dropout) {
  dims <- c(input_dim, hidden_layers, latent_dim)
  enc_layers <- list()
  for (i in seq_len(length(dims) - 1)) {
    enc_layers <- c(enc_layers, list(torch::nn_linear(dims[i], dims[i + 1])))
    if (i < length(dims) - 1) {
      enc_layers <- c(enc_layers, list(.activation_layer(activation)))
      if (dropout > 0) {
        enc_layers <- c(enc_layers, list(torch::nn_dropout(dropout)))
      }
    }
  }
  encoder <- do.call(torch::nn_sequential, enc_layers)

  dec_dims <- c(latent_dim, rev(hidden_layers), input_dim)
  dec_layers <- list()
  for (i in seq_len(length(dec_dims) - 1)) {
    dec_layers <- c(dec_layers, list(torch::nn_linear(dec_dims[i], dec_dims[i + 1])))
    if (i < length(dec_dims) - 1) {
      dec_layers <- c(dec_layers, list(.activation_layer(activation)))
      if (dropout > 0) {
        dec_layers <- c(dec_layers, list(torch::nn_dropout(dropout)))
      }
    }
  }
  decoder <- do.call(torch::nn_sequential, dec_layers)

  Autoencoder <- torch::nn_module(
    initialize = function(encoder, decoder) {
      self$encoder <- encoder
      self$decoder <- decoder
    },
    forward = function(x) {
      z <- self$encoder(x)
      recon <- self$decoder(z)
      list(recon = recon, z = z)
    }
  )

  Autoencoder(encoder = encoder, decoder = decoder)
}


.tensor_to_matrix <- function(x) {
  as.matrix(x$to(device = "cpu")$detach())
}


#' @title Fit Autoencoder
#'
#' @description
#' Fits a torch autoencoder and returns an OmicAutoencoder object.
#'
#' @param x Numeric matrix or data.frame (samples x features).
#' @param latent_dim Integer, size of latent representation.
#' @param hidden_layers Integer vector of hidden layer sizes.
#' @param activation Activation function: "relu", "tanh", "sigmoid", "gelu".
#' @param dropout Dropout rate between 0 and 1.
#' @param epochs Number of training epochs.
#' @param batch_size Batch size.
#' @param lr Learning rate.
#' @param weight_decay L2 weight decay.
#' @param loss Loss function: "mse" or "mae".
#' @param validation_split Fraction of data for validation.
#' @param early_stopping Logical, enable early stopping on validation loss.
#' @param patience Early stopping patience.
#' @param min_delta Minimum improvement to reset patience.
#' @param device "cpu" or "cuda".
#' @param seed Optional random seed.
#' @param pretrained Optional path to a saved autoencoder state (torch).
#' @param freeze_encoder Logical, freeze encoder weights during training.
#'
#' @return An OmicAutoencoder object.
#' @export
autoencoder_fit <- function(x,
                            latent_dim = 32L,
                            hidden_layers = c(128L, 64L),
                            activation = c("relu", "tanh", "sigmoid", "gelu"),
                            dropout = 0.1,
                            epochs = 100L,
                            batch_size = 32L,
                            lr = 1e-3,
                            weight_decay = 0,
                            loss = c("mse", "mae"),
                            validation_split = 0.1,
                            early_stopping = TRUE,
                            patience = 10L,
                            min_delta = 0,
                            device = "cpu",
                            seed = NULL,
                            pretrained = NULL,
                            freeze_encoder = FALSE) {

  check_torch_available()

  activation <- match.arg(activation)
  loss <- match.arg(loss)

  x <- as.matrix(x)
  if (anyNA(x)) {
    stop("autoencoder_fit: input contains NA. Impute before training.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
    torch::torch_manual_seed(seed)
  }

  input_dim <- ncol(x)
  latent_dim <- as.integer(latent_dim)
  if (latent_dim <= 0) {
    stop("latent_dim must be > 0", call. = FALSE)
  }
  if (latent_dim >= input_dim) {
    warning("latent_dim >= input_dim; reducing latent_dim to input_dim - 1.", call. = FALSE)
    latent_dim <- max(1L, input_dim - 1L)
  }

  if (!is.null(pretrained)) {
    loaded <- autoencoder_load(pretrained, device = device)
    model <- loaded$model
    config <- loaded$config
    if (!identical(as.integer(config$latent_dim), latent_dim)) {
      stop("pretrained latent_dim does not match requested latent_dim.", call. = FALSE)
    }
    if (!identical(config$hidden_layers, as.integer(hidden_layers))) {
      stop("pretrained hidden_layers do not match requested hidden_layers.", call. = FALSE)
    }
  } else {
    model <- .build_autoencoder(
      input_dim = input_dim,
      hidden_layers = as.integer(hidden_layers),
      latent_dim = latent_dim,
      activation = activation,
      dropout = dropout
    )
    config <- list(
      input_dim = input_dim,
      hidden_layers = as.integer(hidden_layers),
      latent_dim = latent_dim,
      activation = activation,
      dropout = dropout,
      loss = loss
    )
  }

  device <- .torch_device(device)
  model$to(device = device)

  if (freeze_encoder) {
    for (param in model$encoder$parameters) {
      param$requires_grad_(FALSE)
    }
  }

  x_tensor <- torch::torch_tensor(x, dtype = torch::torch_float())
  x_tensor <- x_tensor$to(device = device)

  n <- nrow(x)
  val_n <- floor(n * validation_split)
  if (val_n < 1) {
    val_n <- 0
  }
  idx <- seq_len(n)
  if (val_n > 0) {
    val_idx <- sample(idx, val_n)
    train_idx <- setdiff(idx, val_idx)
  } else {
    train_idx <- idx
    val_idx <- integer(0)
  }

  criterion <- switch(
    loss,
    mse = torch::nn_mse_loss(),
    mae = torch::nn_l1_loss()
  )

  opt <- torch::optim_adam(
    params = model$parameters,
    lr = lr,
    weight_decay = weight_decay
  )

  history <- list(train_loss = numeric(), val_loss = numeric())
  best_loss <- Inf
  best_state <- NULL
  patience_left <- patience

  for (epoch in seq_len(epochs)) {
    model$train()
    epoch_loss <- 0

    train_idx <- sample(train_idx)
    for (i in seq(1, length(train_idx), by = batch_size)) {
      batch_ids <- train_idx[i:min(i + batch_size - 1, length(train_idx))]
      xb <- x_tensor[batch_ids, ]
      opt$zero_grad()
      out <- model(xb)
      loss_val <- criterion(out$recon, xb)
      loss_val$backward()
      opt$step()
      epoch_loss <- epoch_loss + loss_val$item()
    }

    epoch_loss <- epoch_loss / max(1, ceiling(length(train_idx) / batch_size))
    history$train_loss <- c(history$train_loss, epoch_loss)

    if (val_n > 0) {
      model$eval()
      torch::with_no_grad({
        xb <- x_tensor[val_idx, ]
        out <- model(xb)
        val_loss <- criterion(out$recon, xb)$item()
        history$val_loss <- c(history$val_loss, val_loss)
      })
    } else {
      val_loss <- epoch_loss
      history$val_loss <- c(history$val_loss, val_loss)
    }

    if (isTRUE(early_stopping)) {
      if ((best_loss - val_loss) > min_delta) {
        best_loss <- val_loss
        best_state <- model$state_dict()
        patience_left <- patience
      } else {
        patience_left <- patience_left - 1
        if (patience_left <= 0) {
          break
        }
      }
    }
  }

  if (!is.null(best_state)) {
    model$load_state_dict(best_state)
  }

  result <- list(
    model = model,
    config = config,
    history = history
  )
  class(result) <- c("OmicAutoencoder", "list")
  result
}


#' @title Encode Features with Autoencoder
#'
#' @description
#' Applies a fitted autoencoder to encode features into latent space.
#'
#' @param model OmicAutoencoder object.
#' @param x Numeric matrix or data.frame (samples x features).
#' @param concat Logical, if TRUE returns original + latent features.
#' @param prefix Prefix for latent feature names.
#'
#' @return A data.table with encoded features.
#' @export
autoencoder_encode <- function(model,
                               x,
                               concat = FALSE,
                               prefix = "ae_") {
  check_torch_available()
  if (!inherits(model, "OmicAutoencoder")) {
    stop("model must be an OmicAutoencoder object.", call. = FALSE)
  }

  x <- as.matrix(x)
  if (anyNA(x)) {
    stop("autoencoder_encode: input contains NA. Impute before encoding.", call. = FALSE)
  }

  device <- model$model$parameters[[1]]$device
  x_tensor <- torch::torch_tensor(x, dtype = torch::torch_float())$to(device = device)

  model$model$eval()
  torch::with_no_grad({
    out <- model$model(x_tensor)
    z <- out$z
  })

  z_mat <- .tensor_to_matrix(z)
  z_dt <- data.table::as.data.table(z_mat)
  colnames(z_dt) <- paste0(prefix, seq_len(ncol(z_dt)))

  if (isTRUE(concat)) {
    x_dt <- data.table::as.data.table(x)
    x_names <- colnames(x)
    if (is.null(x_names)) {
      x_names <- paste0("x", seq_len(ncol(x_dt)))
    }
    data.table::setnames(x_dt, x_names)
    return(cbind(x_dt, z_dt))
  }

  z_dt
}


#' @title Save Autoencoder State
#'
#' @param model OmicAutoencoder object.
#' @param path File path to save state (.pt recommended).
#'
#' @export
autoencoder_save <- function(model, path) {
  check_torch_available()
  if (!inherits(model, "OmicAutoencoder")) {
    stop("model must be an OmicAutoencoder object.", call. = FALSE)
  }

  state <- list(
    encoder_state = model$model$encoder$state_dict(),
    decoder_state = model$model$decoder$state_dict(),
    config = model$config
  )
  torch::torch_save(state, path)
  invisible(path)
}


#' @title Load Autoencoder State
#'
#' @param path Path to a saved autoencoder state.
#' @param device "cpu" or "cuda".
#'
#' @return An OmicAutoencoder object.
#' @export
autoencoder_load <- function(path, device = "cpu") {
  check_torch_available()
  if (!file.exists(path)) {
    stop("autoencoder_load: file not found.", call. = FALSE)
  }

  state <- torch::torch_load(path)
  config <- state$config

  model <- .build_autoencoder(
    input_dim = config$input_dim,
    hidden_layers = config$hidden_layers,
    latent_dim = config$latent_dim,
    activation = config$activation,
    dropout = config$dropout
  )
  model$encoder$load_state_dict(state$encoder_state)
  model$decoder$load_state_dict(state$decoder_state)
  model$to(device = .torch_device(device))

  result <- list(
    model = model,
    config = config,
    history = NULL
  )
  class(result) <- c("OmicAutoencoder", "list")
  result
}


#' @title Create Autoencoder PipeOp
#'
#' @description
#' Creates an mlr3pipelines PipeOp that trains a torch autoencoder on the
#' training data and replaces features with latent representations.
#'
#' @param latent_dim Integer, size of latent representation.
#' @param hidden_layers Integer vector of hidden layer sizes.
#' @param activation Activation function: "relu", "tanh", "sigmoid", "gelu".
#' @param dropout Dropout rate between 0 and 1.
#' @param epochs Number of training epochs.
#' @param batch_size Batch size.
#' @param lr Learning rate.
#' @param weight_decay L2 weight decay.
#' @param loss Loss function: "mse" or "mae".
#' @param validation_split Fraction of data for validation.
#' @param early_stopping Logical, enable early stopping.
#' @param patience Early stopping patience.
#' @param min_delta Minimum improvement to reset patience.
#' @param device "cpu" or "cuda".
#' @param seed Optional random seed.
#' @param pretrained Optional path to a saved autoencoder state.
#' @param freeze_encoder Logical, freeze encoder weights during training.
#' @param concat Logical, if TRUE returns original + latent features.
#' @param prefix Prefix for latent feature names.
#'
#' @return A PipeOp object.
#' @export
create_autoencoder_pipeop <- function(latent_dim = 32L,
                                      hidden_layers = c(128L, 64L),
                                      activation = "relu",
                                      dropout = 0.1,
                                      epochs = 100L,
                                      batch_size = 32L,
                                      lr = 1e-3,
                                      weight_decay = 0,
                                      loss = "mse",
                                      validation_split = 0.1,
                                      early_stopping = TRUE,
                                      patience = 10L,
                                      min_delta = 0,
                                      device = "cpu",
                                      seed = NULL,
                                      pretrained = NULL,
                                      freeze_encoder = FALSE,
                                      concat = FALSE,
                                      prefix = "ae_") {

  if (!requireNamespace("mlr3pipelines", quietly = TRUE)) {
    stop("mlr3pipelines is required for create_autoencoder_pipeop().",
         call. = FALSE)
  }

  PipeOpAutoencoder <- R6::R6Class(
    "PipeOpAutoencoder",
    inherit = mlr3pipelines::PipeOpTaskPreproc,

    public = list(
      initialize = function(id = "autoencoder", param_vals = list()) {
        ps <- paradox::ps(
          latent_dim = paradox::p_int(lower = 1L, tags = "train"),
          hidden_layers = paradox::p_uty(tags = "train"),
          activation = paradox::p_fct(levels = c("relu", "tanh", "sigmoid", "gelu"), tags = "train"),
          dropout = paradox::p_dbl(lower = 0, upper = 1, tags = "train"),
          epochs = paradox::p_int(lower = 1L, tags = "train"),
          batch_size = paradox::p_int(lower = 1L, tags = "train"),
          lr = paradox::p_dbl(lower = 1e-6, upper = 1, tags = "train"),
          weight_decay = paradox::p_dbl(lower = 0, tags = "train"),
          loss = paradox::p_fct(levels = c("mse", "mae"), tags = "train"),
          validation_split = paradox::p_dbl(lower = 0, upper = 0.5, tags = "train"),
          early_stopping = paradox::p_lgl(default = TRUE, tags = "train"),
          patience = paradox::p_int(lower = 1L, tags = "train"),
          min_delta = paradox::p_dbl(lower = 0, tags = "train"),
          device = paradox::p_fct(levels = c("cpu", "cuda"), tags = "train"),
          seed = paradox::p_uty(tags = "train"),
          pretrained = paradox::p_uty(tags = "train"),
          freeze_encoder = paradox::p_lgl(default = FALSE, tags = "train"),
          concat = paradox::p_lgl(default = FALSE, tags = "train"),
          prefix = paradox::p_uty(tags = "train")
        )

        ps$values <- list(
          latent_dim = latent_dim,
          hidden_layers = hidden_layers,
          activation = activation,
          dropout = dropout,
          epochs = epochs,
          batch_size = batch_size,
          lr = lr,
          weight_decay = weight_decay,
          loss = loss,
          validation_split = validation_split,
          early_stopping = early_stopping,
          patience = patience,
          min_delta = min_delta,
          device = device,
          seed = seed,
          pretrained = pretrained,
          freeze_encoder = freeze_encoder,
          concat = concat,
          prefix = prefix
        )

        super$initialize(
          id = id,
          param_set = ps,
          param_vals = param_vals,
          feature_types = c("numeric", "integer")
        )
      },

      save_state = function(path) {
        if (is.null(private$.model)) {
          stop("No trained autoencoder in PipeOp. Train first.", call. = FALSE)
        }
        autoencoder_save(private$.model, path)
      }
    ),

    private = list(
      .model = NULL,

      .train_dt = function(dt, levels, target) {
        model <- autoencoder_fit(
          x = as.matrix(dt),
          latent_dim = self$param_set$values$latent_dim,
          hidden_layers = self$param_set$values$hidden_layers,
          activation = self$param_set$values$activation,
          dropout = self$param_set$values$dropout,
          epochs = self$param_set$values$epochs,
          batch_size = self$param_set$values$batch_size,
          lr = self$param_set$values$lr,
          weight_decay = self$param_set$values$weight_decay,
          loss = self$param_set$values$loss,
          validation_split = self$param_set$values$validation_split,
          early_stopping = self$param_set$values$early_stopping,
          patience = self$param_set$values$patience,
          min_delta = self$param_set$values$min_delta,
          device = self$param_set$values$device,
          seed = self$param_set$values$seed,
          pretrained = self$param_set$values$pretrained,
          freeze_encoder = self$param_set$values$freeze_encoder
        )
        private$.model <- model
        autoencoder_encode(
          model = model,
          x = as.matrix(dt),
          concat = self$param_set$values$concat,
          prefix = self$param_set$values$prefix
        )
      },

      .predict_dt = function(dt, levels) {
        if (is.null(private$.model)) {
          return(dt)
        }
        autoencoder_encode(
          model = private$.model,
          x = as.matrix(dt),
          concat = self$param_set$values$concat,
          prefix = self$param_set$values$prefix
        )
      }
    )
  )

  PipeOpAutoencoder$new()
}
