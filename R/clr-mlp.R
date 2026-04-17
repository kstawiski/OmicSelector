#' @title CLR + MLP Models for Within-Sample Classification
#'
#' @description
#' Implements centered log-ratio (CLR) transforms and dense neural-network
#' classifiers for within-sample biomarker panel classification. On log-scale
#' inputs the CLR transform is a simple within-sample centering operation,
#' making it invariant to uniform additive shifts applied to all features of a
#' sample.
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C. (2003).
#' Isometric logratio transformations for compositional data analysis.
#' \emph{Mathematical Geology}, 35(3), 279-300.
#'
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017).
#' Microbiome datasets are compositional: and this is not optional.
#' \emph{Frontiers in Microbiology}, 8, 2224.
#'
#' Quinn TP, Erb I, Richardson MF, Crowley TM. (2020).
#' Understanding sequencing data as compositions: an outlook and review.
#' \emph{Bioinformatics}, 36(16), 4424-4432.
#'
#' @name clr-mlp
NULL


#' @title Centered Log-Ratio Transformation
#'
#' @description
#' Applies the centered log-ratio transform row-wise. For log-scale inputs this
#' subtracts the within-sample mean from each feature. For raw positive inputs
#' the function applies `log(x + pseudocount)` before centering.
#'
#' @param x Numeric vector or matrix with samples in rows and features in
#'   columns.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param pseudocount Small positive constant added to raw inputs before
#'   logging. Ignored when `input_scale = "log"`.
#'
#' @return A numeric vector or matrix of the same dimensions as `x`.
#'
#' @examples
#' clr_transform(c(5, 3, 1))
#' clr_transform(matrix(c(5, 3, 1, 4, 2, 0), nrow = 2, byrow = TRUE))
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C. (2003).
#' Isometric logratio transformations for compositional data analysis.
#' \emph{Mathematical Geology}, 35(3), 279-300.
#'
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017).
#' Microbiome datasets are compositional: and this is not optional.
#' \emph{Frontiers in Microbiology}, 8, 2224.
#'
#' @export
clr_transform <- function(x, input_scale = c("log", "raw"), pseudocount = 1e-8) {
  input_scale <- match.arg(input_scale)
  was_vector <- is.vector(x)
  x <- .paper1_require_matrix(x, arg = "x")

  if (identical(input_scale, "raw")) {
    if (any(x <= 0, na.rm = TRUE)) {
      stop("Raw CLR inputs must be strictly positive.", call. = FALSE)
    }
    x <- log(x + pseudocount)
  }

  centered <- x - rowMeans(x)
  if (was_vector) {
    return(as.numeric(centered[1, ]))
  }
  centered
}


#' @title Train a CLR + MLP Classifier
#'
#' @description
#' Fits the Paper 1 v2.2 dense CoDA baseline: a centered log-ratio transform
#' followed by a multilayer perceptron. The default backend uses torch when it
#' is available, and otherwise falls back to `nnet::nnet()` or logistic
#' regression for environments without torch.
#'
#' @param X_train Numeric training matrix with samples in rows and features in
#'   columns.
#' @param y_train Binary outcome vector.
#' @param epochs Number of training epochs. Defaults to `30`.
#' @param lr Learning rate. Defaults to `0.001`.
#' @param batch_size Mini-batch size. Defaults to `256`.
#' @param class_weight Logical; if `TRUE`, use inverse-frequency weighting.
#' @param device One of `"cpu"`, `"cuda"`, or `NULL`/`"auto"`.
#' @param verbose Logical; print training progress. Defaults to `TRUE`.
#' @param seed Integer random seed. Defaults to `42`.
#' @param positive Optional positive-class label for non-numeric outcomes.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param backend One of `"auto"`, `"torch"`, `"nnet"`, or `"glm"`.
#' @param hidden_units Integer vector with the two hidden-layer sizes used by
#'   the torch backend. Defaults to `c(152, 136)`.
#' @param dropout Dropout rate used by the torch backend. Defaults to `0.3`.
#'
#' @return An object of class `"omicselector_clr_mlp"` suitable for
#'   [predict_clr_mlp()].
#'
#' @examples
#' X_train <- matrix(rnorm(14 * 20), nrow = 20)
#' y_train <- factor(rep(c("control", "case"), each = 10))
#' fit <- train_clr_mlp(X_train, y_train, backend = "glm", verbose = FALSE)
#' fit$backend
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C. (2003).
#' Isometric logratio transformations for compositional data analysis.
#' \emph{Mathematical Geology}, 35(3), 279-300.
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
train_clr_mlp <- function(X_train, y_train,
                          epochs = 30L,
                          lr = 0.001,
                          batch_size = 256L,
                          class_weight = TRUE,
                          device = NULL,
                          verbose = TRUE,
                          seed = 42L,
                          positive = NULL,
                          input_scale = c("log", "raw"),
                          backend = c("auto", "torch", "nnet", "glm"),
                          hidden_units = c(152L, 136L),
                          dropout = 0.3) {
  input_scale <- match.arg(input_scale)
  backend <- match.arg(backend)
  train <- .paper1_prepare_train_inputs(X_train, y_train, positive = positive)
  device <- .paper1_resolve_device(device)

  if (identical(backend, "auto")) {
    if (requireNamespace("torch", quietly = TRUE)) {
      backend <- "torch"
    } else if (requireNamespace("nnet", quietly = TRUE)) {
      backend <- "nnet"
    } else {
      backend <- "glm"
    }
  }

  clr_train <- clr_transform(train$x, input_scale = input_scale)
  normalized <- .paper1_normalize_train_test(clr_train)

  fit <- switch(
    backend,
    torch = .fit_clr_mlp_torch(
      x_train = normalized$train,
      y_train = train$y,
      epochs = epochs,
      lr = lr,
      batch_size = batch_size,
      class_weight = class_weight,
      device = device,
      verbose = verbose,
      seed = seed,
      hidden_units = hidden_units,
      dropout = dropout
    ),
    nnet = .fit_clr_mlp_nnet(
      x_train = normalized$train,
      y_train = train$y,
      class_names = train$class_names,
      positive = train$positive,
      epochs = epochs,
      class_weight = class_weight,
      verbose = verbose
    ),
    glm = .fit_clr_mlp_glm(
      x_train = normalized$train,
      y_train = train$y,
      class_weight = class_weight
    )
  )

  structure(
    c(
      fit,
      list(
        class_names = train$class_names,
        positive = train$positive,
        feature_names = colnames(train$x),
        clr_stats = list(mean = normalized$mean, sd = normalized$sd),
        input_scale = input_scale,
        seed = as.integer(seed),
        hidden_units = as.integer(hidden_units),
        dropout = dropout
      )
    ),
    class = "omicselector_clr_mlp"
  )
}


#' @title Predict with a CLR + MLP Classifier
#'
#' @description
#' Generates positive-class probabilities for new samples from a fit created by
#' [train_clr_mlp()].
#'
#' @param fit A fitted `"omicselector_clr_mlp"` object.
#' @param X_new Numeric matrix of new samples.
#'
#' @return A numeric vector of positive-class probabilities.
#'
#' @examples
#' X_train <- matrix(rnorm(14 * 20), nrow = 20)
#' y_train <- factor(rep(c("control", "case"), each = 10))
#' fit <- train_clr_mlp(X_train, y_train, backend = "glm", verbose = FALSE)
#' predict_clr_mlp(fit, X_train[1:3, , drop = FALSE])
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C. (2003).
#' Isometric logratio transformations for compositional data analysis.
#' \emph{Mathematical Geology}, 35(3), 279-300.
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
predict_clr_mlp <- function(fit, X_new) {
  if (!inherits(fit, "omicselector_clr_mlp")) {
    stop("fit must be an object returned by train_clr_mlp().", call. = FALSE)
  }

  X_new <- .paper1_prepare_new_data(
    X_new,
    n_features = length(fit$feature_names),
    feature_names = fit$feature_names
  )
  clr_new <- clr_transform(X_new, input_scale = fit$input_scale)
  if (!is.na(fit$clr_stats$sd) && fit$clr_stats$sd > 0) {
    clr_new <- (clr_new - fit$clr_stats$mean) / fit$clr_stats$sd
  }

  switch(
    fit$backend,
    torch = .predict_clr_mlp_torch(fit, clr_new),
    nnet = .predict_clr_mlp_nnet(fit, clr_new),
    glm = .predict_clr_mlp_glm(fit, clr_new)
  )
}


#' @title mlr3 CLR + MLP Learner
#'
#' @description
#' `mlr3` learner wrapping [train_clr_mlp()] and [predict_clr_mlp()]. It is
#' registered on package load as `classif.clr_mlp`.
#'
#' @examples
#' \dontrun{
#' learner <- mlr3::lrn("classif.clr_mlp", backend = "glm", verbose = FALSE)
#' learner
#' }
#'
#' @references
#' Quinn TP, Erb I, Richardson MF, Crowley TM. (2020).
#' Understanding sequencing data as compositions: an outlook and review.
#' \emph{Bioinformatics}, 36(16), 4424-4432.
#'
#' @keywords internal
#' @noRd
LearnerClassifCLRMLP <- R6::R6Class(
  "LearnerClassifCLRMLP",
  inherit = mlr3::LearnerClassif,
  public = list(
    initialize = function() {
      ps <- paradox::ps(
        epochs = paradox::p_int(lower = 1L, default = 30L, tags = "train"),
        lr = paradox::p_dbl(lower = 1e-5, upper = 1, default = 0.001, tags = "train"),
        batch_size = paradox::p_int(lower = 1L, default = 256L, tags = "train"),
        class_weight = paradox::p_lgl(default = TRUE, tags = "train"),
        device = paradox::p_fct(
          levels = c("auto", "cpu", "cuda"),
          default = "auto",
          tags = "train"
        ),
        verbose = paradox::p_lgl(default = FALSE, tags = "train"),
        seed = paradox::p_int(lower = 1L, default = 42L, tags = "train"),
        input_scale = paradox::p_fct(
          levels = c("log", "raw"),
          default = "log",
          tags = "train"
        ),
        backend = paradox::p_fct(
          levels = c("auto", "torch", "nnet", "glm"),
          default = "auto",
          tags = "train"
        ),
        dropout = paradox::p_dbl(lower = 0, upper = 1, default = 0.3, tags = "train")
      )

      super$initialize(
        id = "classif.clr_mlp",
        param_set = ps,
        predict_types = c("response", "prob"),
        feature_types = c("integer", "numeric"),
        properties = "twoclass",
        packages = "stats",
        label = "CLR + MLP"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pars <- self$param_set$get_values(tags = "train")
      x_train <- as.matrix(task$data(cols = task$feature_names))
      storage.mode(x_train) <- "numeric"
      self$model <- train_clr_mlp(
        X_train = x_train,
        y_train = task$truth(),
        epochs = pars$epochs %||% 30L,
        lr = pars$lr %||% 0.001,
        batch_size = pars$batch_size %||% 256L,
        class_weight = pars$class_weight %||% TRUE,
        device = pars$device %||% "auto",
        verbose = pars$verbose %||% FALSE,
        seed = pars$seed %||% 42L,
        positive = task$positive,
        input_scale = pars$input_scale %||% "log",
        backend = pars$backend %||% "auto",
        hidden_units = c(152L, 136L),
        dropout = pars$dropout %||% 0.3
      )
      invisible(self$model)
    },
    .predict = function(task) {
      x_new <- as.matrix(task$data(cols = task$feature_names))
      storage.mode(x_new) <- "numeric"
      prob_positive <- predict_clr_mlp(self$model, x_new)
      prob <- .paper1_prob_matrix(
        prob_positive = prob_positive,
        class_names = self$model$class_names,
        positive = self$model$positive
      )
      response <- .paper1_response_from_prob(
        prob_positive = prob_positive,
        class_names = self$model$class_names,
        positive = self$model$positive
      )
      mlr3::PredictionClassif$new(
        task = task,
        row_ids = task$row_ids,
        truth = task$truth(),
        prob = prob,
        response = factor(response, levels = self$model$class_names)
      )
    }
  )
)


#' @keywords internal
.clr_mlp_torch_module <- function(n_features,
                                  hidden_units = c(152L, 136L),
                                  dropout = 0.3) {
  torch::nn_module(
    initialize = function() {
      self$fc1 <- torch::nn_linear(n_features, hidden_units[[1L]])
      self$fc2 <- torch::nn_linear(hidden_units[[1L]], hidden_units[[2L]])
      self$fc3 <- torch::nn_linear(hidden_units[[2L]], 1)
      self$relu <- torch::nn_relu()
      self$dropout <- torch::nn_dropout(dropout)
    },
    forward = function(x) {
      x <- self$dropout(self$relu(self$fc1(x)))
      x <- self$dropout(self$relu(self$fc2(x)))
      self$fc3(x)
    }
  )
}


#' @keywords internal
.fit_clr_mlp_torch <- function(x_train, y_train, epochs, lr, batch_size,
                               class_weight, device, verbose, seed,
                               hidden_units, dropout) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("torch backend requested but package 'torch' is not installed.", call. = FALSE)
  }
  if (!requireNamespace("coro", quietly = TRUE)) {
    stop("torch backend requested but package 'coro' is not installed.", call. = FALSE)
  }

  .paper1_set_seed(seed)
  x_train_t <- torch::torch_tensor(x_train, dtype = torch::torch_float())
  y_train_t <- torch::torch_tensor(y_train, dtype = torch::torch_float())

  if (identical(device, "cuda")) {
    x_train_t <- x_train_t$cuda()
    y_train_t <- y_train_t$cuda()
  }

  pos_weight <- if (isTRUE(class_weight)) {
    n_pos <- sum(y_train == 1L)
    n_neg <- sum(y_train == 0L)
    torch::torch_tensor(n_neg / max(n_pos, 1L), dtype = torch::torch_float())
  } else {
    torch::torch_tensor(1, dtype = torch::torch_float())
  }
  if (identical(device, "cuda")) {
    pos_weight <- pos_weight$cuda()
  }

  dataset <- torch::dataset(
    initialize = function(X, y) {
      self$X <- X
      self$y <- y
    },
    .getitem = function(i) {
      list(x = self$X[i, ], y = self$y[i])
    },
    .length = function() {
      self$X$size()[[1L]]
    }
  )(x_train_t, y_train_t)

  loader <- torch::dataloader(dataset, batch_size = as.integer(batch_size), shuffle = TRUE)
  network <- .clr_mlp_torch_module(
    n_features = ncol(x_train),
    hidden_units = hidden_units,
    dropout = dropout
  )()
  if (identical(device, "cuda")) {
    network <- network$cuda()
  }

  optimizer <- torch::optim_adam(network$parameters, lr = lr, weight_decay = 1e-4)
  criterion <- torch::nn_bce_with_logits_loss(pos_weight = pos_weight)

  network$train()
  for (epoch in seq_len(as.integer(epochs))) {
    epoch_loss <- 0
    coro::loop(for (batch in loader) {
      optimizer$zero_grad()
      logits <- network(batch$x)$squeeze()
      loss <- criterion(logits, batch$y)
      loss$backward()
      optimizer$step()
      epoch_loss <- epoch_loss + loss$item()
    })
    if (isTRUE(verbose) && (epoch == 1L || epoch == epochs || epoch %% 10L == 0L)) {
      message(
        sprintf(
          "clr-mlp epoch %d/%d loss=%.4f",
          epoch,
          as.integer(epochs),
          epoch_loss
        )
      )
    }
  }

  list(
    backend = "torch",
    model = network,
    device = device
  )
}


#' @keywords internal
.fit_clr_mlp_nnet <- function(x_train, y_train, class_names, positive, epochs,
                              class_weight, verbose) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("nnet backend requested but package 'nnet' is not installed.", call. = FALSE)
  }

  negative <- setdiff(class_names, positive)[[1L]]
  outcome <- factor(
    ifelse(y_train == 1L, positive, negative),
    levels = class_names
  )
  train_df <- as.data.frame(x_train)
  train_df$.target <- outcome

  obs_weights <- rep(1, length(y_train))
  if (isTRUE(class_weight)) {
    n_pos <- sum(y_train == 1L)
    n_neg <- sum(y_train == 0L)
    obs_weights[y_train == 1L] <- n_neg / max(n_pos, 1L)
  }

  fit <- nnet::nnet(
    .target ~ .,
    data = train_df,
    size = min(152L, max(2L, ncol(x_train))),
    decay = 1e-4,
    maxit = as.integer(epochs),
    trace = isTRUE(verbose),
    weights = obs_weights,
    MaxNWts = max(10000L, 10L * ncol(x_train))
  )

  list(
    backend = "nnet",
    model = fit
  )
}


#' @keywords internal
.fit_clr_mlp_glm <- function(x_train, y_train, class_weight) {
  train_df <- as.data.frame(x_train)
  train_df$.target <- y_train
  obs_weights <- rep(1, length(y_train))
  if (isTRUE(class_weight)) {
    n_pos <- sum(y_train == 1L)
    n_neg <- sum(y_train == 0L)
    obs_weights[y_train == 1L] <- n_neg / max(n_pos, 1L)
  }

  fit <- stats::glm(
    .target ~ .,
    data = train_df,
    family = stats::binomial(),
    weights = obs_weights
  )

  list(
    backend = "glm",
    model = fit
  )
}


#' @keywords internal
.predict_clr_mlp_torch <- function(fit, x_new) {
  x_new_t <- torch::torch_tensor(x_new, dtype = torch::torch_float())
  if (identical(fit$device, "cuda")) {
    x_new_t <- x_new_t$cuda()
  }
  fit$model$eval()
  preds <- torch::with_no_grad({
    logits <- fit$model(x_new_t)$squeeze()
    torch::torch_sigmoid(logits)$cpu()
  })
  as.numeric(preds)
}


#' @keywords internal
.predict_clr_mlp_nnet <- function(fit, x_new) {
  raw <- stats::predict(fit$model, newdata = as.data.frame(x_new), type = "raw")
  if (is.matrix(raw)) {
    return(as.numeric(raw[, fit$positive]))
  }
  as.numeric(raw)
}


#' @keywords internal
.predict_clr_mlp_glm <- function(fit, x_new) {
  as.numeric(stats::predict(fit$model, newdata = as.data.frame(x_new), type = "response"))
}
