#' @title Ratio Image CNN for Biomarker Panel Classification
#'
#' @description
#' Converts biomarker expression profiles into pairwise log-ratio images and
#' classifies them using a convolutional neural network. The encoding is
#' additive-shift-invariant in log-space: any per-sample additive constant
#' applied uniformly across all features cancels exactly in the pairwise
#' difference pixel(i,j) = log2(x_i) - log2(x_j). This yields a
#' reference-cohort-free representation that preserves within-sample
#' relationships while removing uniform sample-level offsets.
#'
#' @details
#' For a panel of \eqn{p} biomarkers, each sample becomes a \eqn{p \times p}
#' image where pixel(i,j) = log2(biomarker_i / biomarker_j). On log-scale data
#' this is simply the difference \eqn{x_i - x_j}. The default CNN matches the
#' Paper 1 v2.2 prototype: two 3x3 convolutions, adaptive average pooling to
#' 3x3, a 64-unit dense layer with dropout 0.3, and a sigmoid output head.
#'
#' The torch backend is used when available. The learner wrapper stores the
#' fitted torch module and applies the same training-fold normalization to new
#' samples during prediction.
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
#' @name ratio-image-cnn
NULL


#' @title Create Pairwise Log-Ratio Image
#'
#' @description
#' Converts a single sample's expression vector into a pairwise log-ratio image.
#' For log-scale inputs, each pixel is the within-sample difference
#' \eqn{x_i - x_j}.
#'
#' @param x Numeric vector of feature values on log scale.
#'
#' @return A square numeric matrix with one row and one column per feature.
#'
#' @examples
#' x <- c(5, 3, 1)
#' encode_ratio_image(x)
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
encode_ratio_image <- function(x) {
  x <- as.numeric(x)
  outer(x, x, "-")
}


#' @title Create Pairwise Ratio Image from Expression Vector
#'
#' @description
#' Alias for [encode_ratio_image()] retained for backwards compatibility.
#'
#' @inheritParams encode_ratio_image
#'
#' @return A square numeric matrix with one row and one column per feature.
#'
#' @examples
#' x <- c(5, 3, 1)
#' make_ratio_image(x)
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' @export
make_ratio_image <- function(x) {
  encode_ratio_image(x)
}


#' @title Batch-Create Ratio Images from an Expression Matrix
#'
#' @description
#' Applies [encode_ratio_image()] row-wise to a log-scale expression matrix.
#'
#' @param mat Numeric matrix with samples in rows and features in columns.
#'
#' @return A 3D array with dimensions `samples x features x features`.
#'
#' @examples
#' mat <- matrix(c(5, 3, 1, 4, 2, 0), nrow = 2, byrow = TRUE)
#' imgs <- make_ratio_images(mat)
#' dim(imgs)
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' @export
make_ratio_images <- function(mat) {
  mat <- .paper1_require_matrix(mat, arg = "mat")
  n <- nrow(mat)
  p <- ncol(mat)
  imgs <- array(0, dim = c(n, p, p))
  for (i in seq_len(n)) {
    imgs[i, , ] <- encode_ratio_image(mat[i, ])
  }
  imgs
}


#' @title Train Ratio Image CNN
#'
#' @description
#' Fits the Paper 1 v2.2 pairwise-ratio CNN on log-scale biomarker data and
#' returns predictions for new samples.
#'
#' @param X_train Numeric matrix with training samples in rows and features in
#'   columns.
#' @param y_train Binary outcome vector. Factors, characters, logical values,
#'   and integer `0/1` labels are supported.
#' @param X_test Numeric matrix with test samples in rows and features in
#'   columns.
#' @param epochs Number of training epochs. Defaults to `30`.
#' @param lr Learning rate. Defaults to `0.001`.
#' @param batch_size Mini-batch size. Defaults to `256`.
#' @param class_weight Logical; if `TRUE`, use inverse-frequency positive-class
#'   weighting. Defaults to `TRUE`.
#' @param device One of `"cpu"`, `"cuda"`, or `NULL`/`"auto"` for automatic
#'   selection.
#' @param verbose Logical; print training progress. Defaults to `TRUE`.
#' @param seed Integer random seed used for reproducible torch initialization.
#'   Defaults to `42`.
#' @param positive Optional positive-class label for non-numeric outcomes.
#'
#' @return A list with elements `predictions`, `fit`, `model`, `train_images`,
#'   and `image_stats`.
#'
#' @examples
#' \dontrun{
#' X_train <- matrix(rnorm(14 * 20), nrow = 20)
#' y_train <- factor(rep(c("control", "case"), each = 10))
#' X_test <- matrix(rnorm(14 * 5), nrow = 5)
#' fit <- train_ratio_cnn(X_train, y_train, X_test, epochs = 5, verbose = FALSE)
#' fit$predictions
#' }
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
train_ratio_cnn <- function(X_train, y_train, X_test,
                            epochs = 30L,
                            lr = 0.001,
                            batch_size = 256L,
                            class_weight = TRUE,
                            device = NULL,
                            verbose = TRUE,
                            seed = 42L,
                            positive = NULL) {
  fit <- .fit_ratio_cnn(
    X_train = X_train,
    y_train = y_train,
    epochs = epochs,
    lr = lr,
    batch_size = batch_size,
    class_weight = class_weight,
    device = device,
    verbose = verbose,
    seed = seed,
    positive = positive
  )
  predictions <- .predict_ratio_cnn_fit(fit, X_test)
  list(
    predictions = predictions,
    fit = fit,
    model = fit$model,
    train_images = fit$train_images,
    image_stats = fit$image_stats
  )
}


#' @title Train Ratio CNN Across Three Default Seeds
#'
#' @description
#' Repeats [train_ratio_cnn()] across the default Paper 1 v2.2 seeds
#' `c(42, 7, 2026)` and aggregates the predicted probabilities.
#'
#' @inheritParams train_ratio_cnn
#' @param seeds Integer vector of seeds. Defaults to the Paper 1 v2.2
#'   reproducibility triplet `c(42, 7, 2026)`.
#' @param aggregate How to aggregate per-seed predictions; either `"mean"` or
#'   `"median"`.
#'
#' @return A list with aggregated predictions, per-seed predictions, and the
#'   fitted objects.
#'
#' @examples
#' \dontrun{
#' X_train <- matrix(rnorm(14 * 20), nrow = 20)
#' y_train <- factor(rep(c("control", "case"), each = 10))
#' X_test <- matrix(rnorm(14 * 5), nrow = 5)
#' fit <- train_ratio_cnn_multiseed(X_train, y_train, X_test, verbose = FALSE)
#' fit$predictions
#' }
#'
#' @references
#' Quinn TP, Erb I, Richardson MF, Crowley TM. (2020).
#' Understanding sequencing data as compositions: an outlook and review.
#' \emph{Bioinformatics}, 36(16), 4424-4432.
#'
#' @export
train_ratio_cnn_multiseed <- function(X_train, y_train, X_test,
                                      seeds = .paper1_default_seeds(),
                                      aggregate = c("mean", "median"),
                                      epochs = 30L,
                                      lr = 0.001,
                                      batch_size = 256L,
                                      class_weight = TRUE,
                                      device = NULL,
                                      verbose = TRUE,
                                      positive = NULL) {
  aggregate <- match.arg(aggregate)
  seeds <- as.integer(seeds)
  fits <- lapply(
    seeds,
    function(seed) {
      train_ratio_cnn(
        X_train = X_train,
        y_train = y_train,
        X_test = X_test,
        epochs = epochs,
        lr = lr,
        batch_size = batch_size,
        class_weight = class_weight,
        device = device,
        verbose = verbose,
        seed = seed,
        positive = positive
      )
    }
  )

  seed_predictions <- do.call(
    cbind,
    lapply(fits, function(x) x$predictions)
  )
  colnames(seed_predictions) <- paste0("seed_", seeds)
  aggregated <- switch(
    aggregate,
    mean = rowMeans(seed_predictions),
    median = apply(seed_predictions, 1L, stats::median)
  )

  list(
    predictions = aggregated,
    seed_predictions = seed_predictions,
    fits = fits,
    seeds = seeds,
    aggregate = aggregate
  )
}


#' @title mlr3 Ratio CNN Learner
#'
#' @description
#' `mlr3` classification learner wrapping [train_ratio_cnn()] with the Paper 1
#' v2.2 CNN architecture.
#'
#' @details
#' The learner is registered on package load as `classif.ratio_cnn` so it can
#' be instantiated with `mlr3::lrn("classif.ratio_cnn")`.
#'
#' @examples
#' \dontrun{
#' learner <- mlr3::lrn("classif.ratio_cnn", epochs = 5, verbose = FALSE)
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
LearnerClassifRatioCNN <- R6::R6Class(
  "LearnerClassifRatioCNN",
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
        seed = paradox::p_int(lower = 1L, default = 42L, tags = "train")
      )

      super$initialize(
        id = "classif.ratio_cnn",
        param_set = ps,
        predict_types = c("response", "prob"),
        feature_types = c("integer", "numeric"),
        properties = "twoclass",
        packages = c("torch", "coro"),
        label = "Pairwise Ratio CNN"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pars <- self$param_set$get_values(tags = "train")
      x_train <- as.matrix(task$data(cols = task$feature_names))
      storage.mode(x_train) <- "numeric"
      self$model <- .fit_ratio_cnn(
        X_train = x_train,
        y_train = task$truth(),
        epochs = pars$epochs %||% 30L,
        lr = pars$lr %||% 0.001,
        batch_size = pars$batch_size %||% 256L,
        class_weight = pars$class_weight %||% TRUE,
        device = pars$device %||% "auto",
        verbose = pars$verbose %||% FALSE,
        seed = pars$seed %||% 42L,
        positive = task$positive
      )
      self$model$feature_names <- task$feature_names
      self$model$class_names <- task$class_names
      self$model$positive <- task$positive %||% tail(task$class_names, 1L)
      invisible(self$model)
    },
    .predict = function(task) {
      x_new <- as.matrix(task$data(cols = task$feature_names))
      storage.mode(x_new) <- "numeric"
      prob_positive <- .predict_ratio_cnn_fit(self$model, x_new)
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
.ratio_cnn_module <- function() {
  torch::nn_module(
    initialize = function() {
      self$conv1 <- torch::nn_conv2d(1, 16, kernel_size = 3, padding = 1)
      self$conv2 <- torch::nn_conv2d(16, 32, kernel_size = 3, padding = 1)
      self$pool <- torch::nn_adaptive_avg_pool2d(3)
      self$fc1 <- torch::nn_linear(32 * 9, 64)
      self$fc2 <- torch::nn_linear(64, 1)
      self$relu <- torch::nn_relu()
      self$dropout <- torch::nn_dropout(0.3)
    },
    forward = function(x) {
      x <- self$relu(self$conv1(x))
      x <- self$relu(self$conv2(x))
      x <- self$pool(x)
      x <- x$view(c(x$size(1), -1))
      x <- self$dropout(self$relu(self$fc1(x)))
      self$fc2(x)
    }
  )
}


#' @keywords internal
.fit_ratio_cnn <- function(X_train, y_train,
                           epochs = 30L,
                           lr = 0.001,
                           batch_size = 256L,
                           class_weight = TRUE,
                           device = NULL,
                           verbose = TRUE,
                           seed = 42L,
                           positive = NULL) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(
      "Package 'torch' is required for train_ratio_cnn(). ",
      "Install it with install.packages('torch') and torch::install_torch().",
      call. = FALSE
    )
  }
  if (!requireNamespace("coro", quietly = TRUE)) {
    stop(
      "Package 'coro' is required for train_ratio_cnn(). ",
      "Install it with install.packages('coro').",
      call. = FALSE
    )
  }

  train <- .paper1_prepare_train_inputs(X_train, y_train, positive = positive)
  device <- .paper1_resolve_device(device)
  .paper1_set_seed(seed)

  train_imgs_raw <- make_ratio_images(train$x)
  normalized <- .paper1_normalize_train_test(train_imgs_raw)
  train_imgs <- normalized$train

  x_train_t <- torch::torch_tensor(train_imgs, dtype = torch::torch_float())$unsqueeze(2)
  y_train_t <- torch::torch_tensor(train$y, dtype = torch::torch_float())

  if (identical(device, "cuda")) {
    x_train_t <- x_train_t$cuda()
    y_train_t <- y_train_t$cuda()
  }

  pos_weight <- if (isTRUE(class_weight)) {
    n_pos <- sum(train$y == 1L)
    n_neg <- sum(train$y == 0L)
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
      list(x = self$X[i, , , ], y = self$y[i])
    },
    .length = function() {
      self$X$size()[[1L]]
    }
  )(x_train_t, y_train_t)

  loader <- torch::dataloader(dataset, batch_size = as.integer(batch_size), shuffle = TRUE)
  network <- .ratio_cnn_module()()
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
          "ratio-cnn epoch %d/%d loss=%.4f",
          epoch,
          as.integer(epochs),
          epoch_loss
        )
      )
    }
  }

  structure(
    list(
      model = network,
      train_images = train_imgs,
      image_stats = list(mean = normalized$mean, sd = normalized$sd),
      class_names = train$class_names,
      positive = train$positive,
      feature_names = colnames(train$x),
      device = device,
      seed = as.integer(seed),
      backend = "torch"
    ),
    class = "omicselector_ratio_cnn"
  )
}


#' @keywords internal
.predict_ratio_cnn_fit <- function(fit, X_new) {
  if (!inherits(fit, "omicselector_ratio_cnn")) {
    stop("fit must be an object returned by train_ratio_cnn().", call. = FALSE)
  }
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("Package 'torch' is required for prediction with a ratio CNN fit.", call. = FALSE)
  }

  X_new <- .paper1_prepare_new_data(
    X_new,
    n_features = length(fit$feature_names %||% seq_len(dim(fit$train_images)[[2L]])),
    feature_names = fit$feature_names
  )
  test_imgs <- make_ratio_images(X_new)
  if (!is.na(fit$image_stats$sd) && fit$image_stats$sd > 0) {
    test_imgs <- (test_imgs - fit$image_stats$mean) / fit$image_stats$sd
  }

  x_test_t <- torch::torch_tensor(test_imgs, dtype = torch::torch_float())$unsqueeze(2)
  if (identical(fit$device, "cuda")) {
    x_test_t <- x_test_t$cuda()
  }

  fit$model$eval()
  preds <- torch::with_no_grad({
    logits <- fit$model(x_test_t)$squeeze()
    torch::torch_sigmoid(logits)$cpu()
  })
  as.numeric(preds)
}
