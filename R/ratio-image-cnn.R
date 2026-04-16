#' @title Ratio Image CNN for Biomarker Panel Classification
#'
#' @description
#' Converts biomarker expression profiles into pairwise log-ratio images and
#' classifies them using a CNN. This approach is inherently within-sample
#' normalized because pairwise ratios cancel multiplicative batch effects.
#'
#' @details
#' The analogy to histopathology image analysis is exact:
#' \itemize{
#'   \item Batch effects in biomarkers ≈ brightness/contrast variation in images
#'   \item Pairwise log-ratios ≈ color normalization (relative channel intensities)
#'   \item CNN on ratio images ≈ learning spatial patterns invariant to staining
#' }
#'
#' For a panel of p biomarkers, each sample becomes a p×p image where
#' pixel(i,j) = log2(biomarker_i / biomarker_j). On log-scale data, this
#' is simply the difference. The CNN then learns non-linear patterns in the
#' pairwise relationship structure.
#'
#' @references
#' Stawiski K et al. (2026). Pairwise ratio images for batch-effect-free
#' biomarker classification. (in preparation)
#'
#' Sharma A et al. (2019). DeepInsight: A methodology to transform a
#' non-image data to an image for convolution neural network architecture.
#' Scientific Reports, 9(1), 11399.
#'
#' @name ratio-image-cnn
NULL


#' @title Create Pairwise Ratio Image from Expression Vector
#'
#' @description
#' Converts a single sample's expression vector into a p×p pairwise
#' log-ratio matrix. If the input is on log scale, pixel(i,j) = x[i] - x[j]
#' which equals log(expr_i / expr_j).
#'
#' @param x Numeric vector of length p (one sample's expression values, log-scale)
#'
#' @return A p×p numeric matrix (the ratio image)
#'
#' @export
make_ratio_image <- function(x) {
  p <- length(x)
  outer(x, x, "-")
}


#' @title Batch-Create Ratio Images from Expression Matrix
#'
#' @description
#' Creates ratio images for all samples in an expression matrix.
#'
#' @param mat Numeric matrix (samples × features), log-scale
#'
#' @return A 3D array (samples × features × features)
#'
#' @export
make_ratio_images <- function(mat) {
  n <- nrow(mat)
  p <- ncol(mat)
  imgs <- array(0, dim = c(n, p, p))
  for (i in seq_len(n)) {
    imgs[i, , ] <- outer(mat[i, ], mat[i, ], "-")
  }
  imgs
}


#' @title Train Ratio Image CNN for Binary Classification
#'
#' @description
#' Trains a lightweight CNN on pairwise ratio images for binary classification
#' (e.g., cancer vs healthy). Uses torch for GPU acceleration.
#'
#' @param X_train Numeric matrix (N_train × p features), log-scale
#' @param y_train Integer/numeric vector (0/1 labels)
#' @param X_test Numeric matrix (N_test × p features), log-scale
#' @param epochs Number of training epochs (default: 30)
#' @param lr Learning rate (default: 0.001)
#' @param batch_size Batch size (default: 256)
#' @param class_weight Whether to apply inverse-frequency class weighting (default: TRUE)
#' @param device "cuda" or "cpu" (default: auto-detect)
#' @param verbose Print progress (default: TRUE)
#'
#' @return A list with:
#'   \itemize{
#'     \item predictions: numeric vector of predicted probabilities for test set
#'     \item model: trained torch model
#'     \item train_images: ratio images used for training (3D array)
#'     \item image_stats: mean and sd used for image normalization
#'   }
#'
#' @examples
#' \dontrun{
#' result <- train_ratio_cnn(X_train, y_train, X_test, epochs = 30)
#' auc <- pROC::auc(pROC::roc(y_test, result$predictions))
#' }
#'
#' @export
train_ratio_cnn <- function(X_train, y_train, X_test,
                             epochs = 30L, lr = 0.001,
                             batch_size = 256L,
                             class_weight = TRUE,
                             device = NULL,
                             verbose = TRUE) {

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("Package 'torch' is required for CNN training. Install with: install.packages('torch')")
  }

  # Auto-detect device
  if (is.null(device)) {
    device <- if (torch::cuda_is_available()) "cuda" else "cpu"
  }
  if (verbose) message("Training ratio-image CNN on: ", device)

  # Create ratio images
  train_imgs <- make_ratio_images(X_train)
  test_imgs <- make_ratio_images(X_test)

  # Normalize images using training statistics
  img_mean <- mean(train_imgs)
  img_sd <- sd(as.vector(train_imgs))
  if (img_sd > 0) {
    train_imgs <- (train_imgs - img_mean) / img_sd
    test_imgs <- (test_imgs - img_mean) / img_sd
  }

  p <- ncol(X_train)

  # Convert to torch tensors (N, 1, H, W)
  X_tr <- torch::torch_tensor(train_imgs, dtype = torch::torch_float())$unsqueeze(2)
  y_tr <- torch::torch_tensor(as.numeric(y_train), dtype = torch::torch_float())
  X_te <- torch::torch_tensor(test_imgs, dtype = torch::torch_float())$unsqueeze(2)

  if (device == "cuda") {
    X_tr <- X_tr$cuda()
    y_tr <- y_tr$cuda()
    X_te <- X_te$cuda()
  }

  # Class weights
  pos_weight <- if (class_weight) {
    n_pos <- sum(y_train == 1)
    n_neg <- sum(y_train == 0)
    torch::torch_tensor(n_neg / max(n_pos, 1), dtype = torch::torch_float())
  } else {
    torch::torch_tensor(1.0)
  }
  if (device == "cuda") pos_weight <- pos_weight$cuda()

  # Define CNN model
  model <- torch::nn_module(
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

  net <- model()
  if (device == "cuda") net <- net$cuda()

  optimizer <- torch::optim_adam(net$parameters, lr = lr, weight_decay = 1e-4)
  criterion <- torch::nn_bce_with_logits_loss(pos_weight = pos_weight)

  # Training loop
  dataset <- torch::dataset(
    initialize = function(X, y) {
      self$X <- X
      self$y <- y
    },
    .getitem = function(i) {
      list(x = self$X[i, , , ], y = self$y[i])
    },
    .length = function() {
      self$X$size(1)
    }
  )(X_tr, y_tr)

  loader <- torch::dataloader(dataset, batch_size = batch_size, shuffle = TRUE)

  net$train()
  for (epoch in seq_len(epochs)) {
    epoch_loss <- 0
    coro::loop(for (batch in loader) {
      optimizer$zero_grad()
      output <- net(batch$x)$squeeze()
      loss <- criterion(output, batch$y)
      loss$backward()
      optimizer$step()
      epoch_loss <- epoch_loss + loss$item()
    })
    if (verbose && epoch %% 10 == 0) {
      message(sprintf("  Epoch %d/%d, loss: %.4f", epoch, epochs, epoch_loss))
    }
  }

  # Prediction
  net$eval()
  torch::with_no_grad({
    logits <- net(X_te)$squeeze()
    preds <- torch::torch_sigmoid(logits)$cpu()$to(dtype = torch::torch_float())
  })

  list(
    predictions = as.numeric(preds),
    model = net,
    train_images = train_imgs,
    image_stats = list(mean = img_mean, sd = img_sd)
  )
}
