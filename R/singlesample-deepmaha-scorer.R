#' @title Deep-feature class-conditional Mahalanobis LRT (learned frozen embedding)
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing the
#' class-conditional Mahalanobis confidence score of Lee et al. (2018), adapted to
#' compositional specimens. A small BatchNorm-free MLP embedding of the
#' per-sample-rCLR specimen is learned on the TRAINING set (via reticulate-python
#' torch), then \strong{FROZEN}; per-class Gaussian means and a single TIED
#' covariance are estimated in the embedding space, and a new specimen is scored by
#' the class-conditional Mahalanobis log-likelihood ratio in that frozen embedding.
#' Larger score = more case-like.
#'
#' \strong{Torch is confined to FIT.} The MLP is trained with torch (Adam, CE loss)
#' on the rCLR'd, universe-aligned training matrix. After training, the weight
#' matrices and biases up to and INCLUDING the penultimate embedding layer are
#' EXPORTED to R as plain numeric matrices/vectors; the final 2-logit
#' classification head (not needed for scoring) and the python module are
#' DISCARDED. The fitted model therefore holds NO external pointer and NO python
#' dependency at score time: the embedding forward, the Gaussian statistics, and
#' the Mahalanobis LRT are all pure base-R. This makes the default
#' \code{model_digest} a stable deterministic snapshot and lets the package suite
#' run the scorer on any host without the venv.
#'
#' \strong{Fit/score consistency despite float32.} Torch trains in float32; an R
#' float64 forward on the exported (float32-valued) weights differs from the torch
#' penultimate embedding by \eqn{\sim 10^{-6}}. To guarantee that the frozen
#' Gaussian statistics live in EXACTLY the same embedding the scorer evaluates, the
#' class means \eqn{\mu_{\mathrm{case}}, \mu_{\mathrm{ctrl}}} and the tied
#' covariance \eqn{\Sigma} are computed from the PURE-R forward of the training rows
#' (\code{.lrt_deepmaha_forward}), NOT from the torch embeddings. Fit and score thus
#' use the identical computation (the same lesson as nc-ecod-copod /
#' nc-sinkhorn-single: fit and score must share one code path).
#'
#' \strong{The embedding forward.} The exported net is
#' \eqn{\mathrm{input}(p) \to [\mathrm{Linear} \to \mathrm{ReLU}]_{1}
#'   \to \dots \to [\mathrm{Linear} \to \mathrm{ReLU}]_{L-1}
#'   \to \mathrm{Linear} \to \mathrm{embedding}(d)}. There are \code{length(hidden)}
#' linear layers up to the embedding (the last \code{hidden} width is the embedding
#' dim \eqn{d}); the activation (\code{"relu"} default, or \code{"tanh"}) is applied
#' AFTER every linear layer EXCEPT the final embedding layer (no activation follows
#' the embedding -- it is the raw penultimate representation, matching the torch
#' net). \code{.lrt_deepmaha_forward(M, weights)} reimplements this exactly:
#' successive \eqn{M \gets \phi(M W_l^\top + b_l)} for the hidden layers and a final
#' linear map to the \eqn{d}-dim embedding.
#'
#' \strong{Class-conditional Mahalanobis LRT.} With the tied pooled within-class
#' covariance \eqn{\Sigma} (ridge-regularised) and class means
#' \eqn{\mu_{\mathrm{case}}, \mu_{\mathrm{ctrl}}}, the score of an embedded specimen
#' \eqn{z} is
#' \deqn{s(z) = -\tfrac{1}{2}(z-\mu_{\mathrm{case}})^\top \Sigma^{-1}
#'   (z-\mu_{\mathrm{case}}) + \tfrac{1}{2}(z-\mu_{\mathrm{ctrl}})^\top \Sigma^{-1}
#'   (z-\mu_{\mathrm{ctrl}}),}
#' the class-conditional Gaussian log-likelihood ratio (case minus control). With a
#' single SHARED covariance this is exactly LINEAR in \eqn{z} (LDA in embedding
#' space) -- the faithful Lee-2018 estimator, which also sidesteps \eqn{p \gg n}
#' since \eqn{d} is small. Larger = more case-like.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly
#' invariant to per-specimen scaling and uses no cross-row statistic. Universe
#' features absent from a specimen carry the neutral rCLR value \code{0}. The SAME
#' rCLR transform is applied to the frozen training rows at fit and to every query
#' at score, so a specimen's score depends only on that specimen and the frozen
#' model -- it passes \code{\link{singlesample_assert_row_equivariant}}.
#'
#' \strong{Degenerate-neutral.} A query returns the neutral score \code{0} if fewer
#' than \code{hp$min_features} universe features overlap the query columns (a
#' column-overlap floor, batch-independent), or it has empty positive support over
#' the frozen universe tested on the (pre-rCLR) aligned ORIGINAL abundances
#' (\code{!any(X_use[i, ] > 0)}, where \code{X_use} is the query aligned to the
#' frozen universe with absent features set to 0). A FLAT all-equal-positive
#' composition maps to the rCLR origin but is a VALID specimen and is scored
#' normally (NOT floored) -- the rCLR-origin trap (the empty-support test is on the
#' ORIGINAL abundances, NOT on \code{all(z == 0)}).
#'
#' @references
#' Lee K, Lee K, Lee H, Shin J. (2018) A Simple Unified Framework for Detecting
#' Out-of-Distribution Samples and Adversarial Attacks. \emph{Advances in Neural
#' Information Processing Systems (NeurIPS)} 31. arXiv:1807.03888.
#'
#' @name singlesample-deepmaha
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .sinkhorn_single_rclr /
# .tabdpt_rclr.
.lrt_deepmaha_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Row-wise rCLR of a matrix. rbind (not t(vapply)) keeps the n x p shape even
# when n == 1 or p == 1.
.lrt_deepmaha_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .lrt_deepmaha_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0. Returns
# an n x p matrix over `feat_order`. Mirrors the sinkhorn/tabdpt present-subset
# alignment.
.lrt_deepmaha_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R embedding forward (exact reimplementation of the torch net up to d)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported MLP up to (and including) the penultimate
# embedding layer. `weights` is a list of per-layer lists list(W = out x in, b =
# out): layers 1..(L-1) are hidden (Linear then `activation`), layer L is the
# embedding Linear with NO activation (matches the torch net exactly). For a layer
# with weight W (out x in) and bias b, torch computes M %*% t(W) + b (row = sample);
# the activation is applied to the hidden layers only. Returns an n x d embedding.
.lrt_deepmaha_forward <- function(M, weights, activation) {
  act <- if (identical(activation, "tanh")) tanh else function(x) pmax(x, 0)
  L <- length(weights)
  H <- M
  for (l in seq_len(L)) {
    W <- weights[[l]]$W                       # out x in
    b <- weights[[l]]$b                       # length out
    H <- H %*% t(W)                           # n x out
    H <- sweep(H, 2L, b, `+`)
    if (l < L) H <- act(H)                    # no activation after the embedding
  }
  H
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.lrt_deepmaha_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_lrt_deepmaha: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_lrt_deepmaha: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_lrt_deepmaha: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("hidden", "activation", "epochs", "lr", "weight_decay",
               "batch_size", "ridge", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_lrt_deepmaha: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  hidden <- hp$hidden
  if (is.null(hidden)) hidden <- c(64L, 32L)
  if (!is.numeric(hidden) || length(hidden) < 1L || any(!is.finite(hidden)) ||
      any(hidden < 1L) || any(hidden > .Machine$integer.max) ||
      any(hidden != as.integer(hidden))) {
    stop("fit_lrt_deepmaha: hp$hidden must be a vector of positive integers")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_lrt_deepmaha: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_lrt_deepmaha: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_lrt_deepmaha: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_lrt_deepmaha: hp$weight_decay must be a non-negative finite number")
  }

  batch_size <- hp$batch_size
  if (is.null(batch_size)) batch_size <- 32L
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      !is.finite(batch_size) || batch_size < 1L ||
      batch_size > .Machine$integer.max ||
      batch_size != as.integer(batch_size)) {
    stop("fit_lrt_deepmaha: hp$batch_size must be a positive integer")
  }

  ridge <- hp$ridge
  if (is.null(ridge)) ridge <- 1e-3
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) ||
      ridge <= 0) {
    stop("fit_lrt_deepmaha: hp$ridge must be a positive finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_lrt_deepmaha: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_lrt_deepmaha: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_lrt_deepmaha: hp$seed must be a single integer")
  }

  list(
    hidden = as.integer(hidden),
    activation = activation,
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    batch_size = as.integer(batch_size),
    ridge = as.numeric(ridge),
    min_features = as.integer(min_features),
    device = device,
    seed = as.integer(seed)
  )
}


# ----------------------------------------------------------------------------
# Python / torch availability + device resolution
# ----------------------------------------------------------------------------

# Clear, actionable error if reticulate or the python torch module is
# unavailable. Plain torch training needs NO brew-OpenSSL bridge and NO HF
# download (no checkpoint, no _ssl), so only reticulate + `import torch` are
# probed -- unlike the foundation-model scorers (tabdpt/tabpfn) that pull _ssl.
.lrt_deepmaha_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_lrt_deepmaha: package 'reticulate' is required to train the ",
         "embedding. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_lrt_deepmaha: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable, so the method never errors
# on a CPU-only host. The DEFAULT is 'cpu' (exact + bit-reproducible for a tiny
# MLP).
.lrt_deepmaha_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}

# Train the BN-free MLP via torch and return the exported R-side weights (a list
# of per-layer list(W, b), only up to and INCLUDING the embedding layer; the final
# 2-logit head is discarded). The python module is used here only; nothing python
# is returned. The torch RNG (CPU + guarded CUDA) and the global R .Random.seed
# are saved and restored, so fit leaves all global RNG byte-unchanged. `Z` is the
# n x p rCLR'd, universe-aligned training matrix; `y` the integer 0/1 labels.
.lrt_deepmaha_train_export <- function(Z, y, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn

  # ---- save the torch RNG (CPU + guarded CUDA), restore on exit -------------
  # The global R .Random.seed is saved/restored by the caller fit_lrt_deepmaha
  # (captured BEFORE any reticulate import, since the first python init itself
  # seeds R's RNG). Here we only guard the torch generators that THIS fit mutates.
  torch_state <- tryCatch(torch$random$get_rng_state(), error = function(e) NULL)
  cuda_ok <- tryCatch(isTRUE(torch$cuda$is_available()), error = function(e) FALSE)
  cuda_state <- if (cuda_ok) {
    tryCatch(torch$cuda$get_rng_state_all(), error = function(e) NULL)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(torch_state)) {
      tryCatch(torch$random$set_rng_state(torch_state), error = function(e) NULL)
    }
    if (!is.null(cuda_state)) {
      tryCatch(torch$cuda$set_rng_state_all(cuda_state), error = function(e) NULL)
    }
  }, add = TRUE)

  torch$manual_seed(as.integer(hp$seed))

  np <- reticulate::import("numpy", delay_load = FALSE)
  # Cast explicitly to float32 / int64 on the tensor (np$asarray -> as_tensor can
  # round-trip the array through R doubles, so pin the dtype with $float()/$long()).
  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)
  y_py <- torch$as_tensor(np$asarray(as.integer(y), dtype = "int64"))$long()$to(device)

  p <- ncol(Z)
  hidden <- hp$hidden
  # Build input(p) -> [Linear -> act]xlen(hidden) -> Linear -> 2 logits. The last
  # `hidden` width is the embedding dim d; the embedding layer is the final hidden
  # Linear (NO activation after it). The 2-logit head trains the embedding but is
  # discarded after fit.
  act_mod <- if (identical(hp$activation, "tanh")) nn$Tanh else nn$ReLU
  layers <- list()
  in_dim <- p
  for (h in hidden) {
    layers[[length(layers) + 1L]] <- nn$Linear(in_dim, as.integer(h))
    layers[[length(layers) + 1L]] <- act_mod()
    in_dim <- as.integer(h)
  }
  # Drop the activation that follows the embedding (last hidden) layer, so the
  # embedding is the raw penultimate Linear output (matches .lrt_deepmaha_forward).
  layers[[length(layers)]] <- NULL
  d <- in_dim
  head <- nn$Linear(d, 2L)
  layers[[length(layers) + 1L]] <- head
  net <- do.call(nn$Sequential, layers)$to(device)

  opt <- torch$optim$Adam(net$parameters(), lr = hp$lr,
                          weight_decay = hp$weight_decay)
  lossf <- nn$CrossEntropyLoss()
  n <- nrow(Z)
  bs <- min(hp$batch_size, n)

  net$train()
  for (ep in seq_len(hp$epochs)) {
    # Deterministic full-pass minibatching: a fixed contiguous partition (no
    # shuffling) so training is bit-reproducible on cpu for a given seed.
    starts <- seq.int(0L, n - 1L, by = bs)
    for (s0 in starts) {
      idx <- torch$arange(as.integer(s0),
                          as.integer(min(s0 + bs, n)),
                          dtype = torch$long)$to(device)
      xb <- X_py$index_select(0L, idx)
      yb <- y_py$index_select(0L, idx)
      opt$zero_grad()
      out <- net(xb)
      loss <- lossf(out, yb)
      loss$backward()
      opt$step()
    }
  }
  net$eval()

  # ---- export weights up to and INCLUDING the embedding layer ---------------
  # net is a Sequential of Linear/activation modules ending in the 2-logit head.
  # The embedding layers are the Linear modules EXCEPT the last (the head). Export
  # each as W (out x in) and b (length out) in float64.
  modules <- reticulate::iterate(net$modules())
  weights <- list()
  for (m in modules) {
    is_linear <- tryCatch(inherits(m, "torch.nn.modules.linear.Linear") ||
                            !is.null(m$weight), error = function(e) FALSE)
    if (!isTRUE(is_linear)) next
    has_w <- tryCatch(!is.null(m$weight), error = function(e) FALSE)
    if (!isTRUE(has_w)) next
    W <- reticulate::py_to_r(m$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(m$bias$detach()$cpu()$to(torch$float64)$numpy())
    weights[[length(weights) + 1L]] <- list(W = W, b = as.numeric(b))
  }
  # Drop the final (head) Linear: keep only up to the embedding layer.
  if (length(weights) < 2L) {
    stop("fit_lrt_deepmaha: exported net has too few linear layers", call. = FALSE)
  }
  weights[[length(weights)]] <- NULL
  weights
}


# ----------------------------------------------------------------------------
# fit_lrt_deepmaha
# ----------------------------------------------------------------------------

#' @title Fit the deep-feature class-conditional Mahalanobis LRT discriminator
#'
#' @description
#' Trains a BatchNorm-free MLP embedding on the per-sample-rCLR, universe-aligned
#' TRAINING matrix (via reticulate-python torch; Adam, cross-entropy), EXPORTS the
#' weights up to the penultimate embedding layer to R as plain numeric
#' matrices/vectors, DISCARDS the python module and the 2-logit head, then
#' estimates the frozen per-class Gaussian means and a single TIED ridge-regularised
#' covariance in the embedding space using the PURE-R forward the scorer evaluates
#' (guaranteeing fit/score consistency despite float32 training). The fitted model
#' holds no external pointer.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{hidden}
#'   (positive-integer vector of layer widths up to the embedding; the LAST entry
#'   is the embedding dim \eqn{d}; default \code{c(64L, 32L)}), \code{activation}
#'   (\code{"relu"} (default) or \code{"tanh"}), \code{epochs} (positive integer;
#'   default \code{200L}), \code{lr} (positive Adam learning rate; default
#'   \code{1e-3}), \code{weight_decay} (non-negative Adam L2; default \code{1e-4}),
#'   \code{batch_size} (positive integer minibatch; default \code{32L}),
#'   \code{ridge} (positive covariance ridge factor; default \code{1e-3}),
#'   \code{min_features} (feature-overlap floor at scoring, positive integer;
#'   default \code{3L}), \code{device} (\code{"cpu"} (default), \code{"cuda"}, or
#'   \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU with no GPU), and
#'   \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{lrt_deepmaha_model}: a list with
#'   \code{feature_universe}, \code{weights} (exported per-layer \code{list(W, b)}
#'   up to the embedding), \code{activation}, \code{mu_case}, \code{mu_ctrl},
#'   \code{Sigma_reg_inv} (inverse of the ridge-regularised tied covariance),
#'   \code{embedding_dim}, \code{device} (resolved), \code{seed}, and \code{hp}.
#'
#' @details
#' The score is the class-conditional Gaussian log-likelihood ratio with a shared
#' covariance: \eqn{s(z) = -\tfrac12 (z-\mu_{\mathrm{case}})^\top \Sigma^{-1}
#' (z-\mu_{\mathrm{case}}) + \tfrac12 (z-\mu_{\mathrm{ctrl}})^\top \Sigma^{-1}
#' (z-\mu_{\mathrm{ctrl}})}, exactly linear in the embedding \eqn{z} (LDA in
#' embedding space). The tied covariance is the pooled within-class estimator
#' \eqn{\Sigma = ((n_1-1)S_1 + (n_0-1)S_0)/(n_1+n_0-2)}, ridge-regularised as
#' \eqn{\Sigma + \lambda\,\overline{\mathrm{diag}(\Sigma)}\,I} with
#' \eqn{\lambda = } \code{hp$ridge}, with a positive-definite floor if needed. This
#' uses the unbiased pooled denominator; Lee 2018 displays the \eqn{1/N} MLE form,
#' which differs only by a constant scaling of the tied covariance and so leaves the
#' per-specimen score ranking (and AUC) unchanged.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 80; p <- 30; k <- 8
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' model <- fit_lrt_deepmaha(X, y)
#' score_lrt_deepmaha(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Lee K, Lee K, Lee H, Shin J. (2018) A Simple Unified Framework for Detecting
#' Out-of-Distribution Samples and Adversarial Attacks. \emph{NeurIPS} 31.
#' arXiv:1807.03888.
#'
#' @export
fit_lrt_deepmaha <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_lrt_deepmaha", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_lrt_deepmaha", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_lrt_deepmaha")
  hp <- .lrt_deepmaha_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_lrt_deepmaha: X_train must contain at least hp$min_features ",
         "features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's
  # first python initialisation seeds R's RNG as a side effect, so the no-prior-
  # seed state must be captured here (not after the import) to be restored exactly
  # (including removing a leaked .Random.seed when none pre-existed). Fit draws no
  # R randomness of its own; this purely neutralises that init side effect.
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)

  .lrt_deepmaha_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .lrt_deepmaha_rclr_matrix(X_train)        # n x p rCLR over full universe

  device <- .lrt_deepmaha_resolve_device(hp$device)
  weights <- .lrt_deepmaha_train_export(Z_in, y, hp, device)
  activation <- hp$activation

  # Frozen Gaussian stats computed from the PURE-R forward the scorer uses (NOT the
  # torch embeddings), so fit and score share one embedding (float32->float64 safe).
  Z_emb <- .lrt_deepmaha_forward(Z_in, weights, activation)   # n x d
  d <- ncol(Z_emb)
  case <- which(y == 1L)
  ctrl <- which(y == 0L)
  mu_case <- colMeans(Z_emb[case, , drop = FALSE])
  mu_ctrl <- colMeans(Z_emb[ctrl, , drop = FALSE])

  # TIED pooled within-class covariance (Lee 2018 single shared covariance ->
  # linear/LDA score). cov() needs >= 2 rows per class; guard tiny classes by a
  # zero contribution.
  n1 <- length(case); n0 <- length(ctrl)
  S1 <- if (n1 >= 2L) stats::cov(Z_emb[case, , drop = FALSE]) else matrix(0, d, d)
  S0 <- if (n0 >= 2L) stats::cov(Z_emb[ctrl, , drop = FALSE]) else matrix(0, d, d)
  denom <- max(n1 + n0 - 2L, 1L)
  Sigma <- ((n1 - 1L) * S1 + (n0 - 1L) * S0) / denom

  # Ridge regularisation + PD floor on the diagonal.
  ridge <- hp$ridge * mean(diag(Sigma))
  if (!is.finite(ridge) || ridge <= 0) ridge <- hp$ridge
  Sigma_reg <- Sigma + ridge * diag(d)
  # Ensure positive-definite (symmetric eigen floor) before inversion.
  Sigma_reg <- (Sigma_reg + t(Sigma_reg)) / 2
  ev <- eigen(Sigma_reg, symmetric = TRUE)
  floor_val <- max(ev$values) * 1e-10
  if (floor_val <= 0 || !is.finite(floor_val)) floor_val <- 1e-10
  vals <- pmax(ev$values, floor_val)
  Sigma_reg <- ev$vectors %*% (vals * t(ev$vectors))
  Sigma_reg <- (Sigma_reg + t(Sigma_reg)) / 2
  Sigma_reg_inv <- chol2inv(chol(Sigma_reg))
  Sigma_reg_inv <- (Sigma_reg_inv + t(Sigma_reg_inv)) / 2

  model <- list(
    feature_universe = feat_order,
    weights = weights,
    activation = activation,
    mu_case = as.numeric(mu_case),
    mu_ctrl = as.numeric(mu_ctrl),
    Sigma_reg_inv = Sigma_reg_inv,
    embedding_dim = d,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "lrt_deepmaha_model"
  model
}


# Class-conditional Mahalanobis LRT score for one embedded specimen z (length d):
#   -0.5 (z-mu_case)' Sinv (z-mu_case) + 0.5 (z-mu_ctrl)' Sinv (z-mu_ctrl).
# Larger = more case-like.
.lrt_deepmaha_score_one <- function(z, mu_case, mu_ctrl, Sinv) {
  dc <- z - mu_case
  dk <- z - mu_ctrl
  -0.5 * as.numeric(crossprod(dc, Sinv %*% dc)) +
    0.5 * as.numeric(crossprod(dk, Sinv %*% dk))
}


# ----------------------------------------------------------------------------
# score_lrt_deepmaha
# ----------------------------------------------------------------------------

#' @title Score the deep-feature class-conditional Mahalanobis LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_lrt_deepmaha}}, in PURE R (no python). Each query is mapped to
#' the per-sample robust CLR over the frozen feature universe (absent universe
#' features carry the neutral rCLR value \code{0}), embedded by the exported MLP
#' forward, and scored by the class-conditional Mahalanobis log-likelihood ratio
#' (case minus control) with the frozen tied covariance. Larger = more case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on
#' that row and the frozen model and is exactly invariant to per-specimen positive
#' scaling.
#'
#' @param model A \code{lrt_deepmaha_model} object from
#'   \code{\link{fit_lrt_deepmaha}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values are
#'   more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_lrt_deepmaha(X, y)
#' score_lrt_deepmaha(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Lee K, Lee K, Lee H, Shin J. (2018) A Simple Unified Framework for Detecting
#' Out-of-Distribution Samples and Adversarial Attacks. \emph{NeurIPS} 31.
#' arXiv:1807.03888.
#'
#' @export
score_lrt_deepmaha <- function(model, X, meta = NULL) {
  if (!inherits(model, "lrt_deepmaha_model")) {
    stop("score_lrt_deepmaha: model must have class lrt_deepmaha_model")
  }
  X <- .reo_check_matrix(X, "score_lrt_deepmaha", "X")
  .reo_check_meta(meta, nrow(X), "score_lrt_deepmaha", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .lrt_deepmaha_align(X, feat_order)
  weights <- model$weights
  activation <- model$activation
  mu_case <- model$mu_case
  mu_ctrl <- model$mu_ctrl
  Sinv <- model$Sigma_reg_inv

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (rCLR
    # origin) but is a valid specimen that MUST be embedded and scored, not floored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .lrt_deepmaha_rclr(X_use[i, ])
    emb <- .lrt_deepmaha_forward(matrix(z, nrow = 1L), weights, activation)
    out[i] <- .lrt_deepmaha_score_one(as.numeric(emb), mu_case, mu_ctrl, Sinv)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_lrt_deepmaha: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen, score-determining R-side state. The model holds NO live
# external-pointer object (the python module is discarded at fit), so every field
# is a plain R object and this is a deterministic snapshot suitable as the
# `model_digest` for singlesample_assert_row_equivariant() (clause (d): the bytes
# must be identical before and after scoring). Two seed-matched CPU fits on
# identical data produce identical digests.
.lrt_deepmaha_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    weights = model$weights,
    activation = model$activation,
    mu_case = model$mu_case,
    mu_ctrl = model$mu_ctrl,
    Sigma_reg_inv = model$Sigma_reg_inv,
    embedding_dim = model$embedding_dim,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
