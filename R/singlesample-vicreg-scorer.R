#' @title VICReg self-supervised embedding + frozen linear-probe discriminator
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing a
#' self-supervised VICReg representation in the sense of Bardes, Ponce & LeCun
#' (2022), adapted to compositional specimens. A small BatchNorm-free MLP encoder
#' of the per-sample-rCLR specimen is learned on the TRAINING set by SELF-SUPERVISED
#' \strong{VICReg} regularization (variance-invariance-covariance loss on two
#' augmented views, \strong{no labels}) via reticulate-python torch, then
#' \strong{FROZEN}; a frozen linear-probe head (the SSL linear-probe evaluation
#' protocol) is fit on the labeled training embeddings, and a new specimen is scored
#' by its frozen-head logit on its pure-R encoder embedding. Larger score = more
#' case-like.
#'
#' \strong{VICReg self-supervised training (Bardes 2022).} Unlike the
#' cross-entropy-trained \code{lrt-deepmaha} sibling or the prototypical-loss
#' \code{proto-net} sibling, the encoder here is trained with NO labels by the VICReg
#' loss. Each step, for a batch of (rCLR'd, universe-aligned) training rows, TWO
#' augmented views are formed (random feature masking + additive Gaussian jitter in
#' rCLR space, independent draws per view), embedded through the encoder and a small
#' expander/projector to \eqn{Z, Z'}, and the VICReg loss (Bardes 2022 eq. 6)
#' \deqn{\mathcal{L} = \lambda\,s(Z,Z') + \mu\,[v(Z)+v(Z')] + \nu\,[c(Z)+c(Z')]}
#' is minimised, where the invariance term \eqn{s} is the MSE between \eqn{Z} and
#' \eqn{Z'}; the variance term \eqn{v(Z) = \frac1d \sum_j \max(0, \gamma -
#' \sqrt{\mathrm{Var}(Z_j) + \epsilon})} is the hinge (at \eqn{\gamma = 1}) on the
#' per-dimension embedding std; and the covariance term \eqn{c(Z) = \frac1d
#' \sum_{i\neq j} [\mathrm{Cov}(Z)_{ij}]^2} penalises the squared off-diagonal of the
#' embedding covariance. Defaults \eqn{\lambda = 25, \mu = 25, \nu = 1} (the paper's
#' standard). Adam optimises encoder + expander over \code{epochs} steps. The
#' DOWNSTREAM embedding is the ENCODER output (pre-expander) — standard VICReg
#' practice; the expander is discarded and only the encoder is exported.
#'
#' \strong{Torch is confined to FIT.} The encoder is trained with torch on the
#' rCLR'd, universe-aligned training matrix. After SSL training, the encoder weight
#' matrices and biases up to and INCLUDING the embedding layer are EXPORTED to R as
#' plain numeric matrices/vectors; the expander, the python module, and the head's
#' torch object are all DISCARDED (only the head's exported \eqn{w, b} survive). The
#' fitted model therefore holds NO external pointer and NO python dependency at score
#' time: the encoder forward and the linear-head logit are pure base-R. This makes
#' the default \code{model_digest} a stable deterministic snapshot and lets the
#' package suite run the scorer on any host without the venv.
#'
#' \strong{Fit/score consistency despite float32.} Torch trains in float32; an R
#' float64 forward on the exported (float32-valued) weights differs from the torch
#' embedding by \eqn{\sim 10^{-6}}. To guarantee the frozen linear head is fit in
#' EXACTLY the same embedding the scorer evaluates, the head trains on the embeddings
#' produced by the PURE-R forward over ALL training rows
#' (\code{.ssl_vicreg_forward}), NOT the torch embeddings. Fit and score thus use the
#' identical computation (the deepmaha / proto-net fit/score-consistency lesson: fit
#' and score must share one code path).
#'
#' \strong{The embedding forward.} The exported encoder is
#' \eqn{\mathrm{input}(p) \to [\mathrm{Linear} \to \phi]_{1}
#'   \to \dots \to [\mathrm{Linear} \to \phi]_{L-1}
#'   \to \mathrm{Linear} \to \mathrm{embedding}(d)}. There are \code{length(hidden)}
#' linear layers up to the embedding (the last \code{hidden} width is the embedding
#' dim \eqn{d}); the activation (\code{"relu"} default, or \code{"tanh"}) is applied
#' AFTER every linear layer EXCEPT the final embedding layer (no activation follows
#' the embedding -- it is the raw representation, matching the torch encoder).
#' \code{.ssl_vicreg_forward(M, weights)} reimplements this exactly.
#'
#' \strong{Frozen linear-probe logit score.} With the frozen head weight \eqn{w}
#' (length \eqn{d}) and bias \eqn{b}, the score of an embedded specimen \eqn{z} is the
#' linear logit
#' \deqn{s(z) = w^\top z + b,}
#' the SSL linear-probe evaluation read-out. The head is a single BN-free
#' \eqn{\mathrm{Linear}(d \to 1)} trained by seeded, fixed-epoch,
#' \code{head_weight_decay}-regularised optimisation on the FROZEN encoder embeddings
#' (regularisation handles separation on small \eqn{n}). Larger = more case-like.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly
#' invariant to per-specimen scaling and uses no cross-row statistic. Universe
#' features absent from a specimen carry the neutral rCLR value \code{0}. The SAME
#' rCLR transform is applied to the frozen training rows at fit and to every query at
#' score, so a specimen's score depends only on that specimen and the frozen model --
#' it passes \code{\link{singlesample_assert_row_equivariant}}. The augmentation
#' (masking + jitter) is applied ONLY at fit to form the SSL views; scoring uses the
#' clean rCLR with no augmentation and no RNG.
#'
#' \strong{Degenerate-neutral.} A query returns the neutral score \code{0} if fewer
#' than \code{hp$min_features} universe features overlap the query columns (a
#' column-overlap floor, batch-independent), or it has empty positive support over the
#' frozen universe tested on the (pre-rCLR) aligned ORIGINAL abundances
#' (\code{!any(X_use[i, ] > 0)}). A FLAT all-equal-positive composition maps to the
#' rCLR origin but is a VALID specimen and is scored normally (NOT floored) -- the
#' empty-support test is on the ORIGINAL abundances, NOT on \code{all(z == 0)}.
#'
#' @references
#' Bardes A, Ponce J, LeCun Y. (2022) VICReg: Variance-Invariance-Covariance
#' Regularization for Self-Supervised Learning. \emph{International Conference on
#' Learning Representations (ICLR)}. arXiv:2105.04906.
#'
#' @name singlesample-vicreg
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .proto_net_rclr /
# .lrt_deepmaha_rclr.
.ssl_vicreg_rclr <- function(v) {
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
.ssl_vicreg_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .ssl_vicreg_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.ssl_vicreg_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R encoder forward (exact reimplementation of the torch encoder up to d)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported encoder up to (and including) the embedding
# layer. `weights` is a list of per-layer lists list(W = out x in, b = out):
# layers 1..(L-1) are hidden (Linear then `activation`), layer L is the embedding
# Linear with NO activation (matches the torch encoder exactly). For a layer with
# weight W (out x in) and bias b, torch computes M %*% t(W) + b (row = sample).
# Returns an n x d embedding.
.ssl_vicreg_forward <- function(M, weights, activation) {
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

.ssl_vicreg_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ssl_vicreg: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_ssl_vicreg: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_ssl_vicreg: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("hidden", "activation", "epochs", "lr", "weight_decay",
               "mask_p", "noise_sd", "expander_dim", "head_epochs",
               "head_weight_decay", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ssl_vicreg: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  hidden <- hp$hidden
  if (is.null(hidden)) hidden <- c(64L, 32L)
  if (!is.numeric(hidden) || length(hidden) < 1L || any(!is.finite(hidden)) ||
      any(hidden < 1L) || any(hidden > .Machine$integer.max) ||
      any(hidden != as.integer(hidden))) {
    stop("fit_ssl_vicreg: hp$hidden must be a vector of positive integers")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_ssl_vicreg: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_ssl_vicreg: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_ssl_vicreg: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_ssl_vicreg: hp$weight_decay must be a non-negative finite number")
  }

  mask_p <- hp$mask_p
  if (is.null(mask_p)) mask_p <- 0.3
  if (!is.numeric(mask_p) || length(mask_p) != 1L || !is.finite(mask_p) ||
      mask_p <= 0 || mask_p >= 1) {
    stop("fit_ssl_vicreg: hp$mask_p must be a number in (0, 1)")
  }

  noise_sd <- hp$noise_sd
  if (is.null(noise_sd)) noise_sd <- 0.1
  if (!is.numeric(noise_sd) || length(noise_sd) != 1L || !is.finite(noise_sd) ||
      noise_sd < 0) {
    stop("fit_ssl_vicreg: hp$noise_sd must be a non-negative finite number")
  }

  expander_dim <- hp$expander_dim
  if (is.null(expander_dim)) expander_dim <- 64L
  if (!is.numeric(expander_dim) || length(expander_dim) != 1L ||
      !is.finite(expander_dim) || expander_dim < 1L ||
      expander_dim > .Machine$integer.max ||
      expander_dim != as.integer(expander_dim)) {
    stop("fit_ssl_vicreg: hp$expander_dim must be a positive integer")
  }

  head_epochs <- hp$head_epochs
  if (is.null(head_epochs)) head_epochs <- 200L
  if (!is.numeric(head_epochs) || length(head_epochs) != 1L ||
      !is.finite(head_epochs) || head_epochs < 1L ||
      head_epochs > .Machine$integer.max ||
      head_epochs != as.integer(head_epochs)) {
    stop("fit_ssl_vicreg: hp$head_epochs must be a positive integer")
  }

  head_weight_decay <- hp$head_weight_decay
  if (is.null(head_weight_decay)) head_weight_decay <- 1e-2
  if (!is.numeric(head_weight_decay) || length(head_weight_decay) != 1L ||
      !is.finite(head_weight_decay) || head_weight_decay < 0) {
    stop("fit_ssl_vicreg: hp$head_weight_decay must be a non-negative finite ",
         "number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_ssl_vicreg: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_ssl_vicreg: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_ssl_vicreg: hp$seed must be a single integer")
  }

  list(
    hidden = as.integer(hidden),
    activation = activation,
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    mask_p = as.numeric(mask_p),
    noise_sd = as.numeric(noise_sd),
    expander_dim = as.integer(expander_dim),
    head_epochs = as.integer(head_epochs),
    head_weight_decay = as.numeric(head_weight_decay),
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
# download, so only reticulate + `import torch` are probed.
.ssl_vicreg_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_ssl_vicreg: package 'reticulate' is required to train the ",
         "embedding. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_ssl_vicreg: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.ssl_vicreg_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# VICReg self-supervised training + frozen-head training + weight export
# ----------------------------------------------------------------------------

# Train the BN-free MLP encoder by SELF-SUPERVISED VICReg regularization (Bardes
# 2022; NO labels) via torch, then fit a frozen BN-free Linear(d -> 1) head on the
# encoder embeddings of the LABELED training rows. Returns the exported R-side
# encoder weights (a list of per-layer list(W, b), up to and INCLUDING the
# embedding layer) and the torch-trained head's exported (w, b). The python module
# is used here only; nothing python is returned. The torch RNG (CPU + guarded
# CUDA) is saved/restored, so fit leaves the torch generators byte-unchanged; the
# global R .Random.seed is saved/restored by the caller. `Z` is the n x p rCLR'd,
# universe-aligned training matrix; `y` the integer 0/1 labels.
#
# NOTE on the head: the head returned here is the torch float32 head, exported so
# the caller can re-fit/verify it lives in the SAME (pure-R) embedding. The caller
# (fit_ssl_vicreg) re-trains the head deterministically on the PURE-R encoder
# embeddings (float32->float64 consistency); this function exports both the encoder
# (the durable artefact) and the head it trained from the torch embedding (kept for
# the caller's optional use but NOT used for scoring -- see fit_ssl_vicreg).
.ssl_vicreg_train_export <- function(Z, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn

  # ---- save the torch RNG (CPU + guarded CUDA), restore on exit -------------
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

  # Seed the torch generators (weight init + augmentation draws) so two
  # seed-matched CPU fits initialise AND augment identically.
  torch$manual_seed(as.integer(hp$seed))

  np <- reticulate::import("numpy", delay_load = FALSE)
  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)

  p <- ncol(Z)
  hidden <- hp$hidden
  # Encoder: input(p) -> [Linear -> act]xlen(hidden) but DROP the activation after
  # the final (embedding) layer, so the net IS the encoder up to the embedding.
  act_mod <- if (identical(hp$activation, "tanh")) nn$Tanh else nn$ReLU
  enc_layers <- list()
  in_dim <- p
  for (h in hidden) {
    enc_layers[[length(enc_layers) + 1L]] <- nn$Linear(in_dim, as.integer(h))
    enc_layers[[length(enc_layers) + 1L]] <- act_mod()
    in_dim <- as.integer(h)
  }
  enc_layers[[length(enc_layers)]] <- NULL    # drop activation after the embedding
  d <- in_dim
  encoder <- do.call(nn$Sequential, enc_layers)$to(device)

  # Expander/projector: embedding(d) -> [Linear -> act] -> Linear(expander_dim).
  # Standard VICReg expander (2-layer here for a small tabular setting). Discarded
  # after training; only the encoder is exported.
  e_dim <- hp$expander_dim
  expander <- nn$Sequential(
    nn$Linear(d, e_dim), act_mod(), nn$Linear(e_dim, e_dim)
  )$to(device)

  params <- c(reticulate::iterate(encoder$parameters()),
              reticulate::iterate(expander$parameters()))
  opt <- torch$optim$Adam(params, lr = hp$lr, weight_decay = hp$weight_decay)

  n <- nrow(Z)
  lambda <- 25.0       # invariance coefficient (Bardes 2022 standard)
  mu     <- 25.0       # variance   coefficient
  nu     <- 1.0        # covariance coefficient
  gamma  <- 1.0        # variance hinge target (per-dim std)
  eps    <- 1e-4       # variance numerical floor (Bardes 2022)
  mask_p <- hp$mask_p
  noise_sd <- hp$noise_sd

  # One augmented view of the full batch X_py: independent per-element feature
  # mask (Bernoulli keep-prob 1-mask_p, zeroing dropped features in rCLR space =
  # the per-sample geometric-mean reference) + additive Gaussian jitter.
  make_view <- function(Xt) {
    keep <- torch$bernoulli(
      torch$full(Xt$shape, 1.0 - mask_p, device = device))
    noise <- torch$randn(Xt$shape, device = device)$mul(noise_sd)
    Xt$mul(keep)$add(noise)
  }

  # VICReg variance + covariance terms on an n x e embedding tensor Zt. Both use
  # the UNBIASED (n-1) empirical (co)variance over the batch, matching torch's
  # `.var(dim=0)` default and the official VICReg implementation. `e` is the
  # projection (expander) width, known R-side, so no Python tuple indexing.
  denom <- torch$sub(torch$as_tensor(as.numeric(n)), 1.0)   # (n - 1)
  e_dim_int <- as.integer(e_dim)
  diag_mask <- torch$eye(e_dim_int, device = device)
  var_cov <- function(Zt) {
    Zc <- Zt$sub(Zt$mean(0L, keepdim = TRUE))            # centre per dimension
    # variance term: hinge on per-dim std (unbiased variance + eps under sqrt)
    var_j <- Zc$pow(2)$sum(0L)$div(denom)                # length e, unbiased var
    std <- var_j$add(eps)$sqrt()
    v <- nn$functional$relu(torch$sub(gamma, std))$mean()
    # covariance term: squared off-diagonal of the (n-1) sample covariance / e
    cov <- torch$matmul(Zc$t(), Zc)$div(denom)           # e x e covariance
    off <- cov$mul(torch$sub(1.0, diag_mask))            # zero the diagonal
    c_term <- off$pow(2)$sum()$div(torch$as_tensor(as.numeric(e_dim_int)))
    list(v = v, c = c_term)
  }

  encoder$train(); expander$train()
  for (ep in seq_len(hp$epochs)) {
    v1 <- make_view(X_py)
    v2 <- make_view(X_py)
    z1 <- expander(encoder(v1))                          # n x e
    z2 <- expander(encoder(v2))                          # n x e

    # invariance: MSE between the two views' projections
    s_term <- nn$functional$mse_loss(z1, z2)
    vc1 <- var_cov(z1)
    vc2 <- var_cov(z2)
    v_term <- torch$add(vc1$v, vc2$v)
    c_term <- torch$add(vc1$c, vc2$c)

    loss <- torch$add(
      torch$add(s_term$mul(lambda), v_term$mul(mu)),
      c_term$mul(nu))

    opt$zero_grad()
    loss$backward()
    opt$step()
  }
  encoder$eval(); expander$eval()

  # ---- export encoder weights up to and INCLUDING the embedding layer -------
  modules <- reticulate::iterate(encoder$modules())
  weights <- list()
  for (m in modules) {
    has_w <- tryCatch(!is.null(m$weight), error = function(e) FALSE)
    if (!isTRUE(has_w)) next
    W <- reticulate::py_to_r(m$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(m$bias$detach()$cpu()$to(torch$float64)$numpy())
    weights[[length(weights) + 1L]] <- list(W = W, b = as.numeric(b))
  }
  if (length(weights) < 1L) {
    stop("fit_ssl_vicreg: exported encoder has no linear layers", call. = FALSE)
  }
  weights
}

# Train a frozen BN-free Linear(d -> 1) head on the FROZEN (pure-R) encoder
# embeddings `E` (n x d, float64) of the labeled training rows, by seeded,
# fixed-epoch, weight_decay-regularised full-batch Adam, then export (w, b) to R.
# This is the SSL linear-probe evaluation protocol; the embeddings are detached
# constants (the encoder is frozen). Training on the PURE-R embeddings the scorer
# uses guarantees the head lives in EXACTLY that embedding (float32->float64 safe).
# The torch RNG is saved/restored. `y` is the integer 0/1 label vector.
.ssl_vicreg_train_head <- function(E, y, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn

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

  # Seed the head init independently of the SSL stage (seed + 1) so the head's
  # weight init is deterministic AND decoupled from the SSL generator state. Both
  # stages are seeded, so two seed=42 fits replay both stages identically.
  torch$manual_seed(as.integer(hp$seed) + 1L)

  np <- reticulate::import("numpy", delay_load = FALSE)
  d <- ncol(E)
  E_py <- torch$as_tensor(np$asarray(E, dtype = "float32"))$float()$to(device)
  # BCEWithLogits target in {0,1}, shape n x 1.
  y_py <- torch$as_tensor(np$asarray(as.numeric(y), dtype = "float32")
                          )$float()$reshape(c(-1L, 1L))$to(device)

  head <- nn$Linear(as.integer(d), 1L)$to(device)
  opt <- torch$optim$Adam(head$parameters(), lr = hp$lr,
                          weight_decay = hp$head_weight_decay)
  lossf <- nn$BCEWithLogitsLoss()

  head$train()
  for (ep in seq_len(hp$head_epochs)) {
    opt$zero_grad()
    logit <- head(E_py)
    loss <- lossf(logit, y_py)
    loss$backward()
    opt$step()
  }
  head$eval()

  w <- reticulate::py_to_r(
    head$weight$detach()$cpu()$to(torch$float64)$numpy())  # 1 x d
  b <- reticulate::py_to_r(
    head$bias$detach()$cpu()$to(torch$float64)$numpy())    # length 1
  list(w = as.numeric(w), b = as.numeric(b))
}


# ----------------------------------------------------------------------------
# fit_ssl_vicreg
# ----------------------------------------------------------------------------

#' @title Fit the VICReg self-supervised single-sample discriminator
#'
#' @description
#' Trains a BatchNorm-free MLP encoder on the per-sample-rCLR, universe-aligned
#' TRAINING matrix by SELF-SUPERVISED VICReg regularization (Bardes 2022; via
#' reticulate-python torch, Adam, variance-invariance-covariance loss on two
#' augmented views -- NO labels), EXPORTS the encoder weights up to and including the
#' embedding layer to R as plain numeric matrices/vectors, DISCARDS the expander and
#' the python module, then fits a FROZEN linear-probe head (\eqn{w, b}) on the
#' labeled training embeddings computed by the PURE-R forward the scorer uses
#' (guaranteeing fit/score consistency despite float32 training). The fitted model
#' holds no external pointer.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control. Labels are used
#'   ONLY for the frozen linear-probe head; the encoder is trained without labels.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{hidden}
#'   (positive-integer vector of encoder layer widths up to the embedding; the LAST
#'   entry is the embedding dim \eqn{d}; default \code{c(64L, 32L)}),
#'   \code{activation} (\code{"relu"} (default) or \code{"tanh"}), \code{epochs}
#'   (number of SSL training steps, positive integer; default \code{200L}),
#'   \code{lr} (positive Adam learning rate; default \code{1e-3}),
#'   \code{weight_decay} (non-negative encoder/expander Adam L2; default
#'   \code{1e-4}), \code{mask_p} (per-view feature-mask fraction in (0, 1); default
#'   \code{0.3}), \code{noise_sd} (non-negative rCLR Gaussian-jitter sd; default
#'   \code{0.1}), \code{expander_dim} (positive-integer expander/projector width;
#'   default \code{64L}), \code{head_epochs} (linear-probe fixed epochs, positive
#'   integer; default \code{200L}), \code{head_weight_decay} (non-negative head Adam
#'   L2; default \code{1e-2}), \code{min_features} (feature-overlap floor at scoring,
#'   positive integer; default \code{3L}), \code{device} (\code{"cpu"} (default),
#'   \code{"cuda"}, or \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU
#'   with no GPU), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{ssl_vicreg_model}: a list with
#'   \code{feature_universe}, \code{weights} (exported per-layer \code{list(W, b)}
#'   up to the embedding), \code{activation}, \code{head_w} (length-\eqn{d} frozen
#'   head weight), \code{head_b} (scalar frozen head bias), \code{embedding_dim},
#'   \code{device} (resolved), \code{seed}, and \code{hp}.
#'
#' @details
#' The score is the frozen linear-probe logit \eqn{s(z) = w^\top z + b} on the
#' encoder embedding \eqn{z}. The distinction from the cross-entropy
#' \code{lrt-deepmaha} and prototypical \code{proto-net} siblings is the TRAINING:
#' the encoder is learned with NO labels by the VICReg loss (Bardes 2022 eq. 6,
#' \eqn{\lambda = 25, \mu = 25, \nu = 1}) on two augmented views, and only a small
#' linear head reads the labels.
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
#' model <- fit_ssl_vicreg(X, y)
#' score_ssl_vicreg(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Bardes A, Ponce J, LeCun Y. (2022) VICReg: Variance-Invariance-Covariance
#' Regularization for Self-Supervised Learning. \emph{ICLR}. arXiv:2105.04906.
#'
#' @export
fit_ssl_vicreg <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ssl_vicreg", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ssl_vicreg", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ssl_vicreg")
  hp <- .ssl_vicreg_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_ssl_vicreg: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's
  # first python init seeds R's RNG as a side effect. Fit draws no R randomness of
  # its own (augmentation + head init are torch-seeded), so this purely neutralises
  # the init side effect. Restoring on exit (incl. the no-prior-seed removal case)
  # keeps fit globally RNG-neutral.
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

  .ssl_vicreg_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .ssl_vicreg_rclr_matrix(X_train)          # n x p rCLR over full universe

  device <- .ssl_vicreg_resolve_device(hp$device)
  weights <- .ssl_vicreg_train_export(Z_in, hp, device)
  activation <- hp$activation

  # Frozen encoder embeddings from the PURE-R forward over ALL training rows (NOT
  # the torch embeddings) -- fit and score share one embedding code path
  # (float32 -> float64 safe). The linear-probe head trains on THESE embeddings.
  Z_emb <- .ssl_vicreg_forward(Z_in, weights, activation)   # n x d
  d <- ncol(Z_emb)
  head <- .ssl_vicreg_train_head(Z_emb, y, hp, device)

  model <- list(
    feature_universe = feat_order,
    weights = weights,
    activation = activation,
    head_w = as.numeric(head$w),
    head_b = as.numeric(head$b),
    embedding_dim = d,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "ssl_vicreg_model"
  model
}


# Frozen linear-probe logit for one embedded specimen z (length d): w . z + b.
# Larger = more case-like.
.ssl_vicreg_score_one <- function(z, head_w, head_b) {
  sum(head_w * z) + head_b
}


# ----------------------------------------------------------------------------
# score_ssl_vicreg
# ----------------------------------------------------------------------------

#' @title Score the VICReg self-supervised single-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_ssl_vicreg}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), embedded by the exported encoder forward,
#' and scored by the frozen linear-probe logit \eqn{w^\top z + b}. Larger = more
#' case-like. No augmentation and no RNG are used at scoring.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on that
#' row and the frozen model and is exactly invariant to per-specimen positive scaling.
#'
#' @param model A \code{ssl_vicreg_model} object from \code{\link{fit_ssl_vicreg}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity and
#'   ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values are
#'   more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_ssl_vicreg(X, y)
#' score_ssl_vicreg(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Bardes A, Ponce J, LeCun Y. (2022) VICReg: Variance-Invariance-Covariance
#' Regularization for Self-Supervised Learning. \emph{ICLR}. arXiv:2105.04906.
#'
#' @export
score_ssl_vicreg <- function(model, X, meta = NULL) {
  if (!inherits(model, "ssl_vicreg_model")) {
    stop("score_ssl_vicreg: model must have class ssl_vicreg_model")
  }
  X <- .reo_check_matrix(X, "score_ssl_vicreg", "X")
  .reo_check_meta(meta, nrow(X), "score_ssl_vicreg", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .ssl_vicreg_align(X, feat_order)
  weights <- model$weights
  activation <- model$activation
  head_w <- model$head_w
  head_b <- model$head_b

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (rCLR
    # origin) but is a valid specimen that MUST be embedded and scored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .ssl_vicreg_rclr(X_use[i, ])
    emb <- .ssl_vicreg_forward(matrix(z, nrow = 1L), weights, activation)
    out[i] <- .ssl_vicreg_score_one(as.numeric(emb), head_w, head_b)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ssl_vicreg: scorer produced non-finite or wrong-length output")
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
.ssl_vicreg_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    weights = model$weights,
    activation = model$activation,
    head_w = model$head_w,
    head_b = model$head_b,
    embedding_dim = model$embedding_dim,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
