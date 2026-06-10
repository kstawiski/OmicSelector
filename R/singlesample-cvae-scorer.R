#' @title Counterfactual class-conditional VAE single-sample discriminator
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing a
#' counterfactual CLASS-conditional variational autoencoder used as a generative
#' classifier. A class-conditional VAE (encoder \eqn{q(z|x,y)} + decoder
#' \eqn{p(x|z,y)}, both plain BatchNorm-free MLPs conditioned on a one-hot of the
#' class label \eqn{y}) is trained on the TRAINING set by the reparameterized ELBO
#' via reticulate-python torch, then \strong{FROZEN}; a new specimen is scored by the
#' COUNTERFACTUAL evidence log-likelihood-RATIO
#' \deqn{s(x) = \mathrm{ELBO}(x, y=\mathrm{case}) - \mathrm{ELBO}(x, y=\mathrm{control}),}
#' computed in PURE base-R with the latent MEAN \eqn{z = \mu} (DETERMINISTIC, no
#' sampling, no RNG at score). Larger score = more case-like.
#'
#' \strong{PRESPEC DEPARTURE (kit-conditional -> class-conditional).} The present
#' \code{fit_kit_conditional_vae} (R/singlesample-02-learned-kit-aware-vae.R)
#' conditions a VAE on KIT metadata -- an external cross-sample covariate that
#' VIOLATES single-sample deployability (the kit of a held-out specimen need not be
#' known/seen at score). The single-sample \code{cvae} instead conditions on the
#' CLASS label \eqn{y} and is used as a COUNTERFACTUAL generative classifier: it
#' models \eqn{p(x|y)} for each class and scores a query by the counterfactual
#' evidence ratio \eqn{\mathrm{ELBO}(x,\mathrm{case}) -
#' \mathrm{ELBO}(x,\mathrm{control})}. This needs NO kit/cohort input at score (only
#' the query \eqn{x}; the class conditioning is internal -- both classes are
#' evaluated and differenced), so it is single-sample-faithful. This reinterpretation
#' is flagged for the gate.
#'
#' \strong{Torch is confined to FIT.} The class-conditional VAE is trained with torch
#' on the rCLR'd, universe-aligned training matrix. After training, the encoder,
#' mu-head, logvar-head, decoder, and xhat-head weight matrices and biases, plus the
#' scalar decoder variance \eqn{\sigma^2}, are EXPORTED to R as plain numeric
#' matrices/vectors; the python module is DISCARDED. The fitted model therefore holds
#' NO external pointer and NO python dependency at score time: the encoder forward,
#' the decoder forward, and the per-class ELBO are all pure base-R. This makes the
#' default \code{model_digest} a stable deterministic snapshot and lets the package
#' suite run the scorer on any host without the venv.
#'
#' \strong{Fit/score consistency despite float32.} Torch trains in float32; an R
#' float64 forward on the exported (float32-valued) weights differs from the torch
#' forward by \eqn{\sim 10^{-6}}. The score uses ONLY the exported weights through the
#' SAME pure-R forwards (\code{.cvae_encode} / \code{.cvae_decode}) at both fit-time
#' verification and score, so there is no second code path to drift (the
#' deepmaha / proto-net / ssl-vicreg fit/score-consistency lesson).
#'
#' \strong{The forwards.} The encoder is
#' \eqn{\mathrm{concat}[x, e_y](p+2) \to [\mathrm{Linear} \to \phi]_{1..H_e}
#'   \to (\mathrm{Linear}\to\mu,\ \mathrm{Linear}\to\log\sigma^2_z)}; the activation
#' (\code{"relu"} default, or \code{"tanh"}) follows every hidden Linear, and the two
#' output heads (mu, logvar) are bare Linear (no activation). The decoder is
#' \eqn{\mathrm{concat}[z, e_y](L+2) \to [\mathrm{Linear} \to \phi]_{1..H_d}
#'   \to \mathrm{Linear} \to \hat{x}(p)}; activation follows every hidden Linear, the
#' xhat-head is bare Linear. \eqn{e_y = \mathrm{onehot}(y)} (length 2).
#'
#' \strong{The counterfactual score (pure R, deterministic).} For a query rCLR vector
#' \eqn{x}, for EACH class \eqn{c \in \{0,1\}}:
#' \deqn{(\mu_c, \log\sigma^2_{z,c}) = \mathrm{encode}(x, e_c);\quad z_c = \mu_c\ \mathrm{(latent\ mean)};\quad \hat{x}_c = \mathrm{decode}(\mu_c, e_c);}
#' \deqn{\log N_c = -\tfrac12\big(p\,\log(2\pi\sigma^2) + \tfrac{1}{\sigma^2}\sum_j (x_j - \hat{x}_{c,j})^2\big);}
#' \deqn{\mathrm{KL}_c = \tfrac12 \sum_l \big(\mu_{c,l}^2 + e^{\log\sigma^2_{z,c,l}} - \log\sigma^2_{z,c,l} - 1\big);}
#' \deqn{\mathrm{ELBO}_c = \log N_c - \beta\,\mathrm{KL}_c;\qquad s(x) = \mathrm{ELBO}_{\mathrm{case}} - \mathrm{ELBO}_{\mathrm{ctrl}}.}
#' Using \eqn{z = \mu} (NOT a sample) makes the score DETERMINISTIC and RNG-free, and
#' it depends only on \eqn{x} and the frozen model -> single-sample. Per-row ->
#' single-row == batch EXACTLY 0. \eqn{\beta = \mathrm{kl\_beta}} (default 1.0).
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly invariant
#' to per-specimen scaling and uses no cross-row statistic. Universe features absent
#' from a specimen carry the neutral rCLR value \code{0}. The SAME rCLR transform is
#' applied to the frozen training rows at fit and to every query at score, so a
#' specimen's score depends only on that specimen and the frozen model -- it passes
#' \code{\link{singlesample_assert_row_equivariant}}.
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
#' Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes. \emph{International
#' Conference on Learning Representations (ICLR)}. arXiv:1312.6114.
#'
#' Sohn K, Lee H, Yan X. (2015) Learning Structured Output Representation using Deep
#' Conditional Generative Models. \emph{Advances in Neural Information Processing
#' Systems (NeurIPS)} 28.
#'
#' @name singlesample-cvae
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .proto_net_rclr /
# .ssl_vicreg_rclr.
.cvae_rclr <- function(v) {
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
.cvae_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .cvae_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.cvae_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R MLP forwards (exact reimplementations of the torch encoder + decoder)
# ----------------------------------------------------------------------------

# Generic pure base-R forward of an MLP trunk + a set of bare-Linear output heads.
# `trunk` is a list of per-layer list(W = out x in, b = out), each followed by the
# activation. `heads` is a NAMED list of head list(W, b), each a bare Linear (NO
# activation) applied to the trunk output. For a Linear with weight W (out x in)
# and bias b, torch computes M %*% t(W) + b (row = sample). `M` is n x in_dim.
# Returns a named list of n x out matrices, one per head.
.cvae_mlp_forward <- function(M, trunk, heads, activation) {
  act <- if (identical(activation, "tanh")) tanh else function(x) pmax(x, 0)
  H <- M
  for (l in seq_along(trunk)) {
    W <- trunk[[l]]$W                         # out x in
    b <- trunk[[l]]$b                         # length out
    H <- H %*% t(W)                           # n x out
    H <- sweep(H, 2L, b, `+`)
    H <- act(H)                               # activation after EVERY trunk layer
  }
  out <- vector("list", length(heads))
  names(out) <- names(heads)
  for (nm in names(heads)) {
    W <- heads[[nm]]$W                        # out x in
    b <- heads[[nm]]$b                        # length out
    Hh <- H %*% t(W)
    out[[nm]] <- sweep(Hh, 2L, b, `+`)        # bare Linear (no activation)
  }
  out
}

# Encoder forward q(z|x,y): input concat[rclr(x) (p), onehot(y) (2)] -> mu, logvar
# (each n x L). `M` is the n x (p+2) input matrix (rCLR rows with the one-hot
# appended). Returns list(mu = n x L, logvar = n x L).
.cvae_encode <- function(M, weights, activation) {
  res <- .cvae_mlp_forward(M, weights$enc_trunk,
                           list(mu = weights$enc_mu, logvar = weights$enc_logvar),
                           activation)
  list(mu = res$mu, logvar = res$logvar)
}

# Decoder forward p(x|z,y): input concat[z (L), onehot(y) (2)] -> xhat (n x p).
# `M` is the n x (L+2) input matrix. Returns the n x p reconstruction mean.
.cvae_decode <- function(M, weights, activation) {
  res <- .cvae_mlp_forward(M, weights$dec_trunk,
                           list(xhat = weights$dec_xhat), activation)
  res$xhat
}


# ----------------------------------------------------------------------------
# Per-class ELBO + counterfactual score (pure base-R, deterministic z = mu)
# ----------------------------------------------------------------------------

# ELBO of one rCLR query x (length p) under class c (0/1), using z = mu
# (DETERMINISTIC; no sampling). `weights` carries the exported nets, `sigma2` the
# scalar decoder variance, `kl_beta` the KL coefficient, `activation` the trunk
# activation, `p` the feature count. The one-hot of c is appended to x for the
# encoder and to mu for the decoder.
#   logN_c = -0.5 * ( p*log(2*pi*sigma2) + sum((x - xhat_c)^2)/sigma2 )
#   KL_c   = 0.5 * sum( mu_c^2 + exp(logvar_c) - logvar_c - 1 )
#   ELBO_c = logN_c - kl_beta * KL_c
.cvae_elbo_class <- function(x, c, weights, sigma2, kl_beta, activation, p) {
  oh <- if (identical(as.integer(c), 0L)) c(1, 0) else c(0, 1)   # onehot(c), len 2
  enc_in <- matrix(c(x, oh), nrow = 1L)                          # 1 x (p+2)
  enc <- .cvae_encode(enc_in, weights, activation)
  mu_c <- as.numeric(enc$mu)                                     # length L
  logvar_c <- as.numeric(enc$logvar)                            # length L
  dec_in <- matrix(c(mu_c, oh), nrow = 1L)                       # 1 x (L+2)
  xhat_c <- as.numeric(.cvae_decode(dec_in, weights, activation))  # length p
  logN_c <- -0.5 * (p * log(2 * pi * sigma2) +
                      sum((x - xhat_c)^2) / sigma2)
  kl_c <- 0.5 * sum(mu_c^2 + exp(logvar_c) - logvar_c - 1)
  logN_c - kl_beta * kl_c
}

# Counterfactual evidence log-likelihood ratio for one rCLR query x:
#   s(x) = ELBO(x, case) - ELBO(x, control).  Larger = more case-like.
.cvae_score_one <- function(x, weights, sigma2, kl_beta, activation, p) {
  elbo_case <- .cvae_elbo_class(x, 1L, weights, sigma2, kl_beta, activation, p)
  elbo_ctrl <- .cvae_elbo_class(x, 0L, weights, sigma2, kl_beta, activation, p)
  elbo_case - elbo_ctrl
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.cvae_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_cvae: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_cvae: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_cvae: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("enc_hidden", "dec_hidden", "latent_dim", "activation", "epochs",
               "lr", "weight_decay", "kl_beta", "decoder_sigma2", "min_features",
               "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_cvae: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  enc_hidden <- hp$enc_hidden
  if (is.null(enc_hidden)) enc_hidden <- c(64L, 32L)
  if (!is.numeric(enc_hidden) || length(enc_hidden) < 1L ||
      any(!is.finite(enc_hidden)) || any(enc_hidden < 1L) ||
      any(enc_hidden > .Machine$integer.max) ||
      any(enc_hidden != as.integer(enc_hidden))) {
    stop("fit_cvae: hp$enc_hidden must be a vector of positive integers")
  }

  dec_hidden <- hp$dec_hidden
  if (is.null(dec_hidden)) dec_hidden <- c(32L, 64L)
  if (!is.numeric(dec_hidden) || length(dec_hidden) < 1L ||
      any(!is.finite(dec_hidden)) || any(dec_hidden < 1L) ||
      any(dec_hidden > .Machine$integer.max) ||
      any(dec_hidden != as.integer(dec_hidden))) {
    stop("fit_cvae: hp$dec_hidden must be a vector of positive integers")
  }

  latent_dim <- hp$latent_dim
  if (is.null(latent_dim)) latent_dim <- 16L
  if (!is.numeric(latent_dim) || length(latent_dim) != 1L ||
      !is.finite(latent_dim) || latent_dim < 1L ||
      latent_dim > .Machine$integer.max ||
      latent_dim != as.integer(latent_dim)) {
    stop("fit_cvae: hp$latent_dim must be a positive integer")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_cvae: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_cvae: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_cvae: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_cvae: hp$weight_decay must be a non-negative finite number")
  }

  kl_beta <- hp$kl_beta
  if (is.null(kl_beta)) kl_beta <- 1.0
  if (!is.numeric(kl_beta) || length(kl_beta) != 1L || !is.finite(kl_beta) ||
      kl_beta < 0) {
    stop("fit_cvae: hp$kl_beta must be a non-negative finite number")
  }

  # decoder_sigma2: NULL (default) -> learned (softplus of a free param); a positive
  # finite scalar -> fixed at that value (validate + expose).
  decoder_sigma2 <- hp$decoder_sigma2
  if (!is.null(decoder_sigma2)) {
    if (!is.numeric(decoder_sigma2) || length(decoder_sigma2) != 1L ||
        !is.finite(decoder_sigma2) || decoder_sigma2 <= 0) {
      stop("fit_cvae: hp$decoder_sigma2 must be NULL (learned) or a positive ",
           "finite number")
    }
    decoder_sigma2 <- as.numeric(decoder_sigma2)
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_cvae: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_cvae: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_cvae: hp$seed must be a single integer")
  }

  list(
    enc_hidden = as.integer(enc_hidden),
    dec_hidden = as.integer(dec_hidden),
    latent_dim = as.integer(latent_dim),
    activation = activation,
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    kl_beta = as.numeric(kl_beta),
    decoder_sigma2 = decoder_sigma2,       # NULL = learned
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
.cvae_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_cvae: package 'reticulate' is required to train the VAE. ",
         "Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_cvae: Python module 'torch' is not available in the active ",
         "reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.cvae_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# Class-conditional VAE training (reparameterized ELBO) + weight export
# ----------------------------------------------------------------------------

# Train the BN-free class-conditional VAE (encoder q(z|x,y) + decoder p(x|z,y),
# conditioned on a one-hot of the TRUE label y) by maximizing the reparameterized
# ELBO (minimize -ELBO) via torch full-batch Adam, then export every weight to R as
# plain float64 matrices/vectors plus the scalar decoder variance sigma2. The python
# module is used here only; nothing python is returned. The torch RNG (CPU + guarded
# CUDA) is saved/restored, so fit leaves the torch generators byte-unchanged; the
# global R .Random.seed is saved/restored by the caller. `Z` is the n x p rCLR'd,
# universe-aligned training matrix; `y` the integer 0/1 labels.
.cvae_train_export <- function(Z, y, hp, device) {
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

  # SEED BEFORE BUILD. manual_seed BEFORE constructing ANY stochastic module
  # (encoder trunk + mu-head + logvar-head, decoder trunk + xhat-head, and the
  # sigma2 free parameter) so two seed-matched CPU fits initialise AND draw the
  # reparameterization eps identically -- and so the model is cross-PROCESS
  # reproducible (a build-before-seed regression is masked by the on.exit RNG
  # restore within a single process; only a fresh process or the inter-fit RNG
  # perturbation in §7.8 catches it).
  torch$manual_seed(as.integer(hp$seed))

  np <- reticulate::import("numpy", delay_load = FALSE)
  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)
  # One-hot of the TRUE label y (n x 2), columns = (control, case).
  oh <- cbind(as.numeric(y == 0L), as.numeric(y == 1L))
  Y_py <- torch$as_tensor(np$asarray(oh, dtype = "float32"))$float()$to(device)

  p <- ncol(Z)
  L <- hp$latent_dim
  act_mod <- if (identical(hp$activation, "tanh")) nn$Tanh else nn$ReLU

  # ---- ENCODER trunk: concat[x, e_y] (p + 2) -> enc_hidden, act after each -----
  enc_layers <- list()
  in_dim <- p + 2L
  for (h in hp$enc_hidden) {
    enc_layers[[length(enc_layers) + 1L]] <- nn$Linear(in_dim, as.integer(h))
    enc_layers[[length(enc_layers) + 1L]] <- act_mod()
    in_dim <- as.integer(h)
  }
  enc_trunk <- do.call(nn$Sequential, enc_layers)$to(device)
  enc_mu <- nn$Linear(in_dim, as.integer(L))$to(device)          # mu-head
  enc_logvar <- nn$Linear(in_dim, as.integer(L))$to(device)      # logvar-head

  # ---- DECODER trunk: concat[z, e_y] (L + 2) -> dec_hidden, act after each ------
  dec_layers <- list()
  in_dim_d <- as.integer(L) + 2L
  for (h in hp$dec_hidden) {
    dec_layers[[length(dec_layers) + 1L]] <- nn$Linear(in_dim_d, as.integer(h))
    dec_layers[[length(dec_layers) + 1L]] <- act_mod()
    in_dim_d <- as.integer(h)
  }
  dec_trunk <- do.call(nn$Sequential, dec_layers)$to(device)
  dec_xhat <- nn$Linear(in_dim_d, as.integer(p))$to(device)      # xhat-head

  # ---- decoder variance sigma2: learned (softplus of a free param) or fixed -----
  # A learned raw scalar parameter mapped through softplus to a positive sigma2.
  # When decoder_sigma2 is fixed, sigma2 is a constant tensor (no parameter).
  sigma2_learned <- is.null(hp$decoder_sigma2)
  if (sigma2_learned) {
    # softplus(raw0) = log(1+exp(raw0)) ~ 1.0 -> a sensible rCLR-scale init.
    raw0 <- log(exp(1.0) - 1.0)
    sigma2_raw <- nn$Parameter(
      torch$as_tensor(np$asarray(raw0, dtype = "float32"))$float()$to(device))
    softplus <- nn$functional$softplus
  } else {
    sigma2_fixed_t <- torch$as_tensor(
      np$asarray(hp$decoder_sigma2, dtype = "float32"))$float()$to(device)
  }

  # ---- optimizer over ALL trainable parameters ---------------------------------
  params <- c(reticulate::iterate(enc_trunk$parameters()),
              reticulate::iterate(enc_mu$parameters()),
              reticulate::iterate(enc_logvar$parameters()),
              reticulate::iterate(dec_trunk$parameters()),
              reticulate::iterate(dec_xhat$parameters()))
  if (sigma2_learned) params <- c(params, list(sigma2_raw))
  opt <- torch$optim$Adam(params, lr = hp$lr, weight_decay = hp$weight_decay)

  n <- nrow(Z)
  half_log_2pi <- 0.5 * log(2 * pi)
  beta <- hp$kl_beta

  enc_trunk$train(); enc_mu$train(); enc_logvar$train()
  dec_trunk$train(); dec_xhat$train()
  for (ep in seq_len(hp$epochs)) {
    # Encode q(z|x,y): concat[x, e_y] -> mu, logvar.
    enc_h <- enc_trunk(torch$cat(list(X_py, Y_py), 1L))          # n x last_enc
    mu <- enc_mu(enc_h)                                          # n x L
    logvar <- enc_logvar(enc_h)                                 # n x L
    # Reparameterize z = mu + exp(0.5*logvar) * eps, eps ~ N(0,I) (ONE draw/row/step).
    std <- logvar$mul(0.5)$exp()
    eps <- torch$randn(list(n, as.integer(L)), device = device)
    z <- mu$add(std$mul(eps))                                   # n x L
    # Decode p(x|z,y): concat[z, e_y] -> xhat.
    xhat <- dec_xhat(dec_trunk(torch$cat(list(z, Y_py), 1L)))   # n x p

    sigma2 <- if (sigma2_learned) softplus(sigma2_raw) else sigma2_fixed_t
    # Reconstruction log-lik (diagonal Gaussian, shared scalar sigma2), summed over
    # features, meaned over rows:
    #   logN = -0.5*( p*log(2*pi*sigma2) + sum((x - xhat)^2)/sigma2 )
    sq <- xhat$sub(X_py)$pow(2)$sum(1L)                         # n, sum over features
    log_sigma2 <- sigma2$log()
    logN <- torch$mul(
      torch$add(
        torch$mul(log_sigma2, as.numeric(p) * 0.5),            # 0.5*p*log(sigma2)
        sq$div(sigma2)$mul(0.5)),                              # 0.5*sum sq / sigma2
      -1.0)$sub(as.numeric(p) * half_log_2pi)                  # - 0.5*p*log(2*pi)
    # KL( q(z|x,y) || N(0,I) ) per row, summed over latent dims:
    #   KL = 0.5*sum( mu^2 + exp(logvar) - logvar - 1 )
    kl <- torch$mul(
      torch$sub(torch$sub(torch$add(mu$pow(2), logvar$exp()), logvar), 1.0)$sum(1L),
      0.5)                                                      # n
    elbo <- logN$sub(kl$mul(beta))                             # n
    loss <- elbo$mean()$mul(-1.0)                              # minimize -ELBO

    opt$zero_grad()
    loss$backward()
    opt$step()
  }
  enc_trunk$eval(); enc_mu$eval(); enc_logvar$eval()
  dec_trunk$eval(); dec_xhat$eval()

  # ---- export every weight as plain float64 R objects --------------------------
  to_layers <- function(seq_mod) {
    mods <- reticulate::iterate(seq_mod$modules())
    w <- list()
    for (m in mods) {
      has_w <- tryCatch(!is.null(m$weight), error = function(e) FALSE)
      if (!isTRUE(has_w)) next
      W <- reticulate::py_to_r(m$weight$detach()$cpu()$to(torch$float64)$numpy())
      b <- reticulate::py_to_r(m$bias$detach()$cpu()$to(torch$float64)$numpy())
      w[[length(w) + 1L]] <- list(W = W, b = as.numeric(b))
    }
    w
  }
  one_linear <- function(lin) {
    W <- reticulate::py_to_r(lin$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(lin$bias$detach()$cpu()$to(torch$float64)$numpy())
    list(W = W, b = as.numeric(b))
  }

  weights <- list(
    enc_trunk = to_layers(enc_trunk),
    enc_mu = one_linear(enc_mu),
    enc_logvar = one_linear(enc_logvar),
    dec_trunk = to_layers(dec_trunk),
    dec_xhat = one_linear(dec_xhat)
  )
  if (length(weights$enc_trunk) < 1L || length(weights$dec_trunk) < 1L) {
    stop("fit_cvae: exported encoder/decoder has no trunk layers", call. = FALSE)
  }

  # sigma2: read in EVAL mode from the pure (softplus) mapping for the learned case;
  # the fixed case returns the hyperparameter. Computed as a float64 scalar.
  sigma2_out <- if (sigma2_learned) {
    raw_r <- reticulate::py_to_r(
      sigma2_raw$detach()$cpu()$to(torch$float64)$numpy())
    raw_r <- as.numeric(raw_r)
    log1p(exp(raw_r))                                          # softplus in float64
  } else {
    as.numeric(hp$decoder_sigma2)
  }

  list(weights = weights, sigma2 = sigma2_out)
}


# ----------------------------------------------------------------------------
# fit_cvae
# ----------------------------------------------------------------------------

#' @title Fit the counterfactual class-conditional VAE single-sample discriminator
#'
#' @description
#' Trains a BatchNorm-free CLASS-conditional VAE (encoder \eqn{q(z|x,y)} + decoder
#' \eqn{p(x|z,y)}, both plain MLPs conditioned on a one-hot of the class label
#' \eqn{y}) on the per-sample-rCLR, universe-aligned TRAINING matrix by the
#' reparameterized ELBO (Kingma & Welling 2014; Sohn et al. 2015; via
#' reticulate-python torch, full-batch Adam, minimizing -ELBO), EXPORTS every weight
#' plus the scalar decoder variance \eqn{\sigma^2} to R as plain numeric
#' matrices/vectors, and DISCARDS the python module. The fitted model holds no
#' external pointer. The score is the counterfactual evidence ratio
#' \eqn{\mathrm{ELBO}(x,\mathrm{case}) - \mathrm{ELBO}(x,\mathrm{control})} with the
#' latent MEAN \eqn{z = \mu} (deterministic, no RNG), in pure base-R.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control. Labels condition
#'   the VAE at fit (the one-hot \eqn{e_y}).
#' @param meta_train Optional per-sample metadata. Accepted for interface uniformity
#'   and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{enc_hidden}
#'   (positive-integer vector of encoder trunk widths; default \code{c(64L, 32L)}),
#'   \code{dec_hidden} (positive-integer vector of decoder trunk widths; default
#'   \code{c(32L, 64L)}), \code{latent_dim} (positive integer latent size \eqn{L};
#'   default \code{16L}), \code{activation} (\code{"relu"} (default) or
#'   \code{"tanh"}), \code{epochs} (number of full-batch training steps, positive
#'   integer; default \code{200L}), \code{lr} (positive Adam learning rate; default
#'   \code{1e-3}), \code{weight_decay} (non-negative Adam L2; default \code{1e-4}),
#'   \code{kl_beta} (non-negative KL coefficient \eqn{\beta}; default \code{1.0}),
#'   \code{decoder_sigma2} (\code{NULL} = learned via softplus of a free parameter
#'   (default), or a positive scalar to fix \eqn{\sigma^2}), \code{min_features}
#'   (feature-overlap floor at scoring, positive integer; default \code{3L}),
#'   \code{device} (\code{"cpu"} (default), \code{"cuda"}, or \code{"auto"};
#'   \code{"cuda"}/\code{"auto"} fall back to CPU with no GPU), and \code{seed}
#'   (integer; default \code{42L}).
#'
#' @return Object of class \code{cvae_model}: a list with \code{feature_universe},
#'   \code{weights} (exported \code{enc_trunk}/\code{enc_mu}/\code{enc_logvar}/
#'   \code{dec_trunk}/\code{dec_xhat}), \code{sigma2} (scalar decoder variance),
#'   \code{kl_beta}, \code{latent_dim}, \code{activation}, \code{device} (resolved),
#'   \code{seed}, and \code{hp}.
#'
#' @details
#' The score is the counterfactual evidence log-likelihood ratio
#' \eqn{s(x) = \mathrm{ELBO}(x,\mathrm{case}) - \mathrm{ELBO}(x,\mathrm{control})}.
#' Departing from the kit-conditional VAE (which conditions on an external
#' cross-sample covariate), \code{cvae} conditions on the CLASS label and evaluates
#' BOTH classes internally, so it needs no kit/cohort input at score and is
#' single-sample-faithful.
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
#' model <- fit_cvae(X, y)
#' score_cvae(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes. \emph{ICLR}.
#' arXiv:1312.6114.
#'
#' Sohn K, Lee H, Yan X. (2015) Learning Structured Output Representation using Deep
#' Conditional Generative Models. \emph{NeurIPS} 28.
#'
#' @export
fit_cvae <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_cvae", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_cvae", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_cvae")
  hp <- .cvae_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_cvae: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's first
  # python init seeds R's RNG as a side effect. Fit draws no R randomness of its own
  # (weight init + reparameterization eps are torch-seeded), so this purely
  # neutralises the init side effect. Restoring on exit (incl. the no-prior-seed
  # removal case) keeps fit globally RNG-neutral.
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

  .cvae_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .cvae_rclr_matrix(X_train)                 # n x p rCLR over full universe

  device <- .cvae_resolve_device(hp$device)
  trained <- .cvae_train_export(Z_in, y, hp, device)

  model <- list(
    feature_universe = feat_order,
    weights = trained$weights,
    sigma2 = trained$sigma2,
    kl_beta = hp$kl_beta,
    latent_dim = hp$latent_dim,
    activation = hp$activation,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "cvae_model"
  model
}


# ----------------------------------------------------------------------------
# score_cvae
# ----------------------------------------------------------------------------

#' @title Score the counterfactual class-conditional VAE single-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_cvae}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), then scored by the counterfactual evidence
#' log-likelihood ratio \eqn{s(x) = \mathrm{ELBO}(x,\mathrm{case}) -
#' \mathrm{ELBO}(x,\mathrm{control})} using the latent MEAN \eqn{z = \mu}
#' (DETERMINISTIC -- no sampling, no RNG). Larger = more case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on that
#' row and the frozen model and is exactly invariant to per-specimen positive scaling.
#'
#' @param model A \code{cvae_model} object from \code{\link{fit_cvae}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative abundances
#'   with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity and
#'   ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values are
#'   more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_cvae(X, y)
#' score_cvae(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes. \emph{ICLR}.
#' arXiv:1312.6114.
#'
#' @export
score_cvae <- function(model, X, meta = NULL) {
  if (!inherits(model, "cvae_model")) {
    stop("score_cvae: model must have class cvae_model")
  }
  X <- .reo_check_matrix(X, "score_cvae", "X")
  .reo_check_meta(meta, nrow(X), "score_cvae", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .cvae_align(X, feat_order)
  weights <- model$weights
  sigma2 <- model$sigma2
  kl_beta <- model$kl_beta
  activation <- model$activation
  p <- length(feat_order)

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0. Test
    # it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition (full
    # positive support, equal abundances) also has an all-zero rCLR (rCLR origin)
    # but is a valid specimen that MUST be encoded/decoded and scored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .cvae_rclr(X_use[i, ])
    out[i] <- .cvae_score_one(z, weights, sigma2, kl_beta, activation, p)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_cvae: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen, score-determining R-side state. The model holds NO live
# external-pointer object (the python module is discarded at fit), so every field is
# a plain R object and this is a deterministic snapshot suitable as the
# `model_digest` for singlesample_assert_row_equivariant() (clause (d): the bytes
# must be identical before and after scoring). Two seed-matched CPU fits on identical
# data produce identical digests.
.cvae_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    weights = model$weights,
    sigma2 = model$sigma2,
    kl_beta = model$kl_beta,
    latent_dim = model$latent_dim,
    activation = model$activation,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
