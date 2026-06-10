#' @title Self-gated mixture of frozen experts (MoE-gated), single-sample
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing a
#' \strong{self-gated mixture of experts} over the specimen's OWN rCLR embedding,
#' adapted to compositional specimens. A small \strong{LayerNorm} MLP encoder of
#' the per-sample-rCLR specimen is trained on the TRAINING set (via
#' reticulate-python torch) to a pre-activation (penultimate) embedding \eqn{z},
#' then \strong{FROZEN}; \eqn{E} linear logistic expert heads each give an expert
#' logit \eqn{h_e(z)}, and a linear gating network with a softmax over the \eqn{E}
#' experts gives mixture weights \eqn{g_e(z)}; the score is the mixture logit
#' \eqn{s(z) = \sum_e g_e(z)\,h_e(z)}. Larger score = more case-like.
#'
#' \strong{Single-sample faithfulness: SELF-gated, NOT kit-gated.} The manifest
#' names this method "kit-gated mixture of frozen experts," where the gate is
#' driven by sequencing-kit / technology metadata. Kit/technology is an EXTERNAL,
#' cross-sample covariate; conditioning the gate on it at inference would require
#' an external label per query and thereby VIOLATE single-sample deployability (no
#' cross-cohort reference, no external covariate at score). The single-sample-
#' faithful design here is therefore \strong{self-gated}: the gate is a function of
#' the specimen's OWN frozen embedding \eqn{z} alone, so a new specimen is routed
#' purely from its own profile with no external signal. This is the design
#' departure flagged for the gate; the MoE structure (experts + a learned soft
#' router) is preserved, only the routing SIGNAL is made endogenous.
#'
#' \strong{Torch is confined to FIT.} The encoder, the \eqn{E} experts, and the
#' gate are trained JOINTLY with torch (Adam, cross-entropy on the mixture logit)
#' on the rCLR'd, universe-aligned training matrix. After training, the encoder
#' Linear weights/biases, the per-layer LayerNorm \eqn{\gamma,\beta,\varepsilon},
#' the \eqn{E} expert weight/bias and the gate weight/bias are EXPORTED to R as
#' plain numeric matrices/vectors; the python module is DISCARDED. The fitted model
#' therefore holds NO external pointer and NO python dependency at score time: the
#' encoder forward (incl. LayerNorm), the expert logits, the gate softmax, and the
#' mixture are all pure base-R. This makes the default \code{model_digest} a stable
#' deterministic snapshot and lets the package suite run the scorer on any host
#' without the venv.
#'
#' \strong{Pre-activation embedding + LayerNorm (single-sample-safe).} The encoder
#' embeds to the PRE-ACTIVATION (penultimate) representation \eqn{z} (the endorsed
#' trained-torch lesson): activation is applied after every LayerNorm EXCEPT the
#' final (embedding) layer, so \eqn{z} is the raw representation the experts and the
#' gate read. \strong{LayerNorm} normalises ACROSS FEATURES WITHIN A ROW, so it is
#' per-row and single-sample-exact; \strong{BatchNorm}, which couples rows, is
#' forbidden here. torch's LayerNorm uses the BIASED variance (1/N), reproduced
#' R-side with \code{rowMeans} of squared deviations.
#'
#' \strong{Fit/score consistency despite float32.} Torch trains in float32; an R
#' float64 forward on the exported (float32-valued) weights differs from the torch
#' forward by \eqn{\sim 10^{-6}}. The exported experts and gate are torch-trained
#' jointly with the encoder, so the mixture is read out from the SAME exported
#' weights the scorer evaluates; the pure-R forward, expert logits, and gate softmax
#' are the deployed computation (the deepmaha / proto-net lesson: fit and score
#' share one code path).
#'
#' \strong{The mixture forward.} With the frozen embedding
#' \eqn{z = \mathrm{enc}(\mathrm{rclr}(x))}, the \eqn{E} expert logits are
#' \eqn{H = z\,W_e^\top + b_e} (length \eqn{E}), the gate logits are
#' \eqn{z\,W_g^\top + b_g}, the mixture weights are the softmax
#' \eqn{g_e(z) = \mathrm{softmax}(z\,W_g^\top + b_g)_e} (a valid distribution over
#' the experts, rows sum to 1, all \eqn{\geq 0}), and the score is the convex
#' combination \eqn{s(z) = \sum_e g_e(z)\,h_e(z)}. Because the gate is a softmax,
#' the mixture logit lies within \eqn{[\min_e h_e, \max_e h_e]} per row. Everything
#' is per-row, so single-row == batch EXACTLY.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly
#' invariant to per-specimen scaling and uses no cross-row statistic. Universe
#' features absent from a specimen carry the neutral rCLR value \code{0}. The SAME
#' rCLR transform is applied to the frozen training rows at fit and to every query
#' at score; combined with LayerNorm (per-row), the frozen experts, and the frozen
#' gate, a specimen's score depends only on that specimen and the frozen model --
#' single-row == batch EXACTLY, and it passes
#' \code{\link{singlesample_assert_row_equivariant}}.
#'
#' \strong{Degenerate-neutral.} A query returns the neutral score \code{0} if fewer
#' than \code{hp$min_features} universe features overlap the query columns (a
#' column-overlap floor, batch-independent), or it has empty positive support over
#' the frozen universe tested on the (pre-rCLR) aligned ORIGINAL abundances
#' (\code{!any(X_use[i, ] > 0)}). A FLAT all-equal-positive composition maps to the
#' rCLR origin but is a VALID specimen and is scored normally (NOT floored) -- the
#' empty-support test is on the ORIGINAL abundances, NOT on \code{all(z == 0)}.
#'
#' @references
#' Jacobs RA, Jordan MI, Nowlan SJ, Hinton GE. (1991) Adaptive Mixtures of Local
#' Experts. \emph{Neural Computation} 3(1):79-87.
#' Shazeer N, Mirhoseini A, Maziarz K, Davis A, Le Q, Hinton G, Dean J. (2017)
#' Outrageously Large Neural Networks: The Sparsely-Gated Mixture-of-Experts Layer.
#' \emph{International Conference on Learning Representations (ICLR)}.
#' arXiv:1701.06538.
#'
#' @name singlesample-moe-gated
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .ssl_vicreg_rclr /
# .unc_sngp_rclr.
.moe_gated_rclr <- function(v) {
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
.moe_gated_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .moe_gated_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.moe_gated_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R encoder forward (exact reimplementation of the torch LayerNorm MLP)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported LayerNorm MLP encoder to the PRE-ACTIVATION
# (penultimate) embedding z. `weights` is a per-layer list(list(W = out x in,
# b = out)); `layernorm` is a per-layer list(list(gamma = out, beta = out,
# eps = scalar)). Each layer l is
#   H <- (H %*% t(W_l) + b_l)                                   # Linear
#   H <- gamma_l * (H - mean(H)) / sqrt(var(H) + eps_l) + beta_l   # LayerNorm
# applied row-wise (mean/var ACROSS FEATURES WITHIN A ROW -- per-row, batch-free),
# and `act` is applied after every LayerNorm EXCEPT the final (embedding) layer
# (matching the torch encoder; the embedding is the raw pre-activation z). torch's
# LayerNorm uses the BIASED variance (1/N), so we use rowMeans of squared
# deviations (NOT stats::var's 1/(N-1)). Returns the n x d embedding. Per-row exact,
# so single-row == batch.
.moe_gated_forward <- function(M, weights, layernorm, activation) {
  act <- if (identical(activation, "tanh")) tanh else function(x) pmax(x, 0)
  L <- length(weights)
  H <- M
  for (l in seq_len(L)) {
    W <- weights[[l]]$W                       # out x in
    b <- weights[[l]]$b                       # length out
    H <- H %*% t(W)                           # n x out
    H <- sweep(H, 2L, b, `+`)
    ln <- layernorm[[l]]
    g <- ln$gamma; be <- ln$beta; eps <- ln$eps
    mu <- rowMeans(H)                         # mean ACROSS FEATURES per row
    Hc <- H - mu                              # centred (recycles mu down columns)
    v  <- rowMeans(Hc * Hc)                   # biased variance (1/N), per row
    Hn <- Hc / sqrt(v + eps)                  # standardised per row
    H  <- sweep(Hn, 2L, g, `*`)               # gamma * Hn
    H  <- sweep(H, 2L, be, `+`)               # + beta
    if (l < L) H <- act(H)                    # no activation after the embedding
  }
  H
}


# Pure-R self-gated mixture forward on an n x d embedding Z. `W_e` is E x d, `b_e`
# length E (expert logits H = Z W_e^T + b_e, n x E); `W_g` is E x d, `b_g` length E
# (gate logits Z W_g^T + b_g, n x E). The gate weights are the row-wise softmax of
# the gate logits (stable max-shift), and the score is the convex combination
# rowSums(G * H). Returns a list(score = length-n, G = n x E gate weights,
# H = n x E expert logits) so callers (and tests) can inspect the gate/experts.
.moe_gated_mixture <- function(Z, W_e, b_e, W_g, b_g) {
  H <- Z %*% t(W_e)                           # n x E expert logits
  H <- sweep(H, 2L, b_e, `+`)
  Glogit <- Z %*% t(W_g)                      # n x E gate logits
  Glogit <- sweep(Glogit, 2L, b_g, `+`)
  m <- apply(Glogit, 1L, max)                 # per-row max for numeric stability
  Eexp <- exp(Glogit - m)                     # n x E
  G <- Eexp / rowSums(Eexp)                   # row-wise softmax (rows sum to 1)
  score <- rowSums(G * H)                     # convex combination
  list(score = as.numeric(score), G = G, H = H)
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.moe_gated_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_moe_gated: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_moe_gated: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_moe_gated: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("hidden", "activation", "n_experts", "epochs", "lr",
               "weight_decay", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_moe_gated: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  hidden <- hp$hidden
  if (is.null(hidden)) hidden <- c(64L, 32L)
  if (!is.numeric(hidden) || length(hidden) < 1L || any(!is.finite(hidden)) ||
      any(hidden < 1L) || any(hidden > .Machine$integer.max) ||
      any(hidden != as.integer(hidden))) {
    stop("fit_moe_gated: hp$hidden must be a vector of positive integers")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_moe_gated: hp$activation must be one of 'relu', 'tanh'")
  }

  n_experts <- hp$n_experts
  if (is.null(n_experts)) n_experts <- 4L
  if (!is.numeric(n_experts) || length(n_experts) != 1L ||
      !is.finite(n_experts) || n_experts < 2L ||
      n_experts > .Machine$integer.max ||
      n_experts != as.integer(n_experts)) {
    stop("fit_moe_gated: hp$n_experts must be an integer >= 2")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_moe_gated: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_moe_gated: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_moe_gated: hp$weight_decay must be a non-negative finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_moe_gated: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_moe_gated: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_moe_gated: hp$seed must be a single integer")
  }

  list(
    hidden = as.integer(hidden),
    activation = activation,
    n_experts = as.integer(n_experts),
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
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
.moe_gated_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_moe_gated: package 'reticulate' is required to train the ",
         "encoder. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_moe_gated: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.moe_gated_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# Self-gated MoE training (encoder + experts + gate) + weight export
# ----------------------------------------------------------------------------

# Train the LayerNorm MLP encoder + E linear expert heads + a linear gate JOINTLY
# via torch (Adam, cross-entropy on the mixture logit) and return the exported
# R-side encoder weights (per-layer list(W, b)), the per-layer LayerNorm parameters
# (list(gamma, beta, eps)), the E expert weight/bias (W_e E x d, b_e length E), and
# the gate weight/bias (W_g E x d, b_g length E). The python module is used here
# only; nothing python is returned. The torch CPU/CUDA RNG (for weight init +
# training) is saved/restored, so fit leaves the torch generators byte-unchanged;
# the global R .Random.seed is saved/restored by the caller. `Z` is the n x p
# rCLR'd, universe-aligned training matrix; `y` the integer 0/1 labels.
#
# SEED BEFORE BUILD: torch$manual_seed(hp$seed) is called BEFORE constructing ANY
# stochastic module (encoder Linears, LayerNorms, experts, gate). Constructing a
# module before seeding makes the model cross-SESSION non-reproducible even though
# a within-process two-fit digest test passes (the on.exit RNG restore masks it).
.moe_gated_train_export <- function(Z, y, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn
  np <- reticulate::import("numpy", delay_load = FALSE)

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

  # SEED BEFORE BUILD: seed the global torch generator BEFORE constructing any
  # stochastic module so weight init (encoder + LayerNorm + experts + gate) is
  # deterministic from hp$seed alone (cross-session reproducible, not just
  # cross-fit-within-process).
  torch$manual_seed(as.integer(hp$seed))

  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)
  y_py <- torch$as_tensor(np$asarray(as.integer(y), dtype = "int64"))$long()$to(device)

  p <- ncol(Z)
  hidden <- hp$hidden
  E <- as.integer(hp$n_experts)

  # Encoder: per layer a Linear + a LayerNorm; activation after every layer EXCEPT
  # the final (embedding) layer, so the embedding is the raw PRE-ACTIVATION z.
  lins <- list()
  lns <- list()
  in_dim <- p
  for (h in hidden) {
    ho <- as.integer(h)
    lins[[length(lins) + 1L]] <- nn$Linear(in_dim, ho)
    lns[[length(lns) + 1L]] <- nn$LayerNorm(ho)
    in_dim <- ho
  }
  d <- in_dim

  # E linear logistic experts (each Linear(d -> 1)) and a linear gate (Linear(d ->
  # E)); the score is the gate-softmax-weighted mixture of the expert logits. A
  # single 2-logit shape is avoided: we keep each expert a SCALAR logit so the
  # mixture logit feeds a sigmoid-style CE. To use CrossEntropyLoss (2-logit) we
  # form a 2-column [-s/2, +s/2]-style output below.
  experts <- nn$Linear(d, E)                  # E expert logits in one Linear
  gate    <- nn$Linear(d, E)                  # E gate logits

  # Manual forward so the LayerNorm placement, the pre-activation embedding, and
  # the gate-softmax mixture match the pure-R .moe_gated_forward / .moe_gated_mixture
  # EXACTLY.
  L <- length(lins)
  act_fun <- if (identical(hp$activation, "tanh")) torch$tanh else torch$relu
  enc_fwd <- function(xb) {
    H <- xb
    for (l in seq_len(L)) {
      H <- lins[[l]](H)                       # Linear
      H <- lns[[l]](H)                        # LayerNorm
      if (l < L) H <- act_fun(H)              # no activation after the embedding
    }
    H                                         # pre-activation embedding z
  }
  mix_fwd <- function(xb) {
    z <- enc_fwd(xb)
    Hlog <- experts(z)                        # n x E expert logits
    Glog <- gate(z)                           # n x E gate logits
    Gw <- nn$functional$softmax(Glog, dim = 1L)   # row-wise softmax over experts
    torch$sum(torch$mul(Gw, Hlog), dim = 1L)  # length-n mixture logit
  }

  # Register every module's parameters for joint optimisation.
  params <- list()
  for (lin in lins) params <- c(params, reticulate::iterate(lin$parameters()))
  for (ln in lns)   params <- c(params, reticulate::iterate(ln$parameters()))
  params <- c(params, reticulate::iterate(experts$parameters()))
  params <- c(params, reticulate::iterate(gate$parameters()))
  opt <- torch$optim$Adam(params, lr = hp$lr, weight_decay = hp$weight_decay)
  lossf <- nn$BCEWithLogitsLoss()
  y_f <- torch$as_tensor(np$asarray(as.numeric(y), dtype = "float32"))$float()$to(device)

  for (lin in lins) lin$train()
  for (ln in lns)   ln$train()
  experts$train(); gate$train()
  n <- nrow(Z)
  bs <- min(32L, n)
  for (ep in seq_len(hp$epochs)) {
    starts <- seq.int(0L, n - 1L, by = bs)
    for (s0 in starts) {
      idx <- torch$arange(as.integer(s0),
                          as.integer(min(s0 + bs, n)),
                          dtype = torch$long)$to(device)
      xb <- X_py$index_select(0L, idx)
      yb <- y_f$index_select(0L, idx)
      opt$zero_grad()
      s_mix <- mix_fwd(xb)                     # length-bs mixture logit
      loss <- lossf(s_mix, yb)
      loss$backward()
      opt$step()
    }
  }
  for (lin in lins) lin$eval()
  for (ln in lns)   ln$eval()
  experts$eval(); gate$eval()

  # ---- export encoder weights + LayerNorm params ----------------------------
  weights <- vector("list", L)
  layernorm <- vector("list", L)
  for (l in seq_len(L)) {
    lin <- lins[[l]]
    W <- reticulate::py_to_r(lin$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(lin$bias$detach()$cpu()$to(torch$float64)$numpy())
    weights[[l]] <- list(W = W, b = as.numeric(b))
    ln <- lns[[l]]
    g  <- reticulate::py_to_r(ln$weight$detach()$cpu()$to(torch$float64)$numpy())
    be <- reticulate::py_to_r(ln$bias$detach()$cpu()$to(torch$float64)$numpy())
    eps <- as.numeric(reticulate::py_to_r(ln$eps))
    layernorm[[l]] <- list(gamma = as.numeric(g), beta = as.numeric(be),
                           eps = eps)
  }

  # ---- export the E experts and the gate (W E x d, b length E) ---------------
  W_e <- reticulate::py_to_r(experts$weight$detach()$cpu()$to(torch$float64)$numpy())
  b_e <- as.numeric(reticulate::py_to_r(experts$bias$detach()$cpu()$to(torch$float64)$numpy()))
  W_g <- reticulate::py_to_r(gate$weight$detach()$cpu()$to(torch$float64)$numpy())
  b_g <- as.numeric(reticulate::py_to_r(gate$bias$detach()$cpu()$to(torch$float64)$numpy()))
  if (!is.matrix(W_e)) W_e <- matrix(W_e, nrow = E)
  if (!is.matrix(W_g)) W_g <- matrix(W_g, nrow = E)

  list(weights = weights, layernorm = layernorm, embedding_dim = d,
       W_e = W_e, b_e = b_e, W_g = W_g, b_g = b_g)
}


# ----------------------------------------------------------------------------
# fit_moe_gated
# ----------------------------------------------------------------------------

#' @title Fit the self-gated mixture-of-experts single-sample discriminator
#'
#' @description
#' Trains a LayerNorm MLP encoder, \eqn{E} linear logistic expert heads, and a
#' linear gating network JOINTLY on the per-sample-rCLR, universe-aligned TRAINING
#' matrix (via reticulate-python torch; Adam, cross-entropy on the gate-softmax
#' mixture logit), EXPORTS the encoder weights, the LayerNorm parameters, the
#' \eqn{E} expert weight/bias, and the gate weight/bias to R as plain numeric
#' matrices/vectors, and DISCARDS the python module. The fitted model holds no
#' external pointer; the score is the pure-R mixture logit
#' \eqn{s(z) = \sum_e g_e(z)\,h_e(z)} with \eqn{z = \mathrm{enc}(\mathrm{rclr}(x))}.
#'
#' The gate is driven by the specimen's OWN embedding (SELF-gated), NOT by external
#' kit/technology metadata, so the method is single-sample-deployable (no external
#' covariate at score). See the method description for this design departure from
#' the manifest's "kit-gated" name.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored (the gate is self-driven, not kit-driven).
#' @param hp Optional list of hyperparameters. Allowed fields: \code{hidden}
#'   (positive-integer vector of encoder layer widths; the LAST entry is the
#'   embedding dim \eqn{d}; default \code{c(64L, 32L)}), \code{activation}
#'   (\code{"relu"} (default) or \code{"tanh"}), \code{n_experts} (number of experts
#'   \eqn{E}, integer \eqn{\geq 2}; default \code{4L}), \code{epochs} (positive
#'   integer; default \code{200L}), \code{lr} (positive Adam learning rate; default
#'   \code{1e-3}), \code{weight_decay} (non-negative Adam L2; default \code{1e-4}),
#'   \code{min_features} (feature-overlap floor at scoring, positive integer;
#'   default \code{3L}), \code{device} (\code{"cpu"} (default), \code{"cuda"}, or
#'   \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU with no GPU), and
#'   \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{moe_gated_model}: a list with
#'   \code{feature_universe}, \code{weights} (exported per-layer \code{list(W, b)}
#'   encoder), \code{layernorm} (per-layer \code{list(gamma, beta, eps)}),
#'   \code{activation}, \code{W_e} (\eqn{E \times d} expert weights), \code{b_e}
#'   (length-\eqn{E} expert biases), \code{W_g} (\eqn{E \times d} gate weights),
#'   \code{b_g} (length-\eqn{E} gate biases), \code{n_experts}, \code{embedding_dim},
#'   \code{device} (resolved), \code{seed}, and \code{hp}.
#'
#' @details
#' The score is the convex combination \eqn{s(z) = \sum_e g_e(z)\,h_e(z)} where the
#' expert logits are \eqn{h_e(z) = (z W_e^\top + b_e)_e} and the gate weights are
#' the softmax \eqn{g_e(z) = \mathrm{softmax}(z W_g^\top + b_g)_e}. Because the gate
#' is a softmax, the mixture logit lies within \eqn{[\min_e h_e, \max_e h_e]} per
#' row. Larger = more case-like.
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
#' model <- fit_moe_gated(X, y)
#' score_moe_gated(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Jacobs RA, Jordan MI, Nowlan SJ, Hinton GE. (1991) Adaptive Mixtures of Local
#' Experts. \emph{Neural Computation} 3(1):79-87.
#' Shazeer N, et al. (2017) Outrageously Large Neural Networks: The Sparsely-Gated
#' Mixture-of-Experts Layer. \emph{ICLR}. arXiv:1701.06538.
#'
#' @export
fit_moe_gated <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_moe_gated", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_moe_gated", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_moe_gated")
  hp <- .moe_gated_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_moe_gated: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's
  # first python init seeds R's RNG as a side effect. Fit draws no R randomness of
  # its own (all torch stochasticity is torch-seeded), so restoring on exit (incl.
  # the no-prior-seed case) keeps fit globally RNG-neutral.
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

  .moe_gated_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .moe_gated_rclr_matrix(X_train)           # n x p rCLR over full universe

  device <- .moe_gated_resolve_device(hp$device)
  ex <- .moe_gated_train_export(Z_in, y, hp, device)

  model <- list(
    feature_universe = feat_order,
    weights = ex$weights,
    layernorm = ex$layernorm,
    activation = hp$activation,
    W_e = ex$W_e,
    b_e = ex$b_e,
    W_g = ex$W_g,
    b_g = ex$b_g,
    n_experts = as.integer(hp$n_experts),
    embedding_dim = ex$embedding_dim,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "moe_gated_model"
  model
}


# Self-gated mixture logit for one embedding z (length d): rowSums(softmax(z W_g^T +
# b_g) * (z W_e^T + b_e)). Larger = more case-like. Thin wrapper over the matrix
# mixture for the single-row scoring path.
.moe_gated_score_one <- function(z, W_e, b_e, W_g, b_g) {
  .moe_gated_mixture(matrix(z, nrow = 1L), W_e, b_e, W_g, b_g)$score
}


# ----------------------------------------------------------------------------
# score_moe_gated
# ----------------------------------------------------------------------------

#' @title Score the self-gated mixture-of-experts single-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_moe_gated}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), embedded by the exported LayerNorm
#' encoder forward to the pre-activation embedding \eqn{z}, and scored by the
#' self-gated mixture logit \eqn{s(z) = \sum_e g_e(z)\,h_e(z)}. Larger = more
#' case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on
#' that row and the frozen model and is exactly invariant to per-specimen positive
#' scaling.
#'
#' @param model A \code{moe_gated_model} object from \code{\link{fit_moe_gated}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity and
#'   ignored by this method (the gate is self-driven).
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values are
#'   more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_moe_gated(X, y)
#' score_moe_gated(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Jacobs RA, et al. (1991) Adaptive Mixtures of Local Experts. \emph{Neural
#' Computation} 3(1):79-87. Shazeer N, et al. (2017) Outrageously Large Neural
#' Networks. \emph{ICLR}. arXiv:1701.06538.
#'
#' @export
score_moe_gated <- function(model, X, meta = NULL) {
  if (!inherits(model, "moe_gated_model")) {
    stop("score_moe_gated: model must have class moe_gated_model")
  }
  X <- .reo_check_matrix(X, "score_moe_gated", "X")
  .reo_check_meta(meta, nrow(X), "score_moe_gated", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .moe_gated_align(X, feat_order)
  weights <- model$weights
  layernorm <- model$layernorm
  activation <- model$activation
  W_e <- model$W_e; b_e <- model$b_e
  W_g <- model$W_g; b_g <- model$b_g

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (rCLR
    # origin) but is a valid specimen that MUST be embedded and scored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .moe_gated_rclr(X_use[i, ])
    emb <- .moe_gated_forward(matrix(z, nrow = 1L), weights, layernorm, activation)
    out[i] <- .moe_gated_score_one(as.numeric(emb), W_e, b_e, W_g, b_g)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_moe_gated: scorer produced non-finite or wrong-length output")
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
# identical data produce identical digests (encoder + LayerNorm + experts + gate
# all seeded/deterministic).
.moe_gated_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    weights = model$weights,
    layernorm = model$layernorm,
    activation = model$activation,
    W_e = model$W_e,
    b_e = model$b_e,
    W_g = model$W_g,
    b_g = model$b_g,
    n_experts = model$n_experts,
    embedding_dim = model$embedding_dim,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
