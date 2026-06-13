#' @title SCARF self-supervised contrastive single-sample discriminator (frozen linear head)
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing the SCARF
#' self-supervised contrastive representation of Bahri, Jiang, Tay & Metzler
#' (2022), adapted to compositional specimens. A small BatchNorm-free MLP encoder
#' of the per-sample-rCLR specimen is learned on the TRAINING set by
#' SELF-SUPERVISED CONTRASTIVE training (InfoNCE / NT-Xent, NO labels) via
#' reticulate-python torch, then \strong{FROZEN}; a FROZEN linear-probe head
#' (the SCARF evaluation protocol) is fit on the labeled training embeddings, and
#' a new specimen is scored by its frozen linear-head logit on its pure-R encoder
#' embedding. Larger score = more case-like.
#'
#' \strong{SCARF contrastive training (Bahri 2022).} Unlike the cross-entropy
#' \code{lrt-deepmaha} or the prototypical \code{proto-net} siblings, the encoder
#' here is trained by the SCARF InfoNCE loss WITHOUT labels. For a batch of the
#' (rCLR'd, universe-aligned) training rows, a CORRUPTED view \eqn{\tilde x} of
#' each row is built (SCARF eq. 1): independently for each feature, with
#' probability \code{corrupt_p} the entry is replaced by a value sampled from that
#' feature's EMPIRICAL TRAINING MARGINAL (sampled uniformly from the training
#' column's observed rCLR values). The positive pair is the row and its corrupted
#' view, both embedded through encoder then projection head
#' (\eqn{z_i = g(f(x_i))}, \eqn{\tilde z_i = g(f(\tilde x_i))}); the NT-Xent
#' contrastive loss with cosine similarity and temperature \eqn{\tau} has the
#' positive pair \eqn{(z_i, \tilde z_i)} as the numerator and the other batch
#' items as negatives, with the partition function summing over ALL \eqn{k}
#' (including the positive \eqn{k = i}, per SCARF Algorithm 1):
#' \deqn{\ell = \frac1B \sum_i -\log
#'   \frac{\exp(\mathrm{sim}(z_i,\tilde z_i)/\tau)}
#'        {\sum_{k}\exp(\mathrm{sim}(z_i,\tilde z_k)/\tau)}.}
#' This is exactly softmax cross-entropy with diagonal labels. The corruption uses
#' an INDEPENDENT \code{Bernoulli(corrupt_p)} mask per feature (expected
#' \eqn{\lfloor corrupt\_p \cdot M\rfloor} corrupted features per row), a common
#' variant of SCARF Algorithm 1's fixed-count \eqn{q = \lfloor corrupt\_p \cdot M\rfloor}
#' per-sample selection. Adam optimises encoder + projection head over
#' \code{epochs} steps. The
#' projection head is discarded after training (standard contrastive practice);
#' only the ENCODER (the pre-projection representation) is exported.
#'
#' \strong{The corruption marginal is FIT-ONLY.} The per-feature empirical
#' marginal (the frozen training column) is used ONLY at FIT, to build training
#' views. It is NEVER stored in the model and NEVER consulted at score: the score
#' does a CLEAN encoder forward on the query, with no corruption. Single-sample
#' deployability therefore holds -- a query's score depends only on that query and
#' the frozen encoder + frozen linear head, with no cross-row statistic. Corruption
#' randomness (the Bernoulli mask and the marginal-resample row indices) is drawn
#' via torch's SEEDED generator, so two seed-matched CPU fits build bit-identical
#' training views.
#'
#' \strong{Frozen linear head (SCARF evaluation protocol).} After SSL training the
#' encoder is FROZEN. A single BN-free torch \code{Linear(d -> 1)} head is then fit
#' by seeded, fixed-epoch, \code{head_weight_decay}-regularised
#' binary-cross-entropy optimisation on the frozen encoder embeddings of the
#' labeled training rows (the L2 regularisation handles separation on small
#' \eqn{n}); its weight \eqn{w} (length \eqn{d}) and bias \eqn{b} (scalar) are
#' EXPORTED to R. CRITICAL fit/score consistency: the embeddings the head trains on
#' are computed from the SAME pure-R \code{.ai_scarf_forward} the scorer uses (NOT
#' the float32 torch embedding), so the head lives in exactly the embedding the
#' scorer evaluates (the deepmaha / proto-net lesson: fit and score share one
#' embedding code path).
#'
#' \strong{Torch is confined to FIT.} The encoder and head are trained with torch
#' on the rCLR'd, universe-aligned training matrix. After training, the encoder
#' weight matrices/biases and the linear-head \eqn{w, b} are EXPORTED to R as plain
#' numeric matrices/vectors; the projection head and the python module are
#' DISCARDED. The fitted model therefore holds NO external pointer and NO python
#' dependency at score time: the encoder forward and the linear-head logit are pure
#' base-R. This makes the default \code{model_digest} a stable deterministic
#' snapshot and lets the package suite run the scorer on any host without the venv.
#'
#' \strong{The encoder forward.} The exported net is
#' \eqn{\mathrm{input}(p) \to [\mathrm{Linear} \to \phi]_{1}
#'   \to \dots \to [\mathrm{Linear} \to \phi]_{L-1}
#'   \to \mathrm{Linear} \to \mathrm{embedding}(d)}. There are \code{length(hidden)}
#' linear layers up to the embedding (the last \code{hidden} width is the embedding
#' dim \eqn{d}); the activation (\code{"relu"} default, or \code{"tanh"}) is applied
#' AFTER every linear layer EXCEPT the final embedding layer (no activation follows
#' the embedding -- it is the raw representation, matching the torch encoder).
#' \code{.ai_scarf_forward(M, weights)} reimplements this exactly.
#'
#' \strong{Linear-head logit score.} With the frozen head the score of an embedded
#' specimen \eqn{z} is \eqn{s(z) = w \cdot z + b}, the linear-probe logit. Larger =
#' more case-like. Encoder forward + linear head are per-row with no cross-row
#' coupling, so the single-row score equals the batch score EXACTLY.
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
#' \strong{eval-frozen-BN.} The encoder is BN-FREE (the cleanest, exact
#' single-sample forward -- same as deepmaha / proto-net). The manifest's
#' "eval-frozen-BN" alternative is acceptable only if BN running stats are frozen
#' and exported so the forward is per-row independent; BN-free is simpler and
#' exactly per-row, so this method uses BN-free.
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
#' Bahri D, Jiang H, Tay Y, Metzler D. (2022) SCARF: Self-Supervised Contrastive
#' Learning using Random Feature Corruption. \emph{International Conference on
#' Learning Representations (ICLR)}. arXiv:2106.15147.
#'
#' @name singlesample-scarf
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .proto_net_rclr.
.ai_scarf_rclr <- function(v) {
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
.ai_scarf_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .ai_scarf_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.ai_scarf_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R encoder forward (exact reimplementation of the torch encoder up to d)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported MLP encoder up to (and including) the
# embedding layer. `weights` is a list of per-layer lists list(W = out x in, b =
# out): layers 1..(L-1) are hidden (Linear then `activation`), layer L is the
# embedding Linear with NO activation (matches the torch encoder exactly). For a
# layer with weight W (out x in) and bias b, torch computes M %*% t(W) + b
# (row = sample). Returns an n x d embedding.
.ai_scarf_forward <- function(M, weights, activation) {
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

.ai_scarf_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ai_scarf: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_ai_scarf: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_ai_scarf: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("hidden", "activation", "epochs", "lr", "weight_decay",
               "corrupt_p", "temperature", "proj_dim", "head_epochs",
               "head_weight_decay", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ai_scarf: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  hidden <- hp$hidden
  if (is.null(hidden)) hidden <- c(64L, 32L)
  if (!is.numeric(hidden) || length(hidden) < 1L || any(!is.finite(hidden)) ||
      any(hidden < 1L) || any(hidden > .Machine$integer.max) ||
      any(hidden != as.integer(hidden))) {
    stop("fit_ai_scarf: hp$hidden must be a vector of positive integers")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_ai_scarf: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_ai_scarf: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_ai_scarf: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_ai_scarf: hp$weight_decay must be a non-negative finite number")
  }

  corrupt_p <- hp$corrupt_p
  if (is.null(corrupt_p)) corrupt_p <- 0.6
  if (!is.numeric(corrupt_p) || length(corrupt_p) != 1L ||
      !is.finite(corrupt_p) || corrupt_p <= 0 || corrupt_p >= 1) {
    stop("fit_ai_scarf: hp$corrupt_p must be a single number in (0, 1)")
  }

  temperature <- hp$temperature
  if (is.null(temperature)) temperature <- 1.0
  if (!is.numeric(temperature) || length(temperature) != 1L ||
      !is.finite(temperature) || temperature <= 0) {
    stop("fit_ai_scarf: hp$temperature must be a positive finite number")
  }

  proj_dim <- hp$proj_dim
  if (is.null(proj_dim)) proj_dim <- 64L
  if (!is.numeric(proj_dim) || length(proj_dim) != 1L || !is.finite(proj_dim) ||
      proj_dim < 1L || proj_dim > .Machine$integer.max ||
      proj_dim != as.integer(proj_dim)) {
    stop("fit_ai_scarf: hp$proj_dim must be a positive integer")
  }

  head_epochs <- hp$head_epochs
  if (is.null(head_epochs)) head_epochs <- 200L
  if (!is.numeric(head_epochs) || length(head_epochs) != 1L ||
      !is.finite(head_epochs) || head_epochs < 1L ||
      head_epochs > .Machine$integer.max ||
      head_epochs != as.integer(head_epochs)) {
    stop("fit_ai_scarf: hp$head_epochs must be a positive integer")
  }

  head_weight_decay <- hp$head_weight_decay
  if (is.null(head_weight_decay)) head_weight_decay <- 1e-2
  if (!is.numeric(head_weight_decay) || length(head_weight_decay) != 1L ||
      !is.finite(head_weight_decay) || head_weight_decay < 0) {
    stop("fit_ai_scarf: hp$head_weight_decay must be a non-negative finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_ai_scarf: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_ai_scarf: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_ai_scarf: hp$seed must be a single integer")
  }

  list(
    hidden = as.integer(hidden),
    activation = activation,
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    corrupt_p = as.numeric(corrupt_p),
    temperature = as.numeric(temperature),
    proj_dim = as.integer(proj_dim),
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
# unavailable. SCARF contrastive training needs NO brew-OpenSSL bridge and NO HF
# download, so only reticulate + `import torch` are probed.
.ai_scarf_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_ai_scarf: package 'reticulate' is required to train the ",
         "encoder. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_ai_scarf: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.ai_scarf_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# SCARF contrastive training + frozen linear head + weight export
# ----------------------------------------------------------------------------

# Train the BN-free MLP encoder by SCARF SELF-SUPERVISED CONTRASTIVE training
# (Bahri 2022; InfoNCE / NT-Xent, NO labels) via torch, FREEZE it, fit a frozen
# linear-probe head on the labeled training embeddings, and return the exported
# R-side encoder weights plus the head (w, b). The python module is used here
# only; nothing python is returned. The torch RNG (CPU + guarded CUDA) is
# saved/restored, so fit leaves the torch generators byte-unchanged; the global R
# .Random.seed is saved/restored by the caller. `Z` is the n x p rCLR'd,
# universe-aligned training matrix; `y` the integer 0/1 labels; `forward_fn` the
# pure-R encoder forward used to recompute head-training embeddings (fit/score
# consistency). All corruption randomness is drawn via torch's seeded generator.

# InfoNCE / NT-Xent contrastive loss for an n x n clean-vs-corrupted cosine-
# similarity matrix `logits` (already divided by the temperature). The positive
# pair is the diagonal sim(zc_i, zk_i); the partition function sums over ALL k
# INCLUDING the positive k = i (SCARF Algorithm 1, Bahri 2022), so the diagonal
# stays in the denominator. This is exactly softmax cross-entropy with diagonal
# labels: -log( exp(sim_ii) / sum_k exp(sim_ik) ). Exposed as a helper so the
# contrastive loss itself is directly unit-testable (vs nnf cross_entropy).
.ai_scarf_infonce_loss <- function(logits) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  pos <- logits$diagonal()                           # length n: sim_ii / tau
  denom <- torch$logsumexp(logits, dim = 1L)         # log sum_k exp(sim_ik), positive INCLUDED
  torch$mean(torch$sub(denom, pos))                  # mean_i -log( e^{sim_ii} / sum_k e^{sim_ik} )
}

.ai_scarf_train_export <- function(Z, y, hp, device, forward_fn) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn
  F <- nn$functional

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

  np <- reticulate::import("numpy", delay_load = FALSE)
  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)
  n <- nrow(Z)
  p <- ncol(Z)
  hidden <- hp$hidden

  # ---- SSL STAGE seed (BEFORE any module construction) ----------------------
  # Seed the torch generators HERE, before the encoder Linears are created, so
  # the encoder init AND the projection-head init AND the corruption draws (mask
  # + marginal-resample indices) all replay identically for a given seed --
  # including ACROSS processes. Seeding after the encoder is built would leave
  # the encoder weights initialised from the unseeded process-entry RNG state,
  # which is not reproducible across sessions (two in-process fits would still
  # match via the on.exit RNG restore, masking the defect in a single session).
  torch$manual_seed(as.integer(hp$seed))

  # ---- build the BN-free encoder: input(p) -> [Linear -> act] x (L-1) -> ----
  # ---- Linear -> embedding(d). The last `hidden` width is the embedding dim. -
  act_mod <- if (identical(hp$activation, "tanh")) nn$Tanh else nn$ReLU
  enc_layers <- list()
  in_dim <- p
  for (h in hidden) {
    enc_layers[[length(enc_layers) + 1L]] <- nn$Linear(in_dim, as.integer(h))
    enc_layers[[length(enc_layers) + 1L]] <- act_mod()
    in_dim <- as.integer(h)
  }
  enc_layers[[length(enc_layers)]] <- NULL    # drop activation after embedding
  d <- in_dim

  # ---- SSL STAGE: build encoder + projection head, contrastive train --------
  encoder <- do.call(nn$Sequential, enc_layers)$to(device)
  # Projection head g: Linear(d -> proj) -> ReLU -> Linear(proj -> proj). It maps
  # the embedding to the contrastive space; it is DISCARDED after SSL (only the
  # encoder is exported).
  proj <- nn$Sequential(
    nn$Linear(d, as.integer(hp$proj_dim)),
    nn$ReLU(),
    nn$Linear(as.integer(hp$proj_dim), as.integer(hp$proj_dim))
  )$to(device)

  ssl_params <- c(reticulate::iterate(encoder$parameters()),
                  reticulate::iterate(proj$parameters()))
  opt <- torch$optim$Adam(ssl_params, lr = hp$lr,
                          weight_decay = hp$weight_decay)

  tau <- hp$temperature
  cp <- hp$corrupt_p

  encoder$train(); proj$train()
  for (ep in seq_len(hp$epochs)) {
    # ---- SCARF corruption (eq. 1), torch-seeded for reproducibility ---------
    # Bernoulli(corrupt_p) mask over the full n x p batch.
    mask <- torch$bernoulli(torch$full(list(n, p), cp))$to(device)
    # Marginal resample: for each (i, j) draw a row index uniformly in 0..n-1
    # and take that column's value -> a draw from feature j's empirical TRAINING
    # marginal. randint uses the seeded generator, so the draw is reproducible.
    ridx <- torch$randint(0L, as.integer(n), list(n, p))$to(device)$long()
    # gather column j of X_py at the drawn rows: resampled[i, j] = X[ridx[i,j], j].
    resampled <- torch$gather(X_py, 0L, ridx)        # n x p
    x_corr <- torch$add(torch$mul(mask, resampled),
                        torch$mul(torch$sub(torch$ones_like(mask), mask), X_py))

    # ---- embed clean + corrupted view through encoder -> projection head ----
    z_clean <- proj(encoder(X_py))                   # n x proj
    z_corr  <- proj(encoder(x_corr))                 # n x proj
    zc <- F$normalize(z_clean, dim = 1L)             # L2-normalised (cosine)
    zk <- F$normalize(z_corr,  dim = 1L)

    # ---- NT-Xent / InfoNCE (SCARF Algorithm 1): positive = (zc_i, zk_i), negatives
    # = the off-diagonal zk_{k != i}; the partition function sums over ALL k including
    # the positive (= softmax cross-entropy with diagonal labels).
    logits <- torch$matmul(zc, zk$t())$div(tau)      # n x n cosine sims / tau
    loss <- .ai_scarf_infonce_loss(logits)

    opt$zero_grad()
    loss$backward()
    opt$step()
    # Release the per-epoch reticulate tensor WRAPPERS so the Python-side tensors they
    # hold are freed. reticulate releases a Python object only when R garbage-collects
    # its R-side external-pointer wrapper; without a periodic gc the n x n InfoNCE
    # intermediates (logits/zc/zk/x_corr/...) accumulate across all epochs on a large
    # training fold, leaking tens-to-hundreds of GB on the biggest cohorts (observed
    # ~210 GB on one cohort -> OOM). gc() is computation-neutral -- it touches no live
    # training state (net/opt/current tensors) -- so the fitted weights are bit-identical
    # (the §7.8 determinism test confirms an unchanged model digest).
    if (ep %% 10L == 0L) gc(FALSE)
  }
  encoder$eval(); proj$eval()

  # ---- export encoder weights up to and INCLUDING the embedding layer -------
  enc_modules <- reticulate::iterate(encoder$modules())
  weights <- list()
  for (m in enc_modules) {
    has_w <- tryCatch(!is.null(m$weight), error = function(e) FALSE)
    if (!isTRUE(has_w)) next
    W <- reticulate::py_to_r(m$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(m$bias$detach()$cpu()$to(torch$float64)$numpy())
    weights[[length(weights) + 1L]] <- list(W = W, b = as.numeric(b))
  }
  if (length(weights) < 1L) {
    stop("fit_ai_scarf: exported encoder has no linear layers", call. = FALSE)
  }

  # ---- HEAD STAGE: frozen linear probe on the PURE-R encoder embeddings ------
  # Embeddings are recomputed with the pure-R forward the scorer uses (NOT the
  # float32 torch embedding), so the head lives in exactly the scorer's embedding
  # (fit/score consistency). Re-seed torch so the head init + optimisation replay
  # identically for a given seed (two-stage determinism).
  Z_emb <- forward_fn(Z, weights, hp$activation)     # n x d (pure-R, float64)
  torch$manual_seed(as.integer(hp$seed))
  emb_py <- torch$as_tensor(np$asarray(Z_emb, dtype = "float32"))$float()$to(device)
  yt <- torch$as_tensor(np$asarray(as.numeric(y), dtype = "float32"))$float()$to(device)
  head <- nn$Linear(as.integer(d), 1L)$to(device)
  head$train()
  opt_h <- torch$optim$Adam(head$parameters(), lr = hp$lr,
                            weight_decay = hp$head_weight_decay)
  bce <- nn$BCEWithLogitsLoss()
  for (ep in seq_len(hp$head_epochs)) {
    opt_h$zero_grad()
    logit <- head(emb_py)$view(-1L)
    loss_h <- bce(logit, yt)
    loss_h$backward()
    opt_h$step()
    if (ep %% 10L == 0L) gc(FALSE)   # release per-epoch wrappers (computation-neutral; see SSL loop)
  }
  head$eval()
  head_w <- as.numeric(reticulate::py_to_r(
    head$weight$detach()$cpu()$to(torch$float64)$numpy()))   # length d
  head_b <- as.numeric(reticulate::py_to_r(
    head$bias$detach()$cpu()$to(torch$float64)$numpy()))     # scalar

  list(weights = weights, head_w = head_w, head_b = head_b, embedding_dim = d)
}


# ----------------------------------------------------------------------------
# fit_ai_scarf
# ----------------------------------------------------------------------------

#' @title Fit the SCARF self-supervised contrastive single-sample discriminator
#'
#' @description
#' Trains a BatchNorm-free MLP encoder on the per-sample-rCLR, universe-aligned
#' TRAINING matrix by SELF-SUPERVISED CONTRASTIVE training (SCARF InfoNCE /
#' NT-Xent, NO labels; via reticulate-python torch, Adam), FREEZES the encoder,
#' fits a FROZEN linear-probe head on the labeled training embeddings, EXPORTS the
#' encoder weights and the head \eqn{w, b} to R as plain numeric matrices/vectors,
#' and DISCARDS the projection head and the python module. The head-training
#' embeddings are computed with the SAME pure-R forward the scorer uses
#' (guaranteeing fit/score consistency despite float32 training). The fitted model
#' holds no external pointer and the corruption marginal is NOT carried into it.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control. Used only to fit
#'   the linear-probe head; the encoder is trained WITHOUT labels.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{hidden}
#'   (positive-integer vector of encoder layer widths up to the embedding; the LAST
#'   entry is the embedding dim \eqn{d}; default \code{c(64L, 32L)}),
#'   \code{activation} (\code{"relu"} (default) or \code{"tanh"}), \code{epochs}
#'   (number of SSL contrastive steps, positive integer; default \code{200L}),
#'   \code{lr} (positive Adam learning rate; default \code{1e-3}),
#'   \code{weight_decay} (non-negative Adam L2 for the SSL stage; default
#'   \code{1e-4}), \code{corrupt_p} (per-feature corruption probability, in
#'   \eqn{(0, 1)}; default \code{0.6}, the SCARF paper value), \code{temperature}
#'   (InfoNCE \eqn{\tau}, positive; default \code{1.0}), \code{proj_dim}
#'   (projection-head width, positive integer; default \code{64L}),
#'   \code{head_epochs} (linear-probe fixed epochs, positive integer; default
#'   \code{200L}), \code{head_weight_decay} (non-negative Adam L2 for the head;
#'   default \code{1e-2}), \code{min_features} (feature-overlap floor at scoring,
#'   positive integer; default \code{3L}), \code{device} (\code{"cpu"} (default),
#'   \code{"cuda"}, or \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU
#'   with no GPU), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{ai_scarf_model}: a list with
#'   \code{feature_universe}, \code{weights} (exported per-layer \code{list(W, b)}
#'   for the encoder up to the embedding), \code{activation}, \code{head_w}
#'   (length-\eqn{d} frozen head weight), \code{head_b} (scalar frozen head bias),
#'   \code{embedding_dim}, \code{device} (resolved), \code{seed}, and \code{hp}.
#'   The corruption marginal is NOT stored (it is fit-only).
#'
#' @details
#' The SSL stage learns the encoder by the SCARF InfoNCE loss: each step builds a
#' random-feature-corrupted view of every training row (corrupted entries resampled
#' from each feature's empirical training marginal), embeds the clean row and its
#' corrupted view through encoder + projection head, and applies the NT-Xent loss
#' (cosine similarity, temperature \eqn{\tau}) with the positive pair being the row
#' and its corrupted view and the negatives being the other batch items. The
#' projection head is discarded; only the encoder is exported. The head stage fits
#' a single \code{Linear(d -> 1)} probe by \code{head_weight_decay}-regularised BCE
#' on the frozen encoder embeddings. The score is the linear-probe logit
#' \eqn{s(z) = w \cdot z + b}.
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
#' model <- fit_ai_scarf(X, y)
#' score_ai_scarf(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Bahri D, Jiang H, Tay Y, Metzler D. (2022) SCARF: Self-Supervised Contrastive
#' Learning using Random Feature Corruption. \emph{ICLR}. arXiv:2106.15147.
#'
#' @export
fit_ai_scarf <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ai_scarf", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ai_scarf", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ai_scarf")
  hp <- .ai_scarf_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_ai_scarf: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's
  # first python initialisation seeds R's RNG as a side effect, so the no-prior-
  # seed state must be captured here (not after the import) to be restored exactly
  # (including removing a leaked .Random.seed when none pre-existed). Fit draws no
  # R randomness of its own (corruption randomness is torch-side); this purely
  # neutralises that init side effect.
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

  .ai_scarf_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .ai_scarf_rclr_matrix(X_train)            # n x p rCLR over full universe

  device <- .ai_scarf_resolve_device(hp$device)
  trained <- .ai_scarf_train_export(Z_in, y, hp, device, .ai_scarf_forward)
  activation <- hp$activation

  model <- list(
    feature_universe = feat_order,
    weights = trained$weights,
    activation = activation,
    head_w = as.numeric(trained$head_w),
    head_b = as.numeric(trained$head_b),
    embedding_dim = trained$embedding_dim,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "ai_scarf_model"
  model
}


# Linear-head logit score for one embedded specimen z (length d):
#   w . z + b.  Larger = more case-like.
.ai_scarf_score_one <- function(z, head_w, head_b) {
  sum(head_w * z) + head_b
}


# ----------------------------------------------------------------------------
# score_ai_scarf
# ----------------------------------------------------------------------------

#' @title Score the SCARF self-supervised contrastive single-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_ai_scarf}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), embedded by the exported encoder forward
#' (a CLEAN forward -- NO corruption; the corruption marginal is not consulted at
#' score), and scored by the frozen linear-probe logit \eqn{w \cdot z + b}. Larger
#' = more case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on
#' that row and the frozen model and is exactly invariant to per-specimen positive
#' scaling.
#'
#' @param model An \code{ai_scarf_model} object from \code{\link{fit_ai_scarf}}.
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
#' model <- fit_ai_scarf(X, y)
#' score_ai_scarf(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Bahri D, Jiang H, Tay Y, Metzler D. (2022) SCARF: Self-Supervised Contrastive
#' Learning using Random Feature Corruption. \emph{ICLR}. arXiv:2106.15147.
#'
#' @export
score_ai_scarf <- function(model, X, meta = NULL) {
  if (!inherits(model, "ai_scarf_model")) {
    stop("score_ai_scarf: model must have class ai_scarf_model")
  }
  X <- .reo_check_matrix(X, "score_ai_scarf", "X")
  .reo_check_meta(meta, nrow(X), "score_ai_scarf", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .ai_scarf_align(X, feat_order)
  weights <- model$weights
  activation <- model$activation
  head_w <- model$head_w
  head_b <- model$head_b

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (rCLR
    # origin) but is a valid specimen that MUST be embedded and scored, not floored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .ai_scarf_rclr(X_use[i, ])
    emb <- .ai_scarf_forward(matrix(z, nrow = 1L), weights, activation)
    out[i] <- .ai_scarf_score_one(as.numeric(emb), head_w, head_b)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ai_scarf: scorer produced non-finite or wrong-length output")
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
.ai_scarf_model_digest <- function(model) {
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
