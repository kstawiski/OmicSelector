#' @title Prototypical-network single-sample discriminator (frozen class prototypes)
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing a
#' prototypical network in the sense of Snell, Swersky & Zemel (2017), adapted to
#' compositional specimens. A small BatchNorm-free MLP embedding of the
#' per-sample-rCLR specimen is learned on the TRAINING set by EPISODIC
#' PROTOTYPICAL training (NOT cross-entropy) via reticulate-python torch, then
#' \strong{FROZEN}; per-class PROTOTYPES (class-mean embeddings) are computed in
#' the embedding space, and a new specimen is scored by the difference of its
#' squared-Euclidean distances to the frozen control and case prototypes. Larger
#' score = more case-like.
#'
#' \strong{Episodic prototypical training (Snell 2017).} Unlike the
#' cross-entropy-trained \code{lrt-deepmaha} sibling, the embedding here is trained
#' with the prototypical loss. Each episode draws a disjoint SUPPORT set
#' (\code{n_support} per class) and QUERY set (\code{n_query} per class) from the
#' (rCLR'd, universe-aligned) training rows; the prototype of class \eqn{c} is the
#' mean SUPPORT embedding of that class; each query is assigned
#' \eqn{\log\,\mathrm{softmax}(-d_c)} over its squared-Euclidean distances
#' \eqn{d_c = \lVert z_q - \mathrm{proto}_c\rVert^2} to the prototypes, and the loss
#' is the mean negative log-probability of the query's TRUE class. This is exactly
#' Snell 2017 eq. (1)-(2) with the squared-Euclidean distance (the Bregman-divergence
#' choice the paper shows is principled). Adam optimises the embedding net over
#' \code{epochs} episodes.
#'
#' \strong{Torch is confined to FIT.} The MLP is trained with torch on the rCLR'd,
#' universe-aligned training matrix. After training, the weight matrices and biases
#' up to and INCLUDING the embedding layer are EXPORTED to R as plain numeric
#' matrices/vectors; the python module is DISCARDED. The fitted model therefore holds
#' NO external pointer and NO python dependency at score time: the embedding forward,
#' the frozen prototypes, and the Euclidean-to-prototype score are all pure base-R.
#' This makes the default \code{model_digest} a stable deterministic snapshot and lets
#' the package suite run the scorer on any host without the venv.
#'
#' \strong{Fit/score consistency despite float32.} Torch trains in float32; an R
#' float64 forward on the exported (float32-valued) weights differs from the torch
#' embedding by \eqn{\sim 10^{-6}}. To guarantee that the frozen prototypes live in
#' EXACTLY the same embedding the scorer evaluates, the prototypes
#' \eqn{\mu_{\mathrm{case}}, \mu_{\mathrm{ctrl}}} are computed from the PURE-R forward
#' over ALL training rows (\code{.proto_net_forward}), NOT from the torch (or the
#' episodic support-mean) embeddings. Fit and score thus use the identical
#' computation (the deepmaha / nc-ecod-copod / nc-sinkhorn-single lesson: fit and
#' score must share one code path).
#'
#' \strong{The embedding forward.} The exported net is
#' \eqn{\mathrm{input}(p) \to [\mathrm{Linear} \to \phi]_{1}
#'   \to \dots \to [\mathrm{Linear} \to \phi]_{L-1}
#'   \to \mathrm{Linear} \to \mathrm{embedding}(d)}. There are \code{length(hidden)}
#' linear layers up to the embedding (the last \code{hidden} width is the embedding
#' dim \eqn{d}); the activation (\code{"relu"} default, or \code{"tanh"}) is applied
#' AFTER every linear layer EXCEPT the final embedding layer (no activation follows
#' the embedding -- it is the raw representation, matching the torch net).
#' \code{.proto_net_forward(M, weights)} reimplements this exactly.
#'
#' \strong{Euclidean-to-prototype score.} With the frozen prototypes the score of an
#' embedded specimen \eqn{z} is
#' \deqn{s(z) = -\lVert z - \mu_{\mathrm{case}}\rVert^2
#'   + \lVert z - \mu_{\mathrm{ctrl}}\rVert^2,}
#' the squared-Euclidean log-likelihood ratio to the prototypes. This is exactly the
#' deepmaha score with \eqn{\Sigma = I} (no covariance / no inversion) -- proto-net is
#' the identity-metric, prototype-mean special case, but with an embedding trained by
#' the PROTOTYPICAL loss, which is what makes it a distinct method
#' ("Euclidean-to-prototype metric learning"). Larger = more case-like.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly
#' invariant to per-specimen scaling and uses no cross-row statistic. Universe
#' features absent from a specimen carry the neutral rCLR value \code{0}. The SAME
#' rCLR transform is applied to the frozen training rows at fit and to every query at
#' score, so a specimen's score depends only on that specimen and the frozen model --
#' it passes \code{\link{singlesample_assert_row_equivariant}}.
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
#' Snell J, Swersky K, Zemel R. (2017) Prototypical Networks for Few-shot Learning.
#' \emph{Advances in Neural Information Processing Systems (NeurIPS)} 30.
#' arXiv:1703.05175.
#'
#' @name singlesample-protonet
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .lrt_deepmaha_rclr.
.proto_net_rclr <- function(v) {
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
.proto_net_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .proto_net_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.proto_net_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R embedding forward (exact reimplementation of the torch net up to d)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported MLP up to (and including) the embedding
# layer. `weights` is a list of per-layer lists list(W = out x in, b = out):
# layers 1..(L-1) are hidden (Linear then `activation`), layer L is the embedding
# Linear with NO activation (matches the torch net exactly). For a layer with
# weight W (out x in) and bias b, torch computes M %*% t(W) + b (row = sample).
# Returns an n x d embedding.
.proto_net_forward <- function(M, weights, activation) {
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

.proto_net_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_proto_net: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_proto_net: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_proto_net: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("hidden", "activation", "epochs", "lr", "weight_decay",
               "n_support", "n_query", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_proto_net: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  hidden <- hp$hidden
  if (is.null(hidden)) hidden <- c(64L, 32L)
  if (!is.numeric(hidden) || length(hidden) < 1L || any(!is.finite(hidden)) ||
      any(hidden < 1L) || any(hidden > .Machine$integer.max) ||
      any(hidden != as.integer(hidden))) {
    stop("fit_proto_net: hp$hidden must be a vector of positive integers")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_proto_net: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_proto_net: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_proto_net: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_proto_net: hp$weight_decay must be a non-negative finite number")
  }

  n_support <- hp$n_support
  if (is.null(n_support)) n_support <- 5L
  if (!is.numeric(n_support) || length(n_support) != 1L ||
      !is.finite(n_support) || n_support < 1L ||
      n_support > .Machine$integer.max ||
      n_support != as.integer(n_support)) {
    stop("fit_proto_net: hp$n_support must be a positive integer")
  }

  n_query <- hp$n_query
  if (is.null(n_query)) n_query <- 5L
  if (!is.numeric(n_query) || length(n_query) != 1L ||
      !is.finite(n_query) || n_query < 1L ||
      n_query > .Machine$integer.max ||
      n_query != as.integer(n_query)) {
    stop("fit_proto_net: hp$n_query must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_proto_net: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_proto_net: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_proto_net: hp$seed must be a single integer")
  }

  list(
    hidden = as.integer(hidden),
    activation = activation,
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    n_support = as.integer(n_support),
    n_query = as.integer(n_query),
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
.proto_net_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_proto_net: package 'reticulate' is required to train the ",
         "embedding. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_proto_net: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.proto_net_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# Episode index plan (deterministic, R-side seeded)
# ----------------------------------------------------------------------------

# Pre-draw, in PURE R under the fit's .Random.seed guard, the per-episode support
# and query row indices for each class. Doing the episode sampling in R (seeded
# with set.seed(hp$seed)) keeps the draw deterministic and reproducible: two
# seed-matched CPU fits replay the identical episode plan, and the global
# .Random.seed is saved/restored by the caller so fit stays RNG-neutral. The
# plan is a list of `epochs` episodes; each episode is a list(support = int vec,
# query = int vec, query_label = int vec of 0/1 aligned to `query`).
#
# Feasibility: per class we need >= 1 support and >= 1 query, and support/query
# are DISJOINT within the episode. n_support / n_query are clamped per class to
# leave at least 1 row for the query (so a class with m rows uses up to
# min(n_support, m - 1) support and up to min(n_query, m - support) query, with a
# minimum of 1 support and 1 query when m >= 2). A class with a single row
# contributes that row to BOTH support and query for that episode (degenerate but
# valid: prototype == the query, distance 0).
.proto_net_episode_plan <- function(idx_case, idx_ctrl, hp) {
  classes <- list(`0` = idx_ctrl, `1` = idx_case)
  n_sup_req <- hp$n_support
  n_qry_req <- hp$n_query
  plan <- vector("list", hp$epochs)
  for (ep in seq_len(hp$epochs)) {
    sup <- integer(0L); qry <- integer(0L); qry_lab <- integer(0L)
    for (cl in c(0L, 1L)) {
      pool <- classes[[as.character(cl)]]
      m <- length(pool)
      if (m >= 2L) {
        n_sup <- min(n_sup_req, m - 1L)              # leave >= 1 for query
        n_sup <- max(1L, n_sup)
        sup_i <- pool[sample.int(m, n_sup)]
        rest <- setdiff(pool, sup_i)
        n_qry <- min(n_qry_req, length(rest))
        n_qry <- max(1L, n_qry)
        qry_i <- rest[sample.int(length(rest), n_qry)]
      } else {
        # single-row class: same row as support and query (proto == query).
        sup_i <- pool
        qry_i <- pool
      }
      sup <- c(sup, sup_i)
      qry <- c(qry, qry_i)
      qry_lab <- c(qry_lab, rep.int(cl, length(qry_i)))
    }
    plan[[ep]] <- list(support = sup, query = qry, query_label = qry_lab)
  }
  plan
}


# ----------------------------------------------------------------------------
# Episodic prototypical training + weight export
# ----------------------------------------------------------------------------

# Train the BN-free MLP embedding by EPISODIC PROTOTYPICAL training (Snell 2017)
# via torch and return the exported R-side weights (a list of per-layer
# list(W, b), up to and INCLUDING the embedding layer). The python module is used
# here only; nothing python is returned. The torch RNG (CPU + guarded CUDA) is
# saved/restored, so fit leaves the torch generators byte-unchanged; the global R
# .Random.seed is saved/restored by the caller. `Z` is the n x p rCLR'd,
# universe-aligned training matrix; `y` the integer 0/1 labels; `plan` the
# pre-drawn deterministic episode index plan (1-based R indices).
.proto_net_train_export <- function(Z, y, plan, hp, device) {
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

  # Seed the torch generators (weight init + any torch-side randomness) so two
  # seed-matched CPU fits initialise identically. Episode sampling is R-side and
  # already deterministic via `plan`.
  torch$manual_seed(as.integer(hp$seed))

  np <- reticulate::import("numpy", delay_load = FALSE)
  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)

  p <- ncol(Z)
  hidden <- hp$hidden
  # Build input(p) -> [Linear -> act]xlen(hidden) but DROP the activation after the
  # final (embedding) layer, so the net IS the embedding net (no classification
  # head -- prototypical training needs only the embedding).
  act_mod <- if (identical(hp$activation, "tanh")) nn$Tanh else nn$ReLU
  layers <- list()
  in_dim <- p
  for (h in hidden) {
    layers[[length(layers) + 1L]] <- nn$Linear(in_dim, as.integer(h))
    layers[[length(layers) + 1L]] <- act_mod()
    in_dim <- as.integer(h)
  }
  layers[[length(layers)]] <- NULL        # drop activation after the embedding
  d <- in_dim
  net <- do.call(nn$Sequential, layers)$to(device)

  opt <- torch$optim$Adam(net$parameters(), lr = hp$lr,
                          weight_decay = hp$weight_decay)

  net$train()
  for (ep in seq_along(plan)) {
    epi <- plan[[ep]]
    # 0-based torch indices for support + query rows of this episode.
    sup_idx <- torch$as_tensor(np$asarray(as.integer(epi$support - 1L),
                                          dtype = "int64"))$long()$to(device)
    qry_idx <- torch$as_tensor(np$asarray(as.integer(epi$query - 1L),
                                          dtype = "int64"))$long()$to(device)
    sup_emb <- net(X_py$index_select(0L, sup_idx))      # n_sup x d
    qry_emb <- net(X_py$index_select(0L, qry_idx))      # n_qry x d

    # Prototype per class = mean support embedding of that class.
    sl <- .proto_net_support_label(epi, y)
    sup_ctrl_emb <- sup_emb$index_select(
      0L, torch$as_tensor(np$asarray(as.integer(which(sl == 0L) - 1L),
                                     dtype = "int64"))$long()$to(device))
    sup_case_emb <- sup_emb$index_select(
      0L, torch$as_tensor(np$asarray(as.integer(which(sl == 1L) - 1L),
                                     dtype = "int64"))$long()$to(device))
    proto_ctrl <- sup_ctrl_emb$mean(0L, keepdim = TRUE)     # 1 x d
    proto_case <- sup_case_emb$mean(0L, keepdim = TRUE)     # 1 x d
    protos <- torch$cat(list(proto_ctrl, proto_case), 0L)   # 2 x d (ctrl,case)

    # Squared-Euclidean distances query x 2: ||q - proto_c||^2.
    # cdist gives Euclidean; square it for the Bregman (squared) distance.
    dists <- torch$cdist(qry_emb, protos)$pow(2)            # n_qry x 2
    logp <- nn$functional$log_softmax(dists$mul(-1), dim = 1L)  # -d -> log_softmax
    # True-class column index (0 = ctrl, 1 = case) per query.
    tgt <- torch$as_tensor(np$asarray(as.integer(epi$query_label),
                                      dtype = "int64"))$long()$to(device)
    loss <- nn$functional$nll_loss(logp, tgt)               # mean -log p_true

    opt$zero_grad()
    loss$backward()
    opt$step()
  }
  net$eval()

  # ---- export weights up to and INCLUDING the embedding layer ---------------
  modules <- reticulate::iterate(net$modules())
  weights <- list()
  for (m in modules) {
    has_w <- tryCatch(!is.null(m$weight), error = function(e) FALSE)
    if (!isTRUE(has_w)) next
    W <- reticulate::py_to_r(m$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(m$bias$detach()$cpu()$to(torch$float64)$numpy())
    weights[[length(weights) + 1L]] <- list(W = W, b = as.numeric(b))
  }
  # No classification head was added: every Linear up to the embedding is kept.
  if (length(weights) < 1L) {
    stop("fit_proto_net: exported net has no linear layers", call. = FALSE)
  }
  weights
}

# Map an episode's `support` row indices (into the full training matrix) to their
# 0/1 class labels, in the SAME order they appear in `epi$support`. `y` is the
# integer 0/1 label vector over all training rows.
.proto_net_support_label <- function(epi, y) {
  as.integer(y[epi$support])
}


# ----------------------------------------------------------------------------
# fit_proto_net
# ----------------------------------------------------------------------------

#' @title Fit the prototypical-network single-sample discriminator
#'
#' @description
#' Trains a BatchNorm-free MLP embedding on the per-sample-rCLR, universe-aligned
#' TRAINING matrix by EPISODIC PROTOTYPICAL training (Snell 2017; via
#' reticulate-python torch, Adam, prototypical loss -- NOT cross-entropy), EXPORTS
#' the weights up to and including the embedding layer to R as plain numeric
#' matrices/vectors, DISCARDS the python module, then computes the frozen per-class
#' PROTOTYPES (class-mean embeddings) using the PURE-R forward over ALL training rows
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
#'   (\code{"relu"} (default) or \code{"tanh"}), \code{epochs} (number of training
#'   episodes, positive integer; default \code{200L}), \code{lr} (positive Adam
#'   learning rate; default \code{1e-3}), \code{weight_decay} (non-negative Adam L2;
#'   default \code{1e-4}), \code{n_support} (per-class support size per episode,
#'   positive integer; default \code{5L}), \code{n_query} (per-class query size per
#'   episode, positive integer; default \code{5L}), \code{min_features}
#'   (feature-overlap floor at scoring, positive integer; default \code{3L}),
#'   \code{device} (\code{"cpu"} (default), \code{"cuda"}, or \code{"auto"};
#'   \code{"cuda"}/\code{"auto"} fall back to CPU with no GPU), and \code{seed}
#'   (integer; default \code{42L}). \code{n_support}/\code{n_query} are clamped per
#'   class so each episode is feasible (disjoint support/query, \eqn{\geq 1} of each
#'   per class).
#'
#' @return Object of class \code{proto_net_model}: a list with
#'   \code{feature_universe}, \code{weights} (exported per-layer \code{list(W, b)}
#'   up to the embedding), \code{activation}, \code{proto_case}, \code{proto_ctrl}
#'   (length-\eqn{d} frozen prototypes), \code{embedding_dim}, \code{device}
#'   (resolved), \code{seed}, and \code{hp}.
#'
#' @details
#' The score is the squared-Euclidean log-likelihood ratio to the prototypes,
#' \eqn{s(z) = -\lVert z - \mu_{\mathrm{case}}\rVert^2
#' + \lVert z - \mu_{\mathrm{ctrl}}\rVert^2}, exactly the deepmaha score with
#' \eqn{\Sigma = I}. The distinction is the embedding: it is trained by the
#' prototypical loss (CE over the softmax of negative squared-Euclidean distances to
#' episodic class-mean prototypes), not a cross-entropy classification head.
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
#' model <- fit_proto_net(X, y)
#' score_proto_net(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Snell J, Swersky K, Zemel R. (2017) Prototypical Networks for Few-shot Learning.
#' \emph{NeurIPS} 30. arXiv:1703.05175.
#'
#' @export
fit_proto_net <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_proto_net", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_proto_net", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_proto_net")
  hp <- .proto_net_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_proto_net: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import and BEFORE the
  # R-side episode sampling: reticulate's first python init seeds R's RNG, and the
  # episode plan draws R randomness (set.seed below). Restoring on exit (incl. the
  # no-prior-seed case) keeps fit globally RNG-neutral and seed-deterministic.
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

  .proto_net_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .proto_net_rclr_matrix(X_train)           # n x p rCLR over full universe

  # Deterministic episode plan (R-side, seeded). set.seed here is local: the outer
  # on.exit restores the global .Random.seed, so this leaves no global RNG trace.
  set.seed(as.integer(hp$seed))
  idx_case <- which(y == 1L)
  idx_ctrl <- which(y == 0L)
  plan <- .proto_net_episode_plan(idx_case, idx_ctrl, hp)

  device <- .proto_net_resolve_device(hp$device)
  weights <- .proto_net_train_export(Z_in, y, plan, hp, device)
  activation <- hp$activation

  # Frozen DEPLOYMENT prototypes from the PURE-R forward over ALL training rows
  # (NOT the episodic support means, NOT the torch embeddings) -- fit/score share
  # one embedding code path (float32 -> float64 safe).
  Z_emb <- .proto_net_forward(Z_in, weights, activation)   # n x d
  d <- ncol(Z_emb)
  proto_case <- colMeans(Z_emb[idx_case, , drop = FALSE])
  proto_ctrl <- colMeans(Z_emb[idx_ctrl, , drop = FALSE])

  model <- list(
    feature_universe = feat_order,
    weights = weights,
    activation = activation,
    proto_case = as.numeric(proto_case),
    proto_ctrl = as.numeric(proto_ctrl),
    embedding_dim = d,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "proto_net_model"
  model
}


# Euclidean-to-prototype score for one embedded specimen z (length d):
#   -||z - proto_case||^2 + ||z - proto_ctrl||^2.  Larger = more case-like.
# This is the deepmaha score with Sigma = I (no covariance / no inverse).
.proto_net_score_one <- function(z, proto_case, proto_ctrl) {
  -sum((z - proto_case)^2) + sum((z - proto_ctrl)^2)
}


# ----------------------------------------------------------------------------
# score_proto_net
# ----------------------------------------------------------------------------

#' @title Score the prototypical-network single-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_proto_net}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), embedded by the exported MLP forward, and
#' scored by the difference of squared-Euclidean distances to the frozen control and
#' case prototypes. Larger = more case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on that
#' row and the frozen model and is exactly invariant to per-specimen positive scaling.
#'
#' @param model A \code{proto_net_model} object from \code{\link{fit_proto_net}}.
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
#' model <- fit_proto_net(X, y)
#' score_proto_net(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Snell J, Swersky K, Zemel R. (2017) Prototypical Networks for Few-shot Learning.
#' \emph{NeurIPS} 30. arXiv:1703.05175.
#'
#' @export
score_proto_net <- function(model, X, meta = NULL) {
  if (!inherits(model, "proto_net_model")) {
    stop("score_proto_net: model must have class proto_net_model")
  }
  X <- .reo_check_matrix(X, "score_proto_net", "X")
  .reo_check_meta(meta, nrow(X), "score_proto_net", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .proto_net_align(X, feat_order)
  weights <- model$weights
  activation <- model$activation
  proto_case <- model$proto_case
  proto_ctrl <- model$proto_ctrl

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (rCLR
    # origin) but is a valid specimen that MUST be embedded and scored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .proto_net_rclr(X_use[i, ])
    emb <- .proto_net_forward(matrix(z, nrow = 1L), weights, activation)
    out[i] <- .proto_net_score_one(as.numeric(emb), proto_case, proto_ctrl)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_proto_net: scorer produced non-finite or wrong-length output")
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
.proto_net_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    weights = model$weights,
    activation = model$activation,
    proto_case = model$proto_case,
    proto_ctrl = model$proto_ctrl,
    embedding_dim = model$embedding_dim,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
