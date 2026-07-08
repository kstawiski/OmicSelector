#' @title IB-IRM (Information-Bottleneck Invariant Risk Minimization) cross-cohort transfer discriminator
#'
#' @description
#' Single-sample TRANSFER-estimand discriminator (family N) implementing
#' \strong{IB-IRM} (Ahuja, Caballero, Zhang et al. 2021 NeurIPS;
#' arXiv:2106.06607) for cross-cohort out-of-distribution generalization. A small
#' BatchNorm-free MLP encoder of the per-sample-rCLR specimen and a linear
#' classification head \eqn{\mathrm{Linear}(d \to 1)} are trained JOINTLY across
#' multiple training ENVIRONMENTS (cohorts) by the IB-IRM objective via
#' reticulate-python torch, then \strong{FROZEN}; a new specimen is scored by its
#' frozen logit on its pure-R encoder embedding. Larger score = more case-like.
#'
#' \strong{Transfer is in the FIT, not the score.} The cross-cohort aspect of this
#' method lives ENTIRELY in training: the encoder + head are optimised so the
#' per-environment risks are low (mean term), the per-environment classifier is
#' SIMULTANEOUSLY OPTIMAL across cohorts (the IRMv1 gradient-penalty term), and the
#' representation is COMPRESSED (the Information-Bottleneck embedding-variance term),
#' which is IB-IRM's mechanism for distributional robustness across cohorts. The
#' SCORE path is identical to the within-estimand trained-torch siblings
#' (\code{lrt-deepmaha}, \code{proto-net}, \code{ssl-vicreg}, \code{dro-vrex},
#' \code{dg-fishr}): a frozen linear logit on a pure-R encoder embedding of the single
#' specimen, so single-row scoring equals batch scoring EXACTLY and the score is
#' exactly invariant to per-specimen positive scaling.
#'
#' \strong{The IB-IRM objective (Ahuja 2021).} Let \eqn{z = \mathrm{encoder}(x)} be
#' the embedding and the head logit be \eqn{l = w_0^\top z + b_0}. IB-IRM minimises
#' \deqn{\mathcal{L} = \frac{1}{|\mathcal{E}|}\sum_e R_e
#'        \;+\; \lambda_{\mathrm{IRM}} \cdot \Omega_{\mathrm{IRM}}
#'        \;+\; \lambda_{\mathrm{IB}} \cdot \Omega_{\mathrm{IB}},}
#' where \eqn{R_e} is the mean binary-cross-entropy-with-logits risk on environment
#' \eqn{e}. The two penalties are:
#' \enumerate{
#'   \item \strong{IRMv1 gradient penalty} (Arjovsky 2019, the IRM term IB-IRM
#'     inherits). A FIXED dummy scalar classifier multiplier \eqn{w = 1.0} is applied
#'     to the head logit (\eqn{l_e = w \cdot (w_0^\top z_e + b_0)}); for each VALID
#'     environment \eqn{e}, \eqn{R_e(w)} is the mean BCE-with-logits at \eqn{w}, and
#'     \deqn{\Omega_{\mathrm{IRM}} = \frac{1}{|\mathcal{E}|}\sum_e
#'           \left( \left. \frac{\partial R_e}{\partial w} \right|_{w=1} \right)^2.}
#'     The gradient \eqn{\partial R_e / \partial w} is taken by torch
#'     \code{autograd.grad(..., create_graph = TRUE)} so the penalty itself stays
#'     differentiable w.r.t. the encoder + head. Because \eqn{w} is a scalar,
#'     \eqn{\lVert \partial R_e / \partial w \rVert^2} is just the SQUARE of that
#'     scalar gradient. This is the canonical IRMv1 penalty.
#'   \item \strong{Information-Bottleneck (IB) variance penalty} (the term that
#'     distinguishes IB-IRM from plain IRM). The representation is compressed by
#'     penalising the per-dimension variance of the embedding over the POOLED retained
#'     rows:
#'     \deqn{\Omega_{\mathrm{IB}} = \frac{1}{d}\sum_{j=1}^{d}
#'           \mathrm{Var}_i(z_{ij}) = \frac{1}{d}\sum_{j=1}^{d}
#'           \frac{1}{n}\sum_{i=1}^{n}\left(z_{ij} - \bar z_{\cdot j}\right)^2,}
#'     the mean over embedding dimensions of the BIASED (population, \eqn{1/n})
#'     variance of each embedding feature. Minimising representation variance is the
#'     standard tractable IB surrogate (lower feature variance \eqn{\to} lower
#'     representation entropy \eqn{\to} a more compressed, invariant code). The biased
#'     \eqn{1/n} variance is computed EXPLICITLY (not \code{torch$var}, whose default
#'     is the \eqn{N-1} unbiased estimator).
#' }
#' Both \eqn{\lambda_{\mathrm{IRM}}} (\code{irm_lambda}, default \code{1.0}) and
#' \eqn{\lambda_{\mathrm{IB}}} (\code{ib_lambda}, default \code{1.0}) are non-negative.
#' Encoder + head are optimised jointly with Adam over \code{epochs} full-batch steps.
#'
#' \strong{Environments from \code{meta_train}.} Environments are the distinct
#' non-missing labels in \code{meta_train[[cohort_col]]}. \code{cohort_col} defaults
#' to \code{NULL} and is auto-detected from common names
#' (\code{"accession"}, \code{"cohort"}, \code{"batch"}, \code{"study"},
#' \code{"dataset"}, \code{"group"}); set it explicitly to override. An environment
#' is VALID for the IRM penalty only if it carries at least one case AND one
#' control (a single-class environment has a degenerate, label-uninformative risk
#' and is DROPPED from the per-environment gradient set). If fewer than two valid
#' environments remain (including \code{NULL}/absent \code{meta_train}, a missing
#' \code{cohort_col}, or a single distinct cohort), fitting FALLS BACK to plain ERM
#' (a single pooled risk over all rows, BOTH penalties inert) -- mirroring
#' \code{dro-vrex}'s and \code{dg-fishr}'s NULL-meta single-environment ERM
#' degeneration. The realised number of environments is recorded as
#' \code{model$n_environments}; \code{n_environments == 1} means the fit ran as ERM and
#' the transfer estimand was NOT engaged.
#'
#' \strong{Torch is confined to FIT.} After training, the encoder weight matrices and
#' biases up to and INCLUDING the embedding layer, and the head \eqn{(w, b)}, are
#' EXPORTED to R as plain numeric matrices/vectors; the python module and the torch
#' objects are DISCARDED. The fitted model holds NO external pointer and NO python
#' dependency at score time: the encoder forward and the linear-head logit are pure
#' base-R. This makes the default \code{model_digest} a stable deterministic snapshot
#' and lets the package suite run the scorer on any host without the venv.
#'
#' \strong{The embedding forward.} The exported encoder is
#' \eqn{\mathrm{input}(p) \to [\mathrm{Linear} \to \phi]_{1} \to \dots \to
#' [\mathrm{Linear} \to \phi]_{L-1} \to \mathrm{Linear} \to \mathrm{embedding}(d)}.
#' There are \code{length(hidden)} linear layers up to the embedding (the last
#' \code{hidden} width is the embedding dim \eqn{d}); the activation (\code{"relu"}
#' default, or \code{"tanh"}) is applied AFTER every linear layer EXCEPT the final
#' embedding layer. \code{.dg_ibirm_forward(M, weights, activation)} reimplements this
#' exactly. The score of an embedded specimen \eqn{z} is the frozen linear logit
#' \eqn{s(z) = w^\top z + b}; larger = more case-like.
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
#' Ahuja K, Caballero E, Zhang D, Gagnon-Audet J-C, Bengio Y, Mitliagkas I,
#' Rish I. (2021) Invariance Principle Meets Information Bottleneck for
#' Out-of-Distribution Generalization. \emph{Advances in Neural Information
#' Processing Systems (NeurIPS) 34}. arXiv:2106.06607.
#'
#' Arjovsky M, Bottou L, Gulrajani I, Lopez-Paz D. (2019) Invariant Risk
#' Minimization. arXiv:1907.02893.
#'
#' @name singlesample-dg-ibirm
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .dg_fishr_rclr.
.dg_ibirm_rclr <- function(v) {
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
.dg_ibirm_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .dg_ibirm_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.dg_ibirm_align <- function(X, feat_order) {
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
.dg_ibirm_forward <- function(M, weights, activation) {
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

.dg_ibirm_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_dg_ibirm: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_dg_ibirm: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_dg_ibirm: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("hidden", "activation", "epochs", "lr", "weight_decay",
               "irm_lambda", "ib_lambda", "cohort_col", "min_features",
               "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_dg_ibirm: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  hidden <- hp[["hidden"]]
  if (is.null(hidden)) hidden <- c(64L, 32L)
  if (!is.numeric(hidden) || length(hidden) < 1L || any(!is.finite(hidden)) ||
      any(hidden < 1L) || any(hidden > .Machine$integer.max) ||
      any(hidden != as.integer(hidden))) {
    stop("fit_dg_ibirm: hp$hidden must be a vector of positive integers")
  }

  activation <- hp[["activation"]]
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_dg_ibirm: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp[["epochs"]]
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_dg_ibirm: hp$epochs must be a positive integer")
  }

  lr <- hp[["lr"]]
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_dg_ibirm: hp$lr must be a positive finite number")
  }

  weight_decay <- hp[["weight_decay"]]
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_dg_ibirm: hp$weight_decay must be a non-negative finite number")
  }

  irm_lambda <- hp[["irm_lambda"]]
  if (is.null(irm_lambda)) irm_lambda <- 1.0
  if (!is.numeric(irm_lambda) || length(irm_lambda) != 1L ||
      !is.finite(irm_lambda) || irm_lambda < 0) {
    stop("fit_dg_ibirm: hp$irm_lambda must be a non-negative finite number")
  }

  ib_lambda <- hp[["ib_lambda"]]
  if (is.null(ib_lambda)) ib_lambda <- 1.0
  if (!is.numeric(ib_lambda) || length(ib_lambda) != 1L ||
      !is.finite(ib_lambda) || ib_lambda < 0) {
    stop("fit_dg_ibirm: hp$ib_lambda must be a non-negative finite number")
  }

  # cohort_col: NULL (auto-detect at fit) or a single non-empty character string.
  cohort_col <- hp[["cohort_col"]]
  if (!is.null(cohort_col)) {
    if (!is.character(cohort_col) || length(cohort_col) != 1L ||
        is.na(cohort_col) || !nzchar(cohort_col)) {
      stop("fit_dg_ibirm: hp$cohort_col must be NULL or a single non-empty ",
           "character string")
    }
  }

  min_features <- hp[["min_features"]]
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_dg_ibirm: hp$min_features must be a positive integer")
  }

  device <- hp[["device"]]
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_dg_ibirm: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_dg_ibirm: hp$seed must be a single integer")
  }

  list(
    hidden = as.integer(hidden),
    activation = activation,
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    irm_lambda = as.numeric(irm_lambda),
    ib_lambda = as.numeric(ib_lambda),
    cohort_col = cohort_col,            # NULL or a single character string
    min_features = as.integer(min_features),
    device = device,
    seed = as.integer(seed)
  )
}


# ----------------------------------------------------------------------------
# Environment resolution from meta_train
# ----------------------------------------------------------------------------

# Resolve the per-row environment label vector from meta_train. Returns a
# character vector of length n (one label per training row) OR NULL if no usable
# cohort column is available. When hp$cohort_col is NULL, the column is
# auto-detected from a fixed list of common names (first match wins). A supplied
# hp$cohort_col that is absent from meta_train resolves to NULL (ERM fallback).
.dg_ibirm_resolve_environments <- function(meta_train, n, cohort_col) {
  if (is.null(meta_train)) return(NULL)
  meta_train <- as.data.frame(meta_train, stringsAsFactors = FALSE)
  if (nrow(meta_train) != n) {
    stop("fit_dg_ibirm: meta_train must have one row per row of X_train")
  }
  col <- cohort_col
  if (is.null(col)) {
    candidates <- c("accession", "cohort", "batch", "study", "dataset", "group")
    hit <- candidates[candidates %in% names(meta_train)]
    if (length(hit) == 0L) return(NULL)
    col <- hit[1L]
  } else if (!col %in% names(meta_train)) {
    return(NULL)
  }
  as.character(meta_train[[col]])
}


# ----------------------------------------------------------------------------
# Python / torch availability + device resolution
# ----------------------------------------------------------------------------

# Clear, actionable error if reticulate or the python torch module is
# unavailable. Plain torch training needs NO brew-OpenSSL bridge and NO HF
# download, so only reticulate + `import torch` are probed.
.dg_ibirm_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_dg_ibirm: package 'reticulate' is required to train the ",
         "embedding. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_dg_ibirm: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.dg_ibirm_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# IB-IRM joint training of encoder + head + weight export
# ----------------------------------------------------------------------------

# Train the BN-free MLP encoder + a Linear(d -> 1) head JOINTLY by the IB-IRM
# objective (Ahuja 2021) across the supplied environments, then export the encoder
# weights (up to and INCLUDING the embedding layer) AND the head (w, b) to R as
# plain numeric matrices/vectors. The python module is used here only; nothing
# python is returned. The torch RNG (CPU + guarded CUDA) is saved/restored, so fit
# leaves the torch generators byte-unchanged; the global R .Random.seed is
# saved/restored by the caller.
#
# `Z` is the n x p rCLR'd, universe-aligned training matrix; `y` the integer 0/1
# labels; `env_idx` an integer vector of length n giving each row's VALID-environment
# index in 1..n_env (n_env >= 1). When n_env == 1, the loss is plain ERM (a single
# pooled BCE risk); when n_env >= 2, the loss is mean BCE + irm_lambda * IRM penalty
# + ib_lambda * IB penalty.
#
# THE IB-IRM PENALTIES:
#   l0_i = w0 . z_i + b0                head logit (z_i = encoder embedding of row i)
#   --- IRMv1 gradient penalty (Arjovsky 2019) ---
#     A FIXED dummy scalar multiplier  w = 1.0  is applied to the head logit:
#       l_e(w) = w * l0   over env e's rows
#     R_e(w) = mean BCE-with-logits( l_e(w) , y_e )
#     gw_e   = d R_e(w) / d w  |_{w=1}     (torch autograd.grad, create_graph=TRUE)
#     Omega_IRM = mean_e ( gw_e ^ 2 )      (w is scalar, so ||.||^2 == square)
#   --- IB (Information-Bottleneck) embedding-variance penalty ---
#     z is the POOLED retained-row embedding (n x d).
#     Omega_IB = mean_j Var_i(z_{ij})      biased 1/n variance, computed explicitly,
#                                          averaged over embedding dims j.
# Total loss = mean_BCE + irm_lambda * Omega_IRM + ib_lambda * Omega_IB. Both
# penalties are differentiable w.r.t. the encoder + head (z carries gradient and the
# IRMv1 grad-of-grad uses create_graph=TRUE), so they shape both at FIT. The EXPORTED
# frozen encoder+head and the PURE-R score are identical in form to dro-vrex/dg-fishr.
.dg_ibirm_train_export <- function(Z, y, env_idx, n_env, hp, device) {
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

  # CRITICAL determinism lesson: SEED BEFORE constructing ANY stochastic torch
  # module. If the encoder/head are built before the seed, their weights
  # initialise from the unseeded process-entry RNG -> reproducible only WITHIN a
  # process (via the on.exit RNG restore), NOT across fresh sessions. Seed first,
  # then build, so two seed=42 fits in two fresh processes init identically.
  torch$manual_seed(as.integer(hp$seed))

  np <- reticulate::import("numpy", delay_load = FALSE)
  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)
  y_py <- torch$as_tensor(np$asarray(as.numeric(y), dtype = "float32")
                          )$float()$reshape(c(-1L, 1L))$to(device)

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
  head <- nn$Linear(as.integer(d), 1L)$to(device)

  params <- c(reticulate::iterate(encoder$parameters()),
              reticulate::iterate(head$parameters()))
  opt <- torch$optim$Adam(params, lr = hp$lr, weight_decay = hp$weight_decay)
  # Per-row BCE-with-logits (reduction='none'); the mean risk is the pooled mean
  # over the retained valid-environment rows. The IRMv1 penalty is per-ENVIRONMENT
  # (the gradient of each env's risk w.r.t. the dummy scale); the IB penalty is over
  # the pooled embedding (representation compression).
  lossf <- nn$BCEWithLogitsLoss(reduction = "none")
  irm_lambda <- hp$irm_lambda
  ib_lambda <- hp$ib_lambda

  # Precompute per-environment row-index tensors (0-based for python indexing).
  # Build as an explicit torch int64 (long) tensor: a numpy int64 array round-trips
  # through reticulate as R numeric and torch$as_tensor would make it a float
  # tensor, which index_select rejects. $long() pins the integer dtype.
  env_tensors <- lapply(seq_len(n_env), function(e) {
    rows0 <- as.integer(which(env_idx == e) - 1L)
    torch$as_tensor(np$asarray(rows0, dtype = "int64"))$long()$to(device)
  })

  encoder$train(); head$train()
  for (ep in seq_len(hp$epochs)) {
    opt$zero_grad()
    emb <- encoder(X_py)                       # n x d
    logit <- head(emb)                         # n x 1, l0 = w0 . z + b0
    per_row <- lossf(logit, y_py)              # n x 1, per-row BCE
    mean_bce <- per_row$mean()                 # pooled mean risk

    if (n_env >= 2L) {
      # ---- IRMv1 gradient penalty (Arjovsky 2019) ---------------------------
      # A FIXED dummy scalar multiplier w = 1.0 applied to the head logit. For each
      # env e: R_e(w) = mean BCE(w * l0_e, y_e); the IRMv1 penalty is the squared
      # gradient of R_e w.r.t. w at w = 1, averaged over environments. autograd.grad
      # with create_graph=TRUE keeps the grad-of-grad in the graph so the penalty is
      # differentiable through the encoder + head.
      dummy_w <- torch$ones(list(1L), dtype = logit$dtype,
                            requires_grad = TRUE)$to(device)
      irm_terms <- lapply(env_tensors, function(rows0) {
        logit_e <- logit$index_select(0L, rows0)        # n_e x 1
        y_e <- y_py$index_select(0L, rows0)             # n_e x 1
        scaled_e <- torch$mul(logit_e, dummy_w)         # w * l0 over env e
        r_e <- lossf(scaled_e, y_e)$mean()              # R_e(w)
        gw <- torch$autograd$grad(r_e, dummy_w,
                                  create_graph = TRUE)[[1L]]  # dR_e/dw, length-1
        gw$pow(2)$sum()                                 # ||dR_e/dw||^2 (scalar)
      })
      omega_irm <- torch$stack(irm_terms)$mean()        # mean_e ||dR_e/dw||^2

      # ---- IB (Information-Bottleneck) embedding-variance penalty ------------
      # Compress the representation: penalise the per-dimension BIASED (1/n)
      # variance of the pooled embedding, averaged over embedding dims. Computed
      # EXPLICITLY (torch$var default is the N-1 unbiased estimator).
      z_mean <- emb$mean(dim = 0L, keepdim = TRUE)      # 1 x d
      z_var <- emb$sub(z_mean)$pow(2)$mean(dim = 0L)    # length d, biased 1/n
      omega_ib <- z_var$mean()                          # mean over embedding dims

      loss <- torch$add(mean_bce,
                        torch$add(omega_irm$mul(irm_lambda),
                                  omega_ib$mul(ib_lambda)))
    } else {
      # ERM fallback: a single pooled risk; both lambdas are inert.
      loss <- mean_bce
    }

    loss$backward()
    opt$step()
  }
  encoder$eval(); head$eval()

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
    stop("fit_dg_ibirm: exported encoder has no linear layers", call. = FALSE)
  }

  # ---- export the jointly-trained head (w, b) -------------------------------
  head_w <- reticulate::py_to_r(
    head$weight$detach()$cpu()$to(torch$float64)$numpy())  # 1 x d
  head_b <- reticulate::py_to_r(
    head$bias$detach()$cpu()$to(torch$float64)$numpy())    # length 1

  list(weights = weights, head_w = as.numeric(head_w),
       head_b = as.numeric(head_b))
}


# ----------------------------------------------------------------------------
# fit_dg_ibirm
# ----------------------------------------------------------------------------

#' @title Fit the IB-IRM cross-cohort transfer single-sample discriminator
#'
#' @description
#' Trains a BatchNorm-free MLP encoder + a linear head \eqn{\mathrm{Linear}(d \to
#' 1)} JOINTLY across training environments (cohorts) by the IB-IRM objective
#' (Ahuja 2021; via reticulate-python torch, Adam, mean per-row risk plus an
#' \eqn{\lambda_{\mathrm{IRM}}}-weighted IRMv1 gradient penalty and an
#' \eqn{\lambda_{\mathrm{IB}}}-weighted Information-Bottleneck embedding-variance
#' penalty), EXPORTS the encoder weights up to and including the embedding layer AND
#' the head \eqn{(w, b)} to R as plain numeric arrays, and DISCARDS the python module.
#' The fitted model holds no external pointer; scoring is pure base-R.
#'
#' Environments are the distinct non-missing labels in
#' \code{meta_train[[cohort_col]]} (\code{cohort_col} auto-detected when
#' \code{NULL}). Only environments carrying both a case and a control are retained;
#' if fewer than two valid environments remain (including \code{NULL}/absent
#' \code{meta_train}), fitting falls back to plain ERM (a single pooled risk, both
#' penalties inert) and \code{model$n_environments} is set to \code{1}.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata carrying a cohort/environment
#'   column. When it contains the resolved \code{cohort_col} with at least two
#'   both-class cohorts, IB-IRM trains across those environments; otherwise the fit
#'   degenerates to ERM. Must have one row per row of \code{X_train}.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{hidden}
#'   (positive-integer vector of encoder layer widths up to the embedding; the LAST
#'   entry is the embedding dim \eqn{d}; default \code{c(64L, 32L)}),
#'   \code{activation} (\code{"relu"} (default) or \code{"tanh"}), \code{epochs}
#'   (training steps, positive integer; default \code{200L}), \code{lr} (positive
#'   Adam learning rate; default \code{1e-3}), \code{weight_decay} (non-negative
#'   Adam L2; default \code{1e-4}), \code{irm_lambda} (non-negative IRMv1
#'   gradient-penalty coefficient; default \code{1.0}), \code{ib_lambda}
#'   (non-negative Information-Bottleneck embedding-variance coefficient; default
#'   \code{1.0}), \code{cohort_col} (\code{NULL} (default, auto-detect) or a single
#'   character string naming the environment column in \code{meta_train}),
#'   \code{min_features} (feature-overlap floor at scoring, positive integer; default
#'   \code{3L}), \code{device} (\code{"cpu"} (default), \code{"cuda"}, or
#'   \code{"auto"}), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{dg_ibirm_model}: a list with
#'   \code{feature_universe}, \code{weights} (exported per-layer \code{list(W, b)}
#'   up to the embedding), \code{activation}, \code{head_w} (length-\eqn{d} frozen
#'   head weight), \code{head_b} (scalar frozen head bias), \code{embedding_dim},
#'   \code{n_environments} (number of valid both-class environments engaged; \code{1}
#'   means ERM fallback), \code{device} (resolved), \code{seed}, and \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 30; k <- 8
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' meta <- data.frame(accession = rep(paste0("GSE", 1:3), length.out = n))
#' model <- fit_dg_ibirm(X, y, meta_train = meta)
#' score_dg_ibirm(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Ahuja K, et al. (2021) Invariance Principle Meets Information Bottleneck for
#' Out-of-Distribution Generalization. \emph{NeurIPS 34}. arXiv:2106.06607.
#'
#' @export
fit_dg_ibirm <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_dg_ibirm", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_dg_ibirm", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_dg_ibirm")
  hp <- .dg_ibirm_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_dg_ibirm: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's
  # first python init seeds R's RNG as a side effect. Fit draws no R randomness of
  # its own (weight init + optimisation are torch-seeded), so this purely
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

  # ---- resolve environments + filter to valid (both-class) ones -------------
  cohort <- .dg_ibirm_resolve_environments(meta_train, nrow(X_train),
                                           hp$cohort_col)
  valid_levels <- character(0)
  if (!is.null(cohort)) {
    non_missing <- !is.na(cohort)
    levels_all <- unique(cohort[non_missing])
    valid_levels <- levels_all[vapply(levels_all, function(lv) {
      rows <- non_missing & cohort == lv
      any(y[rows] == 1L) && any(y[rows] == 0L)
    }, logical(1L))]
  }

  if (length(valid_levels) >= 2L) {
    # IB-IRM path: env_idx maps each row to its valid-environment index; rows whose
    # cohort is missing or single-class are NOT in any valid environment and are
    # therefore excluded from the per-environment gradient set. We restrict training
    # to the valid-env rows so every training row carries a gradient; this keeps the
    # objective a clean mean/penalty over the retained environments.
    keep <- !is.na(cohort) & cohort %in% valid_levels
    env_idx_full <- match(cohort, valid_levels)            # NA outside valid envs
    env_idx <- env_idx_full[keep]
    n_env <- length(valid_levels)
    X_use <- X_train[keep, , drop = FALSE]
    y_use <- y[keep]
  } else {
    # ERM fallback: single pooled environment over all rows; both lambdas inert.
    env_idx <- rep.int(1L, nrow(X_train))
    n_env <- 1L
    X_use <- X_train
    y_use <- y
  }

  .dg_ibirm_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .dg_ibirm_rclr_matrix(X_use)              # rCLR over full universe

  device <- .dg_ibirm_resolve_device(hp$device)
  trained <- .dg_ibirm_train_export(Z_in, y_use, env_idx, n_env, hp, device)
  weights <- trained$weights
  activation <- hp$activation
  # the embedding dim is the OUTPUT of the last linear layer (nrow of its W).
  d <- nrow(weights[[length(weights)]]$W)

  model <- list(
    feature_universe = feat_order,
    weights = weights,
    activation = activation,
    head_w = as.numeric(trained$head_w),
    head_b = as.numeric(trained$head_b),
    embedding_dim = d,
    n_environments = n_env,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "dg_ibirm_model"
  model
}


# Frozen linear logit for one embedded specimen z (length d): w . z + b.
# Larger = more case-like.
.dg_ibirm_score_one <- function(z, head_w, head_b) {
  sum(head_w * z) + head_b
}


# ----------------------------------------------------------------------------
# score_dg_ibirm
# ----------------------------------------------------------------------------

#' @title Score the IB-IRM cross-cohort transfer single-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_dg_ibirm}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), embedded by the exported encoder forward,
#' and scored by the frozen linear logit \eqn{w^\top z + b}. Larger = more
#' case-like. No RNG is used at scoring.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on that
#' row and the frozen model and is exactly invariant to per-specimen positive scaling.
#'
#' @param model A \code{dg_ibirm_model} object from \code{\link{fit_dg_ibirm}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity and
#'   ignored by this method (the transfer aspect is in the fit, not the score).
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values are
#'   more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_dg_ibirm(X, y, meta_train = meta)
#' score_dg_ibirm(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Ahuja K, et al. (2021) Invariance Principle Meets Information Bottleneck for
#' Out-of-Distribution Generalization. \emph{NeurIPS 34}. arXiv:2106.06607.
#'
#' @export
score_dg_ibirm <- function(model, X, meta = NULL) {
  if (!inherits(model, "dg_ibirm_model")) {
    stop("score_dg_ibirm: model must have class dg_ibirm_model")
  }
  X <- .reo_check_matrix(X, "score_dg_ibirm", "X")
  .reo_check_meta(meta, nrow(X), "score_dg_ibirm", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .dg_ibirm_align(X, feat_order)
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
    z <- .dg_ibirm_rclr(X_use[i, ])
    emb <- .dg_ibirm_forward(matrix(z, nrow = 1L), weights, activation)
    out[i] <- .dg_ibirm_score_one(as.numeric(emb), head_w, head_b)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_dg_ibirm: scorer produced non-finite or wrong-length output")
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
.dg_ibirm_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    weights = model$weights,
    activation = model$activation,
    head_w = model$head_w,
    head_b = model$head_b,
    embedding_dim = model$embedding_dim,
    n_environments = model$n_environments,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
