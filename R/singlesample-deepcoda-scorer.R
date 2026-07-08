#' @title DeepCoDA zero-sum log-contrast bottleneck + self-explaining head (single-sample)
#'
#' @description
#' Single-sample within-cohort discriminator (family O, compositional) implementing
#' DeepCoDA in the sense of Quinn, Nguyen, Rana, Gupta & Venkatesh (2020), adapted to
#' the single-specimen scoring contract. A small set of ZERO-SUM LOG-CONTRAST
#' bottlenecks (scale-invariant balances / log-ratios of the per-sample-rCLR
#' specimen) and a SELF-EXPLAINING per-sample predictor over those bottlenecks are
#' learned end-to-end on the TRAINING set by binary cross-entropy (via
#' reticulate-python torch), then \strong{FROZEN}; a new specimen is scored by its
#' frozen DeepCoDA logit. Larger score = more case-like.
#'
#' \strong{Zero-sum log-contrast bottleneck (DeepCoDA eq. 1).} The first layer maps
#' the (per-sample rCLR'd, universe-aligned) input \eqn{x} of length \eqn{p} to
#' \eqn{B} bottleneck "balances"
#' \deqn{b_k(x) = \sum_j W_{kj}\, \mathrm{rclr}(x)_j,}
#' under the HARD CONSTRAINT \eqn{\sum_j W_{kj} = 0} for every bottleneck \eqn{k}
#' (the log-contrast constraint -- each bottleneck is a valid log-ratio, invariant
#' to the specimen total). The constraint is enforced at training by parameterising
#' \eqn{W = W^{\mathrm{raw}} - \mathrm{rowMeans}(W^{\mathrm{raw}})} (each row is
#' re-centred to sum to zero); the centred (zero-sum) \eqn{W} is EXPORTED. There is
#' NO activation on the bottleneck (it is a linear log-contrast), which also gives
#' the pre-activation-embedding consistency of the deepmaha / proto-net / vicreg
#' siblings. Because a zero-sum contrast annihilates any per-row additive constant,
#' applying it to the rCLR is identical to applying it to the raw log-abundance, so
#' the bottleneck -- and therefore the whole score -- is EXACTLY per-sample
#' scale-invariant (verified maxdiff \eqn{\sim 0}, up to float).
#'
#' \strong{Self-explaining per-sample predictor (Quinn 2020).} DeepCoDA's
#' "personalised interpretability" head is a self-explaining predictor: a small
#' BatchNorm-free MLP \eqn{\theta(b(x))} produces, FROM the bottleneck vector,
#' per-sample coefficients \eqn{\theta_k(x)} (length \eqn{B}), and the logit is
#' \deqn{s(x) = \sum_k \theta_k(x)\, b_k(x) + c,}
#' with a learned scalar bias \eqn{c}. The faithful self-explaining form
#' \eqn{\theta(x)\cdot b(x)} is implemented here (the paper's contribution; a global
#' linear head \eqn{s(x)=w\cdot b(x)+c} with constant coefficients is the simpler
#' defensible variant -- see SEVEN_evidence.md for the flagged design choice). The
#' \eqn{\theta}-net uses \code{relu} / \code{tanh} hidden activation and a
#' \eqn{\mathrm{Linear}(\cdot \to B)} output, with NO BatchNorm (every operation is a function of
#' \eqn{x} alone, so the forward is per-row and single-sample-safe).
#'
#' \strong{Torch is confined to FIT.} The bottleneck + \eqn{\theta}-net + bias are
#' trained end-to-end with torch (Adam, \code{BCEWithLogitsLoss}) on the rCLR'd,
#' universe-aligned training matrix. After training, the centred (zero-sum)
#' bottleneck \eqn{W}, the \eqn{\theta}-net weight matrices/biases, and the scalar
#' bias are EXPORTED to R as plain numeric matrices/vectors; the python module is
#' DISCARDED. The fitted model holds NO external pointer and NO python dependency at
#' score time: the log-contrast forward, the \eqn{\theta}-net forward, and the logit
#' are all pure base-R. This makes the default \code{model_digest} a stable
#' deterministic snapshot and lets the package suite run the scorer on any host
#' without the venv.
#'
#' \strong{Fit/score consistency despite float32.} Torch trains in float32; an R
#' float64 forward on the exported (float32-valued) weights differs from the torch
#' forward by \eqn{\sim 10^{-6}}. There are no frozen training statistics beyond the
#' exported weights themselves (the logit is a deterministic function of the exported
#' parameters), so fit and score share one code path by construction -- the scorer
#' evaluates exactly the exported parameters with \code{.coda_deepcoda_logit}, which
#' is also what fit uses to verify the held-out forward (the deepmaha lesson: fit and
#' score must share one code path; compute any frozen quantity from the pure-R
#' forward the scorer uses, NOT from torch float32 embeddings).
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
#' empty-support test is on the ORIGINAL abundances, NOT on \code{all(z == 0)}. (For a
#' flat composition every \eqn{b_k = W_k \cdot 0 = 0}, so the logit reduces to the
#' learned bias \eqn{c}, a genuine computed value, not a floored 0.)
#'
#' @references
#' Quinn TP, Nguyen D, Rana S, Gupta S, Venkatesh S. (2020) DeepCoDA: personalized
#' interpretability for compositional health data. \emph{Proceedings of the 37th
#' International Conference on Machine Learning (ICML)}, PMLR 119. arXiv:2006.01392.
#'
#' @name singlesample-deepcoda
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .proto_net_rclr /
# .ssl_vicreg_rclr / .lrt_deepmaha_rclr.
.coda_deepcoda_rclr <- function(v) {
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
.coda_deepcoda_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .coda_deepcoda_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.coda_deepcoda_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R DeepCoDA forward (exact reimplementation of the torch model)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported theta-net (the self-explaining coefficient
# network). `theta_weights` is a list of per-layer lists list(W = out x in, b =
# out): layers 1..(L-1) are hidden (Linear then `activation`), layer L is the
# coefficient-output Linear (NO activation -- the coefficients are raw, matching
# the torch net). For a layer with weight W (out x in) and bias b, torch computes
# M %*% t(W) + b (row = sample). `B` is the bottleneck matrix (n x B); returns the
# n x B per-sample coefficient matrix theta(b).
.coda_deepcoda_theta_forward <- function(B, theta_weights, activation) {
  act <- if (identical(activation, "tanh")) tanh else function(x) pmax(x, 0)
  L <- length(theta_weights)
  H <- B
  for (l in seq_len(L)) {
    W <- theta_weights[[l]]$W                 # out x in
    b <- theta_weights[[l]]$b                 # length out
    H <- H %*% t(W)                           # n x out
    H <- sweep(H, 2L, b, `+`)
    if (l < L) H <- act(H)                    # no activation after the output
  }
  H
}

# Full DeepCoDA forward on an n x p rCLR matrix Z (universe-aligned). Returns the
# logit vector s(x) = sum_k theta_k(x) * b_k(x) + bias. The bottleneck is the
# zero-sum log-contrast b(x) = Z %*% t(W_bottleneck) (NO activation, NO bias); the
# self-explaining coefficients theta(b) are the theta-net forward; the logit is the
# elementwise product summed over bottlenecks plus the learned scalar bias. This is
# the SINGLE code path used by both fit (held-out verification) and score.
.coda_deepcoda_logit <- function(Z, W_bottleneck, theta_weights, bias, activation) {
  Bmat <- Z %*% t(W_bottleneck)                          # n x B bottlenecks
  Theta <- .coda_deepcoda_theta_forward(Bmat, theta_weights, activation)  # n x B
  rowSums(Theta * Bmat) + bias                           # n logits
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.coda_deepcoda_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_coda_deepcoda: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_coda_deepcoda: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_coda_deepcoda: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("bottleneck_dim", "theta_hidden", "activation", "epochs", "lr",
               "weight_decay", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_coda_deepcoda: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  bottleneck_dim <- hp$bottleneck_dim
  if (is.null(bottleneck_dim)) bottleneck_dim <- 8L
  if (!is.numeric(bottleneck_dim) || length(bottleneck_dim) != 1L ||
      !is.finite(bottleneck_dim) || bottleneck_dim < 1L ||
      bottleneck_dim > .Machine$integer.max ||
      bottleneck_dim != as.integer(bottleneck_dim)) {
    stop("fit_coda_deepcoda: hp$bottleneck_dim must be a positive integer")
  }

  theta_hidden <- hp$theta_hidden
  if (is.null(theta_hidden)) theta_hidden <- c(16L)
  if (!is.numeric(theta_hidden) || length(theta_hidden) < 1L ||
      any(!is.finite(theta_hidden)) || any(theta_hidden < 1L) ||
      any(theta_hidden > .Machine$integer.max) ||
      any(theta_hidden != as.integer(theta_hidden))) {
    stop("fit_coda_deepcoda: hp$theta_hidden must be a vector of positive integers")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_coda_deepcoda: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_coda_deepcoda: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_coda_deepcoda: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_coda_deepcoda: hp$weight_decay must be a non-negative finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_coda_deepcoda: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_coda_deepcoda: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_coda_deepcoda: hp$seed must be a single integer")
  }

  list(
    bottleneck_dim = as.integer(bottleneck_dim),
    theta_hidden = as.integer(theta_hidden),
    activation = activation,
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
# download (no checkpoint, no _ssl), so only reticulate + `import torch` are
# probed.
.coda_deepcoda_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_coda_deepcoda: package 'reticulate' is required to train the ",
         "model. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_coda_deepcoda: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'
# (exact + bit-reproducible for a tiny model).
.coda_deepcoda_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# End-to-end DeepCoDA training + weight export
# ----------------------------------------------------------------------------

# Train the zero-sum log-contrast bottleneck + self-explaining theta-net + scalar
# bias end-to-end via torch (Adam, BCEWithLogits) and return the exported R-side
# parameters: the CENTRED (zero-sum) bottleneck W (B x p), the theta-net per-layer
# list(W, b), and the scalar bias. The python module is used here only; nothing
# python is returned. The torch RNG (CPU + guarded CUDA) is saved/restored, so fit
# leaves the torch generators byte-unchanged; the global R .Random.seed is
# saved/restored by the caller. `Z` is the n x p rCLR'd, universe-aligned training
# matrix; `y` the integer 0/1 labels.
.coda_deepcoda_train_export <- function(Z, y, hp, device) {
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

  # Seed the torch generators (weight init) so two seed-matched CPU fits
  # initialise identically. Training is full-batch (no minibatch shuffling), so
  # the whole fit is bit-reproducible on cpu for a given seed.
  torch$manual_seed(as.integer(hp$seed))

  np <- reticulate::import("numpy", delay_load = FALSE)
  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)
  # BCEWithLogits target in {0,1}, shape n x 1.
  y_py <- torch$as_tensor(np$asarray(as.numeric(y), dtype = "float32")
                          )$float()$reshape(c(-1L, 1L))$to(device)

  p <- ncol(Z)
  B <- as.integer(hp$bottleneck_dim)

  # ---- the DeepCoDA module --------------------------------------------------
  # Bottleneck: a single Linear(p -> B, bias = FALSE). The zero-sum constraint is
  # enforced in the forward by re-centring each row of the weight to sum to zero
  # (W = Wraw - rowMeans(Wraw)), so EVERY forward (and the exported weight) uses a
  # genuine zero-sum log-contrast. Self-explaining theta-net: input(B) -> [Linear
  # -> act]xlen(theta_hidden) -> Linear(-> B) (NO activation on the output). Logit:
  # sum_k theta_k(b) * b_k + bias, with a learned scalar bias.
  act_mod <- if (identical(hp$activation, "tanh")) nn$Tanh else nn$ReLU
  theta_hidden <- hp$theta_hidden

  reticulate::py_run_string("
import torch
import torch.nn as nn

class _DeepCoDA(nn.Module):
    def __init__(self, p, B, theta_hidden, act_name):
        super().__init__()
        self.bottleneck = nn.Linear(p, B, bias=False)
        act = nn.Tanh if act_name == 'tanh' else nn.ReLU
        layers = []
        in_dim = B
        for h in theta_hidden:
            layers.append(nn.Linear(in_dim, int(h)))
            layers.append(act())
            in_dim = int(h)
        layers.append(nn.Linear(in_dim, B))   # coefficient output, no activation
        self.theta_net = nn.Sequential(*layers)
        self.bias = nn.Parameter(torch.zeros(1))

    def zero_sum_W(self):
        W = self.bottleneck.weight                 # B x p
        return W - W.mean(dim=1, keepdim=True)      # re-centre each row to sum 0

    def bottlenecks(self, x):
        return torch.matmul(x, self.zero_sum_W().t())   # n x B (no activation)

    def forward(self, x):
        b = self.bottlenecks(x)                    # n x B
        theta = self.theta_net(b)                  # n x B coefficients
        logit = (theta * b).sum(dim=1, keepdim=True) + self.bias
        return logit
", convert = FALSE)
  main <- reticulate::import_main()
  net <- main$`_DeepCoDA`(as.integer(p), B,
                          reticulate::r_to_py(as.list(as.integer(theta_hidden))),
                          hp$activation)$to(device)

  opt <- torch$optim$Adam(net$parameters(), lr = hp$lr,
                          weight_decay = hp$weight_decay)
  lossf <- nn$BCEWithLogitsLoss()

  net$train()
  for (ep in seq_len(hp$epochs)) {
    opt$zero_grad()
    logit <- net(X_py)
    loss <- lossf(logit, y_py)
    loss$backward()
    opt$step()
  }
  net$eval()

  # ---- export the CENTRED (zero-sum) bottleneck + theta-net + bias ----------
  # Export the zero-sum W (NOT the raw bottleneck.weight): the re-centred weight is
  # the genuine log-contrast the forward uses, so its exported rows sum to zero.
  W_zs <- net$zero_sum_W()$detach()$cpu()$to(torch$float64)$numpy()
  W_bottleneck <- reticulate::py_to_r(W_zs)                     # B x p
  storage.mode(W_bottleneck) <- "double"

  # theta-net: every Linear in order, exported as W (out x in) + b (length out).
  theta_modules <- reticulate::iterate(net$theta_net$modules())
  theta_weights <- list()
  for (m in theta_modules) {
    has_w <- tryCatch(!is.null(m$weight), error = function(e) FALSE)
    if (!isTRUE(has_w)) next
    W <- reticulate::py_to_r(m$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(m$bias$detach()$cpu()$to(torch$float64)$numpy())
    theta_weights[[length(theta_weights) + 1L]] <- list(W = W, b = as.numeric(b))
  }
  if (length(theta_weights) < 1L) {
    stop("fit_coda_deepcoda: exported theta-net has no linear layers", call. = FALSE)
  }

  bias <- as.numeric(
    reticulate::py_to_r(net$bias$detach()$cpu()$to(torch$float64)$numpy()))

  list(W_bottleneck = W_bottleneck, theta_weights = theta_weights, bias = bias)
}


# ----------------------------------------------------------------------------
# fit_coda_deepcoda
# ----------------------------------------------------------------------------

#' @title Fit the DeepCoDA zero-sum log-contrast single-sample discriminator
#'
#' @description
#' Trains a zero-sum log-contrast bottleneck (DeepCoDA eq. 1) plus a self-explaining
#' per-sample predictor end-to-end on the per-sample-rCLR, universe-aligned TRAINING
#' matrix by binary cross-entropy (via reticulate-python torch, Adam,
#' \code{BCEWithLogitsLoss}), EXPORTS the centred (zero-sum) bottleneck weight, the
#' \eqn{\theta}-net weights, and the scalar bias to R as plain numeric
#' matrices/vectors, and DISCARDS the python module. The fitted model holds no
#' external pointer; scoring is pure R.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{bottleneck_dim}
#'   (number of zero-sum log-contrast bottlenecks \eqn{B}, positive integer; default
#'   \code{8L}), \code{theta_hidden} (positive-integer vector of self-explaining
#'   \eqn{\theta}-net hidden widths; default \code{c(16L)}), \code{activation}
#'   (\code{"relu"} (default) or \code{"tanh"}), \code{epochs} (full-batch training
#'   epochs, positive integer; default \code{200L}), \code{lr} (positive Adam
#'   learning rate; default \code{1e-3}), \code{weight_decay} (non-negative Adam L2;
#'   default \code{1e-4}), \code{min_features} (feature-overlap floor at scoring,
#'   positive integer; default \code{3L}), \code{device} (\code{"cpu"} (default),
#'   \code{"cuda"}, or \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU
#'   with no GPU), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{coda_deepcoda_model}: a list with
#'   \code{feature_universe}, \code{W_bottleneck} (the \eqn{B \times p} zero-sum
#'   bottleneck weight), \code{theta_weights} (exported per-layer \code{list(W, b)}
#'   of the self-explaining \eqn{\theta}-net), \code{bias} (scalar), \code{activation},
#'   \code{bottleneck_dim} (\eqn{B}), \code{device} (resolved), \code{seed}, and
#'   \code{hp}.
#'
#' @details
#' The score is the DeepCoDA logit \eqn{s(x) = \sum_k \theta_k(x)\, b_k(x) + c} where
#' \eqn{b_k(x) = \sum_j W_{kj}\,\mathrm{rclr}(x)_j} is the \eqn{k}-th zero-sum
#' log-contrast bottleneck (\eqn{\sum_j W_{kj} = 0}) and \eqn{\theta_k(x)} the
#' self-explaining per-sample coefficient produced by the \eqn{\theta}-net. The
#' zero-sum constraint makes each bottleneck a scale-invariant log-ratio, so the
#' whole score is exactly per-sample scale-invariant.
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
#' model <- fit_coda_deepcoda(X, y)
#' score_coda_deepcoda(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Quinn TP, Nguyen D, Rana S, Gupta S, Venkatesh S. (2020) DeepCoDA: personalized
#' interpretability for compositional health data. \emph{ICML}, PMLR 119.
#' arXiv:2006.01392.
#'
#' @export
fit_coda_deepcoda <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_coda_deepcoda", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_coda_deepcoda", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_coda_deepcoda")
  hp <- .coda_deepcoda_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_coda_deepcoda: X_train must contain at least hp$min_features ",
         "features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's
  # first python init seeds R's RNG as a side effect. Fit draws no R randomness of
  # its own (weight init is torch-seeded), so this purely neutralises the init side
  # effect. Restoring on exit (incl. the no-prior-seed removal case) keeps fit
  # globally RNG-neutral.
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

  .coda_deepcoda_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .coda_deepcoda_rclr_matrix(X_train)       # n x p rCLR over full universe

  device <- .coda_deepcoda_resolve_device(hp$device)
  exp <- .coda_deepcoda_train_export(Z_in, y, hp, device)
  activation <- hp$activation

  model <- list(
    feature_universe = feat_order,
    W_bottleneck = exp$W_bottleneck,
    theta_weights = exp$theta_weights,
    bias = as.numeric(exp$bias),
    activation = activation,
    bottleneck_dim = nrow(exp$W_bottleneck),
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "coda_deepcoda_model"
  model
}


# ----------------------------------------------------------------------------
# score_coda_deepcoda
# ----------------------------------------------------------------------------

#' @title Score the DeepCoDA zero-sum log-contrast single-sample discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_coda_deepcoda}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), projected through the frozen zero-sum
#' log-contrast bottleneck, passed through the self-explaining \eqn{\theta}-net, and
#' scored by the DeepCoDA logit \eqn{\sum_k \theta_k(x)\,b_k(x) + c}. Larger = more
#' case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally (its bottlenecks are all 0, so the logit
#' reduces to the learned bias -- a genuine computed value, not a floored 0). The
#' score of a row depends only on that row and the frozen model and is exactly
#' invariant to per-specimen positive scaling.
#'
#' @param model A \code{coda_deepcoda_model} object from
#'   \code{\link{fit_coda_deepcoda}}.
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
#' model <- fit_coda_deepcoda(X, y)
#' score_coda_deepcoda(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Quinn TP, Nguyen D, Rana S, Gupta S, Venkatesh S. (2020) DeepCoDA: personalized
#' interpretability for compositional health data. \emph{ICML}, PMLR 119.
#' arXiv:2006.01392.
#'
#' @export
score_coda_deepcoda <- function(model, X, meta = NULL) {
  if (!inherits(model, "coda_deepcoda_model")) {
    stop("score_coda_deepcoda: model must have class coda_deepcoda_model")
  }
  X <- .reo_check_matrix(X, "score_coda_deepcoda", "X")
  .reo_check_meta(meta, nrow(X), "score_coda_deepcoda", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .coda_deepcoda_align(X, feat_order)
  W_bottleneck <- model$W_bottleneck
  theta_weights <- model$theta_weights
  bias <- model$bias
  activation <- model$activation

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (rCLR
    # origin) but is a valid specimen that MUST be projected + scored (its score
    # reduces to the learned bias, not a floored 0).
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .coda_deepcoda_rclr(X_use[i, ])
    out[i] <- .coda_deepcoda_logit(matrix(z, nrow = 1L), W_bottleneck,
                                   theta_weights, bias, activation)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_coda_deepcoda: scorer produced non-finite or wrong-length output")
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
.coda_deepcoda_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    W_bottleneck = model$W_bottleneck,
    theta_weights = model$theta_weights,
    bias = model$bias,
    activation = model$activation,
    bottleneck_dim = model$bottleneck_dim,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
