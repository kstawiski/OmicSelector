#' @title Spectral-normalized neural Gaussian process, distance-aware logit (SNGP)
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing the
#' Spectral-normalized Neural Gaussian Process (SNGP) of Liu et al. (2020), adapted
#' to compositional specimens. A small \strong{LayerNorm} MLP encoder with
#' \strong{spectral-normalized} hidden Linear layers is trained on the TRAINING set
#' (via reticulate-python torch) to embed the per-sample-rCLR specimen, then
#' \strong{FROZEN}; a FIXED random-Fourier-feature (RFF) map approximates an RBF
#' Gaussian process on the embedding, and a closed-form LAPLACE posterior on the GP
#' output weights gives the distance-aware posterior-mean logit. A new specimen is
#' scored by its mean-field-adjusted GP logit. Larger score = more case-like.
#'
#' \strong{Distance awareness (Liu 2020).} SNGP combines two ingredients to make a
#' deterministic deep net distance-aware. (i) The hidden Linear layers are
#' \strong{spectral-normalized} so the encoder is approximately bi-Lipschitz and
#' preserves input distances in the embedding space; combined with
#' \strong{LayerNorm} (which normalises ACROSS FEATURES WITHIN A ROW, so it is
#' per-row and single-sample-exact -- BatchNorm, which couples rows, is forbidden
#' here). (ii) A Gaussian-process output layer, approximated by random Fourier
#' features with a closed-form Laplace posterior, whose predictive variance grows
#' with distance from the training data, so far-from-data specimens get a shrunken
#' (less confident) logit. This is exactly the SNGP recipe.
#'
#' \strong{Torch is confined to FIT.} The encoder is trained with torch (Adam, CE
#' loss) on the rCLR'd, universe-aligned training matrix with
#' \code{torch.nn.utils.parametrizations.spectral_norm} on the hidden Linear
#' layers. After training, the \emph{effective} (post-spectral-norm, coefficient-
#' scaled) Linear weights \eqn{W_{\mathrm{eff}} = c\,W_{\mathrm{SN}}} -- the exact
#' weights torch used in its forward -- the biases, and the LayerNorm
#' \eqn{\gamma,\beta,\varepsilon} are EXPORTED to R as plain numeric
#' matrices/vectors; the classification head and the python module are DISCARDED.
#' The fixed RFF projection \eqn{W_{\mathrm{rff}}, b_{\mathrm{rff}}} (drawn ONCE
#' under a dedicated seeded torch generator) and the Laplace posterior
#' (\eqn{\beta}, \eqn{\Sigma}) are likewise frozen R-side. The fitted model
#' therefore holds NO external pointer and NO python dependency at score time: the
#' encoder forward (incl. LayerNorm), the RFF cosine map, and the mean-field GP
#' logit are all pure base-R. This makes the default \code{model_digest} a stable
#' deterministic snapshot and lets the package suite run the scorer on any host
#' without the venv.
#'
#' \strong{Spectral-norm bake-at-export fidelity.} PyTorch's \code{spectral_norm}
#' parametrisation materialises \code{layer$weight} to the spectrally-normalised
#' weight \eqn{W_{\mathrm{SN}}} (spectral norm \eqn{1}) at every forward; the
#' encoder's forward is \eqn{x\,(c\,W_{\mathrm{SN}})^\top + b} with the bound
#' coefficient \eqn{c=}\code{sn_coeff}. Since the training forward scales the whole
#' linear output by \eqn{c} (\eqn{c\,(H W_{\mathrm{SN}}^\top + b) = H\,(c
#' W_{\mathrm{SN}})^\top + c\,b}), we export BOTH the effective weight
#' \eqn{W_{\mathrm{eff}} = c\,W_{\mathrm{SN}}} and the effective bias
#' \eqn{b_{\mathrm{eff}} = c\,b}, so the pure-R forward is a plain matmul of the
#' SAME weights AND bias torch used for ANY \code{sn_coeff} (not only \eqn{c=1}):
#' the R float64 forward matches the torch float32 forward to \eqn{\sim 10^{-7}}
#' (the float32 round-off only).
#'
#' \strong{Fit/score consistency.} The Laplace posterior (\eqn{\beta},
#' \eqn{\Sigma}) is computed from the PURE-R RFF features
#' \eqn{\Phi(\,.\mathrm{unc\_sngp\_forward}(Z_{\mathrm{train}})\,)} -- the IDENTICAL
#' code path the scorer evaluates -- NOT from a torch GP head. Fit and score thus
#' share one computation (the deepmaha / proto-net lesson: fit and score must share
#' one code path).
#'
#' \strong{The RFF Gaussian-process head.} With the frozen embedding
#' \eqn{z = \mathrm{enc}(\mathrm{rclr}(x))} the random-Fourier feature map is
#' \deqn{\Phi(z) = \sqrt{2/D}\,\cos(z\,W_{\mathrm{rff}}^\top + b_{\mathrm{rff}}),}
#' \eqn{W_{\mathrm{rff}}\sim\mathcal N(0, 1/\ell^2)} (length-scale \eqn{\ell}=
#' \code{rff_ls}), \eqn{b_{\mathrm{rff}}\sim\mathrm{Unif}(0, 2\pi)},
#' \eqn{D=}\code{rff_dim}; this approximates an RBF kernel GP. A ridge-regularised
#' logistic regression of the binary label on \eqn{\Phi(Z_{\mathrm{train}})} (a few
#' Newton/IRLS steps) gives the Laplace posterior MEAN \eqn{\beta} and PRECISION
#' \eqn{\Sigma^{-1} = \lambda I + \sum_i p_i(1-p_i)\,\phi_i\phi_i^\top}
#' (\eqn{\lambda=}\code{ridge}), and \eqn{\Sigma = (\Sigma^{-1})^{-1}}.
#'
#' \strong{Distance-aware mean-field score.} The score of a new specimen is the
#' mean-field-approximated posterior logit (Liu 2020):
#' \deqn{s(x) = \frac{\Phi(z)\cdot\beta}
#'   {\sqrt{1 + (\pi/8)\,\Phi(z)^\top\Sigma\,\Phi(z)}}.}
#' The denominator IS the distance awareness: a specimen far from the training
#' manifold has a large RFF-GP predictive variance \eqn{\Phi^\top\Sigma\Phi} and so
#' a shrunken (less confident) logit. We use the FULL mean-field denominator (not
#' the posterior mean alone): it is the defining SNGP predictive and is the
#' distance-aware mechanism; within-cohort ranking is essentially preserved
#' (Spearman \eqn{\approx 1} vs the mean alone) while the score is the faithful
#' distance-aware logit. Larger = more case-like.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly
#' invariant to per-specimen scaling and uses no cross-row statistic. Universe
#' features absent from a specimen carry the neutral rCLR value \code{0}. The SAME
#' rCLR transform is applied to the frozen training rows at fit and to every query
#' at score; combined with LayerNorm (per-row), the fixed RFF map, and the frozen
#' posterior, a specimen's score depends only on that specimen and the frozen model
#' -- single-row == batch EXACTLY, and it passes
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
#' Liu JZ, Lin Z, Padhy S, Tran D, Bedrax-Weiss T, Lakshminarayanan B. (2020)
#' Simple and Principled Uncertainty Estimation with Deterministic Deep Learning
#' via Distance Awareness. \emph{Advances in Neural Information Processing Systems
#' (NeurIPS)} 33. arXiv:2006.10108.
#'
#' @name singlesample-sngp
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .lrt_deepmaha_rclr /
# .proto_net_rclr.
.unc_sngp_rclr <- function(v) {
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
.unc_sngp_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .unc_sngp_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.unc_sngp_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R encoder forward (exact reimplementation of the torch LayerNorm SN net)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported spectral-normalized, LayerNorm MLP encoder.
# `weights` is a per-layer list(list(W = out x in, b = out)) where W is the
# EFFECTIVE post-spectral-norm, coefficient-scaled weight W_eff = c * W_SN (the
# exact weight torch used in its forward). `layernorm` is a per-layer
# list(list(gamma = out, beta = out, eps = scalar)). Each layer l is
#   H <- (H %*% t(W_l) + b_l)                  # spectral-norm-baked Linear
#   H <- gamma_l * (H - mean(H)) / sqrt(var(H) + eps_l) + beta_l   # LayerNorm
# applied row-wise (mean/var ACROSS FEATURES WITHIN A ROW -- per-row, batch-free),
# and `act` is applied after every LayerNorm EXCEPT the final (embedding) layer
# (matching the torch encoder). torch's LayerNorm uses the BIASED variance (1/N),
# so we use rowMeans of squared deviations (NOT stats::var's 1/(N-1)). Returns the
# n x d embedding. This is per-row exact, so single-row == batch.
.unc_sngp_forward <- function(M, weights, layernorm, activation) {
  act <- if (identical(activation, "tanh")) tanh else function(x) pmax(x, 0)
  L <- length(weights)
  H <- M
  for (l in seq_len(L)) {
    W <- weights[[l]]$W                       # out x in (W_eff)
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


# Pure-R random-Fourier-feature map: Phi(Z) = sqrt(2/D) * cos(Z W_rff^T + b_rff).
# `W_rff` is D x d, `b_rff` length D, `rff_scale` = sqrt(2/D). Z is n x d (the
# encoder embedding); returns n x D. Identical at fit and score (frozen draw).
.unc_sngp_rff <- function(Z, W_rff, b_rff, rff_scale) {
  P <- Z %*% t(W_rff)                         # n x D
  P <- sweep(P, 2L, b_rff, `+`)
  rff_scale * cos(P)
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.unc_sngp_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_unc_sngp: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_unc_sngp: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_unc_sngp: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("hidden", "activation", "epochs", "lr", "weight_decay",
               "sn_coeff", "rff_dim", "rff_ls", "ridge", "min_features",
               "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_unc_sngp: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  hidden <- hp$hidden
  if (is.null(hidden)) hidden <- c(64L, 32L)
  if (!is.numeric(hidden) || length(hidden) < 1L || any(!is.finite(hidden)) ||
      any(hidden < 1L) || any(hidden > .Machine$integer.max) ||
      any(hidden != as.integer(hidden))) {
    stop("fit_unc_sngp: hp$hidden must be a vector of positive integers")
  }

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L ||
      is.na(activation) || !activation %in% c("relu", "tanh")) {
    stop("fit_unc_sngp: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_unc_sngp: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_unc_sngp: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_unc_sngp: hp$weight_decay must be a non-negative finite number")
  }

  sn_coeff <- hp$sn_coeff
  if (is.null(sn_coeff)) sn_coeff <- 1.0
  if (!is.numeric(sn_coeff) || length(sn_coeff) != 1L || !is.finite(sn_coeff) ||
      sn_coeff <= 0) {
    stop("fit_unc_sngp: hp$sn_coeff must be a positive finite number")
  }

  rff_dim <- hp$rff_dim
  if (is.null(rff_dim)) rff_dim <- 128L
  if (!is.numeric(rff_dim) || length(rff_dim) != 1L || !is.finite(rff_dim) ||
      rff_dim < 1L || rff_dim > .Machine$integer.max ||
      rff_dim != as.integer(rff_dim)) {
    stop("fit_unc_sngp: hp$rff_dim must be a positive integer")
  }

  rff_ls <- hp$rff_ls
  if (is.null(rff_ls)) rff_ls <- 1.0
  if (!is.numeric(rff_ls) || length(rff_ls) != 1L || !is.finite(rff_ls) ||
      rff_ls <= 0) {
    stop("fit_unc_sngp: hp$rff_ls must be a positive finite number")
  }

  ridge <- hp$ridge
  if (is.null(ridge)) ridge <- 1.0
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) ||
      ridge <= 0) {
    stop("fit_unc_sngp: hp$ridge must be a positive finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_unc_sngp: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_unc_sngp: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_unc_sngp: hp$seed must be a single integer")
  }

  list(
    hidden = as.integer(hidden),
    activation = activation,
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    sn_coeff = as.numeric(sn_coeff),
    rff_dim = as.integer(rff_dim),
    rff_ls = as.numeric(rff_ls),
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
# download, so only reticulate + `import torch` are probed.
.unc_sngp_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_unc_sngp: package 'reticulate' is required to train the ",
         "encoder. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_unc_sngp: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.unc_sngp_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# Spectral-normalized LayerNorm encoder training + weight/RFF export
# ----------------------------------------------------------------------------

# Train the spectral-normalized, LayerNorm MLP encoder + a transient 2-logit head
# via torch (Adam, CE loss) and return (a) the exported R-side encoder weights (a
# per-layer list(W = c * W_SN effective post-spectral-norm weight, b)), (b) the
# per-layer LayerNorm parameters list(gamma, beta, eps), and (c) the FIXED RFF
# projection (W_rff D x d, b_rff length D) drawn ONCE under a DEDICATED seeded
# torch generator. The python module is used here only; nothing python is
# returned. The torch CPU/CUDA RNG (for weight init + training) is saved/restored,
# so fit leaves the torch generators byte-unchanged; the global R .Random.seed is
# saved/restored by the caller. The RFF draw uses a dedicated Generator (seeded
# from hp$seed) so it (i) is bit-reproducible across fits and (ii) does NOT perturb
# the global torch RNG. `Z` is the n x p rCLR'd, universe-aligned training matrix;
# `y` the integer 0/1 labels.
.unc_sngp_train_export <- function(Z, y, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn
  np <- reticulate::import("numpy", delay_load = FALSE)
  sn_fun <- torch$nn$utils$parametrizations$spectral_norm

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

  # Seed the global torch generator: weight init + any torch-side randomness so
  # two seed-matched CPU fits initialise identically.
  torch$manual_seed(as.integer(hp$seed))

  X_py <- torch$as_tensor(np$asarray(Z, dtype = "float32"))$float()$to(device)
  y_py <- torch$as_tensor(np$asarray(as.integer(y), dtype = "int64"))$long()$to(device)

  p <- ncol(Z)
  hidden <- hp$hidden
  c_sn <- as.numeric(hp$sn_coeff)

  # Build a module list: per layer a spectral-normalized Linear + a LayerNorm.
  # The encoder applies, per layer, H <- c * SN(Linear)(H) ; H <- LayerNorm(H) ;
  # and activation after every layer EXCEPT the final (embedding) layer. A
  # transient Linear(d -> 2) classification head trains the encoder (CE loss) and
  # is discarded after fit. spectral_norm replaces lin$weight with the spectrally
  # normalised weight at every forward (spectral norm 1); we multiply the linear
  # output by c so the effective weight is c * W_SN (the SNGP coefficient bound).
  lins <- list()
  lns <- list()
  in_dim <- p
  for (h in hidden) {
    ho <- as.integer(h)
    lin <- sn_fun(nn$Linear(in_dim, ho), name = "weight",
                  n_power_iterations = 1L)
    lins[[length(lins) + 1L]] <- lin
    lns[[length(lns) + 1L]] <- nn$LayerNorm(ho)
    in_dim <- ho
  }
  d <- in_dim
  head <- nn$Linear(d, 2L)

  # Register all modules so their parameters are optimised. We keep python
  # references in R lists and run the forward manually (so the c-scaling and the
  # LayerNorm placement match .unc_sngp_forward exactly).
  params <- list()
  for (lin in lins) params <- c(params, reticulate::iterate(lin$parameters()))
  for (ln in lns)   params <- c(params, reticulate::iterate(ln$parameters()))
  params <- c(params, reticulate::iterate(head$parameters()))
  opt <- torch$optim$Adam(params, lr = hp$lr, weight_decay = hp$weight_decay)
  lossf <- nn$CrossEntropyLoss()
  L <- length(lins)

  fwd <- function(xb) {
    H <- xb
    for (l in seq_len(L)) {
      lin <- lins[[l]]
      H <- torch$mul(lin(H), c_sn)            # c * SN(Linear)(H)
      H <- lns[[l]](H)                        # LayerNorm
      if (l < L) {
        H <- if (identical(hp$activation, "tanh")) torch$tanh(H) else torch$relu(H)
      }
    }
    H
  }

  for (lin in lins) lin$train()
  for (ln in lns)   ln$train()
  head$train()
  n <- nrow(Z)
  bs <- min(32L, n)
  for (ep in seq_len(hp$epochs)) {
    starts <- seq.int(0L, n - 1L, by = bs)
    for (s0 in starts) {
      idx <- torch$arange(as.integer(s0),
                          as.integer(min(s0 + bs, n)),
                          dtype = torch$long)$to(device)
      xb <- X_py$index_select(0L, idx)
      yb <- y_py$index_select(0L, idx)
      opt$zero_grad()
      emb <- fwd(xb)
      out <- head(emb)
      loss <- lossf(out, yb)
      loss$backward()
      opt$step()
    }
  }
  for (lin in lins) lin$eval()
  for (ln in lns)   ln$eval()
  head$eval()

  # ---- export the EFFECTIVE post-spectral-norm encoder weights --------------
  # spectral_norm materialises lin$weight to the spectrally-normalised weight at
  # access. The training forward scales the WHOLE linear output by c:
  #   c * (H W_SN^T + b) = H (c W_SN)^T + c b,
  # so BOTH the weight and the bias are scaled by c. Export the effective weight
  # W_eff = c W_SN AND the effective bias b_eff = c b, so the pure-R forward
  # `H W_eff^T + b_eff` reproduces the exact trained forward for ANY sn_coeff
  # (not only the default c = 1, where c b = b). Export as W (out x in) and b
  # (out) in float64.
  weights <- vector("list", L)
  layernorm <- vector("list", L)
  for (l in seq_len(L)) {
    lin <- lins[[l]]
    W_sn <- reticulate::py_to_r(lin$weight$detach()$cpu()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(lin$bias$detach()$cpu()$to(torch$float64)$numpy())
    weights[[l]] <- list(W = c_sn * W_sn, b = c_sn * as.numeric(b))
    ln <- lns[[l]]
    g  <- reticulate::py_to_r(ln$weight$detach()$cpu()$to(torch$float64)$numpy())
    be <- reticulate::py_to_r(ln$bias$detach()$cpu()$to(torch$float64)$numpy())
    eps <- as.numeric(reticulate::py_to_r(ln$eps))
    layernorm[[l]] <- list(gamma = as.numeric(g), beta = as.numeric(be),
                           eps = eps)
  }

  # ---- draw the FIXED RFF projection under a DEDICATED seeded generator ------
  # W_rff ~ N(0, 1/ls^2)  (so std = 1/ls);  b_rff ~ Unif(0, 2*pi). A dedicated
  # Generator seeded from hp$seed makes the draw bit-reproducible AND leaves the
  # global torch RNG untouched. Exported in float64 (frozen forever after).
  Dr <- as.integer(hp$rff_dim)
  gen <- torch$Generator()
  gen$manual_seed(as.integer(hp$seed))
  W_rff <- torch$randn(c(Dr, as.integer(d)), generator = gen,
                       dtype = torch$float64)$mul(1.0 / hp$rff_ls)
  b_rff <- torch$rand(Dr, generator = gen, dtype = torch$float64)$mul(2 * pi)
  W_rff <- reticulate::py_to_r(W_rff$cpu()$numpy())
  b_rff <- as.numeric(reticulate::py_to_r(b_rff$cpu()$numpy()))
  if (!is.matrix(W_rff)) W_rff <- matrix(W_rff, nrow = Dr)

  list(weights = weights, layernorm = layernorm, embedding_dim = d,
       W_rff = W_rff, b_rff = b_rff)
}


# ----------------------------------------------------------------------------
# Closed-form Laplace posterior on the RFF-GP output weights (pure R, IRLS)
# ----------------------------------------------------------------------------

# Ridge-regularised logistic regression of the 0/1 label on the RFF features Phi
# (n x D) by IRLS / Newton steps, giving the Laplace posterior MEAN beta (length
# D) and PRECISION Sigma_inv = ridge*I + Phi^T diag(p(1-p)) Phi, then the
# covariance Sigma = solve(Sigma_inv). The ridge IS the Gaussian prior precision
# lambda (SNGP's closed-form last-layer Laplace). No intercept term: the RFF map
# already carries an offset-like component via b_rff, and a centred cos basis is
# the standard RFF-GP parameterisation. Returns list(beta, Sigma, Sigma_inv).
.unc_sngp_laplace <- function(Phi, y, ridge, n_iter = 50L) {
  D <- ncol(Phi)
  lam <- as.numeric(ridge)
  yv <- as.numeric(y)
  beta <- rep.int(0, D)
  Ridge <- lam * diag(D)
  for (it in seq_len(n_iter)) {
    eta <- as.numeric(Phi %*% beta)
    pr <- 1 / (1 + exp(-eta))
    w <- pr * (1 - pr)
    w <- pmax(w, 1e-10)                       # PD floor on the IRLS weights
    grad <- as.numeric(crossprod(Phi, pr - yv)) + lam * beta
    H <- crossprod(Phi, Phi * w) + Ridge      # Phi^T diag(w) Phi + lambda I
    step <- tryCatch(solve(H, grad), error = function(e) {
      solve(H + 1e-8 * diag(D), grad)
    })
    beta_new <- beta - step
    if (max(abs(beta_new - beta)) < 1e-10) { beta <- beta_new; break }
    beta <- beta_new
  }
  # Laplace precision at the posterior mode.
  eta <- as.numeric(Phi %*% beta)
  pr <- 1 / (1 + exp(-eta))
  w <- pmax(pr * (1 - pr), 1e-10)
  Sigma_inv <- crossprod(Phi, Phi * w) + Ridge
  Sigma_inv <- (Sigma_inv + t(Sigma_inv)) / 2
  Sigma <- chol2inv(chol(Sigma_inv))
  Sigma <- (Sigma + t(Sigma)) / 2
  list(beta = as.numeric(beta), Sigma = Sigma, Sigma_inv = Sigma_inv)
}


# ----------------------------------------------------------------------------
# fit_unc_sngp
# ----------------------------------------------------------------------------

#' @title Fit the spectral-normalized neural Gaussian process discriminator (SNGP)
#'
#' @description
#' Trains a spectral-normalized, LayerNorm MLP encoder on the per-sample-rCLR,
#' universe-aligned TRAINING matrix (via reticulate-python torch; Adam,
#' cross-entropy on a transient head), EXPORTS the effective
#' (post-spectral-norm, coefficient-scaled) encoder weights, the LayerNorm
#' parameters, and a FIXED random-Fourier-feature projection to R as plain numeric
#' matrices/vectors, DISCARDS the python module and the head, then fits a
#' closed-form ridge-Laplace logistic posterior (\eqn{\beta}, \eqn{\Sigma}) on the
#' PURE-R RFF features \eqn{\Phi(\mathrm{enc}(Z_{\mathrm{train}}))} the scorer
#' evaluates (guaranteeing fit/score consistency). The fitted model holds no
#' external pointer.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{hidden}
#'   (positive-integer vector of encoder layer widths; the LAST entry is the
#'   embedding dim \eqn{d}; default \code{c(64L, 32L)}), \code{activation}
#'   (\code{"relu"} (default) or \code{"tanh"}), \code{epochs} (positive integer;
#'   default \code{200L}), \code{lr} (positive Adam learning rate; default
#'   \code{1e-3}), \code{weight_decay} (non-negative Adam L2; default \code{1e-4}),
#'   \code{sn_coeff} (positive spectral-norm bound \eqn{c}; default \code{1.0}),
#'   \code{rff_dim} (positive integer number of random Fourier features \eqn{D};
#'   default \code{128L}), \code{rff_ls} (positive RBF length-scale \eqn{\ell};
#'   default \code{1.0}), \code{ridge} (positive Laplace prior precision
#'   \eqn{\lambda}; default \code{1.0}), \code{min_features} (feature-overlap floor
#'   at scoring, positive integer; default \code{3L}), \code{device} (\code{"cpu"}
#'   (default), \code{"cuda"}, or \code{"auto"}; \code{"cuda"}/\code{"auto"} fall
#'   back to CPU with no GPU), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{unc_sngp_model}: a list with
#'   \code{feature_universe}, \code{weights} (exported per-layer \code{list(W, b)}
#'   effective post-spectral-norm encoder), \code{layernorm} (per-layer
#'   \code{list(gamma, beta, eps)}), \code{activation}, \code{W_rff}, \code{b_rff},
#'   \code{rff_scale} (\eqn{\sqrt{2/D}}), \code{beta} (Laplace posterior mean),
#'   \code{Sigma} (Laplace posterior covariance), \code{Sigma_inv} (precision),
#'   \code{embedding_dim}, \code{rff_dim}, \code{device} (resolved), \code{seed},
#'   and \code{hp}.
#'
#' @details
#' The score is the mean-field-adjusted RFF-GP posterior logit
#' \eqn{s(x) = \Phi(z)\cdot\beta / \sqrt{1 + (\pi/8)\,\Phi(z)^\top\Sigma\,\Phi(z)}},
#' with \eqn{z = \mathrm{enc}(\mathrm{rclr}(x))} and
#' \eqn{\Phi(z) = \sqrt{2/D}\cos(z W_{\mathrm{rff}}^\top + b_{\mathrm{rff}})}. The
#' denominator is the SNGP distance awareness: the RFF-GP predictive variance
#' \eqn{\Phi^\top\Sigma\Phi} grows with distance from the training data, shrinking
#' the logit of far-from-data specimens. The Laplace posterior uses
#' \eqn{\Sigma^{-1} = \lambda I + \sum_i p_i(1-p_i)\phi_i\phi_i^\top}.
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
#' model <- fit_unc_sngp(X, y)
#' score_unc_sngp(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Liu JZ, Lin Z, Padhy S, Tran D, Bedrax-Weiss T, Lakshminarayanan B. (2020)
#' Simple and Principled Uncertainty Estimation with Deterministic Deep Learning
#' via Distance Awareness. \emph{NeurIPS} 33. arXiv:2006.10108.
#'
#' @export
fit_unc_sngp <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_unc_sngp", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_unc_sngp", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_unc_sngp")
  hp <- .unc_sngp_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_unc_sngp: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's
  # first python init seeds R's RNG as a side effect. Fit draws no R randomness of
  # its own (the RFF draw is on a dedicated torch generator); restoring on exit
  # (incl. the no-prior-seed case) keeps fit globally RNG-neutral.
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

  .unc_sngp_require_torch()

  feat_order <- colnames(X_train)
  Z_in <- .unc_sngp_rclr_matrix(X_train)            # n x p rCLR over full universe

  device <- .unc_sngp_resolve_device(hp$device)
  ex <- .unc_sngp_train_export(Z_in, y, hp, device)
  weights <- ex$weights
  layernorm <- ex$layernorm
  d <- ex$embedding_dim
  W_rff <- ex$W_rff
  b_rff <- ex$b_rff
  activation <- hp$activation
  rff_scale <- sqrt(2 / hp$rff_dim)

  # Laplace posterior computed from the PURE-R RFF features the scorer uses (NOT a
  # torch GP head): fit and score share one code path (float32 -> float64 safe).
  Z_emb <- .unc_sngp_forward(Z_in, weights, layernorm, activation)   # n x d
  Phi <- .unc_sngp_rff(Z_emb, W_rff, b_rff, rff_scale)               # n x D
  lap <- .unc_sngp_laplace(Phi, y, hp$ridge)

  model <- list(
    feature_universe = feat_order,
    weights = weights,
    layernorm = layernorm,
    activation = activation,
    W_rff = W_rff,
    b_rff = b_rff,
    rff_scale = rff_scale,
    beta = lap$beta,
    Sigma = lap$Sigma,
    Sigma_inv = lap$Sigma_inv,
    embedding_dim = d,
    rff_dim = as.integer(hp$rff_dim),
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "unc_sngp_model"
  model
}


# Distance-aware mean-field GP logit for one RFF feature vector phi (length D):
#   s = (phi . beta) / sqrt(1 + (pi/8) * phi' Sigma phi).  Larger = more case-like.
.unc_sngp_score_one <- function(phi, beta, Sigma) {
  mean_logit <- as.numeric(crossprod(phi, beta))
  var_term <- as.numeric(crossprod(phi, Sigma %*% phi))
  if (!is.finite(var_term) || var_term < 0) var_term <- 0
  mean_logit / sqrt(1 + (pi / 8) * var_term)
}


# ----------------------------------------------------------------------------
# score_unc_sngp
# ----------------------------------------------------------------------------

#' @title Score the spectral-normalized neural Gaussian process discriminator (SNGP)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_unc_sngp}}, in PURE R (no python). Each query is mapped to the
#' per-sample robust CLR over the frozen feature universe (absent universe features
#' carry the neutral rCLR value \code{0}), embedded by the exported
#' spectral-normalized LayerNorm encoder forward, projected by the fixed
#' random-Fourier-feature map, and scored by the mean-field-adjusted RFF-GP
#' posterior logit. Larger = more case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on
#' that row and the frozen model and is exactly invariant to per-specimen positive
#' scaling.
#'
#' @param model A \code{unc_sngp_model} object from \code{\link{fit_unc_sngp}}.
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
#' model <- fit_unc_sngp(X, y)
#' score_unc_sngp(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Liu JZ, Lin Z, Padhy S, Tran D, Bedrax-Weiss T, Lakshminarayanan B. (2020)
#' Simple and Principled Uncertainty Estimation with Deterministic Deep Learning
#' via Distance Awareness. \emph{NeurIPS} 33. arXiv:2006.10108.
#'
#' @export
score_unc_sngp <- function(model, X, meta = NULL) {
  if (!inherits(model, "unc_sngp_model")) {
    stop("score_unc_sngp: model must have class unc_sngp_model")
  }
  X <- .reo_check_matrix(X, "score_unc_sngp", "X")
  .reo_check_meta(meta, nrow(X), "score_unc_sngp", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  X_use <- .unc_sngp_align(X, feat_order)
  weights <- model$weights
  layernorm <- model$layernorm
  activation <- model$activation
  W_rff <- model$W_rff
  b_rff <- model$b_rff
  rff_scale <- model$rff_scale
  beta <- model$beta
  Sigma <- model$Sigma

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition
    # (full positive support, equal abundances) also has an all-zero rCLR (rCLR
    # origin) but is a valid specimen that MUST be embedded and scored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .unc_sngp_rclr(X_use[i, ])
    emb <- .unc_sngp_forward(matrix(z, nrow = 1L), weights, layernorm, activation)
    phi <- as.numeric(.unc_sngp_rff(emb, W_rff, b_rff, rff_scale))
    out[i] <- .unc_sngp_score_one(phi, beta, Sigma)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_unc_sngp: scorer produced non-finite or wrong-length output")
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
# identical data produce identical digests (encoder weights + LayerNorm + the
# seeded RFF draw + the Laplace posterior all seeded/deterministic).
.unc_sngp_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    weights = model$weights,
    layernorm = model$layernorm,
    activation = model$activation,
    W_rff = model$W_rff,
    b_rff = model$b_rff,
    rff_scale = model$rff_scale,
    beta = model$beta,
    Sigma = model$Sigma,
    Sigma_inv = model$Sigma_inv,
    embedding_dim = model$embedding_dim,
    rff_dim = model$rff_dim,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
