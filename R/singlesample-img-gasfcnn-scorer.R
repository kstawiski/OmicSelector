#' @title GASF-image CNN single-sample discriminator (frozen-BatchNorm, python-at-score)
#'
#' @description
#' Single-sample within-cohort discriminator (family M, image/invariance) that turns
#' each specimen's per-sample robust-CLR (rCLR) profile into a Gramian Angular
#' Summation Field (GASF) image (Wang & Oates, 2015) and scores it through a small
#' convolutional network with FROZEN BatchNorm. A specimen of \eqn{p} features maps
#' to a single-channel \eqn{p \times p} GASF image; a small Conv-BN-ReLU-MaxPool CNN
#' with a global-average-pool logit head is trained end-to-end on the TRAINING images
#' by binary cross-entropy, then FROZEN. A new specimen is scored by the logit of the
#' frozen, eval-mode CNN on its GASF image. Larger = more case-like.
#'
#' \strong{Python-at-score contract (no live pointer).} Unlike the BatchNorm-free
#' trained-torch siblings (\code{proto-net}, \code{lrt-deepmaha}) which export weights
#' to a pure-R forward, a CNN with Conv2d/MaxPool/BatchNorm has no compact pure-R
#' forward, so this scorer follows the PROVEN PYTHON-AT-SCORE pattern of
#' \code{\link{fit_tabicl}} / \code{\link{fit_tabpfn}}: the fitted model stores ONLY
#' plain R-side state -- the CNN \code{state_dict} as named float64 numeric arrays
#' (every conv/BN/linear parameter, INCLUDING the BN \code{running_mean},
#' \code{running_var}, \code{num_batches_tracked}), the architecture config, the
#' frozen per-feature GASF min/max bounds, the feature universe, and the seed/hp.
#' There is NO live python pointer on the model. At SCORE, the identical CNN is rebuilt
#' in python from that R state, \code{load_state_dict} loads the exported weights,
#' \code{net.eval()} puts BatchNorm in FROZEN mode (running stats, not per-batch
#' stats), and each query's GASF image is forwarded ONE ROW AT A TIME. Because the
#' R-side state is fully serialisable, the default \code{model_digest} is a stable
#' deterministic snapshot suitable for \code{\link{singlesample_assert_row_equivariant}}.
#'
#' \strong{Frozen BatchNorm + forced row-by-row = exact single-sample equivariance.}
#' Two ingredients make the canonical gate hold EXACTLY. (1) \code{net.eval()} freezes
#' BatchNorm to its training running statistics, so a query's normalisation depends
#' only on the FROZEN stats -- never on the other rows in the scored batch. In
#' \code{net.train()} mode BatchNorm would compute per-batch statistics and COUPLE the
#' rows, breaking single-sample deployability; eval mode is therefore mandatory and is
#' asserted by a dedicated test (§7.10). (2) On a TRAINED net the float32 conv
#' accumulation across the batch dimension makes a 1-row forward differ from the same
#' row inside an \eqn{n>1} batch by \eqn{\sim 2\times 10^{-7}} (a tiling/accumulation
#' artefact, not a cross-row coupling). To remove even that, this scorer forwards in
#' \strong{float64} and ONE ROW AT A TIME, so the \eqn{n=1} path is used uniformly and
#' a row's logit is bit-identical whether scored alone or in any batch (single-row
#' score == batch score EXACTLY 0). float64 also makes the same image in a batch of 1
#' vs a batch of 32 agree to \eqn{\sim 4\times 10^{-16}} (the frozen-BN eval-mode proof).
#'
#' \strong{GASF transform (deterministic, single-sample-safe via FROZEN bounds).} For a
#' specimen's universe-aligned per-sample rCLR vector \eqn{v} (length \eqn{p}, frozen
#' feature order): (1) min-max scale to \eqn{[-1, 1]} using FROZEN TRAINING per-feature
#' bounds \eqn{(\mathrm{lo}_j, \mathrm{hi}_j)} computed once at fit,
#' \eqn{\tilde v_j = \mathrm{clamp}(2(v_j - \mathrm{lo}_j)/(\mathrm{hi}_j -
#' \mathrm{lo}_j) - 1, -1, 1)}; (2) angular \eqn{\phi_j = \arccos(\tilde v_j)};
#' (3) \eqn{\mathrm{GASF}[i, j] = \cos(\phi_i + \phi_j) = \tilde v_i \tilde v_j -
#' \sqrt{1 - \tilde v_i^2}\sqrt{1 - \tilde v_j^2}}, a single-channel \eqn{p \times p}
#' image. Because the bounds are FROZEN (training-derived), a query is scaled by
#' training bounds -- never its own min/max -- so the image depends ONLY on the query
#' and the frozen bounds, making it single-sample. The transform is implemented in
#' python (numpy) at fit and score; a pure-R \code{.img_gasfcnn_gasf} reimplements the
#' closed form and is used only by a §7 test to independently verify the python GASF.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the self-contained
#' per-sample robust CLR over its OWN strictly-positive support (geometric-mean
#' centring on \code{v > 0}) in the FROZEN feature universe; this is exactly invariant
#' to per-specimen scaling and uses no cross-row statistic. Universe features absent
#' from a specimen carry the neutral rCLR value \code{0}. The SAME rCLR transform feeds
#' the GASF at fit and at score.
#'
#' \strong{Degenerate-neutral.} A query returns the neutral score \code{0} if fewer
#' than \code{hp$min_features} universe features overlap the query columns (a
#' column-overlap floor, batch-independent), or it has empty positive support over the
#' frozen universe on the (pre-rCLR) aligned ORIGINAL abundances
#' (\code{!any(X_use[i, ] > 0)}). A FLAT all-equal-positive composition maps to the rCLR
#' origin but is a VALID specimen: it yields a valid GASF and is scored normally (NOT
#' floored) -- the empty-support test is on the ORIGINAL abundances, NOT on
#' \code{all(z == 0)}.
#'
#' @references
#' Wang Z, Oates T. (2015) Imaging Time-Series to Improve Classification and
#' Imputation. \emph{Proceedings of the 24th International Joint Conference on
#' Artificial Intelligence (IJCAI)}, 3939-3945. arXiv:1506.00327.
#'
#' Ioffe S, Szegedy C. (2015) Batch Normalization: Accelerating Deep Network Training
#' by Reducing Internal Covariate Shift. \emph{ICML}.
#'
#' @name singlesample-img-gasfcnn
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Mirrors .proto_net_rclr / .tabicl_rclr.
.img_gasfcnn_rclr <- function(v) {
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
.img_gasfcnn_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .img_gasfcnn_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order, dropping absent universe features to the rCLR-neutral 0.
.img_gasfcnn_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R GASF (closed form) -- used ONLY by a §7 test to verify the python GASF
# ----------------------------------------------------------------------------

# Pure base-R GASF of one rCLR vector `v` (length p) under FROZEN per-feature
# bounds `lo`, `hi` (each length p). vt_j = clamp(2*(v_j-lo_j)/(hi_j-lo_j)-1,-1,1);
# GASF[i,j] = cos(phi_i + phi_j) = vt_i*vt_j - sqrt(1-vt_i^2)*sqrt(1-vt_j^2),
# evaluated via the algebraic closed form (no arccos/cos), exactly the python form.
# A zero-width bound (hi_j == lo_j) maps that coordinate to the midpoint vt_j = 0
# (matching the python guard). Returns a p x p matrix.
.img_gasfcnn_gasf <- function(v, lo, hi) {
  span <- hi - lo
  vt <- ifelse(span > 0, 2 * (v - lo) / span - 1, 0)
  vt <- pmin(pmax(vt, -1), 1)
  s <- sqrt(pmax(1 - vt^2, 0))
  outer(vt, vt) - outer(s, s)
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.img_gasfcnn_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_img_gasfcnn: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_img_gasfcnn: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_img_gasfcnn: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("channels", "epochs", "lr", "weight_decay", "min_features",
               "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_img_gasfcnn: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  channels <- hp$channels
  if (is.null(channels)) channels <- c(8L, 16L)
  if (!is.numeric(channels) || length(channels) < 1L || any(!is.finite(channels)) ||
      any(channels < 1L) || any(channels > .Machine$integer.max) ||
      any(channels != as.integer(channels))) {
    stop("fit_img_gasfcnn: hp$channels must be a vector of positive integers")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_img_gasfcnn: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_img_gasfcnn: hp$lr must be a positive finite number")
  }

  weight_decay <- hp$weight_decay
  if (is.null(weight_decay)) weight_decay <- 1e-4
  if (!is.numeric(weight_decay) || length(weight_decay) != 1L ||
      !is.finite(weight_decay) || weight_decay < 0) {
    stop("fit_img_gasfcnn: hp$weight_decay must be a non-negative finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_img_gasfcnn: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_img_gasfcnn: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_img_gasfcnn: hp$seed must be a single integer")
  }

  list(
    channels = as.integer(channels),
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    weight_decay = as.numeric(weight_decay),
    min_features = as.integer(min_features),
    device = device,
    seed = as.integer(seed)
  )
}


# ----------------------------------------------------------------------------
# Python / torch + torchvision availability + device resolution
# ----------------------------------------------------------------------------

# Clear, actionable error if reticulate, torch, or torchvision is unavailable.
# torchvision is probed because the method's dep_route is
# reticulate-torch-torchvision (a 2D image pipeline), matching the brief.
.img_gasfcnn_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_img_gasfcnn: package 'reticulate' is required to train the CNN. ",
         "Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_img_gasfcnn: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  if (!reticulate::py_module_available("torchvision")) {
    stop("fit_img_gasfcnn: Python module 'torchvision' is not available in the ",
         "active reticulate environment. Install torchvision in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu'.
.img_gasfcnn_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# CNN construction (python) -- one place builds the identical module at fit + score
# ----------------------------------------------------------------------------

# Build the python nn.Module on `device`: K blocks of
# Conv2d(cin->c,3,pad1) -> BatchNorm2d(c) -> ReLU -> MaxPool2d(2) for c in
# `channels`, then AdaptiveAvgPool2d(1) -> Flatten -> Linear(c_last -> 1). Returns
# the module (NOT seeded here; the caller seeds torch BEFORE this call so weight
# init is reproducible). `channels` is an R integer vector.
.img_gasfcnn_build_net <- function(channels, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn
  layers <- list()
  cin <- 1L
  for (c in channels) {
    cc <- as.integer(c)
    layers[[length(layers) + 1L]] <- nn$Conv2d(as.integer(cin), cc,
                                               kernel_size = 3L, padding = 1L)
    layers[[length(layers) + 1L]] <- nn$BatchNorm2d(cc)
    layers[[length(layers) + 1L]] <- nn$ReLU()
    layers[[length(layers) + 1L]] <- nn$MaxPool2d(2L)
    cin <- cc
  }
  layers[[length(layers) + 1L]] <- nn$AdaptiveAvgPool2d(1L)
  layers[[length(layers) + 1L]] <- nn$Flatten()
  layers[[length(layers) + 1L]] <- nn$Linear(as.integer(cin), 1L)
  do.call(nn$Sequential, layers)$to(device)
}

# Export a python module's state_dict to a named R list of plain numeric arrays.
# Float parameters/buffers become float64 R arrays with their torch dims preserved;
# the integer BN counter `num_batches_tracked` becomes a length-1 R integer. This
# is the FULL serialisable snapshot the scorer rebuilds from -- no python pointer.
.img_gasfcnn_export_state_dict <- function(net) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  # reticulate converts the python OrderedDict to a NAMED R list; the names are
  # the state_dict keys in their python (insertion) order.
  sd <- net$state_dict()
  key_list <- names(sd)
  out <- vector("list", length(key_list))
  names(out) <- key_list
  for (k in key_list) {
    t <- sd[[k]]
    if (grepl("num_batches_tracked", k, fixed = TRUE)) {
      out[[k]] <- as.integer(reticulate::py_to_r(t$item()))
    } else {
      arr <- reticulate::py_to_r(t$detach()$cpu()$to(torch$float64)$numpy())
      # Copy to drop the non-writable numpy view and pin dims.
      out[[k]] <- array(as.numeric(arr), dim = dim(as.array(arr)))
    }
  }
  out
}

# Load an exported R state_dict (named list of arrays / counters) into a freshly
# built python module `net` and put it in eval() float64 mode. The module is
# rebuilt with the SAME arch (channels) so every key matches by name. eval()
# freezes BatchNorm to its running stats; float64 makes the forward deterministic
# and batch-position-invariant.
.img_gasfcnn_load_state_dict <- function(net, state, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  net <- net$to(torch$float64)
  sd <- net$state_dict()
  key_list <- names(sd)
  new_sd <- reticulate::dict()
  for (k in key_list) {
    if (grepl("num_batches_tracked", k, fixed = TRUE)) {
      new_sd[[k]] <- torch$tensor(as.integer(state[[k]]),
                                  dtype = torch$long)$to(device)
    } else {
      a <- state[[k]]
      arr <- np$asarray(a, dtype = "float64")
      new_sd[[k]] <- torch$as_tensor(arr)$to(torch$float64)$to(device)
    }
  }
  net$load_state_dict(new_sd)
  net$eval()
  net
}

# Compute the FROZEN per-feature GASF bounds from the n x p training rCLR matrix.
# Per-feature min/max (column-wise) so a query feature is scaled by training-derived
# bounds, never its own. A degenerate zero-width column (constant feature) is left
# as lo == hi; the GASF guard then maps that coordinate to the midpoint.
.img_gasfcnn_bounds <- function(Z) {
  lo <- apply(Z, 2L, min)
  hi <- apply(Z, 2L, max)
  list(lo = as.numeric(lo), hi = as.numeric(hi))
}

# Build the n x 1 x p x p GASF image tensor (numpy -> torch) from the n x p rCLR
# matrix under frozen bounds, in float64. Implemented in python (numpy) per the
# brief; the closed form matches .img_gasfcnn_gasf exactly (verified in §7.1).
.img_gasfcnn_images <- function(Z, lo, hi) {
  np <- reticulate::import("numpy", delay_load = FALSE)
  Zp <- np$asarray(Z, dtype = "float64")
  lop <- np$asarray(lo, dtype = "float64")
  hip <- np$asarray(hi, dtype = "float64")
  span <- np$subtract(hip, lop)
  # vt = clamp(2*(v-lo)/span - 1, -1, 1); zero-width span -> midpoint 0.
  safe_span <- np$where(np$greater(span, 0), span, np$ones_like(span))
  vt <- np$subtract(np$multiply(2, np$divide(np$subtract(Zp, lop), safe_span)), 1)
  vt <- np$where(np$greater(span, 0), vt, np$zeros_like(vt))
  vt <- np$clip(vt, -1, 1)                                  # n x p
  s <- np$sqrt(np$clip(np$subtract(1, np$square(vt)), 0, NULL))
  # GASF[n,i,j] = vt_i*vt_j - s_i*s_j via batched outer products.
  vt_i <- np$expand_dims(vt, axis = 2L)                    # n x p x 1
  vt_j <- np$expand_dims(vt, axis = 1L)                    # n x 1 x p
  s_i <- np$expand_dims(s, axis = 2L)
  s_j <- np$expand_dims(s, axis = 1L)
  G <- np$subtract(np$multiply(vt_i, vt_j), np$multiply(s_i, s_j))  # n x p x p
  G <- np$expand_dims(G, axis = 1L)                        # n x 1 x p x p
  G
}


# ----------------------------------------------------------------------------
# Train the CNN (python) + export the state_dict to R
# ----------------------------------------------------------------------------

# Train the GASF-image CNN by full-batch BCEWithLogitsLoss (Adam, no shuffle,
# seed-before-build) on the training images, call net.eval(), and return the
# exported R-side state_dict. The torch RNG (CPU + guarded CUDA) is saved/restored
# so fit leaves the torch generators byte-unchanged; the global R .Random.seed is
# saved/restored by the caller. `Z` is the n x p rCLR'd training matrix, `y` the
# integer 0/1 labels, `bounds` the frozen GASF bounds.
.img_gasfcnn_train_export <- function(Z, y, bounds, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  nn <- torch$nn

  # ---- save the torch RNG (CPU + guarded CUDA), restore on exit -------------
  torch_state <- tryCatch(torch$random$get_rng_state(), error = function(e) NULL)
  cuda_ok <- tryCatch(isTRUE(torch$cuda$is_available()), error = function(e) FALSE)
  cuda_state <- if (cuda_ok) {
    tryCatch(torch$cuda$get_rng_state_all(), error = function(e) NULL)
  } else {
    NULL
  }
  # Force cuDNN into DETERMINISTIC mode for the duration of the fit so the GPU
  # convolution path is bit-reproducible (the default cuDNN conv backward is
  # nondeterministic; on CPU these flags are inert). Without this, two
  # seed-matched GPU fits diverge by ~1e-5. The previous flag values are saved and
  # restored on exit so the global torch backend config is left unchanged.
  old_det <- tryCatch(torch$backends$cudnn$deterministic, error = function(e) NULL)
  old_bench <- tryCatch(torch$backends$cudnn$benchmark, error = function(e) NULL)
  on.exit({
    if (!is.null(torch_state)) {
      tryCatch(torch$random$set_rng_state(torch_state), error = function(e) NULL)
    }
    if (!is.null(cuda_state)) {
      tryCatch(torch$cuda$set_rng_state_all(cuda_state), error = function(e) NULL)
    }
    if (!is.null(old_det)) {
      tryCatch(torch$backends$cudnn$deterministic <- old_det,
               error = function(e) NULL)
    }
    if (!is.null(old_bench)) {
      tryCatch(torch$backends$cudnn$benchmark <- old_bench,
               error = function(e) NULL)
    }
  }, add = TRUE)
  tryCatch(torch$backends$cudnn$deterministic <- TRUE, error = function(e) NULL)
  tryCatch(torch$backends$cudnn$benchmark <- FALSE, error = function(e) NULL)

  # Seed torch BEFORE building the net so weight init is reproducible (two
  # seed-matched fits initialise identically). Full-batch, no shuffle -> no other
  # torch randomness, so the whole fit is deterministic.
  torch$manual_seed(as.integer(hp$seed))

  net <- .img_gasfcnn_build_net(hp$channels, device)        # float32 training net

  # Training images (float32 for training speed) + labels.
  G <- .img_gasfcnn_images(Z, bounds$lo, bounds$hi)         # n x 1 x p x p float64
  X_py <- torch$as_tensor(G)$to(torch$float32)$to(device)
  y_py <- torch$as_tensor(np$asarray(as.numeric(y), dtype = "float32"))$to(device)

  opt <- torch$optim$Adam(net$parameters(), lr = hp$lr,
                          weight_decay = hp$weight_decay)
  loss_fn <- nn$BCEWithLogitsLoss()

  net$train()
  for (ep in seq_len(hp$epochs)) {
    opt$zero_grad()
    logit <- net(X_py)$reshape(list(-1L))                   # n
    loss <- loss_fn(logit, y_py)
    loss$backward()
    opt$step()
  }
  net$eval()

  list(
    state_dict = .img_gasfcnn_export_state_dict(net),
    channels = as.integer(hp$channels)
  )
}


# ----------------------------------------------------------------------------
# fit_img_gasfcnn
# ----------------------------------------------------------------------------

#' @title Fit the GASF-image CNN single-sample discriminator
#'
#' @description
#' Maps the TRAINING matrix to the per-sample robust CLR over the frozen feature
#' universe (\code{colnames(X_train)}), computes the FROZEN per-feature GASF bounds,
#' renders each training specimen to a GASF image, and trains a small
#' Conv-BN-ReLU-MaxPool CNN end-to-end by full-batch \code{BCEWithLogitsLoss}
#' (Adam, no shuffle, seed-before-build) via reticulate-python torch + torchvision.
#' After training the net is put in \code{eval()} mode (frozen BatchNorm) and its
#' FULL \code{state_dict} is EXPORTED to R as named float64 numeric arrays (incl. BN
#' running stats); the python module is DISCARDED. The fitted model holds NO external
#' pointer -- score rebuilds the CNN from the R state.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{channels}
#'   (positive-integer vector of per-block conv channel widths; default
#'   \code{c(8L, 16L)} -- two blocks), \code{epochs} (full-batch training epochs,
#'   positive integer; default \code{200L}), \code{lr} (positive Adam learning rate;
#'   default \code{1e-3}), \code{weight_decay} (non-negative Adam L2; default
#'   \code{1e-4}), \code{min_features} (feature-overlap floor at scoring, positive
#'   integer; default \code{3L}), \code{device} (\code{"cpu"} (default),
#'   \code{"cuda"}, or \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU
#'   with no GPU), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{img_gasfcnn_model}: a list with
#'   \code{feature_universe}, \code{state_dict} (named float64 arrays + BN counters),
#'   \code{channels}, \code{gasf_lo}, \code{gasf_hi} (frozen per-feature bounds),
#'   \code{device} (resolved), \code{seed}, and \code{hp}. No python pointer.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 80; p <- 16; k <- 6
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' model <- fit_img_gasfcnn(X, y)
#' score_img_gasfcnn(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Wang Z, Oates T. (2015) Imaging Time-Series to Improve Classification and
#' Imputation. \emph{IJCAI}; arXiv:1506.00327.
#'
#' @export
fit_img_gasfcnn <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_img_gasfcnn", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_img_gasfcnn", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_img_gasfcnn")
  hp <- .img_gasfcnn_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_img_gasfcnn: X_train must contain at least hp$min_features features")
  }

  # Save the global R .Random.seed BEFORE any reticulate import: reticulate's first
  # python init seeds R's RNG. Restoring on exit (incl. the no-prior-seed case)
  # keeps fit globally RNG-neutral. The fit draws no R randomness itself.
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

  .img_gasfcnn_require_torch()

  feat_order <- colnames(X_train)
  Z <- .img_gasfcnn_rclr_matrix(X_train)            # n x p rCLR over full universe
  bounds <- .img_gasfcnn_bounds(Z)                  # frozen per-feature GASF bounds

  device <- .img_gasfcnn_resolve_device(hp$device)
  trained <- .img_gasfcnn_train_export(Z, y, bounds, hp, device)

  model <- list(
    feature_universe = feat_order,
    state_dict = trained$state_dict,
    channels = trained$channels,
    gasf_lo = bounds$lo,
    gasf_hi = bounds$hi,
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "img_gasfcnn_model"
  model
}


# ----------------------------------------------------------------------------
# score_img_gasfcnn
# ----------------------------------------------------------------------------

#' @title Score the GASF-image CNN single-sample discriminator (forced row-by-row)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_img_gasfcnn}}. The python CNN is REBUILT from the exported R-side
#' \code{state_dict} (\code{load_state_dict}, \code{eval()} -> FROZEN BatchNorm,
#' float64), then each query is mapped to the per-sample robust CLR over the frozen
#' feature universe (absent universe features carry the neutral rCLR value \code{0}),
#' rendered to its GASF image under the FROZEN bounds, and forwarded through the CNN
#' ONE ROW AT A TIME. The score is the eval-mode logit; larger = more case-like.
#' Forcing the \eqn{n=1} forward path uniformly (float64) makes a row's logit
#' bit-identical whether scored alone or in any batch.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' (pre-rCLR) aligned abundances (\code{!any(X_use[i, ] > 0)}), return the neutral
#' score \code{0}. A FLAT all-equal-positive composition maps to the rCLR origin but
#' is a VALID specimen and is scored normally. The score of a row depends only on that
#' row and the frozen model, and is exactly invariant to per-specimen positive scaling.
#'
#' @param model An \code{img_gasfcnn_model} object from \code{\link{fit_img_gasfcnn}}.
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
#' model <- fit_img_gasfcnn(X, y)
#' score_img_gasfcnn(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Wang Z, Oates T. (2015) Imaging Time-Series to Improve Classification and
#' Imputation. \emph{IJCAI}; arXiv:1506.00327.
#'
#' @export
score_img_gasfcnn <- function(model, X, meta = NULL) {
  if (!inherits(model, "img_gasfcnn_model")) {
    stop("score_img_gasfcnn: model must have class img_gasfcnn_model")
  }
  X <- .reo_check_matrix(X, "score_img_gasfcnn", "X")
  .reo_check_meta(meta, nrow(X), "score_img_gasfcnn", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  .img_gasfcnn_require_torch()
  torch <- reticulate::import("torch", delay_load = FALSE)
  device <- .img_gasfcnn_resolve_device(model$hp$device)

  X_use <- .img_gasfcnn_align(X, feat_order)
  lo <- model$gasf_lo
  hi <- model$gasf_hi

  # Rebuild the identical CNN, load the exported weights, eval() (frozen BN),
  # float64. The torch RNG is saved/restored around the whole loop so scoring
  # leaves the torch generators byte-unchanged (eval-mode forward consumes no RNG,
  # but the guard makes RNG-safety robust). No python pointer survives the call.
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

  net <- .img_gasfcnn_build_net(model$channels, device)
  net <- .img_gasfcnn_load_state_dict(net, model$state_dict, device)

  no_grad <- torch$no_grad()
  no_grad$`__enter__`()
  on.exit(tryCatch(no_grad$`__exit__`(NULL, NULL, NULL), error = function(e) NULL),
          add = TRUE)

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test on the ORIGINAL abundances, NOT on all(z == 0): a FLAT composition has
    # an all-zero rCLR but is a valid specimen that MUST be imaged and scored.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    z <- .img_gasfcnn_rclr(X_use[i, ])
    # One row at a time: the n = 1 float64 forward path is used uniformly, so the
    # logit is bit-identical alone or in any batch.
    G <- .img_gasfcnn_images(matrix(z, nrow = 1L), lo, hi)   # 1 x 1 x p x p
    img <- torch$as_tensor(G)$to(torch$float64)$to(device)
    logit <- net(img)$reshape(list(-1L))
    out[i] <- as.numeric(reticulate::py_to_r(logit$item()))
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_img_gasfcnn: scorer produced non-finite or wrong-length output")
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
# must be identical before and after scoring). Two seed-matched fits produce
# identical digests.
.img_gasfcnn_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    state_dict = model$state_dict,
    channels = model$channels,
    gasf_lo = model$gasf_lo,
    gasf_hi = model$gasf_hi,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
