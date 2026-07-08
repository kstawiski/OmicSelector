#' @title CoDaCoRe stagewise log-ratio-balance single-sample discriminator
#'   (torch relaxation at fit, frozen discrete balances, PURE-R score)
#'
#' @description
#' Single-sample within-cohort discriminator (family O, compositional) implementing
#' CoDaCoRe (Gordon-Rodriguez, Quinn & Cunningham 2021, "Learning sparse log-ratios
#' for high-throughput sequencing data"): a STAGEWISE (boosting) sequence of sparse
#' log-ratio BALANCES, each learned by a CONTINUOUS RELAXATION, scored through a
#' FROZEN logistic head. Larger score = more case-like.
#'
#' \strong{The CoDaCoRe algorithm.} For \eqn{b = 1, \dots, B} (B = \code{hp$max_balances})
#' a stage introduces a learnable per-feature assignment over three states --
#' numerator (+), denominator (-), and neither -- via a temperature-\code{tau}
#' softmax of a logit matrix \eqn{W_b \in \mathbb{R}^{p \times 3}}. The SOFT balance
#' of specimen \eqn{i} is the relaxed log-contrast
#' \deqn{SB_b(i) = \frac{\sum_j s^{+}_{j}\,\log x_{ij}}{\sum_j s^{+}_{j}}
#'   - \frac{\sum_j s^{-}_{j}\,\log x_{ij}}{\sum_j s^{-}_{j}},}
#' with \eqn{s^{+}, s^{-}} the soft (+)/(-) memberships. \eqn{W_b} and a scalar stage
#' slope are trained by Adam to minimise the binary cross-entropy of the boosting
#' RESIDUAL (\eqn{y - \sigma(\text{running score})}). The trained soft memberships are
#' then DISCRETIZED to two disjoint feature sets \eqn{A_b} (numerator, argmax \eqn{=+})
#' and \eqn{B_b} (denominator, argmax \eqn{=-}), with explicit floors \eqn{|A_b| \ge 1}
#' and \eqn{|B_b| \ge 1} and a deterministic tie-break; \eqn{A_b, B_b} are FROZEN. The
#' discrete SBP/ILR balance is
#' \deqn{\mathrm{bal}_b(x) = \sqrt{\frac{|A_b||B_b|}{|A_b|+|B_b|}}\,
#'   \bigl(\overline{\log x_{A_b}} - \overline{\log x_{B_b}}\bigr),}
#' computed over the specimen's OWN strictly-positive support within \eqn{A_b / B_b}.
#' The running score is re-fit (a base-R IRLS ridge logistic on the discrete-balance
#' matrix), the residual recomputed, and the next stage proceeds; a stage is ACCEPTED
#' only if it lowers the training-only deviance, otherwise stagewise growth stops.
#'
#' \strong{Distinct from selbal.} \code{bal-selbal} learns ONE balance by GREEDY
#' forward feature selection; \code{coda-codacore} learns a SEQUENCE of balances by
#' CONTINUOUS-RELAXATION stagewise boosting -- the relaxation + discretization +
#' boosting residual are what make it CoDaCoRe, not selbal.
#'
#' \strong{DEP-ROUTE PIVOT: codacore-tensorflow -> reticulate-torch-plus-R.} The
#' prespecified backend was TensorFlow/Keras, but TensorFlow has NO py3.14 wheel and
#' the R \pkg{codacore} package is not installed. The CoDaCoRe ALGORITHM is preserved
#' EXACTLY; only the autodiff backend changes TF/Keras -> the already-installed torch
#' (the relaxation is a tiny length-\eqn{3p} logit model). Torch is confined to FIT:
#' the relaxation is trained, the balances DISCRETIZED + FROZEN, and the discrete sets
#' + logistic head EXPORTED to R; the python module is DISCARDED. The fitted model
#' holds NO external pointer and NO python dependency at score time.
#'
#' \strong{The score is PURE BASE R.} For a query row \eqn{x}: align to the frozen
#' universe; for each frozen balance \eqn{b} compute \eqn{\mathrm{bal}_b(x)} over
#' \eqn{x}'s own strictly-positive support in \eqn{A_b / B_b}; if a side has NO positive
#' feature for this specimen that balance contributes its neutral value \code{0}. The
#' score is \eqn{\mathrm{intercept} + \sum_b w_b\,\mathrm{bal}_b(x)}. Each balance is a
#' per-specimen difference of mean-logs over that specimen's own features -- it uses no
#' cross-sample statistic and no python -- so a row's score is bit-identical scored alone
#' or in any batch (single-row == batch EXACTLY 0), and is exactly invariant to per-sample
#' positive scaling (a positive rescale \eqn{c} cancels:
#' \eqn{\overline{\log(c\,x_A)} - \overline{\log(c\,x_B)} = \overline{\log x_A} -
#' \overline{\log x_B}}).
#'
#' \strong{Reproducibility: the relaxation is fit on CPU by default.} The DISCRETE
#' balances (the \eqn{A_b, B_b} sets) must be identical across two seed-matched fits for
#' the determinism gate and the \code{model_digest}. GPU reductions are nondeterministic
#' and could flip a near-threshold feature's argmax membership, so the relaxation is
#' trained on CPU by DEFAULT (the per-stage model is tiny, so CPU is fast and EXACTLY
#' reproducible), seeded before build, full-batch with no shuffle, with a deterministic
#' argmax discretization and an explicit lowest-index tie-break. (Tests may probe
#' \code{device = "cuda"}/\code{"auto"} for the FIT; the SCORE is pure-R and therefore
#' device-irrelevant.)
#'
#' \strong{Degenerate-neutral / rCLR-origin.} A query returns the neutral score \code{0}
#' only if fewer than \code{hp$min_features} universe features overlap the query columns
#' (a column-overlap floor, batch-independent) OR it has empty positive support over the
#' frozen universe on the (pre-log) aligned ORIGINAL abundances
#' (\code{!any(X_use[i, ] > 0)}). A FLAT all-equal-positive composition gives every
#' balance the value \code{0} (\eqn{\overline{\log x_A} - \overline{\log x_B} = 0}) but is
#' a VALID specimen, scored to the intercept (NOT floored). The score consumes NO RNG.
#'
#' @references
#' Gordon-Rodriguez E, Quinn TP, Cunningham JP. (2021) Learning sparse log-ratios for
#' high-throughput sequencing data. \emph{Bioinformatics} 38(1):157-163.
#'
#' Egozcue JJ, Pawlowsky-Glahn V. (2005) Groups of parts and their balances in
#' compositional data analysis. \emph{Mathematical Geology} 37(7):795-828.
#'
#' @name singlesample-coda-codacore
NULL


# ----------------------------------------------------------------------------
# Single-sample input transform (per-sample log over the frozen universe)
# ----------------------------------------------------------------------------

# Align an n x q input matrix `X` (its columns named) to the frozen feature
# universe order. Absent universe features are filled with 0 (a sentinel for
# "absent"; the balance helpers exclude non-positive entries from a specimen's
# support, so a 0 here means the feature is simply not in that specimen's support).
.coda_codacore_align <- function(X, feat_order) {
  present <- intersect(feat_order, colnames(X))
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]
  X_use
}


# ----------------------------------------------------------------------------
# Pure-R discrete ILR balance over a specimen's own positive support
# ----------------------------------------------------------------------------

# Discrete SBP/ILR balance for ONE specimen over the frozen numerator (A) /
# denominator (B) index sets, using only features strictly POSITIVE in this
# specimen. `xrow` is the universe-aligned abundance vector (length p), `a_idx`
# and `b_idx` are 1-based positions into the frozen universe. The balance is
#   sqrt(r s / (r + s)) * (mean_{A+} log x - mean_{B+} log x)
# where r, s are the per-specimen counts of present-positive A / B features. If
# EITHER side has no positive feature for this specimen the balance is the neutral
# value 0 (documented single-sample behaviour: a missing side contributes nothing).
.coda_codacore_balance_one <- function(xrow, a_idx, b_idx) {
  va <- xrow[a_idx]; va <- va[is.finite(va) & va > 0]
  vb <- xrow[b_idx]; vb <- vb[is.finite(vb) & vb > 0]
  r <- length(va); s <- length(vb)
  if (r < 1L || s < 1L) return(0)
  k <- sqrt((r * s) / (r + s))
  k * (mean(log(va)) - mean(log(vb)))
}

# n x K matrix of discrete balances for the n x p universe-aligned matrix `X_use`
# over the K frozen balances `balances` (each a list(A = int vec, B = int vec)).
.coda_codacore_balance_matrix <- function(X_use, balances) {
  K <- length(balances)
  n <- nrow(X_use)
  B <- matrix(0, nrow = n, ncol = K)
  if (K == 0L) return(B)
  for (i in seq_len(n)) {
    xi <- X_use[i, ]
    for (b in seq_len(K)) {
      B[i, b] <- .coda_codacore_balance_one(xi, balances[[b]]$A, balances[[b]]$B)
    }
  }
  B
}


# ----------------------------------------------------------------------------
# Base-R ridge IRLS logistic head (frozen) over the discrete-balance matrix
# ----------------------------------------------------------------------------

# FROZEN logistic head over the K discrete balances, fit by IRLS with a small L2
# ridge on the WEIGHTS (intercept unpenalized) for numerical stability. Returns
# list(intercept, weights). Deterministic, no RNG. `lambda` is a small ridge.
.coda_codacore_ridge_logit <- function(Bmat, y, lambda, max_iter = 200L,
                                       tol = 1e-10) {
  Xa <- cbind(1, Bmat)                         # n x (K+1), col 1 = intercept
  D1 <- ncol(Xa)
  pen <- rep(lambda, D1); pen[1L] <- 0         # do not penalize the intercept
  b <- rep(0, D1)
  for (it in seq_len(max_iter)) {
    eta <- as.vector(Xa %*% b)
    mu <- 1 / (1 + exp(-eta))
    w <- mu * (1 - mu)
    w[w < 1e-9] <- 1e-9                         # IRLS weight floor (stability)
    zwork <- eta + (y - mu) / w
    XtW <- t(Xa * w)
    H <- XtW %*% Xa
    diag(H) <- diag(H) + pen
    g <- XtW %*% zwork
    bn <- tryCatch(solve(H, g), error = function(e) {
      solve(H + diag(1e-8, D1), g)             # tiny jitter if ill-conditioned
    })
    bn <- as.vector(bn)
    if (max(abs(bn - b)) < tol) { b <- bn; break }
    b <- bn
  }
  list(intercept = b[1L], weights = b[-1L])
}

# Mean binomial deviance of the frozen head on (Bmat, y). Lower = better fit.
# Used as the stagewise ACCEPT criterion (training-only).
.coda_codacore_deviance <- function(Bmat, y, head) {
  eta <- head$intercept + as.vector(Bmat %*% head$weights)
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-12
  mu <- pmin(pmax(mu, eps), 1 - eps)
  -2 * mean(y * log(mu) + (1 - y) * log(1 - mu))
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.coda_codacore_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_coda_codacore: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_coda_codacore: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_coda_codacore: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("max_balances", "epochs", "lr", "tau", "ridge", "min_features",
               "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_coda_codacore: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  max_balances <- hp$max_balances
  if (is.null(max_balances)) max_balances <- 3L
  if (!is.numeric(max_balances) || length(max_balances) != 1L ||
      !is.finite(max_balances) || max_balances < 1L ||
      max_balances > .Machine$integer.max ||
      max_balances != as.integer(max_balances)) {
    stop("fit_coda_codacore: hp$max_balances must be a positive integer")
  }

  epochs <- hp$epochs
  if (is.null(epochs)) epochs <- 200L
  if (!is.numeric(epochs) || length(epochs) != 1L || !is.finite(epochs) ||
      epochs < 1L || epochs > .Machine$integer.max ||
      epochs != as.integer(epochs)) {
    stop("fit_coda_codacore: hp$epochs must be a positive integer")
  }

  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-2
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_coda_codacore: hp$lr must be a positive finite number")
  }

  tau <- hp$tau
  if (is.null(tau)) tau <- 1.0
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("fit_coda_codacore: hp$tau must be a positive finite number")
  }

  ridge <- hp$ridge
  if (is.null(ridge)) ridge <- 1e-4
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) ||
      ridge < 0) {
    stop("fit_coda_codacore: hp$ridge must be a non-negative finite number")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_coda_codacore: hp$min_features must be a positive integer")
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_coda_codacore: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_coda_codacore: hp$seed must be a single integer")
  }

  list(
    max_balances = as.integer(max_balances),
    epochs = as.integer(epochs),
    lr = as.numeric(lr),
    tau = as.numeric(tau),
    ridge = as.numeric(ridge),
    min_features = as.integer(min_features),
    device = device,
    seed = as.integer(seed)
  )
}


# ----------------------------------------------------------------------------
# Python / torch availability + device resolution
# ----------------------------------------------------------------------------

.coda_codacore_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_coda_codacore: package 'reticulate' is required to train the ",
         "relaxation. Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_coda_codacore: Python module 'torch' is not available in the ",
         "active reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a 'cuda'
# request) fall back to 'cpu' when CUDA is unavailable. The DEFAULT is 'cpu' --
# the relaxation is fit on CPU for EXACT discretization reproducibility.
.coda_codacore_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# One CoDaCoRe stage: train the relaxation in torch, return the SOFT memberships
# ----------------------------------------------------------------------------

# Deterministic, RNG-free symmetry-breaking init for one stage's relaxation logits.
# The relaxation has an exact +/- symmetry at W = 0 (softmax(0) is uniform 1/3, the
# soft balance is identically 0, and the numerator/denominator columns are
# interchangeable, so the gradient at the symmetric point is zero -- a saddle the
# optimiser never leaves). CoDaCoRe breaks this toward each feature's marginal
# log-ratio association with the target; we do the same DETERMINISTICALLY (no RNG):
# `w0_j = corr(clr-log_j, resid)` -- features positively associated with the
# residual get a higher numerator logit, negatively associated a higher denominator
# logit, and the neither column stays at the baseline 0. A constant `init_scale`
# sets the strength. This is computed in R (deterministic), so two seed-matched CPU
# fits init identically -> EXACTLY reproducible discretization. `Lmat`/`Mmask`/
# `resid` as in .coda_codacore_train_stage.
.coda_codacore_stage_init <- function(Lmat, Mmask, resid, init_scale = 2.0) {
  # Per-sample CLR over each specimen's OWN positive support (mean-log centred),
  # matching the balance geometry; non-positive entries contribute 0.
  n <- nrow(Lmat); p <- ncol(Lmat)
  Lc <- matrix(0, n, p)
  for (i in seq_len(n)) {
    m <- Mmask[i, ] > 0
    if (any(m)) {
      lv <- Lmat[i, m]
      Lc[i, m] <- lv - mean(lv)
    }
  }
  rsd <- resid - mean(resid)
  denom_r <- sqrt(sum(rsd^2))
  w0 <- vapply(seq_len(p), function(j) {
    cj <- Lc[, j] - mean(Lc[, j])
    dj <- sqrt(sum(cj^2))
    if (dj < 1e-12 || denom_r < 1e-12) return(0)
    sum(cj * rsd) / (dj * denom_r)               # Pearson corr, in [-1, 1]
  }, numeric(1L))
  w0 * init_scale
}

# Train one stage's continuous relaxation in torch and return the trained soft
# (+)/(-) memberships as plain R numeric vectors. `Lmat` is the n x p per-sample
# LOG-abundance matrix (already log-transformed over the universe, with absent /
# non-positive entries pre-masked to 0 contribution via `Mmask`); `Mmask` is the
# n x p 0/1 mask of strictly-positive entries (so the soft balance sums only over a
# specimen's own support); `resid` is the boosting residual (y - sigmoid(running
# score)) the stage regresses; `hp`/`device` the config. The relaxation is a logit
# matrix W in R^{p x 3} (numerator / denominator / neither) -> softmax(W / tau)
# memberships s+, s-, s0; the soft balance per sample is
#   ( sum_j s+_j m_ij L_ij / sum_j s+_j m_ij ) - ( sum_j s-_j m_ij L_ij / sum_j s-_j m_ij )
# (the per-sample support mask m_i restricts each mean-log to that specimen's
# positive support, matching the discrete balance). A scalar slope + intercept map
# the soft balance to a logit; Adam minimises the BCE-with-logits of `resid` (the
# boosting target). Full-batch, no shuffle, seeded before build -> deterministic on
# CPU. The torch CPU/CUDA RNG is saved/restored by the CALLER (fit-level guard).
.coda_codacore_train_stage <- function(Lmat, Mmask, resid, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)

  # Seed before building the stage parameters so two seed-matched CPU fits init the
  # generators identically. The W init is the DETERMINISTIC, RNG-free
  # correlation-based symmetry break below (not a random init); seeding keeps the
  # stage robust to any torch-side randomness and matches the trained-torch template.
  torch$manual_seed(as.integer(hp$seed))

  n <- nrow(Lmat); p <- ncol(Lmat)
  L_py <- torch$as_tensor(np$asarray(Lmat, dtype = "float64"))$to(torch$float64)$to(device)
  M_py <- torch$as_tensor(np$asarray(Mmask, dtype = "float64"))$to(torch$float64)$to(device)
  r_py <- torch$as_tensor(np$asarray(matrix(resid, ncol = 1L), dtype = "float64"))$to(torch$float64)$to(device)

  # Relaxation logits W (p x 3): col 1 = numerator(+), col 2 = denominator(-),
  # col 3 = neither. DETERMINISTIC symmetry-breaking init: numerator logit =
  # +w0_j, denominator logit = -w0_j, neither = 0, with w0_j the feature's
  # residual correlation (see .coda_codacore_stage_init). This escapes the W = 0
  # saddle along the data-driven direction, so Adam refines a meaningful +/- split,
  # reproducibly. Scalar slope + bias map the soft balance to a logit for the BCE.
  w0 <- .coda_codacore_stage_init(Lmat, Mmask, resid)
  W_init <- cbind(w0, -w0, rep(0, p))            # p x 3 (num, den, neither)
  W <- torch$as_tensor(np$asarray(W_init, dtype = "float64"))$to(torch$float64)$to(device)
  W$requires_grad_(TRUE)
  slope <- torch$zeros(list(1L), dtype = torch$float64,
                       requires_grad = TRUE)$to(device)
  bias <- torch$zeros(list(1L), dtype = torch$float64,
                      requires_grad = TRUE)$to(device)
  W$requires_grad_(TRUE); slope$requires_grad_(TRUE); bias$requires_grad_(TRUE)

  params <- list(W, slope, bias)
  opt <- torch$optim$Adam(params, lr = hp$lr)
  bce <- torch$nn$functional$binary_cross_entropy_with_logits
  tau <- as.numeric(hp$tau)
  eps <- 1e-12

  for (ep in seq_len(hp$epochs)) {
    opt$zero_grad()
    sm <- torch$softmax(W$div(tau), dim = 1L)                  # p x 3 memberships
    s_plus <- sm$index_select(1L, torch$as_tensor(0L)$long()$to(device))$squeeze(1L)   # p
    s_minus <- sm$index_select(1L, torch$as_tensor(1L)$long()$to(device))$squeeze(1L)  # p
    # Per-sample masked weighted-mean-log over the soft + / - memberships. M_py
    # restricts each sample's sum to its own positive support.
    wp <- M_py$mul(s_plus$unsqueeze(0L))                       # n x p  (m_ij s+_j)
    wm <- M_py$mul(s_minus$unsqueeze(0L))                      # n x p  (m_ij s-_j)
    num_plus <- wp$mul(L_py)$sum(1L, keepdim = TRUE)
    den_plus <- wp$sum(1L, keepdim = TRUE)$add(eps)
    num_minus <- wm$mul(L_py)$sum(1L, keepdim = TRUE)
    den_minus <- wm$sum(1L, keepdim = TRUE)$add(eps)
    sb <- num_plus$div(den_plus)$sub(num_minus$div(den_minus)) # n x 1 soft balance
    logit <- sb$mul(slope)$add(bias)                           # n x 1 stage logit
    loss <- bce(logit, r_py)
    loss$backward()
    opt$step()
  }

  # Export the trained soft memberships (numerator / denominator / neither).
  sm_final <- torch$softmax(W$detach()$div(tau), dim = 1L)
  sm_r <- reticulate::py_to_r(sm_final$cpu()$to(torch$float64)$numpy())
  storage.mode(sm_r) <- "double"
  dim(sm_r) <- c(p, 3L)
  sm_r
}


# Discretize one stage's soft memberships (p x 3: +, -, neither) to disjoint sets
# A (numerator) / B (denominator). A feature joins the side of its largest soft
# membership (argmax over the 3 columns); ties resolve to the LOWEST column index
# (+ before -, - before neither) for a deterministic, device-independent split.
# Floors |A| >= 1 and |B| >= 1: if a side is empty, the single most-confident
# feature for that side (largest soft membership in that column, lowest-index
# tie-break) is FORCED into it. Returns list(A = int vec, B = int vec) of 1-based
# universe positions, or NULL if a non-degenerate disjoint split is impossible
# (p < 2, or the forced floors would overlap).
.coda_codacore_discretize <- function(sm) {
  p <- nrow(sm)
  if (p < 2L) return(NULL)
  # argmax over columns with a lowest-index tie-break (which.max already returns
  # the first maximum, i.e. the lowest column index).
  assign <- apply(sm, 1L, which.max)                 # 1 = +, 2 = -, 3 = neither
  A <- which(assign == 1L)
  B <- which(assign == 2L)

  if (length(A) == 0L) {
    # force the most-confident numerator feature (largest col-1 membership);
    # ties -> lowest index. Exclude any feature already forced to B below.
    A <- which.max(sm[, 1L])
  }
  if (length(B) == 0L) {
    cand <- setdiff(seq_len(p), A)                   # keep A and B disjoint
    if (length(cand) == 0L) return(NULL)
    B <- cand[which.max(sm[cand, 2L])]
  }
  # A feature forced into A could in principle also be the argmax for B; the
  # setdiff above guarantees disjointness. Final guard:
  if (length(intersect(A, B)) > 0L) return(NULL)
  list(A = as.integer(A), B = as.integer(B))
}


# ----------------------------------------------------------------------------
# Stagewise CoDaCoRe boosting (CPU-default) -> frozen discrete balances + head
# ----------------------------------------------------------------------------

# Run the full CoDaCoRe stagewise boosting and return the FROZEN discrete balances
# (list of list(A, B)) and the final frozen logistic head. `X_use` is the n x p
# universe-aligned ORIGINAL abundance matrix; `y` the 0/1 labels. Per stage: train
# the relaxation on the boosting residual, discretize to A/B, recompute the
# discrete-balance matrix INCLUDING the new balance, refit the base-R IRLS ridge
# head, and ACCEPT the new balance only if it strictly lowers the training
# deviance (else stop). The running score for the residual uses the SAME pure-R
# discrete balances + head the scorer uses, so fit and score share one code path.
.coda_codacore_boost <- function(X_use, y, hp, device) {
  # Save/restore the torch RNG (CPU + guarded CUDA) INSIDE this helper so the
  # restore (a reticulate call that can re-create R's .Random.seed as a side
  # effect) COMPLETES before control returns to fit_coda_codacore -- the fit-level
  # R-seed cleanup then runs LAST with nothing re-creating the seed afterward (the
  # proto-net pattern). The stages seed torch internally; this leaves the caller's
  # torch generators byte-unchanged.
  torch <- reticulate::import("torch", delay_load = FALSE)
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

  n <- nrow(X_use); p <- ncol(X_use)
  # Per-sample LOG over the positive support; non-positive entries masked out.
  Mmask <- (X_use > 0) * 1
  Lmat <- matrix(0, n, p)
  pos <- X_use > 0
  Lmat[pos] <- log(X_use[pos])

  balances <- list()
  head <- list(intercept = 0, weights = numeric(0))
  # Initial running deviance: intercept-only logistic (no balances yet).
  Bmat <- matrix(0, n, 0)
  head0 <- .coda_codacore_ridge_logit(Bmat, y, hp$ridge)
  cur_dev <- .coda_codacore_deviance(Bmat, y, head0)
  head <- head0

  for (b in seq_len(hp$max_balances)) {
    # Boosting residual = y - sigmoid(running score) under the CURRENT head.
    if (length(balances) == 0L) {
      eta <- rep(head$intercept, n)
    } else {
      Bcur <- .coda_codacore_balance_matrix(X_use, balances)
      eta <- head$intercept + as.vector(Bcur %*% head$weights)
    }
    resid <- y - 1 / (1 + exp(-eta))

    sm <- .coda_codacore_train_stage(Lmat, Mmask, resid, hp, device)
    new_bal <- .coda_codacore_discretize(sm)
    if (is.null(new_bal)) break

    cand_balances <- c(balances, list(new_bal))
    Bcand <- .coda_codacore_balance_matrix(X_use, cand_balances)
    head_cand <- .coda_codacore_ridge_logit(Bcand, y, hp$ridge)
    dev_cand <- .coda_codacore_deviance(Bcand, y, head_cand)

    # ACCEPT only if the new balance strictly improves training deviance.
    if (dev_cand < cur_dev - 1e-8) {
      balances <- cand_balances
      head <- head_cand
      cur_dev <- dev_cand
    } else {
      break
    }
  }

  list(balances = balances, head = head, deviance = cur_dev)
}


# ----------------------------------------------------------------------------
# fit_coda_codacore
# ----------------------------------------------------------------------------

#' @title Fit the CoDaCoRe stagewise log-ratio-balance discriminator
#'
#' @description
#' Trains the CoDaCoRe continuous relaxation stagewise (boosting) via
#' reticulate-python torch (CPU by default for exact reproducibility), DISCRETIZES
#' each balance to two disjoint feature sets, FREEZES them, fits a base-R IRLS ridge
#' logistic head over the discrete balances, EXPORTS the discrete sets + head to R,
#' and DISCARDS the python module. The fitted model holds no external pointer.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{max_balances}
#'   (max stagewise balances, positive integer; default \code{3L}), \code{epochs}
#'   (Adam epochs per stage, positive integer; default \code{200L}), \code{lr}
#'   (positive Adam learning rate; default \code{1e-2}), \code{tau} (positive
#'   softmax relaxation temperature; default \code{1.0}), \code{ridge}
#'   (non-negative L2 ridge on the logistic head; default \code{1e-4}),
#'   \code{min_features} (feature-overlap floor at scoring, positive integer;
#'   default \code{3L}), \code{device} (\code{"cpu"} (default), \code{"cuda"}, or
#'   \code{"auto"}; \code{"cuda"}/\code{"auto"} fall back to CPU with no GPU -- but
#'   note the default \code{"cpu"} is what guarantees an EXACTLY reproducible
#'   discretization), and \code{seed} (integer; default \code{42L}).
#'
#' @return Object of class \code{coda_codacore_model}: a list with
#'   \code{feature_universe}, \code{balances} (a list of \code{list(A, B)} 1-based
#'   universe-position index sets per frozen discrete balance), \code{intercept},
#'   \code{weights} (the frozen logistic head over the balances), \code{n_balances},
#'   \code{device} (resolved), \code{seed}, and \code{hp}. No python pointer.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 16; k <- 3
#' X <- matrix(stats::rgamma(n * p, shape = 30, rate = 2), nrow = n,
#'             dimnames = list(NULL, paste0("f", sprintf("%02d", seq_len(p)))))
#' y <- rep(c(0, 1), each = n / 2)
#' X[y == 1, 1:k] <- X[y == 1, 1:k] * 1.8
#' X[y == 1, (p - k + 1):p] <- X[y == 1, (p - k + 1):p] / 1.8
#' model <- fit_coda_codacore(X, y)
#' score_coda_codacore(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Gordon-Rodriguez E, Quinn TP, Cunningham JP. (2021) Learning sparse log-ratios
#' for high-throughput sequencing data. \emph{Bioinformatics} 38(1):157-163.
#'
#' @export
fit_coda_codacore <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_coda_codacore", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_coda_codacore", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_coda_codacore")
  hp <- .coda_codacore_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_coda_codacore: X_train must contain at least hp$min_features features")
  }
  if (ncol(X_train) < 2L) {
    stop("fit_coda_codacore: X_train must contain at least two features for a balance")
  }

  # Save the global R .Random.seed BEFORE any reticulate import (reticulate's first
  # python init seeds R's RNG). The boosting + IRLS use NO R RNG; the torch RNG
  # save/restore lives INSIDE .coda_codacore_boost (so any reticulate side effect
  # on R's RNG from the torch restore completes BEFORE this cleanup runs). This
  # on.exit therefore runs LAST and restores the pre-fit R-RNG state exactly --
  # keeping fit globally RNG-neutral, including the no-prior-seed case.
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

  .coda_codacore_require_torch()

  feat_order <- colnames(X_train)
  X_use <- .coda_codacore_align(X_train, feat_order)   # n x p (already in order)
  device <- .coda_codacore_resolve_device(hp$device)

  boosted <- .coda_codacore_boost(X_use, y, hp, device)

  model <- list(
    feature_universe = feat_order,
    balances = boosted$balances,
    intercept = boosted$head$intercept,
    weights = boosted$head$weights,
    n_balances = length(boosted$balances),
    device = device,
    seed = hp$seed,
    hp = hp
  )
  class(model) <- "coda_codacore_model"
  model
}


# ----------------------------------------------------------------------------
# score_coda_codacore (PURE BASE R -- no python)
# ----------------------------------------------------------------------------

#' @title Score the CoDaCoRe stagewise log-ratio-balance discriminator (pure base R)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_coda_codacore}}, in PURE base R (no python). Each query is aligned
#' to the frozen feature universe; for each frozen discrete balance \eqn{b} the
#' specimen's discrete ILR balance is computed over its OWN strictly-positive support
#' in \eqn{A_b / B_b} (a side with no positive feature contributes the neutral value
#' \code{0}); the score is the frozen logistic linear predictor
#' \eqn{\mathrm{intercept} + \sum_b w_b\,\mathrm{bal}_b(x)}. Larger = more case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} universe features present (a
#' column-overlap floor), or empty positive support over the frozen universe on the
#' aligned ORIGINAL abundances (\code{!any(X_use[i, ] > 0)}), return the neutral score
#' \code{0}. A FLAT all-equal-positive composition gives every balance value \code{0}
#' but is a VALID specimen, scored to the intercept (NOT floored). The score of a row
#' depends only on that row and the frozen model, uses no random numbers and no
#' scored-batch statistics, and is exactly invariant to per-specimen positive scaling.
#'
#' @param model A \code{coda_codacore_model} object from \code{\link{fit_coda_codacore}}.
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
#' model <- fit_coda_codacore(X, y)
#' score_coda_codacore(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Gordon-Rodriguez E, Quinn TP, Cunningham JP. (2021) Learning sparse log-ratios
#' for high-throughput sequencing data. \emph{Bioinformatics} 38(1):157-163.
#'
#' @export
score_coda_codacore <- function(model, X, meta = NULL) {
  if (!inherits(model, "coda_codacore_model")) {
    stop("score_coda_codacore: model must have class coda_codacore_model")
  }
  X <- .reo_check_matrix(X, "score_coda_codacore", "X")
  .reo_check_meta(meta, nrow(X), "score_coda_codacore", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  out <- rep(0, nrow(X))                              # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }
  # No accepted balances -> the head is intercept-only over the present specimens.
  if (model$n_balances == 0L || length(model$balances) == 0L) {
    X_use0 <- .coda_codacore_align(X, feat_order)
    active0 <- which(apply(X_use0, 1L, function(r) any(r > 0)))
    out[active0] <- model$intercept
    return(out)
  }

  X_use <- .coda_codacore_align(X, feat_order)

  for (i in seq_len(nrow(X_use))) {
    # Empty positive support is the only degenerate row floored to neutral 0.
    # Test it on the ORIGINAL aligned abundances, NOT on the balances being 0: a
    # FLAT composition (full positive support, equal abundances) yields all-zero
    # balances but is a valid specimen scored to the intercept.
    if (!any(X_use[i, ] > 0)) next                    # empty positive support
    xi <- X_use[i, ]
    bal <- vapply(model$balances, function(bb) {
      .coda_codacore_balance_one(xi, bb$A, bb$B)
    }, numeric(1L))
    out[i] <- model$intercept + sum(model$weights * bal)
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_coda_codacore: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen, score-determining R-side state. The model holds NO live
# external-pointer object (the torch relaxation is discarded after export; only the
# discrete A_b/B_b sets + logistic weights remain), so every field is a plain R
# object and this is a deterministic snapshot suitable as the `model_digest` for
# singlesample_assert_row_equivariant() (clause (d): the bytes must be identical
# before and after scoring). Two seed-matched CPU fits produce identical digests.
.coda_codacore_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe,
    balances = model$balances,
    intercept = model$intercept,
    weights = model$weights,
    n_balances = model$n_balances,
    seed = model$seed,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
