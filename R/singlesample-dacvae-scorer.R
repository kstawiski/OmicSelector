#' @title Domain-adversarial conditional VAE single-sample discriminator (DA-cVAE)
#'
#' @description
#' Single-sample within-cohort discriminator (family N) implementing Stage 2/3 of
#' idea_report_1570797495812882433 ("Structurally Informed ILR + DA-cVAE"). The
#' encoder ingests the \strong{Stage-1 structural-ILR balances} (the genomic-cluster /
#' GC-tertile Sequential Binary Partition of \code{\link{fit_ss_struct_ilr}}, computed
#' per specimen via \code{\link{ws_balance_ilr}} and standardised with train-only
#' statistics) and maps them to a latent \eqn{z}. Training combines three FIT-only
#' regularisers on that latent — a conditional VAE decoder \eqn{p(\mathrm{balances}\mid
#' z, \mathrm{platform})}, a KL prior, and a \strong{gradient-reversal platform
#' adversary} — with a supervised disease LABEL head. The encoder + label head are then
#' \strong{FROZEN} and exported to R; a new specimen is scored in PURE base-R by its
#' frozen disease logit on the deterministic latent \eqn{z = \mu}. Larger = more
#' case-like.
#'
#' \strong{Why this architecture (cross-platform transportability).} The platform/assay
#' label conditions ONLY the decoder (conditional generation) and drives the GRL
#' adversary; the ENCODER never sees the platform, so a held-out specimen is scored
#' with no platform/cohort covariate (single-sample-faithful). The reversed adversary
#' gradient pushes \eqn{z} toward a representation whose platform is not recoverable —
#' the mechanism for \emph{reduced platform imprint / improved cross-platform
#' transportability}. As enforced in the design (v2.1 §0), this is the claim ceiling:
#' the method targets reduced platform imprint, NOT platform "invariance", and every
#' result from it is EXPLORATORY and never merged into the frozen primary headline.
#'
#' \strong{The DA-cVAE objective.} Let \eqn{b} be the standardised structural-ILR
#' balances of a specimen, \eqn{(\mu, \log\sigma^2) = \mathrm{enc}(b)},
#' \eqn{z = \mu + \sigma \odot \epsilon} (\eqn{\epsilon \sim N(0,I)}; ONE draw per row
#' per step), \eqn{e_g = \mathrm{onehot}(\mathrm{platform})}. Four terms are minimised
#' jointly with Adam:
#' \deqn{\mathcal{L} = \underbrace{\overline{\mathrm{BCE}}(w^\top\mu + b_0,\, y)}_{\text{label (on }\mu\text{)}}
#'   + \rho\,\underbrace{\overline{\lVert \mathrm{dec}(z, e_g) - b\rVert^2}}_{\text{recon}}
#'   + \beta\,\underbrace{\overline{\mathrm{KL}(q(z\mid b)\,\Vert\,N(0,I))}}_{\text{prior}}
#'   + \underbrace{\mathrm{CE}\big(\mathrm{adv}(\mathrm{GRL}_\lambda(z)),\, g\big)}_{\text{platform adversary}}.}
#' The label head is trained on the DETERMINISTIC \eqn{\mu} (exactly the latent used at
#' score), so the kept readout sees no train/score representation mismatch. The decoder,
#' logvar head, and adversary are FIT-only and are NOT exported.
#'
#' \strong{Platforms from \code{meta_train}.} The platform/assay domains are the
#' distinct non-missing labels in \code{meta_train[[platform_col]]} (auto-detected from
#' \code{"platform"}, \code{"assay"}, \code{"technology"}, \code{"accession"},
#' \code{"batch"}, \code{"study"}, \code{"dataset"} when \code{NULL}). The conditional
#' decoder uses ALL \eqn{\ge 2} platforms whenever they exist, INDEPENDENT of
#' \code{grl_lambda}; the GRL adversary is engaged SEPARATELY, only when
#' \code{grl_lambda > 0} (recorded as \code{model$adv_engaged}). Setting
#' \code{grl_lambda = 0} with \eqn{\ge 2} platforms is therefore the design's
#' \strong{no-GRL ablation}: the SAME architecture (conditional decoder preserved),
#' adversary off — exactly the baseline the design's non-negotiable no-harm guard
#' compares against. A single platform (or \code{NULL} \code{meta_train}) collapses the
#' decoder to unconditional (\code{n_platforms == 1}) with the adversary necessarily
#' off.
#'
#' \strong{Torch is confined to FIT.} After training, the encoder weights (trunk up to
#' and including the \eqn{\mu} head) and the disease LABEL head \eqn{(w, b_0)} are
#' EXPORTED to R as plain float64 matrices/vectors; the python module (logvar head,
#' decoder, adversary, GRL) is DISCARDED. The fitted model holds NO external pointer and
#' NO python dependency at score time — the encoder forward and the linear logit are
#' pure base-R, so the suite runs the scorer on any host without the venv, and the
#' default \code{model_digest} is a stable deterministic snapshot.
#'
#' \strong{The encoder forward + score.} The exported encoder is
#' \eqn{b \to [\mathrm{Linear}\to\phi]_{1}\to\dots\to[\mathrm{Linear}\to\phi]_{H}\to
#' \mathrm{Linear}\to\mu(L)} (activation after every layer EXCEPT the final \eqn{\mu}
#' layer; \code{"relu"} default or \code{"tanh"}), reimplemented exactly by
#' \code{.dacvae_forward}. The disease score is the frozen logit
#' \eqn{s = w^\top\mu + b_0}.
#'
#' \strong{Novelty (secondary capability).} A frozen kNN-distance novelty score over the
#' TRAIN control (reference) latents is also serialised (design F9: multi-reference, not
#' a single centroid); \code{\link{score_dacvae_novelty}} returns, per query, the mean
#' distance to its \code{novelty_k} nearest reference latents. This supports the
#' cross-platform / out-of-distribution capability and is NOT the disease score.
#'
#' \strong{Single-sample contract.} Each specimen's balances are computed from its OWN
#' values (the \code{\link{ws_balance_ilr}} per-sample pseudocount makes them exactly
#' invariant to per-specimen positive scaling), standardised with the FROZEN train
#' mean/sd (missing balances imputed to the standardised mean \code{0}), and embedded by
#' the frozen encoder. No cross-row / batch statistic is used, so a specimen's score
#' depends only on that specimen and the frozen model — it passes
#' \code{\link{singlesample_assert_row_equivariant}} (single row == batch EXACTLY,
#' exact per-specimen scale invariance). Specimens with empty positive support over the
#' frozen universe, or fewer than \code{hp$min_features} universe features present,
#' return the neutral score \code{0}.
#'
#' @references
#' Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes. \emph{ICLR}.
#' arXiv:1312.6114.
#'
#' Sohn K, Lee H, Yan X. (2015) Learning Structured Output Representation using Deep
#' Conditional Generative Models. \emph{NeurIPS} 28.
#'
#' Ganin Y, et al. (2016) Domain-Adversarial Training of Neural Networks. \emph{JMLR}
#' 17(59):1--35. arXiv:1505.07818.
#'
#' @name singlesample-dacvae
NULL


# ----------------------------------------------------------------------------
# Encoder input: standardised structural-ILR balances (reuses Stage 1)
# ----------------------------------------------------------------------------

# Raw structural-ILR balances of X over the FROZEN SBP, as an n x nB double matrix
# (columns ordered by the frozen balance names). Reuses ws_balance_ilr exactly as the
# Stage-1 scorer does, so the encoder input is the SAME structural representation.
.dacvae_raw_balances <- function(X, sbp, pseudocount, aggregate, min_balance_coverage,
                                 balance_names) {
  B <- ws_balance_ilr(X, balances = sbp, pseudocount = pseudocount,
                      aggregate = aggregate,
                      min_balance_coverage = min_balance_coverage)
  B <- as.matrix(B)
  storage.mode(B) <- "double"
  if (!is.null(colnames(B)) && length(balance_names) == ncol(B)) {
    idx <- match(balance_names, colnames(B))
    if (!any(is.na(idx))) B <- B[, idx, drop = FALSE]
  }
  B
}

# Standardise raw balances with frozen train center/scale; impute any missing balance
# (insufficient coverage -> NA) to the standardised mean 0 (a row-local imputation that
# preserves row-equivariance: the NA pattern depends only on the row's own coverage).
.dacvae_standardize <- function(B, center, scale_) {
  Bs <- sweep(sweep(B, 2L, center, "-"), 2L, scale_, "/")
  Bs[!is.finite(Bs)] <- 0
  Bs
}


# ----------------------------------------------------------------------------
# qPCR Ct -> relative abundance via a FROZEN train-derived per-miRNA reference
# ----------------------------------------------------------------------------
# Pre-registered model (design F3): relative abundance = E^(-dCt), dCt = Ct - ref, with
# ref = the per-miRNA MEDIAN Ct over the TRAINING samples (a frozen per-feature vector,
# NOT a sample-specific housekeeping gene, NOT a global scalar). Censored / undetected
# Ct (NA) maps to the MZR floor (a small positive abundance), never to 0 and never
# dropped. The reference + efficiency + floor are computed ONCE on train, serialized,
# and applied UNCHANGED at score time -> single-sample-faithful (no test-batch stat).

# Fit the frozen Ct parameters from a TRAIN Ct matrix. Returns list(reference [named per
# feature], efficiency E, floor). The floor defaults to half the smallest observed
# converted abundance (so a censored Ct sits below every detected one); a positive
# hp$ct_floor overrides it.
.dacvae_fit_ct <- function(X_ct, E, floor_hp) {
  ref <- apply(X_ct, 2L, function(v) stats::median(v, na.rm = TRUE))
  ref[!is.finite(ref)] <- 0                                  # all-NA feature -> ref 0
  names(ref) <- colnames(X_ct)
  A_obs <- E^(-(sweep(X_ct, 2L, ref, "-")))                  # NA Ct -> NA abundance
  obs <- A_obs[is.finite(A_obs) & A_obs > 0]
  floor <- if (!is.null(floor_hp)) floor_hp
           else if (length(obs)) 0.5 * min(obs) else 1e-6
  list(reference = ref, efficiency = E, floor = floor)
}

# Apply the frozen Ct parameters to a Ct matrix -> relative-abundance matrix (same
# dimnames). Features without a frozen reference, and censored/undetected (NA or
# non-positive) conversions, take the MZR floor. Per-cell with a frozen per-feature
# reference -> row-equivariant.
.dacvae_apply_ct <- function(X_ct, ct) {
  ref <- ct$reference; E <- ct$efficiency; fl <- ct$floor
  A <- matrix(fl, nrow = nrow(X_ct), ncol = ncol(X_ct), dimnames = dimnames(X_ct))
  common <- intersect(colnames(X_ct), names(ref))
  if (length(common) > 0L) {
    sub <- X_ct[, common, drop = FALSE]
    conv <- E^(-(sweep(sub, 2L, ref[common], "-")))
    conv[!is.finite(conv) | conv <= 0] <- fl                 # censored -> MZR floor
    A[, common] <- conv
  }
  A
}

#' @title Convert qPCR Ct values to relative abundance with a fitted DA-cVAE's frozen reference
#'
#' @description
#' Applies a fitted \code{ct}-mode \code{\link{fit_dacvae}} model's FROZEN per-miRNA Ct
#' reference (median train Ct), assay efficiency \eqn{E}, and censored-Ct MZR floor to a
#' matrix of qPCR cycle-threshold values, returning the relative-abundance matrix
#' \eqn{E^{-(Ct - \mathrm{ref})}} (censored / undetected Ct \eqn{\to} the floor). Exposed
#' for the cross-platform harness to convert qPCR cohorts to the abundance simplex with
#' the SAME frozen reference the model was fit under (single-sample-faithful). For an
#' \code{abundance}-mode model this is a no-op identity.
#'
#' @param model A \code{dacvae_model} from \code{\link{fit_dacvae}}.
#' @param X Numeric matrix of qPCR Ct values (samples x features); named feature columns.
#' @return Relative-abundance matrix with the same dimensions/dimnames as \code{X}.
#' @seealso \code{\link{fit_dacvae}}, \code{\link{score_dacvae}}
#' @export
dacvae_ct_to_abundance <- function(model, X) {
  if (!inherits(model, "dacvae_model")) {
    stop("dacvae_ct_to_abundance: model must have class dacvae_model")
  }
  X <- .reo_check_matrix(X, "dacvae_ct_to_abundance", "X")
  if (!identical(model$input_type, "ct") || is.null(model$ct)) return(X)
  .dacvae_apply_ct(X, model$ct)
}


# ----------------------------------------------------------------------------
# Pure-R encoder forward (exact reimplementation of the torch encoder up to mu)
# ----------------------------------------------------------------------------

# Pure base-R forward of the exported encoder up to (and including) the mu head.
# `weights` is a list of per-layer list(W = out x in, b = out): layers 1..(L-1) are
# hidden (Linear then `activation`), layer L is the mu head (bare Linear, NO
# activation) -- matching the torch encoder exactly. torch computes M %*% t(W) + b
# (row = sample). Returns an n x L latent mean matrix.
.dacvae_forward <- function(M, weights, activation) {
  act <- if (identical(activation, "tanh")) tanh else function(x) pmax(x, 0)
  L <- length(weights)
  H <- M
  for (l in seq_len(L)) {
    W <- weights[[l]]$W                        # out x in
    b <- weights[[l]]$b                        # length out
    H <- H %*% t(W)                            # n x out
    H <- sweep(H, 2L, b, `+`)
    if (l < L) H <- act(H)                     # no activation after the mu head
  }
  H
}

# Frozen disease logit for one latent mean z (length L): w . z + b0. Larger = case-like.
.dacvae_score_one <- function(z, head_w, head_b) {
  sum(head_w * z) + head_b
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.dacvae_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_dacvae: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) stop("fit_dacvae: all hp entries must be named")
    if (any(duplicated(nm))) {
      stop("fit_dacvae: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("enc_hidden", "dec_hidden", "latent_dim", "activation", "epochs", "lr",
               "weight_decay", "kl_beta", "recon_weight", "grl_lambda", "domain_hidden",
               "platform_col", "novelty_k", "input_type", "ct_efficiency", "ct_floor",
               "pseudocount", "aggregate",
               "min_balance_coverage", "min_part", "min_features", "device", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_dacvae: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  pos_int_vec <- function(v, what) {
    if (!is.numeric(v) || length(v) < 1L || any(!is.finite(v)) || any(v < 1L) ||
        any(v > .Machine$integer.max) || any(v != as.integer(v))) {
      stop(sprintf("fit_dacvae: hp$%s must be a vector of positive integers", what))
    }
    as.integer(v)
  }
  pos_int1 <- function(v, what) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1L ||
        v > .Machine$integer.max || v != as.integer(v)) {
      stop(sprintf("fit_dacvae: hp$%s must be a positive integer", what))
    }
    as.integer(v)
  }
  nn_num1 <- function(v, what) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0) {
      stop(sprintf("fit_dacvae: hp$%s must be a non-negative finite number", what))
    }
    as.numeric(v)
  }

  enc_hidden <- if (is.null(hp$enc_hidden)) c(64L, 32L) else pos_int_vec(hp$enc_hidden, "enc_hidden")
  dec_hidden <- if (is.null(hp$dec_hidden)) c(32L, 64L) else pos_int_vec(hp$dec_hidden, "dec_hidden")
  latent_dim <- if (is.null(hp$latent_dim)) 16L else pos_int1(hp$latent_dim, "latent_dim")

  activation <- hp$activation
  if (is.null(activation)) activation <- "relu"
  if (!is.character(activation) || length(activation) != 1L || is.na(activation) ||
      !activation %in% c("relu", "tanh")) {
    stop("fit_dacvae: hp$activation must be one of 'relu', 'tanh'")
  }

  epochs <- if (is.null(hp$epochs)) 200L else pos_int1(hp$epochs, "epochs")
  lr <- hp$lr
  if (is.null(lr)) lr <- 1e-3
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("fit_dacvae: hp$lr must be a positive finite number")
  }
  weight_decay <- if (is.null(hp$weight_decay)) 1e-4 else nn_num1(hp$weight_decay, "weight_decay")
  kl_beta <- if (is.null(hp$kl_beta)) 1.0 else nn_num1(hp$kl_beta, "kl_beta")
  recon_weight <- if (is.null(hp$recon_weight)) 1.0 else nn_num1(hp$recon_weight, "recon_weight")
  grl_lambda <- if (is.null(hp$grl_lambda)) 1.0 else nn_num1(hp$grl_lambda, "grl_lambda")

  domain_hidden <- hp$domain_hidden
  if (is.null(domain_hidden)) domain_hidden <- c(32L)
  if (!is.numeric(domain_hidden) || any(!is.finite(domain_hidden)) ||
      any(domain_hidden < 1L) || any(domain_hidden > .Machine$integer.max) ||
      any(domain_hidden != as.integer(domain_hidden))) {
    stop("fit_dacvae: hp$domain_hidden must be a vector of positive integers ",
         "(or integer(0) for a single linear adversary head)")
  }

  platform_col <- hp$platform_col
  if (!is.null(platform_col)) {
    if (!is.character(platform_col) || length(platform_col) != 1L ||
        is.na(platform_col) || !nzchar(platform_col)) {
      stop("fit_dacvae: hp$platform_col must be NULL or a single non-empty string")
    }
  }

  novelty_k <- if (is.null(hp$novelty_k)) 5L else pos_int1(hp$novelty_k, "novelty_k")

  # input_type: "abundance" (default; X is non-negative relative abundance) or "ct"
  # (X is qPCR cycle-threshold values, converted to relative abundance via a FROZEN
  # train-derived per-miRNA reference; see .dacvae_fit_ct / .dacvae_apply_ct).
  input_type <- hp$input_type
  if (is.null(input_type)) input_type <- "abundance"
  if (!is.character(input_type) || length(input_type) != 1L || is.na(input_type) ||
      !input_type %in% c("abundance", "ct")) {
    stop("fit_dacvae: hp$input_type must be one of 'abundance', 'ct'")
  }
  ct_efficiency <- hp$ct_efficiency
  if (is.null(ct_efficiency)) ct_efficiency <- 2.0
  if (!is.numeric(ct_efficiency) || length(ct_efficiency) != 1L ||
      !is.finite(ct_efficiency) || ct_efficiency <= 1) {
    stop("fit_dacvae: hp$ct_efficiency must be a finite number > 1 (e.g. 2 for 2^-dCt)")
  }
  ct_floor <- hp$ct_floor                       # NULL -> computed at fit (frozen)
  if (!is.null(ct_floor) && (!is.numeric(ct_floor) || length(ct_floor) != 1L ||
      !is.finite(ct_floor) || ct_floor <= 0)) {
    stop("fit_dacvae: hp$ct_floor must be NULL (fit-computed) or a positive finite number")
  }

  pseudocount <- hp$pseudocount
  if (!is.null(pseudocount) && (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      !is.finite(pseudocount) || pseudocount <= 0)) {
    stop("fit_dacvae: hp$pseudocount must be NULL or a positive finite number")
  }
  aggregate <- hp$aggregate
  if (is.null(aggregate)) aggregate <- "gmean"
  aggregate <- match.arg(aggregate, c("gmean", "trimmed_gmean"))
  mbc <- hp$min_balance_coverage
  if (is.null(mbc)) mbc <- 0.5
  if (!is.numeric(mbc) || length(mbc) != 1L || !is.finite(mbc) || mbc < 0 || mbc > 1) {
    stop("fit_dacvae: hp$min_balance_coverage must be in [0, 1]")
  }
  min_part <- if (is.null(hp$min_part)) 1L else pos_int1(hp$min_part, "min_part")
  min_features <- if (is.null(hp$min_features)) 3L else pos_int1(hp$min_features, "min_features")

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_dacvae: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }
  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_dacvae: hp$seed must be a single integer")
  }

  list(enc_hidden = enc_hidden, dec_hidden = dec_hidden, latent_dim = latent_dim,
       activation = activation, epochs = epochs, lr = as.numeric(lr),
       weight_decay = weight_decay, kl_beta = kl_beta, recon_weight = recon_weight,
       grl_lambda = grl_lambda, domain_hidden = as.integer(domain_hidden),
       platform_col = platform_col, novelty_k = novelty_k,
       input_type = input_type, ct_efficiency = as.numeric(ct_efficiency),
       ct_floor = ct_floor, pseudocount = pseudocount,
       aggregate = aggregate, min_balance_coverage = mbc, min_part = min_part,
       min_features = min_features, device = device, seed = as.integer(seed))
}


# ----------------------------------------------------------------------------
# Platform resolution from meta_train (mirrors .dann_resolve_environments)
# ----------------------------------------------------------------------------

.dacvae_resolve_platforms <- function(meta_train, n, platform_col) {
  if (is.null(meta_train)) return(NULL)
  meta_train <- as.data.frame(meta_train, stringsAsFactors = FALSE)
  if (nrow(meta_train) != n) {
    stop("fit_dacvae: meta_train must have one row per row of X_train")
  }
  col <- platform_col
  if (is.null(col)) {
    candidates <- c("platform", "assay", "technology", "accession", "batch",
                    "study", "dataset")
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

.dacvae_require_torch <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_dacvae: package 'reticulate' is required to train the DA-cVAE. ",
         "Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("torch")) {
    stop("fit_dacvae: Python module 'torch' is not available in the active ",
         "reticulate environment. Install torch in the configured ",
         "RETICULATE_PYTHON venv.", call. = FALSE)
  }
  invisible(TRUE)
}

.dacvae_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    torch_py <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(torch_py$cuda$is_available())
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}


# ----------------------------------------------------------------------------
# DA-cVAE joint training + weight export
# ----------------------------------------------------------------------------

# Train encoder (-> mu, logvar) + disease LABEL head (on mu) + FIT-only conditional
# decoder p(balances | z, platform) + FIT-only GRL platform adversary, jointly by the
# DA-cVAE objective, then export the encoder weights (trunk up to and INCLUDING the mu
# head) and the LABEL head (w, b0). The decoder, logvar head, adversary, and GRL are
# FIT-only and NOT exported. The torch RNG (CPU + guarded CUDA) is saved/restored so fit
# leaves the generators byte-unchanged; the global R .Random.seed is saved/restored by
# the caller.
#
# `Bstd` is the n x nB standardised structural-ILR training matrix (encoder input);
# `y` the integer 0/1 disease labels; `plat_idx` an integer vector of length n giving
# each row's platform index in 1..n_plat (n_plat >= 1, 0-based for python). The
# adversary is engaged only when n_plat >= 2 AND grl_lambda > 0.
.dacvae_train_export <- function(Bstd, y, plat_idx, n_plat, hp, device) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn

  # ---- save the torch RNG (CPU + guarded CUDA), restore on exit -------------
  torch_state <- tryCatch(torch$random$get_rng_state(), error = function(e) NULL)
  cuda_ok <- tryCatch(isTRUE(torch$cuda$is_available()), error = function(e) FALSE)
  cuda_state <- if (cuda_ok) {
    tryCatch(torch$cuda$get_rng_state_all(), error = function(e) NULL)
  } else NULL
  on.exit({
    if (!is.null(torch_state)) {
      tryCatch(torch$random$set_rng_state(torch_state), error = function(e) NULL)
    }
    if (!is.null(cuda_state)) {
      tryCatch(torch$cuda$set_rng_state_all(cuda_state), error = function(e) NULL)
    }
  }, add = TRUE)

  adv_active <- (n_plat >= 2L) && (hp$grl_lambda > 0)

  # ---- define the gradient-reversal autograd.Function (FIT-only) ------------
  # Defined BEFORE seeding/building (a pure def draws no RNG and registers no
  # parameters, so it cannot perturb weight init). Identical to the dann GRL.
  if (adv_active) {
    reticulate::py_run_string("
import torch
import torch.nn as nn

class _DACVAE_GradReverseFn(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, lamb):
        ctx.lamb = lamb
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        return grad_output.neg() * ctx.lamb, None

def _dacvae_grad_reverse(x, lamb):
    return _DACVAE_GradReverseFn.apply(x, lamb)

class _DACVAE_AdvHead(nn.Module):
    def __init__(self, d, domain_hidden, n_dom, act_name):
        super().__init__()
        act = nn.Tanh if act_name == 'tanh' else nn.ReLU
        layers = []
        in_dim = d
        for h in domain_hidden:
            layers.append(nn.Linear(in_dim, int(h)))
            layers.append(act())
            in_dim = int(h)
        layers.append(nn.Linear(in_dim, int(n_dom)))
        self.net = nn.Sequential(*layers)

    def forward(self, z):
        return self.net(z)
", convert = FALSE)
    main <- reticulate::import_main()
    grad_reverse <- main$`_dacvae_grad_reverse`
  }

  # SEED BEFORE BUILD so two seed-matched fits initialise AND draw the reparam eps
  # identically, and the model is cross-process reproducible.
  torch$manual_seed(as.integer(hp$seed))

  np <- reticulate::import("numpy", delay_load = FALSE)
  B_py <- torch$as_tensor(np$asarray(Bstd, dtype = "float32"))$float()$to(device)
  y_py <- torch$as_tensor(np$asarray(as.numeric(y), dtype = "float32")
                          )$float()$reshape(c(-1L, 1L))$to(device)
  # platform one-hot for the conditional decoder (n x n_plat; a single ones column
  # when n_plat == 1 -> a plain conditional VAE on a constant).
  oh <- matrix(0, nrow = nrow(Bstd), ncol = n_plat)
  oh[cbind(seq_len(nrow(Bstd)), plat_idx)] <- 1
  G_py <- torch$as_tensor(np$asarray(oh, dtype = "float32"))$float()$to(device)

  nB <- ncol(Bstd)
  L <- hp$latent_dim
  act_mod <- if (identical(hp$activation, "tanh")) nn$Tanh else nn$ReLU

  # ---- ENCODER trunk: balances(nB) -> enc_hidden (act after each) -----------
  enc_layers <- list()
  in_dim <- nB
  for (h in hp$enc_hidden) {
    enc_layers[[length(enc_layers) + 1L]] <- nn$Linear(in_dim, as.integer(h))
    enc_layers[[length(enc_layers) + 1L]] <- act_mod()
    in_dim <- as.integer(h)
  }
  enc_trunk <- do.call(nn$Sequential, enc_layers)$to(device)
  enc_mu <- nn$Linear(in_dim, as.integer(L))$to(device)       # mu head (exported)
  enc_logvar <- nn$Linear(in_dim, as.integer(L))$to(device)   # logvar head (FIT-only)
  label_head <- nn$Linear(as.integer(L), 1L)$to(device)       # disease head (exported)

  # ---- DECODER (FIT-only): concat[z, e_g] (L + n_plat) -> dec_hidden -> nB ---
  dec_layers <- list()
  in_dim_d <- as.integer(L) + as.integer(n_plat)
  for (h in hp$dec_hidden) {
    dec_layers[[length(dec_layers) + 1L]] <- nn$Linear(in_dim_d, as.integer(h))
    dec_layers[[length(dec_layers) + 1L]] <- act_mod()
    in_dim_d <- as.integer(h)
  }
  dec_layers[[length(dec_layers) + 1L]] <- nn$Linear(in_dim_d, as.integer(nB))
  decoder <- do.call(nn$Sequential, dec_layers)$to(device)

  # ---- ADVERSARY (FIT-only): built only when active -------------------------
  adv_head <- NULL
  plat_tensor <- NULL
  if (adv_active) {
    adv_head <- main$`_DACVAE_AdvHead`(
      as.integer(L), reticulate::r_to_py(as.list(as.integer(hp$domain_hidden))),
      as.integer(n_plat), hp$activation)$to(device)
    plat_tensor <- torch$as_tensor(
      np$asarray(as.integer(plat_idx - 1L), dtype = "int64"))$long()$to(device)
  }

  # ---- optimizer over ALL trainable parameters ------------------------------
  params <- c(reticulate::iterate(enc_trunk$parameters()),
              reticulate::iterate(enc_mu$parameters()),
              reticulate::iterate(enc_logvar$parameters()),
              reticulate::iterate(label_head$parameters()),
              reticulate::iterate(decoder$parameters()))
  if (adv_active) params <- c(params, reticulate::iterate(adv_head$parameters()))
  opt <- torch$optim$Adam(params, lr = hp$lr, weight_decay = hp$weight_decay)
  label_lossf <- nn$BCEWithLogitsLoss()
  adv_lossf <- nn$CrossEntropyLoss()
  beta <- hp$kl_beta
  rho <- hp$recon_weight
  grl_lambda <- hp$grl_lambda
  n <- nrow(Bstd)

  enc_trunk$train(); enc_mu$train(); enc_logvar$train()
  label_head$train(); decoder$train()
  if (adv_active) adv_head$train()
  for (ep in seq_len(hp$epochs)) {
    opt$zero_grad()
    enc_h <- enc_trunk(B_py)                         # n x last_enc
    mu <- enc_mu(enc_h)                              # n x L
    logvar <- enc_logvar(enc_h)                      # n x L
    # Reparameterize z = mu + exp(0.5*logvar) * eps (ONE eps draw / row / step).
    std <- logvar$mul(0.5)$exp()
    eps <- torch$randn(list(n, as.integer(L)), device = device)
    z <- mu$add(std$mul(eps))                        # n x L

    # Label head on the DETERMINISTIC mu (the latent used at score).
    label_loss <- label_lossf(label_head(mu), y_py)
    # Conditional decoder recon: sum over balances, mean over rows.
    xhat <- decoder(torch$cat(list(z, G_py), 1L))    # n x nB
    recon <- xhat$sub(B_py)$pow(2)$sum(1L)$mean()
    # KL( q(z|b) || N(0,I) ), mean over rows.
    kl <- torch$mul(
      torch$sub(torch$sub(torch$add(mu$pow(2), logvar$exp()), logvar), 1.0)$sum(1L),
      0.5)$mean()
    loss <- torch$add(label_loss, torch$add(recon$mul(rho), kl$mul(beta)))
    if (adv_active) {
      z_rev <- grad_reverse(z, grl_lambda)
      adv_loss <- adv_lossf(adv_head(z_rev), plat_tensor)
      loss <- torch$add(loss, adv_loss)
    }
    loss$backward()
    opt$step()
  }
  enc_trunk$eval(); enc_mu$eval(); enc_logvar$eval()
  label_head$eval(); decoder$eval()
  if (adv_active) adv_head$eval()

  # ---- export encoder (trunk up to + including mu head) + label head --------
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

  enc_weights <- c(to_layers(enc_trunk), list(one_linear(enc_mu)))  # trunk + mu head
  if (length(enc_weights) < 1L) {
    stop("fit_dacvae: exported encoder has no linear layers", call. = FALSE)
  }
  lab <- one_linear(label_head)

  list(weights = enc_weights, head_w = as.numeric(lab$W), head_b = as.numeric(lab$b))
}


# ----------------------------------------------------------------------------
# fit_dacvae
# ----------------------------------------------------------------------------

#' @title Fit the DA-cVAE single-sample discriminator
#'
#' @description
#' Trains a domain-adversarial conditional VAE on the standardised Stage-1
#' structural-ILR balances (encoder \eqn{\to\mu}, FIT-only conditional decoder
#' \eqn{p(\mathrm{balances}\mid z,\mathrm{platform})}, KL prior, FIT-only
#' gradient-reversal platform adversary, supervised disease label head on \eqn{\mu}; via
#' reticulate-python torch, full-batch Adam), EXPORTS the encoder (trunk up to and
#' including the \eqn{\mu} head) and the disease LABEL head to R as plain numeric arrays,
#' fits a frozen kNN-distance novelty reference over the train control latents, and
#' DISCARDS the python module. The fitted model holds no external pointer; scoring is
#' pure base-R with no platform/cohort input.
#'
#' @param X_train Numeric matrix (samples x features) of non-negative abundances with
#'   unique, non-empty feature names (colnames).
#' @param y_train 0/1 disease labels (1 = case), length \code{nrow(X_train)}, with at
#'   least one case and one control.
#' @param meta_train Optional per-sample metadata carrying a platform/assay column (the
#'   conditional-decoder + adversary domain). With \eqn{\ge 2} platforms the conditional
#'   decoder is always used; the adversary is engaged only when \code{grl_lambda > 0}
#'   (so \code{grl_lambda = 0} here is the no-GRL ablation: conditional decoder preserved,
#'   adversary off). Must have one row per row of \code{X_train}.
#' @param annotation data.frame with columns \code{feature} (miRNA id), \code{cluster}
#'   (genomic-cluster id), \code{gc} (GC fraction) — the Stage-1 SBP priors.
#' @param hp Optional list of hyperparameters (all frozen into the model). Fields:
#'   \code{enc_hidden} (encoder trunk widths; default \code{c(64L,32L)}),
#'   \code{dec_hidden} (FIT-only decoder widths; default \code{c(32L,64L)}),
#'   \code{latent_dim} (\eqn{L}; default \code{16L}), \code{activation}
#'   (\code{"relu"}/\code{"tanh"}), \code{epochs} (default \code{200L}), \code{lr}
#'   (default \code{1e-3}), \code{weight_decay} (default \code{1e-4}), \code{kl_beta}
#'   (\eqn{\beta}; default \code{1.0}), \code{recon_weight} (\eqn{\rho}; default
#'   \code{1.0}), \code{grl_lambda} (adversary strength; default \code{1.0};
#'   \code{0} = no-GRL ablation), \code{domain_hidden} (FIT-only adversary widths;
#'   default \code{c(32L)}), \code{platform_col} (\code{NULL} = auto-detect),
#'   \code{novelty_k} (kNN novelty neighbours; default \code{5L}), \code{input_type}
#'   (\code{"abundance"} (default) or \code{"ct"} for qPCR cycle-threshold input,
#'   converted to relative abundance via a frozen per-miRNA train-median reference),
#'   \code{ct_efficiency} (assay efficiency \eqn{E} for \eqn{E^{-\Delta Ct}}; default
#'   \code{2}), \code{ct_floor} (censored-Ct MZR floor abundance; \code{NULL} =
#'   fit-computed), \code{pseudocount} /
#'   \code{aggregate} / \code{min_balance_coverage} / \code{min_part} (Stage-1 balance
#'   controls), \code{min_features} (score-time overlap floor; default \code{3L}),
#'   \code{device} (\code{"cpu"}/\code{"cuda"}/\code{"auto"}), \code{seed}
#'   (default \code{42L}).
#'
#' @return Object of class \code{dacvae_model}: a list with \code{feature_universe},
#'   \code{sbp}, \code{balance_names}, balance controls (\code{pseudocount},
#'   \code{aggregate}, \code{min_balance_coverage}), \code{center}/\code{scale}
#'   (frozen balance standardisation), \code{weights} (exported encoder layers),
#'   \code{activation}, \code{head_w}/\code{head_b} (frozen disease head),
#'   \code{latent_dim}, \code{platform_levels}, \code{n_platforms} (platforms used by the
#'   conditional decoder; \code{1} = unconditional), \code{adv_engaged} (\code{FALSE} =
#'   no-GRL ablation), \code{grl_lambda}, \code{novelty_ref} (control-latent reference
#'   matrix), \code{novelty_k}, \code{device}, \code{seed}, \code{hp}.
#'
#' @seealso \code{\link{score_dacvae}}, \code{\link{score_dacvae_novelty}},
#'   \code{\link{fit_ss_struct_ilr}}
#' @export
fit_dacvae <- function(X_train, y_train, meta_train = NULL, annotation, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_dacvae", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_dacvae", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_dacvae")
  if (missing(annotation) || is.null(annotation)) {
    stop("fit_dacvae: annotation (feature/cluster/gc) is required")
  }
  annotation <- .ss_struct_ilr_check_annotation(annotation)
  hp <- .dacvae_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_dacvae: X_train must contain at least hp$min_features features")
  }

  # ---- qPCR Ct -> relative abundance (frozen per-miRNA reference), if requested ----
  # input_type = "ct": derive the per-miRNA Ct reference + MZR floor from TRAIN once,
  # serialize them, and convert to relative abundance; everything downstream operates on
  # abundance. input_type = "abundance" (default) is byte-identical to the prior path.
  if (identical(hp$input_type, "ct")) {
    ct_params <- .dacvae_fit_ct(X_train, hp$ct_efficiency, hp$ct_floor)
    X_ab <- .dacvae_apply_ct(X_train, ct_params)
  } else {
    ct_params <- NULL
    X_ab <- X_train
  }
  features <- colnames(X_ab)
  sbp <- .ss_struct_ilr_build_sbp(features, annotation, hp$min_part)
  balance_names <- names(sbp)

  model <- list(
    feature_universe = features,
    sbp = sbp,
    balance_names = balance_names,
    input_type = hp$input_type, ct = ct_params,
    pseudocount = hp$pseudocount,
    aggregate = hp$aggregate,
    min_balance_coverage = hp$min_balance_coverage,
    center = numeric(0), scale = numeric(0),
    weights = list(), activation = hp$activation,
    head_w = numeric(0), head_b = 0,
    latent_dim = hp$latent_dim,
    platform_levels = character(0), n_platforms = 1L, adv_engaged = FALSE,
    grl_lambda = hp$grl_lambda,
    novelty_ref = matrix(numeric(0), nrow = 0L, ncol = 0L),
    novelty_k = hp$novelty_k,
    device = "cpu", seed = hp$seed, hp = hp
  )
  if (length(sbp) == 0L) {
    # Degenerate: no structural balances constructible -> neutral scorer (all 0).
    class(model) <- "dacvae_model"
    return(model)
  }

  # ---- balances + frozen train-only standardisation -------------------------
  B <- .dacvae_raw_balances(X_ab, sbp, hp$pseudocount, hp$aggregate,
                            hp$min_balance_coverage, balance_names)
  center <- apply(B, 2L, function(v) { m <- mean(v, na.rm = TRUE); if (is.finite(m)) m else 0 })
  scale_ <- apply(B, 2L, function(v) { s <- stats::sd(v, na.rm = TRUE); if (is.finite(s) && s > 0) s else 1 })
  Bstd <- .dacvae_standardize(B, center, scale_)

  # ---- resolve platforms (conditional-decoder domains) ----------------------
  # The conditional decoder p(balances | z, platform) uses ALL platforms whenever >=2
  # exist, INDEPENDENT of grl_lambda. The GRL adversary is engaged SEPARATELY (only
  # when grl_lambda > 0; see adv_engaged + .dacvae_train_export's adv_active). This
  # decoupling is what makes grl_lambda=0 with >=2 platforms the design's true no-GRL
  # ablation: SAME architecture (conditional decoder preserved), adversary off -- the
  # reference baseline for the non-negotiable no-harm guard. (Collapsing the decoder to
  # unconditional when grl_lambda=0 would confound adversary-off with
  # conditional->unconditional and mis-specify that guard.)
  plat <- .dacvae_resolve_platforms(meta_train, nrow(X_train), hp$platform_col)
  plat_levels <- character(0)
  if (!is.null(plat)) plat_levels <- sort(unique(plat[!is.na(plat)]))
  if (length(plat_levels) >= 2L) {
    plat_idx <- match(plat, plat_levels)
    plat_idx[is.na(plat_idx)] <- 1L                  # missing -> first platform
    n_plat <- length(plat_levels)
  } else {
    plat_idx <- rep.int(1L, nrow(X_train))           # single constant platform
    n_plat <- 1L
    plat_levels <- if (length(plat_levels) >= 1L) plat_levels[1L] else "all"
  }
  adv_engaged <- (n_plat >= 2L) && (hp$grl_lambda > 0)

  # Save the global R .Random.seed BEFORE any reticulate import (reticulate's first
  # python init seeds R's RNG as a side effect; fit draws no R randomness of its own).
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)

  .dacvae_require_torch()
  device <- .dacvae_resolve_device(hp$device)
  trained <- .dacvae_train_export(Bstd, y, plat_idx, n_plat, hp, device)

  # ---- frozen novelty reference: train CONTROL latents (mu) -----------------
  mu_train <- .dacvae_forward(Bstd, trained$weights, hp$activation)   # n x L
  ref_rows <- which(y == 0L)
  if (length(ref_rows) == 0L) ref_rows <- seq_len(nrow(mu_train))      # fallback: all
  novelty_ref <- mu_train[ref_rows, , drop = FALSE]
  dimnames(novelty_ref) <- NULL

  model$center <- center
  model$scale <- scale_
  model$weights <- trained$weights
  model$head_w <- as.numeric(trained$head_w)
  model$head_b <- as.numeric(trained$head_b)
  model$platform_levels <- plat_levels
  model$n_platforms <- n_plat
  model$adv_engaged <- adv_engaged
  model$novelty_ref <- novelty_ref
  model$device <- device
  class(model) <- "dacvae_model"
  model
}


# ----------------------------------------------------------------------------
# score_dacvae (disease) + score_dacvae_novelty (out-of-distribution)
# ----------------------------------------------------------------------------

# Shared per-row latent: align X to the present frozen universe, compute structural-ILR
# balances, standardise with frozen stats, encode to mu. Returns list(mu = n x L matrix
# (NA rows for floored specimens), floored = logical n) so both score functions share
# one forward and one degenerate rule.
.dacvae_latents <- function(model, X) {
  n <- nrow(X)
  L <- model$latent_dim
  out <- list(mu = matrix(NA_real_, nrow = n, ncol = L), floored = rep(TRUE, n))
  if (length(model$sbp) == 0L || length(model$weights) == 0L) return(out)
  # Ct-mode: convert qPCR Ct -> relative abundance with the FROZEN reference first.
  if (identical(model$input_type, "ct") && !is.null(model$ct)) {
    X <- .dacvae_apply_ct(X, model$ct)
  }
  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) return(out)
  Xr <- X[, present, drop = FALSE]

  B <- .dacvae_raw_balances(Xr, model$sbp, model$pseudocount, model$aggregate,
                            model$min_balance_coverage, model$balance_names)
  Bstd <- .dacvae_standardize(B, model$center, model$scale)
  mu <- .dacvae_forward(Bstd, model$weights, model$activation)
  # Floor only specimens with empty positive support over the present universe (NOT
  # all(z==0): a flat composition has all-zero balances but is a valid specimen). The
  # empty-support test is on the ORIGINAL aligned abundances.
  floored <- vapply(seq_len(n), function(i) !any(Xr[i, ] > 0), logical(1L))
  mu[floored, ] <- NA_real_
  list(mu = mu, floored = floored)
}

#' @title Score the DA-cVAE single-sample disease discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen \code{\link{fit_dacvae}}
#' model, in PURE R (no python, no platform/cohort input). Each query is mapped to its
#' standardised Stage-1 structural-ILR balances, embedded by the exported encoder to the
#' deterministic latent \eqn{z = \mu}, and scored by the frozen disease logit
#' \eqn{w^\top\mu + b_0}. Larger = more case-like. Specimens with empty positive support
#' over the frozen universe, fewer than \code{hp$min_features} universe features present,
#' or a degenerate (no-balance) model return the neutral score \code{0}. The score of a
#' row depends only on that row and the frozen model and is exactly invariant to
#' per-specimen positive scaling (\code{\link{singlesample_assert_row_equivariant}}).
#'
#' @param model A \code{dacvae_model} from \code{\link{fit_dacvae}}.
#' @param X Numeric matrix (samples x features) of non-negative abundances; named cols.
#' @param meta Optional per-sample metadata; accepted for interface uniformity, ignored
#'   (the platform adversary is FIT-only).
#' @return Finite numeric vector of length \code{nrow(X)}; higher = more case-like.
#' @seealso \code{\link{fit_dacvae}}, \code{\link{score_dacvae_novelty}}
#' @export
score_dacvae <- function(model, X, meta = NULL) {
  if (!inherits(model, "dacvae_model")) {
    stop("score_dacvae: model must have class dacvae_model")
  }
  X <- .reo_check_matrix(X, "score_dacvae", "X")
  .reo_check_meta(meta, nrow(X), "score_dacvae", "meta")

  lat <- .dacvae_latents(model, X)
  out <- rep(0, nrow(X))                              # neutral default
  for (i in seq_len(nrow(X))) {
    if (lat$floored[i]) next
    out[i] <- .dacvae_score_one(lat$mu[i, ], model$head_w, model$head_b)
  }
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_dacvae: scorer produced non-finite or wrong-length output")
  }
  out
}

#' @title Score the DA-cVAE single-sample novelty (out-of-distribution) capability
#'
#' @description
#' Returns, for each row of \code{X}, the mean Euclidean distance in the frozen latent
#' space from its encoded \eqn{\mu} to its \code{model$novelty_k} nearest TRAIN control
#' (reference) latents — a kNN-distance novelty score (design F9: a multi-reference
#' score, not a single-centroid Mahalanobis). Larger = more out-of-distribution relative
#' to the training reference population. This is a SECONDARY capability (it supports the
#' cross-platform / OOD use case) and is NOT the disease score. Floored / degenerate
#' specimens return \code{0}. Pure base-R, row-equivariant, exactly scale-invariant.
#'
#' @param model A \code{dacvae_model} from \code{\link{fit_dacvae}}.
#' @param X Numeric matrix (samples x features) of non-negative abundances; named cols.
#' @param meta Optional per-sample metadata; accepted for interface uniformity, ignored.
#' @return Finite numeric vector of length \code{nrow(X)}; higher = more novel.
#' @seealso \code{\link{fit_dacvae}}, \code{\link{score_dacvae}}
#' @export
score_dacvae_novelty <- function(model, X, meta = NULL) {
  if (!inherits(model, "dacvae_model")) {
    stop("score_dacvae_novelty: model must have class dacvae_model")
  }
  X <- .reo_check_matrix(X, "score_dacvae_novelty", "X")
  .reo_check_meta(meta, nrow(X), "score_dacvae_novelty", "meta")

  ref <- model$novelty_ref
  out <- rep(0, nrow(X))
  if (is.null(ref) || nrow(ref) == 0L) {
    return(out)
  }
  k <- min(as.integer(model$novelty_k), nrow(ref))
  lat <- .dacvae_latents(model, X)
  for (i in seq_len(nrow(X))) {
    if (lat$floored[i]) next
    d <- sqrt(rowSums(sweep(ref, 2L, lat$mu[i, ], "-")^2))
    out[i] <- mean(sort(d)[seq_len(k)])
  }
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_dacvae_novelty: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen, score-determining R-side state. The model holds NO live
# external-pointer object (the python module is discarded at fit), so every field is a
# plain R object and this is a deterministic snapshot suitable as the model_digest for
# singlesample_assert_row_equivariant() (clause (d): bytes identical before/after score).
.dacvae_model_digest <- function(model) {
  state <- list(
    feature_universe = model$feature_universe, sbp = model$sbp,
    balance_names = model$balance_names,
    input_type = model$input_type, ct = model$ct,
    pseudocount = model$pseudocount,
    aggregate = model$aggregate, min_balance_coverage = model$min_balance_coverage,
    center = model$center, scale = model$scale, weights = model$weights,
    activation = model$activation, head_w = model$head_w, head_b = model$head_b,
    latent_dim = model$latent_dim, platform_levels = model$platform_levels,
    n_platforms = model$n_platforms, adv_engaged = model$adv_engaged,
    grl_lambda = model$grl_lambda,
    novelty_ref = model$novelty_ref, novelty_k = model$novelty_k,
    seed = model$seed, hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
