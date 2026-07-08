#' @title TabDPT in-context tabular-foundation discriminator (conditional single-sample)
#'
#' @description
#' Conditional-single-sample discriminator (family O) built on TabDPT, the
#' pretrained tabular foundation model of Ma et al. (2024). Like TabPFN and
#' TabICL, TabDPT classifies a query row by \emph{in-context learning}: a forward
#' pass of a frozen pretrained transformer conditions on a provided labelled
#' TRAINING context (optionally retrieval-augmented) and emits the query's class
#' posterior. There is no per-dataset gradient training -- "fitting" only loads
#' the training context into the model. This scorer is the direct sibling of
#' \code{\link{fit_tabpfn}} and \code{\link{fit_tabicl}}, swapping the in-context
#' forward calls for TabDPT's.
#'
#' \strong{Why this is single-sample-deployable (conditional-in-N).} A specimen
#' is single-sample-classifiable iff the TRAINING context that conditions the
#' transformer is FROZEN into the model at fit time, so that at inference the
#' query's posterior depends only on that query row and the frozen context --
#' never on which other specimens are scored alongside it. This implementation
#' freezes the per-sample-rCLR-transformed training context (rows + labels) and
#' the feature universe at fit, and at score conditions every query against that
#' one frozen context. The frozen context is the TRAINING rows only, so it is
#' leakage-safe (no scored / test data ever enters the context). The scorer pins
#' the context-reduction mode to \code{"subsample"} (a single context shared
#' across queries, not per-query retrieval) and passes the full frozen context
#' (no per-query subsampling), so the conditioning set is identical for every
#' query and never depends on the scored batch.
#'
#' \strong{Row-independent scoring modes.} With the frozen shared (subsample)
#' context and the full context passed to every query, TabDPT's
#' \code{predict_proba} is invariant to the row ORDER of the scored batch and to
#' the scored subset -- a query's posterior does not depend on the other queries
#' (measured companion maxdiff \eqn{0}). The default deployment path still
#' evaluates \code{predict_proba} ONE QUERY ROW AT A TIME, so the \eqn{n=1}
#' forward path is used uniformly. The optional benchmark-only \code{score_batch}
#' path evaluates balanced chunks of query rows. That batched path is
#' leakage-free and AUC-faithful (\eqn{|dAUC| = 0}; specimen ranks/verdicts
#' unchanged), but it is not bit-identical to the \eqn{n=1} row-by-row scorer:
#' the installed build has a deterministic batch-size numerical effect (observed
#' \code{max|b-r| = 1.97e-3}), not a cross-row coupling.
#'
#' \strong{Single-sample transform.} Each specimen is mapped to the
#' self-contained per-sample robust CLR (rCLR) over its OWN strictly-positive
#' support (geometric-mean centring on \code{v > 0}) in the FROZEN feature
#' universe; this is exactly invariant to per-specimen scaling and uses no
#' cross-row statistic. The SAME rCLR transform is applied to the frozen training
#' context at fit and to each query at score. Universe features absent from a
#' specimen carry the neutral rCLR value \code{0} (centred log-ratio of an absent
#' coordinate), matching the missing-feature rule of the other roster scorers.
#'
#' \strong{Token-free pinned checkpoint.} The default TabDPT checkpoint
#' (\code{tabdpt1_2.safetensors}, the TabDPT v1.2 model of Ma et al. 2024) is
#' hosted in the public Hugging Face repository \code{Layer6/TabDPT} and is
#' fetched non-interactively by \code{hf_hub_download} -- there is no license
#' acceptance or interactive token step (unlike TabPFN's gated V2.6). This method
#' pins that checkpoint version explicitly (via the installed \code{tabdpt}
#' package's pinned \code{_VERSION = "1_2"}) so the model is reproducible and the
#' download path is token-free.
#'
#' \strong{Orientation and output.} The score is the frozen-context posterior
#' \eqn{P(\mathrm{case}) = P(y = 1)} for the query, in \eqn{[0, 1]}; larger =
#' more case-like. Degenerate queries (fewer than \code{hp$min_features} universe
#' features present, or an all-zero / empty specimen) return the neutral
#' probability \code{0.5}.
#'
#' @references
#' Ma J, Thomas V, Hosseinzadeh R, Kamkari H, Labach A, Cresswell JC, Golestan S,
#' Yu G, Volkovs M, Caterini AL. (2024) TabDPT: Scaling Tabular Foundation Models.
#' arXiv:2410.18164.
#'
#' Qu J, Holzmuller D, Varoquaux G, Le Morvan M. (2025) TabICL: A Tabular
#' Foundation Model for In-Context Learning on Large Data.
#' \emph{International Conference on Machine Learning (ICML)}; arXiv:2502.05564.
#'
#' Hollmann N, Muller S, Purucker L, et al. (2025) Accurate predictions on small
#' data with a tabular foundation model. \emph{Nature} 637:319-326.
#'
#' @name singlesample-tabdpt
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Deliberately does NOT reuse a package
# rclr helper that centres on a fixed pseudocount or cross-sample reference
# statistics (either would break exact scale-invariance / single-sample
# deployability). Mirrors .tabicl_rclr / .tabpfn_rclr / .sinkhorn_single_rclr.
.tabdpt_rclr <- function(v) {
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
.tabdpt_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .tabdpt_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

.tabdpt_score_chunks <- function(score_idx, chunk_size = 1024L) {
  n <- length(score_idx)
  if (n == 0L) return(list())
  n_chunks <- ceiling(n / chunk_size)
  sizes <- rep.int(n %/% n_chunks, n_chunks)
  extra <- n %% n_chunks
  if (extra > 0L) sizes[seq_len(extra)] <- sizes[seq_len(extra)] + 1L
  split(score_idx, rep.int(seq_len(n_chunks), sizes))
}


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

.tabdpt_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_tabdpt: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_tabdpt: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_tabdpt: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("device", "temperature", "normalizer", "max_context",
               "min_features", "seed", "score_batch")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_tabdpt: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_tabdpt: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  temperature <- hp$temperature
  if (is.null(temperature)) temperature <- 1.0
  if (!is.numeric(temperature) || length(temperature) != 1L ||
      !is.finite(temperature) || temperature <= 0) {
    stop("fit_tabdpt: hp$temperature must be a positive finite number")
  }

  normalizer <- hp$normalizer
  if (is.null(normalizer)) normalizer <- "standard"
  if (!is.character(normalizer) || length(normalizer) != 1L || is.na(normalizer) ||
      !normalizer %in% c("standard", "minmax", "robust", "power",
                         "quantile-uniform", "quantile-normal", "log1p")) {
    stop("fit_tabdpt: hp$normalizer must be one of 'standard', 'minmax', ",
         "'robust', 'power', 'quantile-uniform', 'quantile-normal', 'log1p'")
  }

  max_context <- hp$max_context
  if (is.null(max_context)) max_context <- 4096L
  if (!is.numeric(max_context) || length(max_context) != 1L ||
      !is.finite(max_context) || max_context < 2L ||
      max_context > .Machine$integer.max ||
      max_context != as.integer(max_context)) {
    stop("fit_tabdpt: hp$max_context must be an integer >= 2")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_tabdpt: hp$min_features must be a positive integer")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_tabdpt: hp$seed must be a single integer")
  }

  score_batch <- hp$score_batch
  if (is.null(score_batch)) score_batch <- FALSE
  if (!is.logical(score_batch) || length(score_batch) != 1L ||
      is.na(score_batch)) {
    stop("fit_tabdpt: hp$score_batch must be TRUE or FALSE")
  }
  list(
    device = device,
    temperature = as.numeric(temperature),
    normalizer = normalizer,
    max_context = as.integer(max_context),
    min_features = as.integer(min_features),
    seed = as.integer(seed),
    score_batch = score_batch
  )
}


# ----------------------------------------------------------------------------
# Python / tabdpt availability + device + classifier construction
# ----------------------------------------------------------------------------

# The pinned TabDPT checkpoint version. The installed tabdpt package pins
# _VERSION = "1_2", which resolves the public token-free HF file
# tabdpt1_2.safetensors in repo Layer6/TabDPT (fetched non-interactively via
# hf_hub_download; no license / token step). Recorded for digest / reproducibility.
.TABDPT_CHECKPOINT <- "tabdpt1_2.safetensors"

# Mirrors .tabicl_require_module: a clear, actionable error if reticulate or the
# tabdpt module is unavailable.
.tabdpt_require_module <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_tabdpt: package 'reticulate' is required to run TabDPT. ",
         "Install with install.packages('reticulate').", call. = FALSE)
  }
  if (!reticulate::py_module_available("tabdpt")) {
    stop("fit_tabdpt: Python module 'tabdpt' is not available in the active ",
         "reticulate environment. Install with: pip install tabdpt ",
         "(in the configured RETICULATE_PYTHON venv).", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a
# 'cuda' request) fall back to 'cpu' when CUDA is unavailable, so the method
# never errors on a CPU-only host. Uses torch via reticulate.
.tabdpt_resolve_device <- function(device) {
  if (identical(device, "cpu")) return("cpu")
  cuda_ok <- tryCatch({
    if (!reticulate::py_module_available("torch")) {
      FALSE
    } else {
      torch_py <- reticulate::import("torch", delay_load = FALSE)
      isTRUE(torch_py$cuda$is_available())
    }
  }, error = function(e) FALSE)
  if (isTRUE(cuda_ok)) "cuda" else "cpu"
}

# Build a TabDPT classifier on the pinned token-free checkpoint with the frozen
# device, normalizer and temperature config, and fit it on the frozen rCLR
# context. Returns the live (non-serialisable) python classifier; it is NEVER
# part of the model digest. `Z_ctx` is the n_ctx x p numeric rCLR context, `y_ctx`
# the integer 0/1 labels. Wraps fit in the torch RNG save/restore guard (TabDPT's
# fit / PCA can advance the global + CUDA torch RNG). compile is forced FALSE so
# the build never triggers torch.compile (non-determinism / slow first call); the
# scorer's exactness does not depend on the compiled path. context_reduction is
# pinned to "subsample" so a single context is shared across queries (not per-query
# retrieval), which keeps scoring row-independent.
.tabdpt_build_and_fit <- function(Z_ctx, y_ctx, device, seed, normalizer) {
  tabdpt <- reticulate::import("tabdpt", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)

  clf <- tabdpt$TabDPTClassifier(
    device = device,
    normalizer = normalizer,
    context_reduction = "subsample",
    feature_reduction = "pca",
    compile = FALSE,
    use_flash = identical(device, "cuda"),
    verbose = FALSE
  )
  x_py <- np$asarray(Z_ctx, dtype = "float64")
  y_py <- np$asarray(as.integer(y_ctx), dtype = "int64")

  .tabdpt_with_torch_rng(seed, function() clf$fit(x_py, y_py))

  # TabDPT has no classes_ attribute: predict_proba returns one column per class
  # in ascending label order (the model is conditioned on the integer labels and
  # emits log-softmax over class ids 0..num_classes-1). With 0/1 labels the
  # columns are [P(y=0), P(y=1)], so the case column (label 1) is the position of
  # 1 in the sorted unique labels.
  classes <- sort(unique(as.integer(y_ctx)))
  case_col <- match(1L, classes)
  if (is.na(case_col)) {
    stop("fit_tabdpt: fitted TabDPT context does not contain the case label 1",
         call. = FALSE)
  }
  list(clf = clf, classes = classes, case_col = case_col)
}

# Run `fn()` with the torch RNG seeded to `seed`, restoring the torch RNG state
# afterwards so the global torch RNG is left byte-identical. TabDPT inference is
# deterministic regardless of the torch seed (verified: predict does not advance
# the torch RNG and re-scoring is bit-identical), but fit / PCA can advance it as
# a side effect; this keeps that side effect contained. Saves/restores BOTH the
# CPU torch RNG and the CUDA RNG state(s) so a GPU fit is CUDA-RNG neutral.
.tabdpt_with_torch_rng <- function(seed, fn) {
  torch_py <- tryCatch(reticulate::import("torch", delay_load = FALSE),
                       error = function(e) NULL)
  if (is.null(torch_py)) return(fn())
  state <- tryCatch(torch_py$random$get_rng_state(), error = function(e) NULL)
  # Also preserve the CUDA RNG state(s): a GPU fit advances the CUDA generator,
  # which the CPU-only get/set_rng_state does NOT cover. No-op on CPU-only hosts.
  cuda_ok <- tryCatch(isTRUE(torch_py$cuda$is_available()), error = function(e) FALSE)
  cuda_state <- if (cuda_ok) {
    tryCatch(torch_py$cuda$get_rng_state_all(), error = function(e) NULL)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(state)) {
      tryCatch(torch_py$random$set_rng_state(state), error = function(e) NULL)
    }
    if (!is.null(cuda_state)) {
      tryCatch(torch_py$cuda$set_rng_state_all(cuda_state), error = function(e) NULL)
    }
  }, add = TRUE)
  tryCatch(torch_py$manual_seed(as.integer(seed)), error = function(e) NULL)
  fn()
}


# ----------------------------------------------------------------------------
# fit_tabdpt
# ----------------------------------------------------------------------------

#' @title Fit the TabDPT in-context discriminator (freeze the training context)
#'
#' @description
#' Freezes the training context for TabDPT in-context classification. The
#' training matrix is mapped to the per-sample robust CLR over the frozen feature
#' universe (\code{colnames(X_train)}); if the context exceeds
#' \code{hp$max_context} rows it is reduced to a deterministic, seeded,
#' class-stratified subsample. A TabDPT classifier (pinned token-free checkpoint)
#' is constructed and \code{fit()} on that frozen rCLR context (which, for TabDPT,
#' simply loads the context into the model and fits the preprocessing scaler /
#' optional PCA -- no gradient training). The frozen R-side state needed to
#' reproduce scores (the rCLR context, labels, feature universe, class order,
#' device/config, hp) is stored on the model; the live python classifier is also
#' kept for scoring but is NEVER part of the model digest.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{device}
#'   (\code{"cpu"} (default), \code{"cuda"}, or \code{"auto"}; \code{"cuda"} /
#'   \code{"auto"} fall back to CPU when no GPU is available), \code{temperature}
#'   (TabDPT softmax temperature at predict, positive finite, default \code{1.0}),
#'   \code{normalizer} (TabDPT input scaler, one of \code{"standard"} (default),
#'   \code{"minmax"}, \code{"robust"}, \code{"power"}, \code{"quantile-uniform"},
#'   \code{"quantile-normal"}, \code{"log1p"}), \code{max_context} (row cap on the
#'   frozen context, integer \eqn{\ge 2}, default \code{4096L}; larger contexts are
#'   seeded class-stratified subsampled), \code{min_features} (feature-overlap
#'   floor at scoring, positive integer, default \code{3L}), \code{seed}
#'   (integer; default \code{42L}), and \code{score_batch} (benchmark-only batched
#'   scoring flag, logical, default \code{FALSE}).
#'
#' @return Object of class \code{tabdpt_model}: a list with \code{context_rclr}
#'   (frozen rCLR context matrix), \code{context_y}, \code{feature_universe},
#'   \code{classes}, \code{case_col}, \code{n_context}, \code{device} (resolved),
#'   \code{checkpoint} (pinned checkpoint version), \code{hp}, and the live python
#'   classifier \code{clf} (excluded from the digest).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 20; k <- 6
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' model <- fit_tabdpt(X, y)
#' score_tabdpt(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Ma J, Thomas V, Hosseinzadeh R, et al. (2024) TabDPT: Scaling Tabular
#' Foundation Models. arXiv:2410.18164.
#'
#' @export
fit_tabdpt <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_tabdpt", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_tabdpt", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_tabdpt")
  hp <- .tabdpt_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_tabdpt: X_train must contain at least hp$min_features features")
  }
  .tabdpt_require_module()

  feat_order <- colnames(X_train)
  Z <- .tabdpt_rclr_matrix(X_train)              # n x p frozen rCLR context

  # Deterministic, seeded, class-stratified context cap. The global .Random.seed
  # is saved and restored so the subsample uses an isolated RNG stream.
  keep <- .tabdpt_context_subsample(nrow(Z), y, hp$max_context, hp$seed)
  Z_ctx <- Z[keep, , drop = FALSE]
  y_ctx <- y[keep]
  if (!any(y_ctx == 1L) || !any(y_ctx == 0L)) {
    stop("fit_tabdpt: capped context lost a class; raise hp$max_context")
  }

  device <- .tabdpt_resolve_device(hp$device)
  fitted <- .tabdpt_build_and_fit(Z_ctx, y_ctx, device, hp$seed, hp$normalizer)

  model <- list(
    context_rclr = Z_ctx,
    context_y = y_ctx,
    feature_universe = feat_order,
    classes = fitted$classes,
    case_col = fitted$case_col,
    n_context = nrow(Z_ctx),
    device = device,
    checkpoint = .TABDPT_CHECKPOINT,
    hp = hp,
    clf = fitted$clf
  )
  class(model) <- "tabdpt_model"
  model
}

# Deterministic seeded class-stratified subsample of row indices to at most
# `cap` rows, preserving the original row order of the kept set. Returns the
# full 1:n when n <= cap. Saves/restores the global .Random.seed so it does not
# perturb the caller's RNG stream.
.tabdpt_context_subsample <- function(n, y, cap, seed) {
  if (n <= cap) return(seq_len(n))
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
  set.seed(seed)

  idx_case <- which(y == 1L)
  idx_ctrl <- which(y == 0L)
  # Proportional allocation, at least 1 per class present, summing to cap.
  n_case <- max(1L, min(length(idx_case),
                        as.integer(round(cap * length(idx_case) / n))))
  # Reserve >= 1 control slot when both classes are present and cap >= 2, so a
  # tiny cap with extreme imbalance never drops a class.
  if (length(idx_ctrl) > 0L && cap >= 2L) n_case <- min(n_case, cap - 1L)
  n_ctrl <- cap - n_case
  if (n_ctrl > length(idx_ctrl)) {
    n_ctrl <- length(idx_ctrl)
    n_case <- min(length(idx_case), cap - n_ctrl)
  }
  if (n_case > length(idx_case)) n_case <- length(idx_case)
  # Sample POSITIONS, never the values: base sample(x, k) treats a length-1 x as
  # 1:x (the classic gotcha), so sample(idx_ctrl, 1) with a single control index
  # would pick a wrong row.
  pick <- function(idx, k) idx[sample.int(length(idx), k)]
  keep <- c(
    if (n_case > 0L) pick(idx_case, n_case) else integer(0),
    if (n_ctrl > 0L) pick(idx_ctrl, n_ctrl) else integer(0)
  )
  sort(keep)
}


# ----------------------------------------------------------------------------
# score_tabdpt
# ----------------------------------------------------------------------------

#' @title Score the TabDPT in-context discriminator (default row-by-row)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_tabdpt}}. Each query is mapped to the per-sample robust CLR
#' over the frozen feature universe (absent universe features carry the neutral
#' rCLR value \code{0}), then classified by TabDPT against the FROZEN training
#' context. By default, queries are passed to \code{predict_proba} ONE ROW AT A
#' TIME (the \eqn{n=1} forward path used uniformly); the score is the case
#' posterior \eqn{P(y = 1)} in \eqn{[0, 1]}, larger = more case-like.
#'
#' Queries with fewer than \code{model$hp$min_features} present universe features,
#' or an all-zero / empty specimen, return the neutral probability \code{0.5}.
#' The score of a row depends only on that row and the frozen model, and is
#' exactly invariant to per-specimen positive scaling.
#'
#' @details A benchmark-only \code{score_batch} option evaluates
#' \code{predict_proba} once over balanced chunks of query rows. This path is
#' leakage-free and AUC-faithful (\eqn{|dAUC| = 0}; specimen ranks/verdicts
#' unchanged), but is not bit-identical to the default \eqn{n = 1} row-by-row
#' path because of a deterministic batch-size numerical effect. The default
#' row-by-row single-specimen deployment path is unchanged.
#'
#' @param model A \code{tabdpt_model} object from \code{\link{fit_tabdpt}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)} in \eqn{[0, 1]};
#'   larger values are more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_tabdpt(X, y)
#' score_tabdpt(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Ma J, Thomas V, Hosseinzadeh R, et al. (2024) TabDPT: Scaling Tabular
#' Foundation Models. arXiv:2410.18164.
#'
#' @export
score_tabdpt <- function(model, X, meta = NULL) {
  if (!inherits(model, "tabdpt_model")) {
    stop("score_tabdpt: model must have class tabdpt_model")
  }
  X <- .reo_check_matrix(X, "score_tabdpt", "X")
  .reo_check_meta(meta, nrow(X), "score_tabdpt", "meta")

  feat_order <- model$feature_universe
  present <- intersect(feat_order, colnames(X))

  # Align X to the frozen universe order; absent features stay 0 (no support).
  X_use <- matrix(0, nrow = nrow(X), ncol = length(feat_order),
                  dimnames = list(rownames(X), feat_order))
  if (length(present) > 0L) X_use[, present] <- X[, present, drop = FALSE]

  out <- rep(0.5, nrow(X))                          # neutral default
  if (length(present) < model$hp$min_features) {
    return(out)
  }

  np <- reticulate::import("numpy", delay_load = FALSE)
  clf <- model$clf
  case_col <- model$case_col
  temperature <- model$hp$temperature

  if (isTRUE(model$hp$score_batch)) {
    Z <- .tabdpt_rclr_matrix(X_use)
    score_idx <- which(rowSums(Z != 0) > 0L)          # all-zero rCLR stays neutral
    if (length(score_idx) > 0L) {
      .tabdpt_with_torch_rng(model$hp$seed, function() {
        for (idx in .tabdpt_score_chunks(score_idx, chunk_size = 1024L)) {
          pr <- reticulate::py_to_r(
            clf$predict_proba(np$asarray(Z[idx, , drop = FALSE], dtype = "float64"),
                              temperature = temperature)
          )
          out[idx] <<- as.numeric(pr[, case_col])
        }
      })
    }
  } else {
    # Forced row-by-row predict_proba: the n = 1 forward path is used uniformly.
    # context_size is left NULL (= infinite -> full frozen context), so every
    # query conditions on the identical full frozen context. The torch RNG is
    # restored around the whole loop (TabDPT's predict does not advance it, but
    # the guard makes RNG-safety robust to build differences and to CUDA).
    .tabdpt_with_torch_rng(model$hp$seed, function() {
      for (i in seq_len(nrow(X_use))) {
        z <- .tabdpt_rclr(X_use[i, ])
        if (all(z == 0)) next                          # all-zero specimen: neutral
        zi <- matrix(z, nrow = 1L)
        pr <- reticulate::py_to_r(
          clf$predict_proba(np$asarray(zi, dtype = "float64"),
                            temperature = temperature)
        )
        out[i] <<- as.numeric(pr[1, case_col])
      }
    })
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out)) ||
      any(out < 0 | out > 1)) {
    stop("score_tabdpt: scorer produced invalid output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen score-determining R-side state (rCLR training context,
# labels, feature universe, class order, device/config, hp). The live python
# TabDPT classifier (model$clf, a non-serialisable external pointer) is
# deliberately EXCLUDED, so the digest is a deterministic snapshot suitable as
# the `model_digest` for singlesample_assert_row_equivariant() (clause (d): the
# bytes must be identical before and after scoring). Scoring never mutates this
# frozen state (TabDPT inference is read-only over the frozen context).
.tabdpt_model_digest <- function(model) {
  hp <- model$hp
  hp$score_batch <- NULL
  state <- list(
    context_rclr = model$context_rclr,
    context_y = model$context_y,
    feature_universe = model$feature_universe,
    classes = model$classes,
    case_col = model$case_col,
    n_context = model$n_context,
    device = model$device,
    checkpoint = model$checkpoint,
    hp = hp
  )
  digest::digest(state, algo = "xxhash64")
}
