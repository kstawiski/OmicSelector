#' @title TabPFN-v2 in-context foundation-model discriminator (conditional single-sample)
#'
#' @description
#' Conditional-single-sample discriminator (family O) built on TabPFN-v2, the
#' pretrained tabular foundation model of Hollmann et al. (Nature 2025). TabPFN
#' classifies a query row by \emph{in-context learning}: a single forward pass of
#' a frozen pretrained transformer conditions on a provided labelled TRAINING
#' context and emits the query's class posterior. There is no per-dataset
#' gradient training -- "fitting" only loads the training context into the model.
#'
#' \strong{Why this is single-sample-deployable (conditional-in-N).} A specimen
#' is single-sample-classifiable iff the TRAINING context that conditions the
#' transformer is FROZEN into the model at fit time, so that at inference the
#' query's posterior depends only on that query row and the frozen context --
#' never on which other specimens are scored alongside it. This implementation
#' freezes the per-sample-rCLR-transformed training context (rows + labels) and
#' the feature universe at fit, and at score conditions every query against that
#' one frozen context. The frozen context is the TRAINING rows only, so it is
#' leakage-safe (no scored / test data ever enters the context).
#'
#' \strong{Row-independent scoring modes.} TabPFN's \code{predict_proba} is
#' invariant to the row ORDER of the scored batch and to the scored subset -- a
#' query's posterior does not depend on the other queries (measured companion
#' maxdiff \eqn{0}). The default deployment path still evaluates
#' \code{predict_proba} ONE QUERY ROW AT A TIME, so the \eqn{n=1} forward path is
#' used uniformly. The optional benchmark-only \code{score_batch} path evaluates
#' balanced chunks of query rows. That batched path is leakage-free and
#' AUC-faithful (observed \eqn{|dAUC| = 7.8\mathrm{e}{-5}} at production scale;
#' specimen ranks/verdicts unchanged), but it is not bit-identical to the
#' \eqn{n=1} row-by-row path: the installed build has a deterministic batch-size
#' numerical effect (raw probability shifts around \eqn{1\mathrm{e}{-4}} to
#' \eqn{2.6\mathrm{e}{-3}}), not a cross-row coupling.
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
#' \strong{Ungated TabPFN-v2 weights.} The default TabPFN model version in
#' current \pkg{tabpfn} releases (V2.6) is license-gated and needs an interactive
#' token. This method pins the model to \code{ModelVersion.V2} -- the original
#' TabPFN-v2 (Hollmann et al., Nature 2025) -- which ships with a direct
#' (token-free) download option, so the method installs and runs non-interactively.
#'
#' \strong{Orientation and output.} The score is the frozen-context posterior
#' \eqn{P(\mathrm{case}) = P(y = 1)} for the query, in \eqn{[0, 1]}; larger =
#' more case-like. Degenerate queries (fewer than \code{hp$min_features} universe
#' features present, or an all-zero / empty specimen) return the neutral
#' probability \code{0.5}.
#'
#' @references
#' Hollmann N, Muller S, Purucker L, Krishnakumar A, Korfer M, Hoo SB, Schirrmeister
#' RT, Hutter F. (2025) Accurate predictions on small data with a tabular
#' foundation model. \emph{Nature} 637:319-326.
#'
#' Hollmann N, Muller S, Eggensperger K, Hutter F. (2023) TabPFN: A Transformer
#' That Solves Small Tabular Classification Problems in a Second.
#' \emph{International Conference on Learning Representations (ICLR)}.
#'
#' @name singlesample-tabpfn
NULL


# ----------------------------------------------------------------------------
# Single-sample transform (self-contained per-sample robust CLR)
# ----------------------------------------------------------------------------

# Robust CLR of one specimen over its OWN positive support. z(c * v) = z(v) for
# any c > 0 (exact per-sample scale-invariance); non-positive parts map to 0 and
# are treated as absent from the support. Deliberately does NOT reuse a package
# rclr helper that centres on a fixed pseudocount or cross-sample reference
# statistics (either would break exact scale-invariance / single-sample
# deployability). Mirrors .ecod_copod_rclr / .sinkhorn_single_rclr.
.tabpfn_rclr <- function(v) {
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
.tabpfn_rclr_matrix <- function(A) {
  Z <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) {
    .tabpfn_rclr(A[i, ])
  }))
  dimnames(Z) <- dimnames(A)
  Z
}

.tabpfn_score_chunks <- function(score_idx, chunk_size = 1024L) {
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

.tabpfn_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_tabpfn: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_tabpfn: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_tabpfn: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("device", "n_estimators", "max_context", "min_features", "seed",
               "score_batch")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_tabpfn: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  device <- hp$device
  if (is.null(device)) device <- "cpu"
  if (!is.character(device) || length(device) != 1L || is.na(device) ||
      !device %in% c("cpu", "cuda", "auto")) {
    stop("fit_tabpfn: hp$device must be one of 'cpu', 'cuda', 'auto'")
  }

  n_estimators <- hp$n_estimators
  if (is.null(n_estimators)) n_estimators <- 8L
  if (!is.numeric(n_estimators) || length(n_estimators) != 1L ||
      !is.finite(n_estimators) || n_estimators < 1L ||
      n_estimators > .Machine$integer.max ||
      n_estimators != as.integer(n_estimators)) {
    stop("fit_tabpfn: hp$n_estimators must be a positive integer")
  }

  max_context <- hp$max_context
  if (is.null(max_context)) max_context <- 4096L
  if (!is.numeric(max_context) || length(max_context) != 1L ||
      !is.finite(max_context) || max_context < 2L ||
      max_context > .Machine$integer.max ||
      max_context != as.integer(max_context)) {
    stop("fit_tabpfn: hp$max_context must be an integer >= 2")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_tabpfn: hp$min_features must be a positive integer")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 42L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("fit_tabpfn: hp$seed must be a single integer")
  }

  score_batch <- hp$score_batch
  if (is.null(score_batch)) score_batch <- FALSE
  if (!is.logical(score_batch) || length(score_batch) != 1L ||
      is.na(score_batch)) {
    stop("fit_tabpfn: hp$score_batch must be TRUE or FALSE")
  }
  list(
    device = device,
    n_estimators = as.integer(n_estimators),
    max_context = as.integer(max_context),
    min_features = as.integer(min_features),
    seed = as.integer(seed),
    score_batch = score_batch
  )
}


# ----------------------------------------------------------------------------
# Python / tabpfn availability + device + classifier construction
# ----------------------------------------------------------------------------

# Mirrors .ecod_copod_require_pyod / .tabdl_require_python_module: a clear,
# actionable error if reticulate or the tabpfn module is unavailable.
.tabpfn_require_module <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("fit_tabpfn: package 'reticulate' is required to run TabPFN. ",
         "Install with install.packages('reticulate').", call. = FALSE)
  }
  # Use the token-free direct-download path for the ungated V2 weights, never an
  # interactive browser login, and no network telemetry. Set only if unset so an
  # explicit caller override is respected.
  if (!nzchar(Sys.getenv("TABPFN_NO_BROWSER"))) {
    Sys.setenv(TABPFN_NO_BROWSER = "1")
  }
  if (!nzchar(Sys.getenv("TABPFN_DISABLE_TELEMETRY"))) {
    Sys.setenv(TABPFN_DISABLE_TELEMETRY = "1")
  }
  if (!reticulate::py_module_available("tabpfn")) {
    stop("fit_tabpfn: Python module 'tabpfn' is not available in the active ",
         "reticulate environment. Install with: pip install tabpfn ",
         "(in the configured RETICULATE_PYTHON venv).", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the requested device to a concrete 'cpu' / 'cuda'. 'auto' (and a
# 'cuda' request) fall back to 'cpu' when CUDA is unavailable, so the method
# never errors on a CPU-only host. Uses torch via reticulate.
.tabpfn_resolve_device <- function(device) {
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

# Build a TabPFN-v2 (ungated ModelVersion.V2) classifier with the frozen device,
# seed, and n_estimators, and fit it on the frozen rCLR context. Returns the
# live (non-serialisable) python classifier; it is NEVER part of the model
# digest. `Z_ctx` is the n_ctx x p numeric rCLR context, `y_ctx` the integer
# 0/1 labels. Saves/restores the torch RNG state around fit so the global torch
# RNG is not perturbed.
.tabpfn_build_and_fit <- function(Z_ctx, y_ctx, device, seed, n_estimators) {
  tabpfn <- reticulate::import("tabpfn", delay_load = FALSE)
  constants <- reticulate::import("tabpfn.constants", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  v2 <- constants$ModelVersion$V2

  clf <- tabpfn$TabPFNClassifier$create_default_for_version(
    v2,
    device = device,
    random_state = as.integer(seed),
    n_estimators = as.integer(n_estimators),
    ignore_pretraining_limits = TRUE
  )
  x_py <- np$asarray(Z_ctx, dtype = "float64")
  y_py <- np$asarray(as.integer(y_ctx), dtype = "int64")

  .tabpfn_with_torch_rng(seed, function() clf$fit(x_py, y_py))

  classes <- as.integer(reticulate::py_to_r(clf$classes_))
  case_col <- match(1L, classes)
  if (is.na(case_col)) {
    stop("fit_tabpfn: fitted TabPFN classes_ does not contain the case label 1",
         call. = FALSE)
  }
  list(clf = clf, classes = classes, case_col = case_col)
}

# Run `fn()` with the torch RNG seeded to `seed`, restoring the torch RNG state
# afterwards so the global torch RNG is left byte-identical. TabPFN inference is
# deterministic regardless of the torch seed (verified), but `predict`/`fit`
# advance the torch RNG as a side effect; this keeps that side effect contained.
.tabpfn_with_torch_rng <- function(seed, fn) {
  torch_py <- tryCatch(reticulate::import("torch", delay_load = FALSE),
                       error = function(e) NULL)
  if (is.null(torch_py)) return(fn())
  state <- tryCatch(torch_py$random$get_rng_state(), error = function(e) NULL)
  # Also preserve the CUDA RNG state(s): on a GPU fit, TabPFN advances the CUDA
  # generator, which the CPU-only get/set_rng_state does NOT cover. Guard on CUDA
  # availability so this is a no-op on CPU-only hosts.
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
# fit_tabpfn
# ----------------------------------------------------------------------------

#' @title Fit the TabPFN-v2 in-context discriminator (freeze the training context)
#'
#' @description
#' Freezes the training context for TabPFN-v2 in-context classification. The
#' training matrix is mapped to the per-sample robust CLR over the frozen feature
#' universe (\code{colnames(X_train)}); if the context exceeds
#' \code{hp$max_context} rows it is reduced to a deterministic, seeded,
#' class-stratified subsample. A TabPFN-v2 (ungated \code{ModelVersion.V2})
#' classifier is constructed and \code{fit()} on that frozen rCLR context (which,
#' for TabPFN, simply loads the context into the model -- no gradient training).
#' The frozen R-side state needed to reproduce scores (the rCLR context, labels,
#' feature universe, class order, device/config, hp) is stored on the model; the
#' live python classifier is also kept for scoring but is NEVER part of the model
#' digest.
#' TabPFN-v2 is applied to the full rCLR feature universe; cohorts with
#' \eqn{>500} features, including approximately 2540-2565-feature high-plex
#' miRNA panels, exceed TabPFN's nominal 500-feature ceiling, so its documented
#' \code{ignore_pretraining_limits} override is enabled. Under the pinned
#' \code{ModelVersion.V2} preprocessor path (tabpfn 8.0.7) the per-estimator
#' feature-subsampling threshold is 1,000,000 (the 500-feature random-subspace
#' reduction is the separate V2.5 preset, which is not used here), so TabPFN-v2
#' ingests the full feature set directly and is run outside its nominal
#' \eqn{\le 500}-feature pretraining regime. No per-estimator feature subsampling
#' is introduced (the ensemble subsample indices are all empty), so scoring stays
#' deterministic and exactly row-equivariant. The override also lifts the CPU
#' \eqn{>1000}-sample guard (a compute guard, not a numeric change). This matches
#' how the sibling in-context foundation models TabICL and TabDPT are run on the
#' full representation and keeps the comparison apples-to-apples with the rCLR
#' baseline; a leakage-safe in-regime top-500 feature reduction is evaluated as
#' a companion sensitivity analysis (the prespecified out-of-regime control).
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with unique, non-empty feature names (colnames).
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with at least one case and one control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and otherwise ignored.
#' @param hp Optional list of hyperparameters. Allowed fields: \code{device}
#'   (\code{"cpu"} (default), \code{"cuda"}, or \code{"auto"}; \code{"cuda"} /
#'   \code{"auto"} fall back to CPU when no GPU is available), \code{n_estimators}
#'   (TabPFN ensemble size, positive integer, default \code{8L}),
#'   \code{max_context} (row cap on the frozen context, integer \eqn{\ge 2},
#'   default \code{4096L}; larger contexts are seeded class-stratified
#'   subsampled), \code{min_features} (feature-overlap floor at scoring, positive
#'   integer, default \code{3L}), \code{seed} (integer; default \code{42L}), and
#'   \code{score_batch} (benchmark-only batched scoring flag, logical, default
#'   \code{FALSE}).
#'
#' @return Object of class \code{tabpfn_model}: a list with \code{context_rclr}
#'   (frozen rCLR context matrix), \code{context_y}, \code{feature_universe},
#'   \code{classes}, \code{case_col}, \code{n_context}, \code{device} (resolved),
#'   \code{model_version} (\code{"v2"}),
#'   \code{ignore_pretraining_limits} (\code{TRUE}), \code{hp}, and the live
#'   python classifier \code{clf} (excluded from the digest).
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
#' model <- fit_tabpfn(X, y)
#' score_tabpfn(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Hollmann N, et al. (2025) Accurate predictions on small data with a tabular
#' foundation model. \emph{Nature} 637:319-326.
#'
#' @export
fit_tabpfn <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_tabpfn", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_tabpfn", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_tabpfn")
  hp <- .tabpfn_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_tabpfn: X_train must contain at least hp$min_features features")
  }
  .tabpfn_require_module()

  feat_order <- colnames(X_train)
  Z <- .tabpfn_rclr_matrix(X_train)              # n x p frozen rCLR context

  # Deterministic, seeded, class-stratified context cap. The global .Random.seed
  # is saved and restored so the subsample uses an isolated RNG stream.
  keep <- .tabpfn_context_subsample(nrow(Z), y, hp$max_context, hp$seed)
  Z_ctx <- Z[keep, , drop = FALSE]
  y_ctx <- y[keep]
  if (!any(y_ctx == 1L) || !any(y_ctx == 0L)) {
    stop("fit_tabpfn: capped context lost a class; raise hp$max_context")
  }

  device <- .tabpfn_resolve_device(hp$device)
  fitted <- .tabpfn_build_and_fit(Z_ctx, y_ctx, device, hp$seed, hp$n_estimators)

  model <- list(
    context_rclr = Z_ctx,
    context_y = y_ctx,
    feature_universe = feat_order,
    classes = fitted$classes,
    case_col = fitted$case_col,
    n_context = nrow(Z_ctx),
    device = device,
    model_version = "v2",
    ignore_pretraining_limits = TRUE,
    hp = hp,
    clf = fitted$clf
  )
  class(model) <- "tabpfn_model"
  model
}

# Deterministic seeded class-stratified subsample of row indices to at most
# `cap` rows, preserving the original row order of the kept set. Returns the
# full 1:n when n <= cap. Saves/restores the global .Random.seed so it does not
# perturb the caller's RNG stream.
.tabpfn_context_subsample <- function(n, y, cap, seed) {
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
  # When both classes are present and cap >= 2, RESERVE at least one slot for the
  # control class so the capped context never drops a class. Without this, a tiny
  # cap with extreme imbalance (e.g. cap = 2 over 99 cases + 1 control) rounds the
  # whole cap to the majority class and loses the control class.
  if (length(idx_ctrl) > 0L && cap >= 2L) n_case <- min(n_case, cap - 1L)
  n_ctrl <- cap - n_case
  if (n_ctrl > length(idx_ctrl)) {
    n_ctrl <- length(idx_ctrl)
    n_case <- min(length(idx_case), cap - n_ctrl)
  }
  if (n_case > length(idx_case)) n_case <- length(idx_case)
  # Sample POSITIONS within each class index vector, never the values directly:
  # base sample(x, k) treats a length-1 x as 1:x (the classic gotcha), so
  # sample(idx_ctrl, 1) with a single control index would pick a wrong row.
  pick <- function(idx, k) idx[sample.int(length(idx), k)]
  keep <- c(
    if (n_case > 0L) pick(idx_case, n_case) else integer(0),
    if (n_ctrl > 0L) pick(idx_ctrl, n_ctrl) else integer(0)
  )
  sort(keep)
}


# ----------------------------------------------------------------------------
# score_tabpfn
# ----------------------------------------------------------------------------

#' @title Score the TabPFN-v2 in-context discriminator (default row-by-row)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_tabpfn}}. Each query is mapped to the per-sample robust CLR
#' over the frozen feature universe (absent universe features carry the neutral
#' rCLR value \code{0}), then classified by TabPFN against the FROZEN training
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
#' leakage-free and AUC-faithful (observed \eqn{|dAUC| = 7.8\mathrm{e}{-5}} at
#' production scale; specimen ranks/verdicts unchanged), but is not bit-identical
#' to the default \eqn{n = 1} row-by-row path because of a deterministic
#' batch-size numerical effect. The default row-by-row single-specimen deployment
#' path is unchanged.
#'
#' @param model A \code{tabpfn_model} object from \code{\link{fit_tabpfn}}.
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
#' model <- fit_tabpfn(X, y)
#' score_tabpfn(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Hollmann N, et al. (2025) Accurate predictions on small data with a tabular
#' foundation model. \emph{Nature} 637:319-326.
#'
#' @export
score_tabpfn <- function(model, X, meta = NULL) {
  if (!inherits(model, "tabpfn_model")) {
    stop("score_tabpfn: model must have class tabpfn_model")
  }
  X <- .reo_check_matrix(X, "score_tabpfn", "X")
  .reo_check_meta(meta, nrow(X), "score_tabpfn", "meta")

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

  if (isTRUE(model$hp$score_batch)) {
    Z <- .tabpfn_rclr_matrix(X_use)
    score_idx <- which(rowSums(Z != 0) > 0L)          # all-zero rCLR stays neutral
    if (length(score_idx) > 0L) {
      .tabpfn_with_torch_rng(model$hp$seed, function() {
        for (idx in .tabpfn_score_chunks(score_idx, chunk_size = 1024L)) {
          pr <- reticulate::py_to_r(
            clf$predict_proba(np$asarray(Z[idx, , drop = FALSE], dtype = "float64"))
          )
          out[idx] <<- as.numeric(pr[, case_col])
        }
      })
    }
  } else {
    # Forced row-by-row predict_proba: the n = 1 forward path is used uniformly.
    # The torch RNG is restored around the whole loop (predict advances it but
    # the result is deterministic regardless of the torch seed).
    .tabpfn_with_torch_rng(model$hp$seed, function() {
      for (i in seq_len(nrow(X_use))) {
        z <- .tabpfn_rclr(X_use[i, ])
        if (all(z == 0)) next                          # all-zero specimen: neutral
        zi <- matrix(z, nrow = 1L)
        pr <- reticulate::py_to_r(clf$predict_proba(np$asarray(zi, dtype = "float64")))
        out[i] <<- as.numeric(pr[1, case_col])
      }
    })
  }

  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out)) ||
      any(out < 0 | out > 1)) {
    stop("score_tabpfn: scorer produced invalid output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen score-determining R-side state (rCLR training context,
# labels, feature universe, class order, device/config, hp). The live python
# TabPFN classifier (model$clf, a non-serialisable external pointer) is
# deliberately EXCLUDED, so the digest is a deterministic snapshot suitable as
# the `model_digest` for singlesample_assert_row_equivariant() (clause (d): the
# bytes must be identical before and after scoring). Scoring never mutates this
# frozen state (TabPFN inference is read-only over the frozen context).
.tabpfn_model_digest <- function(model) {
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
    model_version = model$model_version,
    hp = hp
  )
  digest::digest(state, algo = "xxhash64")
}
