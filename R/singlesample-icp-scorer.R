#' @title Fit Invariant Causal Prediction discriminator (transfer at training,
#'   single-sample at inference)
#'
#' @description
#' Learns a frozen logistic discriminator over the subset of robust-CLR (rCLR)
#' features that are both PREDICTIVE of the label and INVARIANT (cohort-stable)
#' across the training environments (cohorts) named by
#' \code{meta_train[[cohort_col]]}, following the Invariant Causal Prediction
#' (ICP) principle of Peters, Buhlmann & Meinshausen (2016): a predictor whose
#' conditional relationship with the label is INVARIANT across environments is a
#' plausible (batch-robust, transferable) causal signal, whereas a predictor
#' whose effect VARIES by cohort signals confounding / non-causal, batch-driven
#' association and is dropped. This makes \code{fit_icp} a transfer-estimand
#' method at TRAINING (selection uses the cross-environment structure), but its
#' deployment is fully SINGLE-SAMPLE: the frozen model is a linear logistic score
#' over the selected features evaluated on one specimen's own rCLR vector. The
#' transfer aspect is entirely in the SELECTION.
#'
#' icp is the structural twin of \code{\link{fit_sel_stablemate}} with ONE
#' change: the feature-SELECTION step. Where sel-stablemate calls StableMate's
#' predictivity+stability bootstrap ensembles, icp selects features by an
#' explicit, tractable per-feature ICP invariance screen:
#' \enumerate{
#'   \item PREDICTIVITY screen: a pooled marginal logistic
#'     \code{glm(y ~ rclr_j, binomial)} over ALL training rows; feature \eqn{j}
#'     is PREDICTIVE if the Wald p-value of its slope \eqn{< \alpha_{pred}}
#'     (default \code{0.05}).
#'   \item INVARIANCE screen (the ICP core): a per-cohort marginal logistic
#'     \code{glm(y ~ rclr_j, binomial)} WITHIN each VALID cohort (\eqn{\ge 1}
#'     case AND \eqn{\ge 1} control AND \eqn{\ge} \code{min_env_n} rows), giving
#'     per-cohort slopes \eqn{b_j^{(e)}} with standard errors \eqn{SE_e}. The
#'     cross-cohort HOMOGENEITY of \eqn{\{b_j^{(e)}\}} is tested with a Cochran's
#'     fixed-effect heterogeneity statistic \eqn{Q = \sum_e w_e (b_e -
#'     \bar{b})^2}, \eqn{w_e = 1/SE_e^2}, \eqn{\bar{b} = \sum_e w_e b_e / \sum_e
#'     w_e}, with \eqn{Q \sim \chi^2_{E-1}} under homogeneity. Feature \eqn{j} is
#'     INVARIANT if the Q p-value \eqn{> \alpha_{inv}} (default \code{0.10}) ---
#'     i.e. we FAIL to reject homogeneity, so the \eqn{y \mid x_j} relationship
#'     is stable across cohorts.
#'   \item SELECT feature \eqn{j} if PREDICTIVE AND INVARIANT.
#' }
#' A frozen \code{glm(y ~ ., binomial)} is then fit over the rCLR of the SELECTED
#' features, pooled over all training rows, and the intercept + per-feature
#' coefficients are stored. The per-cohort glm sweep can touch the RNG only
#' through deterministic IRLS (it does not), but the whole fit is nevertheless
#' wrapped in a global-\code{.Random.seed} save/restore guard with
#' \code{set.seed(hp$seed)} inside it, so fitting is deterministic and leaves the
#' caller's RNG state byte-unchanged --- mirroring the sel-stablemate contract.
#'
#' DEGENERATION (icp ALWAYS produces a usable, non-empty head):
#' \itemize{
#'   \item If \code{meta_train} is NULL, the cohort column is missing, or there
#'     are \eqn{< 2} VALID cohorts, the invariance screen cannot be evaluated and
#'     selection degrades GRACEFULLY to the PREDICTIVITY-ONLY screen (the pooled
#'     marginal logistic alone). \code{n_environments == 1L} records this so a
#'     caller that intends cross-cohort transfer can audit whether a mis-wired
#'     \code{meta_train} silently demoted this transfer method to pooled
#'     predictivity-only selection. \code{n_environments > 1L} means the ICP
#'     invariance screen ran.
#'   \item If no feature passes both screens, the selection RELAXES (documented in
#'     \code{model$selection_mode}): take the most-homogeneous (highest Q
#'     p-value) among the predictive features, up to \code{top_k} (default
#'     \code{min(5, n_predictive)}); if still empty (no predictive feature at
#'     all), take the single most-predictive feature. This guarantees a usable
#'     non-empty head while keeping the relaxation auditable.
#' }
#' Unlike sel-stablemate (a stability method that honestly returns the EMPTY set
#' when nothing is stable), icp's brief mandates an always-non-empty head; the
#' relaxation taken is therefore recorded in \code{model$selection_mode} rather
#' than left empty.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. When it contains
#'   \code{cohort_col}, non-missing values define the training environments used
#'   for the ICP invariance screen.
#' @param hp List of frozen hyperparameters: \code{cohort_col} (cohort/environment
#'   column in \code{meta_train}, default \code{"accession"}); \code{max_features}
#'   (candidate cap by training variance, default \code{50L}, integer
#'   \eqn{\ge 2}); \code{alpha_pred} (predictivity Wald p-value threshold, default
#'   \code{0.05}, in \code{(0, 1)}); \code{alpha_inv} (invariance Cochran's-Q
#'   p-value threshold; a feature is invariant when its Q p-value EXCEEDS this,
#'   default \code{0.10}, in \code{(0, 1)}); \code{min_env_n} (minimum rows for a
#'   cohort to be a valid environment, default \code{8L}, integer \eqn{\ge 4});
#'   \code{top_k} (relaxation cap on most-homogeneous-among-predictive features
#'   when the strict set is empty, default \code{5L}, integer \eqn{\ge 1});
#'   \code{min_selected} (minimum present selected features in a specimen for a
#'   non-neutral score, default \code{1L}, integer \eqn{\ge 1}); and \code{seed}
#'   (default \code{1L}, used to restore the global RNG state). Resolved once at
#'   fit and not tuned inside \code{fit_icp()}.
#'
#' @return A plain list of class \code{icp_model} containing the frozen
#'   \code{selected_features} (the predictive-and-invariant rCLR feature ids), the
#'   logistic \code{intercept} and named \code{coefficients}, the candidate
#'   \code{feature_universe} (the frozen rCLR universe), \code{n_environments}
#'   (the auditable meta-engagement signal: \code{> 1L} = ICP invariance screen
#'   ran, \code{== 1L} = degenerated predictivity-only), \code{selection_mode}
#'   (one of \code{"invariant_predictive"}, \code{"predictivity_only"},
#'   \code{"relaxed_top_k"}, \code{"relaxed_single_most_predictive"}, recording
#'   any relaxation taken), \code{cohort_col}, \code{min_selected}, a
#'   \code{degenerate} flag (\code{TRUE} only if NO feature could be selected even
#'   after relaxation, in which case scoring returns the neutral \code{0}), and
#'   the resolved \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 180L; p <- 12L
#' env <- rep(c("e1", "e2", "e3"), length.out = n)
#' X <- matrix(stats::rgamma(n * p, shape = 5), n, p,
#'             dimnames = list(NULL, paste0("f", seq_len(p))))
#' lin <- log(X[, 1]) + log(X[, 2]) + log(X[, 3])
#' y <- stats::rbinom(n, 1, plogis(scale(lin)))
#' meta <- data.frame(accession = env, stringsAsFactors = FALSE)
#' model <- fit_icp(X, y, meta_train = meta)
#' score <- score_icp(model, X)
#' }
#'
#' @references
#' Peters J, Buhlmann P, Meinshausen N. (2016) Causal inference by using
#' invariant prediction: identification and confidence intervals. \emph{Journal
#' of the Royal Statistical Society: Series B} 78(5): 947-1012.
#'
#' Cochran WG. (1954) The combination of estimates from different experiments.
#' \emph{Biometrics} 10(1): 101-129.
#'
#' Aitchison J. (1986) \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' @export
fit_icp <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_icp", "X_train")
  if (ncol(X_train) < 2L) {
    stop("fit_icp: X_train must contain at least two features")
  }
  .reo_check_meta(meta_train, nrow(X_train), "fit_icp", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_icp")
  hp <- .icp_resolve_hp(hp, ncol(X_train))

  # Frozen candidate universe: top-max_features by training variance, name-tied
  # deterministically. The rCLR feature matrix is per-sample (within-sample),
  # so this universe is the only training-derived structure besides the model.
  universe <- .icp_universe(X_train, hp$max_features)
  Z_train <- .icp_rclr_matrix(X_train[, universe, drop = FALSE])

  env_info <- .icp_environment(meta_train, y, hp$cohort_col, hp$min_env_n)

  sel <- .icp_select(Z_train, y, env_info, hp$alpha_pred, hp$alpha_inv,
                     hp$top_k, hp$seed)
  selected <- sel$selected
  selection_mode <- sel$selection_mode

  degenerate <- length(selected) == 0L
  intercept <- 0
  coefficients <- numeric(0)
  if (!degenerate) {
    cal <- .icp_calibrate(Z_train, y, selected)
    intercept <- cal$intercept
    coefficients <- cal$coefficients
    selected <- names(coefficients)
    degenerate <- length(selected) == 0L
    if (degenerate) selection_mode <- "degenerate"
  }

  model <- list(
    selected_features = unname(selected),
    intercept = intercept,
    coefficients = coefficients,
    feature_universe = universe,
    n_environments = env_info$n_environments,
    selection_mode = selection_mode,
    cohort_col = hp$cohort_col,
    min_selected = hp$min_selected,
    degenerate = degenerate,
    hp = hp
  )
  class(model) <- "icp_model"
  model
}


#' @title Score Invariant Causal Prediction discriminator (single-sample)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen
#' invariant-and-predictive logistic learned by \code{\link{fit_icp}}. For one
#' specimen its per-sample rCLR is computed over the frozen candidate universe
#' (present features only), restricted to the frozen selected features, and the
#' frozen logistic linear predictor
#' \deqn{\eta = intercept + \sum_{j \in S} \beta_j \, z_j}
#' is returned (NOT \code{plogis()} of it); larger is more case-like. Here
#' \eqn{z_j} is the specimen's own rCLR coordinate for selected feature \eqn{j}.
#' Because rCLR centres each specimen by its own geometric mean over present
#' parts, the score is EXACTLY invariant to per-sample positive scaling and uses
#' only that specimen's own values, so the method is single-sample deployable and
#' row-equivariant. Scoring uses no random numbers and no scored-batch statistics.
#'
#' Missing selected features contribute \code{0} (their rCLR coordinate is absent
#' from the specimen): the linear predictor is summed only over selected features
#' present in \code{X}. If fewer than \code{model$min_selected} selected features
#' are present in a specimen, or the model selected nothing at fit time
#' (\code{model$degenerate}), that specimen receives the documented neutral score
#' \code{0}.
#'
#' @param model An \code{icp_model} object returned by \code{\link{fit_icp}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method (the transfer/environment information is used
#'   only at fit time).
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Scores use only
#'   each row's own values plus the frozen model, with no scored-batch centering,
#'   quantiles, or renormalization.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 180L; p <- 12L
#' env <- rep(c("e1", "e2", "e3"), length.out = n)
#' X <- matrix(stats::rgamma(n * p, shape = 5), n, p,
#'             dimnames = list(NULL, paste0("f", seq_len(p))))
#' lin <- log(X[, 1]) + log(X[, 2]) + log(X[, 3])
#' y <- stats::rbinom(n, 1, plogis(scale(lin)))
#' meta <- data.frame(accession = env, stringsAsFactors = FALSE)
#' model <- fit_icp(X, y, meta_train = meta)
#' score_icp(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Peters J, Buhlmann P, Meinshausen N. (2016) Causal inference by using
#' invariant prediction: identification and confidence intervals. \emph{Journal
#' of the Royal Statistical Society: Series B} 78(5): 947-1012.
#'
#' @export
score_icp <- function(model, X, meta = NULL) {
  if (!inherits(model, "icp_model")) {
    stop("score_icp: model must have class icp_model")
  }
  X <- .reo_check_matrix(X, "score_icp", "X")
  .reo_check_meta(meta, nrow(X), "score_icp", "meta")

  if (isTRUE(model$degenerate) || length(model$selected_features) == 0L) {
    return(rep(0, nrow(X)))
  }
  present <- intersect(model$selected_features, colnames(X))
  if (length(present) < model$min_selected) {
    return(rep(0, nrow(X)))
  }

  out <- vapply(seq_len(nrow(X)), function(i) {
    # Restore feature names: single-row slicing of a 1-column matrix drops them,
    # which would otherwise make the rCLR universe / present-feature selection
    # depend on batch shape.
    x_i <- X[i, ]
    names(x_i) <- colnames(X)
    .icp_score_row(
      x = x_i,
      universe = model$feature_universe,
      selected = model$selected_features,
      coefficients = model$coefficients,
      intercept = model$intercept,
      min_selected = model$min_selected
    )
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_icp: scorer produced non-finite or wrong-length output")
  }
  out
}

.icp_resolve_hp <- function(hp, p) {
  if (!is.list(hp)) stop("fit_icp: hp must be a list")
  allowed <- c("cohort_col", "max_features", "alpha_pred", "alpha_inv",
               "min_env_n", "top_k", "min_selected", "seed")
  nm <- names(hp)
  if (length(hp) > 0L && (is.null(nm) || any(!nzchar(nm)))) {
    stop("fit_icp: all hp fields must be named")
  }
  if (any(duplicated(nm))) {
    stop("fit_icp: duplicate hp field(s): ",
         paste(unique(nm[duplicated(nm)]), collapse = ", "))
  }
  unknown <- setdiff(nm, allowed)
  if (length(unknown) > 0L) {
    stop("fit_icp: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }

  cohort_col <- hp[["cohort_col"]]
  if (is.null(cohort_col)) cohort_col <- "accession"
  if (!is.character(cohort_col) || length(cohort_col) != 1L ||
      is.na(cohort_col) || !nzchar(cohort_col)) {
    stop("fit_icp: hp$cohort_col must be a single non-empty character string")
  }

  max_features <- hp[["max_features"]]
  if (is.null(max_features)) max_features <- 50L
  if (!is.numeric(max_features) || length(max_features) != 1L ||
      !is.finite(max_features) || max_features < 2L ||
      max_features > .Machine$integer.max ||
      max_features != as.integer(max_features)) {
    stop("fit_icp: hp$max_features must be an integer >= 2")
  }
  max_features <- min(as.integer(max_features), p)

  alpha_pred <- hp[["alpha_pred"]]
  if (is.null(alpha_pred)) alpha_pred <- 0.05
  if (!is.numeric(alpha_pred) || length(alpha_pred) != 1L ||
      !is.finite(alpha_pred) || alpha_pred <= 0 || alpha_pred >= 1) {
    stop("fit_icp: hp$alpha_pred must be a number in (0, 1)")
  }

  alpha_inv <- hp[["alpha_inv"]]
  if (is.null(alpha_inv)) alpha_inv <- 0.10
  if (!is.numeric(alpha_inv) || length(alpha_inv) != 1L ||
      !is.finite(alpha_inv) || alpha_inv <= 0 || alpha_inv >= 1) {
    stop("fit_icp: hp$alpha_inv must be a number in (0, 1)")
  }

  min_env_n <- hp[["min_env_n"]]
  if (is.null(min_env_n)) min_env_n <- 8L
  if (!is.numeric(min_env_n) || length(min_env_n) != 1L ||
      !is.finite(min_env_n) || min_env_n < 4L ||
      min_env_n > .Machine$integer.max ||
      min_env_n != as.integer(min_env_n)) {
    stop("fit_icp: hp$min_env_n must be an integer >= 4")
  }

  top_k <- hp[["top_k"]]
  if (is.null(top_k)) top_k <- 5L
  if (!is.numeric(top_k) || length(top_k) != 1L || !is.finite(top_k) ||
      top_k < 1L || top_k > .Machine$integer.max ||
      top_k != as.integer(top_k)) {
    stop("fit_icp: hp$top_k must be a positive integer")
  }

  min_selected <- hp[["min_selected"]]
  if (is.null(min_selected)) min_selected <- 1L
  if (!is.numeric(min_selected) || length(min_selected) != 1L ||
      !is.finite(min_selected) || min_selected < 1L ||
      min_selected > .Machine$integer.max ||
      min_selected != as.integer(min_selected)) {
    stop("fit_icp: hp$min_selected must be a positive integer")
  }

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_icp: hp$seed must be a non-negative integer")
  }

  list(
    cohort_col = cohort_col,
    max_features = max_features,
    alpha_pred = as.numeric(alpha_pred),
    alpha_inv = as.numeric(alpha_inv),
    min_env_n = as.integer(min_env_n),
    top_k = as.integer(top_k),
    min_selected = as.integer(min_selected),
    seed = as.integer(seed)
  )
}

# Frozen candidate universe: the top-max_features features by training variance,
# with feature-name radix tie-break so the universe is reproducible across
# machines/locales.
.icp_universe <- function(X_train, max_features) {
  v <- apply(X_train, 2L, stats::var)
  v[!is.finite(v)] <- -Inf
  ord <- order(-v, colnames(X_train), method = "radix")
  colnames(X_train)[ord][seq_len(max_features)]
}

# Resolve the training environment factor from meta and record the auditable
# n_environments engagement signal. A VALID environment must contain at least
# one case AND one control AND at least min_env_n rows (the per-cohort marginal
# logistic needs both classes and enough rows for a usable slope/SE). When fewer
# than two valid environments exist, selection degrades to predictivity-only and
# n_environments == 1L.
.icp_environment <- function(meta_train, y, cohort_col, min_env_n) {
  pooled <- list(env = NULL, valid_levels = character(0), n_environments = 1L)
  if (is.null(meta_train)) return(pooled)
  meta_train <- as.data.frame(meta_train, stringsAsFactors = FALSE)
  if (!cohort_col %in% names(meta_train)) return(pooled)

  cohort <- as.character(meta_train[[cohort_col]])
  non_missing <- !is.na(cohort)
  levels_all <- unique(cohort[non_missing])
  valid <- levels_all[vapply(levels_all, function(lv) {
    rows <- non_missing & cohort == lv
    sum(rows) >= min_env_n && any(y[rows] == 1L) && any(y[rows] == 0L)
  }, logical(1))]
  if (length(valid) <= 1L) return(pooled)

  list(env = cohort, valid_levels = sort(valid),
       n_environments = length(valid))
}

# ICP feature selection under a global-.Random.seed save/restore guard so the
# selection leaves the caller's RNG state byte-unchanged (the glm IRLS used here
# consumes no RNG, but the guard mirrors the sel-stablemate contract and is
# robust to any future RNG-touching internals). Returns the selected feature
# names and the selection_mode tag.
#
# With >= 2 valid environments the full ICP screen runs (predictivity AND
# Cochran's-Q invariance). With a single environment (degeneration) the
# invariance screen is undefined, so selection gracefully degrades to the
# predictivity-only marginal screen. icp ALWAYS yields a non-empty head: if the
# strict screen selects nothing, it relaxes to the most-homogeneous-among-
# predictive top_k, then to the single most-predictive feature.
.icp_select <- function(Z_train, y, env_info, alpha_pred, alpha_inv,
                        top_k, seed) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)

  feats <- colnames(Z_train)

  # Per-feature pooled predictivity p-value (Wald p of the marginal slope).
  p_pred <- vapply(feats, function(f) {
    .icp_pred_pvalue(Z_train[, f], y)
  }, numeric(1))
  predictive <- feats[is.finite(p_pred) & p_pred < alpha_pred]

  # Single-environment degeneration: invariance is undefined, use predictivity
  # only. Still guarantee a non-empty head via the relaxation ladder.
  if (env_info$n_environments <= 1L) {
    if (length(predictive) > 0L) {
      return(list(selected = predictive, selection_mode = "predictivity_only"))
    }
    return(.icp_relax(feats, p_pred, q_pval = NULL, predictive = predictive,
                      top_k = top_k,
                      mode_topk = "relaxed_top_k",
                      mode_single = "relaxed_single_most_predictive"))
  }

  # Per-feature Cochran's-Q invariance p-value over the valid cohorts. A larger Q
  # p-value = MORE homogeneous (invariant) across cohorts. A feature whose
  # per-cohort slopes cannot be estimated in >= 2 cohorts gets NA (treated as
  # NON-invariant: we cannot demonstrate invariance, so we do not claim it).
  q_pval <- vapply(feats, function(f) {
    .icp_invariance_pvalue(Z_train[, f], y, env_info$env, env_info$valid_levels)
  }, numeric(1))

  invariant <- feats[is.finite(q_pval) & q_pval > alpha_inv]
  selected <- intersect(predictive, invariant)
  if (length(selected) > 0L) {
    return(list(selected = selected,
                selection_mode = "invariant_predictive"))
  }

  # Relaxation: among predictive features, keep the most homogeneous (highest Q
  # p-value) up to top_k; if no predictive feature, the single most-predictive.
  .icp_relax(feats, p_pred, q_pval, predictive, top_k,
             mode_topk = "relaxed_top_k",
             mode_single = "relaxed_single_most_predictive")
}

# Relaxation ladder shared by both the single-environment and multi-environment
# empty-selection cases. With predictive features available, pick the top_k by
# homogeneity (highest Q p-value; when q_pval is NULL/all-NA, fall back to most
# predictive). With no predictive feature, pick the single most-predictive
# feature (lowest predictivity p-value). Ties broken by feature name (radix) for
# reproducibility.
.icp_relax <- function(feats, p_pred, q_pval, predictive, top_k,
                       mode_topk, mode_single) {
  if (length(predictive) > 0L) {
    if (!is.null(q_pval)) {
      qp <- q_pval[predictive]
      # Most homogeneous first: highest Q p-value. NA homogeneity -> least
      # preferred. Tie-break by feature name.
      key <- ifelse(is.finite(qp), qp, -Inf)
      ord <- order(-key, predictive, method = "radix")
    } else {
      pp <- p_pred[predictive]
      key <- ifelse(is.finite(pp), pp, Inf)
      ord <- order(key, predictive, method = "radix")
    }
    k <- min(as.integer(top_k), length(predictive))
    return(list(selected = predictive[ord][seq_len(k)],
                selection_mode = mode_topk))
  }
  # No predictive feature at all: take the single most-predictive (lowest p).
  pp <- ifelse(is.finite(p_pred), p_pred, Inf)
  ord <- order(pp, feats, method = "radix")
  list(selected = feats[ord][1L], selection_mode = mode_single)
}

# Pooled marginal-logistic predictivity p-value: Wald p of the slope of
# glm(y ~ z_j, binomial) over ALL training rows. A constant/degenerate column or
# a glm failure -> NA (not predictive).
.icp_pred_pvalue <- function(z, y) {
  if (length(unique(z[is.finite(z)])) < 2L) return(NA_real_)
  df <- data.frame(.y = y, z = z)
  fit <- tryCatch(
    suppressWarnings(stats::glm(.y ~ z, data = df, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  sm <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(sm) || !("z" %in% rownames(sm))) return(NA_real_)
  p <- sm["z", "Pr(>|z|)"]
  if (!is.finite(p)) return(NA_real_)
  p
}

# ICP invariance p-value for one feature: per-cohort marginal-logistic slopes
# b_e with standard errors SE_e over the VALID cohorts, combined via Cochran's
# fixed-effect heterogeneity Q = sum_e w_e (b_e - b_bar)^2, w_e = 1/SE_e^2,
# b_bar = sum w_e b_e / sum w_e, Q ~ chi^2_{E-1}. Returns the homogeneity
# p-value = pchisq(Q, E-1, lower.tail = FALSE): LARGE p = homogeneous/invariant.
# A cohort whose slope/SE is not finite (separation, constant column, glm
# failure) is dropped; if fewer than 2 cohorts yield a usable slope, invariance
# cannot be assessed -> NA (treated as non-invariant by the caller).
.icp_invariance_pvalue <- function(z, y, env, valid_levels) {
  b <- numeric(0)
  se <- numeric(0)
  for (lv in valid_levels) {
    rows <- !is.na(env) & env == lv
    zl <- z[rows]
    yl <- y[rows]
    if (length(unique(zl[is.finite(zl)])) < 2L ||
        !any(yl == 1L) || !any(yl == 0L)) {
      next
    }
    df <- data.frame(.y = yl, z = zl)
    fit <- tryCatch(
      suppressWarnings(stats::glm(.y ~ z, data = df, family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    sm <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(sm) || !("z" %in% rownames(sm))) next
    bj <- sm["z", "Estimate"]
    sej <- sm["z", "Std. Error"]
    if (!is.finite(bj) || !is.finite(sej) || sej <= 0) next
    b <- c(b, bj)
    se <- c(se, sej)
  }
  E <- length(b)
  if (E < 2L) return(NA_real_)
  w <- 1 / (se^2)
  b_bar <- sum(w * b) / sum(w)
  Q <- sum(w * (b - b_bar)^2)
  if (!is.finite(Q) || Q < 0) return(NA_real_)
  stats::pchisq(Q, df = E - 1L, lower.tail = FALSE)
}

# Frozen logistic calibration over the rCLR of the SELECTED features.
# glm(y ~ ., binomial); store intercept + coefficients keyed by feature name.
# Non-finite/dropped coefficients (rank deficiency, perfect separation) are
# removed so only finite linear terms are frozen. If glm fails or no finite
# coefficient survives, the model degenerates (empty selection -> neutral 0).
# (Identical to the sel-stablemate calibration; the only icp departure is the
# selection step upstream.)
.icp_calibrate <- function(Z_train, y, selected) {
  present <- intersect(selected, colnames(Z_train))
  if (length(present) == 0L) {
    return(list(intercept = 0, coefficients = numeric(0)))
  }
  # Use safe placeholder column names (V1, V2, ...) so glm cannot mangle feature
  # ids that contain non-syntactic characters (e.g. "miR-4" -> "`miR-4`"); map
  # the fitted coefficients back to the original feature ids strictly BY POSITION
  # (glm preserves design-column order), avoiding any name-quoting ambiguity.
  Zsel <- Z_train[, present, drop = FALSE]
  safe <- paste0("V", seq_along(present))
  df <- as.data.frame(Zsel)
  names(df) <- safe
  df$.y <- y
  fit <- tryCatch(
    suppressWarnings(stats::glm(.y ~ ., data = df, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(intercept = 0, coefficients = numeric(0)))
  }
  co <- stats::coef(fit)
  intercept <- unname(co[["(Intercept)"]])
  # Coefficients are named V1..Vk in the same order as `present`; re-key to the
  # original feature ids by position (a dropped/aliased term shows up as NA and
  # is removed below).
  beta <- stats::setNames(unname(co[safe]), present)
  beta <- beta[is.finite(beta) & beta != 0]
  if (!is.finite(intercept) || length(beta) == 0L) {
    return(list(intercept = 0, coefficients = numeric(0)))
  }
  list(intercept = intercept, coefficients = beta)
}

# Self-contained per-sample robust CLR of one specimen. The geometric-mean
# centring uses ONLY this specimen's own present (positive) parts, so it is a
# within-sample transform: z(c*v) = z(v) for any c > 0 (exact per-sample
# scale-invariance), and absent/zero parts map to 0. This deliberately does NOT
# reuse any package rclr helper that centres using cross-sample / reference
# statistics. fit and score share this ONE helper.
.icp_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Per-sample rCLR over a feature matrix (rows = samples), preserving dimnames.
.icp_rclr_matrix <- function(X) {
  Z <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    .icp_rclr(X[i, ])
  }))
  dimnames(Z) <- dimnames(X)
  Z
}

# Single-specimen score: frozen logistic linear predictor over the rCLR of the
# selected features present in this specimen's panel. Missing selected features
# contribute 0 (absent rCLR coordinate). Below min_selected present selected
# features -> neutral 0. The rCLR-origin trap (a specimen with NO positive part
# over the present universe) maps to an all-zero rCLR, so the score is the
# intercept-plus-zero linear predictor, i.e. the flat composition is SCORED (not
# floored) when at least min_selected selected features are present.
.icp_score_row <- function(x, universe, selected, coefficients,
                           intercept, min_selected) {
  univ_present <- intersect(universe, names(x))
  if (length(univ_present) == 0L) return(0)
  v <- x[univ_present]
  z <- .icp_rclr(v)
  names(z) <- univ_present

  use <- intersect(selected, univ_present)
  if (length(use) < min_selected) return(0)
  intercept + sum(coefficients[use] * z[use])
}
