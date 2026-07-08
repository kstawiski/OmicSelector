#' @title Fit StableMate stable-predictor discriminator (transfer at training,
#'   single-sample at inference)
#'
#' @description
#' Learns a frozen logistic discriminator over the subset of robust-CLR (rCLR)
#' features that StableMate identifies as both PREDICTIVE and STABLE across the
#' training environments (cohorts) named by \code{meta_train[[cohort_col]]}. The
#' StableMate idea is that predictors whose label-association is consistent
#' across cohorts are the batch-robust ones, while predictors that are
#' informative in a single cohort only are spurious; the stable set is therefore
#' the transfer-robust feature panel. This makes \code{fit_sel_stablemate} a
#' transfer-estimand method at TRAINING (selection uses the cross-environment
#' structure), but its deployment is fully SINGLE-SAMPLE: the frozen model is a
#' linear logistic score over the selected features evaluated on one specimen's
#' own rCLR vector. The transfer aspect is entirely in the SELECTION.
#'
#' Pipeline: (1) transform \code{X_train} to per-sample rCLR over a frozen
#' candidate universe (the \code{max_features} highest-variance training
#' features, frozen, to keep StableMate's bootstrap of \code{K} selections
#' tractable); (2) call \code{StableMate::stablemate} with the cohort factor as
#' the environment to obtain the predictivity and stability ensembles, then
#' extract the stable-and-predictive set with the package's own documented
#' selection rule (\code{print.stablemate}: predictors significant in the
#' predictivity ensemble AND in the stability ensemble at \code{sigthresh});
#' (3) fit a frozen \code{glm(y ~ ., binomial)} over the rCLR of the SELECTED
#' stable features and store the intercept + coefficients. The entire StableMate
#' call is wrapped in a global-\code{.Random.seed} save/restore guard with
#' \code{set.seed(hp$seed)} inside it (StableMate's stochastic stepwise
#' selection \code{st2e} uses the RNG), so fitting is deterministic and leaves
#' the caller's RNG state untouched. StableMate is run on the in-process
#' SEQUENTIAL path (\code{ncore = NULL}), which executes under our seed in the
#' main process with no parallel RNG streams; the \code{ncore = 1L} path is a
#' \pkg{foreach}/SNOW cluster path with independent RNG streams and a namespace
#' export bug, so it is deliberately NOT used.
#'
#' DEGENERATION: if \code{meta_train} is NULL, the cohort column is missing, or
#' there is a single usable environment, StableMate cannot assess
#' cross-environment stability. Fitting then degrades GRACEFULLY to a
#' predictivity-only selection: a single pooled environment is passed to
#' StableMate and the predictivity-significant set is used (stability is
#' undefined with one environment). \code{n_environments == 1L} records this so a
#' caller that intends cross-cohort transfer can audit whether a mis-wired
#' \code{meta_train} silently demoted this transfer method to pooled
#' predictivity-only selection. \code{n_environments} is the auditable
#' meta-engagement signal: \code{> 1L} means env-aware stable selection ran.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. When it contains
#'   \code{cohort_col}, non-missing values define the training environments used
#'   for env-aware stable selection.
#' @param hp List of frozen hyperparameters: \code{cohort_col} (cohort/environment
#'   column in \code{meta_train}, default \code{"accession"}); \code{K}
#'   (StableMate bootstrap selections, default \code{100L}, integer \eqn{\ge 2});
#'   \code{max_features} (candidate cap by training variance, default \code{50L},
#'   integer \eqn{\ge 2}); \code{min_selected} (minimum present selected features
#'   in a specimen for a non-neutral score, default \code{1L}, integer
#'   \eqn{\ge 1}); \code{sigthresh} (StableMate significance threshold for the
#'   predictivity/stability selection, default \code{0.9}, in \code{(0, 1)}); and
#'   \code{seed} (default \code{1L}, used to make the StableMate selection
#'   reproducible while restoring the global RNG state). Resolved once at fit and
#'   not tuned inside \code{fit_sel_stablemate()}.
#'
#' @return A plain list of class \code{sel_stablemate_model} containing the frozen
#'   \code{selected_features} (the stable-and-predictive rCLR feature ids), the
#'   logistic \code{intercept} and named \code{coefficients}, the candidate
#'   \code{feature_universe} (the frozen rCLR universe), \code{n_environments}
#'   (the auditable meta-engagement signal: \code{> 1L} = env-aware stable
#'   selection, \code{== 1L} = degenerated predictivity-only), \code{cohort_col},
#'   \code{min_selected}, a \code{degenerate} flag (\code{TRUE} when no feature
#'   was selected, in which case scoring returns the neutral \code{0}), and the
#'   resolved \code{hp}.
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
#' model <- fit_sel_stablemate(X, y, meta_train = meta)
#' score <- score_sel_stablemate(model, X)
#' }
#'
#' @references
#' Deng Y, Liao Y, Wong NKP, Bhuva DD, Lin Y, Cao Y. (2024) StableMate: a
#' statistical method to select stable predictors in omics data. \emph{bioRxiv}.
#'
#' Aitchison J. (1986) \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' @export
fit_sel_stablemate <- function(X_train, y_train, meta_train = NULL,
                               hp = list()) {
  if (!requireNamespace("StableMate", quietly = TRUE)) {
    stop("fit_sel_stablemate: package 'StableMate' is required for sel-stablemate",
         call. = FALSE)
  }
  X_train <- .reo_check_matrix(X_train, "fit_sel_stablemate", "X_train")
  if (ncol(X_train) < 2L) {
    stop("fit_sel_stablemate: X_train must contain at least two features")
  }
  .reo_check_meta(meta_train, nrow(X_train), "fit_sel_stablemate", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_sel_stablemate")
  hp <- .sel_stablemate_resolve_hp(hp, ncol(X_train))

  # Frozen candidate universe: top-max_features by training variance, name-tied
  # deterministically. The rCLR feature matrix is per-sample (within-sample),
  # so this universe is the only training-derived structure besides the model.
  universe <- .sel_stablemate_universe(X_train, hp$max_features)
  Z_train <- .sel_stablemate_rclr_matrix(X_train[, universe, drop = FALSE])

  env_info <- .sel_stablemate_environment(meta_train, y, hp$cohort_col)

  selected <- .sel_stablemate_select(Z_train, y, env_info,
                                     hp$K, hp$sigthresh, hp$seed)

  degenerate <- length(selected) == 0L
  intercept <- 0
  coefficients <- numeric(0)
  if (!degenerate) {
    cal <- .sel_stablemate_calibrate(Z_train, y, selected)
    intercept <- cal$intercept
    coefficients <- cal$coefficients
    selected <- names(coefficients)
    degenerate <- length(selected) == 0L
  }

  model <- list(
    selected_features = unname(selected),
    intercept = intercept,
    coefficients = coefficients,
    feature_universe = universe,
    n_environments = env_info$n_environments,
    cohort_col = hp$cohort_col,
    min_selected = hp$min_selected,
    degenerate = degenerate,
    hp = hp
  )
  class(model) <- "sel_stablemate_model"
  model
}


#' @title Score StableMate stable-predictor discriminator (single-sample)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen stable-feature
#' logistic learned by \code{\link{fit_sel_stablemate}}. For one specimen its
#' per-sample rCLR is computed over the frozen candidate universe (present
#' features only), restricted to the frozen selected stable features, and the
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
#' @param model A \code{sel_stablemate_model} object returned by
#'   \code{\link{fit_sel_stablemate}}.
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
#' model <- fit_sel_stablemate(X, y, meta_train = meta)
#' score_sel_stablemate(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Deng Y, Liao Y, Wong NKP, Bhuva DD, Lin Y, Cao Y. (2024) StableMate: a
#' statistical method to select stable predictors in omics data. \emph{bioRxiv}.
#'
#' @export
score_sel_stablemate <- function(model, X, meta = NULL) {
  if (!inherits(model, "sel_stablemate_model")) {
    stop("score_sel_stablemate: model must have class sel_stablemate_model")
  }
  X <- .reo_check_matrix(X, "score_sel_stablemate", "X")
  .reo_check_meta(meta, nrow(X), "score_sel_stablemate", "meta")

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
    .sel_stablemate_score_row(
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
    stop("score_sel_stablemate: scorer produced non-finite or wrong-length output")
  }
  out
}

.sel_stablemate_resolve_hp <- function(hp, p) {
  if (!is.list(hp)) stop("fit_sel_stablemate: hp must be a list")
  allowed <- c("cohort_col", "K", "max_features", "min_selected", "sigthresh",
               "seed")
  nm <- names(hp)
  if (length(hp) > 0L && (is.null(nm) || any(!nzchar(nm)))) {
    stop("fit_sel_stablemate: all hp fields must be named")
  }
  if (any(duplicated(nm))) {
    stop("fit_sel_stablemate: duplicate hp field(s): ",
         paste(unique(nm[duplicated(nm)]), collapse = ", "))
  }
  unknown <- setdiff(nm, allowed)
  if (length(unknown) > 0L) {
    stop("fit_sel_stablemate: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  cohort_col <- hp[["cohort_col"]]
  if (is.null(cohort_col)) cohort_col <- "accession"
  if (!is.character(cohort_col) || length(cohort_col) != 1L ||
      is.na(cohort_col) || !nzchar(cohort_col)) {
    stop("fit_sel_stablemate: hp$cohort_col must be a single non-empty character string")
  }

  K <- hp[["K"]]
  if (is.null(K)) K <- 100L
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) ||
      K < 2L || K > .Machine$integer.max || K != as.integer(K)) {
    stop("fit_sel_stablemate: hp$K must be an integer >= 2")
  }
  K <- as.integer(K)

  max_features <- hp[["max_features"]]
  if (is.null(max_features)) max_features <- 50L
  if (!is.numeric(max_features) || length(max_features) != 1L ||
      !is.finite(max_features) || max_features < 2L ||
      max_features > .Machine$integer.max ||
      max_features != as.integer(max_features)) {
    stop("fit_sel_stablemate: hp$max_features must be an integer >= 2")
  }
  max_features <- min(as.integer(max_features), p)

  min_selected <- hp[["min_selected"]]
  if (is.null(min_selected)) min_selected <- 1L
  if (!is.numeric(min_selected) || length(min_selected) != 1L ||
      !is.finite(min_selected) || min_selected < 1L ||
      min_selected > .Machine$integer.max ||
      min_selected != as.integer(min_selected)) {
    stop("fit_sel_stablemate: hp$min_selected must be a positive integer")
  }

  sigthresh <- hp[["sigthresh"]]
  if (is.null(sigthresh)) sigthresh <- 0.9
  if (!is.numeric(sigthresh) || length(sigthresh) != 1L ||
      !is.finite(sigthresh) || sigthresh <= 0 || sigthresh >= 1) {
    stop("fit_sel_stablemate: hp$sigthresh must be a number in (0, 1)")
  }

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_sel_stablemate: hp$seed must be a non-negative integer")
  }

  list(
    cohort_col = cohort_col,
    K = K,
    max_features = max_features,
    min_selected = as.integer(min_selected),
    sigthresh = as.numeric(sigthresh),
    seed = as.integer(seed)
  )
}

# Frozen candidate universe: the top-max_features features by training variance,
# with feature-name radix tie-break so the universe is reproducible across
# machines/locales.
.sel_stablemate_universe <- function(X_train, max_features) {
  v <- apply(X_train, 2L, stats::var)
  v[!is.finite(v)] <- -Inf
  ord <- order(-v, colnames(X_train), method = "radix")
  colnames(X_train)[ord][seq_len(max_features)]
}

# Resolve the training environment factor from meta and record the auditable
# n_environments engagement signal. A usable environment must contain at least
# one case AND one control (StableMate needs both classes per environment to
# assess stability). When fewer than two usable environments exist, all rows
# collapse to a single pooled environment (degeneration to predictivity-only
# selection) and n_environments == 1L.
.sel_stablemate_environment <- function(meta_train, y, cohort_col) {
  pooled <- list(env = factor(rep("all", length(y))), n_environments = 1L)
  if (is.null(meta_train)) return(pooled)
  meta_train <- as.data.frame(meta_train, stringsAsFactors = FALSE)
  if (!cohort_col %in% names(meta_train)) return(pooled)

  cohort <- as.character(meta_train[[cohort_col]])
  non_missing <- !is.na(cohort)
  levels_all <- unique(cohort[non_missing])
  usable <- levels_all[vapply(levels_all, function(lv) {
    rows <- non_missing & cohort == lv
    any(y[rows] == 1L) && any(y[rows] == 0L)
  }, logical(1))]
  if (length(usable) <= 1L) return(pooled)

  # Keep only rows in usable environments; rows in unusable/missing environments
  # are folded into the nearest-usable structure by relabelling them "all" would
  # change n, so instead we drop them from the env definition by assigning them
  # to their own (non-usable) levels -- but StableMate needs every row to have an
  # environment. Rows outside usable environments are assigned to a single shared
  # pooled level so they still contribute to predictivity without defining a
  # spurious singleton environment.
  env_chr <- ifelse(non_missing & cohort %in% usable, cohort, "__pooled__")
  list(env = factor(env_chr), n_environments = length(usable))
}

# Run StableMate predictor selection under a global-.Random.seed save/restore
# guard so the selection is deterministic (set.seed(seed) inside) and the
# caller's RNG state is left untouched. StableMate is run on the in-process
# SEQUENTIAL path (ncore = NULL): it executes under our seed in the main process
# with no parallel RNG streams. (ncore = 1L is a foreach/SNOW cluster path with
# independent RNG streams and a namespace-export bug, so it is NOT used.)
#
# With >= 2 usable environments the full env-aware StableMate::stablemate runs
# and the stable-and-predictive set is extracted with the package's own
# documented selection rule (print.stablemate). With a single environment
# (degeneration) StableMate's stability objective is undefined -- stablemate()
# errors on a single-level env -- so selection gracefully degrades to a
# predictivity-only st2e ensemble (no env, BIC-logit objective), the same
# StableMate engine without the stability machinery.
.sel_stablemate_select <- function(Z_train, y, env_info, K, sigthresh, seed) {
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

  if (env_info$n_environments > 1L) {
    fit <- tryCatch(
      suppressWarnings(suppressMessages(
        StableMate::stablemate(
          Y = y,
          X = Z_train,
          env = env_info$env,
          K = K,
          ret_mod = TRUE,
          ret_imp = TRUE,
          verbose = FALSE,
          ncore = NULL
        )
      )),
      error = function(e) {
        stop("fit_sel_stablemate: StableMate selection failed: ",
             conditionMessage(e), call. = FALSE)
      }
    )
    return(.sel_stablemate_extract_stable(fit, sigthresh))
  }
  .sel_stablemate_predictivity_only(Z_train, y, K, sigthresh)
}

# Extract the stable-and-predictive predictor names using StableMate's own
# documented selection logic (print.stablemate), capturing its invisible return
# without printing to the console. A predictor is selected when it is significant
# (selected more often than the pseudo-predictor at sigthresh) in BOTH the
# predictivity ensemble (joint importance) AND the stability ensemble
# (conditional importance), intersected -- exactly the package's definition. When
# the stability ensemble yields NO significant predictor at sigthresh, the
# intersection is empty and the model is left degenerate (neutral-0): a
# TRANSFER/stability method must NOT fall back to predictivity-only features,
# because those can be environment-specific (batch-driven) signal -- selecting
# them while reporting n_environments > 1 would freeze exactly the spurious
# cross-cohort instability the method exists to reject. "No stable predictor
# found" is an honest, correct outcome, not a failure to be patched over.
.sel_stablemate_extract_stable <- function(fit, sigthresh) {
  pred_imp <- fit$prediction_ensemble$imp_scores$joint
  stab_imp <- fit$stable_ensemble$imp_scores$conditional
  if (is.null(pred_imp) || is.null(stab_imp)) {
    return(character())
  }
  pred_sig <- pred_imp$significance[-1]
  stab_sig <- stab_imp$significance[-1]
  pred <- names(pred_sig)[pred_sig > sigthresh]
  stab <- names(stab_sig)[stab_sig > sigthresh]
  intersect(pred, stab)
}

# Predictivity-only degeneration: a single st2e ensemble with no environment and
# the BIC-logit objective (the same StableMate stochastic-stepwise engine without
# the cross-environment stability objective, which is undefined for one
# environment). ret_imp = TRUE computes the joint importance scores inside st2e;
# a predictor is selected when it is significant (selected more often than the
# pseudo-predictor) at sigthresh.
.sel_stablemate_predictivity_only <- function(Z_train, y, K, sigthresh) {
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      StableMate::st2e(
        Y = y,
        X = Z_train,
        env = NULL,
        obj_fun = StableMate::obj_bic_logit,
        reg_fun = StableMate::reg_logit,
        K = K,
        ret_mod = TRUE,
        ret_imp = TRUE,
        verbose = FALSE,
        ncore = NULL
      )
    )),
    error = function(e) {
      stop("fit_sel_stablemate: StableMate predictivity-only selection failed: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  sig <- fit$imp_scores$joint$significance
  if (is.null(sig)) return(character(0))
  sig <- sig[setdiff(names(sig), "pseudo_predictor")]
  as.character(names(sig)[sig > sigthresh])
}

# Frozen logistic calibration over the rCLR of the SELECTED stable features.
# glm(y ~ ., binomial); store intercept + coefficients keyed by feature name.
# Non-finite/dropped coefficients (rank deficiency, perfect separation) are
# removed so only finite linear terms are frozen. If glm fails or no finite
# coefficient survives, the model degenerates (empty selection -> neutral 0).
.sel_stablemate_calibrate <- function(Z_train, y, selected) {
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
.sel_stablemate_rclr <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Per-sample rCLR over a feature matrix (rows = samples), preserving dimnames.
.sel_stablemate_rclr_matrix <- function(X) {
  Z <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    .sel_stablemate_rclr(X[i, ])
  }))
  dimnames(Z) <- dimnames(X)
  Z
}

# Single-specimen score: frozen logistic linear predictor over the rCLR of the
# selected features present in this specimen's panel. Missing selected features
# contribute 0 (absent rCLR coordinate). Below min_selected present selected
# features -> neutral 0.
.sel_stablemate_score_row <- function(x, universe, selected, coefficients,
                                      intercept, min_selected) {
  univ_present <- intersect(universe, names(x))
  if (length(univ_present) == 0L) return(0)
  v <- x[univ_present]
  z <- .sel_stablemate_rclr(v)
  names(z) <- univ_present

  use <- intersect(selected, univ_present)
  if (length(use) < min_selected) return(0)
  intercept + sum(coefficients[use] * z[use])
}
