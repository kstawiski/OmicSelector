#' @title Group-DRO kit x biofluid robust discriminator (canonical single-sample wrapper)
#'
#' @description
#' Canonical \code{fit_*}/\code{score_*} wrapper around the existing Group-DRO
#' kit x biofluid scorer primitive (\code{\link{fit_group_dro_scorer}} /
#' \code{\link{score_group_dro_scorer}}, family N, transfer estimand). The
#' primitive fits a train-only robust-CLR logistic classifier whose
#' exponentiated-gradient sample weights are pooled over kit x biofluid groups
#' (Sagawa et al., ICLR 2020), so the trained scorer is robust to the
#' worst-performing provenance group rather than only the pooled average. The
#' SCORE path of the primitive is already single-sample / row-equivariant: a
#' specimen is mapped to the frozen train-only CLR with a per-row pseudocount and
#' per-row pivot centring, scaled by the frozen median/MAD constants, and pushed
#' through the frozen logistic head; it carries no cross-row statistic and depends
#' only on the specimen and the frozen model. This wrapper does NOT reimplement any
#' Group-DRO math; it only adapts the primitive to the canonical contract
#' (\code{fit_dro_group(X_train, y_train, meta_train, hp)} /
#' \code{score_dro_group(model, X, meta)}) so the method is dispatchable via
#' \code{\link{singlesample_score_call}}.
#'
#' \strong{Auditable group engagement.} The Group-DRO objective only differs from
#' plain ERM when there is more than one kit x biofluid group. When the canonical
#' \code{meta_train} is \code{NULL} or lacks the configured kit / biofluid columns,
#' this wrapper DEGENERATES GRACEFULLY to a single-group fit by synthesizing a
#' constant single-group metadata frame, so the primitive runs without crashing and
#' still yields a valid single-sample model -- but the run is then mathematically
#' ERM, not group-robust. To never silently mislabel the estimand, the model records
#' \code{n_groups} (the number of distinct kit x biofluid groups actually engaged)
#' and the logical \code{groups_engaged} (\code{n_groups > 1}); these mirror the
#' \code{n_cohorts} / \code{n_environments} engagement signals exposed by
#' \code{reo-metaktsp} and \code{sel-stablemate}.
#'
#' \strong{Orientation.} \code{score_group_dro_scorer} returns the logistic
#' probability of the case class, so larger already means more case-like; the
#' wrapper preserves this orientation unchanged.
#'
#' \strong{Metadata at score time.} \code{meta} passed to \code{score_dro_group} is
#' forwarded to the primitive only as a diagnostic (it attaches test-group audit
#' attributes); it never enters the numeric score, which is stripped to a plain
#' finite numeric vector before return.
#'
#' @references
#' Sagawa S, Koh PW, Hashimoto TB, Liang P. (2020) Distributionally Robust Neural
#' Networks for Group Shifts. \emph{ICLR}. arXiv:1911.08731.
#'
#' @name singlesample-dro-group
NULL


# ----------------------------------------------------------------------------
# Hyperparameter resolver (strict allow-list)
# ----------------------------------------------------------------------------

# Surfaces the primitive's tunables through the canonical `hp` list. Defaults
# match fit_group_dro_scorer()'s formal defaults exactly so fit_dro_group() with
# hp = list() reproduces the primitive's default behaviour. Unknown / duplicate /
# unnamed hp names are rejected before any are consumed (mirrors the ecod-copod
# resolver convention).
.dro_group_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_dro_group: hp must be a list")
  if (length(hp) > 0L) {
    nm <- names(hp)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("fit_dro_group: all hp entries must be named")
    }
    if (any(duplicated(nm))) {
      stop("fit_dro_group: duplicate hp field(s): ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "))
    }
  }
  allowed <- c("kit_col", "biofluid_col", "panel_size", "panel_features",
               "eta_grid", "l2_grid", "inner_folds", "learning_rate", "epochs",
               "tune", "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_dro_group: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  one_string <- function(v, field) {
    if (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v)) {
      stop("fit_dro_group: hp$", field, " must be a single non-empty string")
    }
    v
  }
  pos_int <- function(v, field) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1L ||
        v > .Machine$integer.max || v != as.integer(v)) {
      stop("fit_dro_group: hp$", field, " must be a positive integer")
    }
    as.integer(v)
  }
  pos_num <- function(v, field) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop("fit_dro_group: hp$", field, " must be a positive finite number")
    }
    as.numeric(v)
  }
  pos_num_grid <- function(v, field) {
    if (!is.numeric(v) || length(v) < 1L || any(!is.finite(v)) || any(v <= 0)) {
      stop("fit_dro_group: hp$", field,
           " must be a non-empty numeric vector of positive finite values")
    }
    as.numeric(v)
  }

  kit_col <- if (is.null(hp$kit_col)) "kit_family" else one_string(hp$kit_col, "kit_col")
  biofluid_col <- if (is.null(hp$biofluid_col)) "biofluid" else one_string(hp$biofluid_col, "biofluid_col")
  panel_size <- if (is.null(hp$panel_size)) 20L else pos_int(hp$panel_size, "panel_size")

  panel_features <- hp$panel_features
  if (!is.null(panel_features)) {
    if (!is.character(panel_features) || length(panel_features) < 1L ||
        any(is.na(panel_features)) || any(!nzchar(panel_features))) {
      stop("fit_dro_group: hp$panel_features must be a non-empty character vector")
    }
    panel_features <- as.character(panel_features)
  }

  eta_grid <- if (is.null(hp$eta_grid)) c(0.01, 0.05, 0.10) else pos_num_grid(hp$eta_grid, "eta_grid")
  l2_grid <- if (is.null(hp$l2_grid)) c(0.001, 0.01, 0.10) else pos_num_grid(hp$l2_grid, "l2_grid")
  inner_folds <- if (is.null(hp$inner_folds)) 3L else pos_int(hp$inner_folds, "inner_folds")
  learning_rate <- if (is.null(hp$learning_rate)) 0.05 else pos_num(hp$learning_rate, "learning_rate")
  epochs <- if (is.null(hp$epochs)) 200L else pos_int(hp$epochs, "epochs")

  tune <- if (is.null(hp$tune)) TRUE else hp$tune
  if (!is.logical(tune) || length(tune) != 1L || is.na(tune)) {
    stop("fit_dro_group: hp$tune must be a single logical")
  }

  seed <- if (is.null(hp$seed)) 42L else pos_int(hp$seed, "seed")

  list(
    kit_col = kit_col,
    biofluid_col = biofluid_col,
    panel_size = panel_size,
    panel_features = panel_features,
    eta_grid = eta_grid,
    l2_grid = l2_grid,
    inner_folds = inner_folds,
    learning_rate = learning_rate,
    epochs = epochs,
    tune = tune,
    seed = seed
  )
}

# Resolve the metadata actually handed to the primitive. The primitive REQUIRES
# a sample_meta carrying the kit / biofluid columns. When the canonical
# meta_train is NULL or lacks EITHER configured column (the primitive needs
# both present), synthesize a constant single-group frame so the primitive runs
# as ERM (one kit x biofluid group). Returns the
# frame plus a logical flag stating whether real (caller-supplied) groups were
# engaged. A meta that supplies the columns but happens to be constant is treated
# as caller-supplied (groups_engaged depends only on the realised n_groups, set
# downstream from the fit).
.dro_group_resolve_meta <- function(meta_train, n, kit_col, biofluid_col) {
  if (!is.null(meta_train)) {
    meta_train <- as.data.frame(meta_train, stringsAsFactors = FALSE)
    if (nrow(meta_train) != n) {
      stop("fit_dro_group: meta_train must have one row per row of X_train")
    }
    if (all(c(kit_col, biofluid_col) %in% names(meta_train))) {
      return(list(meta = meta_train, real_meta = TRUE))
    }
  }
  meta <- data.frame(
    rep("all", n), rep("all", n),
    stringsAsFactors = FALSE
  )
  names(meta) <- c(kit_col, biofluid_col)
  list(meta = meta, real_meta = FALSE)
}


# ----------------------------------------------------------------------------
# fit_dro_group
# ----------------------------------------------------------------------------

#' @title Fit the Group-DRO kit x biofluid robust discriminator
#'
#' @description
#' Canonical-contract fit wrapper for the Group-DRO kit x biofluid scorer
#' primitive. The training matrix, labels, and (optional) kit x biofluid metadata
#' are forwarded to \code{\link{fit_group_dro_scorer}}, which performs all
#' train-only CLR panel selection, optional nested-CV hyperparameter tuning, and
#' exponentiated-gradient Group-DRO logistic optimisation. When \code{meta_train}
#' is \code{NULL} or lacks the configured kit / biofluid columns the fit degenerates
#' gracefully to a single constant group (ERM); the realised number of groups is
#' recorded as an auditable engagement signal.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param y_train Numeric / integer 0/1 labels (1 = case), length
#'   \code{nrow(X_train)}, with both classes present.
#' @param meta_train Optional per-sample metadata carrying the kit and biofluid
#'   columns named by \code{hp$kit_col} / \code{hp$biofluid_col}. When \code{NULL}
#'   or missing those columns, the fit runs as a single-group (ERM) Group-DRO model.
#' @param hp Optional list of hyperparameters (strict allow-list). Fields, with
#'   primitive-matching defaults: \code{kit_col} (\code{"kit_family"}),
#'   \code{biofluid_col} (\code{"biofluid"}), \code{panel_size} (\code{20L}),
#'   \code{panel_features} (\code{NULL}; train-only CLR panel selected when unset),
#'   \code{eta_grid} (\code{c(0.01, 0.05, 0.10)}), \code{l2_grid}
#'   (\code{c(0.001, 0.01, 0.10)}), \code{inner_folds} (\code{3L}),
#'   \code{learning_rate} (\code{0.05}), \code{epochs} (\code{200L}), \code{tune}
#'   (\code{TRUE}), \code{seed} (\code{42L}).
#'
#' @return Object of class \code{dro_group_model}: a list with the frozen
#'   primitive \code{fit} (a \code{group_dro_scorer} object), \code{n_groups} (the
#'   number of distinct kit x biofluid groups engaged), \code{groups_engaged}
#'   (\code{n_groups > 1}), \code{train_groups}, the resolved \code{kit_col} /
#'   \code{biofluid_col}, and \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 80; p <- 25; k <- 6
#' L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
#'             dimnames = list(NULL, paste0("hsa-miR-", seq_len(p))))
#' y <- rep(c(0, 1), each = n / 2)
#' L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
#' X <- exp(L)
#' meta <- data.frame(kit_family = rep(c("A", "B"), length.out = n),
#'                    biofluid = rep(c("plasma", "serum"), each = n / 2))
#' model <- fit_dro_group(X, y, meta)
#' score_dro_group(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Sagawa S, Koh PW, Hashimoto TB, Liang P. (2020) Distributionally Robust Neural
#' Networks for Group Shifts. \emph{ICLR}. arXiv:1911.08731.
#'
#' @seealso \code{\link{score_dro_group}}, \code{\link{fit_group_dro_scorer}}
#' @export
fit_dro_group <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_dro_group", "X_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_dro_group")
  hp <- .dro_group_resolve_hp(hp)

  resolved <- .dro_group_resolve_meta(meta_train, nrow(X_train),
                                      hp$kit_col, hp$biofluid_col)

  # fit_group_dro_scorer() calls set.seed() internally (deterministic given the
  # frozen hp$seed) but thereby mutates the global RNG; save/restore so the fit
  # is RNG-safe for the caller.
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

  fit <- fit_group_dro_scorer(
    X_train, y, resolved$meta,
    kit_col = hp$kit_col,
    biofluid_col = hp$biofluid_col,
    panel_features = hp$panel_features,
    panel_size = hp$panel_size,
    eta_grid = hp$eta_grid,
    l2_grid = hp$l2_grid,
    inner_folds = hp$inner_folds,
    learning_rate = hp$learning_rate,
    epochs = hp$epochs,
    tune = hp$tune,
    seed = hp$seed)

  n_groups <- as.integer(fit$n_groups)
  model <- list(
    fit = fit,
    n_groups = n_groups,
    groups_engaged = isTRUE(n_groups > 1L),
    real_meta = resolved$real_meta,
    train_groups = fit$train_groups,
    kit_col = hp$kit_col,
    biofluid_col = hp$biofluid_col,
    hp = hp
  )
  class(model) <- "dro_group_model"
  model
}


# ----------------------------------------------------------------------------
# score_dro_group
# ----------------------------------------------------------------------------

#' @title Score the Group-DRO kit x biofluid robust discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model from
#' \code{\link{fit_dro_group}} by delegating to \code{\link{score_group_dro_scorer}}
#' on the frozen primitive fit. The primitive maps each specimen to the frozen
#' train-only CLR (per-row pseudocount and pivot centring), scales it with the
#' frozen median / MAD constants, and applies the frozen logistic head; the
#' returned logistic probability is already oriented so larger = more case-like.
#' The score of a row depends only on that row and the frozen model.
#'
#' @param model A \code{dro_group_model} object from \code{\link{fit_dro_group}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundances with named feature columns.
#' @param meta Optional per-sample metadata. Forwarded to the primitive only as a
#'   diagnostic (it attaches test-group audit attributes) and never affects the
#'   numeric score, which is returned as a plain finite numeric vector.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}; larger values are
#'   more case-like.
#'
#' @examples
#' \dontrun{
#' model <- fit_dro_group(X, y, meta)
#' score_dro_group(model, X[1, , drop = FALSE])
#' }
#'
#' @seealso \code{\link{fit_dro_group}}, \code{\link{score_group_dro_scorer}}
#' @export
score_dro_group <- function(model, X, meta = NULL) {
  if (!inherits(model, "dro_group_model")) {
    stop("score_dro_group: model must have class dro_group_model")
  }
  X <- .reo_check_matrix(X, "score_dro_group", "X")
  if (!is.null(meta)) {
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (nrow(meta) != nrow(X)) {
      stop("score_dro_group: meta must have one row per row of X")
    }
  }

  out <- score_group_dro_scorer(model$fit, X, sample_meta = meta)
  out <- as.numeric(out)  # drop any diagnostic attributes; keep only the scores
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_dro_group: scorer produced non-finite or wrong-length output")
  }
  out
}


# ----------------------------------------------------------------------------
# Stable model digest for the equivariance harness
# ----------------------------------------------------------------------------

# Hash ONLY the frozen R-side state needed to reproduce the score: the primitive
# fit's panel features, logistic head (beta, intercept), frozen CLR fit and
# scaler, and the group / orientation metadata, plus the engagement signal and
# hp. The model holds no external pointer, so this is a deterministic snapshot
# suitable as the `model_digest` for singlesample_assert_row_equivariant()
# (clause (d): the bytes must be identical before and after scoring).
.dro_group_model_digest <- function(model) {
  fit <- model$fit
  state <- list(
    panel_features = fit$panel_features,
    beta = fit$beta,
    intercept = fit$intercept,
    clr_fit = fit$clr_fit,
    scaler = fit$scaler,
    feature_directions = fit$feature_directions,
    train_groups = fit$train_groups,
    kit_col = fit$kit_col,
    biofluid_col = fit$biofluid_col,
    n_groups = model$n_groups,
    groups_engaged = model$groups_engaged,
    hp = model$hp
  )
  digest::digest(state, algo = "xxhash64")
}
