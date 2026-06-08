#' @title Fit selbal-selected single-balance discriminator
#'
#' @description
#' Learns the two feature groups (numerator vs denominator) of a SINGLE
#' isometric-log-ratio (ILR) balance from the training data with the
#' \pkg{selbal} forward-selection algorithm, freezes those groups, and stores a
#' frozen univariate logistic calibration of the self-contained per-sample
#' balance. The selbal selection is the only part that touches the training
#' labels and is therefore leakage-safe: it runs once at fit time, and the
#' frozen numerator/denominator groups are then evaluated on each specimen
#' independently. \pkg{selbal} differs from \code{\link{fit_ws_balance_ilr}} in
#' exactly one way: instead of a prespecified sequential binary partition, the
#' two sides of the balance are LEARNED from the training data.
#'
#' At inference (see \code{\link{score_bal_selbal}}) the balance over the frozen
#' groups is a self-contained per-sample log-contrast:
#' \deqn{B = \sqrt{r s / (r + s)} \,\bigl(\overline{\log v_N} -
#'   \overline{\log v_D}\bigr),}
#' where the means are over the features of each frozen group that are both
#' PRESENT in the scored panel and strictly POSITIVE in that specimen, and
#' \eqn{r}, \eqn{s} are the per-specimen counts of those present-positive
#' numerator and denominator features. Because the balance is a difference of
#' means of logs, it is EXACTLY invariant to per-sample positive scaling (the
#' per-sample scale cancels), and it uses only that specimen's own values, so
#' the method is single-sample deployable and row-equivariant. The calibration
#' recomputes the training balance through the SAME self-contained helper used
#' at score time (not selbal's internal \code{$balance.values}, which is
#' computed over the full zero-replaced closure) so fit and score share one
#' balance definition and the stored intercept/slope match exactly what
#' \code{\link{score_bal_selbal}} computes.
#'
#' The entire \pkg{selbal} call (whose internal accuracy evaluation may touch
#' the RNG) is wrapped in a global-\code{.Random.seed} save/restore guard with
#' \code{set.seed(hp$seed)} inside it, so fitting is deterministic and leaves
#' the caller's RNG state untouched.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{max_vars} (selbal
#'   \code{maxV}, the cap on the total number of features in the balance,
#'   default \code{10L}, integer \eqn{\ge 2}); \code{logit_acc} (selbal accuracy
#'   measure driving forward selection, default \code{"Dev"}, one of
#'   \code{c("AUC", "Dev", "Rsq", "Tjur")}; deviance is a continuous objective
#'   that keeps adding informative features after class separation, recovering a
#'   more robust multi-feature balance than rank-saturating \code{"AUC"}, which
#'   tends to stop early at a tiny balance); \code{zero_rep} (selbal zero
#'   replacement, default \code{"bayes"}, one of \code{c("bayes", "one")});
#'   \code{min_group_coverage} (minimum present-positive features required in
#'   EACH frozen group for a non-neutral score, default \code{1L}, integer
#'   \eqn{\ge 1}); and \code{seed} (default \code{1L}, used only to make the
#'   selbal selection reproducible while restoring the global RNG state).
#'
#' @return A plain list of class \code{bal_selbal_model} containing the frozen
#'   \code{numerator} and \code{denominator} feature-name vectors, the logistic
#'   \code{intercept} and \code{slope}, an \code{orientation_sign} fallback
#'   (used only when the logistic slope is degenerate/non-finite),
#'   \code{min_group_coverage}, \code{feature_universe}, a \code{degenerate}
#'   flag, and the resolved \code{hp}. If selbal returns an empty/degenerate
#'   balance (a group with no features, or overlapping groups), \code{degenerate}
#'   is \code{TRUE} and \code{\link{score_bal_selbal}} returns the neutral score
#'   \code{0} for every specimen.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(120 * 15, shape = 30, rate = 2), nrow = 120,
#'             dimnames = list(NULL, paste0("f", sprintf("%02d", seq_len(15)))))
#' y <- rep(c(0, 1), each = 60)
#' X[y == 1, c("f01", "f02", "f03")] <- X[y == 1, c("f01", "f02", "f03")] * 1.8
#' X[y == 1, c("f13", "f14", "f15")] <- X[y == 1, c("f13", "f14", "f15")] / 1.8
#' model <- fit_bal_selbal(X, y, hp = list(logit_acc = "Dev"))
#' score <- score_bal_selbal(model, X)
#' }
#'
#' @references
#' Rivera-Pinto J, Egozcue JJ, Pawlowsky-Glahn V, Paredes R, Noguera-Julian M,
#' Calle ML. (2018) Balances: a new perspective for microbiome analysis.
#' \emph{mSystems} 3: e00053-18.
#'
#' Egozcue JJ, Pawlowsky-Glahn V. (2005) Groups of parts and their balances in
#' compositional data analysis. \emph{Mathematical Geology} 37(7): 795-828.
#'
#' @export
fit_bal_selbal <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  if (!requireNamespace("selbal", quietly = TRUE)) {
    stop("fit_bal_selbal: package 'selbal' is required for bal-selbal",
         call. = FALSE)
  }
  X_train <- .reo_check_matrix(X_train, "fit_bal_selbal", "X_train")
  if (ncol(X_train) < 2L) {
    stop("fit_bal_selbal: X_train must contain at least two features for a balance")
  }
  .reo_check_meta(meta_train, nrow(X_train), "fit_bal_selbal", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_bal_selbal")
  hp <- .bal_selbal_resolve_hp(hp, ncol(X_train))

  sel <- .bal_selbal_run_selection(X_train, y, hp)
  numerator <- as.character(sel$numerator)
  denominator <- as.character(sel$denominator)

  degenerate <- .bal_selbal_is_degenerate(numerator, denominator)
  if (degenerate) {
    numerator <- character(0)
    denominator <- character(0)
  }

  # Calibrate on OUR self-contained balance (NOT selbal's $balance.values, which
  # is over the full zero-replaced closure) so fit and score share one balance
  # definition and the stored intercept/slope match what score computes.
  intercept <- 0
  slope <- 0
  orientation_sign <- 1
  if (!degenerate) {
    B_train <- .bal_selbal_balance_vector(X_train, numerator, denominator,
                                          hp$min_group_coverage)
    cal <- .bal_selbal_calibrate(B_train, y)
    intercept <- cal$intercept
    slope <- cal$slope
    orientation_sign <- cal$orientation_sign
  }

  model <- list(
    numerator = numerator,
    denominator = denominator,
    intercept = intercept,
    slope = slope,
    orientation_sign = orientation_sign,
    min_group_coverage = hp$min_group_coverage,
    feature_universe = colnames(X_train),
    degenerate = degenerate,
    hp = hp
  )
  class(model) <- "bal_selbal_model"
  model
}


#' @title Score selbal-selected single-balance discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen balance learned by
#' \code{\link{fit_bal_selbal}}. For one specimen, the self-contained balance is
#' computed over the frozen numerator/denominator groups using only the features
#' that are PRESENT in \code{X} and strictly POSITIVE in that specimen:
#' \deqn{B = \sqrt{r s / (r + s)} \,\bigl(\overline{\log v_N} -
#'   \overline{\log v_D}\bigr),}
#' with \eqn{r}, \eqn{s} the per-specimen counts of present-positive numerator
#' and denominator features. The returned score is the frozen logistic linear
#' predictor \eqn{intercept + slope \cdot B}; larger is more case-like.
#'
#' If a specimen has fewer than \code{model$min_group_coverage} present-positive
#' features in EITHER group (so the balance is undefined/degenerate), it receives
#' the documented neutral score \code{0}. If selbal selected an empty/degenerate
#' balance at fit time (\code{model$degenerate}), every row receives \code{0}.
#' The balance is exactly invariant to per-sample positive scaling, and scoring
#' uses no random numbers and no scored-batch statistics.
#'
#' @param model A \code{bal_selbal_model} object returned by
#'   \code{\link{fit_bal_selbal}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Scores use only
#'   each row's own values plus the frozen model, with no scored-batch centering,
#'   quantiles, or renormalization.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(120 * 15, shape = 30, rate = 2), nrow = 120,
#'             dimnames = list(NULL, paste0("f", sprintf("%02d", seq_len(15)))))
#' y <- rep(c(0, 1), each = 60)
#' X[y == 1, c("f01", "f02", "f03")] <- X[y == 1, c("f01", "f02", "f03")] * 1.8
#' X[y == 1, c("f13", "f14", "f15")] <- X[y == 1, c("f13", "f14", "f15")] / 1.8
#' model <- fit_bal_selbal(X, y, hp = list(logit_acc = "Dev"))
#' score_bal_selbal(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Rivera-Pinto J, Egozcue JJ, Pawlowsky-Glahn V, Paredes R, Noguera-Julian M,
#' Calle ML. (2018) Balances: a new perspective for microbiome analysis.
#' \emph{mSystems} 3: e00053-18.
#'
#' @export
score_bal_selbal <- function(model, X, meta = NULL) {
  if (!inherits(model, "bal_selbal_model")) {
    stop("score_bal_selbal: model must have class bal_selbal_model")
  }
  X <- .reo_check_matrix(X, "score_bal_selbal", "X")
  .reo_check_meta(meta, nrow(X), "score_bal_selbal", "meta")

  if (isTRUE(model$degenerate) ||
      length(model$numerator) == 0L || length(model$denominator) == 0L) {
    return(rep(0, nrow(X)))
  }

  num_present <- intersect(model$numerator, colnames(X))
  den_present <- intersect(model$denominator, colnames(X))
  if (length(num_present) == 0L || length(den_present) == 0L) {
    return(rep(0, nrow(X)))
  }

  out <- vapply(seq_len(nrow(X)), function(i) {
    # Restore feature names: single-row slicing of a 1-column matrix drops them,
    # which would otherwise make present/positive selection depend on batch shape.
    x_i <- X[i, ]
    names(x_i) <- colnames(X)
    B <- .bal_selbal_balance_row(x_i, num_present, den_present,
                                 model$min_group_coverage)
    if (is.na(B)) return(0)
    model$intercept + model$slope * B
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_bal_selbal: scorer produced non-finite or wrong-length output")
  }
  out
}

.bal_selbal_resolve_hp <- function(hp, p) {
  if (!is.list(hp)) stop("fit_bal_selbal: hp must be a list")
  allowed <- c("max_vars", "logit_acc", "zero_rep", "min_group_coverage",
               "seed")
  nm <- names(hp)
  if (length(hp) > 0L && (is.null(nm) || any(!nzchar(nm)))) {
    stop("fit_bal_selbal: all hp fields must be named")
  }
  if (any(duplicated(nm))) {
    stop("fit_bal_selbal: duplicate hp field(s): ",
         paste(unique(nm[duplicated(nm)]), collapse = ", "))
  }
  unknown <- setdiff(nm, allowed)
  if (length(unknown) > 0L) {
    stop("fit_bal_selbal: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  max_vars <- hp[["max_vars"]]
  if (is.null(max_vars)) max_vars <- 10L
  if (!is.numeric(max_vars) || length(max_vars) != 1L ||
      !is.finite(max_vars) || max_vars < 2L ||
      max_vars > .Machine$integer.max ||
      max_vars != as.integer(max_vars)) {
    stop("fit_bal_selbal: hp$max_vars must be an integer >= 2")
  }
  max_vars <- min(as.integer(max_vars), p)

  logit_acc <- hp[["logit_acc"]]
  if (is.null(logit_acc)) logit_acc <- "Dev"
  if (!is.character(logit_acc) || length(logit_acc) != 1L ||
      is.na(logit_acc) || !logit_acc %in% c("AUC", "Dev", "Rsq", "Tjur")) {
    stop("fit_bal_selbal: hp$logit_acc must be one of 'AUC', 'Dev', 'Rsq', 'Tjur'")
  }

  zero_rep <- hp[["zero_rep"]]
  if (is.null(zero_rep)) zero_rep <- "bayes"
  if (!is.character(zero_rep) || length(zero_rep) != 1L ||
      is.na(zero_rep) || !zero_rep %in% c("bayes", "one")) {
    stop("fit_bal_selbal: hp$zero_rep must be 'bayes' or 'one'")
  }

  min_group_coverage <- hp[["min_group_coverage"]]
  if (is.null(min_group_coverage)) min_group_coverage <- 1L
  if (!is.numeric(min_group_coverage) || length(min_group_coverage) != 1L ||
      !is.finite(min_group_coverage) || min_group_coverage < 1L ||
      min_group_coverage > .Machine$integer.max ||
      min_group_coverage != as.integer(min_group_coverage)) {
    stop("fit_bal_selbal: hp$min_group_coverage must be a positive integer")
  }

  seed <- hp[["seed"]]
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_bal_selbal: hp$seed must be a non-negative integer")
  }

  list(
    max_vars = max_vars,
    logit_acc = logit_acc,
    zero_rep = zero_rep,
    min_group_coverage = as.integer(min_group_coverage),
    seed = as.integer(seed)
  )
}

# Run selbal selection under a global-.Random.seed save/restore guard so the
# selection is deterministic (set.seed(seed) inside) and the caller's RNG state
# is left untouched. selbal emits benign glm separation warnings during its
# greedy search; those are suppressed (they are not method warnings).
.bal_selbal_run_selection <- function(X_train, y, hp) {
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
  set.seed(hp$seed)

  x_df <- as.data.frame(X_train, stringsAsFactors = FALSE, check.names = FALSE)
  y_fac <- factor(y, levels = c(0L, 1L))
  sel <- tryCatch(
    suppressWarnings(suppressMessages(
      selbal::selbal(
        x = x_df,
        y = y_fac,
        maxV = hp$max_vars,
        draw = FALSE,
        logit.acc = hp$logit_acc,
        zero.rep = hp$zero_rep
      )
    )),
    error = function(e) {
      stop("fit_bal_selbal: selbal selection failed: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  sel
}

# A selection is degenerate if either group is empty or the groups overlap
# (selbal should never return overlapping groups, but a frozen balance with a
# shared feature has no well-defined log-contrast, so we guard it).
.bal_selbal_is_degenerate <- function(numerator, denominator) {
  if (length(numerator) == 0L || length(denominator) == 0L) return(TRUE)
  if (any(is.na(numerator)) || any(is.na(denominator))) return(TRUE)
  if (length(intersect(numerator, denominator)) > 0L) return(TRUE)
  FALSE
}

# Self-contained per-sample balance for ONE specimen over the frozen groups,
# using only features present in names(x) and strictly positive in x. Returns NA
# if either group has fewer than min_cov present-positive features (the caller
# maps NA to the neutral score 0). This is the ONE balance definition shared by
# fit (calibration) and score.
.bal_selbal_balance_row <- function(x, numerator, denominator, min_cov) {
  num_feat <- intersect(numerator, names(x))
  den_feat <- intersect(denominator, names(x))
  v_num <- x[num_feat]
  v_den <- x[den_feat]
  v_num <- v_num[is.finite(v_num) & v_num > 0]
  v_den <- v_den[is.finite(v_den) & v_den > 0]
  r <- length(v_num)
  s <- length(v_den)
  if (r < min_cov || s < min_cov) return(NA_real_)
  k <- sqrt((r * s) / (r + s))
  k * (mean(log(v_num)) - mean(log(v_den)))
}

# Vectorized training-balance over all rows; NA rows (below coverage) are kept as
# NA so the calibration can drop them.
.bal_selbal_balance_vector <- function(X, numerator, denominator, min_cov) {
  num_present <- intersect(numerator, colnames(X))
  den_present <- intersect(denominator, colnames(X))
  vapply(seq_len(nrow(X)), function(i) {
    x_i <- X[i, ]
    names(x_i) <- colnames(X)
    .bal_selbal_balance_row(x_i, num_present, den_present, min_cov)
  }, numeric(1))
}

# Frozen univariate logistic calibration of the self-contained balance: fit
# glm(y ~ B, binomial) on rows with a defined (non-NA) balance, store
# intercept+slope. If the balance is constant/degenerate (slope non-finite or
# zero, or too few usable rows), fall back to an orientation sign so larger is
# case-like, with slope as the sign and intercept 0.
.bal_selbal_calibrate <- function(B, y) {
  usable <- is.finite(B)
  orientation_sign <- .bal_selbal_orientation(B[usable], y[usable])
  if (sum(usable) < 3L ||
      length(unique(stats::na.omit(B[usable]))) < 2L ||
      length(unique(y[usable])) < 2L) {
    return(list(intercept = 0, slope = orientation_sign,
                orientation_sign = orientation_sign))
  }
  fit <- tryCatch(
    suppressWarnings(stats::glm(y[usable] ~ B[usable],
                                family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(intercept = 0, slope = orientation_sign,
                orientation_sign = orientation_sign))
  }
  co <- stats::coef(fit)
  intercept <- unname(co[1])
  slope <- unname(co[2])
  if (!is.finite(intercept) || !is.finite(slope) || slope == 0) {
    return(list(intercept = 0, slope = orientation_sign,
                orientation_sign = orientation_sign))
  }
  list(intercept = intercept, slope = slope,
       orientation_sign = orientation_sign)
}

.bal_selbal_orientation <- function(B, y) {
  if (length(B) == 0L || !any(y == 1L) || !any(y == 0L)) return(1)
  m1 <- mean(B[y == 1L])
  m0 <- mean(B[y == 0L])
  if (!is.finite(m1) || !is.finite(m0)) return(1)
  if (m1 >= m0) 1 else -1
}
