#' @title Fit Mondrian-conformal class-conditional LRT discriminator
#'
#' @description
#' Learns a frozen two-class discriminator from training data using
#' class-conditional conformal nonconformity scores in per-sample robust
#' centered log-ratio (rCLR) coordinates. Each class stores raw training anchors,
#' a centroid in within-sample rCLR space, and a calibration set of
#' class-specific nonconformity scores. No test data or test-time batch
#' statistics are used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{nonconformity}, either
#'   \code{"euclid"} (default) or \code{"mahalanobis"} for diagonal
#'   class-variance weighting; \code{max_anchors_per_class} positive integer
#'   anchor cap per class (default \code{200L}); \code{min_features} positive
#'   integer feature-overlap floor at scoring (default \code{3L}); \code{eps}
#'   positive finite variance floor used by the Mahalanobis path (default
#'   \code{1e-6}); and \code{seed} non-negative integer used only for
#'   deterministic anchor subsampling while restoring the global RNG state
#'   (default \code{1L}).
#'
#' @return A plain list of class \code{conf_mondrian_model} containing
#'   \code{feature_universe}, raw \code{case_anchors} and
#'   \code{control_anchors}, full-universe \code{repr} for inspection, and the
#'   resolved \code{hp}.
#'
#' @details
#' The conformal p-value is used here as a discriminative atypicality contrast,
#' not as an exact coverage-calibrated guarantee: the same training anchors
#' define the centroids and the calibration nonconformity scores, so
#' exchangeability is approximate.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
#' X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
#' model <- fit_conf_mondrian(X, y)
#' score <- score_conf_mondrian(model, X)
#' }
#'
#' @references
#' Vovk V, Gammerman A, Shafer G. (2005) \emph{Algorithmic Learning in a
#' Random World}. Springer.
#'
#' Shafer G, Vovk V. (2008) A tutorial on conformal prediction. \emph{Journal
#' of Machine Learning Research} 9: 371-421.
#'
#' @export
fit_conf_mondrian <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_conf_mondrian", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_conf_mondrian", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_conf_mondrian")
  hp <- .conf_mondrian_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_conf_mondrian: X_train must contain at least hp$min_features features")
  }

  anchor_idx <- .conf_mondrian_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  case_anchors <- X_train[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_train[anchor_idx$control, , drop = FALSE]

  # Full-universe representation for inspection/tests. Scoring re-derives the
  # same representation from frozen raw anchors over the present feature subset.
  repr <- .conf_mondrian_class_repr(case_anchors, control_anchors, hp)

  model <- list(
    feature_universe = colnames(X_train),
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    repr = repr,
    hp = hp
  )
  class(model) <- "conf_mondrian_model"
  model
}


#' @title Score Mondrian-conformal class-conditional LRT discriminator
#'
#' @description
#' Scores each row of \code{X} independently with a frozen
#' \code{\link{fit_conf_mondrian}} model. A specimen is transformed to
#' within-sample rCLR coordinates over the model features present in \code{X};
#' its class-specific nonconformity scores are compared with the frozen
#' class-conditional calibration scores; and the returned score is
#' \deqn{\log p_{case}(x) - \log p_{control}(x).}
#'
#' @param model A \code{conf_mondrian_model} object returned by
#'   \code{\link{fit_conf_mondrian}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like. If fewer than \code{model$hp$min_features} model
#'   features are present, the documented neutral score \code{0} is returned for
#'   every row.
#'
#' @details
#' At scoring time the present feature set is
#' \code{intersect(model$feature_universe, colnames(X))}. Class centroids,
#' diagonal variances, and calibration nonconformity scores are re-derived over
#' that present subset from the frozen raw anchors only. Each row's score then
#' uses only that row's own rCLR vector plus the frozen model-derived
#' representation, so no scored-batch statistic or cross-row coupling is used.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
#' X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
#' model <- fit_conf_mondrian(X, y)
#' score_conf_mondrian(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Vovk V, Gammerman A, Shafer G. (2005) \emph{Algorithmic Learning in a
#' Random World}. Springer.
#'
#' Shafer G, Vovk V. (2008) A tutorial on conformal prediction. \emph{Journal
#' of Machine Learning Research} 9: 371-421.
#'
#' @export
score_conf_mondrian <- function(model, X, meta = NULL) {
  if (!inherits(model, "conf_mondrian_model")) {
    stop("score_conf_mondrian: model must have class conf_mondrian_model")
  }
  X <- .reo_check_matrix(X, "score_conf_mondrian", "X")
  .reo_check_meta(meta, nrow(X), "score_conf_mondrian", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  repr <- .conf_mondrian_class_repr(
    model$case_anchors[, present, drop = FALSE],
    model$control_anchors[, present, drop = FALSE],
    model$hp
  )
  X_use <- X[, present, drop = FALSE]

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    x_i <- X_use[i, ]
    names(x_i) <- present
    z <- .conf_mondrian_rclr_row(x_i)
    A_case <- .conf_mondrian_nonconformity(
      z, repr$mu_case, repr$var_case, model$hp
    )
    A_control <- .conf_mondrian_nonconformity(
      z, repr$mu_control, repr$var_control, model$hp
    )
    p_case <- (1 + sum(repr$S_case >= A_case)) / (1 + length(repr$S_case))
    p_control <- (1 + sum(repr$S_control >= A_control)) /
      (1 + length(repr$S_control))
    log(p_case) - log(p_control)
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_conf_mondrian: scorer produced non-finite or wrong-length output")
  }
  out
}

.conf_mondrian_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_conf_mondrian: hp must be a list")
  allowed <- c("nonconformity", "max_anchors_per_class", "min_features",
               "eps", "seed")
  if (length(hp) > 0L && (is.null(names(hp)) || any(!nzchar(names(hp))))) {
    stop("fit_conf_mondrian: hp fields must be named")
  }
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_conf_mondrian: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  nonconformity <- hp$nonconformity
  if (is.null(nonconformity)) nonconformity <- "euclid"
  if (!is.character(nonconformity) || length(nonconformity) != 1L ||
      is.na(nonconformity) ||
      !nonconformity %in% c("euclid", "mahalanobis")) {
    stop("fit_conf_mondrian: hp$nonconformity must be 'euclid' or 'mahalanobis'")
  }

  max_anchors_per_class <- hp$max_anchors_per_class
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_conf_mondrian: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 3L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_conf_mondrian: hp$min_features must be a positive integer")
  }

  eps <- hp$eps
  if (is.null(eps)) eps <- 1e-6
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("fit_conf_mondrian: hp$eps must be a positive finite number")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_conf_mondrian: hp$seed must be a non-negative integer")
  }

  list(
    nonconformity = nonconformity,
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    eps = as.numeric(eps),
    seed = as.integer(seed)
  )
}

.conf_mondrian_anchor_indices <- function(y, max_anchors_per_class, seed) {
  case_idx <- which(y == 1L)
  control_idx <- which(y == 0L)
  needs_subsample <- length(case_idx) > max_anchors_per_class ||
    length(control_idx) > max_anchors_per_class
  if (!needs_subsample) {
    return(list(case = case_idx, control = control_idx))
  }

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
  if (length(case_idx) > max_anchors_per_class) {
    case_idx <- sort(sample(case_idx, max_anchors_per_class, replace = FALSE))
  }
  if (length(control_idx) > max_anchors_per_class) {
    control_idx <- sort(sample(control_idx, max_anchors_per_class,
                               replace = FALSE))
  }
  list(case = case_idx, control = control_idx)
}

.conf_mondrian_rclr_row <- function(v) {
  z <- numeric(length(v))
  pos <- which(v > 0)
  if (length(pos) > 0L) {
    g <- mean(log(v[pos]))
    z[pos] <- log(v[pos]) - g
  }
  names(z) <- names(v)
  z
}

.conf_mondrian_rclr_matrix <- function(X) {
  Z <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    .conf_mondrian_rclr_row(X[i, ])
  }))
  dimnames(Z) <- dimnames(X)
  storage.mode(Z) <- "double"
  Z
}

.conf_mondrian_feature_var <- function(Z, mu) {
  if (nrow(Z) < 2L) {
    out <- rep.int(0, ncol(Z))
    names(out) <- colnames(Z)
    return(out)
  }
  centered <- sweep(Z, 2L, mu, "-")
  out <- colSums(centered^2) / (nrow(Z) - 1L)
  names(out) <- colnames(Z)
  out
}

.conf_mondrian_nonconformity <- function(z, mu, var, hp) {
  delta2 <- (z - mu)^2
  if (identical(hp$nonconformity, "mahalanobis")) {
    return(sum(delta2 / (var + hp$eps)))
  }
  sum(delta2)
}

.conf_mondrian_class_scores <- function(Z, mu, var, hp) {
  vapply(seq_len(nrow(Z)), function(i) {
    .conf_mondrian_nonconformity(Z[i, ], mu, var, hp)
  }, numeric(1))
}

.conf_mondrian_class_repr <- function(case_anchors, control_anchors, hp) {
  Z_case <- .conf_mondrian_rclr_matrix(case_anchors)
  Z_control <- .conf_mondrian_rclr_matrix(control_anchors)

  mu_case <- colMeans(Z_case)
  mu_control <- colMeans(Z_control)
  var_case <- .conf_mondrian_feature_var(Z_case, mu_case)
  var_control <- .conf_mondrian_feature_var(Z_control, mu_control)
  S_case <- .conf_mondrian_class_scores(Z_case, mu_case, var_case, hp)
  S_control <- .conf_mondrian_class_scores(
    Z_control, mu_control, var_control, hp
  )

  list(
    mu_case = mu_case,
    mu_control = mu_control,
    var_case = var_case,
    var_control = var_control,
    S_case = S_case,
    S_control = S_control
  )
}
