#' @title Fit kernel mean-embedding witness-at-a-point discriminator
#'
#' @description
#' Learns a frozen two-class kernel mean-embedding witness function from
#' training data. Each specimen is represented by a within-sample compositional
#' transform: it is first closed to proportions \eqn{p = v / \sum v} and then
#' mapped to pseudocount centred log-ratio coordinates
#' \eqn{\phi(v) = \log(p + pc) - \mathrm{mean}(\log(p + pc))}. Closing to
#' proportions first makes \eqn{\phi} \emph{exactly} invariant to per-sample
#' scaling (library size / input amount) for any \code{pc}, because
#' \eqn{p(c v) = p(v)} for \eqn{c > 0}; a raw-abundance pseudocount log-ratio
#' \eqn{\log(v + pc) - \mathrm{mean}(\log(v + pc))} is only approximately
#' scale-invariant and so would not deliver library-size invariance at the
#' default pseudocount. The RBF bandwidth is either supplied as \code{hp$gamma}
#' or learned by the median heuristic on the stored training anchors:
#' \deqn{gamma = 1 / (2 * median(||phi(x_i) - phi(x_j)||^2)),}
#' using all anchor pairs; if the median squared distance is zero or non-finite,
#' \code{gamma = 1}. This bandwidth is panel-frozen: it is learned once over the
#' full training feature universe and applied unchanged at scoring, including
#' under partial feature overlap, because re-tuning \code{gamma} to the scored
#' overlap would require a scored-batch statistic and break single-sample
#' deployability. Fitting stores raw case/control anchor rows, the frozen
#' training feature universe, \code{pc}, \code{gamma}, and resolved
#' hyperparameters. No test data or test-time batch statistic is used.
#'
#' @param X_train Numeric matrix (samples \eqn{\times} features) of
#'   non-negative abundance values. Columns must be uniquely named feature ids.
#' @param y_train Integer/numeric 0/1 labels aligned to \code{X_train}; 1 is
#'   case/disease and 0 is control.
#' @param meta_train Optional per-sample metadata. Accepted for interface
#'   uniformity and ignored by this method.
#' @param hp List of frozen hyperparameters: \code{pc} positive pseudocount on
#'   the closed proportions (default \code{1e-3}; small so it bounds the
#'   log-ratio for absent/near-absent features without flattening the
#'   composition), \code{gamma} positive RBF bandwidth (default \code{NULL},
#'   median heuristic), \code{max_anchors_per_class} positive integer anchor cap
#'   per class (default \code{200L}), \code{min_features} positive integer
#'   feature-overlap floor at scoring (default \code{2L}), and \code{seed}
#'   non-negative integer used only for deterministic anchor subsampling while
#'   restoring the global RNG state (default \code{1L}).
#'
#' @return A plain list of class \code{kme_witness_model} containing
#'   \code{case_anchors}, \code{control_anchors}, \code{feature_universe},
#'   \code{gamma}, \code{pc}, and resolved \code{hp}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
#' X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
#' model <- fit_kme_witness(X, y)
#' score <- score_kme_witness(model, X)
#' }
#'
#' @references
#' Gretton A, Borgwardt KM, Rasch MJ, Scholkopf B, Smola A. (2012) A kernel
#' two-sample test. \emph{Journal of Machine Learning Research} 13: 723-773.
#'
#' Sutherland DJ, Tung HY, Strathmann H, De S, Ramdas A, Smola A, Gretton A.
#' (2017) Generative models and model criticism via optimized maximum mean
#' discrepancy. \emph{International Conference on Learning Representations}.
#'
#' @export
fit_kme_witness <- function(X_train, y_train, meta_train = NULL,
                            hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_kme_witness", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_kme_witness",
                  "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_kme_witness")
  hp <- .kme_witness_resolve_hp(hp)
  if (ncol(X_train) < hp$min_features) {
    stop("fit_kme_witness: X_train must contain at least hp$min_features features")
  }

  anchor_idx <- .kme_witness_anchor_indices(
    y = y,
    max_anchors_per_class = hp$max_anchors_per_class,
    seed = hp$seed
  )
  case_anchors <- X_train[anchor_idx$case, , drop = FALSE]
  control_anchors <- X_train[anchor_idx$control, , drop = FALSE]

  gamma <- hp$gamma
  if (is.null(gamma)) {
    anchor_clr <- .kme_witness_clr_matrix(
      rbind(case_anchors, control_anchors),
      pc = hp$pc
    )
    gamma <- .kme_witness_median_gamma(anchor_clr)
    hp$gamma <- gamma
  }

  model <- list(
    case_anchors = case_anchors,
    control_anchors = control_anchors,
    feature_universe = colnames(X_train),
    gamma = gamma,
    pc = hp$pc,
    hp = hp
  )
  class(model) <- "kme_witness_model"
  model
}


#' @title Score kernel mean-embedding witness-at-a-point discriminator
#'
#' @description
#' Scores each row of \code{X} independently with the frozen model learned by
#' \code{\link{fit_kme_witness}}. For one specimen \eqn{x}, the score is the
#' witness function evaluated at that point:
#' \deqn{mean_i k(phi(x), phi(a_i^{case})) -
#'       mean_j k(phi(x), phi(a_j^{control})),}
#' where \eqn{k(a,b) = exp(-gamma ||a-b||^2)}. At scoring time, the present
#' feature set is \code{intersect(model$feature_universe, colnames(X))}. Both
#' the specimen and the frozen raw anchors are closed to proportions over
#' exactly that present set and then centred-log-ratio transformed, which keeps
#' partial-overlap distances consistent (and the score exactly invariant to
#' per-sample scaling) without using any scored-batch statistic. The frozen
#' bandwidth \code{model$gamma} is applied unchanged in this reduced
#' present-feature space; it is not re-tuned to the overlap, which is what keeps
#' scoring single-sample. Re-deriving the anchor representation over
#' \code{present} uses only frozen training data, so the score depends only on
#' the specimen and the frozen model. If fewer than
#' \code{model$hp$min_features} features are present, the documented neutral
#' score \code{0} is returned for every row.
#'
#' @param model A \code{kme_witness_model} object returned by
#'   \code{\link{fit_kme_witness}}.
#' @param X Numeric matrix (samples \eqn{\times} features) of non-negative
#'   abundance values. Columns must be named feature ids.
#' @param meta Optional per-sample metadata. Accepted for interface uniformity
#'   and ignored by this method.
#'
#' @return Plain finite numeric vector of length \code{nrow(X)}. Larger values
#'   are more case-like because the witness is case mean embedding minus
#'   control mean embedding. Scoring uses only each row's own values plus the
#'   frozen anchors, \code{gamma}, and \code{pc}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
#'             dimnames = list(NULL, paste0("miR-", seq_len(30))))
#' y <- rep(c(0, 1), each = 30)
#' X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
#' X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
#' model <- fit_kme_witness(X, y)
#' score_kme_witness(model, X[1, , drop = FALSE])
#' }
#'
#' @references
#' Gretton A, Borgwardt KM, Rasch MJ, Scholkopf B, Smola A. (2012) A kernel
#' two-sample test. \emph{Journal of Machine Learning Research} 13: 723-773.
#'
#' Sutherland DJ, Tung HY, Strathmann H, De S, Ramdas A, Smola A, Gretton A.
#' (2017) Generative models and model criticism via optimized maximum mean
#' discrepancy. \emph{International Conference on Learning Representations}.
#'
#' @export
score_kme_witness <- function(model, X, meta = NULL) {
  if (!inherits(model, "kme_witness_model")) {
    stop("score_kme_witness: model must have class kme_witness_model")
  }
  X <- .reo_check_matrix(X, "score_kme_witness", "X")
  .reo_check_meta(meta, nrow(X), "score_kme_witness", "meta")

  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) < model$hp$min_features) {
    return(rep(0, nrow(X)))
  }

  X_use <- X[, present, drop = FALSE]
  case_phi <- .kme_witness_clr_matrix(
    model$case_anchors[, present, drop = FALSE],
    pc = model$pc
  )
  control_phi <- .kme_witness_clr_matrix(
    model$control_anchors[, present, drop = FALSE],
    pc = model$pc
  )

  out <- vapply(seq_len(nrow(X_use)), function(i) {
    # Restore feature names: single-row slicing of a 1-column matrix drops them,
    # which would make future name-based row helpers fragile.
    x_i <- X_use[i, ]
    names(x_i) <- present
    phi_i <- .kme_witness_clr_vector(x_i, pc = model$pc)
    .kme_witness_score_phi(
      phi = phi_i,
      case_phi = case_phi,
      control_phi = control_phi,
      gamma = model$gamma
    )
  }, numeric(1))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_kme_witness: scorer produced non-finite or wrong-length output")
  }
  out
}

.kme_witness_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_kme_witness: hp must be a list")
  allowed <- c("pc", "gamma", "max_anchors_per_class", "min_features",
               "seed")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_kme_witness: unknown hp field(s): ",
         paste(unknown, collapse = ", "))
  }

  pc <- hp$pc
  if (is.null(pc)) pc <- 1e-3
  if (!is.numeric(pc) || length(pc) != 1L || !is.finite(pc) || pc <= 0) {
    stop("fit_kme_witness: hp$pc must be a positive finite number")
  }

  gamma <- hp$gamma
  if (!is.null(gamma) &&
      (!is.numeric(gamma) || length(gamma) != 1L ||
       !is.finite(gamma) || gamma <= 0)) {
    stop("fit_kme_witness: hp$gamma must be NULL or a positive finite number")
  }

  max_anchors_per_class <- hp$max_anchors_per_class
  if (is.null(max_anchors_per_class)) max_anchors_per_class <- 200L
  if (!is.numeric(max_anchors_per_class) ||
      length(max_anchors_per_class) != 1L ||
      !is.finite(max_anchors_per_class) || max_anchors_per_class < 1L ||
      max_anchors_per_class > .Machine$integer.max ||
      max_anchors_per_class != as.integer(max_anchors_per_class)) {
    stop("fit_kme_witness: hp$max_anchors_per_class must be a positive integer")
  }

  min_features <- hp$min_features
  if (is.null(min_features)) min_features <- 2L
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      !is.finite(min_features) || min_features < 1L ||
      min_features > .Machine$integer.max ||
      min_features != as.integer(min_features)) {
    stop("fit_kme_witness: hp$min_features must be a positive integer")
  }

  seed <- hp$seed
  if (is.null(seed)) seed <- 1L
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("fit_kme_witness: hp$seed must be a non-negative integer")
  }

  list(
    pc = as.numeric(pc),
    gamma = if (is.null(gamma)) NULL else as.numeric(gamma),
    max_anchors_per_class = as.integer(max_anchors_per_class),
    min_features = as.integer(min_features),
    seed = as.integer(seed)
  )
}

.kme_witness_anchor_indices <- function(y, max_anchors_per_class, seed) {
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

# Within-sample representation: close each specimen to proportions, then apply a
# pseudocount centred log-ratio. Closing to proportions first (p = v / sum(v))
# makes the representation EXACTLY invariant to per-sample scaling for any pc,
# because p(c*v) = p(v) for c > 0. A raw-abundance pseudocount CLR
# (log(v + pc) - mean(log(v + pc))) is only approximately scale-invariant (it
# needs v >> pc), so it does not deliver library-size invariance at the default
# pseudocount -- exact library-size invariance is the property that makes a
# within-sample compositional score robust to the per-sample normalization
# differences this method bank exists to avoid correcting at inference. An
# all-zero specimen maps to the uniform composition, whose CLR is the zero vector.
.kme_witness_clr_matrix <- function(X, pc) {
  rs <- rowSums(X)
  safe <- ifelse(rs > 0, rs, 1)
  P <- X / safe
  if (any(rs <= 0)) P[rs <= 0, ] <- 1 / ncol(X)
  L <- log(P + pc)
  sweep(L, 1L, rowMeans(L), "-")
}

.kme_witness_clr_vector <- function(x, pc) {
  s <- sum(x)
  p <- if (s > 0) x / s else rep.int(1 / length(x), length(x))
  l <- log(p + pc)
  as.numeric(l - mean(l))
}

.kme_witness_median_gamma <- function(anchor_clr) {
  d2 <- as.numeric(stats::dist(anchor_clr, method = "euclidean"))^2
  med <- stats::median(d2, na.rm = TRUE)
  if (!is.finite(med) || med <= 0) return(1)
  1 / (2 * med)
}

.kme_witness_score_phi <- function(phi, case_phi, control_phi, gamma) {
  case_d2 <- rowSums(sweep(case_phi, 2L, phi, "-")^2)
  control_d2 <- rowSums(sweep(control_phi, 2L, phi, "-")^2)
  mean(exp(-gamma * case_d2)) - mean(exp(-gamma * control_d2))
}
