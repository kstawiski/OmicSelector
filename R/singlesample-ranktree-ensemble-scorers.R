#' External rank-tree competitor registry
#'
#' Returns the two maintained `ranktreeEnsemble` algorithms evaluated as a
#' separately multiplicity-controlled external-competitor sensitivity for the
#' OmicSelector paper. They are deliberately not appended to
#' [singlesample_method_roster()], whose 74 rows and 50-route within-cohort
#' family were frozen before this direct competitor was identified.
#'
#' @return A data frame with one row per external competitor and executable
#'   fit/score routing metadata.
#' @export
singlesample_external_competitor_roster <- function() {
  data.frame(
    method_id = c("reo-rankforest", "reo-rankboost"),
    family = "K",
    estimand = "within",
    role = "external_sensitivity",
    tier = "R2",
    dep_route = "ranktreeEnsemble-0.24",
    fit_fn = c("fit_reo_rankforest", "fit_reo_rankboost"),
    score_fn = c("score_reo_rankforest", "score_reo_rankboost"),
    pkg_status = "new-secondary",
    row_source = "per_cohort_rows",
    lopbo_mechanism = "within_cohort",
    primary_family = FALSE,
    notes = c(
      "random rank forest; post-freeze pre-outcome external sensitivity",
      "boosted rank trees; post-freeze pre-outcome external sensitivity"
    ),
    stringsAsFactors = FALSE
  )
}

.reo_ranktree_require_backend <- function() {
  if (!requireNamespace("ranktreeEnsemble", quietly = TRUE)) {
    stop(
      "ranktreeEnsemble 0.24 is required; install it from CRAN to use the ",
      "rank-tree external competitors.",
      call. = FALSE
    )
  }
  version <- as.character(utils::packageVersion("ranktreeEnsemble"))
  if (!identical(version, "0.24")) {
    stop(
      "The audited rank-tree adapter requires ranktreeEnsemble 0.24; found ",
      version, ". Re-audit the backend API before changing this pin.",
      call. = FALSE
    )
  }
  version
}

.reo_ranktree_check_training <- function(X_train, y_train, meta_train, fn) {
  X_train <- .reo_check_matrix(X_train, fn, "X_train")
  .reo_check_meta(meta_train, nrow(X_train), fn, "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), fn)
  if (ncol(X_train) < 5L) {
    stop(fn, ": at least five features are required for the frozen ",
         "ranktreeEnsemble dimension-reduction contract", call. = FALSE)
  }
  class_n <- table(factor(y, levels = c(0L, 1L)))
  if (any(class_n < 5L)) {
    stop(fn, ": at least five training rows per class are required",
         call. = FALSE)
  }
  list(X = X_train, y = y)
}

.reo_ranktree_safe_training_data <- function(X, y) {
  safe_names <- sprintf("feature_%07d", seq_len(ncol(X)))
  X_safe <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
  names(X_safe) <- safe_names
  X_safe$outcome <- factor(
    ifelse(y == 1L, "case", "control"),
    levels = c("control", "case")
  )
  list(data = X_safe, safe_names = safe_names)
}

.reo_ranktree_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_rf_cores <- getOption("rf.cores")
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
    if (is.null(old_rf_cores)) {
      options(rf.cores = NULL)
    } else {
      options(rf.cores = old_rf_cores)
    }
  }, add = TRUE)
  # ranktreeEnsemble 0.24 does not pass its seed to the preliminary
  # randomForestSRC dimension-reduction fit. Its parallel RNG is not
  # reproducible from R's seed alone, whereas its single-core route is.
  options(rf.cores = 1L)
  set.seed(seed)
  eval.parent(substitute(expr))
}

.reo_ranktree_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != as.integer(x)) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  as.integer(x)
}

.reo_ranktree_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || x > 1) {
    stop(name, " must be in (0, 1]", call. = FALSE)
  }
  as.numeric(x)
}

.reo_ranktree_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  x
}

.reo_rankforest_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_reo_rankforest: hp must be a list", call. = FALSE)
  allowed <- c("ntree", "seed", "dimreduce", "class_balance")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown)) {
    stop("fit_reo_rankforest: unknown hp field(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  list(
    ntree = .reo_ranktree_positive_integer(hp$ntree %||% 500L,
                                           "fit_reo_rankforest: hp$ntree"),
    seed = .reo_ranktree_positive_integer(hp$seed %||% 1L,
                                          "fit_reo_rankforest: hp$seed"),
    dimreduce = .reo_ranktree_flag(hp$dimreduce %||% TRUE,
                                   "fit_reo_rankforest: hp$dimreduce"),
    class_balance = .reo_ranktree_flag(
      hp$class_balance %||% TRUE,
      "fit_reo_rankforest: hp$class_balance"
    ),
    datrank = FALSE,
    rf_cores = 1L
  )
}

.reo_rankboost_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_reo_rankboost: hp must be a list", call. = FALSE)
  allowed <- c("ntree", "seed", "dimreduce", "nodedepth", "nodesize",
               "shrinkage", "bag_fraction")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown)) {
    stop("fit_reo_rankboost: unknown hp field(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  list(
    ntree = .reo_ranktree_positive_integer(hp$ntree %||% 500L,
                                           "fit_reo_rankboost: hp$ntree"),
    seed = .reo_ranktree_positive_integer(hp$seed %||% 1L,
                                          "fit_reo_rankboost: hp$seed"),
    dimreduce = .reo_ranktree_flag(hp$dimreduce %||% TRUE,
                                   "fit_reo_rankboost: hp$dimreduce"),
    nodedepth = .reo_ranktree_positive_integer(
      hp$nodedepth %||% 3L, "fit_reo_rankboost: hp$nodedepth"
    ),
    nodesize = .reo_ranktree_positive_integer(
      hp$nodesize %||% 5L, "fit_reo_rankboost: hp$nodesize"
    ),
    shrinkage = .reo_ranktree_probability(
      hp$shrinkage %||% 0.05, "fit_reo_rankboost: hp$shrinkage"
    ),
    bag_fraction = .reo_ranktree_probability(
      hp$bag_fraction %||% 0.5, "fit_reo_rankboost: hp$bag_fraction"
    ),
    datrank = FALSE
  )
}

#' Fit a random rank forest external competitor
#'
#' Fits the random-forest algorithm from `ranktreeEnsemble` 0.24 entirely on
#' training-fold data. The audited adapter deliberately fixes `datrank=FALSE`:
#' version 0.24 otherwise converts training columns to across-sample ranks but
#' predicts from raw within-row pairs, giving different train and deployment
#' representations. Feature names are replaced by a frozen collision-free map
#' before the dependency constructs pair indicators. The preliminary
#' `randomForestSRC` preliminary screen is unweighted because version 0.24
#' exposes no screen-weight hook; the downstream pair forest uses
#' inverse-frequency case weights, with outcome level `case` (numeric target 1)
#' as the positive class. Screening is confined to `rf.cores=1` because version
#' 0.24 does not forward its seed to that fit; the caller's prior option is
#' restored.
#'
#' @param X_train Numeric non-negative matrix of training samples by features.
#' @param y_train Numeric/integer 0/1 training labels; 1 is the case class.
#' @param meta_train Optional training metadata, checked for row alignment and
#'   otherwise ignored.
#' @param hp Optional fixed hyperparameters: `ntree` (500), `seed` (1),
#'   `dimreduce` (`TRUE`) and `class_balance` (`TRUE`).
#' @return A frozen `reo_rankforest_model` suitable for singleton scoring.
#' @references Lu M, Yin R, Chen XS. Ensemble methods of rank-based trees for
#'   single sample classification with gene expression profiles. Journal of
#'   Translational Medicine. 2024;22:140. doi:10.1186/s12967-024-04940-2.
#' @export
fit_reo_rankforest <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  backend_version <- .reo_ranktree_require_backend()
  checked <- .reo_ranktree_check_training(
    X_train, y_train, meta_train, "fit_reo_rankforest"
  )
  hp <- .reo_rankforest_resolve_hp(hp)
  safe <- .reo_ranktree_safe_training_data(checked$X, checked$y)
  case_wt <- NULL
  if (hp$class_balance) {
    counts <- table(factor(checked$y, levels = c(0L, 1L)))
    case_wt <- nrow(checked$X) / (2 * as.numeric(counts[checked$y + 1L]))
  }
  backend <- .reo_ranktree_with_seed(hp$seed, {
    ranktreeEnsemble::rforest(
      outcome ~ ., data = safe$data, dimreduce = hp$dimreduce,
      datrank = FALSE, ntree = hp$ntree, case.wt = case_wt,
      forest = TRUE, seed = -abs(hp$seed)
    )
  })
  model <- list(
    backend = backend,
    feature_names = colnames(checked$X),
    safe_feature_names = safe$safe_names,
    positive_class = "case",
    backend_package = "ranktreeEnsemble",
    backend_version = backend_version,
    hp = hp
  )
  class(model) <- "reo_rankforest_model"
  model
}

#' Fit boosted rank trees external competitor
#'
#' Fits the boosting algorithm from `ranktreeEnsemble` 0.24 using the same
#' fold-local, collision-free feature mapping and deployable `datrank=FALSE`
#' contract as [fit_reo_rankforest()]. The backend's exposed `weights` argument
#' is not forwarded by version 0.24 and is therefore not represented as active
#' in this adapter.
#'
#' @inheritParams fit_reo_rankforest
#' @param hp Optional fixed hyperparameters: `ntree` (500), `seed` (1),
#'   `dimreduce` (`TRUE`), `nodedepth` (3), `nodesize` (5), `shrinkage`
#'   (0.05), and `bag_fraction` (0.5).
#' @return A frozen `reo_rankboost_model` suitable for singleton scoring.
#' @export
fit_reo_rankboost <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  backend_version <- .reo_ranktree_require_backend()
  checked <- .reo_ranktree_check_training(
    X_train, y_train, meta_train, "fit_reo_rankboost"
  )
  hp <- .reo_rankboost_resolve_hp(hp)
  safe <- .reo_ranktree_safe_training_data(checked$X, checked$y)
  backend <- .reo_ranktree_with_seed(hp$seed, {
    ranktreeEnsemble::rboost(
      outcome ~ ., data = safe$data, dimreduce = hp$dimreduce,
      datrank = FALSE, ntree = hp$ntree, nodedepth = hp$nodedepth,
      nodesize = hp$nodesize, shrinkage = hp$shrinkage,
      bag.fraction = hp$bag_fraction, train.fraction = 1,
      cv.folds = 1, keep.data = TRUE, verbose = FALSE, n.cores = 1L
    )
  })
  model <- list(
    backend = backend,
    feature_names = colnames(checked$X),
    safe_feature_names = safe$safe_names,
    positive_class = "case",
    backend_package = "ranktreeEnsemble",
    backend_version = backend_version,
    hp = hp
  )
  class(model) <- "reo_rankboost_model"
  model
}

.reo_ranktree_score_data <- function(model, X, fn) {
  X <- .reo_check_matrix(X, fn, "X")
  if (!setequal(colnames(X), model$feature_names) ||
      ncol(X) != length(model$feature_names)) {
    stop(fn, ": X must contain the complete frozen training feature universe",
         call. = FALSE)
  }
  X <- X[, model$feature_names, drop = FALSE]
  out <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- model$safe_feature_names
  out
}

.reo_ranktree_extract_case_probability <- function(value, n, positive_class, fn) {
  if (is.array(value) && length(dim(value)) == 3L) {
    value_dim <- dim(value)
    value_names <- dimnames(value)
    if (value_dim[[3L]] != 1L || is.null(value_names[[2L]])) {
      stop(fn, ": backend returned an unsupported probability array",
           call. = FALSE)
    }
    value <- matrix(
      value[, , 1L], nrow = value_dim[[1L]], ncol = value_dim[[2L]],
      dimnames = list(value_names[[1L]], value_names[[2L]])
    )
  } else if (is.array(value) && length(dim(value)) > 3L) {
    stop(fn, ": backend returned an unsupported probability array",
         call. = FALSE)
  }
  if (is.matrix(value) || is.data.frame(value)) {
    value <- as.matrix(value)
    if (nrow(value) != n || is.null(colnames(value)) ||
        !(positive_class %in% colnames(value))) {
      stop(fn, ": backend did not return a named case-probability column",
           call. = FALSE)
    }
    out <- value[, positive_class]
  } else if (is.numeric(value) && length(value) == n) {
    out <- value
  } else {
    stop(fn, ": backend returned an unsupported probability shape",
         call. = FALSE)
  }
  out <- as.numeric(out)
  if (length(out) != n || any(!is.finite(out)) ||
      any(out < 0 | out > 1)) {
    stop(fn, ": backend case probabilities must be finite and in [0, 1]",
         call. = FALSE)
  }
  out
}

#' Score a random rank forest one specimen at a time
#'
#' @param model A `reo_rankforest_model` returned by [fit_reo_rankforest()].
#' @param X Numeric matrix of samples by the complete frozen feature universe.
#' @param meta Optional row-aligned metadata, accepted for canonical dispatch
#'   and otherwise ignored.
#' @return Numeric case probabilities, one per row of `X`.
#' @export
score_reo_rankforest <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_rankforest_model")) {
    stop("score_reo_rankforest: model must have class reo_rankforest_model",
         call. = FALSE)
  }
  .reo_ranktree_require_backend()
  data <- .reo_ranktree_score_data(model, X, "score_reo_rankforest")
  .reo_check_meta(meta, nrow(data), "score_reo_rankforest", "meta")
  pred <- ranktreeEnsemble::predict(model$backend, newdata = data)
  .reo_ranktree_extract_case_probability(
    pred$value, nrow(data), model$positive_class, "score_reo_rankforest"
  )
}

#' Score boosted rank trees one specimen at a time
#'
#' @param model A `reo_rankboost_model` returned by [fit_reo_rankboost()].
#' @inheritParams score_reo_rankforest
#' @return Numeric case probabilities, one per row of `X`.
#' @export
score_reo_rankboost <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_rankboost_model")) {
    stop("score_reo_rankboost: model must have class reo_rankboost_model",
         call. = FALSE)
  }
  .reo_ranktree_require_backend()
  data <- .reo_ranktree_score_data(model, X, "score_reo_rankboost")
  .reo_check_meta(meta, nrow(data), "score_reo_rankboost", "meta")
  pred <- ranktreeEnsemble::predict(
    model$backend, newdata = data, n.trees = model$hp$ntree,
    type = "response"
  )
  .reo_ranktree_extract_case_probability(
    pred$value, nrow(data), model$positive_class, "score_reo_rankboost"
  )
}
