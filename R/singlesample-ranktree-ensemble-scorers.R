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

.reo_ranktree_dimreduce_quantile <- function(n_features, max_screen_features,
                                             enabled) {
  n_features <- .reo_ranktree_positive_integer(
    n_features, "rank-tree training feature count"
  )
  max_screen_features <- .reo_ranktree_positive_integer(
    max_screen_features, "rank-tree maximum screened features"
  )
  if (!.reo_ranktree_flag(enabled, "rank-tree dimension reduction flag")) {
    return(FALSE)
  }
  max(0.75, 1 - max_screen_features / n_features)
}

.reo_ranktree_pair_order_varies <- function(X) {
  if (nrow(X) < 2L) return(FALSE)
  feature_index <- seq_len(ncol(X))
  # ranktreeEnsemble::convert_genepairs emits one bit for each i < j:
  # I(X_i < X_j). Sorting by value and breaking ties by decreasing feature
  # index gives the unique permutation whose pair precedences are exactly those
  # backend bits: for a tie, j precedes i and the i<j bit remains zero.
  reference <- order(X[1L, ], -feature_index, method = "radix")
  for (i in 2:nrow(X)) {
    signature <- order(X[i, ], -feature_index, method = "radix")
    if (!identical(signature, reference)) return(TRUE)
  }
  FALSE
}

.reo_rankboost_effective_bag_fraction <- function(n_train, nodesize,
                                                   requested) {
  n_train <- .reo_ranktree_positive_integer(
    n_train, "fit_reo_rankboost: training row count"
  )
  nodesize <- .reo_ranktree_positive_integer(
    nodesize, "fit_reo_rankboost: hp$nodesize"
  )
  requested <- .reo_ranktree_probability(
    requested, "fit_reo_rankboost: hp$bag_fraction"
  )
  effective <- max(requested, (2 * nodesize + 2) / n_train)
  if (effective > 1) {
    stop(
      "fit_reo_rankboost: training partition is too small for nodesize=",
      nodesize, "; need at least ", 2 * nodesize + 2,
      " rows to satisfy the pinned gbm bagging inequality",
      call. = FALSE
    )
  }
  effective
}

.reo_ranktree_intercept_only_model <- function(checked, hp, backend_version,
                                                model_class) {
  safe_names <- sprintf("feature_%07d", seq_len(ncol(checked$X)))
  model <- list(
    backend = NULL,
    feature_names = colnames(checked$X),
    safe_feature_names = safe_names,
    backend_safe_feature_names = character(0),
    positive_class = "case",
    backend_package = "ranktreeEnsemble",
    backend_version = backend_version,
    constant_probability = mean(checked$y),
    fit_status = "intercept_only_constant_pair_order",
    hp = hp
  )
  class(model) <- model_class
  model
}

.reo_ranktree_intercept_only_score <- function(model, X, meta, fn) {
  if (!identical(model$fit_status, "intercept_only_constant_pair_order")) {
    return(NULL)
  }
  X <- .reo_check_matrix(X, fn, "X")
  .reo_check_meta(meta, nrow(X), fn, "meta")
  probability <- model$constant_probability
  if (!is.numeric(probability) || length(probability) != 1L ||
      !is.finite(probability) || probability < 0 || probability > 1 ||
      !is.null(model$backend) || length(model$backend_safe_feature_names)) {
    stop(fn, ": malformed intercept-only rank-tree state", call. = FALSE)
  }
  rep(as.numeric(probability), nrow(X))
}

.reo_ranktree_backend_pair_names <- function(backend, fn) {
  pair_names <- if (inherits(backend, "gbm")) {
    backend$var.names
  } else if (inherits(backend, "rfsrc")) {
    backend$xvar.names
  } else {
    NULL
  }
  if (!is.character(pair_names) || !length(pair_names) ||
      anyNA(pair_names) || any(!nzchar(pair_names)) ||
      anyDuplicated(pair_names)) {
    stop(fn, ": backend does not expose a unique pair-feature contract",
         call. = FALSE)
  }
  pair_names
}

.reo_ranktree_backend_feature_order <- function(backend, safe_names, fn) {
  pair_names <- .reo_ranktree_backend_pair_names(backend, fn)
  parts <- strsplit(pair_names, "_less_than_", fixed = TRUE)
  if (any(lengths(parts) != 2L)) {
    stop(fn, ": backend pair-feature names cannot be decoded", call. = FALSE)
  }
  edges <- do.call(rbind, parts)
  features <- unique(as.vector(t(edges)))
  if (length(features) < 2L || !all(features %in% safe_names) ||
      nrow(edges) != choose(length(features), 2L) ||
      anyDuplicated(apply(edges, 1L, function(x) paste(sort(x), collapse = "\r")))) {
    stop(fn, ": backend pair-feature set is incomplete or inconsistent",
         call. = FALSE)
  }
  incoming <- setNames(integer(length(features)), features)
  for (feature in edges[, 2L]) incoming[[feature]] <- incoming[[feature]] + 1L
  if (!setequal(unname(incoming), seq.int(0L, length(features) - 1L))) {
    stop(fn, ": backend pair directions do not define one feature order",
         call. = FALSE)
  }
  ordered <- names(sort(incoming, method = "radix"))
  expected <- character(0)
  for (i in seq_len(length(ordered) - 1L)) {
    expected <- c(
      expected,
      paste0(ordered[[i]], "_less_than_", ordered[(i + 1L):length(ordered)])
    )
  }
  if (!setequal(pair_names, expected)) {
    stop(fn, ": backend pair-feature names differ from the recovered order",
         call. = FALSE)
  }
  ordered
}

.reo_rankforest_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_reo_rankforest: hp must be a list", call. = FALSE)
  allowed <- c("ntree", "seed", "dimreduce", "class_balance",
               "max_screen_features")
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
    max_screen_features = .reo_ranktree_positive_integer(
      hp$max_screen_features %||% 128L,
      "fit_reo_rankforest: hp$max_screen_features"
    ),
    datrank = FALSE,
    rf_cores = 1L
  )
}

.reo_rankboost_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_reo_rankboost: hp must be a list", call. = FALSE)
  allowed <- c("ntree", "seed", "dimreduce", "nodedepth", "nodesize",
               "shrinkage", "bag_fraction", "max_screen_features")
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
    max_screen_features = .reo_ranktree_positive_integer(
      hp$max_screen_features %||% 128L,
      "fit_reo_rankboost: hp$max_screen_features"
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
#' `randomForestSRC` screen is unweighted because version 0.24
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
#'   `dimreduce` (`TRUE`), `class_balance` (`TRUE`), and
#'   `max_screen_features` (128). The last value preserves the backend's
#'   top-quartile importance screen unless it would be expected to materialize
#'   more than 128 preliminary features before pair expansion. The realized
#'   count can differ when preliminary importance values are tied. If no
#'   pair-order indicator varies in the training data, fitting returns an exact
#'   intercept-only state instead of invoking a backend that cannot split it.
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
  hp$backend_dimreduce_quantile <- .reo_ranktree_dimreduce_quantile(
    ncol(checked$X), hp$max_screen_features, hp$dimreduce
  )
  if (!.reo_ranktree_pair_order_varies(checked$X)) {
    return(.reo_ranktree_intercept_only_model(
      checked, hp, backend_version, "reo_rankforest_model"
    ))
  }
  safe <- .reo_ranktree_safe_training_data(checked$X, checked$y)
  case_wt <- NULL
  if (hp$class_balance) {
    counts <- table(factor(checked$y, levels = c(0L, 1L)))
    case_wt <- nrow(checked$X) / (2 * as.numeric(counts[checked$y + 1L]))
  }
  backend <- .reo_ranktree_with_seed(hp$seed, {
    ranktreeEnsemble::rforest(
      outcome ~ ., data = safe$data, dimreduce = hp$backend_dimreduce_quantile,
      datrank = FALSE, ntree = hp$ntree, case.wt = case_wt,
      forest = TRUE, seed = -abs(hp$seed)
    )
  })
  backend_safe_feature_names <- .reo_ranktree_backend_feature_order(
    backend, safe$safe_names, "fit_reo_rankforest"
  )
  model <- list(
    backend = backend,
    feature_names = colnames(checked$X),
    safe_feature_names = safe$safe_names,
    backend_safe_feature_names = backend_safe_feature_names,
    positive_class = "case",
    backend_package = "ranktreeEnsemble",
    backend_version = backend_version,
    fit_status = "backend",
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
#'   (0.05), requested `bag_fraction` (0.5), and `max_screen_features` (128).
#'   The effective bag fraction applies a conservative training-size-only
#'   one-observation margin above the pinned `gbm` node-size boundary; both the
#'   requested and effective values are retained in the fitted state. A
#'   training set with no variable pair ordering returns an exact
#'   intercept-only state.
#' @return A frozen `reo_rankboost_model` suitable for singleton scoring.
#' @export
fit_reo_rankboost <- function(X_train, y_train, meta_train = NULL, hp = list()) {
  backend_version <- .reo_ranktree_require_backend()
  checked <- .reo_ranktree_check_training(
    X_train, y_train, meta_train, "fit_reo_rankboost"
  )
  hp <- .reo_rankboost_resolve_hp(hp)
  hp$backend_dimreduce_quantile <- .reo_ranktree_dimreduce_quantile(
    ncol(checked$X), hp$max_screen_features, hp$dimreduce
  )
  hp$requested_bag_fraction <- hp$bag_fraction
  if (!.reo_ranktree_pair_order_varies(checked$X)) {
    hp$effective_bag_fraction <- NA_real_
    return(.reo_ranktree_intercept_only_model(
      checked, hp, backend_version, "reo_rankboost_model"
    ))
  }
  hp$effective_bag_fraction <- .reo_rankboost_effective_bag_fraction(
    nrow(checked$X), hp$nodesize, hp$requested_bag_fraction
  )
  safe <- .reo_ranktree_safe_training_data(checked$X, checked$y)
  backend <- .reo_ranktree_with_seed(hp$seed, {
    ranktreeEnsemble::rboost(
      outcome ~ ., data = safe$data, dimreduce = hp$backend_dimreduce_quantile,
      datrank = FALSE, ntree = hp$ntree, nodedepth = hp$nodedepth,
      nodesize = hp$nodesize, shrinkage = hp$shrinkage,
      bag.fraction = hp$effective_bag_fraction, train.fraction = 1,
      cv.folds = 1, keep.data = TRUE, verbose = FALSE, n.cores = 1L
    )
  })
  backend_safe_feature_names <- .reo_ranktree_backend_feature_order(
    backend, safe$safe_names, "fit_reo_rankboost"
  )
  model <- list(
    backend = backend,
    feature_names = colnames(checked$X),
    safe_feature_names = safe$safe_names,
    backend_safe_feature_names = backend_safe_feature_names,
    positive_class = "case",
    backend_package = "ranktreeEnsemble",
    backend_version = backend_version,
    fit_status = "backend",
    hp = hp
  )
  class(model) <- "reo_rankboost_model"
  model
}

.reo_ranktree_score_data <- function(model, X, fn) {
  X <- .reo_check_matrix(X, fn, "X")
  selected <- model$backend_safe_feature_names
  if (!is.character(model$feature_names) ||
      !is.character(model$safe_feature_names) ||
      length(model$feature_names) != length(model$safe_feature_names) ||
      anyDuplicated(model$feature_names) || anyDuplicated(model$safe_feature_names) ||
      !is.character(selected) || length(selected) < 2L || anyNA(selected) ||
      any(!selected %in% model$safe_feature_names) || anyDuplicated(selected)) {
    stop(fn, ": model lacks a valid frozen backend feature subset",
         call. = FALSE)
  }
  selected_index <- match(selected, model$safe_feature_names)
  required_features <- model$feature_names[selected_index]
  missing <- setdiff(required_features, colnames(X))
  if (length(missing)) {
    stop(fn, ": X lacks frozen backend feature(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- as.data.frame(
    X[, required_features, drop = FALSE],
    check.names = FALSE, stringsAsFactors = FALSE
  )
  names(out) <- selected
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
#' @param X Numeric matrix of samples containing every raw feature selected and
#'   frozen by the fitted backend. Unselected training features are optional.
#'   An intercept-only state has no selected raw features.
#' @param meta Optional row-aligned metadata, accepted for canonical dispatch
#'   and otherwise ignored.
#' @return Numeric case probabilities, one per row of `X`.
#' @export
score_reo_rankforest <- function(model, X, meta = NULL) {
  if (!inherits(model, "reo_rankforest_model")) {
    stop("score_reo_rankforest: model must have class reo_rankforest_model",
         call. = FALSE)
  }
  constant <- .reo_ranktree_intercept_only_score(
    model, X, meta, "score_reo_rankforest"
  )
  if (!is.null(constant)) return(constant)
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
  constant <- .reo_ranktree_intercept_only_score(
    model, X, meta, "score_reo_rankboost"
  )
  if (!is.null(constant)) return(constant)
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
