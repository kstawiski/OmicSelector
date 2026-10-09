#' Build a frozen single-sample deployment scorer
#'
#' @description
#' Fits one rostered single-sample method and returns a frozen deployment object
#' that can score one incoming specimen without a scored batch, co-resident
#' reference cohort, or batch-correction step. This is a deployability API, not a
#' benchmark-ranking claim: the default \code{"ws-balance-ilr"} is a reasonable
#' compositional default, not a certified winner. It does not assert that any
#' single-sample method out-discriminates a batch-corrected pipeline; validate
#' the method choice on your own data.
#'
#' Reticulate methods (roster tier R2; require a Python backend / venv) are
#' intentionally not fitted by this base-R deployment wrapper.
#'
#' @param X_train Numeric matrix or data frame of training specimens
#'   (samples x features).
#' @param y_train Numeric/integer 0/1 labels aligned to \code{X_train}.
#' @param method Method identifier from \code{\link{singlesample_method_roster}}.
#' @param meta_train Optional training metadata, one row per training specimen.
#' @param annotation Optional feature annotation passed only to methods whose fit
#'   signature includes an \code{annotation} argument.
#' @param verify Logical; if \code{TRUE}, run
#'   \code{\link{singlesample_assert_row_equivariant}} on the training rows (or a
#'   deterministic subset for large inputs) so construction fails loudly when the
#'   fitted scorer is not single-sample deployable.
#' @param nondetect_limit Finite non-negative abundance threshold, frozen at
#'   deployment construction (default 0). For an assay mapping undetected Ct 40
#'   to abundance 1, use 1. This threshold must be specified on the input scale.
#'
#' @return An S3 object of class \code{singlesample_deployable} containing the
#'   frozen \code{model}, \code{method_id}, resolved \code{score_fn},
#'   \code{n_train}, \code{fit_time}, and roster \code{meta}. Its value is
#'   deployability: scoring an incoming specimen from that specimen plus frozen
#'   fit-time parameters only.
#' @seealso \code{\link{score_specimen}},
#'   \code{\link{is_singlesample_deployable}},
#'   \code{\link{singlesample_method_roster}}
#' @param group_train Optional biological-group identifiers, one per training
#'   row, without missing values. Whole groups stay together in inner resampling.
#'   `NULL` preserves the original row-based fit. Invalid or infeasible groups
#'   raise `singlesample_input_error` with code `INVALID_GROUPS`.
#' @export
deploy_singlesample <- function(X_train, y_train,
                                method = "ws-balance-ilr",
                                meta_train = NULL,
                                annotation = NULL,
                                verify = TRUE,
                                nondetect_limit = 0,
                                group_train = NULL) {
  if (!is.null(group_train)) {
    group_train <- .ss_validate_groups(group_train, nrow(X_train), "deploy_singlesample")
  }
  if (!is.numeric(nondetect_limit) || length(nondetect_limit) != 1L ||
      !is.finite(nondetect_limit) || nondetect_limit < 0) {
    .ss_input_reject("INVALID_NONDETECT_LIMIT",
                     "nondetect_limit must be one finite non-negative number",
                     "deploy_singlesample")
  }
  roster <- singlesample_method_roster()
  method <- .singlesample_deploy_match_method(method, roster)
  row <- roster[roster$method_id == method, , drop = FALSE]

  dep_route <- row$dep_route[[1]]
  tier <- row$tier[[1]]
  if (identical(tier, "R2")) {
    stop("deploy_singlesample: method '", method, "' uses dep_route '",
         dep_route, "' (roster tier R2). Reticulate methods (roster tier R2; ",
         "require a Python backend / venv) are out of scope for this base-R ",
         "deployment wrapper.", call. = FALSE)
  }

  X_train <- .singlesample_deploy_as_score_matrix(X_train)
  .ss_guard_input(X_train, colnames(X_train), "deploy_singlesample",
                  nondetect_limit = nondetect_limit)
  if (!is.null(meta_train)) {
    meta_train <- .singlesample_deploy_as_meta(
      meta_train, nrow(X_train), "deploy_singlesample", "meta_train"
    )
  }

  score_fn_name <- row$score_fn[[1]]
  score_fn <- .singlesample_deploy_resolve_fn(score_fn_name, "score_fn",
                                              method)
  fit_fn_name <- row$fit_fn[[1]]
  has_fit <- !is.na(fit_fn_name) && nzchar(fit_fn_name)

  model <- NULL
  fit_time <- 0
  if (has_fit) {
    fit_fn <- .singlesample_deploy_resolve_fn(fit_fn_name, "fit_fn", method)
    timing <- system.time({
      model <- .singlesample_deploy_fit(
        fit_fn = fit_fn,
        fit_fn_name = fit_fn_name,
        method = method,
        X_train = X_train,
        y_train = y_train,
        meta_train = meta_train,
        annotation = annotation,
        group_train = group_train
      )
    })
    fit_time <- unname(timing[["elapsed"]])
  }

  deployable <- structure(
    list(
      model = model,
      method_id = method,
      score_fn = score_fn,
      score_fn_name = score_fn_name,
      fit_fn_name = if (has_fit) fit_fn_name else NA_character_,
      score_route = .singlesample_deploy_score_route(has_fit, score_fn),
      input_guard = list(
        version = 1L,
        required_features = .ss_required_features(model, colnames(X_train)),
        nondetect_limit = nondetect_limit
      ),
      n_train = nrow(X_train),
      fit_time = fit_time,
      meta = list(
        family = row$family[[1]],
        role = row$role[[1]],
        tier = row$tier[[1]],
        dep_route = dep_route
      ),
      verified = FALSE,
      verify_n = 0L
    ),
    class = "singlesample_deployable"
  )

  .singlesample_deploy_preflight_score(deployable, X_train, meta_train)

  if (isTRUE(verify)) {
    probe <- .singlesample_deploy_probe(X_train, meta_train)
    score_fun <- function(model, X, meta = NULL) {
      deployable$model <- model
      .singlesample_deploy_score(deployable, X, meta)
    }
    singlesample_assert_row_equivariant(score_fun, model, probe$X, probe$meta)
    deployable$verified <- TRUE
    deployable$verify_n <- nrow(probe$X)
  }

  deployable
}

#' Score new specimens with a frozen single-sample deployment object
#'
#' @description
#' Scores one or more incoming specimens using only \code{x_new} plus the frozen
#' parameters stored in \code{deployable}. No scored batch, co-resident reference
#' cohort, or test-time batch statistic is required.
#'
#' @param deployable A \code{singlesample_deployable} object from
#'   \code{\link{deploy_singlesample}}.
#' @param x_new A single named numeric vector, or a numeric matrix/data frame
#'   with one or more specimen rows and feature columns.
#' @param meta Optional per-specimen metadata, one row per row of \code{x_new}.
#' @param detected Optional logical matrix aligned to \code{x_new} (or a named
#'   logical vector for one specimen). TRUE means detected before imputation.
#'   Supply this mask when preprocessing replaces nondetects by positive values
#'   that cannot be distinguished from measured values. No detection status is
#'   inferred from scored-batch variation or from a constant specimen profile.
#' @param strict Logical scalar; defaults to \code{FALSE}. If \code{TRUE},
#'   k-TSP and pair-ratio require every fitted-pair feature, trimmed rCLR
#'   requires the frozen training feature universe, and balance ILR requires
#'   every frozen balance member. Other routes and custom adapters reject with
#'   \code{STRICT_UNSUPPORTED}. Existing input guards apply in both modes.
#'
#' @details
#' Invalid input raises a \code{singlesample_input_error} with machine-readable
#' \code{code} and affected \code{rows}. Codes include \code{EMPTY_INPUT},
#' \code{INVALID_FEATURE_NAMES}, \code{ALL_REQUIRED_ABSENT},
#' \code{ALL_NONDETECT}, \code{EMPTY_FITTED_SUPPORT},
#' \code{CONSTANT_FITTED_SUPPORT}, and \code{NO_USABLE_SUPPORT}.
#' A batch containing an invalid row is rejected as a whole. Checks use only
#' each specimen and frozen training information. Positive imputation cannot
#' be recognized without \code{detected} or an appropriate frozen
#' \code{nondetect_limit}. Refit legacy models to record training support status.
#'
#' @return Numeric score vector with one value per incoming specimen. For k-TSP,
#'   pair-ratio, trimmed-rCLR and balance-ILR panels, the scorer's
#'   \code{support} attribute is retained; see \code{\link{score_reo_ktsp}},
#'   \code{\link{score_reo_pairratio}}, \code{\link{score_ws_rclr_panel}}
#'   and \code{\link{score_ws_balance_ilr}}. The default preserves each
#'   scorer's missing-feature policy. Use \code{strict = TRUE} here to require
#'   complete support. Invalid strict values raise \code{INVALID_STRICT};
#'   missing features raise \code{MISSING_REQUIRED_FEATURES}, and absent frozen
#'   universe metadata raises \code{MISSING_FITTED_UNIVERSE}.
#' @export
score_specimen <- function(deployable, x_new, meta = NULL, detected = NULL,
                           strict = FALSE) {
  .singlesample_deploy_assert_object(deployable)
  X <- .singlesample_deploy_as_score_matrix(x_new)
  .singlesample_deploy_score(deployable, X, meta, detected = detected,
                             strict = strict)
}

#' Check whether a deployment scorer is single-sample deployable
#'
#' @description
#' User-facing wrapper over \code{\link{singlesample_is_row_equivariant}}. Supply
#' representative probe rows in \code{X_probe}; the check asks whether each
#' specimen's singleton score equals the same row's score inside a batch, with
#' the model state unchanged during scoring.
#'
#' @param deployable A \code{singlesample_deployable} object from
#'   \code{\link{deploy_singlesample}}.
#' @param X_probe Numeric matrix/data frame of probe rows. If \code{NULL}, the
#'   check cannot be run; callers should supply probe rows from the intended
#'   panel or validation surface.
#'
#' @return \code{TRUE} when the probe check passes, otherwise \code{FALSE}.
#' @export
is_singlesample_deployable <- function(deployable, X_probe = NULL) {
  .singlesample_deploy_assert_object(deployable)
  if (is.null(X_probe)) {
    stop("is_singlesample_deployable: X_probe is required; supply probe rows ",
         "to test singleton scoring against batch scoring.", call. = FALSE)
  }
  X_probe <- .singlesample_deploy_as_score_matrix(X_probe)
  score_fun <- function(model, X, meta = NULL) {
    deployable$model <- model
    .singlesample_deploy_score(deployable, X, meta)
  }
  singlesample_is_row_equivariant(score_fun, deployable$model, X_probe)
}

#' Print a single-sample deployment object
#'
#' @param x A \code{singlesample_deployable} object.
#' @param ... Unused.
#'
#' @return Invisibly, \code{x}.
#' @export
print.singlesample_deployable <- function(x, ...) {
  .singlesample_deploy_assert_object(x)
  cat("OmicSelector single-sample deployable\n")
  cat("  method_id: ", x$method_id, "\n", sep = "")
  cat("  n_train: ", x$n_train, "\n", sep = "")
  if (isTRUE(x$verified)) {
    cat("  gate: verified row-equivariant scoring on ", x$verify_n,
        " probe rows (single specimen == same row in batch; frozen model)\n",
        sep = "")
  } else {
    cat("  gate: not run for this object; call is_singlesample_deployable() ",
        "with probe rows\n", sep = "")
  }
  cat("  caveat: deployability only; not a claim of superior discrimination ",
      "over a batch-corrected pipeline.\n", sep = "")
  invisible(x)
}

.singlesample_deploy_match_method <- function(method, roster) {
  available <- roster$method_id
  if (!is.character(method) || length(method) != 1L ||
      is.na(method) || !nzchar(method) || !method %in% available) {
    stop("deploy_singlesample: unknown method '", paste(method, collapse = ", "),
         "'. Available methods: ", paste(available, collapse = ", "),
         call. = FALSE)
  }
  method
}

.singlesample_deploy_resolve_fn <- function(fn_name, field, method) {
  if (is.na(fn_name) || !nzchar(fn_name)) {
    stop("deploy_singlesample: method '", method, "' requires a declared ",
         field, ", but none is listed in singlesample_method_roster().",
         call. = FALSE)
  }
  ns <- asNamespace("OmicSelector")
  if (!exists(fn_name, envir = ns, mode = "function", inherits = FALSE)) {
    stop("deploy_singlesample: method '", method, "' requires ", field,
         " '", fn_name, "', but it is not available in the OmicSelector ",
         "namespace.", call. = FALSE)
  }
  get(fn_name, envir = ns, mode = "function", inherits = FALSE)
}

.singlesample_deploy_fit <- function(fit_fn, fit_fn_name, method, X_train,
                                     y_train, meta_train, annotation,
                                     group_train = NULL, hp = NULL) {
  fml <- names(formals(fit_fn))
  data_arg <- intersect(c("X_train", "x_train", "X", "x"), fml)
  if (length(data_arg) == 0L) {
    stop("deploy_singlesample: method '", method, "' requires fit_fn '",
         fit_fn_name, "' to declare one of X_train, x_train, X, or x.",
         call. = FALSE)
  }

  args <- list()
  args[[data_arg[[1L]]]] <- X_train
  label_arg <- intersect(c("y_train", "y"), fml)
  if (length(label_arg) > 0L) args[[label_arg[[1L]]]] <- y_train
  if ("meta_train" %in% fml) args$meta_train <- meta_train
  if ("annotation" %in% fml) args$annotation <- annotation
  if ("meta" %in% fml) args$meta <- meta_train
  if (!is.null(hp) && "hp" %in% fml) args$hp <- hp
  if (!is.null(group_train) && "group_train" %in% fml) {
    args$group_train <- .ss_validate_groups(group_train, nrow(X_train),
                                           "deploy_singlesample")
  }

  tryCatch(
    do.call(fit_fn, args),
    error = function(e) {
      if (inherits(e, "singlesample_input_error")) stop(e)
      detail <- .singlesample_deploy_safe_error_detail(conditionMessage(e))
      msg <- paste0("deploy_singlesample: method '", method, "' requires ",
                    fit_fn_name, "'s backend which is not available or failed ",
                    "during fit.")
      if (nzchar(detail)) msg <- paste0(msg, " Original error: ", detail)
      stop(msg, call. = FALSE)
    }
  )
}

.singlesample_deploy_score_route <- function(has_fit, score_fn) {
  if (isTRUE(has_fit)) return("score_call")
  fml <- names(formals(score_fn))
  if (length(fml) >= 1L && identical(fml[[1]], "model")) {
    return("direct_model")
  }
  "direct_x"
}

.singlesample_deploy_score <- function(deployable, X, meta = NULL,
                                       caller = "score_specimen", detected = NULL,
                                       strict = FALSE) {
  .ss_check_strict(strict, caller)
  if (strict) {
    .ss_strict_scorer(deployable$method_id, caller)
    if (!identical(deployable$score_route, "score_call")) {
      .ss_input_reject("STRICT_UNSUPPORTED",
                       "strict checking is unavailable for this score route", caller)
    }
  }
  guard <- deployable$input_guard
  required <- if (is.null(guard)) .ss_required_features(deployable$model) else guard$required_features
  limit <- if (is.null(guard)) 0 else guard$nondetect_limit
  .ss_guard_input(X, required, caller, detected, limit)
  if (is.list(deployable$model)) .ss_guard_support(deployable$model, caller)
  out <- .singlesample_deploy_score_raw(deployable, X, meta, strict = strict)
  .singlesample_deploy_validate_scores(out, nrow(X), deployable$method_id,
                                       caller)
}

.singlesample_deploy_score_raw <- function(deployable, X, meta = NULL,
                                           strict = FALSE) {
  if (!is.null(meta)) {
    meta <- .singlesample_deploy_as_meta(meta, nrow(X), "score_specimen",
                                         "meta")
  }

  switch(
    deployable$score_route,
    score_call = singlesample_score_call(deployable$method_id,
                                         deployable$model, X, meta, strict = strict),
    direct_model = deployable$score_fn(deployable$model, X, meta),
    direct_x = deployable$score_fn(X),
    stop("score_specimen: unknown score route '", deployable$score_route, "'.",
         call. = FALSE)
  )
}

.singlesample_deploy_as_score_matrix <- function(x_new) {
  if (is.numeric(x_new) && is.null(dim(x_new))) {
    if (!length(x_new)) {
      .ss_input_reject("EMPTY_INPUT", "x_new must have at least one feature", "score_specimen")
    }
    if (is.null(names(x_new)) || anyNA(names(x_new)) ||
        any(!nzchar(names(x_new))) || anyDuplicated(names(x_new))) {
      .ss_input_reject("INVALID_FEATURE_NAMES",
                       "a numeric vector x_new must have unique, non-empty feature names",
                       "score_specimen")
    }
    X <- matrix(as.numeric(x_new), nrow = 1L,
                dimnames = list(NULL, names(x_new)))
    return(X)
  }
  if (is.data.frame(x_new)) x_new <- as.matrix(x_new)
  if (!is.matrix(x_new) || !is.numeric(x_new)) {
    stop("score_specimen: x_new must be a named numeric vector, numeric ",
         "matrix, or numeric data frame.", call. = FALSE)
  }
  if (nrow(x_new) < 1L || ncol(x_new) < 1L) {
    .ss_input_reject("EMPTY_INPUT", "x_new must have at least one row and one feature",
                     "score_specimen")
  }
  storage.mode(x_new) <- "double"
  x_new
}

.singlesample_deploy_as_meta <- function(meta, n, caller, arg) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (nrow(meta) != n) {
    stop(caller, ": ", arg, " must have one row per row of ",
         if (identical(caller, "deploy_singlesample")) "X_train" else "x_new",
         ".", call. = FALSE)
  }
  meta
}

.singlesample_deploy_probe <- function(X, meta = NULL, max_rows = 50L) {
  X <- .singlesample_deploy_as_score_matrix(X)
  if (!is.null(meta)) {
    meta <- .singlesample_deploy_as_meta(meta, nrow(X), "deploy_singlesample",
                                         "meta_train")
  }
  if (nrow(X) > max_rows) {
    idx <- unique(as.integer(round(seq(1L, nrow(X), length.out = max_rows))))
    X <- X[idx, , drop = FALSE]
    if (!is.null(meta)) {
      meta <- as.data.frame(meta, stringsAsFactors = FALSE)[idx, , drop = FALSE]
    }
  }
  list(X = X, meta = meta)
}

.singlesample_deploy_preflight_score <- function(deployable, X_train,
                                                 meta_train) {
  probe <- .singlesample_deploy_probe(X_train, meta_train, max_rows = 5L)
  tryCatch(
    .singlesample_deploy_score(deployable, probe$X, probe$meta,
                               caller = "deploy_singlesample"),
    error = function(e) {
      if (inherits(e, "singlesample_input_error")) stop(e)
      msg <- conditionMessage(e)
      if (startsWith(msg, "deploy_singlesample:")) {
        stop(msg, call. = FALSE)
      }
      stop("deploy_singlesample: method '", deployable$method_id,
           "' is not directly deployable via deploy_singlesample(); its ",
           "scorer failed on a training probe. See ",
           "singlesample_method_roster().", call. = FALSE)
    }
  )
  invisible(TRUE)
}

.singlesample_deploy_validate_scores <- function(out, n, method_id,
                                                caller = "score_specimen") {
  d <- dim(out)
  if (!is.null(d) && !(length(d) == 2L && d[2] == 1L)) {
    .singlesample_deploy_stop_bad_scores(out, method_id, caller)
  }
  if (!is.numeric(out)) {
    .singlesample_deploy_stop_bad_scores(out, method_id, caller)
  }
  support <- attr(out, "support", exact = TRUE)
  out <- as.numeric(out)
  if (!is.null(support)) attr(out, "support") <- support
  if (length(out) != n) {
    .singlesample_deploy_stop_bad_scores(out, method_id, caller)
  }
  if (any(!is.finite(out))) {
    .singlesample_deploy_stop_bad_scores(out, method_id, caller,
                                         nonfinite = TRUE)
  }
  out
}

.singlesample_deploy_stop_bad_scores <- function(out, method_id, caller,
                                                 nonfinite = FALSE) {
  if (identical(caller, "deploy_singlesample")) {
    stop("deploy_singlesample: method '", method_id,
         "' does not return one score per specimen (its scorer produced ",
         .singlesample_deploy_score_shape(out, nonfinite), "); it is not ",
         "directly deployable via deploy_singlesample() -- see ",
         "singlesample_method_roster().", call. = FALSE)
  }

  if (nonfinite) {
    stop("score_specimen: method '", method_id,
         "' returned non-finite score(s); expected one finite numeric score ",
         "per specimen.", call. = FALSE)
  }
  stop("score_specimen: method '", method_id, "' returned ",
       .singlesample_deploy_score_shape(out),
       "; expected one finite numeric score per specimen.", call. = FALSE)
}

.singlesample_deploy_score_shape <- function(out, nonfinite = FALSE) {
  d <- dim(out)
  if (!is.null(d)) return(paste(d, collapse = "x"))
  shape <- paste0(class(out)[[1L]], " length ", length(out))
  if (nonfinite) shape <- paste0(shape, " with non-finite values")
  shape
}

.singlesample_deploy_safe_error_detail <- function(msg) {
  if (!is.character(msg) || length(msg) == 0L || is.na(msg[[1L]])) {
    return("")
  }
  msg <- msg[[1L]]
  cryptic <- c("unused argument", "could not find function", "array with dim")
  if (any(vapply(cryptic, grepl, logical(1L), x = msg, fixed = TRUE))) {
    return("")
  }
  msg
}

.singlesample_deploy_assert_object <- function(x) {
  if (!inherits(x, "singlesample_deployable")) {
    stop("Expected a singlesample_deployable object.", call. = FALSE)
  }
  invisible(TRUE)
}
