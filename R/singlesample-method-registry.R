# Session-local callbacks. Built-in routes are never stored or replaced here.
.singlesample_method_registry <- new.env(parent = emptyenv())

#' Register a custom single-specimen scorer for this R session
#'
#' Adds a fit/score pair to [singlesample_method_roster()],
#' [deploy_singlesample()] and [singlesample_score_call()]. Registration is
#' local to the loaded package in this R session. Built-in routes cannot be
#' replaced. Deployment objects retain their resolved custom scorer and strict
#' support flag after a route is replaced or removed.
#'
#' @param method_id One syntactically valid R identifier, unused by built-in
#'   routes. Existing user routes require `overwrite = TRUE`.
#' @param fit_fn Function declaring `X_train` and `y_train`. Additional arguments
#'   must have defaults, except optional deployment inputs `meta_train`, `meta`
#'   and `annotation`. `...` is allowed but not required.
#' @param score_fn Function declaring `model`, `X` and `meta = NULL`, returning
#'   one finite numeric score per specimen. Additional arguments must have
#'   defaults; `...` is allowed. All fitted parameters must be stored in `model`.
#'   Scoring must not change the model or depend on other specimens or mutable
#'   external state.
#' @param family,role,estimand Non-empty character scalars describing the route.
#'   Use `estimand = "within"` to include a route in within-cohort benchmark
#'   helpers, or `"transfer"` for transfer routing. The descriptive default
#'   is excluded from within-cohort comparisons. Selector eligibility also
#'   requires role `"baseline"` or `"discriminator"`.
#' @param strict_support Logical scalar. If `TRUE`, `score_fn` must accept
#'   `strict` explicitly or through `...`; scoring passes this flag in both
#'   modes. The callback defines its strict feature policy. Otherwise strict
#'   requests fail with `STRICT_UNSUPPORTED`.
#' @param audit Logical scalar. If `TRUE`, fit the probe and run
#'   [singlesample_assert_row_equivariant()], including its serialized fitted-
#'   state immutability check. Both scoring modes are audited when strict
#'   support is declared. `FALSE` skips only registration auditing;
#'   deployment verification still defaults to `TRUE`.
#' @param probe Optional list with elements `X` (numeric matrix with
#'   unique feature names and at least four rows) and `y` (aligned binary 0/1
#'   labels, both classes present), named or supplied in that order. The default
#'   has 12 positive-valued specimens
#'   and six features named `feature1` through `feature6`. Supply a probe for
#'   callbacks requiring a specific panel. Extra fit settings can be captured
#'   in a function wrapper with fixed configuration.
#' @param overwrite Logical scalar. Replace an existing user route only after
#'   all requested checks pass.
#'
#' @details
#' Registration failures are `singlesample_input_error` conditions with codes
#' `REGISTRATION_INVALID_ID`, `REGISTRATION_ID_EXISTS`,
#' `REGISTRATION_INVALID_FUNCTION`, `REGISTRATION_INVALID_ARGUMENT`, or
#' `REGISTRATION_AUDIT_FAILED`. Audit errors retain the original message.
#' Custom scoring uses the deployment input guards. Callback errors use
#' `CUSTOM_SCORER_FAILED`, unless already coded; invalid or non-finite score
#' output uses `CUSTOM_INVALID_SCORES`.
#'
#' Auditing tests the supplied probe with absolute tolerance 1e-8. It is not a
#' proof for every possible input or evidence of clinical performance. Models
#' with external pointers cannot pass the default serialization audit. Audit
#' fitting preserves the caller's random seed, but callbacks must avoid other
#' external side effects. Re-register routes in each new R session to fit new
#' deployments; saved deployments contain their scorer but still require its
#' dependencies and a self-contained callback environment.
#'
#' @return `register_singlesample_method()` returns invisibly `TRUE`.
#'   `unregister_singlesample_method()` returns invisibly whether a user route
#'   was removed. `registered_singlesample_methods()` returns a data frame in
#'   identifier order with roster columns plus `strict_support` and `audited`.
#' @seealso [deploy_singlesample()], [score_specimen()]
#' @export
register_singlesample_method <- function(
    method_id, fit_fn, score_fn, family = "custom", role = "discriminator",
    estimand = "custom single-specimen score", strict_support = FALSE,
    audit = TRUE, probe = NULL, overwrite = FALSE) {
  caller <- "register_singlesample_method"
  .ss_registration_id(method_id, caller)
  for (arg in c("strict_support", "audit", "overwrite")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      .ss_input_reject("REGISTRATION_INVALID_ARGUMENT",
                       paste0(arg, " must be TRUE or FALSE"), caller)
    }
  }
  if (!is.null(.ss_registered_method(method_id)) && !overwrite) {
    .ss_input_reject("REGISTRATION_ID_EXISTS", "user route already exists", caller)
  }
  for (arg in c("family", "role", "estimand")) {
    value <- get(arg)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(value))) {
      .ss_input_reject("REGISTRATION_INVALID_ARGUMENT",
                       paste0(arg, " must be one non-empty string"), caller)
    }
  }
  .ss_registration_function(fit_fn, "fit_fn", c("X_train", "y_train"),
                             c("meta_train", "meta", "annotation"))
  .ss_registration_function(score_fn, "score_fn", c("model", "X", "meta"),
                             if (strict_support) "strict" else character())
  if (!identical(formals(score_fn)$meta, NULL) ||
      (strict_support && !any(c("strict", "...") %in% names(formals(score_fn))))) {
    .ss_input_reject("REGISTRATION_INVALID_FUNCTION",
                     "score_fn must declare meta = NULL and accept strict when strict_support = TRUE",
                     caller)
  }
  row <- data.frame(
    method_id = method_id, family = family, estimand = estimand, role = role,
    tier = "R1", dep_route = "session_registry",
    fit_fn = paste0("user::", method_id, "::fit"),
    score_fn = paste0("user::", method_id, "::score"),
    pkg_status = "user", notes = "Session-registered custom scorer",
    row_source = "user", lopbo_mechanism = if (estimand == "within")
      "within_cohort" else if (estimand == "transfer") "transfer_fold" else "na",
    single_sample_gate = "required",
    is_discriminator = role %in% c("discriminator", "discriminator-conditional"),
    is_negative_control = role == "negative_control", stringsAsFactors = FALSE
  )
  entry <- list(fit_fn = fit_fn, score_fn = score_fn, row = row,
                strict_support = strict_support, audited = audit)
  if (audit) {
    old_seed <- .singlesample_get_seed()
    on.exit(.singlesample_restore_seed(old_seed), add = TRUE)
    tryCatch({
      set.seed(1L)
      if (is.null(probe)) probe <- .ss_registration_probe()
      if (is.list(probe) && length(probe) == 2L && is.null(names(probe))) {
        names(probe) <- c("X", "y")
      }
      if (!is.list(probe) || !all(c("X", "y") %in% names(probe))) {
        stop("probe must be a list with named X and y elements", call. = FALSE)
      }
      X <- if (is.data.frame(probe$X)) as.matrix(probe$X) else probe$X
      .ss_guard_input(X, colnames(X), caller)
      y <- probe$y
      if (nrow(X) < 4L || !is.numeric(y) || length(y) != nrow(X) ||
          anyNA(y) || !setequal(y, c(0, 1))) {
        stop("probe requires at least four rows and aligned binary 0/1 labels with both classes",
             call. = FALSE)
      }
      model <- .singlesample_deploy_fit(fit_fn, row$fit_fn, method_id,
                                        X, y, NULL, NULL)
      for (strict in if (strict_support) c(FALSE, TRUE) else FALSE) {
        score_fun <- function(model, X, meta = NULL) {
          .ss_registered_score(entry, model, X, meta, strict, method_id)
        }
        singlesample_assert_row_equivariant(score_fun, model, X,
                                            n_singletons = nrow(X))
      }
    }, error = function(e) {
      .ss_input_reject("REGISTRATION_AUDIT_FAILED", conditionMessage(e), caller)
    })
  }
  assign(method_id, entry, envir = .singlesample_method_registry)
  invisible(TRUE)
}

#' @rdname register_singlesample_method
#' @export
unregister_singlesample_method <- function(method_id) {
  .ss_registration_id(method_id, "unregister_singlesample_method")
  present <- !is.null(.ss_registered_method(method_id))
  if (present) rm(list = method_id, envir = .singlesample_method_registry)
  invisible(present)
}

#' @rdname register_singlesample_method
#' @export
registered_singlesample_methods <- function() {
  ids <- ls(.singlesample_method_registry, all.names = TRUE)
  entries <- lapply(ids, .ss_registered_method)
  out <- .singlesample_builtin_method_roster()[FALSE, , drop = FALSE]
  if (length(entries)) out <- do.call(rbind, lapply(entries, `[[`, "row"))
  out$strict_support <- vapply(entries, `[[`, logical(1), "strict_support")
  out$audited <- vapply(entries, `[[`, logical(1), "audited")
  rownames(out) <- NULL
  out
}

.ss_registration_id <- function(method_id, caller) {
  if (!is.character(method_id) || length(method_id) != 1L ||
      is.na(method_id) || !nzchar(method_id)) {
    .ss_input_reject("REGISTRATION_INVALID_ID", "method_id must be one R identifier", caller)
  }
  if (method_id %in% .singlesample_builtin_method_roster()$method_id) {
    .ss_input_reject("REGISTRATION_ID_EXISTS", "built-in routes cannot be replaced or removed", caller)
  }
  if (!identical(make.names(method_id), unname(method_id))) {
    .ss_input_reject("REGISTRATION_INVALID_ID", "method_id must be a syntactically valid R identifier", caller)
  }
  invisible(TRUE)
}

.ss_registration_function <- function(fn, arg, required, allowed = character()) {
  fml <- if (is.function(fn)) formals(fn) else NULL
  missing_default <- vapply(fml, identical, logical(1), quote(expr = ))
  extra_required <- setdiff(names(fml)[missing_default],
                            c(required, allowed, "..."))
  if (is.null(fml) || !all(required %in% names(fml)) || length(extra_required)) {
    .ss_input_reject("REGISTRATION_INVALID_FUNCTION",
                     paste0(arg, " must declare ", paste(required, collapse = ", "),
                            " and give defaults to additional arguments"),
                     "register_singlesample_method")
  }
  invisible(TRUE)
}

.ss_registered_method <- function(method_id) {
  if (!is.character(method_id) || length(method_id) != 1L ||
      is.na(method_id) || !nzchar(method_id)) return(NULL)
  if (!exists(method_id, envir = .singlesample_method_registry, inherits = FALSE)) return(NULL)
  get(method_id, envir = .singlesample_method_registry, inherits = FALSE)
}

.ss_registration_probe <- function() {
  X <- matrix(seq_len(72L) / 10 + 1, nrow = 12L,
               dimnames = list(NULL, paste0("feature", seq_len(6L))))
  y <- rep(c(0, 1), each = 6L)
  X[y == 1, 1:2] <- X[y == 1, 1:2] + 4
  list(X = X, y = y)
}

.ss_registered_strict <- function(entry, strict, caller) {
  if (strict && !isTRUE(entry$strict_support)) {
    .ss_input_reject("STRICT_UNSUPPORTED",
                     "strict checking is unavailable for this custom method", caller)
  }
  invisible(TRUE)
}

.ss_registered_score <- function(entry, model, X, meta, strict, method_id,
                                  caller = "singlesample_score_call") {
  .ss_registered_strict(entry, strict, caller)
  if (is.data.frame(X)) X <- as.matrix(X)
  .ss_guard_input(X, .ss_required_features(model), caller)
  if (is.list(model)) .ss_guard_support(model, caller)
  if (!is.null(meta)) {
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (nrow(meta) != nrow(X)) {
      .ss_input_reject("INVALID_INPUT", "meta must have one row per specimen", caller)
    }
  }
  args <- list(model = model, X = X, meta = meta)
  if (entry$strict_support) args$strict <- strict
  out <- tryCatch(do.call(entry$score_fn, args), error = function(e) {
    if (inherits(e, "singlesample_input_error")) stop(e)
    .ss_input_reject("CUSTOM_SCORER_FAILED", conditionMessage(e), caller)
  })
  tryCatch(.singlesample_deploy_validate_scores(out, nrow(X), method_id, caller),
           error = function(e) {
             .ss_input_reject("CUSTOM_INVALID_SCORES", conditionMessage(e), caller)
           })
}
