.ss_foundation_module <- function() {
  reticulate::import_from_path("foundation_inference",
    path = system.file("python", package = "OmicSelector"), convert = FALSE)
}

.ss_python_state_module <- function() {
  reticulate::import_from_path("singlesample_state",
    path = system.file("python", package = "OmicSelector"), convert = FALSE)
}

.ss_with_r_rng <- function(expr) {
  # Reticulate calls can create .Random.seed even without drawing R randomness.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  state <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", state, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  force(expr)
}

.ss_python_backend_check <- function(dep_route, importer = function(module) {
  reticulate::import(module, delay_load = FALSE)
}) {
  modules <- switch(dep_route,
    `reticulate-tabpfn` = c("torch", "tabpfn"),
    `reticulate-tabicl` = c("torch", "tabicl"),
    `reticulate-tabdpt` = c("torch", "tabdpt"),
    `reticulate-kymatio` = c("torch", "kymatio"),
    `reticulate-pyod` = "pyod",
    `reticulate-torch-torchvision` = c("torch", "torchvision"), "torch")
  if (!requireNamespace("reticulate", quietly = TRUE)) .ss_input_reject(
    "PYTHON_BACKEND_UNAVAILABLE", "requires reticulate and an importable Python backend", "deploy_singlesample")
  for (module in modules) tryCatch(importer(module),
    error = function(e) .ss_input_reject("PYTHON_BACKEND_UNAVAILABLE",
      paste0("requires Python module '", module, "': ", conditionMessage(e)), "deploy_singlesample"))
  invisible(TRUE)
}

.ss_python_objects <- function(x) {
  if (inherits(x, "python.builtin.object")) return(list(x))
  if (is.list(x)) return(unlist(lapply(x, .ss_python_objects), recursive = FALSE))
  list()
}

.ss_python_freeze <- function(model) {
  objects <- .ss_python_objects(model)
  if (length(objects)) {
    mod <- .ss_python_state_module()
    for (obj in objects) mod$freeze(obj)
  }
  invisible(model)
}

.ss_python_copy_scope <- new.env(parent = emptyenv())

.ss_with_python_copies <- function(fn) {
  mod <- .ss_python_state_module()
  scope <- .ss_python_copy_scope$current
  outermost <- is.null(scope)
  if (outermost) {
    scope <- mod$CopyTracker()
    .ss_python_copy_scope$current <- scope
  }
  complete <- FALSE
  if (outermost) on.exit({
    .ss_python_copy_scope$current <- NULL
    # Collect nested copies together. Collecting an inner scorer first would
    # promote the still-live deployment copy and force a full collection.
    gc(verbose = FALSE, full = FALSE)
    if (!complete || reticulate::py_to_r(scope$pending())) {
      gc(verbose = FALSE)
      if (!complete || reticulate::py_to_r(scope$pending())) scope$collect_cycles()
    }
  }, add = TRUE)
  result <- tryCatch(fn(scope), error = function(e) {
    # Keep traceback diagnostics, but discard unwound locals owning copies.
    if (inherits(e, "python.builtin.BaseException")) mod$clear_exception_frames(e)
    stop(e)
  })
  complete <- TRUE
  result
}

.ss_python_clone <- function(x, mod, freeze_state = TRUE, scope = NULL) .ss_with_r_rng({
  if (inherits(x, "python.builtin.object")) {
    cloned <- mod$clone_for_inference(x, freeze_state = freeze_state)
    if (!is.null(scope)) scope$track(cloned)
    return(cloned)
  }
  if (is.list(x)) {
    out <- lapply(x, .ss_python_clone, mod = mod, freeze_state = freeze_state, scope = scope)
    attributes(out) <- attributes(x)
    return(out)
  }
  x
})

.ss_python_score <- function(deployable, X, meta, strict) .ss_with_r_rng({
  # A fixed seed and one forward pass per specimen give the same inference
  # path for singleton and batch inputs, including in-context classifiers.
  mod <- .ss_python_state_module()
  # These scorers create one private copy themselves. Avoid a second complete
  # estimator copy in the deployment wrapper, retaining all input/RNG guards.
  foundation <- inherits(deployable$model, c("tabpfn_model", "tabdpt_model"))
  vapply(seq_len(nrow(X)), function(i) {
    context <- mod$inference(as.integer(deployable$model$hp$seed %||% 1L))
    context$`__enter__`()
    on.exit({
      context$`__exit__`(NULL, NULL, NULL)
      # Reticulate owns the scorer's private estimator through small R external
      # pointers. R's allocation threshold does not see their Python tensor
      # storage, so successive specimens can retain whole training caches until
      # an unrelated R collection. The scorer has returned (also on errors):
      # finalize those unreachable wrappers before accepting another specimen.
      if (foundation) gc(verbose = FALSE)
    }, add = TRUE)
    score <- function(model) singlesample_score_call(deployable$method_id, model,
      X[i, , drop = FALSE], .singlesample_subset_meta(meta, i), strict = strict)[[1L]]
    if (foundation) score(deployable$model) else .ss_with_python_copies(function(scope) {
      score(.ss_python_clone(deployable$model, mod, scope = scope))
    })
  }, numeric(1L))
})

.ss_python_snapshot <- function(x) {
  environments <- list()
  walk <- function(value) {
    if (inherits(value, "python.builtin.object")) return(list(
      python_state_sha256 = reticulate::py_to_r(.ss_python_state_module()$state_fingerprint(value))))
    # Snapshot closure bindings recursively: dependency wrappers retain the
    # original function in a local environment. Serializing that function raw
    # would include session-specific JIT bytecode and source references.
    if (is.function(value)) return(list(formals = deparse(formals(value)),
      body = deparse(body(value)), environment = walk(environment(value))))
    if (is.environment(value)) {
      if (isNamespace(value)) return(list(namespace = getNamespaceName(value)))
      if (identical(value, globalenv())) return(list(environment = ".GlobalEnv"))
      if (identical(value, baseenv())) return(list(environment = "base"))
      if (identical(value, emptyenv())) return(list(environment = "empty"))
      for (i in seq_along(environments)) {
        if (identical(value, environments[[i]])) return(list(environment_reference = i))
      }
      index <- length(environments) + 1L
      environments[[index]] <<- value
      keys <- sort(ls(value, all.names = TRUE))
      if (any(vapply(keys, bindingIsActive, logical(1L), env = value)))
        stop("Cannot fingerprint an environment with active bindings.", call. = FALSE)
      bindings <- lapply(keys, function(key) walk(get(key, envir = value, inherits = FALSE)))
      names(bindings) <- keys
      return(list(environment_index = index, bindings = bindings,
                  parent = walk(parent.env(value))))
    }
    if (is.list(value)) {
      out <- lapply(value, walk)
      attributes(out) <- attributes(value)
      return(out)
    }
    value
  }
  walk(x)
}

#' Fingerprint frozen fitted state, including Python parameters and buffers
#'
#' Computes SHA-256 of R state and deterministic content snapshots of nested
#' Python objects. Python snapshots include estimator preprocessing, training
#' context, neural parameters, buffers and evaluation flags. Snapshots use
#' numerical content rather than shared-object identities or pickle bytes.
#' Autograd ownership, TabICL checkpoint paths and temporary inference workspaces
#' are excluded; frozen context caches and estimator inference settings remain
#' included. Fingerprints from earlier package versions must be regenerated.
#' Unknown opaque
#' state is rejected, rather than hashed as a pointer placeholder.
#' @param model A fitted model or single-sample deployment object.
#' @return One SHA-256 string. Equality tests detect state changes, not validity
#'   of the fitted model or clinical performance.
#' @export
singlesample_state_fingerprint <- function(model) {
  digest::digest(.singlesample_default_model_digest(model), algo = "sha256", serialize = FALSE)
}
