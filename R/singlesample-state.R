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

.ss_python_clone <- function(x, mod, freeze_state = TRUE) .ss_with_r_rng({
  if (inherits(x, "python.builtin.object")) return(mod$clone_for_inference(x, freeze_state = freeze_state))
  if (is.list(x)) {
    out <- lapply(x, .ss_python_clone, mod = mod, freeze_state = freeze_state)
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
    singlesample_score_call(deployable$method_id,
      if (foundation) deployable$model else .ss_python_clone(deployable$model, mod),
      X[i, , drop = FALSE], .singlesample_subset_meta(meta, i), strict = strict)[[1L]]
  }, numeric(1L))
})

.ss_python_snapshot <- function(x) {
  if (inherits(x, "python.builtin.object")) return(list(
    python_state_sha256 = reticulate::py_to_r(.ss_python_state_module()$state_fingerprint(x))))
  # Bytecode can change through R's JIT without changing fitted state. Keep
  # code and closure state; deparse avoids unbound formal-argument sentinels.
  if (is.function(x)) return(list(formals = deparse(formals(x)), body = deparse(body(x)),
                                 environment = environment(x)))
  if (is.list(x)) {
    out <- lapply(x, .ss_python_snapshot)
    attributes(out) <- attributes(x)
    return(out)
  }
  x
}

#' Fingerprint frozen fitted state, including Python parameters and buffers
#'
#' Computes SHA-256 of R state and deterministic content snapshots of nested
#' Python objects. Python snapshots include estimator preprocessing, training
#' context, neural parameters, buffers and evaluation flags. Unknown opaque
#' state is rejected, rather than hashed as a pointer placeholder.
#' @param model A fitted model or single-sample deployment object.
#' @return One SHA-256 string. Equality tests detect state changes, not validity
#'   of the fitted model or clinical performance.
#' @export
singlesample_state_fingerprint <- function(model) {
  digest::digest(.singlesample_default_model_digest(model), algo = "sha256", serialize = FALSE)
}
