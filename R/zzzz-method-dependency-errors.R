# Compose the startup hook without changing any v2.6.9 computation source file.
# Successful calls still evaluate the original closure with the original arguments.
.method_dependency_error <- function(e, method) {
  message <- conditionMessage(e)
  dependency <- grepl("package|module|backend|import|checkpoint|weights|reticulate|Python", message, ignore.case = TRUE) &&
    grepl("required|unavailable|not available|not installed|not found|No module|missing|failed|cannot|not importable", message, ignore.case = TRUE)
  if (dependency && !grepl("install_method_dependencies", message, fixed = TRUE)) {
    e$message <- paste0(message, "\nRun install_method_dependencies(methods = \"", method, "\").")
  }
  stop(e)
}
.method_dependency_wrapper <- function(original, method = NULL) {
  force(original); force(method)
  wrapper <- function() NULL
  formals(wrapper) <- formals(original)
  body(wrapper) <- quote({
    supplied <- names(as.list(match.call(expand.dots = FALSE)))[-1L]
    args <- lapply(supplied, as.name)
    names(args) <- ifelse(supplied == "...", "", supplied)
    call <- as.call(c(list(quote(original)), args))
    id <- method
    tryCatch(eval(call, envir = environment()),
             error = function(e) .method_dependency_error(e, id))
  })
  environment(wrapper) <- environment()
  wrapper
}
.method_install_error_hints <- function(ns) {
  roster <- singlesample_method_roster()
  req <- method_requirements()
  roster <- roster[roster$method_id %in% req$method_id[nzchar(req$r_packages)], , drop = FALSE]
  for (i in seq_len(nrow(roster))) for (nm in c(roster$fit_fn[i], roster$score_fn[i])) {
    if (is.na(nm) || !exists(nm, envir = ns, inherits = FALSE)) next
    original <- get(nm, envir = ns, inherits = FALSE)
    assign(nm, .method_dependency_wrapper(original, roster$method_id[i]), envir = ns)
  }
  original <- get("deploy_singlesample", envir = ns, inherits = FALSE)
  assign("deploy_singlesample", .method_dependency_wrapper(original), envir = ns)
  invisible(NULL)
}
.onLoad <- local({
  previous <- .onLoad
  function(libname, pkgname) {
    previous(libname, pkgname)
    .method_install_error_hints(asNamespace(pkgname))
    invisible(NULL)
  }
})
