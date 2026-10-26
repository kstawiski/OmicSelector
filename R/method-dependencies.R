# Dependency metadata is separate from numerical implementation files.
.method_split <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  unique(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}

.method_select <- function(x, methods = NULL, families = NULL) {
  for (arg in c("methods", "families")) {
    value <- get(arg)
    if (!is.null(value) && (!is.character(value) || anyNA(value)))
      stop(arg, " must be NULL or a character vector without NA.", call. = FALSE)
    column <- if (arg == "methods") "method_id" else "family"
    bad <- setdiff(value, x[[column]])
    if (length(bad)) stop("Unknown ", arg, ": ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (is.null(methods) && is.null(families)) return(x)
  x[x$method_id %in% methods | x$family %in% families, , drop = FALSE]
}

#' Requirements and measured costs for registered methods
#' @param methods Optional method identifiers. NULL selects all unless families is set.
#' @param families Optional roster family letters. Methods and families select their union.
#' @return A data frame, one row per method. Semicolon-separated dependency fields
#'   distinguish required and optional packages. Cost columns are medians over
#'   completed CPU benchmark cells; full cells, hardware and source are attached
#'   as the `measurements` attribute. Missing measurements are NA, never zero.
#' @export
method_requirements <- function(methods = NULL, families = NULL) {
  roster <- singlesample_method_roster()
  path <- system.file("extdata", "method_requirements.tsv", package = "OmicSelector")
  meta <- utils::read.delim(path, stringsAsFactors = FALSE, na.strings = "NA")
  x <- cbind(roster[, c("method_id", "family", "tier", "dep_route")],
             meta[match(roster$method_id, meta$method_id), -1L, drop = FALSE])
  rownames(x) <- NULL
  x$requirements_known <- !is.na(x$requirements_source)
  # Registered callbacks have no dependency declaration in the existing API.
  x$requirements_source[!x$requirements_known] <- "user callback; requirements undeclared"
  for (nm in c("fit_seconds", "score_seconds_per_specimen", "peak_memory_bytes",
               "frozen_model_bytes")) x[[nm]] <- NA_real_
  x$measuring_node <- NA_character_
  x$cost_status <- "PB1 measurements pending"
  path <- system.file("extdata", "method_costs.tsv", package = "OmicSelector")
  costs <- if (nzchar(path)) utils::read.delim(path, stringsAsFactors = FALSE) else data.frame()
  if (nrow(costs)) for (i in seq_len(nrow(x))) {
    z <- costs[costs$method == x$method_id[i] & costs$device == "cpu" &
                 costs$status == "complete", , drop = FALSE]
    if (!nrow(z)) { x$cost_status[i] <- "no completed CPU measurement"; next }
    sources <- c(fit_seconds = "fit_wall_seconds", score_seconds_per_specimen = "score_median_seconds",
                 peak_memory_bytes = "fit_peak_rss_bytes", frozen_model_bytes = "frozen_model_bytes")
    for (nm in names(sources)) {
      v <- as.numeric(z[[sources[[nm]]]])
      if (any(is.finite(v))) x[[nm]][i] <- stats::median(v[is.finite(v)])
    }
    x$measuring_node[i] <- paste(unique(paste(z$host, z$cpu_model, sep = ": ")), collapse = "; ")
    x$cost_status[i] <- "median of completed PB1 CPU cells; see measurements attribute"
  }
  x <- .method_select(x, methods, families)
  attr(x, "measurements") <- if (nrow(costs)) costs[costs$method %in% x$method_id, , drop = FALSE] else costs
  x
}

.method_r_available <- function(package) requireNamespace(package, quietly = TRUE)
.method_python_probe <- function(module) {
  if (!.method_r_available("reticulate")) return("R package reticulate is absent; Python not inspected")
  tryCatch({
    reticulate::import(module, delay_load = FALSE, convert = FALSE)
    ""
  }, error = function(e) paste0(module, ": ", gsub("[\r\n]+", " ", conditionMessage(e))))
}
.method_weight <- function(method, download = FALSE) {
  if (!method %in% c("ai-tabpfn", "tab-tabicl", "tab-tabdpt")) return("")
  tryCatch({
    if (method == "ai-tabpfn") {
      mod <- reticulate::import("tabpfn.model_loading", convert = FALSE)
      path <- reticulate::py_str(mod$get_cache_dir())
      path <- file.path(path, "tabpfn-v2-classifier-finetuned-zk73skhh.ckpt")
      if (!file.exists(path) && download) {
        Sys.setenv(TABPFN_NO_BROWSER = "1", TABPFN_DISABLE_TELEMETRY = "1")
        pathlib <- reticulate::import("pathlib", convert = FALSE)
        version <- reticulate::import("tabpfn.constants", convert = FALSE)$ModelVersion$V2
        mod$download_model(pathlib$Path(path), version = version, which = "classifier")
      }
      if (!file.exists(path) || file.info(path)$size <= 0) stop("missing weights: ", basename(path))
    } else {
      hub <- reticulate::import("huggingface_hub", convert = TRUE)
      repo <- if (method == "tab-tabicl") "jingang/TabICL" else "Layer6/TabDPT"
      filename <- if (method == "tab-tabicl") "tabicl-classifier-v2-20260212.ckpt" else "tabdpt1_2.safetensors"
      path <- hub$hf_hub_download(repo_id = repo, filename = filename, local_files_only = !download)
      if (!file.exists(path) || file.info(path)$size <= 0) stop("missing weights: ", filename)
    }
    ""
  }, error = function(e) paste0("weights for ", method, ": ", gsub("[\r\n]+", " ", conditionMessage(e))))
}

#' Inspect installed method dependencies without downloading
#' @inheritParams method_requirements
#' @param python Check imports in the active reticulate Python environment. FALSE
#'   leaves Python methods unverified (available is NA), unless an R dependency is absent.
#' @param weights Check cached foundation weights without network access.
#' @return Requirements plus available, missing_r, missing_python, missing_weights,
#'   reason and install_command columns. Availability concerns dependencies;
#'   panel, metadata and numerical validity still need the existing audits.
#' @export
available_methods <- function(methods = NULL, families = NULL, python = TRUE, weights = TRUE) {
  stopifnot(is.logical(python), length(python) == 1L, !is.na(python),
            is.logical(weights), length(weights) == 1L, !is.na(weights))
  x <- method_requirements(methods, families)
  x$available <- rep(NA, nrow(x))
  for (nm in c("missing_r", "missing_python", "missing_weights", "reason", "install_command")) x[[nm]] <- rep("", nrow(x))
  probe <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(x))) {
    x$install_command[i] <- sprintf('install_method_dependencies(methods = "%s")', x$method_id[i])
    if (!x$requirements_known[i]) { x$reason[i] <- "custom callback requirements are undeclared"; next }
    r <- .method_split(x$r_packages[i]); py <- .method_split(x$python_modules[i])
    missing <- r[!vapply(r, .method_r_available, logical(1))]
    x$missing_r[i] <- paste(missing, collapse = ";")
    unknown <- length(py) > 0L && !python
    errors <- character()
    if (length(py) && python && !length(missing)) for (module in py) {
      if (!exists(module, probe, inherits = FALSE)) assign(module, .method_python_probe(module), probe)
      error <- get(module, probe, inherits = FALSE)
      if (nzchar(error)) errors <- c(errors, error)
    }
    x$missing_python[i] <- paste(errors, collapse = "; ")
    if (length(py) && length(missing)) unknown <- TRUE
    if (nzchar(x$model_weights[i])) {
      if (!weights) unknown <- TRUE
      else if (!unknown && !length(errors)) x$missing_weights[i] <- .method_weight(x$method_id[i])
    }
    reasons <- c(if (length(missing)) paste0("R packages: ", paste(missing, collapse = ", ")),
                 errors, x$missing_weights[i])
    reasons <- reasons[nzchar(reasons)]
    x$available[i] <- if (length(reasons)) FALSE else if (unknown) NA else TRUE
    x$reason[i] <- if (length(reasons)) paste(reasons, collapse = "; ") else if (unknown) "Python or weights not checked" else "dependencies present"
  }
  x
}

#' Install dependencies for selected scoring methods
#' @inheritParams method_requirements
#' @param python Install Python modules and download required foundation weights.
#' @param dry_run Return the plan without installation, imports, or downloads.
#' @param envname Reticulate virtual environment name or path. A fresh R session
#'   is required if a different Python interpreter has already been initialized.
#' @param lib R library receiving newly installed packages.
#' @param repos CRAN repository for R dependencies.
#' @return Invisibly a dependency plan. Already installed packages and cached
#'   weights are retained. With python = FALSE only R dependencies are installed.
#' @export
install_method_dependencies <- function(methods = NULL, families = NULL, python = TRUE,
                                        dry_run = FALSE, envname = "omicselector",
                                        lib = .libPaths()[1L],
                                        repos = "https://cloud.r-project.org") {
  stopifnot(is.logical(python), length(python) == 1L, !is.na(python),
            is.logical(dry_run), length(dry_run) == 1L, !is.na(dry_run))
  x <- method_requirements(methods, families)
  if (any(!x$requirements_known)) stop("Dependencies of custom callbacks are undeclared: ",
    paste(x$method_id[!x$requirements_known], collapse = ", "), call. = FALSE)
  r <- .method_split(x$r_packages)
  py <- if (python) .method_split(x$python_modules) else character()
  distributions <- sub("\\..*$", "", py)
  distributions[distributions == "huggingface_hub"] <- "huggingface-hub"
  pins <- c(tabpfn = "tabpfn==8.0.8", tabicl = "tabicl==2.1.1", tabdpt = "tabdpt==1.2.0")
  packages <- unique(ifelse(distributions %in% names(pins), pins[distributions], distributions))
  missing <- r[!vapply(r, .method_r_available, logical(1))]
  plan <- list(methods = x$method_id, r_packages = r, missing_r = missing,
               python_packages = packages, weights = if (python) .method_split(x$model_weights) else character(),
               envname = envname, lib = lib, dry_run = dry_run)
  if (dry_run) return(plan)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  old_paths <- .libPaths()
  .libPaths(c(lib, old_paths))
  on.exit(.libPaths(old_paths), add = TRUE)
  missing <- r[!vapply(r, .method_r_available, logical(1))]
  if (length(missing)) {
    cran <- setdiff(missing, "selbal")
    if (length(cran)) utils::install.packages(cran, lib = lib, repos = repos)
    if ("selbal" %in% missing) {
      if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes", lib = lib, repos = repos)
      if (!requireNamespace("remotes", quietly = TRUE)) stop("Install R package remotes in lib, then retry.", call. = FALSE)
      remotes::install_github("malucalle/selbal", lib = lib, dependencies = NA, upgrade = "never")
    }
  }
  if (length(py)) {
    if (!requireNamespace("reticulate", quietly = TRUE)) stop("R package reticulate installation failed.", call. = FALSE)
    if (!reticulate::virtualenv_exists(envname)) {
      starter <- tryCatch(reticulate::virtualenv_starter(), error = function(e) NULL)
      if (length(starter) != 1L || is.na(starter) || !file.exists(starter)) {
        root <- reticulate::miniconda_path()
        starter <- file.path(root, if (.Platform$OS.type == "windows") "python.exe" else "bin/python")
        if (!file.exists(starter)) reticulate::install_miniconda(path = root, update = FALSE)
      }
      reticulate::virtualenv_create(envname, python = starter)
    }
    target <- reticulate::virtualenv_python(envname)
    if (reticulate::py_available(initialize = FALSE) &&
        normalizePath(reticulate::py_config()$python) != normalizePath(target))
      stop("A different Python is already initialized. Restart R and call install_method_dependencies() first; target: ", target, call. = FALSE)
    # reticulate::virtualenv_install adds --upgrade even when ignore_installed
    # is FALSE. Ordinary pip install preserves already satisfied requirements.
    status <- system2(target, c("-m", "pip", "install", "--disable-pip-version-check",
                                "--no-user", shQuote(packages)))
    if (!identical(status, 0L)) stop("Python dependency installation failed (exit ",
      status, "). Retry install_method_dependencies() in a fresh R session.", call. = FALSE)
    Sys.setenv(RETICULATE_PYTHON = target)
    reticulate::use_python(target, required = TRUE)
    for (id in x$method_id[nzchar(x$model_weights)]) {
      error <- .method_weight(id, download = TRUE)
      if (nzchar(error)) stop(error, "; retry install_method_dependencies(methods = \"", id, "\")", call. = FALSE)
    }
  }
  check <- available_methods(x$method_id, python = python, weights = python)
  if (any(check$available %in% FALSE)) stop(paste(check$reason[check$available %in% FALSE], collapse = "\n"), call. = FALSE)
  invisible(plan)
}
