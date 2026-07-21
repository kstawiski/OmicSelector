#!/usr/bin/env Rscript

# Exact, fail-closed runner for the OmicSelector paper's separately registered
# external rank-tree competitor sensitivity. Scientific results are produced
# only through the pre-existing grouped-CV benchmark engine after every source,
# cache, environment, package, and provenance input has been hash-bound.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.ranktree_methods <- function() {
  c("ws-rclr-trimmed", "reo-ktsp", "reo-pairratio", "reo-singscore",
    "reo-ucell", "reo-rankforest", "reo-rankboost")
}

.ranktree_seeds <- function() c(101L, 202L, 303L, 404L, 505L)

.ranktree_arg_value <- function(args, flag) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) NA_character_ else args[[hit + 1L]]
}

.ranktree_runner_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  flags <- c(
    package_root = "--package-root", paper_root = "--paper-root",
    analysis_plan = "--analysis-plan",
    analysis_amendment = "--analysis-amendment",
    analysis_amendment2 = "--analysis-amendment2",
    analysis_amendment3 = "--analysis-amendment3",
    cache = "--cache", output_dir = "--output-dir", methods = "--methods",
    provenance_manifest = "--provenance-manifest",
    provenance_union_inventory = "--provenance-union-inventory",
    runtime_image = "--runtime-image",
    runtime_attestation_manifest = "--runtime-attestation-manifest",
    engine_manifest = "--engine-manifest",
    environment_manifest = "--environment-manifest",
    expected_package_version = "--expected-package-version",
    expected_analysis_plan_sha256 = "--expected-analysis-plan-sha256",
    expected_analysis_amendment_sha256 =
      "--expected-analysis-amendment-sha256",
    expected_analysis_amendment2_sha256 =
      "--expected-analysis-amendment2-sha256",
    expected_analysis_amendment3_sha256 =
      "--expected-analysis-amendment3-sha256",
    expected_package_commit = "--expected-package-commit",
    expected_cache_commit = "--expected-cache-commit",
    expected_cache_sha256 = "--expected-cache-sha256",
    expected_cache_manifest_sha256 = "--expected-cache-manifest-sha256",
    expected_cache_digest_sha256 = "--expected-cache-digest-sha256",
    expected_runtime_image_sha256 = "--expected-runtime-image-sha256",
    expected_runtime_attestation_manifest_sha256 =
      "--expected-runtime-attestation-manifest-sha256",
    expected_container_launcher_sha256 =
      "--expected-container-launcher-sha256",
    expected_r_version = "--expected-r-version",
    expected_engine_manifest_sha256 = "--expected-engine-manifest-sha256",
    expected_environment_manifest_sha256 =
      "--expected-environment-manifest-sha256",
    expected_provenance_script_sha256 =
      "--expected-provenance-script-sha256",
    expected_provenance_manifest_sha256 =
      "--expected-provenance-manifest-sha256",
    expected_provenance_union_inventory_sha256 =
      "--expected-provenance-union-inventory-sha256",
    expected_analysis_code_id = "--expected-analysis-code-id"
  )
  out <- lapply(flags, function(flag) .ranktree_arg_value(args, flag))
  missing <- names(out)[vapply(out, function(x) is.na(x) || !nzchar(x), logical(1L))]
  if (length(missing)) {
    stop("Missing required argument(s): ",
         paste(unname(flags[missing]), collapse = ", "), call. = FALSE)
  }
  hex <- c(
    expected_analysis_plan_sha256 = 64L,
    expected_analysis_amendment_sha256 = 64L,
    expected_analysis_amendment2_sha256 = 64L,
    expected_analysis_amendment3_sha256 = 64L,
    expected_package_commit = 40L, expected_cache_commit = 40L,
    expected_cache_sha256 = 64L, expected_cache_manifest_sha256 = 64L,
    expected_cache_digest_sha256 = 64L,
    expected_runtime_image_sha256 = 64L,
    expected_runtime_attestation_manifest_sha256 = 64L,
    expected_container_launcher_sha256 = 64L,
    expected_engine_manifest_sha256 = 64L,
    expected_environment_manifest_sha256 = 64L,
    expected_provenance_script_sha256 = 64L,
    expected_provenance_manifest_sha256 = 64L,
    expected_provenance_union_inventory_sha256 = 64L,
    expected_analysis_code_id = 16L
  )
  for (name in names(hex)) {
    if (!grepl(sprintf("^[0-9a-f]{%d}$", hex[[name]]), out[[name]])) {
      stop(flags[[name]], " must be an exact lowercase hexadecimal pin of ",
           hex[[name]], " characters.", call. = FALSE)
    }
  }
  selected <- strsplit(out$methods, ",", fixed = TRUE)[[1L]]
  if (!length(selected) || any(!nzchar(selected)) || anyDuplicated(selected) ||
      length(setdiff(selected, .ranktree_methods()))) {
    stop("--methods must be a unique comma-separated subset of the seven ",
         "registered sensitivity methods.", call. = FALSE)
  }
  out$methods <- .ranktree_methods()[.ranktree_methods() %in% selected]
  out
}

.ranktree_sha <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

.ranktree_assert_file_pin <- function(path, expected, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  observed <- .ranktree_sha(path)
  if (!identical(observed, expected)) {
    stop(label, " does not match its expected SHA-256 pin.", call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

.ranktree_validate_analysis_plan <- function(path, paper_root,
                                             expected_sha256) {
  paper_root <- normalizePath(paper_root, mustWork = TRUE)
  canonical <- normalizePath(file.path(
    paper_root, "plan", "ranktree_external_competitor_sensitivity_20260719.md"
  ), mustWork = TRUE)
  path <- .ranktree_assert_file_pin(path, expected_sha256,
                                    "Approved rank-tree analysis plan")
  if (!identical(path, canonical)) {
    stop("Rank-tree analysis plan must be the exact canonical path under the ",
         "declared paper root.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  approved <- grepl(
    "^Status:[[:space:]]+(?:`|\\*\\*)?APPROVED\\b", lines, perl = TRUE
  )
  if (sum(approved) != 1L) {
    stop("Rank-tree analysis plan must contain exactly one explicit ",
         "'Status: APPROVED' marker.", call. = FALSE)
  }
  path
}

.ranktree_validate_analysis_amendment <- function(path, paper_root,
                                                  expected_sha256,
                                                  amendment_number = 1L) {
  amendment_number <- as.integer(amendment_number)
  if (length(amendment_number) != 1L || is.na(amendment_number) ||
      !amendment_number %in% c(1L, 2L, 3L)) {
    stop("Rank-tree amendment number must be 1, 2, or 3.", call. = FALSE)
  }
  paper_root <- normalizePath(paper_root, mustWork = TRUE)
  canonical <- normalizePath(file.path(
    paper_root, "plan", sprintf(
      "ranktree_external_competitor_amendment%d_20260720.md",
      amendment_number
    )
  ), mustWork = TRUE)
  path <- .ranktree_assert_file_pin(
    path, expected_sha256, "Approved rank-tree analysis amendment"
  )
  if (!identical(path, canonical)) {
    stop("Rank-tree analysis amendment must be the exact canonical path under ",
         "the declared paper root.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  approved <- grepl(
    "^Status:[[:space:]]+(?:`|\\*\\*)?APPROVED\\b", lines, perl = TRUE
  )
  if (sum(approved) != 1L) {
    stop("Rank-tree analysis amendment must contain exactly one explicit ",
         "'Status: APPROVED' marker.", call. = FALSE)
  }
  path
}

.ranktree_git_provenance <- function(package_root) {
  package_root <- normalizePath(package_root, mustWork = TRUE)
  git <- Sys.which("git")
  if (!nzchar(git)) stop("A git executable is required.", call. = FALSE)
  commit <- system2(git, c("-C", package_root, "rev-parse", "HEAD"),
                    stdout = TRUE, stderr = TRUE)
  status <- system2(git, c("-C", package_root, "status", "--porcelain=v1",
                           "--untracked-files=all"),
                    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(commit, "status")) || length(commit) != 1L ||
      !grepl("^[0-9a-f]{40}$", commit) || !is.null(attr(status, "status"))) {
    stop("Could not resolve exact OmicSelector git provenance.", call. = FALSE)
  }
  list(root = package_root, commit = unname(commit[[1L]]),
       dirty = length(status) > 0L)
}

.ranktree_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  }
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
}

.ranktree_validate_engine_manifest <- function(path, expected_sha256) {
  path <- .ranktree_assert_file_pin(path, expected_sha256, "Engine manifest")
  manifest <- data.table::fread(path)
  required <- c("role", "path", "sha256")
  roles <- c("benchmark_engine", "paired_auc_engine", "fold_engine",
             "public_design")
  if (!all(required %in% names(manifest)) || nrow(manifest) != length(roles) ||
      anyNA(manifest[, ..required]) || anyDuplicated(manifest$role) ||
      !setequal(manifest$role, roles) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Engine manifest must enumerate exactly the four frozen benchmark ",
         "source roles with complete SHA-256 pins.", call. = FALSE)
  }
  manifest[, path := vapply(seq_len(.N), function(i) {
    .ranktree_assert_file_pin(path[[i]], sha256[[i]],
                              paste0("Benchmark source ", role[[i]]))
  }, character(1L))]
  manifest <- manifest[match(roles, role)]
  manifest
}

.ranktree_installed_db <- function() {
  raw <- utils::installed.packages()
  packages <- unique(as.character(raw[, "Package"]))
  lib_paths <- vapply(unique(raw[, "LibPath"]), normalizePath, character(1L),
                      mustWork = TRUE)
  names(lib_paths) <- unique(raw[, "LibPath"])
  selected <- vapply(packages, function(package) {
    # Restrict resolution to the installed libraries. pkgload::load_all() makes
    # find.package() without lib.loc prefer the clean source checkout, which is
    # intentionally not an installed.packages() record and is pinned separately
    # by Git provenance.
    path <- find.package(package, lib.loc = .libPaths(), quiet = TRUE)
    if (!nzchar(path)) return(NA_integer_)
    library <- normalizePath(dirname(path), mustWork = TRUE)
    matches <- which(raw[, "Package"] == package &
                       unname(lib_paths[raw[, "LibPath"]]) == library)
    if (length(matches) == 1L) matches else NA_integer_
  }, integer(1L))
  if (anyNA(selected)) {
    stop("Could not resolve one active installed-library record per package.",
         call. = FALSE)
  }
  out <- raw[selected, , drop = FALSE]
  rownames(out) <- out[, "Package"]
  out
}

.ranktree_runtime_closure <- function(package = "ranktreeEnsemble") {
  installed <- .ranktree_installed_db()
  if (!package %in% rownames(installed)) {
    stop(package, " is not installed.", call. = FALSE)
  }
  dependencies <- tools::package_dependencies(
    package, db = installed, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  )[[package]]
  dependencies <- setdiff(unique(c(package, dependencies)), "R")
  missing <- setdiff(dependencies, rownames(installed))
  if (length(missing)) {
    stop("The installed ranktreeEnsemble dependency closure is incomplete: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  data.table::data.table(
    package = sort(dependencies),
    version = as.character(installed[sort(dependencies), "Version"])
  )
}

.ranktree_validate_environment_manifest <- function(
    path, expected_sha256, closure = .ranktree_runtime_closure()) {
  path <- .ranktree_assert_file_pin(path, expected_sha256,
                                    "Dependency environment manifest")
  manifest <- data.table::fread(path)
  required <- c("package", "version", "source_path", "source_sha256")
  if (!all(required %in% names(manifest)) || anyNA(manifest[, ..required]) ||
      anyDuplicated(manifest$package) ||
      any(!nzchar(manifest$package) | !nzchar(manifest$version) |
            !nzchar(manifest$source_path)) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$source_sha256))) {
    stop("Dependency environment manifest has an invalid schema or pin.",
         call. = FALSE)
  }
  closure <- data.table::as.data.table(closure)
  if (!all(c("package", "version") %in% names(closure)) ||
      anyDuplicated(closure$package) ||
      !setequal(manifest$package, closure$package)) {
    stop("Dependency environment manifest does not enumerate the exact ",
         "recursive ranktreeEnsemble closure.", call. = FALSE)
  }
  observed <- closure$version[match(manifest$package, closure$package)]
  if (any(manifest$version != observed)) {
    stop("Installed dependency versions differ from the frozen environment ",
         "manifest.", call. = FALSE)
  }
  if (!identical(manifest[package == "ranktreeEnsemble", version], "0.24")) {
    stop("The external sensitivity requires ranktreeEnsemble 0.24.", call. = FALSE)
  }
  manifest[, source_path := vapply(
    seq_len(.N), function(i) .ranktree_assert_file_pin(
      source_path[[i]], source_sha256[[i]],
      paste0("Dependency source archive for ", package[[i]])), character(1L)
  )]
  data.table::setorder(manifest, package)
  manifest
}

.ranktree_attestation_closure <- function(
    roots = c("data.table", "digest", "metafor", "pkgload",
              "ranktreeEnsemble")) {
  installed <- .ranktree_installed_db()
  missing <- setdiff(roots, rownames(installed))
  if (length(missing)) {
    stop("Runtime attestation is missing required package(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  dependencies <- unique(unlist(tools::package_dependencies(
    roots, db = installed, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  ), use.names = FALSE))
  packages <- sort(setdiff(unique(c(roots, dependencies)), "R"))
  data.table::data.table(
    package = packages,
    version = as.character(installed[packages, "Version"]),
    root_path = vapply(packages, function(package) {
      normalizePath(file.path(installed[package, "LibPath"], package),
                    mustWork = TRUE)
    }, character(1L))
  )
}

.ranktree_attestation_files <- function(root, scope, owner, version) {
  paths <- list.files(root, recursive = TRUE, all.files = TRUE,
                      full.names = TRUE, include.dirs = FALSE,
                      no.. = TRUE)
  info <- file.info(paths)
  paths <- paths[!is.na(info$isdir) & !info$isdir]
  info <- file.info(paths)
  if (!length(paths) || anyNA(info$size)) {
    stop("Runtime attestation found an empty or unreadable root: ", root,
         call. = FALSE)
  }
  root_prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
  normalized <- normalizePath(paths, mustWork = TRUE)
  relative <- substring(normalized, nchar(root_prefix) + 1L)
  if (any(!nzchar(relative)) || any(startsWith(relative, "../"))) {
    stop("Runtime attestation could not derive safe relative paths.",
         call. = FALSE)
  }
  data.table::data.table(
    scope = scope, owner = owner, version = version,
    root_path = normalizePath(root, mustWork = TRUE),
    relative_path = relative, bytes = as.numeric(info$size),
    sha256 = vapply(normalized, .ranktree_sha, character(1L))
  )
}

.ranktree_capture_runtime_attestation <- function() {
  image_sha <- Sys.getenv("OMICSELECTOR_RUNTIME_IMAGE_SHA256", unset = "")
  launcher_sha <- Sys.getenv("OMICSELECTOR_CONTAINER_LAUNCHER_SHA256",
                             unset = "")
  container_marker <- Sys.getenv("SINGULARITY_CONTAINER", unset = "")
  container_runtime <- "singularity"
  if (!nzchar(container_marker)) {
    container_marker <- Sys.getenv("APPTAINER_CONTAINER", unset = "")
    container_runtime <- "apptainer"
  }
  if (!grepl("^[0-9a-f]{64}$", image_sha) ||
      !grepl("^[0-9a-f]{64}$", launcher_sha) ||
      !nzchar(container_marker)) {
    stop("Runtime attestation must run inside the package launcher and an ",
         "active Singularity/Apptainer container.", call. = FALSE)
  }
  runtime_roots <- unique(normalizePath(
    c(R.home("bin"), R.home("etc"), R.home("lib"), R.home("modules")),
    mustWork = TRUE
  ))
  runtime <- data.table::rbindlist(lapply(runtime_roots, function(root) {
    .ranktree_attestation_files(
      root, "r_runtime", basename(root), R.version.string
    )
  }))
  closure <- .ranktree_attestation_closure()
  libraries <- data.table::rbindlist(lapply(seq_len(nrow(closure)), function(i) {
    .ranktree_attestation_files(
      closure$root_path[[i]], "r_library", closure$package[[i]],
      closure$version[[i]]
    )
  }))
  out <- data.table::rbindlist(list(runtime, libraries), use.names = TRUE)
  out[, `:=`(
    r_version = R.version.string,
    r_home = normalizePath(R.home(), mustWork = TRUE),
    container_runtime = container_runtime,
    container_marker = container_marker,
    runtime_image_sha256 = image_sha,
    container_launcher_sha256 = launcher_sha
  )]
  data.table::setorder(out, scope, owner, relative_path)
  if (anyDuplicated(out[, .(scope, owner, relative_path)]) ||
      any(!grepl("^[0-9a-f]{64}$", out$sha256))) {
    stop("Runtime attestation contains duplicate or invalid file records.",
         call. = FALSE)
  }
  out
}

.ranktree_write_runtime_attestation <- function(path) {
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  target <- file.path(parent, basename(path))
  if (file.exists(target)) {
    stop("Runtime attestation target already exists: ", target, call. = FALSE)
  }
  stage <- tempfile(".runtime-attestation-", tmpdir = parent)
  on.exit(if (file.exists(stage)) unlink(stage), add = TRUE)
  attestation <- .ranktree_capture_runtime_attestation()
  data.table::fwrite(attestation, stage,
                     sep = "\t", quote = FALSE)
  .ranktree_validate_runtime_attestation(
    stage, .ranktree_sha(stage), unique(attestation$r_version),
    unique(attestation$runtime_image_sha256),
    unique(attestation$container_launcher_sha256), current = attestation
  )
  if (!file.rename(stage, target)) {
    stop("Could not atomically write runtime attestation.", call. = FALSE)
  }
  cat("Created runtime attestation:", target, "\n")
  invisible(target)
}

.ranktree_validate_runtime_attestation <- function(
    path, expected_sha256, expected_r_version, expected_image_sha256,
    expected_launcher_sha256,
    current = .ranktree_capture_runtime_attestation()) {
  path <- .ranktree_assert_file_pin(path, expected_sha256,
                                    "Runtime closure attestation")
  manifest <- data.table::fread(path, colClasses = "character")
  required <- c(
    "scope", "owner", "version", "root_path", "relative_path", "bytes",
    "sha256", "r_version", "r_home", "container_runtime",
    "container_marker", "runtime_image_sha256", "container_launcher_sha256"
  )
  if (!identical(names(manifest), required) || !nrow(manifest) ||
      anyNA(manifest) || any(!nzchar(as.matrix(manifest))) ||
      anyDuplicated(manifest[, .(scope, owner, relative_path)]) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Runtime closure attestation has an invalid fail-closed schema.",
         call. = FALSE)
  }
  scalar <- function(column, expected, label) {
    observed <- unique(manifest[[column]])
    if (length(observed) != 1L || !identical(observed, expected)) {
      stop(label, " differs from the exact runtime attestation.",
           call. = FALSE)
    }
  }
  scalar("r_version", expected_r_version, "R version")
  scalar("runtime_image_sha256", expected_image_sha256,
         "Runtime image SHA-256")
  scalar("container_launcher_sha256", expected_launcher_sha256,
         "Container launcher SHA-256")
  current <- data.table::as.data.table(current)
  current[, bytes := as.character(bytes)]
  manifest[, bytes := as.character(bytes)]
  if (!identical(current, manifest)) {
    stop("In-process R runtime/library bytes differ from the exact frozen ",
         "attestation.", call. = FALSE)
  }
  list(path = path, sha256 = expected_sha256, manifest = manifest)
}

.ranktree_run_bound_provenance <- function(
    cohorts, paper_root, provenance_script, provenance_manifest,
    union_inventory, expected_script_sha256, expected_manifest_sha256,
    expected_union_sha256, log_path, preflight_fn) {
  paper_root <- normalizePath(paper_root, mustWork = TRUE)
  canonical <- c(
    script = file.path(paper_root, "code", "check_provenance.R"),
    manifest = file.path(paper_root, "data_bucket",
                         "_PROVENANCE_MANIFEST.md"),
    union = file.path(paper_root, "knowledge", "union_inventory.tsv")
  )
  supplied <- c(script = provenance_script, manifest = provenance_manifest,
                union = union_inventory)
  supplied <- vapply(supplied, normalizePath, character(1L), mustWork = TRUE)
  canonical <- vapply(canonical, normalizePath, character(1L), mustWork = TRUE)
  if (!identical(supplied, canonical)) {
    stop("Provenance inputs must be the exact canonical paths under the ",
         "declared paper root.", call. = FALSE)
  }
  expected <- c(expected_script_sha256, expected_manifest_sha256,
                expected_union_sha256)
  observed_before <- vapply(supplied, .ranktree_sha, character(1L))
  if (!identical(unname(observed_before), expected)) {
    stop("A canonical provenance input differs from its exact pin.",
         call. = FALSE)
  }
  inventory <- data.table::fread(supplied[["union"]], colClasses = "character")
  if (!"accession" %in% names(inventory) || !nrow(inventory) ||
      anyNA(inventory$accession) || any(!nzchar(inventory$accession)) ||
      anyDuplicated(toupper(inventory$accession))) {
    stop("Canonical provenance union inventory has an invalid accession key.",
         call. = FALSE)
  }
  source_accessions <- unique(toupper(unlist(lapply(cohorts, function(cohort) {
    value <- cohort$source_accessions
    if (is.null(value) || !length(value)) value <- cohort$preflight_accessions
    if (is.null(value) || !length(value)) value <- cohort$accession
    as.character(value)
  }), use.names = FALSE)))
  if (!length(source_accessions) || anyNA(source_accessions) ||
      any(!nzchar(source_accessions)) ||
      length(setdiff(source_accessions, toupper(inventory$accession)))) {
    stop("Every benchmark source accession must resolve in the exact union ",
         "inventory.", call. = FALSE)
  }
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(paper_root)
  result <- preflight_fn(cohorts, paper_root, log_path)
  observed_after <- vapply(supplied, .ranktree_sha, character(1L))
  if (!identical(observed_after, observed_before)) {
    stop("A canonical provenance input mutated during preflight.",
         call. = FALSE)
  }
  result
}

.ranktree_structural_reason <- function() {
  paste0("n_valid_folds=0 < 3; per-fold: ", paste(rep(
    "method/baseline constant or all-NA on test fold", 3L), collapse = "; "))
}

.ranktree_abort_engine <- function(stage, method, cohort, seed, fold,
                                   parent = NULL, detail = NULL) {
  message <- paste0(
    "Strict structural-proof engine failure [", stage, "] for ", method,
    " / ", cohort, " / seed ", seed, " / fold ", fold,
    if (!is.null(detail)) paste0(": ", detail) else if (!is.null(parent)) {
      paste0(": ", conditionMessage(parent))
    } else ""
  )
  condition <- structure(
    list(message = message, call = NULL, stage = stage, method = method,
         cohort = cohort, seed = as.integer(seed), fold = as.integer(fold),
         parent = parent),
    class = c("ranktree_engine_error", "error", "condition")
  )
  stop(condition)
}

.ranktree_strict_engine_call <- function(fn, args, stage, method, cohort, seed,
                                         fold) {
  tryCatch(
    do.call(fn, args),
    error = function(error) {
      .ranktree_abort_engine(stage, method, cohort, seed, fold, parent = error)
    }
  )
}

.ranktree_strict_score_vector <- function(value, n, stage, method, cohort,
                                          seed, fold) {
  if (!is.numeric(value) || length(value) != n || anyNA(value) ||
      any(!is.finite(value))) {
    .ranktree_abort_engine(
      stage, method, cohort, seed, fold,
      detail = "score must be a complete finite numeric held-out vector"
    )
  }
  as.numeric(value)
}

.ranktree_prove_structural_exception <- function(
    method, cohort_id, cohort, roster, seed, outer_k,
    grouped_fold_fn = get(".grouped_folds", envir = .GlobalEnv),
    stratified_fold_fn = get(".stratified_folds", envir = .GlobalEnv),
    baseline_fn = get(".srb_rclr_panel_score", envir = .GlobalEnv),
    fit_dispatch = get(".srb_call_fit", envir = .GlobalEnv),
    score_dispatch = OmicSelector::singlesample_score_call,
    direction_fn = get(".srb_training_direction", envir = .GlobalEnv)) {
  if (!identical(cohort_id, "GSE83977") || !identical(as.integer(outer_k), 3L) ||
      !as.integer(seed) %in% .ranktree_seeds()) {
    stop("Structural proof is restricted to the registered GSE83977/3-fold ",
         "exception.", call. = FALSE)
  }
  X <- cohort$expr_per_sample
  y <- as.integer(cohort$y_bin)
  group_id <- cohort$group_id
  if (!is.matrix(X) || nrow(X) != length(y) ||
      any(!y %in% c(0L, 1L)) || length(unique(y)) != 2L) {
    stop("Structural-proof cohort inputs are invalid.", call. = FALSE)
  }
  has_groups <- !is.null(group_id) && length(group_id) == nrow(X)
  folds <- if (has_groups) {
    grouped_fold_fn(y, group_id, k = outer_k, seed = seed)
  } else {
    stratified_fold_fn(y, k = outer_k, seed = seed)
  }
  if (!is.list(folds) || length(folds) != outer_k) {
    stop("Structural proof did not reproduce the exact outer-fold count.",
         call. = FALSE)
  }
  row <- roster[roster$method_id == method, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Structural proof did not resolve one registered method.", call. = FALSE)
  }
  proof <- vector("list", outer_k)
  for (fold in seq_len(outer_k)) {
    test <- as.integer(folds[[fold]])
    train <- setdiff(seq_len(nrow(X)), test)
    if (!length(test) || !length(train) || anyDuplicated(test) ||
        any(test < 1L | test > nrow(X)) ||
        (has_groups && length(intersect(unique(group_id[test]),
                                        unique(group_id[train]))))) {
      stop("Structural proof found an invalid or group-leaking outer fold.",
           call. = FALSE)
    }
    y_train <- y[train]
    y_test <- y[test]
    if (length(unique(y_train)) != 2L || length(unique(y_test)) != 2L ||
        sum(y_test == 0L) < cohort$min_per_fold ||
        sum(y_test == 1L) < cohort$min_per_fold) {
      stop("Structural exception is not attributable solely to constant ",
           "score support.", call. = FALSE)
    }
    X_train <- X[train, , drop = FALSE]
    X_test <- X[test, , drop = FALSE]
    baseline <- .ranktree_strict_engine_call(
      baseline_fn, list(X_train, y_train, X_test, panel_size = 20L),
      "baseline", method, cohort_id, seed, fold
    )
    baseline <- .ranktree_strict_score_vector(
      baseline, length(test), "baseline", method, cohort_id, seed, fold
    )
    if (identical(method, "ws-rclr-trimmed")) {
      method_score <- baseline
    } else {
      fit_name <- as.character(row$fit_fn[[1L]])
      if (is.na(fit_name) || !nzchar(fit_name)) {
        .ranktree_abort_engine(
          "fit", method, cohort_id, seed, fold,
          detail = "registered fit function is absent"
        )
      }
      model <- .ranktree_strict_engine_call(
        fit_dispatch, list(fit_name, X_train, y_train, NULL, list()),
        "fit", method, cohort_id, seed, fold
      )
      training_score <- .ranktree_strict_engine_call(
        score_dispatch,
        list(method, model, X_train, NULL, roster = roster),
        "training_prediction", method, cohort_id, seed, fold
      )
      training_score <- .ranktree_strict_score_vector(
        training_score, length(train), "training_prediction", method,
        cohort_id, seed, fold
      )
      direction <- .ranktree_strict_engine_call(
        direction_fn, list(training_score, y_train), "direction", method,
        cohort_id, seed, fold
      )
      if (!is.numeric(direction) || length(direction) != 1L ||
          !is.finite(direction) || !direction %in% c(-1, 1)) {
        .ranktree_abort_engine(
          "direction", method, cohort_id, seed, fold,
          detail = "training-only direction is unavailable"
        )
      }
      method_score <- .ranktree_strict_engine_call(
        score_dispatch, list(method, model, X_test, NULL, roster = roster),
        "heldout_prediction", method, cohort_id, seed, fold
      )
      method_score <- direction * .ranktree_strict_score_vector(
        method_score, length(test), "heldout_prediction", method, cohort_id,
        seed, fold
      )
    }
    baseline_constant <- data.table::uniqueN(baseline) < 2L
    method_constant <- data.table::uniqueN(method_score) < 2L
    if (!baseline_constant && !method_constant) {
      stop("Strict replay did not prove constant held-out support in every ",
           "outer fold.", call. = FALSE)
    }
    proof[[fold]] <- data.table::data.table(
      fold = as.integer(fold), n_train = length(train), n_test = length(test),
      n_positive_test = sum(y_test == 1L),
      n_negative_test = sum(y_test == 0L),
      baseline_constant = baseline_constant, method_constant = method_constant,
      test_index_sha256 = digest::digest(test, algo = "sha256", serialize = TRUE),
      baseline_score_sha256 = digest::digest(
        baseline, algo = "sha256", serialize = TRUE),
      method_score_sha256 = digest::digest(
        method_score, algo = "sha256", serialize = TRUE)
    )
  }
  proof <- data.table::rbindlist(proof)
  list(
    proven = TRUE,
    basis = paste0(
      "complete finite replay: method or rCLR baseline constant on every ",
      "exact group-safe held-out fold"
    ),
    proof = proof,
    proof_sha256 = digest::digest(proof, algo = "sha256", serialize = TRUE)
  )
}

.ranktree_validate_cell_result <- function(cell) {
  required <- c("method", "cohort", "seed", "eligible", "ineligible_reason",
                "n_eff", "n_valid_folds", "outer_k",
                "n_group_split_violations", "structural_exception_proven",
                "structural_exception_proof_sha256")
  if (!is.data.frame(cell) || nrow(cell) != 1L ||
      !all(required %in% names(cell)) || !is.logical(cell$eligible) ||
      is.na(cell$eligible)) {
    stop("Benchmark cell result has an invalid schema.", call. = FALSE)
  }
  if (isTRUE(cell$eligible)) {
    if (is.na(cell$n_eff) || cell$n_eff < 1L ||
        is.na(cell$n_valid_folds) || is.na(cell$outer_k) ||
        cell$n_valid_folds != cell$outer_k ||
        is.na(cell$n_group_split_violations) ||
        cell$n_group_split_violations != 0L ||
        (!is.na(cell$ineligible_reason) && nzchar(cell$ineligible_reason)) ||
        !identical(cell$structural_exception_proven, FALSE) ||
        !is.na(cell$structural_exception_proof_sha256)) {
      stop("Eligible benchmark cell violates fold/group accounting.", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  accepted <- identical(as.character(cell$cohort), "GSE83977") &&
    as.integer(cell$seed) %in% .ranktree_seeds() &&
    identical(as.integer(cell$n_valid_folds), 0L) &&
    identical(as.integer(cell$outer_k), 3L) &&
    identical(as.integer(cell$n_group_split_violations), 0L) &&
    identical(as.character(cell$ineligible_reason), .ranktree_structural_reason()) &&
    identical(cell$structural_exception_proven, TRUE) &&
    is.character(cell$structural_exception_proof_sha256) &&
    length(cell$structural_exception_proof_sha256) == 1L &&
    grepl("^[0-9a-f]{64}$", cell$structural_exception_proof_sha256)
  if (!accepted) {
    stop("Unreviewed ineligibility, fit/score failure, or dependency/runtime ",
         "failure: ", as.character(cell$ineligible_reason), call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_cell_row <- function(result, seed, package_pins) {
  `%rt||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
  data.table::data.table(
    method = result$method_id, cohort = result$cohort, seed = as.integer(seed),
    eligible = isTRUE(result$eligible),
    ineligible_reason = result$ineligible_reason %rt||% NA_character_,
    auc_method = result$auc_method %rt||% NA_real_,
    auc_baseline = result$auc_baseline %rt||% NA_real_,
    lift = result$est %rt||% NA_real_, lift_se = result$se %rt||% NA_real_,
    n_eff = result$n_eff %rt||% NA_integer_,
    n_valid_folds = result$n_valid_folds %rt||% NA_integer_,
    min_pos = result$min_pos %rt||% NA_integer_,
    min_neg = result$min_neg %rt||% NA_integer_,
    outer_k = result$outer_k %rt||% NA_integer_,
    grouping_policy = result$grouping_policy %rt||% NA_character_,
    n_group_split_violations =
      result$n_group_split_violations %rt||% NA_integer_,
    structural_exception_proven = FALSE,
    structural_exception_basis = NA_character_,
    structural_exception_proof_sha256 = NA_character_,
    package_version = package_pins$version,
    package_commit = package_pins$commit,
    cache_package_commit = package_pins$cache_commit,
    cache_sha256 = package_pins$cache_sha256,
    engine_manifest_sha256 = package_pins$engine_manifest_sha256,
    environment_manifest_sha256 = package_pins$environment_manifest_sha256,
    runtime_image_sha256 = package_pins$runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      package_pins$runtime_attestation_manifest_sha256,
    container_launcher_sha256 = package_pins$container_launcher_sha256,
    r_version = package_pins$r_version,
    analysis_plan_sha256 = package_pins$analysis_plan_sha256,
    analysis_amendment_sha256 = package_pins$analysis_amendment_sha256,
    analysis_amendment2_sha256 = package_pins$analysis_amendment2_sha256,
    analysis_amendment3_sha256 = package_pins$analysis_amendment3_sha256,
    provenance_manifest_sha256 = package_pins$provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      package_pins$provenance_union_inventory_sha256,
    analysis_code_id = package_pins$analysis_code_id
  )
}

.ranktree_add_split_ids <- function(predictions) {
  predictions <- data.table::as.data.table(predictions)
  if (!nrow(predictions)) return(predictions[, split_id := character()])
  required <- c("method", "cohort", "seed", "fold", "sample_id", "group_id",
                "y", "score_m", "n_train", "n_test")
  if (!all(required %in% names(predictions)) ||
      anyNA(predictions[, ..required]) || any(!is.finite(predictions$score_m)) ||
      any(!predictions$y %in% c(0L, 1L)) ||
      anyDuplicated(predictions[, .(method, cohort, seed, sample_id)])) {
    stop("Held-out predictions are incomplete, invalid, or duplicated.",
         call. = FALSE)
  }
  group_fold <- predictions[, .(n_fold = data.table::uniqueN(fold),
                                n_y = data.table::uniqueN(y)),
                            by = .(method, cohort, seed, group_id)]
  if (any(group_fold$n_fold != 1L) || any(group_fold$n_y != 1L)) {
    stop("A provenance group crosses a held-out fold or outcome label.",
         call. = FALSE)
  }
  split_map <- predictions[, {
    reference <- .SD[, .(fold = as.integer(fold), sample_id = as.character(sample_id),
                         group_id = as.character(group_id), y = as.integer(y),
                         n_train = as.integer(n_train), n_test = as.integer(n_test))]
    data.table::setorder(reference, fold, sample_id, group_id, y)
    .(split_id = digest::digest(reference, algo = "sha256", serialize = TRUE))
  }, by = .(method, cohort, seed)]
  out <- merge(predictions, split_map,
               by = c("method", "cohort", "seed"), all.x = TRUE, sort = FALSE)
  cross_method <- unique(out[, .(method, cohort, seed, split_id)])[
    , .(n_split = data.table::uniqueN(split_id)), by = .(cohort, seed)]
  if (any(cross_method$n_split != 1L)) {
    stop("Selected methods do not share byte-identical held-out splits.",
         call. = FALSE)
  }
  out
}

.ranktree_units <- function(cohorts) {
  data.table::rbindlist(lapply(names(cohorts), function(id) {
    cohort <- cohorts[[id]]
    data.table::data.table(
      unit_id = id, disease = as.character(cohort$disease),
      modality = as.character(cohort$modality),
      biospecimen = as.character(cohort$biospecimen),
      provenance_block = as.character(cohort$provenance_block),
      source_accessions = paste(as.character(cohort$source_accessions),
                                collapse = ";"),
      outer_k = as.integer(cohort$outer_k),
      primary_unit = !grepl("whole[ -]?blood", cohort$biospecimen,
                            ignore.case = TRUE),
      n_deposited_profiles = nrow(cohort$expr_per_sample),
      n_unique_groups = data.table::uniqueN(as.character(cohort$group_id))
    )
  }), fill = TRUE)
}

.ranktree_validate_runtime_bundle <- function(cells, predictions, units,
                                               methods,
                                               runtime_image_sha256 = NULL,
                                               runtime_attestation_sha256 = NULL,
                                               container_launcher_sha256 = NULL,
                                               r_version = NULL) {
  expected <- data.table::CJ(method = methods, cohort = units$unit_id,
                             seed = .ranktree_seeds(), unique = TRUE)
  if (nrow(units) != 34L || sum(units$primary_unit) != 33L ||
      anyDuplicated(units$unit_id) ||
      nrow(cells) != nrow(expected) ||
      anyDuplicated(cells[, .(method, cohort, seed)]) ||
      !setequal(do.call(paste, cells[, .(method, cohort, seed)]),
                do.call(paste, expected))) {
    stop("Runtime bundle is not the exact method-by-34-unit-by-five-seed grid.",
         call. = FALSE)
  }
  invisible(lapply(seq_len(nrow(cells)), function(i) {
    .ranktree_validate_cell_result(cells[i])
  }))
  if (any(!cells$eligible)) {
    stop("Amendment 2 requires every rank-tree sensitivity cell to satisfy ",
         "the complete-finite held-out support contract.", call. = FALSE)
  }
  if (!is.null(runtime_image_sha256) &&
      (!"runtime_image_sha256" %in% names(cells) ||
       !"runtime_image_sha256" %in% names(predictions) ||
       anyNA(cells$runtime_image_sha256) ||
       anyNA(predictions$runtime_image_sha256) ||
       any(cells$runtime_image_sha256 != runtime_image_sha256) ||
       any(predictions$runtime_image_sha256 != runtime_image_sha256))) {
    stop("Runtime image pin is absent or inconsistent in cells/predictions.",
         call. = FALSE)
  }
  runtime_columns <- c(
    runtime_attestation_manifest_sha256 = runtime_attestation_sha256,
    container_launcher_sha256 = container_launcher_sha256,
    r_version = r_version
  )
  for (column in names(runtime_columns)) {
    expected_value <- runtime_columns[[column]]
    if (!is.null(expected_value) &&
        (!column %in% names(cells) || !column %in% names(predictions) ||
         anyNA(cells[[column]]) || anyNA(predictions[[column]]) ||
         any(cells[[column]] != expected_value) ||
         any(predictions[[column]] != expected_value))) {
      stop("Runtime attestation field is absent or inconsistent: ", column,
           call. = FALSE)
    }
  }
  eligible <- cells[eligible == TRUE, .(method, cohort, seed, n_eff, outer_k)]
  pred_cells <- unique(predictions[, .(method, cohort, seed)])
  if (!setequal(do.call(paste, eligible[, .(method, cohort, seed)]),
                do.call(paste, pred_cells))) {
    stop("Eligible cells and held-out prediction cells disagree.", call. = FALSE)
  }
  accounting <- merge(
    eligible, predictions[, .(n_prediction = .N,
                              n_fold = data.table::uniqueN(fold)),
                          by = .(method, cohort, seed)],
    by = c("method", "cohort", "seed"), all.x = TRUE
  )
  if (anyNA(accounting$n_prediction) ||
      any(accounting$n_prediction != accounting$n_eff) ||
      any(accounting$n_fold != accounting$outer_k)) {
    stop("Eligible prediction cardinality/fold accounting is incomplete.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_copy <- function(from, to, label) {
  if (!file.copy(from, to, overwrite = FALSE, copy.mode = TRUE)) {
    stop("Could not copy ", label, ".", call. = FALSE)
  }
  invisible(to)
}

.ranktree_checkpoint_contract <- function(output_dir, methods, units, pins,
                                          runner_sha, r_version) {
  data.table::as.data.table(data.frame(
    key = c(
      "schema_version", "output_dir", "methods", "seeds", "units",
      "package_version", "package_commit", "cache_package_commit",
      "cache_sha256", "cache_manifest_sha256", "cache_digest_sha256",
      "runner_script_sha256", "engine_manifest_sha256",
      "environment_manifest_sha256", "runtime_image_sha256",
      "runtime_attestation_manifest_sha256",
      "container_launcher_sha256",
      "analysis_plan_sha256", "analysis_amendment_sha256",
      "analysis_amendment2_sha256",
      "analysis_amendment3_sha256",
      "provenance_script_sha256",
      "provenance_manifest_sha256", "provenance_union_inventory_sha256",
      "analysis_code_id", "r_version"
    ),
    value = c(
      "omicselector-paper3-ranktree-checkpoint-v2",
      normalizePath(output_dir, mustWork = FALSE),
      paste(methods, collapse = ","), paste(.ranktree_seeds(), collapse = ","),
      paste(units$unit_id, collapse = ","), pins$version, pins$commit,
      pins$cache_commit, pins$cache_sha256, pins$cache_manifest_sha256,
      pins$cache_digest_sha256, runner_sha, pins$engine_manifest_sha256,
      pins$environment_manifest_sha256, pins$runtime_image_sha256,
      pins$runtime_attestation_manifest_sha256,
      pins$container_launcher_sha256,
      pins$analysis_plan_sha256,
      pins$analysis_amendment_sha256,
      pins$analysis_amendment2_sha256,
      pins$analysis_amendment3_sha256,
      pins$provenance_script_sha256,
      pins$provenance_manifest_sha256,
      pins$provenance_union_inventory_sha256,
      pins$analysis_code_id, r_version
    ), stringsAsFactors = FALSE
  ))
}

.ranktree_validate_checkpoint_contract <- function(checkpoint_dir, expected) {
  if (!dir.exists(checkpoint_dir) || nzchar(Sys.readlink(checkpoint_dir))) {
    stop("Persistent checkpoint path is absent, not a directory, or a symlink.",
         call. = FALSE)
  }
  path <- file.path(checkpoint_dir, "checkpoint_contract.tsv")
  if (!file.exists(path)) {
    stop("Persistent checkpoint directory has an incomplete contract.",
         call. = FALSE)
  }
  observed <- tryCatch(data.table::fread(
    path, colClasses = c(key = "character", value = "character")),
    error = identity)
  if (inherits(observed, "error") ||
      !identical(names(observed), c("key", "value")) ||
      nrow(observed) != nrow(expected) || anyNA(observed) ||
      anyDuplicated(observed$key) || !identical(observed, expected)) {
    stop("Persistent checkpoint contract is mismatched or incomplete.",
         call. = FALSE)
  }
  invisible(path)
}

.ranktree_open_checkpoints <- function(output_dir, expected) {
  checkpoint_dir <- paste0(normalizePath(output_dir, mustWork = FALSE),
                           ".checkpoints")
  if (file.exists(checkpoint_dir)) {
    .ranktree_validate_checkpoint_contract(checkpoint_dir, expected)
    return(checkpoint_dir)
  }
  parent <- normalizePath(dirname(checkpoint_dir), mustWork = TRUE)
  stage <- tempfile(".ranktree-checkpoint-init-", tmpdir = parent)
  if (!dir.create(stage, recursive = FALSE)) {
    stop("Could not create checkpoint-contract staging directory.",
         call. = FALSE)
  }
  promoted <- FALSE
  on.exit(if (!promoted && dir.exists(stage)) unlink(stage, recursive = TRUE),
          add = TRUE)
  data.table::fwrite(expected, file.path(stage, "checkpoint_contract.tsv"),
                     sep = "\t", quote = FALSE)
  if (!file.rename(stage, checkpoint_dir)) {
    stop("Could not atomically create the persistent checkpoint directory.",
         call. = FALSE)
  }
  promoted <- TRUE
  .ranktree_validate_checkpoint_contract(checkpoint_dir, expected)
  checkpoint_dir
}

.ranktree_bundle_artifacts <- function() {
  c("analysis_plan.md", "analysis_amendment.md", "analysis_amendment2.md",
    "analysis_amendment3.md",
    "engine_manifest.tsv", "environment_manifest.tsv",
    "methods.tsv",
    "provenance_preflight.log", "provenance_resolution.tsv", "report.md",
    "runtime_attestation.tsv",
    "units.tsv", "within_bundle.rds", "within_cells.tsv",
    "within_predictions.tsv.gz")
}

.ranktree_validate_completed_output <- function(output_dir) {
  if (!dir.exists(output_dir) || nzchar(Sys.readlink(output_dir))) {
    stop("Final output is absent, not a directory, or a symlink.",
         call. = FALSE)
  }
  manifest_path <- file.path(output_dir, "output_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop("Final output lacks its completion manifest.", call. = FALSE)
  }
  manifest <- data.table::fread(manifest_path)
  required <- .ranktree_bundle_artifacts()
  if (!all(c("file", "bytes", "sha256") %in% names(manifest)) ||
      anyNA(manifest[, .(file, bytes, sha256)]) || anyDuplicated(manifest$file) ||
      !setequal(manifest$file, required) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Final output completion manifest is malformed or incomplete.",
         call. = FALSE)
  }
  paths <- file.path(output_dir, manifest$file)
  if (any(!file.exists(paths)) || any(file.info(paths)$size != manifest$bytes) ||
      !identical(unname(vapply(paths, .ranktree_sha, character(1L))),
                 as.character(manifest$sha256))) {
    stop("Final output bytes differ from the completion manifest.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_remove_completed_checkpoints <- function(output_dir, expected) {
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  checkpoint_dir <- paste0(output_dir, ".checkpoints")
  .ranktree_validate_completed_output(output_dir)
  .ranktree_validate_checkpoint_contract(checkpoint_dir, expected)
  status <- unlink(checkpoint_dir, recursive = TRUE, force = FALSE)
  if (!identical(status, 0L) || dir.exists(checkpoint_dir)) {
    stop("Final output is complete but its checkpoint directory could not be ",
         "removed.", call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_strict_exception_cell <- function(method, cohort_id, cohort, roster,
                                            out_dir, seed, outer_k, ...) {
  if (!identical(cohort_id, "GSE83977")) {
    stop("Strict exception cell is restricted to GSE83977.", call. = FALSE)
  }
  checkpoint <- file.path(
    out_dir,
    paste0("strict__", cohort_id, "__", gsub("[^A-Za-z0-9_.-]", "_", method),
           ".rds")
  )
  if (file.exists(checkpoint)) {
    cached <- readRDS(checkpoint)
    proof <- cached$structural_exception_proof
    if (!is.list(cached) || !identical(cached$method_id, method) ||
        !identical(cached$cohort, cohort_id) ||
        !identical(as.integer(cached$seed), as.integer(seed)) ||
        !identical(as.integer(cached$outer_k), as.integer(outer_k)) ||
        !is.data.frame(proof) || nrow(proof) != outer_k ||
        !identical(
          digest::digest(proof, algo = "sha256", serialize = TRUE),
          cached$structural_exception_proof_sha256
        )) {
      stop("Strict structural-exception checkpoint is invalid.", call. = FALSE)
    }
    return(cached)
  }
  replay <- .ranktree_prove_structural_exception(
    method, cohort_id, cohort, roster, seed, outer_k, ...
  )
  result <- list(
    method_id = method, cohort = cohort_id, seed = as.integer(seed),
    eligible = FALSE, ineligible_reason = .ranktree_structural_reason(),
    n_eff = NA_integer_, n_valid_folds = 0L,
    min_pos = NA_integer_, min_neg = NA_integer_, outer_k = as.integer(outer_k),
    grouping_policy = if (!is.null(cohort$group_id)) {
      "groupKFold"
    } else "ungrouped_stratified",
    n_group_split_violations = 0L, fold_scores = NULL,
    structural_exception_proven = replay$proven,
    structural_exception_basis = replay$basis,
    structural_exception_proof = replay$proof,
    structural_exception_proof_sha256 = replay$proof_sha256
  )
  stage <- tempfile(".strict-cell-", tmpdir = out_dir, fileext = ".rds")
  on.exit(if (file.exists(stage)) unlink(stage), add = TRUE)
  saveRDS(result, stage, version = 3L)
  if (!file.rename(stage, checkpoint)) {
    stop("Could not atomically persist strict structural-exception cell.",
         call. = FALSE)
  }
  result
}

.ranktree_call_benchmark_cell <- function(
    benchmark_fn, method, cohort_id, cohort, roster, out_dir, seed, outer_k,
    strict_exception_fn = .ranktree_strict_exception_cell) {
  # Amendment 2 applies one uniform complete-finite score contract. Constant
  # held-out scores are valid midrank-AUC evidence, including for GSE83977; the
  # historical strict-exception helper is retained only for audit replay.
  benchmark_fn(
    method, cohort_id, cohort, roster, out_dir = out_dir, seed = seed,
    outer_k = outer_k, min_valid_folds = outer_k,
    min_pos_per_fold = cohort$min_per_fold,
    min_neg_per_fold = cohort$min_per_fold,
    baseline_mode = "rclr_panel", panel_size = 20L,
    overwrite = FALSE, rerun_present = TRUE, orient = "fixed", hp = list()
  )
}

.ranktree_runner_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .ranktree_runner_args(args)
  code_start <- .ranktree_git_provenance(opt$package_root)
  if (code_start$dirty) {
    stop("Rank-tree sensitivity execution requires a clean OmicSelector ",
         "worktree.", call. = FALSE)
  }
  if (!identical(code_start$commit, opt$expected_package_commit)) {
    stop("OmicSelector checkout does not match the expected package commit.",
         call. = FALSE)
  }
  description <- read.dcf(file.path(code_start$root, "DESCRIPTION"))
  package_version <- unname(description[, "Version"])
  if (!identical(package_version, opt$expected_package_version)) {
    stop("OmicSelector DESCRIPTION version does not match its expected pin.",
         call. = FALSE)
  }
  script_path <- .ranktree_script_path()
  expected_script <- normalizePath(file.path(
    code_start$root, "inst", "scripts", basename(script_path)), mustWork = TRUE)
  if (!identical(script_path, expected_script)) {
    stop("Runner must execute from the pinned OmicSelector checkout.",
         call. = FALSE)
  }
  runner_sha <- .ranktree_sha(script_path)
  launcher_path <- .ranktree_assert_file_pin(
    file.path(code_start$root, "inst", "scripts",
              "paper3_run_ranktree_sensitivity_container.sh"),
    opt$expected_container_launcher_sha256, "Package container launcher"
  )

  engine_manifest <- .ranktree_validate_engine_manifest(
    opt$engine_manifest, opt$expected_engine_manifest_sha256)
  environment_manifest <- .ranktree_validate_environment_manifest(
    opt$environment_manifest, opt$expected_environment_manifest_sha256)
  runtime_image <- .ranktree_assert_file_pin(
    opt$runtime_image, opt$expected_runtime_image_sha256,
    "Singularity runtime image")
  if (!identical(R.version.string, opt$expected_r_version)) {
    stop("In-process R version differs from its exact pre-execution pin.",
         call. = FALSE)
  }
  runtime_attestation <- .ranktree_validate_runtime_attestation(
    opt$runtime_attestation_manifest,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_r_version, opt$expected_runtime_image_sha256,
    opt$expected_container_launcher_sha256
  )

  paper_root <- normalizePath(opt$paper_root, mustWork = TRUE)
  analysis_plan <- .ranktree_validate_analysis_plan(
    opt$analysis_plan, paper_root, opt$expected_analysis_plan_sha256
  )
  analysis_amendment <- .ranktree_validate_analysis_amendment(
    opt$analysis_amendment, paper_root,
    opt$expected_analysis_amendment_sha256, amendment_number = 1L
  )
  analysis_amendment2 <- .ranktree_validate_analysis_amendment(
    opt$analysis_amendment2, paper_root,
    opt$expected_analysis_amendment2_sha256, amendment_number = 2L
  )
  analysis_amendment3 <- .ranktree_validate_analysis_amendment(
    opt$analysis_amendment3, paper_root,
    opt$expected_analysis_amendment3_sha256, amendment_number = 3L
  )
  cache <- .ranktree_assert_file_pin(opt$cache, opt$expected_cache_sha256,
                                    "Cohort cache")
  cache_manifest_path <- file.path(dirname(cache), "cohort_cache_manifest.tsv")
  cache_digest_path <- file.path(dirname(cache), "cohort_cache_digest.tsv")
  .ranktree_assert_file_pin(cache_manifest_path,
                            opt$expected_cache_manifest_sha256,
                            "Cohort cache manifest")
  .ranktree_assert_file_pin(cache_digest_path,
                            opt$expected_cache_digest_sha256,
                            "Cohort cache digest")
  provenance_script <- .ranktree_assert_file_pin(
    file.path(paper_root, "code", "check_provenance.R"),
    opt$expected_provenance_script_sha256, "Provenance preflight script")
  provenance_manifest <- .ranktree_assert_file_pin(
    opt$provenance_manifest,
    opt$expected_provenance_manifest_sha256, "Canonical provenance manifest")
  provenance_union_inventory <- .ranktree_assert_file_pin(
    opt$provenance_union_inventory,
    opt$expected_provenance_union_inventory_sha256,
    "Canonical provenance union inventory")

  analysis_code_id <- substr(digest::digest(paste(
    runner_sha, opt$expected_package_commit, opt$expected_cache_sha256,
    opt$expected_analysis_plan_sha256,
    opt$expected_analysis_amendment_sha256,
    opt$expected_analysis_amendment2_sha256,
    opt$expected_analysis_amendment3_sha256,
    opt$expected_engine_manifest_sha256,
    opt$expected_environment_manifest_sha256,
    opt$expected_runtime_image_sha256,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_container_launcher_sha256, opt$expected_r_version,
    opt$expected_provenance_script_sha256,
    opt$expected_provenance_manifest_sha256,
    opt$expected_provenance_union_inventory_sha256, sep = "|"),
    algo = "sha256"),
    1L, 16L)
  if (!identical(analysis_code_id, opt$expected_analysis_code_id)) {
    stop("Computed analysis code id does not match its pre-execution pin.",
         call. = FALSE)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The exact sensitivity runner requires pkgload.", call. = FALSE)
  }
  suppressMessages(suppressWarnings(
    pkgload::load_all(code_start$root, quiet = TRUE, export_all = FALSE,
                      helpers = FALSE)
  ))
  if (!identical(as.character(utils::packageVersion("OmicSelector")),
                 opt$expected_package_version)) {
    stop("Loaded OmicSelector version differs from the pinned checkout.",
         call. = FALSE)
  }

  old_engine_root <- Sys.getenv("OMICSELECTOR_PAPER_ROOT", unset = NA_character_)
  on.exit({
    if (is.na(old_engine_root)) Sys.unsetenv("OMICSELECTOR_PAPER_ROOT")
    else Sys.setenv(OMICSELECTOR_PAPER_ROOT = old_engine_root)
  }, add = TRUE)
  engine_root <- dirname(dirname(dirname(
    engine_manifest[role == "benchmark_engine", path]
  )))
  Sys.setenv(OMICSELECTOR_PAPER_ROOT = engine_root)
  for (role_value in c("paired_auc_engine", "fold_engine", "benchmark_engine",
                       "public_design")) {
    source_path <- engine_manifest[role == role_value, path]
    if (length(source_path) != 1L) {
      stop("Frozen engine role did not resolve uniquely: ", role_value,
           call. = FALSE)
    }
    sys.source(source_path, envir = .GlobalEnv)
  }
  required_engine_fns <- c("srb_verify_cache_bundle",
                           "srb_run_provenance_preflight",
                           "srb_benchmark_within_cell", "srb_expected_outer_k",
                           "srb_sanitize_character_fields")
  if (!all(vapply(required_engine_fns, exists, logical(1L), mode = "function",
                  envir = .GlobalEnv, inherits = FALSE))) {
    stop("Frozen benchmark sources did not expose their required API.",
         call. = FALSE)
  }
  srb_verify_cache_bundle(cache, cache_manifest_path, cache_digest_path)
  cache_manifest <- data.table::fread(cache_manifest_path)
  if (!"package_commit" %in% names(cache_manifest) ||
      data.table::uniqueN(cache_manifest$package_commit) != 1L ||
      !identical(unique(cache_manifest$package_commit),
                 opt$expected_cache_commit)) {
    stop("Cohort cache manifest does not match the expected cache-building ",
         "commit.", call. = FALSE)
  }
  cohorts <- readRDS(cache)
  if (!is.list(cohorts) || length(cohorts) != 34L ||
      anyNA(names(cohorts)) || any(!nzchar(names(cohorts))) ||
      anyDuplicated(names(cohorts))) {
    stop("Cohort cache is not the exact 34-unit public benchmark corpus.",
         call. = FALSE)
  }
  units <- .ranktree_units(cohorts)
  expected_k <- as.integer(srb_expected_outer_k(units$unit_id))
  if (any(units$outer_k != expected_k) || sum(units$primary_unit) != 33L) {
    stop("Cohort cache violates frozen outer-fold or primary-unit contracts.",
         call. = FALSE)
  }

  base_roster <- OmicSelector::singlesample_method_roster()
  anchor <- base_roster[match(.ranktree_methods()[1:5], base_roster$method_id), ]
  external <- OmicSelector::singlesample_external_competitor_roster()
  roster <- data.table::rbindlist(list(anchor, external), fill = TRUE)
  if (nrow(roster) != 7L || anyNA(roster$method_id) ||
      !identical(roster$method_id, .ranktree_methods()) ||
      any(roster$estimand != "within")) {
    stop("Installed package does not expose the frozen seven-method ",
         "sensitivity registry.", call. = FALSE)
  }

  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) stop("Output bundle already exists: ", output_dir)
  stage <- tempfile(".ranktree-run-", tmpdir = output_parent)
  dir.create(stage)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  provenance_log <- file.path(stage, "provenance_preflight.log")
  provenance <- .ranktree_run_bound_provenance(
    cohorts, paper_root, provenance_script, provenance_manifest,
    provenance_union_inventory, opt$expected_provenance_script_sha256,
    opt$expected_provenance_manifest_sha256,
    opt$expected_provenance_union_inventory_sha256, provenance_log,
    srb_run_provenance_preflight
  )
  if (!file.exists(provenance$resolution_path)) {
    stop("Provenance preflight did not emit its resolution table.",
         call. = FALSE)
  }

  pins <- list(
    version = package_version, commit = code_start$commit,
    cache_commit = opt$expected_cache_commit,
    cache_sha256 = opt$expected_cache_sha256,
    cache_manifest_sha256 = opt$expected_cache_manifest_sha256,
    cache_digest_sha256 = opt$expected_cache_digest_sha256,
    engine_manifest_sha256 = opt$expected_engine_manifest_sha256,
    environment_manifest_sha256 = opt$expected_environment_manifest_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_sha256 =
      opt$expected_analysis_amendment_sha256,
    analysis_amendment2_sha256 =
      opt$expected_analysis_amendment2_sha256,
    analysis_amendment3_sha256 =
      opt$expected_analysis_amendment3_sha256,
    provenance_script_sha256 = opt$expected_provenance_script_sha256,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256,
    analysis_code_id = analysis_code_id
  )
  checkpoint_contract <- .ranktree_checkpoint_contract(
    output_dir, opt$methods, units, pins, runner_sha, opt$expected_r_version)
  checkpoint <- .ranktree_open_checkpoints(output_dir, checkpoint_contract)
  cell_rows <- prediction_rows <- list()
  index <- 0L
  for (method in opt$methods) {
    method_roster <- as.data.frame(roster)
    for (seed in .ranktree_seeds()) {
      seed_dir <- file.path(checkpoint, sprintf("seed_%03d", seed))
      dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)
      for (cohort_id in names(cohorts)) {
        message(sprintf("[%d] %s / %s", seed, method, cohort_id))
        outer_k <- as.integer(srb_expected_outer_k(cohort_id))
        result <- .ranktree_call_benchmark_cell(
          srb_benchmark_within_cell, method, cohort_id,
          cohorts[[cohort_id]], method_roster, seed_dir, seed, outer_k
        )
        result$outer_k <- outer_k
        cell <- .ranktree_cell_row(result, seed, pins)
        if (!isTRUE(result$eligible)) {
          stop(
            "Rank-tree sensitivity cell failed the Amendment 2 complete-finite ",
            "support contract: ", method, "/", cohort_id, "/seed ", seed,
            ": ", result$ineligible_reason, call. = FALSE
          )
        }
        .ranktree_validate_cell_result(cell)
        index <- index + 1L
        cell_rows[[index]] <- cell
        if (isTRUE(result$eligible)) {
          if (is.null(result$fold_scores)) {
            stop("Eligible cell lacks held-out fold scores.", call. = FALSE)
          }
          prediction <- data.table::as.data.table(result$fold_scores)
          prediction[, `:=`(
            method = method, cohort = cohort_id, seed = as.integer(seed),
            package_version = package_version,
            package_commit = code_start$commit,
            cache_package_commit = opt$expected_cache_commit,
            runtime_image_sha256 = opt$expected_runtime_image_sha256,
            runtime_attestation_manifest_sha256 =
              opt$expected_runtime_attestation_manifest_sha256,
            container_launcher_sha256 =
              opt$expected_container_launcher_sha256,
            r_version = opt$expected_r_version,
            analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
            analysis_amendment_sha256 =
              opt$expected_analysis_amendment_sha256,
            analysis_amendment2_sha256 =
              opt$expected_analysis_amendment2_sha256,
            analysis_amendment3_sha256 =
              opt$expected_analysis_amendment3_sha256,
            provenance_manifest_sha256 =
              opt$expected_provenance_manifest_sha256,
            provenance_union_inventory_sha256 =
              opt$expected_provenance_union_inventory_sha256,
            analysis_code_id = analysis_code_id
          )]
          prediction_rows[[length(prediction_rows) + 1L]] <- prediction
        }
      }
    }
  }
  cells <- data.table::rbindlist(cell_rows, fill = TRUE)
  predictions <- data.table::rbindlist(prediction_rows, fill = TRUE)
  cells <- srb_sanitize_character_fields(cells)
  predictions <- srb_sanitize_character_fields(predictions)
  predictions <- .ranktree_add_split_ids(predictions)
  data.table::setorder(cells, method, cohort, seed)
  data.table::setorder(predictions, method, cohort, seed, fold, sample_id)
  .ranktree_validate_runtime_bundle(
    cells, predictions, units, opt$methods,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version)

  bundle <- list(
    cells = cells, predictions = predictions, units = units,
    methods = opt$methods, method_registry = roster[method_id %in% opt$methods],
    seeds = .ranktree_seeds(), cohorts = units$unit_id,
    package_version = package_version, package_commit = code_start$commit,
    cache_package_commit = opt$expected_cache_commit,
    cache_sha256 = opt$expected_cache_sha256,
    cache_manifest_sha256 = opt$expected_cache_manifest_sha256,
    cache_digest_sha256 = opt$expected_cache_digest_sha256,
    runner_script_sha256 = runner_sha,
    engine_manifest_sha256 = opt$expected_engine_manifest_sha256,
    environment_manifest_sha256 = opt$expected_environment_manifest_sha256,
    runtime_image_path = runtime_image,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_path = launcher_path,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    analysis_plan_path = analysis_plan,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_path = analysis_amendment,
    analysis_amendment_sha256 = opt$expected_analysis_amendment_sha256,
    analysis_amendment2_path = analysis_amendment2,
    analysis_amendment2_sha256 = opt$expected_analysis_amendment2_sha256,
    analysis_amendment3_path = analysis_amendment3,
    analysis_amendment3_sha256 = opt$expected_analysis_amendment3_sha256,
    provenance_script_sha256 = opt$expected_provenance_script_sha256,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_path = provenance_union_inventory,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256,
    analysis_code_id = analysis_code_id,
    r_version = opt$expected_r_version,
    ranktree_dependency_closure = environment_manifest[, .(package, version,
                                                           source_sha256)],
    provenance_exit_status = provenance$raw_status,
    provenance_overlap_pairs = provenance$n_overlap_pairs
  )
  saveRDS(bundle, file.path(stage, "within_bundle.rds"), version = 3L)
  data.table::fwrite(cells, file.path(stage, "within_cells.tsv"), sep = "\t",
                     quote = FALSE, na = "")
  data.table::fwrite(predictions,
                     file.path(stage, "within_predictions.tsv.gz"), sep = "\t",
                     quote = FALSE, na = "", compress = "gzip")
  data.table::fwrite(units, file.path(stage, "units.tsv"), sep = "\t",
                     quote = FALSE, na = "")
  data.table::fwrite(roster[method_id %in% opt$methods],
                     file.path(stage, "methods.tsv"), sep = "\t",
                     quote = FALSE, na = "")
  .ranktree_copy(opt$engine_manifest,
                 file.path(stage, "engine_manifest.tsv"), "engine manifest")
  .ranktree_copy(opt$environment_manifest,
                 file.path(stage, "environment_manifest.tsv"),
                 "environment manifest")
  .ranktree_copy(runtime_attestation$path,
                 file.path(stage, "runtime_attestation.tsv"),
                 "runtime closure attestation")
  .ranktree_copy(analysis_plan, file.path(stage, "analysis_plan.md"),
                 "approved analysis plan")
  .ranktree_copy(analysis_amendment, file.path(stage, "analysis_amendment.md"),
                 "approved analysis amendment 1")
  .ranktree_copy(analysis_amendment2,
                 file.path(stage, "analysis_amendment2.md"),
                 "approved analysis amendment 2")
  .ranktree_copy(analysis_amendment3,
                 file.path(stage, "analysis_amendment3.md"),
                 "approved analysis amendment 3")
  resolution_target <- file.path(stage, "provenance_resolution.tsv")
  if (!identical(normalizePath(provenance$resolution_path, mustWork = TRUE),
                 normalizePath(resolution_target, mustWork = FALSE))) {
    .ranktree_copy(provenance$resolution_path, resolution_target,
                   "provenance resolution")
    unlink(provenance$resolution_path)
  }
  writeLines(c(
    "# External rank-tree sensitivity method bundle", "",
    paste0("Methods: ", paste(opt$methods, collapse = ", ")),
    paste0("Analysis units: ", nrow(units), " (",
           sum(units$primary_unit), " primary)"),
    paste0("Deposited profiles represented: ", sum(units$n_deposited_profiles)),
    paste0("Deduplicated biological/provenance groups: ",
           sum(units$n_unique_groups)),
    paste0("Singularity runtime image SHA-256: ",
           opt$expected_runtime_image_sha256),
    paste0("Exact in-process R runtime/library attestation SHA-256: ",
           opt$expected_runtime_attestation_manifest_sha256),
    paste0("Package container launcher SHA-256: ",
           opt$expected_container_launcher_sha256),
    paste0("Canonical provenance union inventory SHA-256: ",
           opt$expected_provenance_union_inventory_sha256),
    paste0("Approved analysis plan SHA-256: ",
           opt$expected_analysis_plan_sha256),
    paste0("Approved analysis amendment 1 SHA-256: ",
           opt$expected_analysis_amendment_sha256),
    paste0("Approved analysis amendment 2 SHA-256: ",
           opt$expected_analysis_amendment2_sha256),
    paste0("Approved analysis amendment 3 SHA-256: ",
           opt$expected_analysis_amendment3_sha256),
    "Rank forest contract: unweighted preliminary screen; inverse-frequency ",
    "pair-forest case weights; numeric y=1 (`case`) is the positive target.",
    "Performance is computed as group-collapsed held-out AUC: one mean score ",
    "per provenance group within each exact test fold.",
    "", "This is an unreviewed execution bundle. It is not reportable until ",
    "assembly, independent results review, plausibility assessment and promotion."
  ), file.path(stage, "report.md"))

  # Recheck all mutable external inputs and package state before promotion.
  .ranktree_validate_runtime_attestation(
    runtime_attestation$path,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_r_version, opt$expected_runtime_image_sha256,
    opt$expected_container_launcher_sha256
  )
  end_pins <- c(
    .ranktree_sha(cache), .ranktree_sha(cache_manifest_path),
    .ranktree_sha(cache_digest_path), .ranktree_sha(opt$engine_manifest),
    .ranktree_sha(analysis_plan),
    .ranktree_sha(analysis_amendment),
    .ranktree_sha(analysis_amendment2),
    .ranktree_sha(analysis_amendment3),
    .ranktree_sha(opt$environment_manifest), .ranktree_sha(runtime_image),
    .ranktree_sha(runtime_attestation$path), .ranktree_sha(launcher_path),
    .ranktree_sha(provenance_script),
    .ranktree_sha(provenance_manifest),
    .ranktree_sha(provenance_union_inventory)
  )
  expected_end <- c(
    opt$expected_cache_sha256, opt$expected_cache_manifest_sha256,
    opt$expected_cache_digest_sha256, opt$expected_engine_manifest_sha256,
    opt$expected_analysis_plan_sha256,
    opt$expected_analysis_amendment_sha256,
    opt$expected_analysis_amendment2_sha256,
    opt$expected_analysis_amendment3_sha256,
    opt$expected_environment_manifest_sha256,
    opt$expected_runtime_image_sha256,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_container_launcher_sha256,
    opt$expected_provenance_script_sha256,
    opt$expected_provenance_manifest_sha256,
    opt$expected_provenance_union_inventory_sha256
  )
  code_end <- .ranktree_git_provenance(opt$package_root)
  if (!identical(end_pins, expected_end) || !identical(code_end, code_start) ||
      !identical(.ranktree_sha(script_path), runner_sha)) {
    stop("A pinned input or OmicSelector code provenance changed during ",
         "execution.", call. = FALSE)
  }
  artifacts <- list.files(stage, full.names = TRUE)
  manifest <- data.table::data.table(
    file = basename(artifacts), bytes = file.info(artifacts)$size,
    sha256 = vapply(artifacts, .ranktree_sha, character(1L)),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = package_version, package_commit = code_start$commit,
    runner_script_sha256 = runner_sha, analysis_code_id = analysis_code_id,
    engine_manifest_sha256 = opt$expected_engine_manifest_sha256,
    environment_manifest_sha256 = opt$expected_environment_manifest_sha256,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_sha256 = opt$expected_analysis_amendment_sha256,
    analysis_amendment2_sha256 = opt$expected_analysis_amendment2_sha256,
    analysis_amendment3_sha256 = opt$expected_analysis_amendment3_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256
  )
  data.table::fwrite(manifest, file.path(stage, "output_manifest.tsv"),
                     sep = "\t", quote = FALSE)
  if (!file.rename(stage, output_dir)) {
    stop("Could not atomically promote rank-tree method bundle.", call. = FALSE)
  }
  .ranktree_remove_completed_checkpoints(output_dir, checkpoint_contract)
  cat("Created external-competitor method bundle:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) {
  trailing <- commandArgs(trailingOnly = TRUE)
  capture_index <- match("--write-runtime-attestation", trailing)
  if (!is.na(capture_index)) {
    if (capture_index == length(trailing) ||
        !nzchar(trailing[[capture_index + 1L]])) {
      stop("--write-runtime-attestation requires a target path.",
           call. = FALSE)
    }
    .ranktree_write_runtime_attestation(trailing[[capture_index + 1L]])
  } else {
    .ranktree_runner_main(trailing)
  }
}
