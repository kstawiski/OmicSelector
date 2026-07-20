# Shared fail-closed runtime receipt capture and validation for paper 3.
# Consumers source this file, validate once before analysis, retain the returned
# guard, and validate again with that guard immediately before promotion.

.PAPER3RT_HELPER_PATH <- local({
  candidates <- unlist(lapply(sys.frames(), function(frame) {
    candidate <- frame$ofile
    if (is.null(candidate)) character() else as.character(candidate)
  }), use.names = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    normalizePath(tail(candidates, 1L), mustWork = TRUE)
  } else {
    NA_character_
  }
})

.paper3rt_sha256 <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

.paper3rt_required_thread_env <- c(
  "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
  "GOTO_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
  "RCPP_PARALLEL_NUM_THREADS", "R_DATATABLE_NUM_THREADS"
)

.paper3rt_expected_rngkind <- c(
  "Mersenne-Twister", "Inversion", "Rejection"
)

paper3_runtime_required_thread_env <- function() {
  stats::setNames(rep("1", length(.paper3rt_required_thread_env)),
                  .paper3rt_required_thread_env)
}

.paper3rt_assert_hex <- function(value, width, label) {
  if (length(value) != 1L || is.na(value) ||
      !grepl(sprintf("^[0-9a-f]{%d}$", width), value)) {
    stop(label, " must be an exact lowercase hexadecimal pin of ", width,
         " characters.", call. = FALSE)
  }
  invisible(value)
}

.paper3rt_assert_text <- function(value, label, allow_empty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      (!allow_empty && !nzchar(value)) || grepl("[\r\n\t]", value)) {
    stop(label, " must be one control-character-free string.", call. = FALSE)
  }
  invisible(value)
}

.paper3rt_read_tsv <- function(path, label) {
  if (!file.exists(path) || dir.exists(path) || nzchar(Sys.readlink(path))) {
    stop(label, " is missing, symbolic, or not a regular file.", call. = FALSE)
  }
  data.table::fread(path, colClasses = "character", na.strings = NULL,
                    showProgress = FALSE)
}

.paper3rt_assert_canonical_file <- function(path, label,
                                            reject_supplied_symlink = FALSE) {
  .paper3rt_assert_text(path, label)
  if (!grepl("^/", path) || !file.exists(path) || dir.exists(path)) {
    stop(label, " must be an existing absolute regular-file path.",
         call. = FALSE)
  }
  if (reject_supplied_symlink && nzchar(Sys.readlink(path))) {
    stop(label, " must not be a symbolic link.", call. = FALSE)
  }
  canonical <- normalizePath(path, mustWork = TRUE)
  if (!grepl("^/", canonical) || !file.exists(canonical) ||
      dir.exists(canonical)) {
    stop(label, " could not be resolved to a real absolute file.",
         call. = FALSE)
  }
  canonical
}

.paper3rt_resolve_executable <- function(name) {
  resolved <- unname(Sys.which(name))
  if (length(resolved) != 1L || !nzchar(resolved)) {
    stop("Required executable is unavailable: ", name, call. = FALSE)
  }
  path <- .paper3rt_assert_canonical_file(
    normalizePath(resolved, mustWork = TRUE), paste0("Executable ", name)
  )
  data.table::data.table(
    asset = paste0("executable_", name), path = path,
    bytes = as.character(unname(file.info(path)$size)),
    sha256 = .paper3rt_sha256(path)
  )
}

.paper3rt_package_tree_sha256 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  files <- sort(list.files(
    root, recursive = TRUE, all.files = TRUE, full.names = TRUE,
    include.dirs = FALSE, no.. = TRUE
  ))
  if (!length(files)) {
    stop("Installed package tree is empty: ", root, call. = FALSE)
  }
  info <- file.info(files)
  if (anyNA(info$isdir) || any(info$isdir) || anyNA(info$size) ||
      any(nzchar(Sys.readlink(files)))) {
    stop("Installed package tree contains unreadable, non-regular, or ",
         "symbolic entries: ", root, call. = FALSE)
  }
  canonical <- normalizePath(files, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (any(!startsWith(canonical, prefix))) {
    stop("Installed package tree escapes its canonical root: ", root,
         call. = FALSE)
  }
  relative <- substring(canonical, nchar(prefix) + 1L)
  if (any(!nzchar(relative)) || anyDuplicated(relative) ||
      any(grepl("[\r\n\t]", relative))) {
    stop("Installed package tree has unsafe or duplicate relative paths.",
         call. = FALSE)
  }
  rows <- paste(
    relative, as.character(unname(info$size)),
    vapply(canonical, .paper3rt_sha256, character(1L)), sep = "\t"
  )
  digest::digest(paste(rows, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

paper3_runtime_installed_package_tree_sha256 <- function(
    package = "OmicSelector") {
  root <- find.package(package, lib.loc = .libPaths(), quiet = TRUE)
  if (!length(root) || !nzchar(root)) {
    stop("Required installed package is unavailable: ", package,
         call. = FALSE)
  }
  .paper3rt_package_tree_sha256(root)
}

.paper3rt_active_installed_db <- function() {
  raw <- utils::installed.packages(noCache = TRUE)
  if (!nrow(raw) || !all(c("Package", "Version", "LibPath") %in%
                         colnames(raw))) {
    stop("Could not read the installed-package database.", call. = FALSE)
  }
  library_order <- unique(vapply(
    .libPaths(), normalizePath, character(1L), mustWork = TRUE
  ))
  packages <- unique(as.character(raw[, "Package"]))
  selected <- vapply(packages, function(package) {
    candidates <- which(raw[, "Package"] == package)
    candidate_libraries <- vapply(
      raw[candidates, "LibPath"], normalizePath, character(1L),
      mustWork = TRUE
    )
    # Resolve from installed.packages() itself, in active library precedence.
    # This excludes pkgload source shadows and does not make the complete
    # receipt depend on find.package() resolving unrelated installed records.
    priority <- match(candidate_libraries, library_order)
    if (anyNA(priority)) return(NA_integer_)
    matches <- candidates[priority == min(priority)]
    if (length(matches) == 1L) matches[[1L]] else NA_integer_
  }, integer(1L))
  if (anyNA(selected)) {
    stop("Could not resolve exactly one active library record per installed ",
         "package.", call. = FALSE)
  }
  out <- raw[selected, , drop = FALSE]
  rownames(out) <- out[, "Package"]
  out
}

.paper3rt_capture_package_closure <- function(
    roots = c("OmicSelector", "data.table", "digest", "metafor")) {
  if (anyDuplicated(roots) || any(!nzchar(roots))) {
    stop("Runtime package roots must be unique non-empty names.",
         call. = FALSE)
  }
  installed <- .paper3rt_active_installed_db()
  missing_roots <- setdiff(roots, rownames(installed))
  if (length(missing_roots)) {
    stop("Runtime package closure is missing root package(s): ",
         paste(missing_roots, collapse = ", "), call. = FALSE)
  }
  dependencies <- unique(unlist(tools::package_dependencies(
    roots, db = installed, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  ), use.names = FALSE))
  dependencies <- dependencies[!is.na(dependencies) & nzchar(dependencies)]
  packages <- sort(setdiff(unique(c(roots, dependencies)), "R"))
  missing <- setdiff(packages, rownames(installed))
  if (length(missing)) {
    stop("Runtime recursive package closure is incomplete: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  roots_path <- vapply(packages, function(package) {
    normalizePath(file.path(installed[package, "LibPath"], package),
                  mustWork = TRUE)
  }, character(1L))
  out <- data.table::data.table(
    package = packages,
    version = as.character(installed[packages, "Version"]),
    lib_path = vapply(dirname(roots_path), normalizePath, character(1L),
                      mustWork = TRUE),
    root_path = roots_path,
    tree_sha256 = vapply(roots_path, .paper3rt_package_tree_sha256,
                         character(1L)),
    is_root = ifelse(packages %in% roots, "TRUE", "FALSE")
  )
  data.table::setorder(out, package)
  out
}

.paper3rt_package_closure_sha256 <- function(closure) {
  columns <- c("package", "version", "lib_path", "root_path",
               "tree_sha256", "is_root")
  closure <- data.table::as.data.table(closure)
  if (!identical(names(closure), columns)) {
    stop("Package closure does not have the exact schema.", call. = FALSE)
  }
  digest::digest(
    paste(do.call(paste, c(closure[, ..columns], sep = "\t")),
          collapse = "\n"),
    algo = "sha256", serialize = FALSE
  )
}

.paper3rt_capture_facts <- function() {
  fact_rows <- function(section, key, value) {
    data.table::as.data.table(list(
      section = section, key = key, value = value
    ))
  }
  thread_values <- Sys.getenv(.paper3rt_required_thread_env, unset = NA_character_)
  if (anyNA(thread_values) || any(thread_values != "1")) {
    bad <- .paper3rt_required_thread_env[
      is.na(thread_values) | thread_values != "1"
    ]
    stop("Required single-thread environment variable(s) are not exactly 1: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  dt_threads <- data.table::getDTthreads(verbose = FALSE)
  if (length(dt_threads) != 1L || is.na(dt_threads) || dt_threads != 1L) {
    stop("data.table must use exactly one thread.", call. = FALSE)
  }
  rng <- RNGkind()
  if (!identical(unname(rng), .paper3rt_expected_rngkind)) {
    stop("RNGkind must be exactly Mersenne-Twister/Inversion/Rejection.",
         call. = FALSE)
  }
  external <- extSoftVersion()
  if (!length(external) || is.null(names(external)) ||
      anyNA(names(external)) || any(!nzchar(names(external))) ||
      anyDuplicated(names(external)) || !"BLAS" %in% names(external)) {
    stop("extSoftVersion() lacks a unique named BLAS-inclusive closure.",
         call. = FALSE)
  }
  lapack <- La_library()
  .paper3rt_assert_text(as.character(lapack), "LAPACK runtime")
  out <- data.table::rbindlist(list(
    fact_rows(
      "r", c("version", "platform"),
      c(R.version.string, R.version$platform)
    ),
    fact_rows(
      "rngkind", c("kind", "normal_kind", "sample_kind"), rng
    ),
    fact_rows(
      "numeric_runtime", c("BLAS", "LAPACK"),
      c(unname(external[["BLAS"]]), as.character(lapack))
    ),
    fact_rows(
      "extSoftVersion", names(external), unname(as.character(external))
    ),
    fact_rows(
      "thread_env", .paper3rt_required_thread_env, unname(thread_values)
    ),
    fact_rows(
      "data.table", "threads", as.character(dt_threads)
    )
  ), use.names = TRUE)
  if (anyNA(out) || anyDuplicated(out[, .(section, key)]) ||
      any(grepl("[\r\n\t]", out$section)) ||
      any(grepl("[\r\n\t]", out$key)) ||
      any(grepl("[\r\n\t]", out$value))) {
    stop("Runtime facts contain missing, duplicate, or unsafe values.",
         call. = FALSE)
  }
  data.table::setorder(out, section, key)
  out
}

.paper3rt_capture_assets <- function(runtime_image_path) {
  image <- .paper3rt_assert_canonical_file(
    runtime_image_path, "Runtime image", reject_supplied_symlink = TRUE
  )
  if (!identical(image, runtime_image_path)) {
    stop("Runtime image path must already be its canonical real path.",
         call. = FALSE)
  }
  out <- data.table::rbindlist(list(
    data.table::data.table(
      asset = "runtime_image", path = image,
      bytes = as.character(unname(file.info(image)$size)),
      sha256 = .paper3rt_sha256(image)
    ),
    .paper3rt_resolve_executable("Rscript"),
    .paper3rt_resolve_executable("gzip"),
    .paper3rt_resolve_executable("git")
  ), use.names = TRUE)
  if (nrow(out) != 4L || anyDuplicated(out$asset) || anyNA(out) ||
      !setequal(out$asset, c("runtime_image", "executable_Rscript",
                            "executable_gzip", "executable_git")) ||
      any(!grepl("^[0-9]+$", out$bytes)) ||
      any(!grepl("^[0-9a-f]{64}$", out$sha256)) ||
      any(!grepl("^/", out$path)) ||
      any(out$path != vapply(out$path, normalizePath, character(1L),
                             mustWork = TRUE))) {
    stop("Runtime asset capture is incomplete or non-canonical.",
         call. = FALSE)
  }
  data.table::setorder(out, asset)
  out
}

.paper3rt_capture <- function(runtime_image_path) {
  list(
    facts = .paper3rt_capture_facts(),
    assets = .paper3rt_capture_assets(runtime_image_path),
    packages = .paper3rt_capture_package_closure()
  )
}

.paper3rt_capture_digest <- function(capture) {
  table_payload <- function(name, x) {
    x <- data.table::as.data.table(x)
    c(paste0("[", name, "]"), paste(names(x), collapse = "\t"),
      do.call(paste, c(lapply(x, as.character), sep = "\t")))
  }
  payload <- c(
    table_payload("facts", capture$facts),
    table_payload("assets", capture$assets),
    table_payload("packages", capture$packages)
  )
  digest::digest(paste(payload, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

.paper3rt_assert_table_identical <- function(observed, expected, label) {
  observed <- data.table::as.data.table(observed)
  expected <- data.table::as.data.table(expected)
  if (!identical(names(observed), names(expected)) ||
      !identical(nrow(observed), nrow(expected)) ||
      !all(vapply(names(expected), function(column) {
        identical(as.character(observed[[column]]),
                  as.character(expected[[column]]))
      }, logical(1L)))) {
    stop(label, " differs from the exact frozen receipt.", call. = FALSE)
  }
  invisible(TRUE)
}

.paper3rt_current_source_hashes <- function() {
  helper <- .PAPER3RT_HELPER_PATH
  if (length(helper) != 1L || is.na(helper) || !file.exists(helper)) {
    stop("Could not resolve the sourced runtime-receipt helper bytes.",
         call. = FALSE)
  }
  helper <- normalizePath(helper, mustWork = TRUE)
  producer <- file.path(dirname(helper), "paper3_build_runtime_receipt.R")
  if (!file.exists(producer)) {
    stop("Could not resolve the sibling runtime-receipt producer bytes.",
         call. = FALSE)
  }
  producer <- normalizePath(producer, mustWork = TRUE)
  list(
    helper_path = helper, helper_sha256 = .paper3rt_sha256(helper),
    producer_path = producer, producer_sha256 = .paper3rt_sha256(producer)
  )
}

.paper3rt_validate_manifest <- function(receipt_dir,
                                        expected_manifest_sha256) {
  .paper3rt_assert_hex(expected_manifest_sha256, 64L,
                       "Runtime receipt manifest SHA-256")
  receipt_dir <- normalizePath(receipt_dir, mustWork = TRUE)
  if (!dir.exists(receipt_dir)) {
    stop("Runtime receipt path is not a directory.", call. = FALSE)
  }
  manifest_path <- file.path(receipt_dir, "manifest.tsv")
  if (!file.exists(manifest_path) ||
      !identical(.paper3rt_sha256(manifest_path),
                 expected_manifest_sha256)) {
    stop("Runtime receipt manifest does not match its exact SHA-256 pin.",
         call. = FALSE)
  }
  manifest <- .paper3rt_read_tsv(manifest_path, "Runtime receipt manifest")
  columns <- c(
    "file", "role", "bytes", "sha256", "package_version",
    "package_commit", "installed_omicselector_tree_sha256",
    "package_closure_sha256", "runtime_image_path",
    "runtime_image_sha256", "producer_script_sha256",
    "helper_script_sha256"
  )
  artifacts <- c("runtime_facts.tsv", "runtime_assets.tsv",
                 "package_closure.tsv", "producer_provenance.tsv")
  roles <- c("runtime_facts", "runtime_assets", "package_closure",
             "producer_provenance")
  if (!identical(names(manifest), columns) || nrow(manifest) != 4L ||
      anyNA(manifest) || anyDuplicated(manifest$file) ||
      !setequal(manifest$file, artifacts) ||
      any(!grepl("^[0-9]+$", manifest$bytes)) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Runtime receipt manifest has an invalid exact-inventory schema.",
         call. = FALSE)
  }
  present <- list.files(receipt_dir, recursive = FALSE, all.files = TRUE,
                        include.dirs = TRUE, no.. = TRUE)
  if (!setequal(present, c(artifacts, "manifest.tsv")) ||
      any(dir.exists(file.path(receipt_dir, present)))) {
    stop("Runtime receipt contains an unmanifested artifact or directory.",
         call. = FALSE)
  }
  paths <- file.path(receipt_dir, manifest$file)
  if (any(!file.exists(paths)) || any(nzchar(Sys.readlink(paths))) ||
      any(as.numeric(manifest$bytes) != unname(file.info(paths)$size)) ||
      any(manifest$sha256 != vapply(paths, .paper3rt_sha256, character(1L)))) {
    stop("Runtime receipt artifact bytes differ from its exact manifest.",
         call. = FALSE)
  }
  manifest <- manifest[match(artifacts, file)]
  if (!identical(manifest$role, roles)) {
    stop("Runtime receipt manifest roles do not match its exact inventory.",
         call. = FALSE)
  }
  list(path = manifest_path, sha256 = expected_manifest_sha256,
       table = manifest, artifacts = artifacts, receipt_dir = receipt_dir)
}

paper3_validate_runtime_receipt <- function(
    receipt_dir, expected_manifest_sha256, expected_package_version,
    expected_package_commit, expected_installed_package_tree_sha256,
    expected_runtime_image_path, expected_runtime_image_sha256,
    phase = c("start", "end"), start_guard = NULL) {
  phase <- match.arg(phase)
  .paper3rt_assert_text(expected_package_version,
                        "Expected OmicSelector package version")
  .paper3rt_assert_hex(expected_package_commit, 40L,
                       "Expected package commit")
  .paper3rt_assert_hex(expected_installed_package_tree_sha256, 64L,
                       "Expected installed OmicSelector tree SHA-256")
  .paper3rt_assert_hex(expected_runtime_image_sha256, 64L,
                       "Expected runtime image SHA-256")
  expected_runtime_image_path <- .paper3rt_assert_canonical_file(
    expected_runtime_image_path, "Expected runtime image",
    reject_supplied_symlink = TRUE
  )
  manifest <- .paper3rt_validate_manifest(
    receipt_dir, expected_manifest_sha256
  )
  facts <- .paper3rt_read_tsv(
    file.path(manifest$receipt_dir, "runtime_facts.tsv"), "Runtime facts"
  )
  assets <- .paper3rt_read_tsv(
    file.path(manifest$receipt_dir, "runtime_assets.tsv"), "Runtime assets"
  )
  packages <- .paper3rt_read_tsv(
    file.path(manifest$receipt_dir, "package_closure.tsv"), "Package closure"
  )
  provenance <- .paper3rt_read_tsv(
    file.path(manifest$receipt_dir, "producer_provenance.tsv"),
    "Producer provenance"
  )
  if (!identical(names(facts), c("section", "key", "value")) ||
      !identical(names(assets), c("asset", "path", "bytes", "sha256")) ||
      !identical(names(packages), c(
        "package", "version", "lib_path", "root_path", "tree_sha256",
        "is_root"
      )) || !identical(names(provenance), c(
        "schema_version", "package_version", "package_commit", "git_clean",
        "producer_script_sha256", "helper_script_sha256",
        "installed_omicselector_tree_sha256", "package_closure_sha256",
        "runtime_image_path", "runtime_image_sha256"
      )) || nrow(provenance) != 1L || anyNA(provenance)) {
    stop("Runtime receipt artifact schema is not exact.", call. = FALSE)
  }
  current_sources <- .paper3rt_current_source_hashes()
  closure_sha <- .paper3rt_package_closure_sha256(packages)
  expected_common <- list(
    package_version = expected_package_version,
    package_commit = expected_package_commit,
    installed_omicselector_tree_sha256 =
      expected_installed_package_tree_sha256,
    package_closure_sha256 = closure_sha,
    runtime_image_path = expected_runtime_image_path,
    runtime_image_sha256 = expected_runtime_image_sha256,
    producer_script_sha256 = current_sources$producer_sha256,
    helper_script_sha256 = current_sources$helper_sha256
  )
  for (field in names(expected_common)) {
    value <- expected_common[[field]]
    if (!identical(unique(manifest$table[[field]]), value) ||
        !identical(provenance[[field]][[1L]], value)) {
      stop("Runtime receipt ", field, " differs from its exact expected pin.",
           call. = FALSE)
    }
  }
  if (!identical(provenance$schema_version[[1L]],
                 "OmicSelector-paper3-runtime-receipt-v1") ||
      !identical(provenance$git_clean[[1L]], "TRUE")) {
    stop("Runtime receipt is not bound to a clean supported producer schema.",
         call. = FALSE)
  }
  current <- .paper3rt_capture(expected_runtime_image_path)
  .paper3rt_assert_table_identical(current$facts, facts, "Runtime facts")
  .paper3rt_assert_table_identical(current$assets, assets, "Runtime assets")
  .paper3rt_assert_table_identical(current$packages, packages,
                                    "Installed package closure")
  image <- current$assets[asset == "runtime_image"]
  omicselector <- current$packages[package == "OmicSelector"]
  roots <- current$packages[is_root == "TRUE", package]
  if (nrow(image) != 1L || image$path[[1L]] != expected_runtime_image_path ||
      image$sha256[[1L]] != expected_runtime_image_sha256 ||
      nrow(omicselector) != 1L ||
      omicselector$version[[1L]] != expected_package_version ||
      omicselector$tree_sha256[[1L]] !=
        expected_installed_package_tree_sha256 ||
      !setequal(roots, c("OmicSelector", "data.table", "digest", "metafor"))) {
    stop("Current runtime image or OmicSelector dependency roots differ from ",
         "their exact receipt pins.", call. = FALSE)
  }
  binding <- list(
    receipt_dir = manifest$receipt_dir,
    manifest_sha256 = expected_manifest_sha256,
    package_version = expected_package_version,
    package_commit = expected_package_commit,
    installed_omicselector_tree_sha256 =
      expected_installed_package_tree_sha256,
    runtime_image_path = expected_runtime_image_path,
    runtime_image_sha256 = expected_runtime_image_sha256
  )
  binding_sha <- digest::digest(serialize(binding, NULL, version = 3L),
                                algo = "sha256")
  current_sha <- .paper3rt_capture_digest(current)
  if (identical(phase, "start")) {
    if (!is.null(start_guard)) {
      stop("start_guard must be NULL during start validation.", call. = FALSE)
    }
    return(structure(
      c(binding, list(binding_sha256 = binding_sha,
                      current_runtime_sha256 = current_sha)),
      class = "paper3_runtime_receipt_guard"
    ))
  }
  if (!inherits(start_guard, "paper3_runtime_receipt_guard") ||
      !identical(start_guard$binding_sha256, binding_sha) ||
      !identical(start_guard$current_runtime_sha256, current_sha)) {
    stop("End validation does not match the exact start runtime guard.",
         call. = FALSE)
  }
  invisible(start_guard)
}
