#!/usr/bin/env Rscript

# Build an immutable paper-3 runtime receipt. This captures runtime facts and
# external assets only; it never fits a model or mutates an analysis input.

.PAPER3RT_PRODUCER_PATH <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    candidate <- sub("^--file=", "", file_arg[[1L]])
    if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
  }
  candidates <- unlist(lapply(sys.frames(), function(frame) {
    candidate <- frame$ofile
    if (is.null(candidate)) character() else as.character(candidate)
  }), use.names = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    return(normalizePath(tail(candidates, 1L), mustWork = TRUE))
  }
  NA_character_
})

if (length(.PAPER3RT_PRODUCER_PATH) != 1L ||
    is.na(.PAPER3RT_PRODUCER_PATH)) {
  stop("Could not resolve the runtime-receipt producer path.", call. = FALSE)
}
.paper3rt_helper_source <- file.path(
  dirname(.PAPER3RT_PRODUCER_PATH), "paper3_runtime_receipt_common.R"
)
if (!file.exists(.paper3rt_helper_source)) {
  stop("Could not resolve the runtime-receipt shared helper.", call. = FALSE)
}
sys.source(.paper3rt_helper_source, envir = environment())
.PAPER3RT_HELPER_PATH <- normalizePath(.paper3rt_helper_source,
                                       mustWork = TRUE)

.paper3rt_parse_producer_args <- function(args) {
  allowed <- c(
    "expected-package-version", "expected-package-commit",
    "expected-installed-package-tree-sha256", "runtime-image",
    "expected-runtime-image-sha256", "output-dir"
  )
  out <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Arguments require the exact --name=value form: ", arg,
           call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    if (!key %in% allowed) {
      stop("Unknown argument: --", key, call. = FALSE)
    }
    if (!is.null(out[[key]])) {
      stop("Duplicate argument: --", key, call. = FALSE)
    }
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  missing <- allowed[!vapply(allowed, function(key) {
    !is.null(out[[key]]) && nzchar(out[[key]])
  }, logical(1L))]
  if (length(missing)) {
    stop("Missing required argument(s): --",
         paste(missing, collapse = ", --"), call. = FALSE)
  }
  .paper3rt_assert_text(out[["expected-package-version"]],
                        "Expected package version")
  .paper3rt_assert_hex(out[["expected-package-commit"]], 40L,
                       "Expected package commit")
  .paper3rt_assert_hex(
    out[["expected-installed-package-tree-sha256"]], 64L,
    "Expected installed package tree SHA-256"
  )
  .paper3rt_assert_hex(out[["expected-runtime-image-sha256"]], 64L,
                       "Expected runtime image SHA-256")
  out
}

.paper3rt_git <- function(git, root, args, allow_failure = FALSE) {
  output <- system2(git, c("-C", root, args), stdout = TRUE, stderr = TRUE)
  failed <- !is.null(attr(output, "status"))
  if (failed && !allow_failure) {
    stop("Git provenance command failed: ", paste(output, collapse = "\n"),
         call. = FALSE)
  }
  output
}

.paper3rt_validate_producer_checkout <- function(
    expected_version, expected_commit) {
  producer <- normalizePath(.PAPER3RT_PRODUCER_PATH, mustWork = TRUE)
  helper <- normalizePath(.PAPER3RT_HELPER_PATH, mustWork = TRUE)
  expected_helper <- normalizePath(file.path(
    dirname(producer), "paper3_runtime_receipt_common.R"
  ), mustWork = TRUE)
  if (!identical(helper, expected_helper)) {
    stop("Producer did not source its exact sibling runtime helper.",
         call. = FALSE)
  }
  root <- normalizePath(file.path(dirname(producer), "..", ".."),
                        mustWork = TRUE)
  description <- file.path(root, "DESCRIPTION")
  if (!file.exists(description)) {
    stop("Producer package DESCRIPTION is missing.", call. = FALSE)
  }
  dcf <- read.dcf(description)
  if (!"Version" %in% colnames(dcf) ||
      !identical(unname(dcf[1L, "Version"]), expected_version)) {
    stop("Producer package version differs from its exact expected pin.",
         call. = FALSE)
  }
  git <- .paper3rt_resolve_executable("git")$path[[1L]]
  top <- .paper3rt_git(git, root, c("rev-parse", "--show-toplevel"))
  commit <- .paper3rt_git(git, root, c("rev-parse", "HEAD"))
  status <- .paper3rt_git(
    git, root, c("status", "--porcelain=v1", "--untracked-files=all")
  )
  if (length(top) != 1L ||
      !identical(normalizePath(top[[1L]], mustWork = TRUE), root) ||
      length(commit) != 1L || !identical(commit[[1L]], expected_commit) ||
      length(status)) {
    stop("Producer must be the exact expected clean Git checkout.",
         call. = FALSE)
  }
  prefix <- paste0(root, .Platform$file.sep)
  tracked <- c(producer, helper, normalizePath(description, mustWork = TRUE))
  if (any(!startsWith(tracked, prefix))) {
    stop("Producer provenance files escape the package checkout.",
         call. = FALSE)
  }
  relative <- substring(tracked, nchar(prefix) + 1L)
  for (path in relative) {
    .paper3rt_git(git, root, c("ls-files", "--error-unmatch", path))
  }
  list(
    root = root, package_version = expected_version,
    package_commit = expected_commit, producer_path = producer,
    helper_path = helper, producer_sha256 = .paper3rt_sha256(producer),
    helper_sha256 = .paper3rt_sha256(helper)
  )
}

.paper3rt_write_table <- function(x, path) {
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "")
  invisible(path)
}

.paper3rt_build_manifest <- function(stage, provenance) {
  files <- c("runtime_facts.tsv", "runtime_assets.tsv",
             "package_closure.tsv", "producer_provenance.tsv")
  roles <- c("runtime_facts", "runtime_assets", "package_closure",
             "producer_provenance")
  paths <- file.path(stage, files)
  if (any(!file.exists(paths))) {
    stop("Cannot manifest an incomplete runtime receipt.", call. = FALSE)
  }
  out <- data.table::data.table(
    file = files, role = roles,
    bytes = as.character(unname(file.info(paths)$size)),
    sha256 = vapply(paths, .paper3rt_sha256, character(1L)),
    package_version = provenance$package_version,
    package_commit = provenance$package_commit,
    installed_omicselector_tree_sha256 =
      provenance$installed_omicselector_tree_sha256,
    package_closure_sha256 = provenance$package_closure_sha256,
    runtime_image_path = provenance$runtime_image_path,
    runtime_image_sha256 = provenance$runtime_image_sha256,
    producer_script_sha256 = provenance$producer_script_sha256,
    helper_script_sha256 = provenance$helper_script_sha256
  )
  out
}

.paper3rt_build_runtime_receipt <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  opt <- .paper3rt_parse_producer_args(args)
  expected_version <- opt[["expected-package-version"]]
  expected_commit <- opt[["expected-package-commit"]]
  expected_tree <- opt[["expected-installed-package-tree-sha256"]]
  expected_image_sha <- opt[["expected-runtime-image-sha256"]]
  runtime_image <- .paper3rt_assert_canonical_file(
    opt[["runtime-image"]], "Runtime image", reject_supplied_symlink = TRUE
  )
  if (!identical(runtime_image, opt[["runtime-image"]]) ||
      !identical(.paper3rt_sha256(runtime_image), expected_image_sha)) {
    stop("Runtime image path/bytes differ from their exact canonical pins.",
         call. = FALSE)
  }
  source_start <- .paper3rt_validate_producer_checkout(
    expected_version, expected_commit
  )
  output_parent <- normalizePath(dirname(opt[["output-dir"]]), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt[["output-dir"]]))
  output_link <- Sys.readlink(output_dir)
  output_is_link <- length(output_link) == 1L && !is.na(output_link) &&
    nzchar(output_link)
  if (basename(output_dir) %in% c("", ".", "..") ||
      file.exists(output_dir) || dir.exists(output_dir) ||
      output_is_link) {
    stop("Output directory already exists or has an unsafe basename; runtime ",
         "receipts are immutable.", call. = FALSE)
  }
  if (startsWith(output_dir, paste0(source_start$root, .Platform$file.sep))) {
    stop("Runtime receipt output must be outside the producer checkout.",
         call. = FALSE)
  }
  capture_start <- .paper3rt_capture(runtime_image)
  omicselector <- capture_start$packages[package == "OmicSelector"]
  if (nrow(omicselector) != 1L ||
      omicselector$version[[1L]] != expected_version ||
      omicselector$tree_sha256[[1L]] != expected_tree) {
    stop("Installed OmicSelector version/tree differ from their exact pins.",
         call. = FALSE)
  }
  closure_sha <- .paper3rt_package_closure_sha256(capture_start$packages)
  image <- capture_start$assets[asset == "runtime_image"]
  provenance <- data.table::data.table(
    schema_version = "OmicSelector-paper3-runtime-receipt-v1",
    package_version = expected_version,
    package_commit = expected_commit, git_clean = "TRUE",
    producer_script_sha256 = source_start$producer_sha256,
    helper_script_sha256 = source_start$helper_sha256,
    installed_omicselector_tree_sha256 = expected_tree,
    package_closure_sha256 = closure_sha,
    runtime_image_path = image$path[[1L]],
    runtime_image_sha256 = image$sha256[[1L]]
  )
  stage <- tempfile(paste0(".", basename(output_dir), "-stage-"),
                    tmpdir = output_parent)
  dir.create(stage, recursive = FALSE)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  .paper3rt_write_table(capture_start$facts,
                        file.path(stage, "runtime_facts.tsv"))
  .paper3rt_write_table(capture_start$assets,
                        file.path(stage, "runtime_assets.tsv"))
  .paper3rt_write_table(capture_start$packages,
                        file.path(stage, "package_closure.tsv"))
  .paper3rt_write_table(provenance,
                        file.path(stage, "producer_provenance.tsv"))
  manifest <- .paper3rt_build_manifest(stage, provenance)
  .paper3rt_write_table(manifest, file.path(stage, "manifest.tsv"))
  manifest_sha <- .paper3rt_sha256(file.path(stage, "manifest.tsv"))
  capture_end <- .paper3rt_capture(runtime_image)
  if (!identical(.paper3rt_capture_digest(capture_start),
                 .paper3rt_capture_digest(capture_end))) {
    stop("Runtime facts/assets/package bytes changed during receipt capture.",
         call. = FALSE)
  }
  source_end <- .paper3rt_validate_producer_checkout(
    expected_version, expected_commit
  )
  if (!identical(source_start$producer_sha256, source_end$producer_sha256) ||
      !identical(source_start$helper_sha256, source_end$helper_sha256)) {
    stop("Producer/helper bytes changed during receipt capture.",
         call. = FALSE)
  }
  paper3_validate_runtime_receipt(
    stage, manifest_sha, expected_version, expected_commit, expected_tree,
    runtime_image, expected_image_sha, phase = "start"
  )
  if (!file.rename(stage, output_dir)) {
    stop("Could not atomically promote the runtime receipt directory.",
         call. = FALSE)
  }
  cat("PASS: atomically promoted runtime receipt\n")
  cat("receipt_dir=", output_dir, "\n", sep = "")
  cat("manifest_sha256=", manifest_sha, "\n", sep = "")
  invisible(list(output_dir = output_dir, manifest_sha256 = manifest_sha))
}

if (sys.nframe() == 0L) .paper3rt_build_runtime_receipt()
