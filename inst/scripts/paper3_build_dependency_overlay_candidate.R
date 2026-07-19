#!/usr/bin/env Rscript

# Assemble the paper-3 optional-dependency results into an immutable candidate.
# This utility only validates and serializes completed result shards. It never
# fits or scores a model and never mutates a source bundle or canonical table.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.OVERLAY_PRODUCER_PATH <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    candidate <- sub("^--file=", "", file_arg[[1L]])
    if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
  }
  candidate <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(candidate) && file.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }
  NA_character_
})

.overlay_within_methods <- c(
  "ai-scarf", "ai-tabpfn", "coda-codacore", "coda-deepcoda",
  "img-gasfcnn", "inv-scatter", "lrt-deepmaha", "nc-ecod-copod",
  "proto-net", "ssl-vicreg", "tab-tabdpt", "tab-tabicl", "unc-sngp"
)
.overlay_transfer_methods <- c(
  "cvae", "dann", "dg-fishr", "dg-ibirm", "dro-vrex", "moe-gated",
  "sel-stablemate"
)
.overlay_seeds <- c(101L, 202L, 303L, 404L, 505L)
.overlay_authorized_witness_keys <- data.table::rbindlist(list(
  data.table::CJ(
    method = "lrt-deepmaha", cohort = "GSE188627",
    seed = .overlay_seeds
  ),
  data.table::data.table(
    method = "lrt-deepmaha", cohort = "GSE270497",
    seed = c(101L, 303L, 404L, 505L)
  )
))

.overlay_sha256 <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

.overlay_scalar <- function(x, expected, label) {
  if (length(x) != 1L || is.na(x) || !identical(as.character(x), expected)) {
    stop(label, " does not match its exact expected pin.", call. = FALSE)
  }
  invisible(TRUE)
}

.overlay_safe_relative <- function(path, label = "path") {
  if (!is.character(path) || anyNA(path) || any(!nzchar(path)) ||
      any(grepl("(^|[/\\\\])\\.\\.($|[/\\\\])", path)) ||
      any(grepl("^[/\\\\]|^[A-Za-z]:[/\\\\]", path))) {
    stop(label, " must contain only non-empty relative paths.", call. = FALSE)
  }
  invisible(TRUE)
}

.overlay_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  expanded <- character()
  for (arg in args) {
    if (grepl("^--[^=]+=", arg)) {
      expanded <- c(expanded, sub("=.*$", "", arg), sub("^[^=]*=", "", arg))
    } else {
      expanded <- c(expanded, arg)
    }
  }
  value <- function(flag) {
    hit <- which(expanded == flag)
    if (length(hit) != 1L || hit == length(expanded) ||
        startsWith(expanded[[hit + 1L]], "--")) {
      return(NA_character_)
    }
    expanded[[hit + 1L]]
  }
  flags <- c(
    base_within_bundle = "--base-within-bundle",
    expected_base_within_bundle_sha256 =
      "--expected-base-within-bundle-sha256",
    base_transfer_bundle = "--base-transfer-bundle",
    expected_base_transfer_bundle_sha256 =
      "--expected-base-transfer-bundle-sha256",
    within_tasks = "--within-tasks",
    expected_within_task_sha256 = "--expected-within-task-sha256",
    within_shard_root = "--within-shard-root",
    within_split_manifest_dir = "--within-split-manifest-dir",
    transfer_tasks = "--transfer-tasks",
    expected_transfer_task_sha256 = "--expected-transfer-task-sha256",
    transfer_shard_root = "--transfer-shard-root",
    within_environment_receipt = "--within-environment-receipt",
    expected_within_environment_receipt_sha256 =
      "--expected-within-environment-receipt-sha256",
    transfer_environment_receipt = "--transfer-environment-receipt",
    expected_transfer_environment_receipt_sha256 =
      "--expected-transfer-environment-receipt-sha256",
    structural_witness = "--structural-witness",
    expected_structural_witness_sha256 =
      "--expected-structural-witness-sha256",
    expected_structural_cache_sha256 =
      "--expected-structural-cache-sha256",
    expected_structural_fold_engine_sha256 =
      "--expected-structural-fold-engine-sha256",
    expected_execution_package_version =
      "--expected-execution-package-version",
    expected_execution_package_commit =
      "--expected-execution-package-commit",
    expected_cache_package_commit = "--expected-cache-package-commit",
    expected_within_code_id = "--expected-within-code-id",
    expected_transfer_code_id = "--expected-transfer-code-id",
    expected_producer_package_version =
      "--expected-producer-package-version",
    expected_producer_package_commit =
      "--expected-producer-package-commit",
    output_dir = "--output-dir", run_mode = "--run-mode"
  )
  out <- lapply(flags, value)
  missing <- names(out)[vapply(out, function(x) is.na(x) || !nzchar(x),
                              logical(1L))]
  if (length(missing)) {
    stop("Missing required argument(s): ",
         paste(unname(flags[missing]), collapse = ", "), call. = FALSE)
  }
  duplicated_flags <- unique(expanded[duplicated(expanded) &
                                        expanded %in% unname(flags)])
  if (length(duplicated_flags)) {
    stop("Duplicate argument(s): ", paste(duplicated_flags, collapse = ", "),
         call. = FALSE)
  }
  if (!identical(out$run_mode, "full")) {
    stop("Only --run-mode=full is accepted for a promotable candidate.",
         call. = FALSE)
  }
  exact_hashes <- c(
    out$expected_base_within_bundle_sha256,
    out$expected_base_transfer_bundle_sha256,
    out$expected_within_task_sha256,
    out$expected_transfer_task_sha256,
    out$expected_within_environment_receipt_sha256,
    out$expected_transfer_environment_receipt_sha256,
    out$expected_structural_witness_sha256,
    out$expected_structural_cache_sha256,
    out$expected_structural_fold_engine_sha256
  )
  if (any(!grepl("^[0-9a-f]{64}$", exact_hashes)) ||
      !grepl("^[0-9a-f]{40}$", out$expected_execution_package_commit) ||
      !grepl("^[0-9a-f]{40}$", out$expected_cache_package_commit) ||
      !grepl("^[0-9a-f]{40}$", out$expected_producer_package_commit) ||
      !grepl("^[0-9a-f]{16}$", out$expected_within_code_id) ||
      !grepl("^[0-9a-f]{16}$", out$expected_transfer_code_id)) {
    stop("Task/witness hashes, commits, and code IDs must be exact lowercase ",
         "hexadecimal pins (64, 40, and 16 characters respectively).",
         call. = FALSE)
  }
  out
}

.overlay_git_output <- function(root, args, label) {
  out <- suppressWarnings(system2(
    "git", c("-C", shQuote(root), args), stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("Could not verify producer checkout ", label, ": ",
         paste(out, collapse = "\n"), call. = FALSE)
  }
  out
}

.overlay_validate_producer_checkout <- function(expected_version,
                                                expected_commit) {
  if (is.na(.OVERLAY_PRODUCER_PATH) || !file.exists(.OVERLAY_PRODUCER_PATH)) {
    stop("Could not resolve the package-owned producer script path.",
         call. = FALSE)
  }
  root <- normalizePath(file.path(dirname(.OVERLAY_PRODUCER_PATH), "..", ".."),
                        mustWork = TRUE)
  top <- trimws(.overlay_git_output(root, c("rev-parse", "--show-toplevel"),
                                    "top-level"))
  commit <- trimws(.overlay_git_output(root, c("rev-parse", "HEAD"), "HEAD"))
  status <- .overlay_git_output(
    root, c("status", "--porcelain=v1", "--untracked-files=all"), "status"
  )
  description <- file.path(root, "DESCRIPTION")
  if (length(top) != 1L || !identical(normalizePath(top), root) ||
      length(commit) != 1L || !identical(commit, expected_commit) ||
      length(status) || !file.exists(description)) {
    stop("Producer must run from the exact clean expected OmicSelector ",
         "checkout.", call. = FALSE)
  }
  version <- as.character(read.dcf(description, fields = "Version")[[1L]])
  if (!identical(version, expected_version)) {
    stop("Producer checkout DESCRIPTION version does not match its exact pin.",
         call. = FALSE)
  }
  relative_script <- substring(
    .OVERLAY_PRODUCER_PATH, nchar(paste0(root, .Platform$file.sep)) + 1L
  )
  if (!identical(relative_script,
                 "inst/scripts/paper3_build_dependency_overlay_candidate.R")) {
    stop("Producer script is not at its package-owned tracked path.",
         call. = FALSE)
  }
  tracked_blob <- trimws(.overlay_git_output(
    root, c("rev-parse", paste0("HEAD:", relative_script)),
    "tracked producer blob"
  ))
  current_blob <- trimws(.overlay_git_output(
    root, c("hash-object", shQuote(.OVERLAY_PRODUCER_PATH)),
    "current producer blob"
  ))
  if (length(tracked_blob) != 1L || length(current_blob) != 1L ||
      !identical(tracked_blob, current_blob)) {
    stop("Producer script bytes do not match the tracked expected checkout.",
         call. = FALSE)
  }
  list(
    root = root, version = version, commit = commit,
    script_sha256 = .overlay_sha256(.OVERLAY_PRODUCER_PATH),
    script_blob = current_blob
  )
}

.overlay_read_within_tasks <- function(path, expected_sha256) {
  if (!identical(.overlay_sha256(path), expected_sha256)) {
    stop("Within task table does not match its expected SHA-256.", call. = FALSE)
  }
  tasks <- data.table::fread(
    path, sep = "\t", header = FALSE, quote = "",
    col.names = c("method", "seed", "cohort", "label"),
    colClasses = c("character", "integer", "character", "character")
  )
  expected <- data.table::CJ(
    method = .overlay_within_methods, seed = .overlay_seeds,
    cohort = sort(unique(tasks$cohort)), sorted = TRUE
  )
  expected[, label := sprintf("%s__seed_%d__%s", method, seed, cohort)]
  data.table::setorder(tasks, method, seed, cohort)
  if (nrow(tasks) != 2210L || uniqueN(tasks$cohort) != 34L || anyNA(tasks) ||
      anyDuplicated(tasks[, .(method, seed, cohort)]) ||
      anyDuplicated(tasks$label) ||
      !identical(as.character(tasks$method), as.character(expected$method)) ||
      !identical(as.integer(tasks$seed), as.integer(expected$seed)) ||
      !identical(as.character(tasks$cohort), as.character(expected$cohort)) ||
      !identical(as.character(tasks$label), as.character(expected$label)) ||
      any(!grepl("^[A-Za-z0-9_.-]+$", tasks$label))) {
    stop("Within tasks are not the exact 13-by-34-by-five frozen grid.",
         call. = FALSE)
  }
  tasks
}

.overlay_read_transfer_tasks <- function(path, expected_sha256) {
  if (!identical(.overlay_sha256(path), expected_sha256)) {
    stop("Transfer task table does not match its expected SHA-256.", call. = FALSE)
  }
  tasks <- data.table::fread(
    path, sep = "\t", header = FALSE, quote = "",
    col.names = c("method", "seed", "label"),
    colClasses = c("character", "integer", "character")
  )
  expected <- data.table::CJ(
    method = .overlay_transfer_methods, seed = .overlay_seeds, sorted = TRUE
  )
  expected[, label := sprintf("%s__seed_%03d", method, seed)]
  data.table::setorder(tasks, method, seed)
  if (nrow(tasks) != 35L || anyNA(tasks) ||
      anyDuplicated(tasks[, .(method, seed)]) || anyDuplicated(tasks$label) ||
      !identical(as.character(tasks$method), as.character(expected$method)) ||
      !identical(as.integer(tasks$seed), as.integer(expected$seed)) ||
      !identical(as.character(tasks$label), as.character(expected$label)) ||
      any(!grepl("^[A-Za-z0-9_.-]+$", tasks$label))) {
    stop("Transfer tasks are not the exact seven-by-five frozen grid.",
         call. = FALSE)
  }
  tasks
}

.overlay_path_under <- function(relative, root) {
  .overlay_safe_relative(relative, "Splitter manifest path")
  root <- normalizePath(root, mustWork = TRUE)
  candidate <- file.path(root, relative)
  resolved <- normalizePath(candidate, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(resolved, prefix)) {
    stop("Splitter manifest path escapes --within-shard-root.", call. = FALSE)
  }
  resolved
}

.overlay_validate_split_manifests <- function(
    manifest_dir, shard_root, tasks, pins) {
  expected_names <- paste0(
    "output_manifest__", .overlay_within_methods, ".tsv"
  )
  observed_names <- list.files(
    manifest_dir, pattern = "^output_manifest__.*\\.tsv$",
    full.names = FALSE, recursive = FALSE, all.files = FALSE
  )
  if (!setequal(observed_names, expected_names) ||
      length(observed_names) != length(expected_names)) {
    stop("Expected exactly 13 method splitter manifests; one is missing, ",
         "extra, or duplicated.", call. = FALSE)
  }
  rows <- vector("list", length(expected_names))
  for (i in seq_along(expected_names)) {
    method <- .overlay_within_methods[[i]]
    manifest_path <- file.path(manifest_dir, expected_names[[i]])
    x <- data.table::fread(manifest_path, sep = "\t", quote = "")
    required <- c(
      "label", "method", "seed", "cohort", "path", "bytes", "sha256",
      "source_bundle_sha256", "task_table_sha256", "package_version",
      "package_commit", "cache_package_commit", "analysis_code_id"
    )
    if (nrow(x) != 170L || any(!required %in% names(x)) ||
        anyNA(x[, ..required]) || anyDuplicated(x$label) ||
        anyDuplicated(x[, .(method, seed, cohort)]) ||
        any(x$method != method) || any(x$task_table_sha256 != pins$task_sha) ||
        any(x$package_version != pins$package_version) ||
        any(x$package_commit != pins$package_commit) ||
        any(x$cache_package_commit != pins$cache_commit) ||
        any(x$analysis_code_id != pins$code_id) ||
        any(!grepl("^[0-9a-f]{64}$", x$source_bundle_sha256)) ||
        uniqueN(x$source_bundle_sha256) != 1L ||
        any(!grepl("^[0-9a-f]{64}$", x$sha256))) {
      stop("Splitter manifest violates row, key, or exact-pin contract: ",
           expected_names[[i]], call. = FALSE)
    }
    method_value <- method
    task_part <- data.table::copy(tasks[method == method_value])
    data.table::setorder(task_part, method, seed, cohort)
    data.table::setorder(x, method, seed, cohort)
    if (nrow(x) != nrow(task_part) ||
        !identical(as.character(x$method), as.character(task_part$method)) ||
        !identical(as.integer(x$seed), as.integer(task_part$seed)) ||
        !identical(as.character(x$cohort), as.character(task_part$cohort)) ||
        !identical(as.character(x$label), as.character(task_part$label))) {
      stop("Splitter manifest keys do not equal the frozen task table: ",
           expected_names[[i]], call. = FALSE)
    }
    paths <- vapply(x$path, .overlay_path_under, character(1L),
                    root = shard_root)
    observed_bytes <- unname(file.info(paths)$size)
    observed_sha <- vapply(paths, .overlay_sha256, character(1L))
    if (!identical(as.numeric(x$bytes), as.numeric(observed_bytes)) ||
        !identical(as.character(x$sha256), unname(observed_sha))) {
      stop("A shard changed after its splitter manifest was generated: ",
           expected_names[[i]], call. = FALSE)
    }
    x[, `:=`(
      manifest_name = expected_names[[i]],
      manifest_sha256 = .overlay_sha256(manifest_path),
      resolved_path = paths
    )]
    if (!"structural_witness_sha256" %in% names(x)) {
      x[, structural_witness_sha256 := NA_character_]
    }
    rows[[i]] <- x
  }
  out <- data.table::rbindlist(rows, fill = TRUE)
  data.table::setorder(out, method, seed, cohort)
  if (nrow(out) != 2210L || anyDuplicated(out[, .(method, seed, cohort)]) ||
      !identical(as.character(out$method), as.character(tasks$method)) ||
      !identical(as.integer(out$seed), as.integer(tasks$seed)) ||
      !identical(as.character(out$cohort), as.character(tasks$cohort)) ||
      !identical(as.character(out$label), as.character(tasks$label))) {
    stop("Combined splitter manifests do not equal all 2,210 frozen tasks.",
         call. = FALSE)
  }
  out
}

.overlay_validate_bundle_pin <- function(bundle, pins, label) {
  required <- c("cells", "predictions", "package_version", "package_commit",
                "cache_package_commit", "analysis_code_id")
  if (!is.list(bundle) || any(!required %in% names(bundle))) {
    stop(label, " lacks required bundle fields.", call. = FALSE)
  }
  .overlay_scalar(bundle$package_version, pins$package_version,
                  paste(label, "package_version"))
  .overlay_scalar(bundle$package_commit, pins$package_commit,
                  paste(label, "package_commit"))
  .overlay_scalar(bundle$cache_package_commit, pins$cache_commit,
                  paste(label, "cache_package_commit"))
  .overlay_scalar(bundle$analysis_code_id, pins$code_id,
                  paste(label, "analysis_code_id"))
  invisible(TRUE)
}

.overlay_splitter_environment <- function() {
  path <- file.path(dirname(.OVERLAY_PRODUCER_PATH),
                    "paper3_split_dependency_within_bundle.R")
  if (!file.exists(path)) {
    stop("Package-owned dependency splitter is missing beside the producer.",
         call. = FALSE)
  }
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

.overlay_read_within_shards <- function(manifests, pins) {
  cells <- vector("list", nrow(manifests))
  predictions <- vector("list", nrow(manifests))
  witnesses <- vector("list", nrow(manifests))
  templates <- list()
  for (i in seq_len(nrow(manifests))) {
    row <- manifests[i]
    bundle <- readRDS(row$resolved_path)
    .overlay_validate_bundle_pin(bundle, pins, "Within shard")
    cell <- data.table::as.data.table(bundle$cells)
    pred <- data.table::as.data.table(bundle$predictions)
    required_cell <- c(
      "method", "seed", "cohort", "eligible", "ineligible_reason", "n_eff",
      "package_version", "package_commit", "cache_package_commit",
      "analysis_code_id"
    )
    if (nrow(cell) != 1L || any(!required_cell %in% names(cell)) ||
        !identical(cell$method, row$method) ||
        !identical(as.integer(cell$seed), as.integer(row$seed)) ||
        !identical(cell$cohort, row$cohort) ||
        any(cell$package_version != pins$package_version) ||
        any(cell$package_commit != pins$package_commit) ||
        any(cell$cache_package_commit != pins$cache_commit) ||
        any(cell$analysis_code_id != pins$code_id)) {
      stop("Within shard key or embedded pins disagree with its manifest: ",
           row$label, call. = FALSE)
    }
    if (nrow(pred)) {
      required_pred <- c("method", "seed", "cohort", "sample_idx",
                         "score_m", "score_b", "package_commit",
                         "analysis_code_id")
      if (any(!required_pred %in% names(pred)) ||
          any(pred$method != row$method) || any(pred$seed != row$seed) ||
          any(pred$cohort != row$cohort) ||
          any(pred$package_commit != pins$package_commit) ||
          any(pred$analysis_code_id != pins$code_id)) {
        stop("Within shard predictions violate key or exact-pin contract: ",
             row$label, call. = FALSE)
      }
    }
    cells[[i]] <- data.table::copy(cell)
    predictions[[i]] <- data.table::copy(pred)
    witnesses[i] <- list(bundle$structural_ineligibility_witness)
    if (is.null(templates[[row$method]])) templates[[row$method]] <- bundle
  }
  list(
    cells = data.table::rbindlist(cells, fill = TRUE),
    predictions = data.table::rbindlist(predictions, fill = TRUE),
    witnesses = witnesses, templates = templates
  )
}

.overlay_validate_reassembled_within <- function(within, tasks, manifests,
                                                 witness, witness_sha, pins) {
  splitter <- .overlay_splitter_environment()
  for (method in .overlay_within_methods) {
    method_value <- method
    bundle <- within$templates[[method]]
    bundle$cells <- data.table::copy(within$cells[method == method_value])
    bundle$predictions <- data.table::copy(
      within$predictions[method == method_value]
    )
    bundle$methods <- method
    bundle$seeds <- .overlay_seeds
    bundle$cohorts <- sort(unique(tasks[method == method_value, cohort]))
    expected <- list(
      method = method, package_version = pins$package_version,
      package_commit = pins$package_commit, cache_commit = pins$cache_commit,
      analysis_code_id = pins$code_id
    )
    source_sha <- unique(manifests[method == method_value, source_bundle_sha256])
    structural <- identical(method, "lrt-deepmaha")
    splitter$.dependency_validate_bundle(
      bundle, tasks, expected,
      structural_witness = if (structural) witness else NULL,
      source_bundle_sha256 = if (structural) source_sha else NULL,
      task_table_sha256 = if (structural) pins$task_sha else NULL,
      expected_cache_sha256 = if (structural) pins$structural_cache_sha else NULL,
      expected_fold_engine_sha256 = if (structural) pins$fold_engine_sha else NULL,
      witness_sha256 = if (structural) witness_sha else NULL
    )
  }
  invisible(TRUE)
}

.overlay_read_transfer_shards <- function(root, tasks, pins) {
  paths <- file.path(root, tasks$label, "transfer_bundle.rds")
  if (any(!file.exists(paths))) {
    stop("Missing transfer shard(s): ",
         paste(tasks$label[!file.exists(paths)], collapse = ", "),
         call. = FALSE)
  }
  cells <- vector("list", length(paths))
  predictions <- vector("list", length(paths))
  rows <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    bundle <- readRDS(paths[[i]])
    .overlay_validate_bundle_pin(bundle, pins, "Transfer shard")
    cell <- data.table::as.data.table(bundle$cells)
    pred <- data.table::as.data.table(bundle$predictions)
    required_cell <- c(
      "method", "seed", "disease", "held_out", "eligible",
      "ineligible_reason", "n_test", "n_training_cohorts", "package_version",
      "package_commit", "cache_package_commit", "analysis_code_id"
    )
    if (nrow(cell) != 18L || any(!required_cell %in% names(cell)) ||
        !identical(unique(cell$method), tasks$method[[i]]) ||
        !identical(as.integer(unique(cell$seed)), tasks$seed[[i]]) ||
        any(cell$package_version != pins$package_version) ||
        any(cell$package_commit != pins$package_commit) ||
        any(cell$cache_package_commit != pins$cache_commit) ||
        any(cell$analysis_code_id != pins$code_id)) {
      stop("Transfer shard violates its exact 18-cell key/pin contract: ",
           tasks$label[[i]], call. = FALSE)
    }
    required_pred <- c("method", "seed", "disease", "held_out", "sample_id",
                       "group_id", "y", "score_m", "score_b", "package_commit",
                       "analysis_code_id")
    if (any(!required_pred %in% names(pred)) ||
        (nrow(pred) &&
         (any(pred$method != tasks$method[[i]]) ||
          any(pred$seed != tasks$seed[[i]]) ||
          any(pred$package_commit != pins$package_commit) ||
          any(pred$analysis_code_id != pins$code_id)))) {
      stop("Transfer predictions violate key or exact-pin contract: ",
           tasks$label[[i]], call. = FALSE)
    }
    eligible_keys <- cell[eligible == TRUE, .(disease, held_out)]
    prediction_keys <- unique(pred[, .(disease, held_out)])
    if (!setequal(
      .overlay_key_string(eligible_keys, c("disease", "held_out")),
      .overlay_key_string(prediction_keys, c("disease", "held_out"))
    )) {
      stop("Transfer prediction directions do not exactly equal the eligible ",
           "18-cell shard surface: ", tasks$label[[i]], call. = FALSE)
    }
    cells[[i]] <- data.table::copy(cell)
    predictions[[i]] <- data.table::copy(pred)
    rows[[i]] <- data.table::data.table(
      label = tasks$label[[i]], method = tasks$method[[i]],
      seed = tasks$seed[[i]], path = file.path(tasks$label[[i]],
                                               "transfer_bundle.rds"),
      bytes = file.info(paths[[i]])$size,
      sha256 = .overlay_sha256(paths[[i]]),
      task_table_sha256 = pins$task_sha,
      package_version = pins$package_version,
      package_commit = pins$package_commit,
      cache_package_commit = pins$cache_commit,
      analysis_code_id = pins$code_id
    )
  }
  out_cells <- data.table::rbindlist(cells, fill = TRUE)
  key <- c("disease", "held_out", "seed", "method")
  if (nrow(out_cells) != 630L || anyDuplicated(out_cells[, ..key]) ||
      uniqueN(out_cells[, .(disease, held_out, n_training_cohorts)]) != 18L ||
      sum(out_cells$n_training_cohorts >= 2L) != 420L ||
      sum(out_cells$n_training_cohorts < 2L) != 210L ||
      any(out_cells[n_training_cohorts < 2L, eligible == TRUE])) {
    stop("Transfer shards do not reproduce the exact 630-cell frozen design.",
         call. = FALSE)
  }
  list(
    cells = out_cells,
    predictions = data.table::rbindlist(predictions, fill = TRUE),
    manifest = data.table::rbindlist(rows)
  )
}

.overlay_validate_transfer_science <- function(cells, predictions) {
  if (!is.logical(cells$eligible) || anyNA(cells$eligible) ||
      any(cells[n_training_cohorts >= 2L, eligible != TRUE]) ||
      any(cells[n_training_cohorts < 2L, eligible != FALSE]) ||
      sum(cells$eligible) != 420L ||
      any(cells[eligible == TRUE,
                !is.na(ineligible_reason) & nzchar(ineligible_reason)]) ||
      any(cells[eligible == FALSE,
                is.na(ineligible_reason) | !nzchar(ineligible_reason)]) ||
      any(grepl("driver[[:space:]]+cell[[:space:]]+error:|cell[[:space:]]+error:",
                cells$ineligible_reason, ignore.case = TRUE), na.rm = TRUE) ||
      ("n_group_split_violations" %in% names(cells) &&
       any(cells$n_group_split_violations != 0L, na.rm = TRUE))) {
    stop("Transfer cells contain invalid eligibility, driver errors, or group ",
         "splits.", call. = FALSE)
  }
  required_prediction <- c(
    "method", "seed", "disease", "held_out", "sample_id", "group_id", "y",
    "score_m", "score_b", "package_commit", "analysis_code_id"
  )
  if (nrow(predictions) &&
      (any(!required_prediction %in% names(predictions)) ||
       anyNA(predictions[, ..required_prediction]) ||
       any(!predictions$y %in% c(0L, 1L)) ||
       any(!nzchar(predictions$sample_id)) ||
       any(!nzchar(predictions$group_id)))) {
    stop("Transfer predictions lack complete sample/group/outcome evidence.",
         call. = FALSE)
  }
  if (nrow(predictions)) {
    eligible_keys <- unique(cells[eligible == TRUE,
                                  .(disease, held_out, method, seed)])
    prediction_keys <- unique(predictions[,
                                          .(disease, held_out, method, seed)])
    if (!setequal(
      .overlay_key_string(eligible_keys,
                          c("disease", "held_out", "method", "seed")),
      .overlay_key_string(prediction_keys,
                          c("disease", "held_out", "method", "seed"))
    )) {
      stop("Transfer prediction keys do not exactly equal eligible cell keys.",
           call. = FALSE)
    }
    group_labels <- predictions[, .(
      n_labels = uniqueN(y), n_held_out = uniqueN(held_out)
    ), by = group_id]
    if (any(group_labels$n_labels != 1L)) {
      stop("Transfer provenance group has inconsistent outcome labels.",
           call. = FALSE)
    }
    if (any(group_labels$n_held_out != 1L)) {
      stop("Transfer provenance group appears in more than one held-out ",
           "cohort.", call. = FALSE)
    }
    heldout_surfaces <- predictions[, {
      surface <- unique(.SD)
      data.table::setorderv(surface, c("sample_id", "group_id", "y"))
      .(surface_sha256 = .overlay_semantic_digest(surface))
    }, by = .(disease, held_out, method, seed),
    .SDcols = c("sample_id", "group_id", "y")]
    discordant <- heldout_surfaces[, .(
      n_surfaces = uniqueN(surface_sha256)
    ), by = .(disease, held_out)]
    if (any(discordant$n_surfaces != 1L)) {
      stop("Transfer methods or seeds use different held-out sample/group/",
           "outcome surfaces.", call. = FALSE)
    }
  }
  environment_methods <- c(
    "dann", "dg-fishr", "dg-ibirm", "dro-vrex", "sel-stablemate"
  )
  required_engagement <- c("engagement_metric", "engagement_ok",
                           "engagement_value")
  if (any(!required_engagement %in% names(cells)) ||
      any(cells[
        eligible == TRUE & method %in% environment_methods,
        engagement_metric != "n_environments" | engagement_ok != TRUE |
          !is.finite(engagement_value) | engagement_value <= 1L
      ]) ||
      any(cells[
        eligible == TRUE & method == "cvae",
        engagement_metric != "not_applicable_class_conditional" |
          engagement_ok != TRUE
      ]) ||
      any(cells[
        eligible == TRUE & method == "moe-gated",
        engagement_metric != "not_applicable_generic_mixture_n_experts" |
          engagement_ok != TRUE | !is.finite(engagement_value) |
          engagement_value <= 1L
      ])) {
    stop("Transfer engagement evidence does not match the frozen method ",
         "contract.", call. = FALSE)
  }
  invisible(TRUE)
}

.overlay_key_string <- function(x, key) {
  do.call(paste, c(x[, ..key], sep = "\r"))
}

.overlay_validate_witness_bindings <- function(
    witness, witness_sha, manifests, within_cells, embedded, pins) {
  required <- c(
    "schema_version", "producer", "producer_script_sha256",
    "source_bundle_sha256", "task_table_sha256",
    "cache_sha256", "fold_engine_sha256", "package_version",
    "package_commit", "cache_package_commit", "analysis_code_id", "cells",
    "predictions", "fold_audit"
  )
  if (!is.list(witness) || any(!required %in% names(witness))) {
    stop("Structural witness lacks required fields.", call. = FALSE)
  }
  .overlay_scalar(witness$schema_version,
                  "OmicSelector-dependency-structural-witness-v1",
                  "Structural witness schema_version")
  .overlay_scalar(witness$task_table_sha256, pins$task_sha,
                  "Structural witness task_table_sha256")
  .overlay_scalar(witness$cache_sha256, pins$structural_cache_sha,
                  "Structural witness cache_sha256")
  .overlay_scalar(witness$fold_engine_sha256, pins$fold_engine_sha,
                  "Structural witness fold_engine_sha256")
  .overlay_scalar(witness$package_version, pins$package_version,
                  "Structural witness package_version")
  .overlay_scalar(witness$package_commit, pins$package_commit,
                  "Structural witness package_commit")
  .overlay_scalar(witness$cache_package_commit, pins$cache_commit,
                  "Structural witness cache_package_commit")
  .overlay_scalar(witness$analysis_code_id, pins$code_id,
                  "Structural witness analysis_code_id")
  deepmaha_source <- unique(
    manifests[method == "lrt-deepmaha", source_bundle_sha256]
  )
  if (length(deepmaha_source) != 1L) {
    stop("DeepMaha splitter manifest lacks one source-bundle SHA.",
         call. = FALSE)
  }
  .overlay_scalar(witness$source_bundle_sha256, deepmaha_source,
                  "Structural witness source_bundle_sha256")
  key <- c("method", "cohort", "seed")
  wc <- data.table::as.data.table(witness$cells)
  wp <- data.table::as.data.table(witness$predictions)
  wa <- data.table::as.data.table(witness$fold_audit)
  expected <- data.table::copy(.overlay_authorized_witness_keys)
  data.table::setorder(wc, method, cohort, seed)
  data.table::setorder(expected, method, cohort, seed)
  if (nrow(wc) != 9L || any(!key %in% names(wc)) ||
      anyDuplicated(wc[, ..key]) ||
      !identical(wc[, ..key], expected[, ..key])) {
    stop("Structural witness must authorize exactly the nine reviewed ",
         "DeepMaha cells.", call. = FALSE)
  }
  for (surface in list(predictions = wp, fold_audit = wa)) {
    if (any(!key %in% names(surface)) || !nrow(surface) ||
        !setequal(.overlay_key_string(unique(surface[, ..key]), key),
                  .overlay_key_string(expected, key))) {
      stop("Structural witness prediction/audit surfaces do not cover its ",
           "exact nine-cell authority.", call. = FALSE)
    }
  }
  manifest_witness <- manifests[
    !is.na(structural_witness_sha256) & nzchar(structural_witness_sha256),
    c(key, "structural_witness_sha256"), with = FALSE
  ]
  data.table::setorder(manifest_witness, method, cohort, seed)
  if (nrow(manifest_witness) != 9L ||
      !identical(manifest_witness[, ..key], expected[, ..key]) ||
      any(manifest_witness$structural_witness_sha256 != witness_sha)) {
    stop("Witness SHA must occur on exactly the nine authorized splitter rows ",
         "and nowhere else.", call. = FALSE)
  }
  cell_keys <- within_cells[eligible == FALSE, ..key]
  if (!all(.overlay_key_string(expected, key) %in%
           .overlay_key_string(cell_keys, key))) {
    stop("A witnessed DeepMaha cell is not ineligible in the accepted shards.",
         call. = FALSE)
  }
  source_surfaces <- list(cells = wc, predictions = wp, fold_audit = wa)
  for (i in seq_len(nrow(manifests))) {
    row <- manifests[i]
    has_witness <- !is.null(embedded[[i]])
    authorized <- .overlay_key_string(row, key) %in%
      .overlay_key_string(expected, key)
    if (!identical(has_witness, authorized)) {
      stop("Embedded structural witness appears on a missing or unauthorized ",
           "cell: ", row$label, call. = FALSE)
    }
    if (!authorized) next
    z <- embedded[[i]]
    scalar_fields <- c(
      "schema_version", "producer", "producer_script_sha256",
      "source_bundle_sha256", "task_table_sha256", "cache_sha256"
    )
    for (field in scalar_fields) {
      if (!field %in% names(z) || !field %in% names(witness) ||
          !identical(z[[field]], witness[[field]])) {
        stop("Embedded witness scalar disagrees with the source witness: ",
             row$label, call. = FALSE)
      }
    }
    .overlay_scalar(z$witness_sha256, witness_sha,
                    "Embedded structural witness SHA")
    for (surface in names(source_surfaces)) {
      if (!surface %in% names(z)) {
        stop("Embedded witness lacks ", surface, ": ", row$label,
             call. = FALSE)
      }
      observed <- data.table::as.data.table(z[[surface]])
      source <- source_surfaces[[surface]][
        method == row$method & cohort == row$cohort & seed == row$seed
      ]
      if (!identical(observed, source)) {
        stop("Embedded witness surface disagrees with source witness: ",
             row$label, call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

.overlay_validate_predictions <- function(cells, predictions, key,
                                          profile_key, n_column, label) {
  required <- c(key, profile_key, "score_m", "score_b")
  if (any(!required %in% names(predictions)) && nrow(predictions)) {
    stop(label, " predictions lack required columns.", call. = FALSE)
  }
  if (nrow(predictions) &&
      (anyNA(predictions[, ..required]) ||
       any(!is.finite(predictions$score_m)) ||
       any(!is.finite(predictions$score_b)) ||
       anyDuplicated(predictions[, c(key, profile_key), with = FALSE]))) {
    stop(label, " predictions are missing, non-finite, or duplicated.",
         call. = FALSE)
  }
  counts <- predictions[, .(n_predictions = .N), by = key]
  expected_keys <- unique(cells[eligible == TRUE, ..key])
  observed_keys <- unique(counts[, ..key])
  if (!setequal(.overlay_key_string(expected_keys, key),
                .overlay_key_string(observed_keys, key))) {
    stop(label, " prediction keys do not exactly equal eligible cell keys.",
         call. = FALSE)
  }
  eligible <- merge(
    cells[eligible == TRUE, c(key, n_column), with = FALSE], counts,
    by = key, all.x = TRUE, sort = FALSE
  )
  data.table::setnames(eligible, n_column, "expected_predictions")
  if (nrow(eligible) &&
      (anyNA(eligible$n_predictions) ||
       any(eligible$n_predictions != eligible$expected_predictions))) {
    stop("Eligible ", label,
         " cells do not have one prediction per effective profile.",
         call. = FALSE)
  }
  if (nrow(merge(cells[eligible == FALSE, ..key], counts, by = key))) {
    stop("Ineligible ", label, " cell has predictions.", call. = FALSE)
  }
  invisible(TRUE)
}

.overlay_sanitize <- function(x) {
  x <- data.table::copy(x)
  char <- names(x)[vapply(x, is.character, logical(1L))]
  for (column in char) x[[column]] <- gsub("[\r\n\t]+", " ", x[[column]])
  x
}

.overlay_semantic_digest <- function(x) {
  payload <- list(
    names = names(x), classes = lapply(x, class),
    columns = lapply(x, function(column) {
      attributes(column) <- attributes(column)[intersect(
        names(attributes(column)), c("class", "levels", "tzone", "units")
      )]
      column
    })
  )
  digest::digest(serialize(payload, NULL, version = 3L), algo = "sha256")
}

.overlay_table_identical <- function(x, y) {
  identical(names(x), names(y)) && identical(nrow(x), nrow(y)) &&
    all(vapply(names(x), function(column) {
      identical(x[[column]], y[[column]])
    }, logical(1L)))
}

.overlay_integrate <- function(base_bundle, overlay_cells, overlay_predictions,
                               methods, key, prediction_order, n_expected,
                               label, pins, execution_route) {
  if (!is.list(base_bundle) ||
      any(!c("cells", "predictions") %in% names(base_bundle))) {
    stop("Base ", label, " bundle lacks cells or predictions.", call. = FALSE)
  }
  base_cells <- data.table::copy(data.table::as.data.table(base_bundle$cells))
  base_predictions <- data.table::copy(
    data.table::as.data.table(base_bundle$predictions)
  )
  selected <- base_cells[method %in% methods]
  if (nrow(selected) != nrow(overlay_cells) || any(selected$eligible) ||
      !setequal(.overlay_key_string(selected, key),
                .overlay_key_string(overlay_cells, key))) {
    stop("Base ", label,
         " bundle no longer exposes the exact reviewed ineligible surface.",
         call. = FALSE)
  }
  overlay_cells <- data.table::copy(overlay_cells)
  overlay_predictions <- data.table::copy(overlay_predictions)
  route_value <- execution_route
  required_route_columns <- c(
    "execution_route", "final_package_version", "final_package_commit"
  )
  if (any(!required_route_columns %in% names(base_cells)) ||
      (nrow(base_predictions) &&
       any(!required_route_columns %in% names(base_predictions)))) {
    stop("Base ", label, " bundle lacks route columns required to preserve ",
         "unrelated rows exactly.", call. = FALSE)
  }
  overlay_cells[, `:=`(
    execution_route = route_value,
    final_package_version = pins$package_version,
    final_package_commit = pins$package_commit
  )]
  if (nrow(overlay_predictions)) overlay_predictions[, `:=`(
    execution_route = route_value,
    final_package_version = pins$package_version,
    final_package_commit = pins$package_commit
  )]
  keep_cells <- base_cells[!method %in% methods]
  keep_predictions <- base_predictions[!method %in% methods]
  data.table::setorderv(keep_cells, key)
  if (nrow(keep_predictions)) data.table::setorderv(
    keep_predictions, prediction_order
  )
  unsafe_character <- function(x) {
    columns <- names(x)[vapply(x, is.character, logical(1L))]
    any(vapply(columns, function(column) {
      any(grepl("[\r\n\t]", x[[column]]), na.rm = TRUE)
    }, logical(1L)))
  }
  if (unsafe_character(keep_cells) || unsafe_character(keep_predictions)) {
    stop("Unrelated base ", label,
         " rows are not safe for lossless tabular serialization.",
         call. = FALSE)
  }
  before_cell_sha <- .overlay_semantic_digest(keep_cells)
  before_prediction_sha <- .overlay_semantic_digest(keep_predictions)
  cells <- data.table::rbindlist(list(keep_cells, overlay_cells), fill = TRUE)
  predictions <- data.table::rbindlist(
    list(keep_predictions, overlay_predictions), fill = TRUE
  )
  data.table::setorderv(cells, key)
  if (nrow(predictions)) data.table::setorderv(predictions, prediction_order)
  if (nrow(cells) != n_expected || anyDuplicated(cells[, ..key])) {
    stop("Integrated ", label, " grid is incomplete or duplicated.",
         call. = FALSE)
  }
  if (any(cells[method %in% methods, execution_route] != route_value) ||
      (nrow(predictions[method %in% methods]) &&
       any(predictions[method %in% methods, execution_route] != route_value))) {
    stop("Integrated ", label, " dependency rows did not receive their exact ",
         "execution route.", call. = FALSE)
  }
  after_cells <- cells[!method %in% methods]
  after_predictions <- predictions[!method %in% methods]
  if (!identical(before_cell_sha, .overlay_semantic_digest(after_cells)) ||
      !identical(before_prediction_sha,
                 .overlay_semantic_digest(after_predictions)) ||
      !.overlay_table_identical(keep_cells, after_cells) ||
      !.overlay_table_identical(keep_predictions, after_predictions)) {
    stop("An unrelated ", label, " row changed during overlay.",
         call. = FALSE)
  }
  out <- base_bundle
  out$cells <- cells
  out$predictions <- predictions
  out$dependency_execution_package_version <- pins$package_version
  out$dependency_execution_package_commit <- pins$package_commit
  out$dependency_analysis_code_id <- pins$code_id
  out$dependency_methods <- methods
  out$dependency_execution_route <- execution_route
  list(bundle = out, cells = cells, predictions = predictions)
}

.overlay_receipt_files <- function(path, route) {
  files <- if (identical(route, "transfer_argos")) {
    file.path(path, c(
      "environment_manifest.tsv", "public_model_manifest.tsv",
      "python_freeze.txt", "r_package_manifest.tsv"
    ))
  } else if (identical(route, "within_s5")) {
    list.files(path, recursive = FALSE, all.files = TRUE,
               full.names = TRUE, include.dirs = FALSE, no.. = TRUE)
  } else {
    stop("Unknown environment receipt route.", call. = FALSE)
  }
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files) || any(!file.exists(files)) ||
      any(nzchar(Sys.readlink(files)))) {
    stop("Environment receipt must contain regular immutable files only.",
         call. = FALSE)
  }
  root <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  relative <- substring(normalizePath(files, mustWork = TRUE),
                        nchar(prefix) + 1L)
  data.table::data.table(
    source_path = files, path = relative,
    bytes = unname(file.info(files)$size),
    sha256 = vapply(files, .overlay_sha256, character(1L))
  )[order(path)]
}

.overlay_receipt_digest <- function(files) {
  digest::digest(paste(files$path, files$bytes, files$sha256,
                       sep = "\t", collapse = "\n"), algo = "sha256")
}

.overlay_validate_environment_receipt <- function(path, route,
                                                   package_commit,
                                                   package_version,
                                                   expected_receipt_sha256) {
  files <- .overlay_receipt_files(path, route)
  names <- files$path
  receipt_sha <- .overlay_receipt_digest(files)
  if (!identical(receipt_sha, expected_receipt_sha256)) {
    stop("Environment receipt does not match its exact route-specific SHA-256 ",
         "pin.", call. = FALSE)
  }
  runtime_assets <- NULL
  if (identical(route, "within_s5")) {
    required <- c(
      "captured_utc.txt", "hostname.txt", "package_commit.txt",
      "package_git_status.txt", "python_conda_explicit.txt",
      "python_conda_list.txt", "python_freeze.txt", "python_runtime.txt",
      "r_package_manifest.tsv", "r_runtime.txt", "receipt_files.tsv",
      "receipt_sha256.txt", "runtime_asset_sha256.txt"
    )
    if (!setequal(required, names)) {
      stop("Within-s5 receipt lacks its complete runtime/package/asset closure.",
           call. = FALSE)
    }
    host <- trimws(readLines(file.path(path, "hostname.txt"), warn = FALSE))
    commit <- trimws(readLines(file.path(path, "package_commit.txt"),
                             warn = FALSE))
    status <- readLines(file.path(path, "package_git_status.txt"), warn = FALSE)
    if (!identical(host, "s5") || !identical(commit, package_commit) ||
        length(status)) {
      stop("Within route is not bound to a clean exact s5 environment receipt.",
           call. = FALSE)
    }
    receipt <- readLines(file.path(path, "receipt_sha256.txt"), warn = FALSE)
    matched <- regexec("^([0-9a-f]{64})[[:space:]]+(.+)$", receipt)
    parsed <- regmatches(receipt, matched)
    if (!length(parsed) || any(lengths(parsed) != 3L)) {
      stop("Within-s5 receipt SHA closure is malformed.", call. = FALSE)
    }
    listed_sha <- vapply(parsed, `[[`, character(1L), 2L)
    listed_name <- basename(vapply(parsed, `[[`, character(1L), 3L))
    listed_paths <- file.path(path, listed_name)
    if (any(!file.exists(listed_paths)) || anyDuplicated(listed_name)) {
      stop("Within-s5 receipt SHA closure names missing or duplicate files.",
           call. = FALSE)
    }
    observed <- vapply(listed_paths, .overlay_sha256, character(1L))
    if (!identical(listed_sha, unname(observed))) {
      stop("Within-s5 receipt SHA closure does not match current bytes.",
           call. = FALSE)
    }
    if (!setequal(listed_name, setdiff(required, "receipt_sha256.txt"))) {
      stop("Within-s5 receipt SHA closure does not enumerate every receipt ",
           "artifact except itself.", call. = FALSE)
    }
    python_runtime <- readLines(file.path(path, "python_runtime.txt"),
                                warn = FALSE)
    required_python <- c(
      "python=3.12.12", "numpy=2.4.2", "pandas=2.3.3",
      "sklearn=1.8.0", "scipy=1.16.3", "torch=2.11.0",
      "torch_cuda_available=False", "torch_cuda_version=None"
    )
    r_runtime <- readLines(file.path(path, "r_runtime.txt"), warn = FALSE)
    r_packages <- data.table::fread(file.path(path, "r_package_manifest.tsv"))
    if (!all(required_python %in% python_runtime) ||
        !any(grepl("^R version 4\\.5\\.2", r_runtime)) ||
        nrow(r_packages[Package == "OmicSelector" & Version ==
                          package_version]) != 1L) {
      stop("Within-s5 receipt does not carry the exact reviewed R/Python ",
           "dependency versions.", call. = FALSE)
    }
    asset_lines <- readLines(file.path(path, "runtime_asset_sha256.txt"),
                             warn = FALSE)
    asset_match <- regexec("^([0-9a-f]{64})[[:space:]]+(.+)$", asset_lines)
    asset_parsed <- regmatches(asset_lines, asset_match)
    if (!length(asset_parsed) || any(lengths(asset_parsed) != 3L)) {
      stop("Within-s5 runtime asset manifest is malformed.", call. = FALSE)
    }
    runtime_assets <- data.table::data.table(
      sha256 = vapply(asset_parsed, `[[`, character(1L), 2L),
      path = vapply(asset_parsed, `[[`, character(1L), 3L)
    )
    if (any(!file.exists(runtime_assets$path)) ||
        anyDuplicated(runtime_assets$path) ||
        any(runtime_assets$sha256 != vapply(
          runtime_assets$path, .overlay_sha256, character(1L)
        ))) {
      stop("Within-s5 runtime asset bytes do not match their receipt hashes.",
           call. = FALSE)
    }
  } else if (identical(route, "transfer_argos")) {
    required <- c("environment_manifest.tsv", "public_model_manifest.tsv",
                  "python_freeze.txt", "r_package_manifest.tsv")
    if (any(!required %in% names)) {
      stop("Transfer-Argos receipt lacks its frozen four-file environment ",
           "closure.", call. = FALSE)
    }
    environment <- data.table::fread(file.path(path, "environment_manifest.tsv"))
    models <- data.table::fread(file.path(path, "public_model_manifest.tsv"))
    python_freeze <- readLines(file.path(path, "python_freeze.txt"),
                               warn = FALSE)
    r_packages <- data.table::fread(file.path(path, "r_package_manifest.tsv"))
    required_environment_keys <- c(
      "reticulate", "StableMate", "python", "python_version",
      "StableMate_RemoteSha",
      paste0("python_module_", c(
        "numpy", "pandas", "sklearn", "torch", "torchvision", "kymatio",
        "pyod", "tabpfn", "tabicl", "tabdpt"
      ))
    )
    if (!all(c("key", "value") %in% names(environment)) ||
        !all(required_environment_keys %in% environment$key) ||
        sum(environment$key == "python") != 1L ||
        !grepl("^/home/kgs24/", environment[key == "python", value][[1L]]) ||
        environment[key == "reticulate", value][[1L]] != "1.45.0" ||
        environment[key == "StableMate", value][[1L]] != "0.1.0" ||
        environment[key == "StableMate_RemoteSha", value][[1L]] !=
          "cf0fe9f344c1756b6f89fd33f9fdae1a7f35b24a" ||
        any(!nzchar(environment[key %in% required_environment_keys, value])) ||
        !all(c("path", "bytes", "sha256") %in% names(models)) ||
        !nrow(models) || anyNA(models[, .(path, bytes, sha256)]) ||
        any(!nzchar(models$path)) || any(models$bytes < 1) ||
        any(!grepl("^[0-9a-f]{64}$", models$sha256)) ||
        anyDuplicated(models$path) || length(python_freeze) < 10L ||
        !any(grepl("^tabpfn==8\\.0\\.7$", python_freeze)) ||
        !any(grepl("^tabicl==2\\.1\\.1$", python_freeze)) ||
        !any(grepl("^tabdpt==1\\.2\\.0$", python_freeze)) ||
        nrow(r_packages[Package == "reticulate" & Version == "1.45.0"]) != 1L ||
        nrow(r_packages[Package == "StableMate" & Version == "0.1.0"]) != 1L) {
      stop("Transfer route is not bound to the reviewed Argos environment.",
           call. = FALSE)
    }
    model_names <- basename(models$path)
    if (any(!nzchar(model_names)) || any(model_names %in% c(".", "..")) ||
        anyDuplicated(model_names)) {
      stop("Transfer public-model manifest has unsafe or ambiguous model ",
           "names.", call. = FALSE)
    }
    model_root_path <- file.path(path, "public_models")
    if (!dir.exists(model_root_path) || nzchar(Sys.readlink(model_root_path))) {
      stop("Transfer public-model payload directory is missing or symbolic.",
           call. = FALSE)
    }
    model_root <- normalizePath(model_root_path, mustWork = TRUE)
    resolved_models <- file.path(model_root, model_names)
    if (any(!file.exists(resolved_models)) ||
        any(file.info(resolved_models)$isdir %in% TRUE) ||
        any(nzchar(Sys.readlink(resolved_models)))) {
      stop("Transfer public-model payload is missing or not a regular ",
           "immutable file.", call. = FALSE)
    }
    resolved_models <- normalizePath(resolved_models, mustWork = TRUE)
    if (any(dirname(resolved_models) != model_root)) {
      stop("Transfer public-model payload escapes its receipt directory.",
           call. = FALSE)
    }
    observed_model_bytes <- unname(file.info(resolved_models)$size)
    observed_model_sha <- vapply(resolved_models, .overlay_sha256,
                                 character(1L))
    if (any(!is.finite(as.numeric(models$bytes))) ||
        any(as.numeric(models$bytes) != as.numeric(observed_model_bytes)) ||
        any(models$sha256 != observed_model_sha)) {
      stop("Transfer public-model payload bytes do not match their manifest ",
           "size and SHA-256.", call. = FALSE)
    }
    runtime_assets <- data.table::data.table(
      path = resolved_models, bytes = observed_model_bytes,
      sha256 = observed_model_sha
    )[order(path)]
  } else {
    stop("Unknown execution route.", call. = FALSE)
  }
  list(files = files, sha256 = receipt_sha, runtime_assets = runtime_assets)
}

.overlay_copy_receipt <- function(receipt, destination) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(receipt$files))) {
    target <- file.path(destination, receipt$files$path[[i]])
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(receipt$files$source_path[[i]], target,
                   overwrite = FALSE, copy.mode = TRUE) ||
        !identical(.overlay_sha256(target), receipt$files$sha256[[i]])) {
      stop("Could not copy an environment receipt byte-for-byte.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.overlay_input_snapshot <- function(paths) {
  files <- character()
  for (path in paths) {
    if (!file.exists(path)) stop("Required input does not exist: ", path,
                                call. = FALSE)
    if (dir.exists(path)) {
      found <- list.files(path, recursive = TRUE, all.files = TRUE,
                          full.names = TRUE, include.dirs = FALSE, no.. = TRUE)
      found <- found[file.info(found)$isdir %in% FALSE]
      files <- c(files, found)
    } else {
      files <- c(files, path)
    }
  }
  files <- sort(unique(normalizePath(files, mustWork = TRUE)))
  out <- data.table::data.table(
    path = files, bytes = unname(file.info(files)$size),
    sha256 = vapply(files, .overlay_sha256, character(1L))
  )
  attr(out, "input_roots") <- paths
  out
}

.overlay_assert_sources_unchanged <- function(snapshot) {
  current <- .overlay_input_snapshot(attr(snapshot, "input_roots"))
  if (!identical(snapshot$path, current$path) ||
      !identical(as.numeric(snapshot$bytes), as.numeric(current$bytes)) ||
      !identical(snapshot$sha256, current$sha256)) {
    stop("One or more source inputs changed during candidate construction.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.overlay_write_table <- function(x, path, gzip = FALSE) {
  x <- .overlay_sanitize(x)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "",
                     compress = if (gzip) "gzip" else "none")
  invisible(path)
}

.overlay_manifest <- function(stage, pins, producer_sha) {
  files <- list.files(stage, recursive = TRUE, all.files = TRUE,
                      full.names = TRUE, include.dirs = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  relative <- substring(normalizePath(files, mustWork = TRUE),
                        nchar(paste0(normalizePath(stage),
                                    .Platform$file.sep)) + 1L)
  .overlay_safe_relative(relative, "Candidate manifest path")
  role <- ifelse(
    relative %in% c("within_bundle.rds", "within_cells.tsv",
                    "within_predictions.tsv.gz"), "integrated_within",
    ifelse(relative %in% c("transfer_bundle.rds", "transfer_cells.tsv",
                           "transfer_predictions.tsv.gz"),
           "integrated_transfer",
    ifelse(startsWith(relative, "splitter_manifests/"), "splitter_manifest",
    ifelse(startsWith(relative, "environment/"), "environment_receipt",
    ifelse(relative == "structural_witness.rds", "structural_witness",
    ifelse(startsWith(relative, "inputs/"), "frozen_task_table",
    ifelse(relative == "transfer_shard_manifest.tsv", "transfer_shard_manifest",
    ifelse(relative == "legacy_compatibility_overlay_input_manifest.tsv",
           "legacy_compatibility_manifest",
    ifelse(relative == "source_input_manifest.tsv", "source_input_manifest",
    ifelse(relative == "route_environment_map.tsv", "route_environment_map",
    ifelse(relative == "acceptance_checks.tsv", "acceptance_checks",
           "candidate_provenance")))))))))))
  manifest <- data.table::data.table(
    file = relative, role = role, bytes = unname(file.info(files)$size),
    sha256 = vapply(files, .overlay_sha256, character(1L)),
    run_mode = "full",
    producer_script_sha256 = producer_sha,
    producer_package_version = pins$producer_package_version,
    producer_package_commit = pins$producer_package_commit,
    base_within_bundle_sha256 = pins$base_within_sha,
    base_transfer_bundle_sha256 = pins$base_transfer_sha,
    within_environment_receipt_sha256 = pins$within_environment_sha,
    transfer_environment_receipt_sha256 = pins$transfer_environment_sha,
    execution_package_version = pins$package_version,
    execution_package_commit = pins$package_commit,
    cache_package_commit = pins$cache_commit,
    within_analysis_code_id = pins$within_code_id,
    transfer_analysis_code_id = pins$transfer_code_id,
    within_task_sha256 = pins$within_task_sha,
    transfer_task_sha256 = pins$transfer_task_sha,
    structural_witness_sha256 = pins$witness_sha,
    structural_cache_sha256 = pins$structural_cache_sha,
    structural_fold_engine_sha256 = pins$fold_engine_sha
  )
  data.table::setorder(manifest, file)
  manifest
}

.overlay_candidate_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .overlay_parse_args(args)
  producer <- .overlay_validate_producer_checkout(
    opt$expected_producer_package_version,
    opt$expected_producer_package_commit
  )
  input_files <- c(
    opt$base_within_bundle, opt$base_transfer_bundle, opt$within_tasks,
    opt$within_shard_root, opt$within_split_manifest_dir, opt$transfer_tasks,
    opt$transfer_shard_root, opt$within_environment_receipt,
    opt$transfer_environment_receipt, opt$structural_witness
  )
  normalized_inputs <- vapply(input_files, normalizePath, character(1L),
                              mustWork = TRUE)
  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) {
    stop("--output-dir already exists; candidate directories are immutable.",
         call. = FALSE)
  }
  output_prefix <- paste0(output_dir, .Platform$file.sep)
  input_prefix <- paste0(normalized_inputs, .Platform$file.sep)
  if (any(startsWith(output_dir, input_prefix)) ||
      any(startsWith(normalized_inputs, output_prefix)) ||
      startsWith(output_dir, paste0(producer$root, .Platform$file.sep))) {
    stop("--output-dir must be separate from every source input.",
         call. = FALSE)
  }
  receipt_input_files <- c(
    .overlay_receipt_files(opt$within_environment_receipt,
                           "within_s5")$source_path,
    .overlay_receipt_files(opt$transfer_environment_receipt,
                           "transfer_argos")$source_path
  )
  source_snapshot <- .overlay_input_snapshot(c(
    normalized_inputs[!normalized_inputs %in% c(
      normalizePath(opt$within_environment_receipt, mustWork = TRUE),
      normalizePath(opt$transfer_environment_receipt, mustWork = TRUE)
    )], receipt_input_files
  ))
  stage <- tempfile(paste0(".", basename(output_dir), "-stage-"),
                    tmpdir = output_parent)
  dir.create(stage, recursive = FALSE)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)

  pins <- list(
    package_version = opt$expected_execution_package_version,
    package_commit = opt$expected_execution_package_commit,
    cache_commit = opt$expected_cache_package_commit,
    within_code_id = opt$expected_within_code_id,
    transfer_code_id = opt$expected_transfer_code_id,
    within_task_sha = opt$expected_within_task_sha256,
    transfer_task_sha = opt$expected_transfer_task_sha256,
    witness_sha = opt$expected_structural_witness_sha256,
    structural_cache_sha = opt$expected_structural_cache_sha256,
    fold_engine_sha = opt$expected_structural_fold_engine_sha256,
    producer_package_version = opt$expected_producer_package_version,
    producer_package_commit = opt$expected_producer_package_commit,
    base_within_sha = opt$expected_base_within_bundle_sha256,
    base_transfer_sha = opt$expected_base_transfer_bundle_sha256,
    within_environment_sha = opt$expected_within_environment_receipt_sha256,
    transfer_environment_sha = opt$expected_transfer_environment_receipt_sha256
  )
  within_tasks <- .overlay_read_within_tasks(
    opt$within_tasks, pins$within_task_sha
  )
  transfer_tasks <- .overlay_read_transfer_tasks(
    opt$transfer_tasks, pins$transfer_task_sha
  )
  within_pins <- list(
    package_version = pins$package_version,
    package_commit = pins$package_commit,
    cache_commit = pins$cache_commit, code_id = pins$within_code_id,
    task_sha = pins$within_task_sha,
    structural_cache_sha = pins$structural_cache_sha,
    fold_engine_sha = pins$fold_engine_sha
  )
  transfer_pins <- list(
    package_version = pins$package_version,
    package_commit = pins$package_commit,
    cache_commit = pins$cache_commit, code_id = pins$transfer_code_id,
    task_sha = pins$transfer_task_sha
  )
  manifests <- .overlay_validate_split_manifests(
    opt$within_split_manifest_dir, opt$within_shard_root, within_tasks,
    within_pins
  )
  within <- .overlay_read_within_shards(manifests, within_pins)
  if (nrow(within$cells) != 2210L ||
      anyDuplicated(within$cells[, .(method, seed, cohort)])) {
    stop("Within shards do not contain exactly 2,210 unique cells.",
         call. = FALSE)
  }
  transfer <- .overlay_read_transfer_shards(
    opt$transfer_shard_root, transfer_tasks, transfer_pins
  )
  .overlay_validate_transfer_science(transfer$cells, transfer$predictions)
  witness_sha <- .overlay_sha256(opt$structural_witness)
  if (!identical(witness_sha, pins$witness_sha)) {
    stop("Structural witness does not match its expected SHA-256.",
         call. = FALSE)
  }
  witness <- readRDS(opt$structural_witness)
  .overlay_validate_reassembled_within(
    within, within_tasks, manifests, witness, witness_sha, within_pins
  )
  .overlay_validate_witness_bindings(
    witness, witness_sha, manifests, within$cells, within$witnesses,
    within_pins
  )
  within_receipt <- .overlay_validate_environment_receipt(
    opt$within_environment_receipt, "within_s5", pins$package_commit,
    pins$package_version, pins$within_environment_sha
  )
  transfer_receipt <- .overlay_validate_environment_receipt(
    opt$transfer_environment_receipt, "transfer_argos", pins$package_commit,
    pins$package_version, pins$transfer_environment_sha
  )
  runtime_asset_snapshot <- .overlay_input_snapshot(c(
    within_receipt$runtime_assets$path,
    transfer_receipt$runtime_assets$path
  ))

  if (!identical(.overlay_sha256(opt$base_within_bundle),
                 pins$base_within_sha) ||
      !identical(.overlay_sha256(opt$base_transfer_bundle),
                 pins$base_transfer_sha)) {
    stop("Base within/transfer bundle bytes do not match their exact SHA-256 ",
         "pins.", call. = FALSE)
  }
  within_base <- readRDS(opt$base_within_bundle)
  transfer_base <- readRDS(opt$base_transfer_bundle)
  integrated_within <- .overlay_integrate(
    within_base, within$cells, within$predictions, .overlay_within_methods,
    c("seed", "method", "cohort"),
    c("seed", "method", "cohort", "fold", "sample_idx"), 8500L,
    "within", within_pins, "2.6.5_dependency_backed_exact"
  )
  integrated_transfer <- .overlay_integrate(
    transfer_base, transfer$cells, transfer$predictions,
    .overlay_transfer_methods, c("disease", "held_out", "seed", "method"),
    c("disease", "held_out", "seed", "method", "sample_id"), 2160L,
    "transfer", transfer_pins, "2.6.5_dependency_backed_exact"
  )
  .overlay_validate_predictions(
    integrated_within$cells, integrated_within$predictions,
    c("seed", "method", "cohort"), "sample_idx", "n_eff", "Within"
  )
  .overlay_validate_predictions(
    integrated_transfer$cells, integrated_transfer$predictions,
    c("disease", "held_out", "seed", "method"), "sample_id", "n_test",
    "Transfer"
  )

  integrated_within$cells <- .overlay_sanitize(integrated_within$cells)
  integrated_within$predictions <- .overlay_sanitize(
    integrated_within$predictions
  )
  integrated_within$bundle$cells <- integrated_within$cells
  integrated_within$bundle$predictions <- integrated_within$predictions
  integrated_transfer$cells <- .overlay_sanitize(integrated_transfer$cells)
  integrated_transfer$predictions <- .overlay_sanitize(
    integrated_transfer$predictions
  )
  integrated_transfer$bundle$cells <- integrated_transfer$cells
  integrated_transfer$bundle$predictions <- integrated_transfer$predictions

  saveRDS(integrated_within$bundle, file.path(stage, "within_bundle.rds"),
          version = 3L)
  .overlay_write_table(integrated_within$cells,
                       file.path(stage, "within_cells.tsv"))
  .overlay_write_table(integrated_within$predictions,
                       file.path(stage, "within_predictions.tsv.gz"), TRUE)
  saveRDS(integrated_transfer$bundle, file.path(stage, "transfer_bundle.rds"),
          version = 3L)
  .overlay_write_table(integrated_transfer$cells,
                       file.path(stage, "transfer_cells.tsv"))
  .overlay_write_table(integrated_transfer$predictions,
                       file.path(stage, "transfer_predictions.tsv.gz"), TRUE)

  dir.create(file.path(stage, "inputs"))
  if (!file.copy(opt$within_tasks, file.path(stage, "inputs", "within_tasks.tsv")) ||
      !file.copy(opt$transfer_tasks,
                 file.path(stage, "inputs", "transfer_tasks.tsv"))) {
    stop("Could not copy frozen task tables.", call. = FALSE)
  }
  if (!identical(.overlay_sha256(file.path(stage, "inputs", "within_tasks.tsv")),
                 pins$within_task_sha) ||
      !identical(.overlay_sha256(file.path(stage, "inputs", "transfer_tasks.tsv")),
                 pins$transfer_task_sha)) {
    stop("Copied task tables do not match their exact SHA-256 pins.",
         call. = FALSE)
  }
  dir.create(file.path(stage, "splitter_manifests"))
  split_names <- unique(manifests$manifest_name)
  for (name in split_names) {
    source <- file.path(opt$within_split_manifest_dir, name)
    if (!file.copy(source, file.path(stage, "splitter_manifests", name)) ||
        !identical(.overlay_sha256(source),
                   .overlay_sha256(file.path(stage, "splitter_manifests", name)))) {
      stop("Could not copy a splitter manifest byte-for-byte.", call. = FALSE)
    }
  }
  if (!file.copy(opt$structural_witness,
                 file.path(stage, "structural_witness.rds")) ||
      !identical(witness_sha,
                 .overlay_sha256(file.path(stage, "structural_witness.rds")))) {
    stop("Could not copy structural witness byte-for-byte.", call. = FALSE)
  }
  .overlay_copy_receipt(within_receipt,
                        file.path(stage, "environment", "within_s5"))
  .overlay_copy_receipt(transfer_receipt,
                        file.path(stage, "environment", "transfer_argos"))
  .overlay_write_table(transfer$manifest,
                       file.path(stage, "transfer_shard_manifest.tsv"))
  compatibility <- data.table::data.table(
    path = c(
      normalizePath(c(opt$within_tasks, opt$transfer_tasks), mustWork = TRUE),
      manifests$resolved_path,
      normalizePath(file.path(
        opt$transfer_shard_root, transfer$manifest$path
      ), mustWork = TRUE),
      transfer_receipt$files$source_path
    )
  )
  compatibility[, sha256 := vapply(path, .overlay_sha256, character(1L))]
  if (nrow(compatibility) != 2251L || anyDuplicated(compatibility$path) ||
      any(!file.exists(compatibility$path))) {
    stop("Legacy compatibility manifest is not the exact 2 + 2210 + 35 + 4 ",
         "input closure.", call. = FALSE)
  }
  .overlay_write_table(
    compatibility,
    file.path(stage, "legacy_compatibility_overlay_input_manifest.tsv")
  )
  source_inputs <- data.table::rbindlist(list(
    data.table::copy(source_snapshot), data.table::copy(runtime_asset_snapshot)
  ))
  source_inputs <- unique(source_inputs, by = "path")
  source_inputs[, input_role := "source_closure"]
  source_inputs[grepl("transfer_bundle\\.rds$", path),
                input_role := "transfer_shard"]
  source_inputs[grepl("within_bundle\\.rds$", path),
                input_role := "within_shard"]
  source_inputs[grepl("output_manifest__", path),
                input_role := "within_splitter_manifest"]
  source_inputs[path %in% transfer_receipt$files$source_path,
                input_role := "transfer_environment_receipt"]
  source_inputs[path %in% within_receipt$files$source_path,
                input_role := "within_environment_receipt"]
  source_inputs[path %in% within_receipt$runtime_assets$path,
                input_role := "within_runtime_asset"]
  source_inputs[path %in% transfer_receipt$runtime_assets$path,
                input_role := "transfer_runtime_asset"]
  source_inputs[path == normalizePath(opt$structural_witness),
                input_role := "structural_witness"]
  source_inputs[path == normalizePath(opt$transfer_tasks),
                input_role := "transfer_task_table"]
  source_inputs[path == normalizePath(opt$within_tasks),
                input_role := "within_task_table"]
  source_inputs[path == normalizePath(opt$base_transfer_bundle),
                input_role := "base_transfer_bundle"]
  source_inputs[path == normalizePath(opt$base_within_bundle),
                input_role := "base_within_bundle"]
  data.table::setcolorder(source_inputs,
                          c("input_role", "path", "bytes", "sha256"))
  data.table::setorder(source_inputs, input_role, path)
  .overlay_write_table(source_inputs,
                       file.path(stage, "source_input_manifest.tsv"))
  route_map <- data.table::data.table(
    execution_route = rep("2.6.5_dependency_backed_exact", 2L),
    analysis_surface = c("within", "transfer"),
    environment_route = c("within_s5", "transfer_argos"),
    environment_receipt_path = c("environment/within_s5",
                                 "environment/transfer_argos"),
    environment_receipt_sha256 = c(within_receipt$sha256,
                                   transfer_receipt$sha256),
    execution_package_version = pins$package_version,
    execution_package_commit = pins$package_commit,
    cache_package_commit = pins$cache_commit,
    analysis_code_id = c(pins$within_code_id, pins$transfer_code_id)
  )
  .overlay_write_table(route_map,
                       file.path(stage, "route_environment_map.tsv"))

  checks <- data.table::data.table(
    check_id = c(
      "within_grid_13x34x5", "within_13_splitter_manifests",
      "within_manifest_shard_hashes", "within_task_key_equality",
      "transfer_exact_630_cells", "integrated_within_8500_cells",
      "integrated_transfer_2160_cells", "unrelated_base_rows_unchanged",
      "eligible_prediction_accounting", "deepmaha_exact_nine_witness_cells",
      "deepmaha_embedded_witness_equality",
      "route_specific_environment_receipts", "source_inputs_immutable",
      "producer_exact_clean_checkout", "base_bundles_exact_sha256",
      "durable_source_input_manifest"
    ),
    pass = TRUE,
    detail = c(
      "2210 cells", "13 manifests; 170 rows each", "all bytes and SHA-256 match",
      "combined keys equal frozen within tasks", "630 cells",
      "8500 cells", "2160 cells", "within and transfer cells/predictions",
      "eligible cells only; exact profile counts", "9 cells only",
      "embedded objects equal source witness", "s5 within; Argos transfer",
      "within and transfer runtime/model assets verified before promotion",
      "tracked script and DESCRIPTION match exact clean source commit",
      "within and transfer base pins matched",
      "all source shards, receipts, assets, witness, and bases recorded"
    )
  )
  .overlay_assert_sources_unchanged(source_snapshot)
  .overlay_assert_sources_unchanged(runtime_asset_snapshot)
  .overlay_write_table(checks, file.path(stage, "acceptance_checks.tsv"))

  producer_end <- .overlay_validate_producer_checkout(
    opt$expected_producer_package_version,
    opt$expected_producer_package_commit
  )
  if (!identical(producer$script_sha256, producer_end$script_sha256) ||
      !identical(producer$script_blob, producer_end$script_blob)) {
    stop("Producer checkout changed during candidate construction.",
         call. = FALSE)
  }
  manifest <- .overlay_manifest(stage, pins, producer$script_sha256)
  .overlay_write_table(manifest, file.path(stage, "manifest.tsv"))
  candidate_files <- setdiff(
    list.files(stage, recursive = TRUE, all.files = TRUE, full.names = FALSE,
               include.dirs = FALSE, no.. = TRUE), "manifest.tsv"
  )
  if (!setequal(candidate_files, manifest$file) ||
      any(manifest$sha256 != vapply(file.path(stage, manifest$file),
                                    .overlay_sha256, character(1L)))) {
    stop("Candidate manifest does not cover every current artifact exactly.",
         call. = FALSE)
  }
  .overlay_assert_sources_unchanged(source_snapshot)
  .overlay_assert_sources_unchanged(runtime_asset_snapshot)
  if (!file.rename(stage, output_dir)) {
    stop("Could not atomically promote the dependency overlay candidate.",
         call. = FALSE)
  }
  cat(sprintf(
    "PASS: atomically promoted dependency overlay candidate (%d within, %d transfer cells).\n",
    nrow(integrated_within$cells), nrow(integrated_transfer$cells)
  ))
  invisible(output_dir)
}

if (sys.nframe() == 0L) .overlay_candidate_main()
