#!/usr/bin/env Rscript

# Build the immutable 34-unit x 5-seed task grid for the paper-3 nested
# single-sample selector analysis. This producer validates the frozen cache but
# never summarizes outcomes or inspects performance.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.nested_task_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Arguments require --name=value: ", arg)
    key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg), fixed = TRUE)
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  required <- c("cache", "cache_sha256", "output_dir", "package_commit")
  missing <- required[!required %in% names(out)]
  if (length(missing)) stop("Missing: --", paste(missing, collapse = ", --"))
  if (!grepl("^[0-9a-f]{64}$", tolower(out$cache_sha256)) ||
      !grepl("^[0-9a-f]{40}$", tolower(out$package_commit))) {
    stop("cache_sha256/package_commit must be exact 64/40-character hex pins.")
  }
  out
}

.nested_expected_units <- function() {
  c(
    "TORAY_PUBLIC_CANONICAL", "GSE109319", "GSE117063", "GSE122488",
    "GSE130654", "GSE137109", "GSE145176", "GSE186595", "GSE188232",
    "GSE188627", "GSE198834", "GSE200489", "GSE204951", "GSE210329",
    "GSE210546", "GSE221088", "GSE227778", "GSE259327_DIAGNOSTIC",
    "GSE266859", "GSE270497", "GSE279209", "GSE304572", "GSE31568",
    "GSE320031", "GSE41922", "GSE50013", "GSE64591", "GSE65071",
    "GSE70080", "GSE71008", "GSE83977", "GSE85589",
    "GSE94533_PUBLIC", "PRJEB30542"
  )
}

.nested_expected_k3 <- function() {
  c("GSE117063", "GSE130654", "GSE200489", "GSE204951", "GSE210546",
    "GSE227778", "GSE266859", "GSE83977")
}

.nested_validate_cache <- function(cache) {
  expected <- .nested_expected_units()
  if (!is.list(cache) || !identical(names(cache), expected)) {
    stop("Cache does not contain the exact frozen 34-unit roster and order.")
  }
  required <- c(
    "expr_per_sample", "y_bin", "group_id", "outer_k",
    "provenance_block", "source_accessions", "modality", "biospecimen"
  )
  rows <- lapply(names(cache), function(unit_id) {
    z <- cache[[unit_id]]
    missing <- setdiff(required, names(z))
    if (length(missing)) {
      stop(unit_id, " is missing: ", paste(missing, collapse = ", "))
    }
    X <- as.matrix(z$expr_per_sample)
    y <- as.integer(z$y_bin)
    group_id <- as.character(z$group_id)
    sample_id <- rownames(X)
    if (!is.numeric(X) || is.null(sample_id) || anyDuplicated(sample_id) ||
        nrow(X) != length(y) || length(y) != length(group_id) ||
        anyNA(X) || anyNA(y) || anyNA(group_id) ||
        any(!nzchar(group_id)) || !identical(sort(unique(y)), 0:1)) {
      stop(unit_id, " has invalid matrix, labels, samples, or groups.")
    }
    group_y <- split(y, group_id)
    if (any(vapply(group_y, function(v) length(unique(v)) != 1L,
                   logical(1L)))) {
      stop(unit_id, " contains a biological group with discordant labels.")
    }
    expected_k <- if (unit_id %in% .nested_expected_k3()) 3L else 5L
    if (!identical(as.integer(z$outer_k), expected_k)) {
      stop(unit_id, " violates the frozen outer-fold contract.")
    }
    source_accessions <- unique(as.character(z$source_accessions))
    modality <- as.character(z$modality)
    biospecimen <- as.character(z$biospecimen)
    if (!length(source_accessions) || anyNA(source_accessions) ||
        any(!nzchar(source_accessions)) ||
        !is.character(z$provenance_block) || length(z$provenance_block) != 1L ||
        is.na(z$provenance_block) || !nzchar(z$provenance_block) ||
        length(modality) != 1L || is.na(modality) || !nzchar(modality) ||
        length(biospecimen) != 1L || is.na(biospecimen) || !nzchar(biospecimen)) {
      stop(unit_id, " has incomplete provenance fields.")
    }
    primary_unit <- unit_id != "GSE31568"
    whole_blood <- grepl("whole[ -]?blood", biospecimen, ignore.case = TRUE)
    if (identical(primary_unit, whole_blood)) {
      stop(unit_id, " violates the frozen primary/whole-blood sensitivity ",
           "contract.")
    }
    data.table(
      unit_id = unit_id, outer_k = expected_k, n_profiles = nrow(X),
      n_groups = uniqueN(group_id), n_cases = sum(y == 1L),
      n_controls = sum(y == 0L), n_features = ncol(X),
      primary_unit = primary_unit, modality = modality,
      biospecimen = biospecimen,
      provenance_block = z$provenance_block,
      source_accessions = paste(source_accessions, collapse = ";")
    )
  })
  units <- rbindlist(rows)
  if (sum(units$primary_unit) != 33L ||
      !identical(units[primary_unit == FALSE, unit_id], "GSE31568")) {
    stop("Cache does not reproduce the frozen 33-primary/one-sensitivity contract.")
  }
  units
}

.nested_task_grid <- function(units) {
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  rows <- lapply(seq_len(nrow(units)), function(i) {
    data.table(
      unit_id = units$unit_id[[i]], seed = seeds,
      label = sprintf("%02d_%s_seed%d", i, units$unit_id[[i]], seeds),
      outer_k = units$outer_k[[i]]
    )
  })
  out <- rbindlist(rows)
  out[, task_id := seq_len(.N)]
  setcolorder(out, c("task_id", "unit_id", "seed", "label", "outer_k"))
  if (nrow(out) != 170L || anyDuplicated(out[, .(unit_id, seed)]) ||
      !identical(out$task_id, seq_len(170L))) {
    stop("Nested-selector task grid is not the exact 34 x 5 contract.")
  }
  out
}

.nested_code_provenance <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
  } else normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
  package_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                                mustWork = TRUE)
  git <- Sys.which("git")
  if (!nzchar(git)) stop("Could not locate git.")
  commit <- system2(git, c("-C", package_root, "rev-parse", "HEAD"),
                    stdout = TRUE, stderr = TRUE)
  status <- system2(git, c("-C", package_root, "status", "--porcelain=v1",
                           "--untracked-files=all"),
                    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(commit, "status")) || length(commit) != 1L ||
      !grepl("^[0-9a-f]{40}$", commit) || !is.null(attr(status, "status"))) {
    stop("Could not resolve task-producer code provenance.")
  }
  list(commit = unname(commit[[1L]]), dirty = length(status) > 0L,
       script_path = script_path,
       script_sha256 = digest(script_path, algo = "sha256", file = TRUE))
}

.nested_task_builder_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .nested_task_args(args)
  code <- .nested_code_provenance()
  if (code$dirty || !identical(tolower(opt$package_commit), code$commit)) {
    stop("Task production requires the exact clean package commit pin.")
  }
  cache_path <- normalizePath(opt$cache, mustWork = TRUE)
  before <- digest(cache_path, algo = "sha256", file = TRUE)
  if (!identical(before, tolower(opt$cache_sha256))) {
    stop("Cache does not match --cache-sha256.")
  }
  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) stop("Output already exists: ", output_dir)
  temp_dir <- tempfile("nested-task-grid-", tmpdir = output_parent)
  dir.create(temp_dir)
  on.exit(if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE), add = TRUE)

  cache <- readRDS(cache_path)
  units <- .nested_validate_cache(cache)
  tasks <- .nested_task_grid(units)
  units[, `:=`(cache_sha256 = before, package_commit = code$commit)]
  fwrite(units, file.path(temp_dir, "nested_selector_units.tsv"), sep = "\t")
  fwrite(tasks, file.path(temp_dir, "nested_selector_tasks.tsv"), sep = "\t")
  after <- digest(cache_path, algo = "sha256", file = TRUE)
  final_code <- .nested_code_provenance()
  if (!identical(before, after) || !identical(code, final_code)) {
    stop("Cache or package code changed during task production.")
  }
  manifest <- data.table(
    cache_path = cache_path, cache_sha256 = before,
    n_units = nrow(units), n_tasks = nrow(tasks), seeds = "101;202;303;404;505",
    package_version = "2.6.5.9000", package_commit = code$commit,
    producer_path = code$script_path, producer_sha256 = code$script_sha256,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  fwrite(manifest, file.path(temp_dir, "nested_selector_input_manifest.tsv"),
         sep = "\t")
  if (!file.rename(temp_dir, output_dir)) {
    stop("Could not atomically promote task-grid output.")
  }
  cat("Created nested-selector grid:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .nested_task_builder_main()
