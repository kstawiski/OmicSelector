#!/usr/bin/env Rscript

# Validate and shard one frozen 170-cell optional-dependency benchmark bundle.
# This is a serialization-only utility: it never fits, scores, or otherwise
# recomputes a scientific result.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.dependency_split_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  value <- function(flag) {
    hit <- match(flag, args)
    if (is.na(hit) || hit == length(args)) NA_character_ else args[[hit + 1L]]
  }
  required <- c(
    bundle = "--bundle", tasks = "--tasks", output_root = "--output-root",
    expected_method = "--expected-method",
    expected_package_version = "--expected-package-version",
    expected_package_commit = "--expected-package-commit",
    expected_cache_commit = "--expected-cache-commit",
    expected_analysis_code_id = "--expected-analysis-code-id",
    expected_task_sha256 = "--expected-task-sha256"
  )
  out <- lapply(required, value)
  missing <- names(out)[vapply(out, function(x) is.na(x) || !nzchar(x), logical(1L))]
  if (length(missing)) {
    stop("Missing required argument(s): ",
         paste(unname(required[missing]), collapse = ", "), call. = FALSE)
  }
  out$allow_identical_overlap <- "--allow-identical-overlap" %in% args
  if (!grepl("^[0-9a-f]{40}$", out$expected_package_commit) ||
      !grepl("^[0-9a-f]{40}$", out$expected_cache_commit) ||
      !grepl("^[0-9a-f]{16}$", out$expected_analysis_code_id) ||
      !grepl("^[0-9a-f]{64}$", out$expected_task_sha256)) {
    stop("Expected package/cache commits, analysis code id, and task hash ",
         "must be exact lowercase hexadecimal pins (40/40/16/64 characters).",
         call. = FALSE)
  }
  out
}

.dependency_assert_scalar <- function(observed, expected, name) {
  if (!is.character(observed) || length(observed) != 1L || is.na(observed) ||
      !identical(observed, expected)) {
    stop("Bundle ", name, " does not match its expected scalar pin.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dependency_read_tasks <- function(path) {
  tasks <- data.table::fread(
    path, sep = "\t", header = FALSE, quote = "",
    col.names = c("method", "seed", "cohort", "label"),
    colClasses = c("character", "integer", "character", "character")
  )
  if (nrow(tasks) != 2210L || ncol(tasks) != 4L || anyNA(tasks) ||
      any(!nzchar(tasks$method)) || any(!nzchar(tasks$cohort)) ||
      any(!nzchar(tasks$label)) || anyDuplicated(tasks[, .(method, seed, cohort)]) ||
      anyDuplicated(tasks$label)) {
    stop("Frozen task table must contain 2,210 unique, complete task rows.",
         call. = FALSE)
  }
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  expected_methods <- c(
    "ai-scarf", "ai-tabpfn", "coda-codacore", "coda-deepcoda",
    "img-gasfcnn", "inv-scatter", "lrt-deepmaha", "nc-ecod-copod",
    "proto-net", "ssl-vicreg", "tab-tabdpt", "tab-tabicl", "unc-sngp"
  )
  methods <- unique(tasks$method)
  cohorts <- unique(tasks$cohort)
  if (!setequal(methods, expected_methods) || length(cohorts) != 34L ||
      !setequal(tasks$seed, seeds) ||
      any(tasks[, .N, by = method]$N != 170L) ||
      any(tasks[, .N, by = .(method, seed)]$N != 34L) ||
      any(tasks[, .N, by = .(method, cohort)]$N != 5L)) {
    stop("Frozen task table is not the complete 13-by-five-by-34 design.",
         call. = FALSE)
  }
  expected_label <- sprintf("%s__seed_%d__%s", tasks$method, tasks$seed,
                            tasks$cohort)
  if (!identical(tasks$label, expected_label) ||
      any(!grepl("^[A-Za-z0-9_.-]+$", tasks$label))) {
    stop("Frozen task labels are non-canonical or unsafe as directory names.",
         call. = FALSE)
  }
  data.table::setorder(tasks, method, seed, cohort)
  tasks
}

.dependency_validate_pin_column <- function(x, column, expected, surface) {
  if (!column %in% names(x) || anyNA(x[[column]]) ||
      any(x[[column]] != expected)) {
    stop(surface, " does not carry the expected ", column, " pin.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dependency_validate_ineligible_cells <- function(cells) {
  ineligible <- cells[eligible == FALSE]
  if (!nrow(ineligible)) return(invisible(TRUE))

  gse83977_constant_reason <- paste0(
    "n_valid_folds=0 < 3; per-fold: ",
    paste(rep("method/baseline constant or all-NA on test fold", 3L),
          collapse = "; ")
  )
  reviewed_structural <-
    ineligible$cohort == "GSE83977" &
    ineligible$seed %in% c(101L, 202L, 303L, 404L, 505L) &
    ineligible$n_valid_folds == 0L &
    ineligible$outer_k == 3L &
    ineligible$n_group_split_violations == 0L &
    ineligible$ineligible_reason == gse83977_constant_reason
  reviewed_structural[is.na(reviewed_structural)] <- FALSE

  if (any(!reviewed_structural)) {
    rejected <- unique(ineligible$ineligible_reason[!reviewed_structural])
    rejected <- rejected[!is.na(rejected) & nzchar(rejected)]
    preview <- if (length(rejected)) rejected[[1L]] else "missing reason"
    stop(
      "Bundle contains an unreviewed ineligibility, fit/score failure, or ",
      "dependency/runtime failure: ", preview,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.dependency_validate_predictions <- function(cells, predictions, expected) {
  required <- c(
    "fold", "n_train", "n_test", "sample_idx", "sample_id", "group_id",
    "y", "score_m", "score_b", "method", "cohort", "seed",
    "package_commit", "analysis_code_id"
  )
  if (!is.data.frame(predictions)) {
    stop("Bundle predictions must be a data.frame.", call. = FALSE)
  }
  if (!nrow(predictions)) {
    if (any(cells$eligible)) {
      stop("Eligible cells have no predictions.", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  missing <- setdiff(required, names(predictions))
  if (length(missing)) {
    stop("Predictions lack required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  predictions <- data.table::as.data.table(predictions)
  .dependency_validate_pin_column(predictions, "package_commit",
                                  expected$package_commit, "Predictions")
  .dependency_validate_pin_column(predictions, "analysis_code_id",
                                  expected$analysis_code_id, "Predictions")
  if (anyNA(predictions[, .(fold, n_train, n_test, sample_idx, sample_id,
                            group_id, y, score_m, score_b, method, cohort, seed)]) ||
      any(!nzchar(predictions$sample_id)) || any(!nzchar(predictions$group_id)) ||
      any(!is.finite(predictions$score_m)) ||
      any(!is.finite(predictions$score_b)) ||
      any(!predictions$y %in% c(0L, 1L)) ||
      any(predictions$n_train < 1L) || any(predictions$n_test < 1L) ||
      any(predictions$fold < 1L)) {
    stop("Predictions contain missing, non-finite, or invalid values.",
         call. = FALSE)
  }
  key <- c("method", "cohort", "seed", "sample_idx")
  if (anyDuplicated(predictions[, ..key])) {
    stop("Predictions contain duplicate held-out profile keys.", call. = FALSE)
  }
  pred_cells <- unique(predictions[, .(method, cohort, seed)])
  eligible_cells <- cells[eligible == TRUE, .(method, cohort, seed)]
  if (!setequal(do.call(paste, pred_cells), do.call(paste, eligible_cells))) {
    stop("Prediction keys do not correspond exactly to eligible cells.",
         call. = FALSE)
  }
  counts <- predictions[, .(n_predictions = .N), by = .(method, cohort, seed)]
  accounting <- merge(
    cells[eligible == TRUE, .(method, cohort, seed, n_eff, n_valid_folds,
                              outer_k)],
    counts, by = c("method", "cohort", "seed"), all.x = TRUE, sort = FALSE
  )
  if (anyNA(accounting$n_predictions) ||
      any(accounting$n_predictions != accounting$n_eff) ||
      any(accounting$n_valid_folds != accounting$outer_k)) {
    stop("Eligible prediction cardinality does not match n_eff/outer-K.",
         call. = FALSE)
  }
  fold_counts <- predictions[, .(
    rows = .N, declared_n_test = list(unique(n_test)),
    declared_n_train = list(unique(n_train))
  ), by = .(method, cohort, seed, fold)]
  if (any(lengths(fold_counts$declared_n_test) != 1L) ||
      any(lengths(fold_counts$declared_n_train) != 1L) ||
      any(vapply(fold_counts$declared_n_test, `[[`, numeric(1L), 1L) !=
          fold_counts$rows)) {
    stop("Fold prediction counts disagree with their declared train/test sizes.",
         call. = FALSE)
  }
  fold_sets <- predictions[, .(observed_folds = list(sort(unique(fold)))),
                           by = .(method, cohort, seed)]
  fold_sets <- merge(fold_sets,
                     cells[, .(method, cohort, seed, outer_k)],
                     by = c("method", "cohort", "seed"), all.x = TRUE,
                     sort = FALSE)
  if (any(!vapply(seq_len(nrow(fold_sets)), function(i) {
    identical(as.integer(fold_sets$observed_folds[[i]]),
              seq_len(as.integer(fold_sets$outer_k[[i]])))
  }, logical(1L)))) {
    stop("Eligible predictions do not cover every prespecified outer fold.",
         call. = FALSE)
  }
  group_folds <- predictions[, .(n_folds = uniqueN(fold)),
                             by = .(method, cohort, seed, group_id)]
  if (any(group_folds$n_folds != 1L)) {
    stop("A provenance group crosses held-out folds.", call. = FALSE)
  }
  group_labels <- predictions[, .(n_labels = uniqueN(y)),
                              by = .(method, cohort, seed, group_id)]
  if (any(group_labels$n_labels != 1L)) {
    stop("A provenance group has inconsistent outcome labels.", call. = FALSE)
  }
  invisible(TRUE)
}

.dependency_validate_bundle <- function(bundle, tasks, expected) {
  required_bundle <- c(
    "cells", "predictions", "seeds", "methods", "cohorts",
    "package_version", "package_commit", "cache_package_commit",
    "analysis_code_id"
  )
  missing_bundle <- setdiff(required_bundle, names(bundle))
  if (!is.list(bundle) || length(missing_bundle)) {
    stop("Bundle lacks required field(s): ", paste(missing_bundle, collapse = ", "),
         call. = FALSE)
  }
  .dependency_assert_scalar(bundle$package_version, expected$package_version,
                            "package_version")
  .dependency_assert_scalar(bundle$package_commit, expected$package_commit,
                            "package_commit")
  .dependency_assert_scalar(bundle$cache_package_commit, expected$cache_commit,
                            "cache_package_commit")
  .dependency_assert_scalar(bundle$analysis_code_id, expected$analysis_code_id,
                            "analysis_code_id")
  if (!is.data.frame(bundle$cells)) {
    stop("Bundle cells must be a data.frame.", call. = FALSE)
  }
  cells <- data.table::as.data.table(bundle$cells)
  required_cells <- c(
    "method", "cohort", "seed", "eligible", "ineligible_reason", "n_eff",
    "n_valid_folds", "outer_k", "n_group_split_violations",
    "package_version", "package_commit", "cache_package_commit",
    "analysis_code_id"
  )
  missing_cells <- setdiff(required_cells, names(cells))
  if (length(missing_cells)) {
    stop("Cells lack required column(s): ", paste(missing_cells, collapse = ", "),
         call. = FALSE)
  }
  method_tasks <- tasks[method == expected$method]
  if (nrow(method_tasks) != 170L) {
    stop("Expected method does not own exactly 170 frozen tasks.", call. = FALSE)
  }
  if (nrow(cells) != 170L ||
      anyDuplicated(cells[, .(method, seed, cohort)]) ||
      !identical(unique(cells$method), expected$method) ||
      !setequal(cells$seed, c(101L, 202L, 303L, 404L, 505L)) ||
      uniqueN(cells$cohort) != 34L ||
      !setequal(do.call(paste, cells[, .(method, seed, cohort)]),
                do.call(paste, method_tasks[, .(method, seed, cohort)]))) {
    stop("Bundle cells are not the exact five-seed-by-34-cohort method grid.",
         call. = FALSE)
  }
  if (!identical(as.character(bundle$methods), expected$method) ||
      !setequal(as.integer(bundle$seeds), c(101L, 202L, 303L, 404L, 505L)) ||
      !setequal(as.character(bundle$cohorts), method_tasks$cohort)) {
    stop("Bundle method/seed/cohort metadata disagree with its cell grid.",
         call. = FALSE)
  }
  .dependency_validate_pin_column(cells, "package_version",
                                  expected$package_version, "Cells")
  .dependency_validate_pin_column(cells, "package_commit",
                                  expected$package_commit, "Cells")
  .dependency_validate_pin_column(cells, "cache_package_commit",
                                  expected$cache_commit, "Cells")
  .dependency_validate_pin_column(cells, "analysis_code_id",
                                  expected$analysis_code_id, "Cells")
  if (!is.logical(cells$eligible) || anyNA(cells$eligible) ||
      any(cells[eligible == TRUE, is.na(n_group_split_violations) |
                  n_group_split_violations != 0L]) ||
      any(cells[eligible == FALSE,
                !is.na(n_group_split_violations) &
                  n_group_split_violations != 0L]) ||
      any(grepl("driver[[:space:]]+cell[[:space:]]+error:|cell[[:space:]]+error:",
                cells$ineligible_reason, ignore.case = TRUE), na.rm = TRUE)) {
    stop("Cells contain a driver/cell error or an unproven/group-split result.",
         call. = FALSE)
  }
  if (any(cells[eligible == TRUE,
                is.na(n_eff) | n_eff < 1L | is.na(n_valid_folds) |
                  is.na(outer_k) | n_valid_folds != outer_k]) ||
      any(cells[eligible == TRUE,
                !is.na(ineligible_reason) & nzchar(ineligible_reason)]) ||
      any(cells[eligible == FALSE,
                is.na(ineligible_reason) | !nzchar(ineligible_reason)])) {
    stop("Eligible/ineligible cell accounting is internally inconsistent.",
         call. = FALSE)
  }
  .dependency_validate_ineligible_cells(cells)
  .dependency_validate_predictions(cells, bundle$predictions, expected)
  invisible(list(cells = cells, predictions = data.table::as.data.table(bundle$predictions),
                 method_tasks = method_tasks))
}

.dependency_shard_object <- function(bundle, cell, predictions, expected) {
  out <- bundle
  out$cells <- data.table::copy(cell)
  out$predictions <- data.table::copy(predictions)
  out$seeds <- as.integer(cell$seed[[1L]])
  out$methods <- expected$method
  out$cohorts <- as.character(cell$cohort[[1L]])
  out$package_version <- expected$package_version
  out$package_commit <- expected$package_commit
  out$cache_package_commit <- expected$cache_commit
  out$analysis_code_id <- expected$analysis_code_id
  out
}

.dependency_identical_rds <- function(existing, staged) {
  identical(digest::digest(existing, algo = "sha256", file = TRUE),
            digest::digest(staged, algo = "sha256", file = TRUE)) &&
    identical(readRDS(existing), readRDS(staged))
}

.dependency_split_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .dependency_split_args(args)
  bundle_path <- normalizePath(opt$bundle, mustWork = TRUE)
  tasks_path <- normalizePath(opt$tasks, mustWork = TRUE)
  output_parent <- normalizePath(dirname(opt$output_root), mustWork = TRUE)
  output_root <- file.path(output_parent, basename(opt$output_root))
  if (file.exists(output_root) && !dir.exists(output_root)) {
    stop("--output-root exists but is not a directory.", call. = FALSE)
  }
  dir.create(output_root, recursive = FALSE, showWarnings = FALSE)

  source_sha <- digest::digest(bundle_path, algo = "sha256", file = TRUE)
  task_sha <- digest::digest(tasks_path, algo = "sha256", file = TRUE)
  if (!identical(task_sha, opt$expected_task_sha256)) {
    stop("Frozen task table does not match --expected-task-sha256.",
         call. = FALSE)
  }
  bundle <- readRDS(bundle_path)
  tasks <- .dependency_read_tasks(tasks_path)
  expected <- list(
    method = opt$expected_method,
    package_version = opt$expected_package_version,
    package_commit = opt$expected_package_commit,
    cache_commit = opt$expected_cache_commit,
    analysis_code_id = opt$expected_analysis_code_id
  )
  validated <- .dependency_validate_bundle(bundle, tasks, expected)
  method_tasks <- validated$method_tasks
  cells <- validated$cells
  predictions <- validated$predictions
  data.table::setkey(cells, method, seed, cohort)
  if (nrow(predictions)) data.table::setkey(predictions, method, seed, cohort)

  if (!grepl("^[A-Za-z0-9_.-]+$", expected$method)) {
    stop("Expected method is unsafe as a manifest identifier.", call. = FALSE)
  }
  manifest_name <- paste0("output_manifest__", expected$method, ".tsv")
  manifest_path <- file.path(output_root, manifest_name)
  destinations <- file.path(output_root, method_tasks$label)
  if (!opt$allow_identical_overlap &&
      (file.exists(manifest_path) || any(file.exists(destinations)))) {
    stop("Output manifest or one or more task directories already exist; ",
         "use --allow-identical-overlap only for an exact restart.", call. = FALSE)
  }

  stage <- tempfile(".dependency-split-stage-", tmpdir = output_root)
  dir.create(stage)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  rows <- vector("list", nrow(method_tasks))
  for (i in seq_len(nrow(method_tasks))) {
    task <- method_tasks[i]
    cell <- cells[list(task$method, task$seed, task$cohort)]
    pred <- if (nrow(predictions)) {
      predictions[list(task$method, task$seed, task$cohort)]
    } else {
      predictions
    }
    if (nrow(cell) != 1L) stop("Internal shard lookup did not resolve one cell.")
    shard <- .dependency_shard_object(bundle, cell, pred, expected)
    shard_dir <- file.path(stage, task$label)
    dir.create(shard_dir)
    shard_path <- file.path(shard_dir, "within_bundle.rds")
    saveRDS(shard, shard_path, version = 3L)
    rows[[i]] <- data.table(
      label = task$label, method = task$method, seed = task$seed,
      cohort = task$cohort,
      path = file.path(task$label, "within_bundle.rds"),
      bytes = file.info(shard_path)$size,
      sha256 = digest::digest(shard_path, algo = "sha256", file = TRUE),
      source_bundle_sha256 = source_sha, task_table_sha256 = task_sha,
      package_version = expected$package_version,
      package_commit = expected$package_commit,
      cache_package_commit = expected$cache_commit,
      analysis_code_id = expected$analysis_code_id
    )
  }
  manifest <- data.table::rbindlist(rows)
  staged_manifest <- file.path(stage, manifest_name)
  data.table::fwrite(manifest, staged_manifest, sep = "\t", quote = FALSE,
                     na = "")

  # Fail before promotion if an explicit-overlap restart is not byte- and
  # object-identical. An existing manifest is accepted only for a complete,
  # already-finished output, never for a partial output tree.
  existing <- file.exists(destinations)
  for (i in which(existing)) {
    present <- list.files(destinations[[i]], all.files = TRUE, no.. = TRUE)
    existing_file <- file.path(destinations[[i]], "within_bundle.rds")
    staged_file <- file.path(stage, method_tasks$label[[i]], "within_bundle.rds")
    if (!opt$allow_identical_overlap ||
        !identical(present, "within_bundle.rds") ||
        !file.exists(existing_file) ||
        !.dependency_identical_rds(existing_file, staged_file)) {
      stop("Existing task output is not an exact identical restart: ",
           method_tasks$label[[i]], call. = FALSE)
    }
  }
  if (file.exists(manifest_path)) {
    if (!opt$allow_identical_overlap || !all(existing) ||
        !identical(digest::digest(manifest_path, algo = "sha256", file = TRUE),
                   digest::digest(staged_manifest, algo = "sha256", file = TRUE))) {
      stop("Existing output manifest is not an exact complete restart.",
           call. = FALSE)
    }
  }
  if (!identical(source_sha,
                 digest::digest(bundle_path, algo = "sha256", file = TRUE)) ||
      !identical(task_sha,
                 digest::digest(tasks_path, algo = "sha256", file = TRUE))) {
    stop("Source bundle or frozen task table changed during sharding.",
         call. = FALSE)
  }

  if (file.exists(manifest_path)) {
    cat("PASS: exact 170-shard output already complete and identical.\n")
    return(invisible(output_root))
  }
  for (i in which(!existing)) {
    from <- file.path(stage, method_tasks$label[[i]])
    if (!file.rename(from, destinations[[i]])) {
      stop("Could not atomically promote task directory: ",
           method_tasks$label[[i]], call. = FALSE)
    }
  }
  # The manifest is the completion marker and is deliberately promoted last.
  if (!file.rename(staged_manifest, manifest_path)) {
    stop("Could not atomically promote the output manifest.", call. = FALSE)
  }
  cat(sprintf("PASS: wrote %d validated one-cell task bundles for %s.\n",
              nrow(method_tasks), expected$method))
  invisible(output_root)
}

if (sys.nframe() == 0L) .dependency_split_main()
