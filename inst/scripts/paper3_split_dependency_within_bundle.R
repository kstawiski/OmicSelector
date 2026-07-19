#!/usr/bin/env Rscript

# Validate and shard one frozen 170-cell optional-dependency benchmark bundle.
# This is a serialization-only utility: it never fits, scores, or otherwise
# recomputes a scientific result.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.DEPENDENCY_SPLITTER_PATH <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    candidate <- sub("^--file=", "", file_arg[[1L]])
    if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
  }
  candidate <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(candidate) && file.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }
  root <- normalizePath(getwd(), mustWork = TRUE)
  for (i in 0:6) {
    candidate <- file.path(root, "inst", "scripts",
                           "paper3_split_dependency_within_bundle.R")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(root)
    if (identical(parent, root)) break
    root <- parent
  }
  candidate <- tryCatch(
    system.file("scripts", "paper3_split_dependency_within_bundle.R",
                package = "OmicSelector"),
    error = function(e) ""
  )
  if (nzchar(candidate) && file.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }
  NA_character_
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
  witness_flags <- c(
    structural_witness = "--structural-witness",
    expected_structural_witness_sha256 =
      "--expected-structural-witness-sha256",
    expected_structural_cache_sha256 =
      "--expected-structural-cache-sha256",
    expected_structural_fold_engine_sha256 =
      "--expected-structural-fold-engine-sha256"
  )
  witness_values <- lapply(witness_flags, value)
  present <- !vapply(witness_values, function(x) is.na(x) || !nzchar(x),
                     logical(1L))
  if (any(present) && !all(present)) {
    stop("Structural witness arguments must be supplied together: ",
         paste(unname(witness_flags), collapse = ", "), call. = FALSE)
  }
  out <- c(out, witness_values)
  if (!grepl("^[0-9a-f]{40}$", out$expected_package_commit) ||
      !grepl("^[0-9a-f]{40}$", out$expected_cache_commit) ||
      !grepl("^[0-9a-f]{16}$", out$expected_analysis_code_id) ||
      !grepl("^[0-9a-f]{64}$", out$expected_task_sha256) ||
      (all(present) &&
       (!grepl("^[0-9a-f]{64}$",
              out$expected_structural_witness_sha256) ||
        !grepl("^[0-9a-f]{64}$",
               out$expected_structural_cache_sha256) ||
        !grepl("^[0-9a-f]{64}$",
               out$expected_structural_fold_engine_sha256)))) {
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

.dependency_validate_ineligible_cells <- function(cells,
                                                   witnessed_keys = NULL) {
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
  if (!is.null(witnessed_keys) && nrow(witnessed_keys)) {
    witness_key <- do.call(paste, witnessed_keys[, .(method, cohort, seed)])
    cell_key <- do.call(paste, ineligible[, .(method, cohort, seed)])
    reviewed_structural <- reviewed_structural | cell_key %in% witness_key
  }

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

.dependency_validate_structural_witness <- function(
    witness, cells, expected, source_bundle_sha256, task_table_sha256,
    expected_cache_sha256, expected_fold_engine_sha256, witness_sha256) {
  required <- c(
    "schema_version", "producer", "producer_script_sha256",
    "source_bundle_sha256", "task_table_sha256", "cache_sha256",
    "fold_engine_sha256",
    "package_version", "package_commit", "cache_package_commit",
    "analysis_code_id", "analysis_sources", "replay_sources", "cells",
    "predictions", "fold_audit", "runtime"
  )
  if (!is.list(witness) || any(!required %in% names(witness))) {
    stop("Structural witness lacks one or more required fields.",
         call. = FALSE)
  }
  scalar <- function(x, value, label) {
    if (!is.character(x) || length(x) != 1L || is.na(x) ||
        !identical(x, value)) {
      stop("Structural witness ", label, " does not match its required pin.",
           call. = FALSE)
    }
  }
  scalar(witness$schema_version,
         "OmicSelector-dependency-structural-witness-v1", "schema_version")
  scalar(witness$producer, "paper3_build_dependency_structural_witness.R",
         "producer")
  if (!is.character(witness$producer_script_sha256) ||
      length(witness$producer_script_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", witness$producer_script_sha256)) {
    stop("Structural witness lacks an exact producer-script SHA-256.",
         call. = FALSE)
  }
  if (is.na(.DEPENDENCY_SPLITTER_PATH)) {
    stop("Could not resolve the package-owned splitter/producer directory.",
         call. = FALSE)
  }
  producer_path <- file.path(dirname(.DEPENDENCY_SPLITTER_PATH),
                             "paper3_build_dependency_structural_witness.R")
  if (!file.exists(producer_path) ||
      !identical(witness$producer_script_sha256,
                 digest::digest(producer_path, algo = "sha256", file = TRUE))) {
    stop("Structural witness producer_script_sha256 does not match the ",
         "package-owned producer bytes.", call. = FALSE)
  }
  scalar(witness$source_bundle_sha256, source_bundle_sha256,
         "source_bundle_sha256")
  scalar(witness$task_table_sha256, task_table_sha256,
         "task_table_sha256")
  scalar(witness$cache_sha256, expected_cache_sha256, "cache_sha256")
  scalar(witness$fold_engine_sha256, expected_fold_engine_sha256,
         "fold_engine_sha256")
  scalar(witness$package_version, expected$package_version, "package_version")
  scalar(witness$package_commit, expected$package_commit, "package_commit")
  scalar(witness$cache_package_commit, expected$cache_commit,
         "cache_package_commit")
  scalar(witness$analysis_code_id, expected$analysis_code_id,
         "analysis_code_id")
  if (!is.character(witness_sha256) || length(witness_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", witness_sha256)) {
    stop("Structural witness bytes are not SHA-256 pinned.", call. = FALSE)
  }

  sources <- data.table::as.data.table(witness$analysis_sources)
  expected_roles <- c("benchmark_runner", "benchmark_engine",
                      "paired_auc_engine", "design_contract")
  expected_paths <- c(
    "code/analyses/run_singlesample_public_benchmark_20260715.R",
    "code/methods/singlesample_roster_benchmark.R",
    "code/methods/paired_auc_diff_se.R",
    "code/methods/singlesample_public_design.R"
  )
  if (!identical(names(sources), c("role", "relative_path", "sha256")) ||
      nrow(sources) != 4L || !identical(sources$role, expected_roles) ||
      !identical(sources$relative_path, expected_paths) ||
      any(!grepl("^[0-9a-f]{64}$", sources$sha256))) {
    stop("Structural witness analysis source closure is malformed.",
         call. = FALSE)
  }
  derived_code_id <- substr(digest::digest(
    paste(sources$sha256, collapse = "|"), algo = "sha256"
  ), 1L, 16L)
  if (!identical(derived_code_id, expected$analysis_code_id)) {
    stop("Structural witness source hashes do not reproduce analysis_code_id.",
         call. = FALSE)
  }
  replay_sources <- data.table::as.data.table(witness$replay_sources)
  if (!all(c("role", "relative_path", "sha256") %in% names(replay_sources)) ||
      nrow(replay_sources) != 5L ||
      !identical(replay_sources[seq_len(4L), .(role, relative_path, sha256)],
                 sources[, .(role, relative_path, sha256)]) ||
      sum(replay_sources$role == "fold_engine" &
            replay_sources$relative_path ==
              "code/methods/matched_null_benchmark.R") != 1L ||
      any(!grepl("^[0-9a-f]{64}$", replay_sources$sha256))) {
    stop("Structural witness replay source closure is incomplete.",
         call. = FALSE)
  }
  fold_row <- replay_sources[role == "fold_engine" & relative_path ==
                               "code/methods/matched_null_benchmark.R"]
  if (nrow(fold_row) != 1L ||
      !identical(fold_row$sha256[[1L]], expected_fold_engine_sha256) ||
      !identical(fold_row$sha256[[1L]], witness$fold_engine_sha256)) {
    stop("Structural witness fold-engine SHA does not match its external pin.",
         call. = FALSE)
  }

  wc <- data.table::as.data.table(witness$cells)
  wp <- data.table::as.data.table(witness$predictions)
  wa <- data.table::as.data.table(witness$fold_audit)
  required_cells <- c(
    "method", "cohort", "seed", "ineligible_reason", "n_valid_folds",
    "outer_k", "n_group_split_violations", "n_profiles"
  )
  required_predictions <- c(
    "method", "cohort", "seed", "fold", "n_train", "n_test",
    "sample_idx", "sample_id", "group_id", "y", "score_m", "score_b",
    "package_commit", "cache_package_commit", "analysis_code_id"
  )
  required_audit <- c(
    "method", "cohort", "seed", "fold", "n_train", "n_test",
    "n_pos_test", "n_neg_test", "n_unique_method", "n_unique_baseline",
    "n_unique_outcome", "n_train_finite", "n_test_finite",
    "training_direction", "n_group_split_violations", "fit_status",
    "baseline_status", "training_score_status", "heldout_score_status",
    "status", "structural_basis", "package_commit",
    "cache_package_commit", "analysis_code_id"
  )
  if (!nrow(wc) || any(!required_cells %in% names(wc)) ||
      !nrow(wp) || any(!required_predictions %in% names(wp)) ||
      !nrow(wa) || any(!required_audit %in% names(wa))) {
    stop("Structural witness cell/fold/prediction surfaces are incomplete.",
         call. = FALSE)
  }
  key <- c("method", "cohort", "seed")
  if (anyDuplicated(wc[, ..key]) || anyNA(wc[, ..required_cells]) ||
      any(wc$cohort == "GSE83977") ||
      any(wc$method != expected$method) ||
      any(!wc$seed %in% c(101L, 202L, 303L, 404L, 505L)) ||
      any(wc$n_valid_folds < 0L | wc$n_valid_folds >= wc$outer_k) ||
      any(wc$n_group_split_violations != 0L) ||
      any(wc$n_profiles < 4L)) {
    stop("Structural witness cell set is invalid or outside its authority.",
         call. = FALSE)
  }
  source_cells <- merge(
    wc, cells[, .(method, cohort, seed, eligible, source_reason = ineligible_reason,
                  source_valid = n_valid_folds, source_outer_k = outer_k,
                  source_group_violations = n_group_split_violations,
                  source_n_eff = n_eff)],
    by = key, all.x = TRUE, sort = FALSE
  )
  if (nrow(source_cells) != nrow(wc) || anyNA(source_cells$eligible) ||
      any(source_cells$eligible) ||
      any(!is.na(source_cells$source_n_eff)) ||
      any(source_cells$ineligible_reason != source_cells$source_reason) ||
      any(source_cells$n_valid_folds != source_cells$source_valid) ||
      any(source_cells$outer_k != source_cells$source_outer_k) ||
      any(source_cells$source_group_violations != 0L)) {
    stop("Structural witness cells do not match source bundle cells/reasons.",
         call. = FALSE)
  }

  for (surface in list(wp, wa)) {
    .dependency_validate_pin_column(surface, "package_commit",
                                    expected$package_commit,
                                    "Structural witness")
    .dependency_validate_pin_column(surface, "cache_package_commit",
                                    expected$cache_commit,
                                    "Structural witness")
    .dependency_validate_pin_column(surface, "analysis_code_id",
                                    expected$analysis_code_id,
                                    "Structural witness")
  }
  if (anyNA(wp[, .(method, cohort, seed, fold, n_train, n_test, sample_idx,
                    sample_id, group_id, y, score_m, score_b)]) ||
      any(!nzchar(wp$sample_id)) || any(!nzchar(wp$group_id)) ||
      any(!is.finite(wp$score_m)) || any(!is.finite(wp$score_b)) ||
      any(!wp$y %in% c(0L, 1L)) || any(wp$fold < 1L) ||
      any(wp$n_train < 1L) || any(wp$n_test < 1L) ||
      anyDuplicated(wp[, .(method, cohort, seed, sample_idx)])) {
    stop("Structural witness predictions are missing, non-finite, or invalid.",
         call. = FALSE)
  }
  witness_keys <- unique(wp[, ..key])
  if (!setequal(do.call(paste, witness_keys), do.call(paste, wc[, ..key])) ||
      !setequal(do.call(paste, unique(wa[, ..key])),
                do.call(paste, wc[, ..key]))) {
    stop("Structural witness predictions/audits do not cover its exact cell set.",
         call. = FALSE)
  }
  group_folds <- wp[, .(n_folds = uniqueN(fold), n_labels = uniqueN(y)),
                    by = .(method, cohort, seed, group_id)]
  if (any(group_folds$n_folds != 1L) || any(group_folds$n_labels != 1L)) {
    stop("Structural witness splits a provenance group across folds/outcomes.",
         call. = FALSE)
  }
  if (anyDuplicated(wa[, .(method, cohort, seed, fold)]) ||
      anyNA(wa[, setdiff(required_audit, "structural_basis"), with = FALSE]) ||
      any(wa$fit_status != "ok") || any(wa$baseline_status != "ok") ||
      any(wa$training_score_status != "ok") ||
      any(wa$heldout_score_status != "ok") ||
      any(!wa$status %in% c("valid", "structural_constant")) ||
      any(!wa$training_direction %in% c(-1, 1)) ||
      any(wa$n_group_split_violations != 0L)) {
    stop("Structural witness audit contains an error, missing fold, or ",
         "unrecognized status.", call. = FALSE)
  }

  for (i in seq_len(nrow(wc))) {
    cell <- wc[i]
    pred <- wp[method == cell$method & cohort == cell$cohort &
                 seed == cell$seed]
    audit <- wa[method == cell$method & cohort == cell$cohort &
                  seed == cell$seed]
    data.table::setorder(audit, fold)
    if (!identical(sort(unique(as.integer(pred$fold))),
                   seq_len(as.integer(cell$outer_k))) ||
        !identical(as.integer(audit$fold),
                   seq_len(as.integer(cell$outer_k))) ||
        nrow(pred) != cell$n_profiles ||
        !setequal(as.integer(pred$sample_idx), seq_len(cell$n_profiles))) {
      stop("Structural witness lacks exact profile/fold coverage.",
           call. = FALSE)
    }
    recomputed_status <- character(cell$outer_k)
    recomputed_basis <- character(cell$outer_k)
    for (fold in seq_len(cell$outer_k)) {
      fold_value <- as.integer(fold)
      d <- pred[fold == fold_value]
      a <- audit[fold == fold_value]
      n_um <- uniqueN(d$score_m)
      n_ub <- uniqueN(d$score_b)
      n_uy <- uniqueN(d$y)
      basis <- paste(c(if (n_um < 2L) "method_constant" else character(),
                       if (n_ub < 2L) "baseline_constant" else character()),
                     collapse = "+")
      status <- if (nzchar(basis)) "structural_constant" else "valid"
      if (!nzchar(basis)) basis <- NA_character_
      recomputed_status[[fold]] <- status
      recomputed_basis[[fold]] <- basis
      if (nrow(a) != 1L || nrow(d) != unique(d$n_test) ||
          length(unique(d$n_test)) != 1L ||
          length(unique(d$n_train)) != 1L ||
          unique(d$n_train) + unique(d$n_test) != cell$n_profiles ||
          a$n_train != unique(d$n_train) || a$n_test != nrow(d) ||
          a$n_pos_test != sum(d$y == 1L) ||
          a$n_neg_test != sum(d$y == 0L) ||
          a$n_unique_method != n_um || a$n_unique_baseline != n_ub ||
          a$n_unique_outcome != n_uy || n_uy != 2L || nrow(d) < 4L ||
          a$n_test_finite != nrow(d) || a$n_train_finite != a$n_train ||
          a$status != status || !identical(a$structural_basis[[1L]], basis)) {
        stop("Structural witness fold audit does not derive from its finite ",
             "predictions.", call. = FALSE)
      }
    }
    n_valid <- sum(recomputed_status == "valid")
    reason <- sprintf(
      "n_valid_folds=%d < %d; per-fold: %s", n_valid, cell$outer_k,
      paste(ifelse(recomputed_status == "valid", "ok",
                   "method/baseline constant or all-NA on test fold"),
            collapse = "; ")
    )
    if (n_valid != cell$n_valid_folds ||
        !identical(reason, cell$ineligible_reason)) {
      stop("Structural witness does not reproduce source fold eligibility/reason.",
           call. = FALSE)
    }
  }
  invisible(list(keys = wc[, ..key], witness = witness,
                 witness_sha256 = witness_sha256))
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

.dependency_validate_bundle <- function(bundle, tasks, expected,
                                        structural_witness = NULL,
                                        source_bundle_sha256 = NULL,
                                        task_table_sha256 = NULL,
                                        expected_cache_sha256 = NULL,
                                        expected_fold_engine_sha256 = NULL,
                                        witness_sha256 = NULL) {
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
  witness_validation <- NULL
  if (!is.null(structural_witness)) {
    if (is.null(source_bundle_sha256) || is.null(task_table_sha256) ||
        is.null(expected_cache_sha256) ||
        is.null(expected_fold_engine_sha256) || is.null(witness_sha256)) {
      stop("Structural witness validation lacks required external SHA pins.",
           call. = FALSE)
    }
    witness_validation <- .dependency_validate_structural_witness(
      structural_witness, cells, expected, source_bundle_sha256,
      task_table_sha256, expected_cache_sha256,
      expected_fold_engine_sha256, witness_sha256
    )
  }
  witnessed_keys <- if (is.null(witness_validation)) NULL else
    witness_validation$keys
  .dependency_validate_ineligible_cells(cells, witnessed_keys)
  .dependency_validate_predictions(cells, bundle$predictions, expected)
  invisible(list(cells = cells, predictions = data.table::as.data.table(bundle$predictions),
                 method_tasks = method_tasks,
                 structural_witness = witness_validation))
}

.dependency_shard_object <- function(bundle, cell, predictions, expected,
                                     structural_witness = NULL) {
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
  if (!is.null(structural_witness)) {
    key_match <- function(x) {
      x$method == cell$method[[1L]] & x$cohort == cell$cohort[[1L]] &
        x$seed == cell$seed[[1L]]
    }
    witness <- structural_witness$witness
    wc <- data.table::as.data.table(witness$cells)
    if (any(key_match(wc))) {
      wp <- data.table::as.data.table(witness$predictions)
      wa <- data.table::as.data.table(witness$fold_audit)
      out$structural_ineligibility_witness <- list(
        schema_version = witness$schema_version,
        witness_sha256 = structural_witness$witness_sha256,
        producer = witness$producer,
        producer_script_sha256 = witness$producer_script_sha256,
        source_bundle_sha256 = witness$source_bundle_sha256,
        task_table_sha256 = witness$task_table_sha256,
        cache_sha256 = witness$cache_sha256,
        cells = data.table::copy(wc[key_match(wc)]),
        predictions = data.table::copy(wp[key_match(wp)]),
        fold_audit = data.table::copy(wa[key_match(wa)])
      )
    }
  }
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
  structural_witness <- NULL
  witness_sha <- NULL
  if (!is.na(opt$structural_witness) && nzchar(opt$structural_witness)) {
    witness_path <- normalizePath(opt$structural_witness, mustWork = TRUE)
    witness_sha <- digest::digest(witness_path, algo = "sha256", file = TRUE)
    if (!identical(witness_sha,
                   opt$expected_structural_witness_sha256)) {
      stop("Structural witness does not match its expected SHA-256 pin.",
           call. = FALSE)
    }
    structural_witness <- readRDS(witness_path)
  }
  tasks <- .dependency_read_tasks(tasks_path)
  expected <- list(
    method = opt$expected_method,
    package_version = opt$expected_package_version,
    package_commit = opt$expected_package_commit,
    cache_commit = opt$expected_cache_commit,
    analysis_code_id = opt$expected_analysis_code_id
  )
  validated <- .dependency_validate_bundle(
    bundle, tasks, expected,
    structural_witness = structural_witness,
    source_bundle_sha256 = source_sha,
    task_table_sha256 = task_sha,
    expected_cache_sha256 = opt$expected_structural_cache_sha256,
    expected_fold_engine_sha256 =
      opt$expected_structural_fold_engine_sha256,
    witness_sha256 = witness_sha
  )
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
    shard <- .dependency_shard_object(
      bundle, cell, pred, expected, validated$structural_witness
    )
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
    if (!is.null(validated$structural_witness)) {
      rows[[i]][, structural_witness_sha256 :=
        if ("structural_ineligibility_witness" %in% names(shard))
          validated$structural_witness$witness_sha256 else NA_character_]
    }
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
                 digest::digest(tasks_path, algo = "sha256", file = TRUE)) ||
      (!is.null(witness_sha) &&
       !identical(witness_sha,
                  digest::digest(witness_path, algo = "sha256", file = TRUE)))) {
    stop("Source bundle, frozen task table, or structural witness changed ",
         "during sharding.",
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
