#!/usr/bin/env Rscript

# Fail-closed assembler and inferential summarizer for seven independently
# executed external rank-tree sensitivity method bundles. The 21 pair
# directions are frozen by registry order; no outcome-dependent orientation,
# method selection, or pair deletion is permitted.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.rta_method_ids <- function() {
  c("ws-rclr-trimmed", "reo-ktsp", "reo-pairratio", "reo-singscore",
    "reo-ucell", "reo-rankforest", "reo-rankboost")
}

.rta_seeds <- function() c(101L, 202L, 303L, 404L, 505L)

.rta_arg_value <- function(args, flag) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) NA_character_ else args[[hit + 1L]]
}

.rta_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  flags <- c(
    package_root = "--package-root", paper_root = "--paper-root",
    analysis_plan = "--analysis-plan",
    analysis_amendment = "--analysis-amendment", bundles = "--bundles",
    output_dir = "--output-dir",
    runtime_image = "--runtime-image",
    runtime_attestation_manifest = "--runtime-attestation-manifest",
    expected_package_version = "--expected-package-version",
    expected_analysis_plan_sha256 = "--expected-analysis-plan-sha256",
    expected_analysis_amendment_sha256 =
      "--expected-analysis-amendment-sha256",
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
    expected_runner_script_sha256 = "--expected-runner-script-sha256",
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
  out <- lapply(flags, function(flag) .rta_arg_value(args, flag))
  missing <- names(out)[vapply(out, function(x) is.na(x) || !nzchar(x), logical(1L))]
  if (length(missing)) {
    stop("Missing required argument(s): ",
         paste(unname(flags[missing]), collapse = ", "), call. = FALSE)
  }
  hex <- c(
    expected_analysis_plan_sha256 = 64L,
    expected_analysis_amendment_sha256 = 64L,
    expected_package_commit = 40L, expected_cache_commit = 40L,
    expected_cache_sha256 = 64L, expected_cache_manifest_sha256 = 64L,
    expected_cache_digest_sha256 = 64L,
    expected_runtime_image_sha256 = 64L,
    expected_runtime_attestation_manifest_sha256 = 64L,
    expected_container_launcher_sha256 = 64L,
    expected_runner_script_sha256 = 64L,
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
  out$bundles <- strsplit(out$bundles, ",", fixed = TRUE)[[1L]]
  if (length(out$bundles) != 7L || any(!nzchar(out$bundles)) ||
      anyDuplicated(out$bundles)) {
    stop("--bundles must name exactly seven unique method-bundle RDS files.",
         call. = FALSE)
  }
  out
}

.rta_sha <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

.rta_assert_file_pin <- function(path, expected, label) {
  if (!file.exists(path) || !identical(.rta_sha(path), expected)) {
    stop(label, " differs from its exact expected pin.", call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

.rta_git_provenance <- function(package_root) {
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

.rta_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  }
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
}

.rta_structural_reason <- function() {
  paste0("n_valid_folds=0 < 3; per-fold: ", paste(rep(
    "method/baseline constant or all-NA on test fold", 3L), collapse = "; "))
}

.rta_assert_pin <- function(observed, expected, label) {
  if (!is.character(observed) || length(observed) != 1L || is.na(observed) ||
      !identical(observed, expected)) {
    stop(label, " does not match its exact expected pin.", call. = FALSE)
  }
  invisible(TRUE)
}

.rta_validate_output_manifest <- function(bundle_path) {
  bundle_path <- normalizePath(bundle_path, mustWork = TRUE)
  root <- dirname(bundle_path)
  manifest_path <- file.path(root, "output_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop("Method bundle lacks its completion manifest: ", bundle_path,
         call. = FALSE)
  }
  manifest <- data.table::fread(manifest_path)
  expected_files <- c(
    "analysis_plan.md", "analysis_amendment.md", "engine_manifest.tsv",
    "environment_manifest.tsv",
    "methods.tsv",
    "provenance_preflight.log", "provenance_resolution.tsv", "report.md",
    "runtime_attestation.tsv",
    "units.tsv", "within_bundle.rds", "within_cells.tsv",
    "within_predictions.tsv.gz"
  )
  if (!all(c("file", "bytes", "sha256") %in% names(manifest)) ||
      anyNA(manifest[, .(file, bytes, sha256)]) || anyDuplicated(manifest$file) ||
      !setequal(manifest$file, expected_files) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Method-bundle completion manifest is malformed.", call. = FALSE)
  }
  paths <- file.path(root, manifest$file)
  if (any(!file.exists(paths)) ||
      any(file.info(paths)$size != manifest$bytes) ||
      !identical(unname(vapply(paths, .rta_sha, character(1L))),
                 as.character(manifest$sha256))) {
    stop("Method-bundle artifact bytes differ from their completion manifest.",
         call. = FALSE)
  }
  list(path = bundle_path, root = root, manifest = manifest,
       manifest_sha256 = .rta_sha(manifest_path),
       bundle_sha256 = .rta_sha(bundle_path))
}

.rta_validate_cell <- function(cell) {
  required <- c("eligible", "cohort", "seed", "ineligible_reason", "n_eff",
                "n_valid_folds", "outer_k", "n_group_split_violations",
                "structural_exception_proven",
                "structural_exception_proof_sha256")
  if (!is.data.frame(cell) || nrow(cell) != 1L ||
      !all(required %in% names(cell)) || !is.logical(cell$eligible) ||
      is.na(cell$eligible)) {
    stop("Sensitivity cell has an invalid structural-proof schema.",
         call. = FALSE)
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
      stop("Eligible sensitivity cell violates fold/group accounting.",
           call. = FALSE)
    }
    return(invisible(TRUE))
  }
  accepted <- identical(as.character(cell$cohort), "GSE83977") &&
    as.integer(cell$seed) %in% .rta_seeds() &&
    identical(as.integer(cell$n_valid_folds), 0L) &&
    identical(as.integer(cell$outer_k), 3L) &&
    identical(as.integer(cell$n_group_split_violations), 0L) &&
    identical(as.character(cell$ineligible_reason), .rta_structural_reason()) &&
    identical(cell$structural_exception_proven, TRUE) &&
    is.character(cell$structural_exception_proof_sha256) &&
    length(cell$structural_exception_proof_sha256) == 1L &&
    grepl("^[0-9a-f]{64}$", cell$structural_exception_proof_sha256)
  if (!accepted) {
    stop("Method bundle contains an unreviewed ineligibility or execution ",
         "failure: ", as.character(cell$ineligible_reason), call. = FALSE)
  }
  invisible(TRUE)
}

.rta_validate_bundle <- function(bundle, expected_method, expected, units = NULL,
                                 seeds = .rta_seeds()) {
  required <- c(
    "cells", "predictions", "units", "methods", "method_registry", "seeds",
    "cohorts", "package_version", "package_commit", "cache_package_commit",
    "cache_sha256", "cache_manifest_sha256", "cache_digest_sha256",
    "runner_script_sha256", "engine_manifest_sha256",
    "environment_manifest_sha256", "provenance_script_sha256",
    "runtime_image_sha256", "runtime_attestation_manifest_sha256",
    "container_launcher_sha256", "provenance_manifest_sha256",
    "provenance_union_inventory_sha256", "analysis_code_id", "r_version",
    "analysis_plan_sha256", "analysis_amendment_sha256",
    "ranktree_dependency_closure"
  )
  if (!is.list(bundle) || length(setdiff(required, names(bundle)))) {
    stop("Method bundle lacks required contract fields.", call. = FALSE)
  }
  scalar_pins <- c(
    package_version = "package_version", package_commit = "package_commit",
    cache_commit = "cache_package_commit", cache_sha256 = "cache_sha256",
    cache_manifest_sha256 = "cache_manifest_sha256",
    cache_digest_sha256 = "cache_digest_sha256",
    runner_script_sha256 = "runner_script_sha256",
    engine_manifest_sha256 = "engine_manifest_sha256",
    environment_manifest_sha256 = "environment_manifest_sha256",
    runtime_image_sha256 = "runtime_image_sha256",
    runtime_attestation_manifest_sha256 =
      "runtime_attestation_manifest_sha256",
    container_launcher_sha256 = "container_launcher_sha256",
    r_version = "r_version",
    provenance_script_sha256 = "provenance_script_sha256",
    provenance_manifest_sha256 = "provenance_manifest_sha256",
    provenance_union_inventory_sha256 =
      "provenance_union_inventory_sha256",
    analysis_plan_sha256 = "analysis_plan_sha256",
    analysis_amendment_sha256 = "analysis_amendment_sha256",
    analysis_code_id = "analysis_code_id"
  )
  for (name in names(scalar_pins)) {
    .rta_assert_pin(bundle[[scalar_pins[[name]]]], expected[[name]],
                    paste("Bundle", scalar_pins[[name]]))
  }
  if (!identical(as.character(bundle$methods), expected_method) ||
      !identical(as.integer(bundle$seeds), as.integer(seeds))) {
    stop("Each input must be one exact method-by-five-seed bundle.",
         call. = FALSE)
  }
  method_registry <- data.table::as.data.table(bundle$method_registry)
  if (nrow(method_registry) != 1L ||
      !identical(as.character(method_registry$method_id), expected_method) ||
      !identical(as.character(method_registry$estimand), "within")) {
    stop("Method bundle registry row does not match its expected method.",
         call. = FALSE)
  }
  observed_units <- data.table::as.data.table(bundle$units)
  required_units <- c("unit_id", "disease", "modality", "biospecimen",
                      "provenance_block", "source_accessions", "outer_k",
                      "primary_unit", "n_deposited_profiles", "n_unique_groups")
  if (!all(required_units %in% names(observed_units)) ||
      nrow(observed_units) < 1L || anyDuplicated(observed_units$unit_id) ||
      !is.logical(observed_units$primary_unit) ||
      anyNA(observed_units[, ..required_units]) ||
      any(observed_units$n_deposited_profiles < observed_units$n_unique_groups) ||
      any(observed_units$n_unique_groups < 1L)) {
    stop("Bundle unit metadata/count accounting is invalid.", call. = FALSE)
  }
  if (!is.null(units) && !identical(observed_units, units)) {
    stop("Method bundles do not carry byte-identical unit metadata/counts.",
         call. = FALSE)
  }
  if (!identical(as.character(bundle$cohorts), observed_units$unit_id)) {
    stop("Bundle cohort vector differs from its unit table.", call. = FALSE)
  }
  cells <- data.table::as.data.table(bundle$cells)
  required_cells <- c(
    "method", "cohort", "seed", "eligible", "ineligible_reason", "n_eff",
    "n_valid_folds", "outer_k", "grouping_policy",
    "n_group_split_violations", "structural_exception_proven",
    "structural_exception_basis", "structural_exception_proof_sha256",
    "package_version", "package_commit",
    "cache_package_commit", "cache_sha256", "engine_manifest_sha256",
    "environment_manifest_sha256", "runtime_image_sha256", "analysis_code_id",
    "runtime_attestation_manifest_sha256", "container_launcher_sha256",
    "r_version", "provenance_manifest_sha256",
    "provenance_union_inventory_sha256", "analysis_plan_sha256",
    "analysis_amendment_sha256"
  )
  if (!all(required_cells %in% names(cells)) ||
      nrow(cells) != length(seeds) * nrow(observed_units) ||
      anyDuplicated(cells[, .(method, cohort, seed)]) ||
      !identical(unique(cells$method), expected_method) ||
      !setequal(cells$cohort, observed_units$unit_id) ||
      !setequal(cells$seed, seeds) || !is.logical(cells$eligible) ||
      anyNA(cells$eligible)) {
    stop("Bundle cells are not its complete method-by-unit-by-seed grid.",
         call. = FALSE)
  }
  for (pin in c("package_version", "package_commit", "cache_package_commit",
                "cache_sha256", "engine_manifest_sha256",
                "environment_manifest_sha256", "runtime_image_sha256",
                "runtime_attestation_manifest_sha256",
                "container_launcher_sha256", "r_version",
                "provenance_manifest_sha256",
                "provenance_union_inventory_sha256", "analysis_plan_sha256",
                "analysis_amendment_sha256",
                "analysis_code_id")) {
    expected_name <- switch(
      pin, package_version = "package_version", package_commit = "package_commit",
      cache_package_commit = "cache_commit", pin)
    if (anyNA(cells[[pin]]) || any(cells[[pin]] != expected[[expected_name]])) {
      stop("Cell table differs from bundle pin: ", pin, call. = FALSE)
    }
  }
  invisible(lapply(seq_len(nrow(cells)), function(i) .rta_validate_cell(cells[i])))
  predictions <- data.table::as.data.table(bundle$predictions)
  required_predictions <- c(
    "method", "cohort", "seed", "fold", "split_id", "n_train", "n_test",
    "sample_id", "group_id", "y", "score_m", "package_version",
    "package_commit", "cache_package_commit", "runtime_image_sha256",
    "runtime_attestation_manifest_sha256", "container_launcher_sha256",
    "r_version", "provenance_manifest_sha256",
    "provenance_union_inventory_sha256", "analysis_plan_sha256",
    "analysis_amendment_sha256",
    "analysis_code_id"
  )
  if (!all(required_predictions %in% names(predictions)) ||
      (!nrow(predictions) && any(cells$eligible)) ||
      anyNA(predictions[, ..required_predictions]) ||
      any(!is.finite(predictions$score_m)) ||
      any(!predictions$y %in% c(0L, 1L)) ||
      any(predictions$n_train < 1L | predictions$n_test < 1L) ||
      anyDuplicated(predictions[, .(method, cohort, seed, sample_id)])) {
    stop("Bundle held-out predictions are incomplete, non-finite, or ",
         "duplicated.", call. = FALSE)
  }
  for (pin in c("package_version", "package_commit", "cache_package_commit",
                "runtime_image_sha256",
                "runtime_attestation_manifest_sha256",
                "container_launcher_sha256", "r_version",
                "provenance_manifest_sha256",
                "provenance_union_inventory_sha256", "analysis_plan_sha256",
                "analysis_amendment_sha256",
                "analysis_code_id")) {
    expected_name <- switch(
      pin, package_version = "package_version", package_commit = "package_commit",
      cache_package_commit = "cache_commit", pin)
    if (any(predictions[[pin]] != expected[[expected_name]])) {
      stop("Prediction table differs from bundle pin: ", pin, call. = FALSE)
    }
  }
  eligible <- cells[eligible == TRUE, .(method, cohort, seed, n_eff, outer_k)]
  pred_cells <- unique(predictions[, .(method, cohort, seed)])
  if (!setequal(do.call(paste, eligible[, .(method, cohort, seed)]),
                do.call(paste, pred_cells))) {
    stop("Prediction cells do not equal eligible cells.", call. = FALSE)
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
    stop("Prediction rows do not satisfy eligible n_eff/outer-K accounting.",
         call. = FALSE)
  }
  group_audit <- predictions[, .(n_fold = data.table::uniqueN(fold),
                                 n_y = data.table::uniqueN(y)),
                             by = .(method, cohort, seed, group_id)]
  if (any(group_audit$n_fold != 1L) || any(group_audit$n_y != 1L)) {
    stop("A provenance group crosses folds or outcome labels.", call. = FALSE)
  }
  fold_accounting <- predictions[, .(
    observed_n_test = .N, declared_n_test = data.table::uniqueN(n_test),
    declared_n_train = data.table::uniqueN(n_train),
    n_test_value = unique(n_test), n_train_value = unique(n_train)),
    by = .(method, cohort, seed, fold)]
  fold_accounting <- merge(
    fold_accounting,
    observed_units[, .(cohort = unit_id, n_deposited_profiles)],
    by = "cohort", all.x = TRUE, sort = FALSE)
  if (any(fold_accounting$declared_n_test != 1L) ||
      any(fold_accounting$declared_n_train != 1L) ||
      any(fold_accounting$observed_n_test != fold_accounting$n_test_value) ||
      any(fold_accounting$n_test_value + fold_accounting$n_train_value !=
            fold_accounting$n_deposited_profiles)) {
    stop("Fold train/test declarations do not match held-out profile counts.",
         call. = FALSE)
  }
  split_check <- predictions[, {
    reference <- .SD[, .(
      fold = as.integer(fold), sample_id = as.character(sample_id),
      group_id = as.character(group_id), y = as.integer(y),
      n_train = as.integer(n_train), n_test = as.integer(n_test))]
    data.table::setorder(reference, fold, sample_id, group_id, y)
    .(declared_split_n = data.table::uniqueN(split_id),
      declared_split_id = unique(split_id),
      recomputed_split_id = digest::digest(
        reference, algo = "sha256", serialize = TRUE))
  }, by = .(method, cohort, seed)]
  if (any(split_check$declared_split_n != 1L) ||
      any(split_check$declared_split_id != split_check$recomputed_split_id)) {
    stop("A declared split id does not reproduce from exact held-out rows.",
         call. = FALSE)
  }
  invisible(list(cells = cells, predictions = predictions,
                 units = observed_units, method_registry = method_registry))
}

.rta_common_split_audit <- function(predictions, cells, methods) {
  rows <- list()
  index <- 0L
  for (unit_value in unique(cells$cohort)) {
    for (seed_value in .rta_seeds()) {
      eligible <- cells[cohort == unit_value & seed == seed_value &
                          eligible == TRUE,
                        method]
      if (!length(eligible)) next
      reference <- predictions[method == eligible[[1L]] &
                                 cohort == unit_value & seed == seed_value,
        .(fold, split_id, n_train, n_test, sample_id, group_id, y)]
      data.table::setorder(reference, fold, sample_id, group_id, y)
      for (method_id in eligible) {
        candidate <- predictions[method == method_id &
                                   cohort == unit_value & seed == seed_value,
          .(fold, split_id, n_train, n_test, sample_id, group_id, y)]
        data.table::setorder(candidate, fold, sample_id, group_id, y)
        index <- index + 1L
        pass <- identical(candidate, reference)
        rows[[index]] <- data.table::data.table(
          method_id = method_id, unit_id = unit_value, seed = seed_value,
          reference_method = eligible[[1L]], pass = pass,
          split_id = if (pass) unique(candidate$split_id) else NA_character_,
          n_predictions = nrow(candidate),
          reason = if (pass) "exact" else "heldout_split_or_label_mismatch"
        )
        if (!pass) {
          stop("Methods do not share exact held-out rows/splits: ", unit_value,
               " seed ", seed_value, " method ", method_id, call. = FALSE)
        }
      }
    }
  }
  audit <- data.table::rbindlist(rows)
  if (!nrow(audit) || any(!audit$pass) ||
      length(setdiff(audit$method_id, methods))) {
    stop("Shared-split audit is incomplete.", call. = FALSE)
  }
  audit
}

.rta_auc <- function(y, score) {
  y <- as.integer(y)
  if (length(unique(y)) != 2L) return(NA_real_)
  pos <- score[y == 1L]
  neg <- score[y == 0L]
  ranks <- rank(c(pos, neg), ties.method = "average")
  (sum(ranks[seq_along(pos)]) - length(pos) * (length(pos) + 1L) / 2) /
    (length(pos) * length(neg))
}

.rta_group_predictions <- function(predictions) {
  out <- predictions[, {
    if (data.table::uniqueN(y) != 1L) {
      stop("A held-out provenance group has discordant labels.", call. = FALSE)
    }
    .(y = unique(y), score = mean(score_m), n_profiles = .N)
  }, by = .(method_id = method, unit_id = cohort, seed, fold, split_id,
            group_id)]
  if (any(!is.finite(out$score)) ||
      anyDuplicated(out[, .(method_id, unit_id, seed, fold, group_id)])) {
    stop("Group-collapsed held-out predictions are invalid.", call. = FALSE)
  }
  out[, stratum_id := paste(seed, fold, sep = "::")]
  out
}

.rta_splits <- function(group_predictions, units) {
  one <- unique(group_predictions[, .(unit_id, seed, fold, split_id, group_id,
                                      y)])
  method_counts <- group_predictions[, .(n_method = data.table::uniqueN(method_id)),
                                     by = .(unit_id, seed, fold, split_id,
                                            group_id, y)]
  if (any(method_counts$n_method < 1L)) stop("Internal split accounting failure.")
  one <- unique(method_counts[, .(unit_id, seed, fold, split_id, group_id, y)])
  out <- one[, .(
    n_test_groups = .N,
    n_positive_groups = sum(y == 1L), n_negative_groups = sum(y == 0L)
  ), by = .(unit_id, seed, fold, split_id)]
  out <- merge(out, units[, .(unit_id, n_unique_groups)], by = "unit_id",
               all.x = TRUE, sort = FALSE)
  out[, n_train_groups := n_unique_groups - n_test_groups]
  out[, stratum_id := paste(seed, fold, sep = "::")]
  duplicate_strata <- out[, anyDuplicated(stratum_id), by = unit_id]
  if (any(out$n_positive_groups < 1L) || any(out$n_negative_groups < 1L) ||
      any(out$n_train_groups < 1L) || any(duplicate_strata$V1 != 0L)) {
    stop("Canonical split table contains a one-class or invalid stratum.",
         call. = FALSE)
  }
  out
}

.rta_eligibility <- function(cells, methods, units) {
  out <- cells[, .(
    eligible = .N == length(.rta_seeds()) && all(eligible),
    eligible_seeds = sum(eligible), expected_seeds = length(.rta_seeds()),
    reason = if (all(eligible)) "eligible" else
      paste(unique(stats::na.omit(ineligible_reason[!eligible])), collapse = " | ")
  ), by = .(method_id = method, unit_id = cohort)]
  expected <- data.table::CJ(method_id = methods, unit_id = units$unit_id,
                             unique = TRUE)
  if (nrow(out) != nrow(expected) ||
      !setequal(do.call(paste, out[, .(method_id, unit_id)]),
                do.call(paste, expected))) {
    stop("Method-by-unit eligibility grid is incomplete.", call. = FALSE)
  }
  out
}

.rta_performance <- function(group_predictions, eligibility, splits, units) {
  strata <- group_predictions[, .(auc = .rta_auc(y, score)),
                              by = .(method_id, unit_id, seed, fold,
                                     split_id, stratum_id)]
  strata <- merge(strata,
                  splits[, .(unit_id, seed, fold, split_id,
                             n_test_groups, n_train_groups)],
                  by = c("unit_id", "seed", "fold", "split_id"),
                  all.x = TRUE, sort = FALSE)
  if (anyNA(strata[, .(auc, n_test_groups, n_train_groups)]) ||
      any(!is.finite(strata$auc))) {
    stop("An eligible group-collapsed CV stratum has undefined performance.",
         call. = FALSE)
  }
  rows <- lapply(seq_len(nrow(eligibility)), function(i) {
    e <- eligibility[i]
    if (!isTRUE(e$eligible)) {
      return(data.table::data.table(
        method_id = e$method_id, unit_id = e$unit_id, eligible = FALSE,
        auc = NA_real_, n_strata = 0L, reason = e$reason))
    }
    d <- strata[method_id == e$method_id & unit_id == e$unit_id]
    expected <- splits[unit_id == e$unit_id, stratum_id]
    if (nrow(d) != length(expected) ||
        !setequal(d$stratum_id, expected)) {
      stop("Eligible method/unit lacks complete expected strata: ",
           e$method_id, "/", e$unit_id, call. = FALSE)
    }
    data.table::data.table(
      method_id = e$method_id, unit_id = e$unit_id, eligible = TRUE,
      auc = mean(d$auc), n_strata = nrow(d), reason = "eligible")
  })
  performance <- data.table::rbindlist(rows)
  performance <- merge(performance, units[, .(
    unit_id, primary_unit, n_deposited_profiles, n_unique_groups)],
    by = "unit_id", all.x = TRUE, sort = FALSE)
  list(strata = strata, performance = performance)
}

.rta_pair_unit <- function(pair, unit, group_predictions, eligibility,
                           splits, units) {
  e <- eligibility[unit_id == unit & method_id %in%
                     c(pair$method_a, pair$method_b)]
  base <- data.table::data.table(
    pair_id = pair$pair_id, method_a = pair$method_a,
    method_b = pair$method_b, effect_direction = "method_a_minus_method_b",
    unit_id = unit, common_support = FALSE, estimate = NA_real_,
    se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_zero = NA_real_,
    p_above_margin = NA_real_, p_below_minus_margin = NA_real_,
    expected_m = NA_integer_, observed_m = NA_integer_, correction = NA_real_,
    mean_test_train_ratio = NA_real_, reason = "not_common_eligible"
  )
  if (nrow(e) != 2L || !all(e$eligible)) return(list(summary = base, strata = NULL))
  a <- group_predictions[method_id == pair$method_a & unit_id == unit]
  b <- group_predictions[method_id == pair$method_b & unit_id == unit]
  keys <- c("unit_id", "seed", "fold", "split_id", "stratum_id", "group_id",
            "y", "n_profiles")
  joined <- merge(a[, c(keys, "score"), with = FALSE],
                  b[, c(keys, "score"), with = FALSE], by = keys, all = FALSE,
                  sort = FALSE, suffixes = c("_a", "_b"))
  if (nrow(joined) != nrow(a) || nrow(joined) != nrow(b)) {
    stop("Pair lacks exact group-collapsed held-out prediction support: ",
         pair$pair_id, "/", unit, call. = FALSE)
  }
  expected_strata <- splits[unit_id == unit, stratum_id]
  matched <- OmicSelector::singlesample_matched_pair_auc(
    y = joined$y, score_a = joined$score_a, score_b = joined$score_b,
    stratum_id = joined$stratum_id,
    sample_id = paste(joined$stratum_id, joined$group_id, sep = "\r"),
    group_id = joined$group_id, expected_strata = expected_strata,
    margin = 0.05, conf_level = 0.95
  )
  s <- matched$summary
  base[, `:=`(
    common_support = TRUE, estimate = s$estimate, se = s$se,
    ci_low = s$ci_low, ci_high = s$ci_high,
    p_zero = s$p_zero_two_sided, p_above_margin = s$p_above_margin,
    p_below_minus_margin = s$p_below_minus_margin,
    expected_m = s$expected_m, observed_m = s$observed_m,
    correction = s$correction,
    mean_test_train_ratio = s$mean_test_train_ratio,
    reason = if (s$inference_available) "evaluable" else "zero_corrected_se"
  )]
  u <- units[unit_id == unit]
  base[, `:=`(primary_unit = u$primary_unit,
              n_deposited_profiles = u$n_deposited_profiles,
              n_unique_groups = u$n_unique_groups)]
  strata <- data.table::as.data.table(matched$strata)
  strata[, `:=`(
    pair_id = pair$pair_id, method_a = pair$method_a,
    method_b = pair$method_b, unit_id = unit,
    effect_direction = "method_a_minus_method_b")]
  list(summary = base, strata = strata)
}

.rta_meta_one <- function(pairwise, pair, population, units) {
  support <- sort(pairwise[common_support == TRUE, unit_id])
  row <- data.table::data.table(
    analysis_population = population, pair_id = pair$pair_id,
    method_a = pair$method_a, method_b = pair$method_b,
    effect_direction = "method_a_minus_method_b",
    k_support = length(support), k_meta = 0L,
    support_units = paste(support, collapse = ";"),
    n_deposited_profiles_support =
      sum(units[unit_id %in% support, n_deposited_profiles]),
    n_unique_groups_support = sum(units[unit_id %in% support, n_unique_groups]),
    estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    prediction_low = NA_real_, prediction_high = NA_real_, tau2 = NA_real_,
    i2 = NA_real_, p_zero = NA_real_, p_above_margin = NA_real_,
    p_below_minus_margin = NA_real_, df = NA_real_,
    pooled_testable = length(support) >= 5L,
    status = if (length(support) >= 5L) "meta_unavailable" else "descriptive_k_lt_5",
    fit_reason = if (length(support) >= 5L) "not_fitted" else "k_support_lt_5"
  )
  use <- pairwise[common_support & is.finite(estimate) & is.finite(se) & se > 0]
  row$k_meta <- nrow(use)
  if (nrow(use) < 5L) {
    if (length(support) >= 5L) row$fit_reason <- "fewer_than_five_finite_unit_effects"
    return(row)
  }
  fit <- tryCatch(withCallingHandlers(
    metafor::rma.uni(yi = use$estimate, sei = use$se, method = "REML",
                     test = "knha"),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)
  ), error = identity)
  if (inherits(fit, "error")) {
    row$fit_reason <- paste0("rma_failed:", conditionMessage(fit))
    return(row)
  }
  prediction <- tryCatch(metafor::predict.rma(fit, level = 95), error = identity)
  if (inherits(prediction, "error") || !is.finite(prediction$pi.lb) ||
      !is.finite(prediction$pi.ub)) {
    row$fit_reason <- "prediction_interval_failed"
    return(row)
  }
  estimate <- as.numeric(fit$b[[1L]])
  se <- as.numeric(fit$se[[1L]])
  df <- nrow(use) - 1L
  row[, `:=`(
    estimate = estimate, se = se, ci_low = as.numeric(fit$ci.lb[[1L]]),
    ci_high = as.numeric(fit$ci.ub[[1L]]),
    prediction_low = as.numeric(prediction$pi.lb[[1L]]),
    prediction_high = as.numeric(prediction$pi.ub[[1L]]),
    tau2 = as.numeric(fit$tau2), i2 = as.numeric(fit$I2),
    p_zero = as.numeric(fit$pval[[1L]]),
    p_above_margin = stats::pt((estimate - 0.05) / se, df = df,
                               lower.tail = FALSE),
    p_below_minus_margin = stats::pt((estimate + 0.05) / se, df = df,
                                     lower.tail = TRUE),
    df = df, status = "REML_KH", fit_reason = "eligible")]
  row
}

.rta_meta_all <- function(pairwise, pairs, units) {
  populations <- list(
    primary_circulating_biofluid = units[primary_unit == TRUE, unit_id],
    all_biofluids_sensitivity = units$unit_id
  )
  rows <- list()
  for (population in names(populations)) {
    d <- pairwise[unit_id %in% populations[[population]]]
    for (i in seq_len(nrow(pairs))) {
      rows[[length(rows) + 1L]] <- .rta_meta_one(
        d[pair_id == pairs$pair_id[[i]]], pairs[i], population, units)
    }
  }
  out <- data.table::rbindlist(rows, fill = TRUE)
  frozen <- pairs$pair_id
  for (population in names(populations)) {
    index <- which(out$analysis_population == population)
    zero <- OmicSelector::singlesample_adjust_zero_by(
      out$pair_id[index], out$p_zero[index], frozen)
    relevance <- OmicSelector::singlesample_adjust_relevance_by(
      out$pair_id[index], out$p_above_margin[index],
      out$p_below_minus_margin[index], frozen)
    out[index, `:=`(
      p_zero_by = zero$p_zero_by,
      p_above_margin_by = relevance$p_above_margin_by,
      p_below_minus_margin_by = relevance$p_below_minus_margin_by,
      zero_family_n = zero$family_n,
      relevance_family_n = relevance$family_n)]
  }
  out
}

.rta_complete_support <- function(eligibility, performance, roster, units) {
  populations <- list(
    primary_circulating_biofluid = units[primary_unit == TRUE, unit_id],
    all_biofluids_sensitivity = units$unit_id
  )
  table_rows <- rectangle_rows <- list()
  for (population in names(populations)) {
    panel <- OmicSelector::singlesample_complete_support_panels(
      eligibility, roster = roster, thresholds = c(0.50, 0.80, 0.95),
      primary_units = populations[[population]], min_methods = 5L,
      min_units = 5L)
    summary <- data.table::as.data.table(panel$summary)
    summary[, `:=`(
      analysis_population = population, item_type = "summary",
      rank_allowed = !degenerate,
      conclusion_allowed = !degenerate,
      n_deposited_profiles = vapply(threshold, function(value) {
        selected <- panel$unit_membership$unit_id[
          panel$unit_membership$threshold == value &
            panel$unit_membership$included_unit]
        sum(units[unit_id %in% selected, n_deposited_profiles])
      }, numeric(1L)),
      n_unique_groups = vapply(threshold, function(value) {
        selected <- panel$unit_membership$unit_id[
          panel$unit_membership$threshold == value &
            panel$unit_membership$included_unit]
        sum(units[unit_id %in% selected, n_unique_groups])
      }, numeric(1L)))]
    methods <- data.table::as.data.table(panel$method_membership)
    methods[, `:=`(analysis_population = population, item_type = "method")]
    panel_units <- data.table::as.data.table(panel$unit_membership)
    panel_units[, `:=`(analysis_population = population, item_type = "unit")]
    table_rows[[length(table_rows) + 1L]] <- data.table::rbindlist(
      list(summary, methods, panel_units), fill = TRUE)
    for (threshold_value in panel$summary$threshold) {
      degenerate <- summary[threshold == threshold_value, degenerate]
      if (length(degenerate) != 1L || is.na(degenerate)) {
        stop("Complete-support panel has an undefined degeneracy state.",
             call. = FALSE)
      }
      selected_methods <- panel$method_membership$method_id[
        panel$method_membership$threshold == threshold_value &
          panel$method_membership$included_method]
      selected_units <- panel$unit_membership$unit_id[
        panel$unit_membership$threshold == threshold_value &
          panel$unit_membership$included_unit]
      rectangle <- performance[method_id %in% selected_methods &
                                 unit_id %in% selected_units]
      if (nrow(rectangle)) {
        if (any(!rectangle$eligible) || any(!is.finite(rectangle$auc)) ||
            nrow(rectangle) != length(selected_methods) * length(selected_units)) {
          stop("Complete-support rectangle contains missing performance.",
               call. = FALSE)
        }
        rectangle[, `:=`(
          analysis_population = population, threshold = threshold_value,
          rank_allowed = !degenerate,
          conclusion_allowed = !degenerate,
          within_unit_rank = if (degenerate) {
            NA_real_
          } else rank(-auc, ties.method = "average")), by = unit_id]
        rectangle_rows[[length(rectangle_rows) + 1L]] <- rectangle
      }
    }
  }
  list(table = data.table::rbindlist(table_rows, fill = TRUE),
       rectangle = data.table::rbindlist(rectangle_rows, fill = TRUE))
}

.rta_coverage <- function(eligibility, performance, roster, units) {
  populations <- list(
    primary_circulating_biofluid = units[primary_unit == TRUE, unit_id],
    all_biofluids_sensitivity = units$unit_id
  )
  data.table::rbindlist(lapply(names(populations), function(population) {
    ids <- populations[[population]]
    out <- eligibility[unit_id %in% ids, .(
      eligible_units = sum(eligible), total_units = .N, coverage = mean(eligible),
      eligibility_units = paste(unit_id[eligible], collapse = ";")),
      by = method_id]
    auc <- performance[unit_id %in% ids & eligible == TRUE, .(
      median_auc = stats::median(auc), mean_auc = mean(auc)), by = method_id]
    out <- merge(out, auc, by = "method_id", all.x = TRUE, sort = FALSE)
    out <- merge(roster[, .(method_id, roster_order = .I, family, role, tier,
                            dep_route)], out, by = "method_id", sort = FALSE)
    out[, `:=`(
      analysis_population = population,
      n_deposited_profiles = sum(units[unit_id %in% ids,
                                       n_deposited_profiles]),
      n_unique_groups = sum(units[unit_id %in% ids, n_unique_groups]))]
    out
  }), fill = TRUE)
}

.rta_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("The sensitivity assembler requires metafor.", call. = FALSE)
  }
  opt <- .rta_args(args)
  code_start <- .rta_git_provenance(opt$package_root)
  if (code_start$dirty || !identical(code_start$commit,
                                     opt$expected_package_commit)) {
    stop("Assembly requires the exact clean OmicSelector package checkout.",
         call. = FALSE)
  }
  description <- read.dcf(file.path(code_start$root, "DESCRIPTION"))
  if (!identical(unname(description[, "Version"]),
                 opt$expected_package_version)) {
    stop("OmicSelector DESCRIPTION version differs from its expected pin.",
         call. = FALSE)
  }
  script_path <- .rta_script_path()
  expected_script <- normalizePath(file.path(
    code_start$root, "inst", "scripts", basename(script_path)), mustWork = TRUE)
  if (!identical(script_path, expected_script)) {
    stop("Assembler must execute from the pinned OmicSelector checkout.",
         call. = FALSE)
  }
  assembler_sha <- .rta_sha(script_path)
  runner_path <- .rta_assert_file_pin(
    file.path(code_start$root, "inst", "scripts",
              "paper3_run_ranktree_sensitivity_bundle.R"),
    opt$expected_runner_script_sha256, "Package sensitivity runner"
  )
  launcher_path <- .rta_assert_file_pin(
    file.path(code_start$root, "inst", "scripts",
              "paper3_run_ranktree_sensitivity_container.sh"),
    opt$expected_container_launcher_sha256, "Package container launcher"
  )
  runtime_image <- .rta_assert_file_pin(
    opt$runtime_image, opt$expected_runtime_image_sha256,
    "Singularity runtime image"
  )
  if (!identical(R.version.string, opt$expected_r_version)) {
    stop("In-process R version differs from its exact pre-execution pin.",
         call. = FALSE)
  }
  runtime_env <- new.env(parent = globalenv())
  sys.source(runner_path, envir = runtime_env)
  paper_root <- normalizePath(opt$paper_root, mustWork = TRUE)
  analysis_plan <- runtime_env$.ranktree_validate_analysis_plan(
    opt$analysis_plan, paper_root, opt$expected_analysis_plan_sha256
  )
  analysis_amendment <- runtime_env$.ranktree_validate_analysis_amendment(
    opt$analysis_amendment, paper_root,
    opt$expected_analysis_amendment_sha256
  )
  runtime_attestation <- runtime_env$.ranktree_validate_runtime_attestation(
    opt$runtime_attestation_manifest,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_r_version, opt$expected_runtime_image_sha256,
    opt$expected_container_launcher_sha256
  )
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The exact sensitivity assembler requires pkgload.", call. = FALSE)
  }
  suppressMessages(suppressWarnings(pkgload::load_all(
    code_start$root, quiet = TRUE, export_all = FALSE, helpers = FALSE)))
  if (!identical(as.character(utils::packageVersion("OmicSelector")),
                 opt$expected_package_version)) {
    stop("Loaded OmicSelector version differs from the pinned checkout.",
         call. = FALSE)
  }
  base <- OmicSelector::singlesample_method_roster()
  anchor <- base[match(.rta_method_ids()[1:5], base$method_id), ]
  external <- OmicSelector::singlesample_external_competitor_roster()
  roster <- data.table::rbindlist(list(anchor, external), fill = TRUE)
  roster[, roster_order := .I]
  if (!identical(roster$method_id, .rta_method_ids()) || nrow(roster) != 7L ||
      any(roster$estimand != "within")) {
    stop("Installed package does not expose the frozen seven-method registry.",
         call. = FALSE)
  }
  pairs <- data.table::as.data.table(
    OmicSelector::singlesample_within_method_pairs(roster))
  if (nrow(pairs) != 21L || anyDuplicated(pairs$pair_id) ||
      !identical(pairs$method_a,
                 roster$method_id[pairs$roster_order_a]) ||
      !identical(pairs$method_b,
                 roster$method_id[pairs$roster_order_b])) {
    stop("Could not freeze the registry-ordered 21-pair family.",
         call. = FALSE)
  }

  expected <- list(
    package_version = opt$expected_package_version,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_sha256 = opt$expected_analysis_amendment_sha256,
    package_commit = opt$expected_package_commit,
    cache_commit = opt$expected_cache_commit,
    cache_sha256 = opt$expected_cache_sha256,
    cache_manifest_sha256 = opt$expected_cache_manifest_sha256,
    cache_digest_sha256 = opt$expected_cache_digest_sha256,
    runner_script_sha256 = opt$expected_runner_script_sha256,
    engine_manifest_sha256 = opt$expected_engine_manifest_sha256,
    environment_manifest_sha256 = opt$expected_environment_manifest_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version,
    provenance_script_sha256 = opt$expected_provenance_script_sha256,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256,
    analysis_code_id = opt$expected_analysis_code_id
  )
  file_receipts <- lapply(opt$bundles, .rta_validate_output_manifest)
  bundles <- lapply(file_receipts, function(x) readRDS(x$path))
  observed_methods <- vapply(bundles, function(x) {
    if (!is.character(x$methods) || length(x$methods) != 1L) return(NA_character_)
    x$methods
  }, character(1L))
  if (anyNA(observed_methods) || !setequal(observed_methods, roster$method_id)) {
    stop("Seven input bundles do not cover the exact registered methods once.",
         call. = FALSE)
  }
  order <- match(roster$method_id, observed_methods)
  bundles <- bundles[order]
  file_receipts <- file_receipts[order]
  units <- NULL
  validated <- vector("list", 7L)
  for (i in seq_len(7L)) {
    validated[[i]] <- .rta_validate_bundle(
      bundles[[i]], roster$method_id[[i]], expected, units = units)
    if (is.null(units)) units <- validated[[i]]$units
    registry_fields <- c("method_id", "family", "estimand", "role", "tier",
                         "dep_route", "fit_fn", "score_fn", "pkg_status",
                         "row_source", "lopbo_mechanism")
    observed_registry <- validated[[i]]$method_registry[, ..registry_fields]
    expected_registry <- roster[i, ..registry_fields]
    if (!identical(observed_registry, expected_registry)) {
      stop("Method bundle registry metadata differs from the current frozen ",
           "package registry: ", roster$method_id[[i]], call. = FALSE)
    }
    for (file in c("analysis_plan.md", "analysis_amendment.md",
                   "engine_manifest.tsv",
                   "environment_manifest.tsv", "runtime_attestation.tsv")) {
      path <- file.path(file_receipts[[i]]$root, file)
      expected_file_sha <- if (file == "analysis_plan.md") {
        opt$expected_analysis_plan_sha256
      } else if (file == "analysis_amendment.md") {
        opt$expected_analysis_amendment_sha256
      } else if (file == "engine_manifest.tsv") {
        opt$expected_engine_manifest_sha256
      } else if (file == "environment_manifest.tsv") {
        opt$expected_environment_manifest_sha256
      } else opt$expected_runtime_attestation_manifest_sha256
      if (!file.exists(path) || !identical(.rta_sha(path), expected_file_sha)) {
        stop("Bundle carries a mismatched pinned input artifact: ", file,
             " for ", roster$method_id[[i]], call. = FALSE)
      }
    }
  }
  if (nrow(units) != 34L || sum(units$primary_unit) != 33L) {
    stop("Sensitivity bundles do not retain the exact 34/33 unit population.",
         call. = FALSE)
  }
  closure <- bundles[[1L]]$ranktree_dependency_closure
  if (any(!vapply(bundles[-1L], function(x) {
    identical(x$ranktree_dependency_closure, closure) &&
      identical(x$r_version, bundles[[1L]]$r_version)
  }, logical(1L)))) {
    stop("Method bundles were not executed under one exact R/dependency ",
         "environment.", call. = FALSE)
  }
  cells <- data.table::rbindlist(lapply(validated, `[[`, "cells"), fill = TRUE)
  predictions <- data.table::rbindlist(
    lapply(validated, `[[`, "predictions"), fill = TRUE)
  data.table::setorder(cells, method, cohort, seed)
  data.table::setorder(predictions, method, cohort, seed, fold, sample_id)
  split_audit <- .rta_common_split_audit(predictions, cells, roster$method_id)
  eligibility <- .rta_eligibility(cells, roster$method_id, units)
  group_predictions <- .rta_group_predictions(predictions)
  splits <- .rta_splits(group_predictions, units)
  performance_result <- .rta_performance(
    group_predictions, eligibility, splits, units)
  performance <- performance_result$performance
  coverage <- .rta_coverage(eligibility, performance, roster, units)

  pair_rows <- stratum_rows <- list()
  for (i in seq_len(nrow(pairs))) {
    for (unit in units$unit_id) {
      result <- .rta_pair_unit(pairs[i], unit, group_predictions,
                               eligibility, splits, units)
      pair_rows[[length(pair_rows) + 1L]] <- result$summary
      if (!is.null(result$strata)) {
        stratum_rows[[length(stratum_rows) + 1L]] <- result$strata
      }
    }
  }
  pairwise <- data.table::rbindlist(pair_rows, fill = TRUE)
  pair_strata <- data.table::rbindlist(stratum_rows, fill = TRUE)
  if (nrow(pairwise) != 21L * 34L) {
    stop("Pair-by-unit output is not the complete 21-by-34 grid.",
         call. = FALSE)
  }
  pair_meta <- .rta_meta_all(pairwise, pairs, units)
  if (nrow(pair_meta) != 42L ||
      any(pair_meta$zero_family_n != 21L) ||
      any(pair_meta$relevance_family_n != 42L)) {
    stop("Random-effects/BY output does not retain the frozen 21-pair family.",
         call. = FALSE)
  }
  panels <- .rta_complete_support(eligibility, performance, roster, units)

  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) stop("Output bundle already exists: ", output_dir)
  stage <- tempfile(".ranktree-assembly-", tmpdir = output_parent)
  dir.create(stage)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  data.table::fwrite(roster, file.path(stage, "methods.tsv"), sep = "\t")
  data.table::fwrite(units, file.path(stage, "units.tsv"), sep = "\t")
  data.table::fwrite(cells, file.path(stage, "cells.tsv"), sep = "\t")
  data.table::fwrite(eligibility, file.path(stage, "eligibility.tsv"), sep = "\t")
  data.table::fwrite(group_predictions,
                     file.path(stage, "heldout_group_predictions.tsv.gz"),
                     sep = "\t", compress = "gzip")
  data.table::fwrite(splits, file.path(stage, "splits.tsv"), sep = "\t")
  data.table::fwrite(split_audit, file.path(stage, "shared_split_audit.tsv"),
                     sep = "\t")
  data.table::fwrite(performance_result$strata,
                     file.path(stage, "method_stratum_performance.tsv"),
                     sep = "\t")
  data.table::fwrite(performance,
                     file.path(stage, "method_unit_performance.tsv"), sep = "\t")
  data.table::fwrite(coverage, file.path(stage, "coverage_performance.tsv"),
                     sep = "\t")
  data.table::fwrite(pairs, file.path(stage, "pair_family.tsv"), sep = "\t")
  data.table::fwrite(pairwise, file.path(stage, "pairwise_unit.tsv"), sep = "\t")
  data.table::fwrite(pair_strata, file.path(stage, "pairwise_strata.tsv"),
                     sep = "\t")
  data.table::fwrite(pair_meta, file.path(stage, "pairwise_meta.tsv"), sep = "\t")
  data.table::fwrite(panels$table,
                     file.path(stage, "complete_support_panels.tsv"), sep = "\t")
  data.table::fwrite(panels$rectangle,
                     file.path(stage, "complete_support_performance.tsv"),
                     sep = "\t")
  input_manifest <- data.table::rbindlist(lapply(seq_len(7L), function(i) {
    data.table::data.table(
      method_id = roster$method_id[[i]], path = file_receipts[[i]]$path,
      bytes = file.info(file_receipts[[i]]$path)$size,
      sha256 = file_receipts[[i]]$bundle_sha256,
      completion_manifest_sha256 = file_receipts[[i]]$manifest_sha256)
  }))
  data.table::fwrite(input_manifest,
                     file.path(stage, "input_bundle_manifest.tsv"), sep = "\t")
  common_files <- c("analysis_plan.md", "analysis_amendment.md",
                    "engine_manifest.tsv",
                    "environment_manifest.tsv", "runtime_attestation.tsv")
  for (file in common_files) {
    source <- file.path(file_receipts[[1L]]$root, file)
    expected_sha <- if (file == "analysis_plan.md") {
      opt$expected_analysis_plan_sha256
    } else if (file == "analysis_amendment.md") {
      opt$expected_analysis_amendment_sha256
    } else if (file == "engine_manifest.tsv") {
      opt$expected_engine_manifest_sha256
    } else if (file == "environment_manifest.tsv") {
      opt$expected_environment_manifest_sha256
    } else opt$expected_runtime_attestation_manifest_sha256
    if (!file.exists(source) || !identical(.rta_sha(source), expected_sha) ||
        !file.copy(source, file.path(stage, file), copy.mode = TRUE)) {
      stop("Could not validate/copy common bundle input: ", file,
           call. = FALSE)
    }
  }
  provenance_artifacts <- data.table::rbindlist(lapply(seq_len(7L), function(i) {
    paths <- file.path(file_receipts[[i]]$root,
                       c("provenance_preflight.log", "provenance_resolution.tsv"))
    if (any(!file.exists(paths))) stop("A bundle lacks provenance artifacts.")
    data.table::data.table(
      method_id = roster$method_id[[i]], file = basename(paths),
      path = normalizePath(paths), bytes = file.info(paths)$size,
      sha256 = vapply(paths, .rta_sha, character(1L)))
  }))
  data.table::fwrite(provenance_artifacts,
                     file.path(stage, "provenance_artifact_manifest.tsv"),
                     sep = "\t")
  writeLines(c(
    "# External rank-tree competitor sensitivity", "",
    paste0("Registered methods: ", nrow(roster)),
    paste0("Frozen ordered comparisons (method A minus method B): ", nrow(pairs)),
    paste0("Primary circulating-biofluid units: ", sum(units$primary_unit)),
    paste0("All-biofluid units: ", nrow(units)),
    paste0("Deposited profiles represented: ", sum(units$n_deposited_profiles)),
    paste0("Deduplicated biological/provenance groups: ",
           sum(units$n_unique_groups)),
    paste0("Singularity runtime image SHA-256: ",
           opt$expected_runtime_image_sha256),
    paste0("Exact in-process R runtime/library attestation SHA-256: ",
           opt$expected_runtime_attestation_manifest_sha256),
    paste0("Package container launcher SHA-256: ",
           opt$expected_container_launcher_sha256),
    paste0("Approved analysis plan SHA-256: ",
           opt$expected_analysis_plan_sha256),
    paste0("Approved analysis amendment SHA-256: ",
           opt$expected_analysis_amendment_sha256),
    "Rank forest contract: unweighted preliminary screen; inverse-frequency ",
    "pair-forest case weights; numeric y=1 (`case`) is the positive target.",
    paste0("Group-collapsed held-out AUC: one mean score per provenance group ",
           "within each exact test fold."),
    paste0("Primary REML/Knapp--Hartung fits: ",
           sum(pair_meta$analysis_population == "primary_circulating_biofluid" &
                 pair_meta$status == "REML_KH")),
    paste0("Non-degenerate complete-support panels: ",
           sum(panels$table$item_type == "summary" &
                 !panels$table$degenerate, na.rm = TRUE)), "",
    "All 21 pair directions were frozen by registry order and adjusted with ",
    "Benjamini--Yekutieli (21 zero-difference tests; 42 joint directional ",
    "relevance tests). This exploratory bundle is not reportable until ",
    "independent results review, plausibility assessment and promotion."
  ), file.path(stage, "report.md"))

  code_end <- .rta_git_provenance(opt$package_root)
  runtime_env$.ranktree_validate_runtime_attestation(
    runtime_attestation$path,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_r_version, opt$expected_runtime_image_sha256,
    opt$expected_container_launcher_sha256
  )
  current_input_sha <- vapply(file_receipts, function(x) .rta_sha(x$path),
                              character(1L))
  if (!identical(code_end, code_start) ||
      !identical(.rta_sha(script_path), assembler_sha) ||
      !identical(.rta_sha(analysis_plan),
                 opt$expected_analysis_plan_sha256) ||
      !identical(.rta_sha(analysis_amendment),
                 opt$expected_analysis_amendment_sha256) ||
      !identical(.rta_sha(runtime_image),
                 opt$expected_runtime_image_sha256) ||
      !identical(.rta_sha(launcher_path),
                 opt$expected_container_launcher_sha256) ||
      !identical(unname(current_input_sha), input_manifest$sha256)) {
    stop("Package code or an input bundle changed during assembly.",
         call. = FALSE)
  }
  artifacts <- list.files(stage, full.names = TRUE)
  output_manifest <- data.table::data.table(
    file = basename(artifacts), bytes = file.info(artifacts)$size,
    sha256 = vapply(artifacts, .rta_sha, character(1L)),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = opt$expected_package_version,
    package_commit = opt$expected_package_commit,
    assembler_script_sha256 = assembler_sha,
    runner_script_sha256 = opt$expected_runner_script_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_sha256 = opt$expected_analysis_amendment_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256,
    analysis_code_id = opt$expected_analysis_code_id,
    pair_family_n = 21L, zero_by_family_n = 21L,
    directional_by_family_n = 42L
  )
  data.table::fwrite(output_manifest,
                     file.path(stage, "output_manifest.tsv"), sep = "\t")
  if (!file.rename(stage, output_dir)) {
    stop("Could not atomically promote sensitivity assembly.", call. = FALSE)
  }
  cat("Created external rank-tree sensitivity assembly:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .rta_main()
