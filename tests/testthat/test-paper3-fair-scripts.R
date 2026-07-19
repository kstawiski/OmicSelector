paper3_script_env <- function(script) {
  env <- new.env(parent = globalenv())
  source_path <- testthat::test_path("..", "..", "inst", "scripts", script)
  if (!file.exists(source_path)) {
    source_path <- system.file("scripts", script, package = "OmicSelector")
  }
  if (!nzchar(source_path) || !file.exists(source_path)) {
    stop("Could not locate installed paper-analysis script: ", script,
         call. = FALSE)
  }
  env$.PAPER3_ACTIVE_SCRIPT_PATH <- normalizePath(source_path, mustWork = TRUE)
  sys.source(source_path, envir = env)
  env
}

test_that("snapshot overlay replaces only declared registered keys", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  base_cells <- data.table::data.table(
    method = c("m1", "m2"), cohort = "u1", seed = 1L,
    eligible = c(TRUE, FALSE), ineligible_reason = c(NA, "old")
  )
  base_predictions <- data.table::data.table(
    method = "m1", cohort = "u1", seed = 1L, fold = 1L,
    sample_id = "s1", group_id = "g1", y = 0L, score_m = 0.1
  )
  td <- tempfile("overlay-test-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  cells_path <- file.path(td, "cells.tsv")
  predictions_path <- file.path(td, "predictions.tsv")
  data.table::fwrite(data.table::data.table(
    method = "m2", cohort = "u1", seed = 1L, eligible = TRUE,
    ineligible_reason = NA_character_
  ), cells_path, sep = "\t")
  data.table::fwrite(data.table::data.table(
    method = "m2", cohort = "u1", seed = 1L, fold = 1L,
    sample_id = "s1", group_id = "g1", y = 0L, score_m = 0.2
  ), predictions_path, sep = "\t")

  out <- env$.snapshot_overlay(
    base_cells, base_predictions, cells_path, predictions_path
  )
  expect_equal(nrow(out$cells), 2L)
  expect_true(all(out$cells$eligible))
  expect_setequal(out$predictions$method, c("m1", "m2"))
  expect_equal(out$predictions[method == "m1", score_m], 0.1)
  expect_equal(out$predictions[method == "m2", score_m], 0.2)
})

test_that("snapshot prediction normalization rejects duplicate held-out rows", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  d <- data.table::data.table(
    method = rep("m1", 2L), cohort = rep("u1", 2L), seed = 1L,
    fold = 1L, sample_id = rep("s1", 2L), group_id = rep("g1", 2L),
    y = 0L, score_m = c(0.1, 0.2)
  )
  expect_error(env$.snapshot_normalize_predictions(d), "duplicate")
})

test_that("snapshot validates fixed-direction score orientation", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  predictions <- data.table::data.table(
    method = "m1", cohort = "u1", seed = 1L,
    fold = rep(1:2, each = 4L), sample_id = paste0("s", 1:8),
    group_id = paste0("g", 1:8), y = rep(c(0L, 0L, 1L, 1L), 2L),
    score_m = rep(c(0.1, 0.2, 0.8, 0.9), 2L)
  )
  cells <- data.table::data.table(
    method = "m1", cohort = "u1", seed = 1L,
    eligible = TRUE, auc_method = 1
  )
  expect_true(env$.snapshot_validate_fixed_auc(cells, predictions))
  cells$auc_method <- 0
  expect_error(env$.snapshot_validate_fixed_auc(cells, predictions),
               "do not reproduce")
})

test_that("paper AUC helpers use exact tie-aware midranks", {
  snapshot <- paper3_script_env("paper3_build_fair_snapshot.R")
  fair <- paper3_script_env("paper3_fair_method_comparison.R")
  y <- c(0L, 0L, 0L, 1L, 1L)
  score <- c(0, 1, 1, 1, 2)
  brute <- mean(outer(score[y == 1L], score[y == 0L],
                      function(a, b) (a > b) + 0.5 * (a == b)))
  expect_equal(snapshot$.snapshot_auc(y, score), brute, tolerance = 0)
  expect_equal(fair$.fair_auc(y, score), brute, tolerance = 0)
})

test_that("paper scripts compute one exact installed-package tree digest", {
  snapshot <- paper3_script_env("paper3_build_fair_snapshot.R")
  fair <- paper3_script_env("paper3_fair_method_comparison.R")
  validator <- paper3_script_env("paper3_validate_fair_method_comparison.R")
  root <- system.file(package = "OmicSelector")
  expect_identical(snapshot$.snapshot_package_tree_sha256(root),
                   fair$.fair_package_tree_sha256(root))
  expect_identical(snapshot$.snapshot_package_tree_sha256(root),
                   validator$.fairqa_package_tree_sha256(root))
})

test_that("snapshot upstream receipt is fail closed", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("upstream-receipt-")
  on.exit(unlink(td), add = TRUE)
  checks <- c(
    "within_auc_recomputed",
    "within_complete_prespecified_fold_evaluability",
    "within_eligible_prediction_accounting",
    "within_no_group_split",
    "within_prediction_all_prespecified_folds_per_eligible_seed",
    "within_prediction_group_labels_consistent",
    "within_prediction_groups_not_split",
    "within_prediction_key_unique",
    "within_prediction_route_scoped_pin",
    "within_prediction_scores_finite",
    "within_route_scoped_pin"
  )
  receipt <- data.table::data.table(check_id = checks, pass = TRUE)
  data.table::fwrite(receipt, td, sep = "\t")
  expect_invisible(env$.snapshot_validate_upstream_receipt(td))
  receipt <- data.table::rbindlist(list(
    receipt,
    data.table::data.table(check_id = "unrelated_but_required", pass = FALSE)
  ))
  data.table::fwrite(receipt, td, sep = "\t")
  expect_error(env$.snapshot_validate_upstream_receipt(td),
               "does not pass")
  receipt <- receipt[check_id != "unrelated_but_required"]
  receipt[check_id == "within_no_group_split", pass := FALSE]
  data.table::fwrite(receipt, td, sep = "\t")
  expect_error(env$.snapshot_validate_upstream_receipt(td),
               "does not pass")
})

test_that("snapshot binds every upstream validation input byte", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("validation-input-manifest-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  input <- file.path(td, "input.tsv")
  writeLines("frozen", input)
  manifest <- file.path(td, "validation_input_manifest.tsv")
  data.table::fwrite(data.table::data.table(
    path = "input.tsv", bytes = file.info(input)$size,
    sha256 = digest::digest(input, algo = "sha256", file = TRUE)
  ), manifest, sep = "\t")
  expect_invisible(env$.snapshot_validate_input_manifest(manifest))
  writeLines("changed", input)
  expect_error(env$.snapshot_validate_input_manifest(manifest),
               "differ from their manifest")
})

test_that("snapshot requires the exact canonical provenance union inventory", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("provenance-sources-")
  dir.create(file.path(td, "data_bucket"), recursive = TRUE)
  dir.create(file.path(td, "knowledge"))
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  manifest <- file.path(td, "data_bucket", "_PROVENANCE_MANIFEST.md")
  union <- file.path(td, "knowledge", "union_inventory.tsv")
  writeLines("manifest", manifest)
  data.table::fwrite(data.table::data.table(accession = c("A", "B")),
                     union, sep = "\t")
  pins <- vapply(c(manifest, union), digest::digest, character(1L),
                 algo = "sha256", file = TRUE)
  expect_length(env$.snapshot_validate_provenance_sources(
    td, pins[[1L]], pins[[2L]]
  ), 2L)
  unlink(union)
  expect_error(env$.snapshot_validate_provenance_sources(
    td, pins[[1L]], pins[[2L]]
  ), "missing or unpinned")
})

test_that("snapshot accepts only a fully passing pinned overlay candidate", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("overlay-candidate-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  splitter_methods <- data.table::as.data.table(singlesample_method_roster())[
    estimand == "within" & tier == "R2", method_id
  ]
  files <- c(
    "within_bundle.rds", "within_cells.tsv", "within_predictions.tsv.gz",
    "transfer_bundle.rds", "transfer_cells.tsv",
    "transfer_predictions.tsv.gz", "acceptance_checks.tsv",
    "route_environment_map.tsv", "structural_witness.rds",
    "inputs/within_tasks.tsv", "inputs/transfer_tasks.tsv",
    "transfer_shard_manifest.tsv",
    "legacy_compatibility_overlay_input_manifest.tsv",
    "source_input_manifest.tsv", env$.snapshot_candidate_environment_files(),
    paste0("splitter_manifests/output_manifest__", splitter_methods, ".tsv")
  )
  for (file in files) {
    dir.create(dirname(file.path(td, file)), recursive = TRUE,
               showWarnings = FALSE)
    writeLines(file, file.path(td, file))
  }
  data.table::fwrite(data.table::data.table(
    check_id = env$.snapshot_candidate_check_ids(), pass = TRUE
  ), file.path(td, "acceptance_checks.tsv"), sep = "\t")
  data.table::fwrite(data.table::data.table(
    execution_route = "exact", analysis_surface = c("within", "transfer"),
    environment_route = c("within_s5", "transfer_argos"),
    environment_receipt_path = c(
      "environment/within_s5", "environment/transfer_argos"
    ),
    environment_receipt_sha256 = c(strrep("d", 64L), strrep("e", 64L))
  ), file.path(td, "route_environment_map.tsv"), sep = "\t")
  within_methods <- data.table::as.data.table(singlesample_method_roster())[
    estimand == "within", method_id
  ]
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  within_cells <- data.table::CJ(
    method = within_methods, cohort = paste0("u", seq_len(34L)),
    seed = seeds, sorted = FALSE
  )
  within_cells[, `:=`(
    eligible = TRUE, ineligible_reason = NA_character_, auc_method = 0.5
  )]
  within_cells[1L, `:=`(eligible = FALSE, ineligible_reason = "structural")]
  within_predictions <- within_cells[eligible == TRUE, .(
    fold = 1L,
    sample_id = paste(cohort, seed, paste0("p", 1:4), sep = "-"),
    sample_idx = 1:4,
    group_id = paste(cohort, seed, paste0("g", 1:4), sep = "-"),
    y = c(0L, 0L, 1L, 1L), score_m = 0, score_b = 0
  ), by = .(method, cohort, seed)]
  within_bundle <- list(cells = data.table::copy(within_cells),
                        predictions = data.table::copy(within_predictions))
  saveRDS(within_bundle, file.path(td, "within_bundle.rds"), version = 3L)
  data.table::fwrite(within_cells, file.path(td, "within_cells.tsv"), sep = "\t")
  data.table::fwrite(within_predictions,
                     file.path(td, "within_predictions.tsv.gz"), sep = "\t",
                     compress = "gzip")

  transfer_methods <- data.table::as.data.table(singlesample_method_roster())[
    estimand == "transfer", method_id
  ]
  targets <- data.table::data.table(
    disease = paste0("d", seq_len(18L)),
    held_out = paste0("h", seq_len(18L))
  )
  transfer_cells <- data.table::rbindlist(lapply(seq_len(nrow(targets)),
    function(i) cbind(
      targets[i], data.table::CJ(seed = seeds, method = transfer_methods,
                                 sorted = FALSE)
    )
  ))
  transfer_cells[, `:=`(eligible = TRUE, ineligible_reason = NA_character_)]
  transfer_cells[1L, `:=`(eligible = FALSE, ineligible_reason = "structural")]
  transfer_predictions <- transfer_cells[eligible == TRUE, .(
    disease, held_out, seed, method,
    sample_id = paste0(held_out, "-sample"),
    group_id = paste0(held_out, "-group"), y = 0L, score_m = 0, score_b = 0
  )]
  transfer_bundle <- list(cells = data.table::copy(transfer_cells),
                          predictions = data.table::copy(transfer_predictions))
  saveRDS(transfer_bundle, file.path(td, "transfer_bundle.rds"), version = 3L)
  data.table::fwrite(transfer_cells, file.path(td, "transfer_cells.tsv"), sep = "\t")
  data.table::fwrite(transfer_predictions,
                     file.path(td, "transfer_predictions.tsv.gz"), sep = "\t",
                     compress = "gzip")

  source_root <- tempfile("overlay-source-")
  dir.create(source_root)
  on.exit(unlink(source_root, recursive = TRUE), add = TRUE)
  make_sources <- function(directory, n) {
    root <- file.path(source_root, directory)
    dir.create(root, recursive = TRUE)
    paths <- file.path(root, sprintf("source-%04d.bin", seq_len(n)))
    for (i in seq_along(paths)) writeLines(paste(directory, i), paths[[i]])
    normalizePath(paths)
  }
  source_rows <- function(role, paths) data.table::data.table(
    input_role = role, path = normalizePath(paths),
    bytes = unname(file.info(paths)$size),
    sha256 = vapply(paths, digest::digest, character(1L),
                    algo = "sha256", file = TRUE)
  )

  within_runtime <- make_sources("within-runtime", 1L)
  writeLines(
    paste(digest::digest(within_runtime, algo = "sha256", file = TRUE),
          within_runtime),
    file.path(td, "environment", "within_s5", "runtime_asset_sha256.txt")
  )
  transfer_runtime <- make_sources("transfer-runtime", 3L)
  data.table::fwrite(data.table::data.table(
    path = transfer_runtime, bytes = unname(file.info(transfer_runtime)$size),
    sha256 = vapply(transfer_runtime, digest::digest, character(1L),
                    algo = "sha256", file = TRUE)
  ), file.path(td, "environment", "transfer_argos",
               "public_model_manifest.tsv"), sep = "\t")

  within_shards <- make_sources("within-shards", 2210L)
  for (i in seq_along(splitter_methods)) {
    index <- ((i - 1L) * 170L + 1L):(i * 170L)
    shard_rows <- source_rows("within_shard", within_shards[index])
    data.table::fwrite(
      shard_rows[, .(bytes, sha256)],
      file.path(td, "splitter_manifests", paste0(
        "output_manifest__", splitter_methods[[i]], ".tsv"
      )), sep = "\t"
    )
  }
  splitter_sources <- file.path(
    source_root, "splitter-manifests",
    paste0("output_manifest__", splitter_methods, ".tsv")
  )
  dir.create(dirname(splitter_sources[[1L]]), recursive = TRUE)
  file.copy(file.path(td, "splitter_manifests", basename(splitter_sources)),
            splitter_sources)
  splitter_sources <- normalizePath(splitter_sources)

  transfer_shards <- make_sources("transfer-shards", 35L)
  transfer_shard_rows <- source_rows("transfer_shard", transfer_shards)
  data.table::fwrite(
    transfer_shard_rows[, .(bytes, sha256)],
    file.path(td, "transfer_shard_manifest.tsv"), sep = "\t"
  )

  base_within <- make_sources("base-within", 1L)
  base_transfer <- make_sources("base-transfer", 1L)
  external_within_task <- file.path(source_root, "within_tasks.tsv")
  external_transfer_task <- file.path(source_root, "transfer_tasks.tsv")
  external_witness <- file.path(source_root, "structural_witness.rds")
  file.copy(file.path(td, "inputs", "within_tasks.tsv"), external_within_task)
  file.copy(file.path(td, "inputs", "transfer_tasks.tsv"), external_transfer_task)
  file.copy(file.path(td, "structural_witness.rds"), external_witness)

  copy_receipt_sources <- function(route) {
    candidate_paths <- file.path(
      td, env$.snapshot_candidate_environment_files()[
        startsWith(env$.snapshot_candidate_environment_files(),
                   paste0("environment/", route, "/"))
      ]
    )
    destination <- file.path(source_root, paste0(route, "-receipt"),
                             basename(candidate_paths))
    dir.create(dirname(destination[[1L]]), recursive = TRUE)
    file.copy(candidate_paths, destination)
    normalizePath(destination)
  }
  within_receipt_sources <- copy_receipt_sources("within_s5")
  transfer_receipt_sources <- copy_receipt_sources("transfer_argos")
  receipt_digest <- function(route) {
    relative <- env$.snapshot_candidate_environment_files()[
      startsWith(env$.snapshot_candidate_environment_files(),
                 paste0("environment/", route, "/"))
    ]
    paths <- file.path(td, relative)
    rows <- data.table::data.table(
      path = basename(relative), bytes = unname(file.info(paths)$size),
      sha256 = vapply(paths, digest::digest, character(1L),
                      algo = "sha256", file = TRUE)
    )
    data.table::setorder(rows, path)
    digest::digest(paste(rows$path, rows$bytes, rows$sha256,
                         sep = "\t", collapse = "\n"), algo = "sha256")
  }
  within_receipt_sha <- receipt_digest("within_s5")
  transfer_receipt_sha <- receipt_digest("transfer_argos")
  data.table::fwrite(data.table::data.table(
    execution_route = "exact", analysis_surface = c("within", "transfer"),
    environment_route = c("within_s5", "transfer_argos"),
    environment_receipt_path = c(
      "environment/within_s5", "environment/transfer_argos"
    ),
    environment_receipt_sha256 = c(within_receipt_sha, transfer_receipt_sha)
  ), file.path(td, "route_environment_map.tsv"), sep = "\t")
  source_manifest <- data.table::rbindlist(list(
    source_rows("base_transfer_bundle", base_transfer),
    source_rows("base_within_bundle", base_within),
    source_rows("structural_witness", external_witness),
    source_rows("transfer_environment_receipt", transfer_receipt_sources),
    source_rows("transfer_runtime_asset", transfer_runtime),
    source_rows("transfer_shard", transfer_shards),
    source_rows("transfer_task_table", external_transfer_task),
    source_rows("within_environment_receipt", within_receipt_sources),
    source_rows("within_runtime_asset", within_runtime),
    source_rows("within_shard", within_shards),
    source_rows("within_splitter_manifest", splitter_sources),
    source_rows("within_task_table", external_within_task)
  ))
  data.table::setorder(source_manifest, input_role, path)
  data.table::fwrite(source_manifest, file.path(td, "source_input_manifest.tsv"),
                     sep = "\t")

  paths <- file.path(td, files)
  role <- rep("candidate_provenance", length(files))
  role[files %in% c("within_bundle.rds", "within_cells.tsv",
                    "within_predictions.tsv.gz")] <- "integrated_within"
  role[files %in% c("transfer_bundle.rds", "transfer_cells.tsv",
                    "transfer_predictions.tsv.gz")] <- "integrated_transfer"
  role[startsWith(files, "splitter_manifests/")] <- "splitter_manifest"
  role[startsWith(files, "environment/")] <- "environment_receipt"
  role[files == "structural_witness.rds"] <- "structural_witness"
  role[startsWith(files, "inputs/")] <- "frozen_task_table"
  role[files == "transfer_shard_manifest.tsv"] <- "transfer_shard_manifest"
  role[files == "legacy_compatibility_overlay_input_manifest.tsv"] <-
    "legacy_compatibility_manifest"
  role[files == "source_input_manifest.tsv"] <- "source_input_manifest"
  role[files == "route_environment_map.tsv"] <- "route_environment_map"
  role[files == "acceptance_checks.tsv"] <- "acceptance_checks"
  package_root <- normalizePath(testthat::test_path("..", ".."))
  producer_script <- file.path(
    package_root, "inst", "scripts",
    "paper3_build_dependency_overlay_candidate.R"
  )
  if (!file.exists(producer_script)) {
    installed_producer <- system.file(
      "scripts", "paper3_build_dependency_overlay_candidate.R",
      package = "OmicSelector"
    )
    expect_true(nzchar(installed_producer) && file.exists(installed_producer))
    package_root <- file.path(source_root, "producer-package-root")
    dir.create(file.path(package_root, "inst", "scripts"), recursive = TRUE)
    producer_script <- file.path(
      package_root, "inst", "scripts",
      "paper3_build_dependency_overlay_candidate.R"
    )
    expect_true(file.copy(installed_producer, producer_script))
  }
  manifest <- data.table::data.table(
    file = files, role = role, bytes = file.info(paths)$size,
    sha256 = vapply(paths, digest::digest, character(1L),
                    algo = "sha256", file = TRUE),
    run_mode = "full",
    producer_script_sha256 = digest::digest(
      producer_script, algo = "sha256", file = TRUE
    ),
    producer_package_version = "2.6.5.9000",
    producer_package_commit = strrep("b", 40L),
    base_within_bundle_sha256 = digest::digest(
      base_within, algo = "sha256", file = TRUE
    ),
    base_transfer_bundle_sha256 = digest::digest(
      base_transfer, algo = "sha256", file = TRUE
    ),
    within_environment_receipt_sha256 = within_receipt_sha,
    transfer_environment_receipt_sha256 = transfer_receipt_sha,
    execution_package_version = "2.6.5",
    execution_package_commit = strrep("3", 40L),
    cache_package_commit = strrep("4", 40L),
    within_analysis_code_id = strrep("5", 16L),
    transfer_analysis_code_id = strrep("6", 16L),
    within_task_sha256 = digest::digest(
      external_within_task, algo = "sha256", file = TRUE
    ),
    transfer_task_sha256 = digest::digest(
      external_transfer_task, algo = "sha256", file = TRUE
    ),
    structural_witness_sha256 = digest::digest(
      external_witness, algo = "sha256", file = TRUE
    ),
    structural_cache_sha256 = strrep("a", 64L),
    structural_fold_engine_sha256 = strrep("c", 64L)
  )
  manifest_path <- file.path(td, "manifest.tsv")
  data.table::fwrite(manifest, manifest_path, sep = "\t")
  manifest_sha <- digest::digest(manifest_path, algo = "sha256", file = TRUE)
  expect_equal(
    env$.snapshot_validate_overlay_candidate(
      td, manifest_sha, "2.6.5.9000", strrep("b", 40L), package_root,
      require_full = TRUE
    )$within_cells,
    file.path(td, "within_cells.tsv")
  )

  writeLines("tampered external shard", within_shards[[1L]])
  expect_error(
    env$.snapshot_validate_overlay_candidate(
      td, manifest_sha, "2.6.5.9000", strrep("b", 40L), package_root,
      require_full = TRUE
    ),
    "source-input bytes"
  )
  writeLines("within-shards 1", within_shards[[1L]])

  checks <- data.table::fread(file.path(td, "acceptance_checks.tsv"))
  checks[, pass := FALSE]
  data.table::fwrite(checks, file.path(td, "acceptance_checks.tsv"), sep = "\t")
  expect_error(
    env$.snapshot_validate_overlay_candidate(
      td, manifest_sha, "2.6.5.9000", strrep("b", 40L), package_root,
      require_full = TRUE
    ),
    "bytes differ"
  )
  acceptance_path <- file.path(td, "acceptance_checks.tsv")
  manifest[file == "acceptance_checks.tsv", `:=`(
    bytes = file.info(acceptance_path)$size,
    sha256 = digest::digest(acceptance_path, algo = "sha256", file = TRUE)
  )]
  data.table::fwrite(manifest, manifest_path, sep = "\t")
  failed_manifest_sha <- digest::digest(
    manifest_path, algo = "sha256", file = TRUE
  )
  expect_error(
    env$.snapshot_validate_overlay_candidate(
      td, failed_manifest_sha, "2.6.5.9000", strrep("b", 40L),
      package_root, require_full = TRUE
    ),
    "failed or invalid acceptance"
  )
})

test_that("snapshot rejects duplicate overlay ownership", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  base_cells <- data.table::data.table(
    method = "m1", cohort = "u1", seed = 1L, eligible = FALSE,
    ineligible_reason = "base"
  )
  base_predictions <- data.table::data.table(
    method = character(), cohort = character(), seed = integer(),
    fold = integer(), sample_id = character(), group_id = character(),
    y = integer(), score_m = numeric()
  )
  td <- tempfile("duplicate-overlay-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  cells_path <- file.path(td, "cells.tsv")
  predictions_path <- file.path(td, "predictions.tsv")
  data.table::fwrite(data.table::data.table(
    method = "m1", cohort = "u1", seed = 1L, eligible = TRUE,
    ineligible_reason = NA_character_
  ), cells_path, sep = "\t")
  data.table::fwrite(data.table::data.table(
    method = "m1", cohort = "u1", seed = 1L, fold = 1L,
    sample_id = "s1", group_id = "g1", y = 0L, score_m = 0.2
  ), predictions_path, sep = "\t")
  expect_error(
    env$.snapshot_overlay(
      base_cells, base_predictions,
      c(cells_path, cells_path), c(predictions_path, predictions_path)
    ),
    "same registered cell key"
  )
})

test_that("canonical snapshot reference rejects groups crossing folds", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  reference <- data.table::data.table(
    fold = c(1L, 1L, 2L, 2L), sample_id = paste0("s", 1:4),
    group_id = c("g0", "g1", "g0", "g2"), y = c(0L, 1L, 0L, 1L)
  )
  expect_error(env$.snapshot_validate_reference(reference), "crosses folds")
  reference$group_id <- paste0("g", 1:4)
  expect_true(env$.snapshot_validate_reference(reference))
})

test_that("snapshot excludes partial-seed predictions without hiding failure", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  predictions <- data.table::data.table(
    method = c(rep("complete", 5L), rep("partial", 3L)),
    cohort = "u1", seed = c(101L, 202L, 303L, 404L, 505L,
                             101L, 202L, 303L),
    fold = 1L, sample_id = paste0("s", seq_len(8L))
  )
  eligibility <- data.table::data.table(
    method_id = c("complete", "partial"), unit_id = "u1",
    eligible = c(TRUE, FALSE), eligible_seeds = c(5L, 3L),
    reason = c("eligible", "seed 404 failed | seed 505 failed")
  )
  eligibility_before <- data.table::copy(eligibility)
  out <- env$.snapshot_filter_aggregate_predictions(predictions, eligibility)
  expect_equal(unique(out$method), "complete")
  expect_equal(sort(out$seed), c(101L, 202L, 303L, 404L, 505L))
  expect_identical(eligibility, eligibility_before)
  expect_false(eligibility[method_id == "partial", eligible])
  expect_match(eligibility[method_id == "partial", reason], "seed 404 failed")
})

test_that("snapshot output cannot contain or mutate immutable roots", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("snapshot-disjoint-")
  immutable <- file.path(td, "candidate")
  dir.create(immutable, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  expect_error(env$.snapshot_assert_disjoint_output(
    file.path(immutable, "output"), immutable
  ), "disjoint")
  expect_error(env$.snapshot_assert_disjoint_output(td, immutable), "disjoint")
  expect_equal(env$.snapshot_assert_disjoint_output(
    file.path(td, "sibling"), immutable
  ), file.path(td, "sibling"))
  project <- file.path(td, "paper-project")
  dir.create(file.path(project, "analysis"), recursive = TRUE)
  dir.create(file.path(project, "metadata"), recursive = TRUE)
  provenance_file <- file.path(project, "metadata", "provenance.tsv")
  writeLines("pinned", provenance_file)
  canonical_output <- file.path(project, "analysis", "fair-snapshot")
  expect_equal(env$.snapshot_assert_disjoint_output(
    canonical_output, provenance_file
  ), canonical_output)
  dangling <- file.path(td, "dangling-output")
  expect_true(file.symlink(file.path(td, "missing-target"), dangling))
  expect_true(env$.snapshot_path_entry_exists(dangling))
})

test_that("fair runner output cannot contain immutable roots or dangling links", {
  env <- paper3_script_env("paper3_fair_method_comparison.R")
  td <- tempfile("fair-disjoint-")
  immutable <- file.path(td, "snapshot")
  dir.create(immutable, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  expect_error(env$.fair_assert_disjoint_output(
    file.path(immutable, "result"), immutable
  ), "disjoint")
  dangling <- file.path(td, "dangling-output")
  expect_true(file.symlink(file.path(td, "missing-target"), dangling))
  expect_true(env$.fair_path_entry_exists(dangling))
})

test_that("provenance resolution is bound to accession-to-unit mapping", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("provenance-map-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  script <- file.path(td, "preflight.R")
  writeLines(c('cat("Overlap hits: 1\\n")', "quit(status = 2L)"), script)
  resolution <- file.path(td, "resolution.tsv")
  data.table::fwrite(data.table::data.table(
    accession_a = "a", accession_b = "b", unit_a = "u1", unit_b = "u1"
  ), resolution, sep = "\t")
  units <- data.table::data.table(
    unit_id = c("u1", "u2"), source_accessions = c("a", "b")
  )
  expect_error(
    env$.snapshot_validate_provenance(
      script, td, units, resolution, file.path(td, "preflight.log")
    ),
    "accession-to-unit map"
  )
})

test_that("provenance resolution matches exact fresh overlap identities", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("provenance-identities-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  units <- data.table::data.table(
    unit_id = c("u1", "u2"), source_accessions = c("a;b", "c;d")
  )
  resolution <- file.path(td, "resolution.tsv")
  data.table::fwrite(data.table::data.table(
    accession_a = "c", accession_b = "d", unit_a = "u2", unit_b = "u2"
  ), resolution, sep = "\t")
  script <- file.path(td, "preflight.R")
  writeLines(c(
    'cat("Overlap hits: 1\\n")', 'cat("  A <-> B\\n")',
    "quit(status = 2L)"
  ), script)
  expect_error(env$.snapshot_validate_provenance(
    script, td, units, resolution, file.path(td, "preflight.log")
  ), "identities")

  data.table::fwrite(data.table::data.table(
    accession_a = c("a", "c"), accession_b = c("b", "d"),
    unit_a = c("u1", "u2"), unit_b = c("u1", "u2")
  ), resolution, sep = "\t")
  writeLines(c(
    'cat("Overlap hits: 2\\n")', 'cat("  A <-> B\\n")',
    'cat("  A <-> B\\n")', "quit(status = 2L)"
  ), script)
  expect_error(env$.snapshot_validate_provenance(
    script, td, units, resolution, file.path(td, "preflight.log")
  ), "identities")
})

test_that("provenance preflight requires exactly one explicit overlap count", {
  env <- paper3_script_env("paper3_build_fair_snapshot.R")
  td <- tempfile("provenance-count-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  units <- data.table::data.table(
    unit_id = c("u1", "u2"), source_accessions = c("a", "b")
  )
  resolution <- file.path(td, "resolution.tsv")
  data.table::fwrite(data.table::data.table(
    accession_a = character(), accession_b = character(),
    unit_a = character(), unit_b = character()
  ), resolution, sep = "\t")
  script <- file.path(td, "preflight.R")
  writeLines('cat("preflight complete\\n")', script)
  expect_error(env$.snapshot_validate_provenance(
    script, td, units, resolution, file.path(td, "preflight.log")
  ), "exactly one")
  writeLines(c('cat("Overlap hits: 0\\n")',
               'cat("Overlap hits: 0\\n")'), script)
  expect_error(env$.snapshot_validate_provenance(
    script, td, units, resolution, file.path(td, "preflight.log")
  ), "exactly one")
})

test_that("fair paper pair helper requires an exact matched join", {
  env <- paper3_script_env("paper3_fair_method_comparison.R")
  strata <- rep(seq_len(5L), each = 8L)
  groups <- rep(paste0("g", seq_len(20L)), each = 2L)
  samples <- paste0("p", seq_len(40L))
  y <- rep(rep(c(0L, 0L, 1L, 1L), each = 2L), times = 5L)
  one <- data.table::data.table(
    unit_id = "u1", seed = 1L, fold = strata,
    split_id = "fixed", sample_id = samples, group_id = groups, y = y
  )
  a <- data.table::copy(one)[, `:=`(method_id = "m1", score = y + 0.01)]
  b <- data.table::copy(one)[, `:=`(method_id = "m2", score = 1 - y + 0.01)]
  predictions <- data.table::rbindlist(list(a, b))
  eligibility <- data.table::data.table(
    method_id = c("m1", "m2"), unit_id = "u1", eligible = TRUE,
    reason = "eligible"
  )
  splits <- data.table::data.table(
    unit_id = "u1", seed = 1L, fold = seq_len(5L), split_id = "fixed",
    n_test_groups = 4L, n_test_profiles = 8L
  )
  pair <- data.table::data.table(
    pair_id = "m1::m2", method_a = "m1", method_b = "m2"
  )

  data.table::setkey(predictions, method_id, unit_id)
  data.table::setkey(eligibility, method_id, unit_id)
  data.table::setkey(splits, unit_id, seed, fold)
  stratum_auc <- env$.fair_precompute_strata(predictions, splits)
  counts <- predictions[, .(n_profiles = .N), by = .(method_id, unit_id)]
  data.table::setkey(counts, method_id, unit_id)
  result <- env$.fair_pair_unit(
    pair, "u1", stratum_auc, eligibility, splits, counts
  )
  expect_true(result$join$pass)
  expect_true(result$summary$common_support)
  expect_equal(result$summary$estimate, 1)
  expect_equal(nrow(result$strata), 5L)

  matched <- data.table::dcast(
    predictions,
    unit_id + seed + fold + split_id + sample_id + group_id + y ~ method_id,
    value.var = "score"
  )
  exported <- singlesample_matched_pair_auc(
    y = matched$y, score_a = matched$m1, score_b = matched$m2,
    stratum_id = paste(matched$seed, matched$fold, sep = "::"),
    sample_id = matched$sample_id, group_id = matched$group_id,
    expected_strata = paste(1L, seq_len(5L), sep = "::"),
    analysis_level = "group"
  )
  expect_equal(result$summary$estimate, exported$summary$estimate,
               tolerance = 0)
  expect_equal(result$summary$se, exported$summary$se, tolerance = 0)
  expect_equal(result$strata$effect, exported$strata$effect, tolerance = 0)

  reversed_pair <- data.table::data.table(
    pair_id = "m2::m1", method_a = "m2", method_b = "m1"
  )
  reversed <- env$.fair_pair_unit(
    reversed_pair, "u1", stratum_auc, eligibility, splits, counts
  )
  expect_equal(reversed$summary$estimate, -result$summary$estimate,
               tolerance = 0)
  zero_strata <- data.table::copy(stratum_auc)
  zero_strata[method_id == "m2", auc :=
                zero_strata[method_id == "m1", auc]]
  zero <- env$.fair_pair_unit(
    pair, "u1", zero_strata, eligibility, splits, counts
  )
  expect_equal(zero$summary$estimate, 0, tolerance = 0)
  expect_equal(zero$summary$se, 0, tolerance = 0)
  expect_true(is.na(zero$summary$p_zero))

  broken <- data.table::copy(counts)
  broken[list("m2", "u1"), n_profiles := n_profiles - 1L]
  expect_error(
    env$.fair_pair_unit(
      pair, "u1", stratum_auc, eligibility, splits, broken
    ),
    "exact join failed"
  )
})

test_that("hindsight accepts base-data-frame panel metadata", {
  env <- paper3_script_env("paper3_fair_method_comparison.R")
  methods <- paste0("m", 1:5)
  groups <- paste0("g", 1:8)
  labels <- rep(c(0L, 0L, 1L, 1L), 2L)
  predictions <- data.table::CJ(
    method_id = methods, stratum = 1:2, group_index = 1:4
  )
  predictions[, `:=`(
    unit_id = "u1", seed = 1L, fold = stratum,
    group_id = groups[(stratum - 1L) * 4L + group_index],
    sample_id = paste0("s", stratum, "_", group_index),
    y = labels[(stratum - 1L) * 4L + group_index],
    score = labels[(stratum - 1L) * 4L + group_index] +
      as.integer(sub("m", "", method_id)) / 100
  )]
  performance <- data.table::data.table(
    method_id = methods, unit_id = "u1", auc = 1
  )
  panels <- list(
    summary = data.frame(threshold = 0.95, degenerate = FALSE),
    method_membership = data.frame(
      threshold = 0.95, method_id = methods, included_method = TRUE
    ),
    unit_membership = data.frame(
      threshold = 0.95, unit_id = "u1", included_unit = TRUE
    )
  )
  out <- env$.fair_hindsight(
    predictions, performance, panels, data.table::data.table(unit_id = "u1"),
    reps = 100L, seed = 19L
  )
  expect_equal(out$status, "evaluable")
  expect_equal(out$valid_permutations, 100L)
})

test_that("rank bootstrap is deterministic on a complete rectangle", {
  env <- paper3_script_env("paper3_fair_method_comparison.R")
  methods <- paste0("m", 1:5)
  units <- paste0("u", 1:5)
  performance <- data.table::CJ(method_id = methods, unit_id = units)
  performance[, auc := 0.5 + as.integer(sub("m", "", method_id)) / 100 +
                as.integer(sub("u", "", unit_id)) / 1000]
  panels <- list(
    summary = data.frame(threshold = 0.95, degenerate = FALSE),
    method_membership = data.frame(
      threshold = 0.95, method_id = methods, included_method = TRUE
    ),
    unit_membership = data.frame(
      threshold = 0.95, unit_id = units, included_unit = TRUE
    )
  )
  a <- env$.fair_rank_bootstrap(performance, panels, 100L, 17L)
  b <- env$.fair_rank_bootstrap(performance, panels, 100L, 17L)
  expect_identical(a, b)
  expect_equal(a[method_id == "m5", p_top1], 1)
  expect_equal(a[method_id == "m1", p_top5], 1)
  expect_equal(a$rank_median, a$observed_average_rank)
})

test_that("rank bootstrap treats tied top methods as inclusively top one", {
  env <- paper3_script_env("paper3_fair_method_comparison.R")
  methods <- paste0("m", 1:5)
  units <- paste0("u", 1:5)
  performance <- data.table::CJ(method_id = methods, unit_id = units)
  performance[, auc := 0.7]
  panels <- list(
    summary = data.frame(threshold = 0.95, degenerate = FALSE),
    method_membership = data.frame(
      threshold = 0.95, method_id = methods, included_method = TRUE
    ),
    unit_membership = data.frame(
      threshold = 0.95, unit_id = units, included_unit = TRUE
    )
  )
  out <- env$.fair_rank_bootstrap(performance, panels, 100L, 17L)
  expect_true(all(out$p_top1 == 1))
  expect_true(all(out$observed_average_rank == 3))
  expect_true(all(out$tie_policy ==
    "inclusive_min_rank_for_topk;average_rank_for_intervals"))
})

test_that("full fair-run arguments cannot weaken registered resampling", {
  env <- paper3_script_env("paper3_fair_method_comparison.R")
  runtime_args <- c(
    "--runtime-receipt=/tmp/runtime-receipt",
    paste0("--expected-runtime-receipt-manifest-sha256=", strrep("d", 64L)),
    "--runtime-image=/tmp/runtime.sif",
    paste0("--expected-runtime-image-sha256=", strrep("e", 64L))
  )
  ok <- env$.fair_args(c(c(
    "--input-dir=/tmp/in", "--output-dir=/tmp/out",
    paste0("--expected-snapshot-manifest-sha256=", strrep("a", 64L)),
    "--expected-package-version=2.6.5.9000",
    paste0("--expected-package-commit=", strrep("b", 40L)),
    paste0("--expected-installed-package-tree-sha256=", strrep("c", 64L)),
    "--bootstrap-reps=10000", "--permutation-reps=1000",
    "--seed=20260718", "--run-mode=full"
  ), runtime_args))
  expect_identical(ok$bootstrap_reps, 10000L)
  expect_error(env$.fair_args(c(c(
    "--input-dir=/tmp/in", "--output-dir=/tmp/out",
    paste0("--expected-snapshot-manifest-sha256=", strrep("a", 64L)),
    "--expected-package-version=2.6.5.9000",
    paste0("--expected-package-commit=", strrep("b", 40L)),
    paste0("--expected-installed-package-tree-sha256=", strrep("c", 64L)),
    "--bootstrap-reps=100", "--permutation-reps=100",
    "--seed=1", "--run-mode=full"
  ), runtime_args)), "Full mode requires")
  expect_error(env$.fair_args(c(c(
    "--input-dir=/tmp/in", "--input-dir=/tmp/other",
    "--output-dir=/tmp/out",
    paste0("--expected-snapshot-manifest-sha256=", strrep("a", 64L)),
    "--expected-package-version=2.6.5.9000",
    paste0("--expected-package-commit=", strrep("b", 40L)),
    paste0("--expected-installed-package-tree-sha256=", strrep("c", 64L))
  ), runtime_args)), "Duplicate argument")
  expect_error(env$.fair_args(c(c(
    "--input-dir=/tmp/in", "--output-dir=/tmp/out",
    paste0("--expected-snapshot-manifest-sha256=", strrep("a", 64L)),
    "--expected-package-version=2.6.5.9000",
    paste0("--expected-package-commit=", strrep("b", 40L)),
    paste0("--expected-installed-package-tree-sha256=", strrep("c", 64L)),
    "--permutations=1000"
  ), runtime_args)), "Unknown argument")
})

test_that("nested-selector task arguments freeze the registered bootstrap", {
  env <- paper3_script_env("paper3_nested_selector_unit.R")
  args <- c(
    "--cache=/tmp/cache.rds",
    paste0("--cache-sha256=", paste(rep("a", 64L), collapse = "")),
    "--unit=u1", "--seed=101", "--output-dir=/tmp/out",
    paste0("--package-commit=", paste(rep("b", 40L), collapse = "")),
    paste0("--runtime-image-sha256=", paste(rep("c", 64L), collapse = "")),
    paste0("--r-library-snapshot-sha256=", paste(rep("d", 64L), collapse = "")),
    "--inner-folds=5", "--bootstrap-reps=1000", "--verify-base=true"
  )
  parsed <- env$.selector_task_args(args)
  expect_identical(parsed$seed, 101L)
  expect_identical(parsed$bootstrap_reps, 1000L)
  expect_true(parsed$verify_base)
  expect_error(
    env$.selector_task_args(sub("1000", "999", args, fixed = TRUE)),
    "Registered tasks require"
  )
  expect_error(
    env$.selector_task_args(sub("verify-base=true", "verify-base=false",
                                args, fixed = TRUE)),
    "Registered tasks require"
  )
})

test_that("nested selector task producer freezes the exact 34 by 5 grid", {
  env <- paper3_script_env("paper3_build_nested_selector_tasks.R")
  expected <- env$.nested_expected_units()
  expect_length(expected, 34L)
  cache <- setNames(lapply(seq_along(expected), function(i) {
    n <- 20L
    y <- rep(0:1, each = n / 2L)
    X <- matrix(seq_len(n * 8L), nrow = n,
                dimnames = list(paste0("s", i, "_", seq_len(n)),
                                paste0("f", seq_len(8L))))
    list(
      expr_per_sample = X, y_bin = y,
      group_id = paste0("g", i, "_", seq_len(n)),
      outer_k = if (expected[[i]] %in% env$.nested_expected_k3()) 3L else 5L,
      provenance_block = paste0("block::", expected[[i]]),
      source_accessions = expected[[i]], modality = "microarray",
      biospecimen = if (expected[[i]] == "GSE31568") "whole blood" else "serum"
    )
  }), expected)
  units <- env$.nested_validate_cache(cache)
  tasks <- env$.nested_task_grid(units)
  expect_equal(nrow(units), 34L)
  expect_equal(nrow(tasks), 170L)
  expect_identical(anyDuplicated(tasks[, c("unit_id", "seed")]), 0L)
  expect_setequal(unique(tasks$seed), c(101L, 202L, 303L, 404L, 505L))
  expect_equal(sum(units$outer_k == 3L), 8L)
  expect_equal(sum(units$primary_unit), 33L)
  expect_identical(units[primary_unit == FALSE, unit_id], "GSE31568")
  expect_true(all(nzchar(units$modality)))
  expect_true(all(nzchar(units$biospecimen)))
})

test_that("nested-selector synthesis arguments have an explicit run mode", {
  env <- paper3_script_env("paper3_summarize_nested_selector.R")
  args <- c(
    "--task-root=/tmp/tasks", "--units=/tmp/units.tsv",
    "--output-dir=/tmp/out", "--run-mode=full"
  )
  expect_identical(env$.selector_summary_args(args)$run_mode, "full")
  expect_error(
    env$.selector_summary_args(sub("full", "partial", args, fixed = TRUE)),
    "full or smoke"
  )
})
