.fairqa_script_env <- function() {
  env <- new.env(parent = globalenv())
  path <- testthat::test_path("..", "..", "inst", "scripts",
                             "paper3_validate_fair_method_comparison.R")
  if (!file.exists(path)) path <- system.file(
    "scripts", "paper3_validate_fair_method_comparison.R",
    package = "OmicSelector"
  )
  env$.PAPER3_ACTIVE_SCRIPT_PATH <- normalizePath(path, mustWork = TRUE)
  sys.source(path, envir = env)
  env
}

.fairqa_runner_env <- function() {
  env <- new.env(parent = globalenv())
  path <- testthat::test_path("..", "..", "inst", "scripts",
                             "paper3_fair_method_comparison.R")
  if (!file.exists(path)) path <- system.file(
    "scripts", "paper3_fair_method_comparison.R", package = "OmicSelector"
  )
  env$.PAPER3_ACTIVE_SCRIPT_PATH <- normalizePath(path, mustWork = TRUE)
  sys.source(path, envir = env)
  env
}

.fairqa_args_fixture <- function() c(
  "--input-dir=/tmp/in", "--result-dir=/tmp/result",
  "--output-dir=/tmp/qa",
  paste0("--expected-snapshot-manifest-sha256=", strrep("a", 64L)),
  "--expected-package-version=2.6.5.9000",
  paste0("--expected-package-commit=", strrep("b", 40L)),
  paste0("--expected-installed-package-tree-sha256=", strrep("c", 64L)),
  "--runtime-receipt=/tmp/runtime-receipt",
  paste0("--expected-runtime-receipt-manifest-sha256=", strrep("d", 64L)),
  "--runtime-image=/tmp/runtime.simg",
  paste0("--expected-runtime-image-sha256=", strrep("e", 64L)),
  "--run-mode=full"
)

test_that("fair QA arguments require every snapshot/package/runtime pin", {
  env <- .fairqa_script_env(); args <- .fairqa_args_fixture()
  expect_identical(env$.fairqa_args(args)$run_mode, "full")
  expect_error(env$.fairqa_args(args[!grepl("runtime-receipt=", args)]),
               "Missing required")
  expect_error(env$.fairqa_args(c(args, "--runtime-image=/duplicate")),
               "Duplicate argument")
  expect_error(env$.fairqa_args(c(args, "--unexpected=true")),
               "Unknown argument")
  expect_error(env$.fairqa_args(sub(strrep("d", 64L), "bad", args,
                                    fixed = TRUE)), "not an exact SHA-256")
})

test_that("result manifest binds the exact current 18-artifact surface", {
  env <- .fairqa_script_env()
  td <- tempfile("fairqa-result-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  required <- env$.fairqa_result_artifacts()
  expect_equal(length(required), 18L)
  expect_setequal(required,
                  .fairqa_runner_env()$.fair_expected_result_files())
  expect_true(all(c("pairwise_profile_sensitivity.tsv",
                    "method_unit_profile_sensitivity.tsv") %in% required))
  for (file in required) writeLines(file, file.path(td, file))
  root <- normalizePath(testthat::test_path("..", ".."))
  runner <- file.path(root, "inst", "scripts",
                      "paper3_fair_method_comparison.R")
  manifest <- data.table::data.table(file = required)
  paths <- file.path(td, required)
  manifest[, `:=`(
    bytes = file.info(paths)$size,
    sha256 = vapply(paths, env$.fairqa_sha, character(1L)),
    run_mode = "full", bootstrap_reps = 10000L,
    permutation_reps = 1000L, seed = 20260718L,
    omicselector_version = "2.6.5.9000",
    package_git_commit = strrep("b", 40L), package_git_dirty = FALSE,
    producer_script_sha256 = env$.fairqa_sha(runner),
    installed_package_tree_sha256 = strrep("c", 64L),
    snapshot_manifest_sha256 = strrep("a", 64L),
    runtime_receipt_manifest_sha256 = strrep("d", 64L),
    runtime_image_sha256 = strrep("e", 64L)
  )]
  data.table::fwrite(manifest, file.path(td, "output_manifest.tsv"), sep = "\t")
  out <- env$.fairqa_validate_bundle_manifest(
    td, strrep("a", 64L), "2.6.5.9000", strrep("b", 40L),
    strrep("c", 64L), strrep("d", 64L), strrep("e", 64L),
    "full", root
  )
  expect_equal(nrow(out$manifest), 18L)
  writeLines("tampered profile sensitivity",
             file.path(td, "pairwise_profile_sensitivity.tsv"))
  expect_error(env$.fairqa_validate_bundle_manifest(
    td, strrep("a", 64L), "2.6.5.9000", strrep("b", 40L),
    strrep("c", 64L), strrep("d", 64L), strrep("e", 64L),
    "full", root), "artifact bytes differ")
})

test_that("outer-fold and path-isolation contracts are exact", {
  env <- .fairqa_script_env()
  units <- data.table::data.table(
    unit_id = sprintf("u%02d", seq_len(34L)),
    outer_k = rep(c(3L, 5L), length.out = 34L)
  )
  splits <- units[, data.table::CJ(
    seed = c(101L, 202L, 303L, 404L, 505L), fold = seq_len(outer_k)
  ), by = unit_id]
  expect_invisible(env$.fairqa_validate_fold_contract(units, splits))
  expect_error(env$.fairqa_validate_fold_contract(units, splits[-1L]),
               "exact 34x5")
  td <- tempfile("fairqa-path-"); dir.create(td)
  source_root <- file.path(td, "source"); dir.create(source_root)
  outside <- file.path(td, "outside")
  expect_invisible(env$.fairqa_assert_disjoint_output(outside, source_root))
  expect_error(env$.fairqa_assert_disjoint_output(
    file.path(source_root, "receipt"), source_root), "not disjoint")
  unlink(td, recursive = TRUE)
})

test_that("result and validation-receipt attestations are byte-bound", {
  env <- .fairqa_script_env()
  td <- tempfile("fairqa-receipt-"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  result_manifest <- file.path(td, "result_manifest.tsv")
  writeLines("f\tsha", result_manifest)
  result_sha <- env$.fairqa_sha(result_manifest)
  expect_invisible(env$.fairqa_assert_result_manifest_unchanged(
    result_manifest, result_sha
  ))
  writeLines(c("f\tsha", "coherent\treplacement"), result_manifest)
  expect_error(env$.fairqa_assert_result_manifest_unchanged(
    result_manifest, result_sha
  ), "changed after scientific reconstruction")

  stage <- file.path(td, "stage"); dir.create(stage)
  checks <- data.table::data.table(
    check_id = env$.fairqa_check_ids(), pass = TRUE
  )
  checks_path <- file.path(stage, "validation_checks.tsv")
  data.table::fwrite(checks, checks_path, sep = "\t")
  pins <- list(
    run_mode = "full", omicselector_version = "2.6.5.9000",
    package_git_commit = strrep("b", 40L),
    installed_package_tree_sha256 = strrep("c", 64L),
    producer_script_sha256 = strrep("f", 64L),
    snapshot_manifest_sha256 = strrep("a", 64L),
    result_manifest_sha256 = strrep("1", 64L),
    runtime_receipt_manifest_sha256 = strrep("d", 64L),
    runtime_image_sha256 = strrep("e", 64L)
  )
  manifest <- data.table::data.table(
    file = "validation_checks.tsv", bytes = file.info(checks_path)$size,
    sha256 = env$.fairqa_sha(checks_path), run_mode = pins$run_mode,
    omicselector_version = pins$omicselector_version,
    package_git_commit = pins$package_git_commit,
    installed_package_tree_sha256 = pins$installed_package_tree_sha256,
    producer_script_sha256 = pins$producer_script_sha256,
    snapshot_manifest_sha256 = pins$snapshot_manifest_sha256,
    result_manifest_sha256 = pins$result_manifest_sha256,
    runtime_receipt_manifest_sha256 = pins$runtime_receipt_manifest_sha256,
    runtime_image_sha256 = pins$runtime_image_sha256
  )
  manifest_path <- file.path(stage, "validation_manifest.tsv")
  data.table::fwrite(manifest, manifest_path, sep = "\t")
  expect_invisible(env$.fairqa_validate_staged_receipt(stage, pins))
  manifest$package_git_commit <- strrep("0", 40L)
  data.table::fwrite(manifest, manifest_path, sep = "\t")
  expect_error(env$.fairqa_validate_staged_receipt(stage, pins),
               "pin differs: package_git_commit")
})

.fairqa_small_pair_fixture <- function(env) {
  methods <- c("m1", "m2"); unit <- "u1"
  splits <- data.table::CJ(seed = c(101L, 202L, 303L, 404L, 505L),
                           fold = 1:2)
  splits[, `:=`(unit_id = unit, split_id = paste0("split-", seed),
                n_test_profiles = 4L, n_test_groups = 3L)]
  data.table::setkey(splits, unit_id, seed, fold)
  predictions <- data.table::CJ(method_id = methods,
                                seed = unique(splits$seed), fold = 1:2)
  predictions <- predictions[, data.table::data.table(
    unit_id = unit, split_id = paste0("split-", seed),
    sample_id = paste0("s", fold, "-", 1:4),
    group_id = paste0("f", fold, "-", c("g0", "g0", "g1", "g2")),
    y = c(0L, 0L, 1L, 1L),
    score = if (method_id == "m1") c(.1, .7, .5, .5) else
      c(.2, .4, .45, .55)
  ), by = .(method_id, seed, fold)]
  eligibility <- data.table::CJ(method_id = methods, unit_id = unit)
  eligibility[, `:=`(eligible = TRUE, reason = "eligible")]
  data.table::setkey(eligibility, method_id, unit_id)
  data.table::setkey(predictions, method_id, unit_id)
  counts <- predictions[, .(n_profiles = .N), by = .(method_id, unit_id)]
  data.table::setkey(counts, method_id, unit_id)
  profile <- env$.fairqa_precompute_strata(predictions, splits, "profile")
  group <- env$.fairqa_precompute_strata(predictions, splits, "group")
  pairs <- data.table::data.table(pair_id = "m1::m2", method_a = "m1",
                                  method_b = "m2")
  list(methods = methods, unit = unit, splits = splits,
       predictions = predictions, eligibility = eligibility, counts = counts,
       profile = profile, group = group, pairs = pairs)
}

test_that("group-primary AUC, corrected pairs, and profile sensitivity reconstruct", {
  env <- .fairqa_script_env(); x <- .fairqa_small_pair_fixture(env)
  expect_invisible(env$.fairqa_validate_prediction_splits(
    x$predictions, x$splits
  ))
  wrong_split <- data.table::copy(x$predictions)
  wrong_split[method_id == "m2" & seed == 101L,
              split_id := "method-specific-split"]
  expect_error(env$.fairqa_validate_prediction_splits(wrong_split, x$splits),
               "registered split contract")
  split_group <- data.table::copy(x$predictions)
  split_group[fold == 2L & grepl("g0$", group_id), group_id := "f1-g0"]
  expect_error(env$.fairqa_validate_prediction_splits(split_group, x$splits),
               "more than one held-out fold")
  performance <- env$.fairqa_performance(x$group, x$eligibility, x$splits)
  expect_true(all(performance$auc == 1))
  profile_performance <- env$.fairqa_performance(
    x$profile, x$eligibility, x$splits
  )
  expect_equal(profile_performance[method_id == "m1", auc], .5)
  expect_equal(profile_performance[method_id == "m2", auc], 1)
  method_sensitivity <- merge(performance[, .(
    method_id, unit_id, eligible, group_collapsed_auc = auc, n_strata
  )], profile_performance[, .(
    method_id, unit_id, profile_auc = auc
  )], by = c("method_id", "unit_id"), sort = FALSE)
  method_sensitivity[, `:=`(
    difference_profile_minus_group = profile_auc - group_collapsed_auc,
    repeated_profiles = TRUE
  )]
  data.table::setcolorder(method_sensitivity, c(
    "method_id", "unit_id", "eligible", "group_collapsed_auc", "n_strata",
    "profile_auc", "difference_profile_minus_group", "repeated_profiles"
  ))
  expect_named(method_sensitivity, c(
    "method_id", "unit_id", "eligible", "group_collapsed_auc", "n_strata",
    "profile_auc", "difference_profile_minus_group", "repeated_profiles"
  ))
  expect_equal(method_sensitivity[method_id == "m1",
                                  difference_profile_minus_group], -.5)
  altered_performance <- data.table::copy(performance)
  altered_performance$auc[[1L]] <- .5
  expect_error(env$.fairqa_compare(altered_performance, performance,
    c("method_id", "unit_id"), "Method-unit performance"), "auc")
  pair <- env$.fairqa_pair_family(x$pairs, x$unit, x$group,
                                  x$eligibility, x$splits, x$counts)
  profile <- env$.fairqa_pair_family(x$pairs, x$unit, x$profile,
                                     x$eligibility, x$splits, x$counts)
  expect_equal(nrow(pair$strata), 10L)
  expect_true(pair$pairwise$common_support)
  expect_equal(pair$pairwise$estimate, 0)
  expect_equal(pair$pairwise$expected_m, 10L)
  expect_equal(pair$pairwise$reason, "zero_corrected_se")
  sensitivity <- merge(pair$pairwise[, .(
    pair_id, method_a, method_b, unit_id, common_support,
    group_collapsed_estimate = estimate, group_collapsed_se = se,
    group_collapsed_ci_low = ci_low, group_collapsed_ci_high = ci_high)],
    profile$pairwise[, .(pair_id, unit_id,
      profile_common_support = common_support,
      profile_estimate = estimate, profile_se = se,
      profile_ci_low = ci_low, profile_ci_high = ci_high)],
    by = c("pair_id", "unit_id"))
  sensitivity[, difference_profile_minus_group :=
                profile_estimate - group_collapsed_estimate]
  data.table::setcolorder(sensitivity, c(
    "pair_id", "method_a", "method_b", "unit_id", "common_support",
    "group_collapsed_estimate", "group_collapsed_se",
    "group_collapsed_ci_low", "group_collapsed_ci_high",
    "profile_common_support", "profile_estimate", "profile_se",
    "profile_ci_low", "profile_ci_high", "difference_profile_minus_group"
  ))
  expect_named(sensitivity, c(
    "pair_id", "method_a", "method_b", "unit_id", "common_support",
    "group_collapsed_estimate", "group_collapsed_se",
    "group_collapsed_ci_low", "group_collapsed_ci_high",
    "profile_common_support", "profile_estimate", "profile_se",
    "profile_ci_low", "profile_ci_high", "difference_profile_minus_group"
  ))
  expect_equal(sensitivity$difference_profile_minus_group, -.5)
  altered <- data.table::copy(sensitivity)
  altered$difference_profile_minus_group <- .1
  expect_error(env$.fairqa_compare(altered, sensitivity,
    c("pair_id", "unit_id"), "Pairwise profile sensitivity"),
    "difference_profile_minus_group")
  altered_strata <- data.table::copy(pair$strata); altered_strata$auc_a[[1L]] <- .5
  expect_error(env$.fairqa_compare(altered_strata, pair$strata,
    c("pair_id", "unit_id", "stratum_id"), "Pairwise strata"), "auc_a")
})

test_that("meta, LOCO, and heterogeneity policies are independently checked", {
  skip_if_not_installed("metafor")
  env <- .fairqa_script_env()
  pairs <- data.table::data.table(pair_id = "a::b", method_a = "a", method_b = "b")
  pairwise <- data.table::data.table(
    pair_id = "a::b", method_a = "a", method_b = "b",
    unit_id = paste0("u", 1:6), common_support = TRUE,
    estimate = c(-.20, -.10, .00, .10, .20, .30), se = rep(.02, 6L)
  )
  frozen <- data.table::data.table(pair_id = "a::b", k_support = 6L,
    support_units = paste0("u", 1:6, collapse = ";"), pooled_testable = TRUE)
  meta <- env$.fairqa_meta_all(pairwise, pairs, frozen)
  expect_identical(meta$status, "REML_KH")
  expect_true(all(is.finite(as.matrix(
    meta[, .(estimate, se, tau2, i2, p_zero)]
  ))))
  altered_meta <- data.table::copy(meta); altered_meta$tau2 <- meta$tau2 + .1
  expect_error(env$.fairqa_compare(altered_meta, meta, "pair_id", "Meta"),
               "tau2")
  loco <- env$.fairqa_loco(pairwise, pairs)
  expect_setequal(loco$omitted_unit, pairwise$unit_id)
  expect_true(all(loco$status == "REML_KH"))
  altered_loco <- data.table::copy(loco); altered_loco$status[[1L]] <- "fit_failed"
  expect_error(env$.fairqa_compare(altered_loco, loco,
    c("pair_id", "omitted_unit"), "LOCO"), "status")
  units <- data.table::data.table(unit_id = pairwise$unit_id,
    modality = rep(c("seq", "array"), each = 3L),
    biospecimen = rep(c("serum", "plasma"), 3L))
  modality <- env$.fairqa_meta_strata(pairwise, pairs, units, "modality", TRUE)
  biospecimen <- env$.fairqa_meta_strata(pairwise, pairs, units,
                                         "biospecimen", FALSE)
  expect_true(all(modality$inference_policy == "REML_KH_when_k_at_least_5"))
  expect_true(all(biospecimen$status == "descriptive_only"))
  expect_true(all(is.na(biospecimen$estimate)))
  altered <- data.table::copy(biospecimen); altered$inference_policy[[1L]] <- "REML"
  expect_error(env$.fairqa_compare(altered, biospecimen,
    c("stratum_variable", "stratum_value", "pair_id"), "Heterogeneity"),
    "inference_policy")
})

.fairqa_panel_fixture <- function(env) {
  ids <- paste0("m", 1:5); units <- paste0("u", 1:5)
  methods <- data.table::data.table(
    method_id = ids, roster_order = 1:5, family = "f", role = "primary",
    tier = "core", dep_route = "base", pkg_status = "present"
  )
  eligibility <- data.table::CJ(method_id = ids, unit_id = units)
  eligibility[, `:=`(eligible = TRUE, reason = "eligible")]
  performance <- data.table::copy(eligibility)
  performance[, `:=`(auc = c(.9, .9, .8, .7, .6)[match(method_id, ids)],
                      n_strata = 1L)]
  splits <- data.table::data.table(unit_id = units, seed = 101L, fold = 1L,
    n_test_profiles = 4L, n_test_groups = 4L)
  predictions <- data.table::CJ(method_id = ids, unit_id = units)
  predictions <- predictions[, data.table::data.table(
    seed = 101L, fold = 1L, split_id = paste0("split-", unit_id),
    sample_id = paste0("s", 1:4), group_id = paste0("g", 1:4),
    y = c(0L, 0L, 1L, 1L),
    score = c(.1, .2, .8, .9) - .01 * match(method_id, ids)
  ), by = .(method_id, unit_id)]
  panels <- env$.fairqa_panels(eligibility, methods, units)
  list(ids = ids, units = units, methods = methods, eligibility = eligibility,
       performance = performance, splits = splits, predictions = predictions,
       panels = panels)
}

test_that("panels, coverage, inclusive ranks, and hindsight reject tampering", {
  env <- .fairqa_script_env(); x <- .fairqa_panel_fixture(env)
  expect_true(all(x$panels[item_type == "summary", !degenerate]))
  altered_panels <- data.table::copy(x$panels)
  altered_panels[which(item_type == "method")[[1L]], included_method := FALSE]
  expect_error(env$.fairqa_compare(altered_panels, x$panels,
    c("threshold", "item_type", "method_id", "unit_id"), "Panels"),
    "included_method")
  rectangles <- env$.fairqa_rectangles(x$performance, x$panels)
  expect_equal(nrow(rectangles), 3L * 25L)
  altered_rectangle <- data.table::copy(rectangles)
  altered_rectangle$within_unit_rank[[1L]] <- 99
  expect_error(env$.fairqa_compare(altered_rectangle, rectangles,
    c("threshold", "method_id", "unit_id"), "Rectangles"),
    "within_unit_rank")
  coverage <- env$.fairqa_coverage(x$methods, x$eligibility, x$performance,
    x$predictions, x$splits, x$units)
  expect_true(all(coverage$prediction_complete))
  expect_true(all(coverage$deployment_status == "package_native"))
  expect_true(any(coverage$pareto_optimal))
  altered_coverage <- data.table::copy(coverage)
  altered_coverage$prediction_complete[[1L]] <- FALSE
  expect_error(env$.fairqa_compare(altered_coverage, coverage, "method_id",
                                   "Coverage"), "prediction_complete")
  ranks <- env$.fairqa_rank_bootstrap(x$performance, x$panels, 100L, 20260718L)
  expect_true(all(ranks$bootstrap_reps == 100L))
  expect_true(all(ranks$tie_policy ==
    "inclusive_min_rank_for_topk;average_rank_for_intervals"))
  expect_true(all(ranks[method_id %in% c("m1", "m2"), p_top1] == 1))
  expect_true(all(is.finite(ranks$observed_average_rank)))
  rank_metric <- data.table::copy(x$performance)
  rank_metric[, auc := data.table::fcase(
    method_id == "m1" & unit_id == "u5", .20,
    method_id == "m1", .61,
    method_id == "m2" & unit_id == "u5", 1.00,
    method_id == "m2", .60,
    method_id == "m3", .10,
    method_id == "m4", .05,
    method_id == "m5", .00
  )]
  metric_ranks_all <- env$.fairqa_rank_bootstrap(
    rank_metric, x$panels, 100L, 20260718L
  )
  runner_panels <- list(
    summary = as.data.frame(x$panels[item_type == "summary"]),
    method_membership = as.data.frame(x$panels[item_type == "method"]),
    unit_membership = as.data.frame(x$panels[item_type == "unit"])
  )
  runner_ranks <- .fairqa_runner_env()$.fair_rank_bootstrap(
    rank_metric, runner_panels, 100L, 20260718L
  )
  expect_invisible(env$.fairqa_compare(
    runner_ranks, metric_ranks_all, c("threshold", "method_id"),
    "Runner/validator rank contract"
  ))
  metric_ranks <- metric_ranks_all[threshold == .95]
  expect_equal(metric_ranks[method_id == "m1", observed_average_rank], 1.2)
  expect_equal(metric_ranks[method_id == "m2", observed_average_rank], 1.8)
  expect_equal(metric_ranks[method_id == "m1", rank_median], 1.2)
  expect_equal(metric_ranks[method_id == "m2", rank_median], 1.8)
  altered_rank <- data.table::copy(ranks); altered_rank$p_top1[[1L]] <- 0
  expect_error(env$.fairqa_compare(altered_rank, ranks,
    c("threshold", "method_id"), "Rank"), "p_top1")
  hindsight <- env$.fairqa_hindsight(x$predictions, x$performance, x$panels,
                                     100L, 20260718L)
  expect_equal(nrow(hindsight), 5L)
  expect_true(all(hindsight$permutation_reps == 100L))
  expect_true(all(hindsight$seed == 20260718L))
  expect_true(all(hindsight$analysis_level == "group"))
  altered_hindsight <- data.table::copy(hindsight)
  altered_hindsight$observed_max[[1L]] <- 0
  expect_error(env$.fairqa_compare(altered_hindsight, hindsight, "unit_id",
                                   "Hindsight"), "observed_max")
})

test_that("skeletal all-ineligible science is rejected before inference", {
  env <- .fairqa_script_env()
  td <- tempfile("fairqa-skeletal-"); input <- file.path(td, "input")
  result <- file.path(td, "result"); dir.create(input, recursive = TRUE)
  dir.create(result)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  methods <- data.table::as.data.table(singlesample_method_roster())[
    estimand == "within"
  ]
  units <- data.table::data.table(
    unit_id = sprintf("u%02d", seq_len(34L)), outer_k = 3L,
    primary_unit = c(rep(TRUE, 33L), FALSE), modality = "seq",
    biospecimen = "serum"
  )
  eligibility <- data.table::CJ(method_id = methods$method_id,
                                unit_id = units$unit_id)
  eligibility[, `:=`(eligible = FALSE, reason = "structural")]
  pair_index <- utils::combn(seq_len(nrow(methods)), 2L)
  pairs <- data.table::data.table(
    pair_id = paste(methods$method_id[pair_index[1L, ]],
                    methods$method_id[pair_index[2L, ]], sep = "::"),
    method_a = methods$method_id[pair_index[1L, ]],
    method_b = methods$method_id[pair_index[2L, ]],
    roster_order_a = pair_index[1L, ], roster_order_b = pair_index[2L, ]
  )
  env$singlesample_within_method_pairs <- function(roster) pairs[, .(
    pair_id, method_a, method_b, roster_order_a, roster_order_b
  )]
  pairs[, `:=`(k_support = 0L, support_units = "", pooled_testable = FALSE)]
  data.table::fwrite(methods, file.path(input, "methods.tsv"), sep = "\t")
  data.table::fwrite(units, file.path(input, "units.tsv"), sep = "\t")
  data.table::fwrite(data.table::data.table(x = seq_len(8500L)),
                     file.path(input, "cells.tsv"), sep = "\t")
  data.table::fwrite(eligibility, file.path(input, "eligibility.tsv"), sep = "\t")
  data.table::fwrite(pairs, file.path(input, "pair_family.tsv"), sep = "\t")
  data.table::fwrite(data.table::data.table(unit_id = "u01", seed = 101L,
    fold = 1L, n_test_profiles = 1L, n_test_groups = 1L),
    file.path(input, "splits.tsv"), sep = "\t")
  data.table::fwrite(data.table::data.table(method_id = "m", unit_id = "u",
    seed = 101L, fold = 1L, split_id = "s", sample_id = "p", group_id = "g",
    y = 0L, score = 0), file.path(input, "predictions.tsv.gz"), sep = "\t",
    compress = "gzip")
  expect_error(env$.fairqa_validate_science(
    input, result, "full", 10000L, 1000L, 20260718L
  ), "non-empty frozen design")
})
