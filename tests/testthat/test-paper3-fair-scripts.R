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
  receipt[check_id == "within_no_group_split", pass := FALSE]
  data.table::fwrite(receipt, td, sep = "\t")
  expect_error(env$.snapshot_validate_upstream_receipt(td),
               "does not pass")
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
    n_test_groups = 4L
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
})

test_that("full fair-run arguments cannot weaken registered resampling", {
  env <- paper3_script_env("paper3_fair_method_comparison.R")
  ok <- env$.fair_args(c(
    "--input-dir=/tmp/in", "--output-dir=/tmp/out",
    "--bootstrap-reps=10000", "--permutation-reps=1000",
    "--seed=20260718", "--run-mode=full"
  ))
  expect_identical(ok$bootstrap_reps, 10000L)
  expect_error(env$.fair_args(c(
    "--input-dir=/tmp/in", "--output-dir=/tmp/out",
    "--bootstrap-reps=100", "--permutation-reps=100",
    "--seed=1", "--run-mode=full"
  )), "Full mode requires")
})

test_that("nested-selector task arguments freeze the registered bootstrap", {
  env <- paper3_script_env("paper3_nested_selector_unit.R")
  args <- c(
    "--cache=/tmp/cache.rds",
    paste0("--cache-sha256=", paste(rep("a", 64L), collapse = "")),
    "--unit=u1", "--seed=101", "--output-dir=/tmp/out",
    paste0("--package-commit=", paste(rep("b", 40L), collapse = "")),
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
      source_accessions = expected[[i]]
    )
  }), expected)
  units <- env$.nested_validate_cache(cache)
  tasks <- env$.nested_task_grid(units)
  expect_equal(nrow(units), 34L)
  expect_equal(nrow(tasks), 170L)
  expect_identical(anyDuplicated(tasks[, c("unit_id", "seed")]), 0L)
  expect_setequal(unique(tasks$seed), c(101L, 202L, 303L, 404L, 505L))
  expect_equal(sum(units$outer_k == 3L), 8L)
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
