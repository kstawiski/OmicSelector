dependency_witness_script_env <- function() {
  env <- new.env(parent = globalenv())
  path <- testthat::test_path(
    "..", "..", "inst", "scripts",
    "paper3_build_dependency_structural_witness.R"
  )
  if (!file.exists(path)) {
    path <- system.file(
      "scripts", "paper3_build_dependency_structural_witness.R",
      package = "OmicSelector"
    )
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("Could not locate dependency structural-witness producer.",
         call. = FALSE)
  }
  sys.source(path, envir = env)
  env
}

dependency_witness_replay_fixture <- function() {
  X <- matrix(seq_len(24L), nrow = 8L, ncol = 3L,
              dimnames = list(sprintf("sample-%d", 1:8),
                              sprintf("feature-%d", 1:3)))
  cohort <- list(
    expr_per_sample = X,
    y_bin = rep(c(0L, 1L, 0L, 1L), 2L),
    group_id = sprintf("group-%d", 1:8),
    outer_k = 2L, min_per_fold = 2L
  )
  roster <- data.frame(method_id = "test-method", fit_fn = "fit_test",
                       stringsAsFactors = FALSE)
  folds <- list(1:4, 5:8)
  grouped <- function(y, group_id, k, seed) folds
  stratified <- function(y, k, seed) folds
  fit <- function(fit_fn, X_train, y_train, meta_train, hp) {
    list(train_ids = rownames(X_train))
  }
  score <- function(method, model, X_new, meta, roster) {
    if (identical(rownames(X_new), model$train_ids)) {
      return(seq_len(nrow(X_new)))
    }
    index <- as.integer(sub("sample-", "", rownames(X_new), fixed = TRUE))
    if (all(index <= 4L)) rep(0.5, length(index)) else seq_along(index)
  }
  baseline <- function(X_train, y_train, X_test, panel_size) {
    seq_len(nrow(X_test))
  }
  list(cohort = cohort, roster = roster, grouped = grouped,
       stratified = stratified, fit = fit, score = score,
       baseline = baseline)
}

test_that("witness producer accepts an exact irregular cell-key set", {
  env <- dependency_witness_script_env()
  args <- c(
    "--bundle", "bundle.rds", "--tasks", "tasks.tsv", "--output", "out.rds",
    "--package-root", "package", "--analysis-root", "analysis",
    "--cache", "cache.rds", "--expected-method", "lrt-deepmaha",
    "--expected-package-version", "2.6.5",
    "--expected-package-commit", paste(rep("a", 40L), collapse = ""),
    "--expected-cache-commit", paste(rep("b", 40L), collapse = ""),
    "--expected-analysis-code-id", paste(rep("c", 16L), collapse = ""),
    "--expected-task-sha256", paste(rep("d", 64L), collapse = ""),
    "--expected-cache-sha256", paste(rep("e", 64L), collapse = ""),
    "--expected-fold-engine-sha256", paste(rep("f", 64L), collapse = ""),
    "--cell-keys", "GSE188627:101,GSE188627:202,GSE270497:101"
  )
  parsed <- env$.dsw_args(args)
  expect_equal(
    parsed$cell_keys,
    data.table::data.table(
      cohort = c("GSE188627", "GSE188627", "GSE270497"),
      seed = c(101L, 202L, 101L)
    )
  )
  cells <- data.table::CJ(
    method = "lrt-deepmaha",
    cohort = c("GSE188627", "GSE270497"),
    seed = c(101L, 202L)
  )
  cells[, `:=`(eligible = FALSE, ineligible_reason = "structural")]
  cells[cohort == "GSE270497" & seed == 202L,
        `:=`(eligible = TRUE, ineligible_reason = NA_character_)]
  selected <- env$.dsw_select_cells(
    cells, "lrt-deepmaha", parsed$cell_keys
  )
  expect_equal(nrow(selected), 3L)
  expect_false(any(selected$cohort == "GSE270497" & selected$seed == 202L))
  expect_error(
    env$.dsw_select_cells(
      cells, "lrt-deepmaha",
      data.table::data.table(cohort = "GSE270497", seed = 202L)
    ),
    "exact non-GSE83977 ineligible"
  )
  duplicate_args <- args
  duplicate_args[[match("--cell-keys", duplicate_args) + 1L]] <-
    "GSE188627:101,GSE188627:101"
  expect_error(env$.dsw_args(duplicate_args), "unique comma-separated")
})

test_that("witness replay preserves every finite fold prediction", {
  env <- dependency_witness_script_env()
  x <- dependency_witness_replay_fixture()
  out <- env$.dsw_replay_cell(
    method = "test-method", cohort_id = "GSETEST", seed = 101L,
    cohort = x$cohort, roster = x$roster,
    grouped_fold_fn = x$grouped, stratified_fold_fn = x$stratified,
    fit_call_fn = x$fit, score_call_fn = x$score,
    baseline_call_fn = x$baseline, outer_k = 2L,
    min_pos_per_fold = 2L, min_neg_per_fold = 2L
  )
  expect_equal(nrow(out$predictions), 8L)
  expect_true(all(is.finite(out$predictions$score_m)))
  expect_true(all(is.finite(out$predictions$score_b)))
  expect_identical(out$fold_audit$status,
                   c("structural_constant", "valid"))
  expect_identical(out$n_valid_folds, 1L)
  expect_identical(
    out$ineligible_reason,
    paste0("n_valid_folds=1 < 2; per-fold: ",
           "method/baseline constant or all-NA on test fold; ok")
  )
})

test_that("witness replay makes fit and score failures fatal", {
  env <- dependency_witness_script_env()
  x <- dependency_witness_replay_fixture()
  replay <- function(fit = x$fit, score = x$score,
                     baseline = x$baseline, grouped = x$grouped) {
    env$.dsw_replay_cell(
      method = "test-method", cohort_id = "GSETEST", seed = 101L,
      cohort = x$cohort, roster = x$roster,
      grouped_fold_fn = grouped, stratified_fold_fn = x$stratified,
      fit_call_fn = fit, score_call_fn = score,
      baseline_call_fn = baseline, outer_k = 2L,
      min_pos_per_fold = 2L, min_neg_per_fold = 2L
    )
  }
  expect_error(
    replay(fit = function(...) stop("optimizer failed")),
    "fit failure.*optimizer failed"
  )
  expect_error(
    replay(baseline = function(...) stop("baseline failed")),
    "baseline scoring failure.*baseline failed"
  )
  expect_error(
    replay(score = function(...) stop("scorer failed")),
    "training scoring failure.*scorer failed"
  )
  expect_error(
    replay(score = function(method, model, X_new, meta, roster) {
      rep(NA_real_, nrow(X_new))
    }),
    "non-finite"
  )
})

test_that("witness replay rejects provenance group leakage", {
  env <- dependency_witness_script_env()
  x <- dependency_witness_replay_fixture()
  x$cohort$group_id[c(1L, 5L)] <- "reused-patient"
  expect_error(
    env$.dsw_replay_cell(
      method = "test-method", cohort_id = "GSETEST", seed = 101L,
      cohort = x$cohort, roster = x$roster,
      grouped_fold_fn = x$grouped, stratified_fold_fn = x$stratified,
      fit_call_fn = x$fit, score_call_fn = x$score,
      baseline_call_fn = x$baseline, outer_k = 2L,
      min_pos_per_fold = 2L, min_neg_per_fold = 2L
    ),
    "provenance group split"
  )
})

test_that("witness replay requires exact groups and unique profile IDs", {
  env <- dependency_witness_script_env()
  x <- dependency_witness_replay_fixture()
  replay <- function(cohort) {
    env$.dsw_replay_cell(
      method = "test-method", cohort_id = "GSETEST", seed = 101L,
      cohort = cohort, roster = x$roster,
      grouped_fold_fn = x$grouped, stratified_fold_fn = x$stratified,
      fit_call_fn = x$fit, score_call_fn = x$score,
      baseline_call_fn = x$baseline, outer_k = 2L,
      min_pos_per_fold = 2L, min_neg_per_fold = 2L
    )
  }

  malformed_groups <- x$cohort
  malformed_groups$group_id <- malformed_groups$group_id[-1L]
  expect_error(replay(malformed_groups), "invalid provenance groups")

  missing_groups <- x$cohort
  missing_groups$group_id <- NULL
  expect_error(replay(missing_groups), "invalid provenance groups")

  duplicated_profiles <- x$cohort
  rownames(duplicated_profiles$expr_per_sample)[[8L]] <-
    rownames(duplicated_profiles$expr_per_sample)[[1L]]
  expect_error(replay(duplicated_profiles), "invalid profiles.*sample IDs")
})
