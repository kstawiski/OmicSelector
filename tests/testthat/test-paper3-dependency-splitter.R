dependency_split_script_env <- function() {
  env <- new.env(parent = globalenv())
  source_path <- testthat::test_path(
    "..", "..", "inst", "scripts",
    "paper3_split_dependency_within_bundle.R"
  )
  if (!file.exists(source_path)) {
    source_path <- system.file(
      "scripts", "paper3_split_dependency_within_bundle.R",
      package = "OmicSelector"
    )
  }
  if (!nzchar(source_path) || !file.exists(source_path)) {
    stop("Could not locate dependency-bundle splitter.", call. = FALSE)
  }
  sys.source(source_path, envir = env)
  env
}

dependency_split_fixture <- function() {
  methods <- c(
    "ai-scarf", "ai-tabpfn", "coda-codacore", "coda-deepcoda",
    "img-gasfcnn", "inv-scatter", "lrt-deepmaha", "nc-ecod-copod",
    "proto-net", "ssl-vicreg", "tab-tabdpt", "tab-tabicl", "unc-sngp"
  )
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  cohorts <- c("GSE83977", sprintf("GSE%06d", seq_len(33L)))
  tasks <- data.table::CJ(method = methods, seed = seeds, cohort = cohorts)
  tasks[, label := sprintf("%s__seed_%d__%s", method, seed, cohort)]
  method <- "ai-scarf"
  cells <- data.table::CJ(method = method, seed = seeds, cohort = cohorts)
  cells[, `:=`(
    eligible = TRUE, ineligible_reason = NA_character_,
    auc_method = 0.75, auc_baseline = 0.70, lift = 0.05,
    lift_se = 0.01, n_eff = 4L, n_valid_folds = 2L,
    min_pos = 1L, min_neg = 1L, outer_k = 2L,
    grouping_policy = "groupKFold", n_group_split_violations = 0L,
    package_version = "2.6.5",
    package_commit = "6455bfe0937a32822ccad003eb8413e47feb2b05",
    cache_package_commit = "cc483baf93c63159a4cf22039f4269c0aaa33bbb",
    analysis_code_id = "3d7ba1b678d44884"
  )]
  predictions <- cells[, {
    idx <- seq_len(4L)
    .(
      fold = rep(1:2, each = 2L), n_train = 2L, n_test = 2L,
      sample_idx = idx,
      sample_id = sprintf("%s-%d-%s-s%d", method, seed, cohort, idx),
      group_id = sprintf("%s-%d-%s-g%d", method, seed, cohort, idx),
      y = rep(c(0L, 1L), 2L), score_m = c(0.1, 0.9, 0.2, 0.8),
      score_b = c(0.2, 0.8, 0.3, 0.7),
      package_commit = package_commit[[1L]],
      analysis_code_id = analysis_code_id[[1L]]
    )
  }, by = .(method, cohort, seed)]
  bundle <- list(
    cells = cells, predictions = predictions, seeds = seeds,
    methods = method, cohorts = cohorts, package_version = "2.6.5",
    package_commit = "6455bfe0937a32822ccad003eb8413e47feb2b05",
    cache_package_commit = "cc483baf93c63159a4cf22039f4269c0aaa33bbb",
    analysis_code_id = "3d7ba1b678d44884"
  )
  expected <- list(
    method = method, package_version = bundle$package_version,
    package_commit = bundle$package_commit,
    cache_commit = bundle$cache_package_commit,
    analysis_code_id = bundle$analysis_code_id
  )
  list(tasks = tasks, bundle = bundle, expected = expected)
}

dependency_structural_witness_fixture <- function() {
  fixture <- dependency_split_fixture()
  source_hashes <- vapply(
    c("runner", "engine", "auc", "design"),
    function(x) digest::digest(x, algo = "sha256"), character(1L)
  )
  analysis_code_id <- substr(digest::digest(
    paste(source_hashes, collapse = "|"), algo = "sha256"
  ), 1L, 16L)
  fixture$bundle$analysis_code_id <- analysis_code_id
  fixture$bundle$cells$analysis_code_id <- analysis_code_id
  fixture$bundle$predictions$analysis_code_id <- analysis_code_id
  fixture$expected$analysis_code_id <- analysis_code_id

  method <- fixture$expected$method
  cohort <- "GSE000001"
  seed <- 101L
  method_value <- method
  cohort_value <- cohort
  seed_value <- seed
  reason <- paste0(
    "n_valid_folds=1 < 2; per-fold: ",
    "method/baseline constant or all-NA on test fold; ok"
  )
  fixture$bundle$cells[
    method == method_value & cohort == cohort_value & seed == seed_value,
    `:=`(eligible = FALSE, ineligible_reason = reason, n_eff = NA_integer_,
         n_valid_folds = 1L, outer_k = 2L,
         n_group_split_violations = 0L)
  ]
  fixture$bundle$predictions <- fixture$bundle$predictions[
    !(method == method_value & cohort == cohort_value & seed == seed_value)
  ]

  predictions <- data.table::data.table(
    method = method, cohort = cohort, seed = seed,
    fold = rep(1:2, each = 4L), n_train = 4L, n_test = 4L,
    sample_idx = 1:8, sample_id = sprintf("witness-s%d", 1:8),
    group_id = sprintf("witness-g%d", 1:8),
    y = rep(c(0L, 1L, 0L, 1L), 2L),
    score_m = c(rep(0.5, 4L), 0.1, 0.9, 0.2, 0.8),
    score_b = c(0.1, 0.9, 0.2, 0.8, 0.2, 0.8, 0.3, 0.7),
    package_commit = fixture$expected$package_commit,
    cache_package_commit = fixture$expected$cache_commit,
    analysis_code_id = fixture$expected$analysis_code_id
  )
  audit <- predictions[, .(
    n_train = unique(n_train), n_test = .N,
    n_pos_test = sum(y == 1L), n_neg_test = sum(y == 0L),
    n_unique_method = data.table::uniqueN(score_m),
    n_unique_baseline = data.table::uniqueN(score_b),
    n_unique_outcome = data.table::uniqueN(y), n_train_finite = 4L,
    n_test_finite = sum(is.finite(score_m) & is.finite(score_b)),
    training_direction = 1,
    n_group_split_violations = 0L,
    fit_status = "ok", baseline_status = "ok",
    training_score_status = "ok", heldout_score_status = "ok",
    status = if (.BY$fold == 1L) "structural_constant" else "valid",
    structural_basis = if (.BY$fold == 1L) "method_constant" else NA_character_,
    package_commit = unique(package_commit),
    cache_package_commit = unique(cache_package_commit),
    analysis_code_id = unique(analysis_code_id)
  ), by = .(method, cohort, seed, fold)]
  paths <- c(
    "code/analyses/run_singlesample_public_benchmark_20260715.R",
    "code/methods/singlesample_roster_benchmark.R",
    "code/methods/paired_auc_diff_se.R",
    "code/methods/singlesample_public_design.R"
  )
  sources <- data.table::data.table(
    role = c("benchmark_runner", "benchmark_engine", "paired_auc_engine",
             "design_contract"),
    relative_path = paths, sha256 = source_hashes
  )
  source_bundle_sha <- paste(rep("a", 64L), collapse = "")
  task_sha <- paste(rep("b", 64L), collapse = "")
  cache_sha <- paste(rep("c", 64L), collapse = "")
  producer_path <- testthat::test_path(
    "..", "..", "inst", "scripts",
    "paper3_build_dependency_structural_witness.R"
  )
  if (!file.exists(producer_path)) {
    producer_path <- system.file(
      "scripts", "paper3_build_dependency_structural_witness.R",
      package = "OmicSelector"
    )
  }
  if (!nzchar(producer_path) || !file.exists(producer_path)) {
    stop("Could not locate structural-witness producer.", call. = FALSE)
  }
  producer_sha <- digest::digest(producer_path, algo = "sha256", file = TRUE)
  fold_engine_sha <- paste(rep("e", 64L), collapse = "")
  witness <- list(
    schema_version = "OmicSelector-dependency-structural-witness-v1",
    producer = "paper3_build_dependency_structural_witness.R",
    producer_script_sha256 = producer_sha,
    source_bundle_sha256 = source_bundle_sha,
    task_table_sha256 = task_sha, cache_sha256 = cache_sha,
    fold_engine_sha256 = fold_engine_sha,
    package_version = fixture$expected$package_version,
    package_commit = fixture$expected$package_commit,
    cache_package_commit = fixture$expected$cache_commit,
    analysis_code_id = fixture$expected$analysis_code_id,
    analysis_sources = sources,
    replay_sources = rbind(
      sources,
      data.table::data.table(
        role = "fold_engine",
        relative_path = "code/methods/matched_null_benchmark.R",
        sha256 = fold_engine_sha
      )
    ),
    cells = data.table::data.table(
      method = method, cohort = cohort, seed = seed,
      ineligible_reason = reason, n_valid_folds = 1L, outer_k = 2L,
      n_group_split_violations = 0L, n_profiles = 8L
    ),
    predictions = predictions, fold_audit = audit,
    runtime = list(R.version = R.version.string)
  )
  list(fixture = fixture, witness = witness,
       source_bundle_sha = source_bundle_sha, task_sha = task_sha,
       cache_sha = cache_sha, fold_engine_sha = fold_engine_sha,
       witness_sha = paste(rep("f", 64L), collapse = ""))
}

test_that("dependency splitter validates and atomically emits 170 shards", {
  env <- dependency_split_script_env()
  fixture <- dependency_split_fixture()
  td <- tempfile("dependency-split-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  bundle_path <- file.path(td, "source.rds")
  tasks_path <- file.path(td, "tasks.tsv")
  output_root <- file.path(td, "shards")
  saveRDS(fixture$bundle, bundle_path, version = 3L)
  data.table::fwrite(fixture$tasks, tasks_path, sep = "\t", col.names = FALSE)
  args <- c(
    "--bundle", bundle_path, "--tasks", tasks_path,
    "--output-root", output_root,
    "--expected-method", fixture$expected$method,
    "--expected-package-version", fixture$expected$package_version,
    "--expected-package-commit", fixture$expected$package_commit,
    "--expected-cache-commit", fixture$expected$cache_commit,
    "--expected-analysis-code-id", fixture$expected$analysis_code_id,
    "--expected-task-sha256",
    digest::digest(tasks_path, algo = "sha256", file = TRUE)
  )

  wrong_hash_args <- args
  hash_index <- match("--expected-task-sha256", wrong_hash_args) + 1L
  wrong_hash_args[[hash_index]] <- paste(rep("0", 64L), collapse = "")
  expect_error(
    env$.dependency_split_main(wrong_hash_args),
    "does not match --expected-task-sha256"
  )
  expect_invisible(env$.dependency_split_main(args))
  manifest <- data.table::fread(file.path(
    output_root, "output_manifest__ai-scarf.tsv"
  ))
  expect_equal(nrow(manifest), 170L)
  expect_true(all(file.exists(file.path(output_root, manifest$path))))
  one <- readRDS(file.path(output_root, manifest$path[[1L]]))
  expect_equal(nrow(one$cells), 1L)
  expect_equal(nrow(one$predictions), 4L)
  expect_identical(one$methods, fixture$expected$method)
  expect_identical(one$seeds, as.integer(one$cells$seed))
  expect_identical(one$cohorts, as.character(one$cells$cohort))
  expect_error(env$.dependency_split_main(args), "already exist")
  expect_invisible(env$.dependency_split_main(
    c(args, "--allow-identical-overlap")
  ))
})

test_that("dependency splitter fails closed on cell and prediction defects", {
  env <- dependency_split_script_env()
  fixture <- dependency_split_fixture()
  tasks <- data.table::copy(fixture$tasks)

  td <- tempfile("dependency-task-contract-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  wrong_tasks <- data.table::copy(tasks)
  wrong_tasks[method == "unc-sngp", `:=`(
    method = "fake-route",
    label = sprintf("fake-route__seed_%d__%s", seed, cohort)
  )]
  wrong_tasks_path <- file.path(td, "wrong_tasks.tsv")
  data.table::fwrite(wrong_tasks, wrong_tasks_path, sep = "\t",
                     col.names = FALSE)
  expect_error(
    env$.dependency_read_tasks(wrong_tasks_path),
    "complete 13-by-five-by-34 design"
  )

  structural <- fixture$bundle
  structural$cells <- data.table::copy(structural$cells)
  first <- structural$cells[cohort == "GSE83977" & seed == 101L,
                            .(method, cohort, seed)]
  structural_reason <- paste0(
    "n_valid_folds=0 < 3; per-fold: ",
    paste(rep("method/baseline constant or all-NA on test fold", 3L),
          collapse = "; ")
  )
  structural$cells[
    method == first$method & cohort == first$cohort & seed == first$seed,
    `:=`(eligible = FALSE, ineligible_reason = structural_reason,
         n_eff = NA_integer_, n_valid_folds = 0L, outer_k = 3L,
         n_group_split_violations = 0L)
  ]
  structural$predictions <- data.table::copy(structural$predictions)[
    !(method == first$method & cohort == first$cohort & seed == first$seed)
  ]
  expect_invisible(
    env$.dependency_validate_bundle(structural, tasks, fixture$expected)
  )

  broken <- structural
  broken$cells <- data.table::copy(structural$cells)
  broken$cells[cohort == "GSE83977" & seed == 101L,
               ineligible_reason := paste0(
                 "n_valid_folds=0 < 3; per-fold: fit failed: ",
                 "huggingface_hub.errors.LocalEntryNotFoundError: ",
                 "checkpoint not found in the local cache"
               )]
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "dependency/runtime failure"
  )

  broken <- structural
  broken$cells <- data.table::copy(structural$cells)
  broken$cells[cohort == "GSE83977" & seed == 101L,
               ineligible_reason := "n_valid_folds=0 < 3; per-fold: fit failed: optimizer error"]
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "unreviewed ineligibility"
  )

  broken <- structural
  broken$cells <- data.table::copy(structural$cells)
  broken$cells[cohort == "GSE000001" & seed == 101L,
               `:=`(eligible = FALSE, ineligible_reason = structural_reason,
                    n_eff = NA_integer_, n_valid_folds = 0L, outer_k = 3L,
                    n_group_split_violations = 0L)]
  broken$predictions <- data.table::copy(structural$predictions)[
    !(cohort == "GSE000001" & seed == 101L)
  ]
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "unreviewed ineligibility"
  )

  broken <- fixture$bundle
  broken$cells <- data.table::copy(broken$cells)
  broken$cells$ineligible_reason[[1L]] <- "driver cell error: hidden failure"
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "driver/cell error"
  )

  broken <- fixture$bundle
  broken$predictions <- data.table::copy(broken$predictions)
  first_cell <- broken$predictions[1L, .(method, cohort, seed)]
  idx <- broken$predictions[
    method == first_cell$method & cohort == first_cell$cohort &
      seed == first_cell$seed, which = TRUE
  ]
  broken$predictions$group_id[idx[c(1L, 3L)]] <- "cross-fold-group"
  broken$predictions$y[idx[c(1L, 3L)]] <- 0L
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "crosses held-out folds"
  )

  broken <- fixture$bundle
  broken$predictions <- data.table::copy(broken$predictions)[-1L]
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "cardinality"
  )

  broken <- fixture$bundle
  broken$cells <- data.table::copy(broken$cells)
  broken$cells[cohort == "GSE83977" & seed == 101L,
               `:=`(eligible = FALSE, ineligible_reason = structural_reason,
                    n_eff = NA_integer_, n_valid_folds = 0L, outer_k = 3L,
                    n_group_split_violations = 0L)]
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "correspond exactly to eligible"
  )

  broken <- fixture$bundle
  broken$predictions <- data.table::copy(broken$predictions)
  broken$predictions$package_commit[[1L]] <- paste(rep("0", 40L), collapse = "")
  expect_error(
    env$.dependency_validate_bundle(broken, tasks, fixture$expected),
    "package_commit pin"
  )
})

test_that("dependency splitter rejects a non-identical overlap restart", {
  env <- dependency_split_script_env()
  fixture <- dependency_split_fixture()
  td <- tempfile("dependency-overlap-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  bundle_path <- file.path(td, "source.rds")
  tasks_path <- file.path(td, "tasks.tsv")
  output_root <- file.path(td, "shards")
  saveRDS(fixture$bundle, bundle_path, version = 3L)
  data.table::fwrite(fixture$tasks, tasks_path, sep = "\t", col.names = FALSE)
  args <- c(
    "--bundle", bundle_path, "--tasks", tasks_path,
    "--output-root", output_root,
    "--expected-method", fixture$expected$method,
    "--expected-package-version", fixture$expected$package_version,
    "--expected-package-commit", fixture$expected$package_commit,
    "--expected-cache-commit", fixture$expected$cache_commit,
    "--expected-analysis-code-id", fixture$expected$analysis_code_id,
    "--expected-task-sha256",
    digest::digest(tasks_path, algo = "sha256", file = TRUE)
  )
  env$.dependency_split_main(args)
  manifest <- data.table::fread(file.path(
    output_root, "output_manifest__ai-scarf.tsv"
  ))
  corrupt_path <- file.path(output_root, manifest$path[[1L]])
  corrupt <- readRDS(corrupt_path)
  corrupt$cells$auc_method <- 0
  saveRDS(corrupt, corrupt_path, version = 3L)
  expect_error(
    env$.dependency_split_main(c(args, "--allow-identical-overlap")),
    "not an exact identical restart"
  )
})

test_that("dependency splitter authorizes only an exact structural witness", {
  env <- dependency_split_script_env()
  x <- dependency_structural_witness_fixture()
  expect_error(
    env$.dependency_validate_bundle(
      x$fixture$bundle, x$fixture$tasks, x$fixture$expected
    ),
    "unreviewed ineligibility"
  )
  validated <- env$.dependency_validate_bundle(
    x$fixture$bundle, x$fixture$tasks, x$fixture$expected,
    structural_witness = x$witness,
    source_bundle_sha256 = x$source_bundle_sha,
    task_table_sha256 = x$task_sha,
    expected_cache_sha256 = x$cache_sha,
    expected_fold_engine_sha256 = x$fold_engine_sha,
    witness_sha256 = x$witness_sha
  )
  expect_equal(nrow(validated$structural_witness$keys), 1L)
  cell <- validated$cells[
    cohort == "GSE000001" & seed == 101L
  ]
  shard <- env$.dependency_shard_object(
    x$fixture$bundle, cell, validated$predictions[0L],
    x$fixture$expected, validated$structural_witness
  )
  expect_true("structural_ineligibility_witness" %in% names(shard))
  expect_equal(nrow(shard$structural_ineligibility_witness$predictions), 8L)
  expect_identical(
    shard$structural_ineligibility_witness$witness_sha256,
    x$witness_sha
  )
})

test_that("structural witness rejects stale, partial, forged, and error cells", {
  env <- dependency_split_script_env()
  x <- dependency_structural_witness_fixture()
  validate <- function(witness = x$witness,
                       source_sha = x$source_bundle_sha,
                       cache_sha = x$cache_sha) {
    env$.dependency_validate_bundle(
      x$fixture$bundle, x$fixture$tasks, x$fixture$expected,
      structural_witness = witness,
      source_bundle_sha256 = source_sha,
      task_table_sha256 = x$task_sha,
      expected_cache_sha256 = cache_sha,
      expected_fold_engine_sha256 = x$fold_engine_sha,
      witness_sha256 = x$witness_sha
    )
  }

  stale <- data.table::copy(x$witness)
  stale$source_bundle_sha256 <- paste(rep("0", 64L), collapse = "")
  expect_error(validate(stale), "source_bundle_sha256")

  stale <- data.table::copy(x$witness)
  stale$cache_sha256 <- paste(rep("0", 64L), collapse = "")
  expect_error(validate(stale), "cache_sha256")

  stale <- data.table::copy(x$witness)
  stale$package_commit <- paste(rep("0", 40L), collapse = "")
  expect_error(validate(stale), "package_commit")

  stale <- data.table::copy(x$witness)
  stale$fold_engine_sha256 <- paste(rep("0", 64L), collapse = "")
  expect_error(validate(stale), "fold_engine_sha256")

  stale <- data.table::copy(x$witness)
  stale$replay_sources[role == "fold_engine", sha256 :=
                         paste(rep("0", 64L), collapse = "")]
  expect_error(validate(stale), "fold-engine SHA")

  forged_closure <- data.table::copy(x$witness)
  forged_closure$replay_sources <- rbind(
    forged_closure$replay_sources,
    forged_closure$replay_sources[1L]
  )
  expect_error(validate(forged_closure), "replay source closure")

  forged_producer <- data.table::copy(x$witness)
  forged_producer$producer_script_sha256 <- paste(rep("0", 64L), collapse = "")
  expect_error(validate(forged_producer), "package-owned producer bytes")

  partial <- data.table::copy(x$witness)
  partial$predictions <- partial$predictions[fold == 1L]
  expect_error(validate(partial), "profile/fold coverage|exact cell set")

  partial <- data.table::copy(x$witness)
  partial$fold_audit <- partial$fold_audit[fold == 1L]
  expect_error(validate(partial), "profile/fold coverage")

  forged <- data.table::copy(x$witness)
  forged$predictions$score_m[[1L]] <- NA_real_
  expect_error(validate(forged), "non-finite")

  forged <- data.table::copy(x$witness)
  forged$predictions[fold == 1L, score_m := seq_len(.N)]
  expect_error(validate(forged), "fold audit does not derive")

  errored <- data.table::copy(x$witness)
  errored$fold_audit$fit_status[[1L]] <- "error: optimizer failed"
  expect_error(validate(errored), "contains an error")

  errored <- data.table::copy(x$witness)
  errored$fold_audit$heldout_score_status[[1L]] <- "error: OOM"
  expect_error(validate(errored), "contains an error")

  leaking <- data.table::copy(x$witness)
  leaking$predictions$group_id[c(1L, 5L)] <- "cross-fold-group"
  leaking$predictions$y[c(1L, 5L)] <- 0L
  expect_error(validate(leaking), "splits a provenance group")

  mismatch <- data.table::copy(x$witness)
  mismatch$cells$ineligible_reason <- paste0(
    "n_valid_folds=0 < 2; per-fold: ",
    paste(rep("method/baseline constant or all-NA on test fold", 2L),
          collapse = "; ")
  )
  expect_error(validate(mismatch), "do not match source bundle")

  wrong_sources <- data.table::copy(x$witness)
  wrong_sources$analysis_sources$sha256[[1L]] <-
    paste(rep("0", 64L), collapse = "")
  wrong_sources$replay_sources$sha256[[1L]] <-
    wrong_sources$analysis_sources$sha256[[1L]]
  expect_error(validate(wrong_sources), "do not reproduce analysis_code_id")
})

test_that("splitter main requires the exact witness byte pin", {
  env <- dependency_split_script_env()
  x <- dependency_structural_witness_fixture()
  td <- tempfile("dependency-witness-pin-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  bundle_path <- file.path(td, "source.rds")
  tasks_path <- file.path(td, "tasks.tsv")
  witness_path <- file.path(td, "witness.rds")
  saveRDS(x$fixture$bundle, bundle_path, version = 3L)
  saveRDS(x$witness, witness_path, version = 3L)
  data.table::fwrite(x$fixture$tasks, tasks_path, sep = "\t",
                     col.names = FALSE)
  args <- c(
    "--bundle", bundle_path, "--tasks", tasks_path,
    "--output-root", file.path(td, "shards"),
    "--expected-method", x$fixture$expected$method,
    "--expected-package-version", x$fixture$expected$package_version,
    "--expected-package-commit", x$fixture$expected$package_commit,
    "--expected-cache-commit", x$fixture$expected$cache_commit,
    "--expected-analysis-code-id", x$fixture$expected$analysis_code_id,
    "--expected-task-sha256",
    digest::digest(tasks_path, algo = "sha256", file = TRUE),
    "--structural-witness", witness_path,
    "--expected-structural-witness-sha256",
    paste(rep("0", 64L), collapse = ""),
    "--expected-structural-cache-sha256", x$cache_sha,
    "--expected-structural-fold-engine-sha256", x$fold_engine_sha
  )
  expect_error(env$.dependency_split_main(args),
               "does not match its expected SHA-256")
})
