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
  cohorts <- sprintf("GSE%06d", seq_len(34L))
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
  first <- structural$cells[1L, .(method, cohort, seed)]
  structural$cells[
    method == first$method & cohort == first$cohort & seed == first$seed,
    `:=`(eligible = FALSE, ineligible_reason = "structurally ineligible",
         n_eff = NA_integer_, n_valid_folds = 0L,
         n_group_split_violations = NA_integer_)
  ]
  structural$predictions <- data.table::copy(structural$predictions)[
    !(method == first$method & cohort == first$cohort & seed == first$seed)
  ]
  expect_invisible(
    env$.dependency_validate_bundle(structural, tasks, fixture$expected)
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
  broken$cells[1L, `:=`(eligible = FALSE,
                        ineligible_reason = "structurally ineligible")]
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
