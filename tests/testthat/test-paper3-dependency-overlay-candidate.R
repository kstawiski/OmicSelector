overlay_candidate_script_env <- function() {
  env <- new.env(parent = globalenv())
  source_path <- testthat::test_path(
    "..", "..", "inst", "scripts",
    "paper3_build_dependency_overlay_candidate.R"
  )
  if (!file.exists(source_path)) {
    source_path <- system.file(
      "scripts", "paper3_build_dependency_overlay_candidate.R",
      package = "OmicSelector"
    )
  }
  if (!nzchar(source_path) || !file.exists(source_path)) {
    stop("Could not locate dependency overlay candidate producer.",
         call. = FALSE)
  }
  sys.source(source_path, envir = env)
  env
}

overlay_candidate_fixture <- function(env) {
  td <- tempfile("dependency-overlay-candidate-")
  dir.create(td)
  shard_root <- file.path(td, "within_shards")
  transfer_root <- file.path(td, "transfer_shards")
  within_receipt <- file.path(td, "within_s5_receipt")
  transfer_receipt <- file.path(td, "transfer_argos_receipt")
  dir.create(shard_root)
  dir.create(transfer_root)
  dir.create(within_receipt)
  dir.create(transfer_receipt)

  producer_root <- file.path(td, "producer_checkout")
  producer_script_dir <- file.path(producer_root, "inst", "scripts")
  dir.create(producer_script_dir, recursive = TRUE)
  producer_files <- c(
    "paper3_build_dependency_overlay_candidate.R",
    "paper3_split_dependency_within_bundle.R",
    "paper3_build_dependency_structural_witness.R"
  )
  source_script_dir <- testthat::test_path("..", "..", "inst", "scripts")
  if (!dir.exists(source_script_dir)) {
    source_script_dir <- dirname(system.file(
      "scripts", "paper3_build_dependency_overlay_candidate.R",
      package = "OmicSelector"
    ))
  }
  stopifnot(all(file.copy(
    file.path(source_script_dir, producer_files),
    file.path(producer_script_dir, producer_files)
  )))
  writeLines(c(
    "Package: OmicSelector", "Version: 2.6.5.9000",
    "Title: Test Producer Checkout", "Description: Test fixture.",
    "License: MIT"
  ), file.path(producer_root, "DESCRIPTION"))
  git <- function(args) {
    out <- system2("git", c("-C", producer_root, args),
                   stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(out, "status"))) stop(paste(out, collapse = "\n"))
    out
  }
  git(c("init", "-q"))
  git(c("config", "user.email", "fixture@example.org"))
  git(c("config", "user.name", "Fixture"))
  git(c("add", "."))
  git(c("commit", "-q", "-m", "fixture-producer"))
  producer_commit <- trimws(git(c("rev-parse", "HEAD")))
  env$.OVERLAY_PRODUCER_PATH <- normalizePath(file.path(
    producer_script_dir, "paper3_build_dependency_overlay_candidate.R"
  ))

  package_version <- "2.6.5"
  package_commit <- "6455bfe0937a32822ccad003eb8413e47feb2b05"
  cache_commit <- "cc483baf93c63159a4cf22039f4269c0aaa33bbb"
  analysis_source_hashes <- vapply(
    c("runner", "engine", "auc", "design"),
    digest::digest, character(1L), algo = "sha256"
  )
  within_code_id <- substr(digest::digest(
    paste(analysis_source_hashes, collapse = "|"), algo = "sha256"
  ), 1L, 16L)
  transfer_code_id <- "a94cf7b24519b360"
  structural_cache_sha <- paste(rep("b", 64L), collapse = "")
  fold_engine_sha <- paste(rep("c", 64L), collapse = "")
  methods <- env$.overlay_within_methods
  transfer_methods <- env$.overlay_transfer_methods
  seeds <- env$.overlay_seeds
  cohorts <- sort(c("GSE188627", "GSE270497",
                    sprintf("GSE%06d", seq_len(32L))))
  within_tasks <- data.table::CJ(
    method = methods, seed = seeds, cohort = cohorts, sorted = TRUE
  )
  within_tasks[, label := sprintf("%s__seed_%d__%s", method, seed, cohort)]
  within_tasks_path <- file.path(td, "within_tasks.tsv")
  data.table::fwrite(within_tasks, within_tasks_path, sep = "\t",
                     col.names = FALSE, quote = FALSE)
  within_task_sha <- env$.overlay_sha256(within_tasks_path)

  witness_keys <- data.table::copy(env$.overlay_authorized_witness_keys)
  structural_reason <- paste0(
    "n_valid_folds=1 < 2; per-fold: ",
    "method/baseline constant or all-NA on test fold; ok"
  )
  witness_cells <- witness_keys[, `:=`(
    ineligible_reason = structural_reason,
    n_valid_folds = 1L, outer_k = 2L,
    n_group_split_violations = 0L, n_profiles = 8L
  )]
  witness_predictions <- witness_keys[, {
    data.table::data.table(
      fold = rep(1:2, each = 4L), n_train = 4L, n_test = 4L,
      sample_idx = seq_len(8L),
      sample_id = sprintf("%s-%s-%d-s%d", method, cohort, seed, seq_len(8L)),
      group_id = sprintf("%s-%s-%d-g%d", method, cohort, seed, seq_len(8L)),
      y = rep(c(0L, 1L, 0L, 1L), 2L),
      score_m = c(rep(0.5, 4L), 0.1, 0.9, 0.2, 0.8),
      score_b = c(0.2, 0.8, 0.3, 0.7, 0.2, 0.8, 0.3, 0.7),
      package_commit = package_commit,
      cache_package_commit = cache_commit,
      analysis_code_id = within_code_id
    )
  }, by = .(method, cohort, seed)]
  witness_audit <- witness_predictions[, .(
    n_train = unique(n_train), n_test = .N,
    n_pos_test = sum(y == 1L), n_neg_test = sum(y == 0L),
    n_unique_method = data.table::uniqueN(score_m),
    n_unique_baseline = data.table::uniqueN(score_b),
    n_unique_outcome = data.table::uniqueN(y), n_train_finite = 4L,
    n_test_finite = .N, training_direction = 1,
    n_group_split_violations = 0L, fit_status = "ok",
    baseline_status = "ok", training_score_status = "ok",
    heldout_score_status = "ok",
    status = if (.BY$fold == 1L) "structural_constant" else "valid",
    structural_basis = if (.BY$fold == 1L) "method_constant" else NA_character_,
    package_commit = unique(package_commit),
    cache_package_commit = unique(cache_package_commit),
    analysis_code_id = unique(analysis_code_id)
  ), by = .(method, cohort, seed, fold)]
  deepmaha_source_sha <- digest::digest("lrt-deepmaha-source", algo = "sha256")
  analysis_sources <- data.table::data.table(
    role = c("benchmark_runner", "benchmark_engine", "paired_auc_engine",
             "design_contract"),
    relative_path = c(
      "code/analyses/run_singlesample_public_benchmark_20260715.R",
      "code/methods/singlesample_roster_benchmark.R",
      "code/methods/paired_auc_diff_se.R",
      "code/methods/singlesample_public_design.R"
    ),
    sha256 = analysis_source_hashes
  )
  witness <- list(
    schema_version = "OmicSelector-dependency-structural-witness-v1",
    producer = "paper3_build_dependency_structural_witness.R",
    producer_script_sha256 = env$.overlay_sha256(file.path(
      producer_script_dir, "paper3_build_dependency_structural_witness.R"
    )),
    source_bundle_sha256 = deepmaha_source_sha,
    task_table_sha256 = within_task_sha,
    cache_sha256 = structural_cache_sha,
    fold_engine_sha256 = fold_engine_sha,
    package_version = package_version, package_commit = package_commit,
    cache_package_commit = cache_commit, analysis_code_id = within_code_id,
    analysis_sources = analysis_sources,
    replay_sources = data.table::rbindlist(list(
      analysis_sources,
      data.table::data.table(
        role = "fold_engine",
        relative_path = "code/methods/matched_null_benchmark.R",
        sha256 = fold_engine_sha
      )
    )),
    cells = witness_cells, predictions = witness_predictions,
    fold_audit = witness_audit, runtime = list(R.version = R.version.string)
  )
  witness_path <- file.path(td, "structural_witness.rds")
  saveRDS(witness, witness_path, version = 3L)
  witness_sha <- env$.overlay_sha256(witness_path)
  witness_key_strings <- env$.overlay_key_string(witness_keys,
                                                 c("method", "cohort", "seed"))

  manifest_rows <- vector("list", length(methods))
  within_predictions <- vector("list", nrow(within_tasks))
  for (i in seq_len(nrow(within_tasks))) {
    task <- within_tasks[i]
    key_string <- env$.overlay_key_string(task,
                                          c("method", "cohort", "seed"))
    witnessed <- key_string %in% witness_key_strings
    cell <- data.table::data.table(
      method = task$method, seed = task$seed, cohort = task$cohort,
      eligible = !witnessed,
      ineligible_reason = if (witnessed) structural_reason else NA_character_,
      n_eff = if (witnessed) NA_integer_ else 4L,
      n_valid_folds = if (witnessed) 1L else 2L, outer_k = 2L,
      n_group_split_violations = 0L,
      package_version = package_version, package_commit = package_commit,
      cache_package_commit = cache_commit, analysis_code_id = within_code_id,
      execution_route = "pre-extension-ineligible",
      final_package_version = package_version,
      final_package_commit = package_commit
    )
    prediction <- if (witnessed) {
      data.table::data.table(
        method = character(), seed = integer(), cohort = character(),
        fold = integer(), n_train = integer(), n_test = integer(),
        sample_idx = integer(), sample_id = character(),
        group_id = character(), y = integer(), score_m = numeric(),
        score_b = numeric(), package_commit = character(),
        analysis_code_id = character(), execution_route = character(),
        final_package_version = character(), final_package_commit = character()
      )
    } else {
      data.table::data.table(
        method = task$method, seed = task$seed, cohort = task$cohort,
        fold = rep(1:2, each = 2L), n_train = 2L, n_test = 2L,
        sample_idx = seq_len(4L),
        sample_id = paste0(task$label, "-s", seq_len(4L)),
        group_id = paste0(task$label, "-g", seq_len(4L)),
        y = rep(c(0L, 1L), 2L), score_m = c(0.1, 0.9, 0.2, 0.8),
        score_b = c(0.2, 0.8, 0.3, 0.7), package_commit = package_commit,
        analysis_code_id = within_code_id,
        execution_route = "pre-extension-ineligible",
        final_package_version = package_version,
        final_package_commit = package_commit
      )
    }
    embedded <- NULL
    if (witnessed) {
      embedded <- list(
        schema_version = witness$schema_version,
        witness_sha256 = witness_sha, producer = witness$producer,
        producer_script_sha256 = witness$producer_script_sha256,
        source_bundle_sha256 = witness$source_bundle_sha256,
        task_table_sha256 = witness$task_table_sha256,
        cache_sha256 = witness$cache_sha256,
        cells = data.table::copy(witness_cells[
          method == task$method & cohort == task$cohort & seed == task$seed
        ]),
        predictions = data.table::copy(witness_predictions[
          method == task$method & cohort == task$cohort & seed == task$seed
        ]),
        fold_audit = data.table::copy(witness_audit[
          method == task$method & cohort == task$cohort & seed == task$seed
        ])
      )
    }
    bundle <- list(
      cells = cell, predictions = prediction,
      seeds = task$seed, methods = task$method, cohorts = task$cohort,
      package_version = package_version, package_commit = package_commit,
      cache_package_commit = cache_commit, analysis_code_id = within_code_id
    )
    if (witnessed) bundle$structural_ineligibility_witness <- embedded
    shard_dir <- file.path(shard_root, task$label)
    dir.create(shard_dir)
    shard_path <- file.path(shard_dir, "within_bundle.rds")
    saveRDS(bundle, shard_path, version = 3L)
    within_predictions[[i]] <- prediction
  }
  for (i in seq_along(methods)) {
    method_value <- methods[[i]]
    tasks <- within_tasks[method == method_value]
    paths <- file.path(shard_root, tasks$label, "within_bundle.rds")
    source_sha <- if (method_value == "lrt-deepmaha") {
      deepmaha_source_sha
    } else {
      digest::digest(paste0(method_value, "-source"), algo = "sha256")
    }
    manifest <- tasks[, .(
      label, method, seed, cohort,
      path = file.path(label, "within_bundle.rds"),
      bytes = unname(file.info(paths)$size),
      sha256 = vapply(paths, env$.overlay_sha256, character(1L)),
      source_bundle_sha256 = source_sha,
      task_table_sha256 = within_task_sha,
      package_version = package_version, package_commit = package_commit,
      cache_package_commit = cache_commit, analysis_code_id = within_code_id
    )]
    if (method_value == "lrt-deepmaha") {
      manifest[, structural_witness_sha256 := ifelse(
        env$.overlay_key_string(manifest, c("method", "cohort", "seed")) %in%
          witness_key_strings, witness_sha, NA_character_
      )]
    }
    manifest_path <- file.path(
      shard_root, paste0("output_manifest__", method_value, ".tsv")
    )
    data.table::fwrite(manifest, manifest_path, sep = "\t", quote = FALSE,
                       na = "")
    manifest_rows[[i]] <- manifest
  }

  transfer_tasks <- data.table::CJ(
    method = transfer_methods, seed = seeds, sorted = TRUE
  )
  transfer_tasks[, label := sprintf("%s__seed_%03d", method, seed)]
  transfer_tasks_path <- file.path(td, "transfer_tasks.tsv")
  data.table::fwrite(transfer_tasks, transfer_tasks_path, sep = "\t",
                     col.names = FALSE, quote = FALSE)
  transfer_task_sha <- env$.overlay_sha256(transfer_tasks_path)
  directions <- data.table::data.table(
    disease = sprintf("D%02d", seq_len(18L)),
    held_out = sprintf("H%02d", seq_len(18L)),
    n_training_cohorts = c(rep(2L, 12L), rep(1L, 6L))
  )
  transfer_overlay_cells <- vector("list", nrow(transfer_tasks))
  transfer_predictions <- vector("list", nrow(transfer_tasks))
  for (i in seq_len(nrow(transfer_tasks))) {
    task <- transfer_tasks[i]
    cell <- data.table::copy(directions)[, `:=`(
      method = task$method, seed = task$seed,
      eligible = n_training_cohorts >= 2L,
      ineligible_reason = ifelse(n_training_cohorts >= 2L, NA_character_,
                                 "one training cohort is structural"),
      n_test = 1L, n_group_split_violations = 0L,
      engagement_metric = if (task$method == "cvae") {
        "not_applicable_class_conditional"
      } else if (task$method == "moe-gated") {
        "not_applicable_generic_mixture_n_experts"
      } else "n_environments",
      engagement_ok = TRUE,
      engagement_value = if (task$method == "cvae") NA_real_ else 2,
      package_version = package_version,
      package_commit = package_commit, cache_package_commit = cache_commit,
      analysis_code_id = transfer_code_id,
      execution_route = "pre-extension-ineligible",
      final_package_version = package_version,
      final_package_commit = package_commit
    )]
    prediction <- cell[eligible == TRUE, .(
      method, seed, disease, held_out,
      sample_id = paste(disease, held_out, "sample", sep = "-"),
      group_id = paste(disease, held_out, "group", sep = "-"),
      y = 1L, score_m = 0.8, score_b = 0.7,
      package_commit, analysis_code_id,
      execution_route, final_package_version, final_package_commit
    )]
    bundle <- list(
      cells = cell, predictions = prediction,
      package_version = package_version, package_commit = package_commit,
      cache_package_commit = cache_commit, analysis_code_id = transfer_code_id
    )
    shard_dir <- file.path(transfer_root, task$label)
    dir.create(shard_dir)
    saveRDS(bundle, file.path(shard_dir, "transfer_bundle.rds"), version = 3L)
    transfer_overlay_cells[[i]] <- cell
    transfer_predictions[[i]] <- prediction
  }

  other_within_methods <- sprintf("base-within-%02d", seq_len(37L))
  base_within_cells <- data.table::CJ(
    method = c(methods, other_within_methods), seed = seeds,
    cohort = cohorts, sorted = TRUE
  )[, `:=`(
    eligible = FALSE, ineligible_reason = "pre-extension ineligible",
    n_eff = NA_integer_, n_valid_folds = 0L, outer_k = 2L,
    n_group_split_violations = 0L, package_version = package_version,
    package_commit = package_commit, cache_package_commit = cache_commit,
    analysis_code_id = within_code_id,
    execution_route = "pre-extension-ineligible",
    final_package_version = package_version,
    final_package_commit = package_commit
  )]
  within_prediction_template <- data.table::rbindlist(within_predictions,
                                                       fill = TRUE)[0L]
  preserved_within <- base_within_cells[
    method == "base-within-01" & seed == 101L & cohort == cohorts[[1L]]
  ]
  base_within_cells[
    method == preserved_within$method & seed == preserved_within$seed &
      cohort == preserved_within$cohort,
    `:=`(eligible = TRUE, ineligible_reason = NA_character_, n_eff = 4L,
         n_valid_folds = 2L)
  ]
  preserved_within_predictions <- data.table::data.table(
    method = preserved_within$method, seed = preserved_within$seed,
    cohort = preserved_within$cohort, fold = rep(1:2, each = 2L),
    n_train = 2L, n_test = 2L, sample_idx = seq_len(4L),
    sample_id = paste0("preserved-within-s", seq_len(4L)),
    group_id = paste0("preserved-within-g", seq_len(4L)),
    y = rep(c(0L, 1L), 2L), score_m = c(0.1, 0.9, 0.2, 0.8),
    score_b = c(0.2, 0.8, 0.3, 0.7), package_commit = package_commit,
    analysis_code_id = within_code_id,
    execution_route = "pre-extension-ineligible",
    final_package_version = package_version,
    final_package_commit = package_commit
  )
  base_within <- list(cells = base_within_cells,
                      predictions = data.table::rbindlist(list(
                        within_prediction_template, preserved_within_predictions
                      ), fill = TRUE))
  base_within_path <- file.path(td, "base_within.rds")
  saveRDS(base_within, base_within_path, version = 3L)

  other_transfer_methods <- sprintf("base-transfer-%02d", seq_len(17L))
  base_transfer_cells <- data.table::CJ(
    method = c(transfer_methods, other_transfer_methods), seed = seeds,
    direction = seq_len(18L), sorted = TRUE
  )
  base_transfer_cells <- merge(base_transfer_cells, directions[
    , direction := .I], by = "direction", sort = FALSE
  )[, `:=`(
    direction = NULL, eligible = FALSE,
    ineligible_reason = "pre-extension ineligible", n_test = 1L,
    n_group_split_violations = 0L, engagement_metric = NA_character_,
    engagement_ok = NA, engagement_value = NA_real_,
    package_version = package_version, package_commit = package_commit,
    cache_package_commit = cache_commit, analysis_code_id = transfer_code_id,
    execution_route = "pre-extension-ineligible",
    final_package_version = package_version,
    final_package_commit = package_commit
  )]
  transfer_prediction_template <- data.table::rbindlist(
    transfer_predictions, fill = TRUE
  )[0L]
  preserved_transfer <- base_transfer_cells[
    method == "base-transfer-01" & seed == 101L & disease == "D01"
  ]
  base_transfer_cells[
    method == preserved_transfer$method & seed == preserved_transfer$seed &
      disease == preserved_transfer$disease & held_out == preserved_transfer$held_out,
    `:=`(eligible = TRUE, ineligible_reason = NA_character_)
  ]
  preserved_transfer_prediction <- data.table::data.table(
    method = preserved_transfer$method, seed = preserved_transfer$seed,
    disease = preserved_transfer$disease, held_out = preserved_transfer$held_out,
    sample_id = "preserved-transfer-s1", group_id = "preserved-transfer-g1",
    y = 1L, score_m = 0.8, score_b = 0.7,
    package_commit = package_commit, analysis_code_id = transfer_code_id,
    execution_route = "pre-extension-ineligible",
    final_package_version = package_version,
    final_package_commit = package_commit
  )
  base_transfer <- list(cells = base_transfer_cells,
                        predictions = data.table::rbindlist(list(
                          transfer_prediction_template,
                          preserved_transfer_prediction
                        ), fill = TRUE))
  base_transfer_path <- file.path(td, "base_transfer.rds")
  saveRDS(base_transfer, base_transfer_path, version = 3L)

  runtime_asset_dir <- file.path(td, "runtime_assets")
  dir.create(runtime_asset_dir)
  runtime_asset_paths <- file.path(runtime_asset_dir,
                                   c("runtime.simg", "model.ckpt"))
  writeLines("runtime image fixture", runtime_asset_paths[[1L]])
  writeLines("runtime model fixture", runtime_asset_paths[[2L]])
  writeLines("2026-07-19T00:00:00Z",
             file.path(within_receipt, "captured_utc.txt"))
  writeLines("s5", file.path(within_receipt, "hostname.txt"))
  writeLines(package_commit, file.path(within_receipt, "package_commit.txt"))
  file.create(file.path(within_receipt, "package_git_status.txt"))
  writeLines("fixture explicit environment",
             file.path(within_receipt, "python_conda_explicit.txt"))
  writeLines("fixture conda list",
             file.path(within_receipt, "python_conda_list.txt"))
  writeLines("fixture==1.0", file.path(within_receipt, "python_freeze.txt"))
  writeLines(c(
    "python=3.12.12", "numpy=2.4.2", "pandas=2.3.3",
    "sklearn=1.8.0", "scipy=1.16.3", "torch=2.11.0",
    "torchvision=0.26.0", "torch_cuda_available=False",
    "torch_cuda_version=None"
  ), file.path(within_receipt, "python_runtime.txt"))
  data.table::fwrite(data.table::data.table(
    Package = c("OmicSelector", "data.table"),
    Version = c(package_version, "fixture"), LibPath = "/fixture"
  ), file.path(within_receipt, "r_package_manifest.tsv"), sep = "\t")
  writeLines("R version 4.5.2 (fixture)",
             file.path(within_receipt, "r_runtime.txt"))
  data.table::fwrite(data.table::data.table(
    file = "fixture", bytes = 1L
  ), file.path(within_receipt, "receipt_files.tsv"), sep = "\t")
  runtime_sha <- vapply(runtime_asset_paths, env$.overlay_sha256, character(1L))
  writeLines(sprintf("%s  %s", runtime_sha, runtime_asset_paths),
             file.path(within_receipt, "runtime_asset_sha256.txt"))
  closure_names <- c(
    "captured_utc.txt", "hostname.txt", "package_commit.txt",
    "package_git_status.txt", "python_conda_explicit.txt",
    "python_conda_list.txt", "python_freeze.txt", "python_runtime.txt",
    "r_package_manifest.tsv", "r_runtime.txt", "receipt_files.tsv",
    "runtime_asset_sha256.txt"
  )
  closure_sha <- vapply(file.path(within_receipt, closure_names),
                        env$.overlay_sha256, character(1L))
  writeLines(sprintf("%s  /fixture/%s", closure_sha, closure_names),
             file.path(within_receipt, "receipt_sha256.txt"))

  argos_modules <- c("numpy", "pandas", "sklearn", "torch", "torchvision",
                     "kymatio", "pyod", "tabpfn", "tabicl", "tabdpt")
  data.table::fwrite(data.table::as.data.table(list(
    key = c("reticulate", "StableMate", "StableMate_RemoteSha", "python",
            "python_version", paste0("python_module_", argos_modules)),
    value = c(
      "1.45.0", "0.1.0", "cf0fe9f344c1756b6f89fd33f9fdae1a7f35b24a",
      "/home/kgs24/omicselector_stage/python/bin/python", "Python 3.12.12",
      rep("fixture", length(argos_modules))
    )
  )), file.path(transfer_receipt, "environment_manifest.tsv"), sep = "\t")
  model_dir <- file.path(transfer_receipt, "public_models")
  dir.create(model_dir)
  model_path <- file.path(model_dir, "fixture.ckpt")
  writeLines("fixture model payload", model_path)
  data.table::fwrite(data.table::data.table(
    path = "/home/kgs24/omicselector_stage/public_models/fixture.ckpt",
    bytes = file.info(model_path)$size,
    sha256 = env$.overlay_sha256(model_path)
  ), file.path(transfer_receipt, "public_model_manifest.tsv"), sep = "\t")
  writeLines(c(
    "tabpfn==8.0.7", "tabicl==2.1.1", "tabdpt==1.2.0",
    sprintf("fixture%d==1.0", seq_len(7L))
  ), file.path(transfer_receipt, "python_freeze.txt"))
  data.table::fwrite(data.table::data.table(
    Package = c("reticulate", "StableMate"),
    Version = c("1.45.0", "0.1.0"), LibPath = "/fixture"
  ), file.path(transfer_receipt, "r_package_manifest.tsv"), sep = "\t")

  within_receipt_sha <- env$.overlay_receipt_digest(
    env$.overlay_receipt_files(within_receipt, "within_s5")
  )
  transfer_receipt_sha <- env$.overlay_receipt_digest(
    env$.overlay_receipt_files(transfer_receipt, "transfer_argos")
  )
  base_within_sha <- env$.overlay_sha256(base_within_path)
  base_transfer_sha <- env$.overlay_sha256(base_transfer_path)

  output_dir <- file.path(td, "candidate")
  args <- c(
    paste0("--base-within-bundle=", base_within_path),
    paste0("--expected-base-within-bundle-sha256=", base_within_sha),
    paste0("--base-transfer-bundle=", base_transfer_path),
    paste0("--expected-base-transfer-bundle-sha256=", base_transfer_sha),
    paste0("--within-tasks=", within_tasks_path),
    paste0("--expected-within-task-sha256=", within_task_sha),
    paste0("--within-shard-root=", shard_root),
    paste0("--within-split-manifest-dir=", shard_root),
    paste0("--transfer-tasks=", transfer_tasks_path),
    paste0("--expected-transfer-task-sha256=", transfer_task_sha),
    paste0("--transfer-shard-root=", transfer_root),
    paste0("--within-environment-receipt=", within_receipt),
    paste0("--expected-within-environment-receipt-sha256=",
           within_receipt_sha),
    paste0("--transfer-environment-receipt=", transfer_receipt),
    paste0("--expected-transfer-environment-receipt-sha256=",
           transfer_receipt_sha),
    paste0("--structural-witness=", witness_path),
    paste0("--expected-structural-witness-sha256=", witness_sha),
    paste0("--expected-structural-cache-sha256=", structural_cache_sha),
    paste0("--expected-structural-fold-engine-sha256=", fold_engine_sha),
    paste0("--expected-execution-package-version=", package_version),
    paste0("--expected-execution-package-commit=", package_commit),
    paste0("--expected-cache-package-commit=", cache_commit),
    paste0("--expected-within-code-id=", within_code_id),
    paste0("--expected-transfer-code-id=", transfer_code_id),
    "--expected-producer-package-version=2.6.5.9000",
    paste0("--expected-producer-package-commit=", producer_commit),
    paste0("--output-dir=", output_dir), "--run-mode=full"
  )
  list(
    root = td, shard_root = shard_root, transfer_root = transfer_root,
    within_receipt = within_receipt, transfer_receipt = transfer_receipt,
    within_tasks = within_tasks, within_tasks_path = within_tasks_path,
    transfer_tasks_path = transfer_tasks_path,
    transfer_task_sha = transfer_task_sha,
    within_task_sha = within_task_sha, witness = witness,
    witness_path = witness_path, witness_sha = witness_sha,
    output_dir = output_dir, args = args, package_commit = package_commit,
    package_version = package_version, cache_commit = cache_commit,
    within_code_id = within_code_id, transfer_code_id = transfer_code_id,
    structural_cache_sha = structural_cache_sha,
    fold_engine_sha = fold_engine_sha, base_within_sha = base_within_sha,
    base_transfer_sha = base_transfer_sha,
    within_receipt_sha = within_receipt_sha,
    transfer_receipt_sha = transfer_receipt_sha,
    producer_commit = producer_commit, runtime_asset_paths = runtime_asset_paths,
    transfer_model_path = model_path,
    producer_root = producer_root, base_within_path = base_within_path,
    base_transfer_path = base_transfer_path
  )
}

test_that("dependency overlay candidate is atomic, complete, and fail closed", {
  env <- overlay_candidate_script_env()
  fixture <- overlay_candidate_fixture(env)
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  pins <- list(
    package_version = fixture$package_version,
    package_commit = fixture$package_commit,
    cache_commit = fixture$cache_commit,
    code_id = fixture$within_code_id,
    task_sha = fixture$within_task_sha,
    structural_cache_sha = fixture$structural_cache_sha,
    fold_engine_sha = fixture$fold_engine_sha
  )

  manifest_path <- file.path(
    fixture$shard_root, "output_manifest__ai-scarf.tsv"
  )
  missing_path <- paste0(manifest_path, ".missing")
  file.rename(manifest_path, missing_path)
  expect_error(
    env$.overlay_validate_split_manifests(
      fixture$shard_root, fixture$shard_root, fixture$within_tasks, pins
    ), "exactly 13"
  )
  file.rename(missing_path, manifest_path)
  file.copy(manifest_path,
            file.path(fixture$shard_root, "output_manifest__extra.tsv"))
  expect_error(
    env$.overlay_validate_split_manifests(
      fixture$shard_root, fixture$shard_root, fixture$within_tasks, pins
    ), "exactly 13"
  )
  unlink(file.path(fixture$shard_root, "output_manifest__extra.tsv"))

  manifest <- data.table::fread(manifest_path)
  shard_path <- file.path(fixture$shard_root, manifest$path[[1L]])
  original <- readBin(shard_path, "raw", n = file.info(shard_path)$size)
  altered <- readRDS(shard_path)
  altered$cells$ineligible_reason <- "altered after manifest"
  saveRDS(altered, shard_path, version = 3L)
  expect_error(
    env$.overlay_validate_split_manifests(
      fixture$shard_root, fixture$shard_root, fixture$within_tasks, pins
    ), "changed after"
  )
  writeBin(original, shard_path)

  validated_manifests <- env$.overlay_validate_split_manifests(
    fixture$shard_root, fixture$shard_root, fixture$within_tasks, pins
  )
  within_data <- env$.overlay_read_within_shards(validated_manifests, pins)
  bad_audit_witness <- unserialize(serialize(fixture$witness, NULL))
  bad_audit_witness$fold_audit$fit_status[[1L]] <- "error: fixture"
  expect_error(
    env$.overlay_validate_reassembled_within(
      within_data, fixture$within_tasks, validated_manifests,
      bad_audit_witness, fixture$witness_sha, pins
    ), "contains an error"
  )
  leaking_within <- within_data
  leaking_within$predictions <- data.table::copy(within_data$predictions)
  leak_key <- leaking_within$predictions[method == "ai-scarf"][1L,
    .(method, cohort, seed)]
  leak_rows <- leaking_within$predictions[
    method == leak_key$method & cohort == leak_key$cohort &
      seed == leak_key$seed, which = TRUE
  ]
  leaking_within$predictions$group_id[leak_rows[c(1L, 3L)]] <- "leaking-group"
  leaking_within$predictions$y[leak_rows[c(1L, 3L)]] <- 0L
  expect_error(
    env$.overlay_validate_reassembled_within(
      leaking_within, fixture$within_tasks, validated_manifests,
      fixture$witness, fixture$witness_sha, pins
    ), "crosses held-out folds"
  )

  transfer_tasks <- env$.overlay_read_transfer_tasks(
    fixture$transfer_tasks_path, fixture$transfer_task_sha
  )
  transfer_pins <- list(
    package_version = fixture$package_version,
    package_commit = fixture$package_commit,
    cache_commit = fixture$cache_commit, code_id = fixture$transfer_code_id,
    task_sha = fixture$transfer_task_sha
  )
  transfer_data <- env$.overlay_read_transfer_shards(
    fixture$transfer_root, transfer_tasks, transfer_pins
  )
  bad_engagement <- data.table::copy(transfer_data$cells)
  bad_engagement[method == "dann" & eligible == TRUE,
                 engagement_value := 1]
  expect_error(
    env$.overlay_validate_transfer_science(
      bad_engagement, transfer_data$predictions
    ), "engagement evidence"
  )
  bad_transfer_predictions <- data.table::rbindlist(list(
    transfer_data$predictions,
    data.table::copy(transfer_data$predictions[1L])[
      , `:=`(sample_id = "conflicting-sample", y = 1L - y)
    ]
  ))
  expect_error(
    env$.overlay_validate_transfer_science(
      transfer_data$cells, bad_transfer_predictions
    ), "inconsistent outcome labels"
  )
  off_grid_predictions <- data.table::rbindlist(list(
    transfer_data$predictions,
    data.table::copy(transfer_data$predictions[1L])[
      , `:=`(disease = "OFF_GRID", held_out = "OFF_GRID",
             sample_id = "off-grid-sample", group_id = "off-grid-group")
    ]
  ))
  expect_error(
    env$.overlay_validate_transfer_science(
      transfer_data$cells, off_grid_predictions
    ), "exactly equal eligible cell keys"
  )
  reused_group_predictions <- data.table::copy(transfer_data$predictions)
  reused_group <- reused_group_predictions[held_out == "H01",
                                             unique(group_id)][[1L]]
  reused_group_predictions[held_out == "H02", group_id := reused_group]
  expect_error(
    env$.overlay_validate_transfer_science(
      transfer_data$cells, reused_group_predictions
    ), "more than one held-out cohort"
  )
  substituted_predictions <- data.table::copy(transfer_data$predictions)
  substituted_predictions[
    method == "cvae" & seed == 101L & disease == "D01" & held_out == "H01",
    sample_id := "method-specific-substitution"
  ]
  expect_error(
    env$.overlay_validate_transfer_science(
      transfer_data$cells, substituted_predictions
    ), "different held-out sample/group/outcome surfaces"
  )
  expect_error(
    env$.overlay_validate_predictions(
      transfer_data$cells, off_grid_predictions,
      c("disease", "held_out", "seed", "method"), "sample_id", "n_test",
      "Transfer"
    ), "exactly equal eligible cell keys"
  )

  expect_error(
    env$.overlay_validate_environment_receipt(
      fixture$transfer_receipt, "within_s5", fixture$package_commit,
      fixture$package_version, fixture$transfer_receipt_sha
    ), "Within-s5 receipt"
  )

  runtime_original <- readBin(
    fixture$runtime_asset_paths[[1L]], "raw",
    n = file.info(fixture$runtime_asset_paths[[1L]])$size
  )
  writeLines("altered runtime", fixture$runtime_asset_paths[[1L]])
  expect_error(
    env$.overlay_validate_environment_receipt(
      fixture$within_receipt, "within_s5", fixture$package_commit,
      fixture$package_version, fixture$within_receipt_sha
    ), "runtime asset bytes"
  )
  writeBin(runtime_original, fixture$runtime_asset_paths[[1L]])

  dirty_path <- file.path(fixture$producer_root, "untracked.txt")
  writeLines("dirty", dirty_path)
  dirty_args <- fixture$args
  dirty_args[grepl("^--output-dir=", dirty_args)] <- paste0(
    "--output-dir=", file.path(fixture$root, "dirty_producer_candidate")
  )
  expect_error(env$.overlay_candidate_main(dirty_args), "exact clean expected")
  unlink(dirty_path)

  stale_base_args <- fixture$args
  stale_base_args[grepl(
    "^--expected-base-within-bundle-sha256=", stale_base_args
  )] <- paste0("--expected-base-within-bundle-sha256=", paste0(rep("0", 64L),
                                                               collapse = ""))
  stale_base_output <- file.path(fixture$root, "stale_base_candidate")
  stale_base_args[grepl("^--output-dir=", stale_base_args)] <-
    paste0("--output-dir=", stale_base_output)
  expect_error(env$.overlay_candidate_main(stale_base_args),
               "Base within/transfer bundle bytes")
  expect_false(file.exists(stale_base_output))
  expect_error(
    env$.overlay_validate_environment_receipt(
      fixture$within_receipt, "transfer_argos", fixture$package_commit,
      fixture$package_version, fixture$within_receipt_sha
    ), "Transfer-Argos receipt|regular immutable files|route-specific SHA-256"
  )
  missing_model_path <- paste0(fixture$transfer_model_path, ".missing")
  expect_true(file.rename(fixture$transfer_model_path, missing_model_path))
  expect_error(
    env$.overlay_validate_environment_receipt(
      fixture$transfer_receipt, "transfer_argos", fixture$package_commit,
      fixture$package_version, fixture$transfer_receipt_sha
    ), "public-model payload"
  )
  expect_true(file.rename(missing_model_path, fixture$transfer_model_path))
  model_original <- readBin(
    fixture$transfer_model_path, "raw",
    n = file.info(fixture$transfer_model_path)$size
  )
  writeLines("tampered model payload", fixture$transfer_model_path)
  expect_error(
    env$.overlay_validate_environment_receipt(
      fixture$transfer_receipt, "transfer_argos", fixture$package_commit,
      fixture$package_version, fixture$transfer_receipt_sha
    ), "payload bytes do not match"
  )
  writeBin(model_original, fixture$transfer_model_path)

  bad_witness <- fixture$witness
  bad_witness$cells <- data.table::rbindlist(list(
    bad_witness$cells,
    bad_witness$cells[1L][, cohort := "GSE_UNAUTHORIZED"]
  ))
  bad_witness_path <- file.path(fixture$root, "bad_witness.rds")
  saveRDS(bad_witness, bad_witness_path, version = 3L)
  bad_args <- fixture$args
  bad_args[grepl("^--structural-witness=", bad_args)] <-
    paste0("--structural-witness=", bad_witness_path)
  bad_args[grepl("^--expected-structural-witness-sha256=", bad_args)] <-
    paste0("--expected-structural-witness-sha256=",
           env$.overlay_sha256(bad_witness_path))
  bad_output <- file.path(fixture$root, "failed_candidate")
  bad_args[grepl("^--output-dir=", bad_args)] <-
    paste0("--output-dir=", bad_output)
  expect_error(env$.overlay_candidate_main(bad_args),
               "exactly the nine|do not match source bundle")
  expect_false(file.exists(bad_output))

  input_paths <- c(
    sub("^--base-within-bundle=", "",
        fixture$args[grepl("^--base-within-bundle=", fixture$args)]),
    sub("^--base-transfer-bundle=", "",
        fixture$args[grepl("^--base-transfer-bundle=", fixture$args)]),
    fixture$within_tasks_path, fixture$shard_root,
    sub("^--transfer-tasks=", "",
        fixture$args[grepl("^--transfer-tasks=", fixture$args)]),
    fixture$transfer_root, fixture$within_receipt, fixture$transfer_receipt,
    fixture$witness_path, fixture$runtime_asset_paths
  )
  before <- env$.overlay_input_snapshot(normalizePath(input_paths))
  expect_invisible(env$.overlay_candidate_main(fixture$args))
  after <- env$.overlay_input_snapshot(normalizePath(input_paths))
  expect_identical(before$path, after$path)
  expect_identical(before$sha256, after$sha256)

  expected_top <- c(
    "within_bundle.rds", "within_cells.tsv", "within_predictions.tsv.gz",
    "transfer_bundle.rds", "transfer_cells.tsv", "transfer_predictions.tsv.gz",
    "acceptance_checks.tsv", "route_environment_map.tsv", "manifest.tsv"
  )
  expect_true(all(file.exists(file.path(fixture$output_dir, expected_top))))
  expect_false(dir.exists(file.path(
    fixture$output_dir, "environment", "transfer_argos", "public_models"
  )))
  checks <- data.table::fread(file.path(fixture$output_dir,
                                        "acceptance_checks.tsv"))
  expect_true(nrow(checks) > 0L && all(checks$pass))
  candidate_within <- data.table::fread(file.path(
    fixture$output_dir, "within_cells.tsv"
  ))
  candidate_transfer <- data.table::fread(file.path(
    fixture$output_dir, "transfer_cells.tsv"
  ))
  expect_equal(nrow(candidate_within), 8500L)
  expect_equal(nrow(candidate_transfer), 2160L)
  expect_true(all(candidate_within[
    method %in% env$.overlay_within_methods,
    execution_route == "2.6.5_dependency_backed_exact"
  ]))
  expect_true(all(candidate_transfer[
    method %in% env$.overlay_transfer_methods,
    execution_route == "2.6.5_dependency_backed_exact"
  ]))
  candidate_within_predictions <- data.table::as.data.table(utils::read.delim(
    gzfile(file.path(fixture$output_dir, "within_predictions.tsv.gz")),
    check.names = FALSE, stringsAsFactors = FALSE
  ))
  candidate_transfer_predictions <- data.table::as.data.table(utils::read.delim(
    gzfile(file.path(fixture$output_dir, "transfer_predictions.tsv.gz")),
    check.names = FALSE, stringsAsFactors = FALSE
  ))
  expect_true(all(candidate_within_predictions[
    method %in% env$.overlay_within_methods,
    execution_route == "2.6.5_dependency_backed_exact"
  ]))
  expect_true(all(candidate_transfer_predictions[
    method %in% env$.overlay_transfer_methods,
    execution_route == "2.6.5_dependency_backed_exact"
  ]))
  candidate_manifest <- data.table::fread(file.path(fixture$output_dir,
                                                     "manifest.tsv"))
  candidate_files <- setdiff(list.files(
    fixture$output_dir, recursive = TRUE, include.dirs = FALSE,
    all.files = TRUE, no.. = TRUE
  ), "manifest.tsv")
  expect_setequal(candidate_manifest$file, candidate_files)
  expect_true(all(!startsWith(candidate_manifest$file, "/")))
  expect_true(all(candidate_manifest$run_mode == "full"))
  expect_true(all(candidate_manifest$producer_package_commit ==
                    fixture$producer_commit))
  expect_true(all(candidate_manifest$sha256 == vapply(
    file.path(fixture$output_dir, candidate_manifest$file),
    env$.overlay_sha256, character(1L)
  )))
  route_map <- data.table::fread(file.path(
    fixture$output_dir, "route_environment_map.tsv"
  ))
  expect_identical(route_map$execution_route,
                   rep("2.6.5_dependency_backed_exact", 2L))
  expect_identical(route_map$environment_route,
                   c("within_s5", "transfer_argos"))
  source_inputs <- data.table::fread(file.path(
    fixture$output_dir, "source_input_manifest.tsv"
  ))
  expect_equal(source_inputs[input_role == "base_within_bundle", sha256],
               fixture$base_within_sha)
  expect_equal(source_inputs[input_role == "base_transfer_bundle", sha256],
               fixture$base_transfer_sha)
  expect_setequal(source_inputs[input_role == "within_runtime_asset", path],
                  fixture$runtime_asset_paths)
  expect_identical(source_inputs[input_role == "transfer_runtime_asset", path],
                   normalizePath(fixture$transfer_model_path))
  expect_identical(
    source_inputs[input_role == "transfer_runtime_asset", sha256],
    env$.overlay_sha256(fixture$transfer_model_path)
  )
  compatibility <- data.table::fread(file.path(
    fixture$output_dir, "legacy_compatibility_overlay_input_manifest.tsv"
  ))
  expect_equal(nrow(compatibility), 2251L)
  expect_equal(anyDuplicated(compatibility$path), 0L)
  expect_true(all(file.exists(compatibility$path)))
  expect_true(all(compatibility$sha256 == vapply(
    compatibility$path, env$.overlay_sha256, character(1L)
  )))
  expect_error(env$.overlay_candidate_main(fixture$args), "already exists")
})
