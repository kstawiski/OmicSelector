ranktree_paper_script_env <- function(script) {
  env <- new.env(parent = globalenv())
  path <- testthat::test_path("..", "..", "inst", "scripts", script)
  if (!file.exists(path)) path <- system.file("scripts", script,
                                               package = "OmicSelector")
  if (!nzchar(path) || !file.exists(path)) stop("Missing paper script: ", script)
  sys.source(path, envir = env)
  env
}

ranktree_pin_fixture <- function() {
  list(
    package_version = "2.6.5.9000", package_commit = strrep("a", 40L),
    cache_commit = strrep("b", 40L), cache_sha256 = strrep("c", 64L),
    cache_manifest_sha256 = strrep("d", 64L),
    cache_digest_sha256 = strrep("e", 64L),
    runtime_image_sha256 = strrep("7", 64L),
    runtime_attestation_manifest_sha256 = strrep("8", 64L),
    container_launcher_sha256 = strrep("9", 64L),
    runner_script_sha256 = strrep("f", 64L),
    engine_manifest_sha256 = strrep("1", 64L),
    environment_manifest_sha256 = strrep("2", 64L),
    provenance_script_sha256 = strrep("3", 64L),
    provenance_manifest_sha256 = strrep("4", 64L),
    provenance_union_inventory_sha256 = strrep("6", 64L),
    analysis_plan_sha256 = strrep("0", 64L),
    r_version = "R version 4.5.2 (2025-10-31)",
    analysis_code_id = strrep("5", 16L)
  )
}

test_that("rank-tree runner accepts only registered unique method subsets", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  pin <- ranktree_pin_fixture()
  args <- c(
    "--package-root", tempdir(), "--paper-root", tempdir(),
    "--analysis-plan", tempfile(),
    "--cache", tempfile(), "--output-dir", tempfile(),
    "--provenance-manifest", tempfile(),
    "--provenance-union-inventory", tempfile(),
    "--methods", "reo-rankboost,ws-rclr-trimmed",
    "--runtime-image", tempfile(),
    "--runtime-attestation-manifest", tempfile(),
    "--engine-manifest", tempfile(), "--environment-manifest", tempfile(),
    "--expected-package-version", pin$package_version,
    "--expected-analysis-plan-sha256", pin$analysis_plan_sha256,
    "--expected-package-commit", pin$package_commit,
    "--expected-cache-commit", pin$cache_commit,
    "--expected-cache-sha256", pin$cache_sha256,
    "--expected-cache-manifest-sha256", pin$cache_manifest_sha256,
    "--expected-cache-digest-sha256", pin$cache_digest_sha256,
    "--expected-runtime-image-sha256", pin$runtime_image_sha256,
    "--expected-runtime-attestation-manifest-sha256",
    pin$runtime_attestation_manifest_sha256,
    "--expected-container-launcher-sha256",
    pin$container_launcher_sha256,
    "--expected-r-version", pin$r_version,
    "--expected-engine-manifest-sha256", pin$engine_manifest_sha256,
    "--expected-environment-manifest-sha256",
    pin$environment_manifest_sha256,
    "--expected-provenance-script-sha256", pin$provenance_script_sha256,
    "--expected-provenance-manifest-sha256", pin$provenance_manifest_sha256,
    "--expected-provenance-union-inventory-sha256",
    pin$provenance_union_inventory_sha256,
    "--expected-analysis-code-id", pin$analysis_code_id)
  parsed <- env$.ranktree_runner_args(args)
  expect_identical(parsed$methods,
                   c("ws-rclr-trimmed", "reo-rankboost"))
  bad <- args
  bad[match("--methods", bad) + 1L] <- "reo-rankboost,not-registered"
  expect_error(env$.ranktree_runner_args(bad), "seven registered")
  duplicate <- args
  duplicate[match("--methods", duplicate) + 1L] <-
    "reo-rankboost,reo-rankboost"
  expect_error(env$.ranktree_runner_args(duplicate), "unique")
  malformed_runtime <- args
  malformed_runtime[match("--expected-runtime-image-sha256",
                          malformed_runtime) + 1L] <- "not-a-sha"
  expect_error(env$.ranktree_runner_args(malformed_runtime),
               "runtime-image-sha256.*64")
})

test_that("assembler requires exact plan runtime and provenance pins", {
  env <- ranktree_paper_script_env(
    "paper3_assemble_ranktree_sensitivity.R")
  pin <- ranktree_pin_fixture()
  values <- c(
    "--package-root" = tempdir(), "--paper-root" = tempdir(),
    "--analysis-plan" = tempfile(),
    "--bundles" = paste(tempfile(rep("bundle", 7L)), collapse = ","),
    "--output-dir" = tempfile(), "--runtime-image" = tempfile(),
    "--runtime-attestation-manifest" = tempfile(),
    "--expected-package-version" = pin$package_version,
    "--expected-analysis-plan-sha256" = pin$analysis_plan_sha256,
    "--expected-package-commit" = pin$package_commit,
    "--expected-cache-commit" = pin$cache_commit,
    "--expected-cache-sha256" = pin$cache_sha256,
    "--expected-cache-manifest-sha256" = pin$cache_manifest_sha256,
    "--expected-cache-digest-sha256" = pin$cache_digest_sha256,
    "--expected-runtime-image-sha256" = pin$runtime_image_sha256,
    "--expected-runtime-attestation-manifest-sha256" =
      pin$runtime_attestation_manifest_sha256,
    "--expected-container-launcher-sha256" =
      pin$container_launcher_sha256,
    "--expected-r-version" = pin$r_version,
    "--expected-runner-script-sha256" = pin$runner_script_sha256,
    "--expected-engine-manifest-sha256" = pin$engine_manifest_sha256,
    "--expected-environment-manifest-sha256" =
      pin$environment_manifest_sha256,
    "--expected-provenance-script-sha256" = pin$provenance_script_sha256,
    "--expected-provenance-manifest-sha256" =
      pin$provenance_manifest_sha256,
    "--expected-provenance-union-inventory-sha256" =
      pin$provenance_union_inventory_sha256,
    "--expected-analysis-code-id" = pin$analysis_code_id
  )
  args <- unlist(Map(c, names(values), unname(values)), use.names = FALSE)
  parsed <- env$.rta_args(args)
  expect_equal(length(parsed$bundles), 7L)
  expect_identical(parsed$expected_analysis_plan_sha256,
                   pin$analysis_plan_sha256)
  plan_index <- match("--analysis-plan", args)
  expect_error(env$.rta_args(args[-c(plan_index, plan_index + 1L)]),
               "--analysis-plan")
})

test_that("runtime image bytes must match their exact launch pin", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  image <- tempfile("ranktree-runtime-", fileext = ".sif")
  on.exit(unlink(image), add = TRUE)
  writeBin(as.raw(c(1L, 2L, 3L, 4L)), image)
  sha <- env$.ranktree_sha(image)
  expect_identical(env$.ranktree_assert_file_pin(
    image, sha, "Singularity runtime image"), normalizePath(image))
  writeBin(as.raw(c(4L, 3L, 2L, 1L)), image)
  expect_error(env$.ranktree_assert_file_pin(
    image, sha, "Singularity runtime image"), "does not match")
})

test_that("approved analysis plan is canonical, pinned, and fail closed", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  root <- tempfile("ranktree-plan-root-")
  dir.create(file.path(root, "plan"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  plan <- file.path(
    root, "plan", "ranktree_external_competitor_sensitivity_20260719.md"
  )
  writeLines(c("# Plan", "Status: APPROVED", "Frozen design."), plan)
  sha <- env$.ranktree_sha(plan)
  expect_identical(env$.ranktree_validate_analysis_plan(plan, root, sha),
                   normalizePath(plan))
  writeLines(c("# Plan", "Status: DRAFT"), plan)
  expect_error(env$.ranktree_validate_analysis_plan(
    plan, root, env$.ranktree_sha(plan)), "Status: APPROVED")
  writeLines(c("# Plan", "Status: APPROVED", "mutated"), plan)
  expect_error(env$.ranktree_validate_analysis_plan(plan, root, sha),
               "does not match")
  other <- file.path(root, "plan", "other.md")
  writeLines(c("# Plan", "Status: APPROVED"), other)
  expect_error(env$.ranktree_validate_analysis_plan(
    other, root, env$.ranktree_sha(other)), "exact canonical path")
})

test_that("in-process runtime and library bytes must match attestation", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  pin <- ranktree_pin_fixture()
  td <- tempfile("ranktree-attestation-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  manifest <- data.table::data.table(
    scope = "r_library", owner = "ranktreeEnsemble", version = "0.24",
    root_path = "/pinned/library/ranktreeEnsemble",
    relative_path = "R/ranktreeEnsemble.rdb", bytes = "12",
    sha256 = strrep("a", 64L), r_version = pin$r_version,
    r_home = "/opt/R/4.5.2/lib/R", container_runtime = "singularity",
    container_marker = "/runtime/pinned.simg",
    runtime_image_sha256 = pin$runtime_image_sha256,
    container_launcher_sha256 = pin$container_launcher_sha256
  )
  path <- file.path(td, "runtime.tsv")
  data.table::fwrite(manifest, path, sep = "\t", quote = FALSE)
  current <- data.table::fread(path, colClasses = "character")
  expect_identical(env$.ranktree_validate_runtime_attestation(
    path, env$.ranktree_sha(path), pin$r_version, pin$runtime_image_sha256,
    pin$container_launcher_sha256, current = current)$path,
    normalizePath(path))
  changed <- data.table::copy(current)
  changed$sha256 <- strrep("b", 64L)
  expect_error(env$.ranktree_validate_runtime_attestation(
    path, env$.ranktree_sha(path), pin$r_version, pin$runtime_image_sha256,
    pin$container_launcher_sha256, current = changed),
    "runtime/library bytes differ")
  expect_error(env$.ranktree_validate_runtime_attestation(
    path, env$.ranktree_sha(path), "wrong R", pin$runtime_image_sha256,
    pin$container_launcher_sha256, current = current), "R version differs")
})

test_that("package launcher forwards clean-container attestation via env", {
  runner_env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  launcher <- testthat::test_path(
    "..", "..", "inst", "scripts",
    "paper3_run_ranktree_sensitivity_container.sh"
  )
  if (!file.exists(launcher)) {
    launcher <- system.file(
      "scripts", "paper3_run_ranktree_sensitivity_container.sh",
      package = "OmicSelector"
    )
  }
  expect_true(nzchar(launcher) && file.exists(launcher))
  td <- tempfile("ranktree-launcher-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  image <- file.path(td, "runtime.simg")
  writeBin(as.raw(1:8), image)
  capture <- file.path(td, "capture.txt")
  fake <- file.path(td, "singularity")
  writeLines(c(
    "#!/usr/bin/env bash", "set -euo pipefail",
    "printf '%s\\n' \"${SINGULARITYENV_OMICSELECTOR_RUNTIME_IMAGE_SHA256-}\" \"${SINGULARITYENV_OMICSELECTOR_CONTAINER_LAUNCHER_SHA256-}\" \"$@\" > \"$RANKTREE_LAUNCH_CAPTURE\""
  ), fake)
  Sys.chmod(fake, "0755")
  image_sha <- runner_env$.ranktree_sha(image)
  launcher_sha <- runner_env$.ranktree_sha(launcher)
  old_path <- Sys.getenv("PATH")
  status <- system2(
    "bash", c(launcher, "--entrypoint", "capture-attestation",
              "--capture-runtime-attestation", file.path(td, "attest.tsv"),
              "--runtime-image", image,
              "--expected-runtime-image-sha256", image_sha,
              "--expected-container-launcher-sha256", launcher_sha),
    env = c(paste0("PATH=", td, .Platform$path.sep, old_path),
            paste0("RANKTREE_LAUNCH_CAPTURE=", capture))
  )
  expect_equal(status, 0L)
  lines <- readLines(capture)
  expect_identical(lines[1:2], c(image_sha, launcher_sha))
  expect_true("--cleanenv" %in% lines)
  expect_false("--env" %in% lines)
})

test_that("provenance preflight uses canonical manifest and union from any cwd", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  root <- tempfile("ranktree-provenance-")
  dir.create(file.path(root, "code"), recursive = TRUE)
  dir.create(file.path(root, "data_bucket"))
  dir.create(file.path(root, "knowledge"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  script <- file.path(root, "code", "check_provenance.R")
  manifest <- file.path(root, "data_bucket", "_PROVENANCE_MANIFEST.md")
  union <- file.path(root, "knowledge", "union_inventory.tsv")
  writeLines("# checker", script)
  writeLines("# manifest", manifest)
  data.table::fwrite(data.table::data.table(accession = "GSE1"), union,
                     sep = "\t")
  cohorts <- list(u1 = list(
    source_accessions = "GSE1", expr_per_sample = matrix(1, 2, 1),
    group_id = c("g1", "g2"), provenance_block = "independent"
  ))
  log <- file.path(root, "preflight.log")
  fake <- function(cohorts, paper_root, log_path) {
    expect_identical(getwd(), normalizePath(root))
    writeLines("ok", log_path)
    resolution <- sub("\\.log$", "_resolution.tsv", log_path)
    writeLines("accession_a\taccession_b", resolution)
    list(raw_status = 0L, n_overlap_pairs = 0L,
         resolution_path = resolution)
  }
  old <- getwd()
  setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  out <- env$.ranktree_run_bound_provenance(
    cohorts, root, script, manifest, union, env$.ranktree_sha(script),
    env$.ranktree_sha(manifest), env$.ranktree_sha(union), log, fake
  )
  expect_equal(out$raw_status, 0L)
  mutating <- function(cohorts, paper_root, log_path) {
    write("changed", file = manifest, append = TRUE)
    fake(cohorts, paper_root, log_path)
  }
  original_manifest <- env$.ranktree_sha(manifest)
  expect_error(env$.ranktree_run_bound_provenance(
    cohorts, root, script, manifest, union, env$.ranktree_sha(script),
    original_manifest, env$.ranktree_sha(union), log, mutating), "mutated")
})

test_that("frozen engine manifest verifies all four source bytes", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  td <- tempfile("ranktree-engine-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  roles <- c("benchmark_engine", "paired_auc_engine", "fold_engine",
             "public_design")
  paths <- file.path(td, paste0(roles, ".R"))
  for (i in seq_along(paths)) writeLines(paste("source", i), paths[[i]])
  manifest <- data.table::data.table(
    role = rev(roles), path = rev(paths),
    sha256 = rev(vapply(paths, env$.ranktree_sha, character(1L))))
  manifest_path <- file.path(td, "engine.tsv")
  data.table::fwrite(manifest, manifest_path, sep = "\t")
  out <- env$.ranktree_validate_engine_manifest(
    manifest_path, env$.ranktree_sha(manifest_path))
  expect_identical(out$role, roles)
  writeLines("changed", paths[[1L]])
  expect_error(env$.ranktree_validate_engine_manifest(
    manifest_path, env$.ranktree_sha(manifest_path)), "does not match")
})

test_that("dependency source manifest is exact over recursive closure", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  td <- tempfile("ranktree-deps-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  paths <- file.path(td, c("ranktree.tar.gz", "dependency.tar.gz"))
  writeLines("ranktree source", paths[[1L]])
  writeLines("dependency source", paths[[2L]])
  closure <- data.table::data.table(
    package = c("ranktreeEnsemble", "treeDependency"),
    version = c("0.24", "1.2.3"))
  manifest <- data.table::data.table(
    package = closure$package, version = closure$version, source_path = paths,
    source_sha256 = vapply(paths, env$.ranktree_sha, character(1L)))
  path <- file.path(td, "environment.tsv")
  data.table::fwrite(manifest, path, sep = "\t")
  expect_equal(nrow(env$.ranktree_validate_environment_manifest(
    path, env$.ranktree_sha(path), closure = closure)), 2L)
  extra <- data.table::rbindlist(list(
    closure, data.table::data.table(package = "undeclared", version = "1")))
  expect_error(env$.ranktree_validate_environment_manifest(
    path, env$.ranktree_sha(path), closure = extra), "exact recursive")
  writeLines("mutated", paths[[2L]])
  expect_error(env$.ranktree_validate_environment_manifest(
    path, env$.ranktree_sha(path), closure = closure), "does not match")
})

test_that("runner rejects every non-reviewed ineligibility", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  eligible <- data.table::data.table(
    method = "reo-rankforest", cohort = "u1", seed = 101L,
    eligible = TRUE, ineligible_reason = NA_character_, n_eff = 20L,
    n_valid_folds = 5L, outer_k = 5L, n_group_split_violations = 0L,
    structural_exception_proven = FALSE,
    structural_exception_proof_sha256 = NA_character_)
  expect_invisible(env$.ranktree_validate_cell_result(eligible))
  structural <- data.table::copy(eligible)
  structural[, `:=`(
    cohort = "GSE83977", eligible = FALSE,
    ineligible_reason = env$.ranktree_structural_reason(), n_eff = NA_integer_,
    n_valid_folds = 0L, outer_k = 3L,
    structural_exception_proven = TRUE,
    structural_exception_proof_sha256 = strrep("a", 64L))]
  expect_invisible(env$.ranktree_validate_cell_result(structural))
  failure <- data.table::copy(structural)
  failure$ineligible_reason <- "fit failed: dependency unavailable"
  expect_error(env$.ranktree_validate_cell_result(failure),
               "Unreviewed ineligibility")
  wrong_unit <- data.table::copy(structural)
  wrong_unit$cohort <- "GSE_OTHER"
  expect_error(env$.ranktree_validate_cell_result(wrong_unit),
               "Unreviewed ineligibility")
})

test_that("GSE83977 actual execution path makes engine errors typed and fatal", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  X <- matrix(seq_len(36), nrow = 12,
              dimnames = list(paste0("s", 1:12), paste0("f", 1:3)))
  y <- rep(c(0L, 1L), 6L)
  cohort <- list(expr_per_sample = X, y_bin = y,
                 group_id = paste0("g", 1:12), min_per_fold = 2L)
  roster <- data.frame(method_id = "reo-rankforest", fit_fn = "fit_fake",
                       stringsAsFactors = FALSE)
  folds <- list(c(1L, 2L, 7L, 8L), c(3L, 4L, 9L, 10L),
                c(5L, 6L, 11L, 12L))
  grouped <- function(y, group_id, k, seed) folds
  baseline_ok <- function(X_train, y_train, X_test, panel_size) {
    rep(0, nrow(X_test))
  }
  fit_ok <- function(fit_fn, X_train, y_train, meta_train, hp) list(ok = TRUE)
  score_ok <- function(method_id, model, X, meta, roster) {
    if (nrow(X) == 8L) seq_len(nrow(X)) else rep(1, nrow(X))
  }
  direction_ok <- function(score, y) 1
  run <- function(baseline = baseline_ok, fit = fit_ok, score = score_ok) {
    td <- tempfile("ranktree-strict-path-")
    dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    strict <- function(method, cohort_id, cohort, roster, out_dir, seed,
                       outer_k) {
      env$.ranktree_strict_exception_cell(
        method, cohort_id, cohort, roster, out_dir, seed, outer_k,
        grouped_fold_fn = grouped,
        stratified_fold_fn = function(...) stop("unexpected"),
        baseline_fn = baseline, fit_dispatch = fit,
        score_dispatch = score, direction_fn = direction_ok
      )
    }
    env$.ranktree_call_benchmark_cell(
      function(...) stop("generic benchmark must not run for GSE83977"),
      "reo-rankforest", "GSE83977", cohort, roster, td, 101L, 3L,
      strict_exception_fn = strict
    )
  }
  result <- run()
  expect_false(result$eligible)
  expect_true(result$structural_exception_proven)
  expect_match(result$structural_exception_proof_sha256,
               "^[0-9a-f]{64}$")
  baseline_error <- function(...) stop("injected baseline")
  expect_error(run(baseline = baseline_error), class = "ranktree_engine_error")
  expect_error(run(fit = function(...) stop("injected fit")),
               class = "ranktree_engine_error")
  calls <- 0L
  training_error <- function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("injected training prediction")
    score_ok(...)
  }
  expect_error(run(score = training_error), class = "ranktree_engine_error")
  calls <- 0L
  heldout_error <- function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("injected heldout prediction")
    score_ok(...)
  }
  error <- tryCatch(run(score = heldout_error), error = identity)
  expect_s3_class(error, "ranktree_engine_error")
  expect_identical(error$stage, "heldout_prediction")
})

test_that("split ids are shared and fail on any held-out mismatch", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  base <- data.table::data.table(
    cohort = "u1", seed = 101L, fold = rep(1:2, each = 4L),
    n_train = 4L, n_test = 4L, sample_id = paste0("s", 1:8),
    group_id = paste0("g", 1:8), y = rep(c(0L, 0L, 1L, 1L), 2L),
    score_m = seq_len(8) / 10)
  predictions <- data.table::rbindlist(list(
    data.table::copy(base)[, method := "ws-rclr-trimmed"],
    data.table::copy(base)[, method := "reo-rankforest"]))
  out <- env$.ranktree_add_split_ids(predictions)
  expect_equal(data.table::uniqueN(out$split_id), 1L)
  broken <- data.table::copy(predictions)
  broken[method == "reo-rankforest" & sample_id == "s1", fold := 2L]
  expect_error(env$.ranktree_add_split_ids(broken),
               "provenance group crosses|byte-identical")
})

test_that("assembler bundle validation binds package cache and provenance pins", {
  env <- ranktree_paper_script_env(
    "paper3_assemble_ranktree_sensitivity.R")
  pin <- ranktree_pin_fixture()
  units <- data.table::data.table(
    unit_id = "u1", disease = "cancer", modality = "array",
    biospecimen = "serum", provenance_block = "independent",
    source_accessions = "GSE1", outer_k = 2L, primary_unit = TRUE,
    n_deposited_profiles = 4L, n_unique_groups = 4L)
  cells <- data.table::data.table(
    method = "reo-rankforest", cohort = "u1", seed = 101L,
    eligible = TRUE, ineligible_reason = NA_character_, n_eff = 4L,
    n_valid_folds = 2L, outer_k = 2L, grouping_policy = "groupKFold",
    n_group_split_violations = 0L,
    structural_exception_proven = FALSE,
    structural_exception_basis = NA_character_,
    structural_exception_proof_sha256 = NA_character_,
    package_version = pin$package_version,
    package_commit = pin$package_commit,
    cache_package_commit = pin$cache_commit,
    cache_sha256 = pin$cache_sha256,
    engine_manifest_sha256 = pin$engine_manifest_sha256,
    environment_manifest_sha256 = pin$environment_manifest_sha256,
    runtime_image_sha256 = pin$runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      pin$runtime_attestation_manifest_sha256,
    container_launcher_sha256 = pin$container_launcher_sha256,
    r_version = pin$r_version,
    provenance_manifest_sha256 = pin$provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      pin$provenance_union_inventory_sha256,
    analysis_plan_sha256 = pin$analysis_plan_sha256,
    analysis_code_id = pin$analysis_code_id)
  predictions <- data.table::data.table(
    method = "reo-rankforest", cohort = "u1", seed = 101L,
    fold = rep(1:2, each = 2L), split_id = "fixed", n_train = 2L,
    n_test = 2L, sample_id = paste0("s", 1:4), group_id = paste0("g", 1:4),
    y = c(0L, 1L, 0L, 1L), score_m = c(0.1, 0.9, 0.2, 0.8),
    package_version = pin$package_version,
    package_commit = pin$package_commit,
    cache_package_commit = pin$cache_commit,
    runtime_image_sha256 = pin$runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      pin$runtime_attestation_manifest_sha256,
    container_launcher_sha256 = pin$container_launcher_sha256,
    r_version = pin$r_version,
    provenance_manifest_sha256 = pin$provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      pin$provenance_union_inventory_sha256,
    analysis_plan_sha256 = pin$analysis_plan_sha256,
    analysis_code_id = pin$analysis_code_id)
  split_reference <- predictions[, .(
    fold = as.integer(fold), sample_id = as.character(sample_id),
    group_id = as.character(group_id), y = as.integer(y),
    n_train = as.integer(n_train), n_test = as.integer(n_test))]
  data.table::setorder(split_reference, fold, sample_id, group_id, y)
  predictions[, split_id := digest::digest(
    split_reference, algo = "sha256", serialize = TRUE)]
  bundle <- c(list(
    cells = cells, predictions = predictions, units = units,
    methods = "reo-rankforest",
    method_registry = data.table::data.table(
      method_id = "reo-rankforest", estimand = "within"),
    seeds = 101L, cohorts = "u1", r_version = pin$r_version,
    ranktree_dependency_closure = data.frame(
      package = "ranktreeEnsemble", version = "0.24",
      source_sha256 = strrep("6", 64L))), list(
        package_version = pin$package_version,
        package_commit = pin$package_commit,
        cache_package_commit = pin$cache_commit,
        cache_sha256 = pin$cache_sha256,
        cache_manifest_sha256 = pin$cache_manifest_sha256,
        cache_digest_sha256 = pin$cache_digest_sha256,
        runner_script_sha256 = pin$runner_script_sha256,
        engine_manifest_sha256 = pin$engine_manifest_sha256,
        environment_manifest_sha256 = pin$environment_manifest_sha256,
        runtime_image_sha256 = pin$runtime_image_sha256,
        runtime_attestation_manifest_sha256 =
          pin$runtime_attestation_manifest_sha256,
        container_launcher_sha256 = pin$container_launcher_sha256,
        provenance_script_sha256 = pin$provenance_script_sha256,
        provenance_manifest_sha256 = pin$provenance_manifest_sha256,
        provenance_union_inventory_sha256 =
          pin$provenance_union_inventory_sha256,
        analysis_plan_sha256 = pin$analysis_plan_sha256,
        analysis_code_id = pin$analysis_code_id))
  expect_invisible(env$.rta_validate_bundle(
    bundle, "reo-rankforest", pin, seeds = 101L))
  broken <- bundle
  broken$provenance_manifest_sha256 <- strrep("9", 64L)
  expect_error(env$.rta_validate_bundle(
    broken, "reo-rankforest", pin, seeds = 101L), "exact expected pin")
  broken <- bundle
  broken$predictions$cache_package_commit <- strrep("9", 40L)
  expect_error(env$.rta_validate_bundle(
    broken, "reo-rankforest", pin, seeds = 101L), "Prediction table")
  broken <- bundle
  broken$runtime_image_sha256 <- strrep("8", 64L)
  expect_error(env$.rta_validate_bundle(
    broken, "reo-rankforest", pin, seeds = 101L),
    "runtime_image_sha256.*exact expected pin")
  broken <- bundle
  broken$analysis_plan_sha256 <- strrep("b", 64L)
  expect_error(env$.rta_validate_bundle(
    broken, "reo-rankforest", pin, seeds = 101L),
    "analysis_plan_sha256.*exact expected pin")
})

test_that("seven-method family freezes 21 registry-ordered comparisons", {
  env <- ranktree_paper_script_env(
    "paper3_assemble_ranktree_sensitivity.R")
  ids <- env$.rta_method_ids()
  roster <- data.frame(method_id = ids, estimand = "within",
                       stringsAsFactors = FALSE)
  pairs <- OmicSelector::singlesample_within_method_pairs(roster)
  expect_equal(nrow(pairs), 21L)
  expect_identical(pairs$pair_id,
                   paste(pairs$method_a, pairs$method_b, sep = "::"))
  expect_true(all(pairs$roster_order_a < pairs$roster_order_b))
  expect_identical(pairs$method_a, ids[pairs$roster_order_a])
  expect_identical(pairs$method_b, ids[pairs$roster_order_b])
})

test_that("meta summaries keep BY denominators at 21 and 42", {
  skip_if_not_installed("metafor")
  env <- ranktree_paper_script_env(
    "paper3_assemble_ranktree_sensitivity.R")
  ids <- env$.rta_method_ids()
  roster <- data.frame(method_id = ids, estimand = "within",
                       stringsAsFactors = FALSE)
  pairs <- data.table::as.data.table(
    OmicSelector::singlesample_within_method_pairs(roster))
  units <- data.table::data.table(
    unit_id = paste0("u", 1:6), primary_unit = TRUE,
    n_deposited_profiles = 20L, n_unique_groups = 20L)
  pairwise <- data.table::CJ(pair_id = pairs$pair_id, unit_id = units$unit_id)
  pairwise <- merge(pairwise, pairs[, .(pair_id, method_a, method_b)],
                    by = "pair_id", sort = FALSE)
  pairwise[, `:=`(
    common_support = TRUE,
    estimate = 0.01 * match(unit_id, units$unit_id) +
      0.0001 * match(pair_id, pairs$pair_id),
    se = 0.02 + 0.001 * match(unit_id, units$unit_id))]
  out <- env$.rta_meta_all(pairwise, pairs, units)
  expect_equal(nrow(out), 42L)
  expect_true(all(out$zero_family_n == 21L))
  expect_true(all(out$relevance_family_n == 42L))
  expect_setequal(out$pair_id, pairs$pair_id)
  expect_setequal(out$analysis_population,
                  c("primary_circulating_biofluid",
                    "all_biofluids_sensitivity"))
})

test_that("complete-support outputs require at least five methods and units", {
  env <- ranktree_paper_script_env(
    "paper3_assemble_ranktree_sensitivity.R")
  ids <- env$.rta_method_ids()
  roster <- data.table::data.table(
    method_id = ids, estimand = "within", family = "K", role = "discriminator",
    tier = "R1", dep_route = "test")
  units <- data.table::data.table(
    unit_id = paste0("u", 1:6), primary_unit = TRUE,
    n_deposited_profiles = 10L, n_unique_groups = 10L)
  eligibility <- data.table::CJ(method_id = ids, unit_id = units$unit_id)
  eligibility[, `:=`(eligible = TRUE, eligible_seeds = 5L,
                     expected_seeds = 5L, reason = "eligible")]
  performance <- data.table::copy(eligibility)[, `:=`(
    auc = 0.6, n_strata = 25L, primary_unit = TRUE,
    n_deposited_profiles = 10L, n_unique_groups = 10L)]
  out <- env$.rta_complete_support(eligibility, performance, roster, units)
  summary <- out$table[item_type == "summary"]
  expect_true(all(!summary$degenerate))
  expect_true(all(summary$n_methods == 7L))
  expect_true(all(summary$n_units == 6L))
  expect_equal(nrow(out$rectangle), 2L * 3L * 7L * 6L)
  expect_true(all(out$rectangle$rank_allowed))
  expect_true(all(out$rectangle$conclusion_allowed))
})

test_that("degenerate complete support suppresses ranks and conclusions", {
  env <- ranktree_paper_script_env(
    "paper3_assemble_ranktree_sensitivity.R")
  ids <- env$.rta_method_ids()[1:4]
  roster <- data.table::data.table(
    method_id = ids, estimand = "within", family = "K", role = "discriminator",
    tier = "R1", dep_route = "test")
  units <- data.table::data.table(
    unit_id = paste0("u", 1:6), primary_unit = TRUE,
    n_deposited_profiles = 10L, n_unique_groups = 10L)
  eligibility <- data.table::CJ(method_id = ids, unit_id = units$unit_id)
  eligibility[, `:=`(eligible = TRUE, eligible_seeds = 5L,
                     expected_seeds = 5L, reason = "eligible")]
  performance <- data.table::copy(eligibility)[, `:=`(
    auc = 0.55 + 0.01 * match(method_id, ids), n_strata = 25L,
    primary_unit = TRUE, n_deposited_profiles = 10L,
    n_unique_groups = 10L)]
  out <- env$.rta_complete_support(eligibility, performance, roster, units)
  summary <- out$table[item_type == "summary"]
  expect_true(all(summary$degenerate))
  expect_true(all(!summary$rank_allowed))
  expect_true(all(!summary$conclusion_allowed))
  expect_true(nrow(out$rectangle) > 0L)
  expect_true(all(is.na(out$rectangle$within_unit_rank)))
  expect_true(all(!out$rectangle$rank_allowed))
  expect_true(all(!out$rectangle$conclusion_allowed))
})

ranktree_checkpoint_fixture <- function(env, output_dir,
                                        methods = "reo-rankforest") {
  pin <- ranktree_pin_fixture()
  pins <- list(
    version = pin$package_version, commit = pin$package_commit,
    cache_commit = pin$cache_commit, cache_sha256 = pin$cache_sha256,
    cache_manifest_sha256 = pin$cache_manifest_sha256,
    cache_digest_sha256 = pin$cache_digest_sha256,
    engine_manifest_sha256 = pin$engine_manifest_sha256,
    environment_manifest_sha256 = pin$environment_manifest_sha256,
    runtime_image_sha256 = pin$runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      pin$runtime_attestation_manifest_sha256,
    container_launcher_sha256 = pin$container_launcher_sha256,
    r_version = pin$r_version,
    analysis_plan_sha256 = pin$analysis_plan_sha256,
    provenance_script_sha256 = pin$provenance_script_sha256,
    provenance_manifest_sha256 = pin$provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      pin$provenance_union_inventory_sha256,
    analysis_code_id = pin$analysis_code_id)
  units <- data.table::data.table(unit_id = "u1")
  list(
    pins = pins,
    contract = env$.ranktree_checkpoint_contract(
      output_dir, methods, units, pins, pin$runner_script_sha256,
      pin$r_version))
}

test_that("runtime image pin propagates into the checkpoint contract", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  output <- file.path(tempdir(), "ranktree-runtime-contract")
  fixture <- ranktree_checkpoint_fixture(env, output)
  expect_identical(
    fixture$contract[key == "runtime_image_sha256", value],
    fixture$pins$runtime_image_sha256)
  expect_identical(
    fixture$contract[key == "runtime_attestation_manifest_sha256", value],
    fixture$pins$runtime_attestation_manifest_sha256)
  expect_identical(
    fixture$contract[key == "analysis_plan_sha256", value],
    fixture$pins$analysis_plan_sha256)
  expect_identical(
    fixture$contract[key == "provenance_union_inventory_sha256", value],
    fixture$pins$provenance_union_inventory_sha256)
  changed <- data.table::copy(fixture$contract)
  changed[key == "runtime_image_sha256", value := strrep("8", 64L)]
  expect_false(identical(changed, fixture$contract))
})

test_that("persistent checkpoint cells are reused with overwrite disabled", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  td <- tempfile("ranktree-resume-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  output <- file.path(td, "method-bundle")
  fixture <- ranktree_checkpoint_fixture(env, output)
  checkpoint <- env$.ranktree_open_checkpoints(output, fixture$contract)
  seed_dir <- file.path(checkpoint, "seed_101")
  dir.create(seed_dir)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  fake_benchmark <- function(method_id, cohort_id, cohort, roster, out_dir,
                             seed, outer_k, min_valid_folds,
                             min_pos_per_fold, min_neg_per_fold,
                             baseline_mode, panel_size, overwrite,
                             rerun_present, orient, hp) {
    expect_false(overwrite)
    expect_true(rerun_present)
    expect_identical(baseline_mode, "rclr_panel")
    path <- file.path(out_dir, paste0(cohort_id, "__", method_id, ".rds"))
    if (!overwrite && file.exists(path)) return(readRDS(path))
    calls$n <- calls$n + 1L
    result <- list(
      method_id = method_id, cohort = cohort_id, seed = seed,
      eligible = TRUE, ineligible_reason = NA_character_, n_eff = 4L,
      n_valid_folds = outer_k, n_group_split_violations = 0L,
      grouping_policy = "groupKFold",
      fold_scores = data.frame(fold = c(1L, 1L, 2L, 2L)))
    saveRDS(result, path)
    result
  }
  cohort <- list(min_per_fold = 1L)
  first <- env$.ranktree_call_benchmark_cell(
    fake_benchmark, "reo-rankforest", "u1", cohort, data.frame(), seed_dir,
    101L, 2L)
  second <- env$.ranktree_call_benchmark_cell(
    fake_benchmark, "reo-rankforest", "u1", cohort, data.frame(), seed_dir,
    101L, 2L)
  expect_equal(calls$n, 1L)
  expect_identical(second, first)
  for (result in list(first, second)) {
    result$outer_k <- 2L
    cell <- env$.ranktree_cell_row(result, 101L, fixture$pins)
    expect_invisible(env$.ranktree_validate_cell_result(cell))
  }
})

test_that("persistent checkpoint contract rejects mismatch and incompleteness", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  td <- tempfile("ranktree-contract-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  output <- file.path(td, "method-bundle")
  fixture <- ranktree_checkpoint_fixture(env, output)
  checkpoint <- env$.ranktree_open_checkpoints(output, fixture$contract)
  expect_true(dir.exists(checkpoint))
  mismatch <- ranktree_checkpoint_fixture(
    env, output, methods = "reo-rankboost")$contract
  expect_error(env$.ranktree_open_checkpoints(output, mismatch),
               "mismatched or incomplete")
  unlink(file.path(checkpoint, "checkpoint_contract.tsv"))
  expect_error(env$.ranktree_open_checkpoints(output, fixture$contract),
               "incomplete contract")
})

test_that("checkpoints persist on incomplete output and clear after promotion", {
  env <- ranktree_paper_script_env(
    "paper3_run_ranktree_sensitivity_bundle.R")
  td <- tempfile("ranktree-finalize-")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  output <- file.path(td, "method-bundle")
  fixture <- ranktree_checkpoint_fixture(env, output)
  checkpoint <- env$.ranktree_open_checkpoints(output, fixture$contract)
  dir.create(output)
  expect_error(
    env$.ranktree_remove_completed_checkpoints(output, fixture$contract),
    "lacks its completion manifest")
  expect_true(dir.exists(checkpoint))
  paths <- file.path(output, env$.ranktree_bundle_artifacts())
  for (i in seq_along(paths)) writeLines(paste("artifact", i), paths[[i]])
  manifest <- data.table::data.table(
    file = basename(paths), bytes = file.info(paths)$size,
    sha256 = vapply(paths, env$.ranktree_sha, character(1L)))
  data.table::fwrite(manifest, file.path(output, "output_manifest.tsv"),
                     sep = "\t")
  expect_invisible(env$.ranktree_remove_completed_checkpoints(
    output, fixture$contract))
  expect_false(dir.exists(checkpoint))
})
