runtime_receipt_source_dir <- function() {
  path <- testthat::test_path("..", "..", "inst", "scripts")
  if (dir.exists(path)) return(path)
  installed <- system.file("scripts", package = "OmicSelector")
  if (!nzchar(installed) || !dir.exists(installed)) {
    stop("Could not locate package scripts for runtime-receipt tests.",
         call. = FALSE)
  }
  installed
}

runtime_receipt_git <- function(root, args) {
  out <- system2("git", c("-C", root, args), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status"))) stop(paste(out, collapse = "\n"))
  out
}

runtime_receipt_set_single_thread <- function(env) {
  variables <- env$paper3_runtime_required_thread_env()
  old <- Sys.getenv(names(variables), unset = NA_character_)
  old_dt <- data.table::getDTthreads(verbose = FALSE)
  old_rng <- RNGkind()
  do.call(Sys.setenv, as.list(variables))
  data.table::setDTthreads(1L)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  function() {
    missing <- names(old)[is.na(old)]
    if (length(missing)) Sys.unsetenv(missing)
    present <- old[!is.na(old)]
    if (length(present)) do.call(Sys.setenv, as.list(present))
    data.table::setDTthreads(old_dt)
    do.call(RNGkind, as.list(old_rng))
  }
}

runtime_receipt_fixture <- function() {
  root <- tempfile("paper3-runtime-receipt-")
  dir.create(root)
  checkout <- file.path(root, "producer_checkout")
  scripts <- file.path(checkout, "inst", "scripts")
  dir.create(scripts, recursive = TRUE)
  script_names <- c("paper3_build_runtime_receipt.R",
                    "paper3_runtime_receipt_common.R")
  stopifnot(all(file.copy(
    file.path(runtime_receipt_source_dir(), script_names),
    file.path(scripts, script_names)
  )))
  installed_path <- find.package(
    "OmicSelector", lib.loc = .libPaths(), quiet = TRUE
  )
  if (!nzchar(installed_path)) {
    stop("Runtime-receipt fixture requires an installed OmicSelector record.",
         call. = FALSE)
  }
  version <- unname(read.dcf(
    file.path(installed_path, "DESCRIPTION"), fields = "Version"
  )[1L, "Version"])
  writeLines(c(
    "Package: OmicSelector", paste0("Version: ", version),
    "Title: Runtime Receipt Test Checkout",
    "Description: Exact runtime receipt fixture.", "License: MIT"
  ), file.path(checkout, "DESCRIPTION"))
  runtime_receipt_git(checkout, c("init", "-q"))
  runtime_receipt_git(checkout,
                      c("config", "user.email", "fixture@example.org"))
  runtime_receipt_git(checkout, c("config", "user.name", "Fixture"))
  runtime_receipt_git(checkout, c("add", "."))
  runtime_receipt_git(checkout,
                      c("commit", "-q", "-m", "runtime-receipt-fixture"))
  commit <- trimws(runtime_receipt_git(checkout, c("rev-parse", "HEAD")))
  env <- new.env(parent = globalenv())
  source(file.path(scripts, "paper3_build_runtime_receipt.R"), local = env)
  restore_threads <- runtime_receipt_set_single_thread(env)
  image <- normalizePath(file.path(root, "runtime.simg"), mustWork = FALSE)
  writeLines("immutable runtime image fixture", image)
  image <- normalizePath(image, mustWork = TRUE)
  tree_sha <- env$paper3_runtime_installed_package_tree_sha256("OmicSelector")
  output <- file.path(root, "receipt")
  args <- c(
    paste0("--expected-package-version=", version),
    paste0("--expected-package-commit=", commit),
    paste0("--expected-installed-package-tree-sha256=", tree_sha),
    paste0("--runtime-image=", image),
    paste0("--expected-runtime-image-sha256=", env$.paper3rt_sha256(image)),
    paste0("--output-dir=", output)
  )
  list(
    root = root, checkout = checkout, scripts = scripts, env = env,
    restore_threads = restore_threads, image = image, output = output,
    version = version, commit = commit, tree_sha = tree_sha, args = args
  )
}

test_that("runtime receipt ignores a pkgload source checkout shadow", {
  env <- new.env(parent = globalenv())
  source(file.path(
    runtime_receipt_source_dir(), "paper3_runtime_receipt_common.R"
  ), local = env)
  installed_path <- find.package(
    "OmicSelector", lib.loc = .libPaths(), quiet = TRUE
  )
  skip_if(!nzchar(installed_path), "no installed OmicSelector record")
  source_or_loaded_path <- find.package("OmicSelector", quiet = TRUE)
  skip_if(identical(normalizePath(source_or_loaded_path),
                    normalizePath(installed_path)),
          "test process is not source-shadowed")

  db <- env$.paper3rt_active_installed_db()
  expect_true("OmicSelector" %in% rownames(db))
  expect_identical(
    normalizePath(file.path(db["OmicSelector", "LibPath"], "OmicSelector")),
    normalizePath(installed_path)
  )
  expect_identical(
    env$paper3_runtime_installed_package_tree_sha256("OmicSelector"),
    env$.paper3rt_package_tree_sha256(installed_path)
  )
})

test_that("runtime receipt is atomic, exact, and guards start/end runtime", {
  skip_if_not_installed("OmicSelector")
  skip_if_not_installed("metafor")
  fixture <- runtime_receipt_fixture()
  on.exit(fixture$restore_threads(), add = TRUE)
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  env <- fixture$env

  expect_error(
    env$.paper3rt_parse_producer_args(c(fixture$args, "--unknown=value")),
    "Unknown argument"
  )
  expect_error(
    env$.paper3rt_parse_producer_args(c(
      fixture$args,
      paste0("--expected-package-version=", fixture$version)
    )), "Duplicate argument"
  )

  dirty <- file.path(fixture$checkout, "untracked.txt")
  writeLines("dirty", dirty)
  dirty_args <- fixture$args
  dirty_args[grepl("^--output-dir=", dirty_args)] <- paste0(
    "--output-dir=", file.path(fixture$root, "dirty-output")
  )
  expect_error(env$.paper3rt_build_runtime_receipt(dirty_args),
               "exact expected clean Git checkout")
  unlink(dirty)

  dangling_output <- file.path(fixture$root, "dangling-output")
  expect_true(file.symlink("missing-runtime-receipt-target", dangling_output))
  dangling_args <- fixture$args
  dangling_args[grepl("^--output-dir=", dangling_args)] <- paste0(
    "--output-dir=", dangling_output
  )
  expect_error(env$.paper3rt_build_runtime_receipt(dangling_args),
               "already exists")
  expect_identical(Sys.readlink(dangling_output),
                   "missing-runtime-receipt-target")
  expect_false(file.exists(file.path(fixture$root,
                                     "missing-runtime-receipt-target")))

  result <- env$.paper3rt_build_runtime_receipt(fixture$args)
  expect_identical(result$output_dir, fixture$output)
  expect_true(grepl("^[0-9a-f]{64}$", result$manifest_sha256))
  expected_files <- c(
    "runtime_facts.tsv", "runtime_assets.tsv", "package_closure.tsv",
    "producer_provenance.tsv", "manifest.tsv"
  )
  expect_setequal(list.files(fixture$output, all.files = TRUE, no.. = TRUE),
                  expected_files)
  manifest <- data.table::fread(file.path(fixture$output, "manifest.tsv"),
                                colClasses = "character")
  expect_setequal(manifest$file, setdiff(expected_files, "manifest.tsv"))
  expect_true(all(manifest$sha256 == vapply(
    file.path(fixture$output, manifest$file), env$.paper3rt_sha256,
    character(1L)
  )))
  expect_true(all(as.numeric(manifest$bytes) == file.info(
    file.path(fixture$output, manifest$file)
  )$size))

  assets <- data.table::fread(file.path(fixture$output, "runtime_assets.tsv"),
                              colClasses = "character")
  expect_setequal(assets$asset, c(
    "runtime_image", "executable_Rscript", "executable_gzip",
    "executable_git"
  ))
  expect_true(all(startsWith(assets$path, "/")))
  expect_true(all(assets$path == vapply(
    assets$path, normalizePath, character(1L), mustWork = TRUE
  )))
  expect_true(all(assets$sha256 == vapply(
    assets$path, env$.paper3rt_sha256, character(1L)
  )))
  expect_false(file.exists(file.path(fixture$output, basename(fixture$image))))

  facts <- data.table::fread(file.path(fixture$output, "runtime_facts.tsv"),
                             colClasses = "character", na.strings = NULL)
  required_threads <- names(env$paper3_runtime_required_thread_env())
  expect_setequal(facts[section == "thread_env", key], required_threads)
  expect_true(all(facts[section == "thread_env", value] == "1"))
  expect_identical(facts[section == "data.table" & key == "threads", value],
                   "1")
  expect_equal(nrow(facts[section == "r" & key %in% c("version", "platform")]),
               2L)
  expect_equal(nrow(facts[section == "rngkind"]), 3L)
  expect_identical(
    facts[section == "rngkind"][match(
      c("kind", "normal_kind", "sample_kind"), key
    ), value],
    c("Mersenne-Twister", "Inversion", "Rejection")
  )
  expect_equal(nrow(facts[section == "numeric_runtime" &
                           key %in% c("BLAS", "LAPACK")]), 2L)

  packages <- data.table::fread(
    file.path(fixture$output, "package_closure.tsv"), colClasses = "character"
  )
  expect_setequal(packages[is_root == "TRUE", package],
                  c("OmicSelector", "data.table", "digest", "metafor"))
  expect_identical(packages[package == "OmicSelector", tree_sha256],
                   fixture$tree_sha)
  expect_true(all(packages$tree_sha256 == vapply(
    packages$root_path, env$.paper3rt_package_tree_sha256, character(1L)
  )))

  validate_args <- list(
    receipt_dir = fixture$output,
    expected_manifest_sha256 = result$manifest_sha256,
    expected_package_version = fixture$version,
    expected_package_commit = fixture$commit,
    expected_installed_package_tree_sha256 = fixture$tree_sha,
    expected_runtime_image_path = fixture$image,
    expected_runtime_image_sha256 = env$.paper3rt_sha256(fixture$image)
  )
  guard <- do.call(env$paper3_validate_runtime_receipt,
                   c(validate_args, list(phase = "start")))
  expect_s3_class(guard, "paper3_runtime_receipt_guard")
  expect_invisible(do.call(
    env$paper3_validate_runtime_receipt,
    c(validate_args, list(phase = "end", start_guard = guard))
  ))

  old_omp <- Sys.getenv("OMP_NUM_THREADS")
  Sys.setenv(OMP_NUM_THREADS = "2")
  expect_error(do.call(
    env$paper3_validate_runtime_receipt,
    c(validate_args, list(phase = "end", start_guard = guard))
  ), "not exactly 1")
  Sys.setenv(OMP_NUM_THREADS = old_omp)

  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  expect_error(do.call(
    env$paper3_validate_runtime_receipt,
    c(validate_args, list(phase = "end", start_guard = guard))
  ), "RNGkind must be exactly")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  image_original <- readBin(fixture$image, "raw",
                            n = file.info(fixture$image)$size)
  writeLines("tampered runtime image", fixture$image)
  expect_error(do.call(
    env$paper3_validate_runtime_receipt,
    c(validate_args, list(phase = "start"))
  ), "Runtime assets differs|Runtime asset")
  writeBin(image_original, fixture$image)

  helper_path <- file.path(fixture$scripts,
                           "paper3_runtime_receipt_common.R")
  helper_original <- readBin(helper_path, "raw", n = file.info(helper_path)$size)
  writeLines(c(readLines(helper_path, warn = FALSE), "# tampered"), helper_path)
  expect_error(do.call(
    env$paper3_validate_runtime_receipt,
    c(validate_args, list(phase = "start"))
  ), "helper_script_sha256")
  writeBin(helper_original, helper_path)

  wrong_tree <- validate_args
  wrong_tree$expected_installed_package_tree_sha256 <-
    paste0(rep("0", 64L), collapse = "")
  expect_error(do.call(
    env$paper3_validate_runtime_receipt,
    c(wrong_tree, list(phase = "start"))
  ), "installed_omicselector_tree_sha256")

  expect_error(env$.paper3rt_build_runtime_receipt(fixture$args),
               "already exists")
})
