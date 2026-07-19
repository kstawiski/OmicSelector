#!/usr/bin/env Rscript

# Build the immutable input snapshot for the OmicSelector paper's fair
# all-method comparison. This script never modifies its source benchmark files.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(OmicSelector)
})

.SNAPSHOT_SCRIPT_PATH <- local({
  if (exists(".PAPER3_ACTIVE_SCRIPT_PATH", inherits = TRUE)) {
    candidate <- get(".PAPER3_ACTIVE_SCRIPT_PATH", inherits = TRUE)
    if (length(candidate) == 1L && file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) && file.exists(sub("^--file=", "", file_arg[[1L]]))) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  }
  candidates <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) character() else as.character(frame$ofile)
  }), use.names = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Could not resolve fair-snapshot script path.")
  normalizePath(tail(candidates, 1L), mustWork = TRUE)
})
.SNAPSHOT_RUNTIME_HELPER <- file.path(
  dirname(.SNAPSHOT_SCRIPT_PATH), "paper3_runtime_receipt_common.R"
)
if (!file.exists(.SNAPSHOT_RUNTIME_HELPER)) {
  stop("Package-owned runtime-receipt helper is missing.")
}
sys.source(.SNAPSHOT_RUNTIME_HELPER, envir = environment())

.snapshot_read <- function(path) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    data.table::fread(cmd = paste("gzip -cd", shQuote(path)),
                      integer64 = "double")
  } else {
    data.table::fread(path, integer64 = "double")
  }
}

.snapshot_args <- function(args) {
  values <- list()
  allowed <- c(
    "dependency-overlay-candidate", "expected-overlay-manifest-sha256",
    "cohorts", "output-dir", "provenance-script", "provenance-root",
    "provenance-resolution", "upstream-validation",
    "upstream-validation-manifest", "expected-package-version",
    "expected-package-commit", "expected-installed-package-tree-sha256",
    "runtime-receipt", "expected-runtime-receipt-manifest-sha256",
    "runtime-image", "expected-runtime-image-sha256",
    "expected-cohort-manifest-sha256",
    "expected-provenance-script-sha256",
    "expected-provenance-manifest-sha256",
    "expected-provenance-union-inventory-sha256",
    "run-mode"
  )
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Arguments must use --name=value syntax: ", arg)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% allowed) stop("Unknown argument: --", key)
    if (key %in% names(values)) stop("Duplicate argument: --", key)
    values[[key]] <- value
  }
  required <- allowed
  missing <- setdiff(required, names(values))
  if (length(missing) || any(!nzchar(unlist(values[required])))) {
    stop("Missing required argument(s): --", paste(missing, collapse = ", --"))
  }
  if (!values[["run-mode"]] %in% c("full", "smoke")) {
    stop("--run-mode must be full or smoke.")
  }
  if (!grepl("^[0-9a-f]{64}$", values[["expected-overlay-manifest-sha256"]])) {
    stop("--expected-overlay-manifest-sha256 must be an exact SHA-256 pin.")
  }
  if (!grepl("^[0-9a-f]{40}$", values[["expected-package-commit"]])) {
    stop("--expected-package-commit must be an exact Git commit pin.")
  }
  if (!grepl("^[0-9a-f]{64}$",
             values[["expected-installed-package-tree-sha256"]])) {
    stop("--expected-installed-package-tree-sha256 must be an exact SHA-256 pin.")
  }
  for (key in c("expected-cohort-manifest-sha256",
                "expected-provenance-script-sha256",
                "expected-provenance-manifest-sha256",
                "expected-provenance-union-inventory-sha256",
                "expected-runtime-receipt-manifest-sha256",
                "expected-runtime-image-sha256")) {
    if (!grepl("^[0-9a-f]{64}$", values[[key]])) {
      stop("--", key, " must be an exact SHA-256 pin.")
    }
  }
  values
}

.snapshot_expected_units <- function() {
  c(
    "TORAY_PUBLIC_CANONICAL", "GSE109319", "GSE117063", "GSE122488",
    "GSE130654", "GSE137109", "GSE145176", "GSE186595", "GSE188232",
    "GSE188627", "GSE198834", "GSE200489", "GSE204951", "GSE210329",
    "GSE210546", "GSE221088", "GSE227778", "GSE259327_DIAGNOSTIC",
    "GSE266859", "GSE270497", "GSE279209", "GSE304572", "GSE31568",
    "GSE320031", "GSE41922", "GSE50013", "GSE64591", "GSE65071",
    "GSE70080", "GSE71008", "GSE83977", "GSE85589",
    "GSE94533_PUBLIC", "PRJEB30542"
  )
}

.snapshot_expected_k3 <- function() {
  c("GSE117063", "GSE130654", "GSE200489", "GSE204951", "GSE210546",
    "GSE227778", "GSE266859", "GSE83977")
}

.snapshot_code_provenance <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
  } else {
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
  }
  package_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                                mustWork = TRUE)
  git_bin <- if (file.exists("/usr/bin/git")) "/usr/bin/git" else Sys.which("git")
  if (!nzchar(git_bin)) stop("Could not locate a git executable.")
  commit <- system2(git_bin, c("-C", package_root, "rev-parse", "HEAD"),
                    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(commit, "status")) || length(commit) != 1L ||
      !grepl("^[0-9a-f]{40}$", commit)) {
    stop("Could not resolve the OmicSelector package git commit.")
  }
  status <- system2(
    git_bin, c("-C", package_root, "status", "--porcelain=v1",
             "--untracked-files=all"), stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(status, "status"))) {
    stop("Could not inspect the OmicSelector package worktree.")
  }
  list(
    script_path = script_path,
    script_sha256 = digest::digest(script_path, algo = "sha256", file = TRUE),
    package_root = package_root,
    package_git_commit = unname(commit[[1L]]),
    package_git_dirty = length(status) > 0L
  )
}

.snapshot_package_tree_sha256 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  files <- sort(list.files(root, recursive = TRUE, all.files = TRUE,
                           full.names = TRUE, include.dirs = FALSE,
                           no.. = TRUE))
  if (!length(files)) stop("Installed OmicSelector package tree is empty.")
  prefix <- paste0(root, .Platform$file.sep)
  relative <- substring(files, nchar(prefix) + 1L)
  if (any(!startsWith(files, prefix)) || any(grepl("[\r\n\t]", relative))) {
    stop("Installed OmicSelector package tree has an invalid path.")
  }
  rows <- paste(
    relative, file.info(files)$size,
    vapply(files, digest::digest, character(1L), algo = "sha256", file = TRUE),
    sep = "\t"
  )
  digest::digest(paste(rows, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

.snapshot_validate_upstream_receipt <- function(path) {
  receipt <- .snapshot_read(path)
  required_columns <- c("check_id", "pass")
  required_checks <- c(
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
  if (!all(required_columns %in% names(receipt)) ||
      anyDuplicated(receipt$check_id)) {
    stop("Upstream validation receipt has an invalid schema.")
  }
  if (!nrow(receipt) || anyNA(receipt$check_id) ||
      any(!nzchar(as.character(receipt$check_id))) ||
      !all(required_checks %in% as.character(receipt$check_id)) ||
      !all(tolower(as.character(receipt$pass)) %in% c("true", "1"))) {
    stop("Upstream validation receipt does not pass every required benchmark check.")
  }
  invisible(receipt)
}

.snapshot_validate_input_manifest <- function(path) {
  manifest <- .snapshot_read(path)
  path_column <- intersect(c("path", "file", "input_path"), names(manifest))
  hash_column <- intersect(c("sha256", "input_sha256", "observed_sha256"),
                           names(manifest))
  if (length(path_column) != 1L || length(hash_column) != 1L ||
      !nrow(manifest) || anyDuplicated(as.character(manifest[[path_column]]))) {
    stop("Upstream validation input manifest has an invalid schema.")
  }
  declared <- as.character(manifest[[path_column]])
  declared <- ifelse(
    grepl("^/", declared), declared,
    file.path(dirname(normalizePath(path, mustWork = TRUE)), declared)
  )
  if (anyNA(declared) || any(!file.exists(declared))) {
    stop("An upstream validation input is missing.")
  }
  observed <- vapply(declared, digest::digest, character(1L),
                     algo = "sha256", file = TRUE)
  expected <- tolower(as.character(manifest[[hash_column]]))
  if (any(!grepl("^[0-9a-f]{64}$", expected)) ||
      !identical(unname(observed), expected)) {
    stop("Upstream validation input bytes differ from their manifest.")
  }
  byte_column <- intersect(c("bytes", "size", "input_bytes"), names(manifest))
  if (length(byte_column) > 1L ||
      (length(byte_column) == 1L &&
       any(file.info(declared)$size != as.numeric(manifest[[byte_column]])))) {
    stop("Upstream validation input sizes differ from their manifest.")
  }
  invisible(manifest)
}

.snapshot_validate_provenance_sources <- function(
    root, expected_manifest_sha256, expected_union_sha256) {
  root <- normalizePath(root, mustWork = TRUE)
  paths <- c(
    manifest = file.path(root, "data_bucket", "_PROVENANCE_MANIFEST.md"),
    union = file.path(root, "knowledge", "union_inventory.tsv")
  )
  if (any(!file.exists(paths)) ||
      !identical(
        unname(vapply(paths, digest::digest, character(1L),
                      algo = "sha256", file = TRUE)),
        c(expected_manifest_sha256, expected_union_sha256)
      )) {
    stop("Canonical provenance manifest/union inventory is missing or unpinned.")
  }
  union <- .snapshot_read(paths[["union"]])
  if (!"accession" %in% names(union) || !nrow(union) ||
      anyNA(union$accession) || any(!nzchar(as.character(union$accession))) ||
      anyDuplicated(toupper(as.character(union$accession)))) {
    stop("Canonical provenance union inventory has an invalid accession key.")
  }
  paths
}

.snapshot_candidate_check_ids <- function() {
  c(
    "within_grid_13x34x5", "within_13_splitter_manifests",
    "within_manifest_shard_hashes", "within_task_key_equality",
    "transfer_exact_630_cells", "integrated_within_8500_cells",
    "integrated_transfer_2160_cells", "unrelated_base_rows_unchanged",
    "eligible_prediction_accounting", "deepmaha_exact_nine_witness_cells",
    "deepmaha_embedded_witness_equality",
    "route_specific_environment_receipts", "source_inputs_immutable",
    "producer_exact_clean_checkout", "base_bundles_exact_sha256",
    "durable_source_input_manifest"
  )
}

.snapshot_candidate_environment_files <- function() {
  c(
    file.path("environment", "within_s5", c(
      "captured_utc.txt", "hostname.txt", "package_commit.txt",
      "package_git_status.txt", "python_conda_explicit.txt",
      "python_conda_list.txt", "python_freeze.txt", "python_runtime.txt",
      "r_package_manifest.tsv", "r_runtime.txt", "receipt_files.tsv",
      "receipt_sha256.txt", "runtime_asset_sha256.txt"
    )),
    file.path("environment", "transfer_argos", c(
      "environment_manifest.tsv", "public_model_manifest.tsv",
      "python_freeze.txt", "r_package_manifest.tsv"
    ))
  )
}

.snapshot_multiset_key <- function(bytes, sha256) {
  sort(paste(as.numeric(bytes), as.character(sha256), sep = "\r"))
}

.snapshot_validate_candidate_source_manifest <- function(root, manifest) {
  path <- file.path(root, "source_input_manifest.tsv")
  source <- .snapshot_read(path)
  required_columns <- c("input_role", "path", "bytes", "sha256")
  expected_roles <- c(
    "base_transfer_bundle", "base_within_bundle", "structural_witness",
    "transfer_environment_receipt", "transfer_runtime_asset",
    "transfer_shard", "transfer_task_table", "within_environment_receipt",
    "within_runtime_asset", "within_shard", "within_splitter_manifest",
    "within_task_table"
  )
  exact_counts <- c(
    base_transfer_bundle = 1L, base_within_bundle = 1L,
    structural_witness = 1L, transfer_environment_receipt = 4L,
    transfer_runtime_asset = 3L, transfer_shard = 35L,
    transfer_task_table = 1L, within_environment_receipt = 13L,
    within_shard = 2210L, within_splitter_manifest = 13L,
    within_task_table = 1L
  )
  if (!identical(names(source), required_columns) || !nrow(source) ||
      anyNA(source) || anyDuplicated(source$path) ||
      !setequal(as.character(source$input_role), expected_roles) ||
      any(!grepl("^/", source$path)) ||
      any(!grepl("^[0-9a-f]{64}$", source$sha256))) {
    stop("Dependency-overlay candidate source-input manifest is malformed.")
  }
  role_counts <- source[, .N, by = input_role]
  role_counts <- setNames(role_counts$N, role_counts$input_role)
  if (any(role_counts[names(exact_counts)] != exact_counts) ||
      role_counts[["within_runtime_asset"]] < 1L) {
    stop("Dependency-overlay candidate source-input role closure is incomplete.")
  }
  declared <- as.character(source$path)
  if (any(!file.exists(declared)) || any(file.info(declared)$isdir %in% TRUE) ||
      !identical(declared, normalizePath(declared, mustWork = TRUE)) ||
      any(startsWith(declared, paste0(root, .Platform$file.sep)))) {
    stop("Dependency-overlay candidate source inputs are missing, mutable, or self-referential.")
  }
  observed_sha <- vapply(declared, digest::digest, character(1L),
                         algo = "sha256", file = TRUE)
  if (!identical(unname(observed_sha), as.character(source$sha256)) ||
      any(as.numeric(file.info(declared)$size) != as.numeric(source$bytes))) {
    stop("Dependency-overlay candidate source-input bytes differ from their manifest.")
  }

  unique_pin <- function(column, pattern = "^[0-9a-f]{64}$") {
    value <- unique(as.character(manifest[[column]]))
    if (length(value) != 1L || !grepl(pattern, value)) {
      stop("Dependency-overlay candidate has a non-scalar pin: ", column)
    }
    value
  }
  source_pin <- function(role) unique(as.character(
    source[input_role == role, sha256]
  ))
  pin_pairs <- list(
    base_within_bundle = c("base_within_bundle_sha256", "within_bundle.rds"),
    base_transfer_bundle = c("base_transfer_bundle_sha256", "transfer_bundle.rds"),
    within_task_table = c("within_task_sha256", "inputs/within_tasks.tsv"),
    transfer_task_table = c("transfer_task_sha256", "inputs/transfer_tasks.tsv"),
    structural_witness = c("structural_witness_sha256", "structural_witness.rds")
  )
  for (role in names(pin_pairs)) {
    declared_pin <- unique_pin(pin_pairs[[role]][[1L]])
    candidate_path <- file.path(root, pin_pairs[[role]][[2L]])
    copied_pin_required <- !role %in% c(
      "base_within_bundle", "base_transfer_bundle"
    )
    if (!identical(source_pin(role), declared_pin) ||
        (copied_pin_required &&
         !identical(digest::digest(candidate_path, algo = "sha256", file = TRUE),
                    declared_pin))) {
      stop("Dependency-overlay candidate source pin is not transitive: ", role)
    }
  }

  split_files <- manifest[role == "splitter_manifest", file]
  split_rows <- data.table::rbindlist(lapply(split_files, function(file) {
    x <- .snapshot_read(file.path(root, file))
    if (nrow(x) != 170L || !all(c("bytes", "sha256") %in% names(x))) {
      stop("A copied within splitter manifest has an invalid shard inventory.")
    }
    x[, .(bytes = as.numeric(bytes), sha256 = as.character(sha256))]
  }))
  within_source <- source[input_role == "within_shard"]
  if (!identical(
    .snapshot_multiset_key(split_rows$bytes, split_rows$sha256),
    .snapshot_multiset_key(within_source$bytes, within_source$sha256)
  )) {
    stop("Within splitter manifests do not transitively bind every source shard.")
  }
  split_source <- source[input_role == "within_splitter_manifest"]
  copied_split <- manifest[role == "splitter_manifest"]
  if (!identical(
    .snapshot_multiset_key(split_source$bytes, split_source$sha256),
    .snapshot_multiset_key(copied_split$bytes, copied_split$sha256)
  )) {
    stop("Copied within splitter manifests differ from their source bytes.")
  }
  transfer_manifest <- .snapshot_read(file.path(root, "transfer_shard_manifest.tsv"))
  transfer_source <- source[input_role == "transfer_shard"]
  if (nrow(transfer_manifest) != 35L ||
      !all(c("bytes", "sha256") %in% names(transfer_manifest)) ||
      !identical(
        .snapshot_multiset_key(transfer_manifest$bytes, transfer_manifest$sha256),
        .snapshot_multiset_key(transfer_source$bytes, transfer_source$sha256)
      )) {
    stop("Transfer shard manifest does not transitively bind every source shard.")
  }

  for (route in c("within_s5", "transfer_argos")) {
    role <- paste0(sub("_.*$", "", route), "_environment_receipt")
    copied <- manifest[
      startsWith(file, paste0("environment/", route, "/")), .(bytes, sha256)
    ]
    external <- source[input_role == role, .(bytes, sha256)]
    if (!identical(.snapshot_multiset_key(copied$bytes, copied$sha256),
                   .snapshot_multiset_key(external$bytes, external$sha256))) {
      stop("Copied route receipt differs from its source bytes: ", route)
    }
    receipt_rows <- manifest[
      startsWith(file, paste0("environment/", route, "/")),
      .(path = substring(file, nchar(paste0("environment/", route, "/")) + 1L),
        bytes, sha256)
    ]
    data.table::setorder(receipt_rows, path)
    observed_receipt_sha <- digest::digest(paste(
      receipt_rows$path, receipt_rows$bytes, receipt_rows$sha256,
      sep = "\t", collapse = "\n"
    ), algo = "sha256")
    pin_column <- if (identical(route, "within_s5")) {
      "within_environment_receipt_sha256"
    } else {
      "transfer_environment_receipt_sha256"
    }
    if (!identical(observed_receipt_sha, unique_pin(pin_column))) {
      stop("Route receipt digest differs from the candidate pin: ", route)
    }
  }

  model_manifest <- .snapshot_read(file.path(
    root, "environment", "transfer_argos", "public_model_manifest.tsv"
  ))
  transfer_assets <- source[input_role == "transfer_runtime_asset"]
  if (!all(c("bytes", "sha256") %in% names(model_manifest)) ||
      !identical(
        .snapshot_multiset_key(model_manifest$bytes, model_manifest$sha256),
        .snapshot_multiset_key(transfer_assets$bytes, transfer_assets$sha256)
      )) {
    stop("Transfer model manifest does not bind the source model payloads.")
  }
  within_asset_lines <- readLines(file.path(
    root, "environment", "within_s5", "runtime_asset_sha256.txt"
  ), warn = FALSE)
  matched <- regexec("^([0-9a-f]{64})[[:space:]]+(.+)$", within_asset_lines)
  parsed <- regmatches(within_asset_lines, matched)
  within_assets <- source[input_role == "within_runtime_asset"]
  if (!length(parsed) || any(lengths(parsed) != 3L) ||
      !identical(sort(vapply(parsed, `[[`, character(1L), 2L)),
                   sort(as.character(within_assets$sha256)))) {
    stop("Within runtime receipt does not bind every source runtime asset.")
  }
  invisible(source)
}

.snapshot_validate_overlay_candidate <- function(
    path, expected_sha256, expected_package_version,
    expected_package_commit, package_root, require_full = FALSE) {
  root <- normalizePath(path, mustWork = TRUE)
  manifest_path <- file.path(root, "manifest.tsv")
  if (!file.exists(manifest_path) ||
      !identical(digest::digest(manifest_path, algo = "sha256", file = TRUE),
                 expected_sha256)) {
    stop("Dependency-overlay candidate manifest differs from its exact pin.")
  }
  manifest <- .snapshot_read(manifest_path)
  required_columns <- c(
    "file", "role", "bytes", "sha256", "run_mode",
    "producer_script_sha256", "producer_package_version",
    "producer_package_commit", "base_within_bundle_sha256",
    "base_transfer_bundle_sha256", "within_environment_receipt_sha256",
    "transfer_environment_receipt_sha256", "execution_package_version",
    "execution_package_commit", "cache_package_commit",
    "within_analysis_code_id", "transfer_analysis_code_id",
    "within_task_sha256", "transfer_task_sha256",
    "structural_witness_sha256", "structural_cache_sha256",
    "structural_fold_engine_sha256"
  )
  required_files <- c(
    "within_bundle.rds", "within_cells.tsv", "within_predictions.tsv.gz",
    "transfer_bundle.rds", "transfer_cells.tsv",
    "transfer_predictions.tsv.gz", "acceptance_checks.tsv",
    "route_environment_map.tsv"
  )
  if (!all(required_columns %in% names(manifest)) || !nrow(manifest) ||
      anyDuplicated(manifest$file) ||
      !all(required_files %in% manifest$file) ||
      any(grepl("^/", manifest$file)) ||
      any(grepl("(^|/)\\.\\.(/|$)", manifest$file))) {
    stop("Dependency-overlay candidate manifest has an invalid schema or inventory.")
  }
  declared <- file.path(root, as.character(manifest$file))
  if (any(!file.exists(declared))) {
    stop("Dependency-overlay candidate is missing a manifest-listed artifact.")
  }
  present <- list.files(root, recursive = TRUE, all.files = TRUE,
                        include.dirs = FALSE, no.. = TRUE)
  if (!setequal(present, c(as.character(manifest$file), "manifest.tsv"))) {
    stop("Dependency-overlay candidate contains an unmanifested artifact.")
  }
  resolved <- normalizePath(declared, mustWork = TRUE)
  root_prefix <- paste0(root, .Platform$file.sep)
  if (any(!startsWith(resolved, root_prefix))) {
    stop("Dependency-overlay candidate manifest escapes its immutable root.")
  }
  hashes <- vapply(declared, digest::digest, character(1L),
                   algo = "sha256", file = TRUE)
  if (!identical(unname(hashes), as.character(manifest$sha256)) ||
      any(file.info(declared)$size != manifest$bytes)) {
    stop("Dependency-overlay candidate bytes differ from its manifest.")
  }
  producer_script <- file.path(
    package_root, "inst", "scripts",
    "paper3_build_dependency_overlay_candidate.R"
  )
  expected_splitters <- paste0(
    "splitter_manifests/output_manifest__",
    data.table::as.data.table(singlesample_method_roster())[
      estimand == "within" & tier == "R2", method_id
    ],
    ".tsv"
  )
  exact_roles <- c(
    integrated_within = 3L, integrated_transfer = 3L,
    splitter_manifest = 13L, structural_witness = 1L,
    frozen_task_table = 2L, transfer_shard_manifest = 1L,
    legacy_compatibility_manifest = 1L, source_input_manifest = 1L,
    route_environment_map = 1L, acceptance_checks = 1L,
    environment_receipt = 17L
  )
  role_counts <- manifest[, .N, by = role]
  observed_role_counts <- setNames(role_counts$N, role_counts$role)
  if (!file.exists(producer_script) ||
      any(as.character(manifest$producer_package_version) !=
            expected_package_version) ||
      any(as.character(manifest$producer_package_commit) !=
            expected_package_commit) ||
      any(as.character(manifest$producer_script_sha256) != digest::digest(
        producer_script, algo = "sha256", file = TRUE
      )) ||
      any(is.na(observed_role_counts[names(exact_roles)])) ||
      any(observed_role_counts[names(exact_roles)] != exact_roles) ||
      !setequal(names(observed_role_counts), names(exact_roles)) ||
      !setequal(manifest[role == "splitter_manifest", file],
                expected_splitters) ||
      !setequal(manifest[role == "environment_receipt", file],
                .snapshot_candidate_environment_files())) {
    stop("Dependency-overlay candidate producer or role closure is invalid.")
  }
  checks <- .snapshot_read(file.path(root, "acceptance_checks.tsv"))
  if (!all(c("check_id", "pass") %in% names(checks)) || !nrow(checks) ||
      anyDuplicated(checks$check_id) ||
      !setequal(as.character(checks$check_id),
                .snapshot_candidate_check_ids()) ||
      !all(tolower(as.character(checks$pass)) %in% c("true", "1"))) {
    stop("Dependency-overlay candidate has a failed or invalid acceptance check.")
  }
  route_map <- .snapshot_read(file.path(root, "route_environment_map.tsv"))
  expected_routes <- data.table::data.table(
    analysis_surface = c("within", "transfer"),
    environment_route = c("within_s5", "transfer_argos"),
    environment_receipt_path = c(
      "environment/within_s5", "environment/transfer_argos"
    )
  )
  route_columns <- names(expected_routes)
  if (!all(c(names(expected_routes), "environment_receipt_sha256") %in%
           names(route_map)) || nrow(route_map) != 2L ||
      anyDuplicated(route_map$analysis_surface) ||
      !identical(route_map[, ..route_columns], expected_routes) ||
      any(!vapply(route_map$environment_receipt_path, function(prefix) {
        any(startsWith(manifest$file, paste0(prefix, "/")))
      }, logical(1L))) ||
      !identical(
        as.character(route_map$environment_receipt_sha256),
        c(unique(manifest$within_environment_receipt_sha256),
          unique(manifest$transfer_environment_receipt_sha256))
      )) {
    stop("Dependency-overlay candidate route receipts are incomplete or mismatched.")
  }
  if (require_full) {
    if (!"run_mode" %in% names(manifest) ||
        any(as.character(manifest$run_mode) != "full")) {
      stop("Full snapshots require a full dependency-overlay candidate.")
    }
  }
  .snapshot_validate_candidate_source_manifest(root, manifest)
  .snapshot_validate_candidate_science(root)
  list(
    root = root, manifest = manifest, manifest_path = manifest_path,
    within_cells = file.path(root, "within_cells.tsv"),
    within_predictions = file.path(root, "within_predictions.tsv.gz")
  )
}

.snapshot_auc <- function(y, score) {
  y <- as.integer(y)
  if (length(unique(y)) != 2L) return(NA_real_)
  positive <- score[y == 1L]
  negative <- score[y == 0L]
  n_positive <- length(positive)
  n_negative <- length(negative)
  ranks <- rank(c(positive, negative), ties.method = "average")
  mann_whitney <- sum(ranks[seq_len(n_positive)]) -
    n_positive * (n_positive + 1L) / 2
  mann_whitney / (n_positive * n_negative)
}

.snapshot_validate_fixed_auc <- function(cells, predictions) {
  if (!"auc_method" %in% names(cells)) {
    stop("Cell table lacks its fixed-direction method AUC audit field.")
  }
  recalculated <- predictions[, .(
    fold_auc = .snapshot_auc(y, score_m)
  ), by = .(method, cohort, seed, fold)][, .(
    auc_recalculated = mean(fold_auc)
  ), by = .(method, cohort, seed)]
  observed <- merge(
    cells[eligible == TRUE, .(
      method, cohort, seed, auc_reported = as.numeric(auc_method)
    )],
    recalculated,
    by = c("method", "cohort", "seed"), all = TRUE, sort = FALSE
  )
  if (nrow(observed) != sum(cells$eligible) ||
      anyNA(observed[, .(auc_reported, auc_recalculated)]) ||
      any(!is.finite(observed$auc_reported)) ||
      any(!is.finite(observed$auc_recalculated)) ||
      max(abs(observed$auc_reported - observed$auc_recalculated)) > 1e-12) {
    stop("Imported scores do not reproduce every fixed-direction cell AUC.")
  }
  invisible(TRUE)
}

.snapshot_validate_reference <- function(reference) {
  required <- c("fold", "sample_id", "group_id", "y")
  if (!all(required %in% names(reference)) || !nrow(reference)) {
    stop("Canonical reference split is empty or incomplete.")
  }
  if (anyDuplicated(reference[, paste(fold, sample_id, sep = "\r")])) {
    stop("Reference split has duplicate held-out profiles.")
  }
  if (any(reference[, data.table::uniqueN(fold), by = sample_id]$V1 != 1L) ||
      any(reference[, data.table::uniqueN(fold), by = group_id]$V1 != 1L) ||
      any(reference[, data.table::uniqueN(group_id), by = sample_id]$V1 != 1L) ||
      any(reference[, data.table::uniqueN(y), by = group_id]$V1 != 1L) ||
      any(reference[, data.table::uniqueN(y), by = fold]$V1 != 2L)) {
    stop("A held-out profile/group crosses folds or carries invalid labels/classes.")
  }
  invisible(TRUE)
}

.snapshot_rename <- function(x, old, new) {
  if (old %in% names(x) && !new %in% names(x)) data.table::setnames(x, old, new)
  x
}

.snapshot_normalize_cells <- function(x) {
  x <- data.table::copy(x)
  x <- .snapshot_rename(x, "method_id", "method")
  x <- .snapshot_rename(x, "accession", "cohort")
  required <- c("method", "cohort", "seed", "eligible")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Cell table is missing: ", paste(missing, collapse = ", "))
  x[, `:=`(method = as.character(method), cohort = as.character(cohort),
           seed = as.integer(seed), eligible = as.logical(eligible))]
  if (anyNA(x[, ..required])) stop("Cell keys and eligibility cannot be missing.")
  if (!"ineligible_reason" %in% names(x)) x[, ineligible_reason := NA_character_]
  x[, ineligible_reason := as.character(ineligible_reason)]
  if (any(x[eligible == TRUE,
            !is.na(ineligible_reason) & nzchar(trimws(ineligible_reason))]) ||
      any(x[eligible == FALSE,
            is.na(ineligible_reason) | !nzchar(trimws(ineligible_reason))])) {
    stop("Every ineligible cell, and no eligible cell, must carry a reason.")
  }
  key <- x[, paste(method, cohort, seed, sep = "\r")]
  if (anyDuplicated(key)) stop("Cell table has duplicate method/cohort/seed keys.")
  x
}

.snapshot_normalize_predictions <- function(x) {
  x <- data.table::copy(x)
  x <- .snapshot_rename(x, "method_id", "method")
  x <- .snapshot_rename(x, "accession", "cohort")
  x <- .snapshot_rename(x, "score", "score_m")
  required <- c("method", "cohort", "seed", "fold", "sample_id",
                "group_id", "y", "score_m")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Prediction table is missing: ", paste(missing, collapse = ", "))
  }
  x[, `:=`(
    method = as.character(method), cohort = as.character(cohort),
    seed = as.integer(seed), fold = as.integer(fold),
    sample_id = as.character(sample_id), group_id = as.character(group_id),
    y = as.integer(y), score_m = as.numeric(score_m)
  )]
  if (anyNA(x[, ..required]) || any(!is.finite(x$score_m))) {
    stop("Prediction keys, labels, and scores must be complete and finite.")
  }
  if (any(!x$y %in% 0:1)) stop("Predictions contain non-binary outcomes.")
  key <- x[, paste(method, cohort, seed, fold, sample_id, sep = "\r")]
  if (anyDuplicated(key)) stop("Prediction table has duplicate held-out rows.")
  x
}

.snapshot_table_equal <- function(observed, expected) {
  observed <- data.table::as.data.table(observed)
  expected <- data.table::as.data.table(expected)
  if (!identical(names(observed), names(expected)) ||
      nrow(observed) != nrow(expected)) return(FALSE)
  all(vapply(names(observed), function(column) {
    left <- observed[[column]]
    right <- expected[[column]]
    if (is.character(left) && is.character(right)) {
      left[is.na(left)] <- ""
      right[is.na(right)] <- ""
    }
    isTRUE(all.equal(left, right, check.attributes = FALSE, tolerance = 0))
  }, logical(1L)))
}

.snapshot_validate_candidate_science <- function(root) {
  within_cells_raw <- .snapshot_read(file.path(root, "within_cells.tsv"))
  within_predictions_raw <- .snapshot_read(
    file.path(root, "within_predictions.tsv.gz")
  )
  within_cells <- .snapshot_normalize_cells(within_cells_raw)
  within_predictions <- .snapshot_normalize_predictions(within_predictions_raw)
  roster <- data.table::as.data.table(singlesample_method_roster())
  within_methods <- roster[estimand == "within", method_id]
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  cohorts <- sort(unique(within_cells$cohort))
  expected_within <- data.table::CJ(
    method = within_methods, cohort = cohorts, seed = seeds, sorted = FALSE
  )
  if (nrow(within_cells) != 8500L || length(cohorts) != 34L ||
      !setequal(within_cells[, paste(method, cohort, seed, sep = "\r")],
                expected_within[, paste(method, cohort, seed, sep = "\r")])) {
    stop("Dependency-overlay candidate within grid is not the exact 50x34x5 family.")
  }
  eligible_within <- within_cells[eligible == TRUE,
                                  .(method, cohort, seed)]
  predicted_within <- unique(within_predictions[, .(method, cohort, seed)])
  if (!setequal(eligible_within[, paste(method, cohort, seed, sep = "\r")],
                predicted_within[, paste(method, cohort, seed, sep = "\r")])) {
    stop("Dependency-overlay candidate within predictions do not equal eligible cells.")
  }
  within_surface <- within_predictions[, {
    x <- unique(.SD)
    data.table::setorder(x, fold, sample_id, group_id, y)
    .(surface_sha256 = digest::digest(x, algo = "sha256", serialize = TRUE))
  }, by = .(method, cohort, seed),
  .SDcols = c("fold", "sample_id", "group_id", "y")]
  if (any(within_predictions[, data.table::uniqueN(y), by = group_id]$V1 != 1L) ||
      any(within_surface[, data.table::uniqueN(surface_sha256),
                         by = .(cohort, seed)]$V1 != 1L)) {
    stop("Dependency-overlay candidate within held-out surfaces differ.")
  }
  .snapshot_validate_fixed_auc(within_cells, within_predictions)

  transfer_cells_raw <- .snapshot_read(file.path(root, "transfer_cells.tsv"))
  transfer_predictions_raw <- .snapshot_read(
    file.path(root, "transfer_predictions.tsv.gz")
  )
  transfer_cells <- data.table::copy(transfer_cells_raw)
  transfer_predictions <- data.table::copy(transfer_predictions_raw)
  transfer_cell_columns <- c(
    "disease", "held_out", "seed", "method", "eligible",
    "ineligible_reason"
  )
  transfer_prediction_columns <- c(
    "disease", "held_out", "seed", "method", "sample_id", "group_id",
    "y", "score_m", "score_b"
  )
  if (!all(transfer_cell_columns %in% names(transfer_cells)) ||
      !all(transfer_prediction_columns %in% names(transfer_predictions))) {
    stop("Dependency-overlay candidate transfer surfaces lack required columns.")
  }
  transfer_cells[, `:=`(
    disease = as.character(disease), held_out = as.character(held_out),
    seed = as.integer(seed), method = as.character(method),
    eligible = as.logical(eligible)
  )]
  transfer_predictions[, `:=`(
    disease = as.character(disease), held_out = as.character(held_out),
    seed = as.integer(seed), method = as.character(method),
    sample_id = as.character(sample_id), group_id = as.character(group_id),
    y = as.integer(y), score_m = as.numeric(score_m),
    score_b = as.numeric(score_b)
  )]
  transfer_key <- c("disease", "held_out", "seed", "method")
  target <- unique(transfer_cells[, .(disease, held_out)])
  transfer_methods <- roster[estimand == "transfer", method_id]
  expected_transfer <- data.table::rbindlist(lapply(seq_len(nrow(target)),
    function(i) cbind(
      target[i],
      data.table::CJ(seed = seeds, method = transfer_methods, sorted = FALSE)
    )
  ))
  transfer_key_string <- function(x) do.call(
    paste, c(x[, ..transfer_key], sep = "\r")
  )
  if (nrow(transfer_cells) != 2160L || nrow(target) != 18L ||
      anyNA(transfer_cells[, ..transfer_key]) || anyNA(transfer_cells$eligible) ||
      anyDuplicated(transfer_key_string(transfer_cells)) ||
      !setequal(transfer_key_string(transfer_cells),
                transfer_key_string(expected_transfer))) {
    stop("Dependency-overlay candidate transfer grid is not the exact 24x18x5 family.")
  }
  if (any(transfer_cells[eligible == TRUE,
                         !is.na(ineligible_reason) & nzchar(ineligible_reason)]) ||
      any(transfer_cells[eligible == FALSE,
                         is.na(ineligible_reason) | !nzchar(ineligible_reason)])) {
    stop("Dependency-overlay candidate transfer eligibility reasons are invalid.")
  }
  if (nrow(transfer_predictions) &&
      (anyNA(transfer_predictions[, ..transfer_prediction_columns]) ||
       any(!transfer_predictions$y %in% 0:1) ||
       any(!is.finite(transfer_predictions$score_m)) ||
       any(!is.finite(transfer_predictions$score_b)) ||
       anyDuplicated(transfer_predictions[, paste(
         disease, held_out, seed, method, sample_id, sep = "\r"
       )]))) {
    stop("Dependency-overlay candidate transfer predictions are invalid.")
  }
  eligible_transfer <- transfer_cells[eligible == TRUE, ..transfer_key]
  predicted_transfer <- unique(transfer_predictions[, ..transfer_key])
  if (!setequal(transfer_key_string(eligible_transfer),
                transfer_key_string(predicted_transfer))) {
    stop("Dependency-overlay candidate transfer predictions do not equal eligible cells.")
  }
  if (nrow(transfer_predictions)) {
    group_audit <- transfer_predictions[, .(
      n_labels = data.table::uniqueN(y),
      n_held_out = data.table::uniqueN(held_out)
    ), by = group_id]
    surfaces <- transfer_predictions[, {
      x <- unique(.SD)
      data.table::setorder(x, sample_id, group_id, y)
      .(surface_sha256 = digest::digest(x, algo = "sha256", serialize = TRUE))
    }, by = .(disease, held_out, method, seed),
    .SDcols = c("sample_id", "group_id", "y")]
    if (any(group_audit$n_labels != 1L) ||
        any(group_audit$n_held_out != 1L) ||
        any(surfaces[, data.table::uniqueN(surface_sha256),
                     by = .(disease, held_out)]$V1 != 1L)) {
      stop("Dependency-overlay candidate transfer held-out surfaces differ.")
    }
  }
  within_bundle <- readRDS(file.path(root, "within_bundle.rds"))
  transfer_bundle <- readRDS(file.path(root, "transfer_bundle.rds"))
  if (!is.list(within_bundle) || !is.list(transfer_bundle) ||
      !all(c("cells", "predictions") %in% names(within_bundle)) ||
      !all(c("cells", "predictions") %in% names(transfer_bundle))) {
    stop("Dependency-overlay candidate RDS bundles lack required surfaces.")
  }
  equality <- c(
    within_cells = .snapshot_table_equal(within_bundle$cells,
                                          within_cells_raw),
    within_predictions = .snapshot_table_equal(
      within_bundle$predictions, within_predictions_raw
    ),
    transfer_cells = .snapshot_table_equal(transfer_bundle$cells,
                                            transfer_cells_raw),
    transfer_predictions = .snapshot_table_equal(
      transfer_bundle$predictions, transfer_predictions_raw
    )
  )
  if (any(!equality)) {
    stop("Dependency-overlay candidate RDS/TSV surfaces are not identical: ",
         paste(names(equality)[!equality], collapse = ", "))
  }
  invisible(TRUE)
}

.snapshot_overlay <- function(base_cells, base_predictions,
                              overlay_cell_paths, overlay_prediction_paths) {
  if (length(overlay_cell_paths) != length(overlay_prediction_paths)) {
    stop("Overlay cell and prediction path counts must match.")
  }
  cells <- base_cells
  predictions <- base_predictions
  if (!length(overlay_cell_paths)) {
    return(list(cells = cells, predictions = predictions,
                overlay_cell_keys = character()))
  }
  seen_overlay_keys <- character()
  for (i in seq_along(overlay_cell_paths)) {
    add_cells <- .snapshot_normalize_cells(.snapshot_read(overlay_cell_paths[[i]]))
    add_predictions <- .snapshot_normalize_predictions(
      .snapshot_read(overlay_prediction_paths[[i]])
    )
    cell_key <- add_cells[, paste(method, cohort, seed, sep = "\r")]
    duplicate_overlay <- intersect(cell_key, seen_overlay_keys)
    if (length(duplicate_overlay)) {
      stop("Multiple overlays declare the same registered cell key.")
    }
    seen_overlay_keys <- c(seen_overlay_keys, cell_key)
    prediction_key <- unique(add_predictions[, paste(method, cohort, seed, sep = "\r")])
    expected_prediction_key <- add_cells[
      eligible == TRUE, paste(method, cohort, seed, sep = "\r")
    ]
    if (!setequal(expected_prediction_key, prediction_key)) {
      stop("Overlay eligible cells and prediction keys disagree: ",
           overlay_cell_paths[[i]])
    }
    base_key <- cells[, paste(method, cohort, seed, sep = "\r")]
    unknown <- setdiff(cell_key, base_key)
    if (length(unknown)) stop("Overlay contains unregistered cell keys.")
    cells <- cells[!base_key %in% cell_key]
    predictions <- predictions[!predictions[, paste(method, cohort, seed, sep = "\r")] %in% cell_key]
    cells <- data.table::rbindlist(list(cells, add_cells), fill = TRUE, use.names = TRUE)
    predictions <- data.table::rbindlist(
      list(predictions, add_predictions), fill = TRUE, use.names = TRUE
    )
  }
  list(cells = cells, predictions = predictions,
       overlay_cell_keys = seen_overlay_keys)
}

.snapshot_validate_provenance <- function(script, root, units, resolution_path,
                                          output_path) {
  if (!file.exists(script)) stop("Provenance script does not exist: ", script)
  if (!dir.exists(root)) stop("Provenance root does not exist: ", root)
  if (!file.exists(resolution_path)) {
    stop("Provenance overlap resolution does not exist: ", resolution_path)
  }
  sources <- unique(unlist(strsplit(
    paste(units$source_accessions, collapse = ";"), ";", fixed = TRUE
  )))
  sources <- sort(unique(trimws(sources[nzchar(trimws(sources))])))
  if (length(sources) < 2L) stop("Fewer than two source accessions were declared.")
  accession_map <- data.table::rbindlist(lapply(seq_len(nrow(units)), function(i) {
    accession <- trimws(strsplit(
      units$source_accessions[[i]], ";", fixed = TRUE
    )[[1L]])
    data.table::data.table(
      accession = accession[nzchar(accession)], unit_id = units$unit_id[[i]]
    )
  }))
  if (anyDuplicated(accession_map$accession) ||
      !setequal(accession_map$accession, sources)) {
    stop("Each source accession must map to exactly one canonical analysis unit.")
  }
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(root)
  rscript <- file.path(R.home("bin"), "Rscript")
  if (!file.exists(rscript)) stop("Pinned R runtime lacks its canonical Rscript.")
  output <- suppressWarnings(system2(
    rscript, c(normalizePath(script), sources),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  writeLines(output, output_path)
  if (status %in% c(1L, 3L) || !status %in% c(0L, 2L)) {
    stop("Provenance preflight failed with exit status ", status, ".")
  }
  resolution <- data.table::fread(resolution_path)
  required <- c("accession_a", "accession_b", "unit_a", "unit_b")
  if (!all(required %in% names(resolution))) {
    stop("Provenance resolution is missing required columns.")
  }
  resolution[, `:=`(
    accession_a = trimws(as.character(accession_a)),
    accession_b = trimws(as.character(accession_b)),
    unit_a = as.character(unit_a), unit_b = as.character(unit_b)
  )]
  resolved_a <- accession_map$unit_id[match(resolution$accession_a,
                                            accession_map$accession)]
  resolved_b <- accession_map$unit_id[match(resolution$accession_b,
                                            accession_map$accession)]
  if (nrow(resolution) &&
      (anyNA(c(resolved_a, resolved_b)) ||
       any(!resolution$unit_a %in% units$unit_id) ||
       any(!resolution$unit_b %in% units$unit_id) ||
       any(resolution$unit_a != resolved_a) ||
       any(resolution$unit_b != resolved_b))) {
    stop("Provenance resolution is not bound to the cohort accession-to-unit map.")
  }
  overlap_line <- grep("Overlap hits:", output, value = TRUE)
  if (length(overlap_line) != 1L ||
      !grepl("^[[:space:]]*Overlap hits:[[:space:]]*[0-9]+[[:space:]]*$",
             overlap_line)) {
    stop("Provenance preflight must emit exactly one parsable Overlap hits line.")
  }
  overlap_n <- as.integer(sub(
    "^[[:space:]]*Overlap hits:[[:space:]]*([0-9]+)[[:space:]]*$",
    "\\1", overlap_line
  ))
  pair_match <- regexec(
    "^[[:space:]]*([^[:space:]]+)[[:space:]]+<->[[:space:]]+([^[:space:]]+)[[:space:]]*$",
    output, perl = TRUE
  )
  pair_tokens <- regmatches(output, pair_match)
  pair_tokens <- pair_tokens[lengths(pair_tokens) == 3L]
  hit_pairs <- if (length(pair_tokens)) {
    data.table::rbindlist(lapply(pair_tokens, function(x) {
      a <- toupper(x[[2L]])
      b <- toupper(x[[3L]])
      data.table::data.table(accession_a = min(a, b),
                             accession_b = max(a, b))
    }))
  } else data.table::data.table(
    accession_a = character(), accession_b = character()
  )
  resolution_pairs <- data.table::data.table(
    accession_a = pmin(toupper(resolution$accession_a),
                       toupper(resolution$accession_b)),
    accession_b = pmax(toupper(resolution$accession_a),
                       toupper(resolution$accession_b))
  )
  hit_key <- hit_pairs[, paste(accession_a, accession_b, sep = "\r")]
  resolution_key <- resolution_pairs[
    , paste(accession_a, accession_b, sep = "\r")
  ]
  invalid_pairs <- any(hit_pairs$accession_a == hit_pairs$accession_b) ||
    any(resolution_pairs$accession_a == resolution_pairs$accession_b) ||
    anyDuplicated(hit_key) || anyDuplicated(resolution_key)
  if (status == 2L) {
    if (!is.finite(overlap_n) || overlap_n < 1L ||
        nrow(hit_pairs) != overlap_n || nrow(resolution) != overlap_n ||
        invalid_pairs || !setequal(hit_key, resolution_key)) {
      stop("Overlap resolution identities do not match the fresh preflight.")
    }
    if (anyNA(resolution[, ..required]) || any(resolution$unit_a != resolution$unit_b)) {
      stop("A known overlap crosses canonical analysis units.")
    }
  } else if (nrow(resolution) != 0L || nrow(hit_pairs) != 0L || invalid_pairs ||
             (is.finite(overlap_n) && overlap_n != 0L)) {
    stop("A no-overlap preflight requires an empty provenance resolution.")
  }
  list(status = status, overlap_n = ifelse(is.na(overlap_n), 0L, overlap_n),
       resolution = resolution)
}

.snapshot_assert_disjoint_output <- function(output, roots) {
  output <- file.path(normalizePath(dirname(output), mustWork = TRUE),
                      basename(output))
  roots <- unique(vapply(roots, normalizePath, character(1L), mustWork = TRUE))
  nested <- vapply(roots, function(root) {
    identical(output, root) ||
      startsWith(output, paste0(root, .Platform$file.sep)) ||
      startsWith(root, paste0(output, .Platform$file.sep))
  }, logical(1L))
  if (any(nested)) {
    stop("Output must be disjoint from immutable input and package roots.")
  }
  invisible(output)
}

.snapshot_path_entry_exists <- function(path) {
  file.exists(path) || dir.exists(path) || nzchar(Sys.readlink(path))
}

.snapshot_filter_aggregate_predictions <- function(predictions, eligibility) {
  required_predictions <- c("method", "cohort")
  required_eligibility <- c("method_id", "unit_id", "eligible")
  if (!all(required_predictions %in% names(predictions)) ||
      !all(required_eligibility %in% names(eligibility)) ||
      anyNA(eligibility[, ..required_eligibility]) ||
      anyDuplicated(eligibility[, paste(method_id, unit_id, sep = "\r")])) {
    stop("Aggregate eligibility/prediction keys are incomplete or duplicated.")
  }
  eligible_keys <- eligibility[eligible == TRUE,
    paste(method_id, unit_id, sep = "\r")
  ]
  prediction_keys <- predictions[, paste(method, cohort, sep = "\r")]
  out <- predictions[prediction_keys %in% eligible_keys]
  observed_keys <- unique(out[, paste(method, cohort, sep = "\r")])
  if (!setequal(observed_keys, eligible_keys)) {
    stop("Aggregate-eligible method units lack prediction rows.")
  }
  out
}

.snapshot_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .snapshot_args(args)
  output_dir <- normalizePath(dirname(opt[["output-dir"]]), mustWork = TRUE)
  output_dir <- file.path(output_dir, basename(opt[["output-dir"]]))
  if (basename(output_dir) %in% c("", ".", "..") ||
      .snapshot_path_entry_exists(output_dir)) {
    stop("Output snapshot already exists or has an unsafe basename: ",
         output_dir)
  }
  code_provenance_start <- .snapshot_code_provenance()
  validation_path <- normalizePath(opt[["upstream-validation"]], mustWork = TRUE)
  validation_manifest_path <- normalizePath(
    opt[["upstream-validation-manifest"]], mustWork = TRUE
  )
  if (!identical(dirname(validation_path), dirname(validation_manifest_path)) ||
      !identical(basename(validation_path), "validation_checks.tsv") ||
      !identical(basename(validation_manifest_path),
                 "validation_input_manifest.tsv")) {
    stop("Upstream validation receipt and input manifest must be one canonical pair.")
  }
  installed_version <- as.character(utils::packageVersion("OmicSelector"))
  installed_root <- system.file(package = "OmicSelector")
  installed_script <- system.file(
    "scripts", basename(code_provenance_start$script_path),
    package = "OmicSelector"
  )
  if (!identical(code_provenance_start$package_git_commit,
                 opt[["expected-package-commit"]]) ||
      !identical(installed_version, opt[["expected-package-version"]]) ||
      !nzchar(installed_root) ||
      !identical(
        .snapshot_package_tree_sha256(installed_root),
        opt[["expected-installed-package-tree-sha256"]]
      ) ||
      !nzchar(installed_script) || !file.exists(installed_script) ||
      !identical(
        digest::digest(installed_script, algo = "sha256", file = TRUE),
        code_provenance_start$script_sha256
      )) {
    stop("Executing source clone and loaded OmicSelector installation do not ",
         "match the exact package pins.")
  }
  if (opt[["run-mode"]] == "full" &&
      isTRUE(code_provenance_start$package_git_dirty)) {
    stop("Full snapshots require a clean OmicSelector package worktree.")
  }
  runtime_receipt <- normalizePath(opt[["runtime-receipt"]], mustWork = TRUE)
  runtime_guard <- paper3_validate_runtime_receipt(
    runtime_receipt,
    opt[["expected-runtime-receipt-manifest-sha256"]],
    opt[["expected-package-version"]], opt[["expected-package-commit"]],
    opt[["expected-installed-package-tree-sha256"]],
    opt[["runtime-image"]], opt[["expected-runtime-image-sha256"]],
    phase = "start"
  )
  runtime_manifest_path <- file.path(runtime_receipt, "manifest.tsv")
  runtime_manifest <- .snapshot_read(runtime_manifest_path)
  runtime_artifacts <- file.path(runtime_receipt, runtime_manifest$file)

  candidate <- .snapshot_validate_overlay_candidate(
    opt[["dependency-overlay-candidate"]],
    opt[["expected-overlay-manifest-sha256"]],
    opt[["expected-package-version"]], opt[["expected-package-commit"]],
    code_provenance_start$package_root,
    require_full = opt[["run-mode"]] == "full"
  )
  provenance_root <- normalizePath(opt[["provenance-root"]], mustWork = TRUE)
  provenance_script <- normalizePath(opt[["provenance-script"]], mustWork = TRUE)
  canonical_provenance_script <- normalizePath(
    file.path(provenance_root, "code", "check_provenance.R"), mustWork = TRUE
  )
  if (!identical(provenance_script, canonical_provenance_script) ||
      !identical(
        digest::digest(provenance_script, algo = "sha256", file = TRUE),
        opt[["expected-provenance-script-sha256"]]
      )) {
    stop("Provenance preflight must use the exact pinned canonical script.")
  }
  cohort_manifest <- normalizePath(opt[["cohorts"]], mustWork = TRUE)
  if (!identical(
    digest::digest(cohort_manifest, algo = "sha256", file = TRUE),
    opt[["expected-cohort-manifest-sha256"]]
  )) {
    stop("Cohort manifest differs from its exact external SHA-256 pin.")
  }
  provenance_sources <- .snapshot_validate_provenance_sources(
    provenance_root, opt[["expected-provenance-manifest-sha256"]],
    opt[["expected-provenance-union-inventory-sha256"]]
  )
  provenance_manifest <- provenance_sources[["manifest"]]
  provenance_union <- provenance_sources[["union"]]
  provenance_dependencies <- c(provenance_manifest, provenance_union)
  .snapshot_assert_disjoint_output(output_dir, c(
    candidate$root, code_provenance_start$package_root, installed_root,
    runtime_receipt, runtime_guard$runtime_image_path, cohort_manifest,
    provenance_script, provenance_dependencies,
    opt[["provenance-resolution"]], validation_path, validation_manifest_path
  ))
  temp_dir <- tempfile("fair-snapshot-", tmpdir = dirname(output_dir))
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE), add = TRUE)
  candidate_artifacts <- file.path(
    candidate$root, as.character(candidate$manifest$file)
  )
  source_inputs <- data.table::data.table(
    role = c(
      "dependency_overlay_manifest",
      rep("dependency_overlay_artifact", length(candidate_artifacts)),
      "cohort_manifest", "provenance_script", "provenance_resolution",
      "upstream_validation", "upstream_validation_input_manifest",
      "provenance_manifest", "provenance_union_inventory",
      "runtime_receipt_manifest",
      rep("runtime_receipt_artifact", length(runtime_artifacts)),
      "runtime_image"
    ),
    path = c(
      candidate$manifest_path, candidate_artifacts, opt[["cohorts"]],
      provenance_script, opt[["provenance-resolution"]],
      validation_path, validation_manifest_path,
      provenance_dependencies, runtime_manifest_path, runtime_artifacts,
      runtime_guard$runtime_image_path
    )
  )
  if (any(!file.exists(source_inputs$path))) {
    stop("A declared snapshot source input does not exist.")
  }
  source_inputs[, path := normalizePath(path, mustWork = TRUE)]
  source_inputs[, `:=`(
    bytes = file.info(path)$size,
    sha256 = vapply(path, digest::digest, character(1L),
                    algo = "sha256", file = TRUE)
  )]
  .snapshot_validate_upstream_receipt(validation_path)
  .snapshot_validate_input_manifest(validation_manifest_path)

  roster <- data.table::as.data.table(singlesample_method_roster())
  roster[, roster_order := .I]
  roster <- roster[estimand == "within"]
  if (nrow(roster) != 50L || anyDuplicated(roster$method_id)) {
    stop("Installed OmicSelector does not expose the frozen 50-route within roster.")
  }
  units <- .snapshot_read(cohort_manifest)[status == "included"]
  if (!"accession" %in% names(units)) stop("Cohort table lacks accession.")
  data.table::setnames(units, "accession", "unit_id")
  if (nrow(units) != 34L || anyDuplicated(units$unit_id) ||
      !setequal(as.character(units$unit_id), .snapshot_expected_units())) {
    stop("Cohort table must contain the exact frozen 34-unit roster.")
  }
  required_unit <- c("unit_id", "disease", "modality", "biospecimen",
                     "provenance_block", "source_accessions", "outer_k")
  if (!all(required_unit %in% names(units)) ||
      any(vapply(required_unit, function(column) {
        anyNA(units[[column]]) ||
          (is.character(units[[column]]) && any(!nzchar(trimws(units[[column]]))))
      }, logical(1L)))) stop("Cohort metadata is incomplete.")
  expected_k <- ifelse(
    units$unit_id %in% .snapshot_expected_k3(), 3L, 5L
  )
  if (any(as.numeric(units$outer_k) != expected_k) ||
      any(!as.numeric(units$outer_k) %in% c(3L, 5L))) {
    stop("Cohort table violates the frozen outer-fold contract.")
  }
  units[, outer_k := as.integer(outer_k)]
  units[, roster_order__ := match(unit_id, .snapshot_expected_units())]
  data.table::setorder(units, roster_order__)
  units[, roster_order__ := NULL]
  union_inventory <- .snapshot_read(provenance_union)
  union_accessions <- toupper(trimws(as.character(union_inventory$accession)))
  source_accessions <- toupper(trimws(unlist(strsplit(
    paste(units$source_accessions, collapse = ";"), ";", fixed = TRUE
  ))))
  source_accessions <- unique(source_accessions[nzchar(source_accessions)])
  if (!length(source_accessions) ||
      any(!source_accessions %in% union_accessions)) {
    stop("A cohort source accession is absent from the pinned union inventory.")
  }
  units[, primary_unit := !grepl("whole[ -]?blood", biospecimen,
                                 ignore.case = TRUE)]
  if (sum(units$primary_unit) != 33L ||
      !identical(units[primary_unit == FALSE, unit_id], "GSE31568")) {
    stop("Exactly 33 circulating-biofluid primary units are required.")
  }

  cells <- .snapshot_normalize_cells(.snapshot_read(candidate$within_cells))
  predictions <- .snapshot_normalize_predictions(
    .snapshot_read(candidate$within_predictions)
  )
  methods <- roster$method_id
  seeds <- sort(unique(cells$seed))
  if (length(seeds) != 5L) stop("Exactly five outer seeds are required.")
  registered <- data.table::CJ(method = methods, cohort = units$unit_id,
                               seed = seeds, unique = TRUE)
  registered_key <- registered[, paste(method, cohort, seed, sep = "\r")]
  cell_key <- cells[, paste(method, cohort, seed, sep = "\r")]
  if (!setequal(registered_key, cell_key) || length(cell_key) != 8500L) {
    stop("Integrated cells do not equal the frozen 50 x 34 x 5 grid.")
  }
  required_cell_audit <- c(
    "grouping_policy", "n_group_split_violations", "package_version",
    "package_commit", "analysis_code_id", "execution_route",
    "final_package_version", "final_package_commit"
  )
  if (!all(required_cell_audit %in% names(cells)) ||
      any(cells$grouping_policy != "groupKFold") ||
      any(cells$n_group_split_violations != 0L) ||
      any(vapply(required_cell_audit[-c(1L, 2L)], function(column) {
        value <- as.character(cells[[column]])
        anyNA(value) || any(!nzchar(value))
      }, logical(1L))) ||
      any(!grepl("^[0-9a-f]{40}$", cells$package_commit)) ||
      any(!grepl("^[0-9a-f]{40}$", cells$final_package_commit)) ||
      any(!grepl("^[0-9a-f]{16}$", cells$analysis_code_id))) {
    stop("Integrated cells lack complete grouping/code/orientation provenance pins.")
  }
  .snapshot_validate_fixed_auc(cells, predictions)
  cells[, score_orientation_contract := "training_frozen_upstream_validated"]
  data.table::setorder(cells, method, cohort, seed)
  # The following split audit queries every eligible method/unit/seed cell.
  # Key once so those exact lookups do not repeatedly scan the full prediction
  # table (millions of held-out rows in the registered corpus).
  data.table::setkey(predictions, cohort, seed, method, fold, sample_id)
  prediction_cell_key <- unique(predictions[, paste(method, cohort, seed, sep = "\r")])
  eligible_key <- cells[eligible == TRUE, paste(method, cohort, seed, sep = "\r")]
  if (!setequal(prediction_cell_key, eligible_key)) {
    stop("Predictions must exist for every and only eligible cells.")
  }

  split_rows <- list()
  split_map <- list()
  map_index <- 0L
  for (unit in units$unit_id) {
    for (seed_value in seeds) {
      available <- cells[
        cohort == unit & seed == seed_value & eligible == TRUE, method
      ]
      if (!length(available)) {
        stop("No registered method exposes the frozen split for ", unit,
             " seed ", seed_value, ".")
      }
      reference <- predictions[
        cohort == unit & seed == seed_value & method == available[[1L]],
        .(fold, sample_id, group_id, y)
      ]
      data.table::setorder(reference, fold, sample_id, group_id, y)
      .snapshot_validate_reference(reference)
      expected_folds <- seq_len(units[unit_id == unit, outer_k][[1L]])
      if (!identical(sort(unique(reference$fold)), expected_folds)) {
        stop("Canonical split does not expose exactly outer_k folds: ", unit,
             " seed ", seed_value)
      }
      for (method_value in available[-1L]) {
        candidate <- predictions[
          cohort == unit & seed == seed_value & method == method_value,
          .(fold, sample_id, group_id, y)
        ]
        data.table::setorder(candidate, fold, sample_id, group_id, y)
        if (!identical(candidate, reference)) {
          stop("Eligible methods do not share one canonical split: ", unit,
               " seed ", seed_value, " method ", method_value)
        }
      }
      split_id <- digest::digest(reference, algo = "sha256", serialize = TRUE)
      map_index <- map_index + 1L
      split_map[[map_index]] <- data.table::data.table(
        cohort = unit, seed = seed_value, split_id = split_id
      )
      by_fold <- reference[, .(
        n_test_profiles = .N,
        n_test_groups = data.table::uniqueN(group_id),
        n_positive_groups = data.table::uniqueN(group_id[y == 1L]),
        n_negative_groups = data.table::uniqueN(group_id[y == 0L])
      ), by = fold]
      if (any(by_fold$n_positive_groups < 1L) ||
          any(by_fold$n_negative_groups < 1L)) {
        stop("An eligible canonical split contains a one-class test fold.")
      }
      by_fold[, `:=`(unit_id = unit, seed = seed_value, split_id = split_id)]
      split_rows[[length(split_rows) + 1L]] <- by_fold
    }
  }
  split_map <- data.table::rbindlist(split_map)
  expected_split_map <- data.table::CJ(
    cohort = units$unit_id, seed = seeds, sorted = FALSE
  )
  if (nrow(split_map) != nrow(expected_split_map) ||
      anyDuplicated(split_map[, .(cohort, seed)]) ||
      !setequal(split_map[, paste(cohort, seed, sep = "\r")],
                expected_split_map[, paste(cohort, seed, sep = "\r")])) {
    stop("Every frozen unit/seed must expose one exact canonical split.")
  }
  predictions <- merge(predictions, split_map,
                       by = c("cohort", "seed"), all.x = TRUE, sort = FALSE)
  if (anyNA(predictions$split_id)) stop("A prediction lacks a canonical split id.")
  splits <- data.table::rbindlist(split_rows)
  split_count <- splits[, .(
    n_folds = .N, folds = paste(sort(unique(fold)), collapse = ";")
  ), by = .(unit_id, seed)]
  split_count <- merge(
    split_count, units[, .(unit_id, outer_k)], by = "unit_id",
    all = TRUE, sort = FALSE
  )
  expected_m <- units[, .(expected_strata = 5L * outer_k), by = unit_id]
  observed_m <- splits[, .(observed_strata = .N), by = unit_id]
  m_audit <- merge(expected_m, observed_m, by = "unit_id", all = TRUE)
  if (nrow(split_count) != 34L * 5L ||
      anyNA(split_count[, .(n_folds, outer_k)]) ||
      any(split_count$n_folds != split_count$outer_k) ||
      any(split_count$folds != vapply(
        split_count$outer_k, function(k) paste(seq_len(k), collapse = ";"),
        character(1L)
      )) || anyNA(m_audit) ||
      any(m_audit$observed_strata != m_audit$expected_strata)) {
    stop("Canonical splits do not satisfy m = 5 * outer_k for every unit.")
  }
  data.table::setcolorder(splits, c("unit_id", "seed", "fold", "split_id",
                                   "n_test_profiles", "n_test_groups",
                                   "n_positive_groups", "n_negative_groups"))

  eligibility <- cells[, .(
    eligible = .N == length(seeds) && all(eligible),
    eligible_seeds = sum(eligible),
    expected_seeds = length(seeds),
    reason = if (all(eligible)) "eligible" else
      paste(unique(na.omit(ineligible_reason[!eligible])), collapse = " | ")
  ), by = .(method_id = method, unit_id = cohort)]
  if (nrow(eligibility) != 50L * 34L) stop("Eligibility grid is incomplete.")
  pairs <- data.table::as.data.table(singlesample_within_method_pairs(roster))
  primary_units <- units[primary_unit == TRUE, unit_id]
  pair_family <- data.table::rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i]
    support <- sort(intersect(
      eligibility[method_id == pair$method_a & unit_id %in% primary_units &
                    eligible == TRUE, unit_id],
      eligibility[method_id == pair$method_b & unit_id %in% primary_units &
                    eligible == TRUE, unit_id]
    ))
    data.table::data.table(
      pair_id = pair$pair_id, method_a = pair$method_a,
      method_b = pair$method_b, k_support = length(support),
      support_units = paste(support, collapse = ";"),
      pooled_testable = length(support) >= 5L
    )
  }))
  if (nrow(pair_family) != 1225L || anyDuplicated(pair_family$pair_id)) {
    stop("Could not freeze the complete 1,225-pair eligibility family.")
  }

  provenance <- .snapshot_validate_provenance(
    provenance_script, provenance_root, units,
    opt[["provenance-resolution"]], file.path(temp_dir, "provenance_preflight.log")
  )
  data.table::fwrite(provenance$resolution,
                     file.path(temp_dir, "provenance_resolution.tsv"), sep = "\t")

  final_bytes <- file.info(source_inputs$path)$size
  final_hashes <- vapply(source_inputs$path, digest::digest, character(1L),
                         algo = "sha256", file = TRUE)
  if (any(final_bytes != source_inputs$bytes) ||
      !identical(unname(final_hashes), as.character(source_inputs$sha256))) {
    stop("A source input changed while the snapshot was being integrated.")
  }
  data.table::fwrite(source_inputs, file.path(temp_dir, "source_inputs.tsv"),
                     sep = "\t")

  methods_out <- roster[, .(method_id, roster_order, family, role, tier,
                            dep_route, fit_fn, score_fn, pkg_status)]
  aggregate_predictions <- .snapshot_filter_aggregate_predictions(
    predictions, eligibility
  )
  predictions_out <- aggregate_predictions[, .(
    method_id = method, unit_id = cohort, seed, fold, split_id,
    sample_id, group_id, y, score = score_m,
    n_train_profiles = if ("n_train" %in% names(predictions)) n_train else NA_integer_,
    n_test_profiles = if ("n_test" %in% names(predictions)) n_test else NA_integer_,
    package_version = if ("package_version" %in% names(predictions)) package_version else NA_character_,
    package_commit = if ("package_commit" %in% names(predictions)) package_commit else NA_character_,
    analysis_code_id = if ("analysis_code_id" %in% names(predictions)) analysis_code_id else NA_character_
  )]
  data.table::setorder(predictions_out, method_id, unit_id, seed, fold, sample_id)
  data.table::fwrite(methods_out, file.path(temp_dir, "methods.tsv"), sep = "\t")
  data.table::fwrite(units[, ..required_unit][, primary_unit := units$primary_unit],
                     file.path(temp_dir, "units.tsv"), sep = "\t")
  data.table::fwrite(cells, file.path(temp_dir, "cells.tsv"), sep = "\t")
  data.table::fwrite(eligibility, file.path(temp_dir, "eligibility.tsv"), sep = "\t")
  data.table::fwrite(pair_family, file.path(temp_dir, "pair_family.tsv"), sep = "\t")
  data.table::fwrite(splits, file.path(temp_dir, "splits.tsv"), sep = "\t")
  data.table::fwrite(predictions_out, file.path(temp_dir, "predictions.tsv.gz"),
                     sep = "\t", compress = "gzip")

  code_provenance <- .snapshot_code_provenance()
  if (!identical(code_provenance$script_sha256,
                 code_provenance_start$script_sha256) ||
      !identical(code_provenance$package_git_commit,
                 code_provenance_start$package_git_commit) ||
      !identical(code_provenance$package_git_dirty,
                 code_provenance_start$package_git_dirty)) {
    stop("OmicSelector code provenance changed while building the snapshot.")
  }
  if (!identical(
    .snapshot_package_tree_sha256(installed_root),
    opt[["expected-installed-package-tree-sha256"]]
  )) {
    stop("Loaded OmicSelector installation changed while building the snapshot.")
  }
  artifact_paths <- list.files(temp_dir, full.names = TRUE)
  manifest <- data.table::data.table(
    file = basename(artifact_paths),
    bytes = file.info(artifact_paths)$size,
    sha256 = vapply(artifact_paths, digest::digest, character(1L),
                    algo = "sha256", file = TRUE),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    omicselector_version = as.character(utils::packageVersion("OmicSelector")),
    package_git_commit = code_provenance$package_git_commit,
    package_git_dirty = code_provenance$package_git_dirty,
    installed_package_tree_sha256 =
      opt[["expected-installed-package-tree-sha256"]],
    producer_script_sha256 = code_provenance$script_sha256,
    dependency_overlay_manifest_sha256 =
      opt[["expected-overlay-manifest-sha256"]],
    upstream_validation_manifest_sha256 = digest::digest(
      validation_manifest_path, algo = "sha256", file = TRUE
    ),
    snapshot_run_mode = opt[["run-mode"]],
    runtime_receipt_manifest_sha256 =
      opt[["expected-runtime-receipt-manifest-sha256"]],
    runtime_image_sha256 = opt[["expected-runtime-image-sha256"]],
    provenance_exit = provenance$status,
    provenance_overlap_pairs = provenance$overlap_n
  )
  data.table::fwrite(manifest, file.path(temp_dir, "manifest.tsv"), sep = "\t")
  final_source_hashes <- vapply(source_inputs$path, digest::digest,
                                character(1L), algo = "sha256", file = TRUE)
  if (!identical(unname(final_source_hashes),
                 as.character(source_inputs$sha256)) ||
      any(file.info(source_inputs$path)$size != source_inputs$bytes)) {
    stop("An immutable source input changed before snapshot promotion.")
  }
  .snapshot_validate_overlay_candidate(
    candidate$root, opt[["expected-overlay-manifest-sha256"]],
    opt[["expected-package-version"]], opt[["expected-package-commit"]],
    code_provenance_start$package_root,
    require_full = opt[["run-mode"]] == "full"
  )
  .snapshot_validate_upstream_receipt(validation_path)
  .snapshot_validate_input_manifest(validation_manifest_path)
  .snapshot_validate_provenance_sources(
    provenance_root, opt[["expected-provenance-manifest-sha256"]],
    opt[["expected-provenance-union-inventory-sha256"]]
  )
  if (!identical(.snapshot_code_provenance(), code_provenance_start) ||
      !identical(.snapshot_package_tree_sha256(installed_root),
                 opt[["expected-installed-package-tree-sha256"]])) {
    stop("Package source or installation changed before snapshot promotion.")
  }
  paper3_validate_runtime_receipt(
    runtime_receipt,
    opt[["expected-runtime-receipt-manifest-sha256"]],
    opt[["expected-package-version"]], opt[["expected-package-commit"]],
    opt[["expected-installed-package-tree-sha256"]],
    opt[["runtime-image"]], opt[["expected-runtime-image-sha256"]],
    phase = "end", start_guard = runtime_guard
  )
  if (!identical(
    unname(vapply(source_inputs$path, digest::digest, character(1L),
                  algo = "sha256", file = TRUE)),
    as.character(source_inputs$sha256)
  ) || any(file.info(source_inputs$path)$size != source_inputs$bytes)) {
    stop("A source input changed immediately before snapshot promotion.")
  }
  if (!file.rename(temp_dir, output_dir)) {
    stop("Could not atomically promote snapshot to: ", output_dir)
  }
  cat("Created immutable fair-comparison snapshot:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .snapshot_main()
