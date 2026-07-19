#!/usr/bin/env Rscript

# Replay explicitly selected ineligible dependency-benchmark cells without the
# benchmark harness' error-to-NA fallbacks.  The output is an audit witness, not
# a replacement result: every fit/score/runtime error is fatal and every held-
# out method/baseline score is retained.  The dependency bundle splitter may
# use a separately SHA-pinned witness to authorize only finite, fold-complete,
# group-safe structural constant-score ineligibility.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.dsw_value <- function(args, flag, default = NA_character_) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) default else args[[hit + 1L]]
}

.dsw_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  flags <- c(
    bundle = "--bundle", tasks = "--tasks", output = "--output",
    package_root = "--package-root", analysis_root = "--analysis-root",
    cache = "--cache", expected_method = "--expected-method",
    expected_package_version = "--expected-package-version",
    expected_package_commit = "--expected-package-commit",
    expected_cache_commit = "--expected-cache-commit",
    expected_analysis_code_id = "--expected-analysis-code-id",
    expected_task_sha256 = "--expected-task-sha256",
    expected_cache_sha256 = "--expected-cache-sha256",
    expected_fold_engine_sha256 = "--expected-fold-engine-sha256",
    cell_keys = "--cell-keys"
  )
  out <- lapply(flags, function(flag) .dsw_value(args, flag))
  missing <- names(out)[vapply(out, function(x) {
    length(x) != 1L || is.na(x) || !nzchar(x)
  }, logical(1L))]
  if (length(missing)) {
    stop("Missing required argument(s): ",
         paste(unname(flags[missing]), collapse = ", "), call. = FALSE)
  }
  if (!grepl("^[0-9a-f]{40}$", out$expected_package_commit) ||
      !grepl("^[0-9a-f]{40}$", out$expected_cache_commit) ||
      !grepl("^[0-9a-f]{16}$", out$expected_analysis_code_id) ||
      !grepl("^[0-9a-f]{64}$", out$expected_task_sha256) ||
      !grepl("^[0-9a-f]{64}$", out$expected_cache_sha256) ||
      !grepl("^[0-9a-f]{64}$", out$expected_fold_engine_sha256)) {
    stop("Expected commit/code/cache/task pins must be lowercase hexadecimal ",
         "values of lengths 40/16/64.", call. = FALSE)
  }
  raw_keys <- strsplit(out$cell_keys, ",", fixed = TRUE)[[1L]]
  if (!length(raw_keys) ||
      any(!grepl("^[A-Za-z0-9_.-]+:(101|202|303|404|505)$", raw_keys)) ||
      anyDuplicated(raw_keys)) {
    stop("--cell-keys must be a unique comma-separated cohort:seed list ",
         "using frozen seeds.", call. = FALSE)
  }
  pieces <- data.table::tstrsplit(raw_keys, ":", fixed = TRUE)
  out$cell_keys <- data.table(cohort = pieces[[1L]],
                              seed = as.integer(pieces[[2L]]))
  out
}

.dsw_select_cells <- function(cells, method, cell_keys) {
  cells <- data.table::as.data.table(cells)
  keys <- data.table::copy(data.table::as.data.table(cell_keys))
  if (!identical(names(keys), c("cohort", "seed")) || !nrow(keys) ||
      anyNA(keys) || anyDuplicated(keys) ||
      any(!keys$seed %in% c(101L, 202L, 303L, 404L, 505L))) {
    stop("Witness cell-key table is malformed.", call. = FALSE)
  }
  method_value <- method
  keys[, method := method_value]
  data.table::setcolorder(keys, c("method", "cohort", "seed"))
  selected <- merge(keys, cells, by = c("method", "cohort", "seed"),
                    all.x = TRUE, sort = FALSE)
  if (nrow(selected) != nrow(keys) || anyDuplicated(selected[, .(method, cohort,
                                                                  seed)]) ||
      !"eligible" %in% names(selected) || anyNA(selected$eligible) ||
      any(selected$eligible) || any(selected$cohort == "GSE83977") ||
      !"ineligible_reason" %in% names(selected) ||
      anyNA(selected$ineligible_reason) ||
      any(!nzchar(selected$ineligible_reason))) {
    stop("Selected witness cells must be exact non-GSE83977 ineligible keys.",
         call. = FALSE)
  }
  selected
}

.dsw_git <- function(repo, args) {
  out <- system2("git", c("-C", repo, args), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("git inspection failed for ", repo, ": ", paste(out, collapse = " "),
         call. = FALSE)
  }
  out
}

.dsw_analysis_sources <- function(root) {
  roles <- c("benchmark_runner", "benchmark_engine", "paired_auc_engine",
             "design_contract")
  rel <- c(
    "code/analyses/run_singlesample_public_benchmark_20260715.R",
    "code/methods/singlesample_roster_benchmark.R",
    "code/methods/paired_auc_diff_se.R",
    "code/methods/singlesample_public_design.R"
  )
  paths <- file.path(root, rel)
  if (any(!file.exists(paths))) {
    stop("The frozen benchmark analysis source set is incomplete.",
         call. = FALSE)
  }
  data.table(
    role = roles, relative_path = rel,
    sha256 = vapply(paths, digest::digest, character(1L),
                    algo = "sha256", file = TRUE)
  )
}

.dsw_read_tasks <- function(path) {
  tasks <- data.table::fread(
    path, sep = "\t", header = FALSE, quote = "",
    col.names = c("method", "seed", "cohort", "label"),
    colClasses = c("character", "integer", "character", "character")
  )
  methods <- c(
    "ai-scarf", "ai-tabpfn", "coda-codacore", "coda-deepcoda",
    "img-gasfcnn", "inv-scatter", "lrt-deepmaha", "nc-ecod-copod",
    "proto-net", "ssl-vicreg", "tab-tabdpt", "tab-tabicl", "unc-sngp"
  )
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  if (nrow(tasks) != 2210L || anyNA(tasks) ||
      anyDuplicated(tasks[, .(method, seed, cohort)]) ||
      anyDuplicated(tasks$label) || !setequal(tasks$method, methods) ||
      !setequal(tasks$seed, seeds) || uniqueN(tasks$cohort) != 34L ||
      any(tasks[, .N, by = method]$N != 170L) ||
      any(tasks[, .N, by = .(method, cohort)]$N != 5L) ||
      !identical(tasks$label,
                 sprintf("%s__seed_%d__%s", tasks$method, tasks$seed,
                         tasks$cohort))) {
    stop("Frozen task table is not the complete canonical dependency design.",
         call. = FALSE)
  }
  tasks
}

.dsw_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    path <- sub("^--file=", "", file_arg[[1L]])
    if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  }
  ofile <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && file.exists(ofile)) {
    return(normalizePath(ofile, mustWork = TRUE))
  }
  stop("Could not resolve the package-owned witness producer path.",
       call. = FALSE)
}

.dsw_stage <- function(stage, method, cohort, seed, fold, expr) {
  tryCatch(
    expr,
    error = function(e) {
      stop("Structural replay ", stage, " failure [", method, " / ",
           cohort, " / seed ", seed, " / fold ", fold, "]: ",
           conditionMessage(e), call. = FALSE)
    }
  )
}

.dsw_numeric_scores <- function(x, n, stage, method, cohort, seed, fold) {
  if (!is.numeric(x) || length(x) != n || anyNA(x) || any(!is.finite(x))) {
    stop("Structural replay ", stage, " produced missing, non-finite, or ",
         "wrong-length scores [", method, " / ", cohort, " / seed ", seed,
         " / fold ", fold, "].", call. = FALSE)
  }
  as.numeric(x)
}

.dsw_training_direction <- function(score, y) {
  if (length(score) != length(y) || anyNA(score) || any(!is.finite(score)) ||
      anyNA(y) || !all(y %in% c(0L, 1L)) || length(unique(y)) != 2L) {
    stop("Training predictions do not provide a finite two-class direction.",
         call. = FALSE)
  }
  d <- stats::median(score[y == 1L]) - stats::median(score[y == 0L])
  if (!is.finite(d)) stop("Training score direction is non-finite.", call. = FALSE)
  if (d < 0) -1 else 1
}

.dsw_fold_reason <- function(status) {
  ifelse(status == "valid", "ok",
         "method/baseline constant or all-NA on test fold")
}

.dsw_replay_cell <- function(method, cohort_id, seed, cohort, roster,
                             grouped_fold_fn, stratified_fold_fn,
                             fit_call_fn, score_call_fn, baseline_call_fn,
                             outer_k = cohort$outer_k,
                             min_pos_per_fold = cohort$min_per_fold,
                             min_neg_per_fold = cohort$min_per_fold,
                             panel_size = 20L, assay_panel_size = Inf) {
  X <- cohort$expr_per_sample
  y <- as.integer(cohort$y_bin)
  group_id <- cohort$group_id
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  if (n < 4L || length(y) != n || anyNA(y) || !all(y %in% c(0L, 1L)) ||
      is.null(rownames(X)) || anyNA(rownames(X)) || any(!nzchar(rownames(X))) ||
      anyDuplicated(rownames(X))) {
    stop("Structural replay cohort has invalid profiles, labels, or sample IDs.",
         call. = FALSE)
  }
  if (is.null(group_id) || length(group_id) != n || anyNA(group_id) ||
      any(!nzchar(as.character(group_id)))) {
    stop("Structural replay cohort has invalid provenance groups.", call. = FALSE)
  }
  has_groups <- TRUE
  if (length(outer_k) != 1L || is.na(outer_k) || outer_k < 2L ||
      length(min_pos_per_fold) != 1L || is.na(min_pos_per_fold) ||
      length(min_neg_per_fold) != 1L || is.na(min_neg_per_fold)) {
    stop("Structural replay lacks the frozen fold/eligibility contract.",
         call. = FALSE)
  }
  row <- roster[roster$method_id == method, , drop = FALSE]
  if (nrow(row) != 1L || is.na(row$fit_fn[[1L]]) ||
      !nzchar(row$fit_fn[[1L]])) {
    stop("Structural replay requires one dependency method with a real fit_fn.",
         call. = FALSE)
  }
  folds <- .dsw_stage(
    "fold construction", method, cohort_id, seed, 0L,
    if (has_groups) grouped_fold_fn(y, group_id, k = outer_k, seed = seed)
    else stratified_fold_fn(y, k = outer_k, seed = seed)
  )
  if (!is.list(folds) || length(folds) != as.integer(outer_k)) {
    stop("Structural replay did not reproduce the exact outer-fold count.",
         call. = FALSE)
  }

  predictions <- vector("list", outer_k)
  audit <- vector("list", outer_k)
  for (fold in seq_len(outer_k)) {
    test <- as.integer(folds[[fold]])
    train <- setdiff(seq_len(n), test)
    if (!length(test) || !length(train) || anyNA(test) ||
        any(test < 1L | test > n) || anyDuplicated(test)) {
      stop("Structural replay found an invalid outer fold.", call. = FALSE)
    }
    group_violations <- if (has_groups) {
      length(intersect(unique(group_id[test]), unique(group_id[train])))
    } else 0L
    if (group_violations != 0L) {
      stop("Structural replay found a provenance group split.", call. = FALSE)
    }
    y_test <- y[test]
    y_train <- y[train]
    n_pos <- sum(y_test == 1L)
    n_neg <- sum(y_test == 0L)
    if (n_pos < min_pos_per_fold || n_neg < min_neg_per_fold ||
        length(unique(y_train)) != 2L) {
      stop("Structural replay encountered a non-score eligibility failure; ",
           "this witness may authorize only constant-score ineligibility.",
           call. = FALSE)
    }
    X_train <- X[train, , drop = FALSE]
    X_test <- X[test, , drop = FALSE]
    if (is.finite(assay_panel_size) && ncol(X_train) > assay_panel_size) {
      prevalence <- colMeans(is.finite(X_train) & X_train > 0)
      dispersion <- apply(log1p(X_train), 2L, stats::mad, na.rm = TRUE)
      dispersion[!is.finite(dispersion)] <- -Inf
      ranked <- order(prevalence, dispersion, colnames(X_train),
                      decreasing = TRUE)
      features <- colnames(X_train)[head(ranked, as.integer(assay_panel_size))]
      X_train <- X_train[, features, drop = FALSE]
      X_test <- X_test[, features, drop = FALSE]
    }
    baseline <- .dsw_stage(
      "baseline scoring", method, cohort_id, seed, fold,
      baseline_call_fn(X_train, y_train, X_test, panel_size)
    )
    baseline <- .dsw_numeric_scores(
      baseline, length(test), "baseline scoring", method, cohort_id, seed, fold
    )
    model <- .dsw_stage(
      "fit", method, cohort_id, seed, fold,
      fit_call_fn(row$fit_fn[[1L]], X_train, y_train, NULL, list())
    )
    if (is.null(model)) {
      stop("Structural replay fit returned NULL [", method, " / ", cohort_id,
           " / seed ", seed, " / fold ", fold, "].", call. = FALSE)
    }
    train_score <- .dsw_stage(
      "training scoring", method, cohort_id, seed, fold,
      score_call_fn(method, model, X_train, NULL, roster)
    )
    train_score <- .dsw_numeric_scores(
      train_score, length(train), "training scoring", method, cohort_id, seed,
      fold
    )
    direction <- .dsw_stage(
      "training direction", method, cohort_id, seed, fold,
      .dsw_training_direction(train_score, y_train)
    )
    method_score <- .dsw_stage(
      "held-out scoring", method, cohort_id, seed, fold,
      score_call_fn(method, model, X_test, NULL, roster)
    )
    method_score <- direction * .dsw_numeric_scores(
      method_score, length(test), "held-out scoring", method, cohort_id, seed,
      fold
    )
    n_unique_method <- data.table::uniqueN(method_score)
    n_unique_baseline <- data.table::uniqueN(baseline)
    n_unique_outcome <- data.table::uniqueN(y_test)
    if (length(test) < 4L || n_unique_outcome != 2L) {
      stop("Structural replay test fold is not AUC-evaluable.", call. = FALSE)
    }
    structural_basis <- paste(c(
      if (n_unique_method < 2L) "method_constant" else character(),
      if (n_unique_baseline < 2L) "baseline_constant" else character()
    ), collapse = "+")
    status <- if (nzchar(structural_basis)) "structural_constant" else "valid"
    if (!nzchar(structural_basis)) structural_basis <- NA_character_

    predictions[[fold]] <- data.table(
      method = method, cohort = cohort_id, seed = as.integer(seed),
      fold = as.integer(fold), n_train = length(train), n_test = length(test),
      sample_idx = test, sample_id = rownames(X)[test],
      group_id = if (has_groups) as.character(group_id[test]) else rownames(X)[test],
      y = y_test, score_m = method_score, score_b = baseline
    )
    audit[[fold]] <- data.table(
      method = method, cohort = cohort_id, seed = as.integer(seed),
      fold = as.integer(fold), n_train = length(train), n_test = length(test),
      n_pos_test = n_pos, n_neg_test = n_neg,
      n_unique_method = n_unique_method,
      n_unique_baseline = n_unique_baseline,
      n_unique_outcome = n_unique_outcome,
      n_train_finite = sum(is.finite(train_score)),
      n_test_finite = sum(is.finite(method_score) & is.finite(baseline)),
      training_direction = direction,
      n_group_split_violations = group_violations,
      fit_status = "ok", baseline_status = "ok",
      training_score_status = "ok", heldout_score_status = "ok",
      status = status, structural_basis = structural_basis
    )
  }
  predictions <- rbindlist(predictions)
  audit <- rbindlist(audit, fill = TRUE)
  if (nrow(predictions) != n || anyDuplicated(predictions$sample_idx) ||
      !setequal(predictions$sample_idx, seq_len(n))) {
    stop("Structural replay does not cover every profile exactly once.",
         call. = FALSE)
  }
  n_valid <- sum(audit$status == "valid")
  if (n_valid >= outer_k) {
    stop("Structural replay cell is fully eligible; no ineligibility witness ",
         "may be produced.", call. = FALSE)
  }
  reason <- sprintf(
    "n_valid_folds=%d < %d; per-fold: %s", n_valid, outer_k,
    paste(.dsw_fold_reason(audit$status), collapse = "; ")
  )
  list(predictions = predictions, fold_audit = audit,
       n_profiles = n, n_valid_folds = n_valid,
       outer_k = as.integer(outer_k), ineligible_reason = reason,
       n_group_split_violations = 0L)
}

.dsw_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .dsw_args(args)
  bundle_path <- normalizePath(opt$bundle, mustWork = TRUE)
  tasks_path <- normalizePath(opt$tasks, mustWork = TRUE)
  package_root <- normalizePath(opt$package_root, mustWork = TRUE)
  analysis_root <- normalizePath(opt$analysis_root, mustWork = TRUE)
  cache_path <- normalizePath(opt$cache, mustWork = TRUE)
  output_parent <- normalizePath(dirname(opt$output), mustWork = TRUE)
  output_path <- file.path(output_parent, basename(opt$output))
  if (file.exists(output_path)) {
    stop("--output already exists; witnesses are immutable.", call. = FALSE)
  }

  source_bundle_sha <- digest(bundle_path, algo = "sha256", file = TRUE)
  task_sha <- digest(tasks_path, algo = "sha256", file = TRUE)
  cache_sha <- digest(cache_path, algo = "sha256", file = TRUE)
  if (!identical(task_sha, opt$expected_task_sha256) ||
      !identical(cache_sha, opt$expected_cache_sha256)) {
    stop("Task table or cohort cache does not match its expected SHA-256 pin.",
         call. = FALSE)
  }
  tasks <- .dsw_read_tasks(tasks_path)
  commit <- trimws(.dsw_git(package_root, c("rev-parse", "HEAD"))[[1L]])
  dirty <- .dsw_git(package_root, c("status", "--porcelain"))
  version <- unname(read.dcf(file.path(package_root, "DESCRIPTION"))[, "Version"])
  if (!identical(commit, opt$expected_package_commit) || length(dirty) ||
      !identical(version, opt$expected_package_version)) {
    stop("Witness replay requires the exact clean frozen package checkout.",
         call. = FALSE)
  }

  bundle <- readRDS(bundle_path)
  required_bundle <- c("cells", "package_version", "package_commit",
                       "cache_package_commit", "analysis_code_id")
  if (!is.list(bundle) || any(!required_bundle %in% names(bundle)) ||
      !identical(bundle$package_version, opt$expected_package_version) ||
      !identical(bundle$package_commit, opt$expected_package_commit) ||
      !identical(bundle$cache_package_commit, opt$expected_cache_commit) ||
      !identical(bundle$analysis_code_id, opt$expected_analysis_code_id)) {
    stop("Source bundle does not match the declared benchmark pins.",
         call. = FALSE)
  }
  cells <- as.data.table(bundle$cells)
  selected <- .dsw_select_cells(cells, opt$expected_method, opt$cell_keys)
  selected_task_keys <- merge(
    selected[, .(method, cohort, seed)],
    tasks[, .(method, cohort, seed, label)],
    by = c("method", "cohort", "seed"), all.x = TRUE, sort = FALSE
  )
  if (nrow(selected_task_keys) != nrow(selected) ||
      anyNA(selected_task_keys$label)) {
    stop("One or more witness cells are absent from the frozen task table.",
         call. = FALSE)
  }

  analysis_sources <- .dsw_analysis_sources(analysis_root)
  analysis_code_id <- substr(digest(
    paste(analysis_sources$sha256, collapse = "|"), algo = "sha256"
  ), 1L, 16L)
  if (!identical(analysis_code_id, opt$expected_analysis_code_id)) {
    stop("Frozen analysis source bytes do not reproduce analysis_code_id.",
         call. = FALSE)
  }
  fold_source <- file.path(analysis_root,
                           "code/methods/matched_null_benchmark.R")
  if (!file.exists(fold_source)) stop("Frozen fold engine is missing.", call. = FALSE)
  fold_engine_sha <- digest(fold_source, algo = "sha256", file = TRUE)
  if (!identical(fold_engine_sha, opt$expected_fold_engine_sha256)) {
    stop("Frozen fold engine does not match --expected-fold-engine-sha256.",
         call. = FALSE)
  }

  suppressPackageStartupMessages(library(devtools))
  suppressMessages(suppressWarnings(devtools::load_all(package_root, quiet = TRUE)))
  analysis_env <- new.env(parent = globalenv())
  old_root <- Sys.getenv("OMICSELECTOR_PAPER_ROOT", unset = NA_character_)
  on.exit({
    if (is.na(old_root)) Sys.unsetenv("OMICSELECTOR_PAPER_ROOT")
    else Sys.setenv(OMICSELECTOR_PAPER_ROOT = old_root)
  }, add = TRUE)
  Sys.setenv(OMICSELECTOR_PAPER_ROOT = analysis_root)
  sys.source(fold_source, envir = analysis_env)
  sys.source(file.path(analysis_root,
                       "code/methods/paired_auc_diff_se.R"), envir = analysis_env)
  sys.source(file.path(analysis_root,
                       "code/methods/singlesample_roster_benchmark.R"),
             envir = analysis_env)
  sys.source(file.path(analysis_root,
                       "code/methods/singlesample_public_design.R"),
             envir = analysis_env)
  required_fns <- c("srb_verify_cache_bundle", ".grouped_folds",
                    ".stratified_folds", ".srb_call_fit",
                    ".srb_rclr_panel_score")
  if (any(!vapply(required_fns, exists, logical(1L), envir = analysis_env,
                  mode = "function", inherits = TRUE))) {
    stop("Frozen benchmark functions did not load completely.", call. = FALSE)
  }
  get("srb_verify_cache_bundle", envir = analysis_env)(cache_path)
  cohorts <- readRDS(cache_path)
  selected_cohorts <- unique(selected$cohort)
  if (!is.list(cohorts) || any(!selected_cohorts %in% names(cohorts))) {
    stop("Cohort cache lacks one or more selected witness cohorts.",
         call. = FALSE)
  }
  roster <- OmicSelector::singlesample_method_roster()
  if (nrow(roster) != 74L) stop("Frozen package roster is not 74 rows.", call. = FALSE)

  fit_call <- get(".srb_call_fit", envir = analysis_env)
  baseline_call <- function(X_train, y_train, X_test, panel_size) {
    get(".srb_rclr_panel_score", envir = analysis_env)(
      X_train, y_train, X_test, panel_size = panel_size
    )
  }
  score_call <- function(method, model, X, meta, roster) {
    OmicSelector::singlesample_score_call(method, model, X, meta, roster = roster)
  }
  results <- vector("list", nrow(selected))
  cell_rows <- vector("list", nrow(selected))
  setorder(selected, method, seed, cohort)
  for (i in seq_len(nrow(selected))) {
    cell <- selected[i]
    cohort <- cohorts[[cell$cohort]]
    replay <- .dsw_replay_cell(
      method = cell$method, cohort_id = cell$cohort, seed = cell$seed,
      cohort = cohort, roster = roster,
      grouped_fold_fn = get(".grouped_folds", envir = analysis_env),
      stratified_fold_fn = get(".stratified_folds", envir = analysis_env),
      fit_call_fn = fit_call, score_call_fn = score_call,
      baseline_call_fn = baseline_call,
      outer_k = cohort$outer_k, min_pos_per_fold = cohort$min_per_fold,
      min_neg_per_fold = cohort$min_per_fold
    )
    if (!identical(as.character(replay$ineligible_reason),
                   as.character(cell$ineligible_reason)) ||
        !identical(as.integer(replay$n_valid_folds),
                   as.integer(cell$n_valid_folds)) ||
        !identical(as.integer(replay$outer_k), as.integer(cell$outer_k)) ||
        !identical(as.integer(cell$n_group_split_violations), 0L)) {
      stop("Structural replay does not reproduce the source cell reason/fold ",
           "accounting: ", cell$method, " / ", cell$cohort, " / seed ",
           cell$seed, call. = FALSE)
    }
    replay$predictions[, `:=`(
      package_commit = opt$expected_package_commit,
      cache_package_commit = opt$expected_cache_commit,
      analysis_code_id = opt$expected_analysis_code_id
    )]
    replay$fold_audit[, `:=`(
      package_commit = opt$expected_package_commit,
      cache_package_commit = opt$expected_cache_commit,
      analysis_code_id = opt$expected_analysis_code_id
    )]
    results[[i]] <- replay
    cell_rows[[i]] <- data.table(
      method = cell$method, cohort = cell$cohort, seed = cell$seed,
      ineligible_reason = cell$ineligible_reason,
      n_valid_folds = cell$n_valid_folds, outer_k = cell$outer_k,
      n_group_split_violations = cell$n_group_split_violations,
      n_profiles = replay$n_profiles
    )
  }

  producer_path <- .dsw_script_path()
  witness <- list(
    schema_version = "OmicSelector-dependency-structural-witness-v1",
    producer = "paper3_build_dependency_structural_witness.R",
    producer_script_sha256 = digest(producer_path, algo = "sha256", file = TRUE),
    source_bundle_sha256 = source_bundle_sha,
    task_table_sha256 = task_sha, cache_sha256 = cache_sha,
    fold_engine_sha256 = fold_engine_sha,
    package_version = opt$expected_package_version,
    package_commit = opt$expected_package_commit,
    cache_package_commit = opt$expected_cache_commit,
    analysis_code_id = opt$expected_analysis_code_id,
    analysis_sources = analysis_sources,
    replay_sources = rbind(
      analysis_sources,
      data.table(role = "fold_engine", relative_path =
                   "code/methods/matched_null_benchmark.R",
                 sha256 = fold_engine_sha)
    ),
    cells = rbindlist(cell_rows),
    predictions = rbindlist(lapply(results, `[[`, "predictions")),
    fold_audit = rbindlist(lapply(results, `[[`, "fold_audit"), fill = TRUE),
    runtime = list(R.version = R.version.string, platform = R.version$platform,
                   session_info = utils::capture.output(utils::sessionInfo()))
  )
  setorder(witness$cells, method, seed, cohort)
  setorder(witness$predictions, method, seed, cohort, fold, sample_idx)
  setorder(witness$fold_audit, method, seed, cohort, fold)

  if (!identical(source_bundle_sha,
                 digest(bundle_path, algo = "sha256", file = TRUE)) ||
      !identical(task_sha, digest(tasks_path, algo = "sha256", file = TRUE)) ||
      !identical(cache_sha, digest(cache_path, algo = "sha256", file = TRUE)) ||
      !identical(fold_engine_sha,
                 digest(fold_source, algo = "sha256", file = TRUE)) ||
      !identical(commit, trimws(.dsw_git(package_root,
                                        c("rev-parse", "HEAD"))[[1L]])) ||
      length(.dsw_git(package_root, c("status", "--porcelain")))) {
    stop("A pinned replay input changed during witness construction.",
         call. = FALSE)
  }
  tmp <- tempfile(".structural-witness-", tmpdir = output_parent,
                  fileext = ".rds")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(witness, tmp, version = 3L)
  if (!file.rename(tmp, output_path)) {
    stop("Could not atomically promote the immutable witness.", call. = FALSE)
  }
  cat(sprintf("PASS: wrote %d-cell structural witness %s (sha256=%s).\n",
              nrow(witness$cells), output_path,
              digest(output_path, algo = "sha256", file = TRUE)))
  invisible(output_path)
}

if (sys.nframe() == 0L) .dsw_main()
