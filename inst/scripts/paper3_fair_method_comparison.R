#!/usr/bin/env Rscript

# Exact analysis runner for the OmicSelector paper's fair all-method comparison.
# It consumes only an immutable snapshot produced by
# paper3_build_fair_snapshot.R and writes one atomic output bundle.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(OmicSelector)
})

.FAIR_SCRIPT_PATH <- local({
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
  if (!length(candidates)) stop("Could not resolve fair-runner script path.")
  normalizePath(tail(candidates, 1L), mustWork = TRUE)
})
.FAIR_RUNTIME_HELPER <- file.path(
  dirname(.FAIR_SCRIPT_PATH), "paper3_runtime_receipt_common.R"
)
if (!file.exists(.FAIR_RUNTIME_HELPER)) {
  stop("Package-owned runtime-receipt helper is missing.")
}
sys.source(.FAIR_RUNTIME_HELPER, envir = environment())

.fair_args <- function(args) {
  out <- list(bootstrap_reps = 10000L, permutation_reps = 1000L,
              seed = 20260718L, run_mode = "full")
  allowed <- c(
    "input_dir", "output_dir", "expected_snapshot_manifest_sha256",
    "expected_package_version", "expected_package_commit",
    "expected_installed_package_tree_sha256",
    "runtime_receipt", "expected_runtime_receipt_manifest_sha256",
    "runtime_image", "expected_runtime_image_sha256",
    "bootstrap_reps", "permutation_reps", "seed", "run_mode"
  )
  seen <- character()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Arguments require --name=value: ", arg)
    key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg), fixed = TRUE)
    if (!key %in% allowed) stop("Unknown argument: --", gsub("_", "-", key))
    if (key %in% seen) stop("Duplicate argument: --", gsub("_", "-", key))
    seen <- c(seen, key)
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  required <- c("input_dir", "output_dir", "expected_snapshot_manifest_sha256",
                "expected_package_version", "expected_package_commit",
                "expected_installed_package_tree_sha256",
                "runtime_receipt", "expected_runtime_receipt_manifest_sha256",
                "runtime_image", "expected_runtime_image_sha256")
  missing <- required[!vapply(required, function(x) {
    !is.null(out[[x]]) && nzchar(out[[x]])
  }, logical(1L))]
  if (length(missing)) {
    stop("Missing required fair-comparison argument(s): ",
         paste(gsub("_", "-", missing, fixed = TRUE), collapse = ", "))
  }
  if (!grepl("^[0-9a-f]{64}$", out$expected_snapshot_manifest_sha256)) {
    stop("expected-snapshot-manifest-sha256 must be an exact SHA-256 pin.")
  }
  if (!grepl("^[0-9a-f]{40}$", out$expected_package_commit)) {
    stop("expected-package-commit must be an exact Git commit pin.")
  }
  if (!grepl("^[0-9a-f]{64}$", out$expected_installed_package_tree_sha256)) {
    stop("expected-installed-package-tree-sha256 must be an exact SHA-256 pin.")
  }
  if (!grepl("^[0-9a-f]{64}$",
             out$expected_runtime_receipt_manifest_sha256) ||
      !grepl("^[0-9a-f]{64}$", out$expected_runtime_image_sha256)) {
    stop("Runtime receipt/image SHA-256 arguments must be exact pins.")
  }
  out$bootstrap_reps <- as.integer(out$bootstrap_reps)
  out$permutation_reps <- as.integer(out$permutation_reps)
  out$seed <- as.integer(out$seed)
  if (anyNA(c(out$bootstrap_reps, out$permutation_reps, out$seed)) ||
      out$bootstrap_reps < 100L || out$permutation_reps < 100L) {
    stop("Resampling counts must be integers >= 100 and seed must be an integer.")
  }
  if (!out$run_mode %in% c("full", "smoke")) stop("run_mode must be full or smoke.")
  if (out$run_mode == "full" &&
      (out$bootstrap_reps != 10000L || out$permutation_reps != 1000L ||
       out$seed != 20260718L)) {
    stop("Full mode requires the registered 10,000/1,000 resamples and seed 20260718.")
  }
  out
}

.fair_code_provenance <- function() {
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
  status <- system2(
    git_bin, c("-C", package_root, "status", "--porcelain=v1",
             "--untracked-files=all"), stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(commit, "status")) || length(commit) != 1L ||
      !grepl("^[0-9a-f]{40}$", commit) || !is.null(attr(status, "status"))) {
    stop("Could not resolve the OmicSelector code provenance.")
  }
  list(
    script_path = script_path,
    package_root = package_root,
    script_sha256 = digest::digest(script_path, algo = "sha256", file = TRUE),
    package_git_commit = unname(commit[[1L]]),
    package_git_dirty = length(status) > 0L
  )
}

.fair_package_tree_sha256 <- function(root) {
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

.fair_read <- function(path) {
  if (!file.exists(path)) stop("Missing snapshot file: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    data.table::fread(cmd = paste("gzip -cd", shQuote(path)),
                      integer64 = "double")
  } else data.table::fread(path, integer64 = "double")
}

.fair_validate_snapshot <- function(input_dir, expected_manifest_sha256,
                                    require_clean = FALSE) {
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  required <- c("methods.tsv", "units.tsv", "cells.tsv", "eligibility.tsv",
                "pair_family.tsv", "splits.tsv", "predictions.tsv.gz",
                "provenance_preflight.log",
                "provenance_resolution.tsv", "source_inputs.tsv",
                "manifest.tsv")
  missing <- required[!file.exists(file.path(input_dir, required))]
  if (length(missing)) stop("Snapshot is incomplete: ", paste(missing, collapse = ", "))
  manifest_path <- file.path(input_dir, "manifest.tsv")
  if (!identical(
    digest::digest(manifest_path, algo = "sha256", file = TRUE),
    expected_manifest_sha256
  )) {
    stop("Snapshot manifest differs from its exact pre-execution pin.")
  }
  manifest <- data.table::fread(manifest_path)
  required_manifest <- c(
    "file", "bytes", "sha256", "omicselector_version",
    "package_git_commit", "package_git_dirty", "producer_script_sha256",
    "installed_package_tree_sha256", "snapshot_run_mode"
  )
  if (!all(required_manifest %in% names(manifest)) ||
      anyDuplicated(manifest$file) || !setequal(manifest$file, setdiff(required, "manifest.tsv"))) {
    stop("Snapshot manifest does not enumerate the exact input artifacts.")
  }
  present <- list.files(input_dir, recursive = FALSE, all.files = TRUE,
                        no.. = TRUE)
  if (!setequal(present, required) ||
      any(dir.exists(file.path(input_dir, present)))) {
    stop("Snapshot directory contains an unmanifested artifact.")
  }
  actual <- vapply(file.path(input_dir, manifest$file), digest::digest,
                   character(1L), algo = "sha256", file = TRUE)
  sizes <- file.info(file.path(input_dir, manifest$file))$size
  if (!identical(unname(actual), as.character(manifest$sha256)) ||
      any(sizes != manifest$bytes)) {
    stop("Snapshot bytes do not match its producer manifest.")
  }
  if (require_clean && any(manifest$package_git_dirty != FALSE)) {
    stop("Full analysis requires a snapshot produced from a clean package worktree.")
  }
  if (require_clean && any(manifest$snapshot_run_mode != "full")) {
    stop("Full analysis requires a dependency-complete full snapshot.")
  }
  list(input_dir = input_dir, manifest = manifest,
       manifest_sha256 = expected_manifest_sha256)
}

.fair_assert_disjoint_output <- function(output, roots) {
  output <- file.path(normalizePath(dirname(output), mustWork = TRUE),
                      basename(output))
  roots <- unique(vapply(roots, normalizePath, character(1L), mustWork = TRUE))
  nested <- vapply(roots, function(root) {
    identical(output, root) ||
      startsWith(output, paste0(root, .Platform$file.sep)) ||
      startsWith(root, paste0(output, .Platform$file.sep))
  }, logical(1L))
  if (any(nested)) {
    stop("Output must be disjoint from the snapshot and package roots.")
  }
  invisible(output)
}

.fair_path_entry_exists <- function(path) {
  file.exists(path) || dir.exists(path) || nzchar(Sys.readlink(path))
}

.fair_expected_result_files <- function() {
  c(
    "pairwise_cohort.tsv", "pairwise_strata.tsv",
    "pairwise_profile_sensitivity.tsv", "pairwise_meta.tsv",
    "pairwise_loco.tsv", "pairwise_heterogeneity.tsv",
    "complete_support_panels.tsv", "complete_support_performance.tsv",
    "rank_bootstrap.tsv", "coverage_performance.tsv",
    "method_unit_performance.tsv", "method_unit_profile_sensitivity.tsv",
    "hindsight_ceiling.tsv", "eligibility_audit.tsv",
    "prediction_join_audit.tsv", "provenance_preflight.log",
    "provenance_resolution.tsv", "report.md"
  )
}

.fair_validate_staged_result <- function(stage) {
  manifest_path <- file.path(stage, "output_manifest.tsv")
  manifest <- .fair_read(manifest_path)
  if (!all(c("file", "bytes", "sha256") %in% names(manifest)) ||
      !nrow(manifest) || anyDuplicated(manifest$file) ||
      !setequal(as.character(manifest$file), .fair_expected_result_files())) {
    stop("Staged result manifest is malformed.")
  }
  present <- list.files(stage, recursive = FALSE, all.files = TRUE,
                        include.dirs = TRUE, no.. = TRUE)
  if (!setequal(present, c(as.character(manifest$file),
                           "output_manifest.tsv")) ||
      any(dir.exists(file.path(stage, present)))) {
    stop("Staged result inventory differs from its exact manifest.")
  }
  paths <- file.path(stage, manifest$file)
  if (any(!file.exists(paths)) ||
      any(file.info(paths)$size != as.numeric(manifest$bytes)) ||
      !identical(
        unname(vapply(paths, digest::digest, character(1L),
                      algo = "sha256", file = TRUE)),
        as.character(manifest$sha256)
      )) {
    stop("Staged result bytes differ from their exact manifest.")
  }
  invisible(manifest)
}

.fair_auc <- function(y, score) {
  y <- as.integer(y)
  if (length(unique(y)) != 2L) return(NA_real_)
  pos <- score[y == 1L]
  neg <- score[y == 0L]
  n_pos <- length(pos)
  n_neg <- length(neg)
  ranks <- rank(c(pos, neg), ties.method = "average")
  mann_whitney <- sum(ranks[seq_len(n_pos)]) - n_pos * (n_pos + 1L) / 2
  mann_whitney / (n_pos * n_neg)
}

.fair_group_stratum_auc <- function(d) {
  if (any(d[, uniqueN(y), by = group_id]$V1 != 1L)) {
    stop("A biological group carries discordant labels.")
  }
  collapsed <- d[, .(y = unique(y), score = mean(score)), by = group_id]
  .fair_auc(collapsed$y, collapsed$score)
}

.fair_precompute_strata <- function(predictions, splits,
                                    analysis_level = c("group", "profile")) {
  analysis_level <- match.arg(analysis_level)
  label_check <- predictions[, data.table::uniqueN(y),
                             by = .(unit_id, seed, fold, group_id)]
  if (any(label_check$V1 != 1L)) {
    stop("A biological group carries discordant labels within a stratum.")
  }
  if (identical(analysis_level, "group")) {
    analysis_rows <- predictions[, .(y = unique(y), score = mean(score)),
      by = .(method_id, unit_id, seed, fold, group_id)]
    strata <- analysis_rows[, .(auc = .fair_auc(y, score), n_test = .N),
                            by = .(method_id, unit_id, seed, fold)]
    split_count_column <- "n_test_groups"
  } else {
    strata <- predictions[, .(auc = .fair_auc(y, score), n_test = .N),
                          by = .(method_id, unit_id, seed, fold)]
    split_count_column <- "n_test_profiles"
  }
  split_counts <- merge(
    splits[, .(unit_id, seed, fold, split_id,
               registered_test = get(split_count_column))],
    splits[, .(total_observations = sum(get(split_count_column))),
           by = .(unit_id, seed)],
    by = c("unit_id", "seed"), all.x = TRUE, sort = FALSE
  )
  strata <- merge(strata, split_counts,
                  by = c("unit_id", "seed", "fold"),
                  all.x = TRUE, sort = FALSE)
  if (anyNA(strata[, .(split_id, registered_test, total_observations)]) ||
      any(strata$n_test != strata$registered_test) ||
      any(!is.finite(strata$auc))) {
    stop("Precomputed strata disagree with canonical splits.")
  }
  strata[, `:=`(
    n_train = total_observations - n_test,
    analysis_level = analysis_level,
    stratum_id = paste(seed, fold, sep = "::")
  )]
  if (any(strata$n_train <= 0L)) {
    stop("A precomputed stratum has no remaining training observations.")
  }
  data.table::setkey(strata, method_id, unit_id, seed, fold)
  strata
}

.fair_method_unit_performance <- function(stratum_auc, eligibility, splits) {
  rows <- vector("list", nrow(eligibility))
  for (i in seq_len(nrow(eligibility))) {
    e <- eligibility[i]
    if (!isTRUE(e$eligible)) {
      rows[[i]] <- data.table(
        method_id = e$method_id, unit_id = e$unit_id, eligible = FALSE,
        auc = NA_real_, n_strata = 0L, reason = e$reason
      )
      next
    }
    d <- stratum_auc[list(e$method_id, e$unit_id), nomatch = 0L]
    expected <- splits[list(e$unit_id),
                       paste(seed, fold, sep = "::")]
    if (!setequal(unique(d$stratum_id), expected)) {
      stop("Eligible method/unit lacks expected strata: ", e$method_id, "/", e$unit_id)
    }
    if (nrow(d) != length(expected) || any(!is.finite(d$auc))) {
      stop("Eligible method/unit has an undefined AUC.")
    }
    rows[[i]] <- data.table(
      method_id = e$method_id, unit_id = e$unit_id, eligible = TRUE,
      auc = mean(d$auc), n_strata = nrow(d), reason = "eligible"
    )
  }
  data.table::rbindlist(rows, fill = TRUE)
}

.fair_pair_unit <- function(pair, unit, stratum_auc, eligibility, splits,
                            prediction_counts) {
  level <- unique(as.character(stratum_auc$analysis_level))
  if (length(level) != 1L || !level %in% c("profile", "group")) {
    stop("Pairwise strata must carry one explicit analysis level.")
  }
  ea <- eligibility[list(pair$method_a, unit), nomatch = 0L]
  eb <- eligibility[list(pair$method_b, unit), nomatch = 0L]
  base <- data.table(
    pair_id = pair$pair_id, method_a = pair$method_a, method_b = pair$method_b,
    unit_id = unit, analysis_level = level, common_support = FALSE,
    estimate = NA_real_, se = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_zero = NA_real_,
    p_above_margin = NA_real_, p_below_minus_margin = NA_real_,
    expected_m = NA_integer_, observed_m = NA_integer_,
    correction = NA_real_, mean_test_train_ratio = NA_real_, reason = ""
  )
  if (!isTRUE(ea$eligible) || !isTRUE(eb$eligible)) {
    base$reason <- paste0("method_a:", ea$reason, " | method_b:", eb$reason)
    return(list(summary = base, strata = NULL, join = data.table(
      pair_id = pair$pair_id, unit_id = unit, pass = FALSE,
      reason = "not_common_eligible", n_a = 0L, n_b = 0L, n_join = 0L
    )))
  }
  a <- stratum_auc[list(pair$method_a, unit), nomatch = 0L]
  b <- stratum_auc[list(pair$method_b, unit), nomatch = 0L]
  counts_a <- prediction_counts[list(pair$method_a, unit), n_profiles]
  counts_b <- prediction_counts[list(pair$method_b, unit), n_profiles]
  join_cols <- c("unit_id", "seed", "fold", "split_id", "stratum_id",
                 "n_test", "n_train", "analysis_level")
  joined <- merge(
    a[, c(join_cols, "auc"), with = FALSE],
    b[, c(join_cols, "auc"), with = FALSE],
    by = join_cols, all = FALSE, sort = FALSE, suffixes = c("_a", "_b")
  )
  expected <- splits[list(unit), paste(seed, fold, sep = "::")]
  exact <- length(counts_a) == 1L && length(counts_b) == 1L &&
    counts_a == counts_b && nrow(joined) == nrow(a) &&
    nrow(joined) == nrow(b) && nrow(joined) == length(expected) &&
    setequal(joined$stratum_id, expected)
  join_audit <- data.table(
    pair_id = pair$pair_id, unit_id = unit, pass = exact,
    reason = if (exact) "exact" else "heldout_rows_or_labels_differ",
    n_a = if (length(counts_a)) counts_a else 0L,
    n_b = if (length(counts_b)) counts_b else 0L,
    n_join = if (exact) counts_a else 0L
  )
  if (!exact) stop("Pairwise exact join failed: ", pair$pair_id, "/", unit)
  joined[, effect := auc_a - auc_b]
  s <- singlesample_corrected_repeated_cv(
    effect = joined$effect,
    n_test = joined$n_test,
    n_train = joined$n_train,
    stratum_id = joined$stratum_id,
    expected_m = length(expected), expected_strata = expected
  )
  base[, `:=`(
    common_support = TRUE, estimate = s$estimate, se = s$se,
    ci_low = s$ci_low, ci_high = s$ci_high,
    p_zero = s$p_zero_two_sided, p_above_margin = s$p_above_margin,
    p_below_minus_margin = s$p_below_minus_margin,
    expected_m = s$expected_m, observed_m = s$observed_m,
    correction = s$correction,
    mean_test_train_ratio = s$mean_test_train_ratio,
    reason = if (s$inference_available) "evaluable" else "zero_corrected_se"
  )]
  strata <- joined[, .(
    stratum_id, auc_a, auc_b, effect, n_test, n_train, analysis_level,
    pair_id = pair$pair_id, method_a = pair$method_a,
    method_b = pair$method_b, unit_id = unit
  )]
  list(summary = base, strata = strata, join = join_audit)
}

.fair_pair_family <- function(pairs, units, stratum_auc, eligibility, splits,
                              prediction_counts) {
  pair_rows <- stratum_rows <- join_rows <- list()
  index <- 0L
  for (p in seq_len(nrow(pairs))) {
    pair <- pairs[p]
    for (unit in units) {
      index <- index + 1L
      result <- .fair_pair_unit(
        pair, unit, stratum_auc, eligibility, splits, prediction_counts
      )
      pair_rows[[index]] <- result$summary
      join_rows[[index]] <- result$join
      if (!is.null(result$strata)) {
        stratum_rows[[length(stratum_rows) + 1L]] <- result$strata
      }
    }
  }
  list(
    pairwise = data.table::rbindlist(pair_rows, fill = TRUE),
    strata = data.table::rbindlist(stratum_rows, fill = TRUE),
    joins = data.table::rbindlist(join_rows, fill = TRUE)
  )
}

.fair_meta_one <- function(d, pair, support_units) {
  row <- data.table(
    pair_id = pair$pair_id, method_a = pair$method_a, method_b = pair$method_b,
    k_support = length(support_units), k_meta = 0L,
    support_units = paste(support_units, collapse = ";"),
    estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    tau2 = NA_real_, i2 = NA_real_, p_zero = NA_real_,
    p_above_margin = NA_real_, p_below_minus_margin = NA_real_,
    df = NA_real_, pooled_testable = length(support_units) >= 5L,
    status = if (length(support_units) >= 5L) "meta_unavailable" else "descriptive_k_lt_5",
    fit_reason = if (length(support_units) >= 5L) "not_fitted" else "k_support_lt_5"
  )
  use <- d[common_support & is.finite(estimate) & is.finite(se) & se > 0]
  row$k_meta <- nrow(use)
  if (nrow(use) < 5L) {
    row$fit_reason <- "fewer_than_five_finite_unit_effects"
    return(row)
  }
  fit <- tryCatch(
    withCallingHandlers(
      metafor::rma.uni(yi = use$estimate, sei = use$se,
                       method = "REML", test = "knha"),
      warning = function(w) stop(conditionMessage(w), call. = FALSE)
    ), error = identity
  )
  if (inherits(fit, "error")) {
    row$fit_reason <- paste0("rma_failed:", conditionMessage(fit))
    return(row)
  }
  estimate_value <- as.numeric(fit$b[[1L]])
  se_value <- as.numeric(fit$se[[1L]])
  df_value <- nrow(use) - 1L
  row[, `:=`(
    estimate = estimate_value, se = se_value,
    ci_low = as.numeric(fit$ci.lb[[1L]]), ci_high = as.numeric(fit$ci.ub[[1L]]),
    tau2 = as.numeric(fit$tau2), i2 = as.numeric(fit$I2),
    p_zero = as.numeric(fit$pval[[1L]]),
    p_above_margin = stats::pt((estimate_value - 0.05) / se_value,
                               df = df_value,
                               lower.tail = FALSE),
    p_below_minus_margin = stats::pt((estimate_value + 0.05) / se_value,
                                     df = df_value,
                                     lower.tail = TRUE),
    df = df_value, status = "REML_KH", fit_reason = "eligible"
  )]
  row
}

.fair_meta_all <- function(pairwise, pairs, frozen_family) {
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i]
    d <- pairwise[pair_id == pair$pair_id]
    support <- sort(d[common_support == TRUE, unit_id])
    .fair_meta_one(d, pair, support)
  })
  out <- data.table::rbindlist(rows)
  frozen_family <- frozen_family[match(out$pair_id, pair_id)]
  if (anyNA(frozen_family$pair_id) ||
      any(out$k_support != frozen_family$k_support) ||
      any(out$support_units != frozen_family$support_units)) {
    stop("Runtime common support differs from the frozen eligibility family.")
  }
  out[, pooled_testable := as.logical(frozen_family$pooled_testable)]
  frozen <- out[pooled_testable == TRUE, pair_id]
  if (length(frozen)) {
    zero <- singlesample_adjust_zero_by(
      out[pooled_testable == TRUE]$pair_id,
      out[pooled_testable == TRUE]$p_zero, frozen
    )
    relevance <- singlesample_adjust_relevance_by(
      out[pooled_testable == TRUE]$pair_id,
      out[pooled_testable == TRUE]$p_above_margin,
      out[pooled_testable == TRUE]$p_below_minus_margin, frozen
    )
    out[pooled_testable == TRUE, `:=`(
      p_zero_by = zero$p_zero_by,
      p_above_margin_by = relevance$p_above_margin_by,
      p_below_minus_margin_by = relevance$p_below_minus_margin_by,
      zero_family_n = zero$family_n,
      relevance_family_n = relevance$family_n
    )]
  } else {
    out[, `:=`(p_zero_by = NA_real_, p_above_margin_by = NA_real_,
               p_below_minus_margin_by = NA_real_, zero_family_n = 0L,
               relevance_family_n = 0L)]
  }
  out
}

.fair_meta_strata <- function(pairwise, pairs, units, field,
                              inferential = TRUE) {
  annotated <- merge(pairwise, units[, c("unit_id", field), with = FALSE],
                     by = "unit_id", all.x = TRUE, sort = FALSE)
  rows <- list()
  values <- sort(unique(annotated[[field]]))
  for (value in values) {
    for (i in seq_len(nrow(pairs))) {
      pair <- pairs[i]
      d <- annotated[pair_id == pair$pair_id & get(field) == value]
      support <- sort(d[common_support == TRUE, unit_id])
      values_supported <- d[
        common_support == TRUE & is.finite(estimate), estimate
      ]
      if (isTRUE(inferential)) {
        row <- .fair_meta_one(d, pair, support)
        row[, `:=`(
          descriptive_mean = if (length(values_supported))
            mean(values_supported) else NA_real_,
          descriptive_median = if (length(values_supported))
            stats::median(values_supported) else NA_real_,
          descriptive_min = if (length(values_supported))
            min(values_supported) else NA_real_,
          descriptive_max = if (length(values_supported))
            max(values_supported) else NA_real_,
          inference_policy = "REML_KH_when_k_at_least_5"
        )]
      } else {
        row <- .fair_meta_one(d[0L], pair, support)
        row[, `:=`(
          k_meta = length(values_supported), estimate = NA_real_, se = NA_real_,
          ci_low = NA_real_, ci_high = NA_real_, tau2 = NA_real_, i2 = NA_real_,
          p_zero = NA_real_, p_above_margin = NA_real_,
          p_below_minus_margin = NA_real_, df = NA_real_,
          status = "descriptive_only",
          fit_reason = "biospecimen_inference_not_prespecified",
          descriptive_mean = if (length(values_supported))
            mean(values_supported) else NA_real_,
          descriptive_median = if (length(values_supported))
            stats::median(values_supported) else NA_real_,
          descriptive_min = if (length(values_supported))
            min(values_supported) else NA_real_,
          descriptive_max = if (length(values_supported))
            max(values_supported) else NA_real_,
          inference_policy = "descriptive_only"
        )]
      }
      row[, `:=`(stratum_variable = field, stratum_value = value)]
      rows[[length(rows) + 1L]] <- row
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

.fair_loco <- function(pairwise, pairs) {
  rows <- list()
  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i]
    support <- pairwise[pair_id == pair$pair_id & common_support == TRUE]
    finite <- support[is.finite(estimate) & is.finite(se) & se > 0]
    for (unit in support$unit_id) {
      keep <- finite[unit_id != unit]
      row <- data.table(
        pair_id = pair$pair_id, method_a = pair$method_a,
        method_b = pair$method_b, omitted_unit = unit,
        k_support = nrow(support), k_finite_before_omission = nrow(finite),
        k = nrow(keep), estimate = NA_real_, se = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_, tau2 = NA_real_,
        i2 = NA_real_, p_zero = NA_real_, status = "not_evaluable",
        fit_reason = "fewer_than_five_finite_effects_after_omission"
      )
      if (nrow(keep) < 5L) {
        rows[[length(rows) + 1L]] <- row
        next
      }
      fit <- tryCatch(
        withCallingHandlers(
          metafor::rma.uni(yi = keep$estimate, sei = keep$se,
                           method = "REML", test = "knha"),
          warning = function(w) stop(conditionMessage(w), call. = FALSE)
        ), error = identity
      )
      if (inherits(fit, "error")) {
        row$status <- "fit_failed"
        row$fit_reason <- paste0("rma_failed:", conditionMessage(fit))
      } else {
        row[, `:=`(
          estimate = as.numeric(fit$b[[1L]]), se = as.numeric(fit$se[[1L]]),
          ci_low = as.numeric(fit$ci.lb[[1L]]),
          ci_high = as.numeric(fit$ci.ub[[1L]]), tau2 = as.numeric(fit$tau2),
          i2 = as.numeric(fit$I2), p_zero = as.numeric(fit$pval[[1L]]),
          status = "REML_KH", fit_reason = "eligible"
        )]
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (length(rows)) data.table::rbindlist(rows) else data.table()
}

.fair_rank_bootstrap <- function(performance, panels, reps, seed) {
  rows <- list()
  thresholds <- panels$summary$threshold[!panels$summary$degenerate]
  for (threshold_value in thresholds) {
    methods <- panels$method_membership[
      panels$method_membership$threshold == threshold_value &
        panels$method_membership$included_method, "method_id"
    ]
    units <- panels$unit_membership[
      panels$unit_membership$threshold == threshold_value &
        panels$unit_membership$included_unit, "unit_id"
    ]
    rectangle <- performance[method_id %in% methods & unit_id %in% units]
    ranked <- data.table::copy(rectangle)[
      , observed_within_unit_rank := rank(-auc, ties.method = "average"),
      by = unit_id
    ]
    wide <- data.table::dcast(
      ranked, unit_id ~ method_id, value.var = "observed_within_unit_rank"
    )
    if (anyNA(wide)) stop("Complete-support rank rectangle contains NA.")
    observed_rank <- ranked[
      , .(observed_average_rank = mean(observed_within_unit_rank)),
      by = method_id
    ]
    set.seed(seed + as.integer(round(threshold_value * 1000)))
    rank_draws <- matrix(NA_real_, nrow = reps, ncol = length(methods),
                         dimnames = list(NULL, methods))
    top1 <- top3 <- top5 <- matrix(
      FALSE, nrow = reps, ncol = length(methods),
      dimnames = list(NULL, methods)
    )
    for (b in seq_len(reps)) {
      sampled <- sample(seq_len(nrow(wide)), nrow(wide), replace = TRUE)
      average_ranks <- colMeans(as.matrix(wide[sampled, ..methods]))
      rank_draws[b, ] <- average_ranks
      inclusive_rank <- rank(average_ranks, ties.method = "min")
      top1[b, ] <- inclusive_rank <= 1L
      top3[b, ] <- inclusive_rank <= 3L
      top5[b, ] <- inclusive_rank <= 5L
    }
    for (method in methods) {
      values <- rank_draws[, method]
      rows[[length(rows) + 1L]] <- data.table(
        threshold = threshold_value, method_id = method,
        rank_median = stats::median(values),
        rank_ci_low = unname(stats::quantile(values, 0.025, type = 6L)),
        rank_ci_high = unname(stats::quantile(values, 0.975, type = 6L)),
        observed_average_rank = observed_rank[
          method_id == method, observed_average_rank
        ],
        p_top1 = mean(top1[, method]), p_top3 = mean(top3[, method]),
        p_top5 = mean(top5[, method]), bootstrap_reps = reps, seed = seed,
        tie_policy = "inclusive_min_rank_for_topk;average_rank_for_intervals"
      )
    }
  }
  if (length(rows)) data.table::rbindlist(rows) else data.table()
}

.fair_hindsight <- function(predictions, performance, panels, units, reps, seed) {
  panel_row <- panels$summary[panels$summary$threshold == 0.95, , drop = FALSE]
  if (!nrow(panel_row) || isTRUE(panel_row$degenerate[[1L]])) return(data.table())
  methods <- as.character(panels$method_membership$method_id[
    panels$method_membership$threshold == 0.95 &
      panels$method_membership$included_method
  ])
  panel_units <- as.character(panels$unit_membership$unit_id[
    panels$unit_membership$threshold == 0.95 &
      panels$unit_membership$included_unit
  ])
  rows <- list()
  for (u in seq_along(panel_units)) {
    unit <- panel_units[[u]]
    d <- predictions[unit_id == unit & method_id %in% methods]
    d[, stratum_id := paste(seed, fold, sep = "::")]
    labels <- unique(d[, .(stratum_id, group_id, y)])
    if (any(labels[, uniqueN(y), by = .(stratum_id, group_id)]$V1 != 1L) ||
        any(labels[, uniqueN(y), by = stratum_id]$V1 != 2L)) {
      stop("Discordant labels in hindsight unit: ", unit)
    }
    observed <- performance[unit_id == unit & method_id %in% methods]
    observed_max <- max(observed$auc)
    set.seed(seed + u * 1009L)
    null <- numeric()
    attempts <- 0L
    max_attempts <- reps * 100L
    while (length(null) < reps && attempts < max_attempts) {
      attempts <- attempts + 1L
      # Permute biological-group labels within each held-out stratum. This
      # preserves the fixed stratified fold's group-level class counts and the
      # joint score correlation across methods.
      perm <- labels[, .(group_id, y_perm = sample(y)), by = stratum_id]
      x <- merge(d[, .(method_id, stratum_id, group_id, sample_id, score)],
                 perm, by = c("stratum_id", "group_id"), sort = FALSE)
      collapsed <- x[, .(y_perm = unique(y_perm), score = mean(score)),
                     by = .(method_id, stratum_id, group_id)]
      aucs <- collapsed[, .(auc = .fair_auc(y_perm, score)),
                by = .(method_id, stratum_id)]
      if (any(!is.finite(aucs$auc))) next
      method_mean <- aucs[, .(auc = mean(auc)), by = method_id]
      if (nrow(method_mean) != length(methods)) next
      null <- c(null, max(method_mean$auc))
    }
    valid_fraction <- length(null) / attempts
    if (length(null) < reps || valid_fraction < 0.80) {
      rows[[u]] <- data.table(
        unit_id = unit, observed_max = observed_max, null_median = NA_real_,
        null_ci_low = NA_real_, null_ci_high = NA_real_, headroom = NA_real_,
        valid_permutations = length(null), attempts = attempts,
        valid_fraction = valid_fraction, status = "permutation_guard_failed",
        permutation_reps = reps, seed = seed, panel_threshold = 0.95,
        n_methods = length(methods), analysis_level = "group"
      )
    } else {
      rows[[u]] <- data.table(
        unit_id = unit, observed_max = observed_max,
        null_median = stats::median(null),
        null_ci_low = unname(stats::quantile(null, 0.025, type = 6L)),
        null_ci_high = unname(stats::quantile(null, 0.975, type = 6L)),
        headroom = observed_max - stats::median(null),
        valid_permutations = length(null), attempts = attempts,
        valid_fraction = valid_fraction, status = "evaluable",
        permutation_reps = reps, seed = seed, panel_threshold = 0.95,
        n_methods = length(methods), analysis_level = "group"
      )
    }
  }
  data.table::rbindlist(rows)
}

.fair_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("The exact fair-comparison runner requires the suggested package 'metafor'.")
  }
  opt <- .fair_args(args)
  snapshot <- .fair_validate_snapshot(
    opt$input_dir, opt$expected_snapshot_manifest_sha256,
    require_clean = opt$run_mode == "full"
  )
  input_dir <- snapshot$input_dir
  code_provenance <- .fair_code_provenance()
  installed_version <- as.character(utils::packageVersion("OmicSelector"))
  installed_root <- system.file(package = "OmicSelector")
  installed_script <- system.file(
    "scripts", basename(code_provenance$script_path), package = "OmicSelector"
  )
  if (!identical(code_provenance$package_git_commit,
                 opt$expected_package_commit) ||
      !identical(installed_version, opt$expected_package_version) ||
      !nzchar(installed_root) ||
      !identical(
        .fair_package_tree_sha256(installed_root),
        opt$expected_installed_package_tree_sha256
      ) ||
      !nzchar(installed_script) || !file.exists(installed_script) ||
      !identical(
        digest::digest(installed_script, algo = "sha256", file = TRUE),
        code_provenance$script_sha256
      )) {
    stop("Executing source clone and loaded OmicSelector installation do not ",
         "match the exact package pins.")
  }
  snapshot_commit <- unique(snapshot$manifest$package_git_commit)
  snapshot_version <- unique(snapshot$manifest$omicselector_version)
  snapshot_installed_tree <- unique(
    snapshot$manifest$installed_package_tree_sha256
  )
  if (length(snapshot_commit) != 1L ||
      !identical(snapshot_commit, opt$expected_package_commit) ||
      length(snapshot_version) != 1L ||
      !identical(snapshot_version, opt$expected_package_version) ||
      length(snapshot_installed_tree) != 1L ||
      !identical(snapshot_installed_tree,
                 opt$expected_installed_package_tree_sha256)) {
    stop("Snapshot and runner do not share one OmicSelector version/commit pin.")
  }
  if (opt$run_mode == "full" && code_provenance$package_git_dirty) {
    stop("Full analysis requires a clean OmicSelector package worktree.")
  }
  snapshot_runtime_receipt <- unique(
    snapshot$manifest$runtime_receipt_manifest_sha256
  )
  snapshot_runtime_image <- unique(snapshot$manifest$runtime_image_sha256)
  if (length(snapshot_runtime_receipt) != 1L ||
      !identical(snapshot_runtime_receipt,
                 opt$expected_runtime_receipt_manifest_sha256) ||
      length(snapshot_runtime_image) != 1L ||
      !identical(snapshot_runtime_image, opt$expected_runtime_image_sha256)) {
    stop("Snapshot and runner do not share one exact runtime receipt/image pin.")
  }
  runtime_receipt <- normalizePath(opt$runtime_receipt, mustWork = TRUE)
  runtime_guard <- paper3_validate_runtime_receipt(
    runtime_receipt, opt$expected_runtime_receipt_manifest_sha256,
    opt$expected_package_version, opt$expected_package_commit,
    opt$expected_installed_package_tree_sha256, opt$runtime_image,
    opt$expected_runtime_image_sha256, phase = "start"
  )
  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (basename(output_dir) %in% c("", ".", "..") ||
      .fair_path_entry_exists(output_dir)) {
    stop("Output bundle already exists or has an unsafe basename: ", output_dir)
  }
  .fair_assert_disjoint_output(
    output_dir, c(input_dir, code_provenance$package_root, installed_root,
                  runtime_receipt, opt$runtime_image)
  )
  temp_dir <- tempfile("fair-analysis-", tmpdir = output_parent)
  dir.create(temp_dir, recursive = TRUE)
  on.exit(if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE), add = TRUE)

  methods <- .fair_read(file.path(input_dir, "methods.tsv"))
  units <- .fair_read(file.path(input_dir, "units.tsv"))
  cells <- .fair_read(file.path(input_dir, "cells.tsv"))
  eligibility <- .fair_read(file.path(input_dir, "eligibility.tsv"))
  frozen_family <- .fair_read(file.path(input_dir, "pair_family.tsv"))
  splits <- .fair_read(file.path(input_dir, "splits.tsv"))
  predictions <- .fair_read(file.path(input_dir, "predictions.tsv.gz"))
  required_prediction <- c(
    "method_id", "unit_id", "seed", "fold", "split_id", "sample_id",
    "group_id", "y", "score"
  )
  if (!all(required_prediction %in% names(predictions))) {
    stop("Snapshot predictions are missing required columns.")
  }
  # fread can conservatively retain a numeric score column as character when
  # a very large heterogeneous file is streamed through gzip. Parse the
  # producer-validated field explicitly and fail if any token is not finite.
  predictions[, score := suppressWarnings(as.numeric(score))]
  roster <- data.table::as.data.table(singlesample_method_roster())[estimand == "within"]
  if (nrow(methods) != 50L || !identical(methods$method_id, roster$method_id) ||
      nrow(units) != 34L || sum(units$primary_unit) != 33L ||
      nrow(cells) != 8500L || nrow(eligibility) != 1700L ||
      any(cells$score_orientation_contract !=
            "training_frozen_upstream_validated")) {
    stop("Snapshot roster/unit/eligibility dimensions are not frozen.")
  }
  split_contract <- merge(
    splits[, .(
      n_folds = .N, fold_ids = paste(sort(unique(fold)), collapse = ";")
    ), by = .(unit_id, seed)],
    units[, .(unit_id, outer_k)], by = "unit_id", all = TRUE,
    sort = FALSE
  )
  if (nrow(split_contract) != 34L * 5L || anyNA(split_contract) ||
      any(split_contract$n_folds != split_contract$outer_k) ||
      any(split_contract$fold_ids != vapply(
        split_contract$outer_k,
        function(k) paste(seq_len(k), collapse = ";"), character(1L)
      ))) {
    stop("Snapshot splits violate the exact m = 5 * outer_k contract.")
  }
  pred_key <- predictions[, paste(method_id, unit_id, seed, fold, sample_id,
                                  sep = "\r")]
  if (anyDuplicated(pred_key) || anyNA(predictions$score) ||
      any(!is.finite(predictions$score))) {
    stop("Snapshot predictions are duplicated or non-finite.")
  }
  split_check <- predictions[, uniqueN(split_id), by = .(unit_id, seed, method_id)]
  if (any(split_check$V1 != 1L)) stop("A method has multiple split ids per unit/seed.")
  canonical_check <- unique(predictions[, .(unit_id, seed, split_id)])[
    , uniqueN(split_id), by = .(unit_id, seed)]
  if (any(canonical_check$V1 != 1L)) stop("Methods do not share one canonical split.")
  data.table::setkey(predictions, method_id, unit_id)
  data.table::setkey(eligibility, method_id, unit_id)
  data.table::setkey(splits, unit_id, seed, fold)
  prediction_counts <- predictions[, .(n_profiles = .N),
                                   by = .(method_id, unit_id)]
  data.table::setkey(prediction_counts, method_id, unit_id)
  stratum_auc <- .fair_precompute_strata(
    predictions, splits, analysis_level = "group"
  )
  stratum_auc_profile <- .fair_precompute_strata(
    predictions, splits, analysis_level = "profile"
  )
  eligible_count_check <- merge(
    eligibility[eligible == TRUE, .(method_id, unit_id)],
    prediction_counts, by = c("method_id", "unit_id"), all.x = TRUE
  )
  if (anyNA(eligible_count_check$n_profiles) ||
      any(eligible_count_check[, data.table::uniqueN(n_profiles),
                               by = unit_id]$V1 != 1L)) {
    stop("Eligible methods do not share identical held-out profile counts.")
  }

  pairs <- data.table::as.data.table(singlesample_within_method_pairs(roster))
  if (nrow(pairs) != 1225L || nrow(frozen_family) != 1225L ||
      !identical(pairs$pair_id, frozen_family$pair_id) ||
      anyDuplicated(frozen_family$pair_id)) {
    stop("The exact frozen 1,225-pair family is absent.")
  }
  primary_units <- units[primary_unit == TRUE, unit_id]
  all_units <- units$unit_id
  primary_pair_family <- .fair_pair_family(
    pairs, all_units, stratum_auc, eligibility, splits, prediction_counts
  )
  profile_pair_family <- .fair_pair_family(
    pairs, all_units, stratum_auc_profile, eligibility, splits, prediction_counts
  )
  pairwise <- primary_pair_family$pairwise
  pair_strata <- primary_pair_family$strata
  join_audit <- primary_pair_family$joins
  pairwise_profile <- profile_pair_family$pairwise
  profile_sensitivity <- merge(
    pairwise[, .(
      pair_id, method_a, method_b, unit_id, common_support,
      group_collapsed_estimate = estimate, group_collapsed_se = se,
      group_collapsed_ci_low = ci_low, group_collapsed_ci_high = ci_high
    )],
    pairwise_profile[, .(
      pair_id, unit_id, profile_common_support = common_support,
      profile_estimate = estimate, profile_se = se,
      profile_ci_low = ci_low, profile_ci_high = ci_high
    )],
    by = c("pair_id", "unit_id"), all = TRUE, sort = FALSE
  )
  profile_sensitivity[, difference_profile_minus_group :=
                        profile_estimate - group_collapsed_estimate]
  data.table::setcolorder(profile_sensitivity, c(
    "pair_id", "method_a", "method_b", "unit_id", "common_support",
    "group_collapsed_estimate", "group_collapsed_se",
    "group_collapsed_ci_low", "group_collapsed_ci_high",
    "profile_common_support", "profile_estimate", "profile_se",
    "profile_ci_low", "profile_ci_high", "difference_profile_minus_group"
  ))
  pairwise <- merge(pairwise, units[, .(unit_id, primary_unit)],
                    by = "unit_id", all.x = TRUE, sort = FALSE)
  if (nrow(pairwise) != 1225L * 34L || any(!join_audit$pass &
                                           join_audit$reason != "not_common_eligible")) {
    stop("Pairwise accounting or exact joins failed.")
  }
  if (nrow(pairwise_profile) != nrow(pairwise) ||
      any(profile_pair_family$joins$pass != join_audit$pass) ||
      nrow(profile_sensitivity) != nrow(pairwise) ||
      any(profile_sensitivity$common_support !=
            profile_sensitivity$profile_common_support)) {
    stop("Group-primary and profile-weighted sensitivity families diverge.")
  }
  pairwise_primary <- pairwise[primary_unit == TRUE]
  meta <- .fair_meta_all(pairwise_primary, pairs, frozen_family)
  loco <- .fair_loco(pairwise_primary, pairs)
  heterogeneity <- data.table::rbindlist(list(
    .fair_meta_strata(pairwise_primary, pairs, units, "modality", TRUE),
    .fair_meta_strata(pairwise_primary, pairs, units, "biospecimen", FALSE)
  ), fill = TRUE)

  performance_all <- .fair_method_unit_performance(
    stratum_auc, eligibility, splits
  )
  performance_profile <- .fair_method_unit_performance(
    stratum_auc_profile, eligibility, splits
  )
  performance_sensitivity <- merge(
    performance_all[, .(
      method_id, unit_id, eligible, group_collapsed_auc = auc, n_strata
    )],
    performance_profile[, .(
      method_id, unit_id, profile_auc = auc
    )],
    by = c("method_id", "unit_id"), all = TRUE, sort = FALSE
  )
  performance_sensitivity[, difference_profile_minus_group :=
                            profile_auc - group_collapsed_auc]
  performance_sensitivity <- merge(
    performance_sensitivity,
    splits[, .(repeated_profiles = any(n_test_profiles > n_test_groups)),
           by = unit_id],
    by = "unit_id", all.x = TRUE, sort = FALSE
  )
  data.table::setcolorder(performance_sensitivity, c(
    "method_id", "unit_id", "eligible", "group_collapsed_auc", "n_strata",
    "profile_auc", "difference_profile_minus_group", "repeated_profiles"
  ))
  performance <- performance_all[unit_id %in% primary_units]
  panels <- singlesample_complete_support_panels(
    eligibility[, .(method_id, unit_id, eligible)], roster = roster,
    thresholds = c(0.50, 0.80, 0.95), primary_units = primary_units
  )
  panel_table <- data.table::rbindlist(list(
    data.table::as.data.table(panels$summary)[, item_type := "summary"],
    data.table::as.data.table(panels$method_membership)[, item_type := "method"],
    data.table::as.data.table(panels$unit_membership)[, item_type := "unit"]
  ), fill = TRUE)
  rectangle_rows <- list()
  for (threshold in panels$summary$threshold) {
    selected_methods <- panels$method_membership[
      panels$method_membership$threshold == threshold &
        panels$method_membership$included_method, "method_id"
    ]
    selected_units <- panels$unit_membership[
      panels$unit_membership$threshold == threshold &
        panels$unit_membership$included_unit, "unit_id"
    ]
    rectangle_rows[[length(rectangle_rows) + 1L]] <- performance[
      method_id %in% selected_methods & unit_id %in% selected_units
    ][, threshold := threshold]
  }
  rectangles <- data.table::rbindlist(rectangle_rows, fill = TRUE)
  if (nrow(rectangles)) {
    rectangles[, within_unit_rank := rank(-auc, ties.method = "average"),
               by = .(threshold, unit_id)]
  }
  rank_bootstrap <- .fair_rank_bootstrap(
    performance, panels, opt$bootstrap_reps, opt$seed
  )
  coverage <- merge(
    eligibility[unit_id %in% primary_units,
                .(eligible_units = sum(eligible),
                  eligibility_coverage = mean(eligible),
                  eligibility_units = paste(unit_id[eligible], collapse = ";")),
                by = method_id],
    performance[eligible == TRUE, .(median_auc = stats::median(auc),
                            performance_units = paste(unit_id, collapse = ";")),
                by = method_id],
    by = "method_id", all.x = TRUE
  )
  expected_prediction_rows <- merge(
    eligibility[eligible == TRUE & unit_id %in% primary_units,
                .(method_id, unit_id)],
    splits[, .(expected_rows = sum(n_test_profiles)), by = unit_id],
    by = "unit_id", all.x = TRUE, sort = FALSE
  )[, .(expected_prediction_rows = sum(expected_rows)), by = method_id]
  actual_prediction_rows <- predictions[unit_id %in% primary_units,
    .(prediction_rows = .N), by = method_id]
  coverage <- merge(coverage, expected_prediction_rows, by = "method_id",
                    all.x = TRUE, sort = FALSE)
  coverage <- merge(coverage, actual_prediction_rows, by = "method_id",
                    all.x = TRUE, sort = FALSE)
  coverage <- merge(methods[, .(method_id, roster_order, family, role, tier,
                                dep_route, pkg_status)], coverage,
                    by = "method_id", sort = FALSE)
  coverage[is.na(expected_prediction_rows), expected_prediction_rows := 0L]
  coverage[is.na(prediction_rows), prediction_rows := 0L]
  coverage[, `:=`(
    prediction_completeness = fifelse(
      expected_prediction_rows > 0,
      prediction_rows / expected_prediction_rows, 1
    ),
    prediction_complete = prediction_rows == expected_prediction_rows,
    deployment_status = fcase(
      pkg_status %in% c("present", "primitive"), "package_native",
      pkg_status == "new", "package_added_primary",
      pkg_status == "new-secondary", "package_added_secondary",
      default = "unknown"
    )
  )]
  if (any(!coverage$prediction_complete) ||
      any(abs(coverage$prediction_completeness - 1) > 1e-12)) {
    stop("Eligible-method prediction accounting is incomplete.")
  }
  coverage[, pareto_optimal := vapply(seq_len(.N), function(i) {
    if (!is.finite(median_auc[[i]]) ||
        !is.finite(eligibility_coverage[[i]])) return(FALSE)
    !any(
      seq_len(.N) != i & is.finite(median_auc) &
        eligibility_coverage >= eligibility_coverage[[i]] &
        median_auc >= median_auc[[i]] &
        (eligibility_coverage > eligibility_coverage[[i]] |
           median_auc > median_auc[[i]])
    )
  }, logical(1L))]
  hindsight <- .fair_hindsight(
    predictions, performance, panels, units, opt$permutation_reps, opt$seed
  )
  code_provenance_final <- .fair_code_provenance()
  if (!identical(code_provenance_final, code_provenance)) {
    stop("OmicSelector code provenance changed while running the comparison.")
  }
  if (!identical(
    .fair_package_tree_sha256(installed_root),
    opt$expected_installed_package_tree_sha256
  )) {
    stop("Loaded OmicSelector installation changed during the comparison.")
  }

  data.table::fwrite(pairwise, file.path(temp_dir, "pairwise_cohort.tsv"), sep = "\t")
  data.table::fwrite(pair_strata, file.path(temp_dir, "pairwise_strata.tsv"), sep = "\t")
  data.table::fwrite(profile_sensitivity,
                     file.path(temp_dir, "pairwise_profile_sensitivity.tsv"),
                     sep = "\t")
  data.table::fwrite(meta, file.path(temp_dir, "pairwise_meta.tsv"), sep = "\t")
  data.table::fwrite(loco, file.path(temp_dir, "pairwise_loco.tsv"), sep = "\t")
  data.table::fwrite(heterogeneity, file.path(temp_dir, "pairwise_heterogeneity.tsv"), sep = "\t")
  data.table::fwrite(panel_table, file.path(temp_dir, "complete_support_panels.tsv"), sep = "\t")
  data.table::fwrite(rectangles, file.path(temp_dir, "complete_support_performance.tsv"), sep = "\t")
  data.table::fwrite(rank_bootstrap, file.path(temp_dir, "rank_bootstrap.tsv"), sep = "\t")
  data.table::fwrite(coverage, file.path(temp_dir, "coverage_performance.tsv"), sep = "\t")
  data.table::fwrite(performance_all, file.path(temp_dir, "method_unit_performance.tsv"), sep = "\t")
  data.table::fwrite(performance_sensitivity,
                     file.path(temp_dir, "method_unit_profile_sensitivity.tsv"),
                     sep = "\t")
  data.table::fwrite(hindsight, file.path(temp_dir, "hindsight_ceiling.tsv"), sep = "\t")
  data.table::fwrite(eligibility, file.path(temp_dir, "eligibility_audit.tsv"), sep = "\t")
  data.table::fwrite(join_audit, file.path(temp_dir, "prediction_join_audit.tsv"), sep = "\t")
  provenance_files <- c("provenance_preflight.log", "provenance_resolution.tsv")
  copied <- file.copy(file.path(input_dir, provenance_files),
                      file.path(temp_dir, provenance_files))
  if (length(copied) != length(provenance_files) || !all(copied)) {
    stop("Could not copy the frozen provenance preflight artifacts.")
  }
  report <- c(
    "# Fair all-method comparison", "",
    paste0("Run mode: ", opt$run_mode),
    paste0("Primary units: ", length(primary_units)),
    paste0("Within methods: ", nrow(methods)),
    paste0("Unordered pairs: ", nrow(pairs)),
    "Primary estimand: group-collapsed AUC under provenance-group-safe folds",
    "Sensitivity estimand: profile-weighted AUC for repeated-profile units",
    paste0("Pooled-testable frozen pair family: ", sum(meta$pooled_testable)),
    paste0("REML/Knapp--Hartung fits: ", sum(meta$status == "REML_KH")),
    paste0("Non-degenerate complete-support panels: ",
           sum(!panels$summary$degenerate)), "",
    "Results are descriptive until independent results review, plausibility ",
    "assessment, and canonical promotion. No universal winner is inferred here."
  )
  writeLines(report, file.path(temp_dir, "report.md"))
  artifacts <- list.files(temp_dir, full.names = TRUE)
  if (!setequal(basename(artifacts), .fair_expected_result_files()) ||
      any(dir.exists(artifacts))) {
    stop("Runner did not create the exact registered result inventory.")
  }
  output_manifest <- data.table(
    file = basename(artifacts), bytes = file.info(artifacts)$size,
    sha256 = vapply(artifacts, digest::digest, character(1L),
                    algo = "sha256", file = TRUE),
    run_mode = opt$run_mode, bootstrap_reps = opt$bootstrap_reps,
    permutation_reps = opt$permutation_reps, seed = opt$seed,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    omicselector_version = as.character(utils::packageVersion("OmicSelector")),
    package_git_commit = code_provenance$package_git_commit,
    package_git_dirty = code_provenance$package_git_dirty,
    installed_package_tree_sha256 =
      opt$expected_installed_package_tree_sha256,
    producer_script_sha256 = code_provenance$script_sha256,
    snapshot_manifest_sha256 = snapshot$manifest_sha256,
    runtime_receipt_manifest_sha256 =
      opt$expected_runtime_receipt_manifest_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256
  )
  data.table::fwrite(output_manifest, file.path(temp_dir, "output_manifest.tsv"), sep = "\t")
  .fair_validate_snapshot(
    input_dir, opt$expected_snapshot_manifest_sha256,
    require_clean = opt$run_mode == "full"
  )
  if (!identical(.fair_code_provenance(), code_provenance) ||
      !identical(.fair_package_tree_sha256(installed_root),
                 opt$expected_installed_package_tree_sha256)) {
    stop("Package source or installation changed before result promotion.")
  }
  paper3_validate_runtime_receipt(
    runtime_receipt, opt$expected_runtime_receipt_manifest_sha256,
    opt$expected_package_version, opt$expected_package_commit,
    opt$expected_installed_package_tree_sha256, opt$runtime_image,
    opt$expected_runtime_image_sha256, phase = "end",
    start_guard = runtime_guard
  )
  .fair_validate_snapshot(
    input_dir, opt$expected_snapshot_manifest_sha256,
    require_clean = opt$run_mode == "full"
  )
  .fair_validate_staged_result(temp_dir)
  if (!file.rename(temp_dir, output_dir)) stop("Could not atomically promote output bundle.")
  cat("Created fair all-method comparison bundle:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .fair_main()
