#!/usr/bin/env Rscript

# Independent, fail-closed validator for the paper-3 fair all-method bundle.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(OmicSelector)
})

.FAIRQA_SCRIPT_PATH <- local({
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
  if (!length(candidates)) stop("Could not resolve validator script path.")
  normalizePath(tail(candidates, 1L), mustWork = TRUE)
})
.FAIRQA_RUNTIME_HELPER <- file.path(
  dirname(.FAIRQA_SCRIPT_PATH), "paper3_runtime_receipt_common.R"
)
if (!file.exists(.FAIRQA_RUNTIME_HELPER)) {
  stop("Package-owned runtime-receipt helper is missing.")
}
sys.source(.FAIRQA_RUNTIME_HELPER, envir = environment())
.PAPER3RT_HELPER_PATH <- normalizePath(.FAIRQA_RUNTIME_HELPER, mustWork = TRUE)

.fairqa_args <- function(args) {
  allowed <- c(
    "input_dir", "result_dir", "output_dir",
    "expected_snapshot_manifest_sha256", "expected_package_version",
    "expected_package_commit", "expected_installed_package_tree_sha256",
    "runtime_receipt", "expected_runtime_receipt_manifest_sha256",
    "runtime_image", "expected_runtime_image_sha256", "run_mode"
  )
  out <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Arguments require --name=value: ", arg)
    key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg), fixed = TRUE)
    if (!key %in% allowed) stop("Unknown argument: --", gsub("_", "-", key))
    if (key %in% names(out)) stop("Duplicate argument: --", gsub("_", "-", key))
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  missing <- allowed[!vapply(allowed, function(key) {
    !is.null(out[[key]]) && nzchar(out[[key]])
  }, logical(1L))]
  if (length(missing)) stop("Missing required validation argument(s): ",
                            paste(missing, collapse = ", "))
  for (key in c("expected_snapshot_manifest_sha256",
                "expected_installed_package_tree_sha256",
                "expected_runtime_receipt_manifest_sha256",
                "expected_runtime_image_sha256")) {
    if (!grepl("^[0-9a-f]{64}$", out[[key]])) {
      stop(key, " is not an exact SHA-256 pin.")
    }
  }
  if (!grepl("^[0-9a-f]{40}$", out$expected_package_commit)) {
    stop("Expected package commit is not an exact Git pin.")
  }
  if (!out$run_mode %in% c("full", "smoke")) {
    stop("--run-mode must be full or smoke.")
  }
  out
}

.fairqa_sha <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

.fairqa_read <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("Required validation input is missing: ", path)
  }
  if (!file.info(path)$size) return(data.table::data.table())
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzip <- normalizePath(Sys.which("gzip"), mustWork = TRUE)
    data.table::fread(cmd = paste(shQuote(gzip), "-cd", shQuote(path)))
  } else data.table::fread(path)
}

.fairqa_package_tree_sha256 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  files <- sort(list.files(root, recursive = TRUE, all.files = TRUE,
                           full.names = TRUE, include.dirs = FALSE,
                           no.. = TRUE))
  if (!length(files)) stop("Installed OmicSelector package tree is empty.")
  prefix <- paste0(root, .Platform$file.sep)
  canonical <- normalizePath(files, mustWork = TRUE)
  relative <- substring(canonical, nchar(prefix) + 1L)
  if (any(!startsWith(canonical, prefix)) || anyDuplicated(relative) ||
      any(nzchar(Sys.readlink(files))) || any(grepl("[\r\n\t]", relative))) {
    stop("Installed OmicSelector package tree has unsafe entries.")
  }
  rows <- paste(relative, file.info(files)$size,
                vapply(canonical, .fairqa_sha, character(1L)), sep = "\t")
  digest::digest(paste(rows, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

.fairqa_code_provenance <- function() {
  script_path <- normalizePath(.FAIRQA_SCRIPT_PATH, mustWork = TRUE)
  package_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                                mustWork = TRUE)
  git <- normalizePath(Sys.which("git"), mustWork = TRUE)
  commit <- system2(git, c("-C", package_root, "rev-parse", "HEAD"),
                    stdout = TRUE, stderr = TRUE)
  status <- system2(git, c("-C", package_root, "status", "--porcelain=v1",
                           "--untracked-files=all"), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(commit, "status")) || length(commit) != 1L ||
      !grepl("^[0-9a-f]{40}$", commit) || !is.null(attr(status, "status"))) {
    stop("Could not resolve exact package provenance.")
  }
  list(script_path = script_path, package_root = package_root,
       script_sha256 = .fairqa_sha(script_path), commit = commit[[1L]],
       dirty = length(status) > 0L)
}

.fairqa_exact_inventory <- function(root, manifest, artifacts, label) {
  if (anyDuplicated(manifest$file) || !setequal(manifest$file, artifacts)) {
    stop(label, " manifest does not enumerate the exact inventory.")
  }
  present <- list.files(root, recursive = FALSE, all.files = TRUE,
                        include.dirs = TRUE, no.. = TRUE)
  manifest_name <- if (label == "Snapshot") "manifest.tsv" else
    "output_manifest.tsv"
  if (!setequal(present, c(artifacts, manifest_name)) ||
      any(dir.exists(file.path(root, present)))) {
    stop(label, " directory contains an unmanifested artifact.")
  }
  paths <- file.path(root, manifest$file)
  if (any(!file.exists(paths)) ||
      any(nzchar(Sys.readlink(paths))) ||
      any(as.numeric(manifest$bytes) != file.info(paths)$size) ||
      any(manifest$sha256 != vapply(paths, .fairqa_sha, character(1L)))) {
    stop(label, " artifact bytes differ from its exact manifest.")
  }
  invisible(TRUE)
}

.fairqa_validate_snapshot <- function(input_dir, expected_sha256) {
  artifacts <- c("methods.tsv", "units.tsv", "cells.tsv", "eligibility.tsv",
                 "pair_family.tsv", "splits.tsv", "predictions.tsv.gz",
                 "provenance_preflight.log", "provenance_resolution.tsv",
                 "source_inputs.tsv")
  path <- file.path(input_dir, "manifest.tsv")
  if (!file.exists(path) || !identical(.fairqa_sha(path), expected_sha256)) {
    stop("Input snapshot manifest differs from its exact pin.")
  }
  manifest <- .fairqa_read(path)
  required <- c("file", "bytes", "sha256", "omicselector_version",
                "package_git_commit", "package_git_dirty",
                "installed_package_tree_sha256", "producer_script_sha256",
                "snapshot_run_mode", "runtime_receipt_manifest_sha256",
                "runtime_image_sha256")
  if (!all(required %in% names(manifest)) || anyNA(manifest[, ..required])) {
    stop("Input snapshot manifest schema is incomplete.")
  }
  .fairqa_exact_inventory(input_dir, manifest, artifacts, "Snapshot")
  if (any(as.logical(manifest$package_git_dirty)) ||
      any(manifest$snapshot_run_mode != "full")) {
    stop("Input snapshot is not a clean full-mode artifact.")
  }
  list(manifest = manifest, path = path)
}

.fairqa_validate_source_inputs <- function(input_dir) {
  x <- .fairqa_read(file.path(input_dir, "source_inputs.tsv"))
  required <- c("path", "bytes", "sha256")
  if (!all(required %in% names(x)) || !nrow(x) || anyNA(x[, ..required]) ||
      anyDuplicated(x$path) || any(!grepl("^[0-9a-f]{64}$", x$sha256)) ||
      any(!file.exists(x$path)) ||
      any(as.numeric(x$bytes) != file.info(x$path)$size) ||
      any(x$sha256 != vapply(x$path, .fairqa_sha, character(1L)))) {
    stop("Snapshot source inputs differ from their exact current bytes.")
  }
  normalizePath(x$path, mustWork = TRUE)
}

.fairqa_result_artifacts <- function() c(
  "pairwise_cohort.tsv", "pairwise_strata.tsv",
  "pairwise_profile_sensitivity.tsv", "pairwise_meta.tsv", "pairwise_loco.tsv",
  "pairwise_heterogeneity.tsv", "complete_support_panels.tsv",
  "complete_support_performance.tsv", "rank_bootstrap.tsv",
  "coverage_performance.tsv", "method_unit_performance.tsv",
  "method_unit_profile_sensitivity.tsv", "hindsight_ceiling.tsv",
  "eligibility_audit.tsv", "prediction_join_audit.tsv",
  "provenance_preflight.log", "provenance_resolution.tsv", "report.md"
)

.fairqa_validate_bundle_manifest <- function(
    result_dir, expected_snapshot_sha, expected_version, expected_commit,
    expected_installed_tree, expected_runtime_receipt,
    expected_runtime_image, run_mode, package_root) {
  path <- file.path(result_dir, "output_manifest.tsv")
  manifest <- .fairqa_read(path)
  required <- c(
    "file", "bytes", "sha256", "run_mode", "bootstrap_reps",
    "permutation_reps", "seed", "omicselector_version",
    "package_git_commit", "package_git_dirty", "producer_script_sha256",
    "installed_package_tree_sha256", "snapshot_manifest_sha256",
    "runtime_receipt_manifest_sha256", "runtime_image_sha256"
  )
  if (!all(required %in% names(manifest)) || anyNA(manifest[, ..required])) {
    stop("Fair-comparison output manifest schema is incomplete.")
  }
  .fairqa_exact_inventory(result_dir, manifest, .fairqa_result_artifacts(),
                          "Result")
  runner <- file.path(package_root, "inst", "scripts",
                      "paper3_fair_method_comparison.R")
  if (any(manifest$run_mode != run_mode) ||
      any(manifest$omicselector_version != expected_version) ||
      any(manifest$package_git_commit != expected_commit) ||
      any(as.logical(manifest$package_git_dirty)) ||
      any(manifest$producer_script_sha256 != .fairqa_sha(runner)) ||
      any(manifest$installed_package_tree_sha256 != expected_installed_tree) ||
      any(manifest$snapshot_manifest_sha256 != expected_snapshot_sha) ||
      any(manifest$runtime_receipt_manifest_sha256 !=
            expected_runtime_receipt) ||
      any(manifest$runtime_image_sha256 != expected_runtime_image)) {
    stop("Result manifest differs from exact code/snapshot/runtime pins.")
  }
  if (run_mode == "full" &&
      (any(as.integer(manifest$bootstrap_reps) != 10000L) ||
       any(as.integer(manifest$permutation_reps) != 1000L) ||
       any(as.integer(manifest$seed) != 20260718L))) {
    stop("Full result lacks the registered resampling contract.")
  }
  if (length(unique(as.integer(manifest$bootstrap_reps))) != 1L ||
      length(unique(as.integer(manifest$permutation_reps))) != 1L ||
      length(unique(as.integer(manifest$seed))) != 1L) {
    stop("Result manifest carries inconsistent resampling pins.")
  }
  list(manifest = manifest, path = path,
       bootstrap_reps = unique(as.integer(manifest$bootstrap_reps)),
       permutation_reps = unique(as.integer(manifest$permutation_reps)),
       seed = unique(as.integer(manifest$seed)))
}

.fairqa_assert_disjoint_output <- function(output, roots) {
  output <- file.path(normalizePath(dirname(output), mustWork = TRUE),
                      basename(output))
  roots <- unique(vapply(roots, normalizePath, character(1L), mustWork = TRUE))
  bad <- vapply(roots, function(root) {
    identical(output, root) || startsWith(output, paste0(root, .Platform$file.sep)) ||
      startsWith(root, paste0(output, .Platform$file.sep))
  }, logical(1L))
  if (any(bad)) stop("Validation output is not disjoint from an immutable root.")
  invisible(output)
}

.fairqa_auc <- function(y, score) {
  y <- as.integer(y)
  if (length(unique(y)) != 2L) return(NA_real_)
  pos <- score[y == 1L]; neg <- score[y == 0L]
  ranks <- rank(c(pos, neg), ties.method = "average")
  (sum(ranks[seq_along(pos)]) - length(pos) * (length(pos) + 1L) / 2) /
    (length(pos) * length(neg))
}

.fairqa_precompute_strata <- function(predictions, splits, level) {
  if (level == "group") {
    d <- predictions[, .(y = unique(y), score = mean(score)),
      by = .(method_id, unit_id, seed, fold, group_id)]
    out <- d[, .(auc = .fairqa_auc(y, score), n_test = .N),
             by = .(method_id, unit_id, seed, fold)]
    count <- "n_test_groups"
  } else {
    out <- predictions[, .(auc = .fairqa_auc(y, score), n_test = .N),
                       by = .(method_id, unit_id, seed, fold)]
    count <- "n_test_profiles"
  }
  registered <- merge(
    splits[, .(unit_id, seed, fold, split_id,
               registered_test = get(count))],
    splits[, .(total_observations = sum(get(count))), by = .(unit_id, seed)],
    by = c("unit_id", "seed"), all.x = TRUE, sort = FALSE
  )
  out <- merge(out, registered, by = c("unit_id", "seed", "fold"),
               all.x = TRUE, sort = FALSE)
  if (anyNA(out) || any(out$n_test != out$registered_test) ||
      any(!is.finite(out$auc))) stop("Snapshot strata do not reproduce AUC.")
  out[, `:=`(n_train = total_observations - n_test,
             analysis_level = level, stratum_id = paste(seed, fold, sep = "::"))]
  data.table::setkey(out, method_id, unit_id, seed, fold)
  out
}

.fairqa_performance <- function(strata, eligibility, splits) {
  rows <- lapply(seq_len(nrow(eligibility)), function(i) {
    e <- eligibility[i]
    if (!isTRUE(e$eligible)) return(data.table(
      method_id = e$method_id, unit_id = e$unit_id, eligible = FALSE,
      auc = NA_real_, n_strata = 0L, reason = e$reason
    ))
    d <- strata[list(e$method_id, e$unit_id), nomatch = 0L]
    expected <- splits[list(e$unit_id), paste(seed, fold, sep = "::")]
    if (nrow(d) != length(expected) || !setequal(d$stratum_id, expected)) {
      stop("Eligible method/unit lacks exact strata.")
    }
    data.table(method_id = e$method_id, unit_id = e$unit_id, eligible = TRUE,
               auc = mean(d$auc), n_strata = nrow(d), reason = "eligible")
  })
  data.table::rbindlist(rows)
}

.fairqa_corrected <- function(effect, n_test, n_train) {
  m <- length(effect); ratio <- n_test / n_train
  correction <- 1 / m + mean(ratio)
  estimate <- mean(effect)
  se <- sqrt(max(0, correction * stats::var(effect)))
  df <- m - 1L; critical <- stats::qt(0.975, df)
  available <- is.finite(se) && se > 0
  list(estimate = estimate, se = se,
       ci_low = estimate - critical * se, ci_high = estimate + critical * se,
       p_zero = if (available) 2 * stats::pt(-abs(estimate / se), df) else NA_real_,
       p_above_margin = if (available)
         stats::pt((estimate - 0.05) / se, df, lower.tail = FALSE) else NA_real_,
       p_below_minus_margin = if (available)
         stats::pt((estimate + 0.05) / se, df, lower.tail = TRUE) else NA_real_,
       expected_m = m, observed_m = m, correction = correction,
       mean_test_train_ratio = mean(ratio), available = available)
}

.fairqa_pair_family <- function(pairs, units, strata, eligibility, splits,
                                prediction_counts) {
  summaries <- strata_rows <- joins <- list(); index <- 0L
  level <- unique(strata$analysis_level)
  for (p in seq_len(nrow(pairs))) for (unit in units) {
    index <- index + 1L; pair <- pairs[p]
    ea <- eligibility[list(pair$method_a, unit), nomatch = 0L]
    eb <- eligibility[list(pair$method_b, unit), nomatch = 0L]
    base <- data.table(
      pair_id = pair$pair_id, method_a = pair$method_a,
      method_b = pair$method_b, unit_id = unit, analysis_level = level,
      common_support = FALSE, estimate = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_zero = NA_real_,
      p_above_margin = NA_real_, p_below_minus_margin = NA_real_,
      expected_m = NA_integer_, observed_m = NA_integer_,
      correction = NA_real_, mean_test_train_ratio = NA_real_, reason = ""
    )
    if (!isTRUE(ea$eligible) || !isTRUE(eb$eligible)) {
      base$reason <- paste0("method_a:", ea$reason, " | method_b:", eb$reason)
      summaries[[index]] <- base
      joins[[index]] <- data.table(pair_id = pair$pair_id, unit_id = unit,
        pass = FALSE, reason = "not_common_eligible", n_a = 0L, n_b = 0L,
        n_join = 0L)
      next
    }
    a <- strata[list(pair$method_a, unit), nomatch = 0L]
    b <- strata[list(pair$method_b, unit), nomatch = 0L]
    cols <- c("unit_id", "seed", "fold", "split_id", "stratum_id",
              "n_test", "n_train", "analysis_level")
    joined <- merge(a[, c(cols, "auc"), with = FALSE],
                    b[, c(cols, "auc"), with = FALSE], by = cols,
                    suffixes = c("_a", "_b"), sort = FALSE)
    expected <- splits[list(unit), paste(seed, fold, sep = "::")]
    ca <- prediction_counts[list(pair$method_a, unit), n_profiles]
    cb <- prediction_counts[list(pair$method_b, unit), n_profiles]
    exact <- length(ca) == 1L && length(cb) == 1L && ca == cb &&
      nrow(joined) == nrow(a) && nrow(joined) == nrow(b) &&
      nrow(joined) == length(expected) && setequal(joined$stratum_id, expected)
    if (!exact) stop("Independent exact pair join failed.")
    joins[[index]] <- data.table(pair_id = pair$pair_id, unit_id = unit,
      pass = TRUE, reason = "exact", n_a = ca, n_b = cb, n_join = ca)
    joined[, effect := auc_a - auc_b]
    s <- .fairqa_corrected(joined$effect, joined$n_test, joined$n_train)
    base[, `:=`(common_support = TRUE, estimate = s$estimate, se = s$se,
      ci_low = s$ci_low, ci_high = s$ci_high, p_zero = s$p_zero,
      p_above_margin = s$p_above_margin,
      p_below_minus_margin = s$p_below_minus_margin,
      expected_m = s$expected_m, observed_m = s$observed_m,
      correction = s$correction,
      mean_test_train_ratio = s$mean_test_train_ratio,
      reason = if (s$available) "evaluable" else "zero_corrected_se")]
    summaries[[index]] <- base
    strata_rows[[length(strata_rows) + 1L]] <- joined[, .(
      stratum_id, auc_a, auc_b, effect, n_test, n_train, analysis_level,
      pair_id = pair$pair_id, method_a = pair$method_a,
      method_b = pair$method_b, unit_id = unit)]
  }
  list(pairwise = data.table::rbindlist(summaries, fill = TRUE),
       strata = data.table::rbindlist(strata_rows, fill = TRUE),
       joins = data.table::rbindlist(joins, fill = TRUE))
}

.fairqa_meta_one <- function(d, pair, support) {
  row <- data.table(pair_id = pair$pair_id, method_a = pair$method_a,
    method_b = pair$method_b, k_support = length(support), k_meta = 0L,
    support_units = paste(support, collapse = ";"), estimate = NA_real_,
    se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, tau2 = NA_real_,
    i2 = NA_real_, p_zero = NA_real_, p_above_margin = NA_real_,
    p_below_minus_margin = NA_real_, df = NA_real_,
    pooled_testable = length(support) >= 5L,
    status = if (length(support) >= 5L) "meta_unavailable" else
      "descriptive_k_lt_5",
    fit_reason = if (length(support) >= 5L) "not_fitted" else "k_support_lt_5")
  use <- d[common_support & is.finite(estimate) & is.finite(se) & se > 0]
  row$k_meta <- nrow(use)
  if (nrow(use) < 5L) {
    row$fit_reason <- "fewer_than_five_finite_unit_effects"; return(row)
  }
  fit <- tryCatch(withCallingHandlers(
    metafor::rma.uni(yi = use$estimate, sei = use$se, method = "REML",
                     test = "knha"),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)),
    error = identity)
  if (inherits(fit, "error")) {
    row$fit_reason <- paste0("rma_failed:", conditionMessage(fit)); return(row)
  }
  estimate_value <- as.numeric(fit$b[[1L]])
  se_value <- as.numeric(fit$se[[1L]])
  df_value <- nrow(use) - 1L
  row[, `:=`(estimate = estimate_value, se = se_value,
    ci_low = as.numeric(fit$ci.lb[[1L]]),
    ci_high = as.numeric(fit$ci.ub[[1L]]), tau2 = as.numeric(fit$tau2),
    i2 = as.numeric(fit$I2), p_zero = as.numeric(fit$pval[[1L]]),
    p_above_margin = stats::pt((estimate_value - .05) / se_value, df_value,
                               lower.tail = FALSE),
    p_below_minus_margin = stats::pt((estimate_value + .05) / se_value,
                                     df_value,
                                     lower.tail = TRUE),
    df = df_value, status = "REML_KH", fit_reason = "eligible")]
  row
}

.fairqa_meta_all <- function(pairwise, pairs, frozen) {
  out <- data.table::rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i]; d <- pairwise[pair_id == pair$pair_id]
    .fairqa_meta_one(d, pair, sort(d[common_support == TRUE, unit_id]))
  }))
  frozen <- frozen[match(out$pair_id, pair_id)]
  if (any(out$k_support != frozen$k_support) ||
      any(out$support_units != frozen$support_units)) {
    stop("Runtime support differs from the frozen pair family.")
  }
  out[, pooled_testable := as.logical(frozen$pooled_testable)]
  ids <- out[pooled_testable == TRUE, pair_id]
  if (length(ids)) {
    p0 <- out[pooled_testable == TRUE, p_zero]
    directional <- c(rbind(out[pooled_testable == TRUE, p_above_margin],
                           out[pooled_testable == TRUE,
                               p_below_minus_margin]))
    adjusted <- stats::p.adjust(directional, "BY", n = 2L * length(ids))
    out[pooled_testable == TRUE, `:=`(
      p_zero_by = stats::p.adjust(p0, "BY", n = length(ids)),
      p_above_margin_by = adjusted[seq(1L, length(adjusted), 2L)],
      p_below_minus_margin_by = adjusted[seq(2L, length(adjusted), 2L)],
      zero_family_n = length(ids), relevance_family_n = 2L * length(ids))]
    out[pooled_testable == FALSE, `:=`(
      p_zero_by = NA_real_, p_above_margin_by = NA_real_,
      p_below_minus_margin_by = NA_real_, zero_family_n = NA_integer_,
      relevance_family_n = NA_integer_)]
  } else out[, `:=`(p_zero_by = NA_real_, p_above_margin_by = NA_real_,
    p_below_minus_margin_by = NA_real_, zero_family_n = 0L,
    relevance_family_n = 0L)]
  out
}

.fairqa_loco <- function(pairwise, pairs) {
  rows <- list()
  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i]; support <- pairwise[pair_id == pair$pair_id & common_support]
    finite <- support[is.finite(estimate) & is.finite(se) & se > 0]
    for (unit in support$unit_id) {
      keep <- finite[unit_id != unit]
      row <- data.table(pair_id = pair$pair_id, method_a = pair$method_a,
        method_b = pair$method_b, omitted_unit = unit,
        k_support = nrow(support), k_finite_before_omission = nrow(finite),
        k = nrow(keep), estimate = NA_real_, se = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, tau2 = NA_real_, i2 = NA_real_, p_zero = NA_real_,
        status = "not_evaluable",
        fit_reason = "fewer_than_five_finite_effects_after_omission")
      if (nrow(keep) >= 5L) {
        fit <- tryCatch(withCallingHandlers(metafor::rma.uni(
          yi = keep$estimate, sei = keep$se, method = "REML", test = "knha"),
          warning = function(w) stop(conditionMessage(w), call. = FALSE)),
          error = identity)
        if (inherits(fit, "error")) {
          row$status <- "fit_failed"
          row$fit_reason <- paste0("rma_failed:", conditionMessage(fit))
        } else row[, `:=`(estimate = as.numeric(fit$b[[1L]]),
          se = as.numeric(fit$se[[1L]]), ci_low = as.numeric(fit$ci.lb[[1L]]),
          ci_high = as.numeric(fit$ci.ub[[1L]]), tau2 = as.numeric(fit$tau2),
          i2 = as.numeric(fit$I2), p_zero = as.numeric(fit$pval[[1L]]),
          status = "REML_KH", fit_reason = "eligible")]
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

.fairqa_meta_strata <- function(pairwise, pairs, units, field, inferential) {
  x <- merge(pairwise, units[, c("unit_id", field), with = FALSE],
             by = "unit_id", all.x = TRUE, sort = FALSE)
  rows <- list()
  for (value in sort(unique(x[[field]]))) for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i]; d <- x[pair_id == pair$pair_id & get(field) == value]
    support <- sort(d[common_support == TRUE, unit_id])
    supported <- d[common_support & is.finite(estimate), estimate]
    row <- if (inferential) .fairqa_meta_one(d, pair, support) else
      .fairqa_meta_one(d[0L], pair, support)
    if (!inferential) row[, `:=`(k_meta = length(supported), estimate = NA_real_,
      se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, tau2 = NA_real_,
      i2 = NA_real_, p_zero = NA_real_, p_above_margin = NA_real_,
      p_below_minus_margin = NA_real_, df = NA_real_, status = "descriptive_only",
      fit_reason = "biospecimen_inference_not_prespecified")]
    row[, `:=`(descriptive_mean = if (length(supported)) mean(supported) else NA_real_,
      descriptive_median = if (length(supported)) median(supported) else NA_real_,
      descriptive_min = if (length(supported)) min(supported) else NA_real_,
      descriptive_max = if (length(supported)) max(supported) else NA_real_,
      inference_policy = if (inferential) "REML_KH_when_k_at_least_5" else
        "descriptive_only", stratum_variable = field, stratum_value = value)]
    rows[[length(rows) + 1L]] <- row
  }
  data.table::rbindlist(rows, fill = TRUE)
}

.fairqa_panels <- function(eligibility, methods, primary_units) {
  ids <- methods$method_id; thresholds <- c(.50, .80, .95)
  grid <- eligibility[unit_id %in% primary_units]
  matrix <- dcast(grid, method_id ~ unit_id, value.var = "eligible")
  matrix <- as.matrix(matrix[match(ids, method_id), ..primary_units])
  coverage_n <- rowSums(matrix)
  summaries <- memberships <- unit_rows <- list()
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[[i]]; required <- ceiling(threshold * length(primary_units))
    included <- coverage_n >= required; retained <- ids[included]
    included_unit <- if (length(retained))
      colSums(matrix[included, , drop = FALSE]) == length(retained) else
        rep(FALSE, length(primary_units))
    degenerate <- sum(included) < 5L || sum(included_unit) < 5L
    summaries[[i]] <- data.table(threshold = threshold,
      required_eligible_units = required, total_primary_units = length(primary_units),
      n_methods = sum(included), n_units = sum(included_unit), min_methods = 5L,
      min_units = 5L, degenerate = degenerate, rank_allowed = !degenerate,
      item_type = "summary")
    memberships[[i]] <- data.table(threshold = threshold, method_id = ids,
      roster_order = seq_along(ids), eligible_units = as.integer(coverage_n),
      total_primary_units = length(primary_units),
      coverage = coverage_n / length(primary_units), included_method = included,
      exclusion_reason = ifelse(included, NA_character_,
                                "below_coverage_threshold"), item_type = "method")
    reason <- if (!length(retained)) rep("no_method_meets_threshold",
                                         length(primary_units)) else
      ifelse(included_unit, NA_character_, "ineligible_for_retained_method")
    unit_rows[[i]] <- data.table(threshold = threshold, unit_id = primary_units,
      included_unit = included_unit, exclusion_reason = reason, item_type = "unit")
  }
  data.table::rbindlist(c(summaries, memberships, unit_rows), fill = TRUE)
}

.fairqa_rank_bootstrap <- function(performance, panels, reps, seed) {
  rows <- list(); summary <- panels[item_type == "summary" & !degenerate]
  for (threshold_value in summary$threshold) {
    methods <- panels[item_type == "method" & threshold == threshold_value &
                        included_method, method_id]
    units <- panels[item_type == "unit" & threshold == threshold_value &
                      included_unit, unit_id]
    rectangle <- performance[method_id %in% methods & unit_id %in% units]
    ranked <- copy(rectangle)[
      , observed_within_unit_rank := rank(-auc, ties.method = "average"),
      by = unit_id
    ]
    wide <- dcast(ranked, unit_id ~ method_id,
                  value.var = "observed_within_unit_rank")
    if (anyNA(wide)) stop("Complete-support rank rectangle contains NA.")
    observed <- ranked[, .(
      observed_average_rank = mean(observed_within_unit_rank)
    ), by = method_id]
    set.seed(seed + as.integer(round(threshold_value * 1000)))
    draws <- matrix(NA_real_, reps, length(methods), dimnames = list(NULL, methods))
    top1 <- top3 <- top5 <- draws == 0
    for (b in seq_len(reps)) {
      sampled <- sample(seq_len(nrow(wide)), nrow(wide), replace = TRUE)
      average_ranks <- colMeans(as.matrix(wide[sampled, ..methods]))
      draws[b, ] <- average_ranks
      inclusive <- rank(average_ranks, ties.method = "min")
      top1[b, ] <- inclusive <= 1L; top3[b, ] <- inclusive <= 3L
      top5[b, ] <- inclusive <= 5L
    }
    for (method in methods) rows[[length(rows) + 1L]] <- data.table(
      threshold = threshold_value, method_id = method,
      rank_median = median(draws[, method]),
      rank_ci_low = unname(quantile(draws[, method], .025, type = 6L)),
      rank_ci_high = unname(quantile(draws[, method], .975, type = 6L)),
      observed_average_rank = observed[method_id == method,
                                       observed_average_rank],
      p_top1 = mean(top1[, method]), p_top3 = mean(top3[, method]),
      p_top5 = mean(top5[, method]), bootstrap_reps = reps, seed = seed,
      tie_policy = "inclusive_min_rank_for_topk;average_rank_for_intervals")
  }
  data.table::rbindlist(rows, fill = TRUE)
}

.fairqa_hindsight <- function(predictions, performance, panels, reps, seed) {
  summary <- panels[item_type == "summary" & threshold == .95]
  if (!nrow(summary) || summary$degenerate[[1L]]) return(data.table())
  methods <- panels[item_type == "method" & threshold == .95 &
                      included_method, method_id]
  units <- panels[item_type == "unit" & threshold == .95 & included_unit,
                  unit_id]
  rows <- list()
  for (u in seq_along(units)) {
    unit <- units[[u]]; d <- predictions[unit_id == unit & method_id %in% methods]
    d[, stratum_id := paste(seed, fold, sep = "::")]
    labels <- unique(d[, .(stratum_id, group_id, y)])
    observed_max <- max(performance[unit_id == unit & method_id %in% methods, auc])
    set.seed(seed + u * 1009L); null <- numeric(); attempts <- 0L
    while (length(null) < reps && attempts < reps * 100L) {
      attempts <- attempts + 1L
      perm <- labels[, .(group_id, y_perm = sample(y)), by = stratum_id]
      x <- merge(d[, .(method_id, stratum_id, group_id, sample_id, score)],
                 perm, by = c("stratum_id", "group_id"), sort = FALSE)
      collapsed <- x[, .(y_perm = unique(y_perm), score = mean(score)),
                     by = .(method_id, stratum_id, group_id)]
      auc <- collapsed[, .(auc = .fairqa_auc(y_perm, score)),
                       by = .(method_id, stratum_id)]
      if (any(!is.finite(auc$auc))) next
      means <- auc[, .(auc = mean(auc)), by = method_id]
      if (nrow(means) == length(methods)) null <- c(null, max(means$auc))
    }
    valid_fraction <- length(null) / attempts
    failed <- length(null) < reps || valid_fraction < .8
    rows[[u]] <- data.table(unit_id = unit, observed_max = observed_max,
      null_median = if (failed) NA_real_ else median(null),
      null_ci_low = if (failed) NA_real_ else unname(quantile(null, .025, type = 6L)),
      null_ci_high = if (failed) NA_real_ else unname(quantile(null, .975, type = 6L)),
      headroom = if (failed) NA_real_ else observed_max - median(null),
      valid_permutations = length(null), attempts = attempts,
      valid_fraction = valid_fraction,
      status = if (failed) "permutation_guard_failed" else "evaluable",
      permutation_reps = reps, seed = seed, panel_threshold = .95,
      n_methods = length(methods), analysis_level = "group")
  }
  rbindlist(rows)
}

.fairqa_validate_fold_contract <- function(units, splits) {
  if (nrow(units) != 34L || anyNA(units$outer_k) ||
      any(units$outer_k != as.integer(units$outer_k)) ||
      any(!units$outer_k %in% c(3L, 5L))) {
    stop("Unit outer_k is not the exact 34-unit frozen contract.")
  }
  contract <- merge(splits[, .(
    n = .N, folds = paste(sort(unique(fold)), collapse = ";")
  ), by = .(unit_id, seed)], units[, .(unit_id, outer_k)],
  by = "unit_id", all = TRUE)
  expected_seeds <- c(101L, 202L, 303L, 404L, 505L)
  if (nrow(contract) != 34L * 5L ||
      !setequal(unique(contract$seed), expected_seeds) || anyNA(contract) ||
      any(contract$n != contract$outer_k) ||
      any(contract$folds != vapply(contract$outer_k, function(k)
        paste(seq_len(k), collapse = ";"), character(1L)))) {
    stop("Snapshot does not contain exact 34x5 outer_k folds.")
  }
  invisible(TRUE)
}

.fairqa_validate_prediction_splits <- function(predictions, splits) {
  split_columns <- c("unit_id", "seed", "fold", "split_id")
  prediction_columns <- c("method_id", split_columns, "group_id")
  if (!all(split_columns %in% names(splits)) ||
      !all(prediction_columns %in% names(predictions)) ||
      anyNA(splits[, ..split_columns]) ||
      anyNA(predictions[, ..prediction_columns])) {
    stop("Prediction/split identity columns are incomplete.")
  }
  registered <- unique(splits[, ..split_columns])
  observed_method_units <- unique(predictions[, .(method_id, unit_id)])
  expected <- merge(observed_method_units, registered, by = "unit_id",
                    all.x = TRUE, allow.cartesian = TRUE, sort = FALSE)
  setcolorder(expected, c("method_id", split_columns))
  observed <- unique(predictions[, c("method_id", split_columns), with = FALSE])
  setcolorder(observed, names(expected))
  if (!data.table::fsetequal(observed, expected)) {
    stop("Prediction split rows differ from the registered split contract.")
  }
  per_method <- predictions[, uniqueN(split_id),
                            by = .(method_id, unit_id, seed)]
  canonical <- unique(predictions[, .(unit_id, seed, split_id)])[,
    uniqueN(split_id), by = .(unit_id, seed)]
  grouped <- unique(predictions[, .(unit_id, seed, group_id, fold)])[,
    uniqueN(fold), by = .(unit_id, seed, group_id)]
  if (any(per_method$V1 != 1L) || any(canonical$V1 != 1L)) {
    stop("Predictions do not share one canonical split ID across methods.")
  }
  if (any(grouped$V1 != 1L)) {
    stop("A provenance group occupies more than one held-out fold.")
  }
  invisible(TRUE)
}

.fairqa_rectangles <- function(performance, panels) {
  out <- rbindlist(lapply(c(.50, .80, .95), function(threshold_value) {
    selected_methods <- panels[item_type == "method" & threshold == threshold_value &
                                 included_method, method_id]
    selected_units <- panels[item_type == "unit" & threshold == threshold_value &
                               included_unit, unit_id]
    performance[method_id %in% selected_methods & unit_id %in% selected_units][
      , threshold := threshold_value]
  }), fill = TRUE)
  if (nrow(out)) out[, within_unit_rank := rank(-auc, ties.method = "average"),
                     by = .(threshold, unit_id)]
  out
}

.fairqa_coverage <- function(methods, eligibility, performance, predictions,
                             splits, primary_units) {
  coverage <- merge(eligibility[unit_id %in% primary_units, .(
    eligible_units = sum(eligible), eligibility_coverage = mean(eligible),
    eligibility_units = paste(unit_id[eligible], collapse = ";")), by = method_id],
    performance[eligible == TRUE, .(median_auc = median(auc),
      performance_units = paste(unit_id, collapse = ";")), by = method_id],
    by = "method_id", all.x = TRUE)
  expected_rows <- merge(eligibility[eligible & unit_id %in% primary_units,
    .(method_id, unit_id)], splits[, .(expected_rows = sum(n_test_profiles)),
                                  by = unit_id], by = "unit_id")[,
    .(expected_prediction_rows = sum(expected_rows)), by = method_id]
  actual_rows <- predictions[unit_id %in% primary_units,
                             .(prediction_rows = .N), by = method_id]
  coverage <- Reduce(function(x, y) merge(x, y, by = "method_id", all.x = TRUE,
    sort = FALSE), list(methods[, .(method_id, roster_order, family, role, tier,
      dep_route, pkg_status)], coverage, expected_rows, actual_rows))
  coverage[is.na(expected_prediction_rows), expected_prediction_rows := 0L]
  coverage[is.na(prediction_rows), prediction_rows := 0L]
  coverage[, `:=`(prediction_completeness = fifelse(expected_prediction_rows > 0,
    prediction_rows / expected_prediction_rows, 1),
    prediction_complete = prediction_rows == expected_prediction_rows,
    deployment_status = fcase(pkg_status %in% c("present", "primitive"),
      "package_native", pkg_status == "new", "package_added_primary",
      pkg_status == "new-secondary", "package_added_secondary",
      default = "unknown"))]
  coverage[, pareto_optimal := vapply(seq_len(.N), function(i) {
    if (!is.finite(median_auc[[i]]) || !is.finite(eligibility_coverage[[i]]))
      return(FALSE)
    !any(seq_len(.N) != i & is.finite(median_auc) &
      eligibility_coverage >= eligibility_coverage[[i]] &
      median_auc >= median_auc[[i]] &
      (eligibility_coverage > eligibility_coverage[[i]] |
       median_auc > median_auc[[i]]))
  }, logical(1L))]
  coverage
}

.fairqa_compare <- function(actual, expected, key, label, tolerance = 1e-12) {
  actual <- copy(as.data.table(actual)); expected <- copy(as.data.table(expected))
  if (!identical(names(actual), names(expected)) ||
      nrow(actual) != nrow(expected)) {
    stop(label, " schema or dimensions differ from reconstruction.")
  }
  if (!nrow(expected)) return(invisible(TRUE))
  if (any(!key %in% names(actual))) {
    stop(label, " lacks its exact comparison key.")
  }
  setorderv(actual, key); setorderv(expected, key)
  for (column in names(expected)) {
    a <- actual[[column]]; e <- expected[[column]]
    if (is.numeric(e) || is.integer(e)) {
      a <- suppressWarnings(as.numeric(a)); e <- as.numeric(e)
      same <- (is.na(a) & is.na(e)) |
        (is.finite(a) & is.finite(e) & abs(a - e) <= tolerance)
    } else same <- (is.na(a) & is.na(e)) | as.character(a) == as.character(e)
    if (length(same) != length(e) || anyNA(same) || any(!same)) {
      stop(label, " column differs from reconstruction: ", column)
    }
  }
  invisible(TRUE)
}

.fairqa_validate_science <- function(input_dir, result_dir, run_mode,
                                     bootstrap_reps, permutation_reps, seed) {
  methods <- .fairqa_read(file.path(input_dir, "methods.tsv"))
  units <- .fairqa_read(file.path(input_dir, "units.tsv"))
  cells <- .fairqa_read(file.path(input_dir, "cells.tsv"))
  eligibility <- .fairqa_read(file.path(input_dir, "eligibility.tsv"))
  frozen <- .fairqa_read(file.path(input_dir, "pair_family.tsv"))
  splits <- .fairqa_read(file.path(input_dir, "splits.tsv"))
  predictions <- .fairqa_read(file.path(input_dir, "predictions.tsv.gz"))
  predictions[, score := suppressWarnings(as.numeric(score))]
  roster <- as.data.table(singlesample_method_roster())[estimand == "within"]
  pairs <- as.data.table(singlesample_within_method_pairs(roster))
  if (nrow(methods) != 50L || !identical(methods$method_id, roster$method_id) ||
      nrow(units) != 34L || sum(units$primary_unit) != 33L ||
      nrow(cells) != 8500L || nrow(eligibility) != 1700L ||
      nrow(frozen) != 1225L || !any(eligibility$eligible) ||
      !identical(frozen$pair_id, pairs$pair_id)) {
    stop("Snapshot dimensions/roster are not the non-empty frozen design.")
  }
  .fairqa_validate_fold_contract(units, splits)
  required_prediction <- c("method_id", "unit_id", "seed", "fold", "split_id",
                           "sample_id", "group_id", "y", "score")
  if (!all(required_prediction %in% names(predictions)) || !nrow(predictions) ||
      anyNA(predictions[, ..required_prediction]) ||
      any(!is.finite(predictions$score)) || any(!predictions$y %in% c(0L, 1L)) ||
      anyDuplicated(predictions[, paste(method_id, unit_id, seed, fold,
                                       sample_id, sep = "\r")])) {
    stop("Snapshot predictions are empty, incomplete, or duplicated.")
  }
  .fairqa_validate_prediction_splits(predictions, splits)
  group_labels <- predictions[, uniqueN(y), by = .(unit_id, group_id)]
  if (any(group_labels$V1 != 1L)) stop("A provenance group has discordant labels.")
  data.table::setkey(predictions, method_id, unit_id)
  data.table::setkey(eligibility, method_id, unit_id)
  data.table::setkey(splits, unit_id, seed, fold)
  counts <- predictions[, .(n_profiles = .N), by = .(method_id, unit_id)]
  data.table::setkey(counts, method_id, unit_id)
  expected_counts <- merge(eligibility[, .(method_id, unit_id, eligible)],
    splits[, .(expected = sum(n_test_profiles)), by = unit_id], by = "unit_id")
  observed_counts <- merge(expected_counts, counts, by = c("method_id", "unit_id"),
                           all.x = TRUE)
  observed_counts[is.na(n_profiles), n_profiles := 0L]
  if (any(observed_counts[eligible == TRUE, n_profiles != expected]) ||
      any(observed_counts[eligible == FALSE, n_profiles != 0L])) {
    stop("Snapshot prediction accounting differs from eligibility/splits.")
  }
  surfaces <- predictions[, {
    surface <- unique(.SD)
    setorderv(surface, c("sample_id", "group_id", "y"))
    .(surface_sha256 = digest::digest(
      paste(do.call(paste, c(surface, sep = "\t")), collapse = "\n"),
      algo = "sha256", serialize = FALSE
    ))
  }, by = .(method_id, unit_id, seed, fold),
  .SDcols = c("sample_id", "group_id", "y")]
  surface_check <- surfaces[, .(
    n_methods = uniqueN(method_id), n_surfaces = uniqueN(surface_sha256)
  ), by = .(unit_id, seed, fold)]
  eligible_methods <- eligibility[eligible == TRUE, .(n_expected = .N),
                                  by = unit_id]
  surface_check <- merge(surface_check, eligible_methods, by = "unit_id",
                         all.x = TRUE)
  if (anyNA(surface_check$n_expected) ||
      any(surface_check$n_methods != surface_check$n_expected) ||
      any(surface_check$n_surfaces != 1L)) {
    stop("Eligible methods do not share one exact held-out profile surface.")
  }
  profile <- .fairqa_precompute_strata(predictions, splits, "profile")
  group <- .fairqa_precompute_strata(predictions, splits, "group")
  performance <- .fairqa_performance(group, eligibility, splits)
  performance_profile <- .fairqa_performance(profile, eligibility, splits)
  actual_performance <- .fairqa_read(file.path(result_dir,
                                               "method_unit_performance.tsv"))
  .fairqa_compare(actual_performance, performance,
                  c("method_id", "unit_id"),
                  "Group-collapsed method-unit performance")
  sensitivity_performance <- merge(performance[, .(method_id, unit_id, eligible,
    group_collapsed_auc = auc, n_strata)],
    performance_profile[, .(method_id, unit_id, profile_auc = auc)],
    by = c("method_id", "unit_id"), all = TRUE,
    sort = FALSE)
  sensitivity_performance[, difference_profile_minus_group :=
                            profile_auc - group_collapsed_auc]
  sensitivity_performance <- merge(sensitivity_performance,
    splits[, .(repeated_profiles = any(n_test_profiles > n_test_groups)),
           by = unit_id], by = "unit_id", all.x = TRUE, sort = FALSE)
  setcolorder(sensitivity_performance, c(
    "method_id", "unit_id", "eligible", "group_collapsed_auc", "n_strata",
    "profile_auc", "difference_profile_minus_group", "repeated_profiles"
  ))
  .fairqa_compare(.fairqa_read(file.path(result_dir,
    "method_unit_profile_sensitivity.tsv")), sensitivity_performance,
    c("method_id", "unit_id"), "Method-unit profile sensitivity")
  all_units <- units$unit_id
  primary_units <- units[primary_unit == TRUE, unit_id]
  pf <- .fairqa_pair_family(pairs, all_units, profile, eligibility, splits, counts)
  gf <- .fairqa_pair_family(pairs, all_units, group, eligibility, splits, counts)
  gf$pairwise <- merge(gf$pairwise, units[, .(unit_id, primary_unit)],
                       by = "unit_id", all.x = TRUE, sort = FALSE)
  if (!any(gf$pairwise$common_support) || !nrow(gf$strata)) {
    stop("Fair comparison has no common-support scientific evidence.")
  }
  .fairqa_compare(.fairqa_read(file.path(result_dir, "pairwise_cohort.tsv")),
    gf$pairwise, c("pair_id", "unit_id"),
    "Pairwise group-collapsed summaries")
  .fairqa_compare(.fairqa_read(file.path(result_dir, "pairwise_strata.tsv")),
    gf$strata, c("pair_id", "unit_id", "stratum_id"),
    "Pairwise group-collapsed strata")
  .fairqa_compare(.fairqa_read(file.path(result_dir,
    "prediction_join_audit.tsv")), gf$joins, c("pair_id", "unit_id"),
    "Prediction join audit")
  profile_sensitivity <- merge(gf$pairwise[, .(
    pair_id, method_a, method_b, unit_id, common_support,
    group_collapsed_estimate = estimate, group_collapsed_se = se,
    group_collapsed_ci_low = ci_low, group_collapsed_ci_high = ci_high)],
    pf$pairwise[, .(pair_id, unit_id, profile_common_support = common_support,
      profile_estimate = estimate, profile_se = se,
      profile_ci_low = ci_low, profile_ci_high = ci_high)],
    by = c("pair_id", "unit_id"), all = TRUE, sort = FALSE)
  profile_sensitivity[, difference_profile_minus_group :=
                        profile_estimate - group_collapsed_estimate]
  setcolorder(profile_sensitivity, c(
    "pair_id", "method_a", "method_b", "unit_id", "common_support",
    "group_collapsed_estimate", "group_collapsed_se",
    "group_collapsed_ci_low", "group_collapsed_ci_high",
    "profile_common_support", "profile_estimate", "profile_se",
    "profile_ci_low", "profile_ci_high", "difference_profile_minus_group"
  ))
  .fairqa_compare(.fairqa_read(file.path(result_dir,
    "pairwise_profile_sensitivity.tsv")), profile_sensitivity,
    c("pair_id", "unit_id"), "Pairwise profile sensitivity")
  primary <- gf$pairwise[primary_unit == TRUE]
  meta <- .fairqa_meta_all(primary, pairs, frozen)
  .fairqa_compare(.fairqa_read(file.path(result_dir, "pairwise_meta.tsv")), meta,
                  "pair_id", "REML/KH meta-analysis", tolerance = 1e-10)
  loco <- .fairqa_loco(primary, pairs)
  .fairqa_compare(.fairqa_read(file.path(result_dir, "pairwise_loco.tsv")), loco,
    c("pair_id", "omitted_unit"), "LOCO meta-analysis", tolerance = 1e-10)
  heterogeneity <- rbindlist(list(
    .fairqa_meta_strata(primary, pairs, units, "modality", TRUE),
    .fairqa_meta_strata(primary, pairs, units, "biospecimen", FALSE)), fill = TRUE)
  .fairqa_compare(.fairqa_read(file.path(result_dir,
    "pairwise_heterogeneity.tsv")), heterogeneity,
    c("stratum_variable", "stratum_value", "pair_id"),
    "Heterogeneity analysis", tolerance = 1e-10)
  panels <- .fairqa_panels(eligibility, methods, primary_units)
  .fairqa_compare(.fairqa_read(file.path(result_dir,
    "complete_support_panels.tsv")), panels,
    c("threshold", "item_type", "method_id", "unit_id"),
    "Complete-support panels")
  primary_performance <- performance[unit_id %in% primary_units]
  rectangles <- .fairqa_rectangles(primary_performance, panels)
  .fairqa_compare(.fairqa_read(file.path(result_dir,
    "complete_support_performance.tsv")), rectangles,
    c("threshold", "method_id", "unit_id"), "Complete-support rectangles")
  eligibility_out <- .fairqa_read(file.path(result_dir, "eligibility_audit.tsv"))
  .fairqa_compare(eligibility_out, eligibility,
                  c("method_id", "unit_id"), "Eligibility audit")
  coverage <- .fairqa_coverage(methods, eligibility, primary_performance,
                               predictions, splits, primary_units)
  .fairqa_compare(.fairqa_read(file.path(result_dir,
    "coverage_performance.tsv")), coverage, "method_id", "Coverage/Pareto")
  ranks <- .fairqa_rank_bootstrap(primary_performance, panels, bootstrap_reps, seed)
  .fairqa_compare(.fairqa_read(file.path(result_dir, "rank_bootstrap.tsv")), ranks,
    c("threshold", "method_id"), "Rank bootstrap")
  hindsight <- .fairqa_hindsight(predictions, primary_performance, panels,
                                 permutation_reps, seed)
  .fairqa_compare(.fairqa_read(file.path(result_dir, "hindsight_ceiling.tsv")),
    hindsight, "unit_id", "Hindsight ceiling")
  invisible(TRUE)
}

.fairqa_check_ids <- function() c(
  "runtime_start_end", "snapshot_and_sources", "result_exact_18",
  "exact_outer_folds", "group_collapsed_auc_and_pair_inference",
  "profile_sensitivity", "frozen_support_meta_by", "loco_complete",
  "heterogeneity_policy", "panels_rectangles", "coverage_pareto",
  "rank_10000_inclusive_ties", "hindsight_registered_reconstruction",
  "path_isolation", "provenance_identity"
)

.fairqa_assert_result_manifest_unchanged <- function(path, expected_sha256) {
  if (!file.exists(path) || .fairqa_sha(path) != expected_sha256) {
    stop("Result manifest changed after scientific reconstruction.")
  }
  invisible(TRUE)
}

.fairqa_validate_staged_receipt <- function(stage, pins) {
  expected <- c("validation_checks.tsv", "validation_manifest.tsv")
  present <- list.files(stage, recursive = FALSE, all.files = TRUE,
                        include.dirs = TRUE, no.. = TRUE)
  paths <- file.path(stage, present)
  if (!setequal(present, expected) || any(dir.exists(paths)) ||
      any(nzchar(Sys.readlink(paths)))) {
    stop("Staged validation receipt inventory is not exact.")
  }
  checks <- .fairqa_read(file.path(stage, "validation_checks.tsv"))
  if (!identical(names(checks), c("check_id", "pass")) ||
      !identical(as.character(checks$check_id), .fairqa_check_ids()) ||
      anyNA(checks) || !all(as.logical(checks$pass))) {
    stop("Staged validation checks are incomplete or failed.")
  }
  manifest <- .fairqa_read(file.path(stage, "validation_manifest.tsv"))
  manifest_columns <- c(
    "file", "bytes", "sha256", "run_mode", "omicselector_version",
    "package_git_commit", "installed_package_tree_sha256",
    "producer_script_sha256", "snapshot_manifest_sha256",
    "result_manifest_sha256", "runtime_receipt_manifest_sha256",
    "runtime_image_sha256"
  )
  pin_columns <- manifest_columns[-seq_len(3L)]
  if (!identical(names(manifest), manifest_columns) || nrow(manifest) != 1L ||
      !identical(names(pins), pin_columns) || anyNA(manifest) ||
      manifest$file[[1L]] != "validation_checks.tsv" ||
      as.numeric(manifest$bytes[[1L]]) !=
        file.info(file.path(stage, "validation_checks.tsv"))$size ||
      manifest$sha256[[1L]] !=
        .fairqa_sha(file.path(stage, "validation_checks.tsv"))) {
    stop("Staged validation receipt bytes differ from its manifest.")
  }
  for (column in pin_columns) {
    if (!identical(as.character(manifest[[column]][[1L]]),
                   as.character(pins[[column]]))) {
      stop("Staged validation receipt pin differs: ", column)
    }
  }
  invisible(TRUE)
}

.fairqa_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor is required.")
  opt <- .fairqa_args(args); code <- .fairqa_code_provenance()
  installed_root <- system.file(package = "OmicSelector")
  installed_version <- as.character(packageVersion("OmicSelector"))
  installed_script <- system.file("scripts", basename(code$script_path),
                                  package = "OmicSelector")
  if (code$commit != opt$expected_package_commit ||
      installed_version != opt$expected_package_version ||
      !nzchar(installed_root) || .fairqa_package_tree_sha256(installed_root) !=
        opt$expected_installed_package_tree_sha256 ||
      !nzchar(installed_script) || .fairqa_sha(installed_script) !=
        code$script_sha256 || (opt$run_mode == "full" && code$dirty)) {
    stop("Validator source/install do not match exact clean package pins.")
  }
  input_dir <- normalizePath(opt$input_dir, mustWork = TRUE)
  result_dir <- normalizePath(opt$result_dir, mustWork = TRUE)
  runtime_receipt <- normalizePath(opt$runtime_receipt, mustWork = TRUE)
  runtime_image <- normalizePath(opt$runtime_image, mustWork = TRUE)
  if (!identical(runtime_image, opt$runtime_image)) {
    stop("Runtime image argument must already be its canonical real path.")
  }
  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  output_link <- Sys.readlink(output_dir)
  if (file.exists(output_dir) || dir.exists(output_dir) ||
      (length(output_link) == 1L && !is.na(output_link) && nzchar(output_link))) {
    stop("Validation output already exists.")
  }
  source_inputs <- .fairqa_validate_source_inputs(input_dir)
  runtime_packages <- .fairqa_read(file.path(runtime_receipt,
                                              "package_closure.tsv"))
  .fairqa_assert_disjoint_output(output_dir, c(
    input_dir, result_dir, code$package_root, installed_root, runtime_receipt,
    runtime_image, R.home(), source_inputs, runtime_packages$root_path,
    runtime_packages$lib_path
  ))
  runtime_guard <- paper3_validate_runtime_receipt(
    runtime_receipt, opt$expected_runtime_receipt_manifest_sha256,
    opt$expected_package_version, opt$expected_package_commit,
    opt$expected_installed_package_tree_sha256, runtime_image,
    opt$expected_runtime_image_sha256, phase = "start"
  )
  snapshot <- .fairqa_validate_snapshot(
    input_dir, opt$expected_snapshot_manifest_sha256)
  for (column in c("omicselector_version", "package_git_commit",
                   "installed_package_tree_sha256",
                   "runtime_receipt_manifest_sha256", "runtime_image_sha256")) {
    expected <- switch(column,
      omicselector_version = opt$expected_package_version,
      package_git_commit = opt$expected_package_commit,
      installed_package_tree_sha256 = opt$expected_installed_package_tree_sha256,
      runtime_receipt_manifest_sha256 =
        opt$expected_runtime_receipt_manifest_sha256,
      runtime_image_sha256 = opt$expected_runtime_image_sha256)
    if (!identical(unique(snapshot$manifest[[column]]), expected)) {
      stop("Snapshot differs from validator pin: ", column)
    }
  }
  bundle <- .fairqa_validate_bundle_manifest(
    result_dir, opt$expected_snapshot_manifest_sha256,
    opt$expected_package_version, opt$expected_package_commit,
    opt$expected_installed_package_tree_sha256,
    opt$expected_runtime_receipt_manifest_sha256,
    opt$expected_runtime_image_sha256, opt$run_mode, code$package_root)
  validated_result_manifest_sha256 <- .fairqa_sha(bundle$path)
  .fairqa_validate_science(input_dir, result_dir, opt$run_mode,
    bundle$bootstrap_reps, bundle$permutation_reps, bundle$seed)
  for (file in c("provenance_preflight.log", "provenance_resolution.tsv")) {
    if (.fairqa_sha(file.path(input_dir, file)) !=
        .fairqa_sha(file.path(result_dir, file))) {
      stop("Result provenance bytes differ from snapshot.")
    }
  }
  stage <- tempfile(".fair-validation-stage-", tmpdir = output_parent)
  dir.create(stage)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  checks <- data.table(check_id = .fairqa_check_ids(), pass = TRUE)
  fwrite(checks, file.path(stage, "validation_checks.tsv"), sep = "\t")
  validation_manifest <- data.table(
    file = "validation_checks.tsv",
    bytes = file.info(file.path(stage, "validation_checks.tsv"))$size,
    sha256 = .fairqa_sha(file.path(stage, "validation_checks.tsv")),
    run_mode = opt$run_mode, omicselector_version = installed_version,
    package_git_commit = code$commit,
    installed_package_tree_sha256 = opt$expected_installed_package_tree_sha256,
    producer_script_sha256 = code$script_sha256,
    snapshot_manifest_sha256 = opt$expected_snapshot_manifest_sha256,
    result_manifest_sha256 = validated_result_manifest_sha256,
    runtime_receipt_manifest_sha256 =
      opt$expected_runtime_receipt_manifest_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256)
  fwrite(validation_manifest, file.path(stage, "validation_manifest.tsv"),
         sep = "\t")
  receipt_pins <- list(
    run_mode = opt$run_mode,
    omicselector_version = installed_version,
    package_git_commit = code$commit,
    installed_package_tree_sha256 = opt$expected_installed_package_tree_sha256,
    producer_script_sha256 = code$script_sha256,
    snapshot_manifest_sha256 = opt$expected_snapshot_manifest_sha256,
    result_manifest_sha256 = validated_result_manifest_sha256,
    runtime_receipt_manifest_sha256 =
      opt$expected_runtime_receipt_manifest_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256
  )
  .fairqa_validate_staged_receipt(stage, receipt_pins)
  .fairqa_validate_snapshot(input_dir, opt$expected_snapshot_manifest_sha256)
  .fairqa_validate_source_inputs(input_dir)
  final_bundle <- .fairqa_validate_bundle_manifest(result_dir,
    opt$expected_snapshot_manifest_sha256, opt$expected_package_version,
    opt$expected_package_commit, opt$expected_installed_package_tree_sha256,
    opt$expected_runtime_receipt_manifest_sha256,
    opt$expected_runtime_image_sha256, opt$run_mode, code$package_root)
  .fairqa_assert_result_manifest_unchanged(
    final_bundle$path, validated_result_manifest_sha256
  )
  if (!identical(.fairqa_code_provenance(), code) ||
      .fairqa_package_tree_sha256(installed_root) !=
        opt$expected_installed_package_tree_sha256) {
    stop("Source/install changed immediately before receipt promotion.")
  }
  paper3_validate_runtime_receipt(
    runtime_receipt, opt$expected_runtime_receipt_manifest_sha256,
    opt$expected_package_version, opt$expected_package_commit,
    opt$expected_installed_package_tree_sha256, runtime_image,
    opt$expected_runtime_image_sha256, phase = "end", start_guard = runtime_guard)
  .fairqa_validate_snapshot(input_dir, opt$expected_snapshot_manifest_sha256)
  .fairqa_validate_source_inputs(input_dir)
  if (!identical(.fairqa_code_provenance(), code) ||
      .fairqa_package_tree_sha256(installed_root) !=
        opt$expected_installed_package_tree_sha256) {
    stop("Source/install changed immediately before receipt promotion.")
  }
  final_bundle <- .fairqa_validate_bundle_manifest(result_dir,
    opt$expected_snapshot_manifest_sha256, opt$expected_package_version,
    opt$expected_package_commit, opt$expected_installed_package_tree_sha256,
    opt$expected_runtime_receipt_manifest_sha256,
    opt$expected_runtime_image_sha256, opt$run_mode, code$package_root)
  .fairqa_assert_result_manifest_unchanged(
    final_bundle$path, validated_result_manifest_sha256
  )
  .fairqa_validate_staged_receipt(stage, receipt_pins)
  if (!file.rename(stage, output_dir)) {
    stop("Could not atomically promote validation receipt.")
  }
  cat("Validated fair all-method comparison bundle:", result_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .fairqa_main()
