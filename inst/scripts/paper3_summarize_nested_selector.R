#!/usr/bin/env Rscript

# Aggregate the 34 x 5 immutable nested-selector tasks, apply corrected
# repeated-CV and cross-cohort inference, and preserve all route failures.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(OmicSelector)
})

.selector_summary_args <- function(args) {
  out <- list(run_mode = "full")
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Arguments require --name=value: ", arg)
    key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg), fixed = TRUE)
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  required <- c("task_root", "units", "output_dir")
  missing <- required[!required %in% names(out)]
  if (length(missing)) stop("Missing: --", paste(missing, collapse = ", --"))
  if (!out$run_mode %in% c("full", "smoke")) stop("run_mode must be full or smoke.")
  out
}

.selector_summary_code_provenance <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
  } else normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
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
    stop("Could not resolve selector-synthesis code provenance.")
  }
  list(
    commit = unname(commit[[1L]]), dirty = length(status) > 0L,
    script_sha256 = digest::digest(script_path, algo = "sha256", file = TRUE)
  )
}

.selector_verify_task <- function(path) {
  required <- c("selector_predictions.tsv.gz", "selector_route_audit.tsv",
                "selector_candidate_audit.tsv", "outer_fold_validation.tsv",
                "selector_states.rds", "task_manifest.tsv")
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) stop("Incomplete selector task: ", path)
  manifest <- data.table::fread(file.path(path, "task_manifest.tsv"))
  if (nrow(manifest) != length(required) - 1L || anyDuplicated(manifest$file) ||
      !setequal(manifest$file, setdiff(required, "task_manifest.tsv"))) {
    stop("Task manifest has the wrong artifact set: ", path)
  }
  actual <- vapply(file.path(path, manifest$file), digest::digest, character(1L),
                   algo = "sha256", file = TRUE)
  if (!identical(unname(actual), as.character(manifest$sha256)) ||
      any(file.info(file.path(path, manifest$file))$size != manifest$bytes)) {
    stop("Task bytes do not match manifest: ", path)
  }
  manifest
}

.selector_read_gz <- function(path) {
  data.table::fread(cmd = paste("gzip -cd", shQuote(path)))
}

.selector_meta <- function(d, route) {
  row <- data.table(
    route = route, k_support = nrow(d), k_meta = 0L,
    support_units = paste(sort(d$unit_id), collapse = ";"),
    estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    tau2 = NA_real_, i2 = NA_real_, p_zero = NA_real_,
    p_above_margin = NA_real_, df = NA_real_, status = "meta_unavailable",
    fit_reason = "not_fitted"
  )
  use <- d[is.finite(estimate) & is.finite(se) & se > 0]
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
  estimate <- as.numeric(fit$b[[1L]])
  se <- as.numeric(fit$se[[1L]])
  df <- nrow(use) - 1L
  row[, `:=`(
    estimate = estimate, se = se,
    ci_low = as.numeric(fit$ci.lb[[1L]]), ci_high = as.numeric(fit$ci.ub[[1L]]),
    tau2 = as.numeric(fit$tau2), i2 = as.numeric(fit$I2),
    p_zero = as.numeric(fit$pval[[1L]]),
    p_above_margin = stats::pt((estimate - 0.05) / se, df = df,
                               lower.tail = FALSE),
    df = df, status = "REML_KH", fit_reason = "eligible"
  )]
  row
}

.selector_summary_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Selector synthesis requires the suggested package 'metafor'.")
  }
  opt <- .selector_summary_args(args)
  code_provenance <- .selector_summary_code_provenance()
  if (opt$run_mode == "full" && code_provenance$dirty) {
    stop("Full selector synthesis requires a clean OmicSelector worktree.")
  }
  task_root <- normalizePath(opt$task_root, mustWork = TRUE)
  units <- data.table::fread(opt$units)
  if (!all(c("unit_id", "primary_unit", "outer_k") %in% names(units)) ||
      nrow(units) != 34L || sum(units$primary_unit) != 33L) {
    stop("Units table must contain the frozen 34/33 unit accounting.")
  }
  task_dirs <- list.dirs(task_root, recursive = FALSE, full.names = TRUE)
  if (!length(task_dirs)) stop("No immutable selector task directories were found.")
  manifests <- lapply(task_dirs, .selector_verify_task)
  manifest <- data.table::rbindlist(manifests, fill = TRUE)
  accounting_columns <- c(
    "inner_candidate_fit_count", "outer_candidate_refit_count",
    "baseline_refit_count", "base_fit_count", "fit_budget"
  )
  if (!all(accounting_columns %in% names(manifest)) ||
      any(manifest$fit_budget != 175000L) ||
      any(manifest$base_fit_count !=
            manifest$inner_candidate_fit_count +
              manifest$outer_candidate_refit_count +
              manifest$baseline_refit_count)) {
    stop("Selector tasks have incomplete or inconsistent fit accounting.")
  }
  task_keys <- unique(manifest[, .(unit_id, seed)])
  expected_seeds <- sort(unique(task_keys$seed))
  if (opt$run_mode == "full" &&
      (nrow(task_keys) != 170L || length(expected_seeds) != 5L ||
       !identical(expected_seeds, c(101L, 202L, 303L, 404L, 505L)) ||
       !setequal(task_keys$unit_id, units$unit_id))) {
    stop("Full selector synthesis requires all 34 x 5 unit/seed tasks.")
  }
  if (sum(unique(manifest[, .(unit_id, seed, base_fit_count)])$base_fit_count) >
      175000L) {
    stop("Selector base-fit budget exceeded 175,000.")
  }
  task_outer <- unique(manifest[, .(unit_id, seed, outer_k)])
  task_outer <- merge(task_outer, units[, .(unit_id, registered_outer_k = outer_k)],
                      by = "unit_id", all.x = TRUE)
  if (anyNA(task_outer$registered_outer_k) ||
      any(task_outer$outer_k != task_outer$registered_outer_k)) {
    stop("Task outer_k values do not match the frozen units table.")
  }
  theoretical_max <- sum(units$outer_k) * 5L * (33L * 5L + 33L + 1L)
  if (theoretical_max > 169150L) {
    stop("Frozen unit fold counts exceed the registered 169,150-fit maximum.")
  }
  pin_columns <- c("cache_sha256", "package_version", "package_commit",
                   "script_sha256")
  if (any(vapply(pin_columns, function(column) uniqueN(manifest[[column]]) != 1L,
                 logical(1L)))) {
    stop("Selector tasks do not share one cache/package pin.")
  }
  if (!identical(unique(manifest$package_commit), code_provenance$commit) ||
      !identical(unique(manifest$package_version),
                 as.character(utils::packageVersion("OmicSelector")))) {
    stop("Selector tasks do not match the executing package version/commit.")
  }

  predictions <- audits <- candidates <- states <- list()
  for (i in seq_along(task_dirs)) {
    path <- task_dirs[[i]]
    predictions[[i]] <- .selector_read_gz(
      file.path(path, "selector_predictions.tsv.gz")
    )
    audits[[i]] <- data.table::fread(file.path(path, "selector_route_audit.tsv"))
    candidates[[i]] <- data.table::fread(
      file.path(path, "selector_candidate_audit.tsv")
    )
    states[[i]] <- readRDS(file.path(path, "selector_states.rds"))
  }
  predictions <- data.table::rbindlist(predictions, fill = TRUE)
  route_audit <- data.table::rbindlist(audits, fill = TRUE)
  candidate_audit <- data.table::rbindlist(candidates, fill = TRUE)
  expected_outer_tasks <- sum(unique(manifest[, .(unit_id, seed, outer_k)])$outer_k)
  if (nrow(route_audit) != expected_outer_tasks * 3L ||
      nrow(candidate_audit) != expected_outer_tasks * 33L ||
      anyDuplicated(route_audit[, paste(unit_id, seed, outer_fold, route,
                                        sep = "\r")]) ||
      anyDuplicated(candidate_audit[, paste(unit_id, seed, outer_fold,
                                            method_id, sep = "\r")])) {
    stop("Aggregated selector route/candidate audits have invalid cardinality.")
  }
  key <- predictions[, paste(route, unit_id, seed, outer_fold, sample_id,
                             sep = "\r")]
  if (anyDuplicated(key)) stop("Aggregated selector predictions are duplicated.")
  expected_routes <- c("best_auc", "lower_auc_ci", "simplex_stack")
  if (!setequal(predictions$route, expected_routes)) {
    stop("Aggregated predictions do not contain the three frozen routes.")
  }

  cohort_rows <- list()
  for (unit in units$unit_id) {
    unit_outer_k <- units[unit_id == unit, outer_k][[1L]]
    expected_strata <- as.vector(outer(expected_seeds, seq_len(unit_outer_k),
                                       paste, sep = "::"))
    for (route_value in expected_routes) {
      d <- predictions[unit_id == unit & route == route_value]
      base <- data.table(
        unit_id = unit, route = route_value, complete = FALSE,
        estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        p_zero = NA_real_, p_above_margin = NA_real_, expected_m = length(expected_strata),
        observed_m = 0L, reason = "incomplete_or_ineligible_route"
      )
      if (nrow(d) && all(d$route_status == "eligible") &&
          all(is.finite(d$score)) && all(is.finite(d$score_rclr))) {
        stratum_id <- paste(d$seed, d$outer_fold, sep = "::")
        result <- singlesample_matched_pair_auc(
          d$y, d$score, d$score_rclr, stratum_id, d$sample_id, d$group_id,
          expected_strata = expected_strata
        )
        s <- result$summary
        base[, `:=`(
          complete = TRUE, estimate = s$estimate, se = s$se,
          ci_low = s$ci_low, ci_high = s$ci_high,
          p_zero = s$p_zero_two_sided, p_above_margin = s$p_above_margin,
          observed_m = s$observed_m,
          reason = if (s$inference_available) "evaluable" else "zero_corrected_se"
        )]
      }
      cohort_rows[[length(cohort_rows) + 1L]] <- base
    }
  }
  cohort <- data.table::rbindlist(cohort_rows)
  cohort <- merge(cohort, units[, .(unit_id, primary_unit, modality, biospecimen)],
                  by = "unit_id", all.x = TRUE, sort = FALSE)
  meta <- data.table::rbindlist(lapply(expected_routes, function(route_name) {
    .selector_meta(cohort[primary_unit == TRUE & route == route_name & complete == TRUE],
                   route_name)
  }))
  zero <- singlesample_adjust_zero_by(meta$route, meta$p_zero, expected_routes)
  meta[, `:=`(p_zero_by = zero$p_zero_by,
              p_above_margin_by = stats::p.adjust(
                p_above_margin, method = "BY", n = length(expected_routes)
              ),
              zero_family_n = zero$family_n,
              relevance_family_n = length(expected_routes))]

  selection_routes <- c("best_auc", "lower_auc_ci")
  route_keys <- route_audit[
    status == "eligible" & route %in% selection_routes,
    .(unit_id, seed, outer_fold, route)
  ]
  availability <- merge(
    candidate_audit[complete == TRUE,
                    .(unit_id, seed, outer_fold, method_id)],
    route_keys, by = c("unit_id", "seed", "outer_fold"),
    allow.cartesian = TRUE
  )[, .(available_folds = .N), by = .(route, method_id)]
  selection_grid <- data.table::CJ(
    route = selection_routes,
    method_id = singlesample_selector_candidates(), unique = TRUE
  )
  selection_frequency <- merge(selection_grid, availability,
                               by = c("route", "method_id"), all.x = TRUE)
  selection_frequency <- merge(selection_frequency,
                               route_audit[
                                 status == "eligible" & route %in% selection_routes,
                                 .(evaluable_route_folds = .N), by = route
                               ], by = "route", all.x = TRUE)
  selection_frequency <- merge(selection_frequency,
                               route_audit[
                                 status == "eligible" & route %in% selection_routes &
                                   nzchar(selected_methods),
                                 .(selected_folds = .N),
                                 by = .(route, method_id = selected_methods)
                               ], by = c("route", "method_id"), all.x = TRUE)
  selection_frequency[is.na(available_folds), available_folds := 0L]
  selection_frequency[is.na(selected_folds), selected_folds := 0L]
  if (any(selection_frequency$selected_folds >
            selection_frequency$available_folds)) {
    stop("A selector method was selected in more folds than it was available.")
  }
  selection_frequency[, selection_fraction := data.table::fifelse(
    available_folds > 0L, selected_folds / available_folds, NA_real_
  )]
  candidate_availability <- candidate_audit[, .(
    available_outer_folds = sum(complete), total_outer_folds = .N,
    availability_fraction = mean(complete)
  ), by = method_id]

  stack_rows <- list()
  for (state_index in seq_along(states)) {
    task <- states[[state_index]]
    task_pin <- unique(manifests[[state_index]][, .(unit_id, seed)])
    if (nrow(task_pin) != 1L) stop("A task manifest has multiple unit/seed keys.")
    for (outer_name in names(task)) {
      selector_set <- task[[outer_name]]
      if (!inherits(selector_set, "singlesample_selector_set")) next
      stack <- selector_set$selectors$simplex_stack
      if (inherits(stack, "singlesample_selector_ineligible")) next
      weights <- stack$stack$weights
      stack_rows[[length(stack_rows) + 1L]] <- data.table(
        unit_id = task_pin$unit_id, seed = task_pin$seed,
        outer_fold = as.integer(sub("^outer_fold_", "", outer_name)),
        method_id = names(weights), weight = as.numeric(weights),
        nonzero = as.numeric(weights) > 1e-12,
        effective_size = 1 / sum(as.numeric(weights) ^ 2)
      )
    }
  }
  stack_weights <- if (length(stack_rows)) data.table::rbindlist(stack_rows) else data.table()

  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) stop("Output already exists: ", output_dir)
  temp_dir <- tempfile("selector-summary-", tmpdir = output_parent)
  dir.create(temp_dir, recursive = TRUE)
  on.exit(if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE), add = TRUE)
  data.table::fwrite(cohort, file.path(temp_dir, "selector_cohort.tsv"), sep = "\t")
  data.table::fwrite(meta, file.path(temp_dir, "selector_meta.tsv"), sep = "\t")
  data.table::fwrite(route_audit, file.path(temp_dir, "selector_route_audit.tsv"), sep = "\t")
  data.table::fwrite(candidate_availability,
                     file.path(temp_dir, "selector_candidate_availability.tsv"), sep = "\t")
  data.table::fwrite(selection_frequency,
                     file.path(temp_dir, "selector_selection_frequency.tsv"), sep = "\t")
  data.table::fwrite(stack_weights, file.path(temp_dir, "selector_stack_weights.tsv"), sep = "\t")
  data.table::fwrite(predictions, file.path(temp_dir, "selector_predictions.tsv.gz"),
                     sep = "\t", compress = "gzip")
  code_provenance_final <- .selector_summary_code_provenance()
  if (!identical(code_provenance_final, code_provenance)) {
    stop("OmicSelector synthesis code provenance changed during execution.")
  }
  report <- c(
    "# Fully nested single-sample selector", "",
    paste0("Run mode: ", opt$run_mode),
    paste0("Unit/seed tasks: ", nrow(task_keys)),
    paste0("Base fits: ", sum(unique(manifest[, .(unit_id, seed, base_fit_count)])$base_fit_count),
           " / 175000"),
    paste0("Registered fold-specific maximum: ", theoretical_max,
           " / 169150"),
    paste0("Complete route/unit cells: ", sum(cohort$complete), " / ", nrow(cohort)),
    paste0("REML/Knapp--Hartung routes: ", sum(meta$status == "REML_KH")), "",
    "No route is recommended until independent results review, plausibility ",
    "assessment, and canonical promotion."
  )
  writeLines(report, file.path(temp_dir, "report.md"))
  paths <- list.files(temp_dir, full.names = TRUE)
  out_manifest <- data.table(
    file = basename(paths), bytes = file.info(paths)$size,
    sha256 = vapply(paths, digest::digest, character(1L),
                    algo = "sha256", file = TRUE),
    run_mode = opt$run_mode,
    package_version = unique(manifest$package_version),
    package_commit = unique(manifest$package_commit),
    cache_sha256 = unique(manifest$cache_sha256),
    units_sha256 = digest::digest(opt$units, algo = "sha256", file = TRUE),
    producer_script_sha256 = code_provenance$script_sha256,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  data.table::fwrite(out_manifest, file.path(temp_dir, "output_manifest.tsv"), sep = "\t")
  if (!file.rename(temp_dir, output_dir)) stop("Could not atomically promote summary.")
  cat("Created nested-selector summary:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .selector_summary_main()
