#!/usr/bin/env Rscript

# One immutable unit/seed task for the paper's fully nested single-sample
# selector analysis. Outer test outcomes are attached only after frozen scoring.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(OmicSelector)
})

.selector_task_args <- function(args) {
  out <- list(inner_folds = 5L, bootstrap_reps = 1000L,
              verify_base = TRUE)
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Arguments require --name=value: ", arg)
    key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg), fixed = TRUE)
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  required <- c("cache", "cache_sha256", "unit", "seed", "output_dir",
                "package_commit")
  missing <- required[!required %in% names(out)]
  if (length(missing)) stop("Missing: --", paste(missing, collapse = ", --"))
  out$seed <- as.integer(out$seed)
  out$inner_folds <- as.integer(out$inner_folds)
  out$bootstrap_reps <- as.integer(out$bootstrap_reps)
  verify_token <- tolower(as.character(out$verify_base))
  if (!verify_token %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("verify_base must be an explicit true/false token.")
  }
  out$verify_base <- verify_token %in% c("true", "1", "yes")
  if (anyNA(c(out$seed, out$inner_folds, out$bootstrap_reps)) ||
      !out$seed %in% c(101L, 202L, 303L, 404L, 505L) ||
      out$inner_folds != 5L || out$bootstrap_reps != 1000L ||
      !out$verify_base) {
    stop("Registered tasks require seeds 101/202/303/404/505, five requested ",
         "inner folds, 1,000 bootstraps, and verify_base=true.")
  }
  if (!grepl("^[0-9a-f]{64}$", tolower(out$cache_sha256)) ||
      !grepl("^[0-9a-f]{40}$", tolower(out$package_commit))) {
    stop("cache_sha256/package_commit must be exact 64/40-character hex pins.")
  }
  out
}

.selector_task_code_provenance <- function() {
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
    stop("Could not resolve the OmicSelector task code provenance.")
  }
  list(
    commit = unname(commit[[1L]]), dirty = length(status) > 0L,
    script_sha256 = digest::digest(script_path, algo = "sha256", file = TRUE)
  )
}

.selector_task_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .selector_task_args(args)
  code_provenance <- .selector_task_code_provenance()
  if (code_provenance$dirty ||
      !identical(tolower(opt$package_commit), code_provenance$commit)) {
    stop("Registered selector tasks require the exact clean package commit pin.")
  }
  cache_path <- normalizePath(opt$cache, mustWork = TRUE)
  actual_cache_hash <- digest::digest(cache_path, algo = "sha256", file = TRUE)
  if (!identical(actual_cache_hash, opt$cache_sha256)) {
    stop("Cohort cache does not match --cache-sha256.")
  }
  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) stop("Task output already exists: ", output_dir)
  temp_dir <- tempfile("selector-task-", tmpdir = output_parent)
  dir.create(temp_dir, recursive = TRUE)
  on.exit(if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE), add = TRUE)

  cache <- readRDS(cache_path)
  if (!opt$unit %in% names(cache)) stop("Unit absent from cache: ", opt$unit)
  cohort <- cache[[opt$unit]]
  required <- c("expr_per_sample", "y_bin", "group_id", "outer_k",
                "provenance_block", "source_accessions")
  missing <- setdiff(required, names(cohort))
  if (length(missing)) stop("Cached unit is missing: ", paste(missing, collapse = ", "))
  X <- as.matrix(cohort$expr_per_sample)
  y <- as.integer(cohort$y_bin)
  group_id <- as.character(cohort$group_id)
  sample_ids <- rownames(X)
  if (is.null(sample_ids)) sample_ids <- paste0(opt$unit, "::", seq_len(nrow(X)))
  if (!is.numeric(X) || nrow(X) != length(y) || length(y) != length(group_id) ||
      anyNA(X) || anyNA(y) || anyNA(group_id) || anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) || anyDuplicated(sample_ids) ||
      !setequal(unique(y), 0:1)) {
    stop("Cached unit has invalid matrix, labels, or group identifiers.")
  }
  group_labels <- split(y, group_id)
  if (any(vapply(group_labels, function(z) length(unique(z)) != 1L,
                 logical(1L)))) {
    stop("A biological/provenance group contains both outcome classes.")
  }
  outer_k <- as.integer(cohort$outer_k)
  if (!outer_k %in% c(3L, 5L)) stop("Cached outer_k must be 3 or 5.")
  outer_folds <- os_make_grouped_stratified_folds(
    y, group_id = group_id, n_folds = outer_k, seed = opt$seed
  )
  outer_validation <- os_validate_folds(
    outer_folds, y = y, group_id = group_id
  )
  if (!all(outer_validation$pass)) stop("Outer grouped folds failed validation.")

  methods <- singlesample_selector_candidates()
  if (length(methods) != 33L || !identical(methods[[1L]], "ws-rclr-trimmed")) {
    stop("The installed package does not expose the frozen 33-route selector bank.")
  }
  prediction_rows <- audit_rows <- route_rows <- states <- list()
  inner_candidate_fit_count <- 0L
  outer_candidate_refit_count <- 0L
  baseline_refit_count <- 0L
  for (outer_fold in seq_along(outer_folds)) {
    test <- outer_folds[[outer_fold]]
    train <- setdiff(seq_len(nrow(X)), test)
    inner_seed <- as.integer((as.double(opt$seed) + outer_fold * 1009) %%
                               (.Machine$integer.max - 1) + 1)
    selector <- fit_singlesample_selector(
      X_train = X[train, , drop = FALSE], y_train = y[train],
      methods = methods, route = "all", group_id = group_id[train],
      inner_folds = opt$inner_folds, seed = inner_seed,
      bootstrap_reps = opt$bootstrap_reps,
      verify_base = opt$verify_base, stack_stabilizer = 1e-6
    )
    if (!inherits(selector, "singlesample_selector_set") &&
        !inherits(selector, "singlesample_selector_ineligible")) {
      stop("Selector returned an unexpected object class.")
    }
    accounting <- selector$fit_accounting
    required_accounting <- c("inner_candidate_fits", "final_candidate_refits",
                             "total_selector_fits")
    if (is.null(accounting) ||
        !all(required_accounting %in% names(accounting)) ||
        !identical(as.integer(accounting$total_selector_fits),
                   as.integer(accounting$inner_candidate_fits +
                                accounting$final_candidate_refits))) {
      stop("Selector returned incomplete base-fit accounting.")
    }
    inner_candidate_fit_count <- inner_candidate_fit_count +
      as.integer(accounting$inner_candidate_fits)
    outer_candidate_refit_count <- outer_candidate_refit_count +
      as.integer(accounting$final_candidate_refits)

    # Frozen objects score the outer rows before held-out outcomes are attached.
    route_scores <- if (inherits(selector, "singlesample_selector_set")) {
      score_singlesample_selector(selector, X[test, , drop = FALSE])
    } else {
      as.data.frame(stats::setNames(
        replicate(3L, rep(NA_real_, length(test)), simplify = FALSE),
        c("best_auc", "lower_auc_ci", "simplex_stack")
      ))
    }
    rclr <- deploy_singlesample(
      X[train, , drop = FALSE], y[train], method = "ws-rclr-trimmed",
      verify = opt$verify_base
    )
    baseline_refit_count <- baseline_refit_count + 1L
    rclr_score <- score_specimen(rclr, X[test, , drop = FALSE])
    if (nrow(route_scores) != length(test) || length(rclr_score) != length(test)) {
      stop("Frozen outer scoring returned the wrong row count.")
    }
    selector_roundtrip <- unserialize(serialize(selector, NULL, version = 3L))
    if (inherits(selector, "singlesample_selector_set")) {
      replay <- score_singlesample_selector(
        selector_roundtrip, X[test, , drop = FALSE]
      )
      if (!isTRUE(all.equal(route_scores, replay, tolerance = 0,
                            check.attributes = FALSE))) {
        stop("Serialized selector does not exactly replay frozen outer scores.")
      }
    }
    selector <- selector_roundtrip
    heldout_y <- y[test]
    fold_predictions <- data.table(
      unit_id = opt$unit, seed = opt$seed, outer_fold = outer_fold,
      sample_id = sample_ids[test], group_id = group_id[test],
      y = heldout_y, score_rclr = as.numeric(rclr_score)
    )
    for (route in names(route_scores)) {
      route_object <- if (inherits(selector, "singlesample_selector_set")) {
        selector$selectors[[route]]
      } else selector
      row <- data.table::copy(fold_predictions)
      row[, `:=`(
        route = route,
        score = as.numeric(route_scores[[route]]),
        route_status = route_object$status,
        selected_methods = if (route_object$status == "eligible")
          paste(route_object$selected_methods, collapse = ";") else "",
        ineligible_reason = if (route_object$status == "eligible") NA_character_
          else route_object$reason
      )]
      prediction_rows[[length(prediction_rows) + 1L]] <- row
      route_rows[[length(route_rows) + 1L]] <- data.table(
        unit_id = opt$unit, seed = opt$seed, outer_fold = outer_fold,
        route = route, status = route_object$status,
        selected_methods = if (route_object$status == "eligible")
          paste(route_object$selected_methods, collapse = ";") else "",
        reason = if (route_object$status == "eligible") "eligible"
          else route_object$reason
      )
    }
    if (!is.null(selector$candidate_audit) && nrow(selector$candidate_audit)) {
      audit <- data.table::as.data.table(selector$candidate_audit)
      if (!setequal(audit$method_id, methods) || nrow(audit) != length(methods)) {
        stop("Selector candidate audit does not cover the frozen bank exactly.")
      }
    } else {
      audit <- data.table(
        method_id = methods, complete = FALSE,
        reason = selector$reason %||% "selector_ineligible_before_candidate_fit"
      )
    }
    audit[, `:=`(unit_id = opt$unit, seed = opt$seed,
                 outer_fold = outer_fold)]
    audit_rows[[length(audit_rows) + 1L]] <- audit
    states[[paste0("outer_fold_", outer_fold)]] <- selector
  }
  predictions <- data.table::rbindlist(prediction_rows, fill = TRUE)
  route_audit <- data.table::rbindlist(route_rows, fill = TRUE)
  candidate_audit <- data.table::rbindlist(audit_rows, fill = TRUE)
  expected_prediction_rows <- nrow(X) * 3L
  if (nrow(predictions) != expected_prediction_rows ||
      anyDuplicated(predictions[, paste(route, sample_id, sep = "\r")])) {
    stop("Outer predictions do not contain exactly one row per route/profile.")
  }
  data.table::fwrite(predictions, file.path(temp_dir, "selector_predictions.tsv.gz"),
                     sep = "\t", compress = "gzip")
  data.table::fwrite(route_audit, file.path(temp_dir, "selector_route_audit.tsv"),
                     sep = "\t")
  data.table::fwrite(candidate_audit,
                     file.path(temp_dir, "selector_candidate_audit.tsv"), sep = "\t")
  data.table::fwrite(outer_validation,
                     file.path(temp_dir, "outer_fold_validation.tsv"), sep = "\t")
  saveRDS(states, file.path(temp_dir, "selector_states.rds"), version = 3L)
  code_provenance_final <- .selector_task_code_provenance()
  if (!identical(code_provenance_final, code_provenance)) {
    stop("OmicSelector task code provenance changed during execution.")
  }
  manifest_paths <- list.files(temp_dir, full.names = TRUE)
  base_fit_count <- inner_candidate_fit_count + outer_candidate_refit_count +
    baseline_refit_count
  manifest <- data.table(
    file = basename(manifest_paths), bytes = file.info(manifest_paths)$size,
    sha256 = vapply(manifest_paths, digest::digest, character(1L),
                    algo = "sha256", file = TRUE),
    unit_id = opt$unit, seed = opt$seed, outer_k = outer_k,
    candidate_bank_n = length(methods),
    inner_candidate_fit_count = inner_candidate_fit_count,
    outer_candidate_refit_count = outer_candidate_refit_count,
    baseline_refit_count = baseline_refit_count,
    base_fit_count = base_fit_count,
    fit_budget = 175000L, cache_sha256 = actual_cache_hash,
    package_version = as.character(utils::packageVersion("OmicSelector")),
    package_commit = opt$package_commit,
    script_sha256 = code_provenance$script_sha256,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  data.table::fwrite(manifest, file.path(temp_dir, "task_manifest.tsv"), sep = "\t")
  if (!file.rename(temp_dir, output_dir)) stop("Could not atomically promote task output.")
  cat("Created nested-selector task:", output_dir, "\n")
  invisible(output_dir)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

if (sys.nframe() == 0L) .selector_task_main()
