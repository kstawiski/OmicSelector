#!/usr/bin/env Rscript

# Exact analysis runner for the OmicSelector paper's fair all-method comparison.
# It consumes only an immutable snapshot produced by
# paper3_build_fair_snapshot.R and writes one atomic output bundle.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(OmicSelector)
})

.fair_args <- function(args) {
  out <- list(bootstrap_reps = 10000L, permutation_reps = 1000L,
              seed = 20260718L, run_mode = "full")
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Arguments require --name=value: ", arg)
    key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg), fixed = TRUE)
    out[[key]] <- sub("^--[^=]+=", "", arg)
  }
  if (is.null(out$input_dir) || is.null(out$output_dir)) {
    stop("Required: --input-dir=PATH --output-dir=PATH")
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
    script_sha256 = digest::digest(script_path, algo = "sha256", file = TRUE),
    package_git_commit = unname(commit[[1L]]),
    package_git_dirty = length(status) > 0L
  )
}

.fair_read <- function(path) {
  if (!file.exists(path)) stop("Missing snapshot file: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    data.table::fread(cmd = paste("gzip -cd", shQuote(path)))
  } else data.table::fread(path)
}

.fair_validate_snapshot <- function(input_dir, require_clean = FALSE) {
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  required <- c("methods.tsv", "units.tsv", "cells.tsv", "eligibility.tsv",
                "pair_family.tsv", "splits.tsv", "predictions.tsv.gz",
                "provenance_preflight.log",
                "provenance_resolution.tsv", "source_inputs.tsv",
                "manifest.tsv")
  missing <- required[!file.exists(file.path(input_dir, required))]
  if (length(missing)) stop("Snapshot is incomplete: ", paste(missing, collapse = ", "))
  manifest <- data.table::fread(file.path(input_dir, "manifest.tsv"))
  required_manifest <- c(
    "file", "bytes", "sha256", "omicselector_version",
    "package_git_commit", "package_git_dirty", "producer_script_sha256",
    "snapshot_run_mode"
  )
  if (!all(required_manifest %in% names(manifest)) ||
      anyDuplicated(manifest$file) || !setequal(manifest$file, setdiff(required, "manifest.tsv"))) {
    stop("Snapshot manifest does not enumerate the exact input artifacts.")
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
  list(input_dir = input_dir, manifest = manifest)
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

.fair_precompute_strata <- function(predictions, splits) {
  label_check <- predictions[, data.table::uniqueN(y),
                             by = .(unit_id, seed, fold, group_id)]
  if (any(label_check$V1 != 1L)) {
    stop("A biological group carries discordant labels within a stratum.")
  }
  collapsed <- predictions[, .(y = unique(y), score = mean(score)),
                           by = .(method_id, unit_id, seed, fold, group_id)]
  strata <- collapsed[, .(
    auc = .fair_auc(y, score), n_test_groups = .N
  ), by = .(method_id, unit_id, seed, fold)]
  split_counts <- merge(
    splits[, .(unit_id, seed, fold, split_id,
               registered_test_groups = n_test_groups)],
    splits[, .(total_groups = sum(n_test_groups)), by = .(unit_id, seed)],
    by = c("unit_id", "seed"), all.x = TRUE, sort = FALSE
  )
  strata <- merge(strata, split_counts,
                  by = c("unit_id", "seed", "fold"),
                  all.x = TRUE, sort = FALSE)
  if (anyNA(strata[, .(split_id, registered_test_groups, total_groups)]) ||
      any(strata$n_test_groups != strata$registered_test_groups) ||
      any(!is.finite(strata$auc))) {
    stop("Precomputed group-primary strata disagree with canonical splits.")
  }
  strata[, `:=`(
    n_train_groups = total_groups - n_test_groups,
    stratum_id = paste(seed, fold, sep = "::")
  )]
  if (any(strata$n_train_groups <= 0L)) {
    stop("A precomputed stratum has no remaining training groups.")
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
      stop("Eligible method/unit has an undefined group-primary AUC.")
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
  ea <- eligibility[list(pair$method_a, unit), nomatch = 0L]
  eb <- eligibility[list(pair$method_b, unit), nomatch = 0L]
  base <- data.table(
    pair_id = pair$pair_id, method_a = pair$method_a, method_b = pair$method_b,
    unit_id = unit, common_support = FALSE, estimate = NA_real_, se = NA_real_,
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
                 "n_test_groups", "n_train_groups")
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
    n_test = joined$n_test_groups,
    n_train = joined$n_train_groups,
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
    stratum_id, auc_a, auc_b, effect, n_test_groups, n_train_groups,
    pair_id = pair$pair_id, method_a = pair$method_a,
    method_b = pair$method_b, unit_id = unit
  )]
  list(summary = base, strata = strata, join = join_audit)
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
    p_below_minus_margin = stats::pt((estimate + 0.05) / se, df = df,
                                     lower.tail = TRUE),
    df = df, status = "REML_KH", fit_reason = "eligible"
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

.fair_meta_strata <- function(pairwise, pairs, units, field) {
  annotated <- merge(pairwise, units[, c("unit_id", field), with = FALSE],
                     by = "unit_id", all.x = TRUE, sort = FALSE)
  rows <- list()
  values <- sort(unique(annotated[[field]]))
  for (value in values) {
    for (i in seq_len(nrow(pairs))) {
      pair <- pairs[i]
      d <- annotated[pair_id == pair$pair_id & get(field) == value]
      support <- sort(d[common_support == TRUE, unit_id])
      row <- .fair_meta_one(d, pair, support)
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
    d <- pairwise[pair_id == pair$pair_id & common_support &
                    is.finite(se) & se > 0]
    if (nrow(d) < 6L) next
    for (unit in d$unit_id) {
      keep <- d[unit_id != unit]
      fit <- tryCatch(
        withCallingHandlers(
          metafor::rma.uni(yi = keep$estimate, sei = keep$se,
                           method = "REML", test = "knha"),
          warning = function(w) stop(conditionMessage(w), call. = FALSE)
        ), error = identity
      )
      if (inherits(fit, "error")) next
      rows[[length(rows) + 1L]] <- data.table(
        pair_id = pair$pair_id, omitted_unit = unit, k = nrow(keep),
        estimate = as.numeric(fit$b[[1L]]),
        ci_low = as.numeric(fit$ci.lb[[1L]]),
        ci_high = as.numeric(fit$ci.ub[[1L]])
      )
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
    wide <- data.table::dcast(rectangle, unit_id ~ method_id, value.var = "auc")
    if (anyNA(wide)) stop("Complete-support performance rectangle contains NA.")
    set.seed(seed + as.integer(round(threshold_value * 1000)))
    rank_draws <- matrix(NA_real_, nrow = reps, ncol = length(methods),
                         dimnames = list(NULL, methods))
    for (b in seq_len(reps)) {
      sampled <- sample(seq_len(nrow(wide)), nrow(wide), replace = TRUE)
      means <- colMeans(as.matrix(wide[sampled, ..methods]))
      rank_draws[b, ] <- rank(-means, ties.method = "average")
    }
    for (method in methods) {
      values <- rank_draws[, method]
      rows[[length(rows) + 1L]] <- data.table(
        threshold = threshold_value, method_id = method,
        rank_median = stats::median(values),
        rank_ci_low = unname(stats::quantile(values, 0.025, type = 6L)),
        rank_ci_high = unname(stats::quantile(values, 0.975, type = 6L)),
        p_top1 = mean(values <= 1), p_top3 = mean(values <= 3),
        p_top5 = mean(values <= 5), bootstrap_reps = reps, seed = seed
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
    collapsed <- d[, .(y = unique(y), score = mean(score)),
                   by = .(method_id, stratum_id, group_id)]
    labels <- unique(collapsed[, .(stratum_id, group_id, y)])
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
      # preserves the fixed stratified fold's class counts and the joint score
      # correlation across methods, avoiding invalid one-class null folds.
      perm <- labels[, .(group_id, y_perm = sample(y)), by = stratum_id]
      x <- merge(collapsed[, .(method_id, stratum_id, group_id, score)],
                 perm, by = c("stratum_id", "group_id"), sort = FALSE)
      aucs <- x[, .(auc = .fair_auc(y_perm, score)),
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
        valid_fraction = valid_fraction, status = "permutation_guard_failed"
      )
    } else {
      rows[[u]] <- data.table(
        unit_id = unit, observed_max = observed_max,
        null_median = stats::median(null),
        null_ci_low = unname(stats::quantile(null, 0.025, type = 6L)),
        null_ci_high = unname(stats::quantile(null, 0.975, type = 6L)),
        headroom = observed_max - stats::median(null),
        valid_permutations = length(null), attempts = attempts,
        valid_fraction = valid_fraction, status = "evaluable"
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
    opt$input_dir, require_clean = opt$run_mode == "full"
  )
  input_dir <- snapshot$input_dir
  code_provenance <- .fair_code_provenance()
  snapshot_commit <- unique(snapshot$manifest$package_git_commit)
  snapshot_version <- unique(snapshot$manifest$omicselector_version)
  if (length(snapshot_commit) != 1L ||
      !identical(snapshot_commit, code_provenance$package_git_commit) ||
      length(snapshot_version) != 1L ||
      !identical(snapshot_version,
                 as.character(utils::packageVersion("OmicSelector")))) {
    stop("Snapshot and runner do not share one OmicSelector version/commit pin.")
  }
  if (opt$run_mode == "full" && code_provenance$package_git_dirty) {
    stop("Full analysis requires a clean OmicSelector package worktree.")
  }
  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) stop("Output bundle already exists: ", output_dir)
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
  stratum_auc <- .fair_precompute_strata(predictions, splits)
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
  pair_rows <- stratum_rows <- join_rows <- list()
  index <- 0L
  for (p in seq_len(nrow(pairs))) {
    pair <- pairs[p]
    for (unit in all_units) {
      index <- index + 1L
      result <- .fair_pair_unit(
        pair, unit, stratum_auc, eligibility, splits, prediction_counts
      )
      pair_rows[[index]] <- result$summary
      join_rows[[index]] <- result$join
      if (!is.null(result$strata)) stratum_rows[[length(stratum_rows) + 1L]] <- result$strata
    }
  }
  pairwise <- data.table::rbindlist(pair_rows, fill = TRUE)
  pair_strata <- data.table::rbindlist(stratum_rows, fill = TRUE)
  join_audit <- data.table::rbindlist(join_rows, fill = TRUE)
  pairwise <- merge(pairwise, units[, .(unit_id, primary_unit)],
                    by = "unit_id", all.x = TRUE, sort = FALSE)
  if (nrow(pairwise) != 1225L * 34L || any(!join_audit$pass &
                                           join_audit$reason != "not_common_eligible")) {
    stop("Pairwise accounting or exact joins failed.")
  }
  pairwise_primary <- pairwise[primary_unit == TRUE]
  meta <- .fair_meta_all(pairwise_primary, pairs, frozen_family)
  loco <- .fair_loco(pairwise_primary, pairs)
  heterogeneity <- data.table::rbindlist(list(
    .fair_meta_strata(pairwise_primary, pairs, units, "modality"),
    .fair_meta_strata(pairwise_primary, pairs, units, "biospecimen")
  ), fill = TRUE)

  performance_all <- .fair_method_unit_performance(
    stratum_auc, eligibility, splits
  )
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
                .(eligible_units = sum(eligible), coverage = mean(eligible),
                  eligibility_units = paste(unit_id[eligible], collapse = ";")),
                by = method_id],
    performance[eligible == TRUE, .(median_auc = stats::median(auc),
                            performance_units = paste(unit_id, collapse = ";")),
                by = method_id],
    by = "method_id", all.x = TRUE
  )
  coverage <- merge(methods[, .(method_id, roster_order, family, role, tier,
                                dep_route)], coverage, by = "method_id", sort = FALSE)
  hindsight <- .fair_hindsight(
    predictions, performance, panels, units, opt$permutation_reps, opt$seed
  )
  code_provenance_final <- .fair_code_provenance()
  if (!identical(code_provenance_final, code_provenance)) {
    stop("OmicSelector code provenance changed while running the comparison.")
  }

  data.table::fwrite(pairwise, file.path(temp_dir, "pairwise_cohort.tsv"), sep = "\t")
  data.table::fwrite(pair_strata, file.path(temp_dir, "pairwise_strata.tsv"), sep = "\t")
  data.table::fwrite(meta, file.path(temp_dir, "pairwise_meta.tsv"), sep = "\t")
  data.table::fwrite(loco, file.path(temp_dir, "pairwise_loco.tsv"), sep = "\t")
  data.table::fwrite(heterogeneity, file.path(temp_dir, "pairwise_heterogeneity.tsv"), sep = "\t")
  data.table::fwrite(panel_table, file.path(temp_dir, "complete_support_panels.tsv"), sep = "\t")
  data.table::fwrite(rectangles, file.path(temp_dir, "complete_support_performance.tsv"), sep = "\t")
  data.table::fwrite(rank_bootstrap, file.path(temp_dir, "rank_bootstrap.tsv"), sep = "\t")
  data.table::fwrite(coverage, file.path(temp_dir, "coverage_performance.tsv"), sep = "\t")
  data.table::fwrite(performance_all, file.path(temp_dir, "method_unit_performance.tsv"), sep = "\t")
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
    paste0("Pooled-testable frozen pair family: ", sum(meta$pooled_testable)),
    paste0("REML/Knapp--Hartung fits: ", sum(meta$status == "REML_KH")),
    paste0("Non-degenerate complete-support panels: ",
           sum(!panels$summary$degenerate)), "",
    "Results are descriptive until independent results review, plausibility ",
    "assessment, and canonical promotion. No universal winner is inferred here."
  )
  writeLines(report, file.path(temp_dir, "report.md"))
  artifacts <- list.files(temp_dir, full.names = TRUE)
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
    producer_script_sha256 = code_provenance$script_sha256,
    snapshot_manifest_sha256 = digest::digest(
      file.path(input_dir, "manifest.tsv"), algo = "sha256", file = TRUE
    )
  )
  data.table::fwrite(output_manifest, file.path(temp_dir, "output_manifest.tsv"), sep = "\t")
  if (!file.rename(temp_dir, output_dir)) stop("Could not atomically promote output bundle.")
  cat("Created fair all-method comparison bundle:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .fair_main()
