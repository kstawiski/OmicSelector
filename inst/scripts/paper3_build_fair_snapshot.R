#!/usr/bin/env Rscript

# Build the immutable input snapshot for the OmicSelector paper's fair
# all-method comparison. This script never modifies its source benchmark files.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(OmicSelector)
})

.snapshot_read <- function(path) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    data.table::fread(cmd = paste("gzip -cd", shQuote(path)))
  } else {
    data.table::fread(path)
  }
}

.snapshot_args <- function(args) {
  values <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Arguments must use --name=value syntax: ", arg)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (key %in% names(values)) stop("Duplicate argument: --", key)
    values[[key]] <- value
  }
  required <- c("base-cells", "base-predictions", "cohorts", "output-dir",
                "provenance-script", "provenance-root",
                "provenance-resolution", "upstream-validation", "run-mode")
  missing <- setdiff(required, names(values))
  if (length(missing)) {
    stop("Missing required argument(s): --", paste(missing, collapse = ", --"))
  }
  if (!values[["run-mode"]] %in% c("full", "smoke")) {
    stop("--run-mode must be full or smoke.")
  }
  values
}

.snapshot_paths <- function(value) {
  if (is.null(value) || !nzchar(value)) character() else strsplit(value, ",", fixed = TRUE)[[1L]]
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

.snapshot_validate_upstream_receipt <- function(path) {
  receipt <- .snapshot_read(path)
  required_columns <- c("check_id", "pass")
  if (!all(required_columns %in% names(receipt)) ||
      anyDuplicated(receipt$check_id)) {
    stop("Upstream validation receipt has an invalid schema.")
  }
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
  rows <- receipt[match(required_checks, receipt$check_id)]
  passed <- tolower(as.character(rows$pass)) %in% c("true", "1")
  if (nrow(rows) != length(required_checks) || anyNA(rows$check_id) ||
      !all(passed)) {
    stop("Upstream validation receipt does not pass every required benchmark check.")
  }
  invisible(receipt)
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
  output <- suppressWarnings(system2(
    "Rscript", c(normalizePath(script), sources),
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
  overlap_n <- if (length(overlap_line)) {
    as.integer(sub(".*Overlap hits:[[:space:]]*([0-9]+).*", "\\1", overlap_line[[1L]]))
  } else NA_integer_
  if (status == 2L) {
    if (!is.finite(overlap_n) || overlap_n < 1L || nrow(resolution) != overlap_n) {
      stop("Overlap resolution row count does not match the fresh preflight.")
    }
    if (anyNA(resolution[, ..required]) || any(resolution$unit_a != resolution$unit_b)) {
      stop("A known overlap crosses canonical analysis units.")
    }
  } else if (nrow(resolution) != 0L ||
             (is.finite(overlap_n) && overlap_n != 0L)) {
    stop("A no-overlap preflight requires an empty provenance resolution.")
  }
  list(status = status, overlap_n = ifelse(is.na(overlap_n), 0L, overlap_n),
       resolution = resolution)
}

.snapshot_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .snapshot_args(args)
  output_dir <- normalizePath(dirname(opt[["output-dir"]]), mustWork = TRUE)
  output_dir <- file.path(output_dir, basename(opt[["output-dir"]]))
  if (file.exists(output_dir)) stop("Output snapshot already exists: ", output_dir)
  temp_dir <- tempfile("fair-snapshot-", tmpdir = dirname(output_dir))
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE), add = TRUE)
  code_provenance_start <- .snapshot_code_provenance()
  if (opt[["run-mode"]] == "full" &&
      isTRUE(code_provenance_start$package_git_dirty)) {
    stop("Full snapshots require a clean OmicSelector package worktree.")
  }

  overlay_cell_paths <- .snapshot_paths(opt[["overlay-cells"]])
  overlay_prediction_paths <- .snapshot_paths(opt[["overlay-predictions"]])
  provenance_root <- normalizePath(opt[["provenance-root"]], mustWork = TRUE)
  provenance_manifest <- file.path(
    provenance_root, "data_bucket", "_PROVENANCE_MANIFEST.md"
  )
  provenance_union <- file.path(provenance_root, "knowledge", "union_inventory.tsv")
  if (!file.exists(provenance_manifest)) {
    stop("The canonical provenance manifest is absent from --provenance-root.")
  }
  provenance_dependencies <- c(
    provenance_manifest,
    if (file.exists(provenance_union)) provenance_union else character()
  )
  source_inputs <- data.table::data.table(
    role = c(
      "base_cells", "base_predictions", "cohort_manifest",
      "provenance_script", "provenance_resolution", "upstream_validation",
      "provenance_manifest",
      if (file.exists(provenance_union)) "provenance_union_inventory" else character(),
      rep("overlay_cells", length(overlay_cell_paths)),
      rep("overlay_predictions", length(overlay_prediction_paths))
    ),
    path = c(
      opt[["base-cells"]], opt[["base-predictions"]], opt[["cohorts"]],
      opt[["provenance-script"]], opt[["provenance-resolution"]],
      opt[["upstream-validation"]], provenance_dependencies, overlay_cell_paths,
      overlay_prediction_paths
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
  .snapshot_validate_upstream_receipt(opt[["upstream-validation"]])

  roster <- data.table::as.data.table(singlesample_method_roster())
  roster[, roster_order := .I]
  roster <- roster[estimand == "within"]
  if (nrow(roster) != 50L || anyDuplicated(roster$method_id)) {
    stop("Installed OmicSelector does not expose the frozen 50-route within roster.")
  }
  units <- .snapshot_read(opt$cohorts)[status == "included"]
  if (!"accession" %in% names(units)) stop("Cohort table lacks accession.")
  data.table::setnames(units, "accession", "unit_id")
  if (nrow(units) != 34L || anyDuplicated(units$unit_id)) {
    stop("Cohort table must contain exactly 34 unique included units.")
  }
  required_unit <- c("unit_id", "disease", "modality", "biospecimen",
                     "provenance_block", "source_accessions", "outer_k")
  if (!all(required_unit %in% names(units))) stop("Cohort metadata is incomplete.")
  units[, primary_unit := !grepl("whole[ -]?blood", biospecimen,
                                 ignore.case = TRUE)]
  if (sum(units$primary_unit) != 33L) {
    stop("Exactly 33 circulating-biofluid primary units are required.")
  }

  base_cells <- .snapshot_normalize_cells(.snapshot_read(opt[["base-cells"]]))
  base_predictions <- .snapshot_normalize_predictions(
    .snapshot_read(opt[["base-predictions"]])
  )
  overlaid <- .snapshot_overlay(
    base_cells, base_predictions,
    overlay_cell_paths,
    overlay_prediction_paths
  )
  cells <- overlaid$cells
  predictions <- overlaid$predictions
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
  required_r2 <- registered[method %in% roster[tier == "R2", method_id]]
  required_r2_key <- required_r2[, paste(method, cohort, seed, sep = "\r")]
  if (opt[["run-mode"]] == "full" &&
      (!length(overlaid$overlay_cell_keys) ||
       !setequal(overlaid$overlay_cell_keys, required_r2_key) ||
       length(overlaid$overlay_cell_keys) != length(required_r2_key))) {
    stop("Full snapshots require the exact 13-method tier-R2 cell grid.")
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
      # An analysis unit can be structurally non-evaluable for every registered
      # method (for example, all routes are constant in every fold). Preserve
      # that negative coverage result in cells/eligibility; there is no split
      # to invent and no prediction row to synthesize.
      if (!length(available)) next
      reference <- predictions[
        cohort == unit & seed == seed_value & method == available[[1L]],
        .(fold, sample_id, group_id, y)
      ]
      data.table::setorder(reference, fold, sample_id, group_id, y)
      .snapshot_validate_reference(reference)
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
  predictions <- merge(predictions, split_map,
                       by = c("cohort", "seed"), all.x = TRUE, sort = FALSE)
  if (anyNA(predictions$split_id)) stop("A prediction lacks a canonical split id.")
  splits <- data.table::rbindlist(split_rows)
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
    opt[["provenance-script"]], opt[["provenance-root"]], units,
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
  predictions_out <- predictions[, .(
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
    producer_script_sha256 = code_provenance$script_sha256,
    snapshot_run_mode = opt[["run-mode"]],
    provenance_exit = provenance$status,
    provenance_overlap_pairs = provenance$overlap_n
  )
  data.table::fwrite(manifest, file.path(temp_dir, "manifest.tsv"), sep = "\t")
  if (!file.rename(temp_dir, output_dir)) {
    stop("Could not atomically promote snapshot to: ", output_dir)
  }
  cat("Created immutable fair-comparison snapshot:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) .snapshot_main()
