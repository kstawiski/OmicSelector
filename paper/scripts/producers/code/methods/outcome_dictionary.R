# outcome_dictionary.R
#
# Loader and applier for the v0.4 outcome dictionary that replaces the
# v0.3.1 regex-based outcome inference inside run_paper3_pipeline.R.
# Codex finding N4 (2026-05-04 HOLD) required a per-accession,
# per-(source_column, raw_label) explicit mapping.
#
# Dictionary location: knowledge/outcome_dictionary_v0_4.tsv
# Schema: accession | tier_role | source_metadata_path | source_column |
#         raw_label | y | include_exclude | include_exclude_reason | n | notes
#
# Key invariants enforced here (codex review 2026-05-05):
#   - Strip footer "# total_*" comment lines before parse.
#   - Trim CRLF on the last column.
#   - Validate documented row counts (258 unique-(accession,raw_label),
#     116 accessions).
#   - Compound key is (accession, source_column, raw_label) — multiple
#     source_columns per accession are permitted (HRA006520).
#   - Apply by sample-ID alignment, never row-position.
#   - Return a structured outcome object so the caller can audit
#     class counts, exclusion reasons, and unresolved samples.

suppressPackageStartupMessages({
  library(data.table)
})

DICT_DOCUMENTED_ROWS <- 258L
DICT_DOCUMENTED_ACCESSIONS <- 116L


# ---------------------------------------------------------------------------
# load_outcome_dictionary_v0_4
# ---------------------------------------------------------------------------

#' Load and validate the v0.4 outcome dictionary.
#'
#' @param path Path to the TSV file. Defaults to the canonical location.
#' @param strict If TRUE (default), error when the documented row counts
#'   don't match. If FALSE, warn and continue.
#' @return A data.table with the parsed rows (no comment/footer lines).
load_outcome_dictionary_v0_4 <- function(
  path = "knowledge/outcome_dictionary_v0_4.tsv",
  strict = TRUE) {

  if (!file.exists(path)) {
    stop("load_outcome_dictionary_v0_4: file not found at ", path)
  }

  raw_lines <- readLines(path, warn = FALSE)
  raw_lines <- gsub("\r$", "", raw_lines, perl = TRUE)
  is_comment <- grepl("^#", raw_lines)
  is_blank   <- !nzchar(raw_lines)
  body_lines <- raw_lines[!is_comment & !is_blank]

  parsed <- tryCatch(
    fread(text = body_lines, sep = "\t", quote = "", header = TRUE,
          na.strings = c("", "NA")),
    error = function(e) stop("load_outcome_dictionary_v0_4: parse error: ", e$message)
  )

  required_cols <- c("accession", "tier_role", "source_metadata_path",
                     "source_column", "raw_label", "y",
                     "include_exclude", "include_exclude_reason", "n", "notes")
  missing_cols <- setdiff(required_cols, colnames(parsed))
  if (length(missing_cols) > 0L) {
    stop("load_outcome_dictionary_v0_4: required columns missing: ",
         paste(missing_cols, collapse = ", "))
  }

  parsed[, y := suppressWarnings(as.integer(y))]
  if ("n" %in% colnames(parsed)) {
    parsed[, n := suppressWarnings(as.integer(n))]
  }

  bad_y <- parsed$y[!is.na(parsed$y) & !(parsed$y %in% c(0L, 1L))]
  if (length(bad_y) > 0L) {
    stop("load_outcome_dictionary_v0_4: y values outside {0,1,NA}: ",
         paste(unique(bad_y), collapse = ", "))
  }

  bad_inc <- parsed$include_exclude[!parsed$include_exclude %in%
                                       c("include", "exclude", NA)]
  if (length(bad_inc) > 0L) {
    stop("load_outcome_dictionary_v0_4: include_exclude values not in ",
         "{'include','exclude'}: ", paste(unique(bad_inc), collapse = ", "))
  }

  unique_pairs <- unique(parsed[, .(accession, raw_label)])
  unique_acc   <- unique(parsed$accession)

  msg_count <- sprintf(
    "loaded %d body rows, %d unique (accession,raw_label) pairs, %d unique accessions",
    nrow(parsed), nrow(unique_pairs), length(unique_acc))

  if (nrow(unique_pairs) != DICT_DOCUMENTED_ROWS ||
      length(unique_acc) != DICT_DOCUMENTED_ACCESSIONS) {
    drift_msg <- sprintf(
      "%s — drift from documented (258 unique pairs, 116 accessions)",
      msg_count)
    if (isTRUE(strict)) stop("load_outcome_dictionary_v0_4: ", drift_msg)
    warning("load_outcome_dictionary_v0_4: ", drift_msg)
  }

  attr(parsed, "checksum_sha1") <- digest_sha1_safe(body_lines)
  attr(parsed, "n_body_rows")   <- nrow(parsed)
  attr(parsed, "n_unique_pairs")<- nrow(unique_pairs)
  attr(parsed, "n_accessions")  <- length(unique_acc)
  attr(parsed, "source_path")   <- path
  attr(parsed, "load_time")     <- Sys.time()

  parsed[]
}


# Lightweight sha1 helper without taking a hard `digest` package dependency.
digest_sha1_safe <- function(lines) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(paste(lines, collapse = "\n"), algo = "sha1"))
  }
  # tools::md5sum-equivalent fallback without the extra package
  if (requireNamespace("tools", quietly = TRUE)) {
    tf <- tempfile()
    on.exit(unlink(tf), add = TRUE)
    writeLines(lines, tf)
    return(unname(tools::md5sum(tf)))
  }
  NA_character_
}


# ---------------------------------------------------------------------------
# select_dictionary_metadata_path
# ---------------------------------------------------------------------------

#' Choose the metadata file path for an accession from the dictionary.
#'
#' is `present_in == "both"`, but the dictionary `source_metadata_path`
#' field is authoritative for outcome assignment. This helper returns the
#' first existing dictionary-specified path; on miss, returns NA so the
select_dictionary_metadata_path <- function(dict_subset) {
  if (!is.data.frame(dict_subset) || nrow(dict_subset) == 0L) return(NA_character_)
  paths <- unique(dict_subset$source_metadata_path)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0L) return(NA_character_)
  for (p in paths) if (file.exists(p)) return(p)
  # Codex review 2026-05-05: when no dictionary-specified path exists, fail
  # closed (return NA) rather than handing back a non-existent path; the
  NA_character_
}


# ---------------------------------------------------------------------------
# apply_outcome_dictionary
# ---------------------------------------------------------------------------

#' Apply the outcome dictionary to a metadata data.frame, aligning by the
#' supplied sample IDs.
#'
#' @param meta            data.frame of cohort metadata.
#' @param expr_sample_ids character vector of sample IDs from the expression
#'                        matrix (typically `colnames(expr)` before transpose).
#' @param dict_subset     subset of the loaded dictionary for this accession.
#' @param sample_id_cols  candidate metadata column names containing the
#'                        sample IDs (priority order). Default is the union
#'                        of the GSM/GEO style and a few common aliases.
#' @return a list with `y_bin` (integer 0/1/NA aligned to expr_sample_ids),
#'         `keep_idx` (logical mask for samples retained), `meta_aligned`
#'         (metadata rows aligned to expr_sample_ids), `audit` (counts and
#'         exclusion reasons).
apply_outcome_dictionary <- function(meta,
                                     expr_sample_ids,
                                     dict_subset,
                                     sample_id_cols = c("gsm_id", "GSM", "sample_id",
                                                        "sampleID", "id", "Sample",
                                                        "sample", "geo_accession",
                                                        "title", "sample_name",
                                                        "Sample.Name"),
                                     allow_row_position_fallback = FALSE) {
  audit <- list(
    n_total          = length(expr_sample_ids),
    n_dict_rows      = nrow(dict_subset),
    source_columns   = unique(dict_subset$source_column),
    matched_by       = NA_character_,
    n_aligned_meta   = 0L,
    n_unaligned      = 0L,
    n_y_eq_0         = 0L,
    n_y_eq_1         = 0L,
    n_excluded       = 0L,
    n_unresolved_samples = 0L,
    n_conflicting_samples = 0L,
    conflict_examples = character(0),
    exclusion_reasons= character(0),
    unresolved_labels = character(0)
  )

  if (!is.data.frame(meta) || nrow(meta) == 0L) {
    return(list(y_bin = rep(NA_integer_, length(expr_sample_ids)),
                keep_idx = rep(FALSE, length(expr_sample_ids)),
                meta_aligned = NULL, audit = audit))
  }
  if (nrow(dict_subset) == 0L) {
    return(list(y_bin = rep(NA_integer_, length(expr_sample_ids)),
                keep_idx = rep(FALSE, length(expr_sample_ids)),
                meta_aligned = NULL, audit = audit))
  }

  # Codex review 2026-05-05 (BLOCKER 2): GSE94533 lost 60 samples to ID
  # mismatch (expression "e1036" vs metadata "E-1036"). Normalize before
  # matching: drop hyphens/dashes/underscores/dots, fold case. Only fall
  # back to a normalized match if the literal match is incomplete.
  .normalize_id <- function(x) toupper(gsub("[\\-_\\. ]", "", as.character(x), perl = TRUE))
  expr_ids_str  <- as.character(expr_sample_ids)

  bridge_dirs <- c(
    file.path("./external/Ext_data", dict_subset$accession[1]),
    file.path("data_bucket", dict_subset$accession[1]))
  for (bd in bridge_dirs) {
    bp <- file.path(bd, "sample_id_bridge.tsv")
    if (file.exists(bp)) {
      bridge <- data.table::fread(bp, sep = "\t")
      bridge_col <- intersect(c("title", "title_base", "chip_id", "expr_id"), names(bridge))
      for (bc in bridge_col) {
        bvals <- as.character(bridge[[bc]])
        n_hit <- sum(expr_ids_str %in% bvals)
        if (n_hit > length(expr_ids_str) * 0.5) {
          gsm_map <- setNames(as.character(bridge$gsm_id), bvals)
          expr_ids_str <- ifelse(expr_ids_str %in% names(gsm_map),
                                  gsm_map[expr_ids_str], expr_ids_str)
          break
        }
      }
      break
    }
  }

  expr_ids_norm <- .normalize_id(expr_ids_str)

  cn_meta <- colnames(meta)
  match_col <- NA_character_
  match_kind <- NA_character_  # "literal" | "normalized"
  best_n <- 0L
  for (cand in c(sample_id_cols, cn_meta)) {
    if (!cand %in% cn_meta) next
    col_str  <- as.character(meta[[cand]])
    col_norm <- .normalize_id(col_str)
    n_lit  <- sum(expr_ids_str  %in% col_str)
    n_norm <- sum(expr_ids_norm %in% col_norm)
    if (n_lit > best_n) {
      best_n <- n_lit; match_col <- cand; match_kind <- "literal"
    } else if (n_norm > best_n) {
      best_n <- n_norm; match_col <- cand; match_kind <- "normalized"
    }
    if (best_n == length(expr_sample_ids)) break
  }

  if (is.na(match_col)) {
    # v0.6 §5.10: row-position fallback is disabled by default (fail-closed).
    # Callers must explicitly pass allow_row_position_fallback = TRUE to opt in.
    if (isTRUE(allow_row_position_fallback) &&
        nrow(meta) == length(expr_sample_ids)) {
      audit$matched_by <- "row_position_fallback"
      meta_aligned <- meta
    } else {
      audit$matched_by <- "none"
      audit$n_unaligned <- length(expr_sample_ids)
      return(list(y_bin = rep(NA_integer_, length(expr_sample_ids)),
                  keep_idx = rep(FALSE, length(expr_sample_ids)),
                  meta_aligned = NULL, audit = audit))
    }
  } else {
    audit$matched_by <- sprintf("%s:%s", match_col, match_kind)
    if (match_kind == "literal") {
      idx <- match(expr_ids_str, as.character(meta[[match_col]]))
    } else {
      idx <- match(expr_ids_norm, .normalize_id(meta[[match_col]]))
    }
    audit$n_unaligned <- sum(is.na(idx))
    meta_aligned <- meta[idx, , drop = FALSE]
  }
  audit$n_aligned_meta <- sum(!is.na(meta_aligned[[1L]]))

  src_cols <- unique(dict_subset$source_column)
  src_cols <- src_cols[!is.na(src_cols) & nzchar(src_cols)]
  src_cols <- intersect(src_cols, cn_meta)

  y_bin <- rep(NA_integer_, length(expr_sample_ids))
  resolved_via <- rep(NA_character_, length(expr_sample_ids))
  per_col_y <- vector("list", length(src_cols))
  names(per_col_y) <- src_cols

  for (sc in src_cols) {
    raw_vals <- as.character(meta_aligned[[sc]])
    sub <- dict_subset[dict_subset$source_column == sc &
                         !is.na(dict_subset$raw_label), ]
    if (nrow(sub) == 0L) next
    lookup <- setNames(sub$y, sub$raw_label)
    inc    <- setNames(sub$include_exclude, sub$raw_label)
    reason <- setNames(sub$include_exclude_reason, sub$raw_label)
    col_y  <- rep(NA_integer_, length(raw_vals))

    for (i in seq_along(raw_vals)) {
      lab <- raw_vals[i]
      if (is.na(lab)) next
      if (!lab %in% names(lookup)) {
        audit$unresolved_labels <- unique(c(audit$unresolved_labels, lab))
        next
      }
      col_y[i] <- as.integer(lookup[[lab]])
      if (is.na(y_bin[i])) {
        y_bin[i] <- col_y[i]
        resolved_via[i] <- sc
        if (!is.na(inc[[lab]]) && inc[[lab]] == "exclude") {
          audit$n_excluded <- audit$n_excluded + 1L
          r <- reason[[lab]]
          if (!is.na(r) && nzchar(r)) {
            audit$exclusion_reasons <- c(audit$exclusion_reasons, r)
          }
        }
      }
    }
    per_col_y[[sc]] <- col_y
  }

  # Codex review 2026-05-05: detect per-sample multi-column conflicts so the
  if (length(src_cols) > 1L) {
    n_conflict <- 0L
    examples   <- character(0)
    for (i in seq_along(expr_sample_ids)) {
      vals <- vapply(per_col_y, function(v) v[i], integer(1L))
      vals_obs <- vals[!is.na(vals)]
      if (length(vals_obs) >= 2L && length(unique(vals_obs)) > 1L) {
        n_conflict <- n_conflict + 1L
        if (length(examples) < 5L) {
          examples <- c(examples,
            sprintf("sample=%s {%s}", expr_sample_ids[i],
              paste(sprintf("%s=%s", names(vals), as.character(vals)), collapse = ",")))
        }
        y_bin[i] <- NA_integer_
        resolved_via[i] <- "CONFLICT"
      }
    }
    audit$n_conflicting_samples <- n_conflict
    audit$conflict_examples     <- examples
  }

  audit$n_unresolved_samples <- sum(is.na(resolved_via))

  audit$n_y_eq_0 <- sum(y_bin == 0L, na.rm = TRUE)
  audit$n_y_eq_1 <- sum(y_bin == 1L, na.rm = TRUE)
  if (length(audit$exclusion_reasons) > 0L) {
    audit$exclusion_reasons <- table(audit$exclusion_reasons)
  } else {
    audit$exclusion_reasons <- table(character(0))
  }

  keep_idx <- !is.na(y_bin)
  list(y_bin = y_bin, keep_idx = keep_idx, meta_aligned = meta_aligned, audit = audit)
}


# ---------------------------------------------------------------------------
# flagged_cohorts_v0_4
# ---------------------------------------------------------------------------

#'
#' Source: outcome_dictionary_v0_4_methods.md Section 8 ("Accessions Requiring
#' Manual Review (FLAGGED for v0.4 pre-flight)"). Excluded from the primary
#' transfer path; may appear in supplementary QC tables with flags.
flagged_cohorts_v0_4 <- function() {
  c("GSE210329", "GSE314577", "GSE262260", "GSE201359",
    "GSE139031", "GSE307729", "GSE85589",  "GSE188232",
    "GSE209838", "GSE300882")
}


# ---------------------------------------------------------------------------
# Self-test (runs only when sourced as a script).
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  cwd <- getwd()
  if (basename(cwd) != "OmicSelector_paper") {
    candidate <- "."
    if (dir.exists(candidate)) setwd(candidate)
  }
  cat("[self-test] loading dictionary...\n")
  dict <- load_outcome_dictionary_v0_4()
  stopifnot(attr(dict, "n_unique_pairs") == DICT_DOCUMENTED_ROWS)
  stopifnot(attr(dict, "n_accessions")  == DICT_DOCUMENTED_ACCESSIONS)
  cat(sprintf("[self-test] OK %d rows / %d accessions / sha1=%s\n",
              attr(dict, "n_body_rows"), attr(dict, "n_accessions"),
              substr(attr(dict, "checksum_sha1"), 1, 12)))

  hra <- dict[accession == "HRA006520"]
  src_cols <- unique(hra$source_column)
  cat(sprintf("[self-test] HRA006520 source_columns: %s (n=%d)\n",
              paste(src_cols, collapse = ","), length(src_cols)))
  stopifnot(length(src_cols) >= 2L)

  fc <- flagged_cohorts_v0_4()
  cat(sprintf("[self-test] %d cohorts flagged for manual review\n", length(fc)))
}
