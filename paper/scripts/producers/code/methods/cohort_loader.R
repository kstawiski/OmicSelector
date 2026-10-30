# cohort_loader.R
#
# Shared expression-loading, miRNA-normalisation, and scoring helpers for
# supplementary validation scripts (tables S49-S55).
# Extracted verbatim from code/figures/loco_external_validation.R so
#
# Source this file AFTER:
#   suppressMessages(devtools::load_all("./software/OmicSelector", quiet=TRUE))
#   source("code/methods/mirna_alias_resolver.R")
#   source("code/methods/outcome_dictionary.R")
#   inv  <- fread("knowledge/union_inventory.tsv", sep="\t", quote="")
#   dict <- load_outcome_dictionary_v0_4()

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# ---------------------------------------------------------------------------
# 1. Expression-file candidate resolution
# ---------------------------------------------------------------------------

.candidates_for <- function(acc, present_in) {
  ext_dir <- file.path("./external/Ext_data", acc)
  db_dir  <- file.path("data_bucket", acc)
  ext_cands <- character(0)
  if (present_in %in% c("Ext_data", "both")) {
    ext_cands <- c(
      file.path(ext_dir, "expression_matrix_smrnaseq.tsv"),
      file.path(ext_dir, "expression_matrix_rawcounts.tsv"),
      file.path(ext_dir, "expression_matrix_rawCq.tsv"),
      file.path(ext_dir, "expression_matrix.tsv"),
      file.path(ext_dir, sprintf("%s_expression.tsv.gz", acc)),
      file.path(ext_dir, sprintf("%s_raw_counts.txt.gz", acc)),
      file.path(ext_dir, sprintf("%s_rawCounts.tsv.gz", acc)),
      file.path(ext_dir, sprintf("%s_miRNA_counts.txt.gz", acc)),
      file.path(ext_dir, sprintf("%s_non_normalized.txt.gz", acc)),
      file.path(ext_dir, sprintf("%s_non-normalized.txt.gz", acc)),
      file.path(ext_dir, sprintf("%s_non-normalized_data.txt.gz", acc)),
      file.path(ext_dir, sprintf("%s_processed_data_matrix.txt.gz", acc)),
      Sys.glob(file.path(ext_dir, sprintf("%s_*NR-gene_row.txt.gz", acc))),
      file.path(ext_dir, sprintf("%s_RAW.tar", acc)),
      file.path(ext_dir, sprintf("%s_series_matrix.txt.gz", acc)),
      file.path(ext_dir, sprintf("%s_family.soft.gz", acc)),
      file.path(ext_dir, sprintf("%s_fold-change.txt.gz", acc))
    )
  }
  db_cands <- c(
    file.path(db_dir, "expression_matrix_smrnaseq.tsv"),
    file.path(db_dir, "expression_matrix_rawcounts.tsv"),
    file.path(db_dir, "expression_matrix_rawCq.tsv"),
    file.path(db_dir, "expression_matrix.tsv"),
    file.path(db_dir, sprintf("%s_expression.tsv.gz", acc)),
    file.path(db_dir, sprintf("%s_raw_counts.txt.gz", acc)),
    file.path(db_dir, sprintf("%s_rawCounts.tsv.gz", acc)),
    file.path(db_dir, sprintf("%s_miRNA_counts.txt.gz", acc)),
    file.path(db_dir, sprintf("%s_non_normalized.txt.gz", acc)),
    file.path(db_dir, sprintf("%s_non-normalized.txt.gz", acc)),
    file.path(db_dir, "suppl", sprintf("%s_miRNA_counts.txt.gz", acc)),
    file.path(db_dir, "suppl", sprintf("%s_non_normalized.txt.gz", acc)),
    file.path(db_dir, "suppl", sprintf("%s_non-normalized.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_miRNA_counts.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_raw_counts.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_rawCounts.tsv.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_non_normalized.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_non-normalized.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_non-normalized_data.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_processed_data_matrix.txt.gz", acc)),
    file.path(db_dir, "raw", "suppl", sprintf("%s_miRNA_counts.txt.gz", acc)),
    file.path(db_dir, "raw", "suppl", sprintf("%s_non_normalized.txt.gz", acc)),
    file.path(db_dir, "raw", "suppl", sprintf("%s_non-normalized.txt.gz", acc)),
    file.path(db_dir, "matrix", sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(db_dir, "raw", "matrix", sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_RAW.tar", acc)),
    file.path(db_dir, "raw", sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_family.soft.gz", acc)),
    file.path(db_dir, "soft", sprintf("%s_family.soft.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_fold-change.txt.gz", acc))
  )
  cands <- unique(c(ext_cands, db_cands))
  cands <- cands[!is.na(cands) & nzchar(cands)]
  cands[file.exists(cands)]
}

.metadata_path_for <- function(acc, present_in, dict_subset) {
  declared <- select_dictionary_metadata_path(dict_subset)
  if (!is.na(declared)) return(declared)
  cands <- c()
  if (present_in %in% c("Ext_data", "both"))
    cands <- c(cands,
               file.path("./external/Ext_data", acc, "sample_metadata.tsv"))
  cands <- c(cands, file.path("data_bucket", acc, "metadata.tsv"))
  cands <- cands[file.exists(cands)]
  if (length(cands) == 0L) return(NA_character_)
  cands[1L]
}

# ---------------------------------------------------------------------------
# 2. Raw-file readers
# ---------------------------------------------------------------------------

.read_lines_any <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

.clean_ids <- function(x) {
  x <- trimws(as.character(x))
  gsub('^"|"$', "", x)
}

.normalize_feature_ids <- function(ids) {
  ids <- .clean_ids(ids)
  vapply(ids, function(id) {
    parts <- unlist(strsplit(id, "[_|;,[:space:]]+", perl = TRUE), use.names = FALSE)
    parts <- parts[nzchar(parts)]
    mir_parts <- parts[grepl("^(hsa[-_])?(miR|mir|let)[-_]", parts, ignore.case = TRUE)]
    if (length(mir_parts) > 0L) {
      arm_parts <- mir_parts[grepl("[-_][35]p$", mir_parts, ignore.case = TRUE)]
      id <- tail(if (length(arm_parts) > 0L) arm_parts else mir_parts, 1L)
    }
    id <- gsub("_", "-", id, fixed = TRUE)
    id <- sub("^((hsa-)?(miR|mir|let)-.+)-[0-9]{5,}$", "\\1", id, ignore.case = TRUE)
    id <- sub("^hsa-mir-", "hsa-miR-", id, ignore.case = TRUE)
    id <- sub("^hsa-let-", "hsa-let-", id, ignore.case = TRUE)
    id
  }, character(1L), USE.NAMES = FALSE)
}

.mirbase_v22_cache <- new.env(parent = emptyenv())

.load_mirbase_v22_mature_map <- function() {
  if (exists("map", envir = .mirbase_v22_cache, inherits = FALSE))
    return(get("map", envir = .mirbase_v22_cache, inherits = FALSE))
  fa <- "./data/reference/mirbase_v22/mature.fa"
  gz <- paste0('',fa,'.gz')
  if (!file.exists(fa) && !file.exists(gz)) {
    dir.create(dirname(fa), showWarnings = FALSE, recursive = TRUE)
    utils::download.file("https://www.mirbase.org/ftp/22/mature.fa.gz",
                         gz, mode = "wb", quiet = TRUE)
  }
  src <- if (file.exists(fa)) fa else gz
  con <- if (grepl("\\.gz$", src, ignore.case = TRUE)) gzfile(src, "rt") else file(src, "rt")
  on.exit(close(con))
  headers <- grep("^>hsa-", readLines(con, warn = FALSE), value = TRUE)
  if (length(headers) == 0L) stop("No human mature miRNA records in ", src)
  fields  <- strsplit(sub("^>", "", headers), "[[:space:]]+", perl = TRUE)
  mature  <- vapply(fields, `[`, character(1L), 1L)
  mimat   <- toupper(vapply(fields, `[`, character(1L), 2L))
  name_keys <- toupper(gsub("_", "-", mature, fixed = TRUE))
  map <- list(by_name  = stats::setNames(mature, name_keys),
              by_mimat = stats::setNames(mature, mimat),
              mature = mature, mature_key = name_keys, source = src)
  assign("map", map, envir = .mirbase_v22_cache)
  map
}

.mirna_candidate_token <- function(id) {
  id <- .clean_ids(id)
  if (!nzchar(id)) return(NA_character_)
  mimat <- regmatches(id, regexpr("MIMAT[0-9]{7}", id, ignore.case = TRUE, perl = TRUE))
  if (length(mimat) == 1L && nzchar(mimat)) return(toupper(mimat))
  precursor <- regmatches(id, regexpr("\\bMI[0-9]{7}\\b", id, ignore.case = TRUE, perl = TRUE))
  if (length(precursor) == 1L && nzchar(precursor)) return(toupper(precursor))
  parts <- unique(c(id, unlist(strsplit(id, "[|;,[:space:]]+", perl = TRUE), use.names = FALSE)))
  for (p in parts) {
    p <- .clean_ids(p)
    p <- gsub("_", "-", p, fixed = TRUE)
    p <- sub("\\.[0-9]+$", "", p)
    p <- sub("^hsa-mir-", "hsa-miR-", p, ignore.case = TRUE)
    p <- sub("^hsa-mir", "hsa-miR", p, ignore.case = TRUE)
    p <- sub("^hsa-let-", "hsa-let-", p, ignore.case = TRUE)
    p <- sub("^mir-", "hsa-miR-", p, ignore.case = TRUE)
    p <- sub("^mir",  "hsa-miR",  p, ignore.case = TRUE)
    p <- sub("^let-", "hsa-let-", p, ignore.case = TRUE)
    p <- sub("-([35])P$", "-\\1p", p, ignore.case = TRUE)
    p <- sub("(-[35]p)[.-][0-9]+$", "\\1", p, ignore.case = TRUE)
    if (grepl("^hsa-(miR|let)-", p)) return(p)
  }
  NA_character_
}

.canonicalize_mirna_ids <- function(ids) {
  map    <- .load_mirbase_v22_mature_map()
  tokens <- vapply(ids, .mirna_candidate_token, character(1L), USE.NAMES = FALSE)
  out    <- rep(NA_character_, length(tokens))
  mimat  <- grepl("^MIMAT[0-9]{7}$", tokens, ignore.case = TRUE)
  out[mimat] <- unname(map$by_mimat[toupper(tokens[mimat])])
  precursor <- grepl("^MI[0-9]{7}$", tokens, ignore.case = TRUE) & !mimat
  out[precursor] <- NA_character_
  name_like <- !mimat & !precursor & !is.na(tokens)
  if (any(name_like)) {
    keys   <- toupper(gsub("_", "-", tokens[name_like], fixed = TRUE))
    direct <- unname(map$by_name[keys])
    unresolved <- is.na(direct)
    if (any(unresolved)) {
      base_keys <- sub("-[0-9]+$", "", keys[unresolved])
      direct2   <- unname(map$by_name[base_keys])
      direct[which(unresolved)[!is.na(direct2)]] <- direct2[!is.na(direct2)]
    }
    out[name_like] <- direct
  }
  out
}

normalize_to_canonical_mirna_namespace <- function(expr_per_sample,
                                                    scale_state = c("log", "abundance", "auto")) {
  scale_state <- match.arg(scale_state)
  if (!is.matrix(expr_per_sample)) expr_per_sample <- as.matrix(expr_per_sample)
  if (is.null(colnames(expr_per_sample)))
    stop("normalize_to_canonical_mirna_namespace: matrix must have column names")
  if (scale_state == "auto") {
    scale_state <- if (isTRUE(any(expr_per_sample < 0, na.rm = TRUE))) "log" else "abundance"
  }
  collapse_method <- if (scale_state == "log") "median_log" else "sum_abundance"
  feats     <- colnames(expr_per_sample)
  canonical <- .canonicalize_mirna_ids(feats)
  keep      <- !is.na(canonical) & nzchar(canonical)
  if (!any(keep)) {
    out <- expr_per_sample[, FALSE, drop = FALSE]
    attr(out, "canonical_namespace_source") <- .load_mirbase_v22_mature_map()$source
    attr(out, "n_canonical_output_features") <- 0L
    return(out)
  }
  groups   <- split(which(keep), canonical[keep])
  out_cols <- vector("list", length(groups))
  collapsed_groups <- 0L
  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    if (length(idx) == 1L) {
      out_cols[[g]] <- expr_per_sample[, idx, drop = TRUE]
    } else {
      sub <- expr_per_sample[, idx, drop = FALSE]
      out_cols[[g]] <- if (scale_state == "log") {
        apply(sub, 1L, function(v) median(v, na.rm = TRUE))
      } else {
        rowSums(sub, na.rm = TRUE)
      }
      collapsed_groups <- collapsed_groups + 1L
    }
  }
  out <- do.call(cbind, out_cols)
  colnames(out) <- names(groups)
  rownames(out) <- rownames(expr_per_sample)
  storage.mode(out) <- "double"
  attr(out, "canonical_namespace_source")  <- .load_mirbase_v22_mature_map()$source
  attr(out, "canonical_collapse_method")   <- collapse_method
  attr(out, "n_canonical_output_features") <- ncol(out)
  attr(out, "n_canonical_collapsed_groups") <- collapsed_groups
  out
}

.mirna_fraction <- function(x) {
  x <- .clean_ids(x)
  mean(grepl("^(hsa[-_]|hsa[-_]?miR|miR[-_]|let[-_]|MIMAT|MI[0-9]{4})",
             x, ignore.case = TRUE), na.rm = TRUE)
}

.plain_skip <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  i <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    i <- i + 1L
    s <- trimws(line)
    if (!nzchar(s)) next
    is_hash_header <- grepl("^#?(miRNA|ID_REF|Gene ID|Gene\\t|Sample\\t|ID\\t)", s, ignore.case = TRUE)
    if (startsWith(s, "!")) next
    if (startsWith(s, "#") && !is_hash_header) next
    return(max(0L, i - 1L))
  }
  0L
}

.read_plain_dt <- function(path, acc = NA_character_) {
  sep_char <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
  skip_n <- .plain_skip(path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  widths <- utils::count.fields(con, sep = sep_char, quote = "", skip = skip_n,
                               blank.lines.skip = TRUE, comment.char = "")
  close(con); on.exit(NULL, add = FALSE)
  if (length(widths) < 2L) stop("Expression matrix has no data rows: ", path)

  header_width <- widths[1L]
  body_widths <- widths[-1L]
  if (all(body_widths == header_width)) {
    return(fread(path, sep = sep_char, quote = "", fill = FALSE, skip = skip_n,
                 header = TRUE, check.names = FALSE,
                 na.strings = c("", "NA", "NaN", "null", "N/A", "Undetermined",
                                "undetermined", "UNDETERMINED")))
  }

  # GSE198834's deposited count table omits only the top-left feature-ID
  # header.  Recognise its complete 155-column sample signature; do not infer
  # derivations for arbitrary ragged tables.
  header_con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(header_con), add = TRUE)
  header_line <- readLines(header_con, n = skip_n + 1L, warn = FALSE)[skip_n + 1L]
  close(header_con); on.exit(NULL, add = FALSE)
  header_fields <- strsplit(header_line, sep_char, fixed = TRUE)[[1L]]
  expected_gse198834 <- sprintf("chip.sample%03d.tsv", seq_len(155L))
  recognised_missing_corner <- identical(as.character(acc), "GSE198834") &&
    identical(header_fields, expected_gse198834) &&
    header_width == 155L && all(body_widths == 156L)
  if (!recognised_missing_corner) {
    stop(sprintf(
      "Unrecognized ragged expression matrix: header has %d fields; data widths are %s (%s)",
      header_width, paste(sort(unique(body_widths)), collapse = ","), path))
  }
  fread(path, sep = sep_char, quote = "", fill = FALSE, skip = skip_n + 1L,
        header = FALSE, col.names = c("feature_id", header_fields),
        check.names = FALSE,
        na.strings = c("", "NA", "NaN", "null", "N/A", "Undetermined",
                       "undetermined", "UNDETERMINED"))
}

.series_paths_for <- function(acc, expr_path) {
  ext_dir <- file.path("./external/Ext_data", acc)
  db_dir  <- file.path("data_bucket", acc)
  unique(c(
    file.path(dirname(expr_path), sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(dirname(dirname(expr_path)), sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(ext_dir, sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(db_dir, "matrix", sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(db_dir, "raw", "matrix", sprintf("%s_series_matrix.txt.gz", acc)),
    file.path(db_dir, "raw", sprintf("%s_series_matrix.txt.gz", acc))
  ))
}

.parse_series_metadata_map <- function(path) {
  if (!file.exists(path)) return(character(0))
  lines    <- .read_lines_any(path)
  geo_line <- grep("^!Sample_geo_accession\\t", lines, value = TRUE)
  if (length(geo_line) == 0L) return(character(0))
  geo <- .clean_ids(strsplit(geo_line[1L], "\t", fixed = TRUE)[[1L]][-1L])
  if (length(geo) == 0L) return(character(0))
  out <- setNames(geo, geo)
  out <- c(out, setNames(geo, as.character(seq_along(geo))))
  title_line <- grep("^!Sample_title\\t", lines, value = TRUE)
  if (length(title_line) > 0L) {
    titles <- .clean_ids(strsplit(title_line[1L], "\t", fixed = TRUE)[[1L]][-1L])
    if (length(titles) == length(geo)) {
      for (key_set in list(titles, trimws(sub(":.*$", "", titles)),
                           sub("^(CONTROL|BENIGN|TUMOR|CANCER|CASE|CTRL|NORMAL|HEALTHY)_", "",
                               titles, ignore.case = TRUE))) {
        ok <- nzchar(key_set) & !duplicated(key_set)
        out <- c(out, setNames(geo[ok], key_set[ok]))
      }
    }
  }
  out[!duplicated(names(out))]
}

.remap_samples_from_series <- function(sample_ids, acc, expr_path) {
  paths <- .series_paths_for(acc, expr_path)
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0L) return(sample_ids)
  for (p in paths) {
    mp  <- .parse_series_metadata_map(p)
    if (length(mp) == 0L) next
    mapped <- unname(mp[sample_ids])
    hit    <- !is.na(mapped)
    if (sum(hit) < length(sample_ids) * 0.5) {
      norm_mp <- setNames(unname(mp), toupper(gsub("[\\-_\\. ]", "", names(mp))))
      mapped2 <- unname(norm_mp[toupper(gsub("[\\-_\\. ]", "", sample_ids))])
      hit2    <- !is.na(mapped2)
      if (sum(hit2) > sum(hit)) { mapped <- mapped2; hit <- hit2 }
    }
    if (sum(hit) >= length(sample_ids) * 0.5) {
      sample_ids[hit] <- mapped[hit]
      return(sample_ids)
    }
  }
  sample_ids
}

.dt_to_sample_feature_matrix <- function(dt, expr_path, acc) {
  if (ncol(dt) < 2L || nrow(dt) < 1L) stop("Expression matrix too small: ", expr_path)
  setnames(dt, .clean_ids(names(dt)))
  id_col  <- names(dt)[1L]
  id_vals <- .clean_ids(dt[[1L]])
  mat_dt  <- dt[, -1L, with = FALSE]
  ann_cols <- intersect(names(mat_dt),
                        c("precursor", "Precursor", "sequence", "Sequence",
                          "description", "Description", "Gene.title", "Gene.symbol"))
  if (length(ann_cols) > 0L)
    mat_dt <- mat_dt[, setdiff(names(mat_dt), ann_cols), with = FALSE]
  for (j in seq_len(ncol(mat_dt)))
    if (!is.numeric(mat_dt[[j]]))
      set(mat_dt, j = j, value = suppressWarnings(as.numeric(mat_dt[[j]])))
  X <- as.matrix(mat_dt); storage.mode(X) <- "numeric"
  feature_id_names <- c("miRNA", "#miRNA", "feature_id", "ID_REF", "miRNA_ID",
                        "probe_id", "Gene", "gene", "Symbol", "symbol",
                        "Probe", "probe", "miRBase_ID", "feature", "Gene ID", "ID", "id")
  sample_id_names  <- c("sample_id", "SampleID", "Sample.Name", "title",
                        "geo_accession", "GSM", "Sample", "sample",
                        "Subject ID", "subject_id", "SubjectID", "library_name")
  id_mir_frac  <- .mirna_fraction(id_vals)
  col_mir_frac <- .mirna_fraction(names(mat_dt))
  id_is_feature <- id_col %in% feature_id_names ||
    (id_mir_frac > 0.3 && !(id_col %in% sample_id_names && col_mir_frac > 0.3))
  id_is_sample  <- id_col %in% sample_id_names && !id_is_feature
  if (id_is_feature) {
    rownames(X) <- id_vals; X <- t(X)
  } else if (id_is_sample) {
    rownames(X) <- id_vals
  } else {
    rownames(X) <- id_vals
    if (id_mir_frac > 0.3 && col_mir_frac < id_mir_frac) X <- t(X)
  }
  rownames(X) <- .remap_samples_from_series(rownames(X), acc, expr_path)
  if (is.null(colnames(X))) stop("No feature names after parsing ", acc)
  X
}

.read_series_matrix <- function(path, acc) {
  lines   <- .read_lines_any(path)
  i_begin <- grep("!series_matrix_table_begin", lines, ignore.case = TRUE)
  i_end   <- grep("!series_matrix_table_end",   lines, ignore.case = TRUE)
  if (length(i_begin) == 0L || length(i_end) == 0L) stop("series matrix: no table markers")
  body <- lines[(i_begin[1L] + 1L):(i_end[1L] - 1L)]
  if (length(body) < 2L) stop("series matrix: empty table")
  dt <- fread(text = paste(body, collapse = "\n"), sep = "\t", quote = "\"",
              fill = TRUE, header = TRUE, check.names = FALSE,
              na.strings = c("", "NA", "NaN", "null", "N/A", "Undetermined"))
  .dt_to_sample_feature_matrix(dt, path, acc)
}

.read_soft_table <- function(path, acc) {
  lines   <- .read_lines_any(path)
  i_begin <- grep("!(dataset|series_matrix)_table_begin", lines, ignore.case = TRUE)
  i_end   <- grep("!(dataset|series_matrix)_table_end",   lines, ignore.case = TRUE)
  if (length(i_begin) == 0L || length(i_end) == 0L) stop("SOFT: no table markers")
  body <- lines[(i_begin[1L] + 1L):(i_end[1L] - 1L)]
  if (length(body) < 2L) stop("SOFT: empty table")
  dt <- fread(text = paste(body, collapse = "\n"), sep = "\t", quote = "\"",
              fill = TRUE, header = TRUE, check.names = FALSE,
              na.strings = c("", "NA", "NaN", "null", "N/A", "Undetermined"))
  .dt_to_sample_feature_matrix(dt, path, acc)
}

.tar_cache <- new.env(parent = emptyenv())
.tar_root  <- tempfile("loader_tar_")
dir.create(.tar_root, showWarnings = FALSE, recursive = TRUE)

.extract_tar_once <- function(path) {
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (exists(key, envir = .tar_cache, inherits = FALSE))
    return(get(key, envir = .tar_cache, inherits = FALSE))
  out_dir <- file.path(.tar_root, tools::file_path_sans_ext(basename(path)))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  utils::untar(path, exdir = out_dir)
  assign(key, out_dir, envir = .tar_cache)
  out_dir
}

.read_tar_expression <- function(path, acc) {
  out_dir <- .extract_tar_once(path)
  files   <- list.files(out_dir, pattern = "\\.(txt|tsv|csv)(\\.gz)?$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0L) stop("RAW.tar: no per-sample tables")
  rows <- list()
  for (f in files) {
    dt <- tryCatch(.read_plain_dt(f, acc = acc), error = function(e) NULL)
    if (is.null(dt) || ncol(dt) < 2L || nrow(dt) < 1L) next
    setnames(dt, .clean_ids(names(dt)))
    feat_col  <- intersect(c("id","ID","ID_REF","miRNA","#miRNA","feature_id","Gene ID","Gene","Probe"), names(dt))[1L]
    if (is.na(feat_col)) feat_col <- names(dt)[1L]
    value_col <- intersect(c("counts","count","VALUE","Value","value","expression","signal","Signal"), names(dt))[1L]
    if (is.na(value_col)) {
      num_cols  <- names(dt)[vapply(dt, is.numeric, logical(1L))]
      num_cols  <- setdiff(num_cols, feat_col)
      if (length(num_cols) == 0L) next
      value_col <- num_cols[1L]
    }
    v <- suppressWarnings(as.numeric(dt[[value_col]]))
    names(v) <- .clean_ids(dt[[feat_col]])
    sid <- sub("\\.(txt|tsv|csv)(\\.gz)?$", "", basename(f), ignore.case = TRUE)
    sid <- sub("^(GSM[0-9]+).*", "\\1", sid)
    rows[[sid]] <- v[!is.na(names(v)) & nzchar(names(v))]
  }
  if (length(rows) == 0L) stop("RAW.tar: no usable tables")
  features <- unique(unlist(lapply(rows, names), use.names = FALSE))
  X <- matrix(NA_real_, nrow = length(rows), ncol = length(features),
              dimnames = list(names(rows), features))
  for (i in seq_along(rows)) X[i, names(rows[[i]])] <- rows[[i]]
  X
}

.read_raw_expression_path <- function(path, acc) {
  if (grepl("_RAW\\.tar$", path, ignore.case = TRUE))
    return(.read_tar_expression(path, acc))
  if (grepl("series_matrix\\.txt\\.gz$", path, ignore.case = TRUE))
    return(.read_series_matrix(path, acc))
  if (grepl("family\\.soft\\.gz$", path, ignore.case = TRUE))
    return(.read_soft_table(path, acc))
  .dt_to_sample_feature_matrix(.read_plain_dt(path, acc = acc), path, acc)
}

.combine_by_sample_rows <- function(mats) {
  features <- unique(unlist(lapply(mats, colnames), use.names = FALSE))
  out <- lapply(mats, function(X) {
    Y <- matrix(NA_real_, nrow = nrow(X), ncol = length(features),
                dimnames = list(rownames(X), features))
    Y[, colnames(X)] <- X; Y
  })
  X <- do.call(rbind, out)
  if (any(duplicated(rownames(X)))) X <- X[!duplicated(rownames(X)), , drop = FALSE]
  X
}

.postprocess_expression <- function(X, inv_row, expr_path) {
  scale_attr    <- tolower(trimws(as.character(inv_row$scale %||% "")))
  log_base_attr <- trimws(as.character(inv_row$log_base %||% "none"))
  td_attr       <- tolower(trimws(as.character(inv_row$transform_direction %||% "none")))
  if (length(scale_attr)    == 0L || is.na(scale_attr))    scale_attr    <- ""
  if (length(log_base_attr) == 0L || is.na(log_base_attr)) log_base_attr <- "none"
  if (length(td_attr)       == 0L || is.na(td_attr))       td_attr       <- "none"
  scale_state <- if (scale_attr %in% c("count", "rpm", "abundance")) "abundance" else "log"
  colnames(X) <- .normalize_feature_ids(colnames(X))
  X <- resolve_and_collapse_features(X, scale_state = scale_state)
  X <- normalize_to_canonical_mirna_namespace(X, scale_state = scale_state)
  col_keep <- rep(TRUE, ncol(X))
  for (j in seq_len(ncol(X))) {
    bad <- !is.finite(X[, j])
    if (all(bad)) { col_keep[j] <- FALSE; next }
    if (any(bad)) {
      med <- suppressWarnings(stats::median(X[, j], na.rm = TRUE))
      X[bad, j] <- if (is.finite(med)) med else 0
    }
  }
  X <- X[, col_keep, drop = FALSE]
  attr(X, "scale")              <- scale_attr
  attr(X, "log_base")           <- log_base_attr
  attr(X, "transform_direction") <- td_attr
  attr(X, "expression_path")    <- expr_path
  X
}

.read_expression <- function(acc) {
  inv_row <- inv[accession == acc][1L]
  if (nrow(inv_row) == 0L) stop("No inventory row for ", acc)
  cands <- .candidates_for(acc, inv_row$present_in)
  if (length(cands) == 0L) stop("No expression matrix on disk for ", acc)
  nr_paths <- cands[grepl("_NR-gene_row\\.txt\\.gz$", basename(cands), ignore.case = TRUE)]
  if (length(nr_paths) > 1L) {
    mats <- lapply(nr_paths, function(p) .read_raw_expression_path(p, acc))
    X    <- .combine_by_sample_rows(mats)
    return(.postprocess_expression(X, inv_row, paste(basename(nr_paths), collapse = ";")))
  }
  errors <- character(0)
  for (expr_path in cands) {
    X <- tryCatch(.read_raw_expression_path(expr_path, acc),
                  error = function(e) {
                    errors <<- c(errors, sprintf("%s: %s", basename(expr_path), e$message))
                    NULL
                  })
    if (is.null(X)) next
    return(.postprocess_expression(X, inv_row, expr_path))
  }
  stop("All expression candidates failed for ", acc, ": ", paste(errors, collapse = " | "))
}

.to_abundance <- function(X) {
  scale <- tolower(trimws(as.character(attr(X, "scale") %||% "")))
  td    <- tolower(trimws(as.character(attr(X, "transform_direction") %||% "none")))
  if (length(scale) == 0L || is.na(scale)) scale <- ""
  if (length(td)    == 0L || is.na(td))    td    <- "none"
  base  <- suppressWarnings(as.numeric(attr(X, "log_base") %||% "none"))
  if (length(base) == 0L || is.na(base) || !is.finite(base) || base <= 0) base <- 2
  if (scale %in% c("count", "rpm", "abundance") || td == "none") {
    X[!is.finite(X)] <- 0; X[X < 0] <- 0; return(X)
  }
  if (scale == "ct" || td == "negate_inverse") {
    C   <- 40
    out <- base ^ (C - X)   # monotone-decreasing in Ct; >= 0 for X <= C
    out[!is.finite(out)] <- 0; out[out < 0] <- 0; return(out)
  }
  out <- base ^ X
  out[!is.finite(out)] <- 0; out[out < 0] <- 0; out
}

# ---------------------------------------------------------------------------
# 3. Per-cohort (X, y) loader
# ---------------------------------------------------------------------------

load_notes <- new.env(parent = emptyenv())
.set_load_note <- function(acc, reason) assign(acc, reason, envir = load_notes)
.get_load_note <- function(acc) {
  if (exists(acc, envir = load_notes, inherits = FALSE))
    return(get(acc, envir = load_notes, inherits = FALSE))
  "expression load failure or no usable label set"
}

.load_cohort <- function(acc, modality) {
  inv_row  <- inv[accession == acc][1L]
  dict_sub <- dict[accession == acc]
  meta_path <- .metadata_path_for(acc, inv_row$present_in, dict_sub)
  if (is.na(meta_path)) {
    .set_load_note(acc, "metadata file missing")
    cat(sprintf("  - %s: metadata file missing\n", acc)); return(NULL)
  }
  meta <- fread(meta_path, sep = "\t", quote = "", fill = TRUE,
                na.strings = c("", "NA", "NaN", "null", "N/A"))
  X <- tryCatch(.read_expression(acc), error = function(e) {
    msg <- sprintf("expression load failed: %s", e$message)
    .set_load_note(acc, msg); cat(sprintf("  - %s: %s\n", acc, msg)); NULL
  })
  if (is.null(X)) return(NULL)
  if (nrow(dict_sub) == 0L) {
    .set_load_note(acc, "no outcome dictionary entry")
    cat(sprintf("  - %s: no outcome dictionary entry\n", acc)); return(NULL)
  }
  res  <- apply_outcome_dictionary(meta, rownames(X), dict_sub, allow_row_position_fallback = FALSE)
  y    <- res$y_bin
  keep <- !is.na(y)
  if (sum(keep) < 20L || length(unique(y[keep])) < 2L) {
    msg <- sprintf("usable n=%d (cancer=%d, control=%d)", sum(keep),
                   sum(y == 1L, na.rm = TRUE), sum(y == 0L, na.rm = TRUE))
    .set_load_note(acc, msg); cat(sprintf("  - %s: %s; skipping\n", acc, msg)); return(NULL)
  }
  x_attrs <- attributes(X)[c("scale","log_base","transform_direction","expression_path")]
  X <- X[keep, , drop = FALSE]
  for (nm in names(x_attrs)) attr(X, nm) <- x_attrs[[nm]]
  y <- as.integer(y[keep])
  X_abund <- .to_abundance(X)
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  cat(sprintf("  - %s: usable n=%d (cancer=%d, control=%d), match=%s\n",
              acc, length(y), n1, n0, res$audit$matched_by))
  list(X = X_abund, y = y, accession = acc, modality = modality,
       n_cancer = n1, n_control = n0,
       expression_path = attr(X, "expression_path"))
}

# ---------------------------------------------------------------------------
# 4. Statistical primitives
# ---------------------------------------------------------------------------

.feature_auc <- function(x, y) {
  if (length(unique(y)) < 2L) return(NA_real_)
  if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) return(NA_real_)
  a <- tryCatch(as.numeric(pROC::auc(y, x, quiet = TRUE, direction = "<",
                                     levels = c(0L, 1L))),
                error = function(e) NA_real_)
  if (!is.finite(a)) return(NA_real_)
  max(a, 1 - a)
}

.pool_train_aucs <- function(train_data, common_features) {
  per_cohort <- lapply(train_data, function(td) {
    Xc   <- td$X[, common_features, drop = FALSE]
    aucs <- vapply(seq_along(common_features),
                   function(j) .feature_auc(Xc[, j], td$y), numeric(1L))
    names(aucs) <- common_features
    aucs
  })
  M <- do.call(cbind, per_cohort)
  apply(M, 1L, function(v) median(v, na.rm = TRUE))
}

.hanley_mcneil_ci <- function(auc, n1, n0, level = 0.95) {
  if (!is.finite(auc) || n1 < 1L || n0 < 1L)
    return(list(lo = NA_real_, hi = NA_real_, se = NA_real_))
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  v  <- (auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) +
         (n0 - 1) * (q2 - auc^2)) / (n1 * n0)
  if (!is.finite(v) || v < 0) v <- 0
  se <- sqrt(v)
  z  <- qnorm(1 - (1 - level) / 2)
  list(lo = max(0, auc - z * se), hi = min(1, auc + z * se), se = se)
}

# Run rCLR + signed panel scoring on a single held-out cohort, given
# training data and a panel. Returns list(score, auc, ci).
.score_held_out <- function(held_data, train_data, panel, PANEL_SIZE = 20L) {
  X_held <- held_data$X[, panel, drop = FALSE]
  Z <- tryCatch(ws_rclr_trimmed(X_held, min_centering_size = 8L),
                error = function(e) NULL)
  if (is.null(Z)) return(list(score = NULL, auc = NA_real_, ci = list(lo=NA,hi=NA,se=NA),
                               notes = "rCLR failed"))
  train_signs <- vapply(panel, function(j) {
    vals <- vapply(train_data, function(td) {
      if (!j %in% colnames(td$X)) return(NA_real_)
      Xt <- td$X[, panel[panel %in% colnames(td$X)], drop = FALSE]
      if (!j %in% colnames(Xt)) return(NA_real_)
      Zt <- tryCatch(ws_rclr_trimmed(Xt, min_centering_size = 8L), error = function(e) NULL)
      if (is.null(Zt) || !j %in% colnames(Zt)) return(NA_real_)
      tryCatch(as.numeric(pROC::auc(td$y, Zt[, j], quiet = TRUE,
                                    direction = "<", levels = c(0L, 1L))),
               error = function(e) NA_real_)
    }, numeric(1L))
    med <- median(vals, na.rm = TRUE)
    if (!is.finite(med)) 1 else if (med < 0.5) -1 else 1
  }, numeric(1L))
  Z_sign <- sweep(Z, 2L, train_signs, `*`)
  score  <- rowSums(Z_sign, na.rm = TRUE)
  auc_v  <- tryCatch(as.numeric(pROC::auc(held_data$y, score, quiet = TRUE,
                                           direction = "<", levels = c(0L, 1L))),
                     error = function(e) NA_real_)
  if (!is.finite(auc_v)) auc_v <- NA_real_
  ci <- .hanley_mcneil_ci(auc_v, held_data$n_cancer, held_data$n_control)
  list(score = score, auc = auc_v, ci = ci, notes = "ok")
}

# Build a panel from training data and common feature intersection.
.build_panel <- function(train_data, common_features, PANEL_SIZE = 20L) {
  if (length(common_features) < PANEL_SIZE * 2L) return(NULL)
  pooled <- .pool_train_aucs(train_data, common_features)
  pooled <- pooled[is.finite(pooled)]
  if (length(pooled) < PANEL_SIZE) return(NULL)
  names(sort(pooled, decreasing = TRUE))[seq_len(PANEL_SIZE)]
}

cat("[cohort_loader.R] Shared loading helpers sourced.\n")
