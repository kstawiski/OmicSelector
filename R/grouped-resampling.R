#' @title Build Grouped Stratified Folds
#'
#' @description
#' Creates validation folds that keep all rows from the same biological specimen
#' or other grouping unit in the same fold. This is intended for public-omics
#' audits where duplicated specimens, repeated aliquots, or reused biological IDs
#' can otherwise leak across resampling splits.
#'
#' @param y Binary outcome vector. The second factor level is treated as positive.
#' @param group_id Optional grouping vector. If omitted, each row is its own group.
#' @param strata Optional stratification vector. Defaults to \code{y}.
#' @param n_folds Number of folds.
#' @param seed Random seed used for within-stratum group ordering.
#' @param mixed_group_action What to do when a group contains both outcome
#'   classes: \code{"stop"}, \code{"warn"}, or \code{"allow"}.
#'
#' @return A list of integer test-row indices with class \code{"os_grouped_folds"}.
#' @export
os_make_grouped_stratified_folds <- function(y, group_id = NULL, strata = NULL,
                                             n_folds = 5L, seed = 1L,
                                             mixed_group_action = c("stop", "warn", "allow")) {
  mixed_group_action <- match.arg(mixed_group_action)
  y_num <- .os_binary_y(y)
  n <- length(y_num)
  if (is.null(group_id)) group_id <- seq_len(n)
  if (is.null(strata)) strata <- ifelse(y_num == 1L, "positive", "negative")
  if (length(group_id) != n || length(strata) != n) {
    stop("`y`, `group_id`, and `strata` must have the same length.", call. = FALSE)
  }
  if (anyNA(group_id) || anyNA(strata)) {
    stop("`group_id` and `strata` cannot contain missing values.", call. = FALSE)
  }
  n_folds <- as.integer(n_folds)
  if (is.na(n_folds) || n_folds < 2L) {
    stop("`n_folds` must be an integer >= 2.", call. = FALSE)
  }

  group_id <- as.character(group_id)
  strata <- as.character(strata)
  groups <- unique(group_id)
  group_rows <- vector("list", length(groups))
  meta <- data.frame(
    group_id = groups,
    n = integer(length(groups)),
    n_positive = integer(length(groups)),
    n_negative = integer(length(groups)),
    stratum = character(length(groups)),
    mixed_outcome = logical(length(groups)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(groups)) {
    idx <- which(group_id == groups[[i]])
    group_rows[[i]] <- idx
    yi <- unique(y_num[idx])
    si <- unique(strata[idx])
    meta$n[i] <- length(idx)
    meta$n_positive[i] <- sum(y_num[idx] == 1L)
    meta$n_negative[i] <- sum(y_num[idx] == 0L)
    meta$mixed_outcome[i] <- length(yi) > 1L
    meta$stratum[i] <- paste(sort(si), collapse = "|")
  }

  if (any(meta$mixed_outcome)) {
    msg <- paste0(sum(meta$mixed_outcome), " group(s) contain both outcome classes.")
    if (mixed_group_action == "stop") stop(msg, call. = FALSE)
    if (mixed_group_action == "warn") warning(msg, call. = FALSE)
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) .Random.seed <<- old_seed
  }, add = TRUE)
  set.seed(seed)

  fold_for_group <- setNames(integer(nrow(meta)), meta$group_id)
  total_counts <- integer(n_folds)
  strata_levels <- sort(unique(meta$stratum))

  for (s in strata_levels) {
    rows <- which(meta$stratum == s)
    rows <- rows[order(-meta$n[rows], sample.int(length(rows)))]
    stratum_counts <- integer(n_folds)
    for (row in rows) {
      candidates <- which(stratum_counts == min(stratum_counts))
      fold <- candidates[which.min(total_counts[candidates])]
      fold_for_group[meta$group_id[row]] <- fold
      stratum_counts[fold] <- stratum_counts[fold] + meta$n[row]
      total_counts[fold] <- total_counts[fold] + meta$n[row]
    }
  }

  folds <- lapply(seq_len(n_folds), function(fold) {
    which(fold_for_group[group_id] == fold)
  })
  names(folds) <- paste0("fold", seq_len(n_folds))
  attr(folds, "group_metadata") <- meta
  attr(folds, "diagnostics") <- .os_fold_diagnostics(folds, y_num, group_id, strata)
  class(folds) <- c("os_grouped_folds", "list")
  folds
}

#' @title Validate Grouped Folds
#'
#' @param folds List of integer test-row indices.
#' @param y Optional outcome vector. When supplied, coverage is checked against
#'   \code{length(y)}.
#' @param group_id Optional grouping vector used to verify that groups are not split.
#' @param strata Optional strata vector used to report fold balance.
#'
#' @return A data frame of validation checks.
#' @export
os_validate_folds <- function(folds, y = NULL, group_id = NULL, strata = NULL) {
  if (!is.list(folds) || length(folds) < 2L) {
    stop("`folds` must be a list with at least two folds.", call. = FALSE)
  }
  idx <- unlist(folds, use.names = FALSE)
  n <- if (!is.null(y)) length(y) else max(idx)
  in_range <- length(idx) > 0L && all(idx >= 1L & idx <= n)
  complete <- in_range && setequal(idx, seq_len(n))
  unique_once <- in_range && !anyDuplicated(idx) && length(idx) == n

  group_pass <- TRUE
  split_groups <- character()
  if (!is.null(group_id)) {
    if (length(group_id) != n) stop("`group_id` must have length `length(y)` or `max(unlist(folds))`.", call. = FALSE)
    fold_id <- rep(seq_along(folds), lengths(folds))
    fold_by_row <- rep(NA_integer_, n)
    fold_by_row[idx] <- fold_id
    by_group <- split(fold_by_row, as.character(group_id))
    split_groups <- names(by_group)[vapply(by_group, function(z) length(unique(stats::na.omit(z))) > 1L, logical(1L))]
    group_pass <- length(split_groups) == 0L
  }

  balance <- ""
  if (!is.null(y)) {
    y_num <- .os_binary_y(y)
    diag <- .os_fold_diagnostics(folds, y_num, if (is.null(group_id)) seq_len(n) else group_id, strata)
    balance <- paste(sprintf("%s:n=%s,pos=%s", diag$fold, diag$n, diag$n_positive), collapse = ";")
  }

  out <- data.frame(
    check = c("indices_in_range", "complete_coverage", "unique_test_assignment", "no_group_split", "fold_balance"),
    pass = c(in_range, complete, unique_once, group_pass, TRUE),
    value = c(
      as.character(in_range),
      paste0(length(unique(idx)), "/", n),
      paste0(length(idx), " assignments"),
      if (group_pass) "none" else paste(utils::head(split_groups, 10L), collapse = ","),
      balance
    ),
    stringsAsFactors = FALSE
  )
  class(out) <- c("os_fold_validation", "data.frame")
  out
}

#' @title Estimate AUC over Grouped Resampling Folds
#'
#' @param y Binary outcome vector.
#' @param score Numeric score vector where larger values indicate the positive class.
#' @param folds List of integer test-row indices.
#'
#' @return A list with pooled and per-fold AUC values.
#' @export
os_grouped_resample_auc <- function(y, score, folds) {
  y_num <- .os_binary_y(y)
  if (length(score) != length(y_num)) stop("`score` and `y` must have the same length.", call. = FALSE)
  score <- as.numeric(score)
  fold_auc <- data.frame(
    fold = names(folds),
    n = lengths(folds),
    auc = vapply(folds, function(idx) .auc_fast(y_num[idx], score[idx]), numeric(1L)),
    stringsAsFactors = FALSE
  )
  out <- list(
    pooled_auc = .auc_fast(y_num[unlist(folds, use.names = FALSE)], score[unlist(folds, use.names = FALSE)]),
    fold_auc = fold_auc,
    validation = os_validate_folds(folds, y = y)
  )
  class(out) <- c("os_grouped_resample_auc", "list")
  out
}

.os_fold_diagnostics <- function(folds, y_num, group_id, strata = NULL) {
  if (is.null(strata)) strata <- ifelse(y_num == 1L, "positive", "negative")
  data.frame(
    fold = names(folds),
    n = lengths(folds),
    n_positive = vapply(folds, function(idx) sum(y_num[idx] == 1L), integer(1L)),
    n_negative = vapply(folds, function(idx) sum(y_num[idx] == 0L), integer(1L)),
    n_groups = vapply(folds, function(idx) length(unique(as.character(group_id)[idx])), integer(1L)),
    n_strata = vapply(folds, function(idx) length(unique(as.character(strata)[idx])), integer(1L)),
    stringsAsFactors = FALSE
  )
}
