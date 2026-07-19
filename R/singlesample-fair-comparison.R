#' Corrected repeated-CV inference for a paired method effect
#'
#' Summarises paired effects observed in the same repeated cross-validation
#' strata using the Nadeau--Bengio correction. The standard error is
#' `sqrt((1 / m + mean(n_test / n_train)) * var(effect))`, with `m - 1`
#' degrees of freedom. The function fails closed unless the observed stratum
#' identifiers match the complete prespecified set exactly.
#'
#' The returned directional tests are fixed before outcomes are inspected:
#' `p_above_margin` tests `H0: effect <= margin`, while
#' `p_below_minus_margin` tests `H0: effect >= -margin`. `p_zero_two_sided`
#' tests a zero mean effect. When the corrected standard error is zero, the
#' confidence interval remains the point estimate but inferential p-values are
#' `NA`, rather than treating a degenerate repeated-CV sample as certain.
#'
#' @param effect Numeric paired effect for every shared CV stratum, ordered as
#'   method A minus method B.
#' @param n_test,n_train Positive test and training counts for each stratum.
#'   Use biological-group counts when profiles repeat within a group.
#' @param stratum_id Unique observed stratum identifiers.
#' @param expected_m Prespecified number of shared strata.
#' @param expected_strata Complete prespecified stratum identifiers. Their order
#'   may differ from `stratum_id`, but membership must match exactly.
#' @param margin Positive relevance margin. The default is `0.05` AUC units.
#' @param conf_level Confidence level for the nominal t interval.
#'
#' @return A named list containing the estimate, corrected uncertainty,
#'   confidence interval, zero and shifted-null p-values, correction factor,
#'   and exact stratum accounting.
#'
#' @references Nadeau C, Bengio Y. Inference for the generalization error.
#'   Machine Learning. 2003;52:239--281.
#' @export
singlesample_corrected_repeated_cv <- function(
    effect, n_test, n_train, stratum_id, expected_m, expected_strata,
    margin = 0.05, conf_level = 0.95) {
  expected_m <- .singlesample_scalar_integer(expected_m, "expected_m", 2L)
  margin <- .singlesample_scalar_number(margin, "margin", lower = 0,
                                         strict = TRUE)
  conf_level <- .singlesample_scalar_number(
    conf_level, "conf_level", lower = 0, upper = 1, strict = TRUE
  )

  effect <- as.numeric(effect)
  n_test <- as.numeric(n_test)
  n_train <- as.numeric(n_train)
  observed <- .singlesample_stratum_ids(stratum_id, "stratum_id")
  expected <- .singlesample_stratum_ids(expected_strata, "expected_strata")

  lengths <- c(length(effect), length(n_test), length(n_train), length(observed))
  if (any(lengths != expected_m)) {
    stop(sprintf(
      "effect, counts, and stratum_id must each contain exactly expected_m=%d rows.",
      expected_m
    ), call. = FALSE)
  }
  if (length(expected) != expected_m) {
    stop("expected_strata must contain exactly expected_m identifiers.",
         call. = FALSE)
  }
  if (anyDuplicated(observed)) {
    stop("stratum_id must not contain duplicates.", call. = FALSE)
  }
  if (anyDuplicated(expected)) {
    stop("expected_strata must not contain duplicates.", call. = FALSE)
  }
  if (!setequal(observed, expected)) {
    stop("stratum_id does not match the complete expected_strata set.",
         call. = FALSE)
  }
  if (any(!is.finite(effect))) {
    stop("effect must be finite in every expected stratum.", call. = FALSE)
  }
  if (any(!is.finite(n_test)) || any(!is.finite(n_train)) ||
      any(n_test <= 0) || any(n_train <= 0)) {
    stop("n_test and n_train must be finite and strictly positive.",
         call. = FALSE)
  }

  ratio <- n_test / n_train
  correction <- 1 / expected_m + mean(ratio)
  estimate <- mean(effect)
  se <- sqrt(max(0, correction * stats::var(effect)))
  df <- expected_m - 1L
  critical <- stats::qt((1 + conf_level) / 2, df = df)
  inference_available <- is.finite(se) && se > 0

  p_zero <- p_above <- p_below <- NA_real_
  if (inference_available) {
    p_zero <- 2 * stats::pt(-abs(estimate / se), df = df)
    p_above <- stats::pt((estimate - margin) / se, df = df,
                         lower.tail = FALSE)
    p_below <- stats::pt((estimate + margin) / se, df = df,
                         lower.tail = TRUE)
  }

  list(
    complete = TRUE,
    observed_m = expected_m,
    expected_m = expected_m,
    estimate = estimate,
    se = se,
    ci_low = estimate - critical * se,
    ci_high = estimate + critical * se,
    conf_level = conf_level,
    p_zero_two_sided = p_zero,
    p_above_margin = p_above,
    p_below_minus_margin = p_below,
    margin = margin,
    df = df,
    correction = correction,
    mean_test_train_ratio = mean(ratio),
    inference_available = inference_available,
    stratum_id = observed,
    expected_strata = expected
  )
}

#' Matched group-primary AUC difference across repeated CV strata
#'
#' Computes fixed-direction AUC for two methods on exactly aligned held-out
#' rows, collapses repeated profiles to their biological/provenance group within
#' each stratum, and applies [singlesample_corrected_repeated_cv()] to the
#' stratum-level paired differences. Larger scores must already indicate the
#' positive class according to a training-frozen orientation.
#'
#' @param y Binary held-out outcome.
#' @param score_a,score_b Finite held-out scores for methods A and B on the same
#'   rows.
#' @param stratum_id Repeated-CV stratum per row, such as `seed::fold`.
#' @param sample_id Held-out profile identifier. Each profile may occur only
#'   once in a stratum.
#' @param group_id Biological/provenance group identifier. Group labels must be
#'   homogeneous.
#' @param expected_strata Complete prespecified repeated-CV strata.
#' @param margin Positive relevance margin in AUC units.
#' @param conf_level Confidence level for the nominal corrected t interval.
#'
#' @return A list with one `strata` data frame and one corrected `summary` list.
#' @export
singlesample_matched_pair_auc <- function(
    y, score_a, score_b, stratum_id, sample_id, group_id, expected_strata,
    margin = 0.05, conf_level = 0.95) {
  lengths <- c(length(y), length(score_a), length(score_b),
               length(stratum_id), length(sample_id), length(group_id))
  if (!length(y) || length(unique(lengths)) != 1L) {
    stop("All matched prediction columns must have the same non-zero length.",
         call. = FALSE)
  }
  y <- .os_binary_y(y)
  score_a <- as.numeric(score_a)
  score_b <- as.numeric(score_b)
  if (any(!is.finite(score_a)) || any(!is.finite(score_b))) {
    stop("score_a and score_b must be finite on every matched row.",
         call. = FALSE)
  }
  stratum_id <- .singlesample_stratum_ids(stratum_id, "stratum_id")
  sample_id <- .singlesample_stratum_ids(sample_id, "sample_id")
  group_id <- .singlesample_stratum_ids(group_id, "group_id")
  expected_strata <- .singlesample_stratum_ids(
    expected_strata, "expected_strata"
  )
  if (anyDuplicated(expected_strata)) {
    stop("expected_strata must not contain duplicates.", call. = FALSE)
  }
  row_key <- paste(stratum_id, sample_id, sep = "\r")
  if (anyDuplicated(row_key)) {
    stop("Each sample_id must occur at most once in a CV stratum.",
         call. = FALSE)
  }
  if (!setequal(unique(stratum_id), expected_strata)) {
    stop("Observed stratum_id does not match the complete expected_strata set.",
         call. = FALSE)
  }
  group_y <- split(y, group_id)
  if (any(vapply(group_y, function(z) length(unique(z)) != 1L, logical(1L)))) {
    stop("A biological/provenance group contains both outcome classes.",
         call. = FALSE)
  }
  total_groups <- length(unique(group_id))

  rows <- lapply(expected_strata, function(stratum) {
    index <- which(stratum_id == stratum)
    groups <- unique(group_id[index])
    collapsed_y <- integer(length(groups))
    collapsed_a <- collapsed_b <- numeric(length(groups))
    for (i in seq_along(groups)) {
      group_index <- index[group_id[index] == groups[[i]]]
      values <- unique(y[group_index])
      if (length(values) != 1L) {
        stop("A group has discordant labels within a CV stratum.",
             call. = FALSE)
      }
      collapsed_y[[i]] <- values[[1L]]
      collapsed_a[[i]] <- mean(score_a[group_index])
      collapsed_b[[i]] <- mean(score_b[group_index])
    }
    if (length(unique(collapsed_y)) != 2L) {
      stop("A required group-collapsed CV stratum lacks an outcome class.",
           call. = FALSE)
    }
    n_test <- length(groups)
    n_train <- total_groups - n_test
    if (n_train <= 0L) {
      stop("A CV stratum has no remaining training groups.", call. = FALSE)
    }
    auc_a <- .auc_fast(collapsed_y, collapsed_a)
    auc_b <- .auc_fast(collapsed_y, collapsed_b)
    data.frame(
      stratum_id = stratum,
      auc_a = auc_a,
      auc_b = auc_b,
      effect = auc_a - auc_b,
      n_test_groups = n_test,
      n_train_groups = n_train,
      stringsAsFactors = FALSE
    )
  })
  strata <- do.call(rbind, rows)
  summary <- singlesample_corrected_repeated_cv(
    effect = strata$effect,
    n_test = strata$n_test_groups,
    n_train = strata$n_train_groups,
    stratum_id = strata$stratum_id,
    expected_m = length(expected_strata),
    expected_strata = expected_strata,
    margin = margin,
    conf_level = conf_level
  )
  list(strata = strata, summary = summary)
}

#' Enumerate frozen within-method pairs
#'
#' Enumerates every unordered pair of methods whose roster estimand is
#' `"within"`. Pair direction follows roster row order and never depends on an
#' outcome, score, eligibility, or performance column.
#'
#' @param roster Method roster containing unique `method_id` and `estimand`
#'   columns. Defaults to [singlesample_method_roster()].
#'
#' @return A data frame with one row per unordered pair and columns `pair_id`,
#'   `method_a`, `method_b`, `roster_order_a`, and `roster_order_b`.
#' @export
singlesample_within_method_pairs <- function(
    roster = singlesample_method_roster()) {
  roster <- as.data.frame(roster, stringsAsFactors = FALSE)
  .singlesample_require_columns(roster, c("method_id", "estimand"), "roster")
  ids <- as.character(roster$method_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("roster$method_id must be non-missing, non-empty, and unique.",
         call. = FALSE)
  }
  estimand <- as.character(roster$estimand)
  if (anyNA(estimand) || any(!estimand %in% c("within", "transfer"))) {
    stop("roster$estimand must contain only 'within' or 'transfer'.",
         call. = FALSE)
  }
  within_order <- which(estimand == "within")
  if (length(within_order) < 2L) {
    return(data.frame(
      pair_id = character(), method_a = character(), method_b = character(),
      roster_order_a = integer(), roster_order_b = integer(),
      stringsAsFactors = FALSE
    ))
  }
  pair_index <- utils::combn(within_order, 2L)
  method_a <- ids[pair_index[1L, ]]
  method_b <- ids[pair_index[2L, ]]
  data.frame(
    pair_id = paste(method_a, method_b, sep = "::"),
    method_a = method_a,
    method_b = method_b,
    roster_order_a = pair_index[1L, ],
    roster_order_b = pair_index[2L, ],
    stringsAsFactors = FALSE
  )
}

#' Construct fixed eligibility-only complete-support panels
#'
#' Builds complete method-by-unit panels without reading outcome or performance
#' columns. At each coverage threshold, methods are retained when eligible in
#' at least the threshold proportion of primary units; units are then retained
#' only when every retained method is eligible. Rank summaries are prohibited
#' when fewer than `min_methods` methods or `min_units` units remain.
#'
#' The eligibility table must contain exactly one Boolean cell for every
#' within-method by primary-unit combination. Additional known transfer methods
#' or non-primary units are ignored.
#'
#' @param eligibility Data frame containing method, unit, and eligibility
#'   columns.
#' @param roster Frozen method roster; defaults to
#'   [singlesample_method_roster()]. Roster row order is preserved.
#' @param thresholds Coverage proportions. Defaults to 0.50, 0.80, and 0.95.
#' @param primary_units Optional ordered vector of primary unit identifiers. If
#'   `NULL`, first-appearance order in `eligibility` is used.
#' @param method_col,unit_col,eligible_col Column names in `eligibility`.
#' @param min_methods,min_units Minimum panel dimensions required to permit rank
#'   and top-k summaries.
#'
#' @return A list with `summary`, `method_membership`, and `unit_membership`
#'   data frames. Panel construction uses no columns beyond the identifiers and
#'   Boolean eligibility flag.
#' @export
singlesample_complete_support_panels <- function(
    eligibility, roster = singlesample_method_roster(),
    thresholds = c(0.50, 0.80, 0.95), primary_units = NULL,
    method_col = "method_id", unit_col = "unit_id",
    eligible_col = "eligible", min_methods = 5L, min_units = 5L) {
  eligibility <- as.data.frame(eligibility, stringsAsFactors = FALSE)
  roster <- as.data.frame(roster, stringsAsFactors = FALSE)
  .singlesample_require_columns(
    eligibility, c(method_col, unit_col, eligible_col), "eligibility"
  )
  .singlesample_require_columns(roster, c("method_id", "estimand"), "roster")
  min_methods <- .singlesample_scalar_integer(min_methods, "min_methods", 1L)
  min_units <- .singlesample_scalar_integer(min_units, "min_units", 1L)
  thresholds <- as.numeric(thresholds)
  if (!length(thresholds) || any(!is.finite(thresholds)) ||
      any(thresholds <= 0 | thresholds > 1) || anyDuplicated(thresholds)) {
    stop("thresholds must be unique finite proportions in (0, 1].",
         call. = FALSE)
  }

  roster_ids <- as.character(roster$method_id)
  if (anyNA(roster_ids) || any(!nzchar(roster_ids)) ||
      anyDuplicated(roster_ids)) {
    stop("roster$method_id must be non-missing, non-empty, and unique.",
         call. = FALSE)
  }
  roster_estimand <- as.character(roster$estimand)
  if (anyNA(roster_estimand) ||
      any(!roster_estimand %in% c("within", "transfer"))) {
    stop("roster$estimand must contain only 'within' or 'transfer'.",
         call. = FALSE)
  }
  within_ids <- roster_ids[roster_estimand == "within"]
  if (!length(within_ids)) {
    stop("roster contains no within-estimand methods.", call. = FALSE)
  }

  eligibility[[method_col]] <- as.character(eligibility[[method_col]])
  eligibility[[unit_col]] <- as.character(eligibility[[unit_col]])
  if (anyNA(eligibility[[method_col]]) ||
      any(!nzchar(eligibility[[method_col]])) ||
      anyNA(eligibility[[unit_col]]) || any(!nzchar(eligibility[[unit_col]]))) {
    stop("eligibility method and unit identifiers must be non-empty.",
         call. = FALSE)
  }
  unknown <- setdiff(unique(eligibility[[method_col]]), roster_ids)
  if (length(unknown)) {
    stop("eligibility contains methods absent from roster: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  eligible_value <- eligibility[[eligible_col]]
  if (!is.logical(eligible_value) || anyNA(eligible_value)) {
    stop("eligibility flag must be logical and non-missing.", call. = FALSE)
  }
  if (is.null(primary_units)) {
    primary_units <- unique(eligibility[[unit_col]])
  } else {
    primary_units <- as.character(primary_units)
  }
  if (!length(primary_units) || anyNA(primary_units) ||
      any(!nzchar(primary_units)) || anyDuplicated(primary_units)) {
    stop("primary_units must be a unique, non-empty identifier vector.",
         call. = FALSE)
  }

  keep <- eligibility[[method_col]] %in% within_ids &
    eligibility[[unit_col]] %in% primary_units
  grid <- eligibility[keep, c(method_col, unit_col, eligible_col), drop = FALSE]
  key <- paste(grid[[method_col]], grid[[unit_col]], sep = "\r")
  if (anyDuplicated(key)) {
    stop("eligibility contains duplicate within-method by primary-unit cells.",
         call. = FALSE)
  }
  expected_key <- as.vector(outer(within_ids, primary_units, paste, sep = "\r"))
  if (!setequal(key, expected_key) || length(key) != length(expected_key)) {
    stop("eligibility must contain the complete within-method by primary-unit grid.",
         call. = FALSE)
  }

  matrix_index <- match(expected_key, key)
  eligible_matrix <- matrix(
    grid[[eligible_col]][matrix_index], nrow = length(within_ids),
    ncol = length(primary_units), byrow = FALSE,
    dimnames = list(within_ids, primary_units)
  )
  coverage_n <- rowSums(eligible_matrix)
  n_primary <- length(primary_units)

  summaries <- methods <- units <- vector("list", length(thresholds))
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[[i]]
    required_units <- ceiling(threshold * n_primary)
    included_method <- coverage_n >= required_units
    retained <- within_ids[included_method]
    included_unit <- if (length(retained)) {
      colSums(eligible_matrix[retained, , drop = FALSE]) == length(retained)
    } else {
      rep(FALSE, n_primary)
    }
    n_method <- sum(included_method)
    n_unit <- sum(included_unit)
    degenerate <- n_method < min_methods || n_unit < min_units

    summaries[[i]] <- data.frame(
      threshold = threshold,
      required_eligible_units = required_units,
      total_primary_units = n_primary,
      n_methods = n_method,
      n_units = n_unit,
      min_methods = min_methods,
      min_units = min_units,
      degenerate = degenerate,
      rank_allowed = !degenerate,
      stringsAsFactors = FALSE
    )
    methods[[i]] <- data.frame(
      threshold = threshold,
      method_id = within_ids,
      roster_order = match(within_ids, roster_ids),
      eligible_units = as.integer(coverage_n),
      total_primary_units = n_primary,
      coverage = as.numeric(coverage_n / n_primary),
      included_method = included_method,
      exclusion_reason = ifelse(included_method, NA_character_,
                                "below_coverage_threshold"),
      stringsAsFactors = FALSE
    )
    unit_reason <- if (!length(retained)) {
      rep("no_method_meets_threshold", n_primary)
    } else {
      ifelse(included_unit, NA_character_, "ineligible_for_retained_method")
    }
    units[[i]] <- data.frame(
      threshold = threshold,
      unit_id = primary_units,
      included_unit = included_unit,
      exclusion_reason = unit_reason,
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = do.call(rbind, summaries),
    method_membership = do.call(rbind, methods),
    unit_membership = do.call(rbind, units)
  )
}

#' Adjust the frozen zero-difference test family with BY
#'
#' Applies Benjamini--Yekutieli adjustment using the complete frozen pair family
#' as the denominator. Exactly one row per frozen pair is required; an
#' unavailable p-value must be represented explicitly as `NA` and is retained
#' in the output.
#'
#' @param pair_id Pair identifiers associated with `p_zero`.
#' @param p_zero Raw two-sided zero-difference p-values.
#' @param frozen_pair_ids Complete, ordered, prespecified pooled-testable pair
#'   family.
#'
#' @return A data frame in frozen family order with raw and BY-adjusted p-values
#'   and the frozen family size.
#' @export
singlesample_adjust_zero_by <- function(pair_id, p_zero, frozen_pair_ids) {
  ordered <- .singlesample_frozen_p_values(
    pair_id, p_zero, frozen_pair_ids, "p_zero"
  )
  data.frame(
    pair_id = ordered$id,
    p_zero = ordered$p,
    p_zero_by = stats::p.adjust(
      ordered$p, method = "BY", n = length(ordered$id)
    ),
    family_n = length(ordered$id),
    stringsAsFactors = FALSE
  )
}

#' Adjust the frozen directional relevance family with BY
#'
#' Jointly adjusts the two prespecified directional shifted-null tests for each
#' frozen pair. The denominator is exactly twice the number of frozen pairs;
#' directions cannot be selected after observing their p-values.
#'
#' @param pair_id Pair identifiers associated with the p-values.
#' @param p_above_margin Raw p-values for `H0: effect <= +margin`.
#' @param p_below_minus_margin Raw p-values for `H0: effect >= -margin`.
#' @param frozen_pair_ids Complete, ordered, prespecified pooled-testable pair
#'   family.
#'
#' @return A data frame in frozen pair order containing both raw and BY-adjusted
#'   directional p-values and the frozen directional-family size.
#' @export
singlesample_adjust_relevance_by <- function(
    pair_id, p_above_margin, p_below_minus_margin, frozen_pair_ids) {
  above <- .singlesample_frozen_p_values(
    pair_id, p_above_margin, frozen_pair_ids, "p_above_margin"
  )
  below <- .singlesample_frozen_p_values(
    pair_id, p_below_minus_margin, frozen_pair_ids, "p_below_minus_margin"
  )
  raw <- c(rbind(above$p, below$p))
  family_n <- 2L * length(above$id)
  adjusted <- stats::p.adjust(raw, method = "BY", n = family_n)
  adjusted <- matrix(adjusted, nrow = 2L)
  data.frame(
    pair_id = above$id,
    p_above_margin = above$p,
    p_below_minus_margin = below$p,
    p_above_margin_by = adjusted[1L, ],
    p_below_minus_margin_by = adjusted[2L, ],
    family_n = family_n,
    stringsAsFactors = FALSE
  )
}

.singlesample_require_columns <- function(x, columns, label) {
  missing_columns <- setdiff(columns, names(x))
  if (length(missing_columns)) {
    stop(label, " is missing required columns: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.singlesample_scalar_integer <- function(x, label, minimum) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x != as.integer(x) ||
      x < minimum) {
    stop(label, " must be one integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.singlesample_scalar_number <- function(x, label, lower = -Inf, upper = Inf,
                                         strict = FALSE) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(label, " must be one finite number.", call. = FALSE)
  }
  inside <- if (strict) x > lower && x < upper else x >= lower && x <= upper
  if (!inside) {
    brackets <- if (strict) "(" else "["
    close <- if (strict) ")" else "]"
    stop(label, " must lie in ", brackets, lower, ", ", upper, close, ".",
         call. = FALSE)
  }
  as.numeric(x)
}

.singlesample_stratum_ids <- function(x, label) {
  if (is.null(x)) stop(label, " is required.", call. = FALSE)
  x <- as.character(x)
  if (anyNA(x) || any(!nzchar(x))) {
    stop(label, " must contain non-missing, non-empty identifiers.",
         call. = FALSE)
  }
  x
}

.singlesample_frozen_p_values <- function(id, p, frozen_ids, p_label) {
  id <- as.character(id)
  frozen_ids <- as.character(frozen_ids)
  if (!length(frozen_ids) || anyNA(frozen_ids) || any(!nzchar(frozen_ids)) ||
      anyDuplicated(frozen_ids)) {
    stop("frozen_pair_ids must be non-empty, non-missing, and unique.",
         call. = FALSE)
  }
  if (length(id) != length(p) || anyNA(id) || any(!nzchar(id)) ||
      anyDuplicated(id)) {
    stop("pair_id must uniquely identify every supplied p-value.",
         call. = FALSE)
  }
  if (!setequal(id, frozen_ids) || length(id) != length(frozen_ids)) {
    stop("pair_id must match the complete frozen_pair_ids family exactly.",
         call. = FALSE)
  }
  p <- as.numeric(p)
  if (any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1))) {
    stop(p_label, " must contain p-values in [0, 1] or explicit NA.",
         call. = FALSE)
  }
  index <- match(frozen_ids, id)
  list(id = frozen_ids, p = p[index])
}
