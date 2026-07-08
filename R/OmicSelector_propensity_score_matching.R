#' OmicSelector_propensity_score_matching
#'
#' Lightweight propensity score matching (nearest neighbor) using base R.
#'
#' @param dataset Original dataset.
#' @param match_by Character vector of covariate names to match on.
#' @param class_col Outcome column name (default: "Class").
#' @param positive Label for the treated/positive class (default: "Case").
#' @param method Matching method (currently only "nearest").
#' @param distance Propensity model type (currently only "logit").
#' @param ratio Number of controls per treated (default: 1).
#' @param replace Logical; allow matching controls with replacement.
#' @param caliper Optional maximum absolute propensity score distance.
#' @param na_action How to handle missing values in match_by: "omit" or "median".
#'
#' @return Matched dataset (treated + matched controls).
#'
#' @export
OmicSelector_propensity_score_matching <- function(dataset,
                                                   match_by = c("age_at_diagnosis", "gender.x"),
                                                   class_col = "Class",
                                                   positive = "Case",
                                                   method = "nearest",
                                                   distance = "logit",
                                                   ratio = 1,
                                                   replace = FALSE,
                                                   caliper = NULL,
                                                   na_action = c("omit", "median")) {
  na_action <- match.arg(na_action)

  if (!method %in% c("nearest")) {
    stop("Only method = 'nearest' is supported.", call. = FALSE)
  }

  if (!distance %in% c("logit")) {
    stop("Only distance = 'logit' is supported.", call. = FALSE)
  }

  if (!all(match_by %in% names(dataset))) {
    missing_cols <- setdiff(match_by, names(dataset))
    stop(sprintf("Missing match_by columns: %s", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  if (!class_col %in% names(dataset)) {
    stop(sprintf("class_col '%s' not found in dataset.", class_col), call. = FALSE)
  }

  dat <- dataset[, c(match_by, class_col), drop = FALSE]
  dat$.row_id <- seq_len(nrow(dataset))

  if (na_action == "median") {
    for (col in match_by) {
      if (is.numeric(dat[[col]])) {
        med <- stats::median(dat[[col]], na.rm = TRUE)
        dat[[col]][is.na(dat[[col]])] <- med
      } else {
        tab <- sort(table(dat[[col]]), decreasing = TRUE)
        if (length(tab) > 0) {
          mode_val <- names(tab)[1]
          dat[[col]][is.na(dat[[col]])] <- mode_val
        }
      }
    }
  } else {
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  }

  treated <- dat[[class_col]] == positive
  if (!any(treated)) {
    stop("No treated (positive) samples found.", call. = FALSE)
  }
  if (all(treated)) {
    stop("No control samples found.", call. = FALSE)
  }

  glm_data <- dat
  glm_data[[class_col]] <- treated

  formula <- stats::as.formula(
    paste(class_col, "~", paste(match_by, collapse = " + "))
  )

  fit <- stats::glm(formula, data = glm_data, family = stats::binomial())
  ps <- stats::predict(fit, type = "response")

  treated_idx <- which(treated)
  control_idx <- which(!treated)

  matched_control_idx <- integer(0)
  matched_treated_idx <- integer(0)

  available_controls <- control_idx

  for (ti in treated_idx) {
    if (length(available_controls) == 0 && !replace) {
      break
    }

    pool <- if (replace) control_idx else available_controls
    distances <- abs(ps[pool] - ps[ti])

    if (!is.null(caliper)) {
      keep <- distances <= caliper
      pool <- pool[keep]
      distances <- distances[keep]
    }

    if (length(pool) == 0) {
      next
    }

    ord <- order(distances)
    pick <- pool[ord][seq_len(min(ratio, length(ord)))]

    matched_treated_idx <- c(matched_treated_idx, rep(ti, length(pick)))
    matched_control_idx <- c(matched_control_idx, pick)

    if (!replace) {
      available_controls <- setdiff(available_controls, pick)
    }
  }

  if (length(matched_control_idx) == 0) {
    stop("No matches found. Consider increasing caliper or using replace = TRUE.",
         call. = FALSE)
  }

  matched_rows <- c(matched_treated_idx, matched_control_idx)
  dataset[dat$.row_id[matched_rows], , drop = FALSE]
}
