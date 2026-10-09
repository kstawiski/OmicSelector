#' Freeze specimen-local preprocessing for counts or Ct measurements
#'
#' @description
#' Records an ordered feature universe and conversion constants. No training
#' means, totals, quantiles or other estimated statistics are retained.
#'
#' @param X_raw Raw specimens-by-features matrix or data frame, with unique,
#'   non-empty column names. Counts must be numeric, finite and non-negative.
#'   Ct input may include character markers, which are converted to NA.
#' @param type Measurement type: \code{"counts"} or \code{"ct"}.
#' @param ct_limit Positive finite Ct detection limit (default 40).
#' @param ct_base Finite exponential base greater than 1 (default 2).
#' @param scale For counts, \code{"cpm"} divides each row by its own total
#'   over available frozen features and multiplies by one million;
#'   \code{"none"} leaves counts unchanged. Ignored for Ct conversion.
#' @param features Ordered character vector of unique feature names present
#'   in \code{X_raw}. Defaults to all training columns.
#'
#' @details
#' Counts of zero are non-detects. An all-zero row remains zero.
#' Ct values that are missing, non-finite, non-numeric, or at or above
#' \code{ct_limit} are non-detects. Their Ct is replaced by \code{ct_limit}
#' before computing \code{ct_base^(ct_limit - Ct)}, so they map to abundance 1.
#' Negative finite Ct values below the limit are converted by the same formula.
#'
#' Columns outside the frozen universe are dropped. Missing frozen columns
#' remain absent, without imputation. For CPM, only available frozen columns
#' contribute to that specimen's total; missing columns can therefore change
#' its conversion. Scorer input guards determine whether support is sufficient.
#'
#' Invalid inputs raise \code{singlesample_input_error}. Codes include
#' \code{INVALID_INPUT}, \code{EMPTY_INPUT}, \code{INVALID_FEATURE_NAMES},
#' \code{INVALID_PREPROCESSOR}, \code{INVALID_COUNT_VALUE} (negative counts),
#' \code{NONFINITE_REQUIRED_VALUE} (non-finite frozen count values), and
#' \code{NONFINITE_ABUNDANCE} (conversion or count-total overflow).
#' Errors involving measurement values report affected \code{rows}.
#'
#' @return An S3 \code{singlesample_preprocessor} holding only \code{type},
#'   \code{features}, \code{ct_limit}, \code{ct_base}, and \code{scale}.
#' @seealso \code{\link{apply_singlesample_preprocessor}},
#'   \code{\link{deploy_singlesample}}
#' @examples
#' counts <- rbind(c(a = 10, b = 30), c(a = 6, b = 2))
#' prep <- fit_singlesample_preprocessor(counts, type = "counts")
#' apply_singlesample_preprocessor(prep, counts[1, ])
#' @export
fit_singlesample_preprocessor <- function(X_raw, type = c("counts", "ct"),
                                          ct_limit = 40, ct_base = 2,
                                          scale = c("cpm", "none"),
                                          features = colnames(X_raw)) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  preprocessor <- structure(list(type = type, features = features,
                                 ct_limit = ct_limit, ct_base = ct_base,
                                 scale = scale),
                            class = "singlesample_preprocessor")
  caller <- "fit_singlesample_preprocessor"
  .ss_check_preprocessor(preprocessor, caller)
  X <- .ss_preprocessor_matrix(X_raw, type, caller)
  if (length(setdiff(features, colnames(X)))) {
    .ss_input_reject("INVALID_FEATURE_NAMES",
                     "features must be present in X_raw at fitting", caller)
  }
  apply_singlesample_preprocessor(preprocessor, X)
  preprocessor
}

#' Apply a frozen specimen-local preprocessor
#'
#' @param preprocessor A \code{singlesample_preprocessor} returned by
#'   \code{\link{fit_singlesample_preprocessor}}.
#' @param X_raw Raw specimens-by-features matrix/data frame or a single named
#'   vector. See \code{\link{fit_singlesample_preprocessor}} for conversion
#'   rules and error codes.
#' @return A list with numeric matrix \code{abundance} and logical matrix
#'   \code{detected}, with matching dimensions and names. Columns follow the
#'   frozen order, restricted to supplied features. If no frozen columns are
#'   supplied, both matrices have zero columns.
#' @export
apply_singlesample_preprocessor <- function(preprocessor, X_raw) {
  caller <- "apply_singlesample_preprocessor"
  .ss_check_preprocessor(preprocessor, caller)
  X <- .ss_preprocessor_matrix(X_raw, preprocessor$type, caller)
  present <- preprocessor$features[preprocessor$features %in% colnames(X)]
  X <- X[, present, drop = FALSE]
  if (identical(preprocessor$type, "counts")) {
    bad <- which(rowSums(!is.finite(X)) > 0L)
    if (length(bad)) {
      .ss_input_reject("NONFINITE_REQUIRED_VALUE",
                       "frozen count features must contain finite values", caller, bad)
    }
    bad <- which(rowSums(X < 0) > 0L)
    if (length(bad)) {
      .ss_input_reject("INVALID_COUNT_VALUE", "counts must be non-negative", caller, bad)
    }
    detected <- X > 0
    if (identical(preprocessor$scale, "cpm")) {
      totals <- rowSums(X)
      bad <- which(!is.finite(totals))
      if (length(bad)) {
        .ss_input_reject("NONFINITE_ABUNDANCE", "count total overflow", caller, bad)
      }
      for (i in seq_len(nrow(X))) {
        if (totals[i] > 0) X[i, ] <- (X[i, ] / totals[i]) * 1e6
      }
    }
  } else {
    detected <- is.finite(X) & X < preprocessor$ct_limit
    X[!detected] <- preprocessor$ct_limit
    X <- preprocessor$ct_base^(preprocessor$ct_limit - X)
    bad <- which(rowSums(!is.finite(X)) > 0L)
    if (length(bad)) {
      .ss_input_reject("NONFINITE_ABUNDANCE", "Ct conversion overflow", caller, bad)
    }
  }
  list(abundance = X, detected = detected)
}

#' Print a frozen single-sample preprocessor
#'
#' @param x A \code{singlesample_preprocessor} object.
#' @param ... Unused.
#' @return Invisibly, \code{x}.
#' @export
print.singlesample_preprocessor <- function(x, ...) {
  .ss_check_preprocessor(x, "print.singlesample_preprocessor")
  cat("OmicSelector frozen single-sample preprocessor\n")
  cat("  type: ", x$type, "\n", sep = "")
  cat("  features: ", length(x$features), "\n", sep = "")
  if (identical(x$type, "counts")) {
    cat("  scale: ", x$scale, "\n", sep = "")
  } else {
    cat("  Ct limit: ", x$ct_limit, "; base: ", x$ct_base, "\n", sep = "")
  }
  invisible(x)
}

.ss_check_preprocessor <- function(x, caller) {
  if (!inherits(x, "singlesample_preprocessor") || !is.list(x) ||
      !identical(sort(names(x)), sort(c("type", "features", "ct_limit", "ct_base", "scale"))) ||
      !is.character(x$type) || length(x$type) != 1L || is.na(x$type) ||
      !x$type %in% c("counts", "ct") ||
      !is.character(x$scale) || length(x$scale) != 1L || is.na(x$scale) ||
      !x$scale %in% c("cpm", "none") ||
      !is.numeric(x$ct_limit) || length(x$ct_limit) != 1L ||
      !is.finite(x$ct_limit) || x$ct_limit <= 0 ||
      !is.numeric(x$ct_base) || length(x$ct_base) != 1L ||
      !is.finite(x$ct_base) || x$ct_base <= 1) {
    .ss_input_reject("INVALID_PREPROCESSOR",
                     "requires a frozen counts or Ct preprocessor with valid constants", caller)
  }
  if (!is.character(x$features) || !length(x$features) || anyNA(x$features) ||
      any(!nzchar(x$features)) || anyDuplicated(x$features)) {
    .ss_input_reject("INVALID_FEATURE_NAMES",
                     "frozen features must be unique, non-empty names", caller)
  }
  invisible(TRUE)
}

.ss_preprocessor_matrix <- function(X, type, caller) {
  if (is.atomic(X) && is.null(dim(X))) {
    X <- matrix(X, nrow = 1L, dimnames = list(NULL, names(X)))
  }
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !(is.numeric(X) ||
      (identical(type, "ct") && (is.character(X) || is.logical(X))))) {
    .ss_input_reject("INVALID_INPUT",
                     "requires a named vector or specimens-by-features matrix of raw measurements",
                     caller)
  }
  if (!nrow(X) || !ncol(X)) {
    .ss_input_reject("EMPTY_INPUT", "requires at least one row and one feature", caller)
  }
  features <- colnames(X)
  if (is.null(features) || anyNA(features) || any(!nzchar(features)) ||
      anyDuplicated(features)) {
    .ss_input_reject("INVALID_FEATURE_NAMES", "requires unique, non-empty feature names", caller)
  }
  # Numeric Ct strings and non-numeric markers are handled before verification;
  # this preserves the raw Ct scale without carrying character data into probes.
  suppressWarnings(storage.mode(X) <- "double")
  X
}

.ss_preprocessor_extreme_transform <- function(preprocessor) {
  if (is.null(preprocessor) || !identical(preprocessor$type, "ct")) return(NULL)
  ct_limit <- preprocessor$ct_limit
  shift <- log(1e6) / log(preprocessor$ct_base)
  function(X) {
    detected <- is.finite(X) & X < ct_limit
    X[detected] <- X[detected] - shift
    X
  }
}
