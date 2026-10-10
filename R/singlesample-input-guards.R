# Frozen deployment input checks. These never estimate a threshold, support, or
# variance from the scored batch. See score_specimen() and NEWS for the codes.
.ss_input_reject <- function(code, detail, caller, rows = integer()) {
  suffix <- if (length(rows)) paste0(" (row(s): ", paste(rows, collapse = ", "), ")") else ""
  stop(structure(
    list(message = paste0(caller, ": [", code, "] ", detail, suffix),
         call = NULL, code = code, rows = as.integer(rows)),
    class = c("singlesample_input_error", "error", "condition")
  ))
}

.ss_required_features <- function(model, fallback = character()) {
  if (inherits(model, "singlesample_legacy_model")) return(model$feature_universe)
  if (inherits(model, "reo_ktsp_model")) {
    return(unique(c(model$ktsp$pairs$feature_a, model$ktsp$pairs$feature_b)))
  }
  if (inherits(model, "reo_pairratio_model")) {
    return(unique(c(model$selected_pairs$feature_a, model$selected_pairs$feature_b)))
  }
  if (inherits(model, "ws_balance_ilr_model")) {
    return(intersect(model$sbp_canonical_features, model$feature_universe))
  }
  if (is.list(model)) {
    if (!is.null(model$panel)) return(model$panel)
    if (!is.null(model$feature_universe)) return(model$feature_universe)
  }
  fallback
}

.ss_guard_input <- function(X, required, caller, detected = NULL,
                            nondetect_limit = 0) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    .ss_input_reject("INVALID_INPUT", "requires a numeric specimens-by-features matrix", caller)
  }
  if (!nrow(X) || !ncol(X)) {
    .ss_input_reject("EMPTY_INPUT", "requires at least one row and one feature", caller)
  }
  features <- colnames(X)
  if (is.null(features) || anyNA(features) || any(!nzchar(features)) ||
      anyDuplicated(features)) {
    .ss_input_reject("INVALID_FEATURE_NAMES", "requires unique, non-empty feature names", caller)
  }
  present <- intersect(required, features)
  if (length(required) && !length(present)) {
    .ss_input_reject("ALL_REQUIRED_ABSENT", "all frozen required features are absent", caller)
  }
  # An unregistered positive column must not rescue an empty required panel.
  if (!length(required)) present <- features
  .ss_guard_finite_required(X, present, caller)
  signal <- is.finite(X[, present, drop = FALSE]) &
    X[, present, drop = FALSE] > nondetect_limit
  if (!is.null(detected)) {
    detected <- .ss_detection_mask(detected, X, caller)
    signal <- signal & detected[, present, drop = FALSE]
  }
  empty <- which(rowSums(signal) == 0L)
  if (length(empty)) {
    .ss_input_reject("ALL_NONDETECT",
                     "no finite detected required feature above the frozen nondetect limit",
                     caller, empty)
  }
  invisible(TRUE)
}

.ss_guard_finite_required <- function(X, required, caller) {
  present <- intersect(required, colnames(X))
  bad <- which(rowSums(!is.finite(X[, present, drop = FALSE])) > 0L)
  if (length(bad)) {
    .ss_input_reject("NONFINITE_REQUIRED_VALUE",
                     "required features must contain finite numeric values", caller, bad)
  }
  invisible(TRUE)
}

.ss_detection_mask <- function(detected, X, caller) {
  features <- colnames(X)
  if (is.logical(detected) && is.null(dim(detected)) && nrow(X) == 1L) {
    detected <- matrix(detected, nrow = 1L, dimnames = list(NULL, names(detected)))
  }
  if (is.data.frame(detected)) detected <- as.matrix(detected)
  if (!is.matrix(detected) || !is.logical(detected) || anyNA(detected) ||
      !identical(dim(detected), dim(X)) ||
      !identical(colnames(detected), features) ||
      (!is.null(rownames(detected)) && !is.null(rownames(X)) &&
       !identical(rownames(detected), rownames(X)))) {
    .ss_input_reject("INVALID_DETECTION_MASK",
                     "detected must be a non-missing logical mask aligned to x_new, including feature names",
                     caller)
  }
  detected
}

.ss_support_status <- function(score) {
  if (!length(score) || !any(is.finite(score))) return("empty")
  if (any(!is.finite(score))) return("empty")
  if (all(score == score[[1L]])) "constant" else "ok"
}

.ss_guard_support <- function(model, caller) {
  status <- model$support_status
  if (identical(status, "constant")) {
    .ss_input_reject("CONSTANT_FITTED_SUPPORT",
                     "requires variable fitted pair/balance support; training support was constant", caller)
  }
  empty <- identical(status, "empty") ||
    (inherits(model, "reo_ktsp_model") && !nrow(model$ktsp$pairs)) ||
    (inherits(model, "reo_pairratio_model") && !nrow(model$selected_pairs))
  if (empty) {
    .ss_input_reject("EMPTY_FITTED_SUPPORT", "requires non-empty fitted pair/balance support", caller)
  }
  invisible(TRUE)
}

.ss_check_strict <- function(strict, caller) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    .ss_input_reject("INVALID_STRICT", "strict must be TRUE or FALSE", caller)
  }
  invisible(TRUE)
}

.ss_strict_scorer <- function(method_id, caller) {
  scorers <- c("reo-ktsp" = "score_reo_ktsp",
               "reo-pairratio" = "score_reo_pairratio",
               "ws-rclr-trimmed" = "score_ws_rclr_panel",
               "ws-balance-ilr" = "score_ws_balance_ilr")
  if (length(method_id) != 1L || is.na(method_id) ||
      !method_id %in% names(scorers)) {
    .ss_input_reject("STRICT_UNSUPPORTED",
                     "strict checking is unavailable for this method", caller)
  }
  unname(scorers[[method_id]])
}

.ss_guard_required <- function(required, features, caller) {
  missing <- setdiff(required, features)
  if (length(missing)) {
    .ss_input_reject("MISSING_REQUIRED_FEATURES",
                     paste0("missing frozen required feature(s): ",
                            paste(missing, collapse = ", ")), caller)
  }
  invisible(TRUE)
}

# Counts describe named-column availability, not detection or assay quality.
# Keep the numeric score API; one support row corresponds to each score.
.ss_attach_support <- function(score, n_available, n_fitted) {
  attr(score, "support") <- data.frame(
    n_available = rep.int(as.integer(n_available), length(score)),
    n_fitted = rep.int(as.integer(n_fitted), length(score)),
    full_support = rep.int(n_available == n_fitted, length(score))
  )
  score
}

.ss_guard_pair_input <- function(model, X, pairs, caller, strict = FALSE) {
  .ss_check_strict(strict, caller)
  .ss_guard_support(model, caller)
  .ss_guard_input(X, .ss_required_features(model), caller)
  if (strict) .ss_guard_required(.ss_required_features(model), colnames(X), caller)
  usable <- pairs$feature_a %in% colnames(X) & pairs$feature_b %in% colnames(X)
  if (!any(usable)) {
    .ss_input_reject("NO_USABLE_SUPPORT", "no fitted pair has both features present", caller)
  }
  invisible(TRUE)
}
