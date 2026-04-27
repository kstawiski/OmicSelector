#' @title Summarize Within-Provenance Case-Control Blocks
#'
#' @description
#' Counts cases and controls inside provenance blocks and marks whether any block
#' is large enough, and documented enough, to support a positive disease-contrast
#' claim. Missing direct-control evidence is treated fail-closed.
#'
#' @param data Data frame containing outcome and provenance columns.
#' @param outcome Column name for the binary outcome.
#' @param provenance_block Column name defining provenance blocks.
#' @param direct_control_evidence Optional logical or character vector, or column
#'   name, indicating direct physical-provenance evidence for control rows.
#' @param positive Optional positive-class label. Defaults to the second factor level.
#' @param min_cases Minimum cases required inside a claim-supporting block.
#' @param min_controls Minimum controls required inside a claim-supporting block.
#'
#' @return A data frame with per-block counts and claim-support flags.
#' @export
os_within_provenance_blocks <- function(data, outcome, provenance_block,
                                        direct_control_evidence = NULL,
                                        positive = NULL,
                                        min_cases = 30L,
                                        min_controls = 30L) {
  data <- as.data.frame(data)
  if (!outcome %in% names(data)) stop("`outcome` column not found.", call. = FALSE)
  if (!provenance_block %in% names(data)) stop("`provenance_block` column not found.", call. = FALSE)
  y <- data[[outcome]]
  y_num <- .os_binary_with_positive(y, positive)
  block <- as.character(data[[provenance_block]])
  evidence <- .os_direct_evidence_vector(data, direct_control_evidence, length(y_num))

  rows <- lapply(sort(unique(block)), function(b) {
    idx <- which(block == b)
    ctrl <- idx[y_num[idx] == 0L]
    direct_pass <- length(ctrl) > 0L && all(.os_evidence_true(evidence[ctrl]))
    n_cases <- sum(y_num[idx] == 1L)
    n_controls <- sum(y_num[idx] == 0L)
    data.frame(
      provenance_block = b,
      n_total = length(idx),
      n_cases = n_cases,
      n_controls = n_controls,
      direct_control_evidence_pass = direct_pass,
      supports_positive_claim = n_cases >= min_cases &&
        n_controls >= min_controls &&
        isTRUE(direct_pass),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @title Estimate a Provenance-Only Prediction Floor
#'
#' @description
#' Fits simple train-fold-only logistic provenance models to estimate whether
#' non-expression source variables alone can predict the outcome. This is a
#' lightweight package-level governance helper, not a replacement for a full
#' prespecified modeling script.
#'
#' @param data Data frame containing outcome and provenance predictors.
#' @param outcome Column name for the binary outcome.
#' @param provenance_predictors Character vector of non-expression predictor columns.
#' @param group_id Optional column name or vector defining grouped folds.
#' @param positive Optional positive-class label.
#' @param n_folds Number of grouped folds.
#' @param seed Random seed.
#'
#' @return An \code{os_provenance_floor} object.
#' @export
os_provenance_floor_suite <- function(data, outcome, provenance_predictors,
                                      group_id = NULL, positive = NULL,
                                      n_folds = 5L, seed = 1L) {
  data <- as.data.frame(data)
  missing <- setdiff(c(outcome, provenance_predictors), names(data))
  if (length(missing) > 0L) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  y_num <- .os_binary_with_positive(data[[outcome]], positive)
  gid <- .os_group_vector(data, group_id, length(y_num))
  folds <- os_make_grouped_stratified_folds(y_num, group_id = gid, n_folds = n_folds, seed = seed)

  model_df <- data.frame(y = y_num, data[, provenance_predictors, drop = FALSE], check.names = FALSE)
  for (nm in provenance_predictors) {
    if (is.character(model_df[[nm]])) model_df[[nm]] <- as.factor(model_df[[nm]])
  }
  all_score <- .os_oof_glm_score(model_df, folds, stats::as.formula("y ~ ."))
  single <- lapply(provenance_predictors, function(nm) {
    f <- stats::as.formula(paste("y ~", sprintf("`%s`", nm)))
    score <- .os_oof_glm_score(model_df[, c("y", nm), drop = FALSE], folds, f)
    data.frame(
      predictor = nm,
      auc = .auc_fast(y_num, score),
      stringsAsFactors = FALSE
    )
  })
  single <- do.call(rbind, single)
  out <- list(
    primary_auc = .auc_fast(y_num, all_score),
    primary_score = all_score,
    single_covariate = single,
    single_covariate_max_auc = max(single$auc, na.rm = TRUE),
    folds = folds,
    predictors = provenance_predictors
  )
  class(out) <- c("os_provenance_floor", "list")
  out
}

#' @title Apply a Fail-Closed Identifiability Gate
#'
#' @description
#' Adjudicates whether a public-omics contrast is identifiable enough to permit
#' downstream positive disease-signal language. Any missing or failing criterion
#' yields \code{"FAIL-DEMOTED"}.
#'
#' @param provenance_floor Optional \code{os_provenance_floor} object.
#' @param provenance_auc Provenance-only AUC.
#' @param provenance_auc_ci_hi Upper confidence bound for provenance-only AUC.
#' @param single_covariate_max_auc Maximum single-provenance-covariate AUC.
#' @param sensitivity_auc Sensitivity-model AUC, for example a random-forest floor.
#' @param within_blocks Data frame from \code{os_within_provenance_blocks}.
#' @param u1_u2_overlap_fraction Fractional overlap of disease biological IDs.
#' @param u1_only_cases Count of disease biological IDs unique to U1.
#' @param u1_only_controls Count of control biological IDs unique to U1.
#' @param physical_qc_pass Logical physical-provenance and QC comparability flag.
#' @param anchor_leakage_pass Logical external-anchor leakage flag.
#' @param floor_auc_ceiling Maximum allowed provenance-floor AUC.
#' @param floor_ci_ceiling Maximum allowed upper CI for provenance-floor AUC.
#' @param single_covariate_ceiling Maximum allowed single-covariate AUC.
#' @param sensitivity_auc_ceiling Maximum allowed sensitivity-model AUC.
#' @param max_overlap_fraction Maximum allowed U1/U2 disease-ID overlap.
#' @param min_u1_only_cases Minimum U1-only cases.
#' @param min_u1_only_controls Minimum U1-only controls.
#'
#' @return An \code{os_identifiability_gate} object.
#' @export
os_identifiability_gate <- function(provenance_floor = NULL,
                                    provenance_auc = NULL,
                                    provenance_auc_ci_hi = NULL,
                                    single_covariate_max_auc = NULL,
                                    sensitivity_auc = NULL,
                                    within_blocks = NULL,
                                    u1_u2_overlap_fraction = NULL,
                                    u1_only_cases = NULL,
                                    u1_only_controls = NULL,
                                    physical_qc_pass = FALSE,
                                    anchor_leakage_pass = FALSE,
                                    floor_auc_ceiling = 0.65,
                                    floor_ci_ceiling = 0.70,
                                    single_covariate_ceiling = 0.70,
                                    sensitivity_auc_ceiling = 0.70,
                                    max_overlap_fraction = 0.80,
                                    min_u1_only_cases = 30L,
                                    min_u1_only_controls = 30L) {
  if (!is.null(provenance_floor)) {
    provenance_auc <- provenance_auc %||% provenance_floor$primary_auc
    single_covariate_max_auc <- single_covariate_max_auc %||% provenance_floor$single_covariate_max_auc
  }

  provenance_vals <- c(provenance_auc, provenance_auc_ci_hi, single_covariate_max_auc, sensitivity_auc)
  c1 <- all(is.finite(provenance_vals)) &&
    provenance_auc < floor_auc_ceiling &&
    provenance_auc_ci_hi < floor_ci_ceiling &&
    single_covariate_max_auc < single_covariate_ceiling &&
    sensitivity_auc < sensitivity_auc_ceiling
  c2 <- !is.null(within_blocks) &&
    "supports_positive_claim" %in% names(within_blocks) &&
    any(isTRUE(within_blocks$supports_positive_claim))
  c3 <- is.finite(u1_u2_overlap_fraction) &&
    is.finite(u1_only_cases) &&
    is.finite(u1_only_controls) &&
    u1_u2_overlap_fraction <= max_overlap_fraction &&
    u1_only_cases >= min_u1_only_cases &&
    u1_only_controls >= min_u1_only_controls
  c4 <- isTRUE(physical_qc_pass)
  c5 <- isTRUE(anchor_leakage_pass)

  criteria <- data.frame(
    criterion = c(
      "Provenance floor below ceiling",
      "Adequate within-provenance contrast",
      "Limited universe overlap",
      "Physical provenance and QC comparable",
      "External anchor leakage controlled"
    ),
    status = ifelse(c(c1, c2, c3, c4, c5), "PASS", "FAIL"),
    value = c(
      sprintf("auc=%s;ci_hi=%s;single=%s;sensitivity=%s",
              .os_fmt_num(provenance_auc), .os_fmt_num(provenance_auc_ci_hi),
              .os_fmt_num(single_covariate_max_auc), .os_fmt_num(sensitivity_auc)),
      if (is.null(within_blocks)) "no block table" else paste0(sum(within_blocks$supports_positive_claim), " supporting block(s)"),
      sprintf("overlap=%s;u1_only_cases=%s;u1_only_controls=%s",
              .os_fmt_num(u1_u2_overlap_fraction), .os_fmt_num(u1_only_cases), .os_fmt_num(u1_only_controls)),
      as.character(isTRUE(physical_qc_pass)),
      as.character(isTRUE(anchor_leakage_pass))
    ),
    threshold = c(
      sprintf("auc<%.2f, ci_hi<%.2f, single<%.2f, sensitivity<%.2f",
              floor_auc_ceiling, floor_ci_ceiling, single_covariate_ceiling, sensitivity_auc_ceiling),
      ">=1 block with min cases, min controls, and direct control evidence",
      sprintf("overlap<=%.2f, u1-only cases>=%s, controls>=%s",
              max_overlap_fraction, min_u1_only_cases, min_u1_only_controls),
      "TRUE",
      "TRUE"
    ),
    stringsAsFactors = FALSE
  )
  terminal_state <- if (all(criteria$status == "PASS")) "PASS" else "FAIL-DEMOTED"
  out <- list(
    terminal_state = terminal_state,
    criteria = criteria,
    failed_criteria = criteria$criterion[criteria$status != "PASS"]
  )
  class(out) <- c("os_identifiability_gate", "list")
  out
}

#' @title Build a Terminal Gate Ledger
#'
#' @param gate_id Character vector of gate identifiers.
#' @param status Character vector of terminal statuses.
#' @param note Optional character vector of notes.
#'
#' @return A data frame with fail-closed claim permissions.
#' @export
os_terminal_gate_ledger <- function(gate_id, status, note = NULL) {
  if (length(gate_id) != length(status)) stop("`gate_id` and `status` must have the same length.", call. = FALSE)
  if (is.null(note)) note <- rep("", length(gate_id))
  if (length(note) != length(gate_id)) stop("`note` must be NULL or the same length as `gate_id`.", call. = FALSE)
  out <- data.frame(
    gate_id = as.character(gate_id),
    status = as.character(status),
    claim_permission = ifelse(status == "PASS", "eligible_if_downstream_gates_pass", "blocked_descriptive_only"),
    note = as.character(note),
    stringsAsFactors = FALSE
  )
  class(out) <- c("os_terminal_gate_ledger", "data.frame")
  out
}

.os_binary_with_positive <- function(y, positive = NULL) {
  if (is.null(positive)) return(.os_binary_y(y))
  y_chr <- as.character(y)
  if (!positive %in% y_chr) stop("`positive` is not present in `y`.", call. = FALSE)
  out <- as.integer(y_chr == as.character(positive))
  if (length(unique(out)) != 2L) stop("`y` must contain both positive and negative classes.", call. = FALSE)
  out
}

.os_direct_evidence_vector <- function(data, evidence, n) {
  if (is.null(evidence)) return(rep(NA, n))
  if (length(evidence) == 1L && is.character(evidence) && evidence %in% names(data)) return(data[[evidence]])
  if (length(evidence) != n) stop("`direct_control_evidence` must be a column name or vector with one value per row.", call. = FALSE)
  evidence
}

.os_evidence_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "pass", "direct", "direct_pass", "1")
}

.os_group_vector <- function(data, group_id, n) {
  if (is.null(group_id)) return(seq_len(n))
  if (length(group_id) == 1L && is.character(group_id) && group_id %in% names(data)) return(data[[group_id]])
  if (length(group_id) != n) stop("`group_id` must be a column name or vector with one value per row.", call. = FALSE)
  group_id
}

.os_oof_glm_score <- function(model_df, folds, formula) {
  score <- rep(NA_real_, nrow(model_df))
  for (idx in folds) {
    train <- setdiff(seq_len(nrow(model_df)), idx)
    pred <- tryCatch({
      if (length(unique(model_df$y[train])) < 2L) stop("single-class train fold", call. = FALSE)
      fit <- suppressWarnings(stats::glm(formula, data = model_df[train, , drop = FALSE], family = stats::binomial()))
      as.numeric(stats::predict(fit, newdata = model_df[idx, , drop = FALSE], type = "response"))
    }, error = function(e) {
      rep(mean(model_df$y[train]), length(idx))
    })
    score[idx] <- pred
  }
  score
}

.os_fmt_num <- function(x) {
  if (is.null(x) || length(x) == 0L || !is.finite(x)) return("NA")
  sprintf("%.4f", as.numeric(x))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
