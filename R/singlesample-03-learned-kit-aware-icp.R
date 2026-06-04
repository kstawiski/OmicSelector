# Inductive conformal predictors fitted separately per kit family.

if (!exists(".lka_as_matrix", mode = "function")) {
  stop("Learned kit-aware helper functions are unavailable.", call. = FALSE)
}

.lka_split_calibration <- function(y, calibration_fraction, seed) {
  y <- .lka_as_binary(y)
  set.seed(seed)
  cal <- logical(length(y))
  for (cl in sort(unique(y))) {
    idx <- which(y == cl)
    n_cal <- max(1L, floor(length(idx) * calibration_fraction))
    if (length(idx) - n_cal < 2L) n_cal <- max(1L, length(idx) - 2L)
    if (n_cal > 0L) cal[sample(idx, n_cal)] <- TRUE
  }
  cal
}

.lka_fit_conformal_predictor <- function(x, y, panel_features, calibration_fraction,
                                         seed, label) {
  y <- .lka_as_binary(y)
  if (length(unique(y)) != 2L || min(table(y)) < 4L) {
    return(list(status = "not_evaluable",
                status_reason = paste0(label, ": fewer than 4 samples per class")))
  }
  cal <- .lka_split_calibration(y, calibration_fraction, seed)
  if (sum(!cal & y == 1L) < 2L || sum(!cal & y == 0L) < 2L ||
      sum(cal & y == 1L) < 1L || sum(cal & y == 0L) < 1L) {
    return(list(status = "not_evaluable",
                status_reason = paste0(label, ": split conformal class counts too small")))
  }
  panel <- intersect(panel_features, colnames(x))
  if (length(panel) < 2L) {
    return(list(status = "not_evaluable",
                status_reason = paste0(label, ": panel overlap <2")))
  }
  w <- .lka_feature_weights(x[!cal, , drop = FALSE], y[!cal], colnames(x))
  raw_train <- .lka_score_projection(x[!cal, , drop = FALSE], w, panel)
  if (length(unique(raw_train[is.finite(raw_train)])) < 2L) {
    return(list(status = "not_evaluable",
                status_reason = paste0(label, ": constant training score")))
  }
  df <- data.frame(y = y[!cal], score = raw_train)
  glm_fit <- tryCatch(stats::glm(y ~ score, data = df, family = stats::binomial()),
                      error = function(e) NULL)
  if (is.null(glm_fit)) {
    return(list(status = "not_evaluable",
                status_reason = paste0(label, ": logistic calibration failed")))
  }
  raw_cal <- .lka_score_projection(x[cal, , drop = FALSE], w, panel)
  prob_cal <- as.numeric(stats::predict(glm_fit, newdata = data.frame(score = raw_cal),
                                        type = "response"))
  prob_cal <- pmin(pmax(prob_cal, 1e-6), 1 - 1e-6)
  y_cal <- y[cal]
  nonconf <- ifelse(y_cal == 1L, 1 - prob_cal, prob_cal)
  nonconf <- nonconf[is.finite(nonconf)]
  alpha <- 0.10
  qhat <- if (length(nonconf) > 0L) {
    stats::quantile(nonconf, probs = min(1, ceiling((length(nonconf) + 1) * (1 - alpha)) /
                                          length(nonconf)),
                    names = FALSE, type = 1, na.rm = TRUE)
  } else NA_real_
  list(status = "evaluable",
       status_reason = "",
       kit_label = label,
       panel_features = panel,
       weights = w,
       glm = glm_fit,
       calibration_nonconformity = nonconf,
       qhat_90 = qhat,
       n_proper = sum(!cal),
       n_calibration = sum(cal),
       n_calibration_pos = sum(y_cal == 1L),
       n_calibration_neg = sum(y_cal == 0L))
}

fit_icp_per_kit <- function(train_expr, train_meta = NULL,
                            train_disease_labels,
                            train_kit_labels,
                            train_cohort_labels = NULL,
                            panel_features = NULL,
                            calibration_fraction = 0.25,
                            min_per_kit = 16L,
                            seed = 42L,
                            ...) {
  x <- .lka_as_matrix(train_expr)
  y <- .lka_as_binary(train_disease_labels)
  if (length(y) != nrow(x)) stop("train_disease_labels length does not match train_expr rows.")
  if (is.null(panel_features)) panel_features <- colnames(x)
  panel_features <- intersect(panel_features, colnames(x))
  if (length(panel_features) == 0L) stop("panel_features do not overlap train_expr.")
  kit <- .lka_factor_state(train_kit_labels)
  scaler <- .lka_standardize_fit(x)
  z <- .lka_standardize_apply(x, scaler)

  predictors <- list()
  audit <- list()
  for (lvl in kit$levels) {
    idx <- which(kit$labels == lvl)
    if (length(idx) < min_per_kit || length(unique(y[idx])) != 2L ||
        min(table(y[idx])) < 4L) {
      audit[[lvl]] <- data.frame(
        kit_label = lvl,
        status = "not_evaluable",
        status_reason = sprintf("kit stratum too small: n=%d pos=%d neg=%d",
                                length(idx), sum(y[idx] == 1L), sum(y[idx] == 0L)),
        stringsAsFactors = FALSE)
      next
    }
    pred <- .lka_fit_conformal_predictor(z[idx, , drop = FALSE], y[idx],
                                         panel_features, calibration_fraction,
                                         seed + match(lvl, kit$levels), lvl)
    audit[[lvl]] <- data.frame(
      kit_label = lvl,
      status = pred$status,
      status_reason = pred$status_reason,
      n_proper = pred$n_proper %||% NA_integer_,
      n_calibration = pred$n_calibration %||% NA_integer_,
      qhat_90 = pred$qhat_90 %||% NA_real_,
      stringsAsFactors = FALSE)
    if (identical(pred$status, "evaluable")) predictors[[lvl]] <- pred
  }

  global <- .lka_fit_conformal_predictor(z, y, panel_features,
                                         calibration_fraction, seed + 999L, "global")
  if (identical(global$status, "evaluable")) predictors[["global"]] <- global
  if (length(predictors) == 0L) stop("No evaluable kit-specific or global conformal predictors.")

  list(method = "icp_per_kit",
       framework = "base_R_split_conformal",
       framework_detail = "Hand-rolled split conformal predictors per kit family; test samples use the kit predictor with smallest confidence nonconformity, falling back to global when kit strata are too small",
       scaler = scaler,
       predictors = predictors,
       kit_levels = kit$levels,
       audit = do.call(rbind, audit),
       panel_features = panel_features,
       calibration_fraction = calibration_fraction,
       min_per_kit = min_per_kit)
}

score_icp_per_kit <- function(model, test_expr, panel_features = NULL) {
  z <- .lka_standardize_apply(test_expr, model$scaler)
  panel_features <- panel_features %||% model$panel_features
  scores <- matrix(NA_real_, nrow = nrow(z), ncol = length(model$predictors))
  nonconf <- matrix(NA_real_, nrow = nrow(z), ncol = length(model$predictors))
  colnames(scores) <- names(model$predictors)
  colnames(nonconf) <- names(model$predictors)
  for (nm in names(model$predictors)) {
    pred <- model$predictors[[nm]]
    panel <- intersect(panel_features, names(pred$weights))
    raw <- .lka_score_projection(z, pred$weights, panel)
    prob <- as.numeric(stats::predict(pred$glm, newdata = data.frame(score = raw),
                                      type = "response"))
    prob <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
    scores[, nm] <- stats::qlogis(prob)
    nonconf[, nm] <- 1 - pmax(prob, 1 - prob)
  }
  chosen <- max.col(-nonconf, ties.method = "first")
  out <- scores[cbind(seq_len(nrow(scores)), chosen)]
  attr(out, "selected_predictor") <- colnames(scores)[chosen]
  attr(out, "selected_nonconformity") <- nonconf[cbind(seq_len(nrow(nonconf)), chosen)]
  out
}

predict_icp_per_kit_details <- function(model, test_expr, panel_features = NULL) {
  score <- score_icp_per_kit(model, test_expr, panel_features)
  data.frame(score = score,
             selected_predictor = attr(score, "selected_predictor"),
             selected_nonconformity = attr(score, "selected_nonconformity"),
             stringsAsFactors = FALSE)
}
