# Common deployment adapters for the legacy transform and metadata interfaces.
# These adapters never call a score wrapper without its training-fitted object.
.ss_legacy_ids <- function() c(
  "ws-alr-pivot", "ws-mad-logratio", "ws-dominance", "frozen-ruv",
  "robust-pca-resid", "hemolysis-resid", "kit-stratified-rclr",
  "kit-fe-adjusted-alr", "kit-orthogonal-ilr", "kit-residual-mad",
  "kit-stable-anchor-rclr", "hemolysis-kit-dual-anchor", "rin-weighted-score",
  "tech-stratified-rclr", "cross-tech-harmonized-rclr", "tech-residualized-alr",
  "biofluid-stratified-rclr", "biofluid-anchor-rclr", "biofluid-residualized-alr",
  "blood-cell-marker-excl-rclr", "nc-mahalanobis", "nc-conformal-anom", "nc-isoforest")

.ss_legacy_meta <- function(method, meta, n, schema = NULL, training = FALSE) {
  caller <- if (training) "deploy_singlesample" else "score_specimen"
  required <- if (!is.null(schema)) names(schema) else if (startsWith(method, "kit-") ||
    method == "hemolysis-kit-dual-anchor") "kit" else if (grepl("tech-", method))
      "technology" else if (startsWith(method, "biofluid-")) "biofluid" else if (
        method == "rin-weighted-score") "rin" else character()
  if (is.null(meta)) meta <- data.frame(row.names = seq_len(n))
  meta <- .singlesample_deploy_as_meta(meta, n, caller, if (training) "meta_train" else "meta")
  # Explicit aliases for the existing legacy metadata conventions.
  aliases <- list(kit = c("kit", "kit_label", "library_kit", "extraction_kit"),
                  rin = c("rin", "RIN", "rna_integrity", "rqi"))
  for (nm in required) {
    if (!nm %in% names(meta)) {
      hit <- intersect(aliases[[nm]], names(meta))
      if (length(hit)) meta[[nm]] <- meta[[hit[[1L]]]]
    }
    if (!nm %in% names(meta)) .ss_input_reject("MISSING_METADATA",
      paste0("method '", method, "' requires metadata column '", nm, "'"), caller)
    value <- trimws(as.character(meta[[nm]]))
    missing <- is.na(value) | !nzchar(value) | tolower(value) %in% c("unknown", "na", "n/a")
    if (any(missing)) .ss_input_reject("MISSING_METADATA",
      paste0("method '", method, "' requires non-missing '", nm, "'"), caller,
      rows = which(missing))
    if (nm == "rin") {
      value <- suppressWarnings(as.numeric(value))
      bad <- !is.finite(value) | value <= 0 | value > 10
      if (any(bad)) .ss_input_reject("INVALID_METADATA", "RIN must be numeric in (0, 10]",
                                   caller, rows = which(bad))
    } else {
      if (nm == "technology") value <- tryCatch(.ta_clean_technology(value),
        error = function(e) .ss_input_reject("UNSEEN_METADATA_LEVEL", conditionMessage(e), caller))
      if (nm == "biofluid") {
        value <- as.character(.bfa_canonical_biofluid(value))
        if (any(value == "unknown")) .ss_input_reject("UNSEEN_METADATA_LEVEL",
          "biofluid must identify a supported fluid", caller, rows = which(value == "unknown"))
      }
      if (!training && any(!value %in% schema[[nm]])) .ss_input_reject(
        "UNSEEN_METADATA_LEVEL", paste0("no frozen training level for '", nm, "'"),
        caller, rows = which(!value %in% schema[[nm]]))
    }
    meta[[nm]] <- value
  }
  if (training) {
    schema <- lapply(required, function(nm) if (nm == "rin") numeric() else sort(unique(meta[[nm]])))
    names(schema) <- required
  }
  list(meta = meta, schema = schema)
}

.ss_legacy_hemo <- function(X, meta, anchors = NULL, training = FALSE) {
  caller <- if (training) "deploy_singlesample" else "score_specimen"
  if (is.null(anchors)) anchors <- c(
    colnames(X)[which(colnames(X) %in% c("hsa-miR-451a", "MIMAT0001631"))[1L]],
    colnames(X)[which(colnames(X) %in% c("hsa-miR-23a-3p", "MIMAT0000078"))[1L]])
  if (!is.null(meta$hemo_score)) {
    h <- meta$hemo_score
    if (!is.numeric(h) || length(h) != nrow(X) || any(!is.finite(h)))
      .ss_input_reject("INVALID_METADATA", "hemo_score must contain one finite number per specimen", caller)
  } else {
    if (length(anchors) != 2L || anyNA(anchors) || any(!anchors %in% colnames(X)))
      .ss_input_reject("MISSING_HEMOLYSIS_ANCHORS",
        "requires hemo_score or the frozen miR-451a and miR-23a-3p anchors", caller)
    L <- log(X + 0.5)
    h <- as.numeric(L[, anchors[1L]] - L[, anchors[2L]])
  }
  list(score = h, anchors = anchors)
}

.ss_legacy_fit <- function(method, X, y, meta, fit_fn, group_train = NULL) {
  md <- .ss_legacy_meta(method, meta, nrow(X), training = TRUE)
  meta <- md$meta
  y <- .reo_check_labels(y, nrow(X), "deploy_singlesample")
  model <- list(method = method, feature_universe = colnames(X), metadata_schema = md$schema)
  # Panel and signs are learned once from training rCLR. Matrix transforms
  # below select their signed panel on the transformed training matrix instead.
  pan <- fit_ws_rclr_panel(X, y)
  roster <- singlesample_method_roster()
  model$score_fn_name <- roster$score_fn[match(method, roster$method_id)]
  model$panel <- pan$panel
  model$weights <- pan$weights
  fml <- if (is.function(fit_fn)) names(formals(fit_fn)) else character()
  args <- list()
  data_arg <- intersect(c("X_train", "x_train", "X_ref", "train_expr", "expr_matrix", "expr"), fml)
  if (length(data_arg)) args[[data_arg[1L]]] <- X
  for (nm in intersect(c("meta", "train_meta", "sample_meta"), fml)) args[[nm]] <- meta
  if ("kit_label_col" %in% fml) args$kit_label_col <- "kit"
  if ("kit_col" %in% fml) args$kit_col <- "kit"
  if ("rin_col" %in% fml) args$rin_col <- "rin"
  if ("panel" %in% fml) args$panel <- model$panel
  if ("train_meta" %in% fml) {
    args$train_meta$disease <- y
    if (!"cohort" %in% names(args$train_meta)) args$train_meta$cohort <- "training"
  }
  if (method == "hemolysis-resid") {
    h <- .ss_legacy_hemo(X, meta, training = TRUE)
    model$hemo_anchors <- h$anchors
    args$hemo_score <- h$score
  }
  if (method == "nc-conformal-anom" && !is.null(group_train)) args$group_train <- group_train
  # Negative controls fit their original reference estimators; no labels enter them.
  if (is.function(fit_fn)) model$fit <- do.call(fit_fn, args)
  if (method == "blood-cell-marker-excl-rclr") {
    Z <- .bfa_global_rclr_matrix(X[, setdiff(colnames(X), .bfa_default_blood_cell_mirnas()), drop = FALSE])
    model$weights <- .ka_select_panel_from_matrix(Z, y)
    model$panel <- names(model$weights)
  }
  transforms <- c("ws-alr-pivot", "ws-mad-logratio", "frozen-ruv", "robust-pca-resid", "hemolysis-resid")
  if (method %in% transforms) {
    Z <- .ss_legacy_transform(model, X, meta)
    model$weights <- .ss_deploy_signed_panel(Z, y)
    model$panel <- names(model$weights)
  }
  class(model) <- "singlesample_legacy_model"
  model
}

.ss_legacy_transform <- function(model, X, meta) {
  switch(model$method,
    `ws-alr-pivot` = ws_alr_pivot(X, allow_global_fallback = TRUE),
    `ws-mad-logratio` = ws_mad_logratio(X),
    `frozen-ruv` = apply_frozen_ruv(X, model$fit),
    `robust-pca-resid` = apply_robust_pca_residual(model$fit, X),
    `hemolysis-resid` = apply_hemolysis_rr(X, model$fit,
      .ss_legacy_hemo(X, meta, model$hemo_anchors)$score))
}

.ss_legacy_score <- function(model, X, meta) {
  meta <- .ss_legacy_meta(model$method, meta, nrow(X), model$metadata_schema)$meta
  if (model$method == "hemolysis-resid") .ss_legacy_hemo(X, meta, model$hemo_anchors)
  .ss_guard_required(model$feature_universe, colnames(X), "score_specimen")
  X <- X[, model$feature_universe, drop = FALSE]
  method <- model$method
  fit <- model$fit
  if (method %in% c("ws-alr-pivot", "ws-mad-logratio", "frozen-ruv", "robust-pca-resid", "hemolysis-resid")) {
    Z <- .ss_legacy_transform(model, X, meta)
    return(as.numeric(Z[, model$panel, drop = FALSE] %*% model$weights))
  }
  if (method == "ws-dominance") return(ws_dominance_score(X, top_k = fit$top_k, pseudocount = fit$pseudocount %||% 0.5))
  if (method == "nc-mahalanobis") return(apply_compositional_mahalanobis(fit, X)$distance)
  if (method == "nc-isoforest") return(apply_isolation_forest_logratio(fit, X)$score)
  if (method == "nc-conformal-anom") return(os_conformal_anomaly(X, fit)$score)
  fn <- .singlesample_deploy_resolve_fn(model$score_fn_name, "score_fn", method)
  fml <- names(formals(fn))
  args <- list()
  args[[intersect(c("expr_matrix", "expr"), fml)[1L]]] <- X
  for (nm in intersect(c("meta", "sample_meta"), fml)) args[[nm]] <- meta
  for (nm in intersect(c("panel", "panel_features"), fml)) args[[nm]] <- model$panel
  for (nm in intersect(c("fit", "anchor_fit", "dual_fit", "rin_fit"), fml)) args[[nm]] <- fit
  if ("feature_weights" %in% fml) args$feature_weights <- model$weights
  if ("kit_label_col" %in% fml) args$kit_label_col <- "kit"
  if ("kit_col" %in% fml) args$kit_col <- "kit"
  if ("rin_col" %in% fml) args$rin_col <- "rin"
  as.numeric(do.call(fn, args))
}

# Deployment copies of the benchmark aggregation contract. Keep the benchmark
# and its package primitives unchanged; learn these quantities only at fit time.
.ss_deploy_signed_panel <- function(Z, y, panel_size = 20L) {
  common <- colnames(Z)
  if (length(common) < 2L) stop("fewer than two common transformed features")
  effect <- vapply(common, function(j) {
    x <- Z[, j]; ok <- is.finite(x) & !is.na(y)
    den <- stats::sd(x[ok])
    if (sum(ok) < 5L || !is.finite(den) || den <= 0) return(0)
    abs(mean(x[ok & y == 1L]) - mean(x[ok & y == 0L])) / den
  }, numeric(1L))
  panel <- common[order(effect, decreasing = TRUE)[seq_len(min(panel_size, length(common)))]]
  vapply(panel, function(j) {
    d <- stats::median(Z[y == 1L, j], na.rm = TRUE) -
      stats::median(Z[y == 0L, j], na.rm = TRUE)
    if (is.finite(d) && d < 0) -1 else 1
  }, numeric(1L))
}

.ss_deploy_training_direction <- function(score, y) {
  y <- .reo_check_labels(y, length(score), "deploy_singlesample")
  ok <- is.finite(score) & !is.na(y)
  d <- if (sum(ok) >= 4L && length(unique(y[ok])) == 2L)
    stats::median(score[ok & y == 1L], na.rm = TRUE) -
      stats::median(score[ok & y == 0L], na.rm = TRUE) else NA_real_
  if (!is.finite(d)) stop("deploy_singlesample: training-only score direction unavailable", call. = FALSE)
  if (d < 0) -1 else 1
}
