#!/usr/bin/env Rscript
# Group-DRO classifier for kit x biofluid robust circulating-miRNA scoring.
#
# Runner utilities below require caller-supplied prepared inputs when used from
# the package; external analysis-workspace paths are intentionally not assumed.

if (!exists(".gdro_as_matrix", mode = "function")) {
  stop("Group-DRO helper functions are unavailable.", call. = FALSE)
}

#' Fit a Group-DRO kit x biofluid logistic scorer.
#'
#' @param expr_matrix Numeric training matrix, samples x miRNA features.
#' @param y Binary 0/1 training labels.
#' @param sample_meta Metadata aligned to rows of `expr_matrix`.
#' @param kit_col,biofluid_col Metadata columns used for train-only group
#'   construction. Groups are `kit x biofluid`.
#' @param panel_features Optional training-selected feature panel. If omitted,
#'   a train-only CLR panel is selected inside this function.
#' @param eta_grid,l2_grid Hyperparameter grids for train-only nested CV.
#' @param panel_size Number of features selected when `panel_features` is NULL.
#' @param inner_folds Number of inner folds for hyperparameter tuning.
#' @param learning_rate Gradient-descent learning rate.
#' @param epochs Number of optimization epochs.
#' @param tune Logical. If TRUE, tune `eta` and `l2` on training data.
#' @param seed Integer random seed.
#' @param ... Reserved for scorer-interface compatibility.
#' @return A frozen Group-DRO fit object consumed by
#'   `score_group_dro_scorer()`.
fit_group_dro_scorer <- function(expr_matrix,
                                 y,
                                 sample_meta,
                                 kit_col = "kit_family",
                                 biofluid_col = "biofluid",
                                 panel_features = NULL,
                                 panel_size = 20L,
                                 eta_grid = c(0.01, 0.05, 0.10),
                                 l2_grid = c(0.001, 0.01, 0.10),
                                 inner_folds = 3L,
                                 learning_rate = 0.05,
                                 epochs = 200L,
                                 tune = TRUE,
                                 seed = 42L,
                                 ...) {
  x <- .gdro_as_matrix(expr_matrix)
  y <- .gdro_as_binary(y)
  if (length(y) != nrow(x)) stop("y length must match expr_matrix rows.")
  meta <- .gdro_validate_meta(sample_meta, nrow(x), kit_col, biofluid_col)

  tuned <- NULL
  eta <- eta_grid[1L]
  l2 <- l2_grid[1L]
  if (isTRUE(tune) && length(eta_grid) * length(l2_grid) > 1L &&
      nrow(x) >= 24L && min(table(y)) >= max(4L, inner_folds)) {
    tuned <- .gdro_tune_grid(
      x, y, meta, kit_col, biofluid_col,
      eta_grid = eta_grid, l2_grid = l2_grid, panel_size = panel_size,
      inner_folds = inner_folds, learning_rate = learning_rate,
      epochs = max(30L, floor(epochs / 2L)), seed = seed)
    eta <- tuned$best$eta[1L]
    l2 <- tuned$best$l2[1L]
  }

  clr_fit <- .gdro_clr_fit(x)
  z_all <- .gdro_clr_apply(x, clr_fit)
  if (is.null(panel_features)) {
    feature_directions <- .gdro_select_panel(z_all, y, panel_size = panel_size)
    panel_features <- names(feature_directions)
  } else {
    panel_features <- intersect(as.character(panel_features), colnames(z_all))
    if (length(panel_features) == 0L) stop("panel_features do not overlap training CLR matrix.")
    feature_directions <- .gdro_select_panel(z_all[, panel_features, drop = FALSE],
                                             y, panel_size = length(panel_features))
    feature_directions <- feature_directions[panel_features]
    feature_directions[!is.finite(feature_directions)] <- 1
  }
  z_panel <- z_all[, panel_features, drop = FALSE]
  scaler <- .gdro_scale_fit(z_panel)
  z_scaled <- .gdro_scale_apply(z_panel, scaler)
  groups <- .gdro_make_groups(meta, kit_col, biofluid_col)

  dro <- .gdro_logistic_group_dro(
    z_scaled, y, groups, eta = eta, l2 = l2,
    learning_rate = learning_rate, epochs = epochs, seed = seed)
  structure(list(
    method = "group_dro_kit_robust",
    reference = "Sagawa et al. ICLR 2020, arXiv:1911.08731",
    direct_scorer = TRUE,
    clr_fit = clr_fit,
    scaler = scaler,
    panel_features = panel_features,
    feature_directions = feature_directions,
    beta = dro$beta,
    intercept = dro$intercept,
    group_weights = dro$group_weights,
    group_losses = dro$group_losses,
    train_groups = sort(unique(groups)),
    kit_col = kit_col,
    biofluid_col = biofluid_col,
    eta = eta,
    l2 = l2,
    eta_grid = eta_grid,
    l2_grid = l2_grid,
    tuning = tuned,
    learning_rate = learning_rate,
    epochs = as.integer(epochs),
    n_train = nrow(x),
    n_features = ncol(x),
    n_panel_features = length(panel_features),
    n_groups = length(unique(groups)),
    train_group_table = sort(table(groups), decreasing = TRUE)
  ), class = c("group_dro_scorer", "direct_scorer"))
}

#' Score held-out samples with a frozen Group-DRO scorer.
#'
#' @param fit Output from `fit_group_dro_scorer()`.
#' @param expr_matrix Numeric held-out matrix.
#' @param sample_meta Optional held-out metadata; used only for audit
#'   attributes, never for fitting.
#' @param ... Reserved for scorer-interface compatibility.
#' @return Numeric probability vector, higher meaning more disease-like.
score_group_dro_scorer <- function(fit, expr_matrix, sample_meta = NULL, ...) {
  if (is.null(fit$beta) || is.null(fit$clr_fit) || is.null(fit$scaler)) {
    stop("fit is not a valid group_dro_scorer object.")
  }
  z <- .gdro_clr_apply(expr_matrix, fit$clr_fit, fit$panel_features)
  z <- .gdro_scale_apply(z, fit$scaler)
  beta <- fit$beta[colnames(z)]
  out <- .gdro_sigmoid(as.numeric(fit$intercept + z %*% beta))
  if (!is.null(sample_meta) && all(c(fit$kit_col, fit$biofluid_col) %in% names(sample_meta))) {
    attr(out, "test_groups") <- .gdro_make_groups(sample_meta, fit$kit_col, fit$biofluid_col)
    attr(out, "unseen_test_groups") <- setdiff(unique(attr(out, "test_groups")), fit$train_groups)
  }
  out
}

.gdro_rclr_baseline_fit <- function(expr_matrix, y, panel_size = 20L) {
  x <- .gdro_as_matrix(expr_matrix)
  y <- .gdro_as_binary(y)
  clr_fit <- .gdro_clr_fit(x)
  z <- .gdro_clr_apply(x, clr_fit)
  directions <- .gdro_select_panel(z, y, panel_size = panel_size)
  list(clr_fit = clr_fit,
       panel_features = names(directions),
       feature_directions = directions)
}

.gdro_rclr_baseline_score <- function(fit, expr_matrix) {
  z <- .gdro_clr_apply(expr_matrix, fit$clr_fit, fit$panel_features)
  keep <- intersect(names(fit$feature_directions), colnames(z))
  as.numeric(z[, keep, drop = FALSE] %*% fit$feature_directions[keep])
}

.gdro_load_prepared_ngs <- function(input_path = NULL, kit_metadata_path = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  input_path <- .singlesample_extdata_path("dann_domain_invariant_input.tsv",
                                     input_path,
                                     "prepared Group-DRO input")
  if (!file.exists(input_path)) stop("Missing prepared input: ", input_path)
  dt <- data.table::fread(input_path, sep = "\t", quote = "")
  meta_cols <- c("sample_id", "accession", "cancer_type", "y", "baseline_rclr_auc")
  feat_cols <- setdiff(names(dt), meta_cols)
  x <- as.matrix(dt[, ..feat_cols])
  storage.mode(x) <- "double"
  prepared_scale <- "nonnegative_as_supplied"
  if (any(x < 0, na.rm = TRUE)) {
    x <- exp(x)
    prepared_scale <- "exp_applied_to_negative_prepared_features"
  }
  colnames(x) <- feat_cols
  rownames(x) <- dt$sample_id

  kit_path <- .singlesample_extdata_path("table_S43_matched_null_kit_metadata.tsv",
                                   kit_metadata_path, "S43 kit metadata",
                                   required = FALSE)
  kit <- data.table::data.table(cohort = unique(dt$accession),
                                biofluid = "unknown",
                                library_kit = "unknown",
                                provenance_block = unique(dt$accession))
  if (!is.null(kit_path) && file.exists(kit_path)) {
    md <- data.table::fread(kit_path, sep = "\t", quote = "")
    kit <- md[cohort %in% unique(dt$accession),
              .(cohort, biofluid, library_kit, provenance_block)]
  }
  meta <- merge(dt[, .(sample_id, accession, cancer_type, y)],
                kit, by.x = "accession", by.y = "cohort", all.x = TRUE, sort = FALSE)
  meta[, kit_family := library_kit]
  meta[is.na(kit_family) | !nzchar(kit_family), kit_family := "unknown"]
  meta[is.na(biofluid) | !nzchar(biofluid), biofluid := "unknown"]
  meta[is.na(provenance_block) | !nzchar(provenance_block), provenance_block := accession]
  meta <- meta[match(dt$sample_id, sample_id)]
  list(x = x, y = as.integer(dt$y), meta = meta,
       cohorts = sort(unique(dt$accession)),
       source = paste(input_path, prepared_scale, sep = ";"),
       universe = "prepared_ngs_5_cohorts")
}

.gdro_load_primary22 <- function(log_fun = message) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  if (!exists(".s45_primary_cohorts", mode = "function") ||
      !exists(".s45_make_metadata", mode = "function") ||
      !exists(".s45_load_data", mode = "function")) {
    stop("Primary-22 benchmark helpers are unavailable in this package build.",
         call. = FALSE)
  }
  primary <- .s45_primary_cohorts()
  meta22 <- .s45_make_metadata(primary)
  data_list <- .s45_load_data(meta22, log_fun)
  common <- Reduce(intersect, lapply(data_list, function(x) colnames(x$X)))
  if (length(common) < 10L) {
    stop("Primary-22 common feature intersection too small for Group-DRO: n=", length(common))
  }
  rows <- list()
  mats <- list()
  for (acc in names(data_list)) {
    dat <- data_list[[acc]]
    mr <- meta22[cohort == acc][1L]
    X <- as.matrix(dat$X[, common, drop = FALSE])
    storage.mode(X) <- "double"
    X <- pmax(X, 0)
    sm <- .s45_sample_meta(dat, mr)
    sm$accession <- acc
    sm$cancer_type <- mr$cancer_type
    sm$provenance_block <- mr$provenance_block
    sm$y <- as.integer(dat$y)
    mats[[acc]] <- X
    rows[[acc]] <- data.table::as.data.table(sm)
  }
  x <- do.call(rbind, mats)
  meta <- data.table::rbindlist(rows, fill = TRUE)
  rownames(x) <- meta$sample_id
  list(x = x, y = meta$y, meta = meta,
       cohorts = sort(unique(meta$accession)),
       source = "E1 primary matched-null loader",
       universe = "primary22_common_feature_LOCO")
}

.gdro_eval_fold <- function(held, obj, panel_size, eta_grid, l2_grid,
                            inner_folds, learning_rate, epochs, seed) {
  test <- obj$meta$accession == held
  train <- !test
  y_train <- obj$y[train]
  y_test <- obj$y[test]
  held_meta <- obj$meta[test]
  held_groups <- .gdro_make_groups(held_meta, "kit_family", "biofluid")
  held_group_label <- if (length(unique(held_groups)) == 1L) unique(held_groups) else "multiple_test_groups"
  if (length(unique(y_test)) != 2L || length(unique(y_train)) != 2L) {
    return(data.table::data.table(
      benchmark_component = "LOCO",
      method = c("rCLR_baseline", "group_dro_kit_robust"),
      held_out = held,
      cancer_type = unique(held_meta$cancer_type)[1L],
      kit_label = unique(held_meta$kit_family)[1L],
      biofluid = unique(held_meta$biofluid)[1L],
      kit_biofluid_group = held_group_label,
      n_train = sum(train),
      n_test = sum(test),
      n_test_pos = sum(y_test == 1L),
      n_test_neg = sum(y_test == 0L),
      observed_auc = NA_real_,
      auc_se = NA_real_,
      auc_ci_low = NA_real_,
      auc_ci_high = NA_real_,
      status = "not_evaluable",
      status_reason = "held-out or training split lacks both classes",
      leakage_flag = FALSE))
  }
  base <- .gdro_rclr_baseline_fit(obj$x[train, , drop = FALSE], y_train,
                                  panel_size = panel_size)
  base_score <- .gdro_rclr_baseline_score(base, obj$x[test, , drop = FALSE])
  dro <- fit_group_dro_scorer(
    obj$x[train, , drop = FALSE], y_train, obj$meta[train],
    kit_col = "kit_family", biofluid_col = "biofluid",
    panel_features = base$panel_features, panel_size = panel_size,
    eta_grid = eta_grid, l2_grid = l2_grid, inner_folds = inner_folds,
    learning_rate = learning_rate, epochs = epochs, tune = TRUE, seed = seed)
  dro_score <- score_group_dro_scorer(dro, obj$x[test, , drop = FALSE], obj$meta[test])

  mk <- function(method, score, extra_status = "") {
    auc <- .gdro_safe_auc(y_test, score)
    ci <- .gdro_hanley_mcneil(auc, sum(y_test == 1L), sum(y_test == 0L))
    data.table::data.table(
      benchmark_component = "LOCO",
      method = method,
      held_out = held,
      cancer_type = unique(held_meta$cancer_type)[1L],
      kit_label = unique(held_meta$kit_family)[1L],
      biofluid = unique(held_meta$biofluid)[1L],
      kit_biofluid_group = held_group_label,
      n_train = sum(train),
      n_test = sum(test),
      n_test_pos = sum(y_test == 1L),
      n_test_neg = sum(y_test == 0L),
      panel_size = length(base$panel_features),
      n_features = ncol(obj$x),
      n_train_groups = if (method == "group_dro_kit_robust") dro$n_groups else NA_integer_,
      eta = if (method == "group_dro_kit_robust") dro$eta else NA_real_,
      l2 = if (method == "group_dro_kit_robust") dro$l2 else NA_real_,
      observed_auc = auc,
      auc_se = ci[["se"]],
      auc_ci_low = ci[["low"]],
      auc_ci_high = ci[["high"]],
      average_auc = auc,
      worst_group_auc = auc,
      worst_group_definition = held_group_label,
      rclr_baseline_auc = if (method == "rCLR_baseline") auc else .gdro_safe_auc(y_test, base_score),
      lift_vs_rCLR_baseline = if (method == "rCLR_baseline") 0 else auc - .gdro_safe_auc(y_test, base_score),
      worst_group_lift_vs_rCLR = if (method == "rCLR_baseline") 0 else auc - .gdro_safe_auc(y_test, base_score),
      framework = if (method == "group_dro_kit_robust") "base_R_full_batch_group_DRO_logistic" else "train_only_rCLR_panel_sum",
      status = if (is.finite(auc)) "evaluable" else "not_evaluable",
      status_reason = if (is.finite(auc)) extra_status else "nonfinite_or_constant_score",
      leakage_flag = is.finite(auc) && auc > 0.95,
      provenance_block = unique(held_meta$provenance_block)[1L])
  }
  rows <- data.table::rbindlist(list(
    mk("rCLR_baseline", base_score),
    mk("group_dro_kit_robust", dro_score)
  ), fill = TRUE)
  attr(rows, "sample_scores") <- data.table::data.table(
    sample_id = held_meta$sample_id,
    held_out = held,
    y = y_test,
    kit_biofluid_group = held_groups,
    rCLR_baseline = base_score,
    group_dro_kit_robust = dro_score)
  rows
}

.gdro_add_summary_rows <- function(out, sample_scores, boot_B, seed) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  pred <- data.table::melt(sample_scores,
                           id.vars = c("sample_id", "held_out", "y", "kit_biofluid_group"),
                           variable.name = "method", value.name = "score")
  group_auc <- pred[, {
    auc <- if (length(unique(y)) == 2L) .gdro_safe_auc(y, score) else NA_real_
    .(group_auc = auc, n_group_samples = .N, n_group_pos = sum(y == 1L), n_group_neg = sum(y == 0L))
  }, by = .(method, kit_biofluid_group)]
  worst <- group_auc[is.finite(group_auc), .SD[which.min(group_auc)], by = method]
  pooled <- .gdro_re_pool_auc(out[benchmark_component == "LOCO"])
  if (nrow(pooled) == 0L) return(out)
  summary <- merge(pooled, worst[, .(method, worst_group_auc = group_auc,
                                     worst_group_definition = kit_biofluid_group,
                                     worst_group_n = n_group_samples,
                                     worst_group_pos = n_group_pos,
                                     worst_group_neg = n_group_neg)],
                   by = "method", all.x = TRUE)
  base_worst <- summary[method == "rCLR_baseline", worst_group_auc][1L]
  summary[, `:=`(
    benchmark_component = "RE_pooling_cluster_bootstrap",
    held_out = "pooled",
    cancer_type = "multi-cancer",
    kit_label = "all",
    biofluid = "all",
    kit_biofluid_group = "all",
    observed_auc = average_auc,
    auc_se = average_auc_se,
    auc_ci_low = average_auc_ci_low,
    auc_ci_high = average_auc_ci_high,
    rclr_baseline_auc = summary[method == "rCLR_baseline", average_auc][1L],
    worst_group_lift_vs_rCLR = worst_group_auc - base_worst,
    status = "exploratory_requires_review",
    status_reason = "Random-effects pooled LOCO AUC with provenance-block collapse; requires triple-consensus and plausibility-check",
    leakage_flag = FALSE
  )]
  diff_ci <- .gdro_cluster_boot_diff(out, B = boot_B, seed = seed)
  summary[method == "group_dro_kit_robust",
          `:=`(lift_vs_rCLR_baseline = diff_ci[["mean"]],
               lift_ci_low = diff_ci[["low"]],
               lift_ci_high = diff_ci[["high"]])]
  summary[method == "rCLR_baseline",
          `:=`(lift_vs_rCLR_baseline = 0,
               lift_ci_low = 0,
               lift_ci_high = 0)]
  data.table::rbindlist(list(out, summary), fill = TRUE)
}

.gdro_write_outputs <- function(out, out_tsv, out_png, out_pdf, out_md,
                                universe, n_cohorts, requested_full22) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(out, out_tsv, sep = "\t", na = "NA", quote = FALSE)

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  plot_dt <- out[benchmark_component == "RE_pooling_cluster_bootstrap" &
                   method %in% c("rCLR_baseline", "group_dro_kit_robust")]
  if (nrow(plot_dt) == 0L) {
    plot_dt <- out[benchmark_component == "LOCO" &
                     method %in% c("rCLR_baseline", "group_dro_kit_robust"),
                   .(average_auc = mean(observed_auc, na.rm = TRUE),
                     worst_group_auc = min(worst_group_auc, na.rm = TRUE),
                     status = paste(unique(status), collapse = ";")),
                   by = method]
  }
  plot_dt[, method_label := factor(method,
    levels = c("rCLR_baseline", "group_dro_kit_robust"),
    labels = c("rCLR baseline", "Group-DRO"))]
  fig <- ggplot2::ggplot(plot_dt,
                         ggplot2::aes(x = average_auc, y = worst_group_auc,
                                      label = method_label, colour = method_label)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                         colour = "grey55") +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::geom_text(nudge_y = 0.025, show.legend = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0.45, 1), name = "Average AUC") +
    ggplot2::scale_y_continuous(limits = c(0.45, 1), name = "Worst kit x biofluid AUC") +
    ggplot2::labs(title = "Supplementary Figure S54. Group-DRO worst-group scoring",
                  colour = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  ggplot2::ggsave(out_png, fig, width = 6.5, height = 4.8, dpi = 300, bg = "white")
  ggplot2::ggsave(out_pdf, fig, width = 6.5, height = 4.8, bg = "white")

  leak <- out[benchmark_component == "LOCO" & leakage_flag == TRUE]
  denom_sentence <- if (isTRUE(requested_full22) && n_cohorts == 22L) {
    "The run used the requested 22-cohort primary matched-null universe."
  } else {
    sprintf("This artifact is a denominator-limited implementation/smoke artifact over %d cohort(s) in `%s`, not the full requested 22-cohort acceptance benchmark.", n_cohorts, universe)
  }
  cap <- c(
    "# Figure S54. Group-DRO kit x biofluid robust scoring",
    "",
    "The Group-DRO scorer fits a train-only CLR logistic classifier with exponentiated-gradient weights over kit x biofluid groups. Hyperparameters (`eta`, L2) are selected by train-only nested CV; held-out labels are used only for AUC evaluation.",
    "",
    denom_sentence,
    "",
    sprintf("Plausibility gate: %d LOCO row(s) had AUC > 0.95 and require targeted leakage review if nonzero.", nrow(leak)),
    "",
    "The table reports both average AUC and worst kit x biofluid AUC. These outputs are not validated evidence until the code, denominators, and leakage screen are independently reviewed."
  )
  writeLines(cap, out_md)
}

run_group_dro_benchmark <- function(
    mode = Sys.getenv("GDRO_BENCHMARK_MODE", "prepared_ngs"),
    input_path = NULL,
    kit_metadata_path = NULL,
    output_dir = tempdir(),
    out_tsv = NULL,
    out_png = NULL,
    out_pdf = NULL,
    out_md = NULL,
    panel_size = as.integer(Sys.getenv("GDRO_PANEL_SIZE", "20")),
    eta_grid = as.numeric(strsplit(Sys.getenv("GDRO_ETA_GRID", "0.01,0.05,0.1"), ",")[[1L]]),
    l2_grid = as.numeric(strsplit(Sys.getenv("GDRO_L2_GRID", "0.001,0.01,0.1"), ",")[[1L]]),
    inner_folds = as.integer(Sys.getenv("GDRO_INNER_FOLDS", "3")),
    learning_rate = as.numeric(Sys.getenv("GDRO_LR", "0.05")),
    epochs = as.integer(Sys.getenv("GDRO_EPOCHS", "160")),
    boot_B = as.integer(Sys.getenv("GDRO_BOOT_B", "500")),
    seed = 42L) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  set.seed(seed)
  if (is.null(out_tsv)) out_tsv <- .singlesample_default_output(output_dir, "table_S54_group_dro_benchmarks.tsv")
  if (is.null(out_png)) out_png <- .singlesample_default_output(output_dir, "figure_S54_group_dro.png")
  if (is.null(out_pdf)) out_pdf <- .singlesample_default_output(output_dir, "figure_S54_group_dro.pdf")
  if (is.null(out_md)) out_md <- .singlesample_default_output(output_dir, "figure_S54_group_dro.md")
  log_fun <- function(msg) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), msg))
  log_fun(sprintf("Group-DRO benchmark start mode=%s seed=%d", mode, seed))
  obj <- switch(mode,
    prepared_ngs = .gdro_load_prepared_ngs(input_path, kit_metadata_path),
    primary22 = .gdro_load_primary22(log_fun),
    stop("Unknown GDRO_BENCHMARK_MODE: ", mode)
  )
  log_fun(sprintf("Universe=%s cohorts=%d samples=%d features=%d",
                  obj$universe, length(obj$cohorts), nrow(obj$x), ncol(obj$x)))
  fold_rows <- list()
  sample_rows <- list()
  for (i in seq_along(obj$cohorts)) {
    held <- obj$cohorts[i]
    log_fun(sprintf("LOCO %02d/%02d held_out=%s", i, length(obj$cohorts), held))
    rows <- .gdro_eval_fold(held, obj, panel_size, eta_grid, l2_grid,
                            inner_folds, learning_rate, epochs, seed + i)
    sample_rows[[held]] <- attr(rows, "sample_scores")
    attr(rows, "sample_scores") <- NULL
    fold_rows[[held]] <- rows
  }
  out <- data.table::rbindlist(fold_rows, fill = TRUE)
  sample_scores <- data.table::rbindlist(sample_rows, fill = TRUE)
  out <- .gdro_add_summary_rows(out, sample_scores, boot_B = boot_B, seed = seed)
  out[, `:=`(benchmark_universe = obj$universe,
             requested_cohorts = if (mode == "primary22") 22L else 22L,
             evaluated_cohorts = length(obj$cohorts),
             full_22_complete = length(obj$cohorts) == 22L && mode == "primary22",
             source_input = obj$source)]
  data.table::setcolorder(out, intersect(c(
    "benchmark_component", "method", "held_out", "cancer_type", "kit_label",
    "biofluid", "kit_biofluid_group", "n_train", "n_test", "n_test_pos",
    "n_test_neg", "panel_size", "n_features", "n_train_groups", "eta", "l2",
    "observed_auc", "auc_se", "auc_ci_low", "auc_ci_high",
    "average_auc", "average_auc_ci_low", "average_auc_ci_high",
    "worst_group_auc", "worst_group_definition", "worst_group_n",
    "worst_group_pos", "worst_group_neg", "rclr_baseline_auc",
    "lift_vs_rCLR_baseline", "lift_ci_low", "lift_ci_high",
    "worst_group_lift_vs_rCLR", "tau2", "i2", "n_provenance_blocks",
    "framework", "status", "status_reason", "leakage_flag",
    "benchmark_universe", "requested_cohorts", "evaluated_cohorts",
    "full_22_complete", "source_input", "provenance_block"
  ), names(out)))
  .gdro_write_outputs(out, out_tsv, out_png, out_pdf, out_md,
                      universe = obj$universe, n_cohorts = length(obj$cohorts),
                      requested_full22 = mode == "primary22")
  log_fun(sprintf("Wrote table: %s rows=%d", out_tsv, nrow(out)))
  log_fun(sprintf("Wrote figure: %s and %s", out_png, out_pdf))
  log_fun(sprintf("Wrote caption: %s", out_md))
  invisible(out)
}

if (sys.nframe() == 0L) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  run_group_dro_benchmark()
}
