# Shared helpers for learned kit-aware within-sample scorers.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.lka_as_matrix <- function(x, name = "expr") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop(name, " must be a numeric matrix/data.frame.")
  storage.mode(x) <- "double"
  if (is.null(colnames(x))) colnames(x) <- paste0("feature_", seq_len(ncol(x)))
  x[!is.finite(x)] <- NA_real_
  x
}

.lka_as_binary <- function(y) {
  if (is.factor(y)) y <- as.integer(y) - 1L
  y <- as.integer(y)
  if (!all(y %in% c(0L, 1L), na.rm = TRUE)) {
    stop("Disease labels must be binary 0/1 or a two-level factor.")
  }
  y
}

.lka_factor_state <- function(x, unknown = "unknown") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- unknown
  levels <- sort(unique(x))
  list(labels = x, levels = levels, index = match(x, levels))
}

.lka_standardize_fit <- function(x) {
  x <- .lka_as_matrix(x)
  center <- apply(x, 2L, median, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  madv <- apply(abs(sweep(x, 2L, center, "-")), 2L, median, na.rm = TRUE)
  scale <- 1.4826 * madv
  fallback <- median(scale[is.finite(scale) & scale > 1e-8], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 1e-8) fallback <- 1
  scale[!is.finite(scale) | scale <= 1e-8] <- fallback
  list(center = center, scale = scale, features = colnames(x))
}

.lka_standardize_apply <- function(x, scaler) {
  x <- .lka_as_matrix(x)
  missing <- setdiff(scaler$features, colnames(x))
  if (length(missing) > 0L) {
    add <- matrix(NA_real_, nrow = nrow(x), ncol = length(missing),
                  dimnames = list(rownames(x), missing))
    x <- cbind(x, add)
  }
  x <- x[, scaler$features, drop = FALSE]
  for (j in seq_len(ncol(x))) {
    bad <- !is.finite(x[, j])
    if (any(bad)) x[bad, j] <- scaler$center[j]
  }
  z <- sweep(sweep(x, 2L, scaler$center, "-"), 2L, scaler$scale, "/")
  z[!is.finite(z)] <- 0
  pmax(pmin(z, 8), -8)
}

.lka_safe_auc <- function(y, score) {
  y <- .lka_as_binary(y)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (length(y) < 4L || length(unique(y)) != 2L || length(unique(score)) < 2L) {
    return(NA_real_)
  }
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

.lka_hanley_mcneil_ci <- function(auc, n_pos, n_neg) {
  if (!is.finite(auc) || n_pos < 1L || n_neg < 1L) {
    return(c(se = NA_real_, low = NA_real_, high = NA_real_))
  }
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  var <- (auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) +
            (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg)
  if (!is.finite(var) || var < 0) return(c(se = NA_real_, low = NA_real_, high = NA_real_))
  se <- sqrt(var)
  c(se = se, low = max(0, auc - 1.96 * se), high = min(1, auc + 1.96 * se))
}

.lka_stratified_folds <- function(y, k = 5L, seed = 42L) {
  y <- .lka_as_binary(y)
  set.seed(seed)
  fold <- integer(length(y))
  for (cl in sort(unique(y))) {
    idx <- which(y == cl)
    fold[idx] <- sample(rep_len(seq_len(k), length(idx)))[seq_along(idx)]
  }
  fold
}

.lka_select_panel_fold_aware <- function(x, y, panel_size = 20L, k_folds = 3L,
                                         seed = 42L, min_detection = 0) {
  x <- .lka_as_matrix(x)
  y <- .lka_as_binary(y)
  keep <- colMeans(is.finite(x) & abs(x) > min_detection, na.rm = TRUE) > 0.05
  if (sum(keep) < panel_size) keep[] <- TRUE
  cand <- which(keep)
  fold <- .lka_stratified_folds(y, k = min(k_folds, max(2L, min(table(y)))), seed = seed)
  cv_auc <- rep(NA_real_, length(cand))
  for (ii in seq_along(cand)) {
    j <- cand[ii]
    fold_auc <- vapply(sort(unique(fold)), function(ff) {
      test <- fold == ff
      .lka_safe_auc(y[test], x[test, j])
    }, numeric(1L))
    cv_auc[ii] <- median(pmax(fold_auc, 1 - fold_auc), na.rm = TRUE)
  }
  cv_auc[!is.finite(cv_auc)] <- 0.5
  ord <- order(cv_auc, decreasing = TRUE)
  panel_idx <- cand[ord[seq_len(min(panel_size, length(ord)))]]
  panel <- colnames(x)[panel_idx]
  directions <- vapply(panel_idx, function(j) {
    mu1 <- mean(x[y == 1L, j], na.rm = TRUE)
    mu0 <- mean(x[y == 0L, j], na.rm = TRUE)
    if (is.finite(mu1) && is.finite(mu0) && mu1 >= mu0) 1 else -1
  }, numeric(1L))
  names(directions) <- panel
  list(panel = panel, directions = directions,
       selection = "training-only stratified inner-CV univariate AUC")
}

.lka_feature_weights <- function(projected_train, y, panel_features = NULL) {
  projected_train <- .lka_as_matrix(projected_train, "projected_train")
  y <- .lka_as_binary(y)
  if (is.null(panel_features)) panel_features <- colnames(projected_train)
  panel_features <- intersect(panel_features, colnames(projected_train))
  weights <- setNames(numeric(ncol(projected_train)), colnames(projected_train))
  for (f in panel_features) {
    d <- mean(projected_train[y == 1L, f], na.rm = TRUE) -
      mean(projected_train[y == 0L, f], na.rm = TRUE)
    weights[f] <- if (is.finite(d) && d >= 0) 1 else -1
  }
  weights
}

.lka_score_projection <- function(projected, weights, panel_features = NULL) {
  projected <- .lka_as_matrix(projected, "projected")
  if (is.null(panel_features)) panel_features <- names(weights)
  keep <- intersect(panel_features, intersect(colnames(projected), names(weights)))
  if (length(keep) == 0L) return(rep(NA_real_, nrow(projected)))
  as.numeric(projected[, keep, drop = FALSE] %*% weights[keep])
}

.lka_matrix_for_features <- function(x, features) {
  x <- .lka_as_matrix(x)
  out <- matrix(0, nrow = nrow(x), ncol = length(features),
                dimnames = list(rownames(x), features))
  hit <- intersect(features, colnames(x))
  if (length(hit) > 0L) out[, hit] <- x[, hit, drop = FALSE]
  out
}

.lka_torch_available <- function() {
  requireNamespace("torch", quietly = TRUE)
}

.lka_one_hot <- function(index, n_levels) {
  m <- matrix(0, nrow = length(index), ncol = n_levels)
  m[cbind(seq_along(index), index)] <- 1
  m
}

.lka_fit_linear_fallback <- function(train_expr, y, kit_labels, cohort_labels,
                                     panel_features = NULL,
                                     method = "linear_residual_fallback") {
  x <- .lka_as_matrix(train_expr)
  y <- .lka_as_binary(y)
  scaler <- .lka_standardize_fit(x)
  z <- .lka_standardize_apply(x, scaler)
  if (is.null(panel_features)) panel_features <- colnames(z)
  weights <- .lka_feature_weights(z, y, panel_features)
  list(method = method, framework = "base_R_linear_fallback",
       scaler = scaler, feature_weights = weights,
       panel_features = intersect(panel_features, colnames(z)),
       kit_levels = sort(unique(as.character(kit_labels))),
       cohort_levels = sort(unique(as.character(cohort_labels))))
}

.lka_score_linear_fallback <- function(model, test_expr, panel_features = NULL) {
  z <- .lka_standardize_apply(test_expr, model$scaler)
  .lka_score_projection(z, model$feature_weights, panel_features %||% model$panel_features)
}

.lka_draw_matched_panel <- function(train_expr, observed_panel, panel_size = length(observed_panel),
                                    seed = 42L) {
  x <- .lka_as_matrix(train_expr)
  observed_panel <- intersect(observed_panel, colnames(x))
  pool <- setdiff(colnames(x), observed_panel)
  if (length(pool) < panel_size) return(character(0))
  det <- colMeans(is.finite(x) & x != 0, na.rm = TRUE)
  abn <- colMeans(abs(x), na.rm = TRUE)
  det_bin <- as.integer(cut(rank(det, ties.method = "average"), 4L, include.lowest = TRUE))
  abn_bin <- as.integer(cut(rank(abn, ties.method = "average"), 4L, include.lowest = TRUE))
  strata <- setNames(paste(det_bin, abn_bin, sep = ":"), colnames(x))
  need <- table(strata[observed_panel])
  set.seed(seed)
  picked <- character(0)
  for (s in names(need)) {
    candidates <- setdiff(names(strata)[strata == s], c(observed_panel, picked))
    n_pick <- min(as.integer(need[[s]]), length(candidates))
    if (n_pick > 0L) picked <- c(picked, sample(candidates, n_pick))
  }
  short <- panel_size - length(picked)
  if (short > 0L) {
    extra <- setdiff(pool, picked)
    if (length(extra) < short) return(character(0))
    picked <- c(picked, sample(extra, short))
  }
  picked[seq_len(panel_size)]
}

.lka_matched_null <- function(model, train_expr, test_expr, y_test, observed_panel,
                              score_fun, K = 100L, seed = 42L) {
  observed_score <- score_fun(model, test_expr, observed_panel)
  observed_auc <- .lka_safe_auc(y_test, observed_score)
  null_auc <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    p_k <- .lka_draw_matched_panel(train_expr, observed_panel,
                                   panel_size = length(observed_panel),
                                   seed = seed + k)
    if (length(p_k) == 0L) next
    null_auc[k] <- .lka_safe_auc(y_test, score_fun(model, test_expr, p_k))
  }
  n_valid <- sum(is.finite(null_auc))
  p_emp <- if (is.finite(observed_auc) && n_valid > 0L)
    (1 + sum(null_auc >= observed_auc, na.rm = TRUE)) / (n_valid + 1) else NA_real_
  list(observed_auc = observed_auc, p_emp = p_emp, n_null_valid = n_valid,
       null_auc_median = median(null_auc, na.rm = TRUE))
}

.lka_re_pool <- function(dt) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  d <- data.table::as.data.table(dt)
  d <- d[is.finite(observed_auc) & is.finite(auc_se) & auc_se > 0]
  if (nrow(d) == 0L) return(data.table::data.table())
  d[, {
    yi <- observed_auc
    vi <- auc_se^2
    wi <- 1 / vi
    fixed <- sum(wi * yi) / sum(wi)
    q <- sum(wi * (yi - fixed)^2)
    c_val <- sum(wi) - sum(wi^2) / sum(wi)
    tau2 <- max(0, (q - (.N - 1)) / c_val)
    wre <- 1 / (vi + tau2)
    pooled <- sum(wre * yi) / sum(wre)
    se <- sqrt(1 / sum(wre))
    list(n_loco = .N, re_auc = pooled, re_se = se,
         re_ci_low = max(0, pooled - 1.96 * se),
         re_ci_high = min(1, pooled + 1.96 * se),
         tau2 = tau2, i2 = if (q > 0) max(0, (q - (.N - 1)) / q) else 0)
  }, by = method]
}

.lka_cluster_bootstrap <- function(dt, B = 500L, seed = 42L) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  d <- data.table::as.data.table(dt)
  d <- d[is.finite(observed_auc)]
  if (nrow(d) == 0L) return(data.table::data.table())
  set.seed(seed)
  d[, {
    cohorts <- unique(held_out)
    vals <- replicate(B, {
      sample_coh <- sample(cohorts, length(cohorts), replace = TRUE)
      mean(observed_auc[match(sample_coh, held_out)], na.rm = TRUE)
    })
    list(cluster_boot_mean = mean(vals, na.rm = TRUE),
         cluster_boot_ci_low = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
         cluster_boot_ci_high = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
         cluster_boot_B = B)
  }, by = method]
}

.lka_load_s48_input <- function(path = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  path <- .singlesample_extdata_path("dann_domain_invariant_input.tsv", path,
                               "prepared DANN/learned-kit-aware input")
  if (!file.exists(path)) stop("Missing prepared DANN input: ", path)
  dt <- data.table::fread(path, sep = "\t", quote = "")
  meta_cols <- c("sample_id", "accession", "cancer_type", "y", "baseline_rclr_auc")
  feat_cols <- setdiff(names(dt), meta_cols)
  x <- as.matrix(dt[, ..feat_cols])
  colnames(x) <- feat_cols
  rownames(x) <- dt$sample_id
  list(dt = dt, x = x, y = as.integer(dt$y), features = feat_cols)
}

.lka_kit_map <- function(accessions, kit_metadata_path = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  s43 <- .singlesample_extdata_path("table_S43_matched_null_kit_metadata.tsv",
                              kit_metadata_path, "S43 kit metadata",
                              required = FALSE)
  out <- data.table::data.table(accession = unique(accessions),
                                kit_label = "unknown",
                                kit_source = "missing_table_S43",
                                biofluid = NA_character_)
  if (!is.null(s43) && file.exists(s43)) {
    md <- data.table::fread(s43, sep = "\t", quote = "")
    md <- md[cohort %in% unique(accessions)]
    if (nrow(md) > 0L) {
      md[, kit_label := ifelse(!is.na(library_kit) & nzchar(library_kit),
                               library_kit, "unknown")]
      md[, kit_source := library_kit_assignment_status]
      out <- merge(out[, !"kit_label"], md[, .(accession = cohort, kit_label,
                                                kit_source, biofluid)],
                   by = "accession", all.x = TRUE)
      out[is.na(kit_label) | !nzchar(kit_label), kit_label := "unknown"]
    }
  }
  out
}

run_learned_kit_aware_benchmark <- function(
    input_path = NULL,
    kit_metadata_path = NULL,
    output_dir = tempdir(),
    out_tsv = NULL,
    out_png = NULL,
    out_pdf = NULL,
    out_caption = NULL,
    seed = 42L,
    epochs = as.integer(Sys.getenv("LKA_EPOCHS", "8")),
    K_null = as.integer(Sys.getenv("LKA_K_NULL", "100")),
    panel_size = 20L) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (is.null(out_tsv)) out_tsv <- .singlesample_default_output(output_dir, "table_S48_learned_kit_aware_benchmarks.tsv")
  if (is.null(out_png)) out_png <- .singlesample_default_output(output_dir, "figure_S48_learned_methods.png")
  if (is.null(out_pdf)) out_pdf <- .singlesample_default_output(output_dir, "figure_S48_learned_methods.pdf")
  if (is.null(out_caption)) out_caption <- .singlesample_default_output(output_dir, "figure_S48_caption.md")
  required <- c("fit_dann_kit_extended", "fit_kit_conditional_vae",
                "fit_icp_per_kit")
  if (!all(vapply(required, exists, logical(1), mode = "function"))) {
    stop("Learned kit-aware scorer functions are unavailable.", call. = FALSE)
  }

  inp <- .lka_load_s48_input(input_path)
  dt <- inp$dt
  x_all <- inp$x
  y_all <- inp$y
  accessions <- dt$accession
  kits <- .lka_kit_map(unique(accessions), kit_metadata_path = kit_metadata_path)
  kit_by_acc <- stats::setNames(kits$kit_label, kits$accession)
  kit_labels <- unname(kit_by_acc[accessions])
  kit_labels[is.na(kit_labels)] <- "unknown"
  cohort_labels <- accessions
  accs <- sort(unique(accessions))

  rows <- list()
  null_rows <- list()
  method_specs <- list(
    dann_kit_extended = list(
      fit = fit_dann_kit_extended,
      score = score_dann_kit_extended,
      args = list(epochs = epochs, seed = seed, lambda_cohort = 0.15,
                  lambda_kit = 0.15, disease_weight = 1.0)
    ),
    kit_conditional_vae = list(
      fit = fit_kit_conditional_vae,
      score = score_kit_conditional_vae,
      args = list(epochs = epochs, seed = seed, beta_kl = 0.03,
                  disease_weight = 0.7)
    ),
    icp_per_kit = list(
      fit = fit_icp_per_kit,
      score = score_icp_per_kit,
      args = list(seed = seed, calibration_fraction = 0.25)
    )
  )

  for (held in accs) {
    message(sprintf("[%s] S48 held-out=%s", Sys.time(), held))
    train <- accessions != held
    test <- accessions == held
    x_train <- x_all[train, , drop = FALSE]
    x_test <- x_all[test, , drop = FALSE]
    y_train <- y_all[train]
    y_test <- y_all[test]
    if (length(unique(y_test)) != 2L) next
    panel <- .lka_select_panel_fold_aware(x_train, y_train, panel_size = panel_size,
                                          seed = seed + match(held, accs))$panel
    held_kit <- unique(kit_labels[test])[1L]
    held_cancer <- unique(dt$cancer_type[test])[1L]
    for (mn in names(method_specs)) {
      spec <- method_specs[[mn]]
      model <- tryCatch({
        do.call(spec$fit, c(list(train_expr = x_train, train_meta = NULL,
                                 train_disease_labels = y_train,
                                 train_kit_labels = kit_labels[train],
                                 train_cohort_labels = cohort_labels[train],
                                 panel_features = panel),
                            spec$args))
      }, error = function(e) {
        structure(list(error = conditionMessage(e)), class = "lka_fit_error")
      })
      if (inherits(model, "lka_fit_error")) {
        rows[[length(rows) + 1L]] <- data.table::data.table(
          benchmark_component = "LOCO",
          method = mn, held_out = held, cancer_type = held_cancer,
          kit_label = held_kit, n_train = sum(train), n_test = sum(test),
          n_test_pos = sum(y_test == 1L), n_test_neg = sum(y_test == 0L),
          panel_size = length(panel), n_features = ncol(x_train),
          observed_auc = NA_real_, auc_se = NA_real_, auc_ci_low = NA_real_,
          auc_ci_high = NA_real_, p_emp_matched_null = NA_real_,
          n_null_valid = 0L, framework = NA_character_,
          status = "not_evaluable", status_reason = model$error,
          leakage_flag = FALSE, calibration_fraction = NA_real_)
        next
      }
      score <- spec$score(model, x_test, panel)
      auc <- .lka_safe_auc(y_test, score)
      ci <- .lka_hanley_mcneil_ci(auc, sum(y_test == 1L), sum(y_test == 0L))
      mnul <- .lka_matched_null(model, x_train, x_test, y_test, panel,
                                spec$score, K = K_null,
                                seed = seed + 1000L * match(held, accs))
      null_rows[[length(null_rows) + 1L]] <- data.table::data.table(
        benchmark_component = "kit_stratified_matched_null",
        method = mn, held_out = held, kit_label = held_kit,
        observed_auc = mnul$observed_auc, p_emp_matched_null = mnul$p_emp,
        n_null_valid = mnul$n_null_valid,
        null_auc_median = mnul$null_auc_median,
        K_requested = K_null)
      rows[[length(rows) + 1L]] <- data.table::data.table(
        benchmark_component = "LOCO",
        method = mn, held_out = held, cancer_type = held_cancer,
        kit_label = held_kit, n_train = sum(train), n_test = sum(test),
        n_test_pos = sum(y_test == 1L), n_test_neg = sum(y_test == 0L),
        panel_size = length(panel), n_features = ncol(x_train),
        observed_auc = auc, auc_se = ci[["se"]], auc_ci_low = ci[["low"]],
        auc_ci_high = ci[["high"]], p_emp_matched_null = mnul$p_emp,
        n_null_valid = mnul$n_null_valid, framework = model$framework %||% NA_character_,
        status = if (is.finite(auc)) "evaluable" else "not_evaluable",
        status_reason = if (is.finite(auc)) "" else "nonfinite_or_constant_score",
        leakage_flag = is.finite(auc) && auc > 0.95,
        calibration_fraction = model$calibration_fraction %||% NA_real_)
    }
  }

  out <- data.table::rbindlist(rows, fill = TRUE)
  if (length(null_rows) > 0L) {
    null_dt <- data.table::rbindlist(null_rows, fill = TRUE)
    out <- data.table::rbindlist(list(out, null_dt), fill = TRUE)
  }
  if (nrow(out) == 0L) stop("No S48 rows were generated.")
  leak_by_method <- out[benchmark_component == "LOCO" & leakage_flag == TRUE,
                        .N, by = method]
  if (nrow(leak_by_method) > 0L) {
    bad <- leak_by_method[N >= 2L, method]
    out[method %in% bad, leakage_flag := TRUE]
    out[method %in% bad & status == "evaluable",
        status_reason := "AUC>0.95 in multiple LOCO strata; leakage investigation required"]
  }

  re <- .lka_re_pool(out[benchmark_component == "LOCO"])
  boot <- .lka_cluster_bootstrap(out[benchmark_component == "LOCO"], B = 500L, seed = seed)
  pooled <- merge(re, boot, by = "method", all = TRUE)
  if (nrow(pooled) > 0L) {
    pooled[, `:=`(
      benchmark_component = "RE_pooling_cluster_bootstrap",
      held_out = "pooled",
      cancer_type = "multi-cancer_NGS",
      kit_label = "all",
      status = "exploratory_requires_review",
      status_reason = "DerSimonian-Laird random-effects plus cohort-cluster bootstrap; requires triple-consensus/plausibility review",
      framework = NA_character_,
      observed_auc = re_auc,
      auc_se = re_se,
      auc_ci_low = re_ci_low,
      auc_ci_high = re_ci_high,
      p_emp_matched_null = NA_real_,
      n_null_valid = NA_integer_,
      leakage_flag = FALSE
    )]
    common_cols <- union(names(out), names(pooled))
    out <- data.table::rbindlist(list(out[, intersect(common_cols, names(out)), with = FALSE],
                                     pooled[, intersect(common_cols, names(pooled)), with = FALSE]),
                                fill = TRUE)
  }
  kit_summary <- out[benchmark_component == "LOCO" & is.finite(observed_auc),
                     .(kit_stratum_n = .N,
                       kit_stratum_mean_auc = mean(observed_auc),
                       kit_stratum_min_auc = min(observed_auc),
                       kit_stratum_max_auc = max(observed_auc)),
                     by = .(method, kit_label)]
  if (nrow(kit_summary) > 0L) {
    kit_summary[, `:=`(
      benchmark_component = "kit_stratum_summary",
      held_out = "kit_summary",
      status = "descriptive",
      status_reason = "Kit strata are sparse and confounded with cohort/platform; descriptive only",
      observed_auc = kit_stratum_mean_auc
    )]
    out <- data.table::rbindlist(list(out, kit_summary), fill = TRUE)
  }

  data.table::setcolorder(out, intersect(c(
    "benchmark_component", "method", "held_out", "cancer_type", "kit_label",
    "n_train", "n_test", "n_test_pos", "n_test_neg", "panel_size",
    "n_features", "observed_auc", "auc_se", "auc_ci_low", "auc_ci_high",
    "p_emp_matched_null", "n_null_valid", "framework", "calibration_fraction",
    "re_auc", "tau2", "i2", "cluster_boot_mean", "cluster_boot_ci_low",
    "cluster_boot_ci_high", "status", "status_reason", "leakage_flag"
  ), names(out)))
  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(out, out_tsv, sep = "\t", na = "NA")

  ev <- out[benchmark_component == "LOCO" & status == "evaluable"]
  if (nrow(ev) > 0L) {
    ev[, method_label := factor(method, levels = c("dann_kit_extended",
                                                   "kit_conditional_vae",
                                                   "icp_per_kit"))]
    fig <- ggplot2::ggplot(ev, ggplot2::aes(x = held_out, y = observed_auc,
                                            color = kit_label, group = method_label)) +
      ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", color = "#9b1c1c") +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", color = "#777777") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = auc_ci_low, ymax = auc_ci_high),
                             width = 0.12, alpha = 0.55) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_wrap(~ method_label, ncol = 1L) +
      ggplot2::scale_y_continuous(limits = c(0, 1), name = "LOCO AUC") +
      ggplot2::labs(x = NULL, color = "Library-kit stratum",
                    title = "Supplementary Figure S6. Learned kit-aware scoring") +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                     plot.background = ggplot2::element_rect(fill = "white", color = NA),
                     panel.background = ggplot2::element_rect(fill = "white", color = NA))
    dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_png, fig, width = 8.5, height = 8.0, dpi = 220)
    ggplot2::ggsave(out_pdf, fig, width = 8.5, height = 8.0)
  }

  n_eval <- out[benchmark_component == "LOCO" & status == "evaluable", .N]
  n_leak <- out[benchmark_component == "LOCO" & leakage_flag == TRUE, .N]
  cap <- c(
    "# Figure S48. Learned kit-aware within-sample scoring",
    "",
    sprintf("Three learned/adversarial scorers were evaluated by training-only panel selection, leave-one-cohort-out testing across %d NGS cohorts, kit-stratified matched-null panels (K=%d), random-effects pooling, and cohort-cluster bootstrap.", length(accs), K_null),
    "",
    sprintf("Evaluable LOCO rows: %d. Leakage screen rows with AUC > 0.95: %d. Any method with AUC > 0.95 in multiple strata is flagged in the table and must be investigated before use.", n_eval, n_leak),
    "",
    "Kit strata are sparse and partly unknown/custom in this NGS subset; all kit-stratified summaries are descriptive and confounded with cohort/cancer type. These outputs are not acceptance-ready scientific evidence until /triple-consensus and /plausibility-check review the code, denominator accounting, and leakage screen."
  )
  writeLines(cap, out_caption)
  invisible(out)
}
