#!/usr/bin/env Rscript
# Sinkhorn OT cohort-barycenter projection for reviewer-package Table/Figure S55.
#
# Runner utilities below require caller-supplied manuscript inputs when used
# from the package; manuscript-workspace paths are intentionally not assumed.
#
# Gate status:
#   Reviewer-package artifact only. The generated code, table, figure, and
#   caption require /triple-consensus and /plausibility-check before manuscript
#   acceptance.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}

.sot_env_int <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(as.integer(default))
  out <- suppressWarnings(as.integer(val))
  if (!is.finite(out) || is.na(out) || out <= 0L) as.integer(default) else out
}

.sot_env_num <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(val))
  if (!is.finite(out) || is.na(out) || out <= 0) as.numeric(default) else out
}

.sot_env_logical <- function(name, default = FALSE) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(isTRUE(default))
  toupper(val) %in% c("1", "TRUE", "T", "YES", "Y")
}

.sot_env_num_vec <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(strsplit(val, ",", fixed = TRUE)[[1L]]))
  out <- out[is.finite(out) & out > 0]
  if (length(out) == 0L) as.numeric(default) else out
}

.sot_matrix <- function(x, arg = "expr_matrix") {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(arg, " must be a matrix/data.frame with samples in rows.")
  }
  out <- as.matrix(x)
  storage.mode(out) <- "double"
  if (nrow(out) == 0L || ncol(out) == 0L) stop(arg, " must be non-empty.")
  if (is.null(colnames(out))) stop(arg, " must have feature column names.")
  if (any(out < 0, na.rm = TRUE)) stop(arg, " must be non-negative.")
  out
}

.sot_meta <- function(meta, n, cohort_col = "cohort") {
  if (!is.data.frame(meta)) meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (nrow(meta) != n) stop("meta must have one row per sample.")
  if (!(cohort_col %in% names(meta))) {
    stop("cohort_col not found in meta: ", cohort_col)
  }
  meta[[cohort_col]] <- as.character(meta[[cohort_col]])
  if (anyNA(meta[[cohort_col]]) || any(!nzchar(meta[[cohort_col]]))) {
    stop("cohort_col contains missing/empty cohort labels.")
  }
  meta
}

.sot_row_pseudocount <- function(x, pseudocount = NULL) {
  if (!is.null(pseudocount)) return(rep(as.numeric(pseudocount), nrow(x)))
  pc <- 1e-6 * rowSums(x, na.rm = TRUE)
  pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  pc
}

.sot_log_proportions <- function(x, pseudocount = NULL) {
  pc <- .sot_row_pseudocount(x, pseudocount)
  x_pos <- sweep(x, 1L, pc, `+`)
  denom <- rowSums(x_pos, na.rm = TRUE)
  denom[!is.finite(denom) | denom <= 0] <- NA_real_
  sweep(log(x_pos), 1L, log(denom), `-`)
}

.sot_col_median <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colMedians(as.matrix(x), na.rm = TRUE)
  } else {
    apply(x, 2L, stats::median, na.rm = TRUE)
  }
}

.sot_select_pivot_features <- function(x,
                                       min_pivot_features = 12L,
                                       max_pivot_features = 80L,
                                       min_detection_rate = 0.30,
                                       trim_lower = 0.05,
                                       trim_upper = 0.10,
                                       exclude_features = c("hsa-miR-451a",
                                                            "hsa-miR-16-5p",
                                                            "hsa-miR-486-5p",
                                                            "hsa-miR-144-3p",
                                                            "hsa-miR-223-3p")) {
  det <- colMeans(x > 0, na.rm = TRUE)
  abundance <- .sot_col_median(x)
  names(abundance) <- colnames(x)
  candidates <- names(det)[is.finite(det) & det >= min_detection_rate]
  candidates <- setdiff(candidates, exclude_features)
  if (length(candidates) < min_pivot_features) {
    candidates <- setdiff(names(sort(det, decreasing = TRUE)), exclude_features)
  }
  candidates <- candidates[is.finite(abundance[candidates])]
  if (length(candidates) == 0L) candidates <- colnames(x)
  ord <- candidates[order(abundance[candidates], decreasing = TRUE)]
  if (length(ord) >= min_pivot_features) {
    lo <- floor(trim_upper * length(ord)) + 1L
    hi <- length(ord) - floor(trim_lower * length(ord))
    ord <- ord[seq.int(lo, hi)]
  }
  head(ord, min(length(ord), max(as.integer(min_pivot_features),
                                as.integer(max_pivot_features))))
}

.sot_clr_transform <- function(x, pivot_features, pseudocount = NULL) {
  logp <- .sot_log_proportions(x, pseudocount = pseudocount)
  piv <- intersect(pivot_features, colnames(logp))
  if (length(piv) == 0L) piv <- colnames(logp)
  ctr <- rowMeans(logp[, piv, drop = FALSE], na.rm = TRUE)
  out <- sweep(logp, 1L, ctr, `-`)
  out[!is.finite(out)] <- NA_real_
  out
}

.sot_fit_scaler <- function(z) {
  center <- .sot_col_median(z)
  scale <- apply(z, 2L, stats::mad, na.rm = TRUE)
  scale[!is.finite(scale) | scale < 1e-6] <- 1
  names(center) <- names(scale) <- colnames(z)
  list(center = center, scale = scale)
}

.sot_apply_scaler <- function(z, scaler) {
  common <- intersect(colnames(z), names(scaler$center))
  out <- z[, common, drop = FALSE]
  out <- sweep(out, 2L, scaler$center[common], `-`)
  out <- sweep(out, 2L, scaler$scale[common], `/`)
  out[!is.finite(out)] <- 0
  out
}

.sot_inverse_scaler <- function(z, scaler, features) {
  common <- intersect(features, colnames(z))
  out <- z[, common, drop = FALSE]
  out <- sweep(out, 2L, scaler$scale[common], `*`)
  out <- sweep(out, 2L, scaler$center[common], `+`)
  out[, features, drop = FALSE]
}

.sot_fast_auc <- function(y, score) {
  y <- as.integer(y)
  s <- as.numeric(score)
  ok <- is.finite(s) & !is.na(y)
  y <- y[ok]
  s <- s[ok]
  if (length(unique(y)) != 2L || length(unique(s)) < 2L) return(NA_real_)
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos < 1L || n_neg < 1L) return(NA_real_)
  ranks <- rank(s, ties.method = "average")
  as.numeric((sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
}

.sot_hm_var <- function(auc, n_pos, n_neg) {
  auc <- as.numeric(auc)
  n_pos <- as.numeric(n_pos)
  n_neg <- as.numeric(n_neg)
  if (!is.finite(auc) || !is.finite(n_pos) || !is.finite(n_neg) ||
      n_pos < 1 || n_neg < 1) return(NA_real_)
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  v <- (auc * (1 - auc) +
          (n_pos - 1) * (q1 - auc^2) +
          (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg)
  if (!is.finite(v) || v < 0) NA_real_ else max(v, 1e-10)
}

.sot_ci <- function(auc, n_pos, n_neg) {
  v <- .sot_hm_var(auc, n_pos, n_neg)
  if (!is.finite(v)) return(c(NA_real_, NA_real_))
  c(max(0, auc - 1.96 * sqrt(v)), min(1, auc + 1.96 * sqrt(v)))
}

.sot_weighted_panel_sum <- function(z, panel_features, feature_weights = NULL) {
  panel <- intersect(as.character(panel_features), colnames(z))
  if (length(panel) == 0L) return(rep(NA_real_, nrow(z)))
  if (is.null(feature_weights)) {
    w <- stats::setNames(rep(1, length(panel)), panel)
  } else {
    w <- as.numeric(feature_weights[panel])
    w[!is.finite(w)] <- 1
    w <- stats::setNames(w, panel)
  }
  as.numeric(z[, panel, drop = FALSE] %*% w[panel])
}

.sot_select_panel_from_clr <- function(z_train, y_train, panel_size = 20L,
                                      min_detection_rate = 0.30) {
  y_train <- as.integer(y_train)
  det <- colMeans(is.finite(z_train) & abs(z_train) > 1e-12, na.rm = TRUE)
  candidates <- names(det)[det >= min_detection_rate]
  if (length(candidates) < panel_size) {
    candidates <- names(sort(det, decreasing = TRUE))[seq_len(min(ncol(z_train),
                                                                  panel_size * 4L))]
  }
  aucs <- vapply(candidates, function(feat) .sot_fast_auc(y_train, z_train[, feat]),
                 numeric(1L))
  strength <- pmax(aucs, 1 - aucs, na.rm = TRUE)
  strength[!is.finite(strength)] <- 0.5
  panel <- candidates[order(strength, decreasing = TRUE)]
  panel <- head(panel, min(length(panel), panel_size))
  weights <- vapply(panel, function(feat) {
    med1 <- stats::median(z_train[y_train == 1L, feat], na.rm = TRUE)
    med0 <- stats::median(z_train[y_train == 0L, feat], na.rm = TRUE)
    if (is.finite(med1) && is.finite(med0) && med1 < med0) -1 else 1
  }, numeric(1L))
  stats::setNames(weights, panel)
}

.sot_sample_by_cohort <- function(z_scaled, meta, cohort_col, max_per_cohort, seed) {
  set.seed(seed)
  cohorts <- sort(unique(meta[[cohort_col]]))
  out <- vector("list", length(cohorts))
  names(out) <- cohorts
  for (coh in cohorts) {
    idx <- which(meta[[cohort_col]] == coh)
    if (length(idx) > max_per_cohort) idx <- sort(sample(idx, max_per_cohort))
    cloud <- z_scaled[idx, , drop = FALSE]
    keep <- stats::complete.cases(cloud)
    out[[coh]] <- cloud[keep, , drop = FALSE]
  }
  out[vapply(out, nrow, integer(1L)) > 0L]
}

.sot_cost <- function(a, b) {
  a2 <- rowSums(a * a)
  b2 <- rowSums(b * b)
  cost <- outer(a2, b2, `+`) - 2 * tcrossprod(a, b)
  cost[cost < 0 & is.finite(cost)] <- 0
  finite <- cost[is.finite(cost) & cost > 0]
  med <- if (length(finite)) stats::median(finite) else 1
  if (!is.finite(med) || med <= 0) med <- 1
  cost / med
}

.sot_sinkhorn_plan_r <- function(a, b, cost, epsilon, max_iter = 250L,
                                 tol = 1e-7) {
  a <- as.numeric(a); a <- a / sum(a)
  b <- as.numeric(b); b <- b / sum(b)
  eps <- max(as.numeric(epsilon), 1e-6)
  k <- exp(-cost / eps)
  k[!is.finite(k) | k <= 0] <- 1e-300
  u <- rep(1, length(a))
  v <- rep(1, length(b))
  for (iter in seq_len(max_iter)) {
    old <- u
    kv <- as.numeric(k %*% v)
    kv[kv <= 0] <- 1e-300
    u <- a / kv
    ktu <- as.numeric(t(k) %*% u)
    ktu[ktu <= 0] <- 1e-300
    v <- b / ktu
    if (max(abs(u - old), na.rm = TRUE) < tol) break
  }
  plan <- (u * k) * rep(v, each = nrow(k))
  plan[!is.finite(plan)] <- 0
  s <- sum(plan)
  if (!is.finite(s) || s <= 0) outer(a, b) else plan / s
}

.sot_free_support_barycenter_r <- function(clouds, n_atoms, epsilon,
                                           reg_lambda, barycenter_iter,
                                           max_sinkhorn_iter, tol, seed) {
  set.seed(seed)
  clouds <- clouds[vapply(clouds, nrow, integer(1L)) > 0L]
  pooled <- do.call(rbind, clouds)
  n_atoms <- min(as.integer(n_atoms), nrow(pooled))
  idx <- sample(seq_len(nrow(pooled)), n_atoms, replace = nrow(pooled) < n_atoms)
  atoms <- pooled[idx, , drop = FALSE]
  atom_w <- rep(1 / nrow(atoms), nrow(atoms))
  cohort_w <- rep(1 / length(clouds), length(clouds))
  damping <- min(max(as.numeric(reg_lambda), 0), 1)
  if (damping <= 0) damping <- 0.5
  for (iter in seq_len(max(1L, barycenter_iter))) {
    mapped_sum <- matrix(0, nrow = nrow(atoms), ncol = ncol(atoms),
                         dimnames = dimnames(atoms))
    for (j in seq_along(clouds)) {
      cloud <- clouds[[j]]
      target_w <- rep(1 / nrow(cloud), nrow(cloud))
      plan <- .sot_sinkhorn_plan_r(atom_w, target_w, .sot_cost(atoms, cloud),
                                   epsilon, max_sinkhorn_iter, tol)
      row_mass <- rowSums(plan)
      row_mass[row_mass <= 0] <- 1
      mapped <- plan %*% cloud
      mapped <- sweep(mapped, 1L, row_mass, `/`)
      mapped_sum <- mapped_sum + cohort_w[j] * mapped
    }
    atoms <- (1 - damping) * atoms + damping * mapped_sum
  }
  list(atoms = atoms, weights = atom_w, backend = "R_internal_sinkhorn",
       pot_available = FALSE, n_atoms = nrow(atoms), n_cohorts = length(clouds))
}

.sot_python_module <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  py_path <- system.file("python", "sinkhorn_ot_scorer.py",
                         package = "OmicSelector")
  if (!file.exists(py_path)) return(NULL)
  tryCatch(reticulate::import_from_path("sinkhorn_ot_scorer",
                                        path = dirname(py_path),
                                        convert = TRUE),
           error = function(e) NULL)
}

.sot_build_barycenter <- function(clouds, n_atoms, epsilon, reg_lambda,
                                  barycenter_iter, max_sinkhorn_iter,
                                  tol, seed, backend = c("auto", "python", "R")) {
  backend <- match.arg(backend)
  if (backend %in% c("auto", "python")) {
    mod <- .sot_python_module()
    if (!is.null(mod)) {
      res <- tryCatch(mod$free_support_sinkhorn_barycenter(
        cohort_clouds = clouds,
        cohort_weights = NULL,
        n_atoms = as.integer(n_atoms),
        epsilon = as.numeric(epsilon),
        reg_lambda = as.numeric(reg_lambda),
        barycenter_iter = as.integer(barycenter_iter),
        max_sinkhorn_iter = as.integer(max_sinkhorn_iter),
        tol = as.numeric(tol),
        seed = as.integer(seed),
        use_pot = TRUE
      ), error = function(e) e)
      if (!inherits(res, "error")) return(res)
      if (identical(backend, "python")) stop(res$message)
    }
  }
  .sot_free_support_barycenter_r(clouds, n_atoms, epsilon, reg_lambda,
                                 barycenter_iter, max_sinkhorn_iter, tol, seed)
}

.sot_project_scaled <- function(z_scaled, fit) {
  backend <- fit$backend_request %||% "auto"
  if (backend %in% c("auto", "python")) {
    mod <- .sot_python_module()
    if (!is.null(mod)) {
      res <- tryCatch(mod$project_to_barycenter(
        source_cloud = z_scaled,
        barycenter_atoms = fit$barycenter_atoms_scaled,
        barycenter_weights = fit$barycenter_weights,
        epsilon = as.numeric(fit$epsilon),
        max_sinkhorn_iter = as.integer(fit$max_sinkhorn_iter),
        tol = as.numeric(fit$sinkhorn_tol),
        use_pot = TRUE,
        single_sample_kernel_fallback = TRUE
      ), error = function(e) e)
      if (!inherits(res, "error")) {
        out <- as.matrix(res$projection)
        attr(out, "projection_backend") <- res$backend
        return(out)
      }
      if (identical(backend, "python")) stop(res$message)
    }
  }
  atoms <- fit$barycenter_atoms_scaled
  b <- fit$barycenter_weights / sum(fit$barycenter_weights)
  cost <- .sot_cost(z_scaled, atoms)
  if (nrow(z_scaled) == 1L) {
    w <- as.numeric(b * exp(-cost[1L, ] / max(fit$epsilon, 1e-6)))
    if (!is.finite(sum(w)) || sum(w) <= 0) w <- b
    w <- w / sum(w)
    out <- matrix(as.numeric(w %*% atoms), nrow = 1L,
                  dimnames = list(rownames(z_scaled), colnames(atoms)))
    attr(out, "projection_backend") <- "R_single_sample_entropic_kernel"
    return(out)
  }
  a <- rep(1 / nrow(z_scaled), nrow(z_scaled))
  plan <- .sot_sinkhorn_plan_r(a, b, cost, fit$epsilon,
                               fit$max_sinkhorn_iter, fit$sinkhorn_tol)
  row_mass <- rowSums(plan)
  row_mass[row_mass <= 0] <- 1
  out <- plan %*% atoms
  out <- sweep(out, 1L, row_mass, `/`)
  rownames(out) <- rownames(z_scaled)
  attr(out, "projection_backend") <- "R_internal_sinkhorn"
  out
}

.sot_fit_core <- function(X_train, y_train, meta_train, cohort_col,
                          epsilon, reg_lambda, n_barycenter_atoms,
                          max_cohort_atoms, barycenter_iter,
                          max_sinkhorn_iter, sinkhorn_tol,
                          min_pivot_features, max_pivot_features,
                          panel_size, panel_features, feature_weights,
                          pseudocount, seed, backend) {
  pivots <- .sot_select_pivot_features(
    X_train,
    min_pivot_features = min_pivot_features,
    max_pivot_features = max_pivot_features)
  z_train <- .sot_clr_transform(X_train, pivots, pseudocount = pseudocount)
  scaler <- .sot_fit_scaler(z_train)
  z_scaled <- .sot_apply_scaler(z_train, scaler)
  if (is.null(panel_features)) {
    if (is.null(y_train)) stop("y_train is required when panel_features is NULL.")
    weights <- .sot_select_panel_from_clr(z_train, y_train, panel_size = panel_size)
    panel_features <- names(weights)
    feature_weights <- feature_weights %||% weights
  }
  clouds <- .sot_sample_by_cohort(z_scaled, meta_train, cohort_col,
                                  max_cohort_atoms, seed)
  if (length(clouds) < 1L) stop("No finite cohort clouds available for OT fit.")
  bary <- .sot_build_barycenter(clouds, n_barycenter_atoms, epsilon,
                                reg_lambda, barycenter_iter,
                                max_sinkhorn_iter, sinkhorn_tol, seed, backend)
  atoms <- as.matrix(bary$atoms)
  colnames(atoms) <- colnames(z_scaled)
  structure(list(
    method = "sinkhorn_barycenter_transport",
    features = colnames(X_train),
    pivot_features = pivots,
    scaler = scaler,
    panel_features = intersect(panel_features, colnames(X_train)),
    feature_weights = feature_weights[intersect(panel_features, names(feature_weights))],
    cohort_col = cohort_col,
    train_cohorts = sort(unique(meta_train[[cohort_col]])),
    epsilon = as.numeric(epsilon),
    lambda = as.numeric(reg_lambda),
    n_barycenter_atoms_requested = as.integer(n_barycenter_atoms),
    n_barycenter_atoms = nrow(atoms),
    max_cohort_atoms = as.integer(max_cohort_atoms),
    barycenter_iter = as.integer(barycenter_iter),
    max_sinkhorn_iter = as.integer(max_sinkhorn_iter),
    sinkhorn_tol = as.numeric(sinkhorn_tol),
    backend = bary$backend %||% "unknown",
    backend_request = backend,
    pot_available = isTRUE(bary$pot_available),
    barycenter_atoms_scaled = atoms,
    barycenter_weights = as.numeric(bary$weights),
    n_train_samples = nrow(X_train),
    n_train_cohorts = length(unique(meta_train[[cohort_col]])),
    fit_notes = "Train-only pivot CLR; free-support entropic Sinkhorn barycenter; batchwise barycentric projection."
  ), class = "sinkhorn_ot_scorer")
}

.sot_tune_train_only <- function(X_train, y_train, meta_train, cohort_col,
                                 epsilon_grid, lambda_grid, panel_size,
                                 tune_max_folds, tune_n_atoms,
                                 tune_max_cohort_atoms, seed, ...) {
  cohorts <- sort(unique(meta_train[[cohort_col]]))
  if (length(cohorts) < 3L || length(unique(y_train)) < 2L) {
    return(list(epsilon = epsilon_grid[1L], lambda = lambda_grid[1L],
                trace = data.frame(status = "not_tuned_insufficient_train_cohorts")))
  }
  set.seed(seed)
  val_cohorts <- head(sample(cohorts), min(length(cohorts), tune_max_folds))
  rows <- list()
  for (eps in epsilon_grid) for (lam in lambda_grid) for (val in val_cohorts) {
    tr <- meta_train[[cohort_col]] != val
    va <- meta_train[[cohort_col]] == val
    if (length(unique(y_train[tr])) < 2L || length(unique(y_train[va])) < 2L) next
    fit <- tryCatch(.sot_fit_core(
      X_train[tr, , drop = FALSE], y_train[tr],
      meta_train[tr, , drop = FALSE], cohort_col,
      epsilon = eps, reg_lambda = lam,
      n_barycenter_atoms = tune_n_atoms,
      max_cohort_atoms = tune_max_cohort_atoms,
      panel_size = panel_size, panel_features = NULL, feature_weights = NULL,
      seed = seed + length(rows) + 1L, ...
    ), error = function(e) e)
    if (inherits(fit, "error")) {
      rows[[length(rows) + 1L]] <- data.frame(
        epsilon = eps, lambda = lam, validation_cohort = val,
        auc = NA_real_, status = paste0("fit_failed:", fit$message))
      next
    }
    sc <- tryCatch(score_sinkhorn_ot_scorer(
      X_train[va, , drop = FALSE], fit,
      panel_features = fit$panel_features,
      feature_weights = fit$feature_weights),
      error = function(e) structure(rep(NA_real_, sum(va)), error = e$message))
    rows[[length(rows) + 1L]] <- data.frame(
      epsilon = eps, lambda = lam, validation_cohort = val,
      auc = .sot_fast_auc(y_train[va], sc),
      status = if (is.null(attr(sc, "error"))) "evaluable" else attr(sc, "error"))
  }
  trace <- if (length(rows)) do.call(rbind, rows) else {
    data.frame(status = "not_tuned_no_evaluable_inner_folds")
  }
  ok <- trace[is.finite(trace$auc), , drop = FALSE]
  if (nrow(ok) == 0L) {
    return(list(epsilon = epsilon_grid[1L], lambda = lambda_grid[1L],
                trace = trace))
  }
  agg <- aggregate(auc ~ epsilon + lambda, ok, mean)
  best <- agg[order(-agg$auc, agg$epsilon, agg$lambda), ][1L, ]
  list(epsilon = best$epsilon, lambda = best$lambda, trace = trace)
}

#' Fit Sinkhorn OT cohort barycenter scorer.
#'
#' @param X_train Samples x features non-negative matrix.
#' @param y_train Optional 0/1 outcome vector. Required when panel selection or
#'   hyperparameter tuning is requested.
#' @param meta_train data.frame with one row per sample and a cohort column.
#' @param cohort_col Column identifying training cohorts.
#' @param panel_size Number of features selected when `panel_features` is NULL.
#' @param panel_features Optional fixed panel features.
#' @param feature_weights Optional named numeric panel weights.
#' @param epsilon Entropic-transport regularization. If NULL, selected from
#'   `epsilon_grid` when tuning is enabled.
#' @param lambda Barycenter shrinkage weight. If NULL, selected from
#'   `lambda_grid` when tuning is enabled.
#' @param epsilon_grid,lambda_grid Candidate hyperparameter grids.
#' @param tune Logical. If TRUE, use train-only nested cohort CV to tune
#'   transport hyperparameters.
#' @param tune_max_folds Maximum training folds used during tuning.
#' @param tune_n_atoms Number of barycenter atoms used during tuning.
#' @param tune_max_cohort_atoms Maximum cohort atoms used during tuning.
#' @param n_barycenter_atoms Number of free-support barycenter atoms.
#' @param max_cohort_atoms Maximum samples retained per cohort cloud.
#' @param barycenter_iter Number of barycenter updates.
#' @param max_sinkhorn_iter Maximum Sinkhorn iterations.
#' @param sinkhorn_tol Sinkhorn convergence tolerance.
#' @param min_pivot_features,max_pivot_features Bounds for CLR pivot features.
#' @param pseudocount Optional additive pseudocount before CLR transform.
#' @param seed Integer random seed.
#' @param backend Backend for barycenter/projection computation.
#' @return A frozen `sinkhorn_ot_scorer` object.
fit_sinkhorn_ot_scorer <- function(X_train,
                                   y_train = NULL,
                                   meta_train,
                                   cohort_col = "cohort",
                                   panel_size = 20L,
                                   panel_features = NULL,
                                   feature_weights = NULL,
                                   epsilon = NULL,
                                   lambda = NULL,
                                   epsilon_grid = c(0.05, 0.10, 0.20),
                                   lambda_grid = c(0.35, 0.60, 0.85),
                                   tune = TRUE,
                                   tune_max_folds = 3L,
                                   tune_n_atoms = 120L,
                                   tune_max_cohort_atoms = 60L,
                                   n_barycenter_atoms = 1000L,
                                   max_cohort_atoms = 250L,
                                   barycenter_iter = 4L,
                                   max_sinkhorn_iter = 250L,
                                   sinkhorn_tol = 1e-7,
                                   min_pivot_features = 12L,
                                   max_pivot_features = 80L,
                                   pseudocount = NULL,
                                   seed = 42L,
                                   backend = c("auto", "python", "R")) {
  backend <- match.arg(backend)
  X_train <- .sot_matrix(X_train, "X_train")
  meta_train <- .sot_meta(meta_train, nrow(X_train), cohort_col)
  if (!is.null(y_train)) {
    y_train <- as.integer(y_train)
    if (length(y_train) != nrow(X_train)) stop("y_train length must match nrow(X_train).")
  }
  if (isTRUE(tune) && (is.null(epsilon) || is.null(lambda))) {
    tuned <- .sot_tune_train_only(
      X_train, y_train, meta_train, cohort_col,
      epsilon_grid = epsilon_grid, lambda_grid = lambda_grid,
      panel_size = panel_size, tune_max_folds = tune_max_folds,
      tune_n_atoms = tune_n_atoms,
      tune_max_cohort_atoms = tune_max_cohort_atoms,
      seed = seed, backend = backend, barycenter_iter = max(1L, min(2L, barycenter_iter)),
      max_sinkhorn_iter = max_sinkhorn_iter, sinkhorn_tol = sinkhorn_tol,
      min_pivot_features = min_pivot_features,
      max_pivot_features = max_pivot_features,
      pseudocount = pseudocount)
    epsilon <- epsilon %||% tuned$epsilon
    lambda <- lambda %||% tuned$lambda
  } else {
    tuned <- list(trace = data.frame(status = "not_tuned_user_fixed_hyperparameters"))
    epsilon <- epsilon %||% epsilon_grid[1L]
    lambda <- lambda %||% lambda_grid[1L]
  }
  fit <- .sot_fit_core(
    X_train, y_train, meta_train, cohort_col,
    epsilon = epsilon, reg_lambda = lambda,
    n_barycenter_atoms = n_barycenter_atoms,
    max_cohort_atoms = max_cohort_atoms,
    barycenter_iter = barycenter_iter,
    max_sinkhorn_iter = max_sinkhorn_iter,
    sinkhorn_tol = sinkhorn_tol,
    min_pivot_features = min_pivot_features,
    max_pivot_features = max_pivot_features,
    panel_size = panel_size,
    panel_features = panel_features,
    feature_weights = feature_weights,
    pseudocount = pseudocount,
    seed = seed,
    backend = backend)
  fit$tuning_trace <- tuned$trace
  fit$tuning_protocol <- if (isTRUE(tune)) {
    "train_only_nested_cohort_cv"
  } else {
    "fixed_hyperparameters"
  }
  fit
}

#' Apply a frozen Sinkhorn OT scorer.
#'
#' @param X_test Samples x features matrix in the fitted feature space.
#' @param fit Object from `fit_sinkhorn_ot_scorer()`.
#' @param meta_test Optional metadata retained for a common scorer signature.
#' @return Projected CLR matrix with samples in rows and fit features in columns.
apply_sinkhorn_ot_scorer <- function(X_test, fit, meta_test = NULL) {
  if (!inherits(fit, "sinkhorn_ot_scorer")) {
    stop("fit must inherit from sinkhorn_ot_scorer.")
  }
  X_test <- .sot_matrix(X_test, "X_test")
  missing <- setdiff(fit$features, colnames(X_test))
  if (length(missing) > 0L) {
    stop("X_test is missing fitted features: ", paste(head(missing, 10), collapse = ", "))
  }
  X_test <- X_test[, fit$features, drop = FALSE]
  z <- .sot_clr_transform(X_test, fit$pivot_features)
  z_scaled <- .sot_apply_scaler(z, fit$scaler)
  proj_scaled <- .sot_project_scaled(z_scaled, fit)
  colnames(proj_scaled) <- fit$features
  out <- .sot_inverse_scaler(proj_scaled, fit$scaler, fit$features)
  attr(out, "projection_backend") <- attr(proj_scaled, "projection_backend")
  attr(out, "single_sample_projection") <- if (nrow(X_test) == 1L) {
    "entropic_kernel_fallback"
  } else {
    "batchwise_sinkhorn_plan"
  }
  out
}

#' Score projected samples by summing panel features.
#'
#' @param X_test Samples x features matrix in the fitted feature space.
#' @param fit Object from `fit_sinkhorn_ot_scorer()`.
#' @param panel_features Panel features to sum after projection.
#' @param feature_weights Optional named numeric panel weights.
#' @param meta_test Optional metadata retained for a common scorer signature.
score_sinkhorn_ot_scorer <- function(X_test, fit, panel_features = fit$panel_features,
                                     feature_weights = fit$feature_weights,
                                     meta_test = NULL) {
  proj <- apply_sinkhorn_ot_scorer(X_test, fit, meta_test = meta_test)
  sc <- .sot_weighted_panel_sum(proj, panel_features, feature_weights)
  attr(sc, "projection_backend") <- attr(proj, "projection_backend")
  attr(sc, "single_sample_projection") <- attr(proj, "single_sample_projection")
  sc
}

# ---- Reviewer-package S55 benchmark -----------------------------------------

.sot_default_paths <- function(output_dir = tempdir()) {
  output_dir <- path.expand(as.character(output_dir)[1L])
  list(
    table = file.path(output_dir, "table_S55_sinkhorn_ot_benchmarks.tsv"),
    png = file.path(output_dir, "figure_S55_sinkhorn_ot.png"),
    pdf = file.path(output_dir, "figure_S55_sinkhorn_ot.pdf"),
    caption = file.path(output_dir, "figure_S55_sinkhorn_ot.md"),
    log = file.path(
      output_dir,
      paste0("sinkhorn_ot_S55_", format(Sys.time(), "%Y%m%dT%H%M%SZ",
                                        tz = "UTC"), ".log")
    )
  )
}

.sot_paths <- .sot_default_paths()

.sot_config <- list(
  seed = 42L,
  n_atoms = .sot_env_int("SINKHORN_OT_N_ATOMS", 180L),
  max_cohort_atoms = .sot_env_int("SINKHORN_OT_MAX_COHORT_ATOMS", 80L),
  barycenter_iter = .sot_env_int("SINKHORN_OT_BARYCENTER_ITER", 3L),
  max_sinkhorn_iter = .sot_env_int("SINKHORN_OT_MAX_ITER", 180L),
  tune = .sot_env_logical("SINKHORN_OT_TUNE", TRUE),
  epsilon_grid = .sot_env_num_vec("SINKHORN_OT_EPSILON_GRID", c(0.08, 0.16)),
  lambda_grid = .sot_env_num_vec("SINKHORN_OT_LAMBDA_GRID", c(0.4, 0.7)),
  tune_max_folds = .sot_env_int("SINKHORN_OT_TUNE_MAX_FOLDS", 2L),
  tune_n_atoms = .sot_env_int("SINKHORN_OT_TUNE_N_ATOMS", 80L),
  tune_max_cohort_atoms = .sot_env_int("SINKHORN_OT_TUNE_MAX_COHORT_ATOMS", 40L),
  panel_size = .sot_env_int("SINKHORN_OT_PANEL_SIZE", 20L),
  min_train_cohorts = .sot_env_int("SINKHORN_OT_MIN_TRAIN_COHORTS", 2L),
  min_common_features = .sot_env_int("SINKHORN_OT_MIN_COMMON_FEATURES", 40L),
  boot_B = .sot_env_int("SINKHORN_OT_BOOT_B", 500L),
  limit_cohorts = .sot_env_int("SINKHORN_OT_LIMIT_COHORTS", 999L),
  backend = Sys.getenv("SINKHORN_OT_BACKEND", "auto")
)

.sot_primary_cohorts <- function(s5_path = NULL, s22_path = NULL) {
  s5_path <- .singlesample_extdata_path("table_S5_null_calibration.tsv", s5_path,
                                  "S5 null-calibration table")
  s22_path <- .singlesample_extdata_path("table_S22_delong_paired_auc.tsv", s22_path,
                                   "S22 DeLong paired-AUC table")
  s5 <- data.table::fread(s5_path)
  s22 <- data.table::fread(s22_path)
  c5 <- sort(unique(s5$Accession))
  c22 <- sort(unique(s22[status == "evaluable", cohort]))
  if (!identical(c5, c22)) stop("Primary cohort set mismatch between S5 and S22.")
  head(c5, min(length(c5), .sot_config$limit_cohorts))
}

.sot_technology <- function(x) {
  raw <- as.character(x)
  out <- rep("Other", length(raw))
  out[grepl("NGS|miRNA-seq", raw, ignore.case = TRUE)] <- "NGS"
  out[grepl("Toray|3D", raw, ignore.case = TRUE)] <- "Toray"
  out
}

.sot_load_metadata <- function(primary, kit_metadata_path = NULL,
                               union_inventory_path = NULL) {
  kit_metadata_path <- .singlesample_extdata_path(
    "table_S43_matched_null_kit_metadata.tsv", kit_metadata_path,
    "S43 kit metadata")
  union_inventory_path <- .singlesample_extdata_path(
    "union_inventory.tsv", union_inventory_path, "union inventory")
  meta <- data.table::fread(
    kit_metadata_path, sep = "\t", quote = "", fill = TRUE)
  meta <- meta[cohort %in% primary]
  inv <- data.table::fread(union_inventory_path, sep = "\t",
                           quote = "", fill = TRUE)
  inv <- inv[accession %in% primary,
             .(cohort = accession,
               cancer_type = cancer_type_canonical %||% cancer_type)]
  meta <- merge(meta, inv, by = "cohort", all.x = TRUE, sort = FALSE)
  meta[, technology := .sot_technology(modality)]
  meta[]
}

.sot_source_prepare_env <- function() {
  stop("The reviewer-package S55 runner is not available from the installed package; ",
       "use the manuscript repository analysis script instead.", call. = FALSE)
}

.sot_collapse_duplicate_features <- function(x) {
  nms <- colnames(x)
  if (!anyDuplicated(nms)) return(x)
  groups <- split(seq_along(nms), nms)
  out <- matrix(NA_real_, nrow = nrow(x), ncol = length(groups),
                dimnames = list(rownames(x), names(groups)))
  for (i in seq_along(groups)) {
    idx <- groups[[i]]
    out[, i] <- if (length(idx) == 1L) x[, idx] else {
      rowMeans(x[, idx, drop = FALSE], na.rm = TRUE)
    }
  }
  out
}

.sot_harmonise_feature_names <- function(x) {
  X <- as.matrix(x)
  storage.mode(X) <- "double"
  if (exists("apply_mirna_aliases", mode = "function")) {
    X <- tryCatch(apply_mirna_aliases(
      X, target_namespace = "mirna_name", keep_unresolved = TRUE,
      verbose = FALSE),
      error = function(e) X)
  }
  colnames(X) <- sub("^hsa-mir-", "hsa-miR-", colnames(X), ignore.case = TRUE)
  colnames(X) <- sub("^hsa-let-", "hsa-let-", colnames(X), ignore.case = TRUE)
  .sot_collapse_duplicate_features(X)
}

.sot_load_data <- function(primary, meta, log_fun = message,
                           union_inventory_path = NULL,
                           sample_filters_path = NULL,
                           matched_set_audit_path = NULL) {
  prep_env <- .sot_source_prepare_env()
  union_inventory_path <- .singlesample_extdata_path(
    "union_inventory.tsv", union_inventory_path, "union inventory")
  sample_filters_path <- .singlesample_extdata_path(
    "v0_6_sample_filters.tsv", sample_filters_path, "v0.6 sample filters")
  matched_set_audit_path <- .singlesample_extdata_path(
    "matched_set_audit.tsv", matched_set_audit_path, "matched-set audit")
  inv <- data.table::fread(union_inventory_path, sep = "\t",
                           quote = "", fill = TRUE)
  rows <- inv[accession %in% primary]
  rows <- rows[match(primary, accession)]
  sample_filters <- data.table::fread(sample_filters_path,
                                      sep = "\t", quote = "",
                                      colClasses = "character", fill = TRUE)
  get_fun <- function(name) {
    if (exists(name, envir = prep_env, inherits = FALSE)) {
      get(name, envir = prep_env, inherits = FALSE)
    } else {
      get(name, envir = globalenv(), inherits = TRUE)
    }
  }
  prepare <- get_fun(".prepare_cohort_data")
  outcome_dict <- get_fun("load_outcome_dictionary_v0_4")(strict = TRUE)
  matched_set_audit <- get_fun("load_matched_set_audit")(
    matched_set_audit_path, n_tier_a = 35L)
  out <- list()
  failures <- list()
  for (i in seq_len(nrow(rows))) {
    acc <- rows$accession[i]
    log_fun(sprintf("load %02d/%02d %s", i, nrow(rows), acc))
    dat <- tryCatch(prepare(rows[i], outcome_dict, sample_filters, matched_set_audit),
      error = function(e) e)
    if (inherits(dat, "error")) {
      failures[[acc]] <- dat$message
      next
    }
    dat$X <- .sot_harmonise_feature_names(dat$X)
    cm <- meta[cohort == acc][1L]
    dat$sample_meta <- data.frame(
      sample_id = rownames(dat$X) %||% paste0(acc, "_", seq_len(nrow(dat$X))),
      cohort = acc,
      technology = cm$technology %||% "Other",
      modality = cm$modality %||% NA_character_,
      cancer_type = cm$cancer_type %||% NA_character_,
      biofluid = cm$biofluid %||% NA_character_,
      library_kit_family = cm$library_kit %||% NA_character_,
      extraction_kit_family = cm$extraction_kit_family %||% NA_character_,
      provenance_block = cm$provenance_block %||% acc,
      stringsAsFactors = FALSE)
    dat$technology <- cm$technology %||% "Other"
    dat$modality <- cm$modality %||% NA_character_
    dat$cancer_type <- cm$cancer_type %||% NA_character_
    dat$provenance_block <- cm$provenance_block %||% acc
    out[[acc]] <- dat
  }
  attr(out, "failures") <- failures
  out
}

.sot_common_features <- function(data, accs) {
  Reduce(intersect, lapply(accs, function(a) colnames(data[[a]]$X)))
}

.sot_select_overlap_subset <- function(data, train_accs, held,
                                       min_train_cohorts, min_common_features) {
  selected <- train_accs
  common <- if (length(selected)) .sot_common_features(data, c(selected, held)) else character(0)
  excluded <- character(0)
  while (length(common) < min_common_features &&
         length(selected) > min_train_cohorts) {
    trial <- lapply(selected, function(drop) {
      keep <- setdiff(selected, drop)
      feat <- .sot_common_features(data, c(keep, held))
      list(drop = drop, keep = keep, common = feat, n = length(feat))
    })
    n_best <- vapply(trial, `[[`, integer(1L), "n")
    best <- trial[[which.max(n_best)]]
    if (best$n <= length(common)) break
    excluded <- c(excluded, best$drop)
    selected <- best$keep
    common <- best$common
  }
  list(train_accs = selected, common = common, excluded = excluded)
}

.sot_pooled <- function(data, accs, features) {
  X <- do.call(rbind, lapply(accs, function(a) data[[a]]$X[, features, drop = FALSE]))
  y <- unlist(lapply(accs, function(a) as.integer(data[[a]]$y)), use.names = FALSE)
  meta <- do.call(rbind, lapply(accs, function(a) data[[a]]$sample_meta))
  list(X = X, y = y, meta = meta)
}

.sot_failure_row <- function(held, comparison_type, reason, data) {
  hd <- data[[held]]
  data.table::data.table(
    row_type = "loco",
    method = "sinkhorn_barycenter_transport",
    held_out = held,
    comparison_type = comparison_type,
    cancer_type = hd$cancer_type,
    technology = hd$technology,
    train_technologies = NA_character_,
    train_cohorts = NA_character_,
    n_train_cohorts = 0L,
    n_train_samples = 0L,
    n_test = nrow(hd$X),
    n_test_pos = sum(hd$y == 1L, na.rm = TRUE),
    n_test_neg = sum(hd$y == 0L, na.rm = TRUE),
    n_common_features = NA_integer_,
    panel_size = NA_integer_,
    epsilon = NA_real_,
    lambda = NA_real_,
    n_barycenter_atoms = NA_integer_,
    max_cohort_atoms = NA_integer_,
    sinkhorn_backend = NA_character_,
    rclr_loco_auc = NA_real_,
    sinkhorn_ot_auc = NA_real_,
    auc_ci_low = NA_real_,
    auc_ci_high = NA_real_,
    cross_tech_loco_auc = NA_real_,
    cross_tech_lift_vs_rCLR = NA_real_,
    within_tech_loco_auc = NA_real_,
    within_tech_lift_vs_rCLR = NA_real_,
    auc_var = NA_real_,
    p_auc_gt_0_5_hm = NA_real_,
    q_block_BH = NA_real_,
    provenance_block = hd$provenance_block,
    status = "not_evaluable",
    status_reason = reason,
    leakage_flag = FALSE,
    plausibility_flag = FALSE,
    notes = reason)
}

.sot_eval_one <- function(held, comparison_type, data, log_fun = message) {
  hd <- data[[held]]
  loaded <- names(data)
  held_tech <- hd$technology
  train_accs <- switch(comparison_type,
    within_tech = loaded[vapply(loaded, function(a) {
      !identical(a, held) && identical(data[[a]]$technology, held_tech)
    }, logical(1L))],
    cross_tech = {
      if (identical(held_tech, "NGS")) {
        loaded[vapply(loaded, function(a) identical(data[[a]]$technology, "Toray"),
                      logical(1L))]
      } else if (identical(held_tech, "Toray")) {
        loaded[vapply(loaded, function(a) identical(data[[a]]$technology, "NGS"),
                      logical(1L))]
      } else {
        character(0)
      }
    },
    character(0))
  if (length(train_accs) < .sot_config$min_train_cohorts) {
    return(.sot_failure_row(held, comparison_type,
                            "insufficient_train_cohorts", data))
  }
  overlap <- .sot_select_overlap_subset(data, train_accs, held,
                                        .sot_config$min_train_cohorts,
                                        .sot_config$min_common_features)
  train_accs <- overlap$train_accs
  common <- overlap$common
  if (length(common) < .sot_config$min_common_features) {
    return(.sot_failure_row(held, comparison_type,
                            paste0("insufficient_common_features_n=", length(common)),
                            data))
  }
  tr <- .sot_pooled(data, train_accs, common)
  X_test <- hd$X[, common, drop = FALSE]
  y_test <- as.integer(hd$y)
  if (length(unique(tr$y)) < 2L || length(unique(y_test)) < 2L) {
    return(.sot_failure_row(held, comparison_type,
                            "class_floor_failed", data))
  }
  log_fun(sprintf("fit held=%s comparison=%s train=%d ntrain=%d p=%d",
                  held, comparison_type, length(train_accs), nrow(tr$X), length(common)))
  fit <- tryCatch(fit_sinkhorn_ot_scorer(
    tr$X, tr$y, tr$meta,
    cohort_col = "cohort",
    panel_size = .sot_config$panel_size,
    epsilon_grid = .sot_config$epsilon_grid,
    lambda_grid = .sot_config$lambda_grid,
    tune = .sot_config$tune,
    tune_max_folds = .sot_config$tune_max_folds,
    tune_n_atoms = .sot_config$tune_n_atoms,
    tune_max_cohort_atoms = .sot_config$tune_max_cohort_atoms,
    n_barycenter_atoms = .sot_config$n_atoms,
    max_cohort_atoms = .sot_config$max_cohort_atoms,
    barycenter_iter = .sot_config$barycenter_iter,
    max_sinkhorn_iter = .sot_config$max_sinkhorn_iter,
    seed = .sot_config$seed,
    backend = .sot_config$backend),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(.sot_failure_row(held, comparison_type,
                            paste0("fit_failed:", fit$message), data))
  }
  base_z <- .sot_clr_transform(X_test, fit$pivot_features)
  base_score <- .sot_weighted_panel_sum(base_z, fit$panel_features,
                                        fit$feature_weights)
  ot_score <- tryCatch(score_sinkhorn_ot_scorer(
    X_test, fit, panel_features = fit$panel_features,
    feature_weights = fit$feature_weights),
    error = function(e) structure(rep(NA_real_, nrow(X_test)), error = e$message))
  auc_base <- .sot_fast_auc(y_test, base_score)
  auc_ot <- .sot_fast_auc(y_test, ot_score)
  n_pos <- sum(y_test == 1L, na.rm = TRUE)
  n_neg <- sum(y_test == 0L, na.rm = TRUE)
  ci <- .sot_ci(auc_ot, n_pos, n_neg)
  v <- .sot_hm_var(auc_ot, n_pos, n_neg)
  se <- if (is.finite(v)) sqrt(v) else NA_real_
  p <- if (is.finite(se) && se > 0 && is.finite(auc_ot)) {
    stats::pnorm((auc_ot - 0.5) / se, lower.tail = FALSE)
  } else {
    NA_real_
  }
  status <- if (is.finite(auc_ot) && is.finite(auc_base)) "evaluable" else "not_evaluable"
  reason <- if (identical(status, "evaluable")) "" else {
    attr(ot_score, "error") %||% "constant_or_all_NA_score"
  }
  plaus <- is.finite(auc_ot) && (auc_ot > 0.95 || auc_ot < 0.05)
  data.table::data.table(
    row_type = "loco",
    method = "sinkhorn_barycenter_transport",
    held_out = held,
    comparison_type = comparison_type,
    cancer_type = hd$cancer_type,
    technology = held_tech,
    train_technologies = paste(sort(unique(vapply(train_accs, function(a)
      data[[a]]$technology, character(1L)))), collapse = ","),
    train_cohorts = paste(train_accs, collapse = ","),
    n_train_cohorts = length(train_accs),
    n_train_samples = nrow(tr$X),
    n_test = nrow(X_test),
    n_test_pos = n_pos,
    n_test_neg = n_neg,
    n_common_features = length(common),
    panel_size = length(fit$panel_features),
    epsilon = fit$epsilon,
    lambda = fit$lambda,
    n_barycenter_atoms = fit$n_barycenter_atoms,
    max_cohort_atoms = fit$max_cohort_atoms,
    sinkhorn_backend = paste(fit$backend,
                             attr(ot_score, "projection_backend") %||% "unknown",
                             sep = "|"),
    rclr_loco_auc = auc_base,
    sinkhorn_ot_auc = auc_ot,
    auc_ci_low = ci[1L],
    auc_ci_high = ci[2L],
    cross_tech_loco_auc = if (comparison_type == "cross_tech") auc_ot else NA_real_,
    cross_tech_lift_vs_rCLR = if (comparison_type == "cross_tech") auc_ot - auc_base else NA_real_,
    within_tech_loco_auc = if (comparison_type == "within_tech") auc_ot else NA_real_,
    within_tech_lift_vs_rCLR = if (comparison_type == "within_tech") auc_ot - auc_base else NA_real_,
    auc_var = v,
    p_auc_gt_0_5_hm = p,
    q_block_BH = NA_real_,
    provenance_block = hd$provenance_block,
    status = status,
    status_reason = reason,
    leakage_flag = any(fit$train_cohorts %in% held),
    plausibility_flag = plaus,
    notes = paste0("barycenter_train_only; tune=", fit$tuning_protocol,
                   "; overlap_excluded_train_cohorts=",
                   if (length(overlap$excluded)) paste(overlap$excluded, collapse = ",") else "none",
                   "; projection=", attr(ot_score, "single_sample_projection") %||%
                     "batchwise_sinkhorn_plan"))
}

.sot_add_q <- function(dt) {
  dt[, q_block_BH := NA_real_]
  for (cmp in unique(dt$comparison_type)) {
    idx <- which(dt$comparison_type == cmp & is.finite(dt$p_auc_gt_0_5_hm))
    if (length(idx)) dt$q_block_BH[idx] <- p.adjust(dt$p_auc_gt_0_5_hm[idx], method = "BH")
  }
  dt
}

.sot_summary_rows <- function(dt) {
  ev <- dt[row_type == "loco" & status == "evaluable" & is.finite(sinkhorn_ot_auc)]
  if (nrow(ev) == 0L) return(dt[0])
  out <- ev[, {
    blocks <- unique(provenance_block)
    set.seed(.sot_config$seed)
    boot <- replicate(.sot_config$boot_B, {
      b <- sample(blocks, length(blocks), replace = TRUE)
      rows <- data.table::rbindlist(lapply(b, function(bb) .SD[provenance_block == bb]),
                                    fill = TRUE)
      mean(rows$sinkhorn_ot_auc, na.rm = TRUE)
    })
    lift_col <- if (.BY$comparison_type == "cross_tech") {
      "cross_tech_lift_vs_rCLR"
    } else {
      "within_tech_lift_vs_rCLR"
    }
    lift <- mean(get(lift_col), na.rm = TRUE)
    .(row_type = "summary",
      method = "sinkhorn_barycenter_transport",
      held_out = "__summary__",
      cancer_type = "multi-cancer",
      technology = paste(sort(unique(technology)), collapse = ","),
      train_technologies = paste(sort(unique(train_technologies)), collapse = ";"),
      train_cohorts = NA_character_,
      n_train_cohorts = NA_integer_,
      n_train_samples = sum(n_train_samples, na.rm = TRUE),
      n_test = sum(n_test, na.rm = TRUE),
      n_test_pos = sum(n_test_pos, na.rm = TRUE),
      n_test_neg = sum(n_test_neg, na.rm = TRUE),
      n_common_features = stats::median(n_common_features, na.rm = TRUE),
      panel_size = stats::median(panel_size, na.rm = TRUE),
      epsilon = stats::median(epsilon, na.rm = TRUE),
      lambda = stats::median(lambda, na.rm = TRUE),
      n_barycenter_atoms = stats::median(n_barycenter_atoms, na.rm = TRUE),
      max_cohort_atoms = stats::median(max_cohort_atoms, na.rm = TRUE),
      sinkhorn_backend = paste(sort(unique(sinkhorn_backend)), collapse = ";"),
      rclr_loco_auc = mean(rclr_loco_auc, na.rm = TRUE),
      sinkhorn_ot_auc = mean(sinkhorn_ot_auc, na.rm = TRUE),
      auc_ci_low = stats::quantile(boot, 0.025, na.rm = TRUE, names = FALSE),
      auc_ci_high = stats::quantile(boot, 0.975, na.rm = TRUE, names = FALSE),
      cross_tech_loco_auc = if (.BY$comparison_type == "cross_tech") mean(sinkhorn_ot_auc, na.rm = TRUE) else NA_real_,
      cross_tech_lift_vs_rCLR = if (.BY$comparison_type == "cross_tech") lift else NA_real_,
      within_tech_loco_auc = if (.BY$comparison_type == "within_tech") mean(sinkhorn_ot_auc, na.rm = TRUE) else NA_real_,
      within_tech_lift_vs_rCLR = if (.BY$comparison_type == "within_tech") lift else NA_real_,
      auc_var = stats::var(sinkhorn_ot_auc, na.rm = TRUE),
      p_auc_gt_0_5_hm = NA_real_,
      q_block_BH = NA_real_,
      provenance_block = paste(blocks, collapse = ","),
      status = if (length(blocks) < 2L) "descriptive_single_provenance_block" else
        "evaluable_cluster_bootstrap_by_provenance_block",
      status_reason = "",
      leakage_flag = any(leakage_flag),
      plausibility_flag = any(plausibility_flag),
      notes = sprintf("cluster_bootstrap_B=%d; n_blocks=%d; provisional reviewer artifact",
                      .sot_config$boot_B, length(blocks)))
  }, by = comparison_type]
  data.table::setcolorder(out, names(dt))
  out
}

.sot_make_figure <- function(dt) {
  ev <- dt[row_type == "loco" & status == "evaluable"]
  if (nrow(ev) == 0L) return(invisible(FALSE))
  ev[, comparison_label := factor(comparison_type,
                                  levels = c("within_tech", "cross_tech"),
                                  labels = c("Within technology", "Cross technology"))]
  p <- ggplot2::ggplot(ev, ggplot2::aes(x = comparison_label, y = sinkhorn_ot_auc,
                                        colour = technology)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.15) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0),
                        size = 2.0, alpha = 0.9) +
    ggplot2::labs(x = NULL, y = "LOCO AUC after Sinkhorn OT projection",
                  colour = "Held-out technology") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
  ggplot2::ggsave(.sot_paths$png, p, width = 6.8, height = 4.8, dpi = 320,
                  bg = "white")
  ggplot2::ggsave(.sot_paths$pdf, p, width = 6.8, height = 4.8, bg = "white")
  invisible(TRUE)
}

.sot_write_caption <- function(dt) {
  ev <- dt[row_type == "loco" & status == "evaluable"]
  cross <- dt[row_type == "summary" & comparison_type == "cross_tech"][1L]
  within <- dt[row_type == "summary" & comparison_type == "within_tech"][1L]
  lead <- if (nrow(cross) == 1L && is.finite(cross$cross_tech_loco_auc)) {
    sprintf("In this provisional first OT-based batch correction benchmark on miRNA, cross-technology NGS<->Toray LOCO AUC averaged %.3f (cluster-bootstrap 95%% CI %.3f to %.3f), with mean lift %.3f over train-frozen rCLR.",
            cross$cross_tech_loco_auc, cross$auc_ci_low, cross$auc_ci_high,
            cross$cross_tech_lift_vs_rCLR)
  } else {
    "This provisional first OT-based batch correction benchmark on miRNA did not yield an evaluable cross-technology summary."
  }
  within_txt <- if (nrow(within) == 1L && is.finite(within$within_tech_loco_auc)) {
    sprintf("Within-technology LOCO AUC averaged %.3f (cluster-bootstrap 95%% CI %.3f to %.3f), with mean lift %.3f over train-frozen rCLR.",
            within$within_tech_loco_auc, within$auc_ci_low, within$auc_ci_high,
            within$within_tech_lift_vs_rCLR)
  } else {
    "Within-technology summary was not evaluable."
  }
  leak_n <- sum(ev$leakage_flag == TRUE, na.rm = TRUE)
  plaus_n <- sum(ev$plausibility_flag == TRUE, na.rm = TRUE)
  flag_txt <- if (leak_n == 0L && plaus_n == 0L) {
    "Plausibility gate: no evaluable LOCO row exceeded the prespecified AUC >0.95/<0.05 threshold and no held-out cohort appeared in its training-barycenter cohort list."
  } else {
    paste0("Plausibility gate: ", plaus_n,
           " evaluable LOCO row(s) crossed the AUC >0.95/<0.05 threshold; leakage flags: ",
           leak_n, ". Flagged rows require targeted review before interpretation.")
  }
  txt <- c(
    "# Figure S55. Sinkhorn OT cohort-barycenter projection",
    "",
    lead,
    within_txt,
    "",
    "The scorer fits the CLR pivot pool, Sinkhorn epsilon/lambda, and entropic barycenter using training cohorts only. Held-out labels are used only after projection for AUC calculation. The plotted rows are leave-one-cohort-out evaluations; intervals in the summary table resample provenance blocks.",
    "",
    flag_txt,
    "",
    "This reviewer-package artifact is not manuscript-accepted evidence until code, result table, figure, and interpretation pass /triple-consensus, /plausibility-check, and /consistency-check."
  )
  writeLines(txt, .sot_paths$caption)
}

.sot_main <- function(output_dir = tempdir(), s5_path = NULL, s22_path = NULL,
                      kit_metadata_path = NULL, union_inventory_path = NULL,
                      sample_filters_path = NULL,
                      matched_set_audit_path = NULL,
                      paths = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  old_paths <- .sot_paths
  on.exit(.sot_paths <<- old_paths, add = TRUE)
  .sot_paths <<- utils::modifyList(.sot_default_paths(output_dir), paths %||% list())
  dir.create(dirname(.sot_paths$table), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(.sot_paths$png), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(.sot_paths$log), recursive = TRUE, showWarnings = FALSE)
  log_lines <- character()
  log_fun <- function(msg) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), msg)
    log_lines <<- c(log_lines, line)
    message(line)
  }
  set.seed(.sot_config$seed)
  log_fun("Sinkhorn OT S55 benchmark start")
  log_fun(sprintf("config atoms=%d max_cohort_atoms=%d tune=%s backend=%s",
                  .sot_config$n_atoms, .sot_config$max_cohort_atoms,
                  .sot_config$tune, .sot_config$backend))
  primary <- .sot_primary_cohorts(s5_path = s5_path, s22_path = s22_path)
  meta <- .sot_load_metadata(primary, kit_metadata_path = kit_metadata_path,
                             union_inventory_path = union_inventory_path)
  data <- .sot_load_data(primary, meta, log_fun,
                         union_inventory_path = union_inventory_path,
                         sample_filters_path = sample_filters_path,
                         matched_set_audit_path = matched_set_audit_path)
  failures <- attr(data, "failures")
  if (length(failures)) {
    for (nm in names(failures)) log_fun(sprintf("load_failed %s: %s", nm, failures[[nm]]))
  }
  rows <- list()
  i <- 0L
  for (held in names(data)) {
    for (cmp in c("within_tech", "cross_tech")) {
      i <- i + 1L
      log_fun(sprintf("LOCO %s %s", held, cmp))
      rows[[i]] <- .sot_eval_one(held, cmp, data, log_fun)
    }
  }
  loco <- data.table::rbindlist(rows, fill = TRUE)
  loco <- .sot_add_q(loco)
  summary <- .sot_summary_rows(loco)
  out <- data.table::rbindlist(list(loco, summary), fill = TRUE)
  data.table::fwrite(out, .sot_paths$table, sep = "\t", quote = FALSE, na = "NA")
  .sot_make_figure(out)
  .sot_write_caption(out)
  log_fun(sprintf("wrote table %s rows=%d", .sot_paths$table, nrow(out)))
  log_fun(sprintf("wrote figure %s and %s", .sot_paths$png, .sot_paths$pdf))
  log_fun(sprintf("wrote caption %s", .sot_paths$caption))
  log_fun("status counts:")
  log_lines <- c(log_lines, capture.output(print(out[, .N, by = .(row_type, comparison_type, status)])))
  log_lines <- c(log_lines, capture.output(sessionInfo()))
  writeLines(log_lines, .sot_paths$log)
  invisible(out)
}

if (sys.nframe() == 0L) {
  .sot_main()
}
