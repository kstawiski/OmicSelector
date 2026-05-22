# Internal helpers for Group-DRO kit x biofluid robust scoring.
# Sourced by code/methods/group_dro_scorer.R and tests.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}

.gdro_unknown_tokens <- c("", "unknown", "unk", "na", "n/a", "nr", "not reported",
                          "paper_not_found", "not_applicable", "not_ngs")

.gdro_as_matrix <- function(x, name = "expr_matrix") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop(name, " must be a numeric matrix/data.frame.")
  storage.mode(x) <- "double"
  if (is.null(colnames(x))) stop(name, " must have feature column names.")
  if (nrow(x) == 0L || ncol(x) == 0L) stop(name, " must be non-empty.")
  if (any(x < 0, na.rm = TRUE)) {
    stop(name, " must contain non-negative counts/abundances before CLR.")
  }
  x[!is.finite(x)] <- NA_real_
  x
}

.gdro_as_binary <- function(y) {
  if (is.factor(y)) y <- as.integer(y) - 1L
  y <- as.integer(y)
  if (!all(y %in% c(0L, 1L), na.rm = TRUE)) {
    stop("y must be binary 0/1 or a two-level factor.")
  }
  if (length(unique(y[!is.na(y)])) != 2L) {
    stop("y must contain both outcome classes.")
  }
  y
}

.gdro_clean_label <- function(x, unknown = "unknown") {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in% .gdro_unknown_tokens] <- unknown
  x
}

.gdro_validate_meta <- function(meta, n, kit_col, biofluid_col) {
  if (!is.data.frame(meta)) meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (nrow(meta) != n) stop("sample_meta must have one row per expression sample.")
  missing <- setdiff(c(kit_col, biofluid_col), names(meta))
  if (length(missing) > 0L) stop("sample_meta missing column(s): ", paste(missing, collapse = ", "))
  meta
}

.gdro_make_groups <- function(sample_meta, kit_col = "kit_family",
                              biofluid_col = "biofluid") {
  kit <- .gdro_clean_label(sample_meta[[kit_col]])
  bio <- .gdro_clean_label(sample_meta[[biofluid_col]])
  paste0("kit=", kit, "|biofluid=", bio)
}

.gdro_row_pseudocount <- function(x, pseudocount = NULL) {
  if (!is.null(pseudocount)) return(rep(as.numeric(pseudocount), nrow(x)))
  pc <- 1e-6 * rowSums(x, na.rm = TRUE)
  pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  pc
}

.gdro_log_proportions <- function(x, pseudocount = NULL, features = NULL) {
  pc <- .gdro_row_pseudocount(x, pseudocount)
  x_pos <- sweep(x, 1L, pc, `+`)
  denom <- rowSums(x_pos, na.rm = TRUE)
  denom[!is.finite(denom) | denom <= 0] <- NA_real_
  if (!is.null(features)) {
    features <- intersect(features, colnames(x_pos))
    x_pos <- x_pos[, features, drop = FALSE]
  }
  out <- sweep(log(x_pos), 1L, log(denom), `-`)
  out[!is.finite(out)] <- NA_real_
  out
}

.gdro_col_medians <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colMedians(as.matrix(x), na.rm = TRUE)
  } else {
    apply(x, 2L, stats::median, na.rm = TRUE)
  }
}

.gdro_select_pivot_pool <- function(x,
                                    trim_upper = 0.10,
                                    trim_lower = 0.05,
                                    min_pivots = 8L,
                                    exclude_features = c("hsa-miR-451a",
                                                         "hsa-miR-16-5p",
                                                         "hsa-miR-486-5p",
                                                         "hsa-miR-144-3p",
                                                         "hsa-miR-223-3p")) {
  abundance <- .gdro_col_medians(x)
  names(abundance) <- colnames(x)
  candidates <- setdiff(names(abundance)[is.finite(abundance)], exclude_features)
  if (length(candidates) < min_pivots) candidates <- names(abundance)[is.finite(abundance)]
  if (length(candidates) < min_pivots) return(candidates)
  ord <- candidates[order(abundance[candidates], decreasing = TRUE)]
  n <- length(ord)
  lo <- floor(trim_upper * n) + 1L
  hi <- n - floor(trim_lower * n)
  keep <- ord[seq.int(lo, hi)]
  if (length(keep) < min_pivots) candidates else keep
}

.gdro_clr_fit <- function(x, pivot_features = NULL, pseudocount = NULL,
                          trim_upper = 0.10, trim_lower = 0.05,
                          min_pivots = 8L) {
  x <- .gdro_as_matrix(x)
  if (is.null(pivot_features)) {
    pivot_features <- .gdro_select_pivot_pool(x, trim_upper = trim_upper,
                                              trim_lower = trim_lower,
                                              min_pivots = min_pivots)
  }
  pivot_features <- intersect(pivot_features, colnames(x))
  if (length(pivot_features) == 0L) stop("No pivot features available for train-only CLR.")
  list(features = colnames(x),
       pivot_features = pivot_features,
       pseudocount = pseudocount,
       trim_upper = trim_upper,
       trim_lower = trim_lower,
       min_pivots = min_pivots)
}

.gdro_clr_apply <- function(x, fit, features = fit$features) {
  x <- .gdro_as_matrix(x)
  missing <- setdiff(fit$features, colnames(x))
  if (length(missing) > 0L) {
    add <- matrix(0, nrow = nrow(x), ncol = length(missing),
                  dimnames = list(rownames(x), missing))
    x <- cbind(x, add)
  }
  x <- x[, fit$features, drop = FALSE]
  need <- unique(c(fit$pivot_features, features))
  lp <- .gdro_log_proportions(x, pseudocount = fit$pseudocount, features = need)
  piv <- intersect(fit$pivot_features, colnames(lp))
  if (length(piv) == 0L) piv <- colnames(lp)
  center <- rowMeans(lp[, piv, drop = FALSE], na.rm = TRUE)
  z <- sweep(lp[, intersect(features, colnames(lp)), drop = FALSE], 1L, center, `-`)
  z[!is.finite(z)] <- NA_real_
  z
}

.gdro_scale_fit <- function(z) {
  z <- as.matrix(z)
  center <- apply(z, 2L, stats::median, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  madv <- apply(abs(sweep(z, 2L, center, `-`)), 2L, stats::median, na.rm = TRUE)
  scale <- 1.4826 * madv
  fallback <- stats::median(scale[is.finite(scale) & scale > 1e-8], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 1e-8) fallback <- 1
  scale[!is.finite(scale) | scale <= 1e-8] <- fallback
  list(center = center, scale = scale, features = colnames(z))
}

.gdro_scale_apply <- function(z, scaler) {
  z <- as.matrix(z)
  missing <- setdiff(scaler$features, colnames(z))
  if (length(missing) > 0L) {
    add <- matrix(NA_real_, nrow = nrow(z), ncol = length(missing),
                  dimnames = list(rownames(z), missing))
    z <- cbind(z, add)
  }
  z <- z[, scaler$features, drop = FALSE]
  for (j in seq_len(ncol(z))) {
    bad <- !is.finite(z[, j])
    if (any(bad)) z[bad, j] <- scaler$center[j]
  }
  out <- sweep(sweep(z, 2L, scaler$center, `-`), 2L, scaler$scale, `/`)
  out[!is.finite(out)] <- 0
  pmax(pmin(out, 8), -8)
}

.gdro_sigmoid <- function(x) {
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  ex <- exp(x[!pos])
  out[!pos] <- ex / (1 + ex)
  out
}

.gdro_logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

.gdro_group_indices <- function(groups) {
  groups <- factor(groups)
  split(seq_along(groups), groups)
}

.gdro_logistic_group_dro <- function(z, y, groups,
                                     eta = 0.05,
                                     l2 = 0.01,
                                     learning_rate = 0.05,
                                     epochs = 200L,
                                     max_grad_norm = 10,
                                     seed = 42L) {
  set.seed(seed)
  z <- as.matrix(z)
  y <- .gdro_as_binary(y)
  groups <- factor(groups)
  if (length(y) != nrow(z) || length(groups) != nrow(z)) {
    stop("z, y, and groups must have aligned sample rows.")
  }
  gidx <- .gdro_group_indices(groups)
  glev <- names(gidx)
  p <- ncol(z)
  beta <- rep(0, p)
  ybar <- mean(y)
  intercept <- stats::qlogis(min(max(ybar, 1e-4), 1 - 1e-4))
  log_q <- rep(-log(length(gidx)), length(gidx))
  names(log_q) <- glev
  last_group_loss <- stats::setNames(rep(NA_real_, length(gidx)), glev)

  for (ep in seq_len(max(1L, as.integer(epochs)))) {
    pred <- .gdro_sigmoid(as.numeric(intercept + z %*% beta))
    pred <- pmin(pmax(pred, 1e-8), 1 - 1e-8)
    loss_i <- -(y * log(pred) + (1 - y) * log(1 - pred))
    group_loss <- vapply(gidx, function(idx) mean(loss_i[idx], na.rm = TRUE),
                         numeric(1L))
    group_loss[!is.finite(group_loss)] <- max(group_loss[is.finite(group_loss)] %||% 0)
    log_q <- log_q + eta * group_loss
    log_q <- log_q - .gdro_logsumexp(log_q)
    q <- exp(log_q)

    w <- numeric(length(y))
    for (g in glev) w[gidx[[g]]] <- q[g] / length(gidx[[g]])
    resid <- w * (pred - y)
    grad_beta <- as.numeric(crossprod(z, resid)) + l2 * beta
    grad_intercept <- sum(resid)
    gnorm <- sqrt(sum(grad_beta^2) + grad_intercept^2)
    if (is.finite(gnorm) && gnorm > max_grad_norm) {
      scale <- max_grad_norm / gnorm
      grad_beta <- grad_beta * scale
      grad_intercept <- grad_intercept * scale
    }
    step <- learning_rate / sqrt(ep)
    beta <- beta - step * grad_beta
    intercept <- intercept - step * grad_intercept
    last_group_loss <- group_loss
  }

  names(beta) <- colnames(z)
  q <- exp(log_q)
  list(beta = beta,
       intercept = intercept,
       group_weights = q,
       group_losses = last_group_loss,
       eta = eta,
       l2 = l2,
       learning_rate = learning_rate,
       epochs = as.integer(epochs),
       n_groups = length(gidx))
}

.gdro_safe_auc <- function(y, score) {
  y <- .gdro_as_binary(y)
  score <- as.numeric(score)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (length(y) < 4L || length(unique(y)) != 2L || length(unique(score)) < 2L) {
    return(NA_real_)
  }
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  as.numeric((sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
}

.gdro_hanley_mcneil <- function(auc, n_pos, n_neg) {
  if (!is.finite(auc) || n_pos < 1L || n_neg < 1L) {
    return(c(se = NA_real_, low = NA_real_, high = NA_real_))
  }
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  var <- (auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) +
            (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg)
  if (!is.finite(var) || var < 0) return(c(se = NA_real_, low = NA_real_, high = NA_real_))
  se <- sqrt(max(var, 1e-10))
  c(se = se, low = max(0, auc - 1.96 * se), high = min(1, auc + 1.96 * se))
}

.gdro_stratified_folds <- function(y, k = 3L, seed = 42L) {
  y <- .gdro_as_binary(y)
  k <- min(as.integer(k), max(2L, min(table(y))))
  set.seed(seed)
  fold <- integer(length(y))
  for (cl in sort(unique(y))) {
    idx <- which(y == cl)
    fold[idx] <- sample(rep_len(seq_len(k), length(idx)))[seq_along(idx)]
  }
  fold
}

.gdro_select_panel <- function(z_train, y_train, panel_size = 20L,
                               min_detection_rate = 0.05) {
  z_train <- as.matrix(z_train)
  y_train <- .gdro_as_binary(y_train)
  det <- colMeans(is.finite(z_train) & abs(z_train) > 1e-12, na.rm = TRUE)
  candidates <- names(det)[det >= min_detection_rate]
  if (length(candidates) < panel_size) {
    candidates <- names(sort(det, decreasing = TRUE))[seq_len(min(ncol(z_train), panel_size * 4L))]
  }
  aucs <- vapply(candidates, function(feat) .gdro_safe_auc(y_train, z_train[, feat]),
                 numeric(1L))
  strength <- pmax(aucs, 1 - aucs, na.rm = TRUE)
  strength[!is.finite(strength)] <- 0.5
  ord <- order(strength, decreasing = TRUE)
  panel <- candidates[ord[seq_len(min(panel_size, length(ord)))]]
  directions <- vapply(panel, function(feat) {
    med1 <- stats::median(z_train[y_train == 1L, feat], na.rm = TRUE)
    med0 <- stats::median(z_train[y_train == 0L, feat], na.rm = TRUE)
    if (is.finite(med1) && is.finite(med0) && med1 < med0) -1 else 1
  }, numeric(1L))
  stats::setNames(directions, panel)
}

.gdro_tune_grid <- function(x, y, meta, kit_col, biofluid_col,
                            eta_grid, l2_grid, panel_size,
                            inner_folds = 3L,
                            learning_rate = 0.05,
                            epochs = 120L,
                            seed = 42L) {
  y <- .gdro_as_binary(y)
  fold <- .gdro_stratified_folds(y, k = inner_folds, seed = seed)
  grid <- expand.grid(eta = eta_grid, l2 = l2_grid)
  grid$mean_worst_group_auc <- NA_real_
  grid$mean_auc <- NA_real_

  for (gi in seq_len(nrow(grid))) {
    fold_worst <- numeric(0)
    fold_auc <- numeric(0)
    for (ff in sort(unique(fold))) {
      train <- fold != ff
      val <- fold == ff
      if (length(unique(y[val])) != 2L || length(unique(y[train])) != 2L) next
      clr <- .gdro_clr_fit(x[train, , drop = FALSE])
      z_tr_all <- .gdro_clr_apply(x[train, , drop = FALSE], clr)
      panel <- names(.gdro_select_panel(z_tr_all, y[train], panel_size = panel_size))
      z_tr <- z_tr_all[, panel, drop = FALSE]
      scaler <- .gdro_scale_fit(z_tr)
      z_tr <- .gdro_scale_apply(z_tr, scaler)
      z_val <- .gdro_scale_apply(.gdro_clr_apply(x[val, , drop = FALSE], clr, panel), scaler)
      g_tr <- .gdro_make_groups(meta[train, , drop = FALSE], kit_col, biofluid_col)
      fit <- .gdro_logistic_group_dro(
        z_tr, y[train], g_tr, eta = grid$eta[gi], l2 = grid$l2[gi],
        learning_rate = learning_rate, epochs = epochs, seed = seed + gi + ff)
      prob <- .gdro_sigmoid(as.numeric(fit$intercept + z_val %*% fit$beta))
      g_val <- .gdro_make_groups(meta[val, , drop = FALSE], kit_col, biofluid_col)
      auc <- .gdro_safe_auc(y[val], prob)
      g_auc <- vapply(split(seq_along(prob), g_val), function(idx) {
        if (length(unique(y[val][idx])) != 2L) return(NA_real_)
        .gdro_safe_auc(y[val][idx], prob[idx])
      }, numeric(1L))
      fold_auc <- c(fold_auc, auc)
      finite_g <- g_auc[is.finite(g_auc)]
      if (length(finite_g) > 0L) fold_worst <- c(fold_worst, min(finite_g))
    }
    grid$mean_auc[gi] <- mean(fold_auc, na.rm = TRUE)
    grid$mean_worst_group_auc[gi] <- mean(fold_worst, na.rm = TRUE)
  }
  grid$score <- ifelse(is.finite(grid$mean_worst_group_auc),
                       grid$mean_worst_group_auc, grid$mean_auc)
  grid$score[!is.finite(grid$score)] <- -Inf
  grid <- grid[order(grid$score, grid$mean_auc, decreasing = TRUE), ]
  list(best = grid[1L, , drop = FALSE], grid = grid)
}

.gdro_re_pool_auc <- function(dt, method_col = "method", auc_col = "observed_auc",
                              se_col = "auc_se", block_col = "provenance_block") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  d <- data.table::as.data.table(dt)
  d <- d[is.finite(get(auc_col)) & is.finite(get(se_col)) & get(se_col) > 0]
  if (nrow(d) == 0L) return(data.table::data.table())
  d[, {
    block <- .SD[, {
      vi <- pmax(get(se_col)^2, 1e-10)
      wi <- 1 / vi
      .(yi = sum(wi * get(auc_col)) / sum(wi), vi = 1 / sum(wi))
    }, by = block_col]
    yi <- block$yi
    vi <- block$vi
    wi <- 1 / vi
    fixed <- sum(wi * yi) / sum(wi)
    q <- sum(wi * (yi - fixed)^2)
    c_val <- sum(wi) - sum(wi^2) / sum(wi)
    tau2 <- if (length(yi) > 1L && c_val > 0) max(0, (q - (length(yi) - 1)) / c_val) else NA_real_
    wre <- 1 / (vi + (tau2 %||% 0))
    pooled <- sum(wre * yi) / sum(wre)
    se <- sqrt(1 / sum(wre))
    .(n_loco = .N,
      n_provenance_blocks = data.table::uniqueN(get(block_col)),
      average_auc = pooled,
      average_auc_se = se,
      average_auc_ci_low = max(0, pooled - 1.96 * se),
      average_auc_ci_high = min(1, pooled + 1.96 * se),
      tau2 = tau2,
      i2 = if (is.finite(q) && q > 0) max(0, (q - (length(yi) - 1)) / q) else NA_real_)
  }, by = method_col]
}

.gdro_cluster_boot_diff <- function(dt, B = 500L, seed = 42L,
                                    method_name = "group_dro_kit_robust") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  d <- data.table::as.data.table(dt)
  d <- d[benchmark_component == "LOCO" & method %in% c("rCLR_baseline", method_name) &
           is.finite(observed_auc)]
  wide <- data.table::dcast(d, held_out + provenance_block ~ method, value.var = "observed_auc")
  if (!all(c("rCLR_baseline", method_name) %in% names(wide))) return(c(NA_real_, NA_real_, NA_real_))
  wide <- wide[is.finite(rCLR_baseline) & is.finite(get(method_name))]
  if (nrow(wide) == 0L) return(c(NA_real_, NA_real_, NA_real_))
  set.seed(seed)
  blocks <- unique(wide$provenance_block)
  vals <- replicate(B, {
    b <- sample(blocks, length(blocks), replace = TRUE)
    samp <- data.table::rbindlist(lapply(b, function(bb) wide[provenance_block == bb]), fill = TRUE)
    mean(samp[[method_name]] - samp$rCLR_baseline, na.rm = TRUE)
  })
  c(mean = mean(vals, na.rm = TRUE),
    low = stats::quantile(vals, 0.025, na.rm = TRUE, names = FALSE),
    high = stats::quantile(vals, 0.975, na.rm = TRUE, names = FALSE))
}
