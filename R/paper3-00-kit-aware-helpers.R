# Internal helpers for kit-aware compositional within-sample scoring.
# Internal helpers used by the package kit-aware compositional scorers and tests.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}

.ka_unknown_tokens <- c("", "unknown", "unk", "na", "n/a", "nr", "not reported",
                        "paper_not_found", "not_applicable", "not_ngs")

.ka_clean_label <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- NA_character_
  x
}

.ka_is_known_label <- function(x) {
  x <- .ka_clean_label(x)
  !is.na(x) & nzchar(x) & !(tolower(x) %in% .ka_unknown_tokens)
}

.ka_validate_inputs <- function(expr_matrix, sample_meta, kit_label_col,
                                panel_features = NULL) {
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("expr_matrix must be a matrix/data.frame with samples in rows and features in columns.")
  }
  x <- as.matrix(expr_matrix)
  storage.mode(x) <- "double"
  if (is.null(colnames(x))) stop("expr_matrix must have feature column names.")
  if (nrow(x) == 0L || ncol(x) == 0L) stop("expr_matrix must be non-empty.")
  if (any(x < 0, na.rm = TRUE)) {
    stop("kit-aware compositional scorers require non-negative counts/abundances.")
  }
  if (!is.data.frame(sample_meta)) {
    sample_meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)
  }
  if (nrow(sample_meta) != nrow(x)) {
    stop("sample_meta must have one row per expr_matrix sample.")
  }
  if (!(kit_label_col %in% names(sample_meta))) {
    stop("kit_label_col not found in sample_meta: ", kit_label_col)
  }
  if (!is.null(panel_features)) {
    panel_features <- intersect(as.character(panel_features), colnames(x))
  }
  list(X = x, sample_meta = sample_meta, panel_features = panel_features)
}

.ka_resolve_strata <- function(sample_meta, kit_label_col,
                               biofluid_col = "biofluid",
                               global_label = "global") {
  kit <- .ka_clean_label(sample_meta[[kit_label_col]])
  known_kit <- .ka_is_known_label(kit)
  out <- rep(global_label, nrow(sample_meta))
  source <- rep("global", nrow(sample_meta))
  out[known_kit] <- paste0("kit:", kit[known_kit])
  source[known_kit] <- "kit"

  if (!is.null(biofluid_col) && biofluid_col %in% names(sample_meta)) {
    bio <- .ka_clean_label(sample_meta[[biofluid_col]])
    known_bio <- .ka_is_known_label(bio)
    use_bio <- !known_kit & known_bio
    out[use_bio] <- paste0("biofluid:", bio[use_bio])
    source[use_bio] <- "biofluid"
  }

  list(label = out, source = source)
}

.ka_row_pseudocount <- function(x, pseudocount = NULL) {
  if (!is.null(pseudocount)) return(rep(as.numeric(pseudocount), nrow(x)))
  pc <- 1e-6 * rowSums(x, na.rm = TRUE)
  pc[!is.finite(pc) | pc <= 0] <- .Machine$double.eps
  pc
}

.ka_log_proportions <- function(x, pseudocount = NULL, features = NULL) {
  pc <- .ka_row_pseudocount(x, pseudocount)
  x_pos_full <- sweep(x, 1L, pc, `+`)
  denom <- rowSums(x_pos_full, na.rm = TRUE)
  denom[!is.finite(denom) | denom <= 0] <- NA_real_
  if (!is.null(features)) {
    x_pos_full <- x_pos_full[, intersect(features, colnames(x_pos_full)), drop = FALSE]
  }
  sweep(log(x_pos_full), 1L, log(denom), `-`)
}

.ka_log_cpm <- function(x, pseudocount = 0.5, features = NULL) {
  x_pos_full <- x + pseudocount
  denom <- rowSums(x_pos_full, na.rm = TRUE)
  denom[!is.finite(denom) | denom <= 0] <- NA_real_
  if (!is.null(features)) {
    x_pos_full <- x_pos_full[, intersect(features, colnames(x_pos_full)), drop = FALSE]
  }
  log(sweep(x_pos_full, 1L, denom, `/`) * 1e6)
}

.ka_trimmed_centering_features <- function(x,
                                           trim_upper = 0.10,
                                           trim_lower = 0.05,
                                           exclude_features = c("hsa-miR-451a",
                                                                "hsa-miR-16-5p",
                                                                "hsa-miR-486-5p",
                                                                "hsa-miR-144-3p",
                                                                "hsa-miR-223-3p"),
                                           min_centering_size = 8L) {
  if (ncol(x) == 0L) return(character(0))
  abundance <- colMedians_or_median(x)
  names(abundance) <- colnames(x)
  candidates <- setdiff(names(abundance)[is.finite(abundance)], exclude_features)
  if (length(candidates) < min_centering_size) {
    candidates <- names(abundance)[is.finite(abundance)]
  }
  if (length(candidates) < min_centering_size) return(candidates)
  ord <- candidates[order(abundance[candidates], decreasing = TRUE)]
  n <- length(ord)
  lo <- floor(trim_upper * n) + 1L
  hi <- n - floor(trim_lower * n)
  keep <- ord[seq.int(lo, hi)]
  if (length(keep) < min_centering_size) candidates else keep
}

colMedians_or_median <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colMedians(as.matrix(x), na.rm = TRUE))
  }
  apply(x, 2L, stats::median, na.rm = TRUE)
}

.ka_panel_and_weights <- function(panel_features, feature_weights, available_features) {
  panel <- intersect(as.character(panel_features), available_features)
  if (length(panel) == 0L) stop("No panel_features are present in expr_matrix.")
  if (is.null(feature_weights)) {
    weights <- stats::setNames(rep(1, length(panel)), panel)
  } else {
    weights <- as.numeric(feature_weights[panel])
    weights[!is.finite(weights)] <- 1
    weights <- stats::setNames(weights, panel)
  }
  list(panel = panel, weights = weights)
}

.ka_weighted_panel_sum <- function(mat, panel, weights) {
  panel <- intersect(panel, colnames(mat))
  if (length(panel) == 0L) return(rep(NA_real_, nrow(mat)))
  w <- as.numeric(weights[panel])
  w[!is.finite(w)] <- 1
  as.numeric(mat[, panel, drop = FALSE] %*% w)
}

.ka_fast_auc <- function(y, score) {
  y <- if (is.factor(y)) as.integer(y) - 1L else as.integer(y)
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

.ka_hm_var <- function(auc, n_pos, n_neg) {
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
  if (!is.finite(v) || v < 0) NA_real_ else max(v, 1e-8)
}

.ka_select_panel_from_matrix <- function(z_train, y_train, panel_size = 20L,
                                         min_detection_rate = 0.30) {
  y_train <- as.integer(y_train)
  det_rate <- colMeans(is.finite(z_train) & z_train != 0, na.rm = TRUE)
  candidates <- names(det_rate)[det_rate >= min_detection_rate]
  if (length(candidates) < panel_size) {
    candidates <- names(sort(det_rate, decreasing = TRUE))[seq_len(min(ncol(z_train), panel_size * 4L))]
  }
  aucs <- vapply(candidates, function(feat) .ka_fast_auc(y_train, z_train[, feat]),
                 numeric(1L))
  strength <- pmax(aucs, 1 - aucs, na.rm = TRUE)
  strength[!is.finite(strength)] <- 0.5
  ord <- order(strength, decreasing = TRUE)
  panel <- candidates[ord[seq_len(min(panel_size, length(ord)))]]
  weights <- vapply(panel, function(feat) {
    med1 <- stats::median(z_train[y_train == 1L, feat], na.rm = TRUE)
    med0 <- stats::median(z_train[y_train == 0L, feat], na.rm = TRUE)
    if (is.finite(med1) && is.finite(med0) && med1 < med0) -1 else 1
  }, numeric(1L))
  stats::setNames(weights, panel)
}

.ka_global_rclr_matrix <- function(expr_matrix, centering_features,
                                   pseudocount = NULL) {
  need <- unique(c(centering_features, colnames(expr_matrix)))
  logp <- .ka_log_proportions(expr_matrix, pseudocount = pseudocount,
                              features = need)
  center <- intersect(centering_features, colnames(logp))
  if (length(center) == 0L) center <- colnames(logp)
  ctr <- rowMeans(logp[, center, drop = FALSE], na.rm = TRUE)
  sweep(logp, 1L, ctr, `-`)
}
