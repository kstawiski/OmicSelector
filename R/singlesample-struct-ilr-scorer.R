# =============================================================================
# ss-struct-ilr : structurally-informed ILR single-sample scorer
# -----------------------------------------------------------------------------
# Stage 1 of idea_report_1570797495812882433 ("Structurally Informed ILR + DA-cVAE").
# A frozen, single-sample, reference-free within-cohort discriminator that builds
# a Sequential Binary Partition (SBP) from miRNA genomic-cluster membership crossed
# with GC-content tiers, computes ILR balances per specimen (reusing ws_balance_ilr),
# standardises them with train-only statistics, and scores each specimen with a frozen
# diagonal-LDA linear rule. Quarantining high-GC features into isolated within-cluster
# balances isolates GC-driven dropout. Score time recomputes NOTHING from the test
# batch: the SBP, the per-coordinate standardisation, and the linear weights are all
# locked at fit. EXPLORATORY (paper 3): occupies the pre-reserved 44th within slot;
# never folded into the frozen 0/43 leakage-free primary.
# =============================================================================

# ---- internal: hyperparameter resolution -----------------------------------
.ss_struct_ilr_resolve_hp <- function(hp) {
  if (!is.list(hp)) stop("fit_ss_struct_ilr: hp must be a list")
  allowed <- c("pseudocount", "aggregate", "min_balance_coverage", "min_part")
  unknown <- setdiff(names(hp), allowed)
  if (length(unknown) > 0L) {
    stop("fit_ss_struct_ilr: unknown hp field(s): ", paste(unknown, collapse = ", "))
  }
  pseudocount <- hp$pseudocount   # NULL => ws_balance_ilr per-sample default
  if (!is.null(pseudocount) &&
      (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
       !is.finite(pseudocount) || pseudocount <= 0)) {
    stop("fit_ss_struct_ilr: hp$pseudocount must be NULL or a positive finite number")
  }
  aggregate <- hp$aggregate
  if (is.null(aggregate)) aggregate <- "gmean"
  aggregate <- match.arg(aggregate, c("gmean", "trimmed_gmean"))
  mbc <- hp$min_balance_coverage
  if (is.null(mbc)) mbc <- 0.5
  if (!is.numeric(mbc) || length(mbc) != 1L || !is.finite(mbc) || mbc < 0 || mbc > 1) {
    stop("fit_ss_struct_ilr: hp$min_balance_coverage must be in [0, 1]")
  }
  min_part <- hp$min_part
  if (is.null(min_part)) min_part <- 1L
  if (!is.numeric(min_part) || length(min_part) != 1L || !is.finite(min_part) ||
      min_part < 1L || min_part != as.integer(min_part)) {
    stop("fit_ss_struct_ilr: hp$min_part must be a positive integer")
  }
  list(pseudocount = pseudocount, aggregate = aggregate,
       min_balance_coverage = mbc, min_part = as.integer(min_part))
}

# ---- internal: validate annotation -----------------------------------------
.ss_struct_ilr_check_annotation <- function(annotation) {
  if (!is.data.frame(annotation)) {
    stop("fit_ss_struct_ilr: annotation must be a data.frame with columns ",
         "'feature', 'cluster', 'gc'")
  }
  req <- c("feature", "cluster", "gc")
  if (!all(req %in% names(annotation))) {
    stop("fit_ss_struct_ilr: annotation must have columns: ",
         paste(req, collapse = ", "))
  }
  annotation$feature <- as.character(annotation$feature)
  annotation$cluster <- as.character(annotation$cluster)
  annotation$gc      <- suppressWarnings(as.numeric(annotation$gc))
  annotation
}

# ---- internal: build the genomic-cluster / GC-tertile SBP -------------------
# Deterministic, frozen at fit from TRAIN features only. Two structural layers:
#   (1) cluster-coarse: a sequential binary partition over clusters (cluster k vs
#       the union of clusters k+1..K), in fixed cluster order;
#   (2) GC-within: high-GC tertile vs low-GC tertile WITHIN each cluster, so high-GC
#       dropout is quarantined to isolated within-cluster balances.
# GC tertile cut-offs are computed on the TRAIN annotated features (frozen). The SBP
# is returned as explicit {numerator, denominator} feature-id lists, so all
# structure (cluster membership, tertile assignment) is baked into the frozen object.
.ss_struct_ilr_build_sbp <- function(features, annotation, min_part = 1L) {
  ann <- annotation[match(features, annotation$feature), , drop = FALSE]
  ann$feature <- features
  keep <- !is.na(ann$cluster) & is.finite(ann$gc)
  ann <- ann[keep, , drop = FALSE]
  if (nrow(ann) < 2L) return(list())

  clusters <- sort(unique(ann$cluster))
  balances <- list()

  # (1) cluster-coarse sequential binary partition
  if (length(clusters) >= 2L) {
    for (k in seq_len(length(clusters) - 1L)) {
      num <- ann$feature[ann$cluster == clusters[k]]
      den <- ann$feature[ann$cluster %in% clusters[(k + 1L):length(clusters)]]
      if (length(num) >= min_part && length(den) >= min_part) {
        balances[[paste0("clust__", clusters[k])]] <-
          list(numerator = num, denominator = den)
      }
    }
  }

  # (2) GC high vs low tertile within each cluster (tertiles frozen from train)
  qs <- stats::quantile(ann$gc, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE)
  lo_cut <- qs[1]; hi_cut <- qs[2]
  for (cl in clusters) {
    sub <- ann[ann$cluster == cl, , drop = FALSE]
    if (nrow(sub) < 2L) next
    hi <- sub$feature[sub$gc >= hi_cut]
    lo <- sub$feature[sub$gc <= lo_cut]
    hi <- setdiff(hi, lo)  # exclusive when hi_cut == lo_cut (degenerate GC)
    if (length(hi) >= min_part && length(lo) >= min_part) {
      balances[[paste0("gc__", cl)]] <- list(numerator = hi, denominator = lo)
    }
  }
  balances
}

#' @title Fit the structurally-informed ILR single-sample scorer (ss-struct-ilr)
#'
#' @description
#' Builds a genomic-cluster / GC-tertile Sequential Binary Partition from the training
#' feature set + a miRNA \code{annotation} table, computes ILR balances per specimen
#' via \code{\link{ws_balance_ilr}}, standardises them with train-only mean/sd, and
#' learns a frozen diagonal-LDA linear rule (per-balance standardised case-minus-control
#' mean difference). The returned object serialises every artefact needed to score one
#' new specimen with no test-batch information.
#'
#' @param X_train Numeric matrix (samples x features) of non-negative abundances;
#'   columns are uniquely named miRNA ids.
#' @param y_train 0/1 labels aligned to \code{X_train} (1 = case).
#' @param meta_train Optional per-sample metadata; accepted for interface uniformity
#'   and ignored.
#' @param annotation data.frame with columns \code{feature} (miRNA id),
#'   \code{cluster} (genomic-cluster id), \code{gc} (GC fraction).
#' @param hp List of frozen hyperparameters: \code{pseudocount} (default \code{NULL} =
#'   the scale-invariant per-sample default of \code{\link{ws_balance_ilr}}),
#'   \code{aggregate} ("gmean"/"trimmed_gmean", default "gmean"),
#'   \code{min_balance_coverage} (default 0.5), \code{min_part} (default 1L).
#'
#' @return A \code{ss_struct_ilr_model} list: frozen \code{sbp}, \code{balance_names},
#'   \code{center}, \code{scale}, \code{weights}, \code{feature_universe}, the resolved
#'   \code{hp}, and an \code{annotation_features} record.
#' @seealso \code{\link{score_ss_struct_ilr}}, \code{\link{ws_balance_ilr}}
#' @export
fit_ss_struct_ilr <- function(X_train, y_train, meta_train = NULL,
                              annotation, hp = list()) {
  X_train <- .reo_check_matrix(X_train, "fit_ss_struct_ilr", "X_train")
  .reo_check_meta(meta_train, nrow(X_train), "fit_ss_struct_ilr", "meta_train")
  y <- .reo_check_labels(y_train, nrow(X_train), "fit_ss_struct_ilr")
  if (missing(annotation) || is.null(annotation)) {
    stop("fit_ss_struct_ilr: annotation (feature/cluster/gc) is required")
  }
  annotation <- .ss_struct_ilr_check_annotation(annotation)
  hp <- .ss_struct_ilr_resolve_hp(hp)

  features <- colnames(X_train)
  sbp <- .ss_struct_ilr_build_sbp(features, annotation, hp$min_part)

  model <- list(
    sbp = sbp,
    pseudocount = hp$pseudocount,
    aggregate = hp$aggregate,
    min_balance_coverage = hp$min_balance_coverage,
    min_part = hp$min_part,
    center = numeric(0),
    scale = numeric(0),
    weights = numeric(0),
    balance_names = names(sbp),
    feature_universe = features,
    annotation_features = intersect(annotation$feature, features),
    hp = hp
  )

  if (length(sbp) == 0L) {
    # Degenerate: no structural balances constructible -> neutral scorer (all 0).
    class(model) <- "ss_struct_ilr_model"
    return(model)
  }

  B <- ws_balance_ilr(X_train, balances = sbp, pseudocount = hp$pseudocount,
                      aggregate = hp$aggregate,
                      min_balance_coverage = hp$min_balance_coverage)
  B <- as.matrix(B)
  storage.mode(B) <- "double"

  center <- apply(B, 2L, function(v) {
    m <- mean(v, na.rm = TRUE); if (is.finite(m)) m else 0
  })
  scale_ <- apply(B, 2L, function(v) {
    s <- stats::sd(v, na.rm = TRUE); if (is.finite(s) && s > 0) s else 1
  })
  Bs <- sweep(sweep(B, 2L, center, "-"), 2L, scale_, "/")

  weights <- vapply(seq_len(ncol(Bs)), function(j) {
    v <- Bs[, j]
    m1 <- mean(v[y == 1L], na.rm = TRUE)
    m0 <- mean(v[y == 0L], na.rm = TRUE)
    d <- m1 - m0
    if (is.finite(d)) d else 0
  }, numeric(1L))

  model$center <- center
  model$scale <- scale_
  model$weights <- weights
  model$balance_names <- colnames(B)
  class(model) <- "ss_struct_ilr_model"
  model
}

#' @title Score the structurally-informed ILR single-sample scorer (ss-struct-ilr)
#'
#' @description
#' Scores each row of \code{X} independently with the frozen
#' \code{\link{fit_ss_struct_ilr}} model. Each specimen is restricted to the frozen
#' training feature universe present in \code{X}, ILR balances are computed from the
#' specimen's own values via \code{\link{ws_balance_ilr}} (the per-sample pseudocount
#' makes balances scale-invariant), standardised with the frozen train mean/sd, and
#' combined with the frozen diagonal-LDA weights. No scored-batch centering, quantiles,
#' or renormalisation is used, so the scorer is row-equivariant
#' (\code{\link{singlesample_assert_row_equivariant}}). Specimens with no constructible
#' balance receive the neutral score 0.
#'
#' @param model A \code{ss_struct_ilr_model} from \code{\link{fit_ss_struct_ilr}}.
#' @param X Numeric matrix (samples x features) of non-negative abundances; named cols.
#' @param meta Optional per-sample metadata; accepted for interface uniformity, ignored.
#' @return Finite numeric vector of length \code{nrow(X)}; higher = more case-like.
#' @seealso \code{\link{fit_ss_struct_ilr}}
#' @export
score_ss_struct_ilr <- function(model, X, meta = NULL) {
  if (!inherits(model, "ss_struct_ilr_model")) {
    stop("score_ss_struct_ilr: model must have class ss_struct_ilr_model")
  }
  X <- .reo_check_matrix(X, "score_ss_struct_ilr", "X")
  .reo_check_meta(meta, nrow(X), "score_ss_struct_ilr", "meta")

  if (length(model$sbp) == 0L || length(model$weights) == 0L) {
    return(rep(0, nrow(X)))
  }
  present <- intersect(model$feature_universe, colnames(X))
  if (length(present) == 0L) return(rep(0, nrow(X)))
  Xr <- X[, present, drop = FALSE]

  B <- ws_balance_ilr(Xr, balances = model$sbp, pseudocount = model$pseudocount,
                      aggregate = model$aggregate,
                      min_balance_coverage = model$min_balance_coverage)
  B <- as.matrix(B)
  storage.mode(B) <- "double"
  # Align columns to the frozen balance order (defensive; same sbp -> same order).
  if (!is.null(colnames(B)) && length(model$balance_names) == ncol(B)) {
    idx <- match(model$balance_names, colnames(B))
    if (!any(is.na(idx))) B <- B[, idx, drop = FALSE]
  }

  Bs <- sweep(sweep(B, 2L, model$center, "-"), 2L, model$scale, "/")
  w <- model$weights
  out <- vapply(seq_len(nrow(Bs)), function(i) {
    sum(w * Bs[i, ], na.rm = TRUE)
  }, numeric(1L))
  out <- as.numeric(out)
  if (length(out) != nrow(X) || any(!is.finite(out))) {
    stop("score_ss_struct_ilr: scorer produced non-finite or wrong-length output")
  }
  out
}
