#' @title Within-Sample Normalization for Biomarker Panels
#'
#' @description
#' Functions for within-sample normalization that operate ONLY on a single
#' sample's feature values. Because these methods use no population-level
#' statistics — each sample is normalized independently using only its own
#' feature values — they cancel specific classes of per-sample nuisance
#' variation WITHOUT a reference cohort. Invariance scope is method-specific
#' (see Details). These methods are NOT a general "batch-effect remover":
#' gene-specific batch effects, non-linear platform response, and
#' measurement-specific bias are NOT cancelled.
#'
#' This is the key insight for clinical deployment: a diagnostic test based on
#' within-sample normalization requires NO reference cohort. Draw blood, measure
#' k biomarkers, compute the relative pattern, classify.
#'
#' @details
#' Within-sample normalization addresses one component of the batch-effect
#' problem in circulating biomarker diagnostics: per-sample nuisance variation
#' (collection-site offset, loading amount, total-library scaling). Traditional
#' normalization methods (quantile, ComBat, z-score) require a reference
#' population, which is impractical for point-of-care diagnostics. Within-sample
#' methods sidestep the reference-cohort requirement but each has a specific,
#' narrower invariance scope:
#' \itemize{
#'   \item \code{ws_logratio} / \code{ws_ratio_image}: additive-shift-invariant
#'     in log-space (equivalently, invariant to a multiplicative per-sample
#'     factor applied UNIFORMLY across all features in raw space). NOT
#'     invariant to gene-specific shifts or non-linear platform response.
#'   \item \code{ws_zscore}: invariant to affine per-sample transforms
#'     (scale + shift) applied uniformly across features.
#'   \item \code{ws_rank}: invariant to any strictly monotonic per-sample
#'     transform of all features.
#'   \item \code{ws_minmax}: invariant to uniform per-sample affine transforms
#'     but sensitive to outlier features.
#' }
#'
#' Available methods:
#' \itemize{
#'   \item \code{ws_minmax}: Scale each sample's features to [0,1]
#'   \item \code{ws_rank}: Replace values with within-sample ranks
#'   \item \code{ws_zscore}: Standardize across features within each sample
#'   \item \code{ws_logratio}: Compute all pairwise log-ratios (self-normalizing)
#'   \item \code{ws_ratio_image}: Create pairwise ratio matrix (for CNN input)
#' }
#'
#' @references
#' Stawiski K et al. (2026). Within-sample normalization enables batch-effect-free
#' pan-cancer detection from circulating miRNA signatures. (in preparation)
#'
#' @name within-sample
NULL


#' @title Within-Sample Min-Max Normalization
#'
#' @description
#' Scales each sample's feature values to [0,1] range using only that sample's
#' minimum and maximum. Completely batch-immune.
#'
#' @param x A numeric matrix (samples x features) or vector (single sample)
#'
#' @return Matrix of same dimensions with values in [0,1]
#'
#' @examples
#' \dontrun{
#' # Single sample
#' ws_minmax(c(5.2, 3.1, 8.7, 2.4))
#' # Matrix (samples x features)
#' ws_minmax(matrix(rnorm(100), nrow=10))
#' }
#'
#' @export
ws_minmax <- function(x) {
  if (is.vector(x) || (is.matrix(x) && nrow(x) == 1)) {
    x <- as.numeric(x)
    rng <- range(x, na.rm = TRUE)
    if (rng[2] == rng[1]) return(rep(0.5, length(x)))
    return((x - rng[1]) / (rng[2] - rng[1]))
  }

  t(apply(x, 1, function(row) {
    rng <- range(row, na.rm = TRUE)
    if (rng[2] == rng[1]) return(rep(0.5, length(row)))
    (row - rng[1]) / (rng[2] - rng[1])
  }))
}


#' @title Within-Sample Rank Normalization
#'
#' @description
#' Replaces each sample's feature values with their within-sample fractional
#' ranks. Completely ordinal and batch-immune.
#'
#' @param x A numeric matrix (samples x features) or vector
#'
#' @return Matrix with values in (0,1] representing fractional ranks
#'
#' @export
ws_rank <- function(x) {
  if (is.vector(x)) {
    return(rank(x, ties.method = "average") / length(x))
  }
  t(apply(x, 1, function(row) rank(row, ties.method = "average") / length(row)))
}


#' @title Within-Sample Z-Score Normalization
#'
#' @description
#' Standardizes each sample across its own features (mean=0, sd=1 within sample).
#' Batch-immune because it uses only within-sample statistics.
#'
#' @param x A numeric matrix (samples x features) or vector
#'
#' @return Matrix with within-sample standardized values
#'
#' @export
ws_zscore <- function(x) {
  if (is.vector(x)) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    return((x - m) / s)
  }
  t(apply(x, 1, function(row) {
    m <- mean(row, na.rm = TRUE)
    s <- sd(row, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(row)))
    (row - m) / s
  }))
}


#' @title Within-Sample Pairwise Log-Ratio Features
#'
#' @description
#' Computes all pairwise log-ratios between features within each sample.
#' For p features, produces p*(p-1)/2 ratio features. Each ratio is
#' inherently self-normalizing — multiplicative batch effects cancel.
#'
#' If the input is already on log scale (e.g., log2 expression), the
#' log-ratio is simply the difference: log(A/B) = log(A) - log(B).
#'
#' @param x A numeric matrix (samples x features). Assumed log-scale.
#' @param feature_names Optional character vector of feature names
#'
#' @return Matrix with p*(p-1)/2 columns (pairwise ratios)
#'
#' @export
ws_logratio <- function(x, feature_names = NULL) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)

  p <- ncol(x)
  n_ratios <- p * (p - 1) / 2
  out <- matrix(NA_real_, nrow = nrow(x), ncol = n_ratios)

  rnames <- character(n_ratios)
  if (is.null(feature_names)) feature_names <- paste0("F", seq_len(p))

  k <- 1
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      out[, k] <- x[, i] - x[, j]
      rnames[k] <- paste0(feature_names[i], "_vs_", feature_names[j])
      k <- k + 1
    }
  }
  colnames(out) <- rnames
  out
}


#' @title Within-Sample Pairwise Ratio Image
#'
#' @description
#' Creates a p x p pairwise ratio matrix for each sample, suitable as input
#' to a CNN. pixel(i,j) = feature_i - feature_j (on log scale = log-ratio).
#'
#' This encoding is analogous to color normalization in histopathology imaging:
#' batch effects (brightness/contrast) cancel in the pairwise ratios, just as
#' stain normalization preserves relative channel intensities.
#'
#' @param x A numeric matrix (samples x features) or vector (single sample).
#'          Assumed log-scale.
#'
#' @return A 3D array (samples x features x features) for matrix input,
#'         or a 2D matrix (features x features) for vector input.
#'
#' @examples
#' \dontrun{
#' # Single sample -> 14x14 ratio image
#' img <- ws_ratio_image(c(5.2, 3.1, 8.7, 2.4, 6.1, 4.3, 7.2,
#'                          1.9, 5.5, 9.0, 3.8, 4.7, 6.8, 2.1))
#' image(img, main = "miRNA ratio image")
#'
#' # Matrix -> array of images
#' imgs <- ws_ratio_image(matrix(rnorm(140), nrow = 10, ncol = 14))
#' dim(imgs)  # 10 x 14 x 14
#' }
#'
#' @export
ws_ratio_image <- function(x) {
  if (is.vector(x)) {
    return(encode_ratio_image(x))
  }

  make_ratio_images(x)
}


#' @title Apply Within-Sample Normalization to OmicPipeline
#'
#' @description
#' Creates a PipeOp for within-sample normalization that can be integrated
#' into the OmicSelector mlr3 pipeline. This ensures within-sample normalization
#' is applied inside CV folds with zero leakage.
#'
#' @param method One of "minmax", "rank", "zscore", "logratio", "ratio_image"
#' @param id PipeOp identifier
#'
#' @return A PipeOpTaskPreproc for within-sample normalization
#'
#' @export
create_ws_pipeop <- function(method = "minmax", id = "ws_norm") {

  ws_fn <- switch(method,
    minmax = ws_minmax,
    rank = ws_rank,
    zscore = ws_zscore,
    logratio = ws_logratio,
    stop("Unknown method: ", method, ". Use minmax, rank, zscore, or logratio.")
  )

  PipeOpWSNorm <- R6::R6Class(
    "PipeOpWSNorm",
    inherit = mlr3pipelines::PipeOpTaskPreproc,
    public = list(
      initialize = function(id = "ws_norm", param_vals = list()) {
        ps <- paradox::ps(
          method = paradox::p_fct(
            levels = c("minmax", "rank", "zscore", "logratio"),
            tags = "train"
          )
        )
        ps$values <- list(method = method)
        super$initialize(
          id = id,
          param_set = ps,
          param_vals = param_vals,
          feature_types = c("numeric", "integer")
        )
      }
    ),
    private = list(
      .train_dt = function(dt, levels, target) {
        fn <- switch(self$param_set$values$method,
          minmax = ws_minmax,
          rank = ws_rank,
          zscore = ws_zscore,
          logratio = ws_logratio
        )
        mat <- as.matrix(dt)
        storage.mode(mat) <- "numeric"
        result <- fn(mat)
        as.data.table(result)
      },
      .predict_dt = function(dt, levels) {
        fn <- switch(self$param_set$values$method,
          minmax = ws_minmax,
          rank = ws_rank,
          zscore = ws_zscore,
          logratio = ws_logratio
        )
        mat <- as.matrix(dt)
        storage.mode(mat) <- "numeric"
        result <- fn(mat)
        as.data.table(result)
      }
    )
  )

  PipeOpWSNorm$new(id = id)
}
