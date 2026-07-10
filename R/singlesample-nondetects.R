#' @title Single-sample qPCR non-detect imputation (Module B add-on)
#'
#' @description
#' Two qPCR non-detect handlers used by the single-sample pipeline:
#' \itemize{
#'   \item \code{\link{qpcr_nondetect_impute}}: Bayesian hierarchical
#'     imputation via the \code{nondetects} Bioconductor package
#'     (McCall et al. 2014).
#'   \item \code{\link{qpcr_nondetect_lod_fallback}}: limit-of-detection
#'     fallback that fills undetermined Cts with a fixed value (default 40).
#' }
#'
#' Both operate on a features-by-samples Ct matrix. NAs in the input
#' represent undetermined Cts. The Bayesian path requires both
#' \code{nondetects} and \code{HTqPCR} (Bioconductor); if either is
#' unavailable, the function falls back to LOD imputation with a logged
#' message. This mirrors the single-sample non-detect handling convention.
#'
#' @references
#' McCall MN, McMurray HR, Land H, Almudevar A. (2014) On non-detects in qPCR
#' data. \emph{Bioinformatics} 30(16): 2310-2316.
#'
#' @name singlesample-nondetects
NULL


#' @title Bayesian hierarchical imputation of qPCR non-detects
#'
#' @description
#' Wraps \code{nondetects::qpcrImpute} for Ct-scale qPCR expression matrices.
#' Falls back to LOD imputation (\code{\link{qpcr_nondetect_lod_fallback}})
#' when either \code{nondetects} or \code{HTqPCR} is unavailable, or when
#' the Bayesian fit fails. Sets \code{nondetect_method} and
#' \code{n_imputed} attributes on the returned matrix.
#'
#' @param expr_mat Numeric matrix (features-by-samples) of Ct values; NAs
#'   represent non-detects.
#' @param log_handler Function used to emit progress messages.
#'   Default \code{message}.
#'
#' @return Matrix of the same shape as \code{expr_mat} with non-detects
#'   imputed. Carries attributes \code{nondetect_method} (one of
#'   \code{"nondetects_qpcrImpute"}, \code{"lod_fallback"},
#'   \code{"none_needed"}) and \code{n_imputed} (integer count of
#'   imputed cells).
#'
#' @examples
#' \dontrun{
#' m <- matrix(c(20, NA, 25, 30, NA, NA), nrow = 2,
#'             dimnames = list(c("miR-1", "miR-2"), c("S1", "S2", "S3")))
#' qpcr_nondetect_impute(m)
#' }
#' @export
qpcr_nondetect_impute <- function(expr_mat, log_handler = message) {
  if (!requireNamespace("nondetects", quietly = TRUE)) {
    log_handler("nondetects package not available; falling back to LOD imputation")
    return(qpcr_nondetect_lod_fallback(expr_mat, log_handler = log_handler))
  }
  if (!requireNamespace("HTqPCR", quietly = TRUE)) {
    log_handler("HTqPCR package not available; falling back to LOD imputation")
    return(qpcr_nondetect_lod_fallback(expr_mat, log_handler = log_handler))
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    log_handler("Biobase package not available; falling back to LOD imputation")
    return(qpcr_nondetect_lod_fallback(expr_mat, log_handler = log_handler))
  }

  n_na_before <- sum(is.na(expr_mat))
  if (n_na_before == 0L) {
    log_handler("no non-detects found; skipping imputation")
    attr(expr_mat, "nondetect_method") <- "none_needed"
    attr(expr_mat, "n_imputed") <- 0L
    return(expr_mat)
  }

  log_handler(sprintf("nondetects: %d NAs in %dx%d matrix (%.1f%%)",
                      n_na_before, nrow(expr_mat), ncol(expr_mat),
                      100 * n_na_before / length(expr_mat)))

  result <- tryCatch({
    qpcr_obj <- methods::new(
      "qPCRset",
      exprs = expr_mat,
      featureCategory = ifelse(rowSums(is.na(expr_mat)) > 0L,
                               "Undetermined", "OK")
    )
    imputed <- nondetects::qpcrImpute(qpcr_obj, groupVars = NULL, iterMax = 100)
    Biobase::exprs(imputed)
  }, error = function(e) {
    log_handler(sprintf("nondetects::qpcrImpute failed: %s -- falling back to LOD", e$message))
    NULL
  })

  if (is.null(result)) {
    return(qpcr_nondetect_lod_fallback(expr_mat, log_handler = log_handler))
  }

  n_na_after <- sum(is.na(result))
  log_handler(sprintf("nondetects imputation: %d -> %d NAs", n_na_before, n_na_after))
  attr(result, "nondetect_method") <- "nondetects_qpcrImpute"
  attr(result, "n_imputed") <- n_na_before - n_na_after
  result
}


#' @title Limit-of-detection fallback imputation for qPCR non-detects
#'
#' @description
#' Replaces every NA in \code{expr_mat} with the supplied \code{lod_ct}
#' value (default 40). Used by the single-sample pipeline when Bayesian
#' hierarchical imputation via \code{nondetects} is unavailable or fails.
#' Sets \code{nondetect_method = "lod_fallback"} and \code{n_imputed}
#' attributes.
#'
#' @param expr_mat Numeric matrix of Ct values; NAs represent non-detects.
#' @param lod_ct Numeric scalar; the Ct value used to impute NAs.
#'   Default 40.
#' @param log_handler Function used to emit progress messages. Default
#'   \code{message}.
#'
#' @return Matrix of the same shape as \code{expr_mat}, with NAs filled by
#'   \code{lod_ct}. Carries attributes \code{nondetect_method} and
#'   \code{n_imputed}.
#'
#' @examples
#' m <- matrix(c(20, NA, 25, 30, NA, NA), nrow = 2,
#'             dimnames = list(c("miR-1", "miR-2"), c("S1", "S2", "S3")))
#' qpcr_nondetect_lod_fallback(m, log_handler = function(...) invisible())
#' @export
qpcr_nondetect_lod_fallback <- function(expr_mat, lod_ct = 40, log_handler = message) {
  n_na <- sum(is.na(expr_mat))
  expr_mat[is.na(expr_mat)] <- lod_ct
  log_handler(sprintf("LOD fallback: imputed %d non-detects with Ct=%g", n_na, lod_ct))
  attr(expr_mat, "nondetect_method") <- "lod_fallback"
  attr(expr_mat, "n_imputed") <- n_na
  expr_mat
}
