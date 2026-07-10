#' @title Single-sample preprocessing utilities
#'
#' @description
#' Preprocessing helpers introduced in the OmicSelector single-sample scoring
#' bank. Currently provides inverse-log
#' preprocessing for microarray deposits whose values are already on a
#' log-abundance scale (e.g., Toray 3D-Gene, Agilent log-ratio exports).
#'
#' @name singlesample-preprocessing
NULL


# ----------------------------------------------------------------------------
# preprocess_inverse_log
# ----------------------------------------------------------------------------

#' @title Inverse-log preprocessor for pre-log-transformed microarray deposits
#'
#' @description
#' Converts a matrix of log-transformed expression intensities back to a
#' positive abundance-like scale via \eqn{base^x}. This is required when
#' a public deposit has already applied a \eqn{\log_2} (or \eqn{\log_{10}})
#' transform before archival (common in Toray 3D-Gene and some Agilent
#' deposits). The resulting matrix is on a positive scale comparable to raw
#' NGS counts and can be fed directly to within-sample Module A methods
#' (rCLR, ILR, ALR-pivot, MAD log-ratio, dominance score) which require
#' positive compositional input.
#'
#' Non-finite values (Inf, -Inf, NaN) are replaced with the per-row median
#' of finite values, with a warning.
#'
#' @param x Numeric matrix (features \eqn{\times} samples, or samples
#'   \eqn{\times} features — the function operates element-wise). Both
#'   orientations are accepted; the transform is purely element-wise.
#' @param base Numeric. The log base used by the deposit. Default 2
#'   (Toray 3D-Gene standard). Use 10 for log10 deposits.
#' @param pseudo Numeric. Small offset added to the output to keep strictly
#'   positive values. Default 0 (output of \eqn{base^x} is always
#'   \eqn{> 0} for finite inputs).
#'
#' @return Numeric matrix of the same dimensions as \code{x}, on a positive
#'   abundance scale. The attribute \code{preprocessing} records the
#'   transform applied.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' log2_matrix <- matrix(rnorm(50 * 20), nrow = 50)
#' colnames(log2_matrix) <- paste0("S", seq_len(20))
#' rownames(log2_matrix) <- paste0("miR-", seq_len(50))
#' abundance <- preprocess_inverse_log(log2_matrix, base = 2)
#' stopifnot(all(abundance > 0))
#' }
#'
#' @references
#' Toray Industries, Inc. (3D-Gene microRNA chip) — data format specification.
#'
#' @export
preprocess_inverse_log <- function(x, base = 2, pseudo = 0) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (any(!is.finite(x))) {
    warning("preprocess_inverse_log: ", sum(!is.finite(x)),
             " non-finite values present; replacing with feature row median")
    for (i in seq_len(nrow(x))) {
      bad <- !is.finite(x[i, ])
      if (any(bad) && any(!bad)) {
        x[i, bad] <- stats::median(x[i, !bad], na.rm = TRUE)
      }
    }
  }
  out <- base ^ x + pseudo
  attr(out, "preprocessing") <- sprintf("inverse_log_base_%g", base)
  out
}
