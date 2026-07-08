#' @title Safe Preprocessing Utilities for Cross-Validation
#'
#' @description
#' Functions that enforce leakage-free preprocessing within cross-validation folds.
#' These utilities prevent the most common data leakage bugs in biomarker validation:
#' global median imputation, global scaling, and test-fold self-imputation.
#'
#' @name safe-preprocessing
NULL


#' @title Within-Fold Median Imputation (Leakage-Free)
#'
#' @description
#' Imputes missing values using training-fold medians only. This prevents
#' the data leakage that occurs when medians are computed on the full dataset
#' (including test-fold samples) before the CV loop.
#'
#' **Common bug this prevents**: Computing `median(X[, j])` on the full matrix
#' before the train/test split, which leaks test-fold information into
#' the imputed training set.
#'
#' @param X_train Training fold feature matrix (samples x features)
#' @param X_test Test fold feature matrix (samples x features)
#'
#' @return A list with:
#'   \item{X_train}{Imputed training matrix (training-only medians)}
#'   \item{X_test}{Imputed test matrix (same training medians)}
#'   \item{impute_medians}{Named vector of medians used for imputation}
#'   \item{n_imputed_train}{Number of values imputed in training set}
#'   \item{n_imputed_test}{Number of values imputed in test set}
#'
#' @examples
#' \dontrun{
#' for (fold in folds) {
#'   X_train <- feature_mat[train_idx, ]
#'   X_test  <- feature_mat[test_idx, ]
#'
#'   # CORRECT: within-fold imputation
#'   imp <- impute_within_fold(X_train, X_test)
#'   X_train <- imp$X_train
#'   X_test  <- imp$X_test
#'
#'   # WRONG (leaks!): imputing before the loop
#'   # for (j in 1:ncol(mat)) mat[is.na(mat[,j]), j] <- median(mat[,j], na.rm=TRUE)
#' }
#' }
#'
#' @export
impute_within_fold <- function(X_train, X_test) {

  if (!is.matrix(X_train)) X_train <- as.matrix(X_train)
  if (!is.matrix(X_test))  X_test  <- as.matrix(X_test)

  if (ncol(X_train) != ncol(X_test)) {
    stop("X_train and X_test must have the same number of columns", call. = FALSE)
  }

  medians <- rep(NA_real_, ncol(X_train))
  names(medians) <- colnames(X_train)
  n_imp_train <- 0L
  n_imp_test  <- 0L

  for (j in seq_len(ncol(X_train))) {
    med <- median(X_train[, j], na.rm = TRUE)
    medians[j] <- med

    na_tr <- is.na(X_train[, j])
    if (any(na_tr)) {
      X_train[na_tr, j] <- med
      n_imp_train <- n_imp_train + sum(na_tr)
    }

    na_te <- is.na(X_test[, j])
    if (any(na_te)) {
      X_test[na_te, j] <- med
      n_imp_test <- n_imp_test + sum(na_te)
    }
  }

  list(
    X_train        = X_train,
    X_test         = X_test,
    impute_medians = medians,
    n_imputed_train = n_imp_train,
    n_imputed_test  = n_imp_test
  )
}


#' @title Detect Global Pre-Loop Imputation (Leakage Check)
#'
#' @description
#' Checks whether a feature matrix has already been imputed (no NAs remaining)
#' when NAs are expected. Call this at the top of a CV loop to detect accidental
#' global imputation.
#'
#' @param X Feature matrix to check
#' @param expected_na_rate Expected proportion of NAs (0 to 1). If the matrix
#'   has fewer NAs than expected, a warning is issued.
#' @param label Label for the matrix in warning messages (e.g., "panel_mat")
#'
#' @return Invisible TRUE if check passes, FALSE with warning if suspicious
#'
#' @export
check_no_premature_imputation <- function(X, expected_na_rate = 0, label = "X") {
  actual_rate <- sum(is.na(X)) / length(X)
  if (expected_na_rate > 0 && actual_rate == 0) {
    warning(
      sprintf(
        "Possible data leakage: '%s' has 0%% NAs but %.1f%% expected. ",
        label, expected_na_rate * 100
      ),
      "Check if global imputation was applied before the CV loop.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
}


#' @title Within-Fold Standardization (Leakage-Free)
#'
#' @description
#' Standardizes features using training-fold statistics only. Test-fold
#' features are centered and scaled using the same training-derived parameters.
#'
#' @param X_train Training fold feature matrix
#' @param X_test Test fold feature matrix
#'
#' @return A list with:
#'   \item{X_train}{Standardized training matrix}
#'   \item{X_test}{Standardized test matrix (using training means/sds)}
#'   \item{means}{Training-derived column means}
#'   \item{sds}{Training-derived column SDs (zeros replaced with 1)}
#'
#' @export
standardize_within_fold <- function(X_train, X_test) {

  if (!is.matrix(X_train)) X_train <- as.matrix(X_train)
  if (!is.matrix(X_test))  X_test  <- as.matrix(X_test)

  means <- colMeans(X_train, na.rm = TRUE)
  sds   <- apply(X_train, 2, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1

  X_train_s <- scale(X_train, center = means, scale = sds)
  X_test_s  <- scale(X_test,  center = means, scale = sds)

  list(
    X_train = X_train_s,
    X_test  = X_test_s,
    means   = means,
    sds     = sds
  )
}
