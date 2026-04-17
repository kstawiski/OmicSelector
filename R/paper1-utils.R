#' @title Paper 1 v2.2 Internal Utilities
#' @name paper1-utils
#' @keywords internal
NULL


#' @keywords internal
.paper1_default_seeds <- function() {
  c(42L, 7L, 2026L)
}


#' @keywords internal
.paper1_require_matrix <- function(x, arg = deparse(substitute(x))) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1L)
  }
  if (!is.matrix(x)) {
    stop(arg, " must be a numeric matrix or vector.", call. = FALSE)
  }
  storage.mode(x) <- "numeric"
  if (nrow(x) < 1L || ncol(x) < 2L) {
    stop(arg, " must contain at least one sample and two features.", call. = FALSE)
  }
  x
}


#' @keywords internal
.paper1_binary_label_info <- function(y, positive = NULL) {
  if (is.factor(y)) {
    class_names <- levels(y)
    y_chr <- as.character(y)
  } else if (is.character(y)) {
    y_chr <- y
    class_names <- unique(y_chr)
  } else if (is.logical(y)) {
    y_chr <- ifelse(y, "TRUE", "FALSE")
    class_names <- c("FALSE", "TRUE")
  } else if (is.numeric(y) || is.integer(y)) {
    uniq <- sort(unique(as.integer(y)))
    if (!identical(uniq, c(0L, 1L))) {
      stop("Numeric labels must be binary 0/1.", call. = FALSE)
    }
    return(list(
      y = as.integer(y),
      class_names = c("0", "1"),
      positive = "1",
      negative = "0"
    ))
  } else {
    stop("Unsupported label type for binary classification.", call. = FALSE)
  }

  class_names <- unique(class_names)
  if (length(class_names) != 2L) {
    stop("Exactly two outcome classes are required.", call. = FALSE)
  }
  positive <- positive %||% class_names[[2L]]
  if (!positive %in% class_names) {
    stop("positive class '", positive, "' is not present in y.", call. = FALSE)
  }
  negative <- setdiff(class_names, positive)[[1L]]
  list(
    y = as.integer(y_chr == positive),
    class_names = class_names,
    positive = positive,
    negative = negative
  )
}


#' @keywords internal
.paper1_prepare_train_inputs <- function(X_train, y_train, positive = NULL) {
  X_train <- .paper1_require_matrix(X_train, arg = "X_train")
  if (is.null(colnames(X_train))) {
    colnames(X_train) <- paste0("feature_", seq_len(ncol(X_train)))
  }
  if (nrow(X_train) != length(y_train)) {
    stop("X_train and y_train must have the same number of samples.", call. = FALSE)
  }
  labels <- .paper1_binary_label_info(y_train, positive = positive)
  list(
    x = X_train,
    y = labels$y,
    class_names = labels$class_names,
    positive = labels$positive,
    negative = labels$negative
  )
}


#' @keywords internal
.paper1_prepare_new_data <- function(X_new, n_features, feature_names = NULL, arg = "X_new") {
  X_new <- .paper1_require_matrix(X_new, arg = arg)
  if (ncol(X_new) != n_features) {
    stop(
      arg, " must contain ", n_features, " features; got ", ncol(X_new), ".",
      call. = FALSE
    )
  }
  if (!is.null(feature_names)) {
    if (is.null(colnames(X_new))) {
      colnames(X_new) <- feature_names
    } else {
      X_new <- X_new[, feature_names, drop = FALSE]
    }
  }
  X_new
}


#' @keywords internal
.paper1_resolve_device <- function(device) {
  if (is.null(device) || identical(device, "auto")) {
    if (requireNamespace("torch", quietly = TRUE) &&
        tryCatch(torch::cuda_is_available(), error = function(...) FALSE)) {
      return("cuda")
    }
    return("cpu")
  }
  match.arg(device, choices = c("cpu", "cuda"))
}


#' @keywords internal
.paper1_set_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  set.seed(as.integer(seed))
  if (requireNamespace("torch", quietly = TRUE)) {
    try(torch::torch_manual_seed(as.integer(seed)), silent = TRUE)
    if (tryCatch(torch::cuda_is_available(), error = function(...) FALSE)) {
      torch_ns <- asNamespace("torch")
      cuda_manual_seed <- get0("cuda_manual_seed", envir = torch_ns, mode = "function")
      cuda_manual_seed_all <- get0("cuda_manual_seed_all", envir = torch_ns, mode = "function")
      if (is.function(cuda_manual_seed)) {
        try(cuda_manual_seed(as.integer(seed)), silent = TRUE)
      }
      if (is.function(cuda_manual_seed_all)) {
        try(cuda_manual_seed_all(as.integer(seed)), silent = TRUE)
      }
    }
  }
  invisible(seed)
}


#' @keywords internal
.paper1_normalize_train_test <- function(train, test = NULL) {
  mean_train <- mean(train)
  sd_train <- stats::sd(as.vector(train))
  if (is.na(sd_train) || sd_train <= 0) {
    return(list(
      train = train,
      test = test,
      mean = mean_train,
      sd = sd_train
    ))
  }
  list(
    train = (train - mean_train) / sd_train,
    test = if (is.null(test)) NULL else (test - mean_train) / sd_train,
    mean = mean_train,
    sd = sd_train
  )
}


#' @keywords internal
.paper1_prob_matrix <- function(prob_positive, class_names, positive) {
  negative <- setdiff(class_names, positive)[[1L]]
  prob <- matrix(
    0,
    nrow = length(prob_positive),
    ncol = 2L,
    dimnames = list(NULL, class_names)
  )
  prob[, positive] <- prob_positive
  prob[, negative] <- 1 - prob_positive
  prob
}


#' @keywords internal
.paper1_response_from_prob <- function(prob_positive, class_names, positive, threshold = 0.5) {
  negative <- setdiff(class_names, positive)[[1L]]
  ifelse(prob_positive >= threshold, positive, negative)
}


#' @keywords internal
.paper1_shift_scale <- function(x,
                                condition = c("baseline", "shift", "scale", "both"),
                                seed = 42L,
                                shift_sd = 2,
                                scale_range = c(0.5, 2.0)) {
  condition <- match.arg(condition)
  x <- .paper1_require_matrix(x, arg = "x")
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  set.seed(as.integer(seed))
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (condition == "baseline") {
    return(x)
  }

  shifted <- x
  if (condition %in% c("shift", "both")) {
    shift <- stats::rnorm(nrow(x), mean = 0, sd = shift_sd)
    shifted <- shifted + shift
  }
  if (condition %in% c("scale", "both")) {
    scale <- stats::runif(nrow(x), min = scale_range[[1L]], max = scale_range[[2L]])
    shifted <- shifted * scale
  }
  shifted
}


#' @keywords internal
.paper1_auc_summary <- function(prediction, measure) {
  tryCatch(
    prediction$score(measure),
    error = function(...) NA_real_
  )
}


#' @keywords internal
.paper1_register_learner <- function(key, constructor) {
  if (!requireNamespace("mlr3", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  dictionary <- mlr3::mlr_learners
  if (key %in% dictionary$keys()) {
    return(invisible(FALSE))
  }
  try(dictionary$add(key, constructor), silent = TRUE)
  invisible(TRUE)
}

