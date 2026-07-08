#' @title CoDaCoRe Interfaces for Sparse Log-Contrast Classification
#'
#' @description
#' Thin wrapper around the official CoDaCoRe implementation when the
#' `codacore` package is available, with a deterministic sparse-balance
#' fallback for environments where the upstream TensorFlow backend is missing
#' or fails at runtime.
#'
#' @references
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017).
#' Microbiome datasets are compositional: and this is not optional.
#' \emph{Frontiers in Microbiology}, 8, 2224.
#'
#' Gordon-Rodriguez E, Susin A, McEwen JD, et al. (2022).
#' CoDaCoRe: learning sparse log-ratios for high-throughput sequencing data.
#' \emph{Bioinformatics}, 38(12), 3179-3186.
#'
#' @name codacore-interface
NULL


#' @title Train a CoDaCoRe Classifier
#'
#' @description
#' Fits an official CoDaCoRe model when the `codacore` package is available.
#' If the package is unavailable or the TensorFlow backend errors and
#' `backend = "auto"`, the function falls back to a deterministic sparse-balance
#' logistic model built from ranked CLR effects.
#'
#' @param X_train Numeric training matrix with samples in rows and features in
#'   columns.
#' @param y_train Binary outcome vector.
#' @param epochs Number of optimization epochs. Defaults to `30`.
#' @param lr Learning rate placeholder retained for API compatibility with the
#'   other Paper 1 learners. The fallback backend does not use it directly.
#' @param batch_size Mini-batch size passed to the official `codacore`
#'   optimizer. Defaults to `256`.
#' @param class_weight Logical; if `TRUE`, use inverse-frequency weighting in
#'   the fallback logistic model.
#' @param device Included for API consistency; currently ignored.
#' @param verbose Logical; print backend messages. Defaults to `TRUE`.
#' @param seed Integer random seed. Defaults to `42`.
#' @param positive Optional positive-class label for non-numeric outcomes.
#' @param input_scale Either `"log"` (default) or `"raw"`.
#' @param backend One of `"auto"`, `"codacore"`, or `"fallback"`.
#' @param max_balances Maximum number of sparse balances for the fallback model
#'   and `maxBaseLearners` for the official backend.
#' @param top_k Initial sparse balance size used by the fallback model.
#' @param cv_folds Number of internal cross-validation folds passed to the
#'   official backend.
#' @param fast Logical; passed through to the official `codacore` call.
#'
#' @return An object of class `"omicselector_codacore"` suitable for
#'   [predict_codacore()].
#'
#' @examples
#' X_train <- matrix(rnorm(14 * 20), nrow = 20)
#' y_train <- factor(rep(c("control", "case"), each = 10))
#' fit <- train_codacore(X_train, y_train, backend = "fallback", verbose = FALSE)
#' fit$backend
#'
#' @references
#' Gordon-Rodriguez E, Susin A, McEwen JD, et al. (2022).
#' CoDaCoRe: learning sparse log-ratios for high-throughput sequencing data.
#' \emph{Bioinformatics}, 38(12), 3179-3186.
#'
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017).
#' Microbiome datasets are compositional: and this is not optional.
#' \emph{Frontiers in Microbiology}, 8, 2224.
#'
#' @export
train_codacore <- function(X_train, y_train,
                           epochs = 30L,
                           lr = 0.001,
                           batch_size = 256L,
                           class_weight = TRUE,
                           device = NULL,
                           verbose = TRUE,
                           seed = 42L,
                           positive = NULL,
                           input_scale = c("log", "raw"),
                           backend = c("auto", "codacore", "fallback"),
                           max_balances = 5L,
                           top_k = 3L,
                           cv_folds = 5L,
                           fast = TRUE) {
  input_scale <- match.arg(input_scale)
  backend <- match.arg(backend)
  train <- .paper1_prepare_train_inputs(X_train, y_train, positive = positive)
  logged_train <- .codacore_logged_matrix(train$x, input_scale = input_scale)

  selected_backend <- backend
  fit <- NULL

  if (backend %in% c("auto", "codacore")) {
    fit <- tryCatch(
      .fit_codacore_official(
        x_train_log = logged_train,
        y_train = train$y,
        epochs = epochs,
        batch_size = batch_size,
        verbose = verbose,
        seed = seed,
        max_balances = max_balances,
        cv_folds = cv_folds,
        fast = fast
      ),
      error = function(e) {
        if (identical(backend, "codacore")) {
          stop(e$message, call. = FALSE)
        }
        if (isTRUE(verbose)) {
          message("Official codacore backend unavailable; using fallback sparse balance model.")
          message("Reason: ", conditionMessage(e))
        }
        NULL
      }
    )
    if (!is.null(fit)) {
      selected_backend <- "codacore"
    }
  }

  if (is.null(fit)) {
    selected_backend <- "fallback"
    fit <- .fit_codacore_fallback(
      x_train_log = logged_train,
      y_train = train$y,
      class_weight = class_weight,
      max_balances = max_balances,
      top_k = top_k
    )
  }

  structure(
    c(
      fit,
      list(
        backend = selected_backend,
        feature_names = colnames(train$x),
        class_names = train$class_names,
        positive = train$positive,
        input_scale = input_scale,
        seed = as.integer(seed),
        lr = lr,
        device = device
      )
    ),
    class = "omicselector_codacore"
  )
}


#' @title Predict with a CoDaCoRe Classifier
#'
#' @description
#' Generates positive-class probabilities for new samples from a fit created by
#' [train_codacore()].
#'
#' @param fit A fitted `"omicselector_codacore"` object.
#' @param X_new Numeric matrix of new samples.
#'
#' @return A numeric vector of positive-class probabilities.
#'
#' @examples
#' X_train <- matrix(rnorm(14 * 20), nrow = 20)
#' y_train <- factor(rep(c("control", "case"), each = 10))
#' fit <- train_codacore(X_train, y_train, backend = "fallback", verbose = FALSE)
#' predict_codacore(fit, X_train[1:3, , drop = FALSE])
#'
#' @references
#' Gordon-Rodriguez E, Susin A, McEwen JD, et al. (2022).
#' CoDaCoRe: learning sparse log-ratios for high-throughput sequencing data.
#' \emph{Bioinformatics}, 38(12), 3179-3186.
#'
#' @export
predict_codacore <- function(fit, X_new) {
  if (!inherits(fit, "omicselector_codacore")) {
    stop("fit must be an object returned by train_codacore().", call. = FALSE)
  }

  X_new <- .paper1_prepare_new_data(
    X_new,
    n_features = length(fit$feature_names),
    feature_names = fit$feature_names
  )
  logged_new <- .codacore_logged_matrix(X_new, input_scale = fit$input_scale)

  switch(
    fit$backend,
    codacore = .predict_codacore_official(fit, logged_new),
    fallback = .predict_codacore_fallback(fit, logged_new)
  )
}


#' @title mlr3 CoDaCoRe Learner
#'
#' @description
#' `mlr3` learner wrapping [train_codacore()] and [predict_codacore()]. It is
#' registered on package load as `classif.codacore`.
#'
#' @examples
#' \dontrun{
#' learner <- mlr3::lrn("classif.codacore", backend = "fallback", verbose = FALSE)
#' learner
#' }
#'
#' @references
#' Gordon-Rodriguez E, Susin A, McEwen JD, et al. (2022).
#' CoDaCoRe: learning sparse log-ratios for high-throughput sequencing data.
#' \emph{Bioinformatics}, 38(12), 3179-3186.
#'
#' @keywords internal
#' @noRd
LearnerClassifCoDaCoRe <- R6::R6Class(
  "LearnerClassifCoDaCoRe",
  inherit = mlr3::LearnerClassif,
  public = list(
    initialize = function() {
      ps <- paradox::ps(
        epochs = paradox::p_int(lower = 1L, default = 30L, tags = "train"),
        lr = paradox::p_dbl(lower = 1e-5, upper = 1, default = 0.001, tags = "train"),
        batch_size = paradox::p_int(lower = 1L, default = 256L, tags = "train"),
        class_weight = paradox::p_lgl(default = TRUE, tags = "train"),
        device = paradox::p_fct(
          levels = c("auto", "cpu", "cuda"),
          default = "auto",
          tags = "train"
        ),
        verbose = paradox::p_lgl(default = FALSE, tags = "train"),
        seed = paradox::p_int(lower = 1L, default = 42L, tags = "train"),
        input_scale = paradox::p_fct(
          levels = c("log", "raw"),
          default = "log",
          tags = "train"
        ),
        backend = paradox::p_fct(
          levels = c("auto", "codacore", "fallback"),
          default = "auto",
          tags = "train"
        ),
        max_balances = paradox::p_int(lower = 1L, default = 5L, tags = "train"),
        top_k = paradox::p_int(lower = 1L, default = 3L, tags = "train"),
        cv_folds = paradox::p_int(lower = 2L, default = 5L, tags = "train"),
        fast = paradox::p_lgl(default = TRUE, tags = "train")
      )

      super$initialize(
        id = "classif.codacore",
        param_set = ps,
        predict_types = c("response", "prob"),
        feature_types = c("integer", "numeric"),
        properties = "twoclass",
        packages = "stats",
        label = "CoDaCoRe"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pars <- self$param_set$get_values(tags = "train")
      x_train <- as.matrix(task$data(cols = task$feature_names))
      storage.mode(x_train) <- "numeric"
      self$model <- train_codacore(
        X_train = x_train,
        y_train = task$truth(),
        epochs = pars$epochs %||% 30L,
        lr = pars$lr %||% 0.001,
        batch_size = pars$batch_size %||% 256L,
        class_weight = pars$class_weight %||% TRUE,
        device = pars$device %||% "auto",
        verbose = pars$verbose %||% FALSE,
        seed = pars$seed %||% 42L,
        positive = task$positive,
        input_scale = pars$input_scale %||% "log",
        backend = pars$backend %||% "auto",
        max_balances = pars$max_balances %||% 5L,
        top_k = pars$top_k %||% 3L,
        cv_folds = pars$cv_folds %||% 5L,
        fast = pars$fast %||% TRUE
      )
      invisible(self$model)
    },
    .predict = function(task) {
      x_new <- as.matrix(task$data(cols = task$feature_names))
      storage.mode(x_new) <- "numeric"
      prob_positive <- predict_codacore(self$model, x_new)
      prob <- .paper1_prob_matrix(
        prob_positive = prob_positive,
        class_names = self$model$class_names,
        positive = self$model$positive
      )
      response <- .paper1_response_from_prob(
        prob_positive = prob_positive,
        class_names = self$model$class_names,
        positive = self$model$positive
      )
      mlr3::PredictionClassif$new(
        task = task,
        row_ids = task$row_ids,
        truth = task$truth(),
        prob = prob,
        response = factor(response, levels = self$model$class_names)
      )
    }
  )
)


#' @keywords internal
.codacore_logged_matrix <- function(x, input_scale) {
  input_scale <- match.arg(input_scale, c("log", "raw"))
  x <- .paper1_require_matrix(x, arg = "x")
  if (identical(input_scale, "log")) {
    return(x)
  }
  if (any(x <= 0, na.rm = TRUE)) {
    stop("Raw CoDaCoRe inputs must be strictly positive.", call. = FALSE)
  }
  log(x)
}


#' @keywords internal
.fit_codacore_official <- function(x_train_log, y_train, epochs, batch_size,
                                   verbose, seed, max_balances, cv_folds, fast) {
  if (!requireNamespace("codacore", quietly = TRUE)) {
    stop("Package 'codacore' is not installed.", call. = FALSE)
  }
  codacore_fun <- get("codacore", envir = asNamespace("codacore"))
  x_train_raw <- exp(x_train_log)

  fit <- do.call(
    codacore_fun,
    list(
      x = x_train_raw,
      y = y_train,
      logRatioType = "balances",
      maxBaseLearners = as.integer(max_balances),
      optParams = list(
        epochs = as.integer(epochs),
        batchSize = as.integer(batch_size)
      ),
      cvParams = list(
        numFolds = as.integer(cv_folds),
        maxCutoffs = ncol(x_train_raw)
      ),
      verbose = isTRUE(verbose),
      fast = isTRUE(fast)
    )
  )

  list(
    model = fit
  )
}


#' @keywords internal
.fit_codacore_fallback <- function(x_train_log, y_train, class_weight,
                                   max_balances, top_k) {
  clr_train <- clr_transform(x_train_log, input_scale = "log")
  effect <- colMeans(clr_train[y_train == 1L, , drop = FALSE]) -
    colMeans(clr_train[y_train == 0L, , drop = FALSE])

  pos_rank <- order(effect, decreasing = TRUE)
  neg_rank <- order(effect, decreasing = FALSE)
  n_balances <- min(
    as.integer(max_balances),
    max(1L, floor(ncol(x_train_log) / 2L))
  )

  balance_spec <- vector("list", n_balances)
  for (i in seq_len(n_balances)) {
    size <- min(as.integer(top_k) + i - 1L, floor(ncol(x_train_log) / 2L))
    num_idx <- pos_rank[seq_len(size)]
    den_idx <- neg_rank[seq_len(size)]
    balance_spec[[i]] <- list(
      numerator = num_idx,
      denominator = den_idx
    )
  }

  balance_matrix <- .build_balance_matrix(x_train_log, balance_spec)
  balance_df <- as.data.frame(balance_matrix)
  balance_df$.target <- y_train
  obs_weights <- rep(1, length(y_train))
  if (isTRUE(class_weight)) {
    n_pos <- sum(y_train == 1L)
    n_neg <- sum(y_train == 0L)
    obs_weights[y_train == 1L] <- n_neg / max(n_pos, 1L)
  }

  glm_fit <- stats::glm(
    .target ~ .,
    data = balance_df,
    family = stats::binomial(),
    weights = obs_weights
  )

  list(
    model = glm_fit,
    balance_spec = balance_spec
  )
}


#' @keywords internal
.build_balance_matrix <- function(log_x, balance_spec) {
  out <- vapply(
    seq_along(balance_spec),
    function(i) {
      spec <- balance_spec[[i]]
      rowMeans(log_x[, spec$numerator, drop = FALSE]) -
        rowMeans(log_x[, spec$denominator, drop = FALSE])
    },
    numeric(nrow(log_x))
  )

  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1L)
  }
  colnames(out) <- paste0("balance_", seq_len(ncol(out)))
  out
}


#' @keywords internal
.predict_codacore_official <- function(fit, logged_new) {
  x_new_raw <- exp(logged_new)
  pred <- tryCatch(
    stats::predict(fit$model, x_new_raw, asLogits = FALSE),
    error = function(...) {
      predict_fun <- get("predict.codacore", envir = asNamespace("codacore"))
      predict_fun(fit$model, x_new_raw, asLogits = FALSE)
    }
  )
  pred <- as.numeric(pred)
  pred
}


#' @keywords internal
.predict_codacore_fallback <- function(fit, logged_new) {
  balance_matrix <- .build_balance_matrix(logged_new, fit$balance_spec)
  as.numeric(
    stats::predict(fit$model, newdata = as.data.frame(balance_matrix), type = "response")
  )
}
