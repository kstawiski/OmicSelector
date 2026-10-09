.input_guard_fixture <- function() {
  set.seed(20261003)
  sbp <- ws_default_sbp()
  features <- .ws_balance_ilr_canonical_features(sbp)
  X <- matrix(rexp(60L * length(features)) + 1, nrow = 60L,
              dimnames = list(paste0("S", 1:60), features))
  y <- rep(0:1, each = 30L)
  X[y == 1L, 1:5] <- X[y == 1L, 1:5] + 20
  X[y == 0L, 6:10] <- X[y == 0L, 6:10] + 20
  list(X = X, y = y)
}

test_that("deployment rejects absent, empty, unnamed and nondetect input", {
  dat <- .input_guard_fixture()
  methods <- c("reo-ktsp", "ws-balance-ilr", "ws-rclr-trimmed")
  if (requireNamespace("glmnet", quietly = TRUE)) methods <- c(methods, "reo-pairratio")
  for (method in methods) {
    dep <- deploy_singlesample(dat$X, dat$y, method, verify = FALSE)
    state <- serialize(dep, NULL)
    rng <- .Random.seed
    x <- dat$X[1:3, , drop = FALSE]
    unknown <- matrix(10, 1, 2, dimnames = list(NULL, c("unknown1", "unknown2")))
    expect_error(score_specimen(dep, unknown), "ALL_REQUIRED_ABSENT")
    expect_error(score_specimen(dep, numeric()), "EMPTY_INPUT")
    expect_error(score_specimen(dep, x[, FALSE, drop = FALSE]), "EMPTY_INPUT")
    expect_error(score_specimen(dep, x[FALSE, , drop = FALSE]), "EMPTY_INPUT")
    expect_error(score_specimen(dep, as.numeric(x[1, ])), "INVALID_FEATURE_NAMES")
    for (blank in c(0, NA_real_, NaN, Inf, -Inf)) {
      bad <- x
      bad[2, ] <- blank
      err <- tryCatch(score_specimen(dep, bad), error = identity)
      expect_s3_class(err, "singlesample_input_error")
      code <- if (is.finite(blank)) "ALL_NONDETECT" else "NONFINITE_REQUIRED_VALUE"
      expect_identical(err$code, code)
      expect_identical(err$rows, 2L)
      expect_error(score_specimen(dep, bad[2, , drop = FALSE]), code)
    }
    bad <- x
    bad[, dep$input_guard$required_features] <- 0
    expect_error(score_specimen(dep, cbind(bad, unknown = 1e6)), "ALL_NONDETECT")
    expect_identical(state, serialize(dep, NULL))
    expect_identical(rng, .Random.seed)
  }
})

test_that("nondetect thresholds are frozen and masks survive positive imputation", {
  dat <- .input_guard_fixture()
  dep <- deploy_singlesample(dat$X, dat$y, "reo-ktsp", nondetect_limit = 1)
  x <- dat$X[1:3, , drop = FALSE]
  bad <- x
  bad[2, ] <- 1
  expect_error(score_specimen(dep, bad), "ALL_NONDETECT")
  mask <- matrix(TRUE, nrow(x), ncol(x), dimnames = dimnames(x))
  mask[2, ] <- FALSE
  expect_error(score_specimen(dep, x, detected = mask), "ALL_NONDETECT")
  expect_error(score_specimen(dep, x[2, ], detected = mask[2, ]), "ALL_NONDETECT")
  mask[,] <- TRUE
  expect_identical(score_specimen(dep, x, detected = mask), score_specimen(dep, x))
  expect_error(score_specimen(dep, x, detected = mask[, rev(seq_len(ncol(mask)))]),
               "INVALID_DETECTION_MASK")
  mask[1, 1] <- NA
  expect_error(score_specimen(dep, x, detected = mask), "INVALID_DETECTION_MASK")
  for (limit in list(NA_real_, -1, Inf, c(0, 1), "zero")) {
    expect_error(deploy_singlesample(dat$X, dat$y, nondetect_limit = limit),
                 "INVALID_NONDETECT_LIMIT")
  }
})

test_that("empty and constant fitted support fails closed", {
  dat <- .input_guard_fixture()
  X <- dat$X
  X[,] <- rep(seq_len(ncol(X)), each = nrow(X))
  for (method in c("reo-ktsp", "ws-balance-ilr")) {
    expect_error(deploy_singlesample(X, dat$y, method, verify = FALSE),
                 "CONSTANT_FITTED_SUPPORT|EMPTY_FITTED_SUPPORT")
  }
  if (requireNamespace("glmnet", quietly = TRUE)) {
    model <- fit_reo_pairratio(X, dat$y)
    expect_error(score_reo_pairratio(model, dat$X), "CONSTANT_FITTED_SUPPORT")
    expect_error(deploy_singlesample(X, dat$y, "reo-pairratio", verify = FALSE),
                 "CONSTANT_FITTED_SUPPORT")
    model <- fit_reo_pairratio(dat$X, dat$y)
    model$selected_pairs <- model$selected_pairs[FALSE, ]
    model$coefficients <- numeric()
    expect_error(score_reo_pairratio(model, dat$X), "EMPTY_FITTED_SUPPORT")
  }
  model <- fit_reo_ktsp(dat$X, dat$y)
  model$ktsp$pairs <- model$ktsp$pairs[FALSE, ]
  expect_error(score_reo_ktsp(model, dat$X), "EMPTY_FITTED_SUPPORT")
  model <- fit_ws_balance_ilr(dat$X, dat$y)
  expect_error(score_ws_balance_ilr(model, dat$X[, 1:3, drop = FALSE]),
               "NO_USABLE_SUPPORT")
  expect_error(score_ws_balance_ilr(model, dat$X[, 1, drop = FALSE]),
               "NO_USABLE_SUPPORT")
  noise <- dat$X
  colnames(noise) <- paste0("noise", seq_len(ncol(noise)))
  expect_error(deploy_singlesample(noise, dat$y, "ws-balance-ilr", verify = FALSE),
               "NO_USABLE_SUPPORT")
})

test_that("valid inputs preserve exact singleton/batch arithmetic and frozen state", {
  dat <- .input_guard_fixture()
  methods <- c("reo-ktsp", "ws-balance-ilr", "ws-rclr-trimmed")
  if (requireNamespace("glmnet", quietly = TRUE)) methods <- c(methods, "reo-pairratio")
  for (method in methods) {
    dep <- deploy_singlesample(dat$X, dat$y, method)
    # Constant held-out scores and uniform positive specimens are valid when
    # the fit had variable support. They must not trigger a batch variance gate.
    x <- dat$X[c(1, 31, 1), , drop = FALSE]
    x[3, ] <- 10
    state <- serialize(dep, NULL)
    rng <- .Random.seed
    batch <- score_specimen(dep, x)
    single <- vapply(seq_len(nrow(x)), function(i) score_specimen(dep, x[i, ]), numeric(1L))
    expect_identical(as.numeric(batch), single)
    expect_identical(batch, score_specimen(dep, as.data.frame(x)))
    expect_identical(batch, score_specimen(dep, x[, rev(seq_len(ncol(x))), drop = FALSE]))
    expect_identical(batch[c(3, 1, 2)], as.numeric(score_specimen(dep, x[c(3, 1, 2), , drop = FALSE])))
    expect_identical(as.numeric(score_specimen(dep, x[rep(1, 3), , drop = FALSE])), rep(batch[1], 3))
    expect_identical(batch, dep$score_fn(dep$model, x))
    expect_identical(state, serialize(dep, NULL))
    expect_identical(rng, .Random.seed)
  }
})
