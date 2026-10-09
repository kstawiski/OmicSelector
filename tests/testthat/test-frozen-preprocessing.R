.frozen_fixture <- function() {
  set.seed(101)
  features <- OmicSelector:::.ws_balance_ilr_canonical_features(ws_default_sbp())[1:30]
  X <- matrix(round(stats::rgamma(40 * 30, shape = 20, rate = 2)), 40,
              dimnames = list(paste0("S", 1:40), features))
  y <- rep(0:1, each = 20)
  X[y == 1, 1:5] <- X[y == 1, 1:5] + 80
  X[y == 0, 6:10] <- X[y == 0, 6:10] + 80
  list(X = X, y = y)
}

.frozen_expect_code <- function(expr, code, rows = integer()) {
  err <- tryCatch(force(expr), error = identity)
  expect_s3_class(err, "singlesample_input_error")
  expect_identical(err$code, code)
  expect_identical(err$rows, rows)
}

test_that("count conversion uses only each row and the ordered frozen universe", {
  X <- rbind(s1 = c(a = 1, b = 3, extra = 100),
             s2 = c(a = 0, b = 0, extra = 200),
             s3 = c(a = 0, b = 2, extra = 300))
  prep <- fit_singlesample_preprocessor(X, features = c("b", "a"))
  expect_s3_class(prep, "singlesample_preprocessor")
  expect_identical(names(prep), c("type", "features", "ct_limit", "ct_base", "scale"))
  expect_identical(prep, fit_singlesample_preprocessor(X * 100,
                                                      features = c("b", "a")))
  state <- serialize(prep, NULL)
  out <- apply_singlesample_preprocessor(prep, X)
  expected <- X[, c("b", "a")]
  expected[1, ] <- c(750000, 250000)
  expected[3, ] <- c(1e6, 0)
  expect_identical(out$abundance, expected)
  expect_identical(out$detected, X[, c("b", "a")] > 0)
  for (i in 1:3) {
    expect_identical(apply_singlesample_preprocessor(prep, X[i, , drop = FALSE]),
                     lapply(out, function(m) m[i, , drop = FALSE]))
  }
  expect_identical(out, apply_singlesample_preprocessor(prep, as.data.frame(X)))
  expect_identical(out, apply_singlesample_preprocessor(prep, X[, 3:1]))
  partial <- apply_singlesample_preprocessor(prep, X[1, c("a", "extra")])
  expect_identical(colnames(partial$abundance), "a")
  expect_identical(as.numeric(partial$abundance), 1e6)
  absent <- apply_singlesample_preprocessor(prep, X[, "extra", drop = FALSE])
  expect_identical(dim(absent$abundance), c(3L, 0L))
  expect_identical(dim(absent$detected), c(3L, 0L))
  raw <- fit_singlesample_preprocessor(X, scale = "none")
  expect_identical(apply_singlesample_preprocessor(raw, X)$abundance, X)
  expect_identical(serialize(prep, NULL), state)
  expect_output(print(prep), "scale: cpm")
  expect_identical(withVisible(print(prep))$visible, FALSE)
})

test_that("Ct conversion handles markers, limits and non-finite measurements", {
  X <- rbind(s1 = c(a = "38", b = "40", c = "Undetermined", d = "Inf"),
             s2 = c(a = "NA", b = "-Inf", c = "NaN", d = "41"),
             s3 = c(a = "39", b = "37", c = NA, d = "36"))
  prep <- fit_singlesample_preprocessor(X, "ct")
  out <- apply_singlesample_preprocessor(prep, X)
  expected <- matrix(c(4, 1, 2, 1, 1, 8, 1, 1, 1, 1, 1, 16), 3,
                     dimnames = dimnames(X))
  expect_identical(out$abundance, expected)
  expect_identical(out$detected, expected > 1)
  for (i in 1:3) {
    expect_identical(apply_singlesample_preprocessor(prep, X[i, , drop = FALSE]),
                     lapply(out, function(m) m[i, , drop = FALSE]))
  }
  factors <- as.data.frame(X, stringsAsFactors = TRUE)
  expect_identical(apply_singlesample_preprocessor(prep, factors), out)
  changed <- fit_singlesample_preprocessor(X, "ct", ct_limit = 39, ct_base = 10,
                                           scale = "none")
  expect_identical(as.numeric(apply_singlesample_preprocessor(changed, X)$abundance[1, ]),
                   c(10, 1, 1, 1))
  expect_output(print(prep), "Ct limit: 40; base: 2")
})

test_that("invalid raw measurements and frozen parameters have coded errors", {
  X <- rbind(c(a = 1, b = 2), c(a = 3, b = 4))
  prep <- fit_singlesample_preprocessor(X)
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- X
    bad[2, 1] <- value
    .frozen_expect_code(apply_singlesample_preprocessor(prep, bad),
                        "NONFINITE_REQUIRED_VALUE", 2L)
  }
  bad <- X
  bad[2, 1] <- -1
  .frozen_expect_code(fit_singlesample_preprocessor(bad), "INVALID_COUNT_VALUE", 2L)
  .frozen_expect_code(apply_singlesample_preprocessor(prep, X * .Machine$double.xmax),
                      "NONFINITE_REQUIRED_VALUE", 1:2)
  bad <- X
  bad[2, ] <- .Machine$double.xmax
  .frozen_expect_code(apply_singlesample_preprocessor(prep, bad), "NONFINITE_ABUNDANCE", 2L)
  ct <- fit_singlesample_preprocessor(X, "ct")
  bad[2, 1] <- -10000
  .frozen_expect_code(apply_singlesample_preprocessor(ct, bad), "NONFINITE_ABUNDANCE", 2L)
  for (value in list(NA_real_, Inf, 0, -1, "40", c(30, 40))) {
    .frozen_expect_code(fit_singlesample_preprocessor(X, ct_limit = value),
                        "INVALID_PREPROCESSOR")
  }
  for (value in list(NA_real_, Inf, 0, 1, "2", c(2, 3))) {
    .frozen_expect_code(fit_singlesample_preprocessor(X, ct_base = value),
                        "INVALID_PREPROCESSOR")
  }
  for (features in list(NULL, c("a", "a"), c("a", ""), NA_character_, "missing")) {
    .frozen_expect_code(fit_singlesample_preprocessor(X, features = features),
                        "INVALID_FEATURE_NAMES")
  }
  .frozen_expect_code(apply_singlesample_preprocessor(list(), X), "INVALID_PREPROCESSOR")
  .frozen_expect_code(apply_singlesample_preprocessor(prep, unname(X)), "INVALID_FEATURE_NAMES")
  .frozen_expect_code(apply_singlesample_preprocessor(prep, X[, FALSE, drop = FALSE]), "EMPTY_INPUT")
  .frozen_expect_code(apply_singlesample_preprocessor(prep, X[FALSE, , drop = FALSE]), "EMPTY_INPUT")
  .frozen_expect_code(apply_singlesample_preprocessor(prep, matrix("bad", 1, 1)), "INVALID_INPUT")
})

for (.frozen_type in c("counts", "ct")) {
  for (.frozen_method in c("reo-ktsp", "reo-pairratio", "ws-rclr-trimmed")) {
    local({
      type <- .frozen_type
      method <- .frozen_method
      test_that(paste(type, method, "is raw-row equivariant and matches manual conversion"), {
        if (method == "reo-pairratio") skip_if_not_installed("glmnet")
        dat <- .frozen_fixture()
        raw <- if (type == "counts") dat$X else 40 - log2(dat$X)
        prep <- fit_singlesample_preprocessor(raw, type)
        dep <- deploy_singlesample(raw, dat$y, method, preprocessor = prep)
        expect_true(dep$verified)
        expect_identical(dep$preprocessor, prep)
        manual <- if (type == "counts") (raw / rowSums(raw)) * 1e6 else 2^(40 - raw)
        plain <- deploy_singlesample(manual, dat$y, method)
        explicit_null <- deploy_singlesample(manual, dat$y, method, preprocessor = NULL)
        expect_identical(score_specimen(dep, raw), score_specimen(plain, manual))
        expect_identical(score_specimen(explicit_null, manual), score_specimen(plain, manual))
        expect_identical(score_specimen(plain, manual), plain$score_fn(plain$model, manual))
        state <- serialize(dep, NULL)
        rng <- .Random.seed
        extreme <- if (type == "ct") function(X) X - log2(1e6) else NULL
        for (strict in c(FALSE, TRUE)) {
          score_fun <- function(model, X, meta = NULL) score_specimen(model, X, meta, strict = strict)
          expect_invisible(singlesample_assert_row_equivariant(
            score_fun, dep, raw, tol = 0, n_singletons = nrow(raw),
            extreme_transform = extreme))
          batch <- score_specimen(dep, raw, strict = strict)
          singles <- vapply(seq_len(nrow(raw)), function(i)
            as.numeric(score_specimen(dep, raw[i, ], strict = strict)), numeric(1))
          expect_identical(as.numeric(batch), singles)
          expect_identical(batch, score_specimen(dep, as.data.frame(raw), strict = strict))
          expect_identical(batch, score_specimen(dep, raw[, ncol(raw):1], strict = strict))
          expect_identical(batch, score_specimen(dep, cbind(raw, unknown = 1e9), strict = strict))
          absent <- matrix(100, 2, 1, dimnames = list(NULL, "unknown"))
          .frozen_expect_code(score_specimen(dep, absent, strict = strict), "ALL_REQUIRED_ABSENT")
          bad <- raw[1:3, , drop = FALSE]
          bad[2, ] <- if (type == "counts") 0 else NA_real_
          .frozen_expect_code(score_specimen(dep, bad, strict = strict), "ALL_NONDETECT", 2L)
        }
        expect_true(is_singlesample_deployable(dep, raw))
        expect_identical(serialize(dep, NULL), state)
        expect_identical(.Random.seed, rng)
      })
    })
  }
}

test_that("raw Ct NA values convert and detection masks combine by logical AND", {
  dat <- .frozen_fixture()
  ct <- 40 - log2(dat$X)
  prep <- fit_singlesample_preprocessor(ct, "ct")
  dep <- deploy_singlesample(ct, dat$y, "reo-ktsp", preprocessor = prep)
  x <- ct[1:3, , drop = FALSE]
  required <- dep$input_guard$required_features
  x[2, required[1]] <- NA
  for (strict in c(FALSE, TRUE)) {
    expect_true(all(is.finite(score_specimen(dep, x, strict = strict))))
    manual <- apply_singlesample_preprocessor(prep, x)
    plain <- dep
    plain$preprocessor <- NULL
    expect_identical(score_specimen(dep, x, strict = strict),
      score_specimen(plain, manual$abundance, detected = manual$detected, strict = strict))
    bad <- x
    bad[2, required] <- NA
    mask <- matrix(TRUE, nrow(bad), ncol(bad), dimnames = dimnames(bad))
    .frozen_expect_code(score_specimen(dep, bad, detected = mask, strict = strict), "ALL_NONDETECT", 2L)
    mask[3, required] <- FALSE
    .frozen_expect_code(score_specimen(dep, x, detected = mask, strict = strict), "ALL_NONDETECT", 3L)
    .frozen_expect_code(score_specimen(dep, x[3, ], detected = mask[3, ], strict = strict),
                        "ALL_NONDETECT", 1L)
  }
  chars <- as.data.frame(x)
  chars[[required[1]]][2] <- "Undetermined"
  expect_identical(score_specimen(dep, chars), score_specimen(dep, x))
  for (value in list(40, Inf, NaN, "Undetermined")) {
    bad <- x
    bad[2, ] <- value
    .frozen_expect_code(score_specimen(dep, bad), "ALL_NONDETECT", 2L)
  }
  bad_train <- ct
  bad_train[2, ] <- NA
  .frozen_expect_code(deploy_singlesample(bad_train, dat$y, "reo-ktsp", preprocessor = prep),
                      "ALL_NONDETECT", 2L)
  ct_train <- ct
  ct_train[2, 1] <- NA
  expect_true(deploy_singlesample(ct_train, dat$y, "reo-ktsp", preprocessor = prep)$verified)
  ct_train <- as.data.frame(ct_train)
  ct_train[2, 1] <- "Undetermined"
  expect_true(deploy_singlesample(ct_train, dat$y, "reo-ktsp", preprocessor = prep)$verified)
  with_extra <- cbind(x, unknown = 100)
  mask <- matrix(TRUE, nrow(with_extra), ncol(with_extra), dimnames = dimnames(with_extra))
  expect_identical(score_specimen(dep, with_extra, detected = mask), score_specimen(dep, x))
  .frozen_expect_code(score_specimen(dep, with_extra, detected = mask[, ncol(mask):1]),
                      "INVALID_DETECTION_MASK")
  mask[1, 1] <- NA
  .frozen_expect_code(score_specimen(dep, with_extra, detected = mask), "INVALID_DETECTION_MASK")
})

test_that("required non-finite abundance is rejected on every deployment route", {
  dat <- .frozen_fixture()
  methods <- c("reo-ktsp", "ws-rclr-trimmed", "ws-balance-ilr")
  if (requireNamespace("glmnet", quietly = TRUE)) methods <- c(methods, "reo-pairratio")
  for (method in methods) {
    dep <- deploy_singlesample(dat$X, dat$y, method, verify = FALSE)
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      bad <- dat$X[1:3, , drop = FALSE]
      bad[c(1, 3), dep$input_guard$required_features[1]] <- value
      for (strict in c(FALSE, TRUE)) {
        .frozen_expect_code(score_specimen(dep, bad, strict = strict), "NONFINITE_REQUIRED_VALUE", c(1L, 3L))
        .frozen_expect_code(dep$score_fn(dep$model, bad, strict = strict), "NONFINITE_REQUIRED_VALUE", c(1L, 3L))
        .frozen_expect_code(singlesample_score_call(method, dep$model, bad, strict = strict),
                            "NONFINITE_REQUIRED_VALUE", c(1L, 3L))
      }
    }
  }
  dep <- deploy_singlesample(dat$X, dat$y, "reo-ktsp", verify = FALSE)
  for (route in c("score_call", "direct_model", "direct_x")) {
    dep$score_route <- route
    dep$score_fn <- function(...) stop("scorer must not be called")
    bad <- dat$X[1:3, , drop = FALSE]
    bad[2, dep$input_guard$required_features[1]] <- NA_real_
    .frozen_expect_code(score_specimen(dep, bad), "NONFINITE_REQUIRED_VALUE", 2L)
  }
  model <- structure(list(panel = "a", weights = 1, feature_universe = c("a", "b")),
                     class = "ws_rclr_panel_model")
  bad <- rbind(c(a = 2, b = Inf), c(a = 1, b = 3))
  .frozen_expect_code(score_ws_rclr_panel(model, bad, strict = TRUE), "NONFINITE_REQUIRED_VALUE", 1L)
})

test_that("domain-specific extreme probes still detect batch coupling", {
  X <- matrix(seq(25, 35, length.out = 24), 6, dimnames = list(NULL, letters[1:4]))
  extreme <- function(X) X - log2(1e6)
  coupled <- function(model, X, meta = NULL) rowSums(X) - mean(X)
  expect_error(singlesample_assert_row_equivariant(coupled, NULL, X,
               extreme_transform = extreme), "row-equivariance violated")
  local_score <- function(model, X, meta = NULL) rowSums(X)
  expect_error(singlesample_assert_row_equivariant(local_score, NULL, X,
               extreme_transform = function(X) X[, 1, drop = FALSE]), "preserve numeric matrix")
  expect_error(singlesample_assert_row_equivariant(local_score, NULL, X,
               extreme_transform = 1), "NULL or a function")
})
