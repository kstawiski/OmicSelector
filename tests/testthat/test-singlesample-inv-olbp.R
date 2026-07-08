library(testthat)

# ---------------------------------------------------------------------------
# Synthetic data: a decreasing-mean baseline so the canonical (training-mean)
# order is well defined, plus a planted within-sample *ordinal-texture* signal.
# Over a contiguous block of features, cases get an alternating-parity boost
# (odd-index features up) and controls get the opposite boost (even-index
# features up). Within the block the local up/down comparisons between adjacent
# canonical neighbours therefore flip between the two classes, producing opposite
# ordinal-LBP codes and well-separated code histograms. The signal is purely a
# rank reshaping (scale-free), exactly what inv-olbp is designed to capture, and
# it is balanced across classes so the canonical order stays the feature index
# order (the zig-zag is preserved at canonical-neighbour granularity).
# ---------------------------------------------------------------------------
.make_inv_olbp_data <- function(n = 120L, p = 40L, seed = 71L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  base_mean <- seq(from = 120, to = 12, length.out = p)
  X <- matrix(0, nrow = n, ncol = p,
              dimnames = list(paste0("S", seq_len(n)), features))
  for (j in seq_len(p)) {
    # gamma with mean base_mean[j] (shape 20) -> strictly positive, mild noise
    X[, j] <- stats::rgamma(n, shape = 20, rate = 20 / base_mean[j])
  }
  block <- 5:24
  odd <- block[seq(1L, length(block), by = 2L)]
  even <- block[seq(2L, length(block), by = 2L)]
  X[y == 1L, odd] <- X[y == 1L, odd] + 80    # cases: odd-parity features up
  X[y == 0L, even] <- X[y == 0L, even] + 80  # controls: even-parity features up
  list(X = X, y = y, block = features[block], features = features)
}

.auc_or_wilcoxon <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# In-package, singlesample_score_call() resolves the manifest score_fn
# (score_inv_olbp) directly from the OmicSelector namespace -- no adapter
# registration is needed; the canonical-dispatch test below exercises that path.

test_that("inv-olbp fit/score roundtrip has the right shape and class", {
  dat <- .make_inv_olbp_data()
  model <- fit_inv_olbp(dat$X, dat$y)
  score <- score_inv_olbp(model, dat$X)

  expect_s3_class(model, "inv_olbp_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_setequal(model$order_perm, colnames(dat$X))
  expect_identical(model$n_codes, 2L^4L)
  expect_length(model$neighbor_offsets, 4L)
  expect_length(model$head$center, model$n_codes)
  expect_length(model$head$coef, model$n_codes)

  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})


test_that("inv-olbp scores are larger for case-like specimens", {
  dat <- .make_inv_olbp_data(seed = 72L)
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_inv_olbp(dat$X[train, ], dat$y[train])
  score <- score_inv_olbp(model, dat$X[test, ])

  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon(dat$y[test], score), 0.7)
})


test_that("inv-olbp passes the canonical row-equivariance gate", {
  dat <- .make_inv_olbp_data(seed = 73L)
  model <- fit_inv_olbp(dat$X[1:90, ], dat$y[1:90])
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_inv_olbp(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})


test_that("inv-olbp scoring is exactly invariant to per-sample positive scaling", {
  dat <- .make_inv_olbp_data(seed = 74L)
  model <- fit_inv_olbp(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_inv_olbp(model, row * 7),
               score_inv_olbp(model, row), tolerance = 1e-8)
  expect_equal(score_inv_olbp(model, row * 1e6),
               score_inv_olbp(model, row), tolerance = 1e-8)

  # Scaling some rows of a batch leaves the other rows' scores untouched.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e5
  expect_equal(score_inv_olbp(model, scaled),
               score_inv_olbp(model, batch), tolerance = 1e-8)
})


test_that("inv-olbp handles partial feature overlap (subset + reorder) consistently", {
  dat <- .make_inv_olbp_data(seed = 75L)
  model <- fit_inv_olbp(dat$X, dat$y)

  keep <- sample(colnames(dat$X), 25L)  # subset
  keep <- rev(keep)                     # and reorder the columns
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_inv_olbp(model, X_partial)

  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # A single row scored alone equals its score within the batch (row-independent,
  # and unaffected by which columns / their order are presented).
  expect_equal(score_inv_olbp(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Column order does not matter: re-permuting present columns gives equal scores.
  X_reperm <- X_partial[, sample(ncol(X_partial)), drop = FALSE]
  expect_equal(score_inv_olbp(model, X_reperm), s_partial, tolerance = 1e-12)
})


test_that("inv-olbp is deployable through the canonical dispatch", {
  skip_if_not(exists("singlesample_score_call"))
  dat <- .make_inv_olbp_data(seed = 76L)
  model <- fit_inv_olbp(dat$X[1:90, ], dat$y[1:90])
  X_test <- dat$X[91:120, ]

  expect_equal(
    singlesample_score_call("inv-olbp", model, X_test),
    score_inv_olbp(model, X_test),
    tolerance = 1e-12
  )
})


test_that("inv-olbp returns neutral 0 below the feature-overlap floor", {
  dat <- .make_inv_olbp_data(seed = 77L)
  model <- fit_inv_olbp(dat$X, dat$y)  # default min_features = max(P+1, 4) = 5

  # Fewer than min_features present -> neutral 0 for every row.
  X_few <- dat$X[1:6, model$order_perm[1:3], drop = FALSE]
  expect_equal(score_inv_olbp(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", seq_len(4L)),
                                   c("other-a", "other-b")))
  expect_equal(score_inv_olbp(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)

  # Exactly P (= 4) present, still <= P -> neutral 0 (circular code undefined).
  model_lowmin <- fit_inv_olbp(dat$X, dat$y, hp = list(min_features = 4L))
  X_eqP <- dat$X[1:5, model_lowmin$order_perm[1:4], drop = FALSE]
  expect_equal(score_inv_olbp(model_lowmin, X_eqP), rep(0, nrow(X_eqP)),
               tolerance = 1e-12)
})


test_that("inv-olbp fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_inv_olbp_data(seed = 78L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_inv_olbp(dat$X, dat$y)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_inv_olbp(dat$X, dat$y)
  expect_identical(model_1, model_2)

  s1 <- score_inv_olbp(model_1, dat$X[1:20, ])
  s2 <- score_inv_olbp(model_2, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
})


test_that("inv-olbp hp validation errors are explicit", {
  dat <- .make_inv_olbp_data(seed = 79L)

  expect_error(fit_inv_olbp(dat$X, dat$y, hp = list(foo = 1)), "unknown hp")
  expect_error(fit_inv_olbp(dat$X, dat$y, hp = list(neighbors = 3L)),
               "neighbors")
  expect_error(fit_inv_olbp(dat$X, dat$y, hp = list(neighbors = 10L)),
               "neighbors")
  expect_error(fit_inv_olbp(dat$X, dat$y, hp = list(neighbors = 4.5)),
               "neighbors")
  # Out-of-integer-range must hit the explicit validation error, not as.integer NA.
  expect_error(fit_inv_olbp(dat$X, dat$y, hp = list(neighbors = 3e9)),
               "neighbors")
  expect_error(fit_inv_olbp(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_inv_olbp(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  # Too few / too many features for the chosen P are explicit fit errors.
  X_small <- dat$X[, 1:3, drop = FALSE]
  expect_error(fit_inv_olbp(X_small, dat$y, hp = list(neighbors = 4L)),
               "neighbors|min_features")
})


test_that("inv-olbp scores a single row identically standalone and in batch", {
  dat <- .make_inv_olbp_data(seed = 80L)
  model <- fit_inv_olbp(dat$X, dat$y)
  s_all <- score_inv_olbp(model, dat$X[1:8, ])
  for (i in c(1L, 4L, 8L)) {
    s_one <- score_inv_olbp(model, dat$X[i, , drop = FALSE])
    expect_length(s_one, 1L)
    expect_equal(s_one, s_all[i], tolerance = 1e-12)
  }
})


test_that("inv-olbp OLBP codes and histogram match a hand-computed D=6, P=4 case", {
  # Auditable bit-packing check on a tiny example. Values in canonical order:
  #   a = (10, 60, 20, 50, 30, 40)  ->  ranks (1, 6, 2, 5, 3, 4)
  # Offsets c(-2,-1,1,2) carry weights (1, 2, 4, 8). Hand-derivation gives codes
  #   (0, 15, 1, 14, 9, 6), all distinct, so each of bins {0,1,6,9,14,15} holds 1/6.
  offsets <- .inv_olbp_offsets(4L)
  expect_identical(offsets, as.integer(c(-2, -1, 1, 2)))

  a <- c(10, 60, 20, 50, 30, 40)
  codes <- .inv_olbp_codes(a, offsets)
  expect_identical(codes, as.integer(c(0, 15, 1, 14, 9, 6)))

  h <- .inv_olbp_hist(a, offsets, 16L)
  expect_length(h, 16L)
  expect_equal(sum(h), 1, tolerance = 1e-12)
  expected <- numeric(16L)
  expected[c(0L, 1L, 6L, 9L, 14L, 15L) + 1L] <- 1 / 6
  expect_equal(h, expected, tolerance = 1e-12)

  # Strictly monotone rescaling of the values leaves codes (hence histogram)
  # identical -> the descriptor is purely ordinal.
  expect_identical(.inv_olbp_codes(a * 1000, offsets), codes)
  expect_identical(.inv_olbp_codes(sqrt(a), offsets), codes)
})


test_that("inv-olbp falls back to a stable diagonal-LDA head under GLM separation", {
  # Helper-level regression (load-bearing, non-vacuous): a perfectly
  # class-separating standardised design must NOT be frozen as a divergent
  # converged GLM fit -- it must signal the diagonal-LDA fallback (return NULL).
  try_glm <- OmicSelector:::.inv_olbp_try_glm
  Za <- matrix(c(-2, -2, -1, -1, 1, 1, 2, 2), ncol = 1)
  ysep <- c(0, 0, 0, 0, 1, 1, 1, 1)
  expect_null(try_glm(Za, ysep))
  # A well-separated-but-not-degenerate design still fits a GLM.
  set.seed(1)
  Zok <- matrix(stats::rnorm(40), ncol = 2)
  yok <- rep(c(0, 1), each = 10)
  expect_type(try_glm(Zok, yok), "list")

  # Integration: a strongly class-structured fit yields a valid head type with
  # finite, correctly oriented scores (whichever branch is taken).
  set.seed(303)
  p <- 24L; n <- 80L
  X <- matrix(stats::rgamma(n * p, shape = 8, rate = 1), n, p,
              dimnames = list(paste0("S", seq_len(n)), paste0("f", seq_len(p))))
  y <- rep(c(0, 1), each = n / 2L)
  X[y == 1L, 1:6] <- X[y == 1L, 1:6] + 300
  X[y == 0L, 7:12] <- X[y == 0L, 7:12] + 300
  model <- fit_inv_olbp(X, y)
  expect_true(model$head$type %in% c("glm", "diag_lda"))
  sc <- score_inv_olbp(model, X)
  expect_true(all(is.finite(sc)))
  expect_gt(mean(sc[y == 1L]) - mean(sc[y == 0L]), 0)
})


test_that("inv-olbp eps is an sd FLOOR (active bins kept), not an activity threshold", {
  dat <- .make_inv_olbp_data(seed = 311L)
  # A large eps must only FLOOR the per-bin sd, NOT deactivate positive-variance
  # bins. Under the old activity-threshold code (active <- raw_sd > eps) every
  # histogram bin (sd << 0.5) would be deactivated -> a constant non-
  # discriminating head with score sd 0. This test fails that old construction.
  model_floor <- fit_inv_olbp(dat$X, dat$y, hp = list(neighbors = 4L, eps = 0.5))
  expect_true(all(model_floor$head$scale >= 0.5 - 1e-12))      # every sd floored
  expect_true(any(model_floor$head$active))                    # positive-var kept
  expect_gt(stats::sd(score_inv_olbp(model_floor, dat$X)), 0)  # still discriminates

  # Default eps: only genuinely zero-variance bins are inactive.
  model_def <- fit_inv_olbp(dat$X, dat$y)
  expect_true(all(model_def$head$scale >= model_def$hp$eps - 1e-15))
  expect_true(any(model_def$head$active))
})
