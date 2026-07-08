# Tests for the ss-struct-ilr structural-ILR single-sample scorer (Stage 1 of
# idea_report_1570797495812882433). Synthetic annotation + compositional data only;
# CPU-only, fast. Verifies: frozen fit/score round-trip, row-equivariance (single==batch),
# no test-batch leakage, exact per-specimen scale invariance, and neutral degenerate paths.

test_that("ss-struct-ilr fits, scores, and is row-equivariant (single == batch)", {
  set.seed(42)
  n <- 40L; p <- 30L
  features <- sprintf("hsa-miR-%03d-3p", seq_len(p))
  # 4 genomic clusters, GC fraction spread across [0.3, 0.7]
  annotation <- data.frame(
    feature = features,
    cluster = paste0("clu", rep(1:4, length.out = p)),
    gc = stats::runif(p, 0.30, 0.70),
    stringsAsFactors = FALSE
  )
  # Non-negative abundances; inject a mild class signal into a few features.
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(NULL, features))
  y <- rep(c(0L, 1L), each = n / 2L)
  X[y == 1L, 1:5] <- X[y == 1L, 1:5] * 1.8

  model <- fit_ss_struct_ilr(X, y, annotation = annotation)
  expect_s3_class(model, "ss_struct_ilr_model")
  expect_true(length(model$sbp) >= 1L)
  # every artefact serialised for frozen scoring
  expect_true(all(c("sbp", "center", "scale", "weights", "feature_universe",
                    "balance_names", "pseudocount", "aggregate") %in% names(model)))

  s_full <- score_ss_struct_ilr(model, X)
  expect_length(s_full, n)
  expect_true(all(is.finite(s_full)))

  # No test-batch leakage: scoring one row alone == its position in the full batch.
  for (i in c(1L, 7L, n)) {
    s_one <- score_ss_struct_ilr(model, X[i, , drop = FALSE])
    expect_equal(s_one, s_full[i], tolerance = 1e-9)
  }

  # The canonical single-sample gate (singletons, permutations, varied/duplicated/
  # adversarial batches, no model-state mutation).
  expect_true(singlesample_is_row_equivariant(score_ss_struct_ilr, model, X))
  expect_silent(singlesample_assert_row_equivariant(score_ss_struct_ilr, model, X))
})

test_that("ss-struct-ilr is exactly scale-invariant per specimen", {
  set.seed(7)
  n <- 24L; p <- 24L
  features <- sprintf("hsa-miR-%03d-5p", seq_len(p))
  annotation <- data.frame(
    feature = features,
    cluster = paste0("clu", rep(1:3, length.out = p)),
    gc = stats::runif(p, 0.30, 0.70),
    stringsAsFactors = FALSE
  )
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(NULL, features))
  y <- rep(c(0L, 1L), each = n / 2L)
  model <- fit_ss_struct_ilr(X, y, annotation = annotation)

  base <- score_ss_struct_ilr(model, X)
  Xs <- X * rep(c(1e-3, 1, 7, 1e4), length.out = n)  # per-specimen rescale
  scaled <- score_ss_struct_ilr(model, Xs)
  expect_equal(scaled, base, tolerance = 1e-6)
})

test_that("ss-struct-ilr scores neutrally when no balances are constructible", {
  set.seed(1)
  n <- 8L; p <- 6L
  features <- sprintf("hsa-miR-%03d-3p", seq_len(p))
  # All features in ONE cluster with identical GC -> no cluster-coarse balance and
  # no GC tertile split -> degenerate SBP -> neutral 0 scorer.
  annotation <- data.frame(
    feature = features, cluster = "clu1", gc = 0.5, stringsAsFactors = FALSE
  )
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(NULL, features))
  y <- rep(c(0L, 1L), each = n / 2L)
  model <- fit_ss_struct_ilr(X, y, annotation = annotation)
  expect_equal(score_ss_struct_ilr(model, X), rep(0, n))
})

test_that("ss-struct-ilr GC SBP is a faithful tertile partition (no tertile dropped)", {
  # One cluster, 9 features with GC spanning all three tertiles -> the within-cluster
  # GC SBP must represent ALL three tertiles via two balances (high vs mid+low; mid vs
  # low), never dropping the middle tertile.
  features <- sprintf("hsa-miR-%03d-3p", seq_len(9L))
  annotation <- data.frame(
    feature = features, cluster = "cl1",
    gc = c(0.20, 0.25, 0.30,  0.45, 0.50, 0.55,  0.75, 0.80, 0.85),
    stringsAsFactors = FALSE
  )
  n <- 12L
  X <- matrix(stats::rgamma(n * 9L, shape = 2, rate = 1), nrow = n,
              dimnames = list(NULL, features))
  y <- rep(c(0L, 1L), each = n / 2L)
  model <- fit_ss_struct_ilr(X, y, annotation = annotation)

  expect_true("gc_hi__cl1" %in% names(model$sbp))   # high tertile vs (mid + low)
  expect_true("gc_lo__cl1" %in% names(model$sbp))   # mid tertile vs low tertile
  hi <- model$sbp[["gc_hi__cl1"]]
  lo <- model$sbp[["gc_lo__cl1"]]
  # All three tertiles are represented; the union of both balances covers every feature
  # (the middle tertile is NOT dropped).
  covered <- unique(c(hi$numerator, hi$denominator, lo$numerator, lo$denominator))
  expect_setequal(covered, features)
  expect_setequal(hi$numerator, features[7:9])             # high tertile isolated
  expect_setequal(lo$numerator, features[4:6])             # middle tertile present
  expect_setequal(lo$denominator, features[1:3])           # low tertile
})

test_that("ss-struct-ilr floors an empty-support specimen to neutral 0", {
  set.seed(5)
  n <- 16L; p <- 18L
  features <- sprintf("hsa-miR-%03d-3p", seq_len(p))
  annotation <- data.frame(
    feature = features, cluster = paste0("clu", rep(1:3, length.out = p)),
    gc = stats::runif(p, 0.30, 0.70), stringsAsFactors = FALSE
  )
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(NULL, features))
  y <- rep(c(0L, 1L), each = n / 2L)
  model <- fit_ss_struct_ilr(X, y, annotation = annotation)

  # An all-zero specimen has empty positive support -> EXACTLY neutral 0 (the bug fix:
  # without the floor it scored a spurious non-zero off-origin value).
  z_row <- matrix(0, nrow = 1L, ncol = p, dimnames = list(NULL, features))
  expect_identical(score_ss_struct_ilr(model, z_row), 0)
  # A FLAT all-equal-positive specimen is a VALID composition and is scored normally
  # (not floored); its score need not be 0.
  flat <- matrix(3.5, nrow = 1L, ncol = p, dimnames = list(NULL, features))
  expect_true(is.finite(score_ss_struct_ilr(model, flat)))
})
