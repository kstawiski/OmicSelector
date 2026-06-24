# Tests for the DA-cVAE single-sample scorer (Stage 2/3 of idea_report_1570797495812882433).
# Synthetic annotation + compositional data + synthetic platform factor; a TINY CPU torch
# train (no GPU, no real data). Verifies: frozen fit/score round-trip, row-equivariance
# (single == batch) and exact per-specimen scale invariance for BOTH the disease score and
# the novelty score, no test-batch leakage, the no-GRL ablation + single-platform paths,
# cross-process determinism (seed-matched digests), and the neutral degenerate path.
# Skips cleanly when the reticulate torch venv is unavailable.

skip_if_no_torch <- function() {
  testthat::skip_if_not_installed("reticulate")
  if (!reticulate::py_module_available("torch")) {
    testthat::skip("python torch not available in the active reticulate env")
  }
}

# Synthetic compositional cohort with a genomic-cluster/GC annotation and a platform
# factor. A mild class signal is injected into a handful of features.
mk_dacvae_data <- function(n = 48L, p = 24L, n_platforms = 3L, seed = 11L) {
  set.seed(seed)
  features <- sprintf("hsa-miR-%03d-3p", seq_len(p))
  annotation <- data.frame(
    feature = features,
    cluster = paste0("clu", rep(seq_len(4L), length.out = p)),
    gc = stats::runif(p, 0.30, 0.70),
    stringsAsFactors = FALSE
  )
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(NULL, features))
  y <- rep(c(0L, 1L), length.out = n)
  X[y == 1L, 1:5] <- X[y == 1L, 1:5] * 1.8
  meta <- data.frame(
    platform = paste0("P", rep(seq_len(n_platforms), length.out = n)),
    stringsAsFactors = FALSE
  )
  list(X = X, y = y, meta = meta, annotation = annotation, features = features)
}

test_that("dacvae fits, scores, and is row-equivariant (disease + novelty)", {
  skip_if_no_torch()
  d <- mk_dacvae_data()
  m <- fit_dacvae(d$X, d$y, meta_train = d$meta, annotation = d$annotation,
                  hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "dacvae_model")
  expect_true(m$n_platforms >= 2L)                         # adversary engaged
  expect_true(all(c("weights", "head_w", "head_b", "center", "scale", "sbp",
                    "novelty_ref", "platform_levels") %in% names(m)))

  s <- score_dacvae(m, d$X)
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
  nov <- score_dacvae_novelty(m, d$X)
  expect_length(nov, nrow(d$X))
  expect_true(all(is.finite(nov)) && all(nov >= 0))

  # No test-batch leakage: scoring one row alone == its position in the full batch.
  for (i in c(1L, 7L, nrow(d$X))) {
    expect_equal(score_dacvae(m, d$X[i, , drop = FALSE]), s[i], tolerance = 1e-9)
    expect_equal(score_dacvae_novelty(m, d$X[i, , drop = FALSE]), nov[i],
                 tolerance = 1e-9)
  }

  # The canonical single-sample gate for BOTH score functions.
  expect_true(singlesample_is_row_equivariant(
    score_dacvae, m, d$X, model_digest = .dacvae_model_digest))
  expect_silent(singlesample_assert_row_equivariant(
    score_dacvae, m, d$X, model_digest = .dacvae_model_digest))
  expect_true(singlesample_is_row_equivariant(
    score_dacvae_novelty, m, d$X, model_digest = .dacvae_model_digest))
  expect_silent(singlesample_assert_row_equivariant(
    score_dacvae_novelty, m, d$X, model_digest = .dacvae_model_digest))
})

test_that("dacvae is exactly scale-invariant per specimen (disease + novelty)", {
  skip_if_no_torch()
  d <- mk_dacvae_data(n = 36L, p = 21L, seed = 3L)
  m <- fit_dacvae(d$X, d$y, meta_train = d$meta, annotation = d$annotation,
                  hp = list(epochs = 50L, device = "cpu"))
  base <- score_dacvae(m, d$X)
  base_nov <- score_dacvae_novelty(m, d$X)
  Xs <- d$X * rep(c(1e-3, 1, 7, 1e4), length.out = nrow(d$X))   # per-specimen rescale
  expect_equal(score_dacvae(m, Xs), base, tolerance = 1e-6)
  expect_equal(score_dacvae_novelty(m, Xs), base_nov, tolerance = 1e-6)
})

test_that("dacvae no-GRL ablation and single-platform paths fit and score", {
  skip_if_no_torch()
  d <- mk_dacvae_data(n = 32L, seed = 9L)
  # grl_lambda = 0 -> adversary inert (the no-harm-guard ablation).
  m0 <- fit_dacvae(d$X, d$y, meta_train = d$meta, annotation = d$annotation,
                   hp = list(epochs = 40L, grl_lambda = 0.0))
  expect_true(all(is.finite(score_dacvae(m0, d$X))))
  # A single platform -> n_platforms == 1 (plain conditional VAE + label head).
  meta1 <- data.frame(platform = rep("only", nrow(d$X)), stringsAsFactors = FALSE)
  m1 <- fit_dacvae(d$X, d$y, meta_train = meta1, annotation = d$annotation,
                   hp = list(epochs = 40L))
  expect_identical(m1$n_platforms, 1L)
  expect_true(all(is.finite(score_dacvae(m1, d$X))))
  # NULL meta -> also single platform.
  mn <- fit_dacvae(d$X, d$y, meta_train = NULL, annotation = d$annotation,
                   hp = list(epochs = 30L))
  expect_identical(mn$n_platforms, 1L)
})

test_that("dacvae is cross-process reproducible (seed-matched digests match)", {
  skip_if_no_torch()
  d <- mk_dacvae_data(n = 32L, seed = 4L)
  m1 <- fit_dacvae(d$X, d$y, meta_train = d$meta, annotation = d$annotation,
                   hp = list(epochs = 40L, device = "cpu", seed = 42L))
  m2 <- fit_dacvae(d$X, d$y, meta_train = d$meta, annotation = d$annotation,
                   hp = list(epochs = 40L, device = "cpu", seed = 42L))
  expect_identical(.dacvae_model_digest(m1), .dacvae_model_digest(m2))
})

test_that("dacvae scores neutrally when no balances are constructible", {
  skip_if_no_torch()
  set.seed(1)
  n <- 8L; p <- 6L
  features <- sprintf("hsa-miR-%03d-3p", seq_len(p))
  # One cluster, identical GC -> no cluster-coarse balance, no GC tertile split ->
  # degenerate SBP -> neutral 0 scorer (fit does not even touch torch).
  annotation <- data.frame(feature = features, cluster = "clu1", gc = 0.5,
                           stringsAsFactors = FALSE)
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(NULL, features))
  y <- rep(c(0L, 1L), each = n / 2L)
  meta <- data.frame(platform = rep(c("P1", "P2"), each = n / 2L),
                     stringsAsFactors = FALSE)
  m <- fit_dacvae(X, y, meta_train = meta, annotation = annotation,
                  hp = list(epochs = 5L))
  expect_equal(score_dacvae(m, X), rep(0, n))
  expect_equal(score_dacvae_novelty(m, X), rep(0, n))
})
