library(testthat)

# ---------------------------------------------------------------------------
# Synthetic data with a planted QF-SHAPE difference in the sorted rCLR profile
# (the substrate inv-fdaqf reads). Both classes share the same overall rCLR
# range; only the SHAPE of the distribution differs:
#   - Controls: per-feature log-means follow a smooth LINEAR ramp -> the sorted
#     rCLR profile (quantile function) is near-linear.
#   - Cases: per-feature log-means follow a convex (cubic) profile -> most
#     features sit low and a few rise steeply, giving a HEAVIER UPPER TAIL in the
#     sorted rCLR quantile function.
# The signal is a pure distributional-shape difference; rCLR removes the
# per-sample level/scale, so the FPCA-then-LDA pipeline must separate the classes
# on QF SHAPE alone. The class log-mean profiles are matched in range (2..5).
# ---------------------------------------------------------------------------
.make_inv_fdaqf_data <- function(n = 200L, p = 60L, seed = 71L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(0, nrow = n, ncol = p,
              dimnames = list(paste0("S", seq_len(n)), features))
  rank_frac <- (seq_len(p) - 1) / (p - 1)
  ctrl_mu <- 2 + 3 * rank_frac                              # near-linear QF
  case_mu <- 2 + 3 * rank_frac^3                            # heavier upper tail
  for (i in seq_len(n)) {
    mu <- if (y[i] == 1L) case_mu else ctrl_mu
    X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
  }
  list(X = X, y = y)
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

# Rebuild the training QF matrix exactly as fit_inv_fdaqf does, for the
# non-circular FPCA-correctness checks.
.fdaqf_train_G <- function(X, L, min_features) {
  qfs <- lapply(seq_len(nrow(X)), function(i) .inv_fdaqf_qf(X[i, ], L, min_features))
  keep <- !vapply(qfs, is.null, logical(1))
  do.call(rbind, qfs[keep])
}


test_that("inv-fdaqf fit/score roundtrip has the right shape and class", {
  dat <- .make_inv_fdaqf_data()
  model <- fit_inv_fdaqf(dat$X, dat$y)
  score <- score_inv_fdaqf(model, dat$X)

  expect_s3_class(model, "inv_fdaqf_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_identical(model$descriptor_dim, model$fpca$K)
  expect_identical(model$fpca$K, 6L)                       # full rank here
  expect_equal(dim(model$fpca$V), c(model$hp$resample_len, model$fpca$K))
  expect_length(model$fpca$qbar, model$hp$resample_len)
  expect_length(model$head$center, model$fpca$K)
  expect_length(model$head$scale, model$fpca$K)
  expect_length(model$head$w, model$fpca$K)

  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_null(names(score))                                # plain, unnamed vector
  expect_true(all(is.finite(score)))
})


test_that("inv-fdaqf separates a planted QF-shape difference (heavier upper tail)", {
  dat <- .make_inv_fdaqf_data(seed = 72L)
  tr <- c(1:70, 101:170)
  te <- c(71:100, 171:200)
  model <- fit_inv_fdaqf(dat$X[tr, ], dat$y[tr])
  score <- score_inv_fdaqf(model, dat$X[te, ])

  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[te] == 1]), mean(score[dat$y[te] == 0]))
  expect_gt(.auc_or_wilcoxon(dat$y[te], score), 0.7)
})


test_that("inv-fdaqf passes the canonical row-equivariance gate", {
  dat <- .make_inv_fdaqf_data(seed = 73L)
  tr <- c(1:70, 101:170)
  te <- c(71:100, 171:200)
  model <- fit_inv_fdaqf(dat$X[tr, ], dat$y[tr])
  X_test <- dat$X[te, ]
  score_fun <- function(model, X, meta) score_inv_fdaqf(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})


test_that("inv-fdaqf scoring is exactly invariant to per-sample positive scaling", {
  dat <- .make_inv_fdaqf_data(seed = 74L)
  model <- fit_inv_fdaqf(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_inv_fdaqf(model, row * 7),
               score_inv_fdaqf(model, row), tolerance = 1e-10)
  expect_equal(score_inv_fdaqf(model, row * 1e6),
               score_inv_fdaqf(model, row), tolerance = 1e-10)

  # Scaling some rows of a batch leaves the other rows' scores untouched.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e5
  expect_equal(score_inv_fdaqf(model, scaled),
               score_inv_fdaqf(model, batch), tolerance = 1e-10)

  # Mixed: each row scaled by its own positive factor (column-major recycling of
  # a length-nrow vector scales row i by factor[i]).
  mixed <- batch * c(2, 0.5, 10, 1, 1e3, 7)
  expect_equal(score_inv_fdaqf(model, mixed),
               score_inv_fdaqf(model, batch), tolerance = 1e-10)
})


test_that("inv-fdaqf handles partial feature overlap (subset + reorder) consistently", {
  dat <- .make_inv_fdaqf_data(seed = 75L)
  model <- fit_inv_fdaqf(dat$X, dat$y)

  keep <- rev(sample(colnames(dat$X), 30L))     # subset and reorder columns
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_inv_fdaqf(model, X_partial)

  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # A single row scored alone equals its score within the batch.
  expect_equal(score_inv_fdaqf(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # The QF sorts present rCLR values, so re-permuting present columns cannot
  # change the score.
  X_reperm <- X_partial[, sample(ncol(X_partial)), drop = FALSE]
  expect_equal(score_inv_fdaqf(model, X_reperm), s_partial, tolerance = 1e-12)
})


test_that("inv-fdaqf is deployable through the canonical dispatch", {
  skip_if_not(exists("singlesample_score_call"))
  dat <- .make_inv_fdaqf_data(seed = 76L)
  tr <- c(1:70, 101:170)
  te <- c(71:100, 171:200)
  model <- fit_inv_fdaqf(dat$X[tr, ], dat$y[tr])
  X_test <- dat$X[te, ]
  direct <- score_inv_fdaqf(model, X_test)

  # Post-integration the manifest carries the inv-fdaqf row and score_inv_fdaqf
  # lives in the OmicSelector namespace, so the canonical adapter resolves it
  # directly. In the pre-integration staged self-test score_inv_fdaqf is not yet
  # in the namespace, so we clone a roster row if absent and register a thin
  # canonical adapter to exercise the same dispatch wiring. Either way the
  # dispatched score must equal the direct call.
  roster <- singlesample_method_roster()
  if (!"inv-fdaqf" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "inv-fdaqf"
    tmpl$fit_fn <- "fit_inv_fdaqf"
    tmpl$score_fn <- "score_inv_fdaqf"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_inv_fdaqf", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "inv-fdaqf",
      function(model, X, meta) score_inv_fdaqf(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("inv-fdaqf", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})


test_that("inv-fdaqf returns neutral 0 below the feature-overlap floor", {
  dat <- .make_inv_fdaqf_data(seed = 77L)
  model <- fit_inv_fdaqf(dat$X, dat$y)            # default min_features = 8

  # Fewer than min_features model-universe features present -> neutral 0.
  X_few <- dat$X[1:6, colnames(dat$X)[1:5], drop = FALSE]
  expect_equal(score_inv_fdaqf(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", seq_len(4L)),
                                   c("other-a", "other-b", "other-c")))
  expect_equal(score_inv_fdaqf(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)

  # Enough universe features present, but ONE specimen has too few nonzero
  # coordinates -> that row alone scores the neutral 0; the others score normally.
  Xs <- dat$X[1:5, , drop = FALSE]
  Xs[2, ] <- 0
  Xs[2, 1:4] <- dat$X[2, 1:4]                     # only 4 nonzero < min_features 8
  s <- score_inv_fdaqf(model, Xs)
  expect_equal(s[2], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s)))
  expect_true(any(s[-2] != 0))
})


test_that("inv-fdaqf fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_inv_fdaqf_data(seed = 78L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_inv_fdaqf(dat$X, dat$y)
  # The method consumes no randomness, so the global RNG is untouched by a fit.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_inv_fdaqf(dat$X, dat$y)
  expect_identical(model_1, model_2)

  s1 <- score_inv_fdaqf(model_1, dat$X[1:20, ])
  s2 <- score_inv_fdaqf(model_2, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
})


test_that("inv-fdaqf hp validation errors are explicit", {
  dat <- .make_inv_fdaqf_data(seed = 79L)

  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(foo = 1)), "unknown hp")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(n_components = 6.5)),
               "n_components")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(n_components = 0L)),
               "n_components")
  # Out-of-integer-range must hit the explicit error, not an as.integer overflow.
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(n_components = 3e9)),
               "n_components")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(resample_len = 4L)),
               "resample_len")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(resample_len = 3e9)),
               "resample_len")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(shrink = -0.1)), "shrink")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(min_features = 2L)),
               "min_features")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(seed = -1)), "hp\\$seed")
  # `$` partial matching would let "min" bind to "min_features"; the exact [[ ]]
  # reads + allowed-list reject the unknown field instead.
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(min = 3L)), "unknown hp")
  # Fully-unnamed hp (names(hp) == NULL) must be rejected, not silently ignored:
  # setdiff(NULL, allowed) is empty, so the allowed-list check alone would pass it.
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(64L)), "must be named")
  # A partially-unnamed entry is rejected too.
  expect_error(fit_inv_fdaqf(dat$X, dat$y, hp = list(n_components = 6L, 99)),
               "must be named")
  # Duplicated field names must be rejected, not silently last-wins.
  expect_error(
    fit_inv_fdaqf(dat$X, dat$y, hp = list(n_components = 4L, n_components = 6L)),
    "duplicated hp")
})


test_that("inv-fdaqf scores a single row identically standalone and in batch", {
  dat <- .make_inv_fdaqf_data(seed = 80L)
  model <- fit_inv_fdaqf(dat$X, dat$y)
  s_all <- score_inv_fdaqf(model, dat$X[1:8, ])
  for (i in c(1L, 4L, 8L)) {
    s_one <- score_inv_fdaqf(model, dat$X[i, , drop = FALSE])
    expect_length(s_one, 1L)
    expect_equal(s_one, s_all[i], tolerance = 1e-12)
  }
})


test_that("inv-fdaqf rCLR keeps present coordinates at the geometric mean", {
  # Structural-zero exclusion must key on v > 0, NOT z != 0: a present coordinate
  # whose centred-log equals the geometric mean (z == 0) must be kept, and a
  # constant positive profile must yield a length-(#positive) point mass at 0 so
  # the resampled QF is well-defined, not an empty vector.
  # v = (1, 2, 4): geomean = 2, so the middle present coordinate has rCLR exactly 0.
  zp <- .inv_fdaqf_rclr_pos(c(1, 2, 4))
  expect_length(zp, 3L)
  expect_equal(sort(zp), c(-log(2), 0, log(2)), tolerance = 1e-12)
  # Structural zeros (v == 0) are excluded; the present middle-0 is kept.
  expect_length(.inv_fdaqf_rclr_pos(c(0, 1, 2, 4, 0)), 3L)
  # A constant positive profile -> point mass at 0 of size = #positive features.
  expect_length(.inv_fdaqf_rclr_pos(rep(2, 7)), 7L)
  expect_true(all(.inv_fdaqf_rclr_pos(rep(2, 7)) == 0))
  # End to end: a specimen constant across all present features scores finite
  # (not dropped to neutral 0 because z != 0 emptied it).
  dat <- .make_inv_fdaqf_data(seed = 81L)
  model <- fit_inv_fdaqf(dat$X, dat$y)
  Xc <- dat$X[1:3, , drop = FALSE]
  Xc[1, ] <- 5
  expect_true(all(is.finite(score_inv_fdaqf(model, Xc))))
})


test_that("inv-fdaqf FPCA basis matches an external SVD ground truth", {
  dat <- .make_inv_fdaqf_data(seed = 82L)
  model <- fit_inv_fdaqf(dat$X, dat$y)
  V <- model$fpca$V
  qbar <- model$fpca$qbar
  K <- model$fpca$K
  L <- model$hp$resample_len

  # (1) Frozen loadings are orthonormal: V^T V = I_K.
  expect_equal(crossprod(V), diag(K), tolerance = 1e-10)
  expect_equal(dim(V), c(L, K))

  # (2) Independent SVD ground truth of the centred training QF matrix.
  G <- .fdaqf_train_G(dat$X, L, model$hp$min_features)
  Gc <- sweep(G, 2L, qbar, "-")

  # (2a) The stored basis V IS the top-K right singular subspace of Gc, verified
  # against a FRESH, independent svd() recompute (NOT the stored V). The rank-K
  # orthogonal projector V V^T is basis-orientation-invariant, so this is robust
  # to SVD sign flips and any within-degenerate-block rotation.
  V_indep <- svd(Gc)$v[, seq_len(K), drop = FALSE]
  expect_equal(crossprod(V_indep), diag(K), tolerance = 1e-10)   # also orthonormal
  expect_equal(tcrossprod(V), tcrossprod(V_indep), tolerance = 1e-8)

  svd_scores <- Gc %*% V                          # external rank-K projection

  # Projecting each training QF through the scorer's own projection helper must
  # recover the SVD scores row-for-row (non-circular: per-specimen crossprod path
  # vs the full matrix product).
  proj <- do.call(rbind, lapply(seq_len(nrow(G)), function(i)
    .inv_fdaqf_project(G[i, ], qbar, V)))
  expect_equal(proj, svd_scores, tolerance = 1e-10, ignore_attr = TRUE)

  # (3) Reconstruction qbar + V V^T (q - qbar) IS the rank-K truncation.
  recon <- t(apply(G, 1L, function(q)
    qbar + as.numeric(V %*% crossprod(V, q - qbar))))
  trunc <- sweep(svd_scores %*% t(V), 2L, qbar, "+")
  expect_equal(recon, trunc, tolerance = 1e-10, ignore_attr = TRUE)

  # (4) Components are ordered by DECREASING training score variance.
  comp_var <- apply(svd_scores, 2L, stats::var)
  expect_false(is.unsorted(rev(comp_var)))        # non-increasing across k
})


test_that("inv-fdaqf degrades gracefully on rank-deficient training QFs", {
  p <- 40L
  n <- 12L
  feats <- paste0("f", seq_len(p))

  # (a) All training QFs identical -> centred QF matrix is exactly rank 0 ->
  #     realised K = 0, every specimen scores the constant neutral 0.
  base <- exp(seq(1, 4, length.out = p))
  X_eq <- matrix(rep(base, each = n), nrow = n,
                 dimnames = list(paste0("R", seq_len(n)), feats))
  y <- rep(c(0, 1), each = n / 2L)
  m_eq <- fit_inv_fdaqf(X_eq, y)
  expect_identical(m_eq$fpca$rank, 0L)
  expect_identical(m_eq$fpca$K, 0L)
  expect_true(m_eq$fpca$K <= m_eq$fpca$rank + 0L) # realised K <= rank
  s_eq <- score_inv_fdaqf(m_eq, X_eq)
  expect_true(all(is.finite(s_eq)))
  expect_equal(s_eq, rep(0, n), tolerance = 1e-12)

  # (b) Exactly two distinct QF shapes -> centred QF matrix is rank 1 ->
  #     realised K = min(n_components, 1) = 1, finite valid scores.
  shapeA <- exp(seq(1, 4, length.out = p))
  shapeB <- exp(c(seq(1, 2.5, length.out = p / 2L),
                  seq(4, 2.6, length.out = p / 2L)))
  X2 <- matrix(0, nrow = n, ncol = p,
               dimnames = list(paste0("R", seq_len(n)), feats))
  for (i in seq_len(n)) X2[i, ] <- if (i <= n / 2L) shapeA else shapeB
  m2 <- fit_inv_fdaqf(X2, y)
  expect_identical(m2$fpca$rank, 1L)
  expect_true(m2$fpca$K <= m2$fpca$rank)
  expect_identical(m2$fpca$K, 1L)
  s2 <- score_inv_fdaqf(m2, X2)
  expect_length(s2, n)
  expect_true(all(is.finite(s2)))
})


test_that("inv-fdaqf head solves the ridge-LDA normal equations", {
  set.seed(5)
  Phi <- matrix(stats::rnorm(40 * 6), 40, 6)
  y <- rep(c(0, 1), each = 20)
  Phi[y == 1, ] <- Phi[y == 1, ] + 0.8         # genuine class separation
  head <- .inv_fdaqf_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)

  # Independent recompute of the standardize + ridge-LDA head.
  center <- colMeans(Phi)
  raw_sd <- apply(Phi, 2L, stats::sd)
  active <- raw_sd > 0
  scale <- pmax(raw_sd, 1e-6)
  Z <- sweep(sweep(Phi, 2L, center, "-"), 2L, scale, "/")
  Z[, !active] <- 0
  mu1 <- colMeans(Z[y == 1, , drop = FALSE])
  mu0 <- colMeans(Z[y == 0, , drop = FALSE])
  df <- nrow(Z) - 2L
  Sw <- (crossprod(sweep(Z[y == 1, , drop = FALSE], 2L, mu1, "-")) +
         crossprod(sweep(Z[y == 0, , drop = FALSE], 2L, mu0, "-"))) / df
  Sw <- (Sw + t(Sw)) / 2
  Swr <- Sw
  diag(Swr) <- diag(Swr) + head$ridge
  w_check <- as.numeric(chol2inv(chol(Swr)) %*% (mu1 - mu0))
  b_check <- -0.5 * sum((mu1 + mu0) * w_check)

  expect_equal(head$w, w_check, tolerance = 1e-8)
  expect_equal(head$b, b_check, tolerance = 1e-8)

  # And the scorer's linear predictor equals the hand-computed LDA value.
  dat <- .make_inv_fdaqf_data(seed = 83L)
  model <- fit_inv_fdaqf(dat$X, dat$y)
  q <- .inv_fdaqf_qf(dat$X[7, ], model$hp$resample_len, model$hp$min_features)
  s <- .inv_fdaqf_project(q, model$fpca$qbar, model$fpca$V)
  z <- (s - model$head$center) / model$head$scale
  z[!model$head$active] <- 0
  manual <- model$head$b + sum(model$head$w * z)
  expect_equal(score_inv_fdaqf(model, dat$X[7, , drop = FALSE]), manual,
               tolerance = 1e-9)
})


test_that("inv-fdaqf head stays positive definite after ridge on collinear scores", {
  set.seed(9)
  # 8 samples, 6 components, tiny variance + an exactly collinear pair: the
  # pooled within-class covariance is (near-)singular, so the ridge-to-PD loop
  # must still produce a finite Cholesky-based head.
  Phi <- matrix(stats::rnorm(8 * 6, sd = 1e-7), 8, 6)
  Phi[, 1] <- Phi[, 2]
  y <- rep(c(0, 1), each = 4)
  head <- .inv_fdaqf_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)
  expect_true(all(is.finite(head$w)))
  expect_true(is.finite(head$b))
  expect_true(head$ridge > 0)

  # Integration: a fit on a p > n regime (few rows after the min_features filter)
  # also stays finite end to end.
  dat <- .make_inv_fdaqf_data(n = 12L, p = 40L, seed = 84L)
  model <- fit_inv_fdaqf(dat$X, dat$y)
  sc <- score_inv_fdaqf(model, dat$X)
  expect_length(sc, nrow(dat$X))
  expect_true(all(is.finite(sc)))
})
