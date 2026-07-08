library(testthat)

# ---------------------------------------------------------------------------
# Synthetic data: cases have a STEP FUNCTION in the canonical (decreasing-mean)
# order -- features 1:30 uniformly elevated → rCLR landscape is flat-high then
# flat-low → GLCM mass concentrated on (high,high) and (low,low) cells →
# LOW contrast, HIGH energy. Controls have NO extra boost → smooth decreasing
# gradient → GLCM mass spread across sub-diagonal → HIGHER contrast, LOWER energy.
# This contrast in texture descriptors is the planted signal for the LDA head.
# Strictly decreasing base means (shape=30) keep the canonical order ≈ index order.
# ---------------------------------------------------------------------------
.make_inv_glcm_data <- function(n = 200L, p = 60L, seed = 71L) {
  set.seed(seed)
  y <- rep(c(0L, 1L), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  base_mean <- seq(from = 80, to = 10, length.out = p)
  X <- matrix(0, nrow = n, ncol = p,
              dimnames = list(paste0("S", seq_len(n)), features))
  for (j in seq_len(p)) {
    X[, j] <- stats::rgamma(n, shape = 30, rate = 30 / base_mean[j])
  }
  # Cases: first 30 features uniformly boosted → step in canonical order.
  # Controls: no boost → smooth-gradient profile in canonical order.
  X[y == 1L, 1:30] <- X[y == 1L, 1:30] + 200
  list(X = X, y = y)
}

.auc_wilcox <- function(y, score) {
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  r  <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2L) / (n1 * n0)
}


# ---------------------------------------------------------------------------
# 1. Fit->score roundtrip: shape, class, field presence, type
# ---------------------------------------------------------------------------
test_that("inv-glcm fit/score roundtrip has the right shape and class", {
  dat   <- .make_inv_glcm_data()
  model <- fit_inv_glcm(dat$X, dat$y)
  score <- score_inv_glcm(model, dat$X)

  expect_s3_class(model, "inv_glcm_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_setequal(model$order_perm, colnames(dat$X))
  expect_identical(model$K, 8L)
  expect_identical(model$descriptor_dim, 5L)
  expect_length(model$edges, model$K - 1L)
  expect_length(model$head$center, 5L)
  expect_length(model$head$w, 5L)
  expect_length(model$head$scale, 5L)

  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
  expect_null(names(score))
})


# ---------------------------------------------------------------------------
# 2. Discrimination direction: cases score higher, AUC > 0.7
# ---------------------------------------------------------------------------
test_that("inv-glcm separates a planted step-function texture signal", {
  dat   <- .make_inv_glcm_data(n = 200L, p = 60L, seed = 72L)
  tr    <- c(seq_len(70L), 101L:170L)
  te    <- c(71L:100L, 171L:200L)
  model <- fit_inv_glcm(dat$X[tr, ], dat$y[tr])
  score <- score_inv_glcm(model, dat$X[te, ])

  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[te] == 1L]), mean(score[dat$y[te] == 0L]))
  expect_gt(.auc_wilcox(dat$y[te], score), 0.7)
})


# ---------------------------------------------------------------------------
# 3. Row-equivariance gate (singlesample_assert / singlesample_is)
# ---------------------------------------------------------------------------
test_that("inv-glcm passes the canonical row-equivariance gate", {
  dat   <- .make_inv_glcm_data(seed = 73L)
  tr    <- c(seq_len(70L), 101L:170L)
  te    <- c(71L:100L, 171L:200L)
  model <- fit_inv_glcm(dat$X[tr, ], dat$y[tr])
  X_te  <- dat$X[te, ]
  score_fun <- function(model, X, meta) score_inv_glcm(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_te))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_te))
})


# ---------------------------------------------------------------------------
# 4. Exact per-sample scale-invariance (headline property)
# ---------------------------------------------------------------------------
test_that("inv-glcm scoring is exactly invariant to per-sample positive scaling", {
  dat   <- .make_inv_glcm_data(seed = 74L)
  model <- fit_inv_glcm(dat$X, dat$y)
  row   <- dat$X[3, , drop = FALSE]

  expect_equal(score_inv_glcm(model, row * 7),
               score_inv_glcm(model, row), tolerance = 1e-10)
  expect_equal(score_inv_glcm(model, row * 1e6),
               score_inv_glcm(model, row), tolerance = 1e-10)

  batch  <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e5
  expect_equal(score_inv_glcm(model, scaled),
               score_inv_glcm(model, batch), tolerance = 1e-10)

  # Each row scaled by its own positive factor (column-major recycling by nrow).
  mixed <- batch * c(2, 0.5, 10, 1, 1e3, 7)
  expect_equal(score_inv_glcm(model, mixed),
               score_inv_glcm(model, batch), tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# 5. Partial feature overlap + column reorder
# ---------------------------------------------------------------------------
test_that("inv-glcm handles partial feature overlap and column reorder", {
  dat   <- .make_inv_glcm_data(seed = 75L)
  model <- fit_inv_glcm(dat$X, dat$y)

  keep      <- rev(sample(colnames(dat$X), 30L))
  X_part    <- dat$X[1:12, keep, drop = FALSE]
  s_part    <- score_inv_glcm(model, X_part)

  expect_length(s_part, nrow(X_part))
  expect_true(all(is.finite(s_part)))
  expect_equal(score_inv_glcm(model, X_part[5, , drop = FALSE]),
               s_part[5], tolerance = 1e-12)
  X_reperm <- X_part[, sample(ncol(X_part)), drop = FALSE]
  expect_equal(score_inv_glcm(model, X_reperm), s_part, tolerance = 1e-12)
})


# ---------------------------------------------------------------------------
# 6. singlesample_score_call dispatch
# ---------------------------------------------------------------------------
test_that("singlesample_score_call('inv-glcm', ...) equals score_inv_glcm", {
  skip_if_not(exists("singlesample_score_call"))
  dat   <- .make_inv_glcm_data(seed = 76L)
  tr    <- c(seq_len(70L), 101L:170L)
  te    <- c(71L:100L, 171L:200L)
  model <- fit_inv_glcm(dat$X[tr, ], dat$y[tr])
  X_te  <- dat$X[te, ]
  direct <- score_inv_glcm(model, X_te)

  roster <- singlesample_method_roster()
  if (!"inv-glcm" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "inv-glcm"
    tmpl$fit_fn    <- "fit_inv_glcm"
    tmpl$score_fn  <- "score_inv_glcm"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_inv_glcm", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "inv-glcm",
      function(model, X, meta) score_inv_glcm(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("inv-glcm", model, X_te, roster = roster),
    direct,
    tolerance = 1e-12
  )
})


# ---------------------------------------------------------------------------
# 7. Neutral 0 below min_features
# ---------------------------------------------------------------------------
test_that("inv-glcm returns neutral 0 below the feature-overlap floor", {
  dat   <- .make_inv_glcm_data(seed = 77L)
  model <- fit_inv_glcm(dat$X, dat$y)   # default min_features = 8

  # Fewer than min_features model-universe features present -> neutral 0.
  X_few <- dat$X[1:6, colnames(dat$X)[1:5], drop = FALSE]
  expect_equal(score_inv_glcm(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  # No shared features -> neutral 0.
  X_none <- matrix(1, 4L, 3L,
                   dimnames = list(paste0("N", seq_len(4L)),
                                   c("other-a", "other-b", "other-c")))
  expect_equal(score_inv_glcm(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)

  # Universe-level overlap OK but one row has too few positive values.
  Xs <- dat$X[1:5, , drop = FALSE]
  Xs[2, ]    <- 0
  Xs[2, 1:4] <- dat$X[2, 1:4]   # 4 nonzero < min_features 8
  s <- score_inv_glcm(model, Xs)
  expect_equal(s[2], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s)))
  expect_true(any(s[-2] != 0))
})


# ---------------------------------------------------------------------------
# 8. Determinism and RNG safety
# ---------------------------------------------------------------------------
test_that("inv-glcm fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_inv_glcm_data(seed = 78L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_inv_glcm(dat$X, dat$y)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_inv_glcm(dat$X, dat$y)
  expect_identical(model_1, model_2)

  s1 <- score_inv_glcm(model_1, dat$X[1:20, ])
  s2 <- score_inv_glcm(model_2, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
})


# ---------------------------------------------------------------------------
# 9. hp validation errors
# ---------------------------------------------------------------------------
test_that("inv-glcm hp validation errors are explicit", {
  dat <- .make_inv_glcm_data(seed = 79L)

  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(foo = 1)), "unknown hp")
  # Unnamed hp: names(hp) == NULL -> must be rejected.
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(8L)),  "must be named")
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(n_levels = 4L, 99)),
               "must be named")
  # Duplicated field.
  expect_error(
    fit_inv_glcm(dat$X, dat$y, hp = list(n_levels = 4L, n_levels = 8L)),
    "duplicated hp")
  # n_levels out of range.
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(n_levels = 1L)), "n_levels")
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(n_levels = 33L)), "n_levels")
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(n_levels = 4.5)), "n_levels")
  # Integer overflow guard.
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(n_levels = 3e9)), "n_levels")
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  # shrink negative.
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(shrink = -0.1)), "shrink")
  # eps <= 0.
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  # seed negative.
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(seed = -1)), "hp\\$seed")
  # Partial match must be rejected.
  expect_error(fit_inv_glcm(dat$X, dat$y, hp = list(min = 3L)), "unknown hp")
})


# ---------------------------------------------------------------------------
# 10. Single-row scoring equals batch
# ---------------------------------------------------------------------------
test_that("inv-glcm scores a single row identically standalone and in batch", {
  dat   <- .make_inv_glcm_data(seed = 80L)
  model <- fit_inv_glcm(dat$X, dat$y)
  s_all <- score_inv_glcm(model, dat$X[1:8, ])
  for (i in c(1L, 4L, 8L)) {
    s1 <- score_inv_glcm(model, dat$X[i, , drop = FALSE])
    expect_length(s1, 1L)
    expect_equal(s1, s_all[i], tolerance = 1e-12)
  }
})


# ---------------------------------------------------------------------------
# 11. GLCM + Haralick CORRECTNESS vs hand-computed ground truth (non-circular)
# ---------------------------------------------------------------------------
test_that("inv-glcm GLCM and Haralick features match hand-computed ground truth", {
  # Hand example: q = c(1,1,2,2,1) with K=2
  # Adjacencies (distance 1): (1,1),(1,2),(2,2),(2,1) -> four pairs.
  # Raw C (before symmetrisation): C[1,1]=1, C[1,2]=1, C[2,2]=1, C[2,1]=1
  # After symmetrise (C + t(C)): C[1,1]=2, C[1,2]=2, C[2,1]=2, C[2,2]=2
  # sum(C) = 8 -> P = C/8 = [[0.25,0.25],[0.25,0.25]]
  q_hand <- c(1L, 1L, 2L, 2L, 1L)
  P_hand <- .inv_glcm_build_P(q_hand, K = 2L)
  expect_equal(P_hand, matrix(0.25, 2, 2), tolerance = 1e-12)
  expect_equal(sum(P_hand), 1, tolerance = 1e-12)

  # Haralick for P = [[0.25,0.25],[0.25,0.25]], K=2, i,j in {1,2}:
  # contrast    = (1-1)^2*0.25 + (1-2)^2*0.25 + (2-1)^2*0.25 + (2-2)^2*0.25 = 0.5
  # homogeneity = 0.25/1 + 0.25/2 + 0.25/2 + 0.25/1 = 0.75
  # energy      = 4 * 0.0625 = 0.25
  # entropy     = -4 * (0.25 * log(0.25)) = log(4)
  # correlation: mu_i=mu_j=1.5, sig_i=sig_j=0.5,
  #              cross_mean = 1*1*0.25+1*2*0.25+2*1*0.25+2*2*0.25 = 2.25
  #              corr = (2.25 - 1.5*1.5)/(0.5*0.5) = 0
  h_hand <- .inv_glcm_haralick(P_hand, K = 2L, eps = 1e-6)
  expect_equal(h_hand[1], 0.5,    tolerance = 1e-12)  # contrast
  expect_equal(h_hand[2], 0.75,   tolerance = 1e-12)  # homogeneity
  expect_equal(h_hand[3], 0.25,   tolerance = 1e-12)  # energy
  expect_equal(h_hand[4], log(4), tolerance = 1e-12)  # entropy
  expect_equal(h_hand[5], 0,      tolerance = 1e-12)  # correlation

  # Constant level sequence: q = c(1,1,1,1), K=2
  # Adjacencies: (1,1),(1,1),(1,1). C[1,1]=3; symmetrised=6; P[1,1]=1.
  # contrast 0, energy 1, entropy 0, homogeneity 1, correlation 0 (guard).
  P_const <- .inv_glcm_build_P(c(1L, 1L, 1L, 1L), K = 2L)
  expect_equal(P_const[1, 1], 1,  tolerance = 1e-12)
  expect_equal(sum(P_const), 1,   tolerance = 1e-12)
  h_const <- .inv_glcm_haralick(P_const, K = 2L, eps = 1e-6)
  expect_equal(h_const[1], 0,     tolerance = 1e-12)  # contrast 0
  expect_equal(h_const[3], 1,     tolerance = 1e-12)  # energy 1
  expect_equal(h_const[4], 0,     tolerance = 1e-12)  # entropy 0
  expect_equal(h_const[5], 0,     tolerance = 1e-12)  # correlation 0 (guard)
})


# ---------------------------------------------------------------------------
# 12. Quantiser correctness and tied-edge robustness
# ---------------------------------------------------------------------------
test_that("inv-glcm quantiser places values correctly and handles tied edges", {
  edges <- c(-1, 0, 1)
  q <- .inv_glcm_quantise(c(-2, -1, -0.5, 0, 0.5, 1, 2), edges, K = 4L)
  expect_identical(q, c(1L, 2L, 2L, 3L, 3L, 4L, 4L))

  # Tied edges: findInterval still returns valid levels in {1..K}.
  edges_tied <- c(1.0, 1.0, 1.0)
  q_tied <- .inv_glcm_quantise(c(0.5, 1.0, 1.5), edges_tied, K = 4L)
  expect_true(all(q_tied >= 1L & q_tied <= 4L))

  # K=1: no edges, every value -> level 1.
  q_k1 <- .inv_glcm_quantise(c(-5, 0, 100), numeric(0), K = 1L)
  expect_identical(q_k1, c(1L, 1L, 1L))
})


# ---------------------------------------------------------------------------
# 13. Ridge-LDA head: normal-equations recompute + PD on near-degenerate Phi
# ---------------------------------------------------------------------------
test_that("inv-glcm head solves ridge-LDA normal equations correctly", {
  set.seed(5)
  Phi <- matrix(stats::rnorm(40 * 5), 40, 5)
  y   <- rep(c(0L, 1L), each = 20L)
  Phi[y == 1L, ] <- Phi[y == 1L, ] + 0.8
  head <- .inv_glcm_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)

  center <- colMeans(Phi)
  raw_sd <- apply(Phi, 2L, stats::sd)
  active <- raw_sd > 0
  scale  <- pmax(raw_sd, 1e-6)
  Z  <- sweep(sweep(Phi, 2L, center, "-"), 2L, scale, "/")
  Z[, !active] <- 0
  mu1 <- colMeans(Z[y == 1L, , drop = FALSE])
  mu0 <- colMeans(Z[y == 0L, , drop = FALSE])
  df  <- nrow(Z) - 2L
  Sw  <- (crossprod(sweep(Z[y == 1L, , drop = FALSE], 2L, mu1, "-")) +
          crossprod(sweep(Z[y == 0L, , drop = FALSE], 2L, mu0, "-"))) / df
  Sw  <- (Sw + t(Sw)) / 2
  Swr <- Sw; diag(Swr) <- diag(Swr) + head$ridge
  w_check <- as.numeric(chol2inv(chol(Swr)) %*% (mu1 - mu0))
  b_check <- -0.5 * sum((mu1 + mu0) * w_check)

  expect_equal(head$w, w_check, tolerance = 1e-8)
  expect_equal(head$b, b_check, tolerance = 1e-8)
})

test_that("inv-glcm head stays PD after ridge on a near-degenerate/collinear Phi", {
  set.seed(9)
  Phi <- matrix(stats::rnorm(10 * 5, sd = 1e-7), 10, 5)
  Phi[, 1] <- Phi[, 2]   # exact collinearity
  y <- rep(c(0L, 1L), each = 5L)
  head <- .inv_glcm_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)
  expect_true(all(is.finite(head$w)))
  expect_true(is.finite(head$b))
  expect_true(head$ridge > 0)
})


# ---------------------------------------------------------------------------
# 14. rCLR landscape: v>0 keying, structural-zero exclusion
# ---------------------------------------------------------------------------
test_that("inv-glcm rCLR landscape keys on v>0, not g!=0", {
  # v = (1, 2, 4): all 3 positive. Middle rCLR = 0.
  g <- .inv_glcm_rclr_landscape(c(1, 2, 4))
  expect_length(g, 3L)
  expect_equal(sort(g), c(-log(2), 0, log(2)), tolerance = 1e-12)
  # Structural zeros excluded; present z==0 kept.
  expect_length(.inv_glcm_rclr_landscape(c(0, 1, 2, 4, 0)), 3L)
  # Constant positive: all rCLR = 0, length = #positive.
  g_const <- .inv_glcm_rclr_landscape(rep(3, 7))
  expect_length(g_const, 7L)
  expect_true(all(g_const == 0))

  # End-to-end: a specimen constant across all present features scores finite.
  dat <- .make_inv_glcm_data(seed = 80L)
  model <- fit_inv_glcm(dat$X, dat$y)
  Xc <- dat$X[1:3, , drop = FALSE]
  Xc[1, ] <- 5
  expect_true(all(is.finite(score_inv_glcm(model, Xc))))
})
