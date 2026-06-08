library(testthat)

# ---------------------------------------------------------------------------
# Synthetic data: controls have a smooth monotone rCLR landscape when read in
# the frozen training-mean order. Cases receive an alternating perturbation
# inside one contiguous mid-block, creating secondary valleys on that same
# canonical path. The signal is within-sample shape, not row scale.
# ---------------------------------------------------------------------------
.make_tda_ph_data <- function(n = 160L, p = 72L, seed = 71L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  base_mu <- seq(10, 2, length.out = p)
  block <- 24:52
  pattern <- rep(c(0.45, -0.80, 0.35, -0.70), length.out = length(block))
  X <- matrix(0, nrow = n, ncol = p,
              dimnames = list(paste0("S", seq_len(n)), features))
  for (i in seq_len(n)) {
    mu <- base_mu
    if (y[i] == 1L) mu[block] <- mu[block] + pattern
    X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.015))
  }
  list(X = X, y = y, block = features[block], features = features)
}

.auc_or_wilcoxon_tda <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.expect_empty_diagram <- function(g) {
  d <- .tda_ph_diagram_0d(g)
  expect_equal(nrow(d), 0L)
  expect_identical(colnames(d), c("birth", "death"))
}

.local_min_count_strict <- function(g) {
  sum(vapply(seq_along(g), function(i) {
    nbr <- c(if (i > 1L) i - 1L, if (i < length(g)) i + 1L)
    all(g[i] < g[nbr])
  }, logical(1)))
}


test_that("tda-ph fit/score roundtrip has the right shape and class", {
  dat <- .make_tda_ph_data()
  model <- fit_tda_ph(dat$X, dat$y)
  score <- score_tda_ph(model, dat$X)

  expect_s3_class(model, "tda_ph_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_setequal(model$order_perm, colnames(dat$X))
  expect_identical(model$descriptor_dim, 8L * 8L)
  expect_identical(model$grid$R, 8L)
  expect_length(model$head$center, model$descriptor_dim)
  expect_length(model$head$scale, model$descriptor_dim)
  expect_length(model$head$w, model$descriptor_dim)

  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_null(names(score))
  expect_true(all(is.finite(score)))
  expect_error(score_tda_ph(list(), dat$X), "tda_ph_model")
})


test_that("tda-ph scores are larger for case-like planted valley landscapes", {
  dat <- .make_tda_ph_data(seed = 72L)
  train <- c(1:60, 81:140)
  test <- setdiff(seq_len(nrow(dat$X)), train)
  model <- fit_tda_ph(dat$X[train, ], dat$y[train])
  score <- score_tda_ph(model, dat$X[test, ])

  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1L]), mean(score[dat$y[test] == 0L]))
  expect_gt(.auc_or_wilcoxon_tda(dat$y[test], score), 0.7)
})


test_that("tda-ph passes the canonical row-equivariance gate", {
  dat <- .make_tda_ph_data(seed = 73L)
  model <- fit_tda_ph(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_tda_ph(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})


test_that("tda-ph scoring is exactly invariant to per-sample positive scaling", {
  dat <- .make_tda_ph_data(seed = 74L)
  model <- fit_tda_ph(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_tda_ph(model, row * 7),
               score_tda_ph(model, row), tolerance = 1e-10)
  expect_equal(score_tda_ph(model, row * 1e6),
               score_tda_ph(model, row), tolerance = 1e-10)

  batch <- dat$X[1:8, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e5
  expect_equal(score_tda_ph(model, scaled),
               score_tda_ph(model, batch), tolerance = 1e-10)

  factors <- c(2, 0.5, 10, 1, 1e3, 7, 3.5, 20)
  mixed <- batch * factors
  expect_equal(score_tda_ph(model, mixed),
               score_tda_ph(model, batch), tolerance = 1e-10)
})


test_that("tda-ph handles partial feature overlap and column reordering consistently", {
  dat <- .make_tda_ph_data(seed = 75L)
  model <- fit_tda_ph(dat$X, dat$y)

  X_full <- dat$X[1:12, , drop = FALSE]
  expect_equal(score_tda_ph(model, X_full[, rev(colnames(X_full)), drop = FALSE]),
               score_tda_ph(model, X_full), tolerance = 1e-12)

  keep <- model$order_perm[1:40]
  X_sparse <- X_full
  X_sparse[, setdiff(colnames(X_sparse), keep)] <- 0
  s_full_sparse <- score_tda_ph(model, X_sparse)
  X_subset_reordered <- X_sparse[, rev(keep), drop = FALSE]
  expect_equal(score_tda_ph(model, X_subset_reordered),
               s_full_sparse, tolerance = 1e-12)
  expect_equal(score_tda_ph(model, X_subset_reordered[5, , drop = FALSE]),
               s_full_sparse[5], tolerance = 1e-12)
})


test_that("tda-ph is deployable through singlesample_score_call", {
  skip_if_not(exists("singlesample_score_call"))
  dat <- .make_tda_ph_data(seed = 76L)
  model <- fit_tda_ph(dat$X[1:120, ], dat$y[1:120])
  X_test <- dat$X[121:160, ]
  direct <- score_tda_ph(model, X_test)
  singlesample_register_score_adapter(
    "tda-ph",
    function(model, X, meta) score_tda_ph(model, X, meta)
  )

  expect_equal(singlesample_score_call("tda-ph", model, X_test),
               direct, tolerance = 1e-12)
})


test_that("tda-ph returns neutral 0 below feature and positive-support floors", {
  dat <- .make_tda_ph_data(seed = 77L)
  model <- fit_tda_ph(dat$X, dat$y)

  X_few <- dat$X[1:6, model$order_perm[1:5], drop = FALSE]
  expect_equal(score_tda_ph(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", seq_len(4L)),
                                   c("other-a", "other-b", "other-c")))
  expect_equal(score_tda_ph(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)

  Xs <- dat$X[1:5, , drop = FALSE]
  Xs[2, ] <- 0
  Xs[2, 1:4] <- dat$X[2, 1:4]
  s <- score_tda_ph(model, Xs)
  expect_equal(s[2], 0, tolerance = 1e-12)
  expect_true(all(is.finite(s)))
  expect_true(any(s[-2] != 0))
})


test_that("tda-ph fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_tda_ph_data(seed = 78L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_tda_ph(dat$X, dat$y)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_tda_ph(dat$X, dat$y)
  expect_identical(model_1, model_2)

  s1 <- score_tda_ph(model_1, dat$X[1:20, ])
  s2 <- score_tda_ph(model_2, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
})


test_that("tda-ph hp validation errors are explicit", {
  dat <- .make_tda_ph_data(seed = 79L)

  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(foo = 1)), "unknown hp")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(8L)), "must be named")
  expect_error(fit_tda_ph(dat$X, dat$y,
                          hp = list(grid_res = 8L, 99)), "must be named")
  dup <- structure(list(8L, 9L), names = c("grid_res", "grid_res"))
  expect_error(fit_tda_ph(dat$X, dat$y, hp = dup), "duplicated hp")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(grid_res = 2L)),
               "grid_res")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(grid_res = 17L)),
               "grid_res")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(grid_res = 3e9)),
               "grid_res")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(grid_res = "8")),
               "grid_res")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(sigma_frac = 0)),
               "sigma_frac")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(grid_pad = -0.1)),
               "grid_pad")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(shrink = -0.1)),
               "shrink")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(min_features = 2L)),
               "min_features")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(seed = -1)), "hp\\$seed")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(seed = 3e9)), "hp\\$seed")
  expect_error(fit_tda_ph(dat$X, dat$y, hp = list(grid = 8L)), "unknown hp")
})


test_that("tda-ph scores a single row identically standalone and in batch", {
  dat <- .make_tda_ph_data(seed = 80L)
  model <- fit_tda_ph(dat$X, dat$y)
  s_all <- score_tda_ph(model, dat$X[1:8, ])
  for (i in c(1L, 4L, 8L)) {
    s_one <- score_tda_ph(model, dat$X[i, , drop = FALSE])
    expect_length(s_one, 1L)
    expect_equal(s_one, s_all[i], tolerance = 1e-12)
  }
})


test_that("tda-ph rCLR keeps coordinates at the geometric mean", {
  z <- .tda_ph_rclr_landscape(c(1, 2, 4))
  expect_length(z, 3L)
  expect_equal(sort(z), c(-log(2), 0, log(2)), tolerance = 1e-12)
  expect_length(.tda_ph_rclr_landscape(c(0, 1, 2, 4, 0)), 3L)
  expect_length(.tda_ph_rclr_landscape(rep(2, 9)), 9L)
  expect_true(all(.tda_ph_rclr_landscape(rep(2, 9)) == 0))
})


test_that("tda-ph persistence diagram matches hand-derived path examples", {
  .expect_empty_diagram(c(0, 2))
  .expect_empty_diagram(c(2, 0, 1))

  d1 <- .tda_ph_diagram_0d(c(0, 2, 1, 3))
  expect_equal(d1, matrix(c(1, 2), nrow = 1L,
                          dimnames = list(NULL, c("birth", "death"))),
               tolerance = 1e-12)

  d2 <- .tda_ph_diagram_0d(c(0, 3, 1, 4, 2, 5))
  expected <- matrix(c(1, 2, 3, 4), ncol = 2L,
                     dimnames = list(NULL, c("birth", "death")))
  expect_equal(d2, expected, tolerance = 1e-12)

  set.seed(12)
  for (j in seq_len(10L)) {
    g <- stats::rnorm(20)
    d <- .tda_ph_diagram_0d(g)
    expect_equal(nrow(d), .local_min_count_strict(g) - 1L)
    expect_true(all(d[, "death"] > d[, "birth"]))
  }
})


test_that("tda-ph persistence image matches independent closed-form integrals", {
  diagram <- matrix(c(1.2, 2.6), nrow = 1L,
                    dimnames = list(NULL, c("birth", "death")))
  grid <- list(
    R = 3L,
    sigma = 0.5,
    birth_edges = seq(0, 3, length.out = 4L),
    pers_edges = seq(0, 3, length.out = 4L)
  )
  phi <- .tda_ph_persistence_image(diagram, grid)
  p <- 2.6 - 1.2
  db <- diff(stats::pnorm((grid$birth_edges - 1.2) / grid$sigma))
  dp <- diff(stats::pnorm((grid$pers_edges - p) / grid$sigma))
  expected <- numeric(9L)
  for (a in seq_len(3L)) {
    expected[((a - 1L) * 3L + 1L):(a * 3L)] <- p * db[a] * dp
  }
  expect_equal(phi, expected, tolerance = 1e-12)

  wide_diagram <- matrix(c(0, 1, 10, 21), ncol = 2L,
                         dimnames = list(NULL, c("birth", "death")))
  wide_grid <- list(
    R = 80L,
    sigma = 0.1,
    birth_edges = seq(-10, 11, length.out = 81L),
    pers_edges = seq(0, 40, length.out = 81L)
  )
  expect_equal(sum(.tda_ph_persistence_image(wide_diagram, wide_grid)),
               sum(wide_diagram[, "death"] - wide_diagram[, "birth"]),
               tolerance = 1e-8)
  expect_equal(.tda_ph_persistence_image(.tda_ph_empty_diagram(), grid),
               rep(0, 9L), tolerance = 1e-12)
})


test_that("tda-ph head solves ridge-LDA normal equations", {
  set.seed(5)
  Phi <- matrix(stats::rnorm(40 * 7), 40, 7)
  y <- rep(c(0, 1), each = 20)
  Phi[y == 1L, ] <- Phi[y == 1L, ] + 0.8
  head <- .tda_ph_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)

  center <- colMeans(Phi)
  raw_sd <- apply(Phi, 2L, stats::sd)
  active <- raw_sd > 0
  scale <- pmax(raw_sd, 1e-6)
  Z <- sweep(sweep(Phi, 2L, center, "-"), 2L, scale, "/")
  Z[, !active] <- 0
  mu1 <- colMeans(Z[y == 1L, , drop = FALSE])
  mu0 <- colMeans(Z[y == 0L, , drop = FALSE])
  df <- nrow(Z) - 2L
  Sw <- (crossprod(sweep(Z[y == 1L, , drop = FALSE], 2L, mu1, "-")) +
         crossprod(sweep(Z[y == 0L, , drop = FALSE], 2L, mu0, "-"))) / df
  Sw <- (Sw + t(Sw)) / 2
  Swr <- Sw
  diag(Swr) <- diag(Swr) + head$ridge
  w_check <- as.numeric(chol2inv(chol(Swr)) %*% (mu1 - mu0))
  b_check <- -0.5 * sum((mu1 + mu0) * w_check)
  expect_equal(head$w, w_check, tolerance = 1e-8)
  expect_equal(head$b, b_check, tolerance = 1e-8)

  dat <- .make_tda_ph_data(seed = 81L)
  model <- fit_tda_ph(dat$X, dat$y)
  phi <- .tda_ph_descriptor(dat$X[7, model$order_perm], model$grid,
                            model$hp$min_features)
  z <- (phi - model$head$center) / model$head$scale
  z[!model$head$active] <- 0
  manual <- model$head$b + sum(model$head$w * z)
  expect_equal(score_tda_ph(model, dat$X[7, , drop = FALSE]), manual,
               tolerance = 1e-9)
})


test_that("tda-ph head stays positive definite after ridge on collinear pixels", {
  set.seed(9)
  Phi <- matrix(stats::rnorm(8 * 6, sd = 1e-7), 8, 6)
  Phi[, 1] <- Phi[, 2]
  y <- rep(c(0, 1), each = 4)
  head <- .tda_ph_fit_head(Phi, y, shrink = 0.1, eps = 1e-6)
  expect_true(all(is.finite(head$w)))
  expect_true(is.finite(head$b))
  expect_true(head$ridge > 0)

  dat <- .make_tda_ph_data(n = 12L, p = 40L, seed = 82L)
  model <- fit_tda_ph(dat$X, dat$y)
  sc <- score_tda_ph(model, dat$X)
  expect_length(sc, nrow(dat$X))
  expect_true(all(is.finite(sc)))
})


test_that("tda-ph scores an empty-diagram (monotone) specimen via the standardized head, not the bare intercept", {
  # Regression lock for the documented design: a strictly monotone landscape has
  # an EMPTY finite diagram -> all-zero persistence image, which is STILL passed
  # through the frozen standardize + head (it is not short-circuited). The score
  # is the constant empty-diagram baseline, generally NOT the bare intercept b,
  # and distinct from the hard neutral 0 of the below-min_features abstention.
  dat <- .make_tda_ph_data(seed = 83L)
  model <- fit_tda_ph(dat$X, dat$y)
  # strictly DECREASING abundances in the frozen canonical order -> strictly
  # monotone rCLR landscape -> exactly one local minimum -> empty finite diagram.
  mono <- matrix(exp(seq(6, 1, length.out = ncol(dat$X))), nrow = 1L,
                 dimnames = list("MONO", model$order_perm))
  g <- .tda_ph_rclr_landscape(mono[1, ])
  expect_equal(nrow(.tda_ph_diagram_0d(g)), 0L)               # empty finite diagram
  phi <- .tda_ph_descriptor(mono[1, ], model$grid, model$hp$min_features)
  expect_true(all(phi == 0))                                  # all-zero image
  s <- score_tda_ph(model, mono[, model$order_perm, drop = FALSE])
  expect_length(s, 1L)
  expect_true(is.finite(s))
  # Exactly the frozen head applied to the standardized all-zero image.
  z <- (phi - model$head$center) / model$head$scale
  z[!model$head$active] <- 0
  manual <- model$head$b + sum(model$head$w * z)
  expect_equal(s, manual, tolerance = 1e-9)
  # The empty-diagram baseline is a CONSTANT shared by every monotone specimen
  # (a second, differently-scaled monotone row scores identically).
  mono2 <- matrix(exp(seq(9, 2, length.out = ncol(dat$X))), nrow = 1L,
                  dimnames = list("MONO2", model$order_perm))
  expect_equal(as.numeric(score_tda_ph(model, mono2[, model$order_perm, drop = FALSE])),
               as.numeric(s), tolerance = 1e-9)
  # ... and it is NOT the hard neutral 0 of the below-floor abstention path:
  below <- mono[, model$order_perm[1:4], drop = FALSE]
  expect_equal(as.numeric(score_tda_ph(model, below)), 0)
})


test_that("tda-ph degenerates to a constant-0 head when every training specimen is monotone", {
  # If NO training specimen has any finite class, all training images are zero,
  # every descriptor is deactivated, and the head collapses to w = 0, b = 0.
  p <- 40L
  n <- 20L
  feats <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(0, n, p, dimnames = list(NULL, feats))
  # every specimen strictly decreasing in feature index -> training-mean order is
  # the feature index order -> every landscape is strictly monotone (empty diagram).
  for (i in seq_len(n)) X[i, ] <- exp(seq(6, 1, length.out = p) + i * 1e-3)
  y <- rep(c(0, 1), each = n / 2L)
  model <- fit_tda_ph(X, y)
  expect_true(all(model$head$w == 0))
  expect_equal(model$head$b, 0)
  expect_true(all(score_tda_ph(model, X) == 0))
})
