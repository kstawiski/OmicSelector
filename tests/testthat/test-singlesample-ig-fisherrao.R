library(testthat)

.make_ig_fisherrao_data <- function(n = 120L, p = 40L, seed = 71L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 30, rate = 1), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case_high <- features[1:8]
  control_high <- features[9:16]
  n0 <- sum(y == 0L)
  n1 <- sum(y == 1L)

  X[y == 1L, case_high] <- 140 +
    matrix(stats::runif(n1 * length(case_high), 0, 3), nrow = n1)
  X[y == 1L, control_high] <- 8 +
    matrix(stats::runif(n1 * length(control_high), 0, 3), nrow = n1)
  X[y == 0L, case_high] <- 8 +
    matrix(stats::runif(n0 * length(case_high), 0, 3), nrow = n0)
  X[y == 0L, control_high] <- 140 +
    matrix(stats::runif(n0 * length(control_high), 0, 3), nrow = n0)

  list(X = X, y = y, case_high = case_high, control_high = control_high)
}

.auc_or_wilcoxon_ig <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("ig-fisherrao fit/score separates planted compositional signal", {
  dat <- .make_ig_fisherrao_data()
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_ig_fisherrao(
    dat$X[train, ], dat$y[train],
    hp = list(max_anchors_per_class = 20L, seed = 17L)
  )
  score <- score_ig_fisherrao(model, dat$X[test, ])

  expect_s3_class(model, "ig_fisherrao_model")
  expect_equal(nrow(model$case_anchors), 20L)
  expect_equal(nrow(model$control_anchors), 20L)
  # Class means are unit vectors on the sphere; dispersions are positive.
  expect_equal(sqrt(sum(model$mu_case^2)), 1, tolerance = 1e-10)
  expect_equal(sqrt(sum(model$mu_control^2)), 1, tolerance = 1e-10)
  expect_gt(model$sigma2_case, 0)
  expect_gt(model$sigma2_control, 0)
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_ig(dat$y[test], score), 0.7)
})

test_that("ig-fisherrao geodesic geometry obeys the Fisher-Rao identities", {
  geo <- OmicSelector:::.ig_fisherrao_geodesic
  emb <- OmicSelector:::.ig_fisherrao_sqrt_embed
  # Distance of a composition to itself is 0; to a disjoint-support composition
  # the inner product (Bhattacharyya) is 0 so the distance is the maximum pi.
  p <- emb(c(a = 2, b = 1, c = 1, d = 0), pc = 0)
  q <- emb(c(a = 0, b = 0, c = 0, d = 5), pc = 0)
  expect_equal(geo(p, p), 0, tolerance = 1e-12)
  expect_equal(geo(p, q), pi, tolerance = 1e-12)
  # Embedding is a unit vector and exactly scale-invariant.
  expect_equal(sum(p^2), 1, tolerance = 1e-12)
  expect_equal(emb(c(a = 4, b = 2, c = 2, d = 0), pc = 0), p, tolerance = 1e-12)
})

test_that("ig-fisherrao passes the canonical row-equivariance gate", {
  dat <- .make_ig_fisherrao_data(seed = 72L)
  model <- fit_ig_fisherrao(
    dat$X[1:90, ], dat$y[1:90],
    hp = list(max_anchors_per_class = 18L, seed = 19L)
  )
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_ig_fisherrao(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("ig-fisherrao shared_dispersion path is also single-sample deployable", {
  dat <- .make_ig_fisherrao_data(seed = 73L)
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)  # both classes held out
  model <- fit_ig_fisherrao(
    dat$X[train, ], dat$y[train],
    hp = list(shared_dispersion = TRUE, max_anchors_per_class = 18L, seed = 19L)
  )
  X_test <- dat$X[test, ]
  score_fun <- function(model, X, meta) score_ig_fisherrao(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  score <- score_ig_fisherrao(model, X_test)
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
})

test_that("ig-fisherrao is deployable through the canonical dispatch", {
  dat <- .make_ig_fisherrao_data(seed = 74L)
  model <- fit_ig_fisherrao(
    dat$X[1:90, ], dat$y[1:90],
    hp = list(max_anchors_per_class = 18L, seed = 19L)
  )
  X_test <- dat$X[91:120, ]

  expect_equal(
    singlesample_score_call("ig-fisherrao", model, X_test),
    score_ig_fisherrao(model, X_test),
    tolerance = 1e-12
  )
})

test_that("ig-fisherrao fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_ig_fisherrao_data(seed = 75L)
  hp <- list(max_anchors_per_class = 10L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_ig_fisherrao(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_ig_fisherrao(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  expect_identical(model_1, model_2)

  s1 <- score_ig_fisherrao(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_ig_fisherrao(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("ig-fisherrao fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_ig_fisherrao_data(seed = 76L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  model <- fit_ig_fisherrao(
    dat$X, dat$y,
    hp = list(max_anchors_per_class = 10L, seed = 29L)
  )

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "ig_fisherrao_model")
})

test_that("ig-fisherrao scoring is exactly scale-invariant at the default pc", {
  # Closing to proportions before the square-root embedding makes the score
  # exactly invariant to per-sample scaling at the default pc=0 (the pc>0 case,
  # where the pseudocount smooths the already-closed proportions, is covered by
  # the next test), at ordinary magnitudes -- not only in a contrived regime.
  dat <- .make_ig_fisherrao_data(seed = 77L)
  model <- fit_ig_fisherrao(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_ig_fisherrao(model, row * 7),
               score_ig_fisherrao(model, row), tolerance = 1e-8)
  expect_equal(score_ig_fisherrao(model, row * 1e6),
               score_ig_fisherrao(model, row), tolerance = 1e-8)
  # The batch is scored row-equivariantly: scaling some rows leaves others fixed.
  batch <- dat$X[1:5, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 4L), ] <- scaled[c(2L, 4L), ] * 1e6
  expect_equal(score_ig_fisherrao(model, scaled),
               score_ig_fisherrao(model, batch), tolerance = 1e-8)

  # An all-zero specimen (closure -> uniform composition) scores finite.
  X_zero <- dat$X[1:3, , drop = FALSE]
  X_zero[1, ] <- 0
  expect_true(all(is.finite(score_ig_fisherrao(model, X_zero))))
})

test_that("ig-fisherrao scoring stays exactly scale-invariant when pc > 0", {
  # The pseudocount smooths the CLOSED proportions (a fixed function of the
  # scale-invariant p), so exact per-sample scale-invariance is preserved for any
  # pc. A pc-before-closure construction would FAIL this test.
  dat <- .make_ig_fisherrao_data(seed = 80L)
  model <- fit_ig_fisherrao(dat$X, dat$y, hp = list(pc = 0.5))
  row <- dat$X[4, , drop = FALSE]
  expect_equal(score_ig_fisherrao(model, row * 7),
               score_ig_fisherrao(model, row), tolerance = 1e-8)
  expect_equal(score_ig_fisherrao(model, row * 1e6),
               score_ig_fisherrao(model, row), tolerance = 1e-8)
  # The partial-overlap path re-derives the smoothed representation consistently
  # and stays scale-invariant too.
  keep <- c(dat$case_high[1:3], dat$control_high[1:3], colnames(dat$X)[20:24])
  Xp <- dat$X[1:6, keep, drop = FALSE]
  expect_equal(score_ig_fisherrao(model, Xp * 50),
               score_ig_fisherrao(model, Xp), tolerance = 1e-8)
})

test_that("ig-fisherrao handles partial feature overlap as finite neutral-aware scores", {
  dat <- .make_ig_fisherrao_data(seed = 78L)
  model <- fit_ig_fisherrao(dat$X, dat$y)

  keep <- c(dat$case_high[1:4], dat$control_high[1:4],
            colnames(dat$X)[25:30])
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_ig_fisherrao(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  expect_equal(score_ig_fisherrao(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)

  X_one <- dat$X[1:6, dat$case_high[1], drop = FALSE]
  expect_equal(score_ig_fisherrao(model, X_one),
               rep(0, nrow(X_one)), tolerance = 1e-12)

  X_no_shared <- matrix(1, nrow = 4L, ncol = 1L,
                        dimnames = list(paste0("N", seq_len(4L)), "other"))
  expect_equal(score_ig_fisherrao(model, X_no_shared),
               rep(0, nrow(X_no_shared)), tolerance = 1e-12)
})

test_that("ig-fisherrao input validation errors are explicit", {
  dat <- .make_ig_fisherrao_data(seed = 79L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_ig_fisherrao(X_unnamed, dat$y), "feature names")

  expect_error(fit_ig_fisherrao(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_ig_fisherrao(dat$X, rep(0, nrow(dat$X))),
               "at least one case")
  expect_error(fit_ig_fisherrao(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")

  expect_error(fit_ig_fisherrao(dat$X, dat$y, hp = list(pc = -1)),
               "hp\\$pc")
  expect_error(fit_ig_fisherrao(dat$X, dat$y, hp = list(eps = 0)),
               "hp\\$eps")
  expect_error(
    fit_ig_fisherrao(dat$X, dat$y, hp = list(shared_dispersion = "yes")),
    "shared_dispersion"
  )
  expect_error(
    fit_ig_fisherrao(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
    "max_anchors_per_class"
  )
  # Out-of-integer-range values must hit the explicit validation error, not a
  # generic as.integer() overflow ("missing value where TRUE/FALSE needed").
  expect_error(
    fit_ig_fisherrao(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
    "max_anchors_per_class"
  )
  expect_error(
    fit_ig_fisherrao(dat$X, dat$y, hp = list(min_features = 3e9)),
    "min_features"
  )
})
