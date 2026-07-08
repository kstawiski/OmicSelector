library(testthat)

.make_ot_slicedlrt_data <- function(n = 120L, p = 36L, seed = 171L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(20 + stats::runif(n * p, 0, 3), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case_module <- features[1:8]
  control_module <- features[9:16]
  n0 <- sum(y == 0L)
  n1 <- sum(y == 1L)

  X[y == 1L, case_module] <- 180 +
    matrix(stats::runif(n1 * length(case_module), 0, 4), nrow = n1)
  X[y == 1L, control_module] <- 3 +
    matrix(stats::runif(n1 * length(control_module), 0, 1), nrow = n1)
  X[y == 0L, case_module] <- 3 +
    matrix(stats::runif(n0 * length(case_module), 0, 1), nrow = n0)
  X[y == 0L, control_module] <- 180 +
    matrix(stats::runif(n0 * length(control_module), 0, 4), nrow = n0)

  list(
    X = X,
    y = y,
    case_module = case_module,
    control_module = control_module
  )
}

.auc_or_wilcoxon_ot <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2) / (n1 * n0)
}

.ot_hp <- list(
  n_slices = 80L,
  max_anchors_per_class = 30L,
  min_features = 3L,
  seed = 17L
)

test_that("ot-slicedlrt fit/score roundtrip returns finite vector and dimensions", {
  dat <- .make_ot_slicedlrt_data()
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_ot_slicedlrt(dat$X[train, ], dat$y[train], hp = .ot_hp)
  score <- score_ot_slicedlrt(model, dat$X[test, ])
  repr_again <- .ot_slicedlrt_repr(
    model$case_anchors, model$control_anchors, model$hp, ncol(dat$X)
  )

  expect_s3_class(model, "ot_slicedlrt_model")
  expect_equal(model$feature_universe, colnames(dat$X))
  expect_equal(dim(model$case_anchors), c(.ot_hp$max_anchors_per_class,
                                          ncol(dat$X)))
  expect_equal(dim(model$control_anchors), c(.ot_hp$max_anchors_per_class,
                                             ncol(dat$X)))
  expect_equal(dim(model$repr$Theta), c(.ot_hp$n_slices, ncol(dat$X)))
  expect_equal(rowSums(model$repr$Theta^2), rep(1, .ot_hp$n_slices),
               tolerance = 1e-12)
  expect_length(model$repr$mu_case, .ot_hp$n_slices)
  expect_length(model$repr$v_case, .ot_hp$n_slices)
  expect_true(all(model$repr$v_case >= model$hp$eps))
  expect_true(all(model$repr$v_control >= model$hp$eps))
  expect_equal(model$repr, repr_again, tolerance = 1e-12)
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
})

test_that("ot-slicedlrt score direction is larger for planted case-like specimens", {
  dat <- .make_ot_slicedlrt_data(seed = 172L)
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)
  model <- fit_ot_slicedlrt(dat$X[train, ], dat$y[train], hp = .ot_hp)
  score <- score_ot_slicedlrt(model, dat$X[test, ])

  expect_gt(mean(score[dat$y[test] == 1L]), mean(score[dat$y[test] == 0L]))
  expect_gt(.auc_or_wilcoxon_ot(dat$y[test], score), 0.9)
})

test_that("ot-slicedlrt passes the canonical row-equivariance gate", {
  dat <- .make_ot_slicedlrt_data(seed = 173L)
  model <- fit_ot_slicedlrt(dat$X[1:90, ], dat$y[1:90], hp = .ot_hp)
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_ot_slicedlrt(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("ot-slicedlrt partial feature overlap is finite and reorder-consistent", {
  dat <- .make_ot_slicedlrt_data(seed = 174L)
  model <- fit_ot_slicedlrt(dat$X, dat$y, hp = .ot_hp)
  keep <- c(dat$case_module[1:5], dat$control_module[1:5],
            colnames(dat$X)[25:31])
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  X_reordered <- X_partial[, rev(keep), drop = FALSE]

  s_partial <- score_ot_slicedlrt(model, X_partial)
  s_reordered <- score_ot_slicedlrt(model, X_reordered)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  expect_equal(s_reordered, s_partial, tolerance = 1e-12)
  expect_equal(score_ot_slicedlrt(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
})

test_that("ot-slicedlrt rCLR scoring is invariant to per-sample scaling", {
  dat <- .make_ot_slicedlrt_data(seed = 175L)
  model <- fit_ot_slicedlrt(dat$X, dat$y, hp = .ot_hp)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_ot_slicedlrt(model, row * 7),
               score_ot_slicedlrt(model, row), tolerance = 1e-8)
  expect_equal(score_ot_slicedlrt(model, row * 1e6),
               score_ot_slicedlrt(model, row), tolerance = 1e-8)

  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e5
  expect_equal(score_ot_slicedlrt(model, scaled),
               score_ot_slicedlrt(model, batch), tolerance = 1e-8)
})

test_that("ot-slicedlrt is deployable through the canonical singlesample_score_call path", {
  dat <- .make_ot_slicedlrt_data(seed = 176L)
  model <- fit_ot_slicedlrt(dat$X[1:90, ], dat$y[1:90], hp = .ot_hp)
  X_test <- dat$X[91:120, ]
  # In-package, singlesample_score_call() resolves the manifest score_fn
  # (score_ot_slicedlrt) directly from the OmicSelector namespace -- no adapter
  # registration needed; this exercises the real roster->package dispatch path.
  expect_equal(
    singlesample_score_call("ot-slicedlrt", model, X_test),
    score_ot_slicedlrt(model, X_test),
    tolerance = 1e-12
  )
})

test_that("ot-slicedlrt returns neutral zero below the minimum feature overlap", {
  dat <- .make_ot_slicedlrt_data(seed = 177L)
  hp <- list(n_slices = 40L, max_anchors_per_class = 25L,
             min_features = 6L, seed = 19L)
  model <- fit_ot_slicedlrt(dat$X, dat$y, hp = hp)

  X_few <- dat$X[1:5, dat$case_module[1:5], drop = FALSE]
  expect_equal(score_ot_slicedlrt(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  X_no_shared <- matrix(1, nrow = 4L, ncol = 1L,
                        dimnames = list(paste0("N", seq_len(4L)), "other"))
  expect_equal(score_ot_slicedlrt(model, X_no_shared),
               rep(0, nrow(X_no_shared)), tolerance = 1e-12)
})

test_that("ot-slicedlrt fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_ot_slicedlrt_data(seed = 178L)
  hp <- list(n_slices = 55L, max_anchors_per_class = 12L,
             min_features = 3L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_ot_slicedlrt(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_ot_slicedlrt(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  expect_identical(model_1, model_2)

  s1 <- score_ot_slicedlrt(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_ot_slicedlrt(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
  expect_true(all(is.finite(s1)))

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  model_3 <- fit_ot_slicedlrt(dat$X, dat$y, hp = hp)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_identical(model_1, model_3)
  expect_true(all(is.finite(score_ot_slicedlrt(model_3, dat$X[1:4, ]))))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("ot-slicedlrt hyperparameter validation rejects malformed fields", {
  dat <- .make_ot_slicedlrt_data(seed = 179L)

  expect_error(fit_ot_slicedlrt(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")
  expect_error(fit_ot_slicedlrt(dat$X, dat$y, hp = list(n_slices = 1.5)),
               "n_slices")
  expect_error(fit_ot_slicedlrt(dat$X, dat$y, hp = list(n_slices = 0L)),
               "n_slices")
  expect_error(fit_ot_slicedlrt(dat$X, dat$y, hp = list(n_slices = 3e9)),
               "n_slices")
  expect_error(fit_ot_slicedlrt(dat$X, dat$y, hp = list(eps = 0)),
               "eps")
  expect_error(fit_ot_slicedlrt(dat$X, dat$y, hp = list(seed = -1L)),
               "seed")
  expect_error(fit_ot_slicedlrt(dat$X, dat$y, hp = list(1)),
               "hp fields must be named")
  expect_error(score_ot_slicedlrt(list(), dat$X[1:3, ]),
               "ot_slicedlrt_model")
})

test_that("ot-slicedlrt single-row scoring equals the corresponding batch score", {
  dat <- .make_ot_slicedlrt_data(seed = 180L)
  model <- fit_ot_slicedlrt(dat$X[1:90, ], dat$y[1:90], hp = .ot_hp)
  X_test <- dat$X[91:120, ]
  score <- score_ot_slicedlrt(model, X_test)
  one <- score_ot_slicedlrt(model, X_test[7, , drop = FALSE])

  expect_type(one, "double")
  expect_length(one, 1L)
  expect_true(is.finite(one))
  expect_equal(one, score[7], tolerance = 1e-12)
})

test_that("ot-slicedlrt Gaussian LRT formula matches manual M=2 calculation", {
  repr <- list(
    mu_case = c(1, -1),
    v_case = c(4, 9),
    mu_control = c(0.5, 0),
    v_control = c(2, 3)
  )
  s <- c(2, -2)
  expected <- mean(
    -((s - repr$mu_case)^2) / (2 * repr$v_case) -
      0.5 * log(repr$v_case) +
      ((s - repr$mu_control)^2) / (2 * repr$v_control) +
      0.5 * log(repr$v_control)
  )

  expect_equal(.ot_slicedlrt_gaussian_llr(s, repr), expected,
               tolerance = 1e-12)
})
