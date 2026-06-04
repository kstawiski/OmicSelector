library(testthat)

.make_man_nystrom_data <- function(n = 120L, p = 36L, seed = 131L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 25, rate = 1), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case_module <- features[1:8]
  control_module <- features[9:16]
  n0 <- sum(y == 0L)
  n1 <- sum(y == 1L)

  X[y == 1L, case_module] <- 180 +
    matrix(stats::runif(n1 * length(case_module), 0, 4), nrow = n1)
  X[y == 1L, control_module] <- 7 +
    matrix(stats::runif(n1 * length(control_module), 0, 2), nrow = n1)
  X[y == 0L, case_module] <- 7 +
    matrix(stats::runif(n0 * length(case_module), 0, 2), nrow = n0)
  X[y == 0L, control_module] <- 180 +
    matrix(stats::runif(n0 * length(control_module), 0, 4), nrow = n0)

  list(
    X = X,
    y = y,
    case_module = case_module,
    control_module = control_module
  )
}

.auc_or_wilcoxon_man <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2) / (n1 * n0)
}

.man_hp <- list(
  n_landmarks = 50L,
  n_components = 6L,
  t = 1L,
  min_features = 4L,
  seed = 17L
)

test_that("man-nystrom fit/score roundtrip returns finite vector and dimensions", {
  dat <- .make_man_nystrom_data()
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_man_nystrom(dat$X[train, ], dat$y[train], hp = .man_hp)
  score <- score_man_nystrom(model, dat$X[test, ])

  expect_s3_class(model, "man_nystrom_model")
  expect_equal(model$feature_universe, colnames(dat$X))
  expect_equal(dim(model$landmarks), c(.man_hp$n_landmarks, ncol(dat$X)))
  expect_gt(model$gamma, 0)
  expect_true(is.finite(model$gamma))
  expect_gt(model$n_components, 0L)
  expect_lte(model$n_components, .man_hp$n_components)
  expect_length(model$geometry$lambda, model$n_components)
  expect_true(model$head$type %in% c("glm", "diag_lda", "constant"))
  expect_length(model$head$coef, model$n_components)
  expect_length(model$head$center, model$n_components)
  expect_length(model$head$scale, model$n_components)
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
})

test_that("man-nystrom score direction is larger for planted case-like specimens", {
  dat <- .make_man_nystrom_data(seed = 132L)
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)
  model <- fit_man_nystrom(dat$X[train, ], dat$y[train], hp = .man_hp)
  score <- score_man_nystrom(model, dat$X[test, ])

  expect_gt(mean(score[dat$y[test] == 1L]), mean(score[dat$y[test] == 0L]))
  expect_gt(.auc_or_wilcoxon_man(dat$y[test], score), 0.8)
})

test_that("man-nystrom passes the canonical row-equivariance gate", {
  dat <- .make_man_nystrom_data(seed = 133L)
  model <- fit_man_nystrom(dat$X[1:90, ], dat$y[1:90], hp = .man_hp)
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_man_nystrom(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("man-nystrom Nystrom coordinates reproduce landmark diffusion coordinates", {
  dat <- .make_man_nystrom_data(seed = 134L)
  hp <- list(n_landmarks = 40L, n_components = 5L, t = 2L,
             min_features = 4L, seed = 23L)
  model <- fit_man_nystrom(dat$X, dat$y, hp = hp)
  geom <- .man_nystrom_geometry(model$landmarks, model$gamma, model$hp)

  expect_gt(length(geom$lambda), 0L)
  for (i in seq_len(nrow(model$landmarks))) {
    emb <- .man_nystrom_embed_row(model$landmarks[i, ], model$landmarks,
                                  model$gamma, geom, model$hp)
    target <- as.numeric((geom$lambda^model$hp$t) * geom$psi[i, ])
    expect_equal(emb, target, tolerance = 1e-8)
  }
})

test_that("man-nystrom partial feature overlap is finite and reorder-consistent", {
  dat <- .make_man_nystrom_data(seed = 135L)
  model <- fit_man_nystrom(dat$X, dat$y, hp = .man_hp)
  keep <- c(dat$case_module[1:5], dat$control_module[1:5],
            colnames(dat$X)[25:30])
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  X_reordered <- X_partial[, rev(keep), drop = FALSE]

  s_partial <- score_man_nystrom(model, X_partial)
  s_reordered <- score_man_nystrom(model, X_reordered)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  expect_equal(s_reordered, s_partial, tolerance = 1e-12)
  expect_equal(score_man_nystrom(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
})

test_that("man-nystrom rCLR scoring is exactly invariant to per-sample scaling", {
  dat <- .make_man_nystrom_data(seed = 136L)
  model <- fit_man_nystrom(dat$X, dat$y, hp = .man_hp)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_man_nystrom(model, row * 7),
               score_man_nystrom(model, row), tolerance = 1e-8)
  expect_equal(score_man_nystrom(model, row * 1e6),
               score_man_nystrom(model, row), tolerance = 1e-8)

  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e5
  expect_equal(score_man_nystrom(model, scaled),
               score_man_nystrom(model, batch), tolerance = 1e-8)
})

test_that("man-nystrom is deployable through the canonical singlesample_score_call path", {
  dat <- .make_man_nystrom_data(seed = 137L)
  model <- fit_man_nystrom(dat$X[1:90, ], dat$y[1:90], hp = .man_hp)
  X_test <- dat$X[91:120, ]

  # In-package, singlesample_score_call() resolves the manifest score_fn
  # (score_man_nystrom) directly from the OmicSelector namespace -- no adapter
  # registration is needed; this exercises the real roster->package path.
  expect_equal(
    singlesample_score_call("man-nystrom", model, X_test),
    score_man_nystrom(model, X_test),
    tolerance = 1e-12
  )
})

test_that("man-nystrom returns neutral zero below the minimum feature overlap", {
  dat <- .make_man_nystrom_data(seed = 138L)
  hp <- list(n_landmarks = 40L, n_components = 5L, min_features = 6L,
             seed = 29L)
  model <- fit_man_nystrom(dat$X, dat$y, hp = hp)

  X_few <- dat$X[1:5, dat$case_module[1:5], drop = FALSE]
  expect_equal(score_man_nystrom(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  X_no_shared <- matrix(1, nrow = 4L, ncol = 1L,
                        dimnames = list(paste0("N", seq_len(4L)), "other"))
  expect_equal(score_man_nystrom(model, X_no_shared),
               rep(0, nrow(X_no_shared)), tolerance = 1e-12)
})

test_that("man-nystrom fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_man_nystrom_data(seed = 139L)
  hp <- list(n_landmarks = 35L, n_components = 5L, min_features = 4L,
             seed = 31L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_man_nystrom(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_man_nystrom(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  expect_identical(model_1, model_2)

  s1 <- score_man_nystrom(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_man_nystrom(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("man-nystrom hyperparameter validation rejects malformed fields", {
  dat <- .make_man_nystrom_data(seed = 140L)

  expect_error(fit_man_nystrom(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")
  expect_error(fit_man_nystrom(dat$X, dat$y, hp = list(gamma = -1)),
               "gamma")
  expect_error(fit_man_nystrom(dat$X, dat$y, hp = list(n_landmarks = 1.5)),
               "n_landmarks")
  expect_error(fit_man_nystrom(dat$X, dat$y, hp = list(n_landmarks = 3e9)),
               "n_landmarks")
  expect_error(fit_man_nystrom(dat$X, dat$y, hp = list(1)),
               "hp fields must be named")
})

test_that("man-nystrom single-row scoring equals the corresponding batch score", {
  dat <- .make_man_nystrom_data(seed = 141L)
  model <- fit_man_nystrom(dat$X[1:90, ], dat$y[1:90], hp = .man_hp)
  X_test <- dat$X[91:120, ]
  score <- score_man_nystrom(model, X_test)
  one <- score_man_nystrom(model, X_test[7, , drop = FALSE])

  expect_type(one, "double")
  expect_length(one, 1L)
  expect_true(is.finite(one))
  expect_equal(one, score[7], tolerance = 1e-12)
})

test_that("man-nystrom score validates the model class", {
  dat <- .make_man_nystrom_data(seed = 142L)
  expect_error(score_man_nystrom(list(), dat$X[1:3, ]), "man_nystrom_model")
})


test_that("man-nystrom degenerates gracefully to a neutral constant head when no diffusion mode survives", {
  dat <- .make_man_nystrom_data(seed = 222L)
  # A high eig_floor retains zero non-trivial diffusion modes. This must NOT
  # crash .man_nystrom_geometry (the 0-column psi path); it must produce the
  # intended constant head and neutral-0 scores, at fit AND at score time.
  model <- fit_man_nystrom(dat$X, dat$y,
                           hp = list(n_landmarks = 20L, n_components = 5L,
                                     eig_floor = 2, min_features = 3L, seed = 1L))
  expect_identical(model$n_components, 0L)
  expect_identical(model$head$type, "constant")
  sc <- score_man_nystrom(model, dat$X[1:10, ])
  expect_length(sc, 10L)
  expect_true(all(is.finite(sc)))
  expect_equal(sc, rep(0, 10L), tolerance = 1e-12)
  # Single-row and dispatch paths also survive the degenerate geometry.
  expect_equal(score_man_nystrom(model, dat$X[3, , drop = FALSE]), 0, tolerance = 1e-12)
  expect_equal(singlesample_score_call("man-nystrom", model, dat$X[1:4, ]),
               rep(0, 4L), tolerance = 1e-12)
})
