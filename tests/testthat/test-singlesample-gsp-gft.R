library(testthat)

.make_gsp_gft_data <- function(n = 120L, p = 36L, seed = 91L) {
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

.auc_or_wilcoxon_gsp <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2) / (n1 * n0)
}

.gsp_hp <- list(graph_k = 6L, n_modes = 8L, min_features = 4L)

test_that("gsp-gft fit/score roundtrip returns finite vector and model dimensions", {
  dat <- .make_gsp_gft_data()
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_gsp_gft(dat$X[train, ], dat$y[train], hp = .gsp_hp)
  score <- score_gsp_gft(model, dat$X[test, ])

  expect_s3_class(model, "gsp_gft_model")
  expect_equal(model$feature_universe, colnames(dat$X))
  expect_equal(dim(model$U_low), c(ncol(dat$X), .gsp_hp$n_modes))
  expect_length(model$eigvals_low, .gsp_hp$n_modes)
  expect_true(all(diff(model$eigvals_low) >= -1e-10))
  expect_true(model$head$type %in% c("glm", "diag_lda"))
  expect_length(model$head$coef, .gsp_hp$n_modes)
  expect_length(model$head$center, .gsp_hp$n_modes)
  expect_length(model$head$scale, .gsp_hp$n_modes)
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
})

test_that("gsp-gft score direction is larger for planted case-like specimens", {
  dat <- .make_gsp_gft_data(seed = 92L)
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)
  model <- fit_gsp_gft(dat$X[train, ], dat$y[train], hp = .gsp_hp)
  score <- score_gsp_gft(model, dat$X[test, ])

  expect_gt(mean(score[dat$y[test] == 1L]), mean(score[dat$y[test] == 0L]))
  expect_gt(.auc_or_wilcoxon_gsp(dat$y[test], score), 0.8)
})

test_that("gsp-gft passes the canonical row-equivariance gate", {
  dat <- .make_gsp_gft_data(seed = 93L)
  model <- fit_gsp_gft(dat$X[1:90, ], dat$y[1:90], hp = .gsp_hp)
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_gsp_gft(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("gsp-gft partial feature overlap is finite and reorder-consistent", {
  dat <- .make_gsp_gft_data(seed = 94L)
  model <- fit_gsp_gft(dat$X, dat$y, hp = .gsp_hp)
  keep <- c(dat$case_module[1:5], dat$control_module[1:5],
            colnames(dat$X)[25:30])
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  X_reordered <- X_partial[, rev(keep), drop = FALSE]

  s_partial <- score_gsp_gft(model, X_partial)
  s_reordered <- score_gsp_gft(model, X_reordered)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  expect_equal(s_reordered, s_partial, tolerance = 1e-12)
  expect_equal(score_gsp_gft(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
})

test_that("gsp-gft rCLR scoring is exactly invariant to per-sample scaling", {
  dat <- .make_gsp_gft_data(seed = 95L)
  model <- fit_gsp_gft(dat$X, dat$y, hp = .gsp_hp)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_gsp_gft(model, row * 7),
               score_gsp_gft(model, row), tolerance = 1e-8)
  expect_equal(score_gsp_gft(model, row * 1e6),
               score_gsp_gft(model, row), tolerance = 1e-8)

  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e5
  expect_equal(score_gsp_gft(model, scaled),
               score_gsp_gft(model, batch), tolerance = 1e-8)
})

test_that("gsp-gft is deployable through the canonical singlesample_score_call path", {
  dat <- .make_gsp_gft_data(seed = 96L)
  model <- fit_gsp_gft(dat$X[1:90, ], dat$y[1:90], hp = .gsp_hp)
  X_test <- dat$X[91:120, ]

  # Canonical dispatch resolves the manifest score_fn (score_gsp_gft) directly
  # from the OmicSelector namespace -- no adapter registration needed once the
  # method is part of the package. This exercises the real roster->package path.
  expect_equal(
    singlesample_score_call("gsp-gft", model, X_test),
    score_gsp_gft(model, X_test),
    tolerance = 1e-12
  )
})

test_that("gsp-gft returns neutral zero below the minimum feature overlap", {
  dat <- .make_gsp_gft_data(seed = 97L)
  hp <- list(graph_k = 6L, n_modes = 8L, min_features = 6L)
  model <- fit_gsp_gft(dat$X, dat$y, hp = hp)

  X_few <- dat$X[1:5, dat$case_module[1:5], drop = FALSE]
  expect_equal(score_gsp_gft(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  X_no_shared <- matrix(1, nrow = 4L, ncol = 1L,
                        dimnames = list(paste0("N", seq_len(4L)), "other"))
  expect_equal(score_gsp_gft(model, X_no_shared),
               rep(0, nrow(X_no_shared)), tolerance = 1e-12)
})

test_that("gsp-gft fitting and scoring are deterministic", {
  dat <- .make_gsp_gft_data(seed = 98L)
  model_1 <- fit_gsp_gft(dat$X, dat$y, hp = .gsp_hp)
  model_2 <- fit_gsp_gft(dat$X, dat$y, hp = .gsp_hp)
  s1 <- score_gsp_gft(model_1, dat$X[1:20, ])
  s2 <- score_gsp_gft(model_1, dat$X[1:20, ])

  expect_identical(model_1, model_2)
  expect_identical(s1, s2)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("gsp-gft hyperparameter validation rejects unknown and malformed fields", {
  dat <- .make_gsp_gft_data(seed = 99L)

  expect_error(fit_gsp_gft(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")
  expect_error(fit_gsp_gft(dat$X, dat$y, hp = list(cor_method = "kendall")),
               "cor_method")
  expect_error(fit_gsp_gft(dat$X, dat$y, hp = list(graph_k = 1.5)),
               "graph_k")
  expect_error(fit_gsp_gft(dat$X, dat$y, hp = list(graph_k = 3e9)),
               "graph_k")
  expect_error(fit_gsp_gft(dat$X, dat$y, hp = list(eps = 0)),
               "eps")
  expect_error(fit_gsp_gft(dat$X, dat$y, hp = list(1)),
               "hp fields must be named")
})

test_that("gsp-gft single-row scoring equals the corresponding batch score", {
  dat <- .make_gsp_gft_data(seed = 100L)
  model <- fit_gsp_gft(dat$X[1:90, ], dat$y[1:90], hp = .gsp_hp)
  X_test <- dat$X[91:120, ]
  score <- score_gsp_gft(model, X_test)
  one <- score_gsp_gft(model, X_test[7, , drop = FALSE])

  expect_type(one, "double")
  expect_length(one, 1L)
  expect_true(is.finite(one))
  expect_equal(one, score[7], tolerance = 1e-12)
})

test_that("gsp-gft score validates the model class", {
  dat <- .make_gsp_gft_data(seed = 101L)
  expect_error(score_gsp_gft(list(), dat$X[1:3, ]), "gsp_gft_model")
})
