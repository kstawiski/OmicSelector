library(testthat)

.make_conf_mondrian_data <- function(n = 120L, p = 36L, seed = 611L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 30, rate = 1), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case_high <- features[1:8]
  control_high <- features[9:16]
  n0 <- sum(y == 0L)
  n1 <- sum(y == 1L)

  X[y == 1L, case_high] <- 160 +
    matrix(stats::runif(n1 * length(case_high), 0, 4), nrow = n1)
  X[y == 1L, control_high] <- 6 +
    matrix(stats::runif(n1 * length(control_high), 0, 2), nrow = n1)
  X[y == 0L, case_high] <- 6 +
    matrix(stats::runif(n0 * length(case_high), 0, 2), nrow = n0)
  X[y == 0L, control_high] <- 160 +
    matrix(stats::runif(n0 * length(control_high), 0, 4), nrow = n0)

  list(
    X = X,
    y = y,
    case_high = case_high,
    control_high = control_high
  )
}

.auc_or_wilcoxon_conf <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1L) / 2) / (n1 * n0)
}

.conf_hp <- list(max_anchors_per_class = 20L, min_features = 4L, seed = 17L)

.register_conf_mondrian_staged_adapter <- function() {
  ns <- tryCatch(asNamespace("OmicSelector"), error = function(e) NULL)
  if (!is.null(ns) &&
      exists("score_conf_mondrian", envir = ns, mode = "function",
             inherits = FALSE)) {
    return(invisible(FALSE))
  }
  singlesample_register_score_adapter(
    "conf-mondrian",
    function(model, X, meta) score_conf_mondrian(model, X, meta)
  )
  invisible(TRUE)
}

test_that("conf-mondrian fit/score roundtrip returns finite vector and model dimensions", {
  dat <- .make_conf_mondrian_data()
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_conf_mondrian(dat$X[train, ], dat$y[train], hp = .conf_hp)
  score <- score_conf_mondrian(model, dat$X[test, ])

  expect_s3_class(model, "conf_mondrian_model")
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_equal(nrow(model$case_anchors), .conf_hp$max_anchors_per_class)
  expect_equal(nrow(model$control_anchors), .conf_hp$max_anchors_per_class)
  expect_length(model$repr$mu_case, ncol(dat$X))
  expect_length(model$repr$mu_control, ncol(dat$X))
  expect_length(model$repr$var_case, ncol(dat$X))
  expect_length(model$repr$var_control, ncol(dat$X))
  expect_length(model$repr$S_case, nrow(model$case_anchors))
  expect_length(model$repr$S_control, nrow(model$control_anchors))
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
})

test_that("conf-mondrian scores are larger for planted case-like specimens", {
  dat <- .make_conf_mondrian_data(seed = 612L)
  train <- c(seq_len(45L), 61L:105L)
  test <- setdiff(seq_len(nrow(dat$X)), train)
  model <- fit_conf_mondrian(dat$X[train, ], dat$y[train], hp = .conf_hp)
  score <- score_conf_mondrian(model, dat$X[test, ])

  expect_gt(mean(score[dat$y[test] == 1L]), mean(score[dat$y[test] == 0L]))
  expect_gt(.auc_or_wilcoxon_conf(dat$y[test], score), 0.8)
})

test_that("conf-mondrian passes the canonical row-equivariance gate", {
  dat <- .make_conf_mondrian_data(seed = 613L)
  model <- fit_conf_mondrian(dat$X[1:90, ], dat$y[1:90], hp = .conf_hp)
  X_test <- dat$X[91:120, ]
  score_fun <- function(model, X, meta) score_conf_mondrian(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("conf-mondrian rCLR scoring is exactly invariant to per-sample scaling", {
  dat <- .make_conf_mondrian_data(seed = 614L)
  model <- fit_conf_mondrian(dat$X, dat$y, hp = .conf_hp)
  row <- dat$X[3, , drop = FALSE]

  expect_equal(score_conf_mondrian(model, row * 7),
               score_conf_mondrian(model, row), tolerance = 1e-8)
  expect_equal(score_conf_mondrian(model, row * 1e6),
               score_conf_mondrian(model, row), tolerance = 1e-8)

  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e5
  expect_equal(score_conf_mondrian(model, scaled),
               score_conf_mondrian(model, batch), tolerance = 1e-8)
})

test_that("conf-mondrian handles partial feature overlap and column reordering", {
  dat <- .make_conf_mondrian_data(seed = 615L)
  model <- fit_conf_mondrian(dat$X, dat$y, hp = .conf_hp)
  keep <- c(dat$case_high[1:5], dat$control_high[1:5],
            colnames(dat$X)[25:30])
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  X_reordered <- X_partial[, rev(keep), drop = FALSE]

  s_partial <- score_conf_mondrian(model, X_partial)
  s_reordered <- score_conf_mondrian(model, X_reordered)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  expect_equal(s_reordered, s_partial, tolerance = 1e-12)
  expect_equal(score_conf_mondrian(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
})

test_that("conf-mondrian is deployable through singlesample_score_call", {
  dat <- .make_conf_mondrian_data(seed = 616L)
  model <- fit_conf_mondrian(dat$X[1:90, ], dat$y[1:90], hp = .conf_hp)
  X_test <- dat$X[91:120, ]

  expect_identical(names(formals(score_conf_mondrian))[1:3],
                   c("model", "X", "meta"))
  .register_conf_mondrian_staged_adapter()
  expect_equal(
    singlesample_score_call("conf-mondrian", model, X_test),
    score_conf_mondrian(model, X_test),
    tolerance = 1e-12
  )
})

test_that("conf-mondrian returns neutral zero below the feature-overlap floor", {
  dat <- .make_conf_mondrian_data(seed = 617L)
  model <- fit_conf_mondrian(
    dat$X, dat$y,
    hp = list(max_anchors_per_class = 20L, min_features = 6L, seed = 17L)
  )

  X_few <- dat$X[1:5, dat$case_high[1:5], drop = FALSE]
  expect_equal(score_conf_mondrian(model, X_few), rep(0, nrow(X_few)),
               tolerance = 1e-12)

  X_no_shared <- matrix(1, nrow = 4L, ncol = 2L,
                        dimnames = list(paste0("N", seq_len(4L)),
                                        c("other-a", "other-b")))
  expect_equal(score_conf_mondrian(model, X_no_shared),
               rep(0, nrow(X_no_shared)), tolerance = 1e-12)
})

test_that("conf-mondrian fitting and scoring are deterministic and RNG-safe", {
  dat <- .make_conf_mondrian_data(seed = 618L)
  hp <- list(max_anchors_per_class = 12L, min_features = 4L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_conf_mondrian(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_conf_mondrian(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  expect_identical(model_1, model_2)

  s1 <- score_conf_mondrian(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_conf_mondrian(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("conf-mondrian hp and model validation errors are explicit", {
  dat <- .make_conf_mondrian_data(seed = 619L)

  expect_error(fit_conf_mondrian(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")
  expect_error(
    fit_conf_mondrian(dat$X, dat$y, hp = list(nonconformity = "rank")),
    "nonconformity"
  )
  expect_error(
    fit_conf_mondrian(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
    "max_anchors_per_class"
  )
  expect_error(
    fit_conf_mondrian(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
    "max_anchors_per_class"
  )
  expect_error(fit_conf_mondrian(dat$X, dat$y, hp = list(eps = 0)), "eps")
  expect_error(fit_conf_mondrian(dat$X, dat$y, hp = list(seed = -1)),
               "seed")
  expect_error(score_conf_mondrian(list(), dat$X[1:3, ]),
               "conf_mondrian_model")
})

test_that("conf-mondrian scores a single row identically standalone and in batch", {
  dat <- .make_conf_mondrian_data(seed = 620L)
  model <- fit_conf_mondrian(dat$X[1:90, ], dat$y[1:90], hp = .conf_hp)
  X_test <- dat$X[91:120, ]
  score <- score_conf_mondrian(model, X_test)
  one <- score_conf_mondrian(model, X_test[7, , drop = FALSE])

  expect_type(one, "double")
  expect_length(one, 1L)
  expect_true(is.finite(one))
  expect_equal(one, score[7], tolerance = 1e-12)
})

test_that("conf-mondrian hand-verifies the conformal p-value formula", {
  X <- rbind(
    case_a = c(a = 12, b = 2, c = 1),
    case_b = c(a = 10, b = 2, c = 1),
    ctrl_a = c(a = 1, b = 2, c = 12),
    ctrl_b = c(a = 1, b = 2, c = 10)
  )
  y <- c(1, 1, 0, 0)
  hp <- list(max_anchors_per_class = 10L, min_features = 3L, seed = 1L)
  model <- fit_conf_mondrian(X, y, hp = hp)
  x <- matrix(c(11, 2, 1), nrow = 1L,
              dimnames = list("new_case", colnames(X)))

  z <- .conf_mondrian_rclr_row(x[1, ])
  A_case <- .conf_mondrian_nonconformity(
    z, model$repr$mu_case, model$repr$var_case, model$hp
  )
  A_control <- .conf_mondrian_nonconformity(
    z, model$repr$mu_control, model$repr$var_control, model$hp
  )
  p_case <- (1 + sum(model$repr$S_case >= A_case)) /
    (1 + length(model$repr$S_case))
  p_control <- (1 + sum(model$repr$S_control >= A_control)) /
    (1 + length(model$repr$S_control))
  expected <- log(p_case) - log(p_control)

  expect_equal(score_conf_mondrian(model, x), expected, tolerance = 1e-12)
  expect_equal(p_case, 1, tolerance = 1e-12)
  expect_equal(p_control, 1 / 3, tolerance = 1e-12)
})
