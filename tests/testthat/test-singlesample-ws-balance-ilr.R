library(testthat)

.make_ws_balance_ilr_data <- function(n = 100L, seed = 41L, flip = FALSE) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  sbp <- ws_default_sbp()
  sentinels <- c("ALL_NON_RBC", "ALL_NON_PLT", "TOP_K_BY_ABUNDANCE",
                 "TAIL_BY_ABUNDANCE")
  canonical <- unique(unlist(lapply(sbp, function(b) {
    c(b$numerator, b$denominator)
  }), use.names = FALSE))
  canonical <- canonical[!canonical %in% sentinels]
  noise <- paste0("noise", sprintf("%02d", seq_len(12L)))
  features <- c(canonical, noise)

  X <- matrix(stats::rgamma(n * length(features), shape = 40, rate = 2),
              nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))

  numerator_signal <- unique(c(
    sbp$rbc_vs_rest$numerator,
    sbp$platelet_vs_rest$numerator,
    sbp$let7_a_vs_let7_g$numerator,
    sbp$let7_canonical_vs_b_d_f$numerator,
    sbp$mir17_vs_mir92$numerator,
    sbp$mir200_vs_mir141$numerator,
    sbp$mir371_373_vs_rest_germcell$numerator
  ))
  denominator_signal <- unique(c(
    sbp$let7_a_vs_let7_g$denominator,
    sbp$let7_canonical_vs_b_d_f$denominator,
    sbp$mir17_vs_mir92$denominator,
    sbp$mir200_vs_mir141$denominator,
    sbp$mir371_373_vs_rest_germcell$denominator
  ))

  case_high <- if (flip) denominator_signal else numerator_signal
  control_high <- if (flip) numerator_signal else denominator_signal
  n0 <- sum(y == 0L)
  n1 <- sum(y == 1L)

  X[y == 1L, case_high] <- 120 +
    matrix(stats::runif(n1 * length(case_high), 0, 2), nrow = n1)
  X[y == 1L, control_high] <- 8 +
    matrix(stats::runif(n1 * length(control_high), 0, 2), nrow = n1)
  X[y == 0L, case_high] <- 8 +
    matrix(stats::runif(n0 * length(case_high), 0, 2), nrow = n0)
  X[y == 0L, control_high] <- 120 +
    matrix(stats::runif(n0 * length(control_high), 0, 2), nrow = n0)

  list(X = X, y = y, canonical = canonical)
}

.auc_or_wilcoxon_balance_ilr <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("ws-balance-ilr fit/score separates planted ILR balance signal", {
  dat <- .make_ws_balance_ilr_data()
  train <- c(seq_len(40L), 51L:90L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_ws_balance_ilr(dat$X[train, ], dat$y[train])
  score <- score_ws_balance_ilr(model, dat$X[test, ])

  B <- ws_balance_ilr(
    dat$X[test, model$feature_universe, drop = FALSE],
    balances = model$sbp,
    pseudocount = model$pseudocount,
    aggregate = model$aggregate,
    min_balance_coverage = model$min_balance_coverage
  )
  valid_cols <- which(colSums(!is.na(B)) > 0L)
  expected <- model$sign * rowSums(B[, valid_cols, drop = FALSE],
                                   na.rm = TRUE)

  expect_s3_class(model, "ws_balance_ilr_model")
  expect_equal(model$sign, 1)
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
  expect_equal(score, as.numeric(expected), tolerance = 1e-12)
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_balance_ilr(dat$y[test], score), 0.7)
})

test_that("ws-balance-ilr passes the canonical row-equivariance gate", {
  dat <- .make_ws_balance_ilr_data(seed = 42L)
  model <- fit_ws_balance_ilr(dat$X[1:80, ], dat$y[1:80])
  X_test <- dat$X[81:100, ]
  score_fun <- function(model, X, meta) score_ws_balance_ilr(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("ws-balance-ilr is deployable through the canonical dispatch", {
  dat <- .make_ws_balance_ilr_data(seed = 42L)
  model <- fit_ws_balance_ilr(dat$X[1:80, ], dat$y[1:80])
  X_test <- dat$X[81:100, ]
  # singlesample_score_call (the interface the benchmark uses) must equal the
  # direct call -- the manifest score_fn must point at score_ws_balance_ilr,
  # not the ws_balance_ilr primitive (which has an incompatible signature).
  expect_equal(
    singlesample_score_call("ws-balance-ilr", model, X_test),
    score_ws_balance_ilr(model, X_test),
    tolerance = 1e-12
  )
})

test_that("ws-balance-ilr scores are deterministic finite numeric vectors", {
  dat <- .make_ws_balance_ilr_data(seed = 43L)
  model_1 <- fit_ws_balance_ilr(dat$X, dat$y)
  model_2 <- fit_ws_balance_ilr(dat$X, dat$y)
  s1 <- score_ws_balance_ilr(model_1, dat$X[1:20, ])
  s2 <- score_ws_balance_ilr(model_1, dat$X[1:20, ])

  expect_identical(model_1, model_2)
  expect_identical(s1, s2)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("ws-balance-ilr handles no-coverage panels as finite neutral scores", {
  dat <- .make_ws_balance_ilr_data(seed = 44L)
  model <- fit_ws_balance_ilr(dat$X, dat$y)

  X_noise <- dat$X[1:6, "noise01", drop = FALSE]
  s_noise <- score_ws_balance_ilr(model, X_noise)
  expect_equal(s_noise, rep(0, nrow(X_noise)))

  X_no_shared <- matrix(1, nrow = 4L, ncol = 1L,
                        dimnames = list(paste0("N", seq_len(4L)), "other"))
  expect_equal(score_ws_balance_ilr(model, X_no_shared),
               rep(0, nrow(X_no_shared)))
})

test_that("ws-balance-ilr handles partial SBP overlap without batch coupling", {
  dat <- .make_ws_balance_ilr_data(seed = 45L)
  model <- fit_ws_balance_ilr(dat$X, dat$y)
  drop_features <- c("hsa-miR-144-3p", "hsa-let-7a-3p",
                     "hsa-miR-302d-3p")
  X_partial <- dat$X[1:12, setdiff(colnames(dat$X), drop_features),
                     drop = FALSE]

  s_partial <- score_ws_balance_ilr(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  expect_true(any(abs(s_partial) > 0))
  expect_equal(score_ws_balance_ilr(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
})

test_that("ws-balance-ilr relearns orientation when raw direction is flipped", {
  dat <- .make_ws_balance_ilr_data(seed = 46L, flip = TRUE)
  train <- c(seq_len(40L), 51L:90L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_ws_balance_ilr(dat$X[train, ], dat$y[train])
  score <- score_ws_balance_ilr(model, dat$X[test, ])

  expect_equal(model$sign, -1)
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_balance_ilr(dat$y[test], score), 0.7)
})

test_that("ws-balance-ilr input validation errors are explicit", {
  dat <- .make_ws_balance_ilr_data(seed = 47L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_ws_balance_ilr(X_unnamed, dat$y), "feature names")

  expect_error(fit_ws_balance_ilr(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_ws_balance_ilr(dat$X, rep(0, nrow(dat$X))),
               "at least one case")
  expect_error(fit_ws_balance_ilr(dat$X, dat$y, hp = list(maxRank = 10L)),
               "unknown hp")

  # hp-field validation
  expect_error(fit_ws_balance_ilr(dat$X, dat$y, hp = list(pseudocount = 0)),
               "pseudocount")
  expect_error(fit_ws_balance_ilr(dat$X, dat$y, hp = list(aggregate = "bad")),
               "aggregate|should be one of")
  expect_error(
    fit_ws_balance_ilr(dat$X, dat$y, hp = list(min_balance_coverage = 1.5)),
    "min_balance_coverage")
  expect_error(
    fit_ws_balance_ilr(dat$X, dat$y, hp = list(min_sbp_features = 2.5)),
    "min_sbp_features")
  expect_error(fit_ws_balance_ilr(dat$X, dat$y, hp = "notalist"),
               "hp must be a list")

  # malformed SBP: empty side and numerator/denominator overlap are rejected.
  bad_empty <- list(b1 = list(numerator = character(0),
                              denominator = "hsa-miR-16-5p"))
  expect_error(fit_ws_balance_ilr(dat$X, dat$y, hp = list(sbp = bad_empty)),
               "non-empty")
  bad_overlap <- list(b1 = list(numerator = c("hsa-miR-16-5p", "hsa-miR-21-5p"),
                                denominator = "hsa-miR-16-5p"))
  expect_error(fit_ws_balance_ilr(dat$X, dat$y, hp = list(sbp = bad_overlap)),
               "overlapping")
  # the default SBP (sentinel denominators) must still validate.
  expect_s3_class(fit_ws_balance_ilr(dat$X, dat$y,
                                     hp = list(sbp = ws_default_sbp())),
                  "ws_balance_ilr_model")

  # an all-zero specimen row scores finite (0), not NaN/NA.
  model <- fit_ws_balance_ilr(dat$X, dat$y)
  X_zero <- dat$X[1:3, , drop = FALSE]
  X_zero[1, ] <- 0
  s_zero <- score_ws_balance_ilr(model, X_zero)
  expect_true(all(is.finite(s_zero)))
})
