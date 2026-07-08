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

.ws_balance_ilr_row_loop_reference <- function(x,
                                               balances = ws_default_sbp(),
                                               pseudocount = NULL,
                                               aggregate = c("gmean", "trimmed_gmean"),
                                               min_balance_coverage = 0.8) {
  aggregate <- match.arg(aggregate)
  rows <- lapply(seq_len(nrow(x)), function(i) {
    ws_balance_ilr(x[i, ], balances = balances, pseudocount = pseudocount,
                   aggregate = aggregate,
                   min_balance_coverage = min_balance_coverage)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- rownames(x)
  coverage <- vapply(rows, function(r) {
    v <- attr(r, "coverage")
    if (is.null(v)) NA_real_ else v
  }, numeric(1L))
  coverage_failed <- vapply(rows, function(r) {
    isTRUE(attr(r, "coverage_failed"))
  }, logical(1L))
  attr(out, "coverage") <- coverage
  attr(out, "coverage_failed") <- coverage_failed
  attr(out, "min_balance_coverage") <- min_balance_coverage
  out
}

test_that("ws_balance_ilr matrix path is byte-identical to scalar row-loop reference", {
  dat <- .make_ws_balance_ilr_data(n = 12L, seed = 101L)
  X <- dat$X
  X_na <- X
  X_na[2L, ] <- NA_real_
  X_na[3L, c("hsa-miR-451a", "hsa-miR-16-5p")] <- NA_real_

  noise_cols <- grep("^noise", colnames(X), value = TRUE)
  X_no_coverage <- X[, noise_cols, drop = FALSE]

  empty_sbp <- list(
    missing_balance = list(numerator = "missing_num", denominator = "missing_den"),
    present_balance = list(numerator = colnames(X)[1L], denominator = colnames(X)[2L])
  )

  struct_features <- sprintf("hsa-miR-%03d-3p", seq_len(30L))
  struct_annotation <- data.frame(
    feature = struct_features,
    cluster = paste0("clu", rep(1:4, length.out = length(struct_features))),
    gc = seq(0.30, 0.70, length.out = length(struct_features)),
    stringsAsFactors = FALSE
  )
  struct_sbp <- .ss_struct_ilr_build_sbp(
    struct_features, struct_annotation, min_part = 1L)
  set.seed(102L)
  X_struct <- matrix(stats::rgamma(9L * length(struct_features), shape = 3, rate = 1),
                     nrow = 9L,
                     dimnames = list(paste0("T", seq_len(9L)), struct_features))
  X_struct[4L, 5L] <- NA_real_

  cases <- list(
    default_gmean_null = list(
      x = X, balances = ws_default_sbp(), pseudocount = NULL,
      aggregate = "gmean", min_balance_coverage = 0.25
    ),
    default_trimmed_explicit = list(
      x = X, balances = ws_default_sbp(), pseudocount = 0.125,
      aggregate = "trimmed_gmean", min_balance_coverage = 1
    ),
    default_na_and_coverage_fail = list(
      x = X_na, balances = ws_default_sbp(), pseudocount = NULL,
      aggregate = "gmean", min_balance_coverage = 0.8
    ),
    no_coverage_panel = list(
      x = X_no_coverage, balances = ws_default_sbp(), pseudocount = NULL,
      aggregate = "trimmed_gmean", min_balance_coverage = 0.8
    ),
    default_single_row = list(
      x = X[5L, , drop = FALSE], balances = ws_default_sbp(), pseudocount = NULL,
      aggregate = "gmean", min_balance_coverage = 0.25
    ),
    empty_feature_balance = list(
      x = X[, seq_len(6L), drop = FALSE], balances = empty_sbp,
      pseudocount = 0.5, aggregate = "gmean", min_balance_coverage = 0
    ),
    struct_ilr_sbp = list(
      x = X_struct, balances = struct_sbp, pseudocount = NULL,
      aggregate = "trimmed_gmean", min_balance_coverage = 0.5
    )
  )

  for (nm in names(cases)) {
    z <- cases[[nm]]
    actual <- ws_balance_ilr(
      z$x, balances = z$balances, pseudocount = z$pseudocount,
      aggregate = z$aggregate,
      min_balance_coverage = z$min_balance_coverage
    )
    reference <- .ws_balance_ilr_row_loop_reference(
      z$x, balances = z$balances, pseudocount = z$pseudocount,
      aggregate = z$aggregate,
      min_balance_coverage = z$min_balance_coverage
    )
    expect_identical(actual, reference, info = nm)
    expect_equal(actual, reference, tolerance = 0, info = nm)

    for (i in unique(c(1L, min(2L, nrow(z$x)), nrow(z$x)))) {
      single_matrix <- ws_balance_ilr(
        z$x[i, , drop = FALSE], balances = z$balances,
        pseudocount = z$pseudocount, aggregate = z$aggregate,
        min_balance_coverage = z$min_balance_coverage
      )
      single_vector <- ws_balance_ilr(
        z$x[i, ], balances = z$balances, pseudocount = z$pseudocount,
        aggregate = z$aggregate,
        min_balance_coverage = z$min_balance_coverage
      )
      single_values <- single_vector
      attr(single_values, "coverage") <- NULL
      attr(single_values, "coverage_failed") <- NULL
      attr(single_values, "min_balance_coverage") <- NULL

      expect_identical(actual[i, ], single_matrix[1L, ], info = nm)
      expect_identical(actual[i, ], single_values, info = nm)
      expect_identical(attr(actual, "coverage")[[i]],
                       attr(single_vector, "coverage"), info = nm)
      expect_identical(attr(actual, "coverage_failed")[[i]],
                       attr(single_vector, "coverage_failed"), info = nm)
      expect_identical(attr(actual, "min_balance_coverage"),
                       attr(single_vector, "min_balance_coverage"), info = nm)
    }
  }
})

test_that("ws_balance_ilr gmean matrix path is byte-identical at large scale", {
  default_sbp <- ws_default_sbp()
  sentinels <- c("ALL_NON_RBC", "ALL_NON_PLT", "TOP_K_BY_ABUNDANCE",
                 "TAIL_BY_ABUNDANCE")
  default_features <- unique(unlist(lapply(default_sbp, function(b) {
    c(b$numerator, b$denominator)
  }), use.names = FALSE))
  default_features <- default_features[!default_features %in% sentinels]

  struct_features <- paste0("struct_miR_", sprintf("%03d", seq_len(96L)))
  noise_features <- paste0(
    "noise_large_", sprintf(
      "%03d", seq_len(300L - length(default_features) - length(struct_features))
    )
  )
  features <- c(default_features, struct_features, noise_features)
  expect_equal(length(features), 300L)
  expect_equal(length(unique(features)), 300L)

  struct_annotation <- data.frame(
    feature = struct_features,
    cluster = paste0("clu", sprintf("%02d",
                                    rep(seq_len(12L), each = 8L))),
    gc = seq(0.25, 0.75, length.out = length(struct_features)),
    stringsAsFactors = FALSE
  )
  struct_sbp <- .ss_struct_ilr_build_sbp(
    struct_features, struct_annotation, min_part = 2L)
  expect_gt(length(struct_sbp), 0L)

  balances <- c(default_sbp[c("rbc_vs_rest", "platelet_vs_rest")],
                struct_sbp)
  rbc_den_size <- length(setdiff(features, default_sbp$rbc_vs_rest$numerator))
  plt_den_size <- length(setdiff(features, default_sbp$platelet_vs_rest$numerator))
  expect_gt(rbc_den_size, 250L)
  expect_gt(plt_den_size, 250L)

  set.seed(501L)
  X <- matrix(
    stats::rgamma(3000L * length(features), shape = 9, rate = 0.7),
    nrow = 3000L,
    dimnames = list(paste0("LS", seq_len(3000L)), features)
  )

  actual <- ws_balance_ilr(
    X, balances = balances, pseudocount = NULL, aggregate = "gmean",
    min_balance_coverage = 0
  )
  reference <- .ws_balance_ilr_row_loop_reference(
    X, balances = balances, pseudocount = NULL, aggregate = "gmean",
    min_balance_coverage = 0
  )

  expect_identical(actual, reference)
  expect_equal(actual, reference, tolerance = 0)
})

test_that("ws_balance_ilr one-column matrix path preserves feature names", {
  X <- matrix(c(2, 5, 11, 17), ncol = 1L,
              dimnames = list(paste0("One", seq_len(4L)), "only_feature"))
  balances <- list(
    only_feature_self = list(numerator = "only_feature",
                             denominator = "only_feature")
  )

  actual <- ws_balance_ilr(
    X, balances = balances, pseudocount = NULL, aggregate = "gmean",
    min_balance_coverage = 0
  )

  # OLD row-loop behavior errored here: x[i, ] dropped the single column name,
  # so the scalar branch could not match the balance to a feature. The matrix
  # path now reads colnames(x) once and returns the named single-column result.
  rows <- lapply(seq_len(nrow(X)), function(i) {
    ws_balance_ilr(
      setNames(X[i, 1L], colnames(X)),
      balances = balances, pseudocount = NULL, aggregate = "gmean",
      min_balance_coverage = 0
    )
  })
  reference <- do.call(rbind, rows)
  rownames(reference) <- rownames(X)
  attr(reference, "coverage") <- vapply(rows, attr, numeric(1L), "coverage")
  attr(reference, "coverage_failed") <- vapply(rows, function(r) {
    isTRUE(attr(r, "coverage_failed"))
  }, logical(1L))
  attr(reference, "min_balance_coverage") <- 0

  expect_identical(actual, reference)
  expect_true(all(is.finite(actual)))
  expect_equal(as.numeric(actual[, "only_feature_self"]), rep(0, nrow(X)),
               tolerance = 0)
})

test_that("ws_balance_ilr data.frame input matches matrix and scalar reference", {
  X <- matrix(c(
    4, 8, 15, 16,
    23, 42, 5, 9,
    12, 7, 14, 28,
    33, 6, 18, 21,
    10, 20, 30, 40
  ), nrow = 5L, byrow = TRUE,
  dimnames = list(paste0("DF", seq_len(5L)), c("a", "b", "c", "d")))
  X_df <- as.data.frame(X)
  rownames(X_df) <- rownames(X)
  balances <- list(
    a_vs_cd = list(numerator = "a", denominator = c("c", "d")),
    bd_vs_a = list(numerator = c("b", "d"), denominator = "a")
  )

  # OLD path infinite-recursed/stack-overflowed for multi-column data.frames:
  # x[i, ] stayed a data.frame and repeatedly re-entered the matrix/data.frame
  # branch. The current path coerces once and must equal matrix/scalar results.
  actual <- ws_balance_ilr(
    X_df, balances = balances, pseudocount = 0.25, aggregate = "gmean",
    min_balance_coverage = 0
  )
  matrix_reference <- ws_balance_ilr(
    X, balances = balances, pseudocount = 0.25, aggregate = "gmean",
    min_balance_coverage = 0
  )
  scalar_reference <- .ws_balance_ilr_row_loop_reference(
    X, balances = balances, pseudocount = 0.25, aggregate = "gmean",
    min_balance_coverage = 0
  )

  expect_identical(actual, matrix_reference)
  expect_identical(actual, scalar_reference)
  expect_equal(actual, scalar_reference, tolerance = 0)
})

test_that("ws_balance_ilr preserves TOP_K_BY_ABUNDANCE p<=k scalar quirk", {
  X <- matrix(c(
    8, 7, 6, 5, 4, 3, 2, 1,
    1, 3, 5, 7, 2, 4, 6, 8,
    9, 9, 2, 2, 5, 5, 1, 1,
    4, 10, 3, 11, 2, 12, 1, 13
  ), nrow = 4L, byrow = TRUE,
  dimnames = list(paste0("TK", seq_len(4L)), paste0("f", seq_len(8L))))
  balances <- list(
    top_k_all = list(numerator = "TOP_K_BY_ABUNDANCE",
                     denominator = "TAIL_BY_ABUNDANCE",
                     k = 8L)
  )

  # KNOWN pre-existing quirk, intentionally preserved for byte-identity:
  # when p <= k, R's (k + 1):p sequence makes den_idx non-empty instead of
  # representing an empty tail. Fixing that changes frozen within-sample
  # results and belongs in a separate, explicitly results-changing patch.
  actual <- ws_balance_ilr(
    X, balances = balances, pseudocount = NULL, aggregate = "gmean",
    min_balance_coverage = 0
  )
  reference <- .ws_balance_ilr_row_loop_reference(
    X, balances = balances, pseudocount = NULL, aggregate = "gmean",
    min_balance_coverage = 0
  )

  expect_identical(actual, reference)
  expect_equal(actual, reference, tolerance = 0)
  expect_true(all(is.finite(actual[, "top_k_all"])))
})

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
