library(testthat)

.make_reo_ucell_data <- function(n = 96L, p = 80L, seed = 1L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 2, rate = 1), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  up <- features[1:8]
  down <- features[9:16]
  X[y == 1, up] <- X[y == 1, up] + 8
  X[y == 0, down] <- X[y == 0, down] + 8
  list(X = X, y = y, up = up, down = down)
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

test_that("reo-ucell fit/score separates planted up-signal", {
  dat <- .make_reo_ucell_data()
  train <- c(seq_len(36L), 49L:84L)
  test <- setdiff(seq_len(nrow(dat$X)), train)

  model <- fit_reo_ucell(dat$X[train, ], dat$y[train], hp = list(k = 8L))
  score <- score_reo_ucell(model, dat$X[test, ])

  expect_s3_class(model, "reo_ucell_model")
  expect_equal(model$k, 8L)
  expect_true(all(dat$up %in% model$up_features))
  expect_true(any(dat$down %in% model$down_features))
  expect_type(score, "double")
  expect_length(score, length(test))
  expect_true(all(is.finite(score)))
  expect_gt(.auc_or_wilcoxon(dat$y[test], score), 0.7)
})

test_that("reo-ucell passes the canonical row-equivariance gate", {
  dat <- .make_reo_ucell_data(seed = 2L)
  model <- fit_reo_ucell(dat$X[1:72, ], dat$y[1:72], hp = list(k = 8L))
  X_test <- dat$X[73:96, ]
  score_fun <- function(model, X, meta) score_reo_ucell(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("reo-ucell scores are deterministic finite numeric vectors", {
  dat <- .make_reo_ucell_data(seed = 3L)
  model <- fit_reo_ucell(dat$X, dat$y, hp = list(k = 8L))

  s1 <- score_reo_ucell(model, dat$X[1:20, ])
  s2 <- score_reo_ucell(model, dat$X[1:20, ])

  expect_identical(s1, s2)
  expect_type(s1, "double")
  expect_length(s1, 20L)
  expect_true(all(is.finite(s1)))
})

test_that("reo-ucell handles missing signature features without batch coupling", {
  dat <- .make_reo_ucell_data(seed = 4L)
  model <- fit_reo_ucell(dat$X, dat$y, hp = list(k = 8L))
  drop_features <- unique(c(model$up_features[1:3], model$down_features[1:2]))
  X_short <- dat$X[1:12, setdiff(colnames(dat$X), drop_features), drop = FALSE]

  score_short <- score_reo_ucell(model, X_short)
  expect_length(score_short, nrow(X_short))
  expect_true(all(is.finite(score_short)))

  row_id <- 5L
  score_batch <- score_reo_ucell(model, X_short)
  score_single <- score_reo_ucell(model, X_short[row_id, , drop = FALSE])
  expect_equal(score_single, score_batch[row_id], tolerance = 1e-12)
})

test_that("reo-ucell handles single-column / single-shared-feature inputs", {
  # (a) fit + score must not crash on a 1-column matrix, and must not silently
  # return 0 from dropped row names on single-column slicing.
  X1 <- matrix(c(10, 5, 2, 8), nrow = 4,
               dimnames = list(c("S1", "S2", "S3", "S4"), "miR-001"))
  y1 <- c(1, 0, 1, 0)
  m1 <- fit_reo_ucell(X1, y1, hp = list(k = 1L))
  s1 <- score_reo_ucell(m1, X1)
  expect_length(s1, 4L)
  expect_true(all(is.finite(s1)))
  # A single-feature universe ranks each sample's lone feature first -> UCell = 1.
  expect_equal(s1, rep(1, 4L))

  # (b) reachable in single-sample deployment: a specimen sharing only ONE
  # feature with a larger frozen model universe. The shared feature must still
  # drive a finite, non-zero score (the pre-fix names bug returned 0 here), and
  # the single row must equal its position inside the batch.
  dat <- .make_reo_ucell_data(seed = 9L)
  model <- fit_reo_ucell(dat$X, dat$y, hp = list(k = 8L))
  shared <- model$up_features[1]
  X_one <- dat$X[1:5, shared, drop = FALSE]
  s_one <- score_reo_ucell(model, X_one)
  expect_length(s_one, 5L)
  expect_true(all(is.finite(s_one)))
  expect_true(all(s_one > 0))
  expect_equal(score_reo_ucell(model, X_one[3, , drop = FALSE]), s_one[3],
               tolerance = 1e-12)
})

test_that("reo-ucell input validation errors are explicit", {
  dat <- .make_reo_ucell_data(seed = 5L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_reo_ucell(X_unnamed, dat$y), "feature names")

  expect_error(fit_reo_ucell(dat$X, dat$y[-1]), "length\\(y_train\\)")
  expect_error(fit_reo_ucell(dat$X, rep(0, nrow(dat$X))), "at least one case")

  model <- fit_reo_ucell(dat$X, dat$y, hp = list(k = 8L))
  X_no_shared <- dat$X[1:4, , drop = FALSE]
  colnames(X_no_shared) <- paste0("other-", seq_len(ncol(X_no_shared)))
  expect_error(score_reo_ucell(model, X_no_shared), "no shared features")
})
