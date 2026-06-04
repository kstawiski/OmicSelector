library(testthat)

.make_reo_metaktsp_data <- function(n_per_cohort = 40L) {
  cohorts <- paste0("C", seq_len(3L))
  y <- unlist(lapply(cohorts, function(cohort) {
    rep(c(0, 1), each = n_per_cohort / 2L)
  }))
  meta <- data.frame(accession = rep(cohorts, each = n_per_cohort),
                     stringsAsFactors = FALSE)
  features <- c("ConsA", "ConsB", "DistA", "DistB", "NoiseA", "NoiseB")
  X <- matrix(5, nrow = length(y), ncol = length(features),
              dimnames = list(paste0("S", seq_along(y)), features))

  X[y == 0L, "ConsA"] <- 12
  X[y == 0L, "ConsB"] <- 1
  X[y == 1L, "ConsA"] <- 1
  X[y == 1L, "ConsB"] <- 12
  X[, "NoiseA"] <- 4
  X[, "NoiseB"] <- 6

  for (cohort in cohorts) {
    rows <- meta$accession == cohort
    if (cohort == "C1") {
      X[rows & y == 0L, "DistA"] <- 12
      X[rows & y == 0L, "DistB"] <- 1
      X[rows & y == 1L, "DistA"] <- 1
      X[rows & y == 1L, "DistB"] <- 12
    } else {
      X[rows, "DistA"] <- 1
      X[rows, "DistB"] <- 12
    }
  }

  list(X = X, y = y, meta = meta)
}

.make_reo_metaktsp_skip_data <- function() {
  y <- c(rep(0, 10L), rep(1, 10L), rep(1, 10L), rep(0, 10L))
  meta <- data.frame(accession = c(rep("usable", 20L),
                                   rep("case_only", 10L),
                                   rep(NA_character_, 10L)),
                     stringsAsFactors = FALSE)
  features <- c("KeepA", "KeepB", "SkipA", "SkipB")
  X <- matrix(5, nrow = length(y), ncol = length(features),
              dimnames = list(paste0("K", seq_along(y)), features))

  usable <- meta$accession == "usable" & !is.na(meta$accession)
  X[usable & y == 0L, "KeepA"] <- 12
  X[usable & y == 0L, "KeepB"] <- 1
  X[usable & y == 1L, "KeepA"] <- 1
  X[usable & y == 1L, "KeepB"] <- 12

  X[meta$accession == "case_only" & !is.na(meta$accession), "SkipA"] <- 1
  X[meta$accession == "case_only" & !is.na(meta$accession), "SkipB"] <- 20
  X[is.na(meta$accession), "SkipA"] <- 20
  X[is.na(meta$accession), "SkipB"] <- 1

  list(X = X, y = y, meta = meta)
}

.auc_or_wilcoxon_metaktsp <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.reo_metaktsp_pair_rank <- function(pairs, feature_a, feature_b) {
  which((pairs$feature_a == feature_a & pairs$feature_b == feature_b) |
          (pairs$feature_a == feature_b & pairs$feature_b == feature_a))
}

test_that("reo-metaktsp prioritizes planted cross-cohort pair-order signal", {
  dat <- .make_reo_metaktsp_data()
  model <- fit_reo_metaktsp(dat$X, dat$y, dat$meta, hp = list(k = 15L))
  score <- score_reo_metaktsp(model, dat$X)

  consistent_rank <- .reo_metaktsp_pair_rank(model$pairs, "ConsA", "ConsB")
  distractor_rank <- .reo_metaktsp_pair_rank(model$pairs, "DistA", "DistB")

  expect_s3_class(model, "reo_metaktsp_model")
  expect_equal(model$hp$k, 15L)
  expect_equal(model$n_cohorts, 3L)
  expect_length(consistent_rank, 1L)
  expect_length(distractor_rank, 1L)
  expect_lt(consistent_rank, distractor_rank)
  expect_gt(model$pairs$score[consistent_rank], model$pairs$score[distractor_rank])
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
  expect_true(all(score >= 0 & score <= 1))
  expect_gt(mean(score[dat$y == 1L]), mean(score[dat$y == 0L]))
  expect_gt(.auc_or_wilcoxon_metaktsp(dat$y, score), 0.7)
})

test_that("reo-metaktsp equal-cohort weighting differs from pooled k-TSP", {
  # The defining property: avgTSP weights each cohort equally, so it is NOT the
  # pooled k-TSP. A balanced fixture cannot show this (equal-weight mean ==
  # pooled delta when cohorts are equal-sized and equally balanced). Here the
  # P<Q pair order is reversed in one LARGE cohort but consistent in two SMALL
  # cohorts. Equal-weight avgTSP = mean(-1, +1, +1) = +1/3 orients P<Q case-high,
  # whereas pooled k-TSP is dominated by the large cohort and orients Q<P. This
  # locks the frozen bank against a future silent reversion to pooling.
  acc <- c(rep("L", 80L), rep("S1", 8L), rep("S2", 8L))
  y <- c(rep(c(0, 1), each = 40L), rep(c(0, 1), each = 4L),
         rep(c(0, 1), each = 4L))
  P <- Q <- numeric(length(y))
  iL <- acc == "L"
  iS <- acc %in% c("S1", "S2")
  P[iL & y == 0] <- 1;  Q[iL & y == 0] <- 12
  P[iL & y == 1] <- 12; Q[iL & y == 1] <- 1
  P[iS & y == 0] <- 12; Q[iS & y == 0] <- 1
  P[iS & y == 1] <- 1;  Q[iS & y == 1] <- 12
  X <- cbind(P = P, Q = Q)
  rownames(X) <- paste0("s", seq_along(y))

  meta_model <- fit_reo_metaktsp(
    X, y, data.frame(accession = acc, stringsAsFactors = FALSE),
    hp = list(k = 1L)
  )
  pooled <- fit_reo_ktsp(X, y, hp = list(k = 1L))$ktsp$pairs

  expect_equal(meta_model$n_cohorts, 3L)
  expect_equal(meta_model$pairs$feature_a[1L], "P")  # equal-weight meta: +1/3
  expect_equal(pooled$feature_a[1L], "Q")            # pooled: large cohort wins
})

test_that("reo-metaktsp degenerates to plain k-TSP when meta_train is NULL", {
  dat <- .make_reo_metaktsp_data()
  model <- fit_reo_metaktsp(dat$X, dat$y, meta_train = NULL, hp = list(k = 5L))
  ktsp <- fit_reo_ktsp(dat$X, dat$y, hp = list(k = 5L))

  expect_equal(model$pairs, ktsp$ktsp$pairs)
  expect_equal(model$k, ktsp$k)
  expect_equal(model$n_cohorts, 1L)
})

test_that("reo-metaktsp degenerates to plain k-TSP for one training cohort", {
  dat <- .make_reo_metaktsp_data()
  one_meta <- data.frame(accession = rep("one", nrow(dat$X)),
                         stringsAsFactors = FALSE)
  model <- fit_reo_metaktsp(dat$X, dat$y, one_meta, hp = list(k = 5L))
  ktsp <- fit_reo_ktsp(dat$X, dat$y, hp = list(k = 5L))

  expect_equal(model$pairs, ktsp$ktsp$pairs)
  expect_equal(model$k, ktsp$k)
  expect_equal(model$n_cohorts, 1L)
})

test_that("reo-metaktsp passes the canonical row-equivariance gate", {
  dat <- .make_reo_metaktsp_data()
  model <- fit_reo_metaktsp(dat$X, dat$y, dat$meta, hp = list(k = 6L))
  X_test <- dat$X[21:80, , drop = FALSE]
  score_fun <- function(model, X, meta) score_reo_metaktsp(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("reo-metaktsp is deployable through the canonical dispatch", {
  dat <- .make_reo_metaktsp_data()
  model <- fit_reo_metaktsp(dat$X, dat$y, dat$meta, hp = list(k = 6L))
  X_test <- dat$X[41:100, , drop = FALSE]

  expect_equal(
    singlesample_score_call("reo-metaktsp", model, X_test),
    score_reo_metaktsp(model, X_test),
    tolerance = 1e-12
  )
})

test_that("reo-metaktsp fitting is deterministic and does not create RNG state", {
  dat <- .make_reo_metaktsp_data()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", seed_before, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }

  model_1 <- fit_reo_metaktsp(dat$X, dat$y, dat$meta, hp = list(k = 6L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  model_2 <- fit_reo_metaktsp(dat$X, dat$y, dat$meta, hp = list(k = 6L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  expect_identical(model_1, model_2)
})

test_that("reo-metaktsp handles partial overlap and neutral no-pair cases", {
  dat <- .make_reo_metaktsp_data()
  model <- fit_reo_metaktsp(dat$X, dat$y, dat$meta, hp = list(k = 6L))
  p1 <- model$pairs[1, ]
  keep <- unique(c(p1$feature_a, p1$feature_b))
  X_one_pair <- dat$X[1:8, keep, drop = FALSE]

  score_one_pair <- score_reo_metaktsp(model, X_one_pair)
  expect_length(score_one_pair, nrow(X_one_pair))
  expect_true(all(is.finite(score_one_pair)))
  expect_equal(
    score_reo_metaktsp(model, X_one_pair[5, , drop = FALSE]),
    score_one_pair[5],
    tolerance = 1e-12
  )

  X_no_shared <- matrix(1, nrow = 4L, ncol = 2L,
                        dimnames = list(paste0("N", seq_len(4L)),
                                        c("other-a", "other-b")))
  expect_equal(score_reo_metaktsp(model, X_no_shared),
               rep(0.5, nrow(X_no_shared)))
})

test_that("reo-metaktsp drops NA cohorts, skips single-class cohorts, and falls back", {
  skip_dat <- .make_reo_metaktsp_skip_data()
  model <- fit_reo_metaktsp(skip_dat$X, skip_dat$y, skip_dat$meta,
                            hp = list(k = 6L))
  keep_rank <- .reo_metaktsp_pair_rank(model$pairs, "KeepA", "KeepB")

  expect_equal(model$n_cohorts, 1L)
  expect_length(keep_rank, 1L)
  expect_equal(keep_rank, 1L)
  expect_true(all(is.finite(score_reo_metaktsp(model, skip_dat$X))))

  dat <- .make_reo_metaktsp_data()
  zero_meta <- data.frame(accession = ifelse(dat$y == 1L, "case_only",
                                             "control_only"),
                          stringsAsFactors = FALSE)
  fallback <- fit_reo_metaktsp(dat$X, dat$y, zero_meta, hp = list(k = 5L))
  ktsp <- fit_reo_ktsp(dat$X, dat$y, hp = list(k = 5L))

  expect_equal(fallback$pairs, ktsp$ktsp$pairs)
  expect_equal(fallback$n_cohorts, 1L)
})

test_that("reo-metaktsp input validation errors are explicit", {
  dat <- .make_reo_metaktsp_data()

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_reo_metaktsp(X_unnamed, dat$y, dat$meta), "feature names")

  expect_error(fit_reo_metaktsp(dat$X[, 1, drop = FALSE], dat$y, dat$meta),
               "at least two features")
  expect_error(fit_reo_metaktsp(dat$X, dat$y[-1], dat$meta),
               "length\\(y_train\\)")
  expect_error(fit_reo_metaktsp(dat$X, dat$y, dat$meta[-1, , drop = FALSE]),
               "one row per row")
  expect_error(fit_reo_metaktsp(dat$X, rep(0, nrow(dat$X)), dat$meta),
               "at least one case")
  expect_error(fit_reo_metaktsp(dat$X, rep(1, nrow(dat$X)), dat$meta),
               "at least one control")
  expect_error(fit_reo_metaktsp(dat$X, dat$y, dat$meta,
                                hp = list(maxRank = 10L)),
               "unknown hp")
  expect_error(fit_reo_metaktsp(dat$X, dat$y, dat$meta,
                                hp = list(k = 1.5)),
               "positive integer")
  expect_error(fit_reo_metaktsp(dat$X, dat$y, dat$meta,
                                hp = list(k = 0L)),
               "positive integer")
  expect_error(fit_reo_metaktsp(dat$X, dat$y, dat$meta,
                                hp = list(cohort_col = 1)),
               "cohort_col")
  expect_error(fit_reo_metaktsp(dat$X, dat$y, dat$meta,
                                hp = list(cohort_col = character(0))),
               "cohort_col")
  expect_error(fit_reo_metaktsp(dat$X, dat$y, dat$meta,
                                hp = list(cohort_col = "")),
               "cohort_col")
})
