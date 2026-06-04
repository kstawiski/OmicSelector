library(testthat)

# Synthetic generator with a CASE-specific MEAN contrast AND a CASE-only feature-
# covariance block -- the design domain of a Gaussian LRT. Both signals are built
# as mean-zero contrasts over the first `k` block features so they survive the
# within-sample rCLR centring (a uniform shift would be removed by it). The mean
# contrast is a fixed alternating-sign offset added to every case (a deterministic
# class mean difference in rCLR space); the covariance block is a shared latent
# factor loading with alternating sign onto the same features for cases only (a
# class-specific dependence). Abundances are exp() of the log-scale construction
# so the input is strictly positive (compositional).
.make_lrt_bw_data <- function(n = 160L, p = 30L, k = 12L, mean_shift = 0.9,
                              load = 1.4, sd = 0.6, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case <- which(y == 1L)
  shift <- rep(c(mean_shift, -mean_shift), length.out = k)   # mean-zero contrast
  L[case, seq_len(k)] <- L[case, seq_len(k)] +
    matrix(shift, nrow = length(case), ncol = k, byrow = TRUE)
  f <- stats::rnorm(length(case))                            # shared latent factor
  loadings <- rep(c(load, -load), length.out = k)            # mean-zero loadings
  L[case, seq_len(k)] <- L[case, seq_len(k)] + outer(f, loadings)
  X <- exp(L)
  list(X = X, y = y, block = features[seq_len(k)])
}

.auc_or_wilcoxon_bw <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Standalone per-sample rCLR (matches the scorer's frozen transform) for the
# manual Gaussian-LLR cross-check.
.rclr0_bw <- function(v) {
  z <- rep.int(0, length(v))
  pos <- which(v > 0)
  if (length(pos) >= 1L) {
    lv <- log(v[pos])
    z[pos] <- lv - mean(lv)
  }
  z
}

# Random SPD matrix for the manual Bures-Wasserstein geodesic checks.
.spd_bw <- function(p, seed) {
  set.seed(seed)
  M <- matrix(stats::rnorm(p * p), p, p)
  crossprod(M) + diag(p)              # PD: A^T A + I
}

test_that("lrt-bw fit/score roundtrip has the right shape and types", {
  dat <- .make_lrt_bw_data()
  model <- fit_lrt_bw(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 40L, seed = 7L))
  expect_s3_class(model, "lrt_bw_model")
  expect_equal(nrow(model$case_anchors), 40L)
  expect_equal(nrow(model$control_anchors), 40L)
  expect_identical(model$feature_universe, colnames(dat$X))
  p <- ncol(dat$X)
  # Frozen full-universe precisions are symmetric, the right size and finite.
  expect_equal(dim(model$repr$Rinv_case), c(p, p))
  expect_equal(dim(model$repr$Rinv_control), c(p, p))
  expect_equal(model$repr$Rinv_case, t(model$repr$Rinv_case), tolerance = 1e-10)
  expect_length(model$repr$mu_case, p)
  expect_true(is.finite(model$repr$logdet_case))

  score <- score_lrt_bw(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("lrt-bw separates a planted class mean+covariance difference", {
  dat <- .make_lrt_bw_data(n = 300L, seed = 202L)
  ctrl <- which(dat$y == 0L)
  case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_lrt_bw(dat$X[train, ], dat$y[train],
                      hp = list(max_anchors_per_class = 112L, seed = 11L))
  score <- score_lrt_bw(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_bw(dat$y[test], score), 0.7)
})

test_that("lrt-bw passes the canonical row-equivariance gate", {
  dat <- .make_lrt_bw_data(seed = 43L)
  model <- fit_lrt_bw(dat$X[1:120, ], dat$y[1:120],
                      hp = list(max_anchors_per_class = 30L, seed = 13L))
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_lrt_bw(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("lrt-bw scoring is exactly invariant to per-sample scaling", {
  dat <- .make_lrt_bw_data(seed = 44L)
  model <- fit_lrt_bw(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]
  expect_equal(score_lrt_bw(model, row * 7),
               score_lrt_bw(model, row), tolerance = 1e-8)
  expect_equal(score_lrt_bw(model, row * 1e6),
               score_lrt_bw(model, row), tolerance = 1e-8)
  # Row-equivariant batch: scaling some rows leaves the others' scores fixed.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e6
  expect_equal(score_lrt_bw(model, scaled),
               score_lrt_bw(model, batch), tolerance = 1e-8)
  # A row with mostly zeros (rCLR over its own nonzero parts) still scores finite.
  X_sparse <- dat$X[1:3, , drop = FALSE]
  X_sparse[1, 3:ncol(X_sparse)] <- 0
  expect_true(all(is.finite(score_lrt_bw(model, X_sparse))))
})

test_that("lrt-bw handles partial feature overlap (consistent + single-row)", {
  dat <- .make_lrt_bw_data(seed = 45L)
  model <- fit_lrt_bw(dat$X, dat$y)
  # Subset and REORDER the columns: present = intersect(universe, colnames) keeps
  # feature-universe order, so the score is invariant to the X column order.
  keep <- rev(c(dat$block[1:5], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_lrt_bw(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_lrt_bw(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_lrt_bw(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("lrt-bw matches the canonical dispatch", {
  dat <- .make_lrt_bw_data(seed = 46L)
  model <- fit_lrt_bw(dat$X[1:120, ], dat$y[1:120],
                      hp = list(max_anchors_per_class = 30L, seed = 17L))
  X_test <- dat$X[121:160, ]
  direct <- score_lrt_bw(model, X_test)

  # Post-integration the manifest carries the lrt-bw row and score_lrt_bw lives in
  # the OmicSelector namespace, so the canonical adapter resolves it directly (the
  # real roster -> package dispatch path). In the pre-integration staged self-test
  # score_lrt_bw is not yet in the namespace, so we clone the roster row if absent
  # and register a thin canonical adapter to exercise the same dispatch wiring.
  # Either way the dispatched score must equal the direct call.
  roster <- singlesample_method_roster()
  if (!"lrt-bw" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "lrt-bw"
    tmpl$fit_fn <- "fit_lrt_bw"
    tmpl$score_fn <- "score_lrt_bw"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_lrt_bw", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "lrt-bw",
      function(model, X, meta) score_lrt_bw(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("lrt-bw", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("lrt-bw returns the neutral 0 below the feature-overlap floor", {
  dat <- .make_lrt_bw_data(seed = 47L)
  model <- fit_lrt_bw(dat$X, dat$y, hp = list(min_features = 3L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, colnames(dat$X)[1:2], drop = FALSE]
  expect_equal(score_lrt_bw(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_lrt_bw(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
})

test_that("lrt-bw fitting/scoring is deterministic and RNG-safe", {
  dat <- .make_lrt_bw_data(seed = 48L)
  hp <- list(max_anchors_per_class = 25L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_lrt_bw(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_lrt_bw(dat$X, dat$y, hp = hp)
  expect_identical(model_1, model_2)

  s1 <- score_lrt_bw(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_lrt_bw(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("lrt-bw fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_lrt_bw_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  # max_anchors_per_class > class size -> no subsampling -> RNG never touched.
  model <- fit_lrt_bw(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 200L, seed = 29L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "lrt_bw_model")
})

test_that("lrt-bw leaves no RNG state when subsampling fits in a fresh session", {
  # Companion to the no-subsampling case: a class EXCEEDS the anchor cap, so
  # set.seed() IS called inside the fit. With no prior global seed, the on.exit
  # restore must rm() the .Random.seed that set.seed() created.
  dat <- .make_lrt_bw_data(seed = 53L)            # this set.seed() seeds the RNG
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())       # wipe to a genuine no-seed state
  }
  model <- fit_lrt_bw(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_lte(nrow(model$case_anchors), 20L)        # subsampling actually happened
  expect_lte(nrow(model$control_anchors), 20L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # And subsampling is deterministic across fits (frozen seed).
  model2 <- fit_lrt_bw(dat$X, dat$y,
                       hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_identical(rownames(model$case_anchors), rownames(model2$case_anchors))
})

test_that("lrt-bw stays finite and positive definite in the p > n regime", {
  # p = 30 features, 12 anchors per class -> rank-deficient class covariance.
  # The eigenvalue floor + BW-geodesic shrink + ridge-to-PD loop must keep the
  # fit and every score finite, and the precisions positive definite.
  dat <- .make_lrt_bw_data(n = 24L, p = 30L, k = 8L, seed = 71L)
  model <- fit_lrt_bw(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 12L, seed = 5L))
  expect_equal(nrow(model$case_anchors), 12L)
  expect_true(all(is.finite(model$repr$Rinv_case)))
  expect_true(all(is.finite(model$repr$Rinv_control)))
  expect_true(is.finite(model$repr$logdet_case))
  expect_true(is.finite(model$repr$logdet_control))
  # Precision PD <=> covariance PD.
  ev_case <- eigen(model$repr$Rinv_case, symmetric = TRUE,
                   only.values = TRUE)$values
  ev_control <- eigen(model$repr$Rinv_control, symmetric = TRUE,
                      only.values = TRUE)$values
  expect_true(all(ev_case > 0))
  expect_true(all(ev_control > 0))
  score <- score_lrt_bw(model, dat$X)
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
  # The pooled anchor variant must also be finite/PD in the p > n regime.
  model_p <- fit_lrt_bw(dat$X, dat$y,
                        hp = list(max_anchors_per_class = 12L, seed = 5L,
                                  anchor = "pooled"))
  expect_true(all(is.finite(score_lrt_bw(model_p, dat$X))))
})

test_that("lrt-bw handles a constant/degenerate training feature without NA", {
  # A feature zero across ALL anchors of a class has a constant (all-zero) rCLR
  # column -> zero-variance covariance row/col -> a zero eigenvalue. The make_pd
  # eigenvalue floor and the ridge-to-PD loop must keep the fit and every score
  # finite.
  dat <- .make_lrt_bw_data(seed = 54L)
  X <- dat$X
  X[, ncol(X)] <- 0                       # feature constant (zero) in BOTH classes
  X[dat$y == 1L, ncol(X) - 1L] <- 0       # feature constant (zero) in cases only
  model <- fit_lrt_bw(X, dat$y, hp = list(max_anchors_per_class = 40L))
  expect_true(all(is.finite(model$repr$Rinv_case)))
  expect_true(all(is.finite(model$repr$Rinv_control)))
  expect_true(is.finite(model$repr$logdet_case))
  expect_true(is.finite(model$repr$logdet_control))
  score <- score_lrt_bw(model, X)
  expect_length(score, nrow(X))
  expect_true(all(is.finite(score)))
})

test_that("lrt-bw hyperparameter validation is strict", {
  dat <- .make_lrt_bw_data(seed = 50L)
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(shrink = 1)),
               "hp\\$shrink")
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(shrink = -0.1)),
               "hp\\$shrink")
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(anchor = "euclidean")),
               "hp\\$anchor")
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
               "max_anchors_per_class")
  # Out-of-integer-range must hit the explicit error, not an as.integer overflow.
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
               "max_anchors_per_class")
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  # `$` partial matching would let "min" bind to "min_features"; the exact [[ ]]
  # reads + allowed-list reject the unknown field instead.
  expect_error(fit_lrt_bw(dat$X, dat$y, hp = list(min = 3L)), "unknown hp")
})

test_that("lrt-bw BW geodesic satisfies its defining identities", {
  eps <- 1e-6
  A <- .spd_bw(4L, seed = 101L)
  B <- .spd_bw(4L, seed = 202L)

  # Endpoints.
  expect_equal(.lrt_bw_geo(A, B, 0, eps), A, tolerance = 1e-9)
  expect_equal(.lrt_bw_geo(A, B, 1, eps), B, tolerance = 1e-9)
  # Constant geodesic: BWgeo(A, A, t) == A through the full path (t in (0,1)).
  expect_equal(.lrt_bw_geo(A, A, 0.5, eps), A, tolerance = 1e-7)
  expect_equal(.lrt_bw_geo(A, A, 0.3, eps), A, tolerance = 1e-7)
  # Equal-weight barycenter is symmetric in its arguments (a non-trivial property
  # for non-commuting A, B): BWgeo(A, B, 0.5) == BWgeo(B, A, 0.5).
  expect_equal(.lrt_bw_geo(A, B, 0.5, eps), .lrt_bw_geo(B, A, 0.5, eps),
               tolerance = 1e-7)
  # The midpoint is symmetric and positive definite.
  mid <- .lrt_bw_geo(A, B, 0.5, eps)
  expect_equal(mid, t(mid), tolerance = 1e-10)
  expect_true(all(eigen(mid, symmetric = TRUE, only.values = TRUE)$values > 0))

  # 2x2 diagonal closed form: for commuting (diagonal) A, B the BW geodesic is
  # BWgeo(A, B, t)_ii = ((1-t) sqrt(a_i) + t sqrt(b_i))^2. With A = diag(4, 1),
  # B = diag(9, 16), t = 0.5 the midpoint is diag(6.25, 6.25).
  Ad <- diag(c(4, 1))
  Bd <- diag(c(9, 16))
  expect_equal(.lrt_bw_geo(Ad, Bd, 0.5, eps), diag(c(6.25, 6.25)),
               tolerance = 1e-9)
})

test_that("lrt-bw score equals the manual Gaussian log-likelihood ratio", {
  dat <- .make_lrt_bw_data(seed = 61L)
  model <- fit_lrt_bw(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 40L, seed = 9L))
  repr <- model$repr     # full-universe frozen mean/precision/logdet
  # Scoring over the full universe re-derives exactly this representation, so the
  # score must equal the hand-computed Gaussian LLR.
  for (i in c(1L, 5L, 12L, 80L, 150L)) {
    z <- .rclr0_bw(dat$X[i, ])
    dc <- z - repr$mu_case
    d0 <- z - repr$mu_control
    l_case <- -0.5 * as.numeric(t(dc) %*% repr$Rinv_case %*% dc) -
      0.5 * repr$logdet_case
    l_control <- -0.5 * as.numeric(t(d0) %*% repr$Rinv_control %*% d0) -
      0.5 * repr$logdet_control
    manual <- l_case - l_control
    got <- score_lrt_bw(model, dat$X[i, , drop = FALSE])
    expect_equal(got, manual, tolerance = 1e-9)
  }
})
