library(testthat)

# Synthetic generator with a CASE-ONLY DENSITY (location) shift in a feature
# block. The first `k` of `p` log-abundance features are elevated by `shift` for
# case samples, so after exp() and the within-sample rCLR centring those features
# sit higher RELATIVE to the per-sample geometric mean for cases than for
# controls. The case and control rCLR clouds therefore occupy different regions of
# feature space -- exactly the marginal/location density difference that a direct
# density-ratio method (uLSIF) is built to detect (unlike the copula methods,
# which match marginals out and read only dependence). Abundances are exp() so the
# input is strictly positive (compositional).
.make_dre_ulsif_data <- function(n = 160L, p = 30L, k = 10L, shift = 1.3,
                                 sd = 0.5, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case <- which(y == 1L)
  L[case, seq_len(k)] <- L[case, seq_len(k)] + shift  # case-only density shift
  X <- exp(L)
  list(X = X, y = y, block = features[seq_len(k)])
}

.auc_or_wilcoxon_ulsif <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Independent (loop-based) Gaussian kernel matrix, used only by the manual-math
# block to cross-check the vectorized package helper .dre_ulsif_kernel.
.manual_ulsif_kernel <- function(Z, C, sigma) {
  K <- matrix(0, nrow(Z), nrow(C))
  for (i in seq_len(nrow(Z))) {
    for (l in seq_len(nrow(C))) {
      K[i, l] <- exp(-sum((Z[i, ] - C[l, ])^2) / (2 * sigma^2))
    }
  }
  K
}

test_that("dre-ulsif fit/score roundtrip has the right shape and types", {
  dat <- .make_dre_ulsif_data()
  model <- fit_dre_ulsif(dat$X, dat$y,
                         hp = list(max_anchors_per_class = 40L,
                                   n_centers = 30L, seed = 7L))
  expect_s3_class(model, "dre_ulsif_model")
  expect_equal(nrow(model$case_anchors), 40L)
  expect_equal(nrow(model$control_anchors), 40L)
  expect_equal(nrow(model$center_anchors), 30L)
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_equal(model$hp$lambda, 1e-3)              # default ridge
  # Frozen full-universe representation: centres in rCLR space (b x p), one
  # bandwidth, one coefficient per centre, all non-negative after the clip.
  p <- ncol(dat$X)
  expect_equal(dim(model$repr$centers), c(30L, p))
  expect_length(model$repr$alpha, 30L)
  expect_true(is.finite(model$repr$sigma) && model$repr$sigma > 0)
  expect_true(all(model$repr$alpha >= 0))

  score <- score_dre_ulsif(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("dre-ulsif separates a planted case-vs-control density shift", {
  # Larger cohort + held-out split so the AUC is a stable estimate of the genuine
  # density-ratio signal, not a small-sample fluctuation.
  dat <- .make_dre_ulsif_data(n = 300L, seed = 202L)
  ctrl <- which(dat$y == 0L)
  case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_dre_ulsif(dat$X[train, ], dat$y[train],
                         hp = list(max_anchors_per_class = 112L, seed = 11L))
  score <- score_dre_ulsif(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_ulsif(dat$y[test], score), 0.7)
})

test_that("dre-ulsif passes the canonical row-equivariance gate", {
  dat <- .make_dre_ulsif_data(seed = 43L)
  model <- fit_dre_ulsif(dat$X[1:120, ], dat$y[1:120],
                         hp = list(max_anchors_per_class = 30L,
                                   n_centers = 20L, seed = 13L))
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_dre_ulsif(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("dre-ulsif scoring is exactly invariant to per-sample scaling", {
  dat <- .make_dre_ulsif_data(seed = 44L)
  model <- fit_dre_ulsif(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]
  expect_equal(score_dre_ulsif(model, row * 7),
               score_dre_ulsif(model, row), tolerance = 1e-8)
  expect_equal(score_dre_ulsif(model, row * 1e6),
               score_dre_ulsif(model, row), tolerance = 1e-8)
  # Row-equivariant batch: scaling some rows leaves the others' scores fixed.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e6
  expect_equal(score_dre_ulsif(model, scaled),
               score_dre_ulsif(model, batch), tolerance = 1e-8)
  # A row with mostly zeros (rCLR over its own nonzero parts) still scores finite.
  X_sparse <- dat$X[1:3, , drop = FALSE]
  X_sparse[1, 3:ncol(X_sparse)] <- 0
  expect_true(all(is.finite(score_dre_ulsif(model, X_sparse))))
})

test_that("dre-ulsif handles partial feature overlap (consistent + single-row)", {
  dat <- .make_dre_ulsif_data(seed = 45L)
  model <- fit_dre_ulsif(dat$X, dat$y)
  # Subset and REORDER the columns: present = intersect(universe, colnames) keeps
  # feature-universe order, so the score is invariant to the X column order.
  keep <- rev(c(dat$block[1:5], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_dre_ulsif(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_dre_ulsif(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_dre_ulsif(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("dre-ulsif matches the canonical dispatch", {
  dat <- .make_dre_ulsif_data(seed = 46L)
  model <- fit_dre_ulsif(dat$X[1:120, ], dat$y[1:120],
                         hp = list(max_anchors_per_class = 30L,
                                   n_centers = 20L, seed = 17L))
  X_test <- dat$X[121:160, ]
  direct <- score_dre_ulsif(model, X_test)

  # Post-integration the manifest carries the dre-ulsif row and score_dre_ulsif
  # lives in the OmicSelector namespace, so the canonical adapter resolves it
  # directly (the real roster -> package dispatch path). In the pre-integration
  # staged self-test the namespace entry does not yet exist (and the manifest row
  # may or may not yet be present), so we clone a one-row roster from a committed
  # within-discriminator row when needed and register a thin adapter to exercise
  # the same dispatch wiring. Either way the dispatched score must equal the
  # direct call.
  roster <- singlesample_method_roster()
  if (!"dre-ulsif" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "dre-ulsif"
    tmpl$fit_fn <- "fit_dre_ulsif"
    tmpl$score_fn <- "score_dre_ulsif"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_dre_ulsif", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "dre-ulsif",
      function(model, X, meta) score_dre_ulsif(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("dre-ulsif", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("dre-ulsif returns the neutral 0 below the feature-overlap floor", {
  dat <- .make_dre_ulsif_data(seed = 47L)
  model <- fit_dre_ulsif(dat$X, dat$y, hp = list(min_features = 3L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, colnames(dat$X)[1:2], drop = FALSE]
  expect_equal(score_dre_ulsif(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_dre_ulsif(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
})

test_that("dre-ulsif fitting/scoring is deterministic and RNG-safe", {
  dat <- .make_dre_ulsif_data(seed = 48L)
  hp <- list(max_anchors_per_class = 25L, n_centers = 15L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_dre_ulsif(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_dre_ulsif(dat$X, dat$y, hp = hp)
  expect_identical(model_1, model_2)

  s1 <- score_dre_ulsif(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_dre_ulsif(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  # Scoring never touches the global RNG.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("dre-ulsif fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_dre_ulsif_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  # max_anchors_per_class > class size (80) and n_centers >= case-anchor count, so
  # neither anchors nor centres are subsampled -> the RNG is never touched.
  model <- fit_dre_ulsif(dat$X, dat$y,
                         hp = list(max_anchors_per_class = 200L,
                                   n_centers = 100L, seed = 29L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "dre_ulsif_model")
  expect_equal(nrow(model$center_anchors), nrow(model$case_anchors))  # all centres
})

test_that("dre-ulsif leaves no RNG state when subsampling fits in a fresh session", {
  # Companion to the no-subsampling case above, exercising the OTHER branch of the
  # anchor/centre subsampling RNG save/restore: here both the anchor cap and the
  # centre cap bite, so set.seed() IS called inside the fit. With no prior global
  # seed, the on.exit restore must rm() the .Random.seed that set.seed() created --
  # not leave it behind. (80 anchors/class > cap 20, and 20 case anchors > 10
  # centres -> both subsampling paths taken.)
  dat <- .make_dre_ulsif_data(seed = 53L)            # this set.seed() seeds the RNG
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())          # wipe to a genuine no-seed state
  }
  model <- fit_dre_ulsif(dat$X, dat$y,
                         hp = list(max_anchors_per_class = 20L,
                                   n_centers = 10L, seed = 31L))
  expect_lte(nrow(model$case_anchors), 20L)           # anchor subsampling happened
  expect_lte(nrow(model$control_anchors), 20L)
  expect_equal(nrow(model$center_anchors), 10L)       # centre subsampling happened
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # And subsampling is deterministic across fits (frozen seed).
  model2 <- fit_dre_ulsif(dat$X, dat$y,
                          hp = list(max_anchors_per_class = 20L,
                                    n_centers = 10L, seed = 31L))
  expect_identical(rownames(model$case_anchors), rownames(model2$case_anchors))
  expect_identical(rownames(model$center_anchors),
                   rownames(model2$center_anchors))
})

test_that("dre-ulsif handles a constant/degenerate training feature without NA", {
  # A feature that is zero across ALL anchors of a class has a constant (all-zero)
  # rCLR column, which only shifts kernel arguments; the Gram matrix H stays PSD
  # and the ridge keeps H + lambda I invertible, so the fit and every score must
  # remain finite (the median-heuristic bandwidth also stays positive because the
  # centres still differ on the other features).
  dat <- .make_dre_ulsif_data(seed = 54L)
  X <- dat$X
  X[, ncol(X)] <- 0                       # feature constant (zero) in BOTH classes
  X[dat$y == 1L, ncol(X) - 1L] <- 0       # feature constant (zero) in cases only
  model <- fit_dre_ulsif(X, dat$y, hp = list(max_anchors_per_class = 40L))
  expect_true(is.finite(model$repr$sigma) && model$repr$sigma > 0)
  expect_true(all(is.finite(model$repr$alpha)))
  expect_true(all(model$repr$alpha >= 0))
  score <- score_dre_ulsif(model, X)
  expect_length(score, nrow(X))
  expect_true(all(is.finite(score)))
  # Scoring a specimen that is itself zero on the degenerate features is finite.
  expect_true(all(is.finite(score_dre_ulsif(model, X[1:4, , drop = FALSE]))))
})

test_that("dre-ulsif hyperparameter validation is strict", {
  dat <- .make_dre_ulsif_data(seed = 50L)
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(lambda = 0)),
               "hp\\$lambda")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(lambda = -1)),
               "hp\\$lambda")
  # sigma must be NULL or a positive finite number.
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(sigma = 0)),
               "hp\\$sigma")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(sigma = -2)),
               "hp\\$sigma")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(sigma = c(1, 2))),
               "hp\\$sigma")
  # n_centers must be a single positive integer.
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(n_centers = 1.5)),
               "n_centers")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(n_centers = 0)),
               "n_centers")
  # Out-of-integer-range must hit the explicit error, not an as.integer overflow.
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(n_centers = 3e9)),
               "n_centers")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
               "max_anchors_per_class")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
               "max_anchors_per_class")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_dre_ulsif(dat$X, dat$y, hp = list(seed = -1)), "hp\\$seed")
})

test_that("dre-ulsif closed-form coefficients and ratio match an independent recompute", {
  # Manual small-example check of the core uLSIF math. The fit's normal equations
  #   (H + lambda I) alpha = h,  H = (1/n_ctrl) Phi_ctrl^T Phi_ctrl,
  #   h = colMeans(Phi_case),  alpha clipped to its non-negative part
  # are re-derived here from the rCLR anchors and the frozen centres/bandwidth
  # using an INDEPENDENT loop-based kernel, and the held-out density ratio
  # w(z) = sum_l alpha_l phi_l(z) (and its log-score) are cross-checked against
  # score_dre_ulsif on the SAME rows (full overlap -> the score path re-derives an
  # identical representation).
  dat <- .make_dre_ulsif_data(n = 80L, p = 12L, k = 5L, seed = 77L)
  model <- fit_dre_ulsif(dat$X, dat$y,
                         hp = list(max_anchors_per_class = 30L,
                                   n_centers = 12L, lambda = 1e-3, seed = 5L))

  Zc <- .dre_ulsif_rclr_matrix(model$case_anchors)
  Z0 <- .dre_ulsif_rclr_matrix(model$control_anchors)
  C <- model$repr$centers
  sigma <- model$repr$sigma
  expect_equal(C, .dre_ulsif_rclr_matrix(model$center_anchors),
               tolerance = 1e-12)

  Phi_case <- .manual_ulsif_kernel(Zc, C, sigma)
  Phi_ctrl <- .manual_ulsif_kernel(Z0, C, sigma)
  H <- crossprod(Phi_ctrl) / nrow(Z0)
  h <- colMeans(Phi_case)
  b <- nrow(C)
  alpha_manual <- pmax(solve(H + model$hp$lambda * diag(b), h), 0)
  expect_equal(as.numeric(alpha_manual), as.numeric(model$repr$alpha),
               tolerance = 1e-9)

  # Held-out density ratio and log-score reproduced from the manual coefficients.
  eps <- model$hp$eps
  held <- dat$X[c(1L, 2L, 41L, 42L), , drop = FALSE]
  w <- vapply(seq_len(nrow(held)), function(i) {
    z <- .dre_ulsif_rclr(held[i, ])
    sum(as.numeric(model$repr$alpha) *
          as.numeric(.manual_ulsif_kernel(matrix(z, nrow = 1L), C, sigma)))
  }, numeric(1))
  # The non-negative kernel model yields a non-negative (finite) density ratio.
  expect_true(all(is.finite(w)))
  expect_true(all(w >= -1e-12))
  manual_score <- log(pmax(w, eps))
  expect_equal(manual_score, score_dre_ulsif(model, held), tolerance = 1e-9)
})
