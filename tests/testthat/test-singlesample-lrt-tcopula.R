library(testthat)

# Synthetic generator with a CASE-ONLY, MATCHED-MARGINAL feature-DEPENDENCE block.
# The first `k` features carry a 2-group contrast correlation (+rho within each
# half, -rho between halves) for case samples only; controls have those features
# independent. Crucially the per-feature MARGINALS are IDENTICAL across classes
# (mean 4, sd `sd`, both classes) -- the only class difference is the dependence.
#
# Why dependence-only and not a mean shift (unlike the lrt-copula test): at finite
# degrees of freedom the Student-t copula with R = I is NOT the independence
# copula -- it assigns high density to JOINTLY EXTREME rank vectors (any
# direction). A mean-shift block makes case specimens marginally extreme under the
# CONTROL marginals, so the control copula (R ~ I) rewards that extremeness and
# can win, flipping the score sign. (The Gaussian copula cancels this exactly via
# its -I term; the t copula does not at finite df.) With matched marginals a case
# specimen is NOT marginally extreme under the control marginals, so the LRT reads
# the pure DEPENDENCE -- exactly what the copula method is built to detect. The
# balanced 2-group contrast is mean-zero across the block, so it is preserved by
# the within-sample rCLR centring. Abundances are exp() so the input is strictly
# positive (compositional).
.make_lrt_tcopula_data <- function(n = 160L, p = 30L, k = 16L, rho = 0.7,
                                   sd = 0.5, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  half <- k %/% 2L
  Rb <- matrix(-rho, k, k)
  Rb[seq_len(half), seq_len(half)] <- rho
  Rb[(half + 1L):k, (half + 1L):k] <- rho
  diag(Rb) <- 1
  Lb <- chol(Rb)                                   # Lb^T Lb = Rb
  B <- matrix(stats::rnorm(n * k), n, k)           # iid N(0,1) block (controls)
  case <- which(y == 1L)
  B[case, ] <- B[case, , drop = FALSE] %*% Lb      # case-only block dependence
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  L[, seq_len(k)] <- 4 + sd * B                    # identical N(4, sd^2) marginals
  X <- exp(L)
  list(X = X, y = y, block = features[seq_len(k)])
}

.auc_or_wilcoxon_tcop <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("lrt-tcopula fit/score roundtrip has the right shape and types", {
  dat <- .make_lrt_tcopula_data()
  model <- fit_lrt_tcopula(dat$X, dat$y,
                           hp = list(max_anchors_per_class = 40L, seed = 7L))
  expect_s3_class(model, "lrt_tcopula_model")
  expect_equal(nrow(model$case_anchors), 40L)
  expect_equal(nrow(model$control_anchors), 40L)
  expect_identical(model$feature_universe, colnames(dat$X))
  expect_equal(model$hp$df, 5)                    # default degrees of freedom
  # Frozen full-universe precision matrices are symmetric and the right size.
  p <- ncol(dat$X)
  expect_equal(dim(model$repr$Rinv_case), c(p, p))
  expect_equal(dim(model$repr$Rinv_control), c(p, p))

  score <- score_lrt_tcopula(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("lrt-tcopula separates the planted class-specific dependence block", {
  # Larger cohort so the held-out AUC is a stable estimate of the (genuine, copula-
  # only) dependence signal rather than a small-sample fluctuation.
  dat <- .make_lrt_tcopula_data(n = 300L, seed = 202L)
  ctrl <- which(dat$y == 0L)
  case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_lrt_tcopula(dat$X[train, ], dat$y[train],
                           hp = list(max_anchors_per_class = 112L, seed = 11L))
  score <- score_lrt_tcopula(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_tcop(dat$y[test], score), 0.7)
})

test_that("lrt-tcopula passes the canonical row-equivariance gate", {
  dat <- .make_lrt_tcopula_data(seed = 43L)
  model <- fit_lrt_tcopula(dat$X[1:120, ], dat$y[1:120],
                           hp = list(max_anchors_per_class = 30L, seed = 13L))
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_lrt_tcopula(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("lrt-tcopula scoring is exactly invariant to per-sample scaling", {
  dat <- .make_lrt_tcopula_data(seed = 44L)
  model <- fit_lrt_tcopula(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]
  expect_equal(score_lrt_tcopula(model, row * 7),
               score_lrt_tcopula(model, row), tolerance = 1e-8)
  expect_equal(score_lrt_tcopula(model, row * 1e6),
               score_lrt_tcopula(model, row), tolerance = 1e-8)
  # Row-equivariant batch: scaling some rows leaves the others' scores fixed.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e6
  expect_equal(score_lrt_tcopula(model, scaled),
               score_lrt_tcopula(model, batch), tolerance = 1e-8)
  # A row with mostly zeros (rCLR over its own nonzero parts) still scores finite.
  X_sparse <- dat$X[1:3, , drop = FALSE]
  X_sparse[1, 3:ncol(X_sparse)] <- 0
  expect_true(all(is.finite(score_lrt_tcopula(model, X_sparse))))
})

test_that("lrt-tcopula handles partial feature overlap (consistent + single-row)", {
  dat <- .make_lrt_tcopula_data(seed = 45L)
  model <- fit_lrt_tcopula(dat$X, dat$y)
  # Subset and REORDER the columns: present = intersect(universe, colnames) keeps
  # feature-universe order, so the score is invariant to the X column order.
  keep <- rev(c(dat$block[1:5], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_lrt_tcopula(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_lrt_tcopula(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_lrt_tcopula(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("lrt-tcopula matches the canonical dispatch", {
  dat <- .make_lrt_tcopula_data(seed = 46L)
  model <- fit_lrt_tcopula(dat$X[1:120, ], dat$y[1:120],
                           hp = list(max_anchors_per_class = 30L, seed = 17L))
  X_test <- dat$X[121:160, ]
  direct <- score_lrt_tcopula(model, X_test)

  # Post-integration the manifest carries the lrt-tcopula row and score_lrt_tcopula
  # lives in the OmicSelector namespace, so the canonical adapter resolves it
  # directly (the real roster -> package dispatch path). In the pre-integration
  # staged self-test the manifest row and namespace entry do not yet exist, so we
  # synthesise a one-row roster (cloned from the committed lrt-copula row) and
  # register a thin adapter to exercise the same dispatch wiring. Either way the
  # dispatched score must equal the direct call.
  roster <- singlesample_method_roster()
  if (!"lrt-tcopula" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "lrt-tcopula"
    tmpl$fit_fn <- "fit_lrt_tcopula"
    tmpl$score_fn <- "score_lrt_tcopula"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_lrt_tcopula", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "lrt-tcopula",
      function(model, X, meta) score_lrt_tcopula(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("lrt-tcopula", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("lrt-tcopula returns the neutral 0 below the feature-overlap floor", {
  dat <- .make_lrt_tcopula_data(seed = 47L)
  model <- fit_lrt_tcopula(dat$X, dat$y, hp = list(min_features = 3L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, colnames(dat$X)[1:2], drop = FALSE]
  expect_equal(score_lrt_tcopula(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_lrt_tcopula(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
})

test_that("lrt-tcopula fitting/scoring is deterministic and RNG-safe", {
  dat <- .make_lrt_tcopula_data(seed = 48L)
  hp <- list(max_anchors_per_class = 25L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_lrt_tcopula(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_lrt_tcopula(dat$X, dat$y, hp = hp)
  expect_identical(model_1, model_2)

  s1 <- score_lrt_tcopula(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_lrt_tcopula(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("lrt-tcopula fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_lrt_tcopula_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  # max_anchors_per_class > class size -> no subsampling -> RNG never touched.
  model <- fit_lrt_tcopula(dat$X, dat$y,
                           hp = list(max_anchors_per_class = 200L, seed = 29L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "lrt_tcopula_model")
})

test_that("lrt-tcopula leaves no RNG state when subsampling fits in a fresh session", {
  # Companion to the no-subsampling case above, exercising the OTHER branch of the
  # anchor-subsampling RNG save/restore: here a class EXCEEDS the anchor cap, so
  # set.seed() IS called inside the fit. With no prior global seed, the on.exit
  # restore must rm() the .Random.seed that set.seed() created -- not leave it
  # behind. (80 anchors/class > cap 20 -> subsampling path taken.)
  dat <- .make_lrt_tcopula_data(seed = 53L)       # this set.seed() seeds the RNG
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())       # wipe to a genuine no-seed state
  }
  model <- fit_lrt_tcopula(dat$X, dat$y,
                           hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_lte(nrow(model$case_anchors), 20L)        # subsampling actually happened
  expect_lte(nrow(model$control_anchors), 20L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # And subsampling is deterministic across fits (frozen seed).
  model2 <- fit_lrt_tcopula(dat$X, dat$y,
                            hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_identical(rownames(model$case_anchors), rownames(model2$case_anchors))
})

test_that("lrt-tcopula handles a constant/degenerate training feature without NA", {
  # A feature that is zero across ALL anchors of a class has a constant (all-zero)
  # rCLR column. Its empirical marginal collapses to a SINGLE knot (no approx call,
  # PIT -> the knot plotting position, q ~ 0) and its Student-t-score column has
  # zero variance, so cor() yields NA for that row/column. The marginal single-knot
  # guard and the correlation NA->identity guard must keep the fit and every score
  # finite.
  dat <- .make_lrt_tcopula_data(seed = 54L)
  X <- dat$X
  X[, ncol(X)] <- 0                       # feature constant (zero) in BOTH classes
  X[dat$y == 1L, ncol(X) - 1L] <- 0       # feature constant (zero) in cases only
  model <- fit_lrt_tcopula(X, dat$y, hp = list(max_anchors_per_class = 40L))
  expect_true(all(is.finite(model$repr$Rinv_case)))
  expect_true(all(is.finite(model$repr$Rinv_control)))
  expect_true(is.finite(model$repr$logdet_case))
  expect_true(is.finite(model$repr$logdet_control))
  score <- score_lrt_tcopula(model, X)
  expect_length(score, nrow(X))
  expect_true(all(is.finite(score)))
  # Scoring a specimen that is itself zero on the degenerate features is still finite.
  expect_true(all(is.finite(score_lrt_tcopula(model, X[1:4, , drop = FALSE]))))
})

test_that("lrt-tcopula hyperparameter validation is strict", {
  dat <- .make_lrt_tcopula_data(seed = 50L)
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(shrink = 1)),
               "hp\\$shrink")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(pit_clamp = 0.7)),
               "hp\\$pit_clamp")
  # Both open-interval endpoints are rejected: pit_clamp = 0 would make
  # qt(0, df) = -Inf and pit_clamp = 0.5 collapses the clamp window.
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(pit_clamp = 0)),
               "hp\\$pit_clamp")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(pit_clamp = 0.5)),
               "hp\\$pit_clamp")
  # df must be a single finite number >= 1: df <= 0 and 0 < df < 1 are rejected.
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(df = 0)), "hp\\$df")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(df = -3)), "hp\\$df")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(df = 0.5)), "hp\\$df")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(df = c(3, 5))), "hp\\$df")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(df = Inf)), "hp\\$df")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
               "max_anchors_per_class")
  # Out-of-integer-range must hit the explicit error, not an as.integer overflow.
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
               "max_anchors_per_class")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_lrt_tcopula(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
})

test_that("lrt-tcopula copula log-density ratio is finite with the correct sign", {
  # Tiny 3-feature MATCHED-MARGINAL example. In cases, f1 and f2 share a common
  # latent (positive dependence, corr = rho); in controls they are independent.
  # Per-feature marginals are identical N(3, 0.5^2) across classes (NO mean
  # shift), so a case specimen is not marginally extreme under the control
  # marginals -- the only class difference is the f1-f2 dependence, which is what
  # the copula LRT reads (see the generator comment for why a mean-shift block
  # would not cleanly separate at finite df). A case-consistent specimen (f1, f2
  # jointly elevated, matching the case + dependence) must score above a case-
  # inconsistent one (f1 up while f2 down, opposing it).
  set.seed(123)
  n <- 200L
  y <- rep(c(0, 1), each = n / 2L)
  feats <- c("f1", "f2", "f3")
  rho <- 0.8
  b1 <- stats::rnorm(n)
  b2 <- stats::rnorm(n)
  case <- which(y == 1L)
  g <- stats::rnorm(length(case))                  # shared latent (cases only)
  b1[case] <- sqrt(rho) * g + sqrt(1 - rho) * b1[case]
  b2[case] <- sqrt(rho) * g + sqrt(1 - rho) * b2[case]  # corr=rho, N(0,1) marginal
  L <- matrix(stats::rnorm(n * 3L, 3, 0.5), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), feats))
  L[, 1] <- 3 + 0.5 * b1
  L[, 2] <- 3 + 0.5 * b2                           # matched marginals across classes
  X <- exp(L)
  model <- fit_lrt_tcopula(X, y, hp = list(df = 5))

  x_caselike <- matrix(exp(c(4.0, 4.0, 3.0)), nrow = 1,
                       dimnames = list(NULL, feats))   # f1,f2 jointly high
  x_ctrllike <- matrix(exp(c(4.0, 2.0, 3.0)), nrow = 1,
                       dimnames = list(NULL, feats))   # f1 high, f2 low
  s_case <- score_lrt_tcopula(model, x_caselike)
  s_ctrl <- score_lrt_tcopula(model, x_ctrllike)
  expect_true(is.finite(s_case) && is.finite(s_ctrl))
  expect_length(s_case, 1L)
  expect_gt(s_case, s_ctrl)

  # On the training batch the case group scores higher on average.
  s_all <- score_lrt_tcopula(model, X)
  expect_true(all(is.finite(s_all)))
  expect_gt(mean(s_all[y == 1]), mean(s_all[y == 0]))
})

test_that("lrt-tcopula log-density helper matches an independent joint-minus-marginal computation", {
  # Manual small-example check of the core math. The scorer's combined formula
  #   l_c(q) = log f_MVt(q; R, nu) - sum_j log f_t(q_j; nu)
  # is re-derived here independently from the multivariate-t and univariate-t
  # densities (NOT the simplified combined formula), and the univariate-t part is
  # additionally cross-checked against base R's dt(). q and R are hand-chosen.
  nu <- 6
  R <- rbind(c(1.0, 0.5, 0.2),
             c(0.5, 1.0, 0.3),
             c(0.2, 0.3, 1.0))
  ch <- chol(R)
  Rinv <- chol2inv(ch)
  logdet <- 2 * sum(log(diag(ch)))
  q <- c(0.7, -1.3, 0.4)
  d <- length(q)

  # Independent multivariate-t log-density (keeps the (d/2) log(nu*pi) term).
  log_fmvt <- lgamma((nu + d) / 2) - lgamma(nu / 2) - (d / 2) * log(nu * pi) -
    0.5 * logdet -
    ((nu + d) / 2) * log1p(as.numeric(t(q) %*% Rinv %*% q) / nu)
  # Independent univariate-t log-density per coordinate, cross-checked vs dt().
  log_ft <- function(x, nu) {
    lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi) -
      ((nu + 1) / 2) * log1p(x^2 / nu)
  }
  expect_equal(vapply(q, log_ft, numeric(1), nu = nu),
               stats::dt(q, df = nu, log = TRUE), tolerance = 1e-12)

  manual_l <- log_fmvt - sum(vapply(q, log_ft, numeric(1), nu = nu))
  expect_equal(.lrt_tcopula_logdens(q, Rinv, logdet, nu), manual_l,
               tolerance = 1e-10)

  # Independence sanity: with R = I and standard-t scores the copula log-density
  # is exactly 0 (joint t with diagonal scale equals the product of marginals
  # only in the Gaussian limit, so for finite nu the value is NONzero in general;
  # here we instead confirm the helper is finite and that a positively dependent
  # specimen scores above an anti-dependent one under the same R).
  q_pos <- c(1.2, 1.1, 0.0)     # aligned with the positive off-diagonals of R
  q_neg <- c(1.2, -1.1, 0.0)    # opposes the positive z1-z2 correlation
  expect_gt(.lrt_tcopula_logdens(q_pos, Rinv, logdet, nu),
            .lrt_tcopula_logdens(q_neg, Rinv, logdet, nu))
})

test_that("lrt-tcopula converges to lrt-copula in the Gaussian (large-df) limit", {
  # As nu -> Inf, qt(.,nu) -> qnorm(.) and the t-copula log-density ratio -> the
  # Gaussian-copula log-density ratio, so lrt-tcopula(df = 1e6) must match the
  # committed in-package lrt-copula on the same data, same anchors, same marginals.
  # max_anchors_per_class exceeds the class size, so neither method subsamples and
  # both freeze identical anchors -- the only difference is qt(.,1e6) vs qnorm.
  dat <- .make_lrt_tcopula_data(seed = 71L)
  hp_shared <- list(shrink = 0.05, pit_clamp = 1e-4,
                    max_anchors_per_class = 200L, min_features = 3L,
                    eps = 1e-6, seed = 1L)
  model_t <- fit_lrt_tcopula(dat$X, dat$y, hp = c(hp_shared, list(df = 1e6)))
  model_g <- fit_lrt_copula(dat$X, dat$y, hp = hp_shared)

  X_eval <- dat$X[1:40, ]
  s_t <- score_lrt_tcopula(model_t, X_eval)
  s_g <- score_lrt_copula(model_g, X_eval)
  expect_equal(length(s_t), length(s_g))
  # Both finite and tightly agreeing in the Gaussian limit.
  expect_true(all(is.finite(s_t)))
  # Tight, meaningful bound: observed max|s_t - s_g| ~ 3.2e-3 on a score scale of
  # ~49 (relative ~7e-5). The 2e-2 threshold is the deterministic-seed guard.
  expect_lt(max(abs(s_t - s_g)), 2e-2)
})

test_that("lrt-tcopula R=I copula is tail-aware but not purely radial", {
  # Regression lock for the documented R=I tail character (Codex R1 P1). The
  # internal log-density .lrt_tcopula_logdens(q, Rinv, logdet, nu) is the copula
  # log-density; at Rinv = I, logdet = 0 it reduces to the joint-t radial term plus
  # the coordinatewise marginal-correction term ((nu+1)/2) sum log1p(q_j^2/nu).
  I2 <- diag(2)
  nu <- 5

  # (a) Sign- and permutation-invariance: depends on q only through q_j^2.
  expect_equal(.lrt_tcopula_logdens(c(2.3, -1.7), I2, 0, nu),
               .lrt_tcopula_logdens(c(-2.3, 1.7), I2, 0, nu), tolerance = 1e-12)
  expect_equal(.lrt_tcopula_logdens(c(2.3, -1.7), I2, 0, nu),
               .lrt_tcopula_logdens(c(1.7, 2.3), I2, 0, nu), tolerance = 1e-12)

  # (b) NOT purely radial: at the SAME radius (||q||^2 = 9), mass spread across two
  #     coordinates scores strictly HIGHER than a single dominant coordinate
  #     (the marginal term is concave in q_j^2). This is the exact Codex example.
  l_concentrated <- .lrt_tcopula_logdens(c(3, 0), I2, 0, nu)            # ~ -0.4154
  l_spread       <- .lrt_tcopula_logdens(c(sqrt(4.5), sqrt(4.5)), I2, 0, nu)  # ~ +0.3468
  expect_equal(sum(c(3, 0)^2), sum(c(sqrt(4.5), sqrt(4.5))^2))          # equal radius
  expect_gt(l_spread, l_concentrated)
  expect_equal(l_concentrated, -0.415448, tolerance = 1e-5)
  expect_equal(l_spread,        0.346817, tolerance = 1e-5)

  # (c) The non-radiality is a finite-nu (tail) effect: it vanishes in the Gaussian
  #     limit, where the copula density at R=I is identically 0 for every q.
  expect_lt(abs(.lrt_tcopula_logdens(c(3, 0), I2, 0, 1e8)), 1e-3)
  expect_lt(abs(.lrt_tcopula_logdens(c(sqrt(4.5), sqrt(4.5)), I2, 0, 1e8)), 1e-3)
})
