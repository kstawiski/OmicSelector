library(testthat)

# Synthetic generator with a CASE-ONLY MARGINAL (location) shift in a feature
# block. The first `k` of `p` log-abundance features are elevated by `shift` for
# case samples. After exp() and the within-sample rCLR centring, those block
# features sit higher RELATIVE to the per-sample geometric mean for cases, and the
# remaining features sit correspondingly lower; both directions are per-feature
# class-conditional MARGINAL location differences. This is exactly the signal the
# naive-Bayes KDE log-likelihood-ratio is built to detect (unlike the copula
# methods, which match marginals out and read only dependence). Abundances are
# exp() so the input is strictly positive (compositional).
.make_lrt_nbkde_data <- function(n = 160L, p = 30L, k = 10L, shift = 1.3,
                              sd = 0.5, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case <- which(y == 1L)
  L[case, seq_len(k)] <- L[case, seq_len(k)] + shift  # case-only marginal shift
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

test_that("lrt-nbkde fit/score roundtrip has the right shape and types", {
  dat <- .make_lrt_nbkde_data()
  model <- fit_lrt_nbkde(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 40L, seed = 7L))
  expect_s3_class(model, "lrt_nbkde_model")
  expect_equal(nrow(model$case_anchors), 40L)
  expect_equal(nrow(model$control_anchors), 40L)
  expect_identical(model$feature_universe, colnames(dat$X))
  # Frozen full-universe representation: per-class rCLR anchor matrices and one
  # bandwidth per feature.
  p <- ncol(dat$X)
  expect_equal(dim(model$repr$z_case), c(40L, p))
  expect_equal(dim(model$repr$z_control), c(40L, p))
  expect_length(model$repr$bw, p)
  expect_true(all(is.finite(model$repr$bw)) && all(model$repr$bw > 0))

  score <- score_lrt_nbkde(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("lrt-nbkde separates a planted case-vs-control marginal shift", {
  # Larger cohort + held-out split so the AUC is a stable estimate of the genuine
  # per-feature marginal signal rather than a small-sample fluctuation.
  dat <- .make_lrt_nbkde_data(n = 300L, seed = 202L)
  ctrl <- which(dat$y == 0L)
  case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_lrt_nbkde(dat$X[train, ], dat$y[train],
                      hp = list(max_anchors_per_class = 112L, seed = 11L))
  score <- score_lrt_nbkde(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_bw(dat$y[test], score), 0.7)
})

test_that("lrt-nbkde passes the canonical row-equivariance gate", {
  dat <- .make_lrt_nbkde_data(seed = 43L)
  model <- fit_lrt_nbkde(dat$X[1:120, ], dat$y[1:120],
                      hp = list(max_anchors_per_class = 30L, seed = 13L))
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_lrt_nbkde(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("lrt-nbkde scoring is exactly invariant to per-sample scaling", {
  dat <- .make_lrt_nbkde_data(seed = 44L)
  model <- fit_lrt_nbkde(dat$X, dat$y)
  row <- dat$X[3, , drop = FALSE]
  expect_equal(score_lrt_nbkde(model, row * 7),
               score_lrt_nbkde(model, row), tolerance = 1e-8)
  expect_equal(score_lrt_nbkde(model, row * 1e6),
               score_lrt_nbkde(model, row), tolerance = 1e-8)
  # Row-equivariant batch: scaling some rows leaves the others' scores fixed.
  batch <- dat$X[1:6, , drop = FALSE]
  scaled <- batch
  scaled[c(2L, 5L), ] <- scaled[c(2L, 5L), ] * 1e6
  expect_equal(score_lrt_nbkde(model, scaled),
               score_lrt_nbkde(model, batch), tolerance = 1e-8)
  # A row with mostly zeros (rCLR over its own nonzero parts) still scores finite.
  X_sparse <- dat$X[1:3, , drop = FALSE]
  X_sparse[1, 3:ncol(X_sparse)] <- 0
  expect_true(all(is.finite(score_lrt_nbkde(model, X_sparse))))
})

test_that("lrt-nbkde handles partial feature overlap (consistent + single-row)", {
  dat <- .make_lrt_nbkde_data(seed = 45L)
  model <- fit_lrt_nbkde(dat$X, dat$y)
  # Subset and REORDER the columns: present = intersect(universe, colnames) keeps
  # feature-universe order, so the score is invariant to the X column order.
  keep <- rev(c(dat$block[1:5], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_lrt_nbkde(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_lrt_nbkde(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores.
  expect_equal(score_lrt_nbkde(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("lrt-nbkde matches the canonical dispatch", {
  dat <- .make_lrt_nbkde_data(seed = 46L)
  model <- fit_lrt_nbkde(dat$X[1:120, ], dat$y[1:120],
                      hp = list(max_anchors_per_class = 30L, seed = 17L))
  X_test <- dat$X[121:160, ]
  direct <- score_lrt_nbkde(model, X_test)

  # Post-integration the manifest carries the lrt-nbkde row and score_lrt_nbkde lives in
  # the OmicSelector namespace, so the canonical adapter resolves it directly (the
  # real roster -> package dispatch path). In the pre-integration staged self-test
  # the namespace entry does not yet exist (and the manifest row may or may not yet
  # be present), so we clone a one-row roster from a committed within-discriminator
  # row when needed and register a thin adapter to exercise the same dispatch
  # wiring. Either way the dispatched score must equal the direct call.
  roster <- singlesample_method_roster()
  if (!"lrt-nbkde" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "lrt-nbkde"
    tmpl$fit_fn <- "fit_lrt_nbkde"
    tmpl$score_fn <- "score_lrt_nbkde"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_lrt_nbkde", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "lrt-nbkde",
      function(model, X, meta) score_lrt_nbkde(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("lrt-nbkde", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("lrt-nbkde returns the neutral 0 below the feature-overlap floor", {
  dat <- .make_lrt_nbkde_data(seed = 47L)
  model <- fit_lrt_nbkde(dat$X, dat$y, hp = list(min_features = 3L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, colnames(dat$X)[1:2], drop = FALSE]
  expect_equal(score_lrt_nbkde(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_lrt_nbkde(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
})

test_that("lrt-nbkde fitting/scoring is deterministic and RNG-safe", {
  dat <- .make_lrt_nbkde_data(seed = 48L)
  hp <- list(max_anchors_per_class = 25L, seed = 23L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_lrt_nbkde(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_lrt_nbkde(dat$X, dat$y, hp = hp)
  expect_identical(model_1, model_2)

  s1 <- score_lrt_nbkde(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_lrt_nbkde(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  # Scoring never touches the global RNG.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("lrt-nbkde fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_lrt_nbkde_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  # max_anchors_per_class > class size (80) -> no subsampling -> RNG never touched.
  model <- fit_lrt_nbkde(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 200L, seed = 29L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "lrt_nbkde_model")
  expect_equal(nrow(model$case_anchors), 80L)         # all anchors kept
})

test_that("lrt-nbkde leaves no RNG state when subsampling fits in a fresh session", {
  # Companion to the no-subsampling case above, exercising the OTHER branch of the
  # anchor-subsampling RNG save/restore: here a class EXCEEDS the anchor cap, so
  # set.seed() IS called inside the fit. With no prior global seed, the on.exit
  # restore must rm() the .Random.seed that set.seed() created -- not leave it
  # behind. (80 anchors/class > cap 20 -> subsampling path taken.)
  dat <- .make_lrt_nbkde_data(seed = 53L)               # this set.seed() seeds the RNG
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())          # wipe to a genuine no-seed state
  }
  model <- fit_lrt_nbkde(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_lte(nrow(model$case_anchors), 20L)           # subsampling actually happened
  expect_lte(nrow(model$control_anchors), 20L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # And subsampling is deterministic across fits (frozen seed).
  model2 <- fit_lrt_nbkde(dat$X, dat$y,
                       hp = list(max_anchors_per_class = 20L, seed = 31L))
  expect_identical(rownames(model$case_anchors), rownames(model2$case_anchors))
})

test_that("lrt-nbkde handles a constant/degenerate training feature without NA", {
  # A feature that is zero across ALL anchors of a class has a constant (all-zero)
  # rCLR column. The pooled bw.nrd0 stays positive even for a constant input (it
  # falls back to a unit scale) and is floored at hp$eps, and the Gaussian kernel
  # sum is evaluated through a stable log-sum-exp, so the per-feature KDE
  # log-density -- and therefore the fit and every score -- stays finite.
  dat <- .make_lrt_nbkde_data(seed = 54L)
  X <- dat$X
  X[, ncol(X)] <- 0                       # feature constant (zero) in BOTH classes
  X[dat$y == 1L, ncol(X) - 1L] <- 0       # feature constant (zero) in cases only
  model <- fit_lrt_nbkde(X, dat$y, hp = list(max_anchors_per_class = 40L))
  expect_true(all(is.finite(model$repr$bw)) && all(model$repr$bw > 0))
  expect_true(all(is.finite(model$repr$z_case)))
  expect_true(all(is.finite(model$repr$z_control)))
  score <- score_lrt_nbkde(model, X)
  expect_length(score, nrow(X))
  expect_true(all(is.finite(score)))
  # Scoring a specimen that is itself zero on the degenerate features is finite.
  expect_true(all(is.finite(score_lrt_nbkde(model, X[1:4, , drop = FALSE]))))
})

test_that("lrt-nbkde hyperparameter validation is strict", {
  dat <- .make_lrt_nbkde_data(seed = 50L)
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  # bw must be NULL or a single positive finite number.
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(bw = 0)), "hp\\$bw")
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(bw = -1)), "hp\\$bw")
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(bw = c(1, 2))), "hp\\$bw")
  # bw_adjust must be a single positive finite number.
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(bw_adjust = 0)),
               "hp\\$bw_adjust")
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(bw_adjust = -2)),
               "hp\\$bw_adjust")
  # Integer fields: non-integer and out-of-integer-range both rejected (the latter
  # must hit the explicit error, not an as.integer overflow).
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
               "max_anchors_per_class")
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
               "max_anchors_per_class")
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_lrt_nbkde(dat$X, dat$y, hp = list(seed = -1)), "hp\\$seed")
})

test_that("lrt-nbkde KDE log-density and score match an independent dnorm recompute", {
  # Manual small-example check of the core naive-Bayes KDE math. The per-feature
  # class KDE log-density log f_cj(t) = log( mean_i dnorm(t, z_cj_i, bw_j) ) is
  # re-derived here with an INDEPENDENT vectorized stats::dnorm sum, and the full
  # score is reproduced as the summed per-feature log-density ratio. At full
  # overlap the score path re-derives a representation identical to model$repr.
  dat <- .make_lrt_nbkde_data(n = 80L, p = 12L, k = 5L, seed = 77L)
  model <- fit_lrt_nbkde(dat$X, dat$y,
                      hp = list(max_anchors_per_class = 30L, seed = 5L))

  # The stored full-universe representation equals a fresh re-derivation.
  repr2 <- .lrt_nbkde_class_repr(model$case_anchors, model$control_anchors, model$hp)
  expect_equal(model$repr$bw, repr2$bw, tolerance = 1e-12)
  expect_equal(model$repr$z_case, repr2$z_case, tolerance = 1e-12)

  repr <- model$repr
  held <- dat$X[c(1L, 2L, 41L, 42L), , drop = FALSE]
  p <- ncol(held)

  manual_score <- vapply(seq_len(nrow(held)), function(i) {
    z <- .lrt_nbkde_rclr(held[i, ])
    # Independent per-feature KDE log-densities via dnorm.
    lc <- vapply(seq_len(p), function(j) {
      log(mean(stats::dnorm(z[j], mean = repr$z_case[, j], sd = repr$bw[j])))
    }, numeric(1))
    l0 <- vapply(seq_len(p), function(j) {
      log(mean(stats::dnorm(z[j], mean = repr$z_control[, j], sd = repr$bw[j])))
    }, numeric(1))
    sum(lc - l0)
  }, numeric(1))

  expect_equal(manual_score, score_lrt_nbkde(model, held), tolerance = 1e-9)

  # And the vectorized helper used inside the scorer matches the dnorm formula for
  # a single specimen / feature.
  z1 <- .lrt_nbkde_rclr(held[1, ])
  helper_case <- .lrt_nbkde_class_logdens(z1, repr$z_case, repr$bw)
  dnorm_case <- vapply(seq_len(p), function(j) {
    log(mean(stats::dnorm(z1[j], mean = repr$z_case[, j], sd = repr$bw[j])))
  }, numeric(1))
  expect_equal(helper_case, dnorm_case, tolerance = 1e-10)
})

test_that("lrt-nbkde bw_adjust is not silently aliased to bw (partial-match guard)", {
  # Regression lock for the hp[["bw"]] exact-extraction fix (Claude R1 P2). `bw` is a
  # unique prefix of `bw_adjust`, and R's `$` does partial matching on lists, so a
  # reader written as hp$bw would silently resolve hp$bw_adjust when only bw_adjust is
  # supplied. A bw_adjust=2 fit (= 2 x pooled bw.nrd0 per feature) must therefore
  # produce DIFFERENT scores than a bw=2 fit (= fixed bandwidth 2 for every feature),
  # not identical ones. This test FAILS if the resolver/bandwidth helper reverts to `$`.
  dat <- .make_lrt_nbkde_data(seed = 61L)
  m_adjust <- fit_lrt_nbkde(dat$X, dat$y, hp = list(bw_adjust = 2))
  m_fixed  <- fit_lrt_nbkde(dat$X, dat$y, hp = list(bw = 2))
  s_adjust <- score_lrt_nbkde(m_adjust, dat$X[1:20, ])
  s_fixed  <- score_lrt_nbkde(m_fixed,  dat$X[1:20, ])
  expect_true(all(is.finite(s_adjust)) && all(is.finite(s_fixed)))
  expect_false(isTRUE(all.equal(s_adjust, s_fixed)))
})
