library(testthat)

# Synthetic generator with a CASE-ONLY, MATCHED-MARGINAL feature-DEPENDENCE block
# CONCENTRATED on the top-abundance features the vine universe selects. The first
# `k` features (the highest-abundance ones, mean `top_mean` > the rest) carry a
# 2-group contrast correlation (+rho within each half, -rho between halves) for
# case samples only; controls have those features independent. The per-feature
# MARGINALS are IDENTICAL across classes (same mean, same sd, both classes) -- the
# only class difference is the dependence, which is exactly what a copula (here a
# vine) reads. The balanced 2-group contrast is mean-zero across the block, so it
# is preserved by the within-sample rCLR centring. Abundances are exp() so the
# input is strictly positive (compositional). The block features carry the larger
# mean so they survive the top-d vine feature selection.
.make_lrt_vine_data <- function(n = 200L, p = 30L, k = 8L, rho = 0.7,
                                sd = 0.5, top_mean = 6, base_mean = 4,
                                seed = 41L) {
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
  # Base panel: low-abundance features at base_mean; block features at top_mean so
  # they are the ones the top-d universe picks.
  L <- matrix(stats::rnorm(n * p, mean = base_mean, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  L[, seq_len(k)] <- top_mean + sd * B             # identical marginals across classes
  X <- exp(L)
  list(X = X, y = y, block = features[seq_len(k)])
}

.auc_or_wilcoxon_vine <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("lrt-vinecopula fit/score roundtrip has the right shape and types", {
  dat <- .make_lrt_vine_data()
  model <- fit_lrt_vinecopula(dat$X, dat$y,
                              hp = list(max_anchors_per_class = 40L, seed = 7L,
                                        n_vine_features = 8L))
  expect_s3_class(model, "lrt_vinecopula_model")
  expect_equal(nrow(model$case_anchors), 40L)
  expect_equal(nrow(model$control_anchors), 40L)
  # The frozen universe is the top-8 by training column mean = the 8 block features.
  expect_equal(length(model$feature_universe), 8L)
  expect_setequal(model$feature_universe, dat$block)
  # Anchors are restricted to the frozen universe.
  expect_equal(ncol(model$case_anchors), 8L)
  expect_true(inherits(model$repr$vine_case, "vinecop"))
  expect_true(inherits(model$repr$vine_control, "vinecop"))

  score <- score_lrt_vinecopula(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
})

test_that("lrt-vinecopula separates a dependence-only contrast (matched marginals)", {
  # The decisive value-add test: case vs control differ ONLY in the dependence
  # structure among the top features (identical per-feature marginals). A naive-
  # Bayes / independent-marginal scorer cannot separate these; the vine copula can.
  dat <- .make_lrt_vine_data(n = 360L, k = 8L, rho = 0.75, seed = 202L)
  ctrl <- which(dat$y == 0L)
  case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 130L), utils::head(case, 130L))
  test <- c(utils::tail(ctrl, 50L), utils::tail(case, 50L))
  model <- fit_lrt_vinecopula(dat$X[train, ], dat$y[train],
                              hp = list(max_anchors_per_class = 130L, seed = 11L,
                                        n_vine_features = 8L))
  score <- score_lrt_vinecopula(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_or_wilcoxon_vine(dat$y[test], score), 0.8)
})

test_that("lrt-vinecopula passes the canonical row-equivariance gate", {
  dat <- .make_lrt_vine_data(seed = 43L)
  model <- fit_lrt_vinecopula(dat$X[1:120, ], dat$y[1:120],
                              hp = list(max_anchors_per_class = 30L, seed = 13L,
                                        n_vine_features = 8L))
  X_test <- dat$X[121:200, ]
  score_fun <- function(model, X, meta) score_lrt_vinecopula(model, X, meta)
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

test_that("lrt-vinecopula scoring is exactly invariant to per-sample scaling", {
  dat <- .make_lrt_vine_data(seed = 44L)
  model <- fit_lrt_vinecopula(dat$X, dat$y, hp = list(n_vine_features = 8L))
  row <- dat$X[3, , drop = FALSE]
  expect_equal(score_lrt_vinecopula(model, row * 7),
               score_lrt_vinecopula(model, row), tolerance = 1e-8)
  expect_equal(score_lrt_vinecopula(model, row * 1e6),
               score_lrt_vinecopula(model, row), tolerance = 1e-8)
  # Per-row random positive scaling leaves every score fixed.
  set.seed(5)
  batch <- dat$X[1:8, , drop = FALSE]
  scaled <- batch * exp(stats::rnorm(nrow(batch), sd = 4))
  expect_equal(score_lrt_vinecopula(model, scaled),
               score_lrt_vinecopula(model, batch), tolerance = 1e-6)
  # Row-equivariant batch: scaling some rows leaves the others' scores fixed.
  scaled2 <- batch
  scaled2[c(2L, 5L), ] <- scaled2[c(2L, 5L), ] * 1e6
  expect_equal(score_lrt_vinecopula(model, scaled2),
               score_lrt_vinecopula(model, batch), tolerance = 1e-6)
  # A row with mostly zeros (rCLR over its own nonzero parts) still scores finite.
  X_sparse <- dat$X[1:3, , drop = FALSE]
  X_sparse[1, 3:ncol(X_sparse)] <- 0
  expect_true(all(is.finite(score_lrt_vinecopula(model, X_sparse))))
})

test_that("lrt-vinecopula handles partial feature overlap (consistent + single-row)", {
  dat <- .make_lrt_vine_data(seed = 45L)
  model <- fit_lrt_vinecopula(dat$X, dat$y, hp = list(n_vine_features = 8L))
  # Present = intersect(universe, colnames) keeps the frozen feature-universe order,
  # so the score is invariant to the X column order. Keep 6 of the 8 universe
  # features (still >= min_features) plus some irrelevant columns, reordered.
  keep <- rev(c(dat$block[1:6], colnames(dat$X)[20:26]))
  X_partial <- dat$X[1:12, keep, drop = FALSE]
  s_partial <- score_lrt_vinecopula(model, X_partial)
  expect_length(s_partial, nrow(X_partial))
  expect_true(all(is.finite(s_partial)))
  # Single-row scoring equals the batch entry for that row.
  expect_equal(score_lrt_vinecopula(model, X_partial[5, , drop = FALSE]),
               s_partial[5], tolerance = 1e-12)
  # Reordering columns does not change the scores (column-permutation invariance).
  expect_equal(score_lrt_vinecopula(model, X_partial[, sample(ncol(X_partial))]),
               s_partial, tolerance = 1e-12)
})

test_that("lrt-vinecopula matches the canonical dispatch", {
  dat <- .make_lrt_vine_data(seed = 46L)
  model <- fit_lrt_vinecopula(dat$X[1:120, ], dat$y[1:120],
                              hp = list(max_anchors_per_class = 30L, seed = 17L,
                                        n_vine_features = 8L))
  X_test <- dat$X[121:200, ]
  direct <- score_lrt_vinecopula(model, X_test)

  # Post-integration the manifest carries the lrt-vinecopula row and
  # score_lrt_vinecopula lives in the OmicSelector namespace, so the canonical
  # adapter resolves it directly. In the pre-integration staged self-test the
  # manifest row and namespace entry do not yet exist, so we synthesise a one-row
  # roster (cloned from the committed lrt-copula row) and register a thin adapter to
  # exercise the same dispatch wiring. Either way the dispatched score must equal
  # the direct call.
  roster <- singlesample_method_roster()
  if (!"lrt-vinecopula" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "lrt-vinecopula"
    tmpl$fit_fn <- "fit_lrt_vinecopula"
    tmpl$score_fn <- "score_lrt_vinecopula"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_lrt_vinecopula", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "lrt-vinecopula",
      function(model, X, meta) score_lrt_vinecopula(model, X, meta))
  }
  expect_equal(
    singlesample_score_call("lrt-vinecopula", model, X_test, roster = roster),
    direct,
    tolerance = 1e-12
  )
})

test_that("lrt-vinecopula returns the neutral 0 below the feature-overlap floor", {
  dat <- .make_lrt_vine_data(seed = 47L)
  model <- fit_lrt_vinecopula(dat$X, dat$y,
                              hp = list(min_features = 3L, n_vine_features = 8L))
  # Only two universe features present -> below min_features=3 -> neutral 0.
  X_two <- dat$X[1:6, model$feature_universe[1:2], drop = FALSE]
  expect_equal(score_lrt_vinecopula(model, X_two), rep(0, nrow(X_two)),
               tolerance = 1e-12)
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2")))
  expect_equal(score_lrt_vinecopula(model, X_none), rep(0, nrow(X_none)),
               tolerance = 1e-12)
})

test_that("lrt-vinecopula fitting/scoring is deterministic and RNG-safe", {
  dat <- .make_lrt_vine_data(seed = 48L)
  hp <- list(max_anchors_per_class = 25L, seed = 23L, n_vine_features = 8L)

  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_lrt_vinecopula(dat$X, dat$y, hp = hp)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_lrt_vinecopula(dat$X, dat$y, hp = hp)
  # Two independent fits give identical scores (vine selection is deterministic).
  s_a <- score_lrt_vinecopula(model_1, dat$X[1:20, ])
  s_b <- score_lrt_vinecopula(model_2, dat$X[1:20, ])
  expect_identical(s_a, s_b)

  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s1 <- score_lrt_vinecopula(model_1, dat$X[1:20, ])
  s2 <- score_lrt_vinecopula(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("lrt-vinecopula fitting leaves the global RNG untouched with no prior seed", {
  dat <- .make_lrt_vine_data(seed = 49L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  # max_anchors_per_class > class size -> no anchor subsampling; the vine fit's own
  # set.seed/on.exit restore must also leave no .Random.seed behind.
  model <- fit_lrt_vinecopula(dat$X, dat$y,
                              hp = list(max_anchors_per_class = 200L, seed = 29L,
                                        n_vine_features = 8L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(model, "lrt_vinecopula_model")
  # Scoring with no prior seed also leaves no RNG state.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  score_lrt_vinecopula(model, dat$X[1:5, ])
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("lrt-vinecopula leaves no RNG state when subsampling fits in a fresh session", {
  # Exercises the anchor-subsampling RNG save/restore: a class EXCEEDS the cap so
  # set.seed() IS called inside the fit. With no prior global seed, the on.exit
  # restore must rm() the .Random.seed that set.seed() created.
  dat <- .make_lrt_vine_data(seed = 53L)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  model <- fit_lrt_vinecopula(dat$X, dat$y,
                              hp = list(max_anchors_per_class = 20L, seed = 31L,
                                        n_vine_features = 8L))
  expect_lte(nrow(model$case_anchors), 20L)
  expect_lte(nrow(model$control_anchors), 20L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # Subsampling is deterministic across fits (frozen seed).
  model2 <- fit_lrt_vinecopula(dat$X, dat$y,
                               hp = list(max_anchors_per_class = 20L, seed = 31L,
                                         n_vine_features = 8L))
  expect_identical(rownames(model$case_anchors), rownames(model2$case_anchors))
})

test_that("lrt-vinecopula handles a constant/degenerate training feature without NA", {
  # A feature that is zero across ALL anchors of a class has a constant (all-zero)
  # rCLR column. Its empirical marginal collapses to a SINGLE knot (PIT -> the knot
  # plotting position) and its near-constant pseudo-observation column carries no
  # dependence; vinecop selects it to the independence family. The fit and every
  # score must stay finite. The degenerate features are forced into the top-d
  # universe by giving them a large base value before zeroing parts of them.
  dat <- .make_lrt_vine_data(seed = 54L)
  X <- dat$X
  deg1 <- dat$block[1]
  deg2 <- dat$block[2]
  X[, deg1] <- 0                          # feature constant (zero) in BOTH classes
  X[dat$y == 1L, deg2] <- 0               # feature constant (zero) in cases only
  model <- fit_lrt_vinecopula(X, dat$y,
                              hp = list(max_anchors_per_class = 40L,
                                        n_vine_features = 8L))
  expect_true(inherits(model$repr$vine_case, "vinecop"))
  expect_true(inherits(model$repr$vine_control, "vinecop"))
  score <- score_lrt_vinecopula(model, X)
  expect_length(score, nrow(X))
  expect_true(all(is.finite(score)))
  # Scoring a specimen that is itself zero on the degenerate features is still finite.
  expect_true(all(is.finite(score_lrt_vinecopula(model, X[1:4, , drop = FALSE]))))
})

test_that("lrt-vinecopula hyperparameter validation is strict", {
  dat <- .make_lrt_vine_data(seed = 50L)
  # Unnamed and duplicate hp fields are rejected BEFORE the allowed-list check.
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(8L)),
               "hp fields must be named")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y,
                                  hp = setNames(list(8L, 9L),
                                                c("n_vine_features", "n_vine_features"))),
               "hp fields must be unique")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(bogus = 1)),
               "unknown hp")
  # n_vine_features must be an integer >= 2.
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(n_vine_features = 1L)),
               "n_vine_features")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(n_vine_features = 2.5)),
               "n_vine_features")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(n_vine_features = 3e9)),
               "n_vine_features")
  # trunc_lvl must be a positive integer.
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(trunc_lvl = 0L)),
               "trunc_lvl")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(trunc_lvl = 2.5)),
               "trunc_lvl")
  # family_set must be a non-empty character subset of accepted family names.
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(family_set = character(0))),
               "family_set")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(family_set = 5)),
               "family_set")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y,
                                  hp = list(family_set = c("gaussian", "notafamily"))),
               "unknown family name")
  # "bb" is a plausible-looking but INVALID rvinecopulib family_set token (the BB
  # copulas are bb1/bb6/bb7/bb8 and the group alias is "archimedean", not "bb");
  # rvinecopulib::vinecop() rejects "bb", so the package validator must reject it
  # FIRST with a package-side message rather than letting it fail inside vinecop().
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(family_set = "bb")),
               "unknown family name")
  # pit_clamp open-interval endpoints rejected.
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(pit_clamp = 0)),
               "pit_clamp")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(pit_clamp = 0.5)),
               "pit_clamp")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(pit_clamp = 0.7)),
               "pit_clamp")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(max_anchors_per_class = 1.5)),
               "max_anchors_per_class")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(max_anchors_per_class = 3e9)),
               "max_anchors_per_class")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(min_features = 3e9)),
               "min_features")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(eps = 0)), "hp\\$eps")
  expect_error(fit_lrt_vinecopula(dat$X, dat$y, hp = list(seed = -1L)), "hp\\$seed")
})

test_that("lrt-vinecopula score is the case-minus-control vine log-density ratio", {
  # The score equals log dvinecop(u_case, vine_case) - log dvinecop(u_ctrl,
  # vine_control), recomputed here independently from the model's stored vines.
  dat <- .make_lrt_vine_data(seed = 55L)
  model <- fit_lrt_vinecopula(dat$X, dat$y,
                              hp = list(max_anchors_per_class = 60L,
                                        n_vine_features = 8L))
  present <- model$feature_universe
  repr <- model$repr
  X_eval <- dat$X[1:6, present, drop = FALSE]
  direct <- score_lrt_vinecopula(model, dat$X[1:6, ])
  manual <- vapply(seq_len(nrow(X_eval)), function(i) {
    v <- X_eval[i, ]
    z <- rep.int(0, length(v)); pos <- which(v > 0)
    z[pos] <- log(v[pos]) - mean(log(v[pos]))
    pit <- function(marg) {
      vapply(seq_along(z), function(j) {
        k <- marg$knots[[j]]
        u <- if (length(k$x) == 1L) k$u[1L] else
          stats::approx(k$x, k$u, xout = z[[j]], rule = 2)$y
        min(max(u, model$hp$pit_clamp), 1 - model$hp$pit_clamp)
      }, numeric(1))
    }
    uc <- pit(repr$marg_case); u0 <- pit(repr$marg_control)
    log(rvinecopulib::dvinecop(matrix(uc, 1), repr$vine_case)) -
      log(rvinecopulib::dvinecop(matrix(u0, 1), repr$vine_control))
  }, numeric(1))
  expect_equal(direct, manual, tolerance = 1e-10)
})

test_that("lrt-vinecopula vine density matches the analytic Gaussian copula (external anchor)", {
  # External ground truth, independent of the scorer internals: a 2-D vine fitted
  # with family_set = "gaussian" must reproduce the ANALYTIC bivariate Gaussian
  # copula density (built from scratch with dnorm/qnorm), and an independence vine
  # must have density identically 1. This locks the correctness of the rvinecopulib
  # copula-density machinery the score relies on (the orchestrator sec.7 anchor).
  skip_if_not_installed("rvinecopulib")
  set.seed(2026L)
  n <- 1000L; rho <- 0.6
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  U <- cbind(pnorm(z1), pnorm(z2))
  fit <- rvinecopulib::vinecop(U, var_types = c("c", "c"),
                               family_set = c("gaussian", "indep"),
                               trunc_lvl = 1L, keep_data = FALSE)
  pc <- rvinecopulib::get_pair_copula(fit, tree = 1, edge = 1)
  rho_hat <- as.numeric(rvinecopulib::get_parameters(pc))
  expect_equal(rho_hat, rho, tolerance = 0.05)
  pts <- matrix(c(0.2, 0.3, 0.5, 0.5, 0.8, 0.7, 0.95, 0.1), ncol = 2, byrow = TRUE)
  q <- qnorm(pts)
  analytic <- 1 / sqrt(1 - rho_hat^2) *
    exp(-(rho_hat^2 * (q[, 1]^2 + q[, 2]^2) - 2 * rho_hat * q[, 1] * q[, 2]) /
          (2 * (1 - rho_hat^2)))
  dv <- rvinecopulib::dvinecop(pts, fit)
  expect_equal(dv, analytic, tolerance = 1e-8)
  fit_i <- rvinecopulib::vinecop(U, var_types = c("c", "c"),
                                 family_set = "indep", trunc_lvl = 1L,
                                 keep_data = FALSE)
  expect_equal(rvinecopulib::dvinecop(pts, fit_i), rep(1, nrow(pts)),
               tolerance = 1e-12)
})
