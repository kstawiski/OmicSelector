# icp implements INVARIANT CAUSAL PREDICTION (Peters et al. 2016): base-R, no
# torch/reticulate. It is the structural twin of sel-stablemate with the
# StableMate selection swapped for an explicit per-feature ICP screen
# (predictivity Wald p < alpha_pred AND Cochran's-Q invariance p > alpha_inv).

# Planted multi-environment data. Stable features (miR-1, miR-2, miR-3) are
# elevated in cases CONSISTENTLY across all environments (multiplicative shift,
# so the signal lives in the per-sample rCLR AND the y|x_j slope is invariant
# across cohorts). Spurious features (miR-4, miR-5, miR-6) are elevated in cases
# in ONE environment each only -- predictive within that cohort but NON-invariant
# across cohorts. The remaining features are noise.
.make_icp_data <- function(seed = 3L, n_per = 70L,
                           envs = c("e1", "e2", "e3"),
                           stable_mult = 2.2, spur_mult = 3.0) {
  set.seed(seed)
  p <- 12L
  n <- n_per * length(envs)
  env <- rep(envs, each = n_per)
  y <- rep(rep(c(0, 1), each = n_per / 2L), length(envs))
  X <- matrix(stats::rgamma(n * p, shape = 8, rate = 1), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("miR-", seq_len(p))))
  for (f in 1:3) X[y == 1L, f] <- X[y == 1L, f] * stable_mult
  spur_idx <- c(4L, 5L, 6L)
  for (j in seq_along(envs)) {
    if (j > length(spur_idx)) break
    rows <- env == envs[j] & y == 1L
    X[rows, spur_idx[j]] <- X[rows, spur_idx[j]] * spur_mult
  }
  list(X = X, y = y, env = env,
       meta = data.frame(accession = env, stringsAsFactors = FALSE),
       stable = paste0("miR-", 1:3), spurious = paste0("miR-", 4:6))
}

.auc_mw <- function(y, s) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, s, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(s, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# §7.1 internal-consistency / forward shape
test_that("icp fit/score returns finite numeric vector of right shape", {
  dat <- .make_icp_data()
  model <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  score <- score_icp(model, dat$X)

  expect_s3_class(model, "icp_model")
  expect_type(model$selected_features, "character")
  expect_type(model$coefficients, "double")
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
  expect_true(length(model$selected_features) >= 1L)  # always non-empty head
})

# §7.11 ICP-SELECTION correctness (NON-tautological): a hand-built 2-cohort toy
# with (a) one truly INVARIANT predictive feature [same y|x_j slope in both
# cohorts] -> KEPT, (b) one NON-invariant feature [predictive in cohort 1, null
# in cohort 2, so its y|x_j slope DIFFERS across cohorts] -> DROPPED *by the
# invariance screen* (it is predictive, but the Cochran's-Q homogeneity test
# rejects), (c) one null feature -> DROPPED *by the predictivity screen*.
#
# The toy uses MANY dilution-noise features so no single planted feature
# dominates the per-sample rCLR geometric-mean closure; the planted case-shifts
# are MODERATE for the same reason. We assert the three brief-mandated, decisive
# claims directly on the per-feature screen verdicts (so the assertions are
# attributable to the RIGHT screen, not just to set membership), then on the
# selected set. We do NOT assert that every dilution-noise feature is excluded:
# under compositional closure a planted case-shift induces a small, cross-cohort-
# CONSISTENT anti-correlation in the other parts, so an occasional noise feature
# is genuinely (closure-)predictive-and-invariant -- that is a property of rCLR
# closure, not a selection defect. The load-bearing claims are inv/varying/null.
test_that("icp selection keeps the invariant predictor and drops the cohort-varying (invariance) and null (predictivity) ones", {
  set.seed(42L)
  n_per <- 200L
  feat <- c("inv", "varying", "null", paste0("noise", sprintf("%02d", 1:30)))
  p <- length(feat)
  mk <- function(env_label, varying_mult) {
    X <- matrix(stats::rlnorm(n_per * p, meanlog = log(50), sdlog = 0.45),
                n_per, p, dimnames = list(NULL, feat))
    y <- rep(c(0L, 1L), each = n_per / 2L)
    # (a) invariant predictive: SAME multiplicative case-shift in BOTH cohorts.
    X[y == 1L, "inv"] <- X[y == 1L, "inv"] * 1.8
    # (b) cohort-varying: case-shift in e1 (varying_mult = 2.2), NONE in e2
    #     (varying_mult = 1.0) -> the y|x slope differs across cohorts.
    X[y == 1L, "varying"] <- X[y == 1L, "varying"] * varying_mult
    # (c) "null" + noise: no planted case association in either cohort.
    list(X = X, y = y, env = rep(env_label, n_per))
  }
  e1 <- mk("e1", varying_mult = 2.2)
  e2 <- mk("e2", varying_mult = 1.0)
  X <- rbind(e1$X, e2$X)
  y <- c(e1$y, e2$y)
  meta <- data.frame(accession = c(e1$env, e2$env), stringsAsFactors = FALSE)

  hp <- list(max_features = p, alpha_pred = 0.05, alpha_inv = 0.10,
             min_env_n = 8L)
  model <- fit_icp(X, y, meta_train = meta, hp = hp)

  # Per-feature screen verdicts (recomputed exactly as the fitter does) so each
  # claim is attributable to the SPECIFIC screen the brief names.
  universe <- .icp_universe(X, p)
  Z <- .icp_rclr_matrix(X[, universe, drop = FALSE])
  env_info <- .icp_environment(meta, y, "accession", hp$min_env_n)
  pp <- function(f) .icp_pred_pvalue(Z[, f], y)
  qp <- function(f) .icp_invariance_pvalue(Z[, f], y, env_info$env,
                                           env_info$valid_levels)

  # (a) inv: PREDICTIVE and INVARIANT -> the conjunction holds -> KEPT.
  expect_lt(pp("inv"), hp$alpha_pred)
  expect_gt(qp("inv"), hp$alpha_inv)
  # (b) varying: PREDICTIVE but NON-invariant -> dropped specifically by the
  #     invariance (Cochran's-Q) screen, NOT by predictivity.
  expect_lt(pp("varying"), hp$alpha_pred)        # it IS predictive
  expect_lt(qp("varying"), hp$alpha_inv)         # but the Q-test rejects homogeneity
  # (c) null: NOT predictive -> dropped by the predictivity screen.
  expect_gte(pp("null"), hp$alpha_pred)

  expect_equal(model$n_environments, 2L)                  # invariance screen ran
  expect_equal(model$selection_mode, "invariant_predictive")
  expect_true("inv" %in% model$selected_features)         # (a) invariant -> KEPT
  expect_false("varying" %in% model$selected_features)    # (b) varies -> DROPPED
  expect_false("null" %in% model$selected_features)       # (c) null -> DROPPED
})

# §7.5 held-out-cohort AUC > 0.7 on a planted shift: train e1+e2, hold out e3.
test_that("icp selects invariant predictors and generalizes to a held-out environment", {
  tr <- .make_icp_data(seed = 3L, envs = c("e1", "e2"))
  te <- .make_icp_data(seed = 99L, envs = c("e3"))
  model <- fit_icp(tr$X, tr$y, meta_train = tr$meta)

  n_stable <- length(intersect(model$selected_features, tr$stable))
  n_spur <- length(intersect(model$selected_features, tr$spurious))
  expect_equal(model$n_environments, 2L)
  expect_gte(n_stable, 2L)            # recovers the invariant (stable) set
  expect_gt(n_stable, n_spur)         # invariant overlap dominates spurious

  score_te <- score_icp(model, te$X)
  expect_gt(mean(score_te[te$y == 1L]), mean(score_te[te$y == 0L]))
  expect_gt(.auc_mw(te$y, score_te), 0.7)   # generalizes to held-out env
})

# §7.10 independent recompute of the linear score from the exported coefficients.
test_that("icp score equals the frozen logistic linear predictor over selected rCLR features", {
  dat <- .make_icp_data(seed = 5L)
  model <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  s <- score_icp(model, dat$X)

  rclr <- function(v) {
    z <- rep(0, length(v)); pos <- which(v > 0)
    if (length(pos)) { lv <- log(v[pos]); z[pos] <- lv - mean(lv) }
    z
  }
  manual <- vapply(seq_len(nrow(dat$X)), function(i) {
    xi <- dat$X[i, ]; up <- intersect(model$feature_universe, names(xi))
    v <- xi[up]; z <- rclr(v); names(z) <- up
    use <- intersect(model$selected_features, up)
    if (length(use) < model$min_selected) return(0)
    model$intercept + sum(model$coefficients[use] * z[use])
  }, numeric(1))
  expect_equal(s, manual, tolerance = 1e-12)
})

# §7.4 per-sample scale-invariance (rCLR).
test_that("icp is exactly invariant to per-sample positive scaling", {
  dat <- .make_icp_data(seed = 11L)
  model <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  X_test <- dat$X[1:60, , drop = FALSE]
  s_base <- score_icp(model, X_test)

  expect_equal(score_icp(model, X_test * 13.7), s_base, tolerance = 1e-6)
  set.seed(5)
  scal <- stats::runif(nrow(X_test), 0.05, 40)
  expect_equal(score_icp(model, X_test * scal), s_base, tolerance = 1e-6)
})

# §7.6 canonical row-equivariance gate + singlesample_score_call dispatch.
test_that("icp passes the canonical row-equivariance gate", {
  dat <- .make_icp_data(seed = 12L)
  model <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  X_test <- dat$X[1:50, , drop = FALSE]
  score_fun <- function(model, X, meta) score_icp(model, X, meta)

  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X_test))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X_test))
})

# §7.2 single==batch EXACTLY 0 ; §7.3 row-perm 0.
test_that("icp is column-permutation invariant and single-row == batch", {
  dat <- .make_icp_data(seed = 13L)
  model <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  X_test <- dat$X[1:50, , drop = FALSE]
  s_base <- score_icp(model, X_test)

  set.seed(9)
  perm <- sample(ncol(X_test))
  expect_equal(score_icp(model, X_test[, perm, drop = FALSE]), s_base,
               tolerance = 1e-10)
  expect_equal(score_icp(model, X_test[7, , drop = FALSE]),
               s_base[7], tolerance = 1e-12)
})

# §7.12 degeneration: NULL / single-cohort / missing-column meta ->
# predictivity-only, n_environments == 1, fit succeeds + finite scores + a
# usable non-empty head.
test_that("icp reports n_environments and degenerates gracefully to predictivity-only", {
  dat <- .make_icp_data(seed = 7L)

  # multi-environment -> ICP invariance screen engaged
  m_multi <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  expect_equal(m_multi$n_environments, 3L)
  expect_equal(m_multi$selection_mode, "invariant_predictive")

  # NULL meta -> predictivity-only selection, n_environments == 1
  m_null <- fit_icp(dat$X, dat$y, meta_train = NULL)
  expect_equal(m_null$n_environments, 1L)
  expect_equal(m_null$selection_mode, "predictivity_only")
  expect_gte(length(intersect(m_null$selected_features, dat$stable)), 2L)
  expect_true(all(is.finite(score_icp(m_null, dat$X))))
  expect_true(length(m_null$selected_features) >= 1L)

  # single cohort -> degenerate to predictivity-only, n_environments == 1
  one_meta <- data.frame(accession = rep("one", nrow(dat$X)),
                         stringsAsFactors = FALSE)
  m_one <- fit_icp(dat$X, dat$y, meta_train = one_meta)
  expect_equal(m_one$n_environments, 1L)
  expect_equal(m_one$selection_mode, "predictivity_only")
  expect_true(length(m_one$selected_features) >= 1L)

  # missing cohort column -> predictivity-only, n_environments == 1
  m_miss <- fit_icp(dat$X, dat$y, meta_train = data.frame(other = dat$env))
  expect_equal(m_miss$n_environments, 1L)
  expect_equal(m_miss$selection_mode, "predictivity_only")
})

# §7.7 degenerate -> 0 ; FLAT positive composition scored not floored (== manual).
test_that("icp returns neutral 0 below min_selected and scores a flat composition not floored", {
  dat <- .make_icp_data(seed = 15L)
  model <- fit_icp(dat$X, dat$y, meta_train = dat$meta)

  # no selected features present at all -> neutral 0
  X_none <- matrix(1, nrow = 4L, ncol = 2L,
                   dimnames = list(NULL, c("zzz-a", "zzz-b")))
  expect_equal(score_icp(model, X_none), rep(0, 4L))

  # below min_selected (require 2, present only 1) -> neutral 0
  keep1 <- model$selected_features[1]
  X1 <- dat$X[1:5, keep1, drop = FALSE]
  model_hi <- fit_icp(dat$X, dat$y, meta_train = dat$meta,
                      hp = list(min_selected = 2L))
  expect_equal(score_icp(model_hi, X1), rep(0, 5L))

  # degenerate (empty) selection -> all 0
  model_degen <- model
  model_degen$degenerate <- TRUE
  expect_equal(score_icp(model_degen, dat$X[1:10, ]), rep(0, 10L))

  # rCLR-ORIGIN TRAP: an all-zero specimen row over the universe has NO positive
  # support, so its rCLR is the flat (all-zero) vector. The score is therefore
  # the intercept-plus-zero linear predictor (the flat composition is SCORED, not
  # floored to a special neutral) whenever >= min_selected selected features are
  # present. It must be finite and EQUAL the intercept (not 0 unless intercept 0).
  X_flat <- dat$X[1:3, , drop = FALSE]
  X_flat[1, ] <- 0
  s_flat <- score_icp(model, X_flat)
  expect_true(all(is.finite(s_flat)))
  # row 1: every present selected feature has rCLR coord 0 -> score == intercept.
  present_sel <- intersect(model$selected_features, colnames(X_flat))
  if (length(present_sel) >= model$min_selected) {
    expect_equal(s_flat[1], model$intercept, tolerance = 1e-12)
  }
})

# §7.8 determinism + §7.9 RNG-safety: two seed-matched fits -> identical digest;
# fit leaves global .Random.seed byte-unchanged; score consumes no RNG.
test_that("icp refit is deterministic and leaves the global RNG untouched", {
  dat <- .make_icp_data(seed = 16L)

  set.seed(2024)
  before <- get(".Random.seed", envir = globalenv())
  m1 <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  after_fit <- get(".Random.seed", envir = globalenv())
  s1 <- score_icp(m1, dat$X[1:20, ])
  after_score <- get(".Random.seed", envir = globalenv())
  m2 <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  s2 <- score_icp(m2, dat$X[1:20, ])

  expect_identical(m1$selected_features, m2$selected_features)
  expect_identical(m1$selection_mode, m2$selection_mode)
  expect_identical(m1$intercept, m2$intercept)
  expect_identical(m1$coefficients, m2$coefficients)
  expect_identical(s1, s2)
  expect_identical(before, after_fit)      # RNG untouched by fit
  expect_identical(after_fit, after_score) # RNG untouched by score

  # fit with NO pre-existing seed leaves none
  rm(list = ".Random.seed", envir = globalenv())
  m3 <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m3, "icp_model")
})

# §7.6 dispatch: deployable through the canonical singlesample_score_call.
test_that("icp is deployable through the canonical dispatch", {
  dat <- .make_icp_data(seed = 17L)
  model <- fit_icp(dat$X, dat$y, meta_train = dat$meta)
  X_test <- dat$X[1:40, , drop = FALSE]

  # Self-contained one-row roster so singlesample_score_call resolves icp without
  # depending on manifest integration order. The dispatched score must equal the
  # direct call. (pkg_status "new" so the canonical namespace adapter is used.)
  roster <- data.frame(
    method_id = "icp", family = "H", estimand = "transfer",
    role = "discriminator", tier = "R2", dep_route = "base-r",
    fit_fn = "fit_icp", score_fn = "score_icp",
    pkg_status = "new", notes = "", row_source = "per_cohort_rows",
    lopbo_mechanism = "within_cohort", stringsAsFactors = FALSE)

  # --- STAGED ISOLATION ANCHOR (stripped at integration) ------------------
  # Pre-integration score_icp is NOT yet in the OmicSelector namespace, so the
  # canonical adapter cannot resolve it. Register a transient adapter so the
  # dispatch path is exercised against the STAGED scorer. At integration
  # score_icp lives in R/ and resolves natively (pkg_status "new" -> canonical
  # adapter), so this registration block is removed.
  if (!exists("score_icp", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "icp", function(model, X, meta) score_icp(model, X, meta))
    on.exit(rm("icp", envir = OmicSelector:::.singlesample_score_adapter_registry),
            add = TRUE)
  }
  # ------------------------------------------------------------------------

  expect_equal(
    singlesample_score_call("icp", model, X_test, roster = roster),
    score_icp(model, X_test),
    tolerance = 1e-12)
})

test_that("icp strict hp resolver rejects malformed hp", {
  dat <- .make_icp_data(seed = 18L)

  expect_error(fit_icp(dat$X, dat$y, hp = "notalist"), "hp must be a list")
  expect_error(fit_icp(dat$X, dat$y, hp = list(foo = 5L)), "unknown hp")
  expect_error(fit_icp(dat$X, dat$y, hp = list(10L)), "must be named")
  expect_error(
    fit_icp(dat$X, dat$y,
            hp = stats::setNames(list(5L, 6L), c("top_k", "top_k"))),
    "duplicate hp")
  expect_error(fit_icp(dat$X, dat$y, hp = list(max_features = 1L)),
               "max_features")
  expect_error(fit_icp(dat$X, dat$y, hp = list(alpha_pred = 0)), "alpha_pred")
  expect_error(fit_icp(dat$X, dat$y, hp = list(alpha_pred = 1)), "alpha_pred")
  expect_error(fit_icp(dat$X, dat$y, hp = list(alpha_inv = 0)), "alpha_inv")
  expect_error(fit_icp(dat$X, dat$y, hp = list(alpha_inv = 1)), "alpha_inv")
  expect_error(fit_icp(dat$X, dat$y, hp = list(min_env_n = 3L)), "min_env_n")
  expect_error(fit_icp(dat$X, dat$y, hp = list(top_k = 0L)), "top_k")
  expect_error(fit_icp(dat$X, dat$y, hp = list(min_selected = 0L)),
               "min_selected")
  expect_error(fit_icp(dat$X, dat$y, hp = list(cohort_col = 1)), "cohort_col")
  expect_error(fit_icp(dat$X, dat$y, hp = list(cohort_col = "")), "cohort_col")
  expect_error(fit_icp(dat$X, dat$y, hp = list(seed = -1L)), "seed")
})

test_that("icp input validation errors are explicit", {
  dat <- .make_icp_data(seed = 19L)

  X_unnamed <- dat$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_icp(X_unnamed, dat$y, dat$meta), "feature names")
  expect_error(fit_icp(dat$X[, 1, drop = FALSE], dat$y, dat$meta),
               "at least two features")
  expect_error(fit_icp(dat$X, dat$y[-1], dat$meta), "length\\(y_train\\)")
  expect_error(fit_icp(dat$X, dat$y, dat$meta[-1, , drop = FALSE]),
               "one row per row")
  expect_error(fit_icp(dat$X, rep(0, nrow(dat$X)), dat$meta),
               "at least one case")
  expect_error(fit_icp(dat$X, rep(1, nrow(dat$X)), dat$meta),
               "at least one control")
  expect_error(score_icp(list(), dat$X), "class icp_model")
})
