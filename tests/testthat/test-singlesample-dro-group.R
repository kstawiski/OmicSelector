library(testthat)

# Runs inside the package test harness (devtools::test / R CMD check), where
# fit_dro_group / score_dro_group, the Group-DRO primitive, and all package
# helpers are already loaded. The wrapper is pure R; no Python is required.

# Synthetic generator with a CASE-ONLY location shift in a feature block, split
# across kit x biofluid groups so the Group-DRO objective sees >1 group. The
# first `k` of `p` log-abundance features are elevated by `shift` for cases, so
# after exp() and the train-only CLR those features carry the case signal.
.make_dro_group_data <- function(n = 160L, p = 25L, k = 6L, shift = 1.2,
                                 sd = 0.5, seed = 41L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  feats <- paste0("hsa-miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), feats))
  case <- which(y == 1L)
  L[case, seq_len(k)] <- L[case, seq_len(k)] + shift
  X <- exp(L)
  # Two kits crossed with two biofluids => 4 kit x biofluid groups, balanced
  # against the case/control split so each group carries both classes. The kit
  # and biofluid factors use DIFFERENT periods (2 vs 4) so they genuinely cross
  # to all four joint groups rather than running in lockstep (only 2 groups).
  meta <- data.frame(
    kit_family = rep(c("kitA", "kitB"), times = n / 2L),
    biofluid   = rep(c("plasma", "serum"), each = 2L, length.out = n),
    stringsAsFactors = FALSE)
  list(X = X, y = y, meta = meta, block = feats[seq_len(k)])
}

.auc_mw_dro <- function(y, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE,
                                          direction = "<"))))
  }
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


test_that("dro-group fit/score roundtrip has the right shape and types", {
  dat <- .make_dro_group_data()
  model <- fit_dro_group(dat$X, dat$y, dat$meta,
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  expect_s3_class(model, "dro_group_model")
  expect_s3_class(model$fit, "group_dro_scorer")
  expect_true(is.numeric(model$n_groups))
  expect_true(is.logical(model$groups_engaged))
  expect_identical(model$kit_col, "kit_family")
  expect_identical(model$biofluid_col, "biofluid")

  score <- score_dro_group(model, dat$X)
  expect_type(score, "double")
  expect_length(score, nrow(dat$X))
  expect_true(all(is.finite(score)))
  expect_null(attributes(score))  # plain numeric vector, no diagnostic attrs
})

test_that("dro-group score equals the underlying primitive scorer (non-circularity)", {
  # The wrapper must not alter the numeric score relative to a direct call on the
  # frozen primitive fit stored in model$fit -- the canonical contract is a pure
  # signature adaptation, nothing more.
  dat <- .make_dro_group_data(seed = 202L)
  model <- fit_dro_group(dat$X, dat$y, dat$meta,
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  direct_primitive <- as.numeric(
    score_group_dro_scorer(model$fit, dat$X))
  expect_equal(score_dro_group(model, dat$X), direct_primitive,
               tolerance = 1e-12)
})

test_that("dro-group separates a planted case-vs-control shift (held-out AUC > 0.5)", {
  dat <- .make_dro_group_data(n = 300L, seed = 303L)
  ctrl <- which(dat$y == 0L); case <- which(dat$y == 1L)
  train <- c(utils::head(ctrl, 112L), utils::head(case, 112L))
  test <- c(utils::tail(ctrl, 38L), utils::tail(case, 38L))
  model <- fit_dro_group(dat$X[train, ], dat$y[train], dat$meta[train, ],
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 150L))
  score <- score_dro_group(model, dat$X[test, ])
  expect_true(all(is.finite(score)))
  expect_gt(mean(score[dat$y[test] == 1]), mean(score[dat$y[test] == 0]))
  expect_gt(.auc_mw_dro(dat$y[test], score), 0.5)
})

test_that("dro-group passes the canonical row-equivariance gate", {
  dat <- .make_dro_group_data(seed = 43L)
  model <- fit_dro_group(dat$X[1:120, ], dat$y[1:120], dat$meta[1:120, ],
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  X_test <- dat$X[121:160, ]
  score_fun <- function(model, X, meta) score_dro_group(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(score_fun, model, X_test,
                                        model_digest = .dro_group_model_digest))
  expect_true(
    singlesample_is_row_equivariant(score_fun, model, X_test,
                                    model_digest = .dro_group_model_digest))
})

test_that("dro-group matches the canonical dispatch", {
  dat <- .make_dro_group_data(seed = 46L)
  model <- fit_dro_group(dat$X[1:120, ], dat$y[1:120], dat$meta[1:120, ],
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  X_test <- dat$X[121:160, ]
  direct <- score_dro_group(model, X_test)

  # Post-integration the manifest carries the dro-group row with score_fn
  # score_dro_group in the OmicSelector namespace, so the canonical adapter
  # resolves it directly. In the pre-integration staged self-test the namespace
  # entry / repointed manifest row do not yet exist; clone a one-row roster from a
  # committed within-discriminator row and register a thin adapter so the same
  # dispatch wiring is exercised. Either way the dispatched score must equal the
  # direct call.
  roster <- singlesample_method_roster()
  resolvable <- exists("score_dro_group", envir = asNamespace("OmicSelector"),
                       mode = "function", inherits = FALSE) &&
    identical(roster$score_fn[roster$method_id == "dro-group"], "score_dro_group")
  if (!resolvable) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "dro-group-canonical-selftest"
    tmpl$fit_fn <- "fit_dro_group"
    tmpl$score_fn <- "score_dro_group"
    roster <- rbind(roster[roster$method_id != "dro-group", , drop = FALSE], tmpl)
    singlesample_register_score_adapter(
      "dro-group-canonical-selftest",
      function(model, X, meta) score_dro_group(model, X, meta))
    dispatched <- singlesample_score_call("dro-group-canonical-selftest", model,
                                          X_test, roster = roster)
  } else {
    dispatched <- singlesample_score_call("dro-group", model, X_test,
                                          roster = roster)
  }
  expect_equal(dispatched, direct, tolerance = 1e-12)
})

test_that("dro-group degenerates to a single group when meta is absent or lacks cols", {
  dat <- .make_dro_group_data(seed = 44L)

  # Real kit x biofluid meta engages multiple groups.
  m_real <- fit_dro_group(dat$X, dat$y, dat$meta,
                          hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  expect_gt(m_real$n_groups, 1L)
  expect_true(m_real$groups_engaged)
  expect_true(m_real$real_meta)

  # NULL meta -> single constant group (ERM); still a valid single-sample model.
  m_null <- fit_dro_group(dat$X, dat$y, meta_train = NULL,
                          hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  expect_equal(m_null$n_groups, 1L)
  expect_false(m_null$groups_engaged)
  expect_false(m_null$real_meta)
  s_null <- score_dro_group(m_null, dat$X)
  expect_length(s_null, nrow(dat$X))
  expect_true(all(is.finite(s_null)))

  # Meta present but missing the configured columns -> degenerate single group.
  bad_meta <- data.frame(other = seq_len(nrow(dat$X)), stringsAsFactors = FALSE)
  m_bad <- fit_dro_group(dat$X, dat$y, bad_meta,
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  expect_equal(m_bad$n_groups, 1L)
  expect_false(m_bad$groups_engaged)
  expect_false(m_bad$real_meta)
})

test_that("dro-group score orientation is larger = more case-like", {
  dat <- .make_dro_group_data(seed = 51L)
  model <- fit_dro_group(dat$X, dat$y, dat$meta,
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 150L))
  score <- score_dro_group(model, dat$X)
  expect_gt(mean(score[dat$y == 1]), mean(score[dat$y == 0]))
  expect_gt(.auc_mw_dro(dat$y, score), 0.5)
})

test_that("dro-group meta at score time is diagnostic-only (never changes the score)", {
  dat <- .make_dro_group_data(seed = 52L)
  model <- fit_dro_group(dat$X, dat$y, dat$meta,
                         hp = list(tune = FALSE, panel_size = 15L, epochs = 120L))
  X_test <- dat$X[1:20, ]
  s_no_meta <- score_dro_group(model, X_test)
  s_with_meta <- score_dro_group(model, X_test, meta = dat$meta[1:20, ])
  expect_equal(s_no_meta, s_with_meta, tolerance = 1e-12)
})

test_that("dro-group fitting is deterministic and RNG-safe; scoring never touches RNG", {
  dat <- .make_dro_group_data(seed = 48L)
  hp <- list(tune = FALSE, panel_size = 15L, epochs = 120L, seed = 7L)
  set.seed(999)
  seed_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  model_1 <- fit_dro_group(dat$X, dat$y, dat$meta, hp = hp)
  # Fit restores the global RNG it temporarily perturbs.
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   seed_before)
  model_2 <- fit_dro_group(dat$X, dat$y, dat$meta, hp = hp)
  expect_identical(model_1$fit$beta, model_2$fit$beta)
  expect_identical(model_1$fit$intercept, model_2$fit$intercept)
  expect_identical(.dro_group_model_digest(model_1),
                   .dro_group_model_digest(model_2))

  s1 <- score_dro_group(model_1, dat$X[1:20, ])
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  s2 <- score_dro_group(model_1, dat$X[1:20, ])
  expect_identical(s1, s2)
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)
})

test_that("dro-group hyperparameter validation is strict", {
  dat <- .make_dro_group_data(seed = 50L)
  base_ok <- list(tune = FALSE, panel_size = 15L, epochs = 60L)
  # Unknown / duplicate / unnamed hp fields are rejected.
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(bogus = 1)),
               "unknown hp")
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(15L)),
               "must be named")
  expect_error(
    fit_dro_group(dat$X, dat$y, dat$meta,
                  hp = structure(list(20L, 20L), names = c("panel_size", "panel_size"))),
    "duplicate hp")
  # Type / range violations per field.
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(panel_size = 0)),
               "panel_size")
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(panel_size = 1.5)),
               "panel_size")
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(kit_col = 1)),
               "kit_col")
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(eta_grid = c(0.1, -1))),
               "eta_grid")
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(learning_rate = 0)),
               "learning_rate")
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(tune = NA)),
               "tune")
  expect_error(fit_dro_group(dat$X, dat$y, dat$meta, hp = list(seed = -1)),
               "seed")
})
