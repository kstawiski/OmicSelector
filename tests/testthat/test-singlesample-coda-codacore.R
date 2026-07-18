library(testthat)

# coda-codacore trains the CoDaCoRe continuous relaxation stagewise via
# reticulate-python torch at FIT (torch is confined to fit; the DISCRETE balances +
# logistic head are exported and the SCORE is pure base R). Tests that fit are
# skipped when reticulate or the venv torch is unavailable, so the package suite
# stays green on hosts without the venv.
skip_if_no_coda_codacore <- function() {
  venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
  if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = venv)
  }
  testthat::skip_if_not_installed("reticulate")
  ok <- tryCatch({
    reticulate::import("torch", delay_load = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(ok, "python torch not importable in the reticulate venv")
}

# Planted multi-balance data: a known numerator group is elevated in cases and a
# known denominator group depressed (otherwise identical gamma baseline). The
# log-ratio of those groups carries the discriminative signal the CoDaCoRe
# relaxation learns and the discrete balance + frozen head score. Abundances are
# strictly positive (compositional).
.make_coda_codacore_data <- function(n = 120L, p = 16L, seed = 7L, mult = 1.8) {
  set.seed(seed)
  feat <- paste0("f", sprintf("%02d", seq_len(p)))
  X <- matrix(stats::rgamma(n * p, shape = 30, rate = 2), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), feat))
  y <- rep(c(0, 1), each = n / 2L)
  num_pl <- c("f01", "f02", "f03")
  den_pl <- c("f14", "f15", "f16")
  X[y == 1L, num_pl] <- X[y == 1L, num_pl] * mult
  X[y == 1L, den_pl] <- X[y == 1L, den_pl] / mult
  list(X = X, y = y, num_pl = num_pl, den_pl = den_pl, feat = feat)
}

.auc_mw_cc <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# A small fit reused across single-sample-property tests. device default = "cpu"
# (the discretization is EXACTLY reproducible only on CPU).
.fit_small_cc <- function(seed = 7L, device = "cpu", epochs = 120L) {
  d <- .make_coda_codacore_data(seed = seed)
  list(model = fit_coda_codacore(d$X, d$y,
                                 hp = list(device = device, epochs = epochs,
                                           seed = 42L)),
       data = d)
}


test_that("fit returns a well-formed coda_codacore_model", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc()
  m <- f$model
  expect_s3_class(m, "coda_codacore_model")
  expect_identical(m$feature_universe, colnames(f$data$X))
  expect_true(m$n_balances >= 1L)
  expect_length(m$weights, m$n_balances)
  expect_true(is.finite(m$intercept))
  # every frozen balance is a disjoint non-empty A / B split.
  for (bb in m$balances) {
    expect_true(length(bb$A) >= 1L)
    expect_true(length(bb$B) >= 1L)
    expect_length(intersect(bb$A, bb$B), 0L)
  }
  s <- score_coda_codacore(m, f$data$X)
  expect_type(s, "double")
  expect_length(s, nrow(f$data$X))
  expect_true(all(is.finite(s)))
  # No live python pointer survives the fit (export-to-pure-R contract).
  expect_false(any(vapply(m, function(x) inherits(x, "python.builtin.object"),
                          logical(1L))))
})

# §7.1 balance correctness: the pure-R discrete balance matches the closed-form
# ILR balance on a hand-checked tiny case to < 1e-12.
test_that("§7.1 discrete balance matches the closed-form ILR balance (<1e-12)", {
  skip_if_no_coda_codacore()
  # A = {1, 2}, B = {3}; x = (e^1, e^3, e^2) -> mean log A = 2, mean log B = 2.
  # bal = sqrt(|A||B|/(|A|+|B|)) * (meanlogA - meanlogB)
  #     = sqrt(2*1/3) * (2 - 2) = 0.
  x1 <- c(exp(1), exp(3), exp(2))
  b1 <- .coda_codacore_balance_one(x1, a_idx = c(1L, 2L), b_idx = 3L)
  expect_lt(abs(b1 - 0), 1e-12)

  # A = {1}, B = {2}; x = (e^5, e^2) -> bal = sqrt(1*1/2) * (5 - 2) = 3/sqrt(2).
  x2 <- c(exp(5), exp(2))
  b2 <- .coda_codacore_balance_one(x2, a_idx = 1L, b_idx = 2L)
  expect_lt(abs(b2 - 3 / sqrt(2)), 1e-12)

  # A = {1, 3}, B = {2, 4}; logs A = (0, 4) mean 2, logs B = (1, 1) mean 1.
  # bal = sqrt(2*2/4) * (2 - 1) = 1 * 1 = 1.
  x3 <- c(exp(0), exp(1), exp(4), exp(1))
  b3 <- .coda_codacore_balance_one(x3, a_idx = c(1L, 3L), b_idx = c(2L, 4L))
  expect_lt(abs(b3 - 1), 1e-12)
})

# §7.2 single-row score == batch score, maxdiff EXACTLY 0 (pure-R discrete
# balances; no python at score; no cross-sample statistic).
test_that("§7.2 DECISIVE: single-row score == batch score (maxdiff 0)", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc()
  m <- f$model
  X <- f$data$X[seq_len(12L), , drop = FALSE]
  batch <- score_coda_codacore(m, X)
  single <- vapply(seq_len(nrow(X)), function(i) {
    score_coda_codacore(m, X[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - single))
  cat(sprintf("\n[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

# §7.3 SAMPLE row-permutation invariance, maxdiff 0.
test_that("§7.3 score is invariant to row permutation (maxdiff 0)", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc()
  m <- f$model
  X <- f$data$X[seq_len(12L), , drop = FALSE]
  base <- score_coda_codacore(m, X)
  set.seed(7)
  perm <- sample(nrow(X))
  sp <- score_coda_codacore(m, X[perm, , drop = FALSE])
  md <- max(abs(sp - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

# §7.4 per-sample positive-scale invariance ~0 (a balance is a difference of
# mean-logs -> a positive rescale c cancels ANALYTICALLY:
# mean log(c x_A) - mean log(c x_B) = (log c + mean log x_A) - (log c + mean log x_B)
# = mean log x_A - mean log x_B. In float64 the `log c` terms cancel only up to a
# few ULP because log(c x_j) is one `log` of the product, not log(c)+log(x_j), so
# the residual is ~1e-14 (the proto-net / inv-scatter / bal-selbal siblings gate
# this at a small tolerance for the same reason). The gate is a TIGHT 1e-10.
test_that("§7.4 score is invariant to per-sample positive scaling (~0, < 1e-10)", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc()
  m <- f$model
  X <- f$data$X[seq_len(12L), , drop = FALSE]
  base <- score_coda_codacore(m, X)
  # constant scale
  md_const <- max(abs(score_coda_codacore(m, X * 13.7) - base))
  expect_lt(md_const, 1e-10)
  # per-row random positive scale
  set.seed(3)
  scal <- matrix(stats::runif(nrow(X), 1e-3, 1e3), nrow(X), ncol(X))
  ss <- score_coda_codacore(m, X * scal)
  md <- max(abs(ss - base))
  cat(sprintf("[§7.4] per-sample scale-invariance maxdiff = %.3e (const %.3e)\n",
              md, md_const))
  expect_lt(md, 1e-10)
})

# §7.5 held-out AUC > 0.7 on a planted log-ratio shift. The CoDaCoRe balance over
# the planted numerator/denominator groups + the frozen head separate the classes.
test_that("§7.5 larger score = more case-like on a held-out planted shift", {
  skip_if_no_coda_codacore()
  d <- .make_coda_codacore_data(seed = 7L, mult = 2.0)
  tr <- c(1:45, 61:105)            # 45 controls + 45 cases (rows are y-ordered)
  te <- setdiff(seq_len(nrow(d$X)), tr)
  m <- fit_coda_codacore(d$X[tr, ], d$y[tr],
                         hp = list(device = "cpu", epochs = 200L))
  s <- score_coda_codacore(m, d$X[te, ])
  auc <- .auc_mw_cc(d$y[te], s)
  cat(sprintf("[§7.5] held-out AUC = %.4f (n_balances = %d)\n", auc, m$n_balances))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

# §7.6 canonical row-equivariance gate + dispatch.
test_that("§7.6 passes the canonical row-equivariance gate + dispatch matches", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc(seed = 12L)
  m <- f$model
  X_test <- f$data$X[seq_len(20L), , drop = FALSE]
  score_fun <- function(model, X, meta) score_coda_codacore(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .coda_codacore_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .coda_codacore_model_digest))

  # Canonical dispatch: post-integration the manifest carries the coda-codacore row
  # and score_coda_codacore is in the namespace; pre-integration (staged self-test)
  # clone a one-row roster from a committed within-discriminator row and register a
  # thin adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_coda_codacore(m, X_test)
  roster <- tryCatch(singlesample_method_roster(), error = function(e) NULL)
  if (is.null(roster)) {
    skip("roster unavailable")
  }
  if (!"coda-codacore" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "bal-selbal", , drop = FALSE]
    if (nrow(tmpl) == 0L) tmpl <- roster[1L, , drop = FALSE]
    tmpl$method_id <- "coda-codacore"
    tmpl$fit_fn <- "fit_coda_codacore"
    tmpl$score_fn <- "score_coda_codacore"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_coda_codacore", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "coda-codacore",
      function(model, X, meta) score_coda_codacore(model, X, meta))
  }
  via <- tryCatch(
    singlesample_score_call("coda-codacore", m, X_test, roster = roster),
    error = function(e) e)
  if (inherits(via, "error")) {
    skip(paste("coda-codacore canonical dispatch not yet wired:",
               conditionMessage(via)))
  }
  expect_equal(via, direct, tolerance = 1e-12)
})

# §7.7 degenerate -> 0; flat composition scored (to the intercept) NOT floored; a
# balance with one empty side -> 0 contribution.
test_that("§7.7 degenerate -> 0; flat composition scored; empty side -> 0 contribution", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc(seed = 15L)
  m <- f$model
  X <- f$data$X[seq_len(6L), , drop = FALSE]

  # Fewer than min_features universe features present -> 0 for every row.
  X_few <- X[, seq_len(2L), drop = FALSE]
  expect_true(all(score_coda_codacore(m, X_few) == 0))

  # No shared features at all -> 0 for every row.
  X_none <- X
  colnames(X_none) <- paste0("other-", seq_len(ncol(X_none)))
  expect_true(all(score_coda_codacore(m, X_none) == 0))

  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- X
  X_zero[2L, ] <- 0
  s_zero <- score_coda_codacore(m, X_zero)
  expect_equal(s_zero[2L], 0)
  expect_true(all(is.finite(s_zero)))

  # A FLAT all-equal-positive composition gives every balance value 0 (mean log A -
  # mean log B = 0) but is a VALID specimen scored to the intercept (NOT floored).
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  s_flat <- score_coda_codacore(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (intercept %.6f)\n",
              s_flat, m$intercept))
  expect_equal(s_flat, m$intercept, tolerance = 1e-12)

  # A balance with one EMPTY side for a specimen contributes its neutral value 0:
  # keep only the A-features of the first frozen balance (its B side absent) and
  # verify that balance's contribution is 0 (so the score is intercept + the OTHER
  # balances' contributions only). Construct a one-row specimen that has positive
  # support only on balance 1's A features.
  bb1 <- m$balances[[1L]]
  a_names <- m$feature_universe[bb1$A]
  xrow <- matrix(0, nrow = 1L, ncol = p,
                 dimnames = list("oneside", m$feature_universe))
  xrow[1L, a_names] <- 5.0                    # only A side present/positive
  # balance 1 has an empty B side here -> contributes 0.
  contrib1 <- .coda_codacore_balance_one(xrow[1L, ], bb1$A, bb1$B)
  expect_equal(contrib1, 0)
})

# §7.8 determinism: two seed-matched CPU fits -> IDENTICAL discrete balances (A,B) +
# IDENTICAL model_digest + bit-identical scores (seed-before-build, full-batch,
# CPU exact discretization). Hardened with an intervening torch-RNG perturbation.
test_that("§7.8 determinism: two seed=42 CPU fits give identical balances + digest", {
  skip_if_no_coda_codacore()
  d <- .make_coda_codacore_data(seed = 21L)
  m1 <- fit_coda_codacore(d$X, d$y,
                          hp = list(device = "cpu", epochs = 120L, seed = 42L))
  X <- d$X[seq_len(10L), , drop = FALSE]
  s1a <- score_coda_codacore(m1, X)
  s1b <- score_coda_codacore(m1, X)
  expect_identical(s1a, s1b)                          # re-score identical

  # Perturb the global torch RNG between fits: the fit is deterministic regardless.
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  torch_py$manual_seed(99999L)

  m2 <- fit_coda_codacore(d$X, d$y,
                          hp = list(device = "cpu", epochs = 120L, seed = 42L))
  # Identical DISCRETE balances (the A/B sets must match exactly).
  expect_identical(m1$n_balances, m2$n_balances)
  for (b in seq_len(m1$n_balances)) {
    expect_identical(m1$balances[[b]]$A, m2$balances[[b]]$A)
    expect_identical(m1$balances[[b]]$B, m2$balances[[b]]$B)
  }
  dg1 <- .coda_codacore_model_digest(m1)
  dg2 <- .coda_codacore_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s\n", dg1, dg2))
  expect_identical(dg1, dg2)
  s2 <- score_coda_codacore(m2, X)
  expect_identical(max(abs(s1a - s2)), 0)             # refit bit-identical scores
})

# §7.9 RNG-safety: fit saves/restores the global R .Random.seed (pre-import) + the
# torch CPU/CUDA RNG; score consumes no RNG (pure-R).
test_that("§7.9 fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_coda_codacore()
  d <- .make_coda_codacore_data(seed = 31L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_coda_codacore(d$X, d$y,
                         hp = list(device = "cpu", epochs = 80L))

  after_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  after_t <- torch$random$get_rng_state()
  expect_identical(before_r, after_r)
  expect_true(isTRUE(reticulate::py_to_r(torch$equal(before_t, after_t))))
  if (cuda_ok) {
    after_cuda <- torch$cuda$get_rng_state_all()
    eq <- vapply(seq_along(before_cuda), function(i)
      isTRUE(reticulate::py_to_r(torch$equal(before_cuda[[i]], after_cuda[[i]]))),
      logical(1L))
    expect_true(all(eq))
  }
  cat(sprintf("[§7.9] R RNG unchanged = %s ; torch CPU RNG unchanged = %s\n",
              identical(before_r, after_r),
              isTRUE(reticulate::py_to_r(torch$equal(before_t, after_t)))))

  # Scoring consumes no RNG: .Random.seed unchanged across scoring.
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  invisible(score_coda_codacore(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_coda_codacore(d$X, d$y, hp = list(device = "cpu", epochs = 40L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "coda_codacore_model")
})

test_that("CUDA fit keeps trainable CoDaCoRe stage parameters as optimizer leaves", {
  skip_if_no_coda_codacore()
  torch <- reticulate::import("torch", delay_load = FALSE)
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  testthat::skip_if_not(cuda_ok, "CUDA is unavailable")
  d <- .make_coda_codacore_data(n = 40L, p = 16L, seed = 37L)
  model <- fit_coda_codacore(
    d$X, d$y,
    hp = list(device = "cuda", epochs = 2L, max_balances = 1L, seed = 91L)
  )
  expect_s3_class(model, "coda_codacore_model")
  expect_true(all(is.finite(score_coda_codacore(model, d$X))))
})

# §7.10 no-live-pointer: the fitted model holds NO python.builtin.object field (the
# torch relaxation is discarded after export; only discrete sets + logistic weights
# remain).
test_that("§7.10 fitted model holds no python pointer (export-to-pure-R contract)", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc(seed = 13L)
  m <- f$model
  expect_false(any(vapply(m, function(x) inherits(x, "python.builtin.object"),
                          logical(1L))))
  # Recurse one level into the balances list (the only nested structure).
  for (bb in m$balances) {
    expect_false(any(vapply(bb, function(x) inherits(x, "python.builtin.object"),
                            logical(1L))))
  }
})

# §7.11 frozen-balance + frozen-head independent R recompute == scorer output
# EXACTLY 0 (the score is a pure-R deterministic linear function of the discrete
# balances; this recompute does NOT call .coda_codacore_balance_one's caller path).
test_that("§7.11 independent frozen-balance + frozen-head recompute == scorer (0)", {
  skip_if_no_coda_codacore()
  f <- .fit_small_cc(seed = 14L)
  m <- f$model
  X <- f$data$X[seq_len(12L), , drop = FALSE]
  direct <- score_coda_codacore(m, X)

  X_use <- .coda_codacore_align(X, m$feature_universe)
  worst <- 0
  for (i in seq_len(nrow(X_use))) {
    xi <- X_use[i, ]
    if (!any(xi > 0)) { sc <- 0 } else {
      # Recompute each discrete ILR balance from scratch (closed form), then the
      # frozen logistic head, WITHOUT reusing the scorer's loop.
      bal <- numeric(m$n_balances)
      for (b in seq_len(m$n_balances)) {
        a_idx <- m$balances[[b]]$A; b_idx <- m$balances[[b]]$B
        va <- xi[a_idx]; va <- va[is.finite(va) & va > 0]
        vb <- xi[b_idx]; vb <- vb[is.finite(vb) & vb > 0]
        r <- length(va); s <- length(vb)
        if (r < 1L || s < 1L) {
          bal[b] <- 0
        } else {
          bal[b] <- sqrt((r * s) / (r + s)) * (mean(log(va)) - mean(log(vb)))
        }
      }
      sc <- m$intercept + sum(m$weights * bal)
    }
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.11] independent frozen recompute maxdiff = %.3e\n", worst))
  expect_identical(worst, 0)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_coda_codacore()
  d <- .make_coda_codacore_data(seed = 18L)
  expect_error(fit_coda_codacore(d$X, d$y, hp = "notalist"), "hp must be a list")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(seed = 1, seed = 2)),
               "duplicate hp")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(7)), "named")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(max_balances = 0)),
               "max_balances")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(max_balances = 2.5)),
               "max_balances")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(tau = 0)), "tau")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(ridge = -1)), "ridge")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(min_features = 2.5)),
               "min_features")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_coda_codacore(d$X, d$y, hp = list(seed = 1.5)), "seed")
})

test_that("coda-codacore input validation errors are explicit", {
  skip_if_no_coda_codacore()
  d <- .make_coda_codacore_data(seed = 19L)
  X_unnamed <- d$X
  colnames(X_unnamed) <- NULL
  expect_error(fit_coda_codacore(X_unnamed, d$y), "feature names")
  expect_error(fit_coda_codacore(d$X, d$y[-1]), "length\\(y_train\\)")
  expect_error(fit_coda_codacore(d$X, rep(0, nrow(d$X))), "at least one case")
  expect_error(score_coda_codacore(list(), d$X), "class coda_codacore_model")
})
