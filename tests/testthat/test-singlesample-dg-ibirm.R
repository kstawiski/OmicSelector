# dg-ibirm trains a BN-free MLP encoder + a Linear(d->1) head JOINTLY across
# ENVIRONMENTS (cohorts) by the IB-IRM objective (Ahuja 2021) via reticulate-python
# torch at FIT, then FREEZES and scores in pure R (torch confined to fit). IB-IRM is
# dro-vrex/dg-fishr with ONE change: the invariance penalty becomes the IRMv1
# gradient-penalty PLUS an Information-Bottleneck embedding-variance penalty. Tests
# that fit are skipped when reticulate or the venv torch is unavailable, so the
# package suite stays green on hosts without the venv. Plain torch training needs NO
# brew-OpenSSL bridge and NO HF download, so the skip-guard only probes reticulate +
# `import torch`.
skip_if_no_ibirm <- function() {
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

# Synthetic MULTI-ENVIRONMENT generator. A CASE-ONLY location shift in the first `k`
# of `p` log-abundance features carries the discriminative signal in EVERY cohort
# (the cohort-invariant signal IB-IRM should keep). Each cohort additionally gets a
# cohort-specific additive offset on a DIFFERENT, non-overlapping feature block
# (a nuisance shift that is NOT label-predictive), so the cohorts differ in
# distribution. `n_env` cohorts of equal size, each with both classes.
.make_ibirm_data <- function(n = 120L, p = 30L, k = 8L, shift = 1.4,
                             sd = 0.5, n_env = 3L, seed = 11L) {
  set.seed(seed)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  # Assign cohorts in contiguous blocks, then BALANCE the label WITHIN each cohort
  # (alternating 0/1 inside each block). This guarantees every cohort carries both
  # classes regardless of n_env (env period and label period must not be coupled --
  # an interleaved global y with an even n_env would make each cohort single-class).
  env <- rep(paste0("GSE", seq_len(n_env)), length.out = n)
  env <- sort(env)                               # contiguous cohort blocks
  y <- integer(n)
  for (e in unique(env)) {
    rows <- which(env == e)
    y[rows] <- rep(c(0L, 1L), length.out = length(rows))
  }
  case <- which(y == 1L)
  L[case, seq_len(k)] <- L[case, seq_len(k)] + shift   # cohort-invariant case signal
  # cohort-specific nuisance offset on a disjoint feature block (not label-linked)
  for (e in seq_len(n_env)) {
    rows <- which(env == paste0("GSE", e))
    blk <- (k + 1L + (e - 1L) * 2L)
    blk <- blk[blk <= p]
    if (length(blk) > 0L) L[rows, blk] <- L[rows, blk] + 0.8 * e
  }
  X <- exp(L)
  list(X = X, y = y, env = env,
       meta = data.frame(accession = env, stringsAsFactors = FALSE),
       block = features[seq_len(k)])
}

.auc_mw_ib <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}


test_that("fit returns a well-formed dg_ibirm_model (IB-IRM path engaged)", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data()
  m <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                    hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "dg_ibirm_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$embedding_dim, 32L)             # last hidden = embedding dim
  expect_length(m$weights, 2L)                   # 2 linear layers up to embedding
  expect_equal(dim(m$weights[[1]]$W), c(64L, ncol(d$X)))
  expect_equal(dim(m$weights[[2]]$W), c(32L, 64L))
  expect_length(m$head_w, 32L)
  expect_length(m$head_b, 1L)
  expect_equal(m$n_environments, 3L)             # 3 both-class cohorts -> IB-IRM
  s <- score_dg_ibirm(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("§7.1 R forward matches the torch encoder embedding (~1e-6)", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 21L)
  m <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                    hp = list(epochs = 50L, device = "cpu", seed = 7L))
  # Embed the training rows with the EXPORTED R-side encoder weights (pure-R
  # forward), and with a torch Sequential rebuilt from those SAME exported weights;
  # compare. The float64 R forward must match the float32 torch forward on the
  # exported (float32-valued) weights to ~1e-6.
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  Z_in <- .dg_ibirm_rclr_matrix(d$X)
  layers <- list()
  L <- length(m$weights)
  for (l in seq_len(L)) {
    W <- m$weights[[l]]$W; b <- m$weights[[l]]$b
    lin <- nn$Linear(ncol(W), nrow(W))
    lin$weight$data <- torch$as_tensor(np$asarray(W, dtype = "float32"))$float()
    lin$bias$data   <- torch$as_tensor(np$asarray(b, dtype = "float32"))$float()
    layers[[length(layers) + 1L]] <- lin
    if (l < L) layers[[length(layers) + 1L]] <- nn$ReLU()
  }
  net <- do.call(nn$Sequential, layers); net$eval()
  X_py <- torch$as_tensor(np$asarray(Z_in, dtype = "float32"))$float()
  torch_emb <- reticulate::py_to_r(net(X_py)$detach()$to(torch$float64)$numpy())
  r_emb <- .dg_ibirm_forward(Z_in, m$weights, "relu")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("\n[§7.1] R-vs-torch encoder embedding maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                            # float32 weights: ~1e-6
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 31L)
  m <- fit_dg_ibirm(d$X[1:90, ], d$y[1:90], meta_train = d$meta[1:90, , drop = FALSE],
                    hp = list(epochs = 60L))
  X_test <- d$X[91:120, ]
  batch <- score_dg_ibirm(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_dg_ibirm(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 41L)
  m <- fit_dg_ibirm(d$X[1:90, ], d$y[1:90], meta_train = d$meta[1:90, , drop = FALSE],
                    hp = list(epochs = 60L))
  X_test <- d$X[91:120, ]
  base <- score_dg_ibirm(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_dg_ibirm(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (maxdiff ~0)", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 44L)
  m <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 61L, 100L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_dg_ibirm(m, row * s) -
                                score_dg_ibirm(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_dg_ibirm(m, scaled) -
                                score_dg_ibirm(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-7)
})

test_that("§7.5 held-out AUC on a planted cross-cohort shift is high (> 0.7)", {
  skip_if_no_ibirm()
  # Train on cohorts GSE1/GSE2/GSE3; evaluate on a HELD-OUT cohort GSE4 that shares
  # the cohort-invariant case signal but has its own nuisance offset. IB-IRM should
  # generalise the case signal to the unseen cohort.
  d <- .make_ibirm_data(n = 200L, n_env = 4L, seed = 51L)
  tr <- which(d$env %in% c("GSE1", "GSE2", "GSE3"))
  te <- which(d$env == "GSE4")
  m <- fit_dg_ibirm(d$X[tr, ], d$y[tr], meta_train = d$meta[tr, , drop = FALSE],
                    hp = list(epochs = 300L, irm_lambda = 1.0, ib_lambda = 1.0))
  expect_gt(m$n_environments, 1L)
  s <- score_dg_ibirm(m, d$X[te, ])
  auc <- .auc_mw_ib(d$y[te], s)
  cat(sprintf("[§7.5] held-out-cohort AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.6 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 61L)
  m <- fit_dg_ibirm(d$X[1:90, ], d$y[1:90], meta_train = d$meta[1:90, , drop = FALSE],
                    hp = list(epochs = 60L))
  X_test <- d$X[91:120, ]
  score_fun <- function(model, X, meta) score_dg_ibirm(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .dg_ibirm_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .dg_ibirm_model_digest))

  # Canonical dispatch: post-integration the manifest carries the dg-ibirm row and
  # score_dg_ibirm is in the namespace; pre-integration (staged self-test) clone a
  # one-row roster from a committed within-discriminator row and register a thin
  # adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_dg_ibirm(m, X_test)
  roster <- singlesample_method_roster()
  if (!"dg-ibirm" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "dg-ibirm"
    tmpl$fit_fn <- "fit_dg_ibirm"
    tmpl$score_fn <- "score_dg_ibirm"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_dg_ibirm", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "dg-ibirm",
      function(model, X, meta) score_dg_ibirm(model, X, meta))
  }
  via <- singlesample_score_call("dg-ibirm", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.7 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 71L)
  m <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                    hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_dg_ibirm(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_dg_ibirm(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_dg_ibirm(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be embedded and scored, NOT floored.
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .dg_ibirm_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  emb <- .dg_ibirm_forward(matrix(z_flat, nrow = 1L), m$weights, m$activation)
  manual <- .dg_ibirm_score_one(as.numeric(emb), m$head_w, m$head_b)
  s_flat <- score_dg_ibirm(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (manual %.6f)\n",
              s_flat, manual))
  expect_equal(s_flat, manual, tolerance = 1e-12)

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_dg_ibirm(m, empty), 0)
})

test_that("§7.8 determinism: two seed=42 CPU fits give identical digests (hardened)", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 81L)
  m1 <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                     hp = list(epochs = 60L, device = "cpu", seed = 42L))

  # HARDEN: perturb the GLOBAL torch RNG BETWEEN the two seed=42 fits. A correctly
  # seeded fit (manual_seed BEFORE module construction) is INVARIANT to the inter-fit
  # generator state; the seed-after-build bug would not be. The fit also restores the
  # torch RNG on exit, so this perturbation must not leak into m2 either way -- the
  # point is that even the entry-state torch RNG is irrelevant to a correct fit.
  torch <- reticulate::import("torch", delay_load = FALSE)
  torch$manual_seed(99L)
  invisible(torch$randn(list(1000L, 50L)))       # advance the global torch generator

  m2 <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                     hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .dg_ibirm_model_digest(m1)
  dg2 <- .dg_ibirm_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s (inter-fit torch-RNG perturbed)\n",
              dg1, dg2))
  expect_identical(dg1, dg2)
  # Scores also bit-identical across the two fits.
  s1 <- score_dg_ibirm(m1, d$X[1:10, ])
  s2 <- score_dg_ibirm(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.9 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                    hp = list(epochs = 60L, device = "cpu"))

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
  invisible(score_dg_ibirm(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                     hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "dg_ibirm_model")
})

test_that("§7.10 independent linear-head recompute bit-matches the scorer", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 101L)
  m <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_dg_ibirm(m, X_test)
  # From-scratch recompute: w . z + b via crossprod (non-circular: does not call
  # .dg_ibirm_score_one).
  X_use <- .dg_ibirm_align(X_test, m$feature_universe)
  worst <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .dg_ibirm_rclr(X_use[i, ])
    emb <- as.numeric(.dg_ibirm_forward(matrix(z, nrow = 1L), m$weights,
                                        m$activation))
    sc <- as.numeric(crossprod(m$head_w, emb)) + m$head_b
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.10] independent linear-head recompute maxdiff = %.3e\n",
              worst))
  expect_lt(worst, 1e-9)
})

test_that("§7.11 IB-IRM penalties flow; ERM fallback triggers on NULL/single-cohort meta", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 113L)

  # (a) The penalties are ACTIVE: with >= 2 environments, irm_lambda = ib_lambda = 0
  # (pure ERM-style mean risk) vs LARGE lambdas give measurably different trained
  # weights. This proves the IRMv1 gradient-penalty + IB variance penalty actually
  # contribute gradient. Same seed, data, meta, epochs -> the ONLY difference is the
  # lambdas, so any weight difference is the IB-IRM penalties.
  m_l0 <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                       hp = list(epochs = 120L, seed = 5L,
                                 irm_lambda = 0.0, ib_lambda = 0.0))
  m_lL <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                       hp = list(epochs = 120L, seed = 5L,
                                 irm_lambda = 50.0, ib_lambda = 50.0))
  expect_equal(m_l0$n_environments, 3L)
  expect_equal(m_lL$n_environments, 3L)
  wd <- max(abs(m_l0$head_w - m_lL$head_w),
            abs(m_l0$weights[[1]]$W - m_lL$weights[[1]]$W))
  cat(sprintf("[§7.11] lambdas=0 vs lambdas=50 max weight diff = %.4e\n", wd))
  expect_gt(wd, 1e-3)                            # penalties measurably change the fit

  # (a2) IRM and IB penalties are INDEPENDENTLY active: each alone moves the fit off
  # the pure-ERM solution (so neither penalty is a no-op that the other masks).
  m_irm <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                        hp = list(epochs = 120L, seed = 5L,
                                  irm_lambda = 50.0, ib_lambda = 0.0))
  m_ib  <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                        hp = list(epochs = 120L, seed = 5L,
                                  irm_lambda = 0.0, ib_lambda = 50.0))
  wd_irm <- max(abs(m_l0$head_w - m_irm$head_w),
                abs(m_l0$weights[[1]]$W - m_irm$weights[[1]]$W))
  wd_ib  <- max(abs(m_l0$head_w - m_ib$head_w),
                abs(m_l0$weights[[1]]$W - m_ib$weights[[1]]$W))
  cat(sprintf("[§7.11] IRM-only weight diff = %.4e ; IB-only weight diff = %.4e\n",
              wd_irm, wd_ib))
  expect_gt(wd_irm, 1e-3)                        # IRMv1 gradient penalty flows
  expect_gt(wd_ib, 1e-3)                         # IB variance penalty flows

  # (b) ERM fallback on NULL meta -> n_environments == 1.
  m_null <- fit_dg_ibirm(d$X, d$y, meta_train = NULL, hp = list(epochs = 30L))
  expect_equal(m_null$n_environments, 1L)
  s_null <- score_dg_ibirm(m_null, d$X[1:5, ])
  expect_true(all(is.finite(s_null)))

  # (c) ERM fallback on a SINGLE-cohort meta -> n_environments == 1.
  meta1 <- data.frame(accession = rep("GSE_only", nrow(d$X)),
                      stringsAsFactors = FALSE)
  m_one <- fit_dg_ibirm(d$X, d$y, meta_train = meta1, hp = list(epochs = 30L))
  expect_equal(m_one$n_environments, 1L)

  # (d) ERM fallback when fewer than 2 cohorts carry both classes: one valid
  # (both-class) cohort + one control-only cohort -> only 1 valid environment ->
  # n_environments == 1. (single-class env DROPPED -> falls below 2 valid envs.)
  env2 <- rep("GSE_valid", nrow(d$X))
  ctrl_idx <- which(d$y == 0L)
  env2[utils::head(ctrl_idx, 5L)] <- "GSE_ctrlonly"   # a control-only second cohort
  meta2 <- data.frame(accession = env2, stringsAsFactors = FALSE)
  m_sc <- fit_dg_ibirm(d$X, d$y, meta_train = meta2, hp = list(epochs = 30L))
  expect_equal(m_sc$n_environments, 1L)          # <2 both-class envs -> ERM

  # (e) MIXED case: >= 2 VALID both-class cohorts PLUS extra rows carrying a missing
  # (NA) or single-class cohort label. Those extra rows are DROPPED from training
  # (they belong to no valid environment), while the valid cohorts are retained, so
  # the IB-IRM penalty is still applied over the kept rows and n_environments stays
  # at the number of valid cohorts (NOT collapsed to 1, NOT inflated by the dropped
  # rows). This exercises the keep/drop branch that fires only when >= 2 valid envs
  # coexist with droppable rows -- the branch left without a dedicated test before.
  # Reproduce the generator's contiguous 3-cohort assignment, then carve out a
  # control-only cohort (5 controls) and NA-out 5 cases. Each of GSE1/2/3 loses at
  # most 5 of its ~20-per-class rows, so all three stay both-class (valid).
  base_env <- sort(rep(paste0("GSE", seq_len(3L)), length.out = nrow(d$X)))
  acc_mix <- base_env
  acc_mix[utils::head(which(d$y == 0L), 5L)] <- "GSE_ctrlonly"  # single-class -> dropped
  acc_mix[utils::tail(which(d$y == 1L), 5L)] <- NA              # missing      -> dropped
  meta_mix <- data.frame(accession = acc_mix, stringsAsFactors = FALSE)
  m_mix <- fit_dg_ibirm(d$X, d$y, meta_train = meta_mix,
                        hp = list(epochs = 40L, seed = 7L))
  cat(sprintf("[§7.11] mixed (>=2 valid + dropped rows): n_environments = %d\n",
              m_mix$n_environments))
  expect_equal(m_mix$n_environments, 3L)         # GSE1/2/3 retained; ctrlonly + NA dropped
  s_mix <- score_dg_ibirm(m_mix, d$X[1:5, ])
  expect_true(all(is.finite(s_mix)))             # fit on the kept rows succeeded + scores
})

test_that("§7.12 tanh activation: R forward matches torch + single-row == batch", {
  skip_if_no_ibirm()
  # Parity test for the `tanh` activation branch of .dg_ibirm_forward, mirroring
  # §7.1 (relu). Additive coverage over a byte-unchanged scorer; tanh is a real
  # branch selected by hp$activation = "tanh".
  d <- .make_ibirm_data(seed = 23L)
  m <- fit_dg_ibirm(d$X, d$y, meta_train = d$meta,
                    hp = list(epochs = 50L, device = "cpu", seed = 7L,
                              activation = "tanh"))
  expect_identical(m$activation, "tanh")
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  Z_in <- .dg_ibirm_rclr_matrix(d$X)
  layers <- list()
  L <- length(m$weights)
  for (l in seq_len(L)) {
    W <- m$weights[[l]]$W; b <- m$weights[[l]]$b
    lin <- nn$Linear(ncol(W), nrow(W))
    lin$weight$data <- torch$as_tensor(np$asarray(W, dtype = "float32"))$float()
    lin$bias$data   <- torch$as_tensor(np$asarray(b, dtype = "float32"))$float()
    layers[[length(layers) + 1L]] <- lin
    if (l < L) layers[[length(layers) + 1L]] <- nn$Tanh()  # tanh between layers
  }
  net <- do.call(nn$Sequential, layers); net$eval()
  X_py <- torch$as_tensor(np$asarray(Z_in, dtype = "float32"))$float()
  torch_emb <- reticulate::py_to_r(net(X_py)$detach()$to(torch$float64)$numpy())
  r_emb <- .dg_ibirm_forward(Z_in, m$weights, "tanh")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("[§7.12] tanh R-vs-torch embedding maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                              # float32 weights: ~1e-7
  # the decisive single-sample property holds under tanh too
  batch <- score_dg_ibirm(m, d$X)
  singles <- vapply(seq_len(nrow(d$X)), function(i) {
    score_dg_ibirm(m, d$X[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.12] tanh single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.13 IB-IRM penalty math is correct on a hand-built 2-env toy (< 1e-6)", {
  skip_if_no_ibirm()
  # NON-tautological check of the EXACT IB-IRM penalties the trainer computes. Build a
  # tiny 2-environment toy with KNOWN embeddings z, head (w0, b0), labels y; compute
  # the reference IRMv1 gradient penalty and the IB embedding-variance penalty in
  # PLAIN R from first principles, then reproduce each in torch using the SAME tensor
  # ops the trainer uses (autograd.grad with create_graph for IRM; explicit biased
  # variance for IB), and require agreement < 1e-6.
  torch <- reticulate::import("torch", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)

  # Toy: d = 2 embedding, 2 environments. Hand-pick embeddings z, head (w0,b0), labels.
  z <- matrix(c( 1.0,  2.0,    # env1 row1
                -0.5,  0.3,    # env1 row2
                 0.7, -1.1,    # env2 row1
                 2.0,  0.0,    # env2 row2
                -1.3,  0.6),   # env2 row3
              ncol = 2L, byrow = TRUE)
  w0 <- c(0.4, -0.7)                                 # head weight
  b0 <- 0.15                                         # head bias
  y  <- c(1,    0,    1,    0,    1)                  # labels
  env <- c(1L, 1L, 2L, 2L, 2L)                       # 2-row env1, 3-row env2
  l0 <- as.numeric(z %*% w0) + b0                    # head logits l0 = w0.z + b0

  # --- reference penalties in PLAIN R from first principles ------------------
  # IRMv1 gradient penalty: with a dummy scalar w applied to the logit, R_e(w) =
  # mean BCE-with-logits(w*l0_e, y_e). The derivative of BCE-with-logits per row
  # w.r.t. w at w=1 is (sigmoid(w*l0) - y) * l0 evaluated at w=1 = (sigmoid(l0)-y)*l0;
  # R_e is the mean over env e's rows, so dR_e/dw|_{w=1} = mean_{i in e}((sig(l0_i)-y_i)*l0_i).
  sig <- function(x) 1 / (1 + exp(-x))
  drow <- (sig(l0) - y) * l0                          # d/dw of per-row BCE at w=1
  per_env_grad <- function(rows) mean(drow[rows])
  gw1 <- per_env_grad(which(env == 1L))
  gw2 <- per_env_grad(which(env == 2L))
  omega_irm_ref <- mean(c(gw1^2, gw2^2))              # mean_e ||dR_e/dw||^2

  # IB embedding-variance penalty: mean over embedding dims of the biased (1/n)
  # variance of the pooled embedding.
  n <- nrow(z)
  zmean <- colMeans(z)
  zvar <- colMeans(sweep(z, 2L, zmean, `-`)^2)        # biased 1/n per-dim variance
  omega_ib_ref <- mean(zvar)                          # mean over embedding dims

  # --- penalties via the SAME torch ops the trainer uses --------------------
  emb_t   <- torch$as_tensor(np$asarray(z, dtype = "float64"))$double()
  w0_t    <- torch$as_tensor(np$asarray(matrix(w0, ncol = 1L),
                                        dtype = "float64"))$double()
  y_t     <- torch$as_tensor(np$asarray(matrix(y, ncol = 1L),
                                        dtype = "float64"))$double()
  logit_t <- torch$add(torch$matmul(emb_t, w0_t), b0) # n x 1, l0 = w0.z + b0
  lossf   <- torch$nn$BCEWithLogitsLoss(reduction = "none")
  env_tensors <- list(
    torch$as_tensor(np$asarray(as.integer(which(env == 1L) - 1L),
                               dtype = "int64"))$long(),
    torch$as_tensor(np$asarray(as.integer(which(env == 2L) - 1L),
                               dtype = "int64"))$long())

  # IRMv1: dummy scalar w = 1.0, autograd.grad with create_graph (exactly the trainer)
  dummy_w <- torch$ones(list(1L), dtype = logit_t$dtype, requires_grad = TRUE)
  irm_terms <- lapply(env_tensors, function(rows0) {
    logit_e <- logit_t$index_select(0L, rows0)
    y_e <- y_t$index_select(0L, rows0)
    scaled_e <- torch$mul(logit_e, dummy_w)
    r_e <- lossf(scaled_e, y_e)$mean()
    gw <- torch$autograd$grad(r_e, dummy_w, create_graph = TRUE)[[1L]]
    gw$pow(2)$sum()
  })
  omega_irm_t <- reticulate::py_to_r(torch$stack(irm_terms)$mean()$item())

  # IB: explicit biased 1/n variance per dim, mean over dims (exactly the trainer)
  z_mean_t <- emb_t$mean(dim = 0L, keepdim = TRUE)
  z_var_t  <- emb_t$sub(z_mean_t)$pow(2)$mean(dim = 0L)
  omega_ib_t <- reticulate::py_to_r(z_var_t$mean()$item())

  diff_irm <- abs(omega_irm_ref - omega_irm_t)
  diff_ib  <- abs(omega_ib_ref - omega_ib_t)
  cat(sprintf("[§7.13] IRM penalty: R-ref = %.10f ; torch = %.10f ; diff = %.3e\n",
              omega_irm_ref, omega_irm_t, diff_irm))
  cat(sprintf("[§7.13] IB  penalty: R-ref = %.10f ; torch = %.10f ; diff = %.3e\n",
              omega_ib_ref, omega_ib_t, diff_ib))
  # Non-tautology guard: both penalties are genuinely non-zero on this toy (the two
  # envs have DIFFERENT per-env risk gradients, and the embedding has spread), so the
  # < 1e-6 match is meaningful.
  expect_gt(omega_irm_ref, 1e-6)
  expect_gt(omega_ib_ref, 1e-6)
  expect_lt(diff_irm, 1e-6)
  expect_lt(diff_ib, 1e-6)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_ibirm()
  d <- .make_ibirm_data(seed = 111L)
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(weight_decay = -1)),
               "weight_decay")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(hidden = c(0L, 4L))), "hidden")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(activation = "gelu")),
               "activation")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(irm_lambda = -1)),
               "irm_lambda")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(ib_lambda = -1)),
               "ib_lambda")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(cohort_col = "")), "cohort_col")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(cohort_col = c("a", "b"))),
               "cohort_col")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_dg_ibirm(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_dg_ibirm(d$X, d$y,
                 hp = structure(list(1e-3, 1e-3), names = c("lr", "lr"))),
    "duplicate hp")
})
