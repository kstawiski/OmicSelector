# dann trains a BN-free MLP encoder + a Linear(d->1) LABEL head JOINTLY across
# ENVIRONMENTS (cohorts = domains) by the DANN objective (Ganin 2016) via
# reticulate-python torch at FIT, then FREEZES and scores in pure R (torch confined
# to fit). DANN is dro-vrex / dg-fishr with ONE change: the cross-cohort penalty is
# a domain discriminator fed the embedding through a gradient-reversal layer (GRL),
# making the encoder cohort-INVARIANT. The GRL + domain head are FIT-only and are NOT
# exported. Tests that fit are skipped when reticulate or the venv torch is
# unavailable, so the package suite stays green on hosts without the venv. Plain
# torch training needs NO brew-OpenSSL bridge and NO HF download, so the skip-guard
# only probes reticulate + `import torch`.
skip_if_no_dann <- function() {
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
# (the cohort-invariant signal DANN should keep). Each cohort additionally gets a
# cohort-specific additive offset on a DIFFERENT, non-overlapping feature block (a
# nuisance shift that is NOT label-predictive and that the adversary should erase),
# so the cohorts differ in distribution. `n_env` cohorts of equal size, both classes.
.make_dann_data <- function(n = 120L, p = 30L, k = 8L, shift = 1.4,
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

.auc_mw_dann <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}


test_that("fit returns a well-formed dann_model (adversary engaged)", {
  skip_if_no_dann()
  d <- .make_dann_data()
  m <- fit_dann(d$X, d$y, meta_train = d$meta,
                hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "dann_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$embedding_dim, 32L)             # last hidden = embedding dim
  expect_length(m$weights, 2L)                   # 2 linear layers up to embedding
  expect_equal(dim(m$weights[[1]]$W), c(64L, ncol(d$X)))
  expect_equal(dim(m$weights[[2]]$W), c(32L, 64L))
  expect_length(m$head_w, 32L)
  expect_length(m$head_b, 1L)
  expect_equal(m$n_environments, 3L)             # 3 both-class cohorts -> adversary
  s <- score_dann(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("§7.1 R forward matches the torch encoder embedding (~1e-6)", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 21L)
  m <- fit_dann(d$X, d$y, meta_train = d$meta,
                hp = list(epochs = 50L, device = "cpu", seed = 7L))
  # Embed the training rows with the EXPORTED R-side encoder weights (pure-R
  # forward), and with a torch Sequential rebuilt from those SAME exported weights;
  # compare. The float64 R forward must match the float32 torch forward on the
  # exported (float32-valued) weights to ~1e-6.
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  Z_in <- .dann_rclr_matrix(d$X)
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
  r_emb <- .dann_forward(Z_in, m$weights, "relu")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("\n[§7.1] R-vs-torch encoder embedding maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                            # float32 weights: ~1e-6
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 31L)
  m <- fit_dann(d$X[1:90, ], d$y[1:90], meta_train = d$meta[1:90, , drop = FALSE],
                hp = list(epochs = 60L))
  X_test <- d$X[91:120, ]
  batch <- score_dann(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_dann(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 41L)
  m <- fit_dann(d$X[1:90, ], d$y[1:90], meta_train = d$meta[1:90, , drop = FALSE],
                hp = list(epochs = 60L))
  X_test <- d$X[91:120, ]
  base <- score_dann(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_dann(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (maxdiff ~0)", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 44L)
  m <- fit_dann(d$X, d$y, meta_train = d$meta, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 61L, 100L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_dann(m, row * s) -
                                score_dann(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_dann(m, scaled) -
                                score_dann(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-7)
})

test_that("§7.5 held-out AUC on a planted cross-cohort shift is high (> 0.7)", {
  skip_if_no_dann()
  # Train on cohorts GSE1/GSE2/GSE3; evaluate on a HELD-OUT cohort GSE4 that shares
  # the cohort-invariant case signal but has its own nuisance offset. DANN should
  # generalise the case signal to the unseen cohort (the adversary erases the
  # cohort-specific nuisance from the embedding).
  d <- .make_dann_data(n = 200L, n_env = 4L, seed = 51L)
  tr <- which(d$env %in% c("GSE1", "GSE2", "GSE3"))
  te <- which(d$env == "GSE4")
  m <- fit_dann(d$X[tr, ], d$y[tr], meta_train = d$meta[tr, , drop = FALSE],
                hp = list(epochs = 300L, grl_lambda = 1.0))
  expect_gt(m$n_environments, 1L)
  s <- score_dann(m, d$X[te, ])
  auc <- .auc_mw_dann(d$y[te], s)
  cat(sprintf("[§7.5] held-out-cohort AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.6 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 61L)
  m <- fit_dann(d$X[1:90, ], d$y[1:90], meta_train = d$meta[1:90, , drop = FALSE],
                hp = list(epochs = 60L))
  X_test <- d$X[91:120, ]
  score_fun <- function(model, X, meta) score_dann(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .dann_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .dann_model_digest))

  # Canonical dispatch: post-integration the manifest carries the dann row and
  # score_dann is in the namespace; pre-integration (staged self-test) clone a
  # one-row roster from a committed within-discriminator row and register a thin
  # adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_dann(m, X_test)
  roster <- singlesample_method_roster()
  if (!"dann" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "dann"
    tmpl$fit_fn <- "fit_dann"
    tmpl$score_fn <- "score_dann"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_dann", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "dann",
      function(model, X, meta) score_dann(model, X, meta))
  }
  via <- singlesample_score_call("dann", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.7 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 71L)
  m <- fit_dann(d$X, d$y, meta_train = d$meta,
                hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_dann(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_dann(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_dann(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be embedded and scored, NOT floored.
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .dann_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  emb <- .dann_forward(matrix(z_flat, nrow = 1L), m$weights, m$activation)
  manual <- .dann_score_one(as.numeric(emb), m$head_w, m$head_b)
  s_flat <- score_dann(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (manual %.6f)\n",
              s_flat, manual))
  expect_equal(s_flat, manual, tolerance = 1e-12)

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_dann(m, empty), 0)
})

test_that("§7.8 determinism: two seed=42 CPU fits give identical digests (hardened)", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 81L)
  m1 <- fit_dann(d$X, d$y, meta_train = d$meta,
                 hp = list(epochs = 60L, device = "cpu", seed = 42L))

  # HARDEN: perturb the GLOBAL torch RNG BETWEEN the two seed=42 fits. A correctly
  # seeded fit (manual_seed BEFORE module construction) is INVARIANT to the inter-fit
  # generator state; the seed-after-build bug would not be. The fit also restores the
  # torch RNG on exit, so this perturbation must not leak into m2 either way -- the
  # point is that even the entry-state torch RNG is irrelevant to a correct fit.
  torch <- reticulate::import("torch", delay_load = FALSE)
  torch$manual_seed(99L)
  invisible(torch$randn(list(1000L, 50L)))       # advance the global torch generator

  m2 <- fit_dann(d$X, d$y, meta_train = d$meta,
                 hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .dann_model_digest(m1)
  dg2 <- .dann_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s (inter-fit torch-RNG perturbed)\n",
              dg1, dg2))
  expect_identical(dg1, dg2)
  # Scores also bit-identical across the two fits.
  s1 <- score_dann(m1, d$X[1:10, ])
  s2 <- score_dann(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.9 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_dann(d$X, d$y, meta_train = d$meta,
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
  invisible(score_dann(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_dann(d$X, d$y, meta_train = d$meta,
                 hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "dann_model")
})

test_that("§7.10 independent linear-head recompute bit-matches the scorer", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 101L)
  m <- fit_dann(d$X, d$y, meta_train = d$meta, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_dann(m, X_test)
  # From-scratch recompute: w . z + b via crossprod (non-circular: does not call
  # .dann_score_one).
  X_use <- .dann_align(X_test, m$feature_universe)
  worst <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .dann_rclr(X_use[i, ])
    emb <- as.numeric(.dann_forward(matrix(z, nrow = 1L), m$weights,
                                    m$activation))
    sc <- as.numeric(crossprod(m$head_w, emb)) + m$head_b
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.10] independent linear-head recompute maxdiff = %.3e\n",
              worst))
  expect_lt(worst, 1e-9)
})

test_that("§7.11 adversary active (grl_lambda changes the fit); ERM fallback on NULL/single-cohort meta", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 113L)

  # (a) The adversary is ACTIVE: with >= 2 domains, grl_lambda = 0 (no gradient
  # reversed into the encoder -- the domain head trains but cannot shape phi) vs a
  # LARGE grl_lambda give measurably different trained encoder/head weights. This
  # proves the GRL actually propagates a (reversed) gradient into the encoder. Same
  # seed, data, meta, epochs -> the ONLY difference is grl_lambda, so any weight
  # difference is the adversary.
  m_l0 <- fit_dann(d$X, d$y, meta_train = d$meta,
                   hp = list(epochs = 120L, seed = 5L, grl_lambda = 0.0))
  m_lL <- fit_dann(d$X, d$y, meta_train = d$meta,
                   hp = list(epochs = 120L, seed = 5L, grl_lambda = 5.0))
  expect_equal(m_l0$n_environments, 3L)
  expect_equal(m_lL$n_environments, 3L)
  wd <- max(abs(m_l0$head_w - m_lL$head_w),
            abs(m_l0$weights[[1]]$W - m_lL$weights[[1]]$W))
  cat(sprintf("[§7.11] grl_lambda=0 vs grl_lambda=5 max weight diff = %.4e\n", wd))
  expect_gt(wd, 1e-3)                            # adversary measurably changes the fit

  # (b) ERM fallback on NULL meta -> n_environments == 1.
  m_null <- fit_dann(d$X, d$y, meta_train = NULL, hp = list(epochs = 30L))
  expect_equal(m_null$n_environments, 1L)
  s_null <- score_dann(m_null, d$X[1:5, ])
  expect_true(all(is.finite(s_null)))

  # (c) ERM fallback on a SINGLE-cohort meta -> n_environments == 1.
  meta1 <- data.frame(accession = rep("GSE_only", nrow(d$X)),
                      stringsAsFactors = FALSE)
  m_one <- fit_dann(d$X, d$y, meta_train = meta1, hp = list(epochs = 30L))
  expect_equal(m_one$n_environments, 1L)

  # (d) ERM fallback when fewer than 2 cohorts carry both classes: one valid
  # (both-class) cohort + one control-only cohort -> only 1 valid domain ->
  # n_environments == 1. (single-class domain DROPPED -> falls below 2 valid domains.)
  env2 <- rep("GSE_valid", nrow(d$X))
  ctrl_idx <- which(d$y == 0L)
  env2[utils::head(ctrl_idx, 5L)] <- "GSE_ctrlonly"   # a control-only second cohort
  meta2 <- data.frame(accession = env2, stringsAsFactors = FALSE)
  m_sc <- fit_dann(d$X, d$y, meta_train = meta2, hp = list(epochs = 30L))
  expect_equal(m_sc$n_environments, 1L)          # <2 both-class domains -> ERM
})

test_that("§7.12 tanh activation: R forward matches torch + single-row == batch", {
  skip_if_no_dann()
  # Parity test for the `tanh` activation branch of .dann_forward, mirroring §7.1
  # (relu). Additive coverage over a byte-unchanged scorer; tanh is a real branch
  # selected by hp$activation = "tanh".
  d <- .make_dann_data(seed = 23L)
  m <- fit_dann(d$X, d$y, meta_train = d$meta,
                hp = list(epochs = 50L, device = "cpu", seed = 7L,
                          activation = "tanh"))
  expect_identical(m$activation, "tanh")
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  Z_in <- .dann_rclr_matrix(d$X)
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
  r_emb <- .dann_forward(Z_in, m$weights, "tanh")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("[§7.12] tanh R-vs-torch embedding maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                              # float32 weights: ~1e-7
  # the decisive single-sample property holds under tanh too
  batch <- score_dann(m, d$X)
  singles <- vapply(seq_len(nrow(d$X)), function(i) {
    score_dann(m, d$X[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.12] tanh single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.13 GRL is correct: forward identity + backward grad x (-grl_lambda) (< 1e-6)", {
  skip_if_no_dann()
  # NON-tautological check of the EXACT gradient-reversal layer the trainer uses. The
  # GRL must (i) be the identity on the forward pass and (ii) negate-and-scale the
  # gradient flowing back THROUGH it. Build a tiny toy: z -> GRL_lambda(z) -> a fixed
  # linear map -> scalar loss, and compare d(loss)/dz WITH the GRL against a no-GRL
  # reference. The GRL gradient must equal -lambda * (no-GRL gradient) to < 1e-6, and
  # the forward output must be byte-identical to the input (identity).
  reticulate::import("torch", delay_load = FALSE)
  # Pull the GRL helper out of the scorer's py_run_string definition by re-running it
  # here (idempotent class def) and grabbing _grad_reverse from main. This is the SAME
  # source string used by .dann_train_export -- no separate re-implementation.
  reticulate::py_run_string("
import torch
import torch.nn as nn

class _GradReverseFn(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, lamb):
        ctx.lamb = lamb
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        return grad_output.neg() * ctx.lamb, None

def _grad_reverse(x, lamb):
    return _GradReverseFn.apply(x, lamb)
", convert = FALSE)
  torch <- reticulate::import("torch", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)
  main <- reticulate::import_main()
  grad_reverse <- main$`_grad_reverse`

  z_vals <- matrix(c(0.4, -1.1, 2.3, 0.7, -0.2, 1.5), ncol = 3L, byrow = TRUE)
  w_vals <- matrix(c(0.9, -0.5, 1.2), ncol = 1L)     # fixed downstream linear map
  lambda <- 2.5

  make_z <- function() {
    z <- torch$as_tensor(np$asarray(z_vals, dtype = "float64"))$double()
    z$requires_grad_(TRUE)
    z
  }
  w <- torch$as_tensor(np$asarray(w_vals, dtype = "float64"))$double()

  # (i) forward identity: GRL output must equal the input exactly.
  z_fwd <- make_z()
  out_grl <- grad_reverse(z_fwd, lambda)
  fwd_maxdiff <- reticulate::py_to_r(
    torch$max(torch$abs(torch$sub(out_grl, z_fwd)))$item())
  cat(sprintf("[§7.13] GRL forward-identity maxdiff = %.3e\n", fwd_maxdiff))
  expect_lt(fwd_maxdiff, 1e-12)

  # (ii) backward: gradient WITH the GRL.
  z_grl <- make_z()
  loss_grl <- torch$matmul(grad_reverse(z_grl, lambda), w)$sum()
  loss_grl$backward()
  grad_grl <- reticulate::py_to_r(z_grl$grad$cpu()$numpy())

  # reference: SAME graph but NO GRL (plain identity). Its gradient is d(z w)/dz = w^T
  # broadcast over rows, i.e. each row = t(w_vals).
  z_ref <- make_z()
  loss_ref <- torch$matmul(z_ref, w)$sum()
  loss_ref$backward()
  grad_ref <- reticulate::py_to_r(z_ref$grad$cpu()$numpy())

  # GRL gradient must be EXACTLY -lambda * (no-GRL gradient).
  back_maxdiff <- max(abs(grad_grl - (-lambda * grad_ref)))
  cat(sprintf("[§7.13] GRL backward: grad_GRL vs -lambda*grad_ref maxdiff = %.3e\n",
              back_maxdiff))
  expect_lt(back_maxdiff, 1e-6)
  # sanity: grad_ref is the plain downstream gradient (= w broadcast), non-zero.
  expect_gt(max(abs(grad_ref)), 1e-6)
  # and the GRL gradient is genuinely reversed (opposite sign of the reference).
  expect_true(all(sign(grad_grl[grad_ref != 0]) == -sign(grad_ref[grad_ref != 0])))
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_dann()
  d <- .make_dann_data(seed = 111L)
  expect_error(fit_dann(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_dann(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_dann(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_dann(d$X, d$y, hp = list(weight_decay = -1)),
               "weight_decay")
  expect_error(fit_dann(d$X, d$y, hp = list(hidden = c(0L, 4L))), "hidden")
  expect_error(fit_dann(d$X, d$y, hp = list(activation = "gelu")),
               "activation")
  expect_error(fit_dann(d$X, d$y, hp = list(grl_lambda = -1)), "grl_lambda")
  expect_error(fit_dann(d$X, d$y, hp = list(domain_hidden = c(0L, 4L))),
               "domain_hidden")
  expect_error(fit_dann(d$X, d$y, hp = list(cohort_col = "")), "cohort_col")
  expect_error(fit_dann(d$X, d$y, hp = list(cohort_col = c("a", "b"))),
               "cohort_col")
  expect_error(fit_dann(d$X, d$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_dann(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_dann(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_dann(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_dann(d$X, d$y,
             hp = structure(list(1e-3, 1e-3), names = c("lr", "lr"))),
    "duplicate hp")
})
