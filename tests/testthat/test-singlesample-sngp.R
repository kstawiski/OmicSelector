# unc-sngp trains a small spectral-normalized, LayerNorm MLP encoder via
# reticulate-python torch at FIT (torch is confined to fit; the fixed RFF map, the
# Laplace posterior, and scoring are all pure R). Tests that fit are skipped when
# reticulate or the venv torch is unavailable, so the package suite stays green on
# hosts without the venv. Plain torch training needs NO brew-OpenSSL bridge and NO
# HF download, so the skip-guard only probes reticulate + `import torch`.
skip_if_no_sngp <- function() {
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

# Synthetic generator with a CASE-ONLY location shift in a feature block. The
# first `k` of `p` log-abundance features are elevated by `shift` for case
# samples; after exp() and per-sample rCLR centring those features carry the
# discriminative signal the encoder learns. Abundances are exp() so the input is
# strictly positive (compositional).
.make_sngp_data <- function(n = 80L, p = 30L, k = 8L, shift = 1.4,
                            sd = 0.5, seed = 11L) {
  set.seed(seed)
  y <- rep(c(0, 1), each = n / 2L)
  features <- paste0("miR-", sprintf("%03d", seq_len(p)))
  L <- matrix(stats::rnorm(n * p, mean = 4, sd = sd), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), features))
  case <- which(y == 1L)
  L[case, seq_len(k)] <- L[case, seq_len(k)] + shift
  X <- exp(L)
  list(X = X, y = y, block = features[seq_len(k)])
}

.auc_mw_sngp <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# Rebuild a torch forward of the encoder from the EXPORTED (effective post-SN)
# weights + LayerNorm params, to compare against the pure-R forward. Uses plain
# nn$Linear with weight set to W_eff (so torch reproduces x %*% t(W_eff) + b
# exactly, the same forward the spectral-norm module performed), plus nn$LayerNorm
# with the exported gamma/beta and activation after every layer except the last.
.torch_encoder_forward <- function(m, Z_in, activation) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  L <- length(m$weights)
  lins <- list(); lns <- list()
  for (l in seq_len(L)) {
    W <- m$weights[[l]]$W; b <- m$weights[[l]]$b
    lin <- nn$Linear(ncol(W), nrow(W))
    lin$weight$data <- torch$as_tensor(np$asarray(W, dtype = "float32"))$float()
    lin$bias$data   <- torch$as_tensor(np$asarray(b, dtype = "float32"))$float()
    lin$eval()
    lins[[l]] <- lin
    ln <- nn$LayerNorm(nrow(W))
    ln$weight$data <- torch$as_tensor(
      np$asarray(m$layernorm[[l]]$gamma, dtype = "float32"))$float()
    ln$bias$data <- torch$as_tensor(
      np$asarray(m$layernorm[[l]]$beta, dtype = "float32"))$float()
    ln$eval()
    lns[[l]] <- ln
  }
  X_py <- torch$as_tensor(np$asarray(Z_in, dtype = "float32"))$float()
  H <- X_py
  for (l in seq_len(L)) {
    H <- lins[[l]](H)
    H <- lns[[l]](H)
    if (l < L) {
      H <- if (identical(activation, "tanh")) torch$tanh(H) else torch$relu(H)
    }
  }
  reticulate::py_to_r(H$detach()$to(torch$float64)$numpy())
}


test_that("fit returns a well-formed unc_sngp_model", {
  skip_if_no_sngp()
  d <- .make_sngp_data()
  m <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "unc_sngp_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$embedding_dim, 32L)             # last hidden = embedding dim
  expect_length(m$weights, 2L)                   # 2 SN-Linear layers (encoder)
  expect_equal(dim(m$weights[[1]]$W), c(64L, ncol(d$X)))
  expect_equal(dim(m$weights[[2]]$W), c(32L, 64L))
  expect_length(m$layernorm, 2L)
  expect_length(m$layernorm[[1]]$gamma, 64L)
  expect_length(m$layernorm[[2]]$gamma, 32L)
  expect_equal(dim(m$W_rff), c(128L, 32L))       # rff_dim x embedding_dim
  expect_length(m$b_rff, 128L)
  expect_length(m$beta, 128L)
  expect_equal(dim(m$Sigma), c(128L, 128L))
  expect_equal(m$rff_scale, sqrt(2 / 128L))
  s <- score_unc_sngp(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("§7.1 R encoder forward matches the torch forward (LayerNorm + SN-baked weights, ~1e-6)", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 21L)
  m <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 50L, device = "cpu", seed = 7L,
                                        sn_coeff = 0.95))
  Z_in <- .unc_sngp_rclr_matrix(d$X)
  torch_emb <- .torch_encoder_forward(m, Z_in, "relu")
  r_emb <- .unc_sngp_forward(Z_in, m$weights, m$layernorm, "relu")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("\n[§7.1] R-vs-torch encoder (LayerNorm+SN) maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                            # float32 weights: ~1e-7
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 31L)
  m <- fit_unc_sngp(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  batch <- score_unc_sngp(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_unc_sngp(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 41L)
  m <- fit_unc_sngp(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  base <- score_unc_sngp(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_unc_sngp(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (maxdiff ~0)", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 44L)
  m <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 41L, 70L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_unc_sngp(m, row * s) -
                                score_unc_sngp(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_unc_sngp(m, scaled) -
                                score_unc_sngp(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-7)
})

test_that("§7.5 held-out AUC on a strongly planted shift is high (> 0.7)", {
  skip_if_no_sngp()
  d <- .make_sngp_data(n = 160L, seed = 51L)
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  tr <- c(utils::head(ctrl, 60L), utils::head(case, 60L))
  te <- c(utils::tail(ctrl, 20L), utils::tail(case, 20L))
  m <- fit_unc_sngp(d$X[tr, ], d$y[tr], hp = list(epochs = 200L))
  s <- score_unc_sngp(m, d$X[te, ])
  auc <- .auc_mw_sngp(d$y[te], s)
  cat(sprintf("[§7.5] held-out AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.6 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 61L)
  m <- fit_unc_sngp(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  score_fun <- function(model, X, meta) score_unc_sngp(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .unc_sngp_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .unc_sngp_model_digest))

  # Canonical dispatch: post-integration the manifest carries the unc-sngp row and
  # score_unc_sngp is in the namespace; pre-integration (staged self-test) clone a
  # one-row roster from a committed within-discriminator row and register a thin
  # adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_unc_sngp(m, X_test)
  roster <- singlesample_method_roster()
  if (!"unc-sngp" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "unc-sngp"
    tmpl$fit_fn <- "fit_unc_sngp"
    tmpl$score_fn <- "score_unc_sngp"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_unc_sngp", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "unc-sngp",
      function(model, X, meta) score_unc_sngp(model, X, meta))
  }
  via <- singlesample_score_call("unc-sngp", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.7 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 71L)
  m <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_unc_sngp(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_unc_sngp(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_unc_sngp(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be embedded and scored, NOT floored.
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .unc_sngp_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  emb <- .unc_sngp_forward(matrix(z_flat, nrow = 1L), m$weights, m$layernorm,
                           m$activation)
  phi <- as.numeric(.unc_sngp_rff(emb, m$W_rff, m$b_rff, m$rff_scale))
  manual <- .unc_sngp_score_one(phi, m$beta, m$Sigma)
  s_flat <- score_unc_sngp(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (manual %.6f)\n",
              s_flat, manual))
  expect_equal(s_flat, manual, tolerance = 1e-12)

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_unc_sngp(m, empty), 0)
})

test_that("§7.8 determinism: two seed=42 CPU fits give identical digests", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 81L)
  m1 <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  # Perturb the GLOBAL torch RNG to a different state before the second fit. A
  # correctly seeded fit re-seeds every stochastic step (encoder init INCLUDED,
  # plus the dedicated RFF generator) from hp$seed, so m2 must still be
  # bit-identical to m1 despite the perturbed entry state. This makes the
  # in-process test catch the seed-after-module-construction class of bug (init
  # reading the unseeded process-entry RNG), which a plain back-to-back two-fit
  # test masks via the fit's on.exit RNG restore -- and which would otherwise
  # surface only as cross-session non-reproducibility.
  torch_pert <- reticulate::import("torch", delay_load = FALSE)
  torch_pert$manual_seed(20260610L)
  invisible(torch_pert$randn(list(64L, 64L)))
  m2 <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .unc_sngp_model_digest(m1)
  dg2 <- .unc_sngp_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s\n", dg1, dg2))
  expect_identical(dg1, dg2)
  # Scores also bit-identical across the two fits.
  s1 <- score_unc_sngp(m1, d$X[1:10, ])
  s2 <- score_unc_sngp(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.9 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))

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
  invisible(score_unc_sngp(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "unc_sngp_model")
})

test_that("§7.10 independent GP-logit recompute (incl. LayerNorm + RFF cos) bit-matches the scorer", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 101L)
  m <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_unc_sngp(m, X_test)
  X_use <- .unc_sngp_align(X_test, m$feature_universe)
  beta <- m$beta; Sigma <- m$Sigma; pi8 <- pi / 8
  L <- length(m$weights)
  act <- if (identical(m$activation, "tanh")) tanh else function(x) pmax(x, 0)
  worst <- 0; worst_ln <- 0; worst_rff <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .unc_sngp_rclr(X_use[i, ])
    # --- from-scratch encoder forward (manual LayerNorm, non-circular) ----------
    H <- matrix(z, nrow = 1L)
    for (l in seq_len(L)) {
      W <- m$weights[[l]]$W; b <- m$weights[[l]]$b
      h <- as.numeric(H %*% t(W)) + b
      g <- m$layernorm[[l]]$gamma; be <- m$layernorm[[l]]$beta
      eps <- m$layernorm[[l]]$eps
      mu <- mean(h); v <- mean((h - mu)^2)         # biased variance (1/N)
      hln <- g * (h - mu) / sqrt(v + eps) + be
      if (l < L) hln <- act(hln)
      H <- matrix(hln, nrow = 1L)
    }
    emb <- as.numeric(H)
    # cross-check the encoder forward against .unc_sngp_forward
    emb_ref <- as.numeric(.unc_sngp_forward(matrix(z, nrow = 1L), m$weights,
                                            m$layernorm, m$activation))
    worst_ln <- max(worst_ln, max(abs(emb - emb_ref)))
    # --- from-scratch RFF cos map (non-circular) -------------------------------
    phi <- sqrt(2 / m$rff_dim) * cos(as.numeric(m$W_rff %*% emb) + m$b_rff)
    phi_ref <- as.numeric(.unc_sngp_rff(matrix(emb, nrow = 1L), m$W_rff, m$b_rff,
                                        m$rff_scale))
    worst_rff <- max(worst_rff, max(abs(phi - phi_ref)))
    # --- from-scratch mean-field GP logit (non-circular) -----------------------
    mean_logit <- sum(phi * beta)
    var_term <- as.numeric(t(phi) %*% (Sigma %*% phi))
    sc <- mean_logit / sqrt(1 + pi8 * var_term)
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.10] independent GP-logit recompute maxdiff = %.3e ; LayerNorm-forward maxdiff = %.3e ; RFF-cos maxdiff = %.3e\n",
              worst, worst_ln, worst_rff))
  expect_lt(worst, 1e-9)
  expect_lt(worst_ln, 1e-12)
  expect_lt(worst_rff, 1e-12)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_sngp()
  d <- .make_sngp_data(seed = 111L)
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(weight_decay = -1)),
               "weight_decay")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(hidden = c(0L, 4L))), "hidden")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(activation = "gelu")),
               "activation")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(sn_coeff = 0)), "sn_coeff")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(rff_dim = 0)), "rff_dim")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(rff_dim = 1.5)), "rff_dim")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(rff_ls = 0)), "rff_ls")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(ridge = 0)), "ridge")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_unc_sngp(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_unc_sngp(d$X, d$y,
                 hp = structure(list(1e-3, 1e-3),
                                names = c("lr", "lr"))),
    "duplicate hp")
})

test_that("§7.11 activation = 'tanh': R encoder forward matches torch + invariants hold", {
  skip_if_no_sngp()
  d <- .make_sngp_data(n = 160L, seed = 131L)
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  tr <- c(utils::head(ctrl, 60L), utils::head(case, 60L))
  te <- c(utils::tail(ctrl, 20L), utils::tail(case, 20L))
  m <- fit_unc_sngp(d$X[tr, ], d$y[tr],
                    hp = list(activation = "tanh", epochs = 200L,
                              device = "cpu", seed = 7L))
  expect_identical(m$activation, "tanh")
  Z_in <- .unc_sngp_rclr_matrix(d$X[tr, ])
  torch_emb <- .torch_encoder_forward(m, Z_in, "tanh")
  r_emb <- .unc_sngp_forward(Z_in, m$weights, m$layernorm, "tanh")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("\n[§7.11] tanh R-vs-torch encoder maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)
  # single-row == batch (decisive single-sample, exact) on the tanh forward path
  s_batch <- score_unc_sngp(m, d$X[te, ])
  singles <- vapply(seq_along(te), function(i)
    score_unc_sngp(m, d$X[te[i], , drop = FALSE]), numeric(1L))
  expect_equal(max(abs(s_batch - singles)), 0, tolerance = 1e-6)
  expect_true(all(is.finite(s_batch)))
  expect_gt(.auc_mw_sngp(d$y[te], s_batch), 0.7)        # tanh encoder discriminates
})

test_that("§7.12 mean-field denominator is distance-aware (far-from-data specimens are less certain)", {
  skip_if_no_sngp()
  # The mean-field denominator sqrt(1 + (pi/8) phi' Sigma phi) is the SNGP
  # distance awareness: a far-from-training-manifold specimen has a larger RFF-GP
  # predictive variance phi' Sigma phi (and so a |logit| shrunk vs the posterior
  # mean). The faithful, low-variance statement is over the TRAINING-ROW variance
  # DISTRIBUTION (a single train row can sit at either tail), so we compare a
  # genuinely out-of-distribution specimen against the train-variance mean / max.
  d <- .make_sngp_data(n = 120L, seed = 151L)
  m <- fit_unc_sngp(d$X, d$y, hp = list(epochs = 120L))
  var_of <- function(row) {
    z <- .unc_sngp_rclr(.unc_sngp_align(row, m$feature_universe)[1, ])
    emb <- .unc_sngp_forward(matrix(z, nrow = 1L), m$weights, m$layernorm,
                             m$activation)
    phi <- as.numeric(.unc_sngp_rff(emb, m$W_rff, m$b_rff, m$rff_scale))
    as.numeric(t(phi) %*% (m$Sigma %*% phi))
  }
  v_train <- vapply(seq_len(nrow(d$X)), function(i)
    var_of(d$X[i, , drop = FALSE]), numeric(1L))
  # A genuinely out-of-distribution specimen: an extreme one-feature-dominant
  # composition (deterministic, no RNG) that is far from every block-structured
  # training row. Its predictive variance exceeds the WHOLE training distribution.
  far <- matrix(c(1e6, rep(1, length(m$feature_universe) - 1L)), nrow = 1L,
                dimnames = list("far", m$feature_universe))
  v_far <- var_of(far)
  cat(sprintf("[§7.12] RFF-GP predictive variance: OOD far = %.4f ; train mean = %.4f ; train max = %.4f\n",
              v_far, mean(v_train), max(v_train)))
  expect_gt(v_far, mean(v_train))                      # OOD is less certain than typical train
  expect_gt(v_far, max(v_train))                       # OOD exceeds the whole train distribution
  # Averaged over many random far draws, the predictive variance is systematically
  # higher than the training mean (the distance-awareness signal, not a single-point
  # fluke).
  set.seed(11)
  v_far_many <- vapply(seq_len(50L), function(i)
    var_of(matrix(exp(stats::rnorm(length(m$feature_universe), 0, 6)), nrow = 1L,
                  dimnames = list("z", m$feature_universe))), numeric(1L))
  cat(sprintf("[§7.12] mean OOD variance over 50 random far draws = %.4f (> train mean %.4f)\n",
              mean(v_far_many), mean(v_train)))
  expect_gt(mean(v_far_many), mean(v_train))
  # The score still respects the variance denominator (finite, shrunk).
  s_far <- score_unc_sngp(m, far)
  expect_true(is.finite(s_far))
})

test_that("§7.13 sn_coeff != 1: exported forward reproduces the TRAINING forward (bias scaled by c)", {
  skip_if_no_sngp()
  # Closes the §7.1 blind spot (§7.1 compares R against a reconstruction from the
  # EXPORTED params, so it cannot detect whether the export reproduces TRAINING).
  # The training forward scales the WHOLE linear output by c: c*(H W_SN^T + b), so
  # the effective bias is c*b -- the export must bake b_eff = c*b, not b. Here we
  # build a spectral-normed layer, take its ACTUAL c-scaled training forward, and
  # check the pure-R `.unc_sngp_forward` on the export convention (W_eff=c*W_SN,
  # b_eff=c*b) reproduces it; and that the UNSCALED bias does NOT (regression lock).
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  torch$manual_seed(5L)
  c_sn <- 0.5; p <- 16L; h <- 8L; n <- 6L
  sn <- nn$utils$parametrizations$spectral_norm
  lin <- sn(nn$Linear(p, h)); ln <- nn$LayerNorm(h)
  Xpy <- torch$randn(list(n, p))
  # ACTUAL training forward (mirrors the fit forward: c * lin(X) then LayerNorm):
  H_train <- ln(torch$mul(lin(Xpy), c_sn))
  py_out <- reticulate::py_to_r(H_train$detach()$to(torch$float64)$numpy())
  X  <- reticulate::py_to_r(Xpy$to(torch$float64)$numpy())
  Wsn <- reticulate::py_to_r(lin$weight$detach()$to(torch$float64)$numpy())
  b   <- as.numeric(reticulate::py_to_r(lin$bias$detach()$to(torch$float64)$numpy()))
  g   <- as.numeric(reticulate::py_to_r(ln$weight$detach()$to(torch$float64)$numpy()))
  be  <- as.numeric(reticulate::py_to_r(ln$bias$detach()$to(torch$float64)$numpy()))
  epsv <- as.numeric(reticulate::py_to_r(ln$eps))
  lnp <- list(list(gamma = g, beta = be, eps = epsv))
  # export convention used by fit_unc_sngp: W_eff = c*W_SN, b_eff = c*b
  r_scaled <- .unc_sngp_forward(X, list(list(W = c_sn * Wsn, b = c_sn * b)), lnp, "relu")
  d_scaled <- max(abs(py_out - r_scaled))
  # the buggy unscaled-bias export: W_eff = c*W_SN, b_eff = b
  r_unscaled <- .unc_sngp_forward(X, list(list(W = c_sn * Wsn, b = b)), lnp, "relu")
  d_unscaled <- max(abs(py_out - r_unscaled))
  cat(sprintf("[§7.13] sn_coeff=0.5  scaled-bias maxdiff = %.3e ; unscaled-bias maxdiff = %.3e\n",
              d_scaled, d_unscaled))
  # scaled bias reproduces training to float32 round-off (looser than the trained
  # §7.1 ~1e-7 because this UNTRAINED random spectral-norm+LayerNorm layer can sit
  # near the LayerNorm variance floor, amplifying float32 error); the unscaled bias
  # is off by O(0.1-1) -- a 3+ order-of-magnitude separation either way.
  expect_lt(d_scaled, 1e-3)                            # scaled bias reproduces training
  expect_gt(d_unscaled, 1e-2)                          # unscaled bias is a different net (the bug)
})
