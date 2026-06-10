
# cvae trains a small BN-free CLASS-conditional VAE (encoder q(z|x,y) + decoder
# p(x|z,y)) via reticulate-python torch by the reparameterized ELBO at FIT (torch is
# confined to fit; scoring is pure R, deterministic z = mu). Tests that fit are
# skipped when reticulate or the venv torch is unavailable, so the package suite
# stays green on hosts without the venv. Plain torch training needs NO brew-OpenSSL
# bridge and NO HF download, so the skip-guard only probes reticulate + `import
# torch`.
skip_if_no_cvae <- function() {
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

# Synthetic generator with a CASE-ONLY location shift in a feature block. The first
# `k` of `p` log-abundance features are elevated by `shift` for case samples; after
# exp() and per-sample rCLR centring those features carry the discriminative signal
# the class-conditional VAE learns. Abundances are exp() so the input is strictly
# positive (compositional).
.make_cvae_data <- function(n = 80L, p = 30L, k = 8L, shift = 1.4,
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

.auc_mw_cvae <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# Rebuild a torch Sequential trunk + bare-Linear heads from EXPORTED float32 weights
# so §7.1 can compare the pure-R forwards/ELBO against torch on the SAME frozen
# weights the scorer uses. `act` is the torch activation constructor.
.cvae_build_torch <- function(torch, nn, np, weights, act) {
  mk_lin <- function(layer) {
    W <- layer$W; b <- layer$b
    lin <- nn$Linear(ncol(W), nrow(W))
    lin$weight$data <- torch$as_tensor(np$asarray(W, dtype = "float32"))$float()
    lin$bias$data   <- torch$as_tensor(np$asarray(b, dtype = "float32"))$float()
    lin
  }
  mk_trunk <- function(trunk) {
    layers <- list()
    for (l in seq_along(trunk)) {
      layers[[length(layers) + 1L]] <- mk_lin(trunk[[l]])
      layers[[length(layers) + 1L]] <- act()
    }
    do.call(nn$Sequential, layers)
  }
  enc_trunk <- mk_trunk(weights$enc_trunk); enc_trunk$eval()
  enc_mu <- mk_lin(weights$enc_mu); enc_mu$eval()
  enc_logvar <- mk_lin(weights$enc_logvar); enc_logvar$eval()
  dec_trunk <- mk_trunk(weights$dec_trunk); dec_trunk$eval()
  dec_xhat <- mk_lin(weights$dec_xhat); dec_xhat$eval()
  list(enc_trunk = enc_trunk, enc_mu = enc_mu, enc_logvar = enc_logvar,
       dec_trunk = dec_trunk, dec_xhat = dec_xhat)
}


test_that("fit returns a well-formed cvae_model", {
  skip_if_no_cvae()
  d <- .make_cvae_data()
  m <- fit_cvae(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "cvae_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$latent_dim, 16L)
  expect_length(m$weights$enc_trunk, 2L)                 # 2 encoder trunk layers
  expect_length(m$weights$dec_trunk, 2L)                 # 2 decoder trunk layers
  # encoder trunk first layer: in = p + 2 (one-hot), out = 64.
  expect_equal(dim(m$weights$enc_trunk[[1]]$W), c(64L, ncol(d$X) + 2L))
  expect_equal(dim(m$weights$enc_trunk[[2]]$W), c(32L, 64L))
  expect_equal(dim(m$weights$enc_mu$W), c(16L, 32L))     # mu-head: 32 -> L
  expect_equal(dim(m$weights$enc_logvar$W), c(16L, 32L)) # logvar-head: 32 -> L
  # decoder trunk first layer: in = L + 2 (one-hot), out = 32.
  expect_equal(dim(m$weights$dec_trunk[[1]]$W), c(32L, 16L + 2L))
  expect_equal(dim(m$weights$dec_trunk[[2]]$W), c(64L, 32L))
  expect_equal(dim(m$weights$dec_xhat$W), c(ncol(d$X), 64L))  # xhat-head: 64 -> p
  expect_true(is.numeric(m$sigma2) && length(m$sigma2) == 1L && m$sigma2 > 0)
  s <- score_cvae(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("§7.1 R encode/decode/ELBO match the torch forwards (~1e-6)", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 21L)
  m <- fit_cvae(d$X, d$y, hp = list(epochs = 50L, device = "cpu", seed = 7L))
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  act <- if (identical(m$activation, "tanh")) nn$Tanh else nn$ReLU
  tm <- .cvae_build_torch(torch, nn, np, m$weights, act)

  Z_in <- .cvae_rclr_matrix(d$X)
  p <- ncol(Z_in); L <- m$latent_dim; n <- nrow(Z_in)
  # Evaluate under BOTH classes and compare encoder (mu, logvar), decoder (xhat),
  # and the per-class ELBO (with z = mu, deterministic) against torch.
  worst_mu <- 0; worst_lv <- 0; worst_xh <- 0; worst_elbo <- 0
  half_log_2pi <- 0.5 * log(2 * pi)
  for (cl in c(0L, 1L)) {
    oh <- if (cl == 0L) c(1, 0) else c(0, 1)
    OH <- matrix(rep(oh, each = n), nrow = n)
    # ---- torch forwards ----
    Xt <- torch$as_tensor(np$asarray(Z_in, dtype = "float32"))$float()
    Yt <- torch$as_tensor(np$asarray(OH, dtype = "float32"))$float()
    enc_h <- tm$enc_trunk(torch$cat(list(Xt, Yt), 1L))
    mu_t <- tm$enc_mu(enc_h)
    lv_t <- tm$enc_logvar(enc_h)
    xhat_t <- tm$dec_xhat(tm$dec_trunk(torch$cat(list(mu_t, Yt), 1L)))   # z = mu
    mu_tr <- reticulate::py_to_r(mu_t$detach()$to(torch$float64)$numpy())
    lv_tr <- reticulate::py_to_r(lv_t$detach()$to(torch$float64)$numpy())
    xhat_tr <- reticulate::py_to_r(xhat_t$detach()$to(torch$float64)$numpy())
    # torch ELBO per row (z = mu): logN - beta*KL
    sig2 <- m$sigma2
    logN_t <- -0.5 * (p * log(2 * pi * sig2) +
                        rowSums((Z_in - xhat_tr)^2) / sig2)
    kl_t <- 0.5 * rowSums(mu_tr^2 + exp(lv_tr) - lv_tr - 1)
    elbo_t <- logN_t - m$kl_beta * kl_t
    # ---- pure-R forwards (the scorer's helpers) ----
    M_enc <- cbind(Z_in, OH)
    enc_r <- .cvae_encode(M_enc, m$weights, m$activation)
    M_dec <- cbind(enc_r$mu, OH)
    xhat_r <- .cvae_decode(M_dec, m$weights, m$activation)
    elbo_r <- vapply(seq_len(n), function(i)
      .cvae_elbo_class(Z_in[i, ], cl, m$weights, m$sigma2, m$kl_beta,
                       m$activation, p), numeric(1L))
    worst_mu <- max(worst_mu, max(abs(enc_r$mu - mu_tr)))
    worst_lv <- max(worst_lv, max(abs(enc_r$logvar - lv_tr)))
    worst_xh <- max(worst_xh, max(abs(xhat_r - xhat_tr)))
    worst_elbo <- max(worst_elbo, max(abs(elbo_r - elbo_t)))
  }
  cat(sprintf("\n[§7.1] R-vs-torch maxdiff: mu=%.3e logvar=%.3e xhat=%.3e ELBO=%.3e\n",
              worst_mu, worst_lv, worst_xh, worst_elbo))
  expect_lt(worst_mu, 1e-4)        # float32 weights: ~1e-6
  expect_lt(worst_lv, 1e-4)
  expect_lt(worst_xh, 1e-4)
  expect_lt(worst_elbo, 1e-4)
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 31L)
  m <- fit_cvae(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  batch <- score_cvae(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_cvae(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 41L)
  m <- fit_cvae(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  base <- score_cvae(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_cvae(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (maxdiff ~0)", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 44L)
  m <- fit_cvae(d$X, d$y, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 41L, 70L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_cvae(m, row * s) - score_cvae(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_cvae(m, scaled) - score_cvae(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-7)
})

test_that("§7.5 held-out AUC on a strongly planted shift is high (> 0.7)", {
  skip_if_no_cvae()
  d <- .make_cvae_data(n = 160L, seed = 51L)
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  tr <- c(utils::head(ctrl, 60L), utils::head(case, 60L))
  te <- c(utils::tail(ctrl, 20L), utils::tail(case, 20L))
  m <- fit_cvae(d$X[tr, ], d$y[tr], hp = list(epochs = 400L))
  s <- score_cvae(m, d$X[te, ])
  auc <- .auc_mw_cvae(d$y[te], s)
  cat(sprintf("[§7.5] held-out AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.6 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 61L)
  m <- fit_cvae(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  score_fun <- function(model, X, meta) score_cvae(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .cvae_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .cvae_model_digest))

  # Canonical dispatch: post-integration the manifest carries the cvae row and
  # score_cvae is in the namespace; pre-integration (staged self-test) clone a
  # one-row roster from a committed within-discriminator row and register a thin
  # adapter so the same dispatch wiring is exercised. Either way the dispatched score
  # must equal the direct call.
  direct <- score_cvae(m, X_test)
  roster <- singlesample_method_roster()
  if (!"cvae" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "cvae"
    tmpl$fit_fn <- "fit_cvae"
    tmpl$score_fn <- "score_cvae"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_cvae", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "cvae", function(model, X, meta) score_cvae(model, X, meta))
  }
  via <- singlesample_score_call("cvae", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.7 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 71L)
  m <- fit_cvae(d$X, d$y, hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_cvae(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_cvae(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_cvae(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be encoded/decoded and scored, NOT floored.
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .cvae_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  manual <- .cvae_score_one(z_flat, m$weights, m$sigma2, m$kl_beta,
                            m$activation, p)
  s_flat <- score_cvae(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (manual %.6f)\n",
              s_flat, manual))
  expect_equal(s_flat, manual, tolerance = 1e-12)
  # The flat rCLR (origin) is generically encoded/decoded to a nonzero ELBO
  # difference, so the flat score is genuinely computed (not floored 0).

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_cvae(m, empty), 0)
})

test_that("§7.8 determinism: two seed=42 CPU fits give identical digests (inter-fit RNG perturbed)", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 81L)
  m1 <- fit_cvae(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))

  # Perturb the GLOBAL torch RNG BETWEEN the two seed=42 fits so a seed-before-build
  # regression (constructing a stochastic module before manual_seed) is actually
  # caught: without seed-before-build, the second fit's weight init would depend on
  # this perturbed generator state and the digests would differ.
  torch_pert <- reticulate::import("torch", delay_load = FALSE)
  torch_pert$manual_seed(99L)
  invisible(torch_pert$randn(list(64L, 64L)))

  m2 <- fit_cvae(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .cvae_model_digest(m1)
  dg2 <- .cvae_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s\n", dg1, dg2))
  expect_identical(dg1, dg2)
  # Scores also bit-identical across the two fits.
  s1 <- score_cvae(m1, d$X[1:10, ])
  s2 <- score_cvae(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.9 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_cvae(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))

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

  # Scoring consumes no RNG (z = mu is deterministic): .Random.seed unchanged.
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  invisible(score_cvae(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_cvae(d$X, d$y, hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "cvae_model")
})

test_that("§7.10 independent from-scratch ELBO_case-ELBO_control recompute bit-matches the scorer", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 101L)
  m <- fit_cvae(d$X, d$y, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_cvae(m, X_test)

  # From-scratch recompute: a SELF-CONTAINED base-R encoder/decoder forward + ELBO,
  # NOT calling any of the scorer's private helpers (.cvae_encode/.cvae_decode/
  # .cvae_elbo_class/.cvae_score_one). Reimplements the MLP forwards inline.
  act_r <- if (identical(m$activation, "tanh")) tanh else function(x) pmax(x, 0)
  fwd_trunk <- function(v, trunk) {
    h <- v
    for (l in seq_along(trunk)) {
      h <- as.numeric(trunk[[l]]$W %*% h) + trunk[[l]]$b
      h <- act_r(h)
    }
    h
  }
  lin_head <- function(h, head) as.numeric(head$W %*% h) + head$b
  p <- length(m$feature_universe); sig2 <- m$sigma2; beta <- m$kl_beta
  elbo_c <- function(x, oh) {
    enc_h <- fwd_trunk(c(x, oh), m$weights$enc_trunk)
    mu_c <- lin_head(enc_h, m$weights$enc_mu)
    lv_c <- lin_head(enc_h, m$weights$enc_logvar)
    dec_h <- fwd_trunk(c(mu_c, oh), m$weights$dec_trunk)
    xhat_c <- lin_head(dec_h, m$weights$dec_xhat)
    logN <- -0.5 * (p * log(2 * pi * sig2) + sum((x - xhat_c)^2) / sig2)
    kl <- 0.5 * sum(mu_c^2 + exp(lv_c) - lv_c - 1)
    logN - beta * kl
  }
  X_use <- .cvae_align(X_test, m$feature_universe)
  worst <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .cvae_rclr(X_use[i, ])
    sc <- elbo_c(z, c(0, 1)) - elbo_c(z, c(1, 0))     # case - control
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.10] independent ELBO-ratio recompute maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-12)
})

test_that("§7.11 COUNTERFACTUAL correctness: toy ELBO vs first principles + planted sign", {
  skip_if_no_cvae()

  # (a) Hand-built tiny toy with KNOWN mu/logvar/xhat/sigma2 -> KNOWN ELBO, verifying
  # the recon-loglik + KL formulas from first principles (does NOT depend on torch
  # or on a fit). We construct a 1-hidden-layer encoder + decoder with hand-set
  # weights so that for a chosen input we can read off mu, logvar, xhat exactly, then
  # check .cvae_elbo_class against an independent closed-form ELBO.
  p <- 3L; L <- 2L
  # Identity-ish trunks: a single Linear layer (weight = I block, bias chosen) then
  # the bare heads. We use tanh-free relu activation; to make the trunk output
  # predictable we feed an input whose post-Linear value is positive so ReLU is the
  # identity on it.
  # ENCODER: input dim p+2 = 5, one trunk layer 5 -> 4 (relu), heads 4 -> L.
  # Choose trunk W = first-4-rows selector with +10 bias so relu = identity (>0).
  enc_trunk_W <- matrix(0, nrow = 4L, ncol = p + 2L)
  for (j in 1:4) enc_trunk_W[j, j] <- 1
  enc_trunk_b <- rep(10, 4L)                              # push pre-act positive
  # mu-head: read trunk dims 1,2 (minus the +10 bias) -> mu = input[1:2].
  enc_mu_W <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  enc_mu_b <- c(-10, -10)
  # logvar-head: constant logvar = c(log(2), 0) regardless of input.
  enc_logvar_W <- matrix(0, nrow = L, ncol = 4L)
  enc_logvar_b <- c(log(2), 0)
  # DECODER: input dim L+2 = 4, one trunk layer 4 -> 4 (relu, +10 bias = identity),
  # xhat-head 4 -> p giving a fixed xhat regardless of z (weights 0, bias = xhat0).
  dec_trunk_W <- diag(4L); dec_trunk_b <- rep(10, 4L)
  dec_xhat_W <- matrix(0, nrow = p, ncol = 4L)
  xhat0 <- c(0.5, -0.2, 0.1)
  dec_xhat_b <- xhat0 + 10 * rowSums(matrix(0, nrow = p, ncol = 4L))  # = xhat0
  dec_xhat_b <- xhat0                                     # decoder bias = xhat0 (W=0)
  toy_w <- list(
    enc_trunk = list(list(W = enc_trunk_W, b = enc_trunk_b)),
    enc_mu = list(W = enc_mu_W, b = enc_mu_b),
    enc_logvar = list(W = enc_logvar_W, b = enc_logvar_b),
    dec_trunk = list(list(W = dec_trunk_W, b = dec_trunk_b)),
    dec_xhat = list(W = dec_xhat_W, b = dec_xhat_b)
  )
  sigma2 <- 0.7; beta <- 1.3
  x <- c(0.3, -0.4, 0.9)                                  # rCLR-like query (length p)
  # For class c = case (onehot c(0,1)): trunk-in = c(x, 0, 1); after selector trunk
  # (dims 1:4 of the 5-vector) + bias 10 + relu = c(x1,x2,x3,0)+10 (all > 0). Then
  # mu = (trunk[1],trunk[2]) - 10 = (x1, x2); logvar = (log2, 0); xhat = xhat0.
  mu_known <- c(x[1], x[2])
  logvar_known <- c(log(2), 0)
  xhat_known <- xhat0
  logN_known <- -0.5 * (p * log(2 * pi * sigma2) +
                          sum((x - xhat_known)^2) / sigma2)
  kl_known <- 0.5 * sum(mu_known^2 + exp(logvar_known) - logvar_known - 1)
  elbo_known <- logN_known - beta * kl_known

  elbo_fn <- .cvae_elbo_class(x, 1L, toy_w, sigma2, beta, "relu", p)
  cat(sprintf("[§7.11a] toy ELBO: closed-form %.10f vs scorer %.10f (diff %.2e)\n",
              elbo_known, elbo_fn, abs(elbo_known - elbo_fn)))
  expect_lt(abs(elbo_known - elbo_fn), 1e-10)
  # Cross-check the intermediate forwards too (mu, logvar, xhat) via the helpers.
  oh <- c(0, 1)
  enc <- .cvae_encode(matrix(c(x, oh), nrow = 1L), toy_w, "relu")
  expect_lt(max(abs(as.numeric(enc$mu) - mu_known)), 1e-12)
  expect_lt(max(abs(as.numeric(enc$logvar) - logvar_known)), 1e-12)
  xhat_r <- as.numeric(.cvae_decode(matrix(c(mu_known, oh), nrow = 1L), toy_w, "relu"))
  expect_lt(max(abs(xhat_r - xhat_known)), 1e-12)

  # (b) Planted-sign non-tautology: a model trained on a clear case-shift scores a
  # clearly case-like query > 0 and a clearly control-like query < 0.
  d <- .make_cvae_data(n = 160L, k = 10L, shift = 2.0, seed = 121L)
  m <- fit_cvae(d$X, d$y, hp = list(epochs = 400L))
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  # Prototype-like queries: per-feature class means (a clean case / control specimen).
  case_proto <- matrix(colMeans(d$X[case, , drop = FALSE]), nrow = 1L,
                       dimnames = list("case", colnames(d$X)))
  ctrl_proto <- matrix(colMeans(d$X[ctrl, , drop = FALSE]), nrow = 1L,
                       dimnames = list("ctrl", colnames(d$X)))
  s_case <- score_cvae(m, case_proto)
  s_ctrl <- score_cvae(m, ctrl_proto)
  cat(sprintf("[§7.11b] planted case-like score = %.4f (>0); control-like = %.4f (<0)\n",
              s_case, s_ctrl))
  expect_gt(s_case, 0)
  expect_lt(s_ctrl, 0)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_cvae()
  d <- .make_cvae_data(seed = 111L)
  expect_error(fit_cvae(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_cvae(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_cvae(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_cvae(d$X, d$y, hp = list(weight_decay = -1)), "weight_decay")
  expect_error(fit_cvae(d$X, d$y, hp = list(enc_hidden = c(0L, 4L))), "enc_hidden")
  expect_error(fit_cvae(d$X, d$y, hp = list(dec_hidden = c(0L, 4L))), "dec_hidden")
  expect_error(fit_cvae(d$X, d$y, hp = list(latent_dim = 0)), "latent_dim")
  expect_error(fit_cvae(d$X, d$y, hp = list(activation = "gelu")), "activation")
  expect_error(fit_cvae(d$X, d$y, hp = list(kl_beta = -1)), "kl_beta")
  expect_error(fit_cvae(d$X, d$y, hp = list(decoder_sigma2 = 0)), "decoder_sigma2")
  expect_error(fit_cvae(d$X, d$y, hp = list(decoder_sigma2 = -1)), "decoder_sigma2")
  expect_error(fit_cvae(d$X, d$y, hp = list(min_features = 0)), "min_features")
  expect_error(fit_cvae(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_cvae(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_cvae(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_cvae(d$X, d$y, hp = structure(list(1e-3, 1e-3), names = c("lr", "lr"))),
    "duplicate hp")
})
