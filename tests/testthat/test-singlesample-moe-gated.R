
# moe-gated trains a LayerNorm MLP encoder + E linear experts + a linear gate
# JOINTLY via reticulate-python torch at FIT (torch is confined to fit; the
# encoder forward, expert logits, gate softmax, and mixture are all pure R at
# score). Tests that fit are skipped when reticulate or the venv torch is
# unavailable, so the package suite stays green on hosts without the venv. Plain
# torch training needs NO brew-OpenSSL bridge and NO HF download, so the skip-guard
# only probes reticulate + `import torch`.
skip_if_no_moe_gated <- function() {
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
# discriminative signal the encoder + experts + gate learn. Abundances are exp()
# so the input is strictly positive (compositional).
.make_moe_data <- function(n = 80L, p = 30L, k = 8L, shift = 1.4,
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

.auc_mw_moe <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# Rebuild a torch forward of the full mixture (encoder LayerNorm MLP -> pre-act
# embedding -> gate-softmax-weighted expert mixture) from the EXPORTED weights, to
# compare against the pure-R mixture. Uses plain nn$Linear / nn$LayerNorm with the
# exported params; activation after every encoder layer except the last; gate
# softmax over the E experts (dim=1, per row); mixture = sum_E(g * h).
.torch_mixture_forward <- function(m, Z_in) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  L <- length(m$weights)
  lins <- list(); lns <- list()
  for (l in seq_len(L)) {
    W <- m$weights[[l]]$W; b <- m$weights[[l]]$b
    lin <- nn$Linear(ncol(W), nrow(W))
    lin$weight$data <- torch$as_tensor(np$asarray(W, dtype = "float32"))$float()
    lin$bias$data   <- torch$as_tensor(np$asarray(b, dtype = "float32"))$float()
    lin$eval(); lins[[l]] <- lin
    ln <- nn$LayerNorm(nrow(W))
    ln$weight$data <- torch$as_tensor(
      np$asarray(m$layernorm[[l]]$gamma, dtype = "float32"))$float()
    ln$bias$data <- torch$as_tensor(
      np$asarray(m$layernorm[[l]]$beta, dtype = "float32"))$float()
    ln$eval(); lns[[l]] <- ln
  }
  E <- m$n_experts
  experts <- nn$Linear(ncol(m$W_e), nrow(m$W_e))
  experts$weight$data <- torch$as_tensor(np$asarray(m$W_e, dtype = "float32"))$float()
  experts$bias$data   <- torch$as_tensor(np$asarray(m$b_e, dtype = "float32"))$float()
  experts$eval()
  gate <- nn$Linear(ncol(m$W_g), nrow(m$W_g))
  gate$weight$data <- torch$as_tensor(np$asarray(m$W_g, dtype = "float32"))$float()
  gate$bias$data   <- torch$as_tensor(np$asarray(m$b_g, dtype = "float32"))$float()
  gate$eval()
  act_fun <- if (identical(m$activation, "tanh")) torch$tanh else torch$relu

  X_py <- torch$as_tensor(np$asarray(Z_in, dtype = "float32"))$float()
  H <- X_py
  for (l in seq_len(L)) {
    H <- lins[[l]](H)
    H <- lns[[l]](H)
    if (l < L) H <- act_fun(H)
  }
  Hlog <- experts(H)
  Glog <- gate(H)
  Gw <- nn$functional$softmax(Glog, dim = 1L)
  s <- torch$sum(torch$mul(Gw, Hlog), dim = 1L)
  reticulate::py_to_r(s$detach()$to(torch$float64)$numpy())
}


test_that("fit returns a well-formed moe_gated_model", {
  skip_if_no_moe_gated()
  d <- .make_moe_data()
  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "moe_gated_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$embedding_dim, 32L)             # last hidden = embedding dim
  expect_length(m$weights, 2L)                   # 2 Linear layers (encoder)
  expect_equal(dim(m$weights[[1]]$W), c(64L, ncol(d$X)))
  expect_equal(dim(m$weights[[2]]$W), c(32L, 64L))
  expect_length(m$layernorm, 2L)
  expect_length(m$layernorm[[1]]$gamma, 64L)
  expect_length(m$layernorm[[2]]$gamma, 32L)
  expect_equal(m$n_experts, 4L)
  expect_equal(dim(m$W_e), c(4L, 32L))           # E x d expert weights
  expect_length(m$b_e, 4L)
  expect_equal(dim(m$W_g), c(4L, 32L))           # E x d gate weights
  expect_length(m$b_g, 4L)
  s <- score_moe_gated(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("CUDA fit places encoder, experts, and gate on the requested device", {
  skip_if_no_moe_gated()
  torch <- reticulate::import("torch", delay_load = FALSE)
  testthat::skip_if_not(
    isTRUE(reticulate::py_to_r(torch$cuda$is_available())),
    "CUDA is unavailable"
  )
  d <- .make_moe_data(n = 20L, p = 8L, k = 3L, seed = 19L)
  m <- fit_moe_gated(
    d$X, d$y,
    hp = list(hidden = c(8L, 4L), epochs = 2L, device = "cuda", seed = 19L)
  )
  expect_identical(m$device, "cuda")
  expect_true(all(is.finite(score_moe_gated(m, d$X))))
})

test_that("§7.1 R encoder+MoE forward matches the torch mixture (LayerNorm + weights baked, ~1e-6)", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 21L)
  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 50L, device = "cpu", seed = 7L))
  Z_in <- .moe_gated_rclr_matrix(d$X)
  torch_s <- .torch_mixture_forward(m, Z_in)
  emb <- .moe_gated_forward(Z_in, m$weights, m$layernorm, m$activation)
  r_s <- .moe_gated_mixture(emb, m$W_e, m$b_e, m$W_g, m$b_g)$score
  maxdiff <- max(abs(torch_s - r_s))
  cat(sprintf("\n[§7.1] R-vs-torch mixture (LayerNorm+experts+gate) maxdiff = %.3e\n",
              maxdiff))
  expect_lt(maxdiff, 1e-4)                            # float32 weights: ~1e-6
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 31L)
  m <- fit_moe_gated(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  batch <- score_moe_gated(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_moe_gated(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 41L)
  m <- fit_moe_gated(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  base <- score_moe_gated(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_moe_gated(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (maxdiff ~0)", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 44L)
  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 41L, 70L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_moe_gated(m, row * s) -
                                score_moe_gated(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_moe_gated(m, scaled) -
                                score_moe_gated(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-12)
})

test_that("§7.5 held-out AUC on a strongly planted shift is high (> 0.7)", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(n = 160L, seed = 51L)
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  tr <- c(utils::head(ctrl, 60L), utils::head(case, 60L))
  te <- c(utils::tail(ctrl, 20L), utils::tail(case, 20L))
  m <- fit_moe_gated(d$X[tr, ], d$y[tr], hp = list(epochs = 200L))
  s <- score_moe_gated(m, d$X[te, ])
  auc <- .auc_mw_moe(d$y[te], s)
  cat(sprintf("[§7.5] held-out AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.6 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 61L)
  m <- fit_moe_gated(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  score_fun <- function(model, X, meta) score_moe_gated(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .moe_gated_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .moe_gated_model_digest))

  # Canonical dispatch: post-integration the manifest carries the moe-gated row and
  # score_moe_gated is in the namespace; pre-integration (staged self-test) clone a
  # one-row roster from a committed within-discriminator row and register a thin
  # adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_moe_gated(m, X_test)
  roster <- singlesample_method_roster()
  if (!"moe-gated" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "moe-gated"
    tmpl$fit_fn <- "fit_moe_gated"
    tmpl$score_fn <- "score_moe_gated"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_moe_gated", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "moe-gated",
      function(model, X, meta) score_moe_gated(model, X, meta))
  }
  via <- singlesample_score_call("moe-gated", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.7 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 71L)
  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_moe_gated(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_moe_gated(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_moe_gated(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be embedded and scored, NOT floored.
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .moe_gated_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  emb <- .moe_gated_forward(matrix(z_flat, nrow = 1L), m$weights, m$layernorm,
                            m$activation)
  manual <- .moe_gated_score_one(as.numeric(emb), m$W_e, m$b_e, m$W_g, m$b_g)
  s_flat <- score_moe_gated(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (manual %.6f)\n",
              s_flat, manual))
  expect_equal(s_flat, manual, tolerance = 1e-12)
  # The flat embedding is the bias-propagated LayerNorm forward; the mixture logit
  # on it is generically nonzero, so the flat score is genuinely computed (not 0).

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_moe_gated(m, empty), 0)
})

test_that("§7.8 determinism: two seed=42 CPU fits give identical digests (inter-fit RNG perturbed)", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 81L)
  m1 <- fit_moe_gated(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  # Perturb the GLOBAL torch RNG to a different state before the second fit. A
  # correctly SEED-BEFORE-BUILD fit re-seeds the global torch generator from
  # hp$seed BEFORE constructing any stochastic module (encoder + LayerNorm +
  # experts + gate), so m2 must still be bit-identical to m1 despite the perturbed
  # entry state. This makes the in-process test catch the seed-after-module-
  # construction class of bug (init reading the unseeded process-entry RNG), which
  # a plain back-to-back two-fit test masks via the fit's on.exit RNG restore --
  # and which would otherwise surface only as cross-session non-reproducibility.
  torch_pert <- reticulate::import("torch", delay_load = FALSE)
  torch_pert$manual_seed(20260610L)
  invisible(torch_pert$randn(list(64L, 64L)))
  m2 <- fit_moe_gated(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .moe_gated_model_digest(m1)
  dg2 <- .moe_gated_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s\n", dg1, dg2))
  expect_identical(dg1, dg2)
  # Scores also bit-identical across the two fits.
  s1 <- score_moe_gated(m1, d$X[1:10, ])
  s2 <- score_moe_gated(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.9 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))

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
  invisible(score_moe_gated(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_moe_gated(d$X, d$y, hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "moe_gated_model")
})

test_that("§7.10 independent from-scratch mixture recompute bit-matches the scorer", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 101L)
  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_moe_gated(m, X_test)
  X_use <- .moe_gated_align(X_test, m$feature_universe)
  L <- length(m$weights)
  act <- if (identical(m$activation, "tanh")) tanh else function(x) pmax(x, 0)
  worst <- 0; worst_ln <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .moe_gated_rclr(X_use[i, ])
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
    # cross-check the encoder forward against .moe_gated_forward
    emb_ref <- as.numeric(.moe_gated_forward(matrix(z, nrow = 1L), m$weights,
                                             m$layernorm, m$activation))
    worst_ln <- max(worst_ln, max(abs(emb - emb_ref)))
    # --- from-scratch self-gated mixture (manual softmax + convex mix) ----------
    h_e <- as.numeric(m$W_e %*% emb) + m$b_e        # E expert logits
    g_l <- as.numeric(m$W_g %*% emb) + m$b_g        # E gate logits
    gw <- exp(g_l - max(g_l)); gw <- gw / sum(gw)   # softmax over experts
    sc <- sum(gw * h_e)                             # convex combination
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.10] independent mixture recompute maxdiff = %.3e ; LayerNorm-forward maxdiff = %.3e\n",
              worst, worst_ln))
  expect_lt(worst, 1e-12)
  expect_lt(worst_ln, 1e-12)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(seed = 111L)
  expect_error(fit_moe_gated(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(weight_decay = -1)),
               "weight_decay")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(hidden = c(0L, 4L))), "hidden")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(activation = "gelu")),
               "activation")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(n_experts = 1L)), "n_experts")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(n_experts = 2.5)), "n_experts")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_moe_gated(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_moe_gated(d$X, d$y,
                  hp = structure(list(1e-3, 1e-3),
                                 names = c("lr", "lr"))),
    "duplicate hp")
})

test_that("§7.11 GATING-SPECIFIC: gate is a valid softmax and the mixture is a convex combination", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(n = 120L, seed = 121L)
  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 120L, n_experts = 4L))
  X_use <- .moe_gated_align(d$X, m$feature_universe)
  Z <- .moe_gated_rclr_matrix(X_use)
  emb <- .moe_gated_forward(Z, m$weights, m$layernorm, m$activation)
  mix <- .moe_gated_mixture(emb, m$W_e, m$b_e, m$W_g, m$b_g)
  G <- mix$G; H <- mix$H; sc <- mix$score
  # Gate weights are a valid softmax: rows sum to 1, all >= 0.
  rs <- rowSums(G)
  cat(sprintf("[§7.11] gate rowSums in [%.6f, %.6f]; min gate weight = %.3e\n",
              min(rs), max(rs), min(G)))
  expect_lt(max(abs(rs - 1)), 1e-12)
  expect_true(all(G >= 0))
  expect_equal(ncol(G), m$n_experts)
  # The mixture score is a convex combination -> lies within [min_e h_e, max_e h_e]
  # per row (a softmax-weighted average of the expert logits).
  hmin <- apply(H, 1L, min); hmax <- apply(H, 1L, max)
  tol <- 1e-9
  expect_true(all(sc >= hmin - tol))
  expect_true(all(sc <= hmax + tol))
  cat(sprintf("[§7.11] convex-combination slack: min(sc-hmin) = %.3e ; min(hmax-sc) = %.3e\n",
              min(sc - hmin), min(hmax - sc)))
  # The scorer's score equals the mixture-helper score (one code path).
  expect_equal(score_moe_gated(m, d$X), sc, tolerance = 1e-12)
})

test_that("§7.12 EXPERT-DIVERSITY: experts are not all identical; a hand-built 2-expert toy reproduces the mixture", {
  skip_if_no_moe_gated()
  d <- .make_moe_data(n = 120L, seed = 131L)
  m <- fit_moe_gated(d$X, d$y, hp = list(epochs = 150L, n_experts = 4L))
  # Non-tautology: the exported expert weight rows are NOT all identical (the MoE
  # is not a degenerate single expert). Compare each expert row against the first.
  E <- m$n_experts
  max_pair_diff <- 0
  for (e in 2:E) {
    max_pair_diff <- max(max_pair_diff,
                         max(abs(m$W_e[e, ] - m$W_e[1, ])),
                         abs(m$b_e[e] - m$b_e[1]))
  }
  cat(sprintf("[§7.12] max expert-vs-expert1 weight/bias diff = %.4f (experts differ)\n",
              max_pair_diff))
  expect_gt(max_pair_diff, 1e-3)                  # experts genuinely differ

  # Hand-built 2-expert, 2-feature toy with KNOWN weights reproduces the pure-R
  # mixture score to < 1e-10 (independent of the trained model). z = (1, -2).
  z <- c(1, -2)
  W_e <- matrix(c(0.5, -1.0,    # expert 1: h1 = 0.5*z1 - 1.0*z2 + 0.1
                  2.0,  0.3),   # expert 2: h2 = 2.0*z1 + 0.3*z2 - 0.4
                nrow = 2L, byrow = TRUE)
  b_e <- c(0.1, -0.4)
  W_g <- matrix(c(1.0,  0.0,    # gate logit 1 = 1.0*z1
                  0.0,  1.0),   # gate logit 2 = 1.0*z2
                nrow = 2L, byrow = TRUE)
  b_g <- c(0.0, 0.0)
  h1 <- sum(W_e[1, ] * z) + b_e[1]               # 0.5 - (-2) + 0.1 = 2.6 -> 0.5*1 -1*(-2)+0.1
  h2 <- sum(W_e[2, ] * z) + b_e[2]
  gl <- c(sum(W_g[1, ] * z) + b_g[1], sum(W_g[2, ] * z) + b_g[2])   # (1, -2)
  gw <- exp(gl - max(gl)); gw <- gw / sum(gw)
  expected <- gw[1] * h1 + gw[2] * h2
  got <- .moe_gated_mixture(matrix(z, nrow = 1L), W_e, b_e, W_g, b_g)$score
  cat(sprintf("[§7.12] hand-built 2-expert toy: expected = %.10f ; got = %.10f\n",
              expected, got))
  expect_equal(got, expected, tolerance = 1e-10)
  # Gate softmax of the toy is a valid distribution.
  toy_G <- .moe_gated_mixture(matrix(z, nrow = 1L), W_e, b_e, W_g, b_g)$G
  expect_equal(sum(toy_G), 1, tolerance = 1e-12)
  expect_true(all(toy_G >= 0))
})
