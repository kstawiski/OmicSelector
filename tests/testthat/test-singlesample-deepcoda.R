
# coda-deepcoda trains a zero-sum log-contrast bottleneck + a self-explaining
# theta-net + scalar bias end-to-end via reticulate-python torch (BCEWithLogits) at
# FIT, then freezes (torch is confined to fit; scoring is pure R). Tests that fit are
# skipped when reticulate or the venv torch is unavailable, so the package suite
# stays green on hosts without the venv. Plain torch training needs NO brew-OpenSSL
# bridge and NO HF download, so the skip-guard only probes reticulate + `import
# torch`.
skip_if_no_deepcoda <- function() {
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
# the bottleneck + theta-net learn. Abundances are exp() so the input is strictly
# positive (compositional).
.make_deepcoda_data <- function(n = 80L, p = 30L, k = 8L, shift = 1.4,
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

.auc_mw_dc <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# Rebuild the DeepCoDA forward in torch from the EXPORTED R-side weights and return
# the logit on the training rows. The bottleneck uses the exported (already
# zero-sum) W via a plain matmul; the theta-net is a rebuilt Sequential; the logit
# is sum_k theta_k * b_k + bias. Used by the R-vs-torch parity test (§7.1).
.torch_deepcoda_logit <- function(m, Z_in) {
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  act_mod <- if (identical(m$activation, "tanh")) nn$Tanh else nn$ReLU
  X_py <- torch$as_tensor(np$asarray(Z_in, dtype = "float32"))$float()
  W_py <- torch$as_tensor(np$asarray(m$W_bottleneck, dtype = "float32"))$float()
  Bmat <- torch$matmul(X_py, W_py$t())                 # n x B (no activation)
  # rebuild the theta-net Sequential from the exported per-layer weights
  layers <- list()
  L <- length(m$theta_weights)
  for (l in seq_len(L)) {
    W <- m$theta_weights[[l]]$W; b <- m$theta_weights[[l]]$b
    lin <- nn$Linear(ncol(W), nrow(W))
    lin$weight$data <- torch$as_tensor(np$asarray(W, dtype = "float32"))$float()
    lin$bias$data   <- torch$as_tensor(np$asarray(b, dtype = "float32"))$float()
    layers[[length(layers) + 1L]] <- lin
    if (l < L) layers[[length(layers) + 1L]] <- act_mod()
  }
  net <- do.call(nn$Sequential, layers); net$eval()
  Theta <- net(Bmat)                                   # n x B
  logit <- Theta$mul(Bmat)$sum(1L)$add(m$bias)         # n
  reticulate::py_to_r(logit$detach()$to(torch$float64)$numpy())
}


test_that("fit returns a well-formed coda_deepcoda_model", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data()
  m <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "coda_deepcoda_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$bottleneck_dim, 8L)                   # default B
  expect_equal(dim(m$W_bottleneck), c(8L, ncol(d$X)))  # B x p
  # theta-net: Linear(B->16) -> Linear(16->B); 2 linear layers exported
  expect_length(m$theta_weights, 2L)
  expect_equal(dim(m$theta_weights[[1]]$W), c(16L, 8L))
  expect_equal(dim(m$theta_weights[[2]]$W), c(8L, 16L))
  expect_length(m$bias, 1L)
  s <- score_coda_deepcoda(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("§7.1 R forward matches the torch DeepCoDA logit (~1e-6)", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 21L)
  m <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 50L, device = "cpu", seed = 7L))
  Z_in <- .coda_deepcoda_rclr_matrix(d$X)
  # pure-R logit (the scorer's single code path) vs the torch model rebuilt from the
  # SAME exported weights. The float64 R forward must match the float32 torch forward
  # on the exported (float32-valued) weights to ~1e-6.
  r_logit <- .coda_deepcoda_logit(Z_in, m$W_bottleneck, m$theta_weights,
                                  m$bias, m$activation)
  t_logit <- .torch_deepcoda_logit(m, Z_in)
  maxdiff <- max(abs(r_logit - t_logit))
  cat(sprintf("\n[§7.1] R-vs-torch DeepCoDA logit maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                              # float32 weights: ~1e-6
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 31L)
  m <- fit_coda_deepcoda(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  batch <- score_coda_deepcoda(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_coda_deepcoda(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 41L)
  m <- fit_coda_deepcoda(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  base <- score_coda_deepcoda(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_coda_deepcoda(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (zero-sum contrast -> exact, maxdiff ~0)", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 44L)
  m <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 41L, 70L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_coda_deepcoda(m, row * s) -
                                score_coda_deepcoda(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_coda_deepcoda(m, scaled) -
                                score_coda_deepcoda(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-7)
})

test_that("§7.5 zero-sum constraint: exported bottleneck rows sum to zero", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 47L)
  m <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L))
  rs <- rowSums(m$W_bottleneck)
  worst <- max(abs(rs))
  cat(sprintf("[§7.5] max|rowSums(W_bottleneck)| = %.3e\n", worst))
  expect_lt(worst, 1e-6)
  # different bottleneck_dim also stays zero-sum
  m2 <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 30L, bottleneck_dim = 4L))
  expect_equal(nrow(m2$W_bottleneck), 4L)
  expect_lt(max(abs(rowSums(m2$W_bottleneck))), 1e-6)
})

test_that("§7.6 held-out AUC on a strongly planted shift is high (> 0.7)", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(n = 160L, seed = 51L)
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  tr <- c(utils::head(ctrl, 60L), utils::head(case, 60L))
  te <- c(utils::tail(ctrl, 20L), utils::tail(case, 20L))
  m <- fit_coda_deepcoda(d$X[tr, ], d$y[tr], hp = list(epochs = 300L))
  s <- score_coda_deepcoda(m, d$X[te, ])
  auc <- .auc_mw_dc(d$y[te], s)
  cat(sprintf("[§7.6] held-out AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.7 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 61L)
  m <- fit_coda_deepcoda(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  score_fun <- function(model, X, meta) score_coda_deepcoda(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .coda_deepcoda_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .coda_deepcoda_model_digest))

  # Canonical dispatch: post-integration the manifest carries the coda-deepcoda row
  # and score_coda_deepcoda is in the namespace; pre-integration (staged self-test)
  # clone a one-row roster from a committed within-discriminator row and register a
  # thin adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_coda_deepcoda(m, X_test)
  roster <- singlesample_method_roster()
  if (!"coda-deepcoda" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "coda-deepcoda"
    tmpl$fit_fn <- "fit_coda_deepcoda"
    tmpl$score_fn <- "score_coda_deepcoda"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_coda_deepcoda", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "coda-deepcoda",
      function(model, X, meta) score_coda_deepcoda(model, X, meta))
  }
  via <- singlesample_score_call("coda-deepcoda", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.8 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 71L)
  m <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_coda_deepcoda(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_coda_deepcoda(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_coda_deepcoda(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be projected + scored, NOT floored. Its
  # bottlenecks are all 0, so the logit reduces to the learned bias (a genuine
  # computed value, distinct from a floored 0 in general).
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .coda_deepcoda_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  manual <- .coda_deepcoda_logit(matrix(z_flat, nrow = 1L), m$W_bottleneck,
                                 m$theta_weights, m$bias, m$activation)
  s_flat <- score_coda_deepcoda(m, flat)
  cat(sprintf("[§7.8] flat-composition score = %.6f (manual %.6f ; bias %.6f)\n",
              s_flat, manual, m$bias))
  expect_equal(s_flat, as.numeric(manual), tolerance = 1e-12)
  # For all-zero bottlenecks the logit equals the learned bias exactly.
  expect_equal(s_flat, m$bias, tolerance = 1e-12)

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_coda_deepcoda(m, empty), 0)
})

test_that("§7.9 determinism: two seed=42 CPU fits give identical digests", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 81L)
  m1 <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  # Perturb the GLOBAL torch RNG between the two fits. A correct seed-before-build
  # fit calls manual_seed(seed) BEFORE constructing the (stochastic) module, so this
  # perturbation is overwritten and the two digests stay identical. But if a future
  # edit ever moved module construction BEFORE the seed, the two fits would diverge
  # here and this test would catch the cross-session non-reproducibility (the
  # ai-scarf bug class that a back-to-back same-seed digest check alone would mask).
  torch_pert <- reticulate::import("torch", delay_load = FALSE)
  torch_pert$manual_seed(99L)
  invisible(torch_pert$randn(list(64L, 64L)))
  m2 <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .coda_deepcoda_model_digest(m1)
  dg2 <- .coda_deepcoda_model_digest(m2)
  cat(sprintf("[§7.9] digest1 = %s ; digest2 = %s\n", dg1, dg2))
  expect_identical(dg1, dg2)
  # Scores also bit-identical across the two fits.
  s1 <- score_coda_deepcoda(m1, d$X[1:10, ])
  s2 <- score_coda_deepcoda(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.10 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))

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
  cat(sprintf("[§7.10] R RNG unchanged = %s ; torch CPU RNG unchanged = %s\n",
              identical(before_r, after_r),
              isTRUE(reticulate::py_to_r(torch$equal(before_t, after_t)))))

  # Scoring consumes no RNG: .Random.seed unchanged across scoring.
  score_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  invisible(score_coda_deepcoda(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "coda_deepcoda_model")
})

test_that("§7.11 independent non-circular logit recompute bit-matches the scorer", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 101L)
  m <- fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_coda_deepcoda(m, X_test)
  # From-scratch recompute that does NOT call .coda_deepcoda_logit /
  # .coda_deepcoda_theta_forward: bottleneck via explicit per-row inner products,
  # theta-net via an independent layer loop, logit via sum(theta * b) + bias.
  act <- if (identical(m$activation, "tanh")) tanh else function(x) pmax(x, 0)
  X_use <- .coda_deepcoda_align(X_test, m$feature_universe)
  worst <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .coda_deepcoda_rclr(X_use[i, ])
    # bottleneck b_k = sum_j W_kj * z_j (explicit), zero-sum W already exported
    b <- as.numeric(m$W_bottleneck %*% z)            # length B
    # theta-net forward, layer by layer (independent of the scorer helper)
    h <- b
    Lw <- length(m$theta_weights)
    for (l in seq_len(Lw)) {
      W <- m$theta_weights[[l]]$W; bb <- m$theta_weights[[l]]$b
      h <- as.numeric(W %*% h) + bb
      if (l < Lw) h <- act(h)
    }
    theta <- h                                       # length B coefficients
    sc <- sum(theta * b) + m$bias
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.11] independent logit recompute maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-9)
})

test_that("§7.12 tanh activation: R forward matches torch + single-row == batch", {
  skip_if_no_deepcoda()
  # Additive coverage over a byte-unchanged scorer for the `tanh` theta-net branch
  # (.coda_deepcoda_theta_forward selects tanh when hp$activation == "tanh"),
  # mirroring §7.1 (relu).
  d <- .make_deepcoda_data(seed = 23L)
  m <- fit_coda_deepcoda(d$X, d$y,
                         hp = list(epochs = 50L, device = "cpu", seed = 7L,
                                   activation = "tanh"))
  expect_identical(m$activation, "tanh")
  Z_in <- .coda_deepcoda_rclr_matrix(d$X)
  r_logit <- .coda_deepcoda_logit(Z_in, m$W_bottleneck, m$theta_weights,
                                  m$bias, m$activation)
  t_logit <- .torch_deepcoda_logit(m, Z_in)
  maxdiff <- max(abs(r_logit - t_logit))
  cat(sprintf("[§7.12] tanh R-vs-torch logit maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                              # float32 weights: ~1e-7
  # the decisive single-sample property holds under tanh too
  batch <- score_coda_deepcoda(m, d$X)
  singles <- vapply(seq_len(nrow(d$X)), function(i) {
    score_coda_deepcoda(m, d$X[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.12] tanh single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_deepcoda()
  d <- .make_deepcoda_data(seed = 111L)
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(bottleneck_dim = 0)),
               "bottleneck_dim")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(bottleneck_dim = 1.5)),
               "bottleneck_dim")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(theta_hidden = c(0L, 4L))),
               "theta_hidden")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(activation = "gelu")),
               "activation")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(weight_decay = -1)),
               "weight_decay")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_coda_deepcoda(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_coda_deepcoda(d$X, d$y,
                      hp = structure(list(1e-3, 1e-3),
                                     names = c("lr", "lr"))),
    "duplicate hp")
})
