# ai-scarf trains a small BN-free MLP encoder via reticulate-python torch by
# SCARF SELF-SUPERVISED CONTRASTIVE training (Bahri 2022; InfoNCE / NT-Xent, NO
# labels) at FIT, FREEZES it, fits a frozen linear-probe head, and scores in pure
# R. Tests that fit are skipped when reticulate or the venv torch is unavailable,
# so the package suite stays green on hosts without the venv. SCARF training needs
# NO brew-OpenSSL bridge and NO HF download, so the skip-guard only probes
# reticulate + `import torch`.
skip_if_no_scarf <- function() {
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
# discriminative signal the encoder + head learn. Abundances are exp() so the
# input is strictly positive (compositional).
.make_scarf_data <- function(n = 80L, p = 30L, k = 8L, shift = 1.4,
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

.auc_mw_sc <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}


test_that("fit returns a well-formed ai_scarf_model", {
  skip_if_no_scarf()
  d <- .make_scarf_data()
  m <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "ai_scarf_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$embedding_dim, 32L)             # last hidden = embedding dim
  expect_length(m$weights, 2L)                   # 2 linear layers up to embedding
  expect_equal(dim(m$weights[[1]]$W), c(64L, ncol(d$X)))
  expect_equal(dim(m$weights[[2]]$W), c(32L, 64L))
  expect_length(m$head_w, 32L)                   # linear head weight, length d
  expect_length(m$head_b, 1L)                    # linear head bias, scalar
  expect_null(m$marginal)                        # corruption marginal NOT stored
  s <- score_ai_scarf(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("§7.1 R forward matches the torch encoder embedding (~1e-6)", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 21L)
  m <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 50L, device = "cpu", seed = 7L))
  # Embed the training rows with the EXPORTED R-side encoder weights (pure-R
  # forward), and with a torch Sequential rebuilt from those SAME exported weights;
  # compare. This is the fit/score-consistency check: the float64 R forward must
  # match the float32 torch forward on the exported (float32-valued) weights.
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  Z_in <- .ai_scarf_rclr_matrix(d$X)
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
  r_emb <- .ai_scarf_forward(Z_in, m$weights, "relu")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("\n[§7.1] R-vs-torch encoder embedding maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                            # float32 weights: ~1e-6
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 31L)
  m <- fit_ai_scarf(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  batch <- score_ai_scarf(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_ai_scarf(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 41L)
  m <- fit_ai_scarf(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  base <- score_ai_scarf(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_ai_scarf(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (maxdiff ~0)", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 44L)
  m <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 41L, 70L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_ai_scarf(m, row * s) -
                                score_ai_scarf(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_ai_scarf(m, scaled) -
                                score_ai_scarf(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-7)
})

test_that("§7.5 held-out AUC on a strongly planted shift is high (> 0.7)", {
  skip_if_no_scarf()
  d <- .make_scarf_data(n = 160L, seed = 51L)
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  tr <- c(utils::head(ctrl, 60L), utils::head(case, 60L))
  te <- c(utils::tail(ctrl, 20L), utils::tail(case, 20L))
  m <- fit_ai_scarf(d$X[tr, ], d$y[tr], hp = list(epochs = 200L))
  s <- score_ai_scarf(m, d$X[te, ])
  auc <- .auc_mw_sc(d$y[te], s)
  cat(sprintf("[§7.5] held-out AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.6 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 61L)
  m <- fit_ai_scarf(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  score_fun <- function(model, X, meta) score_ai_scarf(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .ai_scarf_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .ai_scarf_model_digest))

  # Canonical dispatch: post-integration the manifest carries the ai-scarf row and
  # score_ai_scarf is in the namespace; pre-integration (staged self-test) clone a
  # one-row roster from a committed within-discriminator row and register a thin
  # adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_ai_scarf(m, X_test)
  roster <- singlesample_method_roster()
  if (!"ai-scarf" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "ai-scarf"
    tmpl$fit_fn <- "fit_ai_scarf"
    tmpl$score_fn <- "score_ai_scarf"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_ai_scarf", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "ai-scarf",
      function(model, X, meta) score_ai_scarf(model, X, meta))
  }
  via <- singlesample_score_call("ai-scarf", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.7 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 71L)
  m <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_ai_scarf(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_ai_scarf(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_ai_scarf(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be embedded and scored, NOT floored.
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .ai_scarf_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  emb <- .ai_scarf_forward(matrix(z_flat, nrow = 1L), m$weights, m$activation)
  manual <- .ai_scarf_score_one(as.numeric(emb), m$head_w, m$head_b)
  s_flat <- score_ai_scarf(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (manual %.6f)\n",
              s_flat, manual))
  expect_equal(s_flat, manual, tolerance = 1e-12)
  # The flat embedding is the bias-only forward; the head logit on it is
  # generically nonzero, so the flat score is genuinely computed (not floored 0).

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_ai_scarf(m, empty), 0)
})

test_that("§7.8 determinism: two seed=42 CPU fits give identical digests", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 81L)
  m1 <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  # Perturb the GLOBAL torch RNG to a different state before the second fit. A
  # correctly seeded fit re-seeds every stochastic step (encoder init INCLUDED)
  # from hp$seed, so m2 must still be bit-identical to m1 despite the perturbed
  # entry state. This makes the in-process test catch the seed-after-module-
  # construction class of bug (encoder Linears initialised from the unseeded
  # process-entry RNG), which a plain back-to-back two-fit test masks via the
  # fit's on.exit RNG restore -- and which would otherwise surface only as
  # cross-session non-reproducibility.
  torch_pert <- reticulate::import("torch", delay_load = FALSE)
  torch_pert$manual_seed(20260609L)
  invisible(torch_pert$randn(list(64L, 64L)))
  m2 <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .ai_scarf_model_digest(m1)
  dg2 <- .ai_scarf_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s\n", dg1, dg2))
  expect_identical(dg1, dg2)
  # Both the SSL-stage encoder weights and the head (w, b) are bit-identical.
  expect_identical(m1$weights, m2$weights)
  expect_identical(m1$head_w, m2$head_w)
  expect_identical(m1$head_b, m2$head_b)
  # Scores also bit-identical across the two fits.
  s1 <- score_ai_scarf(m1, d$X[1:10, ])
  s2 <- score_ai_scarf(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.9 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))

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
  invisible(score_ai_scarf(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "ai_scarf_model")
})

test_that("§7.10 corruption marginal is fit-only (not in score path) + independent head recompute", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 101L)
  m <- fit_ai_scarf(d$X, d$y, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_ai_scarf(m, X_test)

  # (a) PROVE the frozen training marginal is NOT consulted at score: the model
  # carries no marginal/training field, and scoring with ONLY the exported model
  # (rebuilt from a stripped copy holding solely the frozen score-determining
  # fields) reproduces the in-fit score EXACTLY. This mimics a fresh process that
  # only ever saw the exported model -- it cannot access the corruption marginal.
  expect_null(m$marginal)
  expect_null(m$X_train); expect_null(m$Z_train); expect_null(m$training_columns)
  model_keys <- c("feature_universe", "weights", "activation", "head_w",
                  "head_b", "embedding_dim", "device", "seed", "hp")
  expect_true(all(names(m) %in% model_keys))     # no extra (marginal-like) fields
  m_export <- m[model_keys]
  class(m_export) <- "ai_scarf_model"
  s_export <- score_ai_scarf(m_export, X_test)
  cat(sprintf("[§7.10] export-only vs in-fit score maxdiff = %.3e\n",
              max(abs(s_export - direct))))
  expect_identical(s_export, direct)

  # (b) Independent (non-circular) linear-head recompute: sum(w * z) + b computed
  # from scratch (does NOT call .ai_scarf_score_one) must bit-match the scorer.
  X_use <- .ai_scarf_align(X_test, m$feature_universe)
  worst <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .ai_scarf_rclr(X_use[i, ])
    emb <- as.numeric(.ai_scarf_forward(matrix(z, nrow = 1L), m$weights,
                                        m$activation))
    sc <- as.numeric(crossprod(m$head_w, emb)) + m$head_b   # w . z + b
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.10] independent linear-head recompute maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-9)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_scarf()
  d <- .make_scarf_data(seed = 111L)
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(weight_decay = -1)),
               "weight_decay")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(hidden = c(0L, 4L))), "hidden")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(activation = "gelu")),
               "activation")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(corrupt_p = 0)), "corrupt_p")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(corrupt_p = 1)), "corrupt_p")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(temperature = 0)), "temperature")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(proj_dim = 0)), "proj_dim")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(head_epochs = 1.5)),
               "head_epochs")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(head_weight_decay = -1)),
               "head_weight_decay")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_ai_scarf(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_ai_scarf(d$X, d$y,
                 hp = structure(list(1e-3, 1e-3),
                                names = c("lr", "lr"))),
    "duplicate hp")
})

test_that("§7.11 InfoNCE loss includes the positive in the partition (== cross-entropy)", {
  skip_if_no_scarf()
  # Directly exercises the SCARF contrastive loss (the method's defining component,
  # otherwise untested). The canonical InfoNCE / NT-Xent loss (Bahri 2022 Alg.1)
  # sums the partition function over ALL k INCLUDING the positive k=i, which is
  # exactly softmax cross-entropy with diagonal labels. Guards against regressing
  # to a positive-EXCLUDED (masked-diagonal) denominator.
  torch <- reticulate::import("torch", delay_load = FALSE)
  F <- torch$nn$functional
  torch$manual_seed(7L)
  n <- 6L
  zc <- F$normalize(torch$randn(list(n, 8L)), dim = 1L)
  zk <- F$normalize(torch$randn(list(n, 8L)), dim = 1L)
  logits <- torch$matmul(zc, zk$t())$div(1.0)
  ours <- as.numeric(reticulate::py_to_r(.ai_scarf_infonce_loss(logits)$item()))
  ce   <- as.numeric(reticulate::py_to_r(
    F$cross_entropy(logits, torch$arange(n)$long())$item()))
  cat(sprintf("[§7.11] InfoNCE loss = %.6f ; cross_entropy = %.6f ; diff = %.2e\n",
              ours, ce, abs(ours - ce)))
  expect_lt(abs(ours - ce), 1e-5)                    # positive IS in the partition
  # the positive-EXCLUDED variant (masked diagonal) must measurably differ -- this
  # is the bug the gate caught; lock it out.
  eye <- torch$eye(n)
  excl <- as.numeric(reticulate::py_to_r(torch$mean(torch$sub(
    torch$logsumexp(torch$add(logits, torch$mul(eye, -1e9)), dim = 1L),
    logits$diagonal()))$item()))
  cat(sprintf("[§7.11] positive-excluded variant = %.6f (must differ)\n", excl))
  expect_gt(abs(ours - excl), 1e-3)
})
