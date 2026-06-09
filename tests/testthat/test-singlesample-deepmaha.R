
# lrt-deepmaha trains a small BN-free MLP embedding via reticulate-python torch at
# FIT (torch is confined to fit; scoring is pure R). Tests that fit are skipped
# when reticulate or the venv torch is unavailable, so the package suite stays
# green on hosts without the venv. Plain torch training needs NO brew-OpenSSL
# bridge and NO HF download (no _ssl / no checkpoint), so the skip-guard only
# probes reticulate + `import torch`.
skip_if_no_deepmaha <- function() {
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
# discriminative signal the MLP learns. Abundances are exp() so the input is
# strictly positive (compositional).
.make_deepmaha_data <- function(n = 80L, p = 30L, k = 8L, shift = 1.4,
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

.auc_mw_dm <- function(y, score) {
  case <- score[y == 1L]; ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}


test_that("fit returns a well-formed lrt_deepmaha_model", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data()
  m <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))
  expect_s3_class(m, "lrt_deepmaha_model")
  expect_identical(m$feature_universe, colnames(d$X))
  expect_equal(m$embedding_dim, 32L)             # last hidden = embedding dim
  expect_length(m$weights, 2L)                   # 2 linear layers up to embedding
  expect_equal(dim(m$weights[[1]]$W), c(64L, ncol(d$X)))
  expect_equal(dim(m$weights[[2]]$W), c(32L, 64L))
  expect_length(m$mu_case, 32L)
  expect_length(m$mu_ctrl, 32L)
  expect_equal(dim(m$Sigma_reg_inv), c(32L, 32L))
  expect_true(all(is.finite(m$Sigma_reg_inv)))
  s <- score_lrt_deepmaha(m, d$X)
  expect_type(s, "double")
  expect_length(s, nrow(d$X))
  expect_true(all(is.finite(s)))
})

test_that("§7.1 R forward matches the torch penultimate embedding (~1e-6)", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 21L)
  m <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 50L, device = "cpu", seed = 7L))
  # Rebuild the SAME torch net (same seed/arch), train identically, and compare the
  # penultimate (embedding) torch output to the pure-R forward on the training rows.
  torch <- reticulate::import("torch", delay_load = FALSE)
  nn <- torch$nn; np <- reticulate::import("numpy", delay_load = FALSE)
  Z_in <- .lrt_deepmaha_rclr_matrix(d$X)
  torch$manual_seed(7L)
  X_py <- torch$as_tensor(np$asarray(Z_in, dtype = "float32"))$float()
  y_py <- torch$as_tensor(np$asarray(as.integer(d$y), dtype = "int64"))$long()
  net <- nn$Sequential(nn$Linear(ncol(d$X), 64L), nn$ReLU(),
                       nn$Linear(64L, 32L), nn$ReLU(), nn$Linear(32L, 2L))
  opt <- torch$optim$Adam(net$parameters(), lr = 1e-3, weight_decay = 1e-4)
  lossf <- nn$CrossEntropyLoss()
  n <- nrow(Z_in); bs <- min(32L, n)
  net$train()
  for (ep in seq_len(50L)) {
    for (s0 in seq.int(0L, n - 1L, by = bs)) {
      idx <- torch$arange(as.integer(s0), as.integer(min(s0 + bs, n)),
                          dtype = torch$long)
      opt$zero_grad()
      loss <- lossf(net(X_py$index_select(0L, idx)), y_py$index_select(0L, idx))
      loss$backward(); opt$step()
    }
  }
  net$eval()
  # Penultimate embedding = forward up to (and incl.) the second Linear, no ReLU.
  # Build it from the trained net's first three modules (Linear, ReLU, Linear).
  mods <- reticulate::iterate(net$children())
  emb_net <- nn$Sequential(mods[[1]], mods[[2]], mods[[3]])
  emb_net$eval()
  torch_emb <- reticulate::py_to_r(
    emb_net(X_py)$detach()$to(torch$float64)$numpy())
  # Embed those SAME trained weights with the pure-R forward (drop the head).
  exp_w <- list()
  for (mm in reticulate::iterate(net$children())) {
    has_w <- tryCatch(!is.null(mm$weight), error = function(e) FALSE)
    if (!isTRUE(has_w)) next
    W <- reticulate::py_to_r(mm$weight$detach()$to(torch$float64)$numpy())
    b <- reticulate::py_to_r(mm$bias$detach()$to(torch$float64)$numpy())
    exp_w[[length(exp_w) + 1L]] <- list(W = W, b = as.numeric(b))
  }
  exp_w[[length(exp_w)]] <- NULL                     # drop the 2-logit head
  r_emb <- .lrt_deepmaha_forward(Z_in, exp_w, "relu")
  maxdiff <- max(abs(torch_emb - r_emb))
  cat(sprintf("\n[§7.1] R-vs-torch embedding maxdiff = %.3e\n", maxdiff))
  expect_lt(maxdiff, 1e-4)                            # float32 weights: ~1e-6
})

test_that("§7.2 single-row == batch (decisive single-sample, maxdiff 0)", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 31L)
  m <- fit_lrt_deepmaha(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  batch <- score_lrt_deepmaha(m, X_test)
  singles <- vapply(seq_len(nrow(X_test)), function(i) {
    score_lrt_deepmaha(m, X_test[i, , drop = FALSE])
  }, numeric(1L))
  md <- max(abs(batch - singles))
  cat(sprintf("[§7.2] single-row vs batch maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.3 row-permutation invariance (maxdiff 0)", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 41L)
  m <- fit_lrt_deepmaha(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  base <- score_lrt_deepmaha(m, X_test)
  set.seed(7); perm <- sample(nrow(X_test))
  permuted <- score_lrt_deepmaha(m, X_test[perm, , drop = FALSE])
  md <- max(abs(permuted - base[perm]))
  cat(sprintf("[§7.3] row-permutation maxdiff = %.3e\n", md))
  expect_identical(md, 0)
})

test_that("§7.4 per-sample scale-invariance (maxdiff ~0)", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 44L)
  m <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 60L))
  worst <- 0
  for (i in c(1L, 3L, 41L, 70L)) {
    row <- d$X[i, , drop = FALSE]
    for (s in c(7, 1e6, 1e-3)) {
      worst <- max(worst, abs(score_lrt_deepmaha(m, row * s) -
                                score_lrt_deepmaha(m, row)))
    }
  }
  batch <- d$X[1:8, , drop = FALSE]
  set.seed(5); scaled <- batch * stats::runif(nrow(batch), 0.01, 100)
  worst <- max(worst, max(abs(score_lrt_deepmaha(m, scaled) -
                                score_lrt_deepmaha(m, batch))))
  cat(sprintf("[§7.4] scale-invariance maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-7)
})

test_that("§7.5 held-out AUC on a strongly planted shift is high (> 0.7)", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(n = 160L, seed = 51L)
  ctrl <- which(d$y == 0L); case <- which(d$y == 1L)
  tr <- c(utils::head(ctrl, 60L), utils::head(case, 60L))
  te <- c(utils::tail(ctrl, 20L), utils::tail(case, 20L))
  m <- fit_lrt_deepmaha(d$X[tr, ], d$y[tr], hp = list(epochs = 200L))
  s <- score_lrt_deepmaha(m, d$X[te, ])
  auc <- .auc_mw_dm(d$y[te], s)
  cat(sprintf("[§7.5] held-out AUC = %.4f\n", auc))
  expect_true(all(is.finite(s)))
  expect_gt(mean(s[d$y[te] == 1]), mean(s[d$y[te] == 0]))
  expect_gt(auc, 0.7)
})

test_that("§7.6 canonical row-equivariance gate PASSES + dispatch matches", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 61L)
  m <- fit_lrt_deepmaha(d$X[1:60, ], d$y[1:60], hp = list(epochs = 60L))
  X_test <- d$X[61:80, ]
  score_fun <- function(model, X, meta) score_lrt_deepmaha(model, X, meta)
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_fun, m, X_test, model_digest = .lrt_deepmaha_model_digest))
  expect_true(
    singlesample_is_row_equivariant(
      score_fun, m, X_test, model_digest = .lrt_deepmaha_model_digest))

  # Canonical dispatch: post-integration the manifest carries the lrt-deepmaha row
  # and score_lrt_deepmaha is in the namespace; pre-integration (staged self-test)
  # clone a one-row roster from a committed within-discriminator row and register a
  # thin adapter so the same dispatch wiring is exercised. Either way the dispatched
  # score must equal the direct call.
  direct <- score_lrt_deepmaha(m, X_test)
  roster <- singlesample_method_roster()
  if (!"lrt-deepmaha" %in% roster$method_id) {
    tmpl <- roster[roster$method_id == "lrt-copula", , drop = FALSE]
    tmpl$method_id <- "lrt-deepmaha"
    tmpl$fit_fn <- "fit_lrt_deepmaha"
    tmpl$score_fn <- "score_lrt_deepmaha"
    roster <- rbind(roster, tmpl)
  }
  if (!exists("score_lrt_deepmaha", envir = asNamespace("OmicSelector"),
              mode = "function", inherits = FALSE)) {
    singlesample_register_score_adapter(
      "lrt-deepmaha",
      function(model, X, meta) score_lrt_deepmaha(model, X, meta))
  }
  via <- singlesample_score_call("lrt-deepmaha", m, X_test, roster = roster)
  expect_equal(via, direct, tolerance = 1e-12)
})

test_that("§7.7 degenerate -> 0; FLAT positive composition is NOT floored", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 71L)
  m <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 60L, min_features = 3L))

  # Fewer than min_features present -> neutral 0 for every row.
  X_two <- d$X[1:6, colnames(d$X)[1:2], drop = FALSE]
  expect_equal(score_lrt_deepmaha(m, X_two), rep(0, nrow(X_two)))
  # No shared features at all -> neutral 0.
  X_none <- matrix(1, nrow = 4L, ncol = 3L,
                   dimnames = list(paste0("N", 1:4), c("zzz1", "zzz2", "zzz3")))
  expect_equal(score_lrt_deepmaha(m, X_none), rep(0, nrow(X_none)))
  # All-zero specimen (empty positive support) -> 0 for that row, finite elsewhere.
  X_zero <- d$X[1:3, , drop = FALSE]; X_zero[1, ] <- 0
  s_zero <- score_lrt_deepmaha(m, X_zero)
  expect_equal(s_zero[1], 0)
  expect_true(all(is.finite(s_zero)))

  # rCLR-origin regression: a FLAT all-equal-positive composition has all(z == 0)
  # but full positive support -> it MUST be embedded and scored, NOT floored.
  p <- length(m$feature_universe)
  flat <- matrix(7.0, nrow = 1L, ncol = p,
                 dimnames = list("flat", m$feature_universe))
  z_flat <- .lrt_deepmaha_rclr(flat[1, ])
  expect_true(all(flat > 0))
  expect_true(all(z_flat == 0))                  # rCLR origin
  emb <- .lrt_deepmaha_forward(matrix(z_flat, nrow = 1L), m$weights, m$activation)
  manual <- .lrt_deepmaha_score_one(as.numeric(emb), m$mu_case, m$mu_ctrl,
                                    m$Sigma_reg_inv)
  s_flat <- score_lrt_deepmaha(m, flat)
  cat(sprintf("[§7.7] flat-composition score = %.6f (manual %.6f)\n",
              s_flat, manual))
  expect_equal(s_flat, manual, tolerance = 1e-12)
  # The flat embedding is the bias-only forward; it is generically nonzero, so the
  # flat score is genuinely computed (not the floored 0). Guard against the rare
  # exact-0 coincidence by checking equality to the manual recompute above only.

  # Contrast: a genuinely empty positive support IS floored to 0.
  empty <- matrix(0, nrow = 1L, ncol = p,
                  dimnames = list("empty", m$feature_universe))
  expect_equal(score_lrt_deepmaha(m, empty), 0)
})

test_that("§7.8 determinism: two seed=42 CPU fits give identical digests", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 81L)
  m1 <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  m2 <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 60L, device = "cpu", seed = 42L))
  dg1 <- .lrt_deepmaha_model_digest(m1)
  dg2 <- .lrt_deepmaha_model_digest(m2)
  cat(sprintf("[§7.8] digest1 = %s ; digest2 = %s\n", dg1, dg2))
  expect_identical(dg1, dg2)
  # Scores also bit-identical across the two fits.
  s1 <- score_lrt_deepmaha(m1, d$X[1:10, ])
  s2 <- score_lrt_deepmaha(m2, d$X[1:10, ])
  expect_identical(s1, s2)
})

test_that("§7.9 RNG-safety: fit leaves global R + torch RNG byte-unchanged; score consumes none", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 91L)
  torch <- reticulate::import("torch", delay_load = FALSE)

  set.seed(999)
  before_r <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  before_t <- torch$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch$cuda$get_rng_state_all() else NULL

  m <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 60L, device = "cpu"))

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
  invisible(score_lrt_deepmaha(m, d$X[1:20, ]))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   score_seed)

  # Fit restores RNG even when no global seed pre-exists.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  m2 <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 30L, device = "cpu"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(m2, "lrt_deepmaha_model")
})

test_that("§7.10 independent Mahalanobis recompute bit-matches the scorer", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 101L)
  m <- fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 60L))
  X_test <- d$X[1:12, ]
  direct <- score_lrt_deepmaha(m, X_test)
  # From-scratch recompute using stats::mahalanobis on the (covariance) matrix
  # S = solve(Sigma_reg_inv): -0.5 maha(z, mu_case, S) + 0.5 maha(z, mu_ctrl, S).
  S <- solve(m$Sigma_reg_inv)
  X_use <- .lrt_deepmaha_align(X_test, m$feature_universe)
  worst <- 0
  for (i in seq_len(nrow(X_use))) {
    z <- .lrt_deepmaha_rclr(X_use[i, ])
    emb <- as.numeric(.lrt_deepmaha_forward(matrix(z, nrow = 1L), m$weights,
                                            m$activation))
    sc <- -0.5 * stats::mahalanobis(matrix(emb, 1L), m$mu_case, S) +
      0.5 * stats::mahalanobis(matrix(emb, 1L), m$mu_ctrl, S)
    worst <- max(worst, abs(sc - direct[i]))
  }
  cat(sprintf("[§7.10] independent mahalanobis recompute maxdiff = %.3e\n", worst))
  expect_lt(worst, 1e-9)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_deepmaha()
  d <- .make_deepmaha_data(seed = 111L)
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(bogus = 1)), "unknown hp")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(epochs = 0)), "epochs")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(lr = 0)), "hp\\$lr")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(weight_decay = -1)),
               "weight_decay")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(hidden = c(0L, 4L))), "hidden")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(activation = "gelu")),
               "activation")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(batch_size = 1.5)),
               "batch_size")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(ridge = 0)), "ridge")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(min_features = 0)),
               "min_features")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(device = "tpu")), "device")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(seed = 1.5)), "seed")
  expect_error(fit_lrt_deepmaha(d$X, d$y, hp = list(3L)), "must be named")
  expect_error(
    fit_lrt_deepmaha(d$X, d$y,
                     hp = structure(list(1e-3, 1e-3),
                                    names = c("lr", "lr"))),
    "duplicate hp")
})
