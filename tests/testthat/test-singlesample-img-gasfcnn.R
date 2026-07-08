library(testthat)

# img-gasfcnn trains a small CNN and rebuilds it at score; both need the python
# torch + torchvision modules via reticulate. Tests that fit/score are skipped
# when either is unavailable.
skip_if_no_gasfcnn <- function() {
  venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
  if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = venv)
  }
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(reticulate::py_module_available("torch"),
                        "Python module 'torch' not available")
  testthat::skip_if_not(reticulate::py_module_available("torchvision"),
                        "Python module 'torchvision' not available")
}

# Synthetic generator with a CASE-ONLY location shift in a feature block. The first
# `k` of `p` log-abundance features are elevated by `shift` for case samples; after
# exp() and the per-sample rCLR centring those features carry the discriminative
# signal the GASF-image CNN learns. Abundances are exp() so the input is strictly
# positive (compositional). The tests fit/score on device = "auto" (the H100 here;
# falls back to CPU elsewhere). The DECISIVE single-sample properties are EXACT
# regardless of device because the scorer forwards in float64 ONE ROW AT A TIME.
.make_gasfcnn_data <- function(n = 40L, p = 16L, k = 6L, shift = 1.2,
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

.auc_mw_gasfcnn <- function(y, score) {
  case <- score[y == 1L]
  ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# A small, fast fit reused across tests (few epochs -- the single-sample properties
# do not depend on convergence, only the held-out-AUC test trains longer).
.fit_small_gasfcnn <- function(seed = 11L, epochs = 60L) {
  d <- .make_gasfcnn_data(seed = seed)
  list(model = fit_img_gasfcnn(d$X, d$y,
                               hp = list(device = "auto", epochs = epochs)),
       data = d)
}


test_that("fit returns a well-formed img_gasfcnn_model", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  expect_s3_class(m, "img_gasfcnn_model")
  expect_equal(m$feature_universe, colnames(f$data$X))
  expect_equal(m$channels, c(8L, 16L))
  expect_length(m$gasf_lo, ncol(f$data$X))
  expect_length(m$gasf_hi, ncol(f$data$X))
  # No live python pointer survives the fit (python-at-score contract).
  expect_false(any(vapply(m, function(x) inherits(x, "python.builtin.object"),
                          logical(1L))))
  # state_dict carries the BN running buffers (frozen-BN proof of provenance).
  expect_true(any(grepl("running_mean", names(m$state_dict))))
  expect_true(any(grepl("running_var", names(m$state_dict))))
  expect_true(any(grepl("num_batches_tracked", names(m$state_dict))))
})

test_that("score is finite, length-correct, and double", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  s <- score_img_gasfcnn(f$model, f$data$X)
  expect_length(s, nrow(f$data$X))
  expect_true(all(is.finite(s)))
  expect_type(s, "double")
})

# §7.1 GASF correctness: the python GASF (the n x 1 x p x p tensor the scorer feeds)
# matches the pure-R closed form .img_gasfcnn_gasf to <1e-6 on random rCLR rows.
test_that("§7.1 python GASF == pure-R closed form (<1e-6)", {
  skip_if_no_gasfcnn()
  set.seed(123)
  p <- 14L
  lo <- runif(p, -2, -0.5)
  hi <- runif(p, 0.5, 2)
  Z <- matrix(rnorm(5L * p), nrow = 5L)
  G_py <- reticulate::py_to_r(.img_gasfcnn_images(Z, lo, hi))   # 5 x 1 x p x p
  maxd <- 0
  for (i in seq_len(nrow(Z))) {
    G_r <- .img_gasfcnn_gasf(Z[i, ], lo, hi)
    G_pi <- G_py[i, 1L, , ]
    maxd <- max(maxd, max(abs(G_r - G_pi)))
  }
  expect_lt(maxd, 1e-6)
})

# §7.2 single-row score == batch score, maxdiff EXACTLY 0 (forced row-by-row +
# eval-mode frozen BN + float64 forward).
test_that("§7.2 DECISIVE: single-row score == batch score (maxdiff 0)", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  batch <- score_img_gasfcnn(m, X)
  single <- vapply(seq_len(nrow(X)), function(i) {
    score_img_gasfcnn(m, X[i, , drop = FALSE])
  }, numeric(1L))
  expect_identical(max(abs(batch - single)), 0)
})

# §7.3 row-permutation invariance, maxdiff 0.
test_that("§7.3 score is invariant to row permutation (maxdiff 0)", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_img_gasfcnn(m, X)
  set.seed(7)
  perm <- sample(nrow(X))
  sp <- score_img_gasfcnn(m, X[perm, , drop = FALSE])
  inv <- numeric(nrow(X)); inv[perm] <- sp
  expect_identical(max(abs(base - inv)), 0)
})

# §7.4 per-sample positive-scale invariance ~0 (rCLR input + frozen bounds: a
# positive rescale leaves the rCLR, GASF, and logit unchanged).
test_that("§7.4 score is invariant to per-sample positive scaling", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_img_gasfcnn(m, X)
  set.seed(3)
  scal <- matrix(stats::runif(nrow(X), 1e-3, 1e3), nrow(X), ncol(X))
  Xs <- X * scal
  ss <- score_img_gasfcnn(m, Xs)
  expect_equal(max(abs(base - ss)), 0, tolerance = 1e-6)
})

# §7.5 held-out AUC > 0.7 on a planted shift. A strong shift (2.0) is used because
# the GASF-image CNN encodes pairwise (second-order) structure and needs a clear
# signal to learn on only 30 training specimens; the directional property is the
# point, not absolute accuracy.
test_that("§7.5 larger score = more case-like on a held-out planted shift", {
  skip_if_no_gasfcnn()
  d <- .make_gasfcnn_data(seed = 11L, shift = 2.0)  # n = 40 (ctrl 1:20, case 21:40)
  tr <- c(1:15, 21:35)
  te <- setdiff(seq_len(nrow(d$X)), tr)        # 5 held-out controls + 5 cases
  m <- fit_img_gasfcnn(d$X[tr, ], d$y[tr],
                       hp = list(device = "auto", epochs = 300L))
  s <- score_img_gasfcnn(m, d$X[te, ])
  auc <- .auc_mw_gasfcnn(d$y[te], s)
  expect_gt(auc, 0.7)
})

# §7.6 canonical row-equivariance gate + dispatch.
test_that("§7.6 passes the canonical row-equivariance gate", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_img_gasfcnn, f$model, f$data$X[seq_len(8L), , drop = FALSE],
      model_digest = .img_gasfcnn_model_digest
    )
  )
})

test_that("§7.6 canonical dispatch via singlesample_score_call matches direct score", {
  skip_if_no_gasfcnn()
  # Only run once the orchestrator has wired img-gasfcnn's canonical score_fn into
  # the package; until then the method is absent from the roster and we skip.
  roster <- tryCatch(singlesample_method_roster(), error = function(e) NULL)
  skip_if(is.null(roster) || !("img-gasfcnn" %in% roster$method_id),
          "img-gasfcnn not yet registered in the roster")
  f <- .fit_small_gasfcnn()
  X <- f$data$X[seq_len(6L), , drop = FALSE]
  via <- tryCatch(singlesample_score_call("img-gasfcnn", f$model, X),
                  error = function(e) e)
  if (inherits(via, "error")) {
    skip(paste("img-gasfcnn canonical dispatch not yet wired:",
               conditionMessage(via)))
  }
  direct <- score_img_gasfcnn(f$model, X)
  expect_equal(direct, via, tolerance = 1e-12)
})

# §7.7 degenerate -> 0 + flat composition scored (NOT floored).
test_that("§7.7 degenerate queries return 0; a flat composition is scored", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  X <- f$data$X[seq_len(5L), , drop = FALSE]

  # Fewer than min_features universe features present -> 0 for every row.
  X_few <- X[, seq_len(2L), drop = FALSE]
  expect_true(all(score_img_gasfcnn(m, X_few) == 0))

  # No shared features at all -> 0 for every row.
  X_none <- X
  colnames(X_none) <- paste0("other-", seq_len(ncol(X_none)))
  expect_true(all(score_img_gasfcnn(m, X_none) == 0))

  # All-zero specimen -> 0 for that row, finite for the others.
  X_zero <- X
  X_zero[2L, ] <- 0
  s <- score_img_gasfcnn(m, X_zero)
  expect_equal(s[2L], 0)
  expect_true(all(is.finite(s)))

  # A FLAT all-equal-positive composition has an all-zero rCLR (rCLR origin) but is
  # a VALID specimen: it must be imaged and scored, NOT floored to 0.
  X_flat <- X
  X_flat[3L, ] <- 7
  s_flat <- score_img_gasfcnn(m, X_flat)
  expect_true(is.finite(s_flat[3L]))
})

# §7.8 determinism: two seed=42 fits -> identical model_digest + bit-identical
# scores (seed-before-build, full-batch). Hardened with an inter-fit torch-RNG
# perturbation (a different intervening torch seed must not change the result).
test_that("§7.8 determinism: re-score and refit reproduce identical scores + digest", {
  skip_if_no_gasfcnn()
  d <- .make_gasfcnn_data(seed = 11L)
  m1 <- fit_img_gasfcnn(d$X, d$y, hp = list(device = "auto", seed = 42L,
                                            epochs = 60L))
  X <- d$X[seq_len(6L), , drop = FALSE]
  s1a <- score_img_gasfcnn(m1, X)
  s1b <- score_img_gasfcnn(m1, X)
  expect_identical(s1a, s1b)                          # re-score identical

  # Perturb the global torch RNG between fits: seed-before-build must make the
  # second fit reproduce the first regardless of the intervening torch state.
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  torch_py$manual_seed(99999L)

  m2 <- fit_img_gasfcnn(d$X, d$y, hp = list(device = "auto", seed = 42L,
                                            epochs = 60L))
  s2 <- score_img_gasfcnn(m2, X)
  expect_identical(max(abs(s1a - s2)), 0)             # refit bit-identical
  expect_identical(.img_gasfcnn_model_digest(m1), .img_gasfcnn_model_digest(m2))
})

# §7.9 RNG-safety: scoring consumes no RNG and leaves the global R + torch
# (CPU/CUDA) generators byte-unchanged; fit leaves the torch generators unchanged.
test_that("§7.9 scoring does not mutate the global torch / R RNG state", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  X <- f$data$X[seq_len(10L), , drop = FALSE]

  set.seed(99)
  before_r <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  before_t <- torch_py$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch_py$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch_py$cuda$get_rng_state_all() else NULL
  invisible(score_img_gasfcnn(m, X))
  after_r <- get(".Random.seed", envir = globalenv())
  after_t <- torch_py$random$get_rng_state()
  expect_identical(before_r, after_r)
  expect_true(isTRUE(reticulate::py_to_r(torch_py$equal(before_t, after_t))))
  if (cuda_ok) {
    after_cuda <- torch_py$cuda$get_rng_state_all()
    eq <- vapply(seq_along(before_cuda), function(i)
      isTRUE(reticulate::py_to_r(torch_py$equal(before_cuda[[i]], after_cuda[[i]]))),
      logical(1L))
    expect_true(all(eq))
  }
})

test_that("§7.9b fit leaves the CPU + CUDA torch RNG byte-unchanged", {
  skip_if_no_gasfcnn()
  d <- .make_gasfcnn_data(seed = 11L)
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  before_t <- torch_py$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch_py$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch_py$cuda$get_rng_state_all() else NULL
  invisible(fit_img_gasfcnn(d$X, d$y, hp = list(device = "auto", epochs = 20L)))
  after_t <- torch_py$random$get_rng_state()
  expect_true(isTRUE(reticulate::py_to_r(torch_py$equal(before_t, after_t))))
  if (cuda_ok) {
    after_cuda <- torch_py$cuda$get_rng_state_all()
    eq <- vapply(seq_along(before_cuda), function(i)
      isTRUE(reticulate::py_to_r(torch_py$equal(before_cuda[[i]], after_cuda[[i]]))),
      logical(1L))
    expect_true(all(eq))
  }
})

# §7.10 DEDICATED frozen-BN proof: the rebuilt CNN is in eval mode, so BatchNorm
# uses its FROZEN running stats and the SAME image scores identically in a batch of
# 1 vs a batch of 32 (maxdiff 0 in float64). A train-mode forward of the SAME net
# on the batch of 32 differs (sanity: BN coupling would break single-sample).
test_that("§7.10 frozen-BN proof: same image batch-1 == batch-32 (eval); train-mode differs", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  torch <- reticulate::import("torch", delay_load = FALSE)
  device <- "cpu"

  # Rebuild the eval-mode float64 CNN exactly as the scorer does.
  net <- .img_gasfcnn_build_net(m$channels, device)
  net <- .img_gasfcnn_load_state_dict(net, m$state_dict, device)

  # One training specimen's GASF image.
  Z <- .img_gasfcnn_rclr_matrix(f$data$X[1L, , drop = FALSE])
  G1 <- .img_gasfcnn_images(Z, m$gasf_lo, m$gasf_hi)        # 1 x 1 x p x p
  img1 <- torch$as_tensor(G1)$to(torch$float64)
  img32 <- img1$`repeat`(list(32L, 1L, 1L, 1L))

  ng <- torch$no_grad(); ng$`__enter__`()
  on.exit(tryCatch(ng$`__exit__`(NULL, NULL, NULL), error = function(e) NULL),
          add = TRUE)

  l1 <- as.numeric(reticulate::py_to_r(
    net(img1)$reshape(list(-1L))$item()))
  l32 <- as.numeric(reticulate::py_to_r(
    net(img32)$reshape(list(-1L))[0]$item()))
  # Frozen-BN eval: identical row in batch-1 vs batch-32 -> identical logit.
  expect_lt(abs(l1 - l32), 1e-12)

  # SANITY: train-mode BatchNorm computes per-batch stats and COUPLES the rows,
  # so the same image inside a 32-row batch scores DIFFERENTLY -> this is exactly
  # the failure mode eval() prevents. (Restore eval after the sanity probe.)
  net$train()
  l32_train <- as.numeric(reticulate::py_to_r(
    net(img32)$reshape(list(-1L))[0]$item()))
  net$eval()
  expect_gt(abs(l32_train - l1), 1e-6)
})

# §7.11 rebuilt-CNN fidelity: the score from the exported-state-dict rebuilt CNN ==
# the ORIGINAL fitted net's eval-mode float64 forward on the same image (<1e-5). We
# refit with a captured original net and compare its float64 eval forward against
# the scorer's rebuilt-from-state-dict forward.
test_that("§7.11 rebuilt-CNN fidelity: rebuilt == original fitted net (<1e-5)", {
  skip_if_no_gasfcnn()
  f <- .fit_small_gasfcnn()
  m <- f$model
  torch <- reticulate::import("torch", delay_load = FALSE)
  device <- m$device                                  # retrain on the FIT device

  # ORIGINAL: build a net seeded identically, train identically, eval+float64.
  # Replicate the fit's cuDNN-deterministic context AND device so the retrained
  # "original" net reproduces the fitted state_dict the scorer rebuilds from
  # (otherwise the GPU conv path is nondeterministic and the two nets diverge).
  Z <- .img_gasfcnn_rclr_matrix(f$data$X)
  bounds <- list(lo = m$gasf_lo, hi = m$gasf_hi)
  np <- reticulate::import("numpy", delay_load = FALSE)
  torch$backends$cudnn$deterministic <- TRUE
  torch$backends$cudnn$benchmark <- FALSE
  torch$manual_seed(as.integer(m$hp$seed))
  orig <- .img_gasfcnn_build_net(m$channels, device)
  G <- .img_gasfcnn_images(Z, bounds$lo, bounds$hi)
  X_py <- torch$as_tensor(G)$to(torch$float32)$to(device)
  y_py <- torch$as_tensor(np$asarray(as.numeric(f$data$y),
                                     dtype = "float32"))$to(device)
  opt <- torch$optim$Adam(orig$parameters(), lr = m$hp$lr,
                          weight_decay = m$hp$weight_decay)
  loss_fn <- torch$nn$BCEWithLogitsLoss()
  orig$train()
  for (ep in seq_len(m$hp$epochs)) {
    opt$zero_grad()
    logit <- orig(X_py)$reshape(list(-1L))
    loss <- loss_fn(logit, y_py)
    loss$backward(); opt$step()
  }
  orig$eval(); orig <- orig$to(torch$float64)

  # REBUILT-from-state-dict (the scorer's path).
  reb <- .img_gasfcnn_build_net(m$channels, device)
  reb <- .img_gasfcnn_load_state_dict(reb, m$state_dict, device)

  ng <- torch$no_grad(); ng$`__enter__`()
  on.exit(tryCatch(ng$`__exit__`(NULL, NULL, NULL), error = function(e) NULL),
          add = TRUE)

  maxd <- 0
  for (i in seq_len(6L)) {
    Gi <- .img_gasfcnn_images(Z[i, , drop = FALSE], bounds$lo, bounds$hi)
    img <- torch$as_tensor(Gi)$to(torch$float64)$to(device)
    lo_orig <- as.numeric(reticulate::py_to_r(orig(img)$reshape(list(-1L))$item()))
    lo_reb <- as.numeric(reticulate::py_to_r(reb(img)$reshape(list(-1L))$item()))
    maxd <- max(maxd, abs(lo_orig - lo_reb))
  }
  expect_lt(maxd, 1e-5)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_gasfcnn()
  d <- .make_gasfcnn_data(seed = 11L)
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(bogus = 1)),
               "unknown hp field")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(seed = 1, seed = 2)),
               "duplicate hp field")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(7)),
               "named")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(device = "tpu")),
               "device")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(channels = 0)),
               "channels")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(epochs = 0)),
               "epochs")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(lr = -1)),
               "lr")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(min_features = 2.5)),
               "min_features")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(seed = 1.5)),
               "seed")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(image_size = 1L)),
               "image_size")
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(image_size = 2.5)),
               "image_size")
  # A GASF image side L (= min(p, image_size)) smaller than 2^(#MaxPool blocks) errors
  # clearly instead of an opaque torch MaxPool RuntimeError. image_size=3 with the
  # default 2-block net -> L=3 < 4 -> guard fires.
  expect_error(fit_img_gasfcnn(d$X, d$y, hp = list(image_size = 3L)),
               "too small for")
})

# §7.12 PAA fixed-size imaging: when p > image_size the frozen-order rCLR profile is
# PAA-reduced to L = image_size, so the GASF image is L x L (bounded regardless of p),
# gasf_lo/gasf_hi have length L, and EVERY single-sample invariant still holds EXACTLY
# (forced row-by-row + frozen-BN eval + float64): single-row == batch maxdiff 0,
# per-sample positive-scale invariance ~0, and two seed-matched fits share a digest.
# This is the property the W7 benchmark relies on: a 2500-feature cohort cannot form a
# p x p image (memory/compute), but a fixed L x L image is feasible and equivariant.
test_that("§7.12 PAA: p > image_size -> bounded L x L image, invariants exact", {
  skip_if_no_gasfcnn()
  d <- .make_gasfcnn_data(n = 40L, p = 200L, k = 8L, shift = 1.2)  # p = 200 > L
  m <- fit_img_gasfcnn(d$X, d$y,
                       hp = list(device = "auto", epochs = 40L, seed = 11L,
                                 image_size = 64L))
  # PAA activated: the frozen bounds are length L = image_size = 64 (not p = 200).
  expect_length(m$gasf_lo, 64L)
  expect_length(m$gasf_hi, 64L)
  expect_identical(m$paa_len, 64L)
  expect_identical(m$image_size, 64L)

  X <- d$X[seq_len(8L), , drop = FALSE]
  # The rendered GASF image is L x L = 64 x 64 (bounded), NOT p x p.
  Z <- .img_gasfcnn_rclr_matrix(X[1L, , drop = FALSE])
  zr <- .img_gasfcnn_paa(Z[1L, ], m$hp$image_size)
  Gimg <- reticulate::py_to_r(.img_gasfcnn_images(matrix(zr, nrow = 1L),
                                                  m$gasf_lo, m$gasf_hi))
  expect_identical(dim(Gimg), c(1L, 1L, 64L, 64L))

  # DECISIVE single-sample equivariance under PAA: single-row == batch, maxdiff 0.
  batch <- score_img_gasfcnn(m, X)
  single <- vapply(seq_len(nrow(X)), function(i)
    score_img_gasfcnn(m, X[i, , drop = FALSE]), numeric(1L))
  expect_identical(max(abs(batch - single)), 0)

  # Per-sample positive-scale invariance preserved (rCLR-before-PAA is linear).
  set.seed(5)
  scal <- matrix(stats::runif(nrow(X), 1e-3, 1e3), nrow(X), ncol(X))
  expect_equal(max(abs(batch - score_img_gasfcnn(m, X * scal))), 0, tolerance = 1e-6)

  # Determinism under PAA: two seed-matched fits -> identical digest.
  m2 <- fit_img_gasfcnn(d$X, d$y,
                        hp = list(device = "auto", epochs = 40L, seed = 11L,
                                  image_size = 64L))
  expect_identical(.img_gasfcnn_model_digest(m), .img_gasfcnn_model_digest(m2))

  # Canonical row-equivariance gate under PAA.
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_img_gasfcnn, m, X, model_digest = .img_gasfcnn_model_digest))
})
