library(testthat)

# Runs inside the package test harness (devtools::test / R CMD check), where
# fit_tabdpt / score_tabdpt / .tabdpt_model_digest and all package helpers are
# already loaded. TabDPT inference needs the Python `tabdpt` module via
# reticulate; tests that fit/score are skipped when tabdpt is unavailable.
skip_if_no_tabdpt <- function() {
  venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
  if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = venv)
  }
  # The pinned checkpoint is fetched once (token-free) then cached; load from
  # cache only, so the test needs no network / TLS (the brew-python TLS stack is
  # not load-path-stable under R) and is hermetic.
  Sys.setenv(HF_HUB_OFFLINE = "1")
  testthat::skip_if_not_installed("reticulate")
  cuda_ok <- tryCatch(suppressWarnings({
    torch <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  }), error = function(e) FALSE)
  testthat::skip_if_not(
    cuda_ok, "CUDA unavailable; TabDPT CPU fallback is too slow for R CMD check")
  # Actually IMPORT tabdpt (not just module-available): tabdpt's import chain
  # pulls in Python's _ssl, whose brew-OpenSSL build needs the brew libcrypto on
  # the loader path BEFORE R loads the system one. On this host run R with
  #   LD_PRELOAD=/home/linuxbrew/.linuxbrew/opt/openssl@3/lib/libcrypto.so.3:.../libssl.so.3
  # (the brew-torch bridge). Where that bridge is inactive the import fails; skip
  # gracefully so the package suite stays green on unbridged hosts/CI.
  # Importing Python's `ssl` eagerly loads `_ssl` (the brew-OpenSSL .so), which is
  # the exact thing that fails without the bridge -- tabdpt's __init__ does NOT
  # pull it, but the checkpoint load at fit time does, so probe `ssl` here.
  ok <- tryCatch({
    reticulate::import("ssl")
    reticulate::import("tabdpt")
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(
    ok, "tabdpt not importable (module missing or brew-OpenSSL _ssl bridge inactive)")
}

# Synthetic generator with a CASE-ONLY location shift in a feature block. The
# first `k` of `p` log-abundance features are elevated by `shift` for case
# samples; after exp() and the per-sample rCLR centring those features carry the
# discriminative signal TabDPT learns in-context. Abundances are exp() so the
# input is strictly positive (compositional).
# n defaults small (40) and the tests fit/score on device = "auto" (the H100 GPU
# here; falls back to CPU elsewhere) because the row-by-row TabDPT forward pass is
# far too slow on CPU for the suite. The DECISIVE single-sample properties
# (single-row == batch, determinism, row-equivariance, scale-invariance) are all
# bit-EXACT on GPU here (identical inputs -> identical n=1 forward); only the
# per-sample scale-invariance test keeps a 1e-5 tolerance to absorb any GPU cuBLAS
# float noise across builds. The scorer's own DEFAULT device stays "cpu" (exact /
# reproducible deployment); these tests opt into the GPU purely for speed. TabDPT
# is designed for in-context data, so a 40-row context with a strong planted shift
# is ample.
.make_tabdpt_data <- function(n = 40L, p = 20L, k = 6L, shift = 1.2,
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

.auc_mw_tabdpt <- function(y, score) {
  case <- score[y == 1L]
  ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# A small, fast GPU fit reused across tests.
.fit_small_tabdpt <- function(seed = 11L) {
  d <- .make_tabdpt_data(seed = seed)
  list(model = fit_tabdpt(d$X, d$y, hp = list(device = "auto")), data = d)
}


test_that("fit returns a well-formed tabdpt_model", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
  m <- f$model
  expect_s3_class(m, "tabdpt_model")
  expect_equal(m$feature_universe, colnames(f$data$X))
  expect_true(m$case_col %in% c(1L, 2L))
  expect_equal(m$classes, c(0L, 1L))
  expect_equal(m$checkpoint, "tabdpt1_2.safetensors")
  expect_equal(nrow(m$context_rclr), m$n_context)
})

test_that("score is finite, length-correct, and in [0, 1]", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
  s <- score_tabdpt(f$model, f$data$X)
  expect_length(s, nrow(f$data$X))
  expect_true(all(is.finite(s)))
  expect_true(all(s >= 0 & s <= 1))
  expect_type(s, "double")
})

test_that("larger score = more case-like on a held-out planted shift", {
  skip_if_no_tabdpt()
  d <- .make_tabdpt_data(seed = 11L)          # n = 40 (controls 1:20, cases 21:40)
  tr <- c(1:15, 21:35)
  te <- setdiff(seq_len(nrow(d$X)), tr)       # 5 held-out controls + 5 cases
  m <- fit_tabdpt(d$X[tr, ], d$y[tr], hp = list(device = "auto"))
  s <- score_tabdpt(m, d$X[te, ])
  auc <- .auc_mw_tabdpt(d$y[te], s)
  expect_gt(auc, 0.5)
  expect_gt(auc, 0.8)   # strong planted signal; clearly directional
})

test_that("DECISIVE: single-row score == batch score (conditional single-sample)", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  batch <- score_tabdpt(m, X)
  single <- vapply(seq_len(nrow(X)), function(i) {
    score_tabdpt(m, X[i, , drop = FALSE])
  }, numeric(1L))
  expect_equal(max(abs(batch - single)), 0, tolerance = 1e-6)
})

test_that("score is invariant to row permutation", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_tabdpt(m, X)
  set.seed(7)
  perm <- sample(nrow(X))
  sp <- score_tabdpt(m, X[perm, , drop = FALSE])
  inv <- numeric(nrow(X)); inv[perm] <- sp
  expect_equal(max(abs(base - inv)), 0, tolerance = 1e-8)
})

test_that("score is invariant to per-sample positive scaling", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_tabdpt(m, X)
  set.seed(3)
  scal <- matrix(stats::runif(nrow(X), 1e-3, 1e3), nrow(X), ncol(X))
  Xs <- X * scal
  ss <- score_tabdpt(m, Xs)
  # rCLR is exactly scale-invariant; any residual is GPU cuBLAS float noise
  # (measured 0 on this build, <=~2e-6 worst-case on cuda), so the tolerance is
  # 1e-5 (still far tighter than any score-meaningful difference).
  expect_equal(max(abs(base - ss)), 0, tolerance = 1e-5)
})

test_that("passes the canonical row-equivariance gate", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_tabdpt, f$model, f$data$X[seq_len(8L), , drop = FALSE],
      model_digest = .tabdpt_model_digest
    )
  )
})

test_that("degenerate queries return the neutral 0.5", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
  m <- f$model
  X <- f$data$X[seq_len(5L), , drop = FALSE]

  # Fewer than min_features universe features present -> 0.5 for every row.
  X_few <- X[, seq_len(2L), drop = FALSE]
  expect_true(all(score_tabdpt(m, X_few) == 0.5))

  # No shared features at all -> 0.5 for every row.
  X_none <- X
  colnames(X_none) <- paste0("other-", seq_len(ncol(X_none)))
  expect_true(all(score_tabdpt(m, X_none) == 0.5))

  # All-zero specimen -> 0.5 for that row, finite for the others.
  X_zero <- X
  X_zero[2L, ] <- 0
  s <- score_tabdpt(m, X_zero)
  expect_equal(s[2L], 0.5)
  expect_true(all(is.finite(s)))
})

test_that("determinism: re-score and refit reproduce identical scores", {
  skip_if_no_tabdpt()
  d <- .make_tabdpt_data(seed = 11L)
  m1 <- fit_tabdpt(d$X, d$y, hp = list(device = "auto", seed = 42L))
  X <- d$X[seq_len(6L), , drop = FALSE]
  s1a <- score_tabdpt(m1, X)
  s1b <- score_tabdpt(m1, X)
  expect_identical(s1a, s1b)                       # re-score identical

  m2 <- fit_tabdpt(d$X, d$y, hp = list(device = "auto", seed = 42L))
  s2 <- score_tabdpt(m2, X)
  expect_equal(max(abs(s1a - s2)), 0, tolerance = 1e-8)  # refit identical

  # Same frozen score-determining state -> identical digest across two fits.
  expect_identical(.tabdpt_model_digest(m1), .tabdpt_model_digest(m2))
})

test_that("scoring does not mutate the global torch / R RNG state", {
  skip_if_no_tabdpt()
  f <- .fit_small_tabdpt()
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
  invisible(score_tabdpt(m, X))
  after_r <- get(".Random.seed", envir = globalenv())
  after_t <- torch_py$random$get_rng_state()
  expect_identical(before_r, after_r)
  expect_true(isTRUE(reticulate::py_to_r(torch_py$equal(before_t, after_t))))
  # Scoring must also leave the CUDA RNG byte-unchanged on a GPU; no-op on CPU.
  if (cuda_ok) {
    after_cuda <- torch_py$cuda$get_rng_state_all()
    eq <- vapply(seq_along(before_cuda), function(i)
      isTRUE(reticulate::py_to_r(torch_py$equal(before_cuda[[i]], after_cuda[[i]]))),
      logical(1L))
    expect_true(all(eq))
  }
})

test_that("canonical dispatch via singlesample_score_call matches direct score", {
  skip_if_no_tabdpt()
  # Only run once the orchestrator has wired tab-tabdpt's canonical score_fn /
  # adapter into the package; until then dispatch raises a "not yet implemented"
  # / "register a canonical adapter" error, which we treat as a skip.
  roster <- tryCatch(singlesample_method_roster(), error = function(e) NULL)
  skip_if(is.null(roster) || !("tab-tabdpt" %in% roster$method_id),
          "tab-tabdpt not yet registered in the roster")
  f <- .fit_small_tabdpt()
  X <- f$data$X[seq_len(6L), , drop = FALSE]
  via <- tryCatch(singlesample_score_call("tab-tabdpt", f$model, X),
                  error = function(e) e)
  if (inherits(via, "error")) {
    skip(paste("tab-tabdpt canonical dispatch not yet wired:",
               conditionMessage(via)))
  }
  direct <- score_tabdpt(f$model, X)
  expect_equal(direct, via, tolerance = 1e-12)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_tabdpt()
  d <- .make_tabdpt_data(seed = 11L)
  expect_error(fit_tabdpt(d$X, d$y, hp = list(bogus = 1)),
               "unknown hp field")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(seed = 1, seed = 2)),
               "duplicate hp field")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(7)),
               "named")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(device = "tpu")),
               "device")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(temperature = 0)),
               "temperature")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(normalizer = "bogus")),
               "normalizer")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(max_context = 1)),
               "max_context")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(min_features = 2.5)),
               "min_features")
  expect_error(fit_tabdpt(d$X, d$y, hp = list(seed = 1.5)),
               "seed")
})

test_that("context cap produces a stratified, seeded, deterministic subsample", {
  skip_if_no_tabdpt()
  d <- .make_tabdpt_data(n = 120L, seed = 11L)
  m <- fit_tabdpt(d$X, d$y, hp = list(device = "auto", max_context = 40L))
  expect_lte(m$n_context, 40L)
  expect_true(any(m$context_y == 1L) && any(m$context_y == 0L))
  # Deterministic: a second capped fit yields the identical frozen context.
  m2 <- fit_tabdpt(d$X, d$y, hp = list(device = "auto", max_context = 40L))
  expect_identical(m$context_rclr, m2$context_rclr)
  expect_identical(m$context_y, m2$context_y)
})

test_that("context cap keeps both classes under extreme imbalance + tiny cap", {
  # Regression (pure R, no tabdpt needed) for two bugs that dropped a class:
  # (a) proportional rounding could give the whole cap to the majority class with
  # no reserved minority slot; (b) sample(idx, 1) on a SINGLE-element class index
  # hit base R's sample(x, 1) -> sample(1:x, 1) gotcha and picked a wrong row.
  sub <- .tabdpt_context_subsample
  combos <- list(c(99L, 1L, 2L), c(1L, 99L, 2L), c(199L, 1L, 2L),
                 c(98L, 2L, 2L), c(2L, 98L, 3L), c(50L, 50L, 2L))
  for (cc in combos) {
    y <- c(rep(1L, cc[1L]), rep(0L, cc[2L]))
    keep <- sub(length(y), y, cc[3L], 42L)
    kept <- y[keep]
    expect_true(0L %in% kept && 1L %in% kept,
                info = sprintf("%dcase:%dctrl cap=%d dropped a class",
                               cc[1L], cc[2L], cc[3L]))
    expect_lte(length(keep), cc[3L])
  }
  y <- c(rep(1L, 99L), 0L)
  expect_identical(sub(100L, y, 2L, 42L), sub(100L, y, 2L, 42L))
})

test_that("fit leaves the CUDA torch RNG byte-unchanged (GPU only)", {
  skip_if_no_tabdpt()
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  skip_if_not(isTRUE(reticulate::py_to_r(torch_py$cuda$is_available())),
              "CUDA not available")
  d <- .make_tabdpt_data(seed = 11L)
  before <- torch_py$cuda$get_rng_state_all()
  invisible(fit_tabdpt(d$X, d$y, hp = list(device = "cuda")))
  after <- torch_py$cuda$get_rng_state_all()
  eq <- vapply(seq_along(before), function(i)
    isTRUE(reticulate::py_to_r(torch_py$equal(before[[i]], after[[i]]))),
    logical(1L))
  expect_true(all(eq))
})

test_that("batch scoring is companion-invariant and digest-neutral", {
  skip_if_no_tabdpt()
  expect_false(.tabdpt_resolve_hp(list())$score_batch)
  expect_error(.tabdpt_resolve_hp(list(score_batch = 1)), "score_batch")
  expect_true(.tabdpt_resolve_hp(list(score_batch = TRUE))$score_batch)

  d <- .make_tabdpt_data(n = 120L, p = 300L, k = 12L, shift = 2.5,
                         sd = 0.35, seed = 2027L)
  hp <- list(device = "auto", max_context = 120L, seed = 42L)
  m <- fit_tabdpt(d$X, d$y, hp = hp)

  query <- 8L
  comp_a <- 9L:24L
  comp_b <- 61L:76L
  X_a <- d$X[c(query, comp_a), , drop = FALSE]
  X_b <- d$X[c(query, comp_b), , drop = FALSE]
  rownames(X_a)[1L] <- "query"
  rownames(X_b)[1L] <- "query"
  zero <- matrix(0, nrow = 1L, ncol = ncol(X_a),
                 dimnames = list("zero", colnames(X_a)))
  X_a_zero <- rbind(X_a, zero)

  digest_before <- .tabdpt_model_digest(m)
  invisible(score_tabdpt(m, X_a))
  expect_identical(digest_before, .tabdpt_model_digest(m))

  mb <- m
  mb$hp$score_batch <- TRUE
  batch_a <- score_tabdpt(mb, X_a_zero)
  batch_b <- score_tabdpt(mb, X_b)
  expect_length(batch_a, nrow(X_a_zero))
  expect_equal(abs(batch_a[1L] - batch_b[1L]), 0, tolerance = 1e-12)
  expect_equal(batch_a[nrow(X_a_zero)], 0.5)
  expect_identical(.tabdpt_model_digest(m), .tabdpt_model_digest(mb))
})

test_that("multi-chunk batch scoring avoids singleton tails", {
  skip_if_no_tabdpt()
  chunks_1025 <- .tabdpt_score_chunks(seq_len(1025L), chunk_size = 1024L)
  chunks_2049 <- .tabdpt_score_chunks(seq_len(2049L), chunk_size = 1024L)
  for (chunks in list(chunks_1025, chunks_2049)) {
    sizes <- lengths(chunks)
    expect_true(all(sizes >= 2L))
    expect_lte(max(sizes) - min(sizes), 1L)
    expect_identical(unname(unlist(chunks)), seq_len(sum(sizes)))
  }

  d <- .make_tabdpt_data(n = 40L, p = 20L, k = 6L, shift = 2.5,
                         sd = 0.35, seed = 2028L)
  m <- fit_tabdpt(d$X, d$y, hp = list(device = "auto", max_context = 40L,
                                      seed = 42L))
  mb <- m
  mb$hp$score_batch <- TRUE
  q_idx <- rep(seq_len(nrow(d$X)), length.out = 1026L)
  Xq <- d$X[q_idx, , drop = FALSE]
  rownames(Xq) <- paste0("Q", seq_len(nrow(Xq)))

  s_1025 <- score_tabdpt(mb, Xq[seq_len(1025L), , drop = FALSE])
  s_1026 <- score_tabdpt(mb, Xq[seq_len(1026L), , drop = FALSE])
  expect_equal(s_1025[1025L], s_1026[1025L], tolerance = 1e-3)
  expect_identical(.tabdpt_model_digest(m), .tabdpt_model_digest(mb))
})
