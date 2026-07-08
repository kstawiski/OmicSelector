library(testthat)

# Runs inside the package test harness (devtools::test / R CMD check), where
# fit_tabicl / score_tabicl / .tabicl_model_digest and all package helpers are
# already loaded. TabICL inference needs the Python `tabicl` module via
# reticulate; tests that fit/score are skipped when tabicl is unavailable.
skip_if_no_tabicl <- function() {
  venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
  if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = venv)
  }
  testthat::skip_if_not_installed("reticulate")
  cuda_ok <- tryCatch(suppressWarnings({
    torch <- reticulate::import("torch", delay_load = FALSE)
    isTRUE(reticulate::py_to_r(torch$cuda$is_available()))
  }), error = function(e) FALSE)
  testthat::skip_if_not(
    cuda_ok, "CUDA unavailable; TabICL CPU fallback is too slow for R CMD check")
  testthat::skip_if_not(reticulate::py_module_available("tabicl"),
                        "Python module 'tabicl' not available")
}

# Synthetic generator with a CASE-ONLY location shift in a feature block. The
# first `k` of `p` log-abundance features are elevated by `shift` for case
# samples; after exp() and the per-sample rCLR centring those features carry the
# discriminative signal TabICL learns in-context. Abundances are exp() so the
# input is strictly positive (compositional).
# n defaults small (40) and the tests fit/score on device = "auto" (the H100 GPU
# here; falls back to CPU elsewhere) because the row-by-row TabICL forward pass is
# far too slow on CPU for the suite. The DECISIVE single-sample properties
# (single-row == batch, determinism, row-equivariance, scale-invariance) are all
# bit-EXACT on GPU here (identical inputs -> identical n=1 forward); only the
# per-sample scale-invariance test keeps a 1e-5 tolerance to absorb any GPU cuBLAS
# float noise across builds. The scorer's own DEFAULT device stays "cpu" (exact /
# reproducible deployment); these tests opt into the GPU purely for speed. TabICL
# is designed for in-context data, so a 40-row context with a strong planted shift
# is ample.
.make_tabicl_data <- function(n = 40L, p = 20L, k = 6L, shift = 1.2,
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

.auc_mw_tabicl <- function(y, score) {
  case <- score[y == 1L]
  ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# A small, fast GPU fit reused across tests.
.fit_small_tabicl <- function(seed = 11L) {
  d <- .make_tabicl_data(seed = seed)
  list(model = fit_tabicl(d$X, d$y, hp = list(device = "auto")), data = d)
}


test_that("fit returns a well-formed tabicl_model", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
  m <- f$model
  expect_s3_class(m, "tabicl_model")
  expect_equal(m$feature_universe, colnames(f$data$X))
  expect_true(m$case_col %in% c(1L, 2L))
  expect_equal(m$classes, c(0L, 1L))
  expect_equal(m$checkpoint, "tabicl-classifier-v2-20260212.ckpt")
  expect_equal(nrow(m$context_rclr), m$n_context)
})

test_that("score is finite, length-correct, and in [0, 1]", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
  s <- score_tabicl(f$model, f$data$X)
  expect_length(s, nrow(f$data$X))
  expect_true(all(is.finite(s)))
  expect_true(all(s >= 0 & s <= 1))
  expect_type(s, "double")
})

test_that("larger score = more case-like on a held-out planted shift", {
  skip_if_no_tabicl()
  d <- .make_tabicl_data(seed = 11L)          # n = 40 (controls 1:20, cases 21:40)
  tr <- c(1:15, 21:35)
  te <- setdiff(seq_len(nrow(d$X)), tr)       # 5 held-out controls + 5 cases
  m <- fit_tabicl(d$X[tr, ], d$y[tr], hp = list(device = "auto"))
  s <- score_tabicl(m, d$X[te, ])
  auc <- .auc_mw_tabicl(d$y[te], s)
  expect_gt(auc, 0.5)
  expect_gt(auc, 0.8)   # strong planted signal; clearly directional
})

test_that("DECISIVE: single-row score == batch score (conditional single-sample)", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  batch <- score_tabicl(m, X)
  single <- vapply(seq_len(nrow(X)), function(i) {
    score_tabicl(m, X[i, , drop = FALSE])
  }, numeric(1L))
  expect_equal(max(abs(batch - single)), 0, tolerance = 1e-6)
})

test_that("score is invariant to row permutation", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_tabicl(m, X)
  set.seed(7)
  perm <- sample(nrow(X))
  sp <- score_tabicl(m, X[perm, , drop = FALSE])
  inv <- numeric(nrow(X)); inv[perm] <- sp
  expect_equal(max(abs(base - inv)), 0, tolerance = 1e-8)
})

test_that("score is invariant to per-sample positive scaling", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_tabicl(m, X)
  set.seed(3)
  scal <- matrix(stats::runif(nrow(X), 1e-3, 1e3), nrow(X), ncol(X))
  Xs <- X * scal
  ss <- score_tabicl(m, Xs)
  # rCLR is exactly scale-invariant; any residual is GPU cuBLAS float noise
  # (measured 0 on this build, <=~2e-6 worst-case on cuda), so the tolerance is
  # 1e-5 (still far tighter than any score-meaningful difference).
  expect_equal(max(abs(base - ss)), 0, tolerance = 1e-5)
})

test_that("passes the canonical row-equivariance gate", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_tabicl, f$model, f$data$X[seq_len(8L), , drop = FALSE],
      model_digest = .tabicl_model_digest
    )
  )
})

test_that("degenerate queries return the neutral 0.5", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
  m <- f$model
  X <- f$data$X[seq_len(5L), , drop = FALSE]

  # Fewer than min_features universe features present -> 0.5 for every row.
  X_few <- X[, seq_len(2L), drop = FALSE]
  expect_true(all(score_tabicl(m, X_few) == 0.5))

  # No shared features at all -> 0.5 for every row.
  X_none <- X
  colnames(X_none) <- paste0("other-", seq_len(ncol(X_none)))
  expect_true(all(score_tabicl(m, X_none) == 0.5))

  # All-zero specimen -> 0.5 for that row, finite for the others.
  X_zero <- X
  X_zero[2L, ] <- 0
  s <- score_tabicl(m, X_zero)
  expect_equal(s[2L], 0.5)
  expect_true(all(is.finite(s)))
})

test_that("determinism: re-score and refit reproduce identical scores", {
  skip_if_no_tabicl()
  d <- .make_tabicl_data(seed = 11L)
  m1 <- fit_tabicl(d$X, d$y, hp = list(device = "auto", seed = 42L))
  X <- d$X[seq_len(6L), , drop = FALSE]
  s1a <- score_tabicl(m1, X)
  s1b <- score_tabicl(m1, X)
  expect_identical(s1a, s1b)                       # re-score identical

  m2 <- fit_tabicl(d$X, d$y, hp = list(device = "auto", seed = 42L))
  s2 <- score_tabicl(m2, X)
  expect_equal(max(abs(s1a - s2)), 0, tolerance = 1e-8)  # refit identical

  # Same frozen score-determining state -> identical digest across two fits.
  expect_identical(.tabicl_model_digest(m1), .tabicl_model_digest(m2))
})

test_that("scoring does not mutate the global torch / R RNG state", {
  skip_if_no_tabicl()
  f <- .fit_small_tabicl()
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
  invisible(score_tabicl(m, X))
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
  skip_if_no_tabicl()
  # Only run once the orchestrator has wired tab-tabicl's canonical score_fn /
  # adapter into the package; until then dispatch raises a "not yet implemented"
  # / "register a canonical adapter" error, which we treat as a skip.
  roster <- tryCatch(singlesample_method_roster(), error = function(e) NULL)
  skip_if(is.null(roster) || !("tab-tabicl" %in% roster$method_id),
          "tab-tabicl not yet registered in the roster")
  f <- .fit_small_tabicl()
  X <- f$data$X[seq_len(6L), , drop = FALSE]
  via <- tryCatch(singlesample_score_call("tab-tabicl", f$model, X),
                  error = function(e) e)
  if (inherits(via, "error")) {
    skip(paste("tab-tabicl canonical dispatch not yet wired:",
               conditionMessage(via)))
  }
  direct <- score_tabicl(f$model, X)
  expect_equal(direct, via, tolerance = 1e-12)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_tabicl()
  d <- .make_tabicl_data(seed = 11L)
  expect_error(fit_tabicl(d$X, d$y, hp = list(bogus = 1)),
               "unknown hp field")
  expect_error(fit_tabicl(d$X, d$y, hp = list(seed = 1, seed = 2)),
               "duplicate hp field")
  expect_error(fit_tabicl(d$X, d$y, hp = list(7)),
               "named")
  expect_error(fit_tabicl(d$X, d$y, hp = list(device = "tpu")),
               "device")
  expect_error(fit_tabicl(d$X, d$y, hp = list(n_estimators = 0)),
               "n_estimators")
  expect_error(fit_tabicl(d$X, d$y, hp = list(max_context = 1)),
               "max_context")
  expect_error(fit_tabicl(d$X, d$y, hp = list(min_features = 2.5)),
               "min_features")
  expect_error(fit_tabicl(d$X, d$y, hp = list(seed = 1.5)),
               "seed")
})

test_that("context cap produces a stratified, seeded, deterministic subsample", {
  skip_if_no_tabicl()
  d <- .make_tabicl_data(n = 120L, seed = 11L)
  m <- fit_tabicl(d$X, d$y, hp = list(device = "auto", max_context = 40L))
  expect_lte(m$n_context, 40L)
  expect_true(any(m$context_y == 1L) && any(m$context_y == 0L))
  # Deterministic: a second capped fit yields the identical frozen context.
  m2 <- fit_tabicl(d$X, d$y, hp = list(device = "auto", max_context = 40L))
  expect_identical(m$context_rclr, m2$context_rclr)
  expect_identical(m$context_y, m2$context_y)
})

test_that("context cap keeps both classes under extreme imbalance + tiny cap", {
  # Regression (pure R, no tabicl needed) for two bugs that dropped a class:
  # (a) proportional rounding could give the whole cap to the majority class with
  # no reserved minority slot; (b) sample(idx, 1) on a SINGLE-element class index
  # hit base R's sample(x, 1) -> sample(1:x, 1) gotcha and picked a wrong row.
  sub <- .tabicl_context_subsample
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
  skip_if_no_tabicl()
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  skip_if_not(isTRUE(reticulate::py_to_r(torch_py$cuda$is_available())),
              "CUDA not available")
  d <- .make_tabicl_data(seed = 11L)
  before <- torch_py$cuda$get_rng_state_all()
  invisible(fit_tabicl(d$X, d$y, hp = list(device = "cuda")))
  after <- torch_py$cuda$get_rng_state_all()
  eq <- vapply(seq_along(before), function(i)
    isTRUE(reticulate::py_to_r(torch_py$equal(before[[i]], after[[i]]))),
    logical(1L))
  expect_true(all(eq))
})

test_that("batch scoring equals row-by-row (>500 features)", {
  skip_if_no_tabicl()
  expect_true(.tabicl_resolve_hp(list(score_batch = TRUE))$score_batch)
  expect_error(.tabicl_resolve_hp(list(score_batch = 1)), "score_batch")

  d <- .make_tabicl_data(n = 64L, p = 700L, k = 12L, shift = 2.5,
                         sd = 0.35, seed = 813L)
  hp <- list(device = "auto", n_estimators = 2L, seed = 42L)
  m <- fit_tabicl(d$X, d$y, hp = hp)

  query <- c(seq_len(12L), 32L + seq_len(12L))
  Xq <- d$X[query, , drop = FALSE]
  yq <- d$y[query]
  Xq <- rbind(
    Xq,
    zero = matrix(0, nrow = 1L, ncol = ncol(Xq), dimnames = list(NULL, colnames(Xq)))
  )

  row <- score_tabicl(m, Xq)
  mb <- m
  mb$hp$score_batch <- TRUE
  batch <- score_tabicl(mb, Xq)

  expect_equal(batch, row, tolerance = 1e-6)
  expect_equal(batch[nrow(Xq)], 0.5)
  expect_gt(.auc_mw_tabicl(yq, batch[seq_along(yq)]), 0.5)
  expect_identical(.tabicl_model_digest(m), .tabicl_model_digest(mb))
})

test_that("multi-chunk batch scoring avoids singleton tails", {
  skip_if_no_tabicl()
  chunks_1025 <- .tabicl_score_chunks(seq_len(1025L), chunk_size = 1024L)
  chunks_2049 <- .tabicl_score_chunks(seq_len(2049L), chunk_size = 1024L)
  for (chunks in list(chunks_1025, chunks_2049)) {
    sizes <- lengths(chunks)
    expect_true(all(sizes >= 2L))
    expect_lte(max(sizes) - min(sizes), 1L)
    expect_identical(unname(unlist(chunks)), seq_len(sum(sizes)))
  }

  d <- .make_tabicl_data(n = 40L, p = 20L, k = 6L, shift = 2.5,
                         sd = 0.35, seed = 814L)
  m <- fit_tabicl(d$X, d$y, hp = list(device = "auto", n_estimators = 2L,
                                      max_context = 40L, seed = 42L))
  mb <- m
  mb$hp$score_batch <- TRUE
  q_idx <- rep(seq_len(nrow(d$X)), length.out = 1026L)
  Xq <- d$X[q_idx, , drop = FALSE]
  rownames(Xq) <- paste0("Q", seq_len(nrow(Xq)))

  s_1025 <- score_tabicl(mb, Xq[seq_len(1025L), , drop = FALSE])
  s_1026 <- score_tabicl(mb, Xq[seq_len(1026L), , drop = FALSE])
  expect_identical(unname(s_1025[1025L]), unname(s_1026[1025L]))
  expect_identical(.tabicl_model_digest(m), .tabicl_model_digest(mb))
})
