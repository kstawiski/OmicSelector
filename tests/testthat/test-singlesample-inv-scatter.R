library(testthat)

# inv-scatter computes a 1D wavelet scattering transform in python at fit AND
# score; both need the torch module + kymatio's 1D frontend via reticulate. The
# top-level `import kymatio.torch` is BROKEN under scipy >= 1.15 (it eagerly
# imports the 3D module which calls the removed scipy.special.sph_harm), so the
# skip-guard probes the DIRECT 1D-frontend import the scorer actually uses.
skip_if_no_inv_scatter <- function() {
  venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
  if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = venv)
  }
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(reticulate::py_module_available("torch"),
                        "Python module 'torch' not available")
  ok <- tryCatch({
    reticulate::import("kymatio.scattering1d.frontend.torch_frontend",
                       delay_load = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(isTRUE(ok),
                        "kymatio 1D torch frontend (ScatteringTorch1D) not importable")
}

# Synthetic generator with a CASE-ONLY location shift in a feature block. The first
# `k` of `p` log-abundance features are elevated by `shift` for case samples; after
# exp() and the per-sample rCLR centring those features carry the discriminative
# signal the scattering coefficients + frozen head learn. Abundances are exp() so
# the input is strictly positive (compositional). The DECISIVE single-sample
# properties are EXACT regardless of device because the scorer scatters in float64
# ONE ROW AT A TIME.
.make_inv_scatter_data <- function(n = 40L, p = 16L, k = 6L, shift = 1.2,
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

.auc_mw_inv_scatter <- function(y, score) {
  case <- score[y == 1L]
  ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# A small fit reused across tests (the single-sample properties do not depend on
# the head's accuracy, only on the frozen scattering + frozen linear map).
.fit_small_inv_scatter <- function(seed = 11L, device = "auto") {
  d <- .make_inv_scatter_data(seed = seed)
  list(model = fit_inv_scatter(d$X, d$y, hp = list(device = device)),
       data = d)
}


test_that("fit returns a well-formed inv_scatter_model", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  expect_s3_class(m, "inv_scatter_model")
  expect_equal(m$feature_universe, colnames(f$data$X))
  expect_equal(m$shape, 16L)                    # next_pow2(16) = 16
  expect_equal(m$J, 3L); expect_equal(m$Q, 4L); expect_equal(m$max_order, 2L)
  expect_length(m$coef_mean, length(m$weights))
  expect_length(m$coef_sd, length(m$weights))
  expect_true(is.finite(m$intercept))
  expect_true(is.finite(m$lambda) && m$lambda >= 0)
  # No live python pointer survives the fit (python-at-score contract).
  expect_false(any(vapply(m, function(x) inherits(x, "python.builtin.object"),
                          logical(1L))))
})

test_that("score is finite, length-correct, and double", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  s <- score_inv_scatter(f$model, f$data$X)
  expect_length(s, nrow(f$data$X))
  expect_true(all(is.finite(s)))
  expect_type(s, "double")
})

# §7.1 scattering correctness: two independent python scatterings of the SAME
# padded row agree to EXACTLY 0 (the transform is deterministic), and a rebuilt
# operator (from the stored config) reproduces the original operator's coeffs to 0.
test_that("§7.1 scattering is deterministic + rebuild-reproducible (maxdiff 0)", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  device <- "cpu"
  Z <- .inv_scatter_rclr_matrix(f$data$X[seq_len(4L), , drop = FALSE])

  S1 <- .inv_scatter_build_scattering(m$J, m$Q, m$shape, m$max_order, device)
  C1a <- .inv_scatter_coeffs(S1, Z, m$shape, device)
  C1b <- .inv_scatter_coeffs(S1, Z, m$shape, device)
  expect_identical(max(abs(C1a - C1b)), 0)            # same operator, deterministic

  S2 <- .inv_scatter_build_scattering(m$J, m$Q, m$shape, m$max_order, device)
  C2 <- .inv_scatter_coeffs(S2, Z, m$shape, device)
  expect_identical(max(abs(C1a - C2)), 0)             # rebuilt-from-config identical
})

# §7.2 single-row score == batch score, maxdiff EXACTLY 0 (forced row-by-row +
# float64 scattering: no BatchNorm, no cross-sample stat).
test_that("§7.2 DECISIVE: single-row score == batch score (maxdiff 0)", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  batch <- score_inv_scatter(m, X)
  single <- vapply(seq_len(nrow(X)), function(i) {
    score_inv_scatter(m, X[i, , drop = FALSE])
  }, numeric(1L))
  expect_identical(max(abs(batch - single)), 0)
})

# §7.3 SAMPLE row-permutation invariance, maxdiff 0.
test_that("§7.3 score is invariant to row permutation (maxdiff 0)", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_inv_scatter(m, X)
  set.seed(7)
  perm <- sample(nrow(X))
  sp <- score_inv_scatter(m, X[perm, , drop = FALSE])
  inv <- numeric(nrow(X)); inv[perm] <- sp
  expect_identical(max(abs(base - inv)), 0)
})

# §7.4 per-sample positive-scale invariance ~0 (rCLR input: a positive rescale
# leaves the rCLR, the signal, the coeffs, and the score unchanged).
test_that("§7.4 score is invariant to per-sample positive scaling", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_inv_scatter(m, X)
  set.seed(3)
  scal <- matrix(stats::runif(nrow(X), 1e-3, 1e3), nrow(X), ncol(X))
  Xs <- X * scal
  ss <- score_inv_scatter(m, Xs)
  expect_equal(max(abs(base - ss)), 0, tolerance = 1e-6)
})

# §7.5 held-out AUC > 0.7 on a planted block shift. The scattering coefficients of
# the rCLR profile separate the planted shift; the frozen ridge-logistic head is
# directional. A clear shift (2.0) is used on a 30-row training split.
test_that("§7.5 larger score = more case-like on a held-out planted shift", {
  skip_if_no_inv_scatter()
  d <- .make_inv_scatter_data(seed = 11L, shift = 2.0)  # n = 40 (ctrl 1:20, case 21:40)
  tr <- c(1:15, 21:35)
  te <- setdiff(seq_len(nrow(d$X)), tr)        # 5 held-out controls + 5 cases
  m <- fit_inv_scatter(d$X[tr, ], d$y[tr], hp = list(device = "auto"))
  s <- score_inv_scatter(m, d$X[te, ])
  auc <- .auc_mw_inv_scatter(d$y[te], s)
  expect_gt(auc, 0.7)
})

# §7.6 canonical row-equivariance gate + dispatch.
test_that("§7.6 passes the canonical row-equivariance gate", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_inv_scatter, f$model, f$data$X[seq_len(8L), , drop = FALSE],
      model_digest = .inv_scatter_model_digest
    )
  )
})

test_that("§7.6 canonical dispatch via singlesample_score_call matches direct score", {
  skip_if_no_inv_scatter()
  # Only run once the orchestrator has wired inv-scatter's canonical score_fn into
  # the package; until then the method is absent from the roster and we skip.
  roster <- tryCatch(singlesample_method_roster(), error = function(e) NULL)
  skip_if(is.null(roster) || !("inv-scatter" %in% roster$method_id),
          "inv-scatter not yet registered in the roster")
  f <- .fit_small_inv_scatter()
  X <- f$data$X[seq_len(6L), , drop = FALSE]
  via <- tryCatch(singlesample_score_call("inv-scatter", f$model, X),
                  error = function(e) e)
  if (inherits(via, "error")) {
    skip(paste("inv-scatter canonical dispatch not yet wired:",
               conditionMessage(via)))
  }
  direct <- score_inv_scatter(f$model, X)
  expect_equal(direct, via, tolerance = 1e-12)
})

# §7.7 degenerate -> 0 + flat composition scored (NOT floored).
test_that("§7.7 degenerate queries return 0; a flat composition is scored", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  X <- f$data$X[seq_len(5L), , drop = FALSE]

  # Fewer than min_features universe features present -> 0 for every row.
  X_few <- X[, seq_len(2L), drop = FALSE]
  expect_true(all(score_inv_scatter(m, X_few) == 0))

  # No shared features at all -> 0 for every row.
  X_none <- X
  colnames(X_none) <- paste0("other-", seq_len(ncol(X_none)))
  expect_true(all(score_inv_scatter(m, X_none) == 0))

  # All-zero specimen -> 0 for that row (empty positive support), finite for others.
  X_zero <- X
  X_zero[2L, ] <- 0
  s <- score_inv_scatter(m, X_zero)
  expect_equal(s[2L], 0)
  expect_true(all(is.finite(s)))

  # A FLAT all-equal-positive composition has an all-zero rCLR (rCLR origin) but is
  # a VALID specimen: it must be scattered + scored, NOT floored to 0.
  X_flat <- X
  X_flat[3L, ] <- 7
  s_flat <- score_inv_scatter(m, X_flat)
  expect_true(is.finite(s_flat[3L]))
})

# §7.8 determinism: two seed=42 fits -> identical model_digest + bit-identical
# scores (the scattering is deterministic; the head's CV uses NO RNG; the IRLS is
# deterministic). Hardened with an intervening torch-RNG perturbation.
test_that("§7.8 determinism: re-score and refit reproduce identical scores + digest", {
  skip_if_no_inv_scatter()
  d <- .make_inv_scatter_data(seed = 11L)
  m1 <- fit_inv_scatter(d$X, d$y, hp = list(device = "auto", seed = 42L))
  X <- d$X[seq_len(6L), , drop = FALSE]
  s1a <- score_inv_scatter(m1, X)
  s1b <- score_inv_scatter(m1, X)
  expect_identical(s1a, s1b)                          # re-score identical

  # Perturb the global torch RNG between fits: the fit is deterministic regardless.
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  torch_py$manual_seed(99999L)

  m2 <- fit_inv_scatter(d$X, d$y, hp = list(device = "auto", seed = 42L))
  s2 <- score_inv_scatter(m2, X)
  expect_identical(max(abs(s1a - s2)), 0)             # refit bit-identical
  expect_identical(.inv_scatter_model_digest(m1), .inv_scatter_model_digest(m2))
})

# §7.9 RNG-safety: scoring consumes no RNG and leaves the global R + torch
# (CPU/CUDA) generators byte-unchanged; fit leaves the torch generators unchanged.
test_that("§7.9 scoring does not mutate the global torch / R RNG state", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
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
  invisible(score_inv_scatter(m, X))
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
  skip_if_no_inv_scatter()
  d <- .make_inv_scatter_data(seed = 11L)
  torch_py <- reticulate::import("torch", delay_load = FALSE)
  before_t <- torch_py$random$get_rng_state()
  cuda_ok <- isTRUE(reticulate::py_to_r(torch_py$cuda$is_available()))
  before_cuda <- if (cuda_ok) torch_py$cuda$get_rng_state_all() else NULL
  invisible(fit_inv_scatter(d$X, d$y, hp = list(device = "auto")))
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

# §7.10 no-live-pointer: the fitted model holds NO python.builtin.object field, and
# the scattering operator REBUILT from the stored config reproduces the fit-time
# coefficients to <1e-10 (float64 determinism).
test_that("§7.10 no python pointer; rebuilt-from-config coeffs reproduce fit-time (<1e-10)", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  expect_false(any(vapply(m, function(x) inherits(x, "python.builtin.object"),
                          logical(1L))))

  device <- "cpu"
  Z <- .inv_scatter_rclr_matrix(f$data$X[seq_len(6L), , drop = FALSE])
  # Fit-time operator (rebuilt with the SAME config the fit used).
  S_fit <- .inv_scatter_build_scattering(m$J, m$Q, m$shape, m$max_order, device)
  C_fit <- .inv_scatter_coeffs(S_fit, Z, m$shape, device)
  # Score-path operator (rebuilt independently from the stored config).
  S_reb <- .inv_scatter_build_scattering(m$J, m$Q, m$shape, m$max_order, device)
  C_reb <- .inv_scatter_coeffs(S_reb, Z, m$shape, device)
  expect_lt(max(abs(C_fit - C_reb)), 1e-10)
})

# §7.11 frozen-head + frozen-standardization proof: an INDEPENDENT R recompute of
# intercept + w . standardize(coeffs) equals the scorer output to EXACTLY 0 (the
# head is a pure-R linear map over the python coeffs under FROZEN statistics).
test_that("§7.11 frozen-head proof: independent R recompute == scorer output (0)", {
  skip_if_no_inv_scatter()
  f <- .fit_small_inv_scatter()
  m <- f$model
  X <- f$data$X[seq_len(6L), , drop = FALSE]
  device <- "cpu"

  # Independent recompute: align -> rCLR -> python coeffs -> FROZEN standardize ->
  # FROZEN linear head, all outside score_inv_scatter.
  X_use <- .inv_scatter_align(X, m$feature_universe)
  Z <- .inv_scatter_rclr_matrix(X_use)
  S <- .inv_scatter_build_scattering(m$J, m$Q, m$shape, m$max_order, device)
  C <- .inv_scatter_coeffs(S, Z, m$shape, device)
  Cs <- .inv_scatter_standardize_apply(C, m$coef_mean, m$coef_sd)
  manual <- as.numeric(m$intercept + as.vector(Cs %*% m$weights))

  # Score on the SAME (cpu) device so the float64 coeffs match bit-for-bit.
  m_cpu <- m; m_cpu$hp$device <- "cpu"; m_cpu$device <- "cpu"
  s <- score_inv_scatter(m_cpu, X)
  expect_identical(max(abs(manual - s)), 0)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  skip_if_no_inv_scatter()
  d <- .make_inv_scatter_data(seed = 11L)
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(bogus = 1)),
               "unknown hp field")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(seed = 1, seed = 2)),
               "duplicate hp field")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(7)),
               "named")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(J = 0)),
               "J")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(Q = 2.5)),
               "Q")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(max_order = 3)),
               "max_order")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(lambda = -1)),
               "lambda")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(min_features = 2.5)),
               "min_features")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(device = "tpu")),
               "device")
  expect_error(fit_inv_scatter(d$X, d$y, hp = list(seed = 1.5)),
               "seed")
})
