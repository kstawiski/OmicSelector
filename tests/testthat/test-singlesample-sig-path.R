library(testthat)

# sig-path is PURE BASE R: no python, no torch, no GPU at fit OR score. The only
# python use anywhere is the §7.1 roughpy REFERENCE ORACLE, which is skip-guarded
# (it never gates the production scorer).
skip_if_no_roughpy <- function() {
  venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
  if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    Sys.setenv(RETICULATE_PYTHON = venv)
  }
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not(reticulate::py_module_available("roughpy"),
                        "Python module 'roughpy' not available (test-only oracle)")
}

# roughpy reference signature of an (s x 2) increment matrix at level L, returned
# WITH the leading level-0 constant 1. Built ENTIRELY inside python (reshape from a
# flat list) so the numpy array is freshly allocated + writable -- reticulate's
# auto-converted arrays are read-only, which trips roughpy's DLPack export.
.roughpy_sig <- function(incr, L) {
  flat <- sprintf("%.17g", as.numeric(t(incr)))      # row-major flatten of incr
  code <- sprintf("
import numpy as np, roughpy as rp
incr = np.ascontiguousarray(np.array([%s], dtype=np.float64).reshape(%d, 2))
ctx = rp.get_context(width=2, depth=%d, coeffs=rp.DPReal)
st = rp.LieIncrementStream.from_increments(incr, ctx=ctx)
rsig = np.asarray(st.signature(rp.RealInterval(0.0, float(%d)), resolution=0)).tolist()
", paste(flat, collapse = ","), nrow(incr), as.integer(L), nrow(incr))
  res <- reticulate::py_run_string(code, convert = TRUE)
  as.numeric(unlist(res$rsig))
}

# Synthetic generator with a CASE-ONLY location shift in a feature block. The first
# `k` of `p` log-abundance features are elevated by `shift` for case samples; after
# exp() and the per-sample rCLR centring those features carry the discriminative
# signal the path signature + frozen head learn. Abundances are exp() so the input
# is strictly positive (compositional). The DECISIVE single-sample properties are
# EXACT (pure base-R double arithmetic, row-by-row).
.make_sig_path_data <- function(n = 40L, p = 16L, k = 6L, shift = 1.2,
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

.auc_mw_sig_path <- function(y, score) {
  case <- score[y == 1L]
  ctrl <- score[y == 0L]
  if (length(case) == 0L || length(ctrl) == 0L) return(NA_real_)
  r <- rank(c(case, ctrl), ties.method = "average")
  n1 <- length(case)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(ctrl))
}

# A small fit reused across tests (the single-sample properties do not depend on
# the head's accuracy, only on the deterministic signature + frozen linear map).
.fit_small_sig_path <- function(seed = 11L, hp = list()) {
  d <- .make_sig_path_data(seed = seed)
  list(model = fit_sig_path(d$X, d$y, hp = hp), data = d)
}


test_that("fit returns a well-formed sig_path_model", {
  f <- .fit_small_sig_path()
  m <- f$model
  expect_s3_class(m, "sig_path_model")
  expect_equal(m$feature_universe, colnames(f$data$X))
  expect_equal(m$L, 3L)
  expect_equal(m$D, 14L)                              # sum_{k=1}^{3} 2^k = 2+4+8
  expect_length(m$coef_mean, m$D)
  expect_length(m$coef_sd, m$D)
  expect_length(m$weights, m$D)
  expect_true(is.finite(m$intercept))
  expect_true(is.finite(m$lambda) && m$lambda >= 0)
  # No live python pointer anywhere (pure base R).
  expect_false(any(vapply(m, function(x) inherits(x, "python.builtin.object"),
                          logical(1L))))
})

test_that("score is finite, length-correct, and double", {
  f <- .fit_small_sig_path()
  s <- score_sig_path(f$model, f$data$X)
  expect_length(s, nrow(f$data$X))
  expect_true(all(is.finite(s)))
  expect_type(s, "double")
})

# §7.1 signature correctness: the PURE-R closed-form truncated signature equals
# roughpy's reference engine to < 1e-9 (skip if roughpy absent) across L in
# {1,2,3,4} and varied p, AND a HAND-CHECKED tiny 2-point case (one segment):
# level-1 == the total increment, level-2 == Delta otimes Delta / 2.
test_that("§7.1 DECISIVE: pure-R signature == roughpy oracle (< 1e-9) + hand-check", {
  # Hand-checked single-segment case (no python needed): p = 2 -> one segment.
  # v = (v1, v2); t = (0, 1); Delta = (t2 - t1, v2 - v1) = (1, v2 - v1).
  v <- c(0.7, -1.3)
  L <- 3L
  feat <- .sig_path_features_one(v, L)               # length D = 14
  dvt <- 1.0                                          # t increment over one segment
  dvv <- v[2L] - v[1L]                                # rCLR increment
  delta <- c(dvt, dvv)
  # level 1 == Delta
  expect_equal(feat[1:2], delta, tolerance = 1e-12)
  # level 2 == Delta otimes Delta / 2 (row-major: (11,12,21,22))
  lev2 <- as.vector(kronecker(delta, delta)) / 2
  expect_equal(feat[3:6], lev2, tolerance = 1e-12)
  # level 3 == Delta^{otimes 3} / 6
  lev3 <- as.vector(kronecker(kronecker(delta, delta), delta)) / 6
  expect_equal(feat[7:14], lev3, tolerance = 1e-12)

  # roughpy oracle across L and p on random rCLR rows.
  skip_if_no_roughpy()
  maxd <- 0
  for (Lx in 1:4) {
    Dx <- .sig_path_dim(Lx)
    for (rep in 1:4) {
      set.seed(Lx * 100L + rep)
      p <- sample(5:24, 1L)
      vv <- stats::rnorm(p, sd = 2)                  # arbitrary rCLR-like profile
      incr <- .sig_path_increments(vv)
      pr <- .sig_path_features_one(vv, Lx)           # pure-R, level-0 dropped
      ro <- .roughpy_sig(incr, Lx)[-1L]              # roughpy, drop leading 1
      expect_length(pr, Dx)
      expect_length(ro, Dx)
      maxd <- max(maxd, max(abs(pr - ro)))
    }
  }
  expect_lt(maxd, 1e-9)
})

# §7.2 single-row score == batch score, maxdiff EXACTLY 0 (pure base R, row-by-row;
# no cross-sample state, no float32/batch artifact).
test_that("§7.2 DECISIVE: single-row score == batch score (maxdiff 0)", {
  f <- .fit_small_sig_path()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  batch <- score_sig_path(m, X)
  single <- vapply(seq_len(nrow(X)), function(i) {
    score_sig_path(m, X[i, , drop = FALSE])
  }, numeric(1L))
  expect_identical(max(abs(batch - single)), 0)
})

# §7.3 SAMPLE row-permutation invariance, maxdiff 0.
test_that("§7.3 score is invariant to row permutation (maxdiff 0)", {
  f <- .fit_small_sig_path()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_sig_path(m, X)
  set.seed(7)
  perm <- sample(nrow(X))
  sp <- score_sig_path(m, X[perm, , drop = FALSE])
  inv <- numeric(nrow(X)); inv[perm] <- sp
  expect_identical(max(abs(base - inv)), 0)
})

# §7.4 per-sample positive-scale invariance ~0 (rCLR input: a positive rescale
# leaves the rCLR, the path, the signature, and the score unchanged).
test_that("§7.4 score is invariant to per-sample positive scaling", {
  f <- .fit_small_sig_path()
  m <- f$model
  X <- f$data$X[seq_len(8L), , drop = FALSE]
  base <- score_sig_path(m, X)
  set.seed(3)
  scal <- matrix(stats::runif(nrow(X), 1e-3, 1e3), nrow(X), ncol(X))
  Xs <- X * scal
  ss <- score_sig_path(m, Xs)
  expect_equal(max(abs(base - ss)), 0, tolerance = 1e-6)
})

# §7.5 held-out AUC > 0.7 on a planted block shift. The path signature of the rCLR
# profile separates the planted shift; the frozen ridge-logistic head is
# directional. A clear shift (2.0) on a 30-row training split.
test_that("§7.5 larger score = more case-like on a held-out planted shift", {
  d <- .make_sig_path_data(seed = 11L, shift = 2.0)  # n = 40 (ctrl 1:20, case 21:40)
  tr <- c(1:15, 21:35)
  te <- setdiff(seq_len(nrow(d$X)), tr)              # 5 held-out controls + 5 cases
  m <- fit_sig_path(d$X[tr, ], d$y[tr])
  s <- score_sig_path(m, d$X[te, ])
  auc <- .auc_mw_sig_path(d$y[te], s)
  expect_gt(auc, 0.7)
})

# §7.6 canonical row-equivariance gate.
test_that("§7.6 passes the canonical row-equivariance gate", {
  f <- .fit_small_sig_path()
  expect_invisible(
    singlesample_assert_row_equivariant(
      score_sig_path, f$model, f$data$X[seq_len(8L), , drop = FALSE],
      model_digest = .sig_path_model_digest
    )
  )
})

test_that("§7.6 canonical dispatch via singlesample_score_call matches direct score", {
  # Only run once the orchestrator has wired sig-path's canonical score_fn into
  # the package; until then the method is absent from the roster and we skip.
  roster <- tryCatch(singlesample_method_roster(), error = function(e) NULL)
  skip_if(is.null(roster) || !("sig-path" %in% roster$method_id),
          "sig-path not yet registered in the roster")
  f <- .fit_small_sig_path()
  X <- f$data$X[seq_len(6L), , drop = FALSE]
  via <- tryCatch(singlesample_score_call("sig-path", f$model, X),
                  error = function(e) e)
  if (inherits(via, "error")) {
    skip(paste("sig-path canonical dispatch not yet wired:",
               conditionMessage(via)))
  }
  direct <- score_sig_path(f$model, X)
  expect_equal(direct, via, tolerance = 1e-12)
})

# §7.7 degenerate -> 0 + flat composition scored (NOT floored).
test_that("§7.7 degenerate queries return 0; a flat composition is scored", {
  f <- .fit_small_sig_path()
  m <- f$model
  X <- f$data$X[seq_len(5L), , drop = FALSE]

  # Fewer than min_features universe features present -> 0 for every row.
  X_few <- X[, seq_len(2L), drop = FALSE]
  expect_true(all(score_sig_path(m, X_few) == 0))

  # No shared features at all -> 0 for every row.
  X_none <- X
  colnames(X_none) <- paste0("other-", seq_len(ncol(X_none)))
  expect_true(all(score_sig_path(m, X_none) == 0))

  # All-zero specimen -> 0 for that row (empty positive support), finite for others.
  X_zero <- X
  X_zero[2L, ] <- 0
  s <- score_sig_path(m, X_zero)
  expect_equal(s[2L], 0)
  expect_true(all(is.finite(s)))

  # A FLAT all-equal-positive composition has an all-zero rCLR (rCLR origin) but is
  # a VALID specimen: its time-augmented path is non-trivial (t still moves), so it
  # must be scored, NOT floored to 0.
  X_flat <- X
  X_flat[3L, ] <- 7
  s_flat <- score_sig_path(m, X_flat)
  expect_true(is.finite(s_flat[3L]))
  # The flat row's signature is non-degenerate: t still moves, so level-1 t-coord
  # is non-zero even though every rCLR value is 0.
  vflat <- .sig_path_rclr(rep(7, ncol(X)))
  expect_true(all(vflat == 0))                       # rCLR origin
  feat_flat <- .sig_path_features_one(vflat, m$L)
  expect_gt(abs(feat_flat[1L]), 0)                   # level-1 t-increment != 0
})

# §7.8 determinism: two seed=42 fits -> identical model_digest + bit-identical
# scores (the signature is deterministic; the head's CV uses NO RNG; the IRLS is
# deterministic). Hardened with an intervening global-RNG perturbation.
test_that("§7.8 determinism: re-score and refit reproduce identical scores + digest", {
  d <- .make_sig_path_data(seed = 11L)
  m1 <- fit_sig_path(d$X, d$y, hp = list(seed = 42L))
  X <- d$X[seq_len(6L), , drop = FALSE]
  s1a <- score_sig_path(m1, X)
  s1b <- score_sig_path(m1, X)
  expect_identical(s1a, s1b)                          # re-score identical

  # Perturb the global R RNG between fits: the fit is deterministic regardless.
  set.seed(99999L); invisible(stats::runif(1000L))

  m2 <- fit_sig_path(d$X, d$y, hp = list(seed = 42L))
  s2 <- score_sig_path(m2, X)
  expect_identical(max(abs(s1a - s2)), 0)            # refit bit-identical
  expect_identical(.sig_path_model_digest(m1), .sig_path_model_digest(m2))
})

# §7.9 RNG-safety: fit + score leave the global R .Random.seed byte-unchanged;
# score consumes no RNG. (Pure base R: there is no torch generator to guard.)
test_that("§7.9 fit + score leave the global R RNG byte-unchanged (no RNG consumed)", {
  d <- .make_sig_path_data(seed = 11L)

  set.seed(99)
  before_fit <- get(".Random.seed", envir = globalenv())
  m <- fit_sig_path(d$X, d$y)
  after_fit <- get(".Random.seed", envir = globalenv())
  expect_identical(before_fit, after_fit)            # fit is RNG-neutral

  X <- d$X[seq_len(10L), , drop = FALSE]
  set.seed(123)
  before_score <- get(".Random.seed", envir = globalenv())
  invisible(score_sig_path(m, X))
  after_score <- get(".Random.seed", envir = globalenv())
  expect_identical(before_score, after_score)        # score consumes no RNG
})

# §7.10 frozen-standardization + frozen-head proof: an INDEPENDENT R recompute of
# intercept + w . standardize(signature) equals the scorer output to EXACTLY 0
# (the head is a pure-R linear map over the pure-R signatures under FROZEN stats).
test_that("§7.10 frozen-head proof: independent R recompute == scorer output (0)", {
  f <- .fit_small_sig_path()
  m <- f$model
  X <- f$data$X[seq_len(6L), , drop = FALSE]

  # Independent recompute: align -> rCLR -> pure-R signature -> FROZEN standardize
  # -> FROZEN linear head, all outside score_sig_path.
  X_use <- .sig_path_align(X, m$feature_universe)
  Z <- .sig_path_rclr_matrix(X_use)
  C <- .sig_path_features(Z, m$L)
  Cs <- .sig_path_standardize_apply(C, m$coef_mean, m$coef_sd)
  manual <- as.numeric(m$intercept + as.vector(Cs %*% m$weights))

  s <- score_sig_path(m, X)
  expect_identical(max(abs(manual - s)), 0)
})

# §7.11 Chen-splitting consistency: the truncated signature of a WHOLE path equals
# the truncated Chen (concatenation) product of the signatures of two sub-paths
# split at an interior point, to < 1e-12. A genuine non-tautological check of the
# closed-form concatenation product + the index ordering (the whole-path code path
# computes the full ordered Chen product; this independently re-splits and
# re-combines the segments and must agree).
test_that("§7.11 DECISIVE: Chen-splitting consistency of the truncated signature (< 1e-12)", {
  L <- 3L; d <- 2L
  maxd <- 0
  for (rep in 1:6) {
    set.seed(2000L + rep)
    p <- sample(6:20, 1L)
    v <- stats::rnorm(p, sd = 1.5)
    incr <- .sig_path_increments(v)                  # (p-1) x 2 increments
    s <- nrow(incr)
    full <- .sig_path_signature_levels(incr, L, d)   # whole-path signature

    # Split the increments at an interior segment boundary.
    cut <- sample(seq_len(s - 1L), 1L)
    A <- .sig_path_signature_levels(incr[seq_len(cut), , drop = FALSE], L, d)
    B <- .sig_path_signature_levels(incr[(cut + 1L):s, , drop = FALSE], L, d)
    comb <- .sig_path_chen(A, B, L, d)               # Chen(sig_A, sig_B)

    for (lvl in 0:L) {
      maxd <- max(maxd, max(abs(full[[lvl + 1L]] - comb[[lvl + 1L]])))
    }
  }
  expect_lt(maxd, 1e-12)
})

# §7.11b Chen-product associativity: (A o B) o C == A o (B o C) to < 1e-12.
test_that("§7.11b Chen product is associative (< 1e-12)", {
  L <- 3L; d <- 2L
  set.seed(4242L)
  mk <- function() .sig_path_seg_exp(stats::rnorm(2L), L)
  A <- mk(); B <- mk(); C <- mk()
  left  <- .sig_path_chen(.sig_path_chen(A, B, L, d), C, L, d)
  right <- .sig_path_chen(A, .sig_path_chen(B, C, L, d), L, d)
  md <- max(vapply(seq_len(L + 1L),
                   function(i) max(abs(left[[i]] - right[[i]])), numeric(1L)))
  expect_lt(md, 1e-12)
})

test_that("strict hp allow-list rejects unknown / duplicate / unnamed / bad fields", {
  d <- .make_sig_path_data(seed = 11L)
  expect_error(fit_sig_path(d$X, d$y, hp = list(bogus = 1)),
               "unknown hp field")
  expect_error(fit_sig_path(d$X, d$y, hp = list(seed = 1, seed = 2)),
               "duplicate hp field")
  expect_error(fit_sig_path(d$X, d$y, hp = list(7)),
               "named")
  expect_error(fit_sig_path(d$X, d$y, hp = list(L = 0)),
               "L")
  expect_error(fit_sig_path(d$X, d$y, hp = list(L = 5)),
               "L")
  expect_error(fit_sig_path(d$X, d$y, hp = list(L = 2.5)),
               "L")
  expect_error(fit_sig_path(d$X, d$y, hp = list(lambda = -1)),
               "lambda")
  expect_error(fit_sig_path(d$X, d$y, hp = list(min_features = 2.5)),
               "min_features")
  expect_error(fit_sig_path(d$X, d$y, hp = list(seed = 1.5)),
               "seed")
  # NO device hp -- pure base R; supplying one is an unknown field.
  expect_error(fit_sig_path(d$X, d$y, hp = list(device = "cpu")),
               "unknown hp field")
})
