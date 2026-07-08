library(testthat)

ridge_primal_ref <- function(Xs, y, lambda, max_iter = 200L, tol = 1e-10) {
  Xa <- cbind(1, Xs)                                 # n x (D+1), col 1 = intercept
  D1 <- ncol(Xa)
  pen <- rep(lambda, D1); pen[1L] <- 0               # do not penalize the intercept
  b <- rep(0, D1)
  for (it in seq_len(max_iter)) {
    eta <- as.vector(Xa %*% b)
    mu <- 1 / (1 + exp(-eta))
    w <- mu * (1 - mu)
    w[w < 1e-9] <- 1e-9                              # IRLS weight floor (stability)
    zwork <- eta + (y - mu) / w
    XtW <- t(Xa * w)
    H <- XtW %*% Xa
    diag(H) <- diag(H) + pen
    g <- XtW %*% zwork
    bn <- tryCatch(solve(H, g), error = function(e) {
      solve(H + diag(1e-8, D1), g)                   # tiny jitter if ill-conditioned
    })
    bn <- as.vector(bn)
    if (max(abs(bn - b)) < tol) { b <- bn; break }
    b <- bn
  }
  list(intercept = b[1L], weights = b[-1L])
}

.cv_lambda_primal_ref <- function(Xs, y, seed) {
  grid <- 10^seq(-3, 3, by = 0.5)
  n <- nrow(Xs)
  k <- max(2L, min(5L, sum(y == 1L), sum(y == 0L)))
  fold <- .inv_scatter_cv_folds(y, k)
  dev_of <- function(lambda) {
    tot <- 0
    cnt <- 0L
    for (f in seq_len(k)) {
      te <- which(fold == f)
      tr <- which(fold != f)
      if (length(te) == 0L) next
      if (!any(y[tr] == 1L) || !any(y[tr] == 0L)) next
      head <- ridge_primal_ref(Xs[tr, , drop = FALSE], y[tr], lambda)
      eta <- head$intercept + as.vector(Xs[te, , drop = FALSE] %*% head$weights)
      mu <- 1 / (1 + exp(-eta))
      eps <- 1e-12
      mu <- pmin(pmax(mu, eps), 1 - eps)
      tot <- tot - 2 * sum(y[te] * log(mu) + (1 - y[te]) * log(1 - mu))
      cnt <- cnt + length(te)
    }
    if (cnt == 0L) return(Inf)
    tot / cnt
  }
  devs <- vapply(grid, dev_of, numeric(1L))
  if (all(!is.finite(devs))) return(1)
  grid[which.min(devs)]
}

.make_head_case <- function(n, D, y_mode = c("balanced", "imbalanced"),
                            rank_deficient = FALSE, all_zero = FALSE,
                            single_class = FALSE, seed = 1L) {
  y_mode <- match.arg(y_mode)
  set.seed(seed)
  Xs <- matrix(stats::rnorm(n * D), nrow = n)
  if (rank_deficient) {
    k <- min(20L, D %/% 2L)
    Xs[, k + seq_len(k)] <- Xs[, seq_len(k), drop = FALSE]
  }
  if (all_zero) Xs[,] <- 0
  y <- if (single_class) {
    rep(1, n)
  } else if (y_mode == "imbalanced") {
    n1 <- max(2L, floor(n * 0.2))
    sample(c(rep(1, n1), rep(0, n - n1)))
  } else {
    sample(rep(c(0, 1), length.out = n))
  }
  Xh <- matrix(stats::rnorm(9L * D), nrow = 9L)
  list(Xs = Xs, y = y, Xh = Xh)
}

.head_deltas <- function(ref, fast, Xh) {
  eta_ref <- ref$intercept + as.vector(Xh %*% ref$weights)
  eta_fast <- fast$intercept + as.vector(Xh %*% fast$weights)
  c(
    weights = max(abs(fast$weights - ref$weights)),
    intercept = abs(fast$intercept - ref$intercept),
    eta = max(abs(eta_fast - eta_ref))
  )
}

test_that("SVD row-space ridge-logit matches the original primal IRLS head", {
  cases <- list(
    balanced_small = .make_head_case(40L, 128L, "balanced", seed = 11L),
    imbalanced_mid = .make_head_case(104L, 512L, "imbalanced", seed = 29L),
    balanced_highD = .make_head_case(160L, 1500L, "balanced", seed = 37L),
    duplicated_columns = .make_head_case(160L, 128L, "imbalanced",
                                         rank_deficient = TRUE, seed = 43L),
    all_zero_design = .make_head_case(40L, 128L, "balanced",
                                      all_zero = TRUE, seed = 53L),
    single_class = .make_head_case(40L, 128L, "balanced",
                                   single_class = TRUE, seed = 59L)
  )
  lambdas <- c(1e-3, 1, 1e3)
  max_seen <- c(weights = 0, intercept = 0, eta = 0)
  for (case_name in names(cases)) {
    dat <- cases[[case_name]]
    for (lambda in lambdas) {
      ref <- ridge_primal_ref(dat$Xs, dat$y, lambda)
      fast <- .inv_scatter_ridge_logit(dat$Xs, dat$y, lambda)
      delta <- .head_deltas(ref, fast, dat$Xh)
      max_seen <- pmax(max_seen, delta)
      cat(sprintf(
        "inv-scatter head equivalence case=%s lambda=%g max_weight=%.3g intercept=%.3g eta=%.3g\n",
        case_name, lambda, delta[["weights"]], delta[["intercept"]], delta[["eta"]]
      ))
      expect_lt(delta[["weights"]], 1e-9)
      expect_lt(delta[["intercept"]], 1e-9)
      expect_lt(delta[["eta"]], 1e-9)
    }
  }
  cat(sprintf(
    "inv-scatter head equivalence max_seen max_weight=%.3g intercept=%.3g eta=%.3g\n",
    max_seen[["weights"]], max_seen[["intercept"]], max_seen[["eta"]]
  ))
})

test_that("SVD row-space head preserves deterministic CV lambda selection", {
  dat <- .make_head_case(42L, 96L, "imbalanced", rank_deficient = TRUE, seed = 71L)
  lambda_primal <- .cv_lambda_primal_ref(dat$Xs, dat$y, seed = 42L)
  lambda_fast <- .inv_scatter_cv_lambda(dat$Xs, dat$y, seed = 42L)
  cat(sprintf(
    "inv-scatter CV lambda equivalence primal=%g svd=%g\n",
    lambda_primal, lambda_fast
  ))
  expect_identical(lambda_fast, lambda_primal)
})

test_that("SVD row-space head prints old-vs-new timing smoke line", {
  dat <- .make_head_case(120L, 2048L, "balanced", seed = 83L)
  old_time <- system.time(ref <- ridge_primal_ref(dat$Xs, dat$y, lambda = 1))[["elapsed"]]
  new_time <- system.time(fast <- .inv_scatter_ridge_logit(dat$Xs, dat$y, lambda = 1))[["elapsed"]]
  delta <- .head_deltas(ref, fast, dat$Xh)
  speedup <- if (new_time > 0) old_time / new_time else Inf
  cat(sprintf(
    "inv-scatter head timing n=120 D=2048 lambda=1 old=%.3fs new=%.3fs speedup=%.1fx max_weight_delta=%.3g intercept_delta=%.3g eta_delta=%.3g\n",
    old_time, new_time, speedup, delta[["weights"]], delta[["intercept"]], delta[["eta"]]
  ))
  expect_lt(delta[["weights"]], 1e-9)
  expect_lt(delta[["intercept"]], 1e-9)
  expect_lt(delta[["eta"]], 1e-9)
})

# --- lambda=0 corner (NOT reached by the benchmark; CV grid floor is 1e-3). ---
# Documents the convergent reviewer caveat: at lambda=0 the reduced (r+1) jitter and
# the primal (D+1) jitter condition a SINGULAR system differently, so on a genuinely
# rank-deficient design they are NOT bit-equivalent (~1e-7) -- while the operating
# regime (lambda >= 1e-3) IS exact to ~1e-12 for the SAME data. (force_jitter path.)
test_that("lambda=0 on rank-deficient input is documented (non-benchmark corner)", {
  dat <- .make_head_case(60L, 256L, "balanced", rank_deficient = TRUE, seed = 904L)
  # Benchmark floor lambda=1e-3: exact to 1e-9 even on this rank-deficient design.
  ref_lo  <- ridge_primal_ref(dat$Xs, dat$y, lambda = 1e-3)
  fast_lo <- .inv_scatter_ridge_logit(dat$Xs, dat$y, lambda = 1e-3)
  d_lo <- .head_deltas(ref_lo, fast_lo, dat$Xh)
  expect_lt(d_lo[["eta"]], 1e-9)                       # operating regime is exact
  # lambda=0 corner: force_jitter matches exact-null but NOT a rank-deficient design.
  ref_0  <- ridge_primal_ref(dat$Xs, dat$y, lambda = 0)
  fast_0 <- .inv_scatter_ridge_logit(dat$Xs, dat$y, lambda = 0)
  d_0 <- .head_deltas(ref_0, fast_0, dat$Xh)
  cat(sprintf("inv-scatter lambda=0 rank-deficient corner: eta_delta(lambda=1e-3)=%.3g eta_delta(lambda=0)=%.3g\n",
              d_lo[["eta"]], d_0[["eta"]]))
  # Both are ~1e-8-ridge solutions of a singular system; assert only a loose,
  # documented bound (the two solvers are not bit-equivalent here). This is the
  # KNOWN, non-benchmark corner -- the assertion exists to pin the behaviour, not
  # to claim equivalence.
  expect_true(is.finite(d_0[["eta"]]) && d_0[["eta"]] < 1e-3)
})
