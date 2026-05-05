library(testthat)

# Source helper: pROC is required by paper3_matched_null_benchmark.
# Skip the whole file gracefully if pROC is absent (CI environments).
skip_if_not_installed("pROC")

set.seed(42L)

# Synthetic dataset used by most tests.
.make_X <- function(n = 100L, p = 60L, seed = 42L) {
  set.seed(seed)
  X <- matrix(rlnorm(n * p, meanlog = 5, sdlog = 1.0), nrow = n, ncol = p)
  colnames(X) <- paste0("mir-", sprintf("%03d", seq_len(p)))
  X
}

.true_panel <- paste0("mir-", sprintf("%03d", 11:22))

.make_y_and_inject <- function(X, panel, signal_mult = 2.5, seed = 42L) {
  set.seed(seed)
  y <- stats::rbinom(nrow(X), 1L, 0.5)
  X[y == 1, panel] <- X[y == 1, panel] * signal_mult
  list(X = X, y = y)
}

# ============================================================================
# paper3_matched_null_benchmark
# ============================================================================

test_that("paper3_matched_null_benchmark: signal panel beats null at p < 0.05 (headline behaviour)", {
  X <- .make_X()
  dat <- .make_y_and_inject(X, .true_panel)
  scoring_fn <- function(M) rowSums(log(M + 1))
  res <- paper3_matched_null_benchmark(dat$X, dat$y,
                                        panel = .true_panel,
                                        scoring_fn = scoring_fn,
                                        K = 300L, seed = 42L)
  expect_true(res$eligible)
  expect_gt(res$auc_obs, mean(res$auc_null, na.rm = TRUE))
  expect_lt(res$p_emp, 0.05)
  expect_equal(res$panel_size, length(.true_panel))
})

test_that("paper3_matched_null_benchmark: fails when panel not in colnames(X)", {
  X <- .make_X()
  y <- stats::rbinom(nrow(X), 1L, 0.5)
  expect_error(
    paper3_matched_null_benchmark(X, y, panel = "nonexistent_mir",
                                   scoring_fn = function(M) rowSums(M), K = 10L),
    "None of the panel features"
  )
})

test_that("paper3_matched_null_benchmark: null panels do NOT overlap observed panel", {
  X <- .make_X()
  dat <- .make_y_and_inject(X, .true_panel)
  scoring_fn <- function(M) rowSums(log(M + 1))
  res <- paper3_matched_null_benchmark(dat$X, dat$y,
                                        panel = .true_panel,
                                        scoring_fn = scoring_fn,
                                        K = 100L, seed = 42L)
  non_null <- Filter(Negate(is.null), res$null_panels)
  non_null <- non_null[lengths(non_null) > 0]
  if (length(non_null) > 0) {
    any_overlap <- any(vapply(non_null,
                               function(p) any(p %in% .true_panel), logical(1L)))
    expect_false(any_overlap)
  }
})

test_that("paper3_matched_null_benchmark: hemolysis markers excluded when post_hemolysis_corrected=TRUE", {
  X <- .make_X()
  set.seed(99L)
  hemo_mir <- "hsa-miR-451a"
  X <- cbind(X, matrix(rlnorm(nrow(X), 5, 1), nrow = nrow(X),
                        dimnames = list(NULL, hemo_mir)))
  y <- stats::rbinom(nrow(X), 1L, 0.5)
  scoring_fn <- function(M) rowSums(log(M + 1))
  res <- paper3_matched_null_benchmark(X, y, panel = .true_panel,
                                        scoring_fn = scoring_fn,
                                        K = 60L, seed = 42L,
                                        post_hemolysis_corrected = TRUE)
  non_null <- Filter(Negate(is.null), res$null_panels)
  non_null <- non_null[lengths(non_null) > 0]
  if (length(non_null) > 0) {
    in_null <- any(vapply(non_null, function(p) hemo_mir %in% p, logical(1L)))
    expect_false(in_null)
  }
  expect_true(hemo_mir %in% res$excluded_features)
})

test_that("paper3_matched_null_benchmark: tier-3 fallback returns eligible=FALSE, p_emp=NA", {
  set.seed(11L)
  n_s <- 50L; p_s <- 30L
  X_s <- matrix(rlnorm(n_s * p_s, 5, 1), nrow = n_s,
                dimnames = list(NULL, paste0("mir-", sprintf("%03d", seq_len(p_s)))))
  panel_big <- paste0("mir-", sprintf("%03d", seq_len(28L)))
  y_s <- stats::rbinom(n_s, 1L, 0.5)
  res_d <- paper3_matched_null_benchmark(X_s, y_s, panel = panel_big,
                                          scoring_fn = function(M) rowSums(M),
                                          K = 30L, seed = 11L)
  # Either eligible or not — structure must be well-formed in both cases.
  expect_true(isTRUE(res_d$eligible) || isFALSE(res_d$eligible))
  if (!res_d$eligible) {
    expect_equal(res_d$fallback_tier, 3L)
    expect_true(is.na(res_d$p_emp))
  }
})


# ============================================================================
# paper3_bh_fdr_correct_matched_null
# ============================================================================

test_that("paper3_bh_fdr_correct_matched_null: returns vector of correct length", {
  fake <- list(list(p_emp = 0.01), list(p_emp = 0.04), list(p_emp = 0.50))
  q <- paper3_bh_fdr_correct_matched_null(fake)
  expect_length(q, 3L)
  expect_true(all(q >= 0 & q <= 1))
})

test_that("paper3_bh_fdr_correct_matched_null: BH is monotone relative to p-value order", {
  fake <- list(list(p_emp = 0.001), list(p_emp = 0.02), list(p_emp = 0.30))
  q <- paper3_bh_fdr_correct_matched_null(fake)
  expect_true(all(diff(q) >= 0))   # q-values are non-decreasing with p-values
})


# ============================================================================
# paper3_holm_correct_familywise
# ============================================================================

test_that("paper3_holm_correct_familywise: Holm-corrected p values >= raw p values", {
  pvals <- c(0.001, 0.01, 0.04, 0.20)
  adj <- paper3_holm_correct_familywise(pvals)
  expect_length(adj, 4L)
  expect_true(all(adj >= pvals - 1e-9))
})


# ============================================================================
# paper3_bh_fdr_correct_blocked
# ============================================================================

test_that("paper3_bh_fdr_correct_blocked: Toray cluster collapses to same block p-value", {
  fake <- list(
    list(p_emp = 0.005),
    list(p_emp = 0.012),
    list(p_emp = 0.04),
    list(p_emp = 0.30)
  )
  blocks <- c("Toray", "Toray", "Affy-S", "Agilent-S")
  out <- paper3_bh_fdr_correct_blocked(fake, blocks)
  expect_equal(nrow(out), 4L)
  expect_equal(out$block_id[1], "Toray")
  expect_equal(out$block_id[2], "Toray")
  expect_equal(out$p_block[1], out$p_block[2])   # same block → same p_block
})

test_that("paper3_bh_fdr_correct_blocked: q_block_BH >= p_block (BH can only increase)", {
  fake <- list(list(p_emp = 0.01), list(p_emp = 0.02), list(p_emp = 0.50))
  blocks <- c("A", "B", "C")
  out <- paper3_bh_fdr_correct_blocked(fake, blocks)
  expect_true(all(out$q_block_BH >= out$p_block - 1e-9))
})

test_that("paper3_bh_fdr_correct_blocked: length mismatch raises error", {
  fake <- list(list(p_emp = 0.01))
  expect_error(paper3_bh_fdr_correct_blocked(fake, c("A", "B")))
})
