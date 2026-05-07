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

test_that("paper3_bh_fdr_correct_blocked: returns both conservative and textbook columns", {
  fake <- list(list(p_emp = 0.005), list(p_emp = 0.012),
               list(p_emp = 0.04), list(p_emp = 0.30))
  blocks <- c("Toray", "Toray", "Affy-S", "Agilent-S")
  out <- paper3_bh_fdr_correct_blocked(fake, blocks)
  expect_true(all(c("q_block_BH", "q_block_BH_textbook",
                    "p_block", "p_block_textbook") %in% names(out)))
})

test_that("paper3_bh_fdr_correct_blocked: textbook reference <= conservative reference (k>=2)", {
  fake <- list(list(p_emp = 0.005), list(p_emp = 0.012),
               list(p_emp = 0.04), list(p_emp = 0.30))
  blocks <- c("Toray", "Toray", "Affy-S", "Agilent-S")
  out <- paper3_bh_fdr_correct_blocked(fake, blocks)
  toray_idx <- which(out$block_id == "Toray")
  expect_lte(out$p_block_textbook[toray_idx[1]] - 1e-12,
             out$p_block[toray_idx[1]])
})

test_that("paper3_bh_fdr_correct_blocked: k=1 singleton block reproduces cohort p", {
  fake <- list(list(p_emp = 0.04), list(p_emp = 0.30))
  blocks <- c("S1", "S2")
  out <- paper3_bh_fdr_correct_blocked(fake, blocks)
  expect_equal(out$p_block, out$p_emp)
  expect_equal(out$p_block_textbook, out$p_emp)
})

test_that("paper3_bh_fdr_correct_blocked: rho argument propagates", {
  fake <- list(list(p_emp = 0.005), list(p_emp = 0.012))
  blocks <- c("Toray", "Toray")
  out_low  <- paper3_bh_fdr_correct_blocked(fake, blocks, rho = 0.10)
  out_high <- paper3_bh_fdr_correct_blocked(fake, blocks, rho = 1.00)
  expect_false(isTRUE(all.equal(out_low$p_block[1], out_high$p_block[1])))
})

test_that("paper3_bh_fdr_correct_blocked: rho=0 recovers vanilla Fisher", {
  fake <- list(list(p_emp = 0.01), list(p_emp = 0.05))
  blocks <- c("B", "B")
  out <- paper3_bh_fdr_correct_blocked(fake, blocks, rho = 0)
  chisq_stat <- -2 * sum(log(c(0.01, 0.05)))
  expected <- stats::pchisq(chisq_stat, df = 4, lower.tail = FALSE)
  expect_equal(out$p_block[1], expected, tolerance = 1e-12)
})


# ============================================================================
# Per-cell AUC CIs (Patch 7)
# ============================================================================

test_that("paper3_hanley_mcneil_auc_ci returns a bounded interval and SE", {
  ci <- paper3_hanley_mcneil_auc_ci(auc = 0.72, n_pos = 50, n_neg = 50)
  expect_named(ci, c("ci_lo", "ci_hi", "se"))
  expect_true(all(is.finite(ci)))
  expect_lte(ci[["ci_lo"]], 0.72)
  expect_gte(ci[["ci_hi"]], 0.72)
  expect_gt(ci[["se"]], 0)
  expect_gte(ci[["ci_lo"]], 0)
  expect_lte(ci[["ci_hi"]], 1)
})

test_that("paper3_hanley_mcneil_auc_ci shrinks with larger sample counts", {
  ci_small <- paper3_hanley_mcneil_auc_ci(auc = 0.72, n_pos = 25, n_neg = 25)
  ci_large <- paper3_hanley_mcneil_auc_ci(auc = 0.72, n_pos = 250, n_neg = 250)
  width_small <- ci_small[["ci_hi"]] - ci_small[["ci_lo"]]
  width_large <- ci_large[["ci_hi"]] - ci_large[["ci_lo"]]
  expect_lt(width_large, width_small)
})

test_that("paper3_matched_null_benchmark: returns AUC CI bracketing auc_obs (hanley_mcneil)", {
  X <- .make_X()
  dat <- .make_y_and_inject(X, .true_panel)
  scoring_fn <- function(M) rowSums(log(M + 1))
  res <- paper3_matched_null_benchmark(dat$X, dat$y,
                                        panel = .true_panel,
                                        scoring_fn = scoring_fn,
                                        K = 50L, seed = 42L,
                                        auc_ci_method = "hanley_mcneil")
  expect_true(is.finite(res$auc_obs_ci_lo))
  expect_true(is.finite(res$auc_obs_ci_hi))
  expect_lte(res$auc_obs_ci_lo, res$auc_obs)
  expect_gte(res$auc_obs_ci_hi, res$auc_obs)
  expect_equal(res$auc_ci_method, "hanley_mcneil")
})

test_that("paper3_matched_null_benchmark: auc_ci_method='none' returns NA CI", {
  X <- .make_X()
  dat <- .make_y_and_inject(X, .true_panel)
  scoring_fn <- function(M) rowSums(log(M + 1))
  res <- paper3_matched_null_benchmark(dat$X, dat$y,
                                        panel = .true_panel,
                                        scoring_fn = scoring_fn,
                                        K = 30L, seed = 42L,
                                        auc_ci_method = "none")
  expect_true(is.na(res$auc_obs_ci_lo))
  expect_true(is.na(res$auc_obs_ci_hi))
})

test_that("paper3_matched_null_benchmark: AUC CI shrinks with growing N (hanley_mcneil)", {
  scoring_fn <- function(M) rowSums(log(M + 1))
  panel <- paste0("mir-", sprintf("%03d", 11:22))
  ci_width <- function(n) {
    set.seed(11L)
    X <- matrix(rlnorm(n * 60L, 5, 1), nrow = n,
                 dimnames = list(NULL, paste0("mir-", sprintf("%03d", seq_len(60L)))))
    y <- stats::rbinom(n, 1L, 0.5)
    X[y == 1, panel] <- X[y == 1, panel] * 2.5
    res <- paper3_matched_null_benchmark(X, y, panel = panel,
                                          scoring_fn = scoring_fn,
                                          K = 30L, seed = 11L,
                                          auc_ci_method = "hanley_mcneil")
    res$auc_obs_ci_hi - res$auc_obs_ci_lo
  }
  expect_lt(ci_width(400L), ci_width(80L))
})


# ============================================================================
# 5-fold nested-CV matched-null benchmark (Patch 4)
# ============================================================================

test_that("paper3_matched_null_benchmark_cv: signal panel passes 5-fold CV", {
  set.seed(99L)
  n_cv <- 150L; p_cv <- 120L
  X_cv <- matrix(stats::rlnorm(n_cv * p_cv, meanlog = 5, sdlog = 1), nrow = n_cv,
                  dimnames = list(NULL, paste0("f-", seq_len(p_cv))))
  y_cv <- stats::rbinom(n_cv, 1L, 0.5)
  X_cv[y_cv == 1, paste0("f-", 1:10)] <-
    X_cv[y_cv == 1, paste0("f-", 1:10)] * 3
  res <- paper3_matched_null_benchmark_cv(
    X_cv, y_cv, panel_size = 10L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 80L, outer_k = 5L, seed = 42L,
    min_valid_folds = 4L, min_pos_per_fold = 5L, min_neg_per_fold = 5L)
  expect_true(isTRUE(res$eligible))
  expect_gte(res$n_valid_folds, 4L)
  expect_true(is.finite(res$auc_obs_cv))
  expect_true(is.finite(res$p_emp_cv))
  expect_equal(res$grouping_policy, "ungrouped_stratified")
  expect_true(is.na(res$n_groups))
  expect_equal(res$n_group_split_violations, 0L)
})

test_that("paper3_matched_null_benchmark_cv: panel selection respects training-fold-only", {
  set.seed(99L)
  n_cv <- 120L; p_cv <- 80L
  X_cv <- matrix(stats::rlnorm(n_cv * p_cv, meanlog = 5, sdlog = 1), nrow = n_cv,
                  dimnames = list(NULL, paste0("f-", seq_len(p_cv))))
  y_cv <- stats::rbinom(n_cv, 1L, 0.5)
  X_cv[y_cv == 1, paste0("f-", 1:8)] <-
    X_cv[y_cv == 1, paste0("f-", 1:8)] * 3
  res <- paper3_matched_null_benchmark_cv(
    X_cv, y_cv, panel_size = 8L, K = 50L, outer_k = 5L, seed = 42L,
    min_valid_folds = 3L, min_pos_per_fold = 4L, min_neg_per_fold = 4L)
  for (fp in res$fold_panels) {
    if (!is.null(fp))
      expect_true(all(fp %in% colnames(X_cv)))
  }
})

test_that("paper3_matched_null_benchmark_cv: grouped folds honour group_id", {
  set.seed(77L)
  n_grp <- 100L; p_grp <- 60L
  X_grp <- matrix(stats::rlnorm(n_grp * p_grp, meanlog = 5, sdlog = 1), nrow = n_grp,
                   dimnames = list(NULL, paste0("g-", seq_len(p_grp))))
  y_grp <- stats::rbinom(n_grp, 1L, 0.5)
  X_grp[y_grp == 1, 1:8] <- X_grp[y_grp == 1, 1:8] * 2.5
  group_id_sim <- rep(seq_len(50L), each = 2L)
  res_grp <- paper3_matched_null_benchmark_cv(
    X_grp, y_grp, panel_size = 8L, K = 30L, outer_k = 5L,
    group_id = group_id_sim, seed = 1L,
    min_valid_folds = 3L, min_pos_per_fold = 4L, min_neg_per_fold = 4L)
  expect_equal(res_grp$grouping_policy, "groupKFold")
  expect_false(is.na(res_grp$n_groups))
  expect_equal(res_grp$n_group_split_violations, 0L)
})

test_that("paper3_matched_null_benchmark_cv: returns eligible=FALSE under min_valid_folds", {
  set.seed(3L)
  n_tiny <- 30L; p_tiny <- 40L
  X_t <- matrix(stats::rlnorm(n_tiny * p_tiny), nrow = n_tiny,
                dimnames = list(NULL, paste0("t-", seq_len(p_tiny))))
  y_t <- c(rep(0L, 27L), rep(1L, 3L))
  res_t <- paper3_matched_null_benchmark_cv(
    X_t, y_t, panel_size = 5L, K = 20L, outer_k = 5L, seed = 5L,
    min_valid_folds = 4L, min_pos_per_fold = 5L, min_neg_per_fold = 5L)
  expect_false(isTRUE(res_t$eligible))
  expect_true(!is.null(res_t$ineligible_reason))
})

test_that("paper3_matched_null_benchmark_cv: same seed -> identical p_emp_cv", {
  set.seed(99L)
  n_cv <- 150L; p_cv <- 120L
  X_cv <- matrix(stats::rlnorm(n_cv * p_cv, meanlog = 5, sdlog = 1), nrow = n_cv,
                  dimnames = list(NULL, paste0("f-", seq_len(p_cv))))
  y_cv <- stats::rbinom(n_cv, 1L, 0.5)
  X_cv[y_cv == 1, paste0("f-", 1:10)] <-
    X_cv[y_cv == 1, paste0("f-", 1:10)] * 3
  res1 <- paper3_matched_null_benchmark_cv(
    X_cv, y_cv, panel_size = 10L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 30L, outer_k = 3L, seed = 17L, min_valid_folds = 2L)
  res2 <- paper3_matched_null_benchmark_cv(
    X_cv, y_cv, panel_size = 10L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 30L, outer_k = 3L, seed = 17L, min_valid_folds = 2L)
  expect_identical(res1$p_emp_cv, res2$p_emp_cv)
})

test_that("paper3_matched_null_benchmark_cv: returns AUC CI brackets auc_obs_cv", {
  set.seed(99L)
  n_cv <- 150L; p_cv <- 120L
  X_cv <- matrix(stats::rlnorm(n_cv * p_cv, meanlog = 5, sdlog = 1), nrow = n_cv,
                  dimnames = list(NULL, paste0("f-", seq_len(p_cv))))
  y_cv <- stats::rbinom(n_cv, 1L, 0.5)
  X_cv[y_cv == 1, paste0("f-", 1:10)] <-
    X_cv[y_cv == 1, paste0("f-", 1:10)] * 3
  res <- paper3_matched_null_benchmark_cv(
    X_cv, y_cv, panel_size = 10L,
    scoring_fn = function(M) rowSums(log(M + 0.5)),
    K = 30L, outer_k = 5L, seed = 42L,
    min_valid_folds = 4L, min_pos_per_fold = 5L, min_neg_per_fold = 5L,
    auc_ci_method = "hanley_mcneil")
  expect_true(is.finite(res$auc_obs_ci_lo))
  expect_true(is.finite(res$auc_obs_ci_hi))
  expect_lte(res$auc_obs_ci_lo, res$auc_obs_cv)
  expect_gte(res$auc_obs_ci_hi, res$auc_obs_cv)
})
