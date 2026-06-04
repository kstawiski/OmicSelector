# Tests for Paper 3 within-sample compositional methods (Module A, P1).
# Source files: R/singlesample-within-sample.R.
# Acceptance criteria mirror the synthetic-data self-tests previously
# embedded in /umed-projekty/JAJNIKI/OmicSelector_paper/code/methods/.

test_that("ws_rclr_trimmed centers around zero on the centering subset", {
  set.seed(42L)
  x <- rlnorm(50L, meanlog = 5, sdlog = 1.2)
  names(x) <- c(paste0("hsa-miR-", sprintf("%03d", seq_len(45))),
                "hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                "hsa-miR-144-3p", "hsa-miR-223-3p")
  z <- ws_rclr_trimmed(x)
  expect_type(z, "double")
  expect_equal(length(z), length(x))
  expect_lt(abs(mean(z[attr(z, "centering_features")])), 1e-9)
  expect_gte(attr(z, "n_centering"), 8L)
  # Hemolysis markers must NOT be in the centering subset.
  expect_false("hsa-miR-451a" %in% attr(z, "centering_features"))
  expect_false("hsa-miR-486-5p" %in% attr(z, "centering_features"))
})

test_that("ws_rclr_trimmed processes matrices row-wise", {
  set.seed(43L)
  X <- matrix(abs(rnorm(20L * 50L, mean = 100, sd = 30)), nrow = 20L,
              dimnames = list(paste0("S", 1:20),
                              c(paste0("hsa-miR-", sprintf("%03d", seq_len(45))),
                                "hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                                "hsa-miR-144-3p", "hsa-miR-223-3p")))
  Z <- ws_rclr_trimmed(X)
  expect_equal(dim(Z), dim(X))
  expect_equal(rownames(Z), rownames(X))
  expect_equal(colnames(Z), colnames(X))
})

test_that("ws_rclr_trimmed refuses degenerate panel sizes", {
  expect_error(ws_rclr_trimmed(c(a = 1, b = 2)),
               "min_centering_size")
})

test_that("ws_balance_ilr returns one balance per partition entry", {
  set.seed(42L)
  feat_names <- c(paste0("hsa-let-7", letters[1:7], "-5p"),
                  paste0("hsa-miR-", c("17","18a","19a","19b","20a","92a"), "-3p"),
                  paste0("hsa-miR-200", letters[1:3], "-3p"),
                  "hsa-miR-141-3p", "hsa-miR-141-5p", "hsa-miR-429",
                  paste0("hsa-miR-37", 1:3, "-3p"),
                  "hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                  "hsa-miR-144-3p", "hsa-miR-223-3p", "hsa-miR-126-3p",
                  paste0("hsa-miR-x", sprintf("%03d", 1:30)))
  x <- rlnorm(length(feat_names), meanlog = 5, sdlog = 1.2)
  names(x) <- feat_names
  b <- ws_balance_ilr(x)
  expect_equal(length(b), length(ws_default_sbp()))
  expect_equal(names(b), names(ws_default_sbp()))
  expect_true(any(is.finite(b)))
})

test_that("ws_balance_ilr returns NA when panel does not intersect partition", {
  x <- rlnorm(10L, meanlog = 5)
  names(x) <- paste0("MIMAT00000", 1:10)
  # Disable the 80%-coverage gate to exercise the per-balance NA path.
  b <- ws_balance_ilr(x, min_balance_coverage = 0)
  # Most balances will be NA because Toray-style MIMAT IDs don't match the
  # biology-keyed default SBP. The dominant_vs_tail balance is always
  # computable on numeric panels, so at least it should be finite.
  expect_true(is.finite(b["dominant_vs_tail"]))
})

test_that("ws_balance_ilr enforces 80% balance-coverage threshold by default", {
  x <- rlnorm(10L, meanlog = 5)
  names(x) <- paste0("MIMAT00000", 1:10)
  b <- ws_balance_ilr(x)  # default min_balance_coverage = 0.8
  expect_true(isTRUE(attr(b, "coverage_failed")))
  expect_true(all(is.na(b)))
  expect_true(is.numeric(attr(b, "coverage")))
  expect_lt(attr(b, "coverage"), 0.8)
})

test_that("ws_balance_ilr exposes coverage attributes when gate passes", {
  set.seed(42L)
  feat_names <- c(paste0("hsa-let-7", letters[1:7], "-5p"),
                  paste0("hsa-miR-", c("17","18a","19a","19b","20a","92a"), "-3p"),
                  paste0("hsa-miR-200", letters[1:3], "-3p"),
                  "hsa-miR-141-3p", "hsa-miR-141-5p", "hsa-miR-429",
                  paste0("hsa-miR-37", 1:3, "-3p"),
                  "hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
                  "hsa-miR-144-3p", "hsa-miR-223-3p", "hsa-miR-126-3p",
                  paste0("hsa-miR-x", sprintf("%03d", 1:30)))
  x <- rlnorm(length(feat_names), meanlog = 5, sdlog = 1.2)
  names(x) <- feat_names
  b <- ws_balance_ilr(x)
  expect_false(isTRUE(attr(b, "coverage_failed")))
  expect_gte(attr(b, "coverage"), 0.8)
  expect_equal(attr(b, "min_balance_coverage"), 0.8)
})

test_that("ws_balance_ilr requires named features", {
  expect_error(ws_balance_ilr(rlnorm(10L)), "feature names")
})

test_that("ws_alr_pivot transforms with default pivot pool", {
  set.seed(42L)
  feat <- c(ws_default_pivot_pool(), paste0("hsa-miR-test-", seq_len(50L)))
  m <- matrix(abs(rnorm(20L * length(feat), mean = 100, sd = 30)),
              nrow = 20L, ncol = length(feat),
              dimnames = list(NULL, feat))
  out <- ws_alr_pivot(m)
  expect_equal(dim(out), dim(m))
  expect_true(all(is.finite(out)))
  expect_equal(length(attr(out, "pivot_features")), 6L)
  expect_false(isTRUE(attr(out, "fallback_used")))
})

test_that("ws_alr_pivot is fail-closed on missing pivots", {
  m <- matrix(abs(rnorm(20L * 30L, mean = 100, sd = 30)),
              nrow = 20L, ncol = 30L,
              dimnames = list(NULL, paste0("unknown-", seq_len(30L))))
  expect_error(ws_alr_pivot(m), "pivot features present")
  out <- ws_alr_pivot(m, allow_global_fallback = TRUE)
  expect_true(isTRUE(attr(out, "fallback_used")))
})

test_that("ws_alr_pivot rejects negative inputs", {
  pool <- ws_default_pivot_pool()
  ncol_m <- length(pool) + 1L
  m <- matrix(abs(rnorm(20L * ncol_m, mean = 100, sd = 30)),
              nrow = 20L, ncol = ncol_m,
              dimnames = list(NULL, c(pool, "extra")))
  m[1L, 1L] <- -1
  expect_error(ws_alr_pivot(m), "negative values")
})

test_that("ws_mad_logratio yields per-sample median ~ 0 and MAD ~ 1", {
  set.seed(42L)
  x <- matrix(abs(rnorm(20L * 50L, mean = 100, sd = 30)), nrow = 20L, ncol = 50L)
  out <- ws_mad_logratio(x)
  expect_equal(dim(out), dim(x))
  expect_true(all(is.finite(out)))
  expect_lt(abs(median(apply(out, 1L, median))), 1e-9)
  expect_lt(abs(median(apply(out, 1L, mad)) - 1), 1e-3)
})

test_that("ws_mad_logratio without MAD scaling preserves median centering only", {
  set.seed(43L)
  x <- matrix(abs(rnorm(10L * 30L, mean = 100, sd = 30)), nrow = 10L, ncol = 30L)
  out <- ws_mad_logratio(x, scale_by_mad = FALSE)
  expect_equal(attr(out, "scaling"), "none")
  expect_lt(abs(median(apply(out, 1L, median))), 1e-9)
})

test_that("ws_dominance_score flags pre-analytical failures", {
  set.seed(42L)
  feat <- paste0("miR-", seq_len(20L))
  x_normal <- matrix(abs(rnorm(10L * 20L, mean = 100, sd = 30)), nrow = 10L,
                      dimnames = list(NULL, feat))
  x_dom <- matrix(abs(rnorm(2L * 20L, mean = 100, sd = 30)), nrow = 2L,
                   dimnames = list(NULL, feat))
  x_dom[1L, 1L] <- 10000
  x_dom[2L, 1:2] <- c(8000, 6000)
  combined <- rbind(x_normal, x_dom)
  scores <- ws_dominance_score(combined, top_k = 2L)
  expect_equal(length(scores), 12L)
  expect_true(all(scores[1:10] < 0.60))
  flags <- ws_dominance_flag(combined, threshold = 0.60)
  expect_true(all(flags[11:12]))
})

test_that("fit_dominance_threshold is between training median and known-pathological scores", {
  set.seed(42L)
  feat <- paste0("miR-", seq_len(20L))
  x_train <- matrix(abs(rnorm(50L * 20L, mean = 100, sd = 30)), nrow = 50L,
                     dimnames = list(NULL, feat))
  cal <- fit_dominance_threshold(x_train, alpha = 0.05)
  expect_named(cal, c("threshold", "alpha", "top_k", "reference_n", "reference_quantiles"))
  expect_equal(cal$reference_n, 50L)
  expect_gt(cal$threshold, median(ws_dominance_score(x_train)))
})

test_that("ws_dominance_score rejects negative inputs", {
  set.seed(42L)
  x <- matrix(abs(rnorm(10L * 20L, mean = 100, sd = 30)), nrow = 10L)
  x[1L, 1L] <- -1
  expect_error(ws_dominance_score(x), "negative values")
})

test_that("ws_default_sbp returns expected named partition list", {
  sbp <- ws_default_sbp()
  expect_named(sbp, c("rbc_vs_rest", "platelet_vs_rest", "let7_a_vs_let7_g",
                      "let7_canonical_vs_b_d_f", "mir17_vs_mir92",
                      "mir200_vs_mir141", "mir371_373_vs_rest_germcell",
                      "dominant_vs_tail"))
  expect_error(ws_default_sbp("circulating_v999"), "Unknown SBP version")
})

test_that("ws_default_pivot_pool excludes miR-92a-3p (erythrocyte-derived)", {
  pool <- ws_default_pivot_pool()
  expect_length(pool, 6L)
  expect_false("hsa-miR-92a-3p" %in% pool)
  expect_true(all(grepl("^hsa-(miR|let-7[a-z]?)-", pool)))
})
