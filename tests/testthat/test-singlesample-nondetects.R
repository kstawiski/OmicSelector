library(testthat)

test_that("qpcr_nondetect_lod_fallback fills NAs with Ct=40", {
  m <- matrix(c(20, NA, 25, 30, NA, NA), nrow = 2,
              dimnames = list(c("miR-1", "miR-2"), c("S1", "S2", "S3")))
  out <- qpcr_nondetect_lod_fallback(m, log_handler = function(...) invisible())
  expect_equal(sum(is.na(out)), 0L)
  expect_equal(sum(out == 40), 3L)
  expect_equal(attr(out, "nondetect_method"), "lod_fallback")
  expect_equal(attr(out, "n_imputed"), 3L)
})

test_that("qpcr_nondetect_lod_fallback no-op when no NAs", {
  m <- matrix(c(20, 25, 30, 35), nrow = 2,
              dimnames = list(c("miR-1", "miR-2"), c("S1", "S2")))
  out <- qpcr_nondetect_lod_fallback(m, log_handler = function(...) invisible())
  expect_equal(unclass(out), m, ignore_attr = TRUE)
  expect_equal(attr(out, "n_imputed"), 0L)
  expect_equal(attr(out, "nondetect_method"), "lod_fallback")
})

test_that("qpcr_nondetect_lod_fallback honours custom lod_ct", {
  m <- matrix(c(20, NA, 25, NA), nrow = 2,
              dimnames = list(c("miR-1", "miR-2"), c("S1", "S2")))
  out <- qpcr_nondetect_lod_fallback(m, lod_ct = 38,
                                      log_handler = function(...) invisible())
  expect_equal(sum(out == 38), 2L)
  expect_equal(attr(out, "n_imputed"), 2L)
})

test_that("qpcr_nondetect_impute falls back to LOD when nondetects absent", {
  skip_if(requireNamespace("nondetects", quietly = TRUE) &&
            requireNamespace("HTqPCR", quietly = TRUE) &&
            requireNamespace("Biobase", quietly = TRUE),
          "Bayesian path available; this test exercises the fallback only")
  m <- matrix(c(20, NA, 25, 30), nrow = 2,
              dimnames = list(c("miR-1", "miR-2"), c("S1", "S2")))
  out <- qpcr_nondetect_impute(m, log_handler = function(...) invisible())
  expect_equal(attr(out, "nondetect_method"), "lod_fallback")
  expect_equal(sum(is.na(out)), 0L)
})

test_that("qpcr_nondetect_impute returns none_needed when no NAs (Bayesian path)", {
  skip_if_not_installed("nondetects")
  skip_if_not_installed("HTqPCR")
  skip_if_not_installed("Biobase")
  m <- matrix(c(20, 25, 30, 35), nrow = 2,
              dimnames = list(c("miR-1", "miR-2"), c("S1", "S2")))
  out <- qpcr_nondetect_impute(m, log_handler = function(...) invisible())
  expect_equal(attr(out, "nondetect_method"), "none_needed")
  expect_equal(attr(out, "n_imputed"), 0L)
})
