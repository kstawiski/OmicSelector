library(OmicSelector)

test_that("pairwise ratio image matches pairwise differences", {
  x <- c(5, 3, 1)
  img <- encode_ratio_image(x)

  expect_equal(dim(img), c(3, 3))
  expect_equal(img, outer(x, x, "-"))
  expect_equal(img[1, 2], x[1] - x[2])
  expect_equal(img[2, 3], x[2] - x[3])
})

test_that("pairwise ratio image is invariant to uniform additive shifts", {
  x <- c(5, 3, 1, -2)
  expect_equal(
    encode_ratio_image(x),
    encode_ratio_image(x + 7.25),
    tolerance = 1e-12
  )
})

test_that("simple, correlation, and DeepInsight encodings produce expected layouts", {
  simple <- encode_simple_grid(1:6, rows = 2, cols = 3)
  expect_equal(simple, matrix(1:6, nrow = 2, byrow = TRUE))

  corr <- encode_corr_grid(1:6, order = c(6, 5, 4, 3, 2, 1), rows = 2, cols = 3)
  expect_equal(as.vector(t(corr)), 6:1)

  skip_if_not_installed("Rtsne")
  X_train <- matrix(rnorm(8 * 6), nrow = 8, ncol = 6)
  layout <- fit_deepinsight(X_train, grid_size = 6, seed = 1)
  expect_equal(nrow(layout$positions), 6)
  expect_equal(length(unique(apply(layout$positions, 1, paste, collapse = "_"))), 6)

  encoded <- encode_deepinsight(1:6, layout)
  expect_equal(dim(encoded), c(6, 6))
  expect_equal(sum(encoded != 0), 6)
})

test_that("encode_batch supports ratio images", {
  X <- matrix(rnorm(12), nrow = 4, ncol = 3)
  imgs <- encode_batch(X, method = "ratio")
  expect_equal(dim(imgs), c(4, 3, 3))
  expect_equal(imgs[1, , ], encode_ratio_image(X[1, ]))
})

