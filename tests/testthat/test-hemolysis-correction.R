library(OmicSelector)

test_that("hemolysis proxy score uses numerator minus denominator", {
  x <- c(
    MIMAT0000092 = 5.5,
    MIMAT0006764 = 3.0,
    MIMAT0005793 = 2.0
  )

  expect_equal(hemolysis_proxy_score(x), 2.5, tolerance = 1e-12)
})

test_that("hemolysis prefilter rejects strongly hemolyzed proxy samples", {
  set.seed(11)
  feature_names <- c(
    "MIMAT0000092", "MIMAT0006764", "MIMAT0005793", "MIMAT0000451",
    "MIMAT0000682", "MIMAT0000617"
  )
  controls <- matrix(
    rnorm(6 * 40, mean = 5, sd = 0.25),
    nrow = 40,
    ncol = length(feature_names),
    dimnames = list(NULL, feature_names)
  )

  fit <- fit_hemolysis_prefilter(controls, seed = 99)
  clean_idx <- which.min(hemolysis_proxy_score(controls))
  clean_sample <- controls[clean_idx, , drop = FALSE]
  hemo_sample <- clean_sample
  hemo_sample[, "MIMAT0000092"] <- hemo_sample[, "MIMAT0000092"] + 3

  pred <- apply_hemolysis_prefilter(fit, rbind(clean_sample, hemo_sample))

  expect_s3_class(fit, "omicselector_hemolysis_prefilter")
  expect_true(pred$keep[[1]])
  expect_true(pred$rejected[[2]])
  expect_gt(pred$hemolysis_score[[2]], pred$threshold[[2]])
})

test_that("weighted CLR remains invariant to uniform additive shifts", {
  x <- matrix(
    c(
      5, 3, 1, 4,
      4, 2, 0, 1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("MIMAT0000092", "MIMAT0006764", "MIMAT0000682", "MIMAT0000617"))
  )

  wclr_a <- weighted_clr_transform(x)
  wclr_b <- weighted_clr_transform(x + 11)

  expect_equal(wclr_a, wclr_b, tolerance = 1e-12)
})

test_that("weighted CLR MLP fallback fit predicts probabilities", {
  set.seed(1)
  feature_names <- c(
    "MIMAT0000682", "MIMAT0000617", "MIMAT0000086", "MIMAT0000092",
    "MIMAT0005793", "MIMAT0006764", "MIMAT0000264", "MIMAT0000765",
    "MIMAT0004909", "MIMAT0005898", "MIMAT0022727", "MIMAT0000451",
    "MIMAT0000418", "MIMAT0000090"
  )
  X <- matrix(rnorm(14 * 24), nrow = 24, ncol = 14, dimnames = list(NULL, feature_names))
  y <- factor(rep(c("control", "case"), each = 12), levels = c("control", "case"))
  X[y == "case", c("MIMAT0000682", "MIMAT0000617", "MIMAT0000086")] <-
    X[y == "case", c("MIMAT0000682", "MIMAT0000617", "MIMAT0000086")] + 0.8

  fit <- suppressWarnings(
    train_weighted_clr_mlp(X, y, backend = "glm", class_weight = FALSE, verbose = FALSE, seed = 1)
  )
  pred <- suppressWarnings(predict_weighted_clr_mlp(fit, X[1:5, , drop = FALSE]))

  expect_s3_class(fit, "omicselector_weighted_clr_mlp")
  expect_length(pred, 5)
  expect_true(all(pred >= 0 & pred <= 1))
})
