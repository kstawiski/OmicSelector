library(OmicSelector)


simulate_coda_signal <- function(n = 60, p = 14, effect = 1) {
  stopifnot(n %% 2 == 0)
  X <- matrix(rnorm(n * p, sd = 0.5), nrow = n, ncol = p)
  y <- factor(rep(c("ctrl", "case"), each = n / 2), levels = c("ctrl", "case"))
  X[y == "case", 1] <- X[y == "case", 1] + effect
  X[y == "case", 2] <- X[y == "case", 2] - effect
  colnames(X) <- paste0("feature_", seq_len(p))
  list(X = X, y = y)
}


run_coda_methods <- function(X, y, n_features = 5) {
  out <- list()

  set.seed(11)
  out$plr <- codaFS_plr_variance(
    X,
    y,
    n_features = n_features,
    max_pair_features = min(ncol(X), 40)
  )

  set.seed(11)
  out$selbal <- codaFS_selbal_wrapper(
    X,
    y,
    n_features = n_features,
    backend = "fallback",
    max_terms = n_features,
    max_pair_features = min(ncol(X), 40)
  )

  set.seed(11)
  out$codacore <- suppressWarnings(
    codaFS_codacore_wrapper(
      X,
      y,
      n_features = n_features,
      backend = "fallback",
      max_balances = 2,
      top_k = 2
    )
  )

  set.seed(11)
  out$logcontrast <- codaFS_logcontrast_lasso(
    X,
    y,
    n_features = n_features,
    nfolds = 3,
    nlambda = 25
  )

  set.seed(11)
  out$stability <- codaFS_stability_logratio(
    X,
    y,
    n_features = n_features,
    n_subsamples = 4,
    nfolds = 2,
    nlambda = 15,
    max_pair_features = min(ncol(X), 20)
  )

  out
}


test_that("CoDA feature selectors are invariant to per-sample additive shifts", {
  skip_if_not_installed("glmnet")

  dat <- simulate_coda_signal(n = 60, p = 14, effect = 1)
  shift <- matrix(rnorm(nrow(dat$X), sd = 2), nrow = nrow(dat$X), ncol = ncol(dat$X))
  X_shifted <- dat$X + shift
  colnames(X_shifted) <- colnames(dat$X)

  base <- run_coda_methods(dat$X, dat$y, n_features = 5)
  shifted <- run_coda_methods(X_shifted, dat$y, n_features = 5)

  expect_identical(base$plr$selected_features, shifted$plr$selected_features)
  expect_identical(base$selbal$selected_features, shifted$selbal$selected_features)
  expect_identical(base$codacore$selected_features, shifted$codacore$selected_features)
  expect_identical(base$logcontrast$selected_features, shifted$logcontrast$selected_features)
  expect_identical(base$stability$selected_features, shifted$stability$selected_features)
})


test_that("log-contrast lasso coefficients satisfy the zero-sum constraint after projection", {
  skip_if_not_installed("glmnet")

  dat <- simulate_coda_signal(n = 60, p = 14, effect = 1)
  fit <- codaFS_logcontrast_lasso(dat$X, dat$y, n_features = 5, nfolds = 3, nlambda = 25)

  expect_lt(abs(fit$zero_sum_error), 1e-8)
  expect_equal(sum(fit$coefficients), 0, tolerance = 1e-8)
})


test_that("CoDA selectors recover the discriminative log-ratio signal", {
  skip_if_not_installed("glmnet")

  dat <- simulate_coda_signal(n = 80, p = 14, effect = 1.2)
  fits <- run_coda_methods(dat$X, dat$y, n_features = 5)

  for (fit in fits) {
    expect_true(all(c("feature_1", "feature_2") %in% fit$selected_features))
  }
})


test_that("CoDA filters are registered in mlr3filters and can score tasks", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3filters")
  skip_if_not_installed("glmnet")

  ids <- c(
    "coda_plr_variance",
    "coda_selbal",
    "coda_codacore",
    "coda_logcontrast_lasso",
    "coda_stability_logratio"
  )
  expect_true(all(ids %in% mlr3filters::mlr_filters$keys()))

  dat <- simulate_coda_signal(n = 40, p = 14, effect = 1)
  task_df <- as.data.frame(dat$X)
  task_df$outcome <- dat$y
  task <- mlr3::TaskClassif$new(
    id = "coda_filters",
    backend = task_df,
    target = "outcome",
    positive = "case"
  )

  filter_specs <- list(
    coda_plr_variance = list(),
    coda_selbal = list(backend = "fallback", max_terms = 5L, max_pair_features = 14L),
    coda_codacore = list(backend = "fallback", max_balances = 2L, top_k = 2L),
    coda_logcontrast_lasso = list(nfolds = 3L, nlambda = 20L),
    coda_stability_logratio = list(
      n_subsamples = 3L,
      nfolds = 2L,
      nlambda = 12L,
      max_pair_features = 10L
    )
  )

  for (id in ids) {
    filt <- do.call(mlr3filters::flt, c(list(.key = id), filter_specs[[id]]))
    suppressWarnings(filt$calculate(task))
    scores <- filt$scores
    expect_type(scores, "double")
    expect_equal(sort(names(scores)), sort(colnames(dat$X)))
  }
})


test_that("CoDA feature selectors run at p = 14, 100, and 500", {
  skip_if_not_installed("glmnet")

  for (p in c(14, 100, 500)) {
    dat <- simulate_coda_signal(n = 40, p = p, effect = 0.9)

    expect_length(
      codaFS_plr_variance(dat$X, dat$y, n_features = 5, max_pair_features = min(p, 40))$selected_features,
      5
    )
    expect_length(
      codaFS_selbal_wrapper(
        dat$X,
        dat$y,
        n_features = 5,
        backend = "fallback",
        max_terms = 5,
        max_pair_features = min(p, 40)
      )$selected_features,
      5
    )
    expect_length(
      suppressWarnings(
        codaFS_codacore_wrapper(
          dat$X,
          dat$y,
          n_features = 5,
          backend = "fallback",
          max_balances = 2,
          top_k = 2
        )
      )$selected_features,
      5
    )
    expect_length(
      codaFS_logcontrast_lasso(
        dat$X,
        dat$y,
        n_features = 5,
        nfolds = 3,
        nlambda = 20
      )$selected_features,
      5
    )
    expect_length(
      codaFS_stability_logratio(
        dat$X,
        dat$y,
        n_features = 5,
        n_subsamples = 3,
        nfolds = 2,
        nlambda = 12,
        max_pair_features = min(p, 20)
      )$selected_features,
      5
    )
  }
})
