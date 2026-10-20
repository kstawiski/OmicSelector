.grouped_inner_data <- function() {
  set.seed(712)
  X <- exp(matrix(rnorm(48 * 16), 48,
                   dimnames = list(paste0("s", 1:48), paste0("f", 1:16))))
  y <- rep(0:1, each = 24)
  X[y == 1, 1:4] <- X[y == 1, 1:4] * 4
  list(X = X, y = y, groups = rep(paste0("g", 1:24), each = 2))
}

.grouped_inner_check <- function(fold, groups) {
  expect_true(all(fold > 0L))
  expect_true(all(vapply(split(fold, groups), function(z) length(unique(z)) == 1L,
                         logical(1L))))
}

test_that("all four CV routes keep groups intact and preserve their NULL folds", {
  d <- .grouped_inner_data()
  routes <- list(
    pairratio = function(g) OmicSelector:::.reo_pairratio_foldid(d$y, 3L, 19L, g),
    dro = function(g) OmicSelector:::.gdro_stratified_folds(d$y, 3L, 19L, g),
    sig = function(g) OmicSelector:::.sig_path_cv_folds(d$y, 3L, g, 19L),
    scatter = function(g) OmicSelector:::.inv_scatter_cv_folds(d$y, 3L, g, 19L)
  )
  for (fn in routes) {
    set.seed(34); before <- .Random.seed
    fold <- fn(d$groups)
    .grouped_inner_check(fold, d$groups)
    expect_identical(.Random.seed, before)
    expect_identical(fold, fn(d$groups))
    expect_true(all(table(fold, d$y) > 0L))
  }
  # Literal pre-change assignment algorithms, including their RNG order.
  set.seed(19)
  old_pair <- integer(48)
  for (cl in 0:1) {
    idx <- which(d$y == cl)
    idx <- sample(idx, length(idx), replace = FALSE)
    old_pair[idx] <- rep(1:3, length.out = length(idx))
  }
  expect_identical(routes$pairratio(NULL), old_pair)
  set.seed(19)
  old_dro <- integer(48)
  for (cl in 0:1) {
    idx <- which(d$y == cl)
    old_dro[idx] <- sample(rep_len(1:3, length(idx)))[seq_along(idx)]
  }
  expect_identical(routes$dro(NULL), old_dro)
  expect_identical(routes$sig(NULL), rep(rep(1:3, 8), 2))
  expect_identical(routes$scatter(NULL), rep(rep(1:3, 8), 2))
  # Mixed outcomes and unequal sizes must not split a person across folds.
  mixed <- rep(paste0("p", 1:12), length.out = 48)
  for (fn in routes) .grouped_inner_check(fn(mixed), mixed)
  unequal <- rep(1:12, times = rep(c(2L, 6L), 6L))
  for (fn in routes) .grouped_inner_check(fn(unequal), unequal)
  arbitrary <- d$groups
  arbitrary[arbitrary == "g1"] <- ""
  for (fn in routes) .grouped_inner_check(fn(arbitrary), arbitrary)
})

test_that("calibration and prototype episodes keep whole groups on each side", {
  d <- .grouped_inner_data()
  set.seed(19)
  cal <- OmicSelector:::.ss_grouped_calibration_indices(d$groups, .3)
  expect_length(intersect(d$groups[cal], d$groups[-cal]), 0L)
  m <- fit_conformal_anomaly(d$X, k = 3L, seed = 19L, group_train = d$groups)
  expect_identical(m$X_ref_manifold, d$X[-cal, , drop = FALSE])
  expect_true(all(is.finite(m$cal_scores)))
  for (seed in 1:5) {
    unequal <- c(rep("large", 40), rep("small", 8))
    m <- fit_conformal_anomaly(d$X, k = 10L, seed = seed, group_train = unequal)
    expect_identical(m$n_ref, 40L)
    expect_identical(m$n_cal, 8L)
  }
  hp <- list(n_support = 3L, n_query = 4L, epochs = 6L, seed = 19L)
  for (groups in list(d$groups, rep(1:12, length.out = 48))) {
    set.seed(19)
    plan <- OmicSelector:::.proto_net_episode_plan(which(d$y == 1L),
      which(d$y == 0L), hp, groups)
    for (ep in plan) {
      expect_length(intersect(groups[ep$support], groups[ep$query]), 0L)
      expect_setequal(unique(d$y[ep$support]), 0:1)
      expect_setequal(unique(d$y[ep$query]), 0:1)
      expect_identical(ep$query_label, d$y[ep$query])
      for (idx in list(ep$support, ep$query)) {
        expect_setequal(idx, which(groups %in% groups[idx]))
      }
    }
  }
})

test_that("group errors have the documented code and adapters tolerate legacy fits", {
  d <- .grouped_inner_data()
  fits <- c("fit_reo_pairratio", "fit_dro_group", "fit_sig_path", "fit_inv_scatter",
            "fit_proto_net", "fit_singlesample_selector",
            "deploy_singlesample")
  for (fn in fits) {
    for (g in list(1:3, c(rep("a", 47), NA), matrix(1:48, 48))) {
      err <- tryCatch(do.call(get(fn, asNamespace("OmicSelector")),
                              list(d$X, d$y, group_train = g)), error = identity)
      expect_s3_class(err, "singlesample_input_error")
      expect_identical(err$code, "INVALID_GROUPS")
    }
  }
  expect_error(fit_conformal_anomaly(d$X, group_train = NA), "INVALID_GROUPS")
  expect_error(fit_singlesample_selector(d$X, d$y, group_id = 1:48,
      group_train = d$groups), "INVALID_GROUPS")
  expect_error(OmicSelector:::fit_group_dro_scorer(d$X, d$y, NULL, group_train = NA),
               "INVALID_GROUPS")
  expect_error(OmicSelector:::.ss_grouped_inner_folds(d$y, rep("one", 48), 3),
               "INVALID_GROUPS")
  expect_error(OmicSelector:::.reo_pairratio_foldid(d$y, 3, 1, rep(1:2, each = 24)),
               "INVALID_GROUPS")
  expect_error(OmicSelector:::.ss_grouped_inner_folds(d$y,
    c(rep("case", 24), rep(1:12, each = 2)), 3), "INVALID_GROUPS")
  legacy <- function(X_train, y_train) list(n = nrow(X_train), y = y_train)
  a <- OmicSelector:::.singlesample_deploy_fit(legacy, "legacy", "legacy",
      d$X, d$y, NULL, NULL, group_train = d$groups)
  expect_identical(a, legacy(d$X, d$y))
  capture <- function(X_train, y_train, group_train = NULL) group_train
  expect_identical(OmicSelector:::.singlesample_deploy_fit(capture, "capture",
      "capture", d$X, d$y, NULL, NULL, group_train = d$groups), d$groups)
  # Older external benchmark adapter: grouped fits use the package adapter.
  never <- function(...) stop("old adapter invoked")
  m <- OmicSelector:::.ss_fit_benchmark_route(never, "fit_sig_path", d$X,
      d$y, NULL, list(L = 2L), d$groups)
  expect_identical(m, fit_sig_path(d$X, d$y, hp = list(L = 2L), group_train = d$groups))
  # A route without inner splits keeps its existing metadata/argument mapping.
  old_adapter <- function(fit_fn_name, X_train, y_train, meta_train, hp) {
    list(route = fit_fn_name, metadata = meta_train, hp = hp)
  }
  meta <- data.frame(cohort = rep("a", 48))
  expect_identical(OmicSelector:::.ss_fit_benchmark_route(old_adapter,
      "fit_frozen_ruv", d$X, d$y, meta, list(k = 2L), d$groups),
      old_adapter("fit_frozen_ruv", d$X, d$y, meta, list(k = 2L)))
})

test_that("grouped deployment passes behavioural verification", {
  d <- .grouped_inner_data()
  for (route in c("reo-ktsp", "sig-path", "reo-pairratio", "dro-group")) {
    if (route == "reo-pairratio") skip_if_not_installed("glmnet")
    m <- suppressWarnings(deploy_singlesample(d$X, d$y, method = route,
                                               group_train = d$groups))
    expect_true(m$verified)
    before <- serialize(m, NULL, version = 3)
    score <- score_specimen(m, d$X)
    expect_true(all(is.finite(score)))
    expect_identical(as.numeric(score_specimen(m, d$X[1, , drop = FALSE])),
                     as.numeric(score[1]))
    expect_identical(serialize(m, NULL, version = 3), before)
  }
})

test_that("selector forwards aligned groups through inner and final fits", {
  d <- .grouped_inner_data()
  seen <- list()
  original <- OmicSelector::deploy_singlesample
  testthat::local_mocked_bindings(deploy_singlesample = function(
      X_train, y_train, method, meta_train, annotation, verify, group_train = NULL) {
    seen[[length(seen) + 1L]] <<- list(rows = rownames(X_train), groups = group_train)
    original(X_train, y_train, method, meta_train, annotation, verify,
              group_train = group_train)
  }, .package = "OmicSelector")
  out <- fit_singlesample_selector(d$X, d$y, methods = "sig-path",
      group_train = d$groups, inner_folds = 2L, bootstrap_reps = 20L)
  expect_false(inherits(out, "singlesample_selector_ineligible"))
  expect_length(seen, 3L)
  for (s in seen) expect_identical(s$groups, d$groups[match(s$rows, rownames(d$X))])
})

# The base snapshot is generated by output/base-fixtures.R against a separately
# installed base commit. It is optional in ordinary package test runs.
for (.route in c("pairratio", "dro", "sig", "scatter", "proto", "stable", "conformal")) {
  local({
    route <- .route
    test_that(paste("installed base models and scores are identical with NULL groups:", route), {
      path <- Sys.getenv("OMICSELECTOR_GROUPED_BASE_FIXTURE")
      skip_if(!nzchar(path) || !file.exists(path), "separate base installation snapshot not supplied")
      d <- readRDS(path); z <- d$cases[[route]]
      skip_if(!is.null(z$error), paste("base backend unavailable:", z$error))
      venv <- "/home/konrad/.virtualenvs/omicselector_torch/bin/python"
      if (file.exists(venv) && !nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
        Sys.setenv(RETICULATE_PYTHON = venv)
      }
      fit <- get(z$fit, asNamespace("OmicSelector"))
      score <- get(z$score, asNamespace("OmicSelector"))
      set.seed(991)
      m <- suppressWarnings(do.call(fit, c(z$args, list(group_train = NULL))))
      expect_identical(m, z$model)
      actual <- if (route == "conformal") score(d$X, m)$score else score(m, d$X)
      expect_identical(actual, z$scores)
      set.seed(991)
      omitted <- suppressWarnings(do.call(fit, z$args))
      expect_identical(omitted, z$model)
      set.seed(991)
      grouped <- suppressWarnings(do.call(fit, c(z$args, list(group_train = d$groups))))
      scorer <- if (route == "conformal") function(model, X, meta = NULL) {
        os_conformal_anomaly(X, model)$score
      } else score
      expect_true(all(is.finite(scorer(grouped, d$X))))
      expect_silent(singlesample_assert_row_equivariant(scorer, grouped, d$X))
    })
  })
}
