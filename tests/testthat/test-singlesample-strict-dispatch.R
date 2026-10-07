.strict_dispatch_fixture <- function() {
  X <- rbind(c(a = 2, b = 1, c = 1, d = 3, background = 10),
             c(a = 1, b = 2, c = 4, d = 2, background = 20),
             c(a = 4, b = 1, c = 5, d = 2, background = 15),
             c(a = 3, b = 2, c = 2, d = 1, background = 12))
  pairs <- data.frame(feature_a = c("a", "c"), feature_b = c("b", "d"))
  models <- list(
    "reo-ktsp" = structure(list(ktsp = list(pairs = pairs),
      feature_universe = colnames(X)), class = "reo_ktsp_model"),
    "reo-pairratio" = structure(list(selected_pairs = pairs,
      coefficients = c(2, -1), intercept = 0.25, pseudocount = 0.5,
      feature_universe = colnames(X)), class = "reo_pairratio_model"),
    "ws-rclr-trimmed" = structure(list(panel = "a", weights = 1,
      feature_universe = colnames(X),
      transform_args = list(pseudocount = 0.5, trim_upper = 0, trim_lower = 0,
        exclude_features = character(), min_centering_size = 1)),
      class = "ws_rclr_panel_model"),
    "ws-balance-ilr" = fit_ws_balance_ilr(X, c(0, 0, 1, 1), hp = list(
      sbp = list(ab_cd = list(numerator = c("a", "b"),
                             denominator = c("c", "d"))),
      min_sbp_features = 2, pseudocount = 0.5)))
  deps <- lapply(names(models), function(id) structure(list(
    method_id = id, model = models[[id]], score_route = "score_call"),
    class = "singlesample_deployable"))
  names(deps) <- names(models)
  list(X = X, deps = deps)
}

.strict_expect_code <- function(expr, code) {
  err <- tryCatch(force(expr), error = identity)
  expect_s3_class(err, "singlesample_input_error")
  expect_identical(err$code, code)
}

test_that("canonical strict scoring accepts full support and rejects partial support", {
  dat <- .strict_dispatch_fixture()
  for (id in names(dat$deps)) {
    dep <- dat$deps[[id]]
    full <- score_specimen(dep, dat$X, strict = TRUE)
    expect_identical(full, score_specimen(dep, dat$X, strict = FALSE))
    expect_identical(full, singlesample_score_call(id, dep$model, dat$X, strict = TRUE))
    expect_identical(full, score_specimen(dep, as.data.frame(dat$X), strict = TRUE))
    missing <- if (id == "ws-rclr-trimmed") "background" else "d"
    partial <- dat$X[, setdiff(colnames(dat$X), missing), drop = FALSE]
    expect_no_error(score_specimen(dep, partial))
    .strict_expect_code(score_specimen(dep, partial, strict = TRUE),
                        "MISSING_REQUIRED_FEATURES")
    .strict_expect_code(singlesample_score_call(id, dep$model, partial, strict = TRUE),
                        "MISSING_REQUIRED_FEATURES")
    expect_true(all(attr(full, "support")$full_support))
    # Fitted pairs and fixed ILR balances do not require unused background.
    if (id != "ws-rclr-trimmed") {
      expect_identical(full, score_specimen(dep, dat$X[, 1:4], strict = TRUE))
    }
  }
})

test_that("strict singleton and batch scores and support agree exactly", {
  dat <- .strict_dispatch_fixture()
  for (dep in dat$deps) {
    batch <- score_specimen(dep, dat$X, strict = TRUE)
    singles <- lapply(seq_len(nrow(dat$X)), function(i)
      score_specimen(dep, dat$X[i, ], strict = TRUE))
    expect_identical(as.numeric(batch), vapply(singles, as.numeric, numeric(1)))
    support <- do.call(rbind, lapply(singles, attr, "support"))
    rownames(support) <- NULL
    expect_identical(attr(batch, "support"), support)
    mask <- dat$X > 0
    expect_identical(batch, score_specimen(dep, dat$X, detected = mask, strict = TRUE))
    mask[2, ] <- FALSE
    err <- tryCatch(score_specimen(dep, dat$X, detected = mask, strict = TRUE),
                    error = identity)
    expect_identical(err$code, "ALL_NONDETECT")
    expect_identical(err$rows, 2L)
  }
})

test_that("strict mode fails closed for unsupported routes and custom adapters", {
  dat <- .strict_dispatch_fixture()
  dep <- dat$deps[[1]]
  dep$method_id <- "ws-alr-pivot"
  .strict_expect_code(score_specimen(dep, dat$X, strict = TRUE), "STRICT_UNSUPPORTED")
  .strict_expect_code(singlesample_score_call("ws-alr-pivot", NULL, dat$X,
                                             strict = TRUE), "STRICT_UNSUPPORTED")
  dep <- dat$deps[[1]]
  for (route in c("direct_model", "direct_x", "unknown")) {
    dep$score_route <- route
    .strict_expect_code(score_specimen(dep, dat$X, strict = TRUE), "STRICT_UNSUPPORTED")
  }
  reg <- OmicSelector:::.singlesample_score_adapter_registry
  id <- "reo-ktsp"
  had <- exists(id, reg, inherits = FALSE)
  old <- if (had) get(id, reg, inherits = FALSE)
  on.exit(if (had) assign(id, old, reg) else rm(list = id, envir = reg))
  singlesample_register_score_adapter(id, function(model, X, meta) rep(0, nrow(X)))
  expect_identical(singlesample_score_call(id, NULL, dat$X), rep(0, nrow(dat$X)))
  .strict_expect_code(score_specimen(dat$deps[[1]], dat$X, strict = TRUE),
                      "STRICT_UNSUPPORTED")
  rm(list = id, envir = reg)
  roster <- singlesample_method_roster()
  roster$score_fn[roster$method_id == id] <- "score_ws_rclr_panel"
  .strict_expect_code(singlesample_score_call(id, NULL, dat$X, roster = roster,
                                             strict = TRUE), "STRICT_UNSUPPORTED")
  if (had) assign(id, old, reg) else assign(id, function(...) NULL, reg)
})

test_that("strict must be one nonmissing logical value at both interfaces", {
  dat <- .strict_dispatch_fixture()
  for (strict in list(NULL, logical(), NA, 1, "TRUE", c(TRUE, FALSE))) {
    .strict_expect_code(score_specimen(dat$deps[[1]], dat$X, strict = strict),
                        "INVALID_STRICT")
    .strict_expect_code(singlesample_score_call("reo-ktsp", dat$deps[[1]]$model,
                                               dat$X, strict = strict), "INVALID_STRICT")
    .strict_expect_code(score_ws_balance_ilr(dat$deps[[4]]$model, dat$X,
                                             strict = strict), "INVALID_STRICT")
  }
})

test_that("ILR support counts unique fitted members and expands sentinel balances", {
  dat <- .strict_dispatch_fixture()
  fixed <- dat$deps[["ws-balance-ilr"]]$model
  expect_identical(attr(score_ws_balance_ilr(fixed, dat$X), "support"),
    data.frame(n_available = rep(4L, 4), n_fitted = rep(4L, 4),
               full_support = rep(TRUE, 4)))
  partial <- dat$X[, -4, drop = FALSE]
  expect_identical(attr(score_ws_balance_ilr(fixed, partial), "support"),
    data.frame(n_available = rep(3L, 4), n_fitted = rep(4L, 4),
               full_support = rep(FALSE, 4)))
  .strict_expect_code(score_ws_balance_ilr(fixed, partial, strict = TRUE),
                      "MISSING_REQUIRED_FEATURES")
  fixed$sbp$repeat_balance <- fixed$sbp[[1]]
  fixed$sbp$repeat_balance$numerator <- c("a", "b", "never_fitted")
  expect_true(all(attr(score_ws_balance_ilr(fixed, dat$X, strict = TRUE),
                       "support")$n_fitted == 4L))
  for (named in c(FALSE, TRUE)) {
    for (sentinel in c("ALL_NON_RBC", "ALL_NON_PLT", "TOP_K_BY_ABUNDANCE")) {
      model <- dat$deps[["ws-balance-ilr"]]$model
      model$sbp$dynamic <- if (sentinel == "TOP_K_BY_ABUNDANCE")
        list(numerator = sentinel, denominator = "TAIL_BY_ABUNDANCE", k = 2L) else
        list(numerator = c("a", "b"), denominator = sentinel)
      if (named) {
        names(model$sbp$dynamic$numerator) <- seq_along(model$sbp$dynamic$numerator)
        names(model$sbp$dynamic$denominator) <- seq_along(model$sbp$dynamic$denominator)
      }
      full <- score_ws_balance_ilr(model, dat$X, strict = TRUE)
      expect_identical(full, score_ws_balance_ilr(model, dat$X))
      expect_true(all(attr(full, "support")$n_fitted == 5L))
      expect_no_error(score_ws_balance_ilr(model, dat$X[, 1:4]))
      .strict_expect_code(score_ws_balance_ilr(model, dat$X[, 1:4], strict = TRUE),
                          "MISSING_REQUIRED_FEATURES")
    }
  }
  model <- dat$deps[["ws-balance-ilr"]]$model
  model$feature_universe <- NULL
  .strict_expect_code(score_ws_balance_ilr(model, dat$X, strict = TRUE),
                      "MISSING_FITTED_UNIVERSE")
})
