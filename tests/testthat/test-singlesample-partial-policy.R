.partial_policy_fixture <- function() {
  pairs <- data.frame(feature_a = c("a", "c"), feature_b = c("b", "d"))
  X <- rbind(c(a = 2, b = 1, c = 1, d = 3, background = 10),
             c(a = 1, b = 2, c = 4, d = 2, background = 20))
  list(
    X = X,
    ktsp = structure(list(ktsp = list(pairs = pairs),
                          feature_universe = colnames(X)), class = "reo_ktsp_model"),
    ratio = structure(list(selected_pairs = pairs, coefficients = c(2, -1),
                           intercept = 0.25, pseudocount = 0.5,
                           feature_universe = colnames(X)), class = "reo_pairratio_model"),
    rclr = structure(list(panel = "a", weights = 1,
                          feature_universe = colnames(X),
                          transform_args = list(pseudocount = 0.5, trim_upper = 0,
                                                trim_lower = 0, exclude_features = character(),
                                                min_centering_size = 1)),
                     class = "ws_rclr_panel_model")
  )
}

.expect_policy_support <- function(score, available, fitted) {
  expect_identical(attr(score, "support"), data.frame(
    n_available = rep.int(as.integer(available), length(score)),
    n_fitted = rep.int(as.integer(fitted), length(score)),
    full_support = rep.int(available == fitted, length(score))))
}

test_that("pair policies expose reduced support and reject every missing required feature in strict mode", {
  dat <- .partial_policy_fixture()
  scorers <- list(ktsp = score_reo_ktsp, ratio = score_reo_pairratio)
  for (id in names(scorers)) {
    score <- scorers[[id]]
    model <- dat[[id]]
    full <- score(model, dat$X)
    partial_X <- dat$X[, c("a", "b"), drop = FALSE]
    partial <- score(model, partial_X)
    expect_type(partial, "double")
    .expect_policy_support(full, 2L, 2L)
    .expect_policy_support(partial, 1L, 2L)
    expect_identical(full, score(model, dat$X, strict = TRUE))
    expect_identical(full, score(model, dat$X[, 1:4], strict = TRUE))
    expect_identical(partial, score(model, as.data.frame(partial_X)))
    for (missing in c("a", "b", "c", "d")) {
      incomplete <- dat$X[, setdiff(colnames(dat$X), missing), drop = FALSE]
      expect_no_error(score(model, incomplete))
      err <- tryCatch(score(model, incomplete, strict = TRUE), error = identity)
      expect_s3_class(err, "singlesample_input_error")
      expect_identical(err$code, "MISSING_REQUIRED_FEATURES")
      expect_match(conditionMessage(err), missing, fixed = TRUE)
    }
    for (strict in c(FALSE, TRUE)) {
      expect_error(score(model, dat$X[, "background", drop = FALSE], strict = strict),
                   "ALL_REQUIRED_ABSENT")
      expect_error(score(model, partial_X * 0, strict = strict), "ALL_NONDETECT")
      expect_error(score(model, dat$X[FALSE, ], strict = strict), "EMPTY_INPUT")
    }
    expect_error(score(model, dat$X[, c("a", "c"), drop = FALSE]), "NO_USABLE_SUPPORT")
  }
  expect_identical(as.numeric(score_reo_ktsp(dat$ktsp, dat$X[, c("a", "b")])), c(0, 1))
  expected <- 0.25 + 2 * (log(dat$X[, "a"] + 0.5) - log(dat$X[, "b"] + 0.5))
  expect_identical(as.numeric(score_reo_pairratio(dat$ratio, dat$X[, c("a", "b")])), expected)
  full_expected <- expected - (log(dat$X[, "c"] + 0.5) - log(dat$X[, "d"] + 0.5))
  expect_equal(as.numeric(score_reo_pairratio(dat$ratio, dat$X)), full_expected)
})

test_that("rCLR support includes background and transforms the supplied matrix before selecting the panel", {
  dat <- .partial_policy_fixture()
  model <- dat$rclr
  full <- score_ws_rclr_panel(model, dat$X, strict = TRUE)
  partial_X <- dat$X[, c("a", "b"), drop = FALSE]
  partial <- score_ws_rclr_panel(model, partial_X)
  expected <- log(partial_X[, "a"] + 0.5) - rowMeans(log(partial_X + 0.5))
  expect_equal(as.numeric(partial), expected)
  expect_gt(max(abs(as.numeric(full) - as.numeric(partial))), 0.1)
  .expect_policy_support(full, 5L, 5L)
  .expect_policy_support(partial, 2L, 5L)
  expect_identical(full, score_ws_rclr_panel(model, dat$X))
  for (missing in colnames(dat$X)) {
    incomplete <- dat$X[, setdiff(colnames(dat$X), missing), drop = FALSE]
    expect_error(score_ws_rclr_panel(model, incomplete, strict = TRUE),
                 "MISSING_REQUIRED_FEATURES|ALL_REQUIRED_ABSENT")
    if (missing != "a") expect_no_error(score_ws_rclr_panel(model, incomplete))
  }
  expect_error(score_ws_rclr_panel(model, dat$X[, -1]), "ALL_REQUIRED_ABSENT")
  extra <- score_ws_rclr_panel(model, cbind(dat$X, extra = 1e6), strict = TRUE)
  .expect_policy_support(extra, 5L, 5L)
  expect_gt(max(abs(as.numeric(full) - as.numeric(extra))), 0.1)
  model$feature_universe <- NULL
  expect_error(score_ws_rclr_panel(model, dat$X, strict = TRUE), "MISSING_FITTED_UNIVERSE")
})

test_that("all three policies reject invalid strict options and preserve specimen-local values and support", {
  dat <- .partial_policy_fixture()
  scorers <- list(ktsp = score_reo_ktsp, ratio = score_reo_pairratio, rclr = score_ws_rclr_panel)
  for (id in names(scorers)) {
    score <- scorers[[id]]
    model <- dat[[id]]
    state <- serialize(model, NULL)
    set.seed(510L)
    rng <- .Random.seed
    for (bad in list(NA, NULL, logical(), c(TRUE, FALSE), 1, "TRUE")) {
      expect_error(score(model, dat$X, strict = bad), "INVALID_STRICT")
    }
    for (strict in c(FALSE, TRUE)) {
      inputs <- list(dat$X)
      if (!strict) {
        inputs <- c(inputs, list(dat$X[, c("a", "b"), drop = FALSE]))
        if (id == "rclr") inputs <- c(inputs, list(dat$X[, "a", drop = FALSE]))
      }
      for (X in inputs) {
        batch <- score(model, X, strict = strict)
        singleton <- lapply(seq_len(nrow(X)), function(i) score(model, X[i, , drop = FALSE], strict = strict))
        expect_identical(as.numeric(batch), vapply(singleton, as.numeric, numeric(1)))
        support <- do.call(rbind, lapply(singleton, attr, which = "support"))
        expect_identical(attr(batch, "support"), support)
        expect_identical(batch, score(model, X[, rev(seq_len(ncol(X))), drop = FALSE], strict = strict))
        perm <- c(2L, 1L, 2L)
        permuted <- score(model, X[perm, , drop = FALSE], strict = strict)
        expect_identical(as.numeric(permuted), as.numeric(batch)[perm])
        expect_equal(attr(permuted, "support")$n_available, attr(batch, "support")$n_available[perm])
      }
    }
    expect_identical(serialize(model, NULL), state)
    expect_identical(.Random.seed, rng)
  }
})

test_that("the published partial-panel examples execute and report their support", {
  for (name in c("score_reo_ktsp", "score_reo_pairratio", "score_ws_rclr_panel")) {
    rd <- tools::parse_Rd(test_path("..", "..", "man", paste0(name, ".Rd")))
    con <- textConnection("code", "w", local = TRUE)
    tools::Rd2ex(rd, out = con)
    close(con)
    env <- new.env(parent = environment())
    expect_no_error(capture.output(eval(parse(text = code), envir = env)))
    .expect_policy_support(env$full, if (name == "score_ws_rclr_panel") 8L else 2L,
                           if (name == "score_ws_rclr_panel") 8L else 2L)
    .expect_policy_support(env$partial, 1L, if (name == "score_ws_rclr_panel") 8L else 2L)
  }
})
