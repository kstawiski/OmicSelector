.registry_data <- function() {
  X <- matrix(seq_len(72) / 7 + 1, 12, 6,
               dimnames = list(NULL, paste0("feature", 1:6)))
  list(X = X, y = rep(0:1, each = 6))
}

.registry_fit <- function(X_train, y_train, ...) {
  list(panel = colnames(X_train)[1L], center = mean(X_train[, 1L]))
}

.registry_score <- function(model, X, meta = NULL, ...) {
  as.numeric(X[, model$panel] - model$center)
}

.registry_expect_code <- function(expr, code) {
  err <- tryCatch(force(expr), error = identity)
  expect_s3_class(err, "singlesample_input_error")
  expect_identical(err$code, code)
  invisible(err)
}

test_that("built-in roster bytes match the base commit", {
  base <- readRDS(test_path("fixtures", "custom-scorer-base-roster.rds"))
  expect_identical(registered_singlesample_methods()$method_id, character())
  expect_identical(singlesample_method_roster(), base)
  expect_identical(serialize(singlesample_method_roster(), NULL),
                    serialize(base, NULL))
  register_singlesample_method("test_local", .registry_fit, .registry_score)
  on.exit(unregister_singlesample_method("test_local"), add = TRUE)
  roster <- singlesample_method_roster()
  expect_identical(roster[seq_len(nrow(base)), ], base)
  expect_identical(roster$row_source[nrow(roster)], "user")
  expect_identical(names(roster), names(base))
})

test_that("a local scorer registers, deploys, verifies and scores alone", {
  dat <- .registry_data()
  set.seed(29)
  before_seed <- .Random.seed
  expect_true(register_singlesample_method("test_local", .registry_fit,
                                            .registry_score, probe = dat))
  on.exit(unregister_singlesample_method("test_local"), add = TRUE)
  expect_identical(.Random.seed, before_seed)
  rows <- registered_singlesample_methods()
  expect_identical(rows$method_id, "test_local")
  expect_true(rows$audited)
  expect_false(rows$strict_support)
  dep <- deploy_singlesample(dat$X, dat$y, method = "test_local")
  expect_s3_class(dep, "singlesample_deployable")
  expect_true(dep$verified)
  before <- serialize(dep$model, NULL)
  batch <- score_specimen(dep, dat$X)
  single <- vapply(seq_len(nrow(dat$X)), function(i) {
    score_specimen(dep, dat$X[i, ])
  }, numeric(1))
  expect_identical(batch, single)
  expect_true(all(is.finite(batch)))
  expect_identical(batch, singlesample_score_call("test_local", dep$model, dat$X))
  expect_identical(batch, singlesample_score_call("test_local", dep$model,
                                                  as.data.frame(dat$X)))
  expect_true(is_singlesample_deployable(dep, dat$X))
  expect_identical(serialize(dep$model, NULL), before)
  restored <- unserialize(serialize(dep, NULL))
  expect_identical(score_specimen(restored, dat$X), batch)
  expect_true(register_singlesample_method("test_local", .registry_fit,
    .registry_score, probe = unname(dat), overwrite = TRUE))
})

test_that("batch centering and fitted-state mutation refuse registration", {
  batch_score <- function(model, X, meta = NULL) X[, 1L] - mean(X[, 1L])
  err <- .registry_expect_code(register_singlesample_method(
    "test_batch", .registry_fit, batch_score), "REGISTRATION_AUDIT_FAILED")
  expect_match(conditionMessage(err), "row-equivariance")
  expect_false("test_batch" %in% registered_singlesample_methods()$method_id)

  mutable_fit <- function(X_train, y_train) {
    model <- new.env(parent = emptyenv())
    model$calls <- 0L
    model
  }
  mutable_score <- function(model, X, meta = NULL) {
    model$calls <- model$calls + 1L
    as.numeric(X[, 1L])
  }
  err <- .registry_expect_code(register_singlesample_method(
    "test_mutating", mutable_fit, mutable_score), "REGISTRATION_AUDIT_FAILED")
  expect_match(conditionMessage(err), "model state changed")
  expect_false("test_mutating" %in% registered_singlesample_methods()$method_id)
})

test_that("invalid identifiers, clashes and signatures have coded errors", {
  for (id in list(NULL, character(), NA_character_, "", "bad-name", "1name",
                  "if", "two words", c("first", "second"), 1)) {
    .registry_expect_code(register_singlesample_method(
      id, .registry_fit, .registry_score), "REGISTRATION_INVALID_ID")
  }
  for (id in c("reo-ktsp", "dann")) {
    for (overwrite in c(FALSE, TRUE)) {
      .registry_expect_code(register_singlesample_method(
        id, .registry_fit, .registry_score, overwrite = overwrite),
        "REGISTRATION_ID_EXISTS")
    }
  }
  for (fn in list(NULL, 1, sum, function(...) NULL,
                  function(X_train) NULL,
                  function(X_train, y_train, required) NULL)) {
    .registry_expect_code(register_singlesample_method(
      "test_invalid", fn, .registry_score, audit = FALSE),
      "REGISTRATION_INVALID_FUNCTION")
  }
  for (fn in list(NULL, 1, sum, function(model, X) NULL,
                  function(model, X, meta) NULL,
                  function(model, X, meta = NULL, required) NULL)) {
    .registry_expect_code(register_singlesample_method(
      "test_invalid", .registry_fit, fn, audit = FALSE),
      "REGISTRATION_INVALID_FUNCTION")
  }
  .registry_expect_code(register_singlesample_method(
    "test_invalid", .registry_fit, function(model, X, meta = NULL) X[, 1],
    strict_support = TRUE), "REGISTRATION_INVALID_FUNCTION")
  expect_false("test_invalid" %in% registered_singlesample_methods()$method_id)
})

test_that("replacement is atomic and existing deployments retain their scorer", {
  dat <- .registry_data()
  register_singlesample_method("test_replace", .registry_fit, .registry_score)
  on.exit(unregister_singlesample_method("test_replace"), add = TRUE)
  dep <- deploy_singlesample(dat$X, dat$y, "test_replace")
  scores <- score_specimen(dep, dat$X)
  .registry_expect_code(register_singlesample_method(
    "test_replace", .registry_fit, .registry_score), "REGISTRATION_ID_EXISTS")
  .registry_expect_code(register_singlesample_method(
    "test_replace", .registry_fit,
    function(model, X, meta = NULL) X[, 1] - mean(X[, 1]), overwrite = TRUE),
    "REGISTRATION_AUDIT_FAILED")
  expect_identical(singlesample_score_call("test_replace", dep$model, dat$X), scores)
  register_singlesample_method("test_replace", .registry_fit,
    function(model, X, meta = NULL) rep(7, nrow(X)), overwrite = TRUE)
  expect_identical(singlesample_score_call("test_replace", dep$model, dat$X), rep(7, 12))
  expect_identical(score_specimen(dep, dat$X), scores)
  expect_true(unregister_singlesample_method("test_replace"))
  expect_false(unregister_singlesample_method("test_replace"))
  expect_false("test_replace" %in% singlesample_method_roster()$method_id)
  expect_error(deploy_singlesample(dat$X, dat$y, "test_replace"), "unknown method")
  expect_error(singlesample_score_call("test_replace", dep$model, dat$X), "not found")
  expect_identical(score_specimen(dep, dat$X), scores)
  expect_true(is_singlesample_deployable(dep, dat$X))
})

test_that("strict requests follow the registered and frozen support flag", {
  dat <- .registry_data()
  register_singlesample_method("test_strict", .registry_fit, .registry_score)
  on.exit(unregister_singlesample_method("test_strict"), add = TRUE)
  unsupported <- deploy_singlesample(dat$X, dat$y, "test_strict")
  .registry_expect_code(score_specimen(unsupported, dat$X, strict = TRUE),
                        "STRICT_UNSUPPORTED")
  .registry_expect_code(singlesample_score_call("test_strict", unsupported$model,
                                               dat$X, strict = TRUE), "STRICT_UNSUPPORTED")
  strict_score <- function(model, X, meta = NULL, strict = FALSE) {
    out <- as.numeric(X[, model$panel] - model$center)
    attr(out, "support") <- data.frame(strict = rep(strict, nrow(X)))
    out
  }
  register_singlesample_method("test_strict", .registry_fit, strict_score,
                                strict_support = TRUE, overwrite = TRUE)
  dep <- deploy_singlesample(dat$X, dat$y, "test_strict")
  expect_true(all(attr(score_specimen(dep, dat$X, strict = TRUE), "support")$strict))
  expect_false(any(attr(score_specimen(dep, dat$X), "support")$strict))
  expect_identical(score_specimen(dep, dat$X, strict = TRUE),
    singlesample_score_call("test_strict", dep$model, dat$X, strict = TRUE))
  .registry_expect_code(score_specimen(unsupported, dat$X, strict = TRUE),
                        "STRICT_UNSUPPORTED")
  .registry_expect_code(score_specimen(dep, dat$X, strict = NA), "INVALID_STRICT")
  unregister_singlesample_method("test_strict")
  expect_true(all(attr(score_specimen(dep, dat$X, strict = TRUE), "support")$strict))
  # A declaration also audits strict-mode batch dependence.
  .registry_expect_code(register_singlesample_method("test_strict", .registry_fit,
    function(model, X, meta = NULL, strict = FALSE) {
      if (strict) X[, 1L] - mean(X[, 1L]) else as.numeric(X[, 1L])
    }, strict_support = TRUE), "REGISTRATION_AUDIT_FAILED")
})

test_that("custom scoring retains the deployment input guards", {
  dat <- .registry_data()
  register_singlesample_method("test_guards", .registry_fit, .registry_score)
  on.exit(unregister_singlesample_method("test_guards"), add = TRUE)
  dep <- deploy_singlesample(dat$X, dat$y, "test_guards", nondetect_limit = 1)
  .registry_expect_code(score_specimen(dep, dat$X[FALSE, ]), "EMPTY_INPUT")
  .registry_expect_code(score_specimen(dep, "invalid"), "INVALID_INPUT")
  .registry_expect_code(score_specimen(dep, dat$X, meta = data.frame(a = 1)),
                        "INVALID_INPUT")
  .registry_expect_code(score_specimen(dep, unname(dat$X)), "INVALID_FEATURE_NAMES")
  .registry_expect_code(score_specimen(dep, dat$X[, -1]), "ALL_REQUIRED_ABSENT")
  empty <- dat$X
  empty[2, 1] <- 1
  err <- .registry_expect_code(score_specimen(dep, empty), "ALL_NONDETECT")
  expect_identical(err$rows, 2L)
  mask <- matrix(TRUE, nrow(dat$X), ncol(dat$X), dimnames = dimnames(dat$X))
  mask[3, 1] <- FALSE
  .registry_expect_code(score_specimen(dep, dat$X, detected = mask), "ALL_NONDETECT")
  .registry_expect_code(singlesample_score_call("test_guards", dep$model,
                                               dat$X[, -1]), "ALL_REQUIRED_ABSENT")
})

test_that("audit failures and opt-out are explicit without registry residue", {
  for (probe in list(list(), list(X = matrix(1, 3, 2), y = 0:1),
                     list(X = .registry_data()$X, y = rep(0, 12)))) {
    .registry_expect_code(register_singlesample_method(
      "test_audit", .registry_fit, .registry_score, probe = probe),
      "REGISTRATION_AUDIT_FAILED")
  }
  for (score in list(function(model, X, meta = NULL) rep(Inf, nrow(X)),
                     function(model, X, meta = NULL) rep("1", nrow(X)))) {
    .registry_expect_code(register_singlesample_method(
      "test_audit", .registry_fit, score), "REGISTRATION_AUDIT_FAILED")
    register_singlesample_method("test_audit", .registry_fit, score, audit = FALSE)
    .registry_expect_code(singlesample_score_call("test_audit", list(),
                                                 .registry_data()$X), "CUSTOM_INVALID_SCORES")
    unregister_singlesample_method("test_audit")
  }
  register_singlesample_method("test_audit", .registry_fit,
    function(model, X, meta = NULL) X[, 1] - mean(X[, 1]), audit = FALSE)
  on.exit(unregister_singlesample_method("test_audit"), add = TRUE)
  expect_false(registered_singlesample_methods()$audited)
  dat <- .registry_data()
  expect_error(deploy_singlesample(dat$X, dat$y, "test_audit"), "row-equivariance")
})

test_that("roster benchmark helpers handle custom method identifiers", {
  base_pairs <- singlesample_within_method_pairs()
  register_singlesample_method("test_descriptive", .registry_fit, .registry_score)
  on.exit(unregister_singlesample_method("test_descriptive"), add = TRUE)
  expect_identical(singlesample_within_method_pairs(), base_pairs)
  register_singlesample_method("test_within", .registry_fit, .registry_score,
                                estimand = "within")
  on.exit(unregister_singlesample_method("test_within"), add = TRUE)
  pairs <- singlesample_within_method_pairs()
  expect_true("test_within" %in% pairs$method_b)
  expect_true("test_within" %in% singlesample_selector_candidates())
  ids <- singlesample_method_roster()$method_id
  eligibility <- expand.grid(method_id = ids, unit_id = c("a", "b"),
                               stringsAsFactors = FALSE)
  eligibility$eligible <- TRUE
  panels <- singlesample_complete_support_panels(eligibility, min_methods = 1,
                                                  min_units = 1)
  expect_true("test_within" %in% panels$method_membership$method_id)
  expect_false("test_descriptive" %in% panels$method_membership$method_id)
  dat <- .registry_data()
  fit <- fit_singlesample_selector(dat$X, dat$y, methods = "test_within",
    inner_folds = 2L, bootstrap_reps = 20L, verify_base = TRUE)
  expect_s3_class(fit, "singlesample_selector")
  scores <- score_singlesample_selector(fit, dat$X)
  expect_true(all(is.finite(scores)))
  expect_equal(score_singlesample_selector(fit, dat$X[1, , drop = FALSE]), scores[1])
})
