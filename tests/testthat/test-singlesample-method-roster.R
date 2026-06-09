# Foundation tests for the Paper 3 Amendment #4 method roster, the canonical
# single-sample score adapter, and the generic row-equivariance harness.

test_that("Amendment #4 roster loads with canonical columns and frozen membership", {
  ros <- singlesample_method_roster()
  expect_equal(nrow(ros), 74L)
  expect_true(all(c("method_id", "estimand", "role", "single_sample_gate",
                    "row_source", "lopbo_mechanism", "is_discriminator",
                    "is_negative_control") %in% names(ros)))
  expect_setequal(ros$method_id, .singlesample_frozen_roster_ids())
  expect_false(any(duplicated(ros$method_id)))
  expect_true(all(ros$single_sample_gate == "required"))
})

test_that("roster per-estimand counts match the frozen prespecification (section 4.4)", {
  ros <- singlesample_method_roster()
  expect_equal(sum(ros$role == "discriminator" & ros$estimand == "within"), 43L)
  expect_equal(sum(ros$role == "discriminator-conditional"), 1L)   # ai-tabpfn
  expect_equal(sum(ros$role == "discriminator" & ros$estimand == "transfer"), 24L)
  expect_equal(sum(ros$role == "negative_control"), 5L)
  expect_equal(sum(ros$role == "baseline"), 1L)
  # Headline discriminator bank = 43 within firm + 24 transfer = 67 (+1 conditional).
  expect_equal(sum(ros$is_discriminator & ros$role != "discriminator-conditional"), 67L)
})

test_that("roster routing columns use only the allowed enums", {
  ros <- singlesample_method_roster()
  expect_true(all(ros$estimand %in% c("within", "transfer")))
  expect_true(all(ros$role %in% c("baseline", "discriminator",
                                  "discriminator-conditional", "negative_control")))
  expect_true(all(ros$tier %in% c("R1", "R2")))
  expect_true(all(ros$pkg_status %in% c("present", "primitive", "new",
                                        "new-secondary")))
  expect_true(all(ros$row_source %in% c("per_cohort_rows",
                                        "stratum_pooled_upstream",
                                        "transfer_fold_rows")))
  expect_true(all(ros$lopbo_mechanism %in% c("within_cohort", "transfer_fold",
                                             "deferred", "m4_decision", "na")))
})

test_that("negative controls are excluded from the discriminator bank and BY family", {
  ros <- singlesample_method_roster()
  nc <- ros[ros$is_negative_control, , drop = FALSE]
  expect_equal(nrow(nc), 5L)
  expect_false(any(nc$is_discriminator))
  expect_true(all(nc$family == "NC"))
})

test_that("present-method primitives resolve in the package namespace (roster<->package binding)", {
  ros <- singlesample_method_roster()
  present <- ros[ros$pkg_status == "present", , drop = FALSE]
  fns <- unique(stats::na.omit(c(present$fit_fn, present$score_fn)))
  ns <- asNamespace("OmicSelector")
  resolved <- vapply(fns, function(f) exists(f, envir = ns, mode = "function"),
                     logical(1))
  expect_identical(fns[!resolved], character(0))
})

test_that("singlesample_method_bank() carries the Amendment #4 route columns without regressing", {
  bank <- singlesample_method_bank()
  expect_true(all(c("method_id", "estimand", "role", "single_sample_gate",
                    "negative_control", "row_source", "lopbo_mechanism") %in%
                    names(bank)))
  # At least the present roster methods that live in the bank are annotated.
  expect_true(any(!is.na(bank$method_id)))
  # Annotation agrees with the roster for every matched row.
  ros <- singlesample_method_roster()
  matched <- bank[!is.na(bank$method_id), , drop = FALSE]
  idx <- match(matched$method_id, ros$method_id)
  expect_true(all(matched$estimand == ros$estimand[idx]))
  expect_true(all(matched$role == ros$role[idx]))
  # Back-compat invariants preserved (existing bank contract).
  expect_true(all(bank$functions_exist))
  expect_true(all(bank$functions_exported))
})

test_that("row-equivariance harness PASSES a single-sample-deployable scorer", {
  set.seed(7)
  model <- list(w = rnorm(6), center = rnorm(6))
  # Each row's score depends only on that row -> row-equivariant.
  score_fun <- function(model, X, meta = NULL) {
    as.numeric(sweep(X, 2L, model$center, "-") %*% model$w)
  }
  X <- matrix(rnorm(15 * 6), nrow = 15,
              dimnames = list(paste0("S", 1:15), paste0("f", 1:6)))
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X))
  expect_true(singlesample_is_row_equivariant(score_fun, model, X))
})

test_that("row-equivariance harness REJECTS a batch-coupled scorer", {
  # Centering on the SCORING matrix's column means couples rows: the same
  # specimen scores differently depending on its batch-mates (fails the gate).
  set.seed(11)
  model <- list(w = rnorm(6))
  score_fun <- function(model, X, meta = NULL) {
    as.numeric(scale(X, center = TRUE, scale = FALSE) %*% model$w)
  }
  X <- matrix(rnorm(15 * 6), nrow = 15, dimnames = list(NULL, paste0("f", 1:6)))
  expect_error(singlesample_assert_row_equivariant(score_fun, model, X),
               regexp = "row-equivariance violated")
  expect_false(singlesample_is_row_equivariant(score_fun, model, X))
})

test_that("row-equivariance harness REJECTS a scorer that mutates model state", {
  # Model is a reference (environment): a scorer that updates running state
  # during scoring (BatchNorm-style) must fail clause (d).
  model <- new.env(parent = emptyenv())
  model$w <- c(1, -1, 0.5, 0.25)
  model$calls <- 0L
  score_fun <- function(model, X, meta = NULL) {
    model$calls <- model$calls + 1L   # mutates model state
    as.numeric(X %*% model$w)
  }
  X <- matrix(rnorm(10 * 4), nrow = 10, dimnames = list(NULL, paste0("f", 1:4)))
  expect_error(singlesample_assert_row_equivariant(score_fun, model, X),
               regexp = "model state changed")
  expect_false(singlesample_is_row_equivariant(score_fun, model, X))
})

test_that("row-equivariance harness REJECTS a duplicate-coupled scorer", {
  # Score shifts by the count of exact duplicates of each row in the batch:
  # a specimen scores differently when a copy of it is also present.
  model <- list(w = rnorm(5))
  score_fun <- function(model, X, meta = NULL) {
    key <- apply(round(X, 8), 1L, paste, collapse = "|")
    dup_count <- ave(seq_along(key), key, FUN = length)
    as.numeric(X %*% model$w) + dup_count
  }
  X <- matrix(rnorm(12 * 5), nrow = 12, dimnames = list(NULL, paste0("f", 1:5)))
  expect_false(singlesample_is_row_equivariant(score_fun, model, X))
})

test_that("default model_digest refuses external-pointer (torch-like) state", {
  digest <- OmicSelector:::.singlesample_default_model_digest
  has_ptr <- OmicSelector:::.singlesample_contains_externalptr
  # Plain / environment models digest fine and are not flagged.
  expect_false(has_ptr(list(w = rnorm(3), center = rnorm(3))))
  expect_type(digest(list(w = rnorm(3))), "raw")
  env_model <- new.env(parent = emptyenv()); env_model$w <- rnorm(2)
  expect_false(has_ptr(env_model))
  # A model carrying an external pointer (as a torch/reticulate backend would)
  # is detected and refuses the blind default digest.
  ptr <- tryCatch(methods::new("externalptr"), error = function(e) NULL)
  skip_if(is.null(ptr) || typeof(ptr) != "externalptr",
          "no external pointer constructible in this R build")
  expect_true(has_ptr(list(weights = rnorm(4), backend = ptr)))
  expect_error(digest(list(weights = rnorm(4), backend = ptr)),
               regexp = "external-pointer state")
  # With a custom digest supplied, an external-pointer model is accepted.
  model <- list(weights = c(1, 2, 3), backend = ptr)
  score_fun <- function(model, X, meta = NULL) as.numeric(X %*% model$weights)
  X <- matrix(rnorm(10 * 3), nrow = 10)
  expect_invisible(singlesample_assert_row_equivariant(
    score_fun, model, X, model_digest = function(m) m$weights))
})

test_that("row-equivariance harness propagates per-specimen metadata correctly", {
  set.seed(3)
  model <- list(w = rnorm(4))
  # Score uses a per-specimen meta offset; still row-independent.
  score_fun <- function(model, X, meta = NULL) {
    off <- if (is.null(meta)) 0 else meta$offset
    as.numeric(X %*% model$w) + off
  }
  X <- matrix(rnorm(12 * 4), nrow = 12, dimnames = list(NULL, paste0("f", 1:4)))
  meta <- data.frame(offset = rnorm(12))
  expect_invisible(singlesample_assert_row_equivariant(score_fun, model, X, meta))
})

test_that("singlesample_score_call dispatches a registered adapter", {
  reg <- OmicSelector:::.singlesample_score_adapter_registry
  on.exit(if (exists("ws-alr-pivot", envir = reg, inherits = FALSE))
    rm("ws-alr-pivot", envir = reg), add = TRUE)
  singlesample_register_score_adapter("ws-alr-pivot", function(model, X, meta) {
    as.numeric(rowSums(X))
  })
  X <- matrix(1:6, nrow = 3)
  expect_equal(singlesample_score_call("ws-alr-pivot", model = NULL, X = X),
               as.numeric(rowSums(X)))
})

test_that("singlesample_score_call errors clearly for a not-yet-implemented method", {
  # `sig-path` is a frozen-roster method whose fit/score functions are not yet
  # implemented (a Tier-R2 reticulate-iisignature method; repointed from `ot-lot`
  # when the linearized-OT/CDT method landed, to keep this error-path test
  # exercising a genuinely unimplemented slot).
  expect_error(
    singlesample_score_call("sig-path", model = NULL, X = matrix(0, 2, 3)),
    regexp = "not yet implemented")
})

test_that("singlesample_score_call validates the score output contract by type", {
  reg <- OmicSelector:::.singlesample_score_adapter_registry
  on.exit(if (exists("ws-alr-pivot", envir = reg, inherits = FALSE))
    rm("ws-alr-pivot", envir = reg), add = TRUE)
  X <- matrix(rnorm(12), nrow = 4)

  # n-by-1 numeric matrix is accepted and normalized to a bare vector.
  singlesample_register_score_adapter("ws-alr-pivot",
                                function(model, X, meta) matrix(rowSums(X), ncol = 1))
  out <- singlesample_score_call("ws-alr-pivot", model = NULL, X = X)
  expect_null(dim(out))
  expect_equal(out, rowSums(X))

  # Character output (even character-numeric) is rejected, not coerced.
  singlesample_register_score_adapter("ws-alr-pivot",
                                function(model, X, meta) as.character(rowSums(X)))
  expect_error(singlesample_score_call("ws-alr-pivot", model = NULL, X = X),
               regexp = "numeric vector")

  # Logical output is rejected.
  singlesample_register_score_adapter("ws-alr-pivot",
                                function(model, X, meta) rowSums(X) > 0)
  expect_error(singlesample_score_call("ws-alr-pivot", model = NULL, X = X),
               regexp = "numeric vector")

  # n-by-2 matrix is rejected as a dimensioned output.
  singlesample_register_score_adapter("ws-alr-pivot",
                                function(model, X, meta) cbind(rowSums(X), rowSums(X)))
  expect_error(singlesample_score_call("ws-alr-pivot", model = NULL, X = X),
               regexp = "plain numeric vector|array with dim")

  # Wrong length is rejected.
  singlesample_register_score_adapter("ws-alr-pivot",
                                function(model, X, meta) rowSums(X)[1])
  expect_error(singlesample_score_call("ws-alr-pivot", model = NULL, X = X),
               regexp = "must return 4 numeric scores")
})

test_that("every implemented non-present roster score_fn is canonically dispatchable", {
  # The canonical adapter calls score_fn(model, X, meta) positionally, so an
  # implemented method on the canonical path MUST take `model` as its first
  # argument -- pointing score_fn at a raw primitive (e.g. ws_balance_ilr(x,...))
  # silently breaks singlesample_score_call() while the direct call still works.
  ros <- singlesample_method_roster()
  ns <- asNamespace("OmicSelector")
  reg <- OmicSelector:::.singlesample_score_adapter_registry
  # Legacy primitives whose roster score_fn still points at a non-canonical
  # function (first arg is not `model`) pending their own implementation gate,
  # where they get a canonical fit_*/score_* wrapper or a registered adapter and
  # the manifest score_fn is repointed. Remove a method here once it is wrapped;
  # a NEWLY-implemented method that mis-points its score_fn will still fail.
  # dro-group was wrapped (fit_dro_group/score_dro_group) and removed here.
  pending_canonical_wrapper <- character(0)
  checked <- 0L
  for (i in seq_len(nrow(ros))) {
    sf <- ros$score_fn[i]
    if (is.na(sf) || !nzchar(sf)) next
    if (!exists(sf, envir = ns, mode = "function")) next   # not yet implemented
    if (identical(ros$pkg_status[i], "present")) next       # adapter-dispatch path
    if (ros$method_id[i] %in% pending_canonical_wrapper) next
    if (exists(ros$method_id[i], envir = reg, inherits = FALSE)) next  # registered
    fn <- get(sf, envir = ns, mode = "function")
    fm <- names(formals(fn))
    expect_identical(
      fm[1], "model",
      info = paste0(ros$method_id[i], ": score_fn '", sf,
                    "' first arg is '", fm[1],
                    "', not 'model' -> not dispatchable via ",
                    "singlesample_score_call()"))
    checked <- checked + 1L
  }
  # At least the methods implemented so far must be covered.
  expect_gte(checked, 5L)
})
