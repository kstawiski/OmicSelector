selector_fixture <- function(seed = 20260718L) {
  set.seed(seed)
  n_groups <- 24L
  group_id <- rep(paste0("patient", seq_len(n_groups)), each = 2L)
  y_group <- rep(0:1, each = n_groups / 2L)
  y <- rep(y_group, each = 2L)
  X <- matrix(
    stats::rgamma(length(y) * 30L, shape = 2, rate = 1),
    nrow = length(y),
    dimnames = list(paste0("profile", seq_along(y)),
                    paste0("hsa-miR-", seq_len(30L)))
  )
  X[y == 1L, 1:6] <- X[y == 1L, 1:6] + 4
  list(X = X, y = y, group_id = group_id)
}

test_that("selector bank is the frozen 33-route R1 deployment family", {
  roster <- singlesample_method_roster()
  bank <- singlesample_selector_candidates(roster)

  expect_length(bank, 33L)
  expect_identical(bank[[1L]], "ws-rclr-trimmed")
  rows <- roster[match(bank, roster$method_id), , drop = FALSE]
  expect_true(all(rows$estimand == "within"))
  expect_true(all(rows$tier == "R1"))
  expect_true(all(rows$role %in% c("baseline", "discriminator")))
  expect_false(any(roster$role[match(bank, roster$method_id)] ==
                     "negative_control"))
})

test_that("selector rejects research-only and negative-control routes", {
  d <- selector_fixture()
  expect_error(
    fit_singlesample_selector(
      d$X, d$y, methods = "ai-scarf", group_id = d$group_id,
      inner_folds = 3L, verify_base = FALSE
    ),
    "tier-R1"
  )
  expect_error(
    fit_singlesample_selector(
      d$X, d$y, methods = "nc-mahalanobis", group_id = d$group_id,
      inner_folds = 3L, verify_base = FALSE
    ),
    "tier-R1"
  )
})

test_that("best-AUC selector uses grouped cross-fitting and scores independently", {
  d <- selector_fixture()
  fit <- fit_singlesample_selector(
    d$X, d$y,
    methods = c("ws-rclr-trimmed", "reo-singscore"),
    route = "best_auc", group_id = d$group_id, inner_folds = 3L,
    bootstrap_reps = 30L, seed = 11L, verify_base = FALSE
  )

  expect_s3_class(fit, "singlesample_selector")
  expect_identical(fit$status, "eligible")
  expect_true(all(fit$fold_validation$pass))
  expect_true(all(fit$candidate_audit$complete))
  expect_length(fit$selected_methods, 1L)
  expect_true(fit$selected_methods %in%
                c("ws-rclr-trimmed", "reo-singscore"))
  expect_identical(fit$contract$score_orientation, "training_only")

  probe <- d$X[1:6, , drop = FALSE]
  batch <- score_singlesample_selector(fit, probe)
  singleton <- vapply(seq_len(nrow(probe)), function(i) {
    score_singlesample_selector(fit, probe[i, , drop = FALSE])
  }, numeric(1L))
  expect_equal(batch, singleton, tolerance = 0)
  expect_equal(
    score_singlesample_selector(fit, probe[c(6, 2, 2, 1), , drop = FALSE]),
    batch[c(6, 2, 2, 1)], tolerance = 0
  )
})

test_that("lower-limit selector is deterministic and records group bootstrap", {
  d <- selector_fixture()
  args <- list(
    X_train = d$X, y_train = d$y,
    methods = c("ws-rclr-trimmed", "reo-singscore"),
    route = "lower_auc_ci", group_id = d$group_id, inner_folds = 3L,
    bootstrap_reps = 30L, seed = 27L, verify_base = FALSE
  )
  fit1 <- do.call(fit_singlesample_selector, args)
  fit2 <- do.call(fit_singlesample_selector, args)

  expect_identical(fit1$selected_methods, fit2$selected_methods)
  expect_equal(fit1$candidate_audit$lower_group_auc_ci,
               fit2$candidate_audit$lower_group_auc_ci, tolerance = 0)
  expect_true(all(is.finite(fit1$candidate_audit$lower_group_auc_ci)))
})

test_that("simplex stack freezes non-negative unit-sum weights", {
  d <- selector_fixture()
  fit <- fit_singlesample_selector(
    d$X, d$y,
    methods = c("ws-rclr-trimmed", "reo-singscore"),
    route = "simplex_stack", group_id = d$group_id, inner_folds = 3L,
    bootstrap_reps = 30L, seed = 91L, verify_base = FALSE
  )

  expect_s3_class(fit, "singlesample_selector")
  expect_identical(fit$route, "simplex_stack")
  expect_gte(min(fit$stack$weights), 0)
  expect_equal(sum(fit$stack$weights), 1, tolerance = 1e-10)
  expect_identical(fit$stack$stabilizer, 1e-6)
  expect_length(fit$selected_methods, 2L)
  expect_false(any(vapply(fit$models, inherits, logical(1L), "python.builtin.object")))

  restored <- unserialize(serialize(fit, NULL))
  probe <- d$X[1:5, , drop = FALSE]
  expect_equal(score_singlesample_selector(restored, probe),
               score_singlesample_selector(fit, probe), tolerance = 0)
  expect_equal(
    score_singlesample_selector(fit, probe[3, , drop = FALSE]),
    score_singlesample_selector(fit, probe)[[3L]], tolerance = 0
  )
})

test_that("large simplex stacks are exactly batch-versus-singleton invariant", {
  set.seed(20260719)
  k <- 33L
  n <- 100L
  X <- matrix(
    stats::rnorm(n * k), nrow = n,
    dimnames = list(NULL, paste0("f", seq_len(k)))
  )
  methods <- paste0("m", seq_len(k))
  models <- stats::setNames(lapply(seq_len(k), function(j) {
    structure(
      list(
        model = j,
        method_id = methods[[j]],
        score_fn = function(model, X, meta = NULL) X[, model],
        score_route = "direct_model"
      ),
      class = "singlesample_deployable"
    )
  }), methods)
  weights <- stats::runif(k)
  weights <- weights / sum(weights)
  selector <- structure(
    list(
      status = "eligible",
      route = "simplex_stack",
      selected_methods = methods,
      models = models,
      directions = stats::setNames(rep(1, k), methods),
      stack = list(
        center = stats::setNames(stats::rnorm(k), methods),
        scale = stats::setNames(stats::runif(k, 0.5, 2), methods),
        intercept = 0.3141592653589793,
        weights = stats::setNames(weights, methods)
      )
    ),
    class = "singlesample_selector"
  )

  batch <- score_singlesample_selector(selector, X)
  singleton <- vapply(seq_len(n), function(i) {
    score_singlesample_selector(selector, X[i, , drop = FALSE])
  }, numeric(1L))
  expect_identical(batch, singleton)
})

test_that("all selector routes reuse one candidate audit and score as a set", {
  d <- selector_fixture()
  fit <- fit_singlesample_selector(
    d$X, d$y,
    methods = c("ws-rclr-trimmed", "reo-singscore"),
    route = "all", group_id = d$group_id, inner_folds = 3L,
    bootstrap_reps = 30L, seed = 123L, verify_base = FALSE
  )

  expect_s3_class(fit, "singlesample_selector_set")
  expect_setequal(names(fit$selectors),
                  c("best_auc", "lower_auc_ci", "simplex_stack"))
  expect_true(all(vapply(fit$selectors, function(x) {
    identical(x$candidate_audit, fit$candidate_audit)
  }, logical(1L))))
  expect_identical(
    fit$fit_accounting$total_selector_fits,
    fit$fit_accounting$inner_candidate_fits +
      fit$fit_accounting$final_candidate_refits
  )
  expect_identical(fit$fit_accounting$inner_candidate_fits, 6L)
  expect_lte(fit$fit_accounting$final_candidate_refits, 2L)
  probe <- d$X[1:4, , drop = FALSE]
  scores <- score_singlesample_selector(fit, probe)
  expect_s3_class(scores, "data.frame")
  expect_identical(names(scores), names(fit$selectors))
  expect_equal(scores$best_auc,
               score_singlesample_selector(fit$selectors$best_auc, probe))
  expect_equal(scores$lower_auc_ci,
               score_singlesample_selector(fit$selectors$lower_auc_ci, probe))
  expect_equal(scores$simplex_stack,
               score_singlesample_selector(fit$selectors$simplex_stack, probe))
})

test_that("stack with fewer than two complete candidates fails explicitly", {
  d <- selector_fixture()
  fit <- fit_singlesample_selector(
    d$X, d$y, methods = "ws-rclr-trimmed", route = "simplex_stack",
    group_id = d$group_id, inner_folds = 3L, bootstrap_reps = 20L,
    seed = 4L, verify_base = FALSE
  )

  expect_s3_class(fit, "singlesample_selector_ineligible")
  expect_match(fit$reason, "fewer_than_two")
  expect_identical(fit$fit_accounting$total_selector_fits, 3L)
  expect_error(score_singlesample_selector(fit, d$X[1, , drop = FALSE]),
               "ineligible")
})

test_that("mixed-outcome biological groups fail closed", {
  d <- selector_fixture()
  d$group_id[1:2] <- c("mixed", "other")
  d$group_id[25] <- "mixed"
  expect_error(
    fit_singlesample_selector(
      d$X, d$y, methods = "ws-rclr-trimmed", group_id = d$group_id,
      inner_folds = 3L, bootstrap_reps = 20L, verify_base = FALSE
    ),
    "both outcome classes"
  )
})
