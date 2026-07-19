test_that("corrected repeated-CV uses the exact Nadeau-Bengio formula", {
  effect <- c(-0.02, 0.01, 0.03, 0.04, 0.00)
  n_test <- c(20, 21, 19, 20, 20)
  n_train <- c(80, 79, 81, 80, 80)
  strata <- paste0("s", 1:5)

  out <- singlesample_corrected_repeated_cv(
    effect, n_test, n_train, strata,
    expected_m = 5L, expected_strata = rev(strata)
  )
  correction <- 1 / 5 + mean(n_test / n_train)

  expect_true(out$complete)
  expect_equal(out$estimate, mean(effect))
  expect_equal(out$correction, correction)
  expect_equal(out$se, sqrt(correction * var(effect)))
  expect_equal(out$df, 4L)
  expect_true(out$p_zero_two_sided >= 0 && out$p_zero_two_sided <= 1)
  expect_true(out$p_above_margin >= 0 && out$p_above_margin <= 1)
  expect_true(out$p_below_minus_margin >= 0 &&
                out$p_below_minus_margin <= 1)
})

test_that("corrected repeated-CV is symmetric under method-order reversal", {
  effect <- c(-0.08, -0.03, 0.01, 0.06, 0.11)
  n_test <- c(16, 18, 17, 16, 18)
  n_train <- c(64, 62, 63, 64, 62)
  strata <- paste0("repeat", 1:5)

  forward <- singlesample_corrected_repeated_cv(
    effect, n_test, n_train, strata, 5L, strata
  )
  reverse <- singlesample_corrected_repeated_cv(
    -effect, n_test, n_train, strata, 5L, strata
  )

  expect_equal(reverse$estimate, -forward$estimate)
  expect_equal(reverse$se, forward$se)
  expect_equal(reverse$ci_low, -forward$ci_high)
  expect_equal(reverse$ci_high, -forward$ci_low)
  expect_equal(reverse$p_zero_two_sided, forward$p_zero_two_sided)
  expect_equal(reverse$p_below_minus_margin, forward$p_above_margin)
  expect_equal(reverse$p_above_margin, forward$p_below_minus_margin)
})

test_that("corrected repeated-CV treats identical zero effects conservatively", {
  strata <- paste0("s", 1:5)
  out <- singlesample_corrected_repeated_cv(
    rep(0, 5), rep(20, 5), rep(80, 5), strata, 5L, strata
  )

  expect_equal(out$estimate, 0)
  expect_equal(out$se, 0)
  expect_equal(out$ci_low, 0)
  expect_equal(out$ci_high, 0)
  expect_false(out$inference_available)
  expect_true(all(is.na(c(out$p_zero_two_sided, out$p_above_margin,
                          out$p_below_minus_margin))))
})

test_that("corrected repeated-CV fails closed on incomplete or wrong strata", {
  strata <- paste0("s", 1:5)
  args <- list(
    effect = seq(-0.02, 0.02, length.out = 5),
    n_test = rep(20, 5), n_train = rep(80, 5),
    stratum_id = strata, expected_m = 5L, expected_strata = strata
  )

  expect_error(do.call(singlesample_corrected_repeated_cv,
                       lapply(args, function(x) if (length(x) == 5L) x[-1] else x)),
               "exactly expected_m")
  wrong <- args
  wrong$stratum_id[5] <- "other"
  expect_error(do.call(singlesample_corrected_repeated_cv, wrong),
               "does not match")
  duplicate <- args
  duplicate$stratum_id[5] <- duplicate$stratum_id[4]
  expect_error(do.call(singlesample_corrected_repeated_cv, duplicate),
               "duplicates")
  nonfinite <- args
  nonfinite$effect[3] <- NA_real_
  expect_error(do.call(singlesample_corrected_repeated_cv, nonfinite),
               "finite in every")
})

test_that("matched pair AUC is group-primary and antisymmetric", {
  strata <- rep(paste0("s", 1:5), each = 8L)
  sample_id <- paste0(strata, "_p", rep(1:8, times = 5L))
  group_id <- rep(paste0("g", 1:20), each = 2L)
  y <- rep(rep(c(0L, 0L, 1L, 1L), each = 2L), times = 5L)
  score_a <- y + rep(seq_len(8L) / 100, times = 5L)
  score_b <- (1L - y) + rep(seq_len(8L) / 100, times = 5L)

  forward <- singlesample_matched_pair_auc(
    y, score_a, score_b, strata, sample_id, group_id, paste0("s", 1:5)
  )
  reverse <- singlesample_matched_pair_auc(
    y, score_b, score_a, strata, sample_id, group_id, paste0("s", 1:5)
  )

  expect_equal(forward$strata$n_test_groups, rep(4L, 5L))
  expect_equal(forward$strata$n_train_groups, rep(16L, 5L))
  expect_equal(reverse$strata$effect, -forward$strata$effect)
  expect_equal(reverse$summary$estimate, -forward$summary$estimate)
  expect_equal(reverse$summary$se, forward$summary$se)
})

test_that("matched pair AUC exposes a profile-weighted estimand", {
  strata <- rep(paste0("s", 1:5), each = 6L)
  sample_id <- paste0(strata, "_p", rep(1:6, times = 5L))
  group_id <- paste0(
    strata, "_", rep(c("g0", "g0", "g0", "g1", "g2", "g3"), times = 5L)
  )
  y <- rep(c(0L, 0L, 0L, 0L, 1L, 1L), times = 5L)
  score_a <- rep(c(0.9, 0.8, 0.7, 0.1, 0.6, 0.5), times = 5L)
  score_b <- rep(c(0.1, 0.2, 0.3, 0.4, 0.8, 0.9), times = 5L)

  profile <- singlesample_matched_pair_auc(
    y, score_a, score_b, strata, sample_id, group_id, paste0("s", 1:5),
    analysis_level = "profile"
  )
  group <- singlesample_matched_pair_auc(
    y, score_a, score_b, strata, sample_id, group_id, paste0("s", 1:5),
    analysis_level = "group"
  )

  expect_equal(profile$strata$n_test, rep(6L, 5L))
  expect_true(all(is.na(profile$strata$n_test_groups)))
  expect_equal(group$strata$n_test_groups, rep(4L, 5L))
  expect_false(isTRUE(all.equal(profile$summary$estimate,
                                group$summary$estimate)))
})

test_that("all 1,225 within-method pairs follow frozen roster order", {
  roster <- singlesample_method_roster()
  within <- roster$method_id[roster$estimand == "within"]
  pairs <- singlesample_within_method_pairs(roster)

  expect_equal(length(within), 50L)
  expect_equal(nrow(pairs), choose(50, 2))
  expect_equal(nrow(pairs), 1225L)
  expect_true(all(pairs$roster_order_a < pairs$roster_order_b))
  expect_equal(pairs$method_a[1], within[1])
  expect_equal(pairs$method_b[1], within[2])
  expect_equal(tail(pairs$method_a, 1), within[49])
  expect_equal(tail(pairs$method_b, 1), within[50])
  expect_equal(anyDuplicated(pairs$pair_id), 0L)
})

test_that("pair enumeration ignores outcome and performance columns", {
  roster <- singlesample_method_roster()
  baseline <- singlesample_within_method_pairs(roster)
  roster$outcome <- rev(seq_len(nrow(roster))) %% 2L
  roster$auc <- seq(1, 0, length.out = nrow(roster))

  expect_identical(singlesample_within_method_pairs(roster), baseline)
  expect_identical(
    singlesample_within_method_pairs(roster[c("method_id", "estimand")]),
    baseline
  )
})

.fair_test_roster <- function(n = 6L) {
  data.frame(
    method_id = c(paste0("m", seq_len(n)), "transfer"),
    estimand = c(rep("within", n), "transfer"),
    stringsAsFactors = FALSE
  )
}

.fair_test_eligibility <- function() {
  grid <- expand.grid(
    method_id = c(paste0("m", 1:6), "transfer"),
    unit_id = paste0("u", 1:6),
    stringsAsFactors = FALSE
  )
  grid$eligible <- TRUE
  grid$eligible[grid$method_id %in% c("m5", "m6") &
                  !grid$unit_id %in% paste0("u", 1:3)] <- FALSE
  grid
}

test_that("complete-support panels depend only on frozen eligibility", {
  roster <- .fair_test_roster()
  eligibility <- .fair_test_eligibility()
  base <- singlesample_complete_support_panels(
    eligibility, roster, thresholds = c(0.5, 0.8),
    primary_units = paste0("u", 1:6)
  )

  with_outcomes <- eligibility
  with_outcomes$outcome <- rep(c(0L, 1L), length.out = nrow(with_outcomes))
  with_outcomes$auc <- rev(seq_len(nrow(with_outcomes))) / nrow(with_outcomes)
  changed <- singlesample_complete_support_panels(
    with_outcomes, roster, thresholds = c(0.5, 0.8),
    primary_units = paste0("u", 1:6)
  )

  expect_identical(changed, base)
  expect_equal(base$summary$n_methods, c(6L, 4L))
  expect_equal(base$summary$n_units, c(3L, 6L))
  expect_true(all(base$summary$degenerate))
  expect_false(any(base$summary$rank_allowed))
})

test_that("complete-support panels permit ranks only above both guards", {
  roster <- .fair_test_roster(5L)
  eligibility <- expand.grid(
    method_id = roster$method_id,
    unit_id = paste0("u", 1:5), stringsAsFactors = FALSE
  )
  eligibility$eligible <- TRUE

  panel <- singlesample_complete_support_panels(
    eligibility, roster, thresholds = 0.95,
    primary_units = paste0("u", 1:5)
  )
  expect_equal(panel$summary$n_methods, 5L)
  expect_equal(panel$summary$n_units, 5L)
  expect_false(panel$summary$degenerate)
  expect_true(panel$summary$rank_allowed)

  incomplete <- eligibility[-1, ]
  expect_error(
    singlesample_complete_support_panels(
      incomplete, roster, thresholds = 0.95,
      primary_units = paste0("u", 1:5)
    ),
    "complete within-method"
  )
})

test_that("zero-difference BY adjustment keeps its frozen denominator", {
  ids <- c("a::b", "a::c", "b::c")
  p <- c(0.01, NA, 0.20)
  out <- singlesample_adjust_zero_by(ids, p, rev(ids))
  expected <- stats::p.adjust(c(0.20, NA, 0.01), method = "BY", n = 3L)

  expect_equal(out$pair_id, rev(ids))
  expect_equal(out$p_zero_by, expected)
  expect_equal(out$family_n, rep(3L, 3L))
  expect_error(singlesample_adjust_zero_by(ids[-1], p[-1], ids),
               "complete frozen")
})

test_that("directional relevance tests share one frozen two-direction family", {
  ids <- c("a::b", "a::c", "b::c")
  above <- c(0.01, 0.03, NA)
  below <- c(0.90, 0.20, 0.04)
  out <- singlesample_adjust_relevance_by(ids, above, below, ids)
  raw_interleaved <- c(rbind(above, below))
  adjusted <- matrix(
    stats::p.adjust(raw_interleaved, method = "BY", n = 6L), nrow = 2L
  )

  expect_equal(out$p_above_margin_by, adjusted[1, ])
  expect_equal(out$p_below_minus_margin_by, adjusted[2, ])
  expect_equal(out$family_n, rep(6L, 3L))
  expect_error(
    singlesample_adjust_relevance_by(ids, above, below, ids[-1]),
    "complete frozen"
  )
})
