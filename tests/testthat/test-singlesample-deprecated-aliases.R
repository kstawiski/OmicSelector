library(testthat)

# Back-compat: every former paper3_* public name must still resolve to an
# exported function that forwards to its singlesample_* counterpart.

.singlesample_alias_pairs <- c(
  "assert_method_bank_exports", "assert_row_equivariant",
  "bh_fdr_correct_blocked", "bh_fdr_correct_matched_null",
  "hanley_mcneil_auc_ci", "holm_correct_familywise", "is_row_equivariant",
  "make_loco_splits", "make_locto_splits", "matched_null_benchmark",
  "matched_null_benchmark_cv", "method_bank", "method_roster",
  "paired_auc_diff_se", "register_score_adapter", "score_call",
  "technology_lift_delong"
)

test_that("all 17 paper3_* aliases exist and are exported functions", {
  ns <- asNamespace("OmicSelector")
  exported <- getNamespaceExports("OmicSelector")
  for (stem in .singlesample_alias_pairs) {
    old <- paste0("paper3_", stem)
    new <- paste0("singlesample_", stem)
    expect_true(exists(old, envir = ns, mode = "function"),
                info = paste("missing alias", old))
    expect_true(exists(new, envir = ns, mode = "function"),
                info = paste("missing function", new))
    expect_true(old %in% exported, info = paste(old, "not exported"))
    expect_true(new %in% exported, info = paste(new, "not exported"))
  }
})

test_that("paper3_* aliases forward to their singlesample_* counterparts", {
  # No-argument / pure-read functions: outputs must be identical.
  expect_identical(paper3_method_roster(), singlesample_method_roster())
  expect_identical(paper3_method_bank(), singlesample_method_bank())

  # A representative computational alias: same result for the same input.
  set.seed(1)
  model <- list(w = rnorm(4))
  score_fun <- function(model, X, meta = NULL) as.numeric(X %*% model$w)
  X <- matrix(rnorm(20), nrow = 5, dimnames = list(NULL, paste0("f", 1:4)))
  expect_identical(
    paper3_is_row_equivariant(score_fun, model, X),
    singlesample_is_row_equivariant(score_fun, model, X)
  )

  # The manuscript pipeline's two callers resolve.
  expect_true(is.function(paper3_bh_fdr_correct_blocked))
  expect_true(is.function(paper3_matched_null_benchmark_cv))
})
