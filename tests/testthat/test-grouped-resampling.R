library(OmicSelector)

test_that("grouped stratified folds keep biological IDs together", {
  y <- rep(c(0, 1), each = 12)
  group_id <- rep(paste0("bio", 1:12), each = 2)
  folds <- os_make_grouped_stratified_folds(y, group_id = group_id, n_folds = 4, seed = 42)
  val <- os_validate_folds(folds, y = y, group_id = group_id)
  expect_s3_class(folds, "os_grouped_folds")
  expect_true(all(val$pass))
  expect_equal(sort(unlist(folds, use.names = FALSE)), seq_along(y))

  fold_for_row <- rep(seq_along(folds), lengths(folds))
  fold_by_row <- rep(NA_integer_, length(y))
  fold_by_row[unlist(folds, use.names = FALSE)] <- fold_for_row
  expect_true(all(vapply(split(fold_by_row, group_id), function(z) length(unique(z)) == 1L, logical(1L))))
})

test_that("grouped resample AUC reports pooled and fold-level summaries", {
  y <- rep(c(0, 1), each = 20)
  group_id <- rep(paste0("g", 1:20), each = 2)
  score <- y + seq_along(y) / 1000
  folds <- os_make_grouped_stratified_folds(y, group_id = group_id, n_folds = 5, seed = 1)
  out <- os_grouped_resample_auc(y, score, folds)
  expect_s3_class(out, "os_grouped_resample_auc")
  expect_gt(out$pooled_auc, 0.95)
  expect_equal(nrow(out$fold_auc), 5L)
})

test_that("mixed outcome groups fail closed by default", {
  y <- c(0, 1, 0, 1)
  group_id <- c("a", "a", "b", "b")
  expect_error(
    os_make_grouped_stratified_folds(y, group_id = group_id, n_folds = 2),
    "contain both outcome classes"
  )
})
