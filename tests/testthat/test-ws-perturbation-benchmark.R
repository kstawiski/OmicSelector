library(OmicSelector)

test_that("within-sample perturbation benchmark returns per-condition AUCs", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  set.seed(1)
  task_df <- data.frame(
    matrix(rnorm(14 * 30), nrow = 30, ncol = 14),
    outcome = factor(rep(c("control", "case"), each = 15), levels = c("control", "case"))
  )
  task_df[task_df$outcome == "case", 1:3] <- task_df[task_df$outcome == "case", 1:3] + 1.2

  task <- mlr3::TaskClassif$new(
    id = "ws_perturbation",
    backend = task_df,
    target = "outcome",
    positive = "case"
  )

  result <- suppressWarnings(
    ws_perturbation_benchmark(
      task = task,
      learners = list(
        mlr3::lrn("classif.rpart", predict_type = "prob"),
        mlr3::lrn(
          "classif.clr_mlp",
          backend = "glm",
          class_weight = FALSE,
          verbose = FALSE,
          seed = 1
        )
      ),
      conditions = c("baseline", "shift"),
      resampling = mlr3::rsmp("cv", folds = 3),
      seed = 1
    )
  )

  expect_s3_class(result, "ws_perturbation_benchmark")
  expect_true(all(c("baseline", "shift") %in% result$aggregate_auc$condition))
  expect_true(all(result$aggregate_auc$auc_folds > 0))
})
