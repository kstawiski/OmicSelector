library(OmicSelector)

test_that("TCGA smoke test runs on tiny subset", {
  skip_on_cran()
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("mlr3tuning")
  skip_if_not_installed("mlr3filters")
  skip_if_not_installed("rpart")

  data("original_TCGA_data", package = "OmicSelector", envir = environment())

  feature_cols <- grep("^hsa\\.", names(original_TCGA_data), value = TRUE)
  feature_cols <- head(feature_cols, 25)

  df_all <- as.data.frame(original_TCGA_data)
  df <- df_all[, c("patient", "sample_type", feature_cols), drop = FALSE]
  df$sample_type <- factor(df$sample_type)

  idx_tumor <- which(df$sample_type == "PrimaryTumor")
  idx_normal <- which(df$sample_type == "SolidTissueNormal")

  set.seed(1)
  df <- df[c(sample(idx_tumor, 10), sample(idx_normal, 10)), , drop = FALSE]
  df <- df[sample(nrow(df)), , drop = FALSE]

  pipeline <- OmicPipeline$new(
    data = df,
    target = "sample_type",
    positive = "PrimaryTumor",
    patient_id = "patient"
  )

  learner <- pipeline$create_graph_learner(
    filter = "anova",
    model = "rpart",
    n_features = 5
  )

  result <- pipeline$benchmark(
    learners = learner,
    outer_folds = 2,
    inner_folds = 2,
    seed = 1,
    parallel = FALSE,
    threads = 1
  )

  expect_s3_class(result, "NestedCVResult")
})
