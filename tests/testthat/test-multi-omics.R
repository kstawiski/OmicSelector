# Tests for Multi-Omics Support

test_that("validate_omics_input handles single data.frame", {
  data <- data.frame(
    gene1 = rnorm(50),
    gene2 = rnorm(50),
    outcome = factor(rep(c("A", "B"), 25))
  )

  result <- validate_omics_input(data, target = "outcome")

  expect_s3_class(result, "OmicsInput")
  expect_false(result$is_multi_omics)
  expect_equal(length(result$modalities), 1)
  expect_equal(names(result$modalities), "main")
  expect_equal(length(result$target_data), 50)
})

test_that("validate_omics_input handles named list of matrices", {
  set.seed(42)
  n <- 30

  rna_data <- data.frame(
    rna_gene1 = rnorm(n),
    rna_gene2 = rnorm(n),
    outcome = factor(sample(c("Case", "Control"), n, replace = TRUE))
  )
  rownames(rna_data) <- paste0("sample", 1:n)

  mirna_data <- data.frame(
    mir21 = rnorm(n),
    mir155 = rnorm(n)
  )
  rownames(mirna_data) <- paste0("sample", 1:n)

  multi_data <- list(
    rna = rna_data,
    mirna = mirna_data
  )

  result <- validate_omics_input(multi_data, target = "outcome")

  expect_s3_class(result, "OmicsInput")
  expect_true(result$is_multi_omics)
  expect_equal(length(result$modalities), 2)
  expect_true("rna" %in% names(result$modalities))
  expect_true("mirna" %in% names(result$modalities))

  # Check namespaced features
  expect_true("rna::rna_gene1" %in% result$feature_map$rna)
  expect_true("mirna::mir21" %in% result$feature_map$mirna)
})

test_that("validate_omics_input rejects unnamed list", {
  data_list <- list(
    data.frame(x = 1:10),
    data.frame(y = 1:10)
  )

  expect_error(
    validate_omics_input(data_list, target = "outcome"),
    "named list"
  )
})

test_that("validate_omics_input rejects target in multiple modalities", {
  data1 <- data.frame(x = 1:10, outcome = factor(rep(c("A", "B"), 5)))
  data2 <- data.frame(y = 1:10, outcome = factor(rep(c("A", "B"), 5)))

  expect_error(
    validate_omics_input(list(mod1 = data1, mod2 = data2), target = "outcome"),
    "multiple modalities"
  )
})

test_that("merge_omics_data creates namespaced data.frame", {
  set.seed(42)
  n <- 20

  rna_data <- data.frame(
    gene1 = rnorm(n),
    gene2 = rnorm(n),
    outcome = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  rownames(rna_data) <- paste0("s", 1:n)

  prot_data <- data.frame(
    prot1 = rnorm(n),
    prot2 = rnorm(n)
  )
  rownames(prot_data) <- paste0("s", 1:n)

  multi_data <- list(rna = rna_data, prot = prot_data)
  omics_input <- validate_omics_input(multi_data, target = "outcome")

  merged <- merge_omics_data(omics_input)

  # Check namespaced columns exist
  expect_true("rna::gene1" %in% names(merged))
  expect_true("rna::gene2" %in% names(merged))
  expect_true("prot::prot1" %in% names(merged))
  expect_true("prot::prot2" %in% names(merged))
})

test_that("get_modality_info returns correct summary", {
  set.seed(42)
  n <- 15

  data <- list(
    rna = data.frame(
      g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n),
      outcome = factor(rep(c("A", "B", "C"), 5))
    ),
    mirna = data.frame(
      m1 = rnorm(n), m2 = rnorm(n)
    )
  )
  rownames(data$rna) <- paste0("s", 1:n)
  rownames(data$mirna) <- paste0("s", 1:n)

  omics_input <- validate_omics_input(data, target = "outcome")
  info <- get_modality_info(omics_input)

  expect_equal(nrow(info), 2)
  expect_true("rna" %in% info$modality)
  expect_true("mirna" %in% info$modality)

  rna_row <- info[info$modality == "rna", ]
  expect_equal(rna_row$n_features, 3)  # g1, g2, g3
  expect_true(rna_row$has_target)

  mirna_row <- info[info$modality == "mirna", ]
  expect_equal(mirna_row$n_features, 2)  # m1, m2
  expect_false(mirna_row$has_target)
})

test_that("OmicPipeline accepts multi-omics input", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 40

  multi_data <- list(
    rna = data.frame(
      gene1 = rnorm(n),
      gene2 = rnorm(n),
      outcome = factor(sample(c("Case", "Control"), n, replace = TRUE))
    ),
    mirna = data.frame(
      mir21 = rnorm(n),
      mir155 = rnorm(n)
    )
  )
  rownames(multi_data$rna) <- paste0("sample", 1:n)
  rownames(multi_data$mirna) <- paste0("sample", 1:n)

  pipeline <- OmicPipeline$new(
    data = multi_data,
    target = "outcome",
    positive = "Case"
  )

  expect_s3_class(pipeline, "OmicPipeline")
  expect_true(pipeline$is_multi_omics())

  # Check namespaced features
  features <- pipeline$get_feature_names()
  expect_true(any(grepl("^rna::", features)))
  expect_true(any(grepl("^mirna::", features)))

  # Check modality info
  mod_info <- pipeline$get_modality_info()
  expect_equal(nrow(mod_info), 2)
})

test_that("OmicPipeline get_modality_features works", {
  skip_if_not_installed("R6")
  skip_if_not_installed("mlr3")

  set.seed(42)
  n <- 30

  multi_data <- list(
    rna = data.frame(
      gene1 = rnorm(n), gene2 = rnorm(n),
      outcome = factor(rep(c("A", "B"), n/2))
    ),
    prot = data.frame(
      prot1 = rnorm(n), prot2 = rnorm(n), prot3 = rnorm(n)
    )
  )
  rownames(multi_data$rna) <- paste0("s", 1:n)
  rownames(multi_data$prot) <- paste0("s", 1:n)

  pipeline <- OmicPipeline$new(multi_data, target = "outcome", positive = "A")

  rna_features <- pipeline$get_modality_features("rna")
  expect_equal(length(rna_features), 2)  # gene1, gene2
  expect_true(all(grepl("^rna::", rna_features)))

  prot_features <- pipeline$get_modality_features("prot")
  expect_equal(length(prot_features), 3)  # prot1, prot2, prot3
  expect_true(all(grepl("^prot::", prot_features)))

  # Invalid modality
  expect_error(pipeline$get_modality_features("invalid"), "not found")
})

test_that("print.OmicsInput shows correct info", {
  data <- list(
    rna = data.frame(g1 = 1:5, g2 = 1:5, y = factor(c("A", "A", "B", "B", "A"))),
    mirna = data.frame(m1 = 1:5)
  )
  rownames(data$rna) <- letters[1:5]
  rownames(data$mirna) <- letters[1:5]

  omics_input <- validate_omics_input(data, target = "y")

  output <- capture.output(print(omics_input))

  expect_true(any(grepl("Multi-Omics", output)))
  expect_true(any(grepl("rna", output)))
  expect_true(any(grepl("mirna", output)))
})
