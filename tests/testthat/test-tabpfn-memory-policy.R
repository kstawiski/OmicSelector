test_that("TabPFN query budget bounds dense attention on either device", {
  skip_if_not_installed("reticulate")
  policy <- reticulate::import_from_path("tabpfn_memory",
    path = system.file("python", package = "OmicSelector"), convert = TRUE)
  # Small contexts retain the previous 1024-row outer scoring partition.
  expect_gte(policy$query_limit(list(c(20L, 32L, 8L, 4L))), 1024)
  # Toray-like cached attention exceeds the 11 GiB GPU for 1024 queries.
  layout <- c(80L, 4096L, 8L, 4L)
  limit <- policy$query_limit(list(layout))
  per_row <- 3 * prod(layout)
  expect_lte(limit * per_row, 4 * 1024^3)
  expect_gt((limit + 1) * per_row, 4 * 1024^3)
  expect_lt(limit, 1024)
  expect_equal(policy$query_limit(list(layout, c(1L, 1L, 1L, 4L))), limit)
  expect_equal(policy$query_limit(list(c(1L, 1L, 1L, 4L)), 12L), 1)
  expect_error(policy$query_limit(list(c(1L, 1L, 1L, 4L)), 11L),
               "one query.*Reduce hp\\$max_context")
  expect_error(policy$query_limit(list(layout), 0L), "budget must be positive")
  expect_error(policy$query_limit(list(c(1L, 0L, 1L, 4L))), "positive integers")
  expect_equal(policy$query_limit(list()), 1024)
  # Balanced chunks preserve every row, including the short final remainder.
  chunks <- .tabpfn_score_chunks(seq_len(15313L), as.integer(limit))
  expect_identical(unname(unlist(chunks)), seq_len(15313L))
  expect_lte(max(lengths(chunks)), limit)
  expect_lte(diff(range(lengths(chunks))), 1L)
})
