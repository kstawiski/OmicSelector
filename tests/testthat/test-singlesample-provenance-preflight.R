library(testthat)

test_that("os_provenance_preflight returns NO_OVERLAP for audited-independent accessions", {
  res <- os_provenance_preflight(c("GSE145176", "GSE94533"))
  expect_equal(res$status, "NO_OVERLAP")
  expect_equal(nrow(res$overlapping_pairs), 0L)
  expect_length(res$unknown_accessions, 0L)
})

test_that("os_provenance_preflight returns KNOWN_OVERLAP for Toray pair", {
  res <- os_provenance_preflight(c("GSE211692", "GSE106817"))
  expect_equal(res$status, "KNOWN_OVERLAP")
  expect_gt(nrow(res$overlapping_pairs), 0L)
  expect_true("Toray-3D-cluster" %in% res$overlapping_pairs$block_id)
})

test_that("os_provenance_preflight returns UNKNOWN_ACCESSION for unknown accs", {
  res <- os_provenance_preflight(c("GSE_FAKE_999"))
  expect_equal(res$status, "UNKNOWN_ACCESSION")
  expect_true("GSE_FAKE_999" %in% res$unknown_accessions)
})

test_that("os_provenance_preflight ships a default manifest in inst/extdata", {
  manifest_path <- system.file("extdata", "provenance_manifest.tsv",
                                package = "OmicSelector")
  expect_true(nzchar(manifest_path))
  expect_true(file.exists(manifest_path))
  m <- utils::read.delim(manifest_path, stringsAsFactors = FALSE)
  expect_true(all(c("acc1", "acc2", "block_id") %in% names(m)))
})

test_that("os_provenance_preflight UNKNOWN takes precedence over KNOWN_OVERLAP", {
  res <- os_provenance_preflight(c("GSE211692", "GSE106817", "GSE_NOT_REGISTERED"))
  expect_equal(res$status, "UNKNOWN_ACCESSION")
  expect_true("GSE_NOT_REGISTERED" %in% res$unknown_accessions)
})

test_that("os_provenance_preflight raises error when manifest missing", {
  tmp <- tempfile(fileext = ".tsv")
  expect_error(os_provenance_preflight("GSE0", manifest_path = tmp),
               "Provenance manifest not found")
})

test_that("os_provenance_preflight respects user-supplied manifest", {
  tmp <- tempfile(fileext = ".tsv")
  utils::write.table(
    data.frame(acc1 = "GSE_A", acc2 = "GSE_B",
               block_id = "test-block", stringsAsFactors = FALSE),
    file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tmp), add = TRUE)
  res <- os_provenance_preflight(c("GSE_A", "GSE_B"), manifest_path = tmp)
  expect_equal(res$status, "KNOWN_OVERLAP")
  expect_equal(res$overlapping_pairs$block_id, "test-block")
})
