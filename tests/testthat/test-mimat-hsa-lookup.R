# Tests for MIMAT <-> hsa-miR namespace lookup (R/mimat_hsa_lookup.R).
# Covers the codex Task A v2 criteria:
#   1) Known-good resolution (positive case per direction)
#   2) NA on unknown input
#   3) Vector handling (order + length preservation)
#   4) Partial-unknown vectors (mix of known + unknown)
#   5) Case-insensitivity (mixed-case input -> canonical output)
#   6) Round-trip identity (MIMAT -> hsa -> MIMAT)
#   7) Idempotency on already-canonical input

# ---------------------------------------------------------------------------
# 1. Known-good resolution
# ---------------------------------------------------------------------------

test_that("resolve_mimat_to_hsa returns canonical hsa-miR name for known MIMAT", {
  expect_equal(resolve_mimat_to_hsa("MIMAT0000062"), "hsa-let-7a-5p")
  expect_equal(resolve_mimat_to_hsa("MIMAT0004481"), "hsa-let-7a-3p")
})

test_that("resolve_hsa_to_mimat returns canonical MIMAT for known hsa-miR name", {
  expect_equal(resolve_hsa_to_mimat("hsa-let-7a-5p"), "MIMAT0000062")
  expect_equal(resolve_hsa_to_mimat("hsa-let-7a-3p"), "MIMAT0004481")
})

# ---------------------------------------------------------------------------
# 2. NA on unknown input
# ---------------------------------------------------------------------------

test_that("resolve_mimat_to_hsa returns NA for unknown MIMAT", {
  expect_true(is.na(resolve_mimat_to_hsa("MIMAT9999999")))
  expect_true(is.na(resolve_mimat_to_hsa("not-a-mimat")))
})

test_that("resolve_hsa_to_mimat returns NA for unknown hsa-miR name", {
  expect_true(is.na(resolve_hsa_to_mimat("hsa-not-real")))
  expect_true(is.na(resolve_hsa_to_mimat("MIMAT0000062"))) # wrong direction
})

# ---------------------------------------------------------------------------
# 3. Vector handling (order + length preservation)
# ---------------------------------------------------------------------------

test_that("resolve_mimat_to_hsa preserves input length and order", {
  ids <- c("MIMAT0000062", "MIMAT0004481", "MIMAT0010195")
  out <- resolve_mimat_to_hsa(ids)
  expect_length(out, length(ids))
  expect_equal(out, c("hsa-let-7a-5p", "hsa-let-7a-3p", "hsa-let-7a-2-3p"))
})

test_that("resolve_hsa_to_mimat preserves input length and order", {
  nms <- c("hsa-let-7a-5p", "hsa-let-7a-3p", "hsa-let-7a-2-3p")
  out <- resolve_hsa_to_mimat(nms)
  expect_length(out, length(nms))
  expect_equal(out, c("MIMAT0000062", "MIMAT0004481", "MIMAT0010195"))
})

# ---------------------------------------------------------------------------
# 4. Partial-unknown vectors (mix of known + unknown)
# ---------------------------------------------------------------------------

test_that("resolve_mimat_to_hsa handles mixed known + unknown input", {
  out <- resolve_mimat_to_hsa(c("MIMAT0000062", "MIMAT9999999", "MIMAT0004481"))
  expect_equal(out[c(1, 3)], c("hsa-let-7a-5p", "hsa-let-7a-3p"))
  expect_true(is.na(out[2]))
  expect_length(out, 3L)
})

test_that("resolve_hsa_to_mimat handles mixed known + unknown input", {
  out <- resolve_hsa_to_mimat(c("hsa-let-7a-5p", "hsa-not-real", "hsa-let-7a-3p"))
  expect_equal(out[c(1, 3)], c("MIMAT0000062", "MIMAT0004481"))
  expect_true(is.na(out[2]))
  expect_length(out, 3L)
})

# ---------------------------------------------------------------------------
# 5. Case-insensitivity (input case is normalised; output is canonical case)
# ---------------------------------------------------------------------------

test_that("resolve_mimat_to_hsa is case-insensitive on input and returns canonical output", {
  variants <- c("MIMAT0000062", "mimat0000062", "Mimat0000062")
  out <- resolve_mimat_to_hsa(variants)
  expect_equal(out, rep("hsa-let-7a-5p", 3L))
})

test_that("resolve_hsa_to_mimat is case-insensitive on input and returns canonical output", {
  variants <- c("hsa-let-7a-5p", "HSA-LET-7A-5P", "Hsa-Let-7a-5p")
  out <- resolve_hsa_to_mimat(variants)
  expect_equal(out, rep("MIMAT0000062", 3L))
})

# ---------------------------------------------------------------------------
# 6. Round-trip identity (MIMAT -> hsa -> MIMAT)
# ---------------------------------------------------------------------------

test_that("MIMAT -> hsa -> MIMAT round trip is identity on known IDs", {
  ids <- c("MIMAT0000062", "MIMAT0004481", "MIMAT0010195")
  expect_equal(resolve_hsa_to_mimat(resolve_mimat_to_hsa(ids)), ids)
})

test_that("hsa -> MIMAT -> hsa round trip is identity on known names", {
  nms <- c("hsa-let-7a-5p", "hsa-let-7a-3p", "hsa-let-7a-2-3p")
  expect_equal(resolve_mimat_to_hsa(resolve_hsa_to_mimat(nms)), nms)
})

# ---------------------------------------------------------------------------
# 7. Idempotency: canonical input returns canonical output unchanged
# ---------------------------------------------------------------------------

test_that("resolve_mimat_to_hsa output is itself canonical (idempotent on already-canonical mapping)", {
  # The output of resolve_mimat_to_hsa is canonical hsa-miR; mapping it back
  # through the resolver in the opposite direction should be stable.
  hsa <- resolve_mimat_to_hsa("MIMAT0000062")
  expect_equal(hsa, "hsa-let-7a-5p")
  expect_equal(resolve_hsa_to_mimat(hsa), "MIMAT0000062")
})

test_that("input type validation: non-character input errors out", {
  expect_error(resolve_mimat_to_hsa(123L))
  expect_error(resolve_hsa_to_mimat(123L))
})
