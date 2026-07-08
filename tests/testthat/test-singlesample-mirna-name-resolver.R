# Tests for Paper 3 miRNA alias resolver (singlesample-mirna-name-resolver.R).
# Run via:
#   Rscript -e 'source("R/singlesample-mirna-name-resolver.R"); \
#               testthat::test_file("tests/testthat/test-singlesample-mirna-name-resolver.R", \
#                                   reporter="minimal")'

# ---------------------------------------------------------------------------
# mirna_alias_table structure
# ---------------------------------------------------------------------------

test_that("mirna_alias_table returns a data.frame with required columns", {
  tbl <- mirna_alias_table()
  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("mirna_name", "mimat", "mimat_pre", "aliases") %in% colnames(tbl)))
})

test_that("mirna_alias_table contains at least 60 unique miRNAs", {
  tbl <- mirna_alias_table()
  expect_gte(nrow(tbl), 60L)
  expect_equal(length(unique(tbl$mirna_name)), nrow(tbl))
})

test_that("mirna_alias_table covers all ws_default_sbp members", {
  tbl <- mirna_alias_table()
  sbp_members <- c(
    # rbc_vs_rest
    "hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p", "hsa-miR-144-3p",
    # platelet_vs_rest
    "hsa-miR-223-3p", "hsa-miR-126-3p",
    # let7 balances
    "hsa-let-7a-5p", "hsa-let-7g-5p", "hsa-let-7c-5p", "hsa-let-7e-5p",
    "hsa-let-7i-5p", "hsa-let-7b-5p", "hsa-let-7d-5p", "hsa-let-7f-5p",
    # miR-17 cluster
    "hsa-miR-17-5p", "hsa-miR-18a-5p", "hsa-miR-19a-3p",
    "hsa-miR-19b-3p", "hsa-miR-20a-5p", "hsa-miR-92a-3p",
    # miR-200 / miR-141
    "hsa-miR-200a-3p", "hsa-miR-200b-3p", "hsa-miR-200c-3p",
    "hsa-miR-141-3p", "hsa-miR-141-5p", "hsa-miR-429",
    # miR-371-373 / miR-302
    "hsa-miR-371a-3p", "hsa-miR-372-3p", "hsa-miR-373-3p",
    "hsa-miR-302a-3p", "hsa-miR-302b-3p", "hsa-miR-302c-3p", "hsa-miR-302d-3p"
  )
  missing <- setdiff(sbp_members, tbl$mirna_name)
  expect_equal(length(missing), 0L,
               info = paste("Missing from alias table:", paste(missing, collapse = ", ")))
})

test_that("mirna_alias_table covers all ws_default_pivot_pool members", {
  tbl <- mirna_alias_table()
  pivots <- c("hsa-miR-103a-3p", "hsa-miR-191-5p", "hsa-miR-26a-5p",
              "hsa-miR-30c-5p", "hsa-let-7g-5p", "hsa-miR-93-5p")
  missing <- setdiff(pivots, tbl$mirna_name)
  expect_equal(length(missing), 0L,
               info = paste("Missing pivots:", paste(missing, collapse = ", ")))
})

test_that("mirna_alias_table covers the five canonical haemolysis markers", {
  tbl <- mirna_alias_table()
  haem <- c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
            "hsa-miR-144-3p", "hsa-miR-92a-3p")
  missing <- setdiff(haem, tbl$mirna_name)
  expect_equal(length(missing), 0L,
               info = paste("Missing haem markers:", paste(missing, collapse = ", ")))
})

# ---------------------------------------------------------------------------
# resolve_mirna_aliases — canonical name lookup
# ---------------------------------------------------------------------------

test_that("resolve_mirna_aliases returns identity for canonical mirna_name inputs", {
  names_in <- c("hsa-miR-451a", "hsa-let-7a-5p", "hsa-miR-200c-3p")
  out <- resolve_mirna_aliases(names_in, target_namespace = "mirna_name")
  expect_equal(out, names_in)
})

# ---------------------------------------------------------------------------
# resolve_mirna_aliases — Toray / MIMAT namespace
# ---------------------------------------------------------------------------

test_that("resolve_mirna_aliases maps Toray MIMAT IDs to canonical names", {
  # MIMAT0001631 = hsa-miR-451a, MIMAT0000062 = hsa-let-7a-5p
  out <- resolve_mirna_aliases(c("MIMAT0001631", "MIMAT0000062"),
                               target_namespace = "mirna_name")
  expect_equal(out[1L], "hsa-miR-451a")
  expect_equal(out[2L], "hsa-let-7a-5p")
})

test_that("resolve_mirna_aliases MIMAT lookup is case-insensitive", {
  out_upper <- resolve_mirna_aliases("MIMAT0001631")
  out_lower <- resolve_mirna_aliases("mimat0001631")
  out_mixed <- resolve_mirna_aliases("Mimat0001631")
  expect_equal(out_upper, out_lower)
  expect_equal(out_upper, out_mixed)
  expect_equal(out_upper, "hsa-miR-451a")
})

test_that("resolve_mirna_aliases strips leading/trailing whitespace", {
  out <- resolve_mirna_aliases(c("  hsa-miR-451a  ", " MIMAT0001631 "),
                               target_namespace = "mirna_name")
  expect_equal(out[1L], "hsa-miR-451a")
  expect_equal(out[2L], "hsa-miR-451a")
})

# ---------------------------------------------------------------------------
# resolve_mirna_aliases — alias / legacy name lookup
# ---------------------------------------------------------------------------

test_that("resolve_mirna_aliases resolves common legacy name variants", {
  # "hsa-miR-451" (no 'a' suffix) is in aliases for hsa-miR-451a
  out <- resolve_mirna_aliases("hsa-miR-451", target_namespace = "mirna_name")
  expect_equal(out, "hsa-miR-451a")
})

test_that("resolve_mirna_aliases resolves abbreviated miRNA name without hsa prefix", {
  out <- resolve_mirna_aliases("miR-122", target_namespace = "mirna_name")
  expect_equal(out, "hsa-miR-122-5p")
})

# ---------------------------------------------------------------------------
# resolve_mirna_aliases — round-trip
# ---------------------------------------------------------------------------

test_that("round-trip mirna_name -> mimat -> mirna_name is identity", {
  tbl <- mirna_alias_table()
  # Use a subset of entries where the MIMAT is unambiguous (not shared)
  dup_mimat <- tbl$mimat[duplicated(tbl$mimat)]
  sub_tbl   <- tbl[!tbl$mimat %in% dup_mimat, ]
  sample_names <- sub_tbl$mirna_name[seq_len(min(20L, nrow(sub_tbl)))]

  step1 <- resolve_mirna_aliases(sample_names, target_namespace = "mimat")
  step2 <- resolve_mirna_aliases(step1,        target_namespace = "mirna_name")
  expect_equal(step2, sample_names)
})

# ---------------------------------------------------------------------------
# resolve_mirna_aliases — unresolved features handling
# ---------------------------------------------------------------------------

test_that("unresolved features return NA when keep_unresolved = TRUE", {
  out <- resolve_mirna_aliases(c("hsa-miR-451a", "not_a_real_mirna"),
                               keep_unresolved = TRUE)
  expect_equal(out[1L], "hsa-miR-451a")
  expect_true(is.na(out[2L]))
})

test_that("unresolved features raise an error when keep_unresolved = FALSE", {
  expect_error(
    resolve_mirna_aliases("not_a_real_mirna", keep_unresolved = FALSE),
    regexp = "could not be resolved"
  )
})

test_that("all-NA input returns all-NA output when keep_unresolved = TRUE", {
  out <- resolve_mirna_aliases(c("fake1", "fake2"), keep_unresolved = TRUE)
  expect_equal(out, c(NA_character_, NA_character_))
})

# ---------------------------------------------------------------------------
# resolve_mirna_aliases — reverse direction (canonical name → MIMAT)
# ---------------------------------------------------------------------------

test_that("resolve_mirna_aliases maps canonical names to MIMAT accessions", {
  out <- resolve_mirna_aliases(c("hsa-miR-451a", "hsa-let-7a-5p"),
                               target_namespace = "mimat")
  expect_equal(out[1L], "MIMAT0001631")
  expect_equal(out[2L], "MIMAT0000062")
})

# ---------------------------------------------------------------------------
# apply_mirna_aliases — named numeric vector
# ---------------------------------------------------------------------------

test_that("apply_mirna_aliases renames elements of a named numeric vector", {
  v <- c(MIMAT0001631 = 12.3, MIMAT0000062 = 4.1, unknown_probe = 0.9)
  out <- apply_mirna_aliases(v, keep_unresolved = TRUE)
  expect_named(out, c("hsa-miR-451a", "hsa-let-7a-5p", "unknown_probe"))
  expect_equal(unname(out), unname(v))
})

test_that("apply_mirna_aliases drops unresolved elements when keep_unresolved = FALSE", {
  v <- c(MIMAT0001631 = 12.3, junk_id = 0.5)
  out <- apply_mirna_aliases(v, keep_unresolved = FALSE)
  expect_equal(length(out), 1L)
  expect_named(out, "hsa-miR-451a")
})

# ---------------------------------------------------------------------------
# apply_mirna_aliases — matrix
# ---------------------------------------------------------------------------

test_that("apply_mirna_aliases renames columns of a samples × features matrix", {
  M <- matrix(runif(6L), nrow = 2L,
              dimnames = list(c("S1", "S2"),
                              c("MIMAT0001631", "MIMAT0000062", "junk")))
  out <- apply_mirna_aliases(M, keep_unresolved = TRUE)
  expect_equal(dim(out), dim(M))
  expect_equal(rownames(out), c("S1", "S2"))
  expect_equal(colnames(out)[1L], "hsa-miR-451a")
  expect_equal(colnames(out)[2L], "hsa-let-7a-5p")
  expect_equal(colnames(out)[3L], "junk")  # unresolved, kept as-is
})

test_that("apply_mirna_aliases drops unresolved columns when keep_unresolved = FALSE", {
  M <- matrix(runif(4L), nrow = 2L,
              dimnames = list(c("S1", "S2"),
                              c("MIMAT0001631", "junk")))
  out <- apply_mirna_aliases(M, keep_unresolved = FALSE)
  expect_equal(ncol(out), 1L)
  expect_equal(colnames(out), "hsa-miR-451a")
})

test_that("apply_mirna_aliases preserves row names of matrix", {
  M <- matrix(runif(4L), nrow = 2L,
              dimnames = list(c("patient_A", "patient_B"),
                              c("MIMAT0000070", "MIMAT0001631")))
  out <- apply_mirna_aliases(M)
  expect_equal(rownames(out), c("patient_A", "patient_B"))
})

test_that("apply_mirna_aliases errors on matrix without column names", {
  M <- matrix(1:4, nrow = 2L)
  expect_error(apply_mirna_aliases(M), regexp = "column names")
})

# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

test_that("resolve_mirna_aliases handles empty character vector", {
  out <- resolve_mirna_aliases(character(0L))
  expect_equal(length(out), 0L)
  expect_type(out, "character")
})

test_that("resolve_mirna_aliases handles NA elements gracefully", {
  out <- resolve_mirna_aliases(c("hsa-miR-451a", NA_character_),
                               keep_unresolved = TRUE)
  expect_equal(out[1L], "hsa-miR-451a")
  expect_true(is.na(out[2L]))
})

test_that("resolve_mirna_aliases errors on missing required table columns", {
  bad_tbl <- data.frame(mirna_name = "hsa-miR-451a", mimat = "MIMAT0001631")
  expect_error(
    resolve_mirna_aliases("MIMAT0001631", table = bad_tbl),
    regexp = "missing columns"
  )
})
