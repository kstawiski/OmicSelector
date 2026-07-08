library(OmicSelector)

test_that("within-provenance blocks require size and direct control evidence", {
  d <- data.frame(
    y = c(rep("OC", 35), rep("CTRL", 35), rep("OC", 40), rep("CTRL", 10)),
    block = c(rep("documented", 70), rep("sparse", 50)),
    direct = c(rep(TRUE, 70), rep(TRUE, 40), rep(FALSE, 10))
  )
  blocks <- os_within_provenance_blocks(
    d,
    outcome = "y",
    provenance_block = "block",
    direct_control_evidence = "direct",
    positive = "OC"
  )
  expect_true(blocks$supports_positive_claim[blocks$provenance_block == "documented"])
  expect_false(blocks$supports_positive_claim[blocks$provenance_block == "sparse"])
})

test_that("identifiability gate fails closed when provenance floor is too high", {
  blocks <- data.frame(supports_positive_claim = TRUE)
  gate <- os_identifiability_gate(
    provenance_auc = 1.00,
    provenance_auc_ci_hi = 1.00,
    single_covariate_max_auc = 1.00,
    sensitivity_auc = 1.00,
    within_blocks = blocks,
    u1_u2_overlap_fraction = 0.10,
    u1_only_cases = 50,
    u1_only_controls = 50,
    physical_qc_pass = TRUE,
    anchor_leakage_pass = TRUE
  )
  expect_s3_class(gate, "os_identifiability_gate")
  expect_equal(gate$terminal_state, "FAIL-DEMOTED")
  expect_true("Provenance floor below ceiling" %in% gate$failed_criteria)
})

test_that("identifiability gate passes only when every criterion passes", {
  blocks <- data.frame(supports_positive_claim = TRUE)
  gate <- os_identifiability_gate(
    provenance_auc = 0.55,
    provenance_auc_ci_hi = 0.62,
    single_covariate_max_auc = 0.60,
    sensitivity_auc = 0.61,
    within_blocks = blocks,
    u1_u2_overlap_fraction = 0.20,
    u1_only_cases = 50,
    u1_only_controls = 50,
    physical_qc_pass = TRUE,
    anchor_leakage_pass = TRUE
  )
  expect_equal(gate$terminal_state, "PASS")
  expect_true(all(gate$criteria$status == "PASS"))
})

test_that("terminal gate ledger blocks non-pass statuses", {
  ledger <- os_terminal_gate_ledger(
    gate_id = c("source_qc", "identifiability", "null_benchmark"),
    status = c("PASS", "FAIL-DEMOTED", "NOT_APPLICABLE_WITH_CONSENSUS")
  )
  expect_s3_class(ledger, "os_terminal_gate_ledger")
  expect_equal(ledger$claim_permission[ledger$gate_id == "source_qc"], "eligible_if_downstream_gates_pass")
  expect_true(all(ledger$claim_permission[ledger$gate_id != "source_qc"] == "blocked_descriptive_only"))
})
