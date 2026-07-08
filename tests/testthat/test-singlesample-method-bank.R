test_that("Paper 3 method bank records exported package functions", {
  bank <- singlesample_method_bank()

  expect_true(nrow(bank) >= 33L)
  expect_true(all(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J") %in%
                    bank$family))
  expect_true(all(bank$functions_exist))
  expect_true(all(bank$functions_exported))
  expect_silent(singlesample_assert_method_bank_exports(bank))
})

.singlesample_family_fixture <- function(seed = 37L, n = 96L, p_extra = 52L) {
  set.seed(seed)
  base_features <- c(
    "hsa-miR-451a", "hsa-miR-23a-3p", "hsa-miR-103a-3p",
    "hsa-miR-191-5p", "hsa-miR-26a-5p", "hsa-miR-30c-5p",
    "hsa-let-7g-5p", "hsa-miR-93-5p", "hsa-miR-486-5p",
    "hsa-miR-16-5p", "hsa-miR-92a-3p", "hsa-miR-144-3p",
    "hsa-miR-223-3p"
  )
  generic <- paste0("hsa-miR-test-", seq_len(p_extra))
  features <- c(base_features, generic)
  y <- rep(c(0L, 1L), length.out = n)
  X <- matrix(rlnorm(n * length(features), meanlog = 4.5, sdlog = 0.35),
              nrow = n, ncol = length(features),
              dimnames = list(paste0("S", seq_len(n)), features))
  signal <- generic[seq_len(10L)]
  X[y == 1L, signal] <- X[y == 1L, signal] * 2.0
  kit <- rep(c("kitA", "kitB", "kitC", "kitD"), length.out = n)
  bio <- rep(c("serum-clot", "plasma-EDTA", "plasma-other", "exosome"),
             length.out = n)
  tech <- rep(c("NGS", "Toray", "qPCR", "NanoString"), length.out = n)
  cohort <- rep(paste0("C", seq_len(8L)), length.out = n)
  block <- ifelse(cohort %in% c("C1", "C2"), "B12", cohort)
  cancer <- rep(c("breast", "ovarian", "lung", "colorectal"),
                length.out = n)
  X[kit == "kitB", generic[11:18]] <- X[kit == "kitB", generic[11:18]] * 1.7
  meta <- data.frame(
    sample_id = rownames(X),
    disease = y,
    cohort = cohort,
    provenance_block = block,
    cancer_type = cancer,
    technology = tech,
    kit_family = kit,
    kit_label = kit,
    library_kit = kit,
    biofluid = bio,
    RIN = rep(c(NA, 6.5, 7.2, 8.1), length.out = n),
    stringsAsFactors = FALSE
  )
  list(X = X, y = y, meta = meta, panel = signal[seq_len(8L)])
}

test_that("kit, provenance, technology, biofluid, and anchor scorers are callable", {
  fx <- .singlesample_family_fixture()

  kit_fit <- fit_kit_residual_mad(fx$X, fx$meta, kit_label_col = "kit_family")
  expect_true(all(is.finite(score_kit_residual_mad(
    fx$X, fx$meta, fx$panel, kit_label_col = "kit_family", fit = kit_fit))))
  expect_true(all(is.finite(score_kit_stratified_rclr(
    fx$X, fx$meta, fx$panel, kit_label_col = "kit_family"))))

  cohort_fit <- fit_cohort_z_then_pool_score(fx$X, fx$meta, fx$panel)
  block_fit <- fit_block_fe_rclr(fx$X, fx$meta, fx$panel)
  mixed_fit <- fit_mixed_effects_scorer(fx$X, fx$meta, fx$panel)
  expect_true(all(is.finite(predict_cohort_z_then_pool_score(
    cohort_fit, fx$X, fx$meta))))
  expect_true(all(is.finite(predict_block_fe_rclr(block_fit, fx$X, fx$meta))))
  expect_true(all(is.finite(predict_mixed_effects_scorer(
    mixed_fit, fx$X, fx$meta))))
  expect_true(all(is.finite(score_cluster_robust_ensemble(
    fx$X, fx$meta, fx$panel))))

  tech_fit <- fit_tech_stratified_rclr(fx$X, fx$meta, fx$panel)
  cross_fit <- fit_cross_tech_harmonized_rclr(
    fx$X, fx$meta, fx$panel,
    required_technologies = c("NGS", "Toray", "qPCR", "NanoString"),
    min_detection_rate = 0.01, n_anchor = 10L)
  tech_alr <- fit_tech_residualized_alr(
    fx$X, fx$meta, fx$panel, outcome_col = "disease")
  expect_true(all(is.finite(predict_tech_stratified_rclr(
    tech_fit, fx$X, fx$meta))))
  expect_true(all(is.finite(predict_cross_tech_harmonized_rclr(
    cross_fit, fx$X, fx$meta, min_anchors_present = 3L))))
  expect_true(all(is.finite(predict_tech_residualized_alr(
    tech_alr, fx$X, fx$meta))))

  bio_fit <- fit_biofluid_stratified_rclr(fx$X, fx$meta, fx$panel)
  bio_anchor <- fit_biofluid_anchor_rclr(fx$X, fx$meta, fx$panel)
  bio_alr <- fit_biofluid_residualized_alr(
    fx$X, fx$meta, fx$panel, outcome_col = "disease")
  expect_true(all(is.finite(predict_biofluid_stratified_rclr(
    bio_fit, fx$X, fx$meta))))
  expect_true(all(is.finite(predict_biofluid_anchor_rclr(
    bio_anchor, fx$X, fx$meta))))
  expect_true(all(is.finite(predict_biofluid_residualized_alr(
    bio_alr, fx$X, fx$meta))))
  expect_true(all(is.finite(score_platelet_exclusion_rclr(
    fx$X, fx$meta, fx$panel))))

  anchor_fit <- fit_kit_stable_anchors(
    fx$X, fx$meta, kit_col = "kit_label", min_anchors = 5L)
  dual_fit <- fit_hemolysis_kit_dual_anchor(
    fx$X, fx$meta, kit_col = "kit_label", min_anchors = 5L)
  rin_fit <- fit_rin_weighted_reference(fx$X, fx$meta, rin_col = "RIN")
  expect_true(all(is.finite(score_kit_stable_anchor_rclr(
    fx$X, fx$meta, fx$panel, anchor_fit = anchor_fit,
    min_eval_anchors = 3L))))
  expect_true(all(is.finite(score_hemolysis_kit_dual_anchor(
    fx$X, fx$meta, fx$panel, dual_fit = dual_fit,
    min_eval_anchors = 3L))))
  expect_true(all(is.finite(score_rin_weighted_score(
    fx$X, fx$meta, fx$panel, rin_fit = rin_fit))))
})

test_that("learned, Group-DRO, Sinkhorn, and external-validation helpers work", {
  fx <- .singlesample_family_fixture(n = 72L, p_extra = 36L)

  dann <- os_fit_dann_kit_extended(
    fx$X, train_disease_labels = fx$y,
    train_kit_labels = fx$meta$kit_family,
    train_cohort_labels = fx$meta$cohort,
    panel_features = fx$panel,
    use_torch = FALSE)
  expect_true(all(is.finite(os_score_dann_kit_extended(dann, fx$X))))

  vae <- os_fit_kit_conditional_vae(
    fx$X, train_disease_labels = fx$y,
    train_kit_labels = fx$meta$kit_family,
    train_cohort_labels = fx$meta$cohort,
    panel_features = fx$panel,
    use_torch = FALSE)
  expect_true(all(is.finite(os_score_kit_conditional_vae(vae, fx$X))))

  icp <- suppressWarnings(os_fit_icp_per_kit(
    fx$X, train_disease_labels = fx$y,
    train_kit_labels = fx$meta$kit_family,
    train_cohort_labels = fx$meta$cohort,
    panel_features = fx$panel,
    min_per_kit = 12L))
  expect_true(all(is.finite(os_score_icp_per_kit(icp, fx$X))))

  gdro <- os_fit_group_dro_scorer(
    fx$X, fx$y, fx$meta, kit_col = "kit_family",
    biofluid_col = "biofluid", panel_features = fx$panel,
    tune = FALSE, epochs = 25L, seed = 11L)
  gdro_score <- os_score_group_dro_scorer(gdro, fx$X, fx$meta)
  expect_true(all(is.finite(gdro_score)))
  expect_true(all(gdro_score >= 0 & gdro_score <= 1))

  sink <- os_fit_sinkhorn_ot_scorer(
    fx$X, fx$y, fx$meta, cohort_col = "cohort",
    panel_features = fx$panel, tune = FALSE,
    epsilon = 0.10, lambda = 0.5,
    n_barycenter_atoms = 8L, max_cohort_atoms = 6L,
    barycenter_iter = 1L, max_sinkhorn_iter = 40L,
    backend = "R", seed = 13L)
  expect_true(all(is.finite(os_score_sinkhorn_ot_scorer(fx$X[1:6, ], sink))))

  loco <- singlesample_make_loco_splits(fx$meta, min_train_cohorts = 1L)
  expect_equal(length(loco), length(unique(fx$meta$cohort)))
  c1 <- loco[[match("C1", vapply(loco, `[[`, character(1), "held_out"))]]
  expect_false(any(fx$meta$provenance_block[c1$train] == "B12"))

  locto <- singlesample_make_locto_splits(fx$meta, min_train_cohorts = 1L)
  expect_equal(length(locto), length(unique(fx$meta$cancer_type)))
})
