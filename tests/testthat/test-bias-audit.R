library(OmicSelector)

# -----------------------------------------------------------------------------
# Setup: synthetic data that reproduces the 3D-Gene bias-floor pattern
# -----------------------------------------------------------------------------
make_provenance_data <- function(seed = 7L) {
  set.seed(seed)
  # 4 cohorts with very asymmetric case/control structure (mimics GSE106817
  # OC-only, GSE211692 mostly-healthy, small OC-contributors)
  cohorts <- c(
    rep("A", 317), rep("B", 3321), rep("C", 63), rep("D", 254)
  )
  y <- c(
    rep(1, 317),                      # A: 317 cases / 0 healthy
    c(rep(1, 397), rep(0, 2924)),     # B: 397 / 2924
    c(rep(1, 15),  rep(0, 48)),       # C: 15 / 48
    c(rep(1, 25),  rep(0, 229))       # D: 25 / 229
  )
  cohorts <- factor(cohorts)

  # 5 features: 2 "real biology" (independent of cohort), 2 "batch-proxy"
  # (carry strong cohort signal even within healthy samples), 1 "pure noise".
  n <- length(y)
  real_bio <- cbind(
    stats::rnorm(n, mean = ifelse(y == 1, 1.0, 0.0), sd = 1.0),
    stats::rnorm(n, mean = ifelse(y == 1, 0.5, 0.0), sd = 1.2)
  )
  # Batch proxy: shifted by cohort, and separately by y (hemolysis-like
  # combined signal). Aligns with the 3D-Gene empirical pattern.
  cohort_shift <- c(A = 3.0, B = 0.0, C = 1.2, D = -0.5)[as.character(cohorts)]
  batch_proxy <- cbind(
    stats::rnorm(n, mean = cohort_shift + ifelse(y == 1, 0.6, 0.0), sd = 0.8),
    stats::rnorm(n, mean = cohort_shift + ifelse(y == 1, 0.4, 0.0), sd = 0.9)
  )
  noise <- stats::rnorm(n, sd = 1.0)
  X <- cbind(real_bio, batch_proxy, noise)
  colnames(X) <- c("bio1", "bio2", "batch1", "batch2", "noise")
  list(X = X, y = y, cohort = cohorts)
}


# -----------------------------------------------------------------------------
# os_bias_floor_auc
# -----------------------------------------------------------------------------
test_that("bias floor detects cohort-identity signal in asymmetric 4-cohort design", {
  d <- make_provenance_data()
  bf <- os_bias_floor_auc(d$y, d$cohort, cv_folds = 5L, seed = 1L)
  expect_s3_class(bf, "os_bias_floor")
  # With A=OC-only and B=mostly-Healthy, apparent AUC must be substantial
  expect_gt(bf$apparent_auc, 0.68)
  expect_false(is.na(bf$cv_auc))
  expect_gt(bf$cv_auc, 0.65)
  expect_true(grepl("HIGH|MODERATE", bf$interpretation))
  expect_equal(nrow(bf$per_cohort), 4L)
})

test_that("bias floor is near 0.5 when cohort is independent of case status", {
  set.seed(3L)
  n <- 1000L
  cohort <- sample(letters[1:3], n, replace = TRUE)
  # case status assigned independently of cohort
  y <- stats::rbinom(n, 1L, 0.4)
  bf <- os_bias_floor_auc(y, cohort, cv_folds = 5L, seed = 1L)
  expect_lt(abs(bf$apparent_auc - 0.5), 0.08)
  expect_true(grepl("LOW", bf$interpretation))
})

test_that("bias floor errors when y has only one level", {
  expect_error(
    os_bias_floor_auc(rep(1, 20), rep(letters[1:2], 10)),
    "two levels"
  )
})


# -----------------------------------------------------------------------------
# os_covariate_only_auc
# -----------------------------------------------------------------------------
test_that("covariate-only AUC recovers age-confounding signal", {
  set.seed(5L)
  n <- 1000L
  age <- c(stats::rnorm(400, 57, 8), stats::rnorm(600, 45, 7))
  y <- c(rep(1, 400), rep(0, 600))  # cases older than controls
  out <- os_covariate_only_auc(y, age, covariate_name = "age")
  expect_gt(out$auc, 0.80)
  expect_equal(out$direction, "higher values = cases")
  expect_gt(out$effect_size, 1.0) # large Cohen's d
})

test_that("covariate-only AUC handles factor covariate with Cramer's V", {
  set.seed(7L)
  sex <- factor(c(rep("M", 150), rep("F", 250), rep("M", 50), rep("F", 150)))
  y   <- c(rep(1, 400), rep(0, 200))
  out <- os_covariate_only_auc(y, sex, covariate_name = "sex")
  expect_true(out$auc >= 0.5)
  expect_equal(out$effect_metric, "Cramer's V")
})


# -----------------------------------------------------------------------------
# os_per_feature_batch_signal
# -----------------------------------------------------------------------------
test_that("per-feature batch signal flags batch proxy miRNAs over real biology", {
  d <- make_provenance_data()
  fs <- os_per_feature_batch_signal(d$X, d$y, d$cohort,
                                    feature_names = colnames(d$X))
  expect_equal(sort(fs$feature), sort(colnames(d$X)))
  # batch1/batch2 must rank above bio1/bio2 on dataset_case_corr
  top2 <- head(fs$feature, 2L)
  expect_true(all(top2 %in% c("batch1", "batch2")))
  # F_dataset_healthy for batch features should exceed bio features
  f_batch <- fs$F_dataset_healthy[fs$feature %in% c("batch1", "batch2")]
  f_bio   <- fs$F_dataset_healthy[fs$feature %in% c("bio1", "bio2")]
  expect_true(min(f_batch, na.rm = TRUE) > max(f_bio, na.rm = TRUE))
})


# -----------------------------------------------------------------------------
# os_clustered_bootstrap_auc
# -----------------------------------------------------------------------------
test_that("clustered bootstrap CI widens vs naive bootstrap under duplication", {
  set.seed(11L)
  n_unique <- 200L
  # Each "true specimen" duplicated 3x
  y_u <- stats::rbinom(n_unique, 1L, 0.4)
  score_u <- y_u + stats::rnorm(n_unique, 0, 0.6)
  cluster_id <- rep(seq_len(n_unique), each = 3)
  y <- rep(y_u, each = 3)
  score <- rep(score_u, each = 3)

  cl <- os_clustered_bootstrap_auc(y, score, cluster_id,
                                   n_boot = 400L, seed = 1L)
  naive <- os_clustered_bootstrap_auc(y, score, seq_along(y),
                                      n_boot = 400L, seed = 1L)
  cl_width <- cl$ci_upper - cl$ci_lower
  naive_width <- naive$ci_upper - naive$ci_lower
  expect_gt(cl_width, naive_width * 1.3)
})

test_that("one-sided-lower bootstrap returns NA upper bound", {
  set.seed(13L)
  n <- 200L
  y <- stats::rbinom(n, 1L, 0.5)
  score <- y + stats::rnorm(n, 0, 0.5)
  out <- os_clustered_bootstrap_auc(y, score, seq_along(y),
                                    n_boot = 200L, seed = 1L,
                                    ci_type = "one-sided-lower")
  expect_true(is.na(out$ci_upper))
  expect_false(is.na(out$ci_lower))
})


# -----------------------------------------------------------------------------
# os_detect_cross_cohort_duplicates
# -----------------------------------------------------------------------------
test_that("cross-cohort duplicates detected by exact feature equality", {
  set.seed(19L)
  n <- 40L
  X <- matrix(stats::rnorm(n * 5), n, 5L)
  colnames(X) <- paste0("F", 1:5)
  # Make sample 3 (cohort A) identical to sample 22 (cohort B)
  X[22, ] <- X[3, ]
  cohort <- c(rep("A", 20), rep("B", 20))
  sample_id <- sprintf("S%03d", seq_len(n))
  dup <- os_detect_cross_cohort_duplicates(X, cohort, sample_id = sample_id)
  expect_gte(nrow(dup), 1L)
  expect_true(any(dup$sample_id_a == "S003" | dup$sample_id_b == "S003"))
})

test_that("cross-cohort duplicate tolerance and feature support are enforced", {
  X <- rbind(
    c(1, 2, 3, 4),
    c(1 + 5e-5, 2 - 5e-5, 99, NA),
    c(1, 2, 3, 4)
  )
  colnames(X) <- paste0("F", 1:4)
  cohort <- c("A", "B", "B")
  sample_id <- c("A1", "B1", "B2")

  exact <- os_detect_cross_cohort_duplicates(X, cohort, sample_id)
  expect_equal(nrow(exact), 1L)
  expect_setequal(c(exact$sample_id_a, exact$sample_id_b), c("A1", "B2"))

  tolerant <- os_detect_cross_cohort_duplicates(
    X[, 1:2], cohort, sample_id, tolerance = 1e-4, min_features = 2L
  )
  expect_equal(nrow(tolerant), 2L)
  expect_true("tolerance_feature_equality" %in% tolerant$match_type)

  insufficient <- os_detect_cross_cohort_duplicates(
    X[1:2, ], cohort[1:2], sample_id[1:2],
    tolerance = 100, min_features = 4L
  )
  expect_equal(nrow(insufficient), 0L)

  expect_error(
    os_detect_cross_cohort_duplicates(X, cohort, sample_id, tolerance = -1),
    "non-negative"
  )
  expect_error(
    os_detect_cross_cohort_duplicates(X, cohort, sample_id, min_features = 5L),
    "1 through ncol"
  )
})

test_that("cross-cohort duplicate matching requires unique sample IDs", {
  X <- matrix(1:8, nrow = 2L)
  expect_error(
    os_detect_cross_cohort_duplicates(X, c("A", "B"), c("same", "same")),
    "sample_id must be unique"
  )
})

test_that("feature-only matches do not borrow a discordant specimen ID", {
  X <- rbind(c(1, 2, 3), c(1, 2, 3))
  out <- os_detect_cross_cohort_duplicates(
    X, c("A", "B"), c("A1", "B1"), specimen_id = c("SP1", "SP2")
  )
  expect_equal(nrow(out), 1L)
  expect_true(is.na(out$specimen_id))
  expect_identical(out$match_type, "exact_feature_equality")
})

test_that("cross-cohort duplicates via explicit specimen_id crosswalk", {
  X <- matrix(stats::rnorm(15 * 3), 15L, 3L)
  cohort <- c(rep("A", 10), rep("B", 5))
  sample_id <- sprintf("S%02d", seq_len(15L))
  specimen_id <- c(
    sprintf("OV%04d", 1:10),
    c("OV0003", "OV0004", "OV0005", "OV9997", "OV9998")  # 3 overlaps w/ A
  )
  dup <- os_detect_cross_cohort_duplicates(X, cohort,
                                           sample_id = sample_id,
                                           specimen_id = specimen_id)
  crosswalk <- dup[dup$match_type == "specimen_id_crosswalk", , drop = FALSE]
  expect_equal(nrow(crosswalk), 3L)
})


# -----------------------------------------------------------------------------
# os_bias_audit (orchestrator)
# -----------------------------------------------------------------------------
test_that("os_bias_audit returns structured report with all components", {
  d <- make_provenance_data()
  # fake model: a classifier that rides mostly on the batch proxies
  fake_preds <- 0.7 * d$X[, "batch1"] + 0.3 * d$X[, "bio1"]
  age <- ifelse(d$y == 1, 57 + stats::rnorm(length(d$y), 0, 4),
                45 + stats::rnorm(length(d$y), 0, 3))
  audit <- os_bias_audit(
    X = d$X, y = d$y, cohort = d$cohort,
    predictions = fake_preds,
    covariates = list(age = age),
    feature_names = colnames(d$X)
  )
  expect_s3_class(audit, "os_bias_audit")
  expect_named(audit, c(
    "bias_floor", "covariate_aucs", "feature_signal", "duplicates",
    "provenance_summary", "classifier_vs_floor", "n", "n_cohorts",
    "audit_version", "audit_date"
  ))
  expect_s3_class(audit$bias_floor, "os_bias_floor")
  expect_true(!is.null(audit$classifier_vs_floor))
  expect_gt(audit$classifier_vs_floor$classifier_auc, audit$classifier_vs_floor$bias_floor_auc)
  expect_true(grepl("CONCERN|MODERATE|ACCEPTABLE", audit$classifier_vs_floor$interpretation))
})

test_that("os_bias_audit prints without error", {
  d <- make_provenance_data()
  audit <- os_bias_audit(X = d$X, y = d$y, cohort = d$cohort)
  expect_output(print(audit), "OmicSelector Cohort-Provenance Bias Audit")
  expect_output(print(audit), "Bias floor")
  expect_output(print(audit), "Top-5 features")
})


# -----------------------------------------------------------------------------
# os_bias_audit_report (plain-text convenience wrapper)
# -----------------------------------------------------------------------------
test_that("os_bias_audit_report returns a character string with the report", {
  d <- make_provenance_data()
  txt <- os_bias_audit_report(X = d$X, y = d$y, cohort = d$cohort)
  expect_s3_class(txt, "os_bias_audit_report")
  expect_type(txt, "character")
  expect_length(txt, 1L)
  expect_match(txt, "OmicSelector Cohort-Provenance Bias Audit", fixed = TRUE)
  expect_match(txt, "Bias floor", fixed = TRUE)
  expect_s3_class(attr(txt, "audit"), "os_bias_audit")
})

test_that("os_bias_audit_report accepts a pre-computed audit object", {
  d <- make_provenance_data()
  audit <- os_bias_audit(X = d$X, y = d$y, cohort = d$cohort)
  txt <- os_bias_audit_report(audit = audit)
  expect_s3_class(txt, "os_bias_audit_report")
  expect_identical(attr(txt, "audit"), audit)
})

test_that("os_bias_audit_report validates its inputs", {
  expect_error(
    os_bias_audit_report(),
    "Provide either `audit` or all of `X`, `y`, `cohort`"
  )
  expect_error(
    os_bias_audit_report(audit = list(foo = 1)),
    "must be an `os_bias_audit` object"
  )
})
