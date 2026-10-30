# =============================================================================
# A2 fix 2a (F6) — paired AUC-difference standard errors for the WITHIN-COHORT
# matched-null lift (method AUC minus the standalone rCLR baseline AUC).
#
# WHY ------------------------------------------------------------------------
# code/analyses/lift_threshold_sensitivity_v19.R currently estimates the
# per-cohort lift SE with a specified paired correlation (paired_lift_se, with
# RHO_PAIRED = 0.5). The method AUC and the rCLR baseline AUC are computed on the
# SAME cohort samples, so the correct paired SE uses the EMPIRICAL covariance
# between the two AUCs rather than an assumed rho. This module derives that
# covariance from the per-sample observed scores recovered by
# reproduce the gated reassay_design per-fold AUCs exactly).
#
# ESTIMAND (decisive for the SE construction) --------------------------------
# reassay_design's "Observed AUC" is the MEAN of per-fold AUCs, not a pooled out-of-fold
# AUC: matched_null_benchmark.R sets
#     auc_obs_cv <- mean(auc_obs_folds[valid_fold_mask]).
# Therefore the per-cohort lift is
#     lift = mean_f ( AUC_method,f - AUC_baseline,f )  over the K valid folds.
# The DeLong-primary SE is built to match THIS estimand: the per-fold DeLong
# variance of the paired AUC difference, aggregated across the K disjoint valid
# folds as the variance of their mean:
#     Var(lift) = (1 / K^2) * sum_f Var_DeLong( AUC_m,f - AUC_b,f ).
#
# ASSUMPTION made explicit for the code-review gate --------------------------
# The K outer TEST folds are disjoint (each sample is tested once), so the
# per-fold EVALUATION noise is independent and the variance-of-the-mean form is
# justified for the within-fold sampling component. It does NOT add a
# between-fold covariance term arising from shared TRAINING data across folds.
# For a single AUC that would understate the CV variance, but the lift is a
# DIFFERENCE of two AUCs trained and evaluated on the SAME folds, so the shared
# training component largely cancels in the difference and the dominant term is
# the within-fold method-vs-baseline covariance that DeLong captures exactly.
# The specified-rho = 0.5 paired SE is retained by the driver as a CONSERVATIVE
# cross-check column, and the driver's existing block-stratified cluster
# bootstrap (table_S74) provides the non-parametric companion at the pooled
# level. The substantive conclusion (no within-cohort method clears the +0.05
# lift boundary) must hold under both the DeLong and the rho = 0.5 SE; the SE
# choice changes CI width, not the near-zero lift point estimates.
#
# Reference: DeLong, DeLong & Clarke-Pearson (1988) Biometrics 44:837-845;
# Sun & Xu (2014) IEEE SP Lett 21:1389-1393 (fast placement values).
# Ties contribute 0.5 throughout, consistent with pROC::auc midrank handling.
# =============================================================================

# --- DeLong placement (structural) components for one ROC on one sample set --
# V10[i] = (1/nn) sum_j psi(pos_i, neg_j); V01[j] = (1/np) sum_i psi(pos_i,neg_j)
# with psi(a,b) = 1 if a>b, 0.5 if a==b, 0 if a<b. Fold test sets are small
# (tens of samples) so the exact O(np*nn) outer comparison is used.
#
# `orient` selects how the ROC is oriented to AUC >= 0.5, to MATCH the engine
# that produced the AUC being instrumented:
#   "median" - reflect when median(controls) > median(cases). This is exactly
#              pROC::auc(direction = "auto"), used by the within-cohort
#              matched-null engine (matched_null_benchmark.R).
#   "auc"    - reflect when AUC < 0.5 (i.e. report max(AUC, 1 - AUC)). This is
#              exactly ta_auc_auto (technology_aware_helpers.R), used by the
#              cross-cohort transfer engine.
#   "specified"  - do not inspect held-out outcomes to choose direction. Use this
#              when score direction was frozen from training data.
# The two rules disagree near the boundary, so the correct one must be passed by
# the caller. Reflection leaves the single-AUC DeLong variance unchanged
# (Var(1 - V) = Var(V)); it fixes the auc value to match the engine and aligns
# the paired covariance with the engine's oriented AUC difference.
.delong_placements <- function(y, score, orient = c("fixed", "median", "auc")) {
  orient <- match.arg(orient)
  pos <- score[y == 1L]; neg <- score[y == 0L]
  np <- length(pos); nn <- length(neg)
  if (np == 0L || nn == 0L) {
    return(list(V10 = numeric(0), V01 = numeric(0), np = np, nn = nn, auc = NA_real_))
  }
  cmp <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
  V10 <- rowMeans(cmp); V01 <- colMeans(cmp); auc <- mean(cmp)
  reflect <- if (orient == "fixed") FALSE else if (orient == "auc") {
    is.finite(auc) && auc < 0.5
  } else {
    stats::median(neg) > stats::median(pos)
  }
  if (reflect) {
    V10 <- 1 - V10; V01 <- 1 - V01; auc <- 1 - auc
  }
  list(V10 = V10, V01 = V01, np = np, nn = nn, auc = auc)
}

# --- single-AUC DeLong variance (exposed for the self-test vs pROC) ----------
delong_var_single <- function(y, score, orient = "fixed") {
  ok <- is.finite(score) & !is.na(y)            # drop non-finite before placements
  y <- y[ok]; score <- score[ok]
  p <- .delong_placements(y, score, orient = orient)
  if (p$np < 2L || p$nn < 2L) return(NA_real_)
  stats::var(p$V10) / p$np + stats::var(p$V01) / p$nn
}

# --- per-fold DeLong variance of the paired AUC difference (method-baseline) --
# Same samples, two scores. Returns the difference variance plus the two AUCs.
# Subsets to the PAIRED-COMPLETE finite set first: an NA in either score would
# otherwise propagate into median()/outer() and error, and the difference
# variance must be computed on the samples where both AUCs are defined.
.delong_var_diff_fold <- function(y, score_m, score_b, orient = "fixed") {
  ok <- is.finite(score_m) & is.finite(score_b) & !is.na(y)
  y <- y[ok]; score_m <- score_m[ok]; score_b <- score_b[ok]
  pm <- .delong_placements(y, score_m, orient = orient)
  pb <- .delong_placements(y, score_b, orient = orient)
  if (pm$np < 2L || pm$nn < 2L) {
    return(list(var = NA_real_, auc_m = pm$auc, auc_b = pb$auc, cov = NA_real_))
  }
  # 2x2 covariance matrices of the placement components (denominator n-1).
  S10 <- stats::cov(cbind(pm$V10, pb$V10))
  S01 <- stats::cov(cbind(pm$V01, pb$V01))
  np <- pm$np; nn <- pm$nn
  var_m  <- S10[1, 1] / np + S01[1, 1] / nn
  var_b  <- S10[2, 2] / np + S01[2, 2] / nn
  cov_mb <- S10[1, 2] / np + S01[1, 2] / nn
  list(var = max(0, var_m + var_b - 2 * cov_mb),
       auc_m = pm$auc, auc_b = pb$auc, cov = cov_mb)
}

# --- cohort lift SE: DeLong primary ------------------------------------------
# fold_dt: data.frame/data.table for ONE (cohort x method) with columns
#   fold (int), y (0/1), score_m (method), score_b (baseline rCLR).
# Returns lift (mean-of-folds AUC difference recomputed here as a sanity check),
# se (DeLong, estimand-consistent), and k (valid folds used).
delong_lift_se <- function(fold_dt, orient = "fixed") {
  folds <- sort(unique(fold_dt$fold))
  vars <- numeric(0); diffs <- numeric(0)
  for (f in folds) {
    sub <- fold_dt[fold_dt$fold == f, , drop = FALSE]
    if (length(unique(sub$y)) < 2L) next
    res <- .delong_var_diff_fold(sub$y, sub$score_m, sub$score_b, orient = orient)
    if (!is.finite(res$var)) next
    vars  <- c(vars, res$var)
    diffs <- c(diffs, res$auc_m - res$auc_b)
  }
  k <- length(vars)
  if (k == 0L) return(list(lift = NA_real_, se = NA_real_, k = 0L))
  list(lift = mean(diffs), se = sqrt(sum(vars) / k^2), k = k)
}
