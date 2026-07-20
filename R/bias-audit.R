#' @title Cohort-Provenance Bias Audit
#'
#' @description
#' Diagnostic tools for detecting and quantifying non-biological signal in
#' biomarker classifier evaluation. These are generic miRNA/omics QC checks
#' that every biomarker study should pass before any AUC claim. They expose
#' confounding from:
#' \itemize{
#'   \item cohort identity (institution, batch, dataset),
#'   \item demographic covariates (age, sex, BMI),
#'   \item pre-analytical asymmetry (storage, collection site, tube type),
#'   \item specimen duplication across accessions.
#' }
#'
#' @details
#' The motivating problem is simple: on the four public 3D-Gene ovarian-cancer
#' cohorts (GSE106817, GSE211692, GSE113486, GSE113740), a logistic regression
#' using \emph{only} dataset-identity dummy variables reaches AUC 0.72 against
#' cancer/healthy labels without touching a single miRNA value. Any classifier
#' trained on the same cohort structure therefore inherits this "bias floor"
#' even if its apparent performance is higher.
#'
#' The functions here make that floor visible before any biological claim is
#' made, and extend to:
#' \itemize{
#'   \item covariate-only AUCs (age, sex, etc.),
#'   \item clustered bootstrap CIs that respect specimen duplication,
#'   \item feature-level batch-signal ANOVA within healthy samples,
#'   \item case/control pre-analytical asymmetry tables.
#' }
#'
#' @references
#' Leek JT, Scharpf RB, Bravo HC, et al. (2010). Tackling the widespread and
#' critical impact of batch effects in high-throughput data.
#' \emph{Nature Reviews Genetics}, 11(10), 733-739.
#'
#' Pritchard CC, Kroh E, Wood B, et al. (2012). Blood cell origin of circulating
#' microRNAs: a cautionary note for cancer biomarker studies.
#' \emph{Cancer Prevention Research}, 5(3), 492-497.
#'
#' Kirschner MB, Kao SC, Edelman JJ, et al. (2011). Haemolysis during sample
#' preparation alters microRNA content of plasma. \emph{PLoS ONE}, 6(9), e24145.
#'
#' @name bias-audit
NULL


#' @title Dataset-Identity AUC ("Bias Floor")
#'
#' @description
#' Fits a logistic regression using ONLY cohort/dataset dummy variables as
#' predictors and reports its AUC against the case/control label. This is
#' the empirical lower bound on how much of an apparent classifier AUC can
#' be explained by cohort identity alone, before any biomarker is used.
#'
#' If the bias floor is high (e.g., > 0.70), any pooled classifier AUC on
#' the same data inherits that share of discrimination from cohort provenance
#' rather than biology. A responsible biomarker claim must show AUC
#' \emph{materially} above the bias floor, on data where the floor itself
#' has been mitigated.
#'
#' @param y Binary outcome vector (0/1 or factor with two levels). Length n.
#' @param cohort Vector or factor of cohort/dataset identity. Length n.
#' @param cv_folds Integer. If non-NULL, also computes K-fold cross-validated
#'   AUC. Default 5. Set to NULL to skip CV.
#' @param seed Random seed for CV folds. Default 42L.
#'
#' @return A list with class \code{"os_bias_floor"}:
#'   \item{apparent_auc}{Logistic regression AUC on the training data
#'     (upper bound — optimistic by design)}
#'   \item{cv_auc}{Cross-validated AUC, if requested}
#'   \item{per_cohort}{Data frame: cohort, n, case rate, mean predicted probability}
#'   \item{interpretation}{Character string summarizing the floor}
#'
#' @examples
#' \dontrun{
#' y <- c(rep(1, 317), rep(0, 200), rep(1, 40), rep(0, 100))
#' cohort <- c(rep("A", 317), rep("A", 200), rep("B", 40), rep("B", 100))
#' bf <- os_bias_floor_auc(y, cohort)
#' print(bf)
#' }
#'
#' @export
os_bias_floor_auc <- function(y, cohort, cv_folds = 5L, seed = 42L) {
  if (length(y) != length(cohort)) {
    stop("`y` and `cohort` must have the same length.", call. = FALSE)
  }
  y_num <- as.numeric(as.factor(y)) - 1L
  if (any(is.na(y_num))) {
    keep <- !is.na(y_num) & !is.na(cohort)
    y_num <- y_num[keep]
    cohort <- cohort[keep]
  }
  cohort <- factor(cohort)

  if (length(unique(y_num)) < 2L) {
    stop("`y` must have two levels.", call. = FALSE)
  }

  df <- data.frame(y = y_num, cohort = cohort)

  fit <- suppressWarnings(
    stats::glm(y ~ cohort, data = df, family = stats::binomial())
  )
  probs <- stats::predict(fit, type = "response")
  apparent_auc <- .auc_fast(y_num, probs)

  cv_auc <- NA_real_
  if (!is.null(cv_folds) && cv_folds > 1L) {
    set.seed(seed)
    idx <- sample(seq_len(nrow(df)))
    fold <- rep(seq_len(cv_folds), length.out = nrow(df))[idx]
    preds <- rep(NA_real_, nrow(df))
    for (k in seq_len(cv_folds)) {
      tr <- fold != k
      te <- fold == k
      if (length(unique(df$y[tr])) < 2L) next
      fit_k <- suppressWarnings(
        stats::glm(y ~ cohort, data = df[tr, ], family = stats::binomial())
      )
      preds[te] <- suppressWarnings(stats::predict(fit_k, newdata = df[te, ], type = "response"))
    }
    ok <- !is.na(preds)
    if (length(unique(df$y[ok])) == 2L) {
      cv_auc <- .auc_fast(df$y[ok], preds[ok])
    }
  }

  per_cohort <- aggregate(
    df$y, by = list(cohort = df$cohort),
    FUN = function(v) c(n = length(v), case_rate = mean(v))
  )
  per_cohort <- data.frame(
    cohort = per_cohort$cohort,
    n = per_cohort$x[, "n"],
    case_rate = per_cohort$x[, "case_rate"]
  )

  verdict <- if (apparent_auc > 0.70) {
    "HIGH bias floor: cohort identity alone achieves substantial discrimination; biomarker claims must show AUC materially above this floor on cohort-balanced data."
  } else if (apparent_auc > 0.60) {
    "MODERATE bias floor: cohort composition is informative. Report bias floor alongside classifier AUC."
  } else {
    "LOW bias floor: cohort composition alone cannot discriminate; biomarker AUC is not obviously confounded by dataset identity."
  }

  out <- list(
    apparent_auc = apparent_auc,
    cv_auc = cv_auc,
    per_cohort = per_cohort,
    interpretation = verdict,
    n = length(y_num),
    n_cohorts = nlevels(cohort)
  )
  class(out) <- c("os_bias_floor", "list")
  out
}


#' @title Covariate-Only AUC
#'
#' @description
#' Fits a logistic regression using a single covariate (typically age or sex)
#' as the only predictor of case status. Quantifies how much of the apparent
#' classifier signal could be explained by demographic confounding alone.
#'
#' @param y Binary outcome vector.
#' @param covariate Numeric or factor covariate (e.g., age in years).
#' @param covariate_name Human-readable label for reporting.
#'
#' @return A list:
#'   \item{auc}{AUC of covariate-only classifier}
#'   \item{direction}{"higher values = cases" or "higher values = controls"}
#'   \item{effect_size}{Cohen's d (numeric) or Cramer's V (factor)}
#'
#' @export
os_covariate_only_auc <- function(y, covariate, covariate_name = "covariate") {
  y_num <- as.numeric(as.factor(y)) - 1L
  keep <- !is.na(y_num) & !is.na(covariate)
  y_num <- y_num[keep]
  covariate <- covariate[keep]

  if (length(unique(y_num)) < 2L) {
    stop("`y` must have two levels after dropping NA.", call. = FALSE)
  }

  if (is.numeric(covariate)) {
    df <- data.frame(y = y_num, x = covariate)
    fit <- suppressWarnings(stats::glm(y ~ x, data = df, family = stats::binomial()))
    probs <- stats::predict(fit, type = "response")
    auc <- .auc_fast(y_num, probs)
    d <- (mean(covariate[y_num == 1L]) - mean(covariate[y_num == 0L])) /
      sqrt((stats::var(covariate[y_num == 1L]) + stats::var(covariate[y_num == 0L])) / 2)
    direction <- if (mean(covariate[y_num == 1L]) > mean(covariate[y_num == 0L])) {
      "higher values = cases"
    } else {
      "higher values = controls"
    }
    list(
      covariate = covariate_name,
      auc = auc,
      direction = direction,
      effect_size = unname(d),
      effect_metric = "Cohen's d"
    )
  } else {
    covariate <- factor(covariate)
    df <- data.frame(y = y_num, x = covariate)
    fit <- suppressWarnings(stats::glm(y ~ x, data = df, family = stats::binomial()))
    probs <- stats::predict(fit, type = "response")
    auc <- .auc_fast(y_num, probs)
    tab <- table(covariate, y_num)
    chi <- suppressWarnings(stats::chisq.test(tab))
    v <- sqrt(chi$statistic / (sum(tab) * (min(dim(tab)) - 1L)))
    list(
      covariate = covariate_name,
      auc = auc,
      direction = "categorical",
      effect_size = unname(as.numeric(v)),
      effect_metric = "Cramer's V"
    )
  }
}


#' @title Per-Feature Batch-vs-Case Signal Partitioning
#'
#' @description
#' For each feature, computes three quantities that jointly distinguish
#' biological case signal from dataset/batch provenance:
#' \itemize{
#'   \item \code{F_dataset_healthy}: one-way ANOVA of feature expression
#'     against dataset identity, \strong{using only healthy samples}
#'     (so case biology is removed);
#'   \item \code{F_case}: one-way ANOVA against case/control status,
#'     using all samples;
#'   \item \code{dataset_case_corr}: Pearson correlation between the
#'     sample-wise dataset-mean signature and the sample-wise case-mean
#'     signature.
#' }
#' Features with high \code{F_dataset_healthy} AND high \code{dataset_case_corr}
#' are the most suspect: their case signal is colinear with dataset identity.
#'
#' @param X Numeric matrix or data frame (n samples x p features).
#' @param y Binary case vector (0/1). Length n.
#' @param cohort Factor or vector of cohort/dataset identity. Length n.
#' @param feature_names Optional character vector of feature labels.
#'
#' @return A data frame with one row per feature, including a
#'   \code{suspect_score} that combines healthy-only dataset ANOVA magnitude
#'   with \code{dataset_case_corr}, ranked from most to least provenance-like.
#'
#' @export
os_per_feature_batch_signal <- function(X, y, cohort, feature_names = NULL) {
  X <- as.matrix(X)
  if (is.null(feature_names)) {
    feature_names <- colnames(X)
    if (is.null(feature_names)) feature_names <- sprintf("F%04d", seq_len(ncol(X)))
  }
  y_num <- as.numeric(as.factor(y)) - 1L
  cohort <- factor(cohort)

  healthy_idx <- which(y_num == 0L)
  dataset_mean <- tapply(seq_len(nrow(X)), cohort, mean) # dataset mean (for corr)

  # Per-sample dataset-mean signature and case-mean signature
  per_sample_dataset_score <- numeric(nrow(X))
  per_sample_case_score <- numeric(nrow(X))

  results <- data.frame(
    feature = feature_names,
    F_dataset_healthy = NA_real_,
    p_dataset_healthy = NA_real_,
    F_case = NA_real_,
    p_case = NA_real_,
    dataset_case_corr = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    # Healthy-only dataset ANOVA
    if (length(unique(cohort[healthy_idx])) >= 2L) {
      a1 <- stats::aov(x[healthy_idx] ~ cohort[healthy_idx])
      s1 <- summary(a1)[[1L]]
      results$F_dataset_healthy[j] <- s1[["F value"]][1L]
      results$p_dataset_healthy[j] <- s1[["Pr(>F)"]][1L]
    }
    # Case ANOVA (all samples)
    a2 <- stats::aov(x ~ factor(y_num))
    s2 <- summary(a2)[[1L]]
    results$F_case[j] <- s2[["F value"]][1L]
    results$p_case[j] <- s2[["Pr(>F)"]][1L]

    # Per-sample dataset-mean score (mean of x within each sample's dataset)
    ds_mean <- tapply(x, cohort, mean)
    case_mean <- tapply(x, factor(y_num), mean)
    per_sample_dataset_score <- ds_mean[as.character(cohort)]
    per_sample_case_score <- case_mean[as.character(y_num)]
    results$dataset_case_corr[j] <- suppressWarnings(
      stats::cor(per_sample_dataset_score, per_sample_case_score)
    )
  }

  # Rank by joint provenance-like signal, not correlation alone. A high
  # correlation is only suspicious when accompanied by substantial
  # healthy-only dataset structure.
  results$suspect_score <- with(
    results,
    log1p(pmax(F_dataset_healthy, 0)) * abs(dataset_case_corr)
  )

  results <- results[
    order(-results$suspect_score, -results$dataset_case_corr, -results$F_dataset_healthy),
    ,
    drop = FALSE
  ]
  rownames(results) <- NULL
  results
}


#' @title Clustered Bootstrap CI for AUC
#'
#' @description
#' Bootstrap the AUC of a set of predictions respecting within-cluster
#' dependence (e.g., duplicated specimens with the same internal ID across
#' accessions; repeated measures from the same patient). Rows with the same
#' cluster ID are resampled as a unit, which preserves the within-cluster
#' correlation structure and prevents anti-conservative CIs.
#'
#' @param y Binary outcome vector.
#' @param predictions Numeric predicted probabilities / scores.
#' @param cluster_id Vector of cluster identifiers (length n).
#' @param n_boot Number of bootstrap replicates. Default 2000.
#' @param stratify_by Optional factor for stratified resampling
#'   (e.g., fold x class). Strata are preserved at the cluster level.
#' @param ci_type One of \code{"two-sided"} or \code{"one-sided-lower"}.
#'   Default \code{"two-sided"}.
#' @param alpha Significance level (default 0.05).
#' @param seed Random seed.
#'
#' @return A list with the point AUC and bootstrap CI bounds.
#'
#' @export
os_clustered_bootstrap_auc <- function(y, predictions, cluster_id,
                                       n_boot = 2000L,
                                       stratify_by = NULL,
                                       ci_type = c("two-sided", "one-sided-lower"),
                                       alpha = 0.05,
                                       seed = 42L) {
  ci_type <- match.arg(ci_type)
  n <- length(y)
  if (length(predictions) != n || length(cluster_id) != n) {
    stop("y, predictions, cluster_id must all have the same length.", call. = FALSE)
  }
  y_num <- as.numeric(as.factor(y)) - 1L

  set.seed(seed)

  # Build per-cluster row index list, optionally stratified
  if (is.null(stratify_by)) {
    cluster_rows <- split(seq_len(n), factor(cluster_id))
    stratum_of_cluster <- rep("all", length(cluster_rows))
  } else {
    stratify_by <- as.factor(stratify_by)
    # assume cluster is nested within stratum; use majority stratum per cluster
    cluster_rows <- split(seq_len(n), factor(cluster_id))
    stratum_of_cluster <- vapply(cluster_rows, function(rows) {
      tt <- table(stratify_by[rows])
      names(tt)[which.max(tt)]
    }, FUN.VALUE = character(1L))
  }

  # Group cluster indices by stratum for stratified resampling
  strata <- split(seq_along(cluster_rows), stratum_of_cluster)

  draws <- numeric(n_boot)
  point_auc <- .auc_fast(y_num, predictions)

  for (b in seq_len(n_boot)) {
    picked_clusters <- unlist(
      lapply(strata, function(ix) sample(ix, length(ix), replace = TRUE)),
      use.names = FALSE
    )
    rows <- unlist(cluster_rows[picked_clusters], use.names = FALSE)
    if (length(unique(y_num[rows])) < 2L) {
      draws[b] <- NA_real_
      next
    }
    draws[b] <- .auc_fast(y_num[rows], predictions[rows])
  }

  draws <- draws[!is.na(draws)]
  if (ci_type == "two-sided") {
    lo <- stats::quantile(draws, alpha / 2)
    hi <- stats::quantile(draws, 1 - alpha / 2)
  } else {
    lo <- stats::quantile(draws, alpha)
    hi <- NA_real_
  }

  list(
    point_auc = point_auc,
    ci_lower = unname(lo),
    ci_upper = unname(hi),
    n_boot_effective = length(draws),
    ci_type = ci_type,
    alpha = alpha
  )
}


#' @title Specimen-Duplication Detection Across Cohorts
#'
#' @description
#' Detects putative specimen duplicates across multiple cohorts/datasets by
#' either (a) exact feature equality on a chosen fingerprint matrix, or
#' (b) an explicit crosswalk of internal specimen IDs.
#'
#' In circulating-miRNA cancer meta-analyses this is critical because the
#' same biological samples are often submitted to GEO under different
#' series accessions with different GSM IDs. Failing to dedup leads to
#' specimen-level train/test leakage and anti-conservative bootstrap CIs
#' (documented in GSE106817/GSE113486/GSE113740 for ovarian cases).
#'
#' @param X Numeric matrix with one row per sample, columns are features.
#' @param cohort Vector of cohort identity per sample.
#' @param sample_id Vector of sample IDs.
#' @param specimen_id Optional explicit internal specimen IDs (preferred).
#' @param tolerance Numeric tolerance for exact equality on feature vectors.
#'   Default 0 (numeric equality without rounding).
#' @param min_features Minimum number of jointly finite features required to
#'   declare a match. All jointly finite features must agree within
#'   `tolerance`. Default: all columns, which therefore requires complete
#'   profiles.
#'
#' @return A data frame with columns:
#'   \code{sample_id_a}, \code{cohort_a}, \code{sample_id_b}, \code{cohort_b},
#'   \code{specimen_id} (if provided), \code{match_type}.
#'
#' @export
os_detect_cross_cohort_duplicates <- function(X, cohort, sample_id,
                                              specimen_id = NULL,
                                              tolerance = 0,
                                              min_features = NULL) {
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 1L) {
    stop("X must be a numeric matrix with at least one feature.", call. = FALSE)
  }
  n <- nrow(X)
  if (length(cohort) != n || length(sample_id) != n) {
    stop("cohort and sample_id must match nrow(X).", call. = FALSE)
  }
  cohort <- as.character(cohort)
  sample_id <- as.character(sample_id)
  if (anyNA(cohort) || any(!nzchar(cohort)) ||
      anyNA(sample_id) || any(!nzchar(sample_id)) || anyDuplicated(sample_id)) {
    stop("cohort and sample_id must be non-missing; sample_id must be unique.",
         call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one finite non-negative number.", call. = FALSE)
  }
  if (is.null(min_features)) min_features <- ncol(X)
  if (!is.numeric(min_features) || length(min_features) != 1L ||
      is.na(min_features) || !is.finite(min_features) ||
      min_features != as.integer(min_features) || min_features < 1L ||
      min_features > ncol(X)) {
    stop("min_features must be one integer from 1 through ncol(X).",
         call. = FALSE)
  }
  min_features <- as.integer(min_features)

  matches <- list()
  matched_pair_keys <- character(0)
  pair_key <- function(i, j) paste(sort(c(sample_id[[i]], sample_id[[j]])),
                                    collapse = "\r")
  add_match <- function(i, j, specimen, match_type) {
    key <- pair_key(i, j)
    if (key %in% matched_pair_keys) return(invisible(FALSE))
    matches[[length(matches) + 1L]] <<- data.frame(
      sample_id_a = sample_id[[i]], cohort_a = cohort[[i]],
      sample_id_b = sample_id[[j]], cohort_b = cohort[[j]],
      specimen_id = specimen, match_type = match_type,
      stringsAsFactors = FALSE
    )
    matched_pair_keys <<- c(matched_pair_keys, key)
    invisible(TRUE)
  }

  # (a) Explicit specimen ID crosswalk
  if (!is.null(specimen_id)) {
    if (length(specimen_id) != n) {
      stop("specimen_id must match nrow(X).", call. = FALSE)
    }
    specimen_id <- trimws(as.character(specimen_id))
    specimen_id[is.na(specimen_id) | !nzchar(specimen_id)] <- NA_character_
    by_spec <- split(seq_len(n), specimen_id, drop = TRUE)
    for (sp in names(by_spec)) {
      idx <- by_spec[[sp]]
      if (length(idx) >= 2L) {
        for (left in seq_len(length(idx) - 1L)) {
          for (right in (left + 1L):length(idx)) {
            i <- idx[[left]]
            j <- idx[[right]]
            if (cohort[[i]] != cohort[[j]]) {
              add_match(i, j, sp, "specimen_id_crosswalk")
            }
          }
        }
      }
    }
  }

  # (b) Feature-equality duplicate detection
  rows_match <- function(i, j) {
    usable <- is.finite(X[i, ]) & is.finite(X[j, ])
    if (sum(usable) < min_features) return(FALSE)
    all(abs(X[i, usable] - X[j, usable]) <= tolerance)
  }

  compare_pairs <- function(idx) {
    if (length(idx) < 2L) return(invisible(NULL))
    for (left in seq_len(length(idx) - 1L)) {
      for (right in (left + 1L):length(idx)) {
        i <- idx[[left]]
        j <- idx[[right]]
        if (cohort[[i]] != cohort[[j]] && rows_match(i, j)) {
          specimen <- if (!is.null(specimen_id) &&
                          !is.na(specimen_id[[i]]) &&
                          identical(specimen_id[[i]], specimen_id[[j]])) {
            specimen_id[[i]]
          } else {
            NA_character_
          }
          add_match(
            i, j, specimen,
            if (tolerance == 0) "exact_feature_equality" else
              "tolerance_feature_equality"
          )
        }
      }
    }
    invisible(NULL)
  }

  # Complete, exact matching is the common large-corpus route. Hash only to
  # create candidates, then verify numerically so formatting collisions cannot
  # produce false matches. Missing/tolerant comparisons take the exhaustive,
  # explicitly requested route.
  complete_exact <- tolerance == 0 && min_features == ncol(X) && all(is.finite(X))
  if (complete_exact) {
    hash_of_row <- apply(X, 1L, function(v) {
      paste(format(v, digits = 17L, scientific = TRUE, trim = TRUE),
            collapse = "|")
    })
    by_hash <- split(seq_len(n), hash_of_row)
    for (idx in by_hash[lengths(by_hash) >= 2L]) compare_pairs(idx)
  } else {
    compare_pairs(seq_len(n))
  }

  if (length(matches) == 0L) {
    return(data.frame(
      sample_id_a = character(0), cohort_a = character(0),
      sample_id_b = character(0), cohort_b = character(0),
      specimen_id = character(0), match_type = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, matches)
}


#' @title One-Shot Bias Audit Report
#'
#' @description
#' Orchestrator that runs all bias-audit diagnostics on a model + data bundle
#' and returns a structured report. Intended to be emitted as a first-class
#' output of every OmicSelector biomarker analysis.
#'
#' @param X Feature matrix (n x p).
#' @param y Binary outcome.
#' @param cohort Cohort/dataset identity per sample.
#' @param predictions Optional numeric vector of model predictions (n).
#' @param covariates Optional named list of numeric or factor covariates
#'   (e.g., list(age = age_years, sex = sex_factor)).
#' @param specimen_id Optional specimen-level IDs for duplicate detection.
#' @param sample_id Optional sample-level IDs used when reporting duplicate
#'   pairs. If omitted, row names of \code{X} are used when available.
#' @param feature_names Optional feature labels.
#' @param seed Random seed.
#'
#' @return An \code{os_bias_audit} object (list) with fields:
#'   \code{bias_floor}, \code{covariate_aucs}, \code{feature_signal},
#'   \code{duplicates}, \code{provenance_summary},
#'   \code{classifier_vs_floor} (if predictions provided).
#'
#' @export
os_bias_audit <- function(X, y, cohort,
                          predictions = NULL,
                          covariates = NULL,
                          specimen_id = NULL,
                          sample_id = NULL,
                          feature_names = NULL,
                          seed = 42L) {
  bf <- os_bias_floor_auc(y, cohort, cv_folds = 5L, seed = seed)

  cov_out <- NULL
  if (!is.null(covariates) && length(covariates) > 0L) {
    cov_out <- lapply(names(covariates), function(nm) {
      os_covariate_only_auc(y, covariates[[nm]], covariate_name = nm)
    })
    names(cov_out) <- names(covariates)
  }

  fs <- os_per_feature_batch_signal(X, y, cohort, feature_names = feature_names)

  dup <- NULL
  if (!is.null(specimen_id) || !is.null(sample_id)) {
    sid <- if (!is.null(sample_id)) sample_id else rownames(X)
    if (is.null(sid)) sid <- as.character(seq_len(nrow(X)))
    dup <- os_detect_cross_cohort_duplicates(
      X, cohort, sample_id = sid, specimen_id = specimen_id
    )
  }

  provenance <- data.frame(
    cohort = names(table(cohort)),
    n = as.integer(table(cohort)),
    stringsAsFactors = FALSE
  )
  provenance$case_rate <- vapply(provenance$cohort, function(c) {
    mean(as.numeric(as.factor(y))[cohort == c] - 1L)
  }, FUN.VALUE = numeric(1L))

  classifier_vs_floor <- NULL
  if (!is.null(predictions)) {
    clf_auc <- .auc_fast(as.numeric(as.factor(y)) - 1L, predictions)
    classifier_vs_floor <- list(
      classifier_auc = clf_auc,
      bias_floor_auc = bf$apparent_auc,
      excess_over_floor = clf_auc - bf$apparent_auc,
      interpretation = if ((clf_auc - bf$apparent_auc) < 0.10) {
        "CONCERN: classifier AUC is within 0.10 of the cohort-identity floor; most apparent signal may be provenance, not biology."
      } else if ((clf_auc - bf$apparent_auc) < 0.20) {
        "MODERATE: classifier AUC exceeds the floor by < 0.20; report excess over floor as headline quantity and validate on cohort-balanced data."
      } else {
        "ACCEPTABLE separation: classifier AUC exceeds floor by >= 0.20; still report floor alongside classifier AUC for transparency."
      }
    )
  }

  out <- list(
    bias_floor = bf,
    covariate_aucs = cov_out,
    feature_signal = fs,
    duplicates = dup,
    provenance_summary = provenance,
    classifier_vs_floor = classifier_vs_floor,
    n = length(y),
    n_cohorts = length(unique(cohort)),
    audit_version = "1.0.0",
    audit_date = format(Sys.time(), "%Y-%m-%d")
  )
  class(out) <- c("os_bias_audit", "list")
  out
}


#' @export
print.os_bias_audit <- function(x, ...) {
  cat("OmicSelector Cohort-Provenance Bias Audit\n")
  cat(sprintf("  Audit version: %s | Date: %s\n", x$audit_version, x$audit_date))
  cat(sprintf("  n = %d samples across %d cohorts\n\n", x$n, x$n_cohorts))

  cat("=== Bias floor (cohort-identity-only classifier) ===\n")
  cat(sprintf("  Apparent AUC: %.3f\n", x$bias_floor$apparent_auc))
  if (!is.na(x$bias_floor$cv_auc)) {
    cat(sprintf("  5-fold CV AUC: %.3f\n", x$bias_floor$cv_auc))
  }
  cat(sprintf("  Interpretation: %s\n\n", x$bias_floor$interpretation))

  if (!is.null(x$classifier_vs_floor)) {
    cat("=== Classifier vs bias floor ===\n")
    cat(sprintf("  Classifier AUC: %.3f\n", x$classifier_vs_floor$classifier_auc))
    cat(sprintf("  Bias floor AUC: %.3f\n", x$classifier_vs_floor$bias_floor_auc))
    cat(sprintf("  Excess over floor: %.3f\n", x$classifier_vs_floor$excess_over_floor))
    cat(sprintf("  %s\n\n", x$classifier_vs_floor$interpretation))
  }

  if (!is.null(x$covariate_aucs) && length(x$covariate_aucs) > 0L) {
    cat("=== Covariate-only AUCs ===\n")
    for (nm in names(x$covariate_aucs)) {
      c0 <- x$covariate_aucs[[nm]]
      cat(sprintf("  %s: AUC=%.3f (%s, %s=%.2f)\n",
                  nm, c0$auc, c0$direction, c0$effect_metric, c0$effect_size))
    }
    cat("\n")
  }

  if (!is.null(x$feature_signal)) {
    cat("=== Top-5 features by dataset-case confounding ===\n")
    head_fs <- utils::head(x$feature_signal, 5L)
    for (i in seq_len(nrow(head_fs))) {
      cat(sprintf("  %s: corr=%.3f | F_dataset_healthy=%.2f | F_case=%.2f\n",
                  head_fs$feature[i], head_fs$dataset_case_corr[i],
                  head_fs$F_dataset_healthy[i], head_fs$F_case[i]))
    }
    cat("\n")
  }

  if (!is.null(x$duplicates) && nrow(x$duplicates) > 0L) {
    cat(sprintf("=== Cross-cohort duplicates detected: %d pairs ===\n",
                nrow(x$duplicates)))
    cat(sprintf("  Match types: %s\n",
                paste(names(table(x$duplicates$match_type)), collapse = ", ")))
    cat("  Recommend specimen-level dedup + clustered bootstrap.\n\n")
  } else if (!is.null(x$duplicates)) {
    cat("=== No cross-cohort duplicates detected. ===\n\n")
  }

  invisible(x)
}


#' @title Plain-Text Bias-Audit Report
#'
#' @description
#' Convenience wrapper that runs \code{\link{os_bias_audit}} on a model/data
#' bundle and returns the formatted report as a single character string.
#' Useful when the report needs to be captured into a log, supplementary
#' table, or pasted into a manuscript figure legend without relying on
#' the console print method.
#'
#' @inheritParams os_bias_audit
#' @param audit An optional pre-computed \code{os_bias_audit} object. If
#'   supplied, the audit is not recomputed and the remaining arguments are
#'   ignored. This makes it cheap to produce the text report alongside the
#'   structured result when both are needed.
#'
#' @return A length-one character vector containing the full text report.
#'   The same content that \code{print.os_bias_audit} would write to the
#'   console is captured, with an additional attribute \code{"audit"} that
#'   holds the underlying \code{os_bias_audit} object for downstream use.
#'
#' @examples
#' \dontrun{
#' txt <- os_bias_audit_report(
#'   X = X, y = y, cohort = cohort,
#'   predictions = preds,
#'   covariates = list(age = age_years)
#' )
#' cat(txt)
#' writeLines(txt, "bias_audit.txt")
#' attr(txt, "audit")$bias_floor$apparent_auc
#' }
#'
#' @export
os_bias_audit_report <- function(X = NULL, y = NULL, cohort = NULL,
                                 predictions = NULL,
                                 covariates = NULL,
                                 specimen_id = NULL,
                                 sample_id = NULL,
                                 feature_names = NULL,
                                 seed = 42L,
                                 audit = NULL) {
  if (is.null(audit)) {
    if (is.null(X) || is.null(y) || is.null(cohort)) {
      stop("Provide either `audit` or all of `X`, `y`, `cohort`.", call. = FALSE)
    }
    audit <- os_bias_audit(
      X = X, y = y, cohort = cohort,
      predictions = predictions,
      covariates = covariates,
      specimen_id = specimen_id,
      sample_id = sample_id,
      feature_names = feature_names,
      seed = seed
    )
  } else if (!inherits(audit, "os_bias_audit")) {
    stop("`audit` must be an `os_bias_audit` object.", call. = FALSE)
  }

  txt <- paste(utils::capture.output(print(audit)), collapse = "\n")
  attr(txt, "audit") <- audit
  class(txt) <- c("os_bias_audit_report", "character")
  txt
}


# Internal AUC: trapezoidal on rank-based approximation (no dependency on pROC).
.auc_fast <- function(y, score) {
  if (length(unique(y)) != 2L) return(NA_real_)
  y <- as.numeric(y)
  pos <- score[y == 1L]
  neg <- score[y == 0L]
  np <- length(pos); nn <- length(neg)
  if (np == 0L || nn == 0L) return(NA_real_)
  # Mann-Whitney U
  r <- rank(c(pos, neg), ties.method = "average")
  sum_r_pos <- sum(r[seq_len(np)])
  U <- sum_r_pos - np * (np + 1L) / 2
  U / (np * nn)
}
