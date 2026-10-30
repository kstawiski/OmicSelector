placements <- function(y, score) {
    pos_idx <- which(y == 1L)
    neg_idx <- which(y == 0L)
    if (length(pos_idx) < 2L || length(neg_idx) < 2L) 
        return(NULL)
    cmp <- outer(score[pos_idx], score[neg_idx], function(a, b) (a > b) + 0.5 * (a == b))
    list(auc = mean(cmp), pos_idx = pos_idx, neg_idx = neg_idx, v10 = rowMeans(cmp), v01 = colMeans(cmp))
}
corrected_repeated_cv <- function(lift, test_train_ratio, expected_m) {
    expected_m <- as.integer(expected_m)
    if (length(expected_m) != 1L || is.na(expected_m) || expected_m < 2L) {
        stop("expected_m must be one integer >= 2")
    }
    ok <- is.finite(lift) & is.finite(test_train_ratio) & test_train_ratio > 0
    lift <- lift[ok]
    test_train_ratio <- test_train_ratio[ok]
    m <- length(lift)
    mean_ratio <- if (m) 
        mean(test_train_ratio)
    else NA_real_
    if (m != expected_m) 
        return(list(complete = FALSE, observed_m = m, expected_m = expected_m, estimate = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p_zero = NA_real_, p_margin = NA_real_, df = expected_m - 1L, correction = NA_real_, mean_ratio = mean_ratio))
    estimate <- mean(lift)
    correction <- 1/expected_m + mean_ratio
    se <- sqrt(max(0, correction * stats::var(lift)))
    df <- expected_m - 1L
    critical <- stats::qt(0.975, df = df)
    list(complete = TRUE, observed_m = m, expected_m = expected_m, estimate = estimate, se = se, ci_lo = estimate - critical * se, ci_hi = estimate + critical * se, p_zero = if (se > 0) stats::pt(estimate/se, df = df, lower.tail = FALSE) else NA_real_, p_margin = if (se > 0) stats::pt((estimate - 0.05)/se, df = df, lower.tail = FALSE) else NA_real_, df = df, correction = correction, mean_ratio = mean(test_train_ratio))
}
cluster_auc_summary <- function(d, stratum_cols, corrected = FALSE, expected_strata = NULL) {
    d <- copy(d[is.finite(score_m) & is.finite(score_b) & !is.na(y)])
    if (!nrow(d)) 
        return(NULL)
    if (any(d[, uniqueN(y), by = group_id]$V1 != 1L)) {
        stop("A biological group carries discordant outcome labels")
    }
    d[, `:=`(stratum, do.call(paste, c(.SD, sep = "::"))), .SDcols = stratum_cols]
    influence <- list()
    auc_rows <- list()
    group_auc_rows <- list()
    valid <- 0L
    total_groups <- uniqueN(d$group_id)
    for (s in sort(unique(d$stratum))) {
        x <- d[stratum == s]
        pm <- placements(x$y, x$score_m)
        pb <- placements(x$y, x$score_b)
        if (is.null(pm) || is.null(pb)) 
            next
        valid <- valid + 1L
        ratio <- if (all(c("n_test", "n_train") %in% names(x))) {
            unique_ratio <- unique(x$n_test/x$n_train)
            if (length(unique_ratio) != 1L) {
                stop("A repeated-CV stratum does not have one test/train ratio")
            }
            unique_ratio
        }
        else NA_real_
        auc_rows[[valid]] <- c(method = pm$auc, baseline = pb$auc, lift = pm$auc - pb$auc, ratio = ratio)
        gx <- x[, .(y = unique(y), score_m = mean(score_m), score_b = mean(score_b)), by = group_id]
        if (any(lengths(gx$y) != 1L)) 
            stop("Discordant group outcome within stratum")
        pgm <- placements(as.integer(gx$y), gx$score_m)
        pgb <- placements(as.integer(gx$y), gx$score_b)
        if (is.null(pgm) || is.null(pgb)) {
            stop("Group-collapsed fold lost an outcome class")
        }
        n_test_groups <- nrow(gx)
        n_train_groups <- total_groups - n_test_groups
        group_auc_rows[[valid]] <- c(method = pgm$auc, baseline = pgb$auc, lift = pgm$auc - pgb$auc, ratio = if (n_train_groups > 0L) n_test_groups/n_train_groups else NA_real_)
        im <- numeric(nrow(x))
        ib <- numeric(nrow(x))
        im[pm$pos_idx] <- (pm$v10 - pm$auc)/length(pm$pos_idx)
        im[pm$neg_idx] <- (pm$v01 - pm$auc)/length(pm$neg_idx)
        ib[pb$pos_idx] <- (pb$v10 - pb$auc)/length(pb$pos_idx)
        ib[pb$neg_idx] <- (pb$v01 - pb$auc)/length(pb$neg_idx)
        influence[[valid]] <- data.table(group_id = x$group_id, method = im, baseline = ib)
    }
    if (!valid) 
        return(NULL)
    aucs <- do.call(rbind, auc_rows)
    group_aucs <- do.call(rbind, group_auc_rows)
    infl <- rbindlist(influence)
    infl[, `:=`(method = method/valid, baseline = baseline/valid)]
    cl <- infl[, .(method = sum(method), baseline = sum(baseline)), by = group_id]
    G <- nrow(cl)
    if (G < 3L) 
        return(NULL)
    correction <- G/(G - 1)
    var_m <- correction * sum(cl$method^2)
    var_b <- correction * sum(cl$baseline^2)
    var_d <- correction * sum((cl$method - cl$baseline)^2)
    auc_m <- mean(aucs[, "method"])
    auc_b <- mean(aucs[, "baseline"])
    lift <- auc_m - auc_b
    se_m <- sqrt(max(0, var_m))
    se_b <- sqrt(max(0, var_b))
    se_d <- sqrt(max(0, var_d))
    rcv <- if (isTRUE(corrected)) {
        corrected_repeated_cv(aucs[, "lift"], aucs[, "ratio"], expected_strata)
    }
    else NULL
    group_rcv <- if (isTRUE(corrected)) {
        corrected_repeated_cv(group_aucs[, "lift"], group_aucs[, "ratio"], expected_strata)
    }
    else NULL
    primary_se <- if (isTRUE(corrected)) 
        rcv$se
    else se_d
    primary_ci_lo <- if (isTRUE(corrected)) 
        rcv$ci_lo
    else lift - 1.96 * se_d
    primary_ci_hi <- if (isTRUE(corrected)) 
        rcv$ci_hi
    else lift + 1.96 * se_d
    primary_p_zero <- if (isTRUE(corrected)) 
        rcv$p_zero
    else if (se_d > 0) {
        pnorm(lift/se_d, lower.tail = FALSE)
    }
    else NA_real_
    primary_p_margin <- if (isTRUE(corrected)) 
        rcv$p_margin
    else if (se_d > 0) {
        pnorm((lift - 0.05)/se_d, lower.tail = FALSE)
    }
    else NA_real_
    data.table(auc_method = auc_m, auc_method_se_cluster = se_m, auc_method_ci_lo = max(0, auc_m - 1.96 * se_m), auc_method_ci_hi = min(1, auc_m + 1.96 * se_m), auc_baseline = auc_b, auc_baseline_se_cluster = se_b, auc_baseline_ci_lo = max(0, auc_b - 1.96 * se_b), auc_baseline_ci_hi = min(1, auc_b + 1.96 * se_b), lift = lift, lift_se = primary_se, lift_ci_lo = primary_ci_lo, lift_ci_hi = primary_ci_hi, p_lift_gt_zero = primary_p_zero, p_lift_gt_margin_0_05 = primary_p_margin, lift_se_corrected_rcv = if (isTRUE(corrected)) 
        rcv$se
    else NA_real_, corrected_rcv_complete = if (isTRUE(corrected)) 
        rcv$complete
    else NA, corrected_rcv_expected_strata = if (isTRUE(corrected)) 
        rcv$expected_m
    else NA_integer_, corrected_rcv_df = if (isTRUE(corrected)) 
        rcv$df
    else NA_integer_, corrected_rcv_factor = if (isTRUE(corrected)) 
        rcv$correction
    else NA_real_, mean_test_train_ratio = if (isTRUE(corrected)) 
        rcv$mean_ratio
    else NA_real_, lift_se_group_cluster_diagnostic = se_d, lift_ci_lo_group_cluster_diagnostic = lift - 1.96 * se_d, lift_ci_hi_group_cluster_diagnostic = lift + 1.96 * se_d, auc_method_group_collapsed = mean(group_aucs[, "method"]), auc_baseline_group_collapsed = mean(group_aucs[, "baseline"]), lift_group_collapsed = mean(group_aucs[, "lift"]), lift_group_collapsed_se_corrected_rcv = if (isTRUE(corrected)) 
        group_rcv$se
    else NA_real_, group_corrected_rcv_complete = if (isTRUE(corrected)) 
        group_rcv$complete
    else NA, lift_group_collapsed_ci_lo = if (isTRUE(corrected)) 
        group_rcv$ci_lo
    else NA_real_, lift_group_collapsed_ci_hi = if (isTRUE(corrected)) 
        group_rcv$ci_hi
    else NA_real_, lift_profile_minus_group_collapsed = lift - mean(group_aucs[, "lift"]), n_groups = G, n_profiles = uniqueN(d$sample_id), n_predictions = nrow(d), n_valid_strata = valid)
}
