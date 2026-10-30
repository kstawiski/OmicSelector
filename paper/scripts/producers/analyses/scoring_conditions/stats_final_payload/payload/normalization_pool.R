meta_fit <- function(x, min_k = 5L) {
    x <- x[is.finite(lift) & is.finite(lift_se) & lift_se > 0]
    if (nrow(x) < min_k) 
        return(data.table(k = nrow(x), pooled_lift = NA_real_, pooled_se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p_lift_gt_zero = NA_real_, p_lift_gt_margin_0_05 = NA_real_, tau2 = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_))
    f <- metafor::rma.uni(yi = x$lift, sei = x$lift_se, method = "REML", test = "knha")
    pr <- predict(f)
    data.table(k = f$k, pooled_lift = as.numeric(f$b), pooled_se = f$se, ci_lo = f$ci.lb, ci_hi = f$ci.ub, p_lift_gt_zero = pt(as.numeric(f$b)/f$se, f$ddf, lower.tail = FALSE), p_lift_gt_margin_0_05 = pt((as.numeric(f$b) - 0.05)/f$se, f$ddf, lower.tail = FALSE), tau2 = f$tau2, I2 = f$I2, pi_lo = pr$pi.lb, pi_hi = pr$pi.ub)
}
pool_safe <- function(d) {
    x <- data.table(lift = d$estimate, lift_se = d$se)
    msg <- character()
    ans <- tryCatch(withCallingHandlers({
        z <- meta_fit(x)
        z[, `:=`(failure, NA_character_)]
        z
    }, warning = function(w) {
        msg <<- c(msg, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) {
        msg <<- c(msg, conditionMessage(e))
        xx <- x[is.finite(lift) & is.finite(lift_se) & lift_se > 0]
        if (nrow(xx) >= 5) {
            f <- tryCatch(metafor::rma.uni(yi = xx$lift, sei = xx$lift_se, method = "REML", test = "knha", control = list(stepadj = 0.5, maxiter = 1000)), error = function(e) NULL)
            if (!is.null(f)) {
                pr <- predict(f)
                return(data.table(k = f$k, pooled_lift = as.numeric(f$b), pooled_se = f$se, ci_lo = f$ci.lb, ci_hi = f$ci.ub, p_lift_gt_zero = pt(as.numeric(f$b)/f$se, f$ddf, lower.tail = FALSE), p_lift_gt_margin_0_05 = pt((as.numeric(f$b) - 0.05)/f$se, f$ddf, lower.tail = FALSE), tau2 = f$tau2, I2 = f$I2, pi_lo = pr$pi.lb, pi_hi = pr$pi.ub, failure = NA_character_))
            }
        }
        z <- meta_fit(data.table(lift = numeric(), lift_se = numeric()))
        z[, `:=`(k = nrow(xx), failure = conditionMessage(e))]
        z
    })
    ans[, `:=`(warnings, paste(unique(msg), collapse = "; "))]
    ans
}
