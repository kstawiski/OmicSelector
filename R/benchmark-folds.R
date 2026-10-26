#' Reproduce the archived benchmark grouped folds
#'
#' Groups are enumerated in first-row order, classified by the rounded mean
#' binary outcome, shuffled within outcome (0 then 1), and assigned round-robin.
#' This is the benchmark harness algorithm, including its change to the global
#' RNG state through set.seed(seed). Mixed-outcome groups use round(mean(y));
#' a mean of 0.5 rounds to 0. Supply rows in the original benchmark order and
#' use the same R RNGkind settings (the archive uses the R default).
#'
#' The existing os_make_grouped_stratified_folds() remains available for
#' balancing row counts and validating mixed groups. Its partitions differ.
#' @param y Numeric binary outcome vector, containing only 0 and 1.
#' @param group_id Nonmissing group identifiers, one per row.
#' @param k Number of folds, at least two.
#' @param seed Random seed.
#' @return An unnamed list of sorted integer test-row indices.
#' @export
os_make_benchmark_folds <- function (y, group_id, k = 5L, seed = 42L)
{
    if (!is.numeric(y) || !length(y) || anyNA(y) || any(!y %in% c(0, 1)))
        stop("y must contain numeric binary outcomes (0 and 1).", call. = FALSE)
    if (length(group_id) != length(y) || anyNA(group_id))
        stop("group_id must be nonmissing and have length(y) entries.", call. = FALSE)
    if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k < 2 || k != as.integer(k))
        stop("k must be an integer >= 2.", call. = FALSE)
    set.seed(seed)
    ugroups <- unique(group_id)
    group_cls <- vapply(ugroups, function(g) {
        yg <- y[group_id == g]
        as.integer(round(mean(yg)))
    }, integer(1L))
    fold_assign <- integer(length(ugroups))
    for (cls in c(0L, 1L)) {
        g_idx <- which(group_cls == cls)
        if (length(g_idx) == 0L)
            next
        g_idx <- g_idx[sample.int(length(g_idx))]
        fold_assign[g_idx] <- rep_len(seq_len(k), length(g_idx))
    }
    folds <- vector("list", k)
    for (fi in seq_len(k)) {
        test_groups <- ugroups[fold_assign == fi]
        folds[[fi]] <- sort(which(group_id %in% test_groups))
    }
    folds
}
