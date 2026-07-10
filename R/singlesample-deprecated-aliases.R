# Back-compatibility aliases for the former `paper3_*` public API.
#
# The single-sample-deployable method-bank machinery was renamed from the
# meaningless internal codename `paper3_*` to the descriptive `singlesample_*`.
# These thin aliases keep existing callers (including the older analysis
# pipeline) working unchanged. Prefer the `singlesample_*` names in new code.

#' Deprecated `paper3_*` aliases
#'
#' Back-compatibility aliases for the renamed `singlesample_*` functions. Each
#' forwards verbatim to its `singlesample_*` counterpart and is retained so that
#' code written against the former `paper3_*` names keeps working. Prefer the
#' `singlesample_*` names; the `paper3_*` aliases may be removed in a future
#' major release.
#'
#' @param ... Arguments passed on to the corresponding `singlesample_*` function.
#' @return The value of the corresponding `singlesample_*` function.
#' @name paper3-deprecated
#' @keywords internal
NULL

#' @rdname paper3-deprecated
#' @export
paper3_assert_method_bank_exports <- function(...) singlesample_assert_method_bank_exports(...)

#' @rdname paper3-deprecated
#' @export
paper3_assert_row_equivariant <- function(...) singlesample_assert_row_equivariant(...)

#' @rdname paper3-deprecated
#' @export
paper3_bh_fdr_correct_blocked <- function(...) singlesample_bh_fdr_correct_blocked(...)

#' @rdname paper3-deprecated
#' @export
paper3_bh_fdr_correct_matched_null <- function(...) singlesample_bh_fdr_correct_matched_null(...)

#' @rdname paper3-deprecated
#' @export
paper3_hanley_mcneil_auc_ci <- function(...) singlesample_hanley_mcneil_auc_ci(...)

#' @rdname paper3-deprecated
#' @export
paper3_holm_correct_familywise <- function(...) singlesample_holm_correct_familywise(...)

#' @rdname paper3-deprecated
#' @export
paper3_is_row_equivariant <- function(...) singlesample_is_row_equivariant(...)

#' @rdname paper3-deprecated
#' @export
paper3_make_loco_splits <- function(...) singlesample_make_loco_splits(...)

#' @rdname paper3-deprecated
#' @export
paper3_make_locto_splits <- function(...) singlesample_make_locto_splits(...)

#' @rdname paper3-deprecated
#' @export
paper3_matched_null_benchmark <- function(...) singlesample_matched_null_benchmark(...)

#' @rdname paper3-deprecated
#' @export
paper3_matched_null_benchmark_cv <- function(...) singlesample_matched_null_benchmark_cv(...)

#' @rdname paper3-deprecated
#' @export
paper3_method_bank <- function(...) singlesample_method_bank(...)

#' @rdname paper3-deprecated
#' @export
paper3_method_roster <- function(...) singlesample_method_roster(...)

#' @rdname paper3-deprecated
#' @export
paper3_paired_auc_diff_se <- function(...) singlesample_paired_auc_diff_se(...)

#' @rdname paper3-deprecated
#' @export
paper3_register_score_adapter <- function(...) singlesample_register_score_adapter(...)

#' @rdname paper3-deprecated
#' @export
paper3_score_call <- function(...) singlesample_score_call(...)

#' @rdname paper3-deprecated
#' @export
paper3_technology_lift_delong <- function(...) singlesample_technology_lift_delong(...)
