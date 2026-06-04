#' Paper 3 Amendment #4 method roster (frozen expansion)
#'
#' Returns the frozen 72-row method roster for the OmicSelector manuscript
#' method-bank expansion (Statistical Analysis Plan Amendment #4). Unlike
#' [paper3_method_bank()], which is the export-check view of functions already
#' callable in the package, the roster enumerates the full pre-specified set of
#' single-sample scoring methods (implemented and not-yet-implemented) together
#' with the routing metadata the benchmark instrument consumes.
#'
#' The roster is shipped as package data
#' (`inst/extdata/paper3_method_manifest.csv`) so the R package and the
#' analysis instrument read one authoritative artifact. Each row carries:
#' \describe{
#'   \item{method_id}{Stable identifier (matches the analysis manifest).}
#'   \item{family}{Human-readable grouping letter (A-O, NC).}
#'   \item{estimand}{`within` or `transfer`; routes the SE/FDR family.}
#'   \item{role}{`baseline`, `discriminator`, `discriminator-conditional`, or
#'     `negative_control`.}
#'   \item{tier}{Execution substrate: `R1` (pure-R) or `R2` (reticulate).}
#'   \item{dep_route}{Dependency / implementation route.}
#'   \item{fit_fn, score_fn}{Package fit / score function names.}
#'   \item{pkg_status}{`present`, `primitive`, `new`, or `new-secondary`.}
#'   \item{row_source}{Instrument assembly path: `per_cohort_rows`,
#'     `stratum_pooled_upstream`, or `transfer_fold_rows`.}
#'   \item{lopbo_mechanism}{`within_cohort`, `transfer_fold`, `deferred`,
#'     `m4_decision` (placeholder resolved at benchmark time), or `na`.}
#' }
#' plus derived `single_sample_gate`, `is_discriminator`, and
#' `is_negative_control` columns.
#'
#' @return A data frame with one row per frozen-roster method.
#' @seealso [paper3_method_bank()], [paper3_score_call()],
#'   [paper3_assert_row_equivariant()]
#' @export
paper3_method_roster <- function() {
  path <- system.file("extdata", "paper3_method_manifest.csv",
                      package = "OmicSelector")
  if (!nzchar(path) || !file.exists(path)) {
    stop("paper3_method_manifest.csv not found in package extdata.",
         call. = FALSE)
  }
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                       na.strings = c("NA", ""))
  req_cols <- c("method_id", "family", "estimand", "role", "tier", "dep_route",
                "fit_fn", "score_fn", "pkg_status", "notes", "row_source",
                "lopbo_mechanism")
  miss <- setdiff(req_cols, names(m))
  if (length(miss) > 0L) {
    stop("paper3 method manifest is missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  m <- m[, req_cols, drop = FALSE]
  .paper3_validate_roster(m)
  m$single_sample_gate <- "required"
  m$is_discriminator <- m$role %in% c("discriminator", "discriminator-conditional")
  m$is_negative_control <- m$role == "negative_control"
  rownames(m) <- NULL
  m
}

#' Frozen Amendment #4 roster identifiers
#'
#' Independent enumeration of the 72 pre-specified method identifiers,
#' transcribed from the prespecification (SAP Amendment #4 sections 4.1-4.4 and
#' the round-2 survey, section B). Kept separate from the shipped manifest so a
#' completeness test can assert the manifest implements exactly this frozen set
#' (`setequal`): editing the CSV without editing this vector (or vice versa)
#' fails the gate.
#'
#' @return Character vector of 72 method identifiers.
#' @keywords internal
.paper3_frozen_roster_ids <- function() {
  baseline <- "ws-rclr-trimmed"
  within_existing <- c(                                  # 6 existing within disc
    "ws-alr-pivot", "ws-mad-logratio", "ws-dominance",
    "frozen-ruv", "robust-pca-resid", "hemolysis-resid")
  transfer_existing <- c(                                # 17 existing transfer disc
    "kit-stratified-rclr", "kit-fe-adjusted-alr", "kit-orthogonal-ilr",
    "kit-residual-mad", "kit-stable-anchor-rclr", "hemolysis-kit-dual-anchor",
    "rin-weighted-score", "tech-stratified-rclr", "cross-tech-harmonized-rclr",
    "tech-residualized-alr", "biofluid-stratified-rclr", "biofluid-anchor-rclr",
    "biofluid-residualized-alr", "blood-cell-marker-excl-rclr",
    "dann", "cvae", "icp")
  within_new <- c(                                       # 35 new within disc
    "reo-ktsp", "reo-pairratio", "reo-singscore", "reo-ucell",
    "bal-selbal", "ws-balance-ilr", "lrt-bw", "lrt-copula", "conf-mondrian",
    "ai-scarf", "inv-scatter", "inv-olbp", "tda-ph", "sig-path", "frac-mfdfa",
    "img-gasfcnn", "ot-lot", "ot-slicedlrt", "kme-witness", "ig-fisherrao",
    "dre-ulsif", "gsp-gft", "man-nystrom", "unc-sngp", "ssl-vicreg",
    "coda-deepcoda", "tab-tabicl", "tab-tabdpt", "inv-glcm", "inv-fdaqf",
    "inv-css", "lrt-vinecopula", "proto-net", "lrt-deepmaha", "coda-codacore")
  within_conditional <- "ai-tabpfn"                      # +1 conditional within
  transfer_new <- c(                                     # 7 new transfer disc
    "reo-metaktsp", "dro-group", "dro-vrex", "moe-gated",
    "dg-fishr", "dg-ibirm", "sel-stablemate")
  negative_controls <- c(                                # 5 negative controls
    "nc-mahalanobis", "nc-conformal-anom", "nc-isoforest",
    "nc-ecod-copod", "nc-sinkhorn-single")
  c(baseline, within_existing, transfer_existing, within_new,
    within_conditional, transfer_new, negative_controls)
}

.paper3_roster_enums <- list(
  estimand        = c("within", "transfer"),
  role            = c("baseline", "discriminator", "discriminator-conditional",
                      "negative_control"),
  tier            = c("R1", "R2"),
  pkg_status      = c("present", "primitive", "new", "new-secondary"),
  row_source      = c("per_cohort_rows", "stratum_pooled_upstream",
                      "transfer_fold_rows"),
  lopbo_mechanism = c("within_cohort", "transfer_fold", "deferred",
                      "m4_decision", "na")
)

.paper3_validate_roster <- function(m) {
  if (anyNA(m$method_id) || any(!nzchar(m$method_id))) {
    stop("paper3 roster has missing method_id values.", call. = FALSE)
  }
  if (any(duplicated(m$method_id))) {
    stop("paper3 roster has duplicate method_id values: ",
         paste(unique(m$method_id[duplicated(m$method_id)]), collapse = ", "),
         call. = FALSE)
  }
  for (col in names(.paper3_roster_enums)) {
    allowed <- .paper3_roster_enums[[col]]
    bad <- setdiff(unique(m[[col]]), allowed)
    if (length(bad) > 0L) {
      stop(sprintf("paper3 roster column '%s' has disallowed values: %s",
                   col, paste(bad, collapse = ", ")), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Canonical single-sample score adapter
# ---------------------------------------------------------------------------

# Per-method adapter registry for non-canonical (legacy) score signatures.
# New methods follow the canonical convention score_*(model, X, meta = NULL)
# and need no adapter. Legacy/present methods whose signature differs register
# a shim via paper3_register_score_adapter() in their own per-method gate.
#
# NOTE for downstream method files: register adapters from the package load hook
# (.onLoad in R/zzz.R), not at the top level of a source file, so registration
# does not depend on R's alphabetical source load order relative to this file.
.paper3_score_adapter_registry <- new.env(parent = emptyenv())

#' Register a non-canonical score adapter for a roster method
#'
#' New roster methods follow the canonical convention
#' `score_*(model, X, meta = NULL)` and need no adapter. Legacy ("present")
#' methods whose score function uses a different argument order register a shim
#' mapping them onto the canonical `function(model, X, meta)` signature.
#'
#' @param method_id Roster method identifier.
#' @param fn A function `function(model, X, meta)` returning one numeric score
#'   per row of `X`.
#' @return Invisibly `TRUE`.
#' @seealso [paper3_score_call()]
#' @export
paper3_register_score_adapter <- function(method_id, fn) {
  if (!is.character(method_id) || length(method_id) != 1L || is.na(method_id)) {
    stop("method_id must be a single non-NA string.", call. = FALSE)
  }
  if (!is.function(fn)) stop("fn must be a function.", call. = FALSE)
  assign(method_id, fn, envir = .paper3_score_adapter_registry)
  invisible(TRUE)
}

#' Score a roster method with the canonical single-sample interface
#'
#' Resolves a roster method to its score function and invokes it through the
#' canonical adapter `function(model, X_rows, meta_rows) -> numeric`. New methods
#' follow `score_*(model, X, meta = NULL)` directly; legacy ("present") methods
#' MUST be dispatched through a shim registered with
#' [paper3_register_score_adapter()] (their argument order predates the canonical
#' convention, so calling them positionally would mis-bind). Methods whose score
#' function is not yet implemented raise an explicit not-yet-implemented error.
#'
#' @param method_id Roster method identifier.
#' @param model A fitted model object produced by the method's `fit_*` function.
#' @param X Numeric matrix (samples x features), or a data frame coercible to one.
#' @param meta Optional per-sample metadata (data frame with `nrow(X)` rows).
#' @param roster The roster to resolve against; defaults to
#'   [paper3_method_roster()].
#' @return Numeric vector of one score per row of `X`.
#' @seealso [paper3_method_roster()], [paper3_assert_row_equivariant()]
#' @export
paper3_score_call <- function(method_id, model, X, meta = NULL,
                              roster = paper3_method_roster()) {
  row <- roster[roster$method_id == method_id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop(sprintf("method_id '%s' not found uniquely in roster.", method_id),
         call. = FALSE)
  }
  X <- .paper3_as_score_matrix(X)
  if (!is.null(meta)) {
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (nrow(meta) != nrow(X)) {
      stop("meta must have one row per row of X.", call. = FALSE)
    }
  }
  if (exists(method_id, envir = .paper3_score_adapter_registry,
             inherits = FALSE)) {
    adapter <- get(method_id, envir = .paper3_score_adapter_registry,
                   inherits = FALSE)
  } else if (identical(row$pkg_status[[1]], "present")) {
    stop(sprintf(paste0("method '%s' is a present/legacy-signature scorer: ",
                        "register a canonical adapter via ",
                        "paper3_register_score_adapter('%s', ",
                        "function(model, X, meta) ...) before scoring it ",
                        "through paper3_score_call()."),
                 method_id, method_id), call. = FALSE)
  } else {
    adapter <- .paper3_canonical_adapter(method_id, row$score_fn[[1]])
  }
  out <- adapter(model, X, meta)
  .paper3_validate_score_output(out, nrow(X), method_id)
}

.paper3_canonical_adapter <- function(method_id, score_fn_name) {
  if (is.na(score_fn_name) || !nzchar(score_fn_name)) {
    stop(sprintf("method '%s' has no score function declared.", method_id),
         call. = FALSE)
  }
  fn <- .paper3_resolve_fn(score_fn_name)
  if (is.null(fn)) {
    stop(sprintf(paste0("method '%s' is not yet implemented: score function ",
                        "'%s' not found in the OmicSelector namespace. ",
                        "Implement fit_*/score_* and, if its signature is not ",
                        "the canonical score(model, X, meta), register an ",
                        "adapter via paper3_register_score_adapter()."),
                 method_id, score_fn_name), call. = FALSE)
  }
  function(model, X, meta) fn(model, X, meta)
}

# Resolve a score function by name in the OmicSelector namespace ONLY
# (inherits = FALSE). A stray object on the global search path must NOT satisfy
# a not-yet-implemented roster method, so the base/global fallback is removed.
.paper3_resolve_fn <- function(name) {
  ns <- tryCatch(asNamespace("OmicSelector"), error = function(e) NULL)
  if (!is.null(ns) &&
      exists(name, envir = ns, mode = "function", inherits = FALSE)) {
    return(get(name, envir = ns, mode = "function", inherits = FALSE))
  }
  NULL
}

.paper3_as_score_matrix <- function(X) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) {
    stop("X must be a matrix or data frame of samples x features.", call. = FALSE)
  }
  if (!is.numeric(X)) storage.mode(X) <- "double"
  X
}

# Enforce the documented contract: a plain numeric vector of length n. An
# n-by-1 numeric matrix is normalized to a vector; a non-numeric result
# (character/factor/logical) or any other dimensioned output is rejected by type
# rather than blindly coerced (coercion would silently accept "1.5" or turn a
# non-numeric string into NA).
.paper3_validate_score_output <- function(out, n, method_id) {
  d <- dim(out)
  if (!is.null(d) && !(length(d) == 2L && d[2] == 1L)) {
    stop(sprintf(paste0("method '%s' score_call must return a plain numeric ",
                        "vector, not an array with dim %s."),
                 method_id, paste(d, collapse = "x")), call. = FALSE)
  }
  if (!is.numeric(out)) {
    stop(sprintf(paste0("method '%s' score_call must return a numeric vector ",
                        "(got %s)."), method_id, class(out)[1]), call. = FALSE)
  }
  out <- as.numeric(out)  # lossless: strips an n-by-1 numeric matrix dim only
  if (length(out) != n) {
    stop(sprintf("method '%s' score_call must return %d numeric scores (got %d).",
                 method_id, n, length(out)), call. = FALSE)
  }
  out
}

# ---------------------------------------------------------------------------
# Generic row-equivariance (single-sample deployability) test harness
# ---------------------------------------------------------------------------

#' Assert that a scorer is row-equivariant (single-sample deployable)
#'
#' The decisive inclusion property for the manuscript method bank: at inference a
#' specimen's score must be computable from that specimen alone given a frozen
#' fitted model, with no dependence on which other specimens are scored alongside
#' it. This harness implements the strengthened single-sample-deployability unit
#' test prespecified in SAP Amendment #4 (section 2): for the same specimen, the
#' batch score must equal its standalone score across ALL of
#' \itemize{
#'   \item (a) the singleton `{x_i}` (the decisive single-sample case);
#'   \item (b) differently-composed batches (varied size/composition) and an
#'     adversarial batch (the specimen among extreme, scaled neighbours);
#'   \item (c) row-order permutations and batches with duplicated rows present;
#'   \item (d) no model-state mutation: the model's bytes are identical before
#'     and after scoring (catches BatchNorm running-stat updates / online
#'     recalibration).
#' }
#' Permutation alone cannot detect mean/quantile-centering coupling (column means
#' are permutation-invariant), which is why the singleton, subset, duplicate, and
#' adversarial cases always run.
#'
#' @param score_fun A function `function(model, X, meta)` returning one numeric
#'   score per row of `X`.
#' @param model A fitted model object.
#' @param X Numeric matrix (samples x features), at least 4 rows.
#' @param meta Optional per-sample metadata (data frame with `nrow(X)` rows).
#' @param n_perm Number of random row permutations to test.
#' @param tol Absolute tolerance for score agreement (default 1e-8).
#' @param n_singletons Number of distinct rows to test as singletons (always
#'   includes the first and last row). Pass `nrow(X)` for an exhaustive sweep.
#' @param model_digest Optional `function(model)` returning a comparable snapshot
#'   of the model's mutable state. Defaults to a serialization snapshot; supply a
#'   custom closure for torch / external-pointer models that do not serialize
#'   deterministically.
#' @param seed Random seed for reproducible permutations/subsets.
#' @return Invisibly `TRUE` on success; otherwise an error naming the failing
#'   case.
#' @seealso [paper3_is_row_equivariant()]
#' @export
paper3_assert_row_equivariant <- function(score_fun, model, X, meta = NULL,
                                          n_perm = 4L, tol = 1e-8,
                                          n_singletons = 4L,
                                          model_digest = NULL, seed = 1L) {
  if (!is.function(score_fun)) stop("score_fun must be a function.", call. = FALSE)
  X <- .paper3_as_score_matrix(X)
  n <- nrow(X)
  if (n < 4L) {
    stop("X must have at least 4 rows for the equivariance test.", call. = FALSE)
  }
  if (!is.null(meta)) {
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (nrow(meta) != n) stop("meta must have one row per row of X.", call. = FALSE)
  }
  if (is.null(model_digest)) model_digest <- .paper3_default_model_digest
  if (!is.function(model_digest)) stop("model_digest must be a function.", call. = FALSE)

  old_seed <- .paper3_get_seed()
  on.exit(.paper3_restore_seed(old_seed), add = TRUE)
  set.seed(seed)

  before <- model_digest(model)
  base <- .paper3_num_scores(score_fun(model, X, meta), n, "full batch")

  # (c) row-order permutations
  for (i in seq_len(n_perm)) {
    .paper3_assert_batch(score_fun, model, X, meta, base,
                         src = sample.int(n), extreme = logical(n),
                         tol = tol, where = sprintf("permutation %d", i))
  }

  # (a) singletons (always incl. first and last) - the decisive single-sample case
  n_s <- max(2L, min(as.integer(n_singletons), n))
  singles <- unique(c(1L, n, sample.int(n, n_s)))
  for (i in singles) {
    .paper3_assert_batch(score_fun, model, X, meta, base,
                         src = i, extreme = FALSE,
                         tol = tol, where = sprintf("singleton {%d}", i))
  }

  # (b) differently-composed batches (varied size) - also breaks mean/quantile
  #     centering coupling that permutations cannot detect
  sizes <- unique(c(2L, max(1L, n %/% 3L), n - 1L))
  sizes <- sizes[sizes >= 1L & sizes <= n]
  for (sz in sizes) {
    src <- sort(sample.int(n, sz))
    .paper3_assert_batch(score_fun, model, X, meta, base,
                         src = src, extreme = logical(length(src)),
                         tol = tol, where = sprintf("subset size %d", sz))
  }

  # (c) duplicated rows present - each occurrence must score identically
  d <- sample.int(n, 1L)
  src_dup <- c(seq_len(n), d, d)
  .paper3_assert_batch(score_fun, model, X, meta, base,
                       src = src_dup, extreme = logical(length(src_dup)),
                       tol = tol, where = sprintf("duplicated row {%d}", d))

  # (b) adversarial batch: target specimen among extreme (scaled) neighbours
  t <- sample.int(n, 1L)
  others <- sample.int(n, min(n, 3L))
  src_adv <- c(t, t, t, others)
  extreme_adv <- c(FALSE, TRUE, TRUE, rep(FALSE, length(others)))
  .paper3_assert_batch(score_fun, model, X, meta, base,
                       src = src_adv, extreme = extreme_adv,
                       tol = tol, where = sprintf("adversarial batch (target %d)", t))

  # (d) no model-state mutation across scoring
  after <- model_digest(model)
  if (!identical(before, after)) {
    stop(paste0("single-sample gate violated: model state changed during ",
                "scoring (e.g. BatchNorm running-stat update or online ",
                "recalibration)."), call. = FALSE)
  }
  invisible(TRUE)
}

#' Test whether a scorer is row-equivariant (non-throwing)
#'
#' Convenience wrapper around [paper3_assert_row_equivariant()] returning a
#' logical instead of erroring. Used both for affirmative checks and for the
#' negative test that a deliberately batch-coupled scorer (e.g. the batchwise
#' Sinkhorn plan) is *not* row-equivariant.
#'
#' @inheritParams paper3_assert_row_equivariant
#' @return `TRUE` if row-equivariant, `FALSE` otherwise.
#' @seealso [paper3_assert_row_equivariant()]
#' @export
paper3_is_row_equivariant <- function(score_fun, model, X, meta = NULL,
                                      n_perm = 4L, tol = 1e-8,
                                      n_singletons = 4L, model_digest = NULL,
                                      seed = 1L) {
  ok <- tryCatch({
    paper3_assert_row_equivariant(score_fun, model, X, meta, n_perm, tol,
                                  n_singletons, model_digest, seed)
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}

# Score the batch X[src, ] (extreme positions scaled to stress-test the target)
# and assert every NON-extreme ("anchored") position equals its standalone base
# score. Extreme positions are noise neighbours whose own scores are not checked.
.paper3_assert_batch <- function(score_fun, model, X, meta, base, src, extreme,
                                 tol, where) {
  if (length(extreme) == 1L) extreme <- rep(extreme, length(src))
  Xb <- X[src, , drop = FALSE]
  if (any(extreme)) Xb[extreme, ] <- Xb[extreme, , drop = FALSE] * 1e6
  mb <- .paper3_subset_meta(meta, src)
  sb <- .paper3_num_scores(score_fun(model, Xb, mb), length(src), where)
  anchored <- !extreme
  .paper3_check_close(sb[anchored], base[src[anchored]], tol, where)
}

.paper3_default_model_digest <- function(model) {
  # serialize() is blind to external-pointer (torch/reticulate) state: it writes
  # a placeholder for an externalptr, so two models differing only in their
  # libtorch C++ buffer (e.g. a BatchNorm running-stat update) would produce
  # identical bytes and silently pass clause (d). Refuse loudly instead, so each
  # external-pointer method must supply a real mutable-state model_digest.
  if (.paper3_contains_externalptr(model)) {
    stop(paste0("default model_digest cannot snapshot external-pointer state ",
                "(torch/reticulate model): pass a custom model_digest closure ",
                "capturing the model's mutable parameters as R objects (e.g. ",
                "the torch state_dict as numeric vectors)."), call. = FALSE)
  }
  tryCatch(serialize(model, connection = NULL),
           error = function(e)
             stop(paste0("default model_digest could not serialize the model (",
                         conditionMessage(e), "); pass a custom model_digest ",
                         "closure capturing the model's mutable state (e.g. ",
                         "torch parameters as R vectors)."), call. = FALSE))
}

# Recursively detect an external pointer anywhere in `model` (the pointee that
# serialize() cannot capture). Lists/pairlists, attributes, and closure/object
# environments are scanned; named system environments (global, base, package,
# namespace) are skipped, and visited environments are tracked to break cycles.
.paper3_contains_externalptr <- function(x, seen = list(), depth = 0L) {
  if (depth > 64L) return(FALSE)
  tx <- typeof(x)
  if (tx == "externalptr") return(TRUE)
  if (tx == "environment") {
    if (identical(x, emptyenv()) || identical(x, globalenv()) ||
        identical(x, baseenv()) || nzchar(environmentName(x))) {
      return(FALSE)
    }
    for (e in seen) if (identical(e, x)) return(FALSE)
    seen <- c(seen, list(x))
    for (nm in ls(x, all.names = TRUE)) {
      val <- tryCatch(get(nm, envir = x, inherits = FALSE),
                      error = function(e) NULL)
      if (.paper3_contains_externalptr(val, seen, depth + 1L)) return(TRUE)
    }
    return(FALSE)
  }
  if (tx == "closure") {
    return(.paper3_contains_externalptr(environment(x), seen, depth + 1L))
  }
  if (is.list(x) || tx == "pairlist") {
    for (el in x) {
      if (.paper3_contains_externalptr(el, seen, depth + 1L)) return(TRUE)
    }
  }
  at <- attributes(x)
  if (!is.null(at)) {
    for (a in at) {
      if (.paper3_contains_externalptr(a, seen, depth + 1L)) return(TRUE)
    }
  }
  FALSE
}

.paper3_num_scores <- function(out, n, where) {
  d <- dim(out)
  if (!is.null(d) && !(length(d) == 2L && d[2] == 1L)) {
    stop(sprintf(paste0("scorer returned an array with dim %s at %s ",
                        "(expected a numeric vector of length %d)."),
                 paste(d, collapse = "x"), where, n), call. = FALSE)
  }
  if (!is.numeric(out)) {
    stop(sprintf("scorer returned a non-numeric result (%s) at %s.",
                 class(out)[1], where), call. = FALSE)
  }
  out <- as.numeric(out)
  if (length(out) != n) {
    stop(sprintf("scorer returned length %d at %s (expected %d).",
                 length(out), where, n), call. = FALSE)
  }
  out
}

.paper3_check_close <- function(a, b, tol, where) {
  d <- max(abs(a - b))
  if (!is.finite(d) || d > tol) {
    stop(sprintf("row-equivariance violated at %s: max |delta| = %.3g > tol %.3g.",
                 where, d, tol), call. = FALSE)
  }
  invisible(TRUE)
}

.paper3_subset_meta <- function(meta, idx) {
  if (is.null(meta)) return(NULL)
  meta[idx, , drop = FALSE]
}

.paper3_get_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

.paper3_restore_seed <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", old, envir = globalenv())
  }
  invisible(NULL)
}
