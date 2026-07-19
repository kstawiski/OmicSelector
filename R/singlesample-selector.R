#' Eligible methods for the single-sample selector
#'
#' Returns the frozen tier-R1 method bank that can be considered by the
#' experimental single-sample selector. The bank contains only within-cohort
#' baselines and discriminators accepted by [deploy_singlesample()]; transfer
#' estimands, conditional routes, negative controls, and reticulate routes are
#' excluded.
#'
#' @param roster Method roster from [singlesample_method_roster()].
#'
#' @return Character vector in frozen roster order.
#' @export
singlesample_selector_candidates <- function(roster = singlesample_method_roster()) {
  required <- c("method_id", "estimand", "role", "tier")
  missing <- setdiff(required, names(roster))
  if (length(missing) > 0L) {
    stop("`roster` is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  keep <- roster$estimand == "within" &
    roster$tier == "R1" &
    roster$role %in% c("baseline", "discriminator")
  as.character(roster$method_id[keep])
}

#' Fit a leakage-resistant single-sample method selector
#'
#' @description
#' Fits an experimental meta-scorer that selects or combines registered
#' single-sample methods using grouped cross-fitted predictions from the
#' supplied training data only. The fitted object scores each new specimen from
#' frozen base models; it never uses a co-resident test batch or test outcomes.
#'
#' Three prespecified routes are available:
#' \describe{
#'   \item{`best_auc`}{Select the method with the largest mean grouped inner-fold
#'     AUC.}
#'   \item{`lower_auc_ci`}{Select the method with the largest lower 95 percent
#'     stratified group-bootstrap limit for pooled cross-fitted AUC.}
#'   \item{`simplex_stack`}{Fit a non-negative logistic stack whose weights sum
#'     to one, using cross-fitted group-level scores and a fixed numerical
#'     stabilizer. No stack hyperparameter is tuned.}
#' }
#'
#' This API establishes a leakage-resistant fitting contract; it is not a claim
#' that selection or stacking improves discrimination. Such a claim requires a
#' separate outer validation.
#'
#' @param X_train Numeric matrix or data frame (specimens x features).
#' @param y_train Binary training outcome.
#' @param methods Character method identifiers. `NULL` uses
#'   [singlesample_selector_candidates()].
#' @param route One of `"best_auc"`, `"lower_auc_ci"`, `"simplex_stack"`, or
#'   `"all"`. The `"all"` route reuses one cross-fitting pass and returns a
#'   `singlesample_selector_set` containing all three prespecified routes.
#' @param group_id Biological/provenance group per row. Rows in one group are
#'   never split across inner folds. Defaults to one group per row.
#' @param meta_train Optional training metadata, one row per specimen.
#' @param annotation Optional feature annotation passed to compatible base
#'   methods.
#' @param inner_folds Requested number of grouped inner folds. It is reduced to
#'   the smaller class-specific group count and must remain at least two.
#' @param seed Integer seed for fold construction and deterministic derived
#'   per-method seeds.
#' @param bootstrap_reps Number of stratified biological-group bootstrap draws
#'   for `lower_auc_ci`.
#' @param verify_base Whether [deploy_singlesample()] should run its row-
#'   equivariance gate for every base fit.
#' @param stack_stabilizer Fixed non-negative L2 numerical stabilizer for the
#'   simplex stack. It is not tuned from outcomes.
#'
#' @return A `singlesample_selector`, `singlesample_selector_set`, or a
#'   `singlesample_selector_ineligible` object with an explicit audit when the
#'   requested route cannot be fitted.
#' @seealso [score_singlesample_selector()], [deploy_singlesample()]
#' @export
fit_singlesample_selector <- function(
    X_train, y_train, methods = NULL,
    route = c("best_auc", "lower_auc_ci", "simplex_stack", "all"),
    group_id = NULL, meta_train = NULL, annotation = NULL,
    inner_folds = 5L, seed = 1L, bootstrap_reps = 1000L,
    verify_base = TRUE, stack_stabilizer = 1e-6) {
  route <- match.arg(route)
  X_train <- .singlesample_deploy_as_score_matrix(X_train)
  y_train <- .os_binary_y(y_train)
  n <- nrow(X_train)
  if (length(y_train) != n) {
    stop("`y_train` must have one value per row of `X_train`.", call. = FALSE)
  }
  if (length(unique(y_train)) != 2L) {
    stop("`y_train` must contain both outcome classes.", call. = FALSE)
  }
  if (is.null(group_id)) group_id <- seq_len(n)
  if (length(group_id) != n || anyNA(group_id)) {
    stop("`group_id` must contain one non-missing value per training row.",
         call. = FALSE)
  }
  group_id <- as.character(group_id)
  if (!is.null(meta_train)) {
    meta_train <- .singlesample_deploy_as_meta(
      meta_train, n, "fit_singlesample_selector", "meta_train"
    )
  }
  inner_folds <- .ss_selector_scalar_integer(inner_folds, "inner_folds", 2L)
  seed <- .ss_selector_scalar_integer(seed, "seed", 0L)
  bootstrap_reps <- .ss_selector_scalar_integer(
    bootstrap_reps, "bootstrap_reps", 20L
  )
  if (!is.numeric(stack_stabilizer) || length(stack_stabilizer) != 1L ||
      is.na(stack_stabilizer) || !is.finite(stack_stabilizer) ||
      stack_stabilizer < 0) {
    stop("`stack_stabilizer` must be one finite non-negative number.",
         call. = FALSE)
  }

  roster <- singlesample_method_roster()
  eligible_bank <- singlesample_selector_candidates(roster)
  if (is.null(methods)) methods <- eligible_bank
  if (!is.character(methods) || length(methods) == 0L || anyNA(methods) ||
      any(!nzchar(methods)) || anyDuplicated(methods)) {
    stop("`methods` must be unique non-empty method identifiers.",
         call. = FALSE)
  }
  invalid <- setdiff(methods, eligible_bank)
  if (length(invalid) > 0L) {
    stop("The selector accepts only tier-R1 within-cohort baselines and ",
         "discriminators. Ineligible method(s): ",
         paste(invalid, collapse = ", "), call. = FALSE)
  }
  roster_order <- match(methods, roster$method_id)
  methods <- methods[order(roster_order)]
  roster_order <- sort(roster_order)

  group_labels <- vapply(split(y_train, group_id), function(z) {
    u <- unique(z)
    if (length(u) != 1L) NA_integer_ else u[[1L]]
  }, integer(1L))
  if (anyNA(group_labels)) {
    stop("A biological/provenance group contains both outcome classes.",
         call. = FALSE)
  }
  class_group_counts <- table(group_labels)
  actual_folds <- min(inner_folds, min(class_group_counts))
  if (!is.finite(actual_folds) || actual_folds < 2L) {
    return(.ss_selector_ineligible(
      route, methods, "fewer_than_two_groups_in_an_outcome_class",
      data.frame(), requested_folds = inner_folds,
      actual_folds = as.integer(actual_folds), seed = seed,
      fit_accounting = list(
        inner_candidate_fits = 0L, final_candidate_refits = 0L,
        total_selector_fits = 0L
      )
    ))
  }
  actual_folds <- as.integer(actual_folds)
  folds <- os_make_grouped_stratified_folds(
    y_train, group_id = group_id, n_folds = actual_folds, seed = seed
  )
  fold_validation <- os_validate_folds(
    folds, y = y_train, group_id = group_id
  )
  if (!all(fold_validation$pass)) {
    stop("Grouped inner folds failed validation.", call. = FALSE)
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  oof <- matrix(NA_real_, nrow = n, ncol = length(methods),
                dimnames = list(rownames(X_train), methods))
  audit_rows <- vector("list", length(methods))
  fold_directions <- vector("list", length(methods))
  names(fold_directions) <- methods

  for (j in seq_along(methods)) {
    method <- methods[[j]]
    failures <- character()
    directions <- numeric(actual_folds)
    fold_auc <- rep(NA_real_, actual_folds)
    elapsed <- system.time({
      for (fold_index in seq_along(folds)) {
        test <- folds[[fold_index]]
        train <- setdiff(seq_len(n), test)
        set.seed(.ss_selector_derived_seed(seed, roster_order[[j]], fold_index))
        fitted <- tryCatch(
          deploy_singlesample(
            X_train[train, , drop = FALSE], y_train[train], method = method,
            meta_train = if (is.null(meta_train)) NULL else
              meta_train[train, , drop = FALSE],
            annotation = annotation, verify = verify_base
          ),
          error = identity
        )
        if (inherits(fitted, "error")) {
          failures <- c(failures, paste0(
            "fold", fold_index, ":fit:",
            .singlesample_deploy_safe_error_detail(conditionMessage(fitted))
          ))
          next
        }
        train_score <- tryCatch(
          score_specimen(
            fitted, X_train[train, , drop = FALSE],
            if (is.null(meta_train)) NULL else meta_train[train, , drop = FALSE]
          ),
          error = identity
        )
        test_score <- tryCatch(
          score_specimen(
            fitted, X_train[test, , drop = FALSE],
            if (is.null(meta_train)) NULL else meta_train[test, , drop = FALSE]
          ),
          error = identity
        )
        if (inherits(train_score, "error") || inherits(test_score, "error")) {
          err <- if (inherits(train_score, "error")) train_score else test_score
          failures <- c(failures, paste0(
            "fold", fold_index, ":score:",
            .singlesample_deploy_safe_error_detail(conditionMessage(err))
          ))
          next
        }
        if (length(train_score) != length(train) ||
            length(test_score) != length(test) ||
            any(!is.finite(train_score)) || any(!is.finite(test_score))) {
          failures <- c(failures, paste0("fold", fold_index,
                                          ":non_finite_or_wrong_length"))
          next
        }
        direction <- .ss_selector_training_direction(
          y_train[train], train_score, group_id[train]
        )
        directions[[fold_index]] <- direction
        oof[test, j] <- direction * as.numeric(test_score)
        fold_auc[[fold_index]] <- .ss_selector_group_auc(
          y_train[test], oof[test, j], group_id[test]
        )
      }
    })[["elapsed"]]
    complete <- length(failures) == 0L && all(is.finite(oof[, j])) &&
      all(is.finite(fold_auc))
    pooled_auc <- if (complete) {
      .ss_selector_group_auc(y_train, oof[, j], group_id)
    } else NA_real_
    lower_ci <- if (complete) {
      .ss_selector_bootstrap_auc_lower(
        y_train, oof[, j], group_id, bootstrap_reps,
        .ss_selector_derived_seed(seed, roster_order[[j]], 100003L)
      )
    } else NA_real_
    audit_rows[[j]] <- data.frame(
      method_id = method,
      roster_order = roster_order[[j]],
      complete = complete,
      expected_inner_rows = n,
      observed_inner_rows = sum(is.finite(oof[, j])),
      mean_fold_auc = if (complete) mean(fold_auc) else NA_real_,
      pooled_group_auc = pooled_auc,
      lower_group_auc_ci = lower_ci,
      elapsed_seconds = unname(elapsed),
      reason = if (complete) "eligible" else paste(unique(failures), collapse = " | "),
      stringsAsFactors = FALSE
    )
    fold_directions[[j]] <- directions
  }
  audit <- do.call(rbind, audit_rows)
  complete_methods <- audit$method_id[audit$complete]
  if (length(complete_methods) == 0L) {
    return(.ss_selector_ineligible(
      route, methods, "no_complete_inner_candidate", audit,
      requested_folds = inner_folds, actual_folds = actual_folds, seed = seed,
      fold_validation = fold_validation,
      fit_accounting = list(
        inner_candidate_fits = length(methods) * actual_folds,
        final_candidate_refits = 0L,
        total_selector_fits = length(methods) * actual_folds
      )
    ))
  }

  candidate_rows <- audit[audit$complete, , drop = FALSE]
  best_method <- candidate_rows$method_id[
    order(-candidate_rows$mean_fold_auc, candidate_rows$roster_order)[[1L]]
  ]
  lower_method <- candidate_rows$method_id[
    order(-candidate_rows$lower_group_auc_ci,
          candidate_rows$roster_order)[[1L]]
  ]
  stack_requested <- route %in% c("simplex_stack", "all")
  stack_methods <- character()
  stack <- NULL
  stack_reason <- "not_requested"
  if (stack_requested) {
    collapsed <- .ss_selector_collapse_matrix(
      y_train, oof[, complete_methods, drop = FALSE], group_id
    )
    centers <- colMeans(collapsed$scores)
    scales <- apply(collapsed$scores, 2L, stats::sd)
    stack_methods <- names(scales)[is.finite(scales) & scales > 0]
    audit$stack_eligible <- audit$method_id %in% stack_methods
    stack_reason <- NULL
    if (length(stack_methods) < 2L) {
      stack_reason <- "fewer_than_two_nonconstant_complete_candidates"
    } else {
      centers <- centers[stack_methods]
      scales <- scales[stack_methods]
      Z <- sweep(collapsed$scores[, stack_methods, drop = FALSE], 2L,
                 centers, "-")
      Z <- sweep(Z, 2L, scales, "/")
      stack_fit <- tryCatch(
        .ss_selector_fit_simplex_stack(
          Z, collapsed$y, stabilizer = stack_stabilizer
        ),
        error = identity
      )
      if (inherits(stack_fit, "error")) {
        stack_reason <- paste0(
          "simplex_stack_fit_failed:",
          .singlesample_deploy_safe_error_detail(conditionMessage(stack_fit))
        )
      } else {
        stack <- stack_fit
        names(stack$weights) <- stack_methods
        stack$center <- centers
        stack$scale <- scales
      }
    }
  } else {
    audit$stack_eligible <- NA
  }
  if (route == "simplex_stack" && !is.null(stack_reason)) {
    return(.ss_selector_ineligible(
      route, methods, stack_reason, audit,
      requested_folds = inner_folds, actual_folds = actual_folds, seed = seed,
      fold_validation = fold_validation,
      fit_accounting = list(
        inner_candidate_fits = length(methods) * actual_folds,
        final_candidate_refits = 0L,
        total_selector_fits = length(methods) * actual_folds
      )
    ))
  }
  specs <- list(
    best_auc = list(methods = best_method, stack = NULL),
    lower_auc_ci = list(methods = lower_method, stack = NULL),
    simplex_stack = if (is.null(stack_reason)) {
      list(methods = stack_methods, stack = stack)
    } else NULL
  )
  requested_routes <- if (route == "all") names(specs) else route
  selected_methods <- unique(unlist(lapply(
    specs[requested_routes], function(x) if (is.null(x)) character() else x$methods
  ), use.names = FALSE))
  selected_methods <- selected_methods[
    order(match(selected_methods, roster$method_id))
  ]
  fit_accounting <- list(
    inner_candidate_fits = length(methods) * actual_folds,
    final_candidate_refits = length(selected_methods),
    total_selector_fits = length(methods) * actual_folds +
      length(selected_methods)
  )

  final_models <- vector("list", length(selected_methods))
  final_directions <- numeric(length(selected_methods))
  names(final_models) <- names(final_directions) <- selected_methods
  final_failures <- stats::setNames(character(), character())
  for (j in seq_along(selected_methods)) {
    method <- selected_methods[[j]]
    method_order <- match(method, roster$method_id)
    set.seed(.ss_selector_derived_seed(seed, method_order, 200003L))
    fitted <- tryCatch(
      deploy_singlesample(
        X_train, y_train, method = method, meta_train = meta_train,
        annotation = annotation, verify = verify_base
      ),
      error = identity
    )
    if (inherits(fitted, "error")) {
      final_failures[[method]] <- paste0(
        "fit:", .singlesample_deploy_safe_error_detail(conditionMessage(fitted))
      )
      next
    }
    train_score <- tryCatch(
      score_specimen(fitted, X_train, meta_train), error = identity
    )
    if (inherits(train_score, "error") || length(train_score) != n ||
        any(!is.finite(train_score))) {
      detail <- if (inherits(train_score, "error")) {
        .singlesample_deploy_safe_error_detail(conditionMessage(train_score))
      } else "non_finite_or_wrong_length"
      final_failures[[method]] <- paste0("score:", detail)
      next
    }
    final_models[[j]] <- fitted
    final_directions[[j]] <- .ss_selector_training_direction(
      y_train, train_score, group_id
    )
  }
  build_one <- function(route_name) {
    spec <- specs[[route_name]]
    if (is.null(spec)) {
      return(.ss_selector_ineligible(
        route_name, methods, stack_reason, audit,
        requested_folds = inner_folds, actual_folds = actual_folds, seed = seed,
        fold_validation = fold_validation,
        fit_accounting = fit_accounting
      ))
    }
    route_failures <- intersect(spec$methods, names(final_failures))
    if (length(route_failures)) {
      details <- paste0(route_failures, ":", final_failures[route_failures])
      return(.ss_selector_ineligible(
        route_name, methods,
        paste0("final_refit_failed:", paste(details, collapse = " | ")),
        audit, requested_folds = inner_folds, actual_folds = actual_folds,
        seed = seed, fold_validation = fold_validation,
        fit_accounting = fit_accounting
      ))
    }
    .ss_selector_build_object(
      route = route_name, methods = methods, complete_methods = complete_methods,
      selected_methods = spec$methods,
      final_models = final_models[spec$methods],
      final_directions = final_directions[spec$methods], stack = spec$stack,
      audit = audit, n = n, group_id = group_id, inner_folds = inner_folds,
      actual_folds = actual_folds, seed = seed,
      bootstrap_reps = bootstrap_reps, fold_validation = fold_validation,
      stack_stabilizer = stack_stabilizer,
      fit_accounting = fit_accounting
    )
  }
  if (route == "all") {
    selectors <- stats::setNames(lapply(names(specs), build_one), names(specs))
    return(structure(
      list(
        status = if (all(vapply(selectors, function(x) x$status == "eligible",
                                logical(1L)))) "eligible" else "partial",
        route = "all",
        selectors = selectors,
        candidate_audit = audit,
        methods_requested = methods,
        methods_complete = complete_methods,
        requested_inner_folds = inner_folds,
        actual_inner_folds = actual_folds,
        seed = seed,
        bootstrap_reps = bootstrap_reps,
        fold_validation = fold_validation,
        fit_accounting = fit_accounting
      ),
      class = "singlesample_selector_set"
    ))
  }
  build_one(route)
}

.ss_selector_build_object <- function(
    route, methods, complete_methods, selected_methods, final_models,
    final_directions, stack, audit, n, group_id, inner_folds, actual_folds,
    seed, bootstrap_reps, fold_validation, stack_stabilizer,
    fit_accounting) {
  structure(
    list(
      status = "eligible",
      route = route,
      methods_requested = methods,
      methods_complete = complete_methods,
      selected_methods = selected_methods,
      models = final_models,
      directions = final_directions,
      stack = stack,
      candidate_audit = audit,
      n_train = n,
      n_groups = length(unique(group_id)),
      requested_inner_folds = inner_folds,
      actual_inner_folds = actual_folds,
      seed = seed,
      bootstrap_reps = bootstrap_reps,
      fold_validation = fold_validation,
      fit_accounting = fit_accounting,
      contract = list(
        candidate_tier = "R1",
        estimand = "within",
        score_orientation = "training_only",
        test_batch_required = FALSE,
        stack_stabilizer = if (route == "simplex_stack")
          stack_stabilizer else NA_real_
      )
    ),
    class = "singlesample_selector"
  )
}

#' Score specimens with a frozen single-sample selector
#'
#' @param selector A fitted `singlesample_selector` from
#'   [fit_singlesample_selector()].
#' @param x_new Named numeric vector, matrix, or data frame of new specimens.
#' @param meta Optional metadata, one row per new specimen.
#'
#' @return Numeric discrimination score per new specimen. For a
#'   `singlesample_selector_set`, a data frame with one column per route is
#'   returned and an explicitly ineligible route is `NA`. Larger values indicate
#'   the positive outcome class according to the training-frozen orientation.
#' @export
score_singlesample_selector <- function(selector, x_new, meta = NULL) {
  if (inherits(selector, "singlesample_selector_set")) {
    scores <- lapply(selector$selectors, function(x) {
      if (inherits(x, "singlesample_selector_ineligible")) {
        X <- .singlesample_deploy_as_score_matrix(x_new)
        rep(NA_real_, nrow(X))
      } else score_singlesample_selector(x, x_new, meta)
    })
    return(as.data.frame(scores, check.names = FALSE,
                         stringsAsFactors = FALSE))
  }
  if (inherits(selector, "singlesample_selector_ineligible")) {
    stop("Cannot score an ineligible single-sample selector: ",
         selector$reason, call. = FALSE)
  }
  if (!inherits(selector, "singlesample_selector")) {
    stop("`selector` must be a fitted `singlesample_selector` object.",
         call. = FALSE)
  }
  X <- .singlesample_deploy_as_score_matrix(x_new)
  if (!is.null(meta)) {
    meta <- .singlesample_deploy_as_meta(
      meta, nrow(X), "score_singlesample_selector", "meta"
    )
  }
  base <- vapply(selector$selected_methods, function(method) {
    as.numeric(score_specimen(selector$models[[method]], X, meta)) *
      selector$directions[[method]]
  }, numeric(nrow(X)))
  base <- matrix(base, nrow = nrow(X), ncol = length(selector$selected_methods),
                 dimnames = list(NULL, selector$selected_methods))
  if (selector$route != "simplex_stack") {
    return(as.numeric(base[, 1L]))
  }
  Z <- sweep(base[, selector$selected_methods, drop = FALSE], 2L,
             selector$stack$center, "-")
  Z <- sweep(Z, 2L, selector$stack$scale, "/")
  as.numeric(selector$stack$intercept + Z %*% selector$stack$weights)
}

#' Print a single-sample selector
#'
#' @param x A `singlesample_selector` or
#'   `singlesample_selector_ineligible` object.
#' @param ... Unused.
#'
#' @return Invisibly, `x`.
#' @export
print.singlesample_selector <- function(x, ...) {
  cat("OmicSelector experimental single-sample selector\n")
  cat("  status: ", x$status, "\n", sep = "")
  cat("  route: ", x$route, "\n", sep = "")
  cat("  selected: ", paste(x$selected_methods, collapse = ", "), "\n",
      sep = "")
  cat("  caveat: fitting contract only; superiority requires external outer ",
      "validation.\n", sep = "")
  invisible(x)
}

#' @param x A `singlesample_selector_set` object.
#' @param ... Unused.
#' @rdname fit_singlesample_selector
#' @export
print.singlesample_selector_set <- function(x, ...) {
  cat("OmicSelector experimental single-sample selector set\n")
  cat("  status: ", x$status, "\n", sep = "")
  for (route in names(x$selectors)) {
    selector <- x$selectors[[route]]
    selected <- if (selector$status == "eligible") {
      paste(selector$selected_methods, collapse = ", ")
    } else paste0("ineligible (", selector$reason, ")")
    cat("  ", route, ": ", selected, "\n", sep = "")
  }
  cat("  caveat: fitting contract only; superiority requires external outer ",
      "validation.\n", sep = "")
  invisible(x)
}

#' @rdname fit_singlesample_selector
#' @export
print.singlesample_selector_ineligible <- function(x, ...) {
  cat("OmicSelector experimental single-sample selector\n")
  cat("  status: ineligible\n")
  cat("  route: ", x$route, "\n", sep = "")
  cat("  reason: ", x$reason, "\n", sep = "")
  invisible(x)
}

.ss_selector_scalar_integer <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < minimum) {
    stop("`", name, "` must be one integer >= ", minimum, ".",
         call. = FALSE)
  }
  as.integer(x)
}

.ss_selector_derived_seed <- function(seed, method_order, fold_index) {
  modulus <- .Machine$integer.max - 1
  value <- (as.double(seed) + 104729 * as.double(method_order) +
              1009 * as.double(fold_index)) %% modulus
  as.integer(value + 1)
}

.ss_selector_group_auc <- function(y, score, group_id) {
  collapsed <- .ss_selector_collapse_matrix(
    y, matrix(score, ncol = 1L), group_id
  )
  .auc_fast(collapsed$y, collapsed$scores[, 1L])
}

.ss_selector_training_direction <- function(y, score, group_id) {
  auc <- .ss_selector_group_auc(y, score, group_id)
  if (!is.finite(auc)) {
    stop("Training data do not support a finite fixed-direction AUC.",
         call. = FALSE)
  }
  if (auc < 0.5) -1 else 1
}

.ss_selector_collapse_matrix <- function(y, scores, group_id) {
  scores <- as.matrix(scores)
  groups <- unique(as.character(group_id))
  y_group <- integer(length(groups))
  score_group <- matrix(NA_real_, nrow = length(groups), ncol = ncol(scores),
                        dimnames = list(groups, colnames(scores)))
  for (i in seq_along(groups)) {
    idx <- which(as.character(group_id) == groups[[i]])
    values <- unique(y[idx])
    if (length(values) != 1L) {
      stop("A biological/provenance group contains both outcome classes.",
           call. = FALSE)
    }
    y_group[[i]] <- values[[1L]]
    score_group[i, ] <- colMeans(scores[idx, , drop = FALSE])
  }
  list(y = y_group, scores = score_group, group_id = groups)
}

.ss_selector_bootstrap_auc_lower <- function(y, score, group_id, reps, seed) {
  collapsed <- .ss_selector_collapse_matrix(
    y, matrix(score, ncol = 1L), group_id
  )
  values <- collapsed$scores[, 1L]
  pos <- which(collapsed$y == 1L)
  neg <- which(collapsed$y == 0L)
  if (length(pos) < 2L || length(neg) < 2L) return(NA_real_)
  set.seed(seed)
  draws <- replicate(reps, {
    idx <- c(sample(pos, length(pos), replace = TRUE),
             sample(neg, length(neg), replace = TRUE))
    .auc_fast(collapsed$y[idx], values[idx])
  })
  unname(stats::quantile(draws, probs = 0.025, names = FALSE,
                         na.rm = TRUE, type = 6L))
}

.ss_selector_fit_simplex_stack <- function(Z, y, stabilizer) {
  Z <- as.matrix(Z)
  k <- ncol(Z)
  if (k < 2L) stop("A simplex stack requires at least two score columns.",
                    call. = FALSE)
  objective <- function(theta) {
    weights <- c(theta[-1L], 1 - sum(theta[-1L]))
    eta <- as.numeric(theta[[1L]] + Z %*% weights)
    logloss <- mean(pmax(eta, 0) + log1p(exp(-abs(eta))) - y * eta)
    logloss + stabilizer * sum(weights ^ 2)
  }
  gradient <- function(theta) {
    weights <- c(theta[-1L], 1 - sum(theta[-1L]))
    eta <- as.numeric(theta[[1L]] + Z %*% weights)
    residual <- stats::plogis(eta) - y
    contrasts <- sweep(Z[, seq_len(k - 1L), drop = FALSE], 1L,
                       Z[, k], "-")
    c(
      mean(residual),
      colMeans(contrasts * residual) +
        2 * stabilizer * (weights[seq_len(k - 1L)] - weights[[k]])
    )
  }
  theta <- c(stats::qlogis(mean(y)), rep(1 / k, k - 1L))
  ui <- matrix(0, nrow = k, ncol = k)
  if (k > 1L) ui[cbind(seq_len(k - 1L), 1L + seq_len(k - 1L))] <- 1
  ui[k, -1L] <- -1
  ci <- c(rep(0, k - 1L), -1)
  fit <- stats::constrOptim(
    theta = theta, f = objective, grad = gradient, ui = ui, ci = ci,
    method = "BFGS", control = list(maxit = 1000L),
    outer.iterations = 100L
  )
  if (!identical(as.integer(fit$convergence), 0L) ||
      any(!is.finite(fit$par)) || !is.finite(fit$value)) {
    stop("The simplex optimizer did not converge to a finite solution.",
         call. = FALSE)
  }
  weights <- c(fit$par[-1L], 1 - sum(fit$par[-1L]))
  weights[abs(weights) < 1e-12] <- 0
  weights <- pmax(weights, 0)
  weights <- weights / sum(weights)
  list(
    intercept = unname(fit$par[[1L]]),
    weights = weights,
    convergence = fit$convergence,
    objective = fit$value,
    stabilizer = stabilizer
  )
}

.ss_selector_ineligible <- function(route, methods, reason, audit,
                                    requested_folds, actual_folds, seed,
                                    fold_validation = NULL,
                                    fit_accounting = NULL) {
  structure(
    list(
      status = "ineligible",
      route = route,
      reason = reason,
      methods_requested = methods,
      methods_complete = if (nrow(audit))
        audit$method_id[audit$complete] else character(),
      selected_methods = character(),
      candidate_audit = audit,
      requested_inner_folds = requested_folds,
      actual_inner_folds = actual_folds,
      seed = seed,
      fold_validation = fold_validation,
      fit_accounting = fit_accounting
    ),
    class = c("singlesample_selector_ineligible", "singlesample_selector")
  )
}
