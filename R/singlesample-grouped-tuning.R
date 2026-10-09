# Training-only group guards and split builders. NULL branches remain in each
# caller so legacy fits retain their exact split and RNG behaviour.
.ss_validate_groups <- function(group_train, n, caller) {
  if (is.null(group_train)) return(NULL)
  if (!is.atomic(group_train) || !is.null(dim(group_train)) ||
      length(n) != 1L || length(group_train) != n || anyNA(group_train)) {
    .ss_input_reject("INVALID_GROUPS",
                     "group_train must contain one non-missing identifier per training row",
                     caller)
  }
  as.character(group_train)
}

# The shipped replay helpers may receive an older external srb fit adapter.
# Keep that adapter for NULL groups; otherwise call the registered package fit
# through the same declared-formals adapter used by deployment.
.ss_fit_benchmark_route <- function(fit_call, fit_fn_name, X_train, y_train,
                                    meta_train, hp, group_train = NULL) {
  if (is.null(group_train)) {
    return(fit_call(fit_fn_name, X_train, y_train, meta_train, hp))
  }
  group_train <- .ss_validate_groups(group_train, nrow(X_train), "benchmark fit")
  if ("group_train" %in% names(formals(fit_call))) {
    return(fit_call(fit_fn_name, X_train, y_train, meta_train, hp,
                    group_train = group_train))
  }
  fit <- get(fit_fn_name, envir = asNamespace("OmicSelector"),
             mode = "function", inherits = FALSE)
  if (!"group_train" %in% names(formals(fit))) {
    return(fit_call(fit_fn_name, X_train, y_train, meta_train, hp))
  }
  .singlesample_deploy_fit(fit, fit_fn_name, fit_fn_name, X_train, y_train,
                           meta_train, annotation = NULL,
                           group_train = group_train, hp = hp)
}

.ss_grouped_inner_folds <- function(y, group_train, k, seed = 1L,
                                    min_folds = 2L) {
  caller <- "grouped inner resampling"
  groups <- .ss_validate_groups(group_train, length(y), caller)
  k <- min(as.integer(k), length(unique(groups)))
  if (k < min_folds) {
    .ss_input_reject("INVALID_GROUPS",
                     paste0("inner resampling requires at least ", min_folds,
                            " distinct groups"), caller)
  }
  folds <- os_make_grouped_stratified_folds(
    y, group_id = match(groups, unique(groups)), n_folds = k, seed = seed,
    mixed_group_action = "allow"
  )
  fold <- integer(length(y))
  for (i in seq_along(folds)) fold[folds[[i]]] <- i
  if (any(lengths(folds) == 0L) || any(vapply(seq_len(k), function(i) {
    length(unique(y[fold != i])) < 2L
  }, logical(1L)))) {
    .ss_input_reject("INVALID_GROUPS",
                     "groups cannot form non-empty folds with both outcomes in every training partition",
                     caller)
  }
  fold
}

# Unstratified calibration holdout. Choose whole groups nearest the requested
# row count; never truncate a group to hit the target exactly.
.ss_grouped_calibration_indices <- function(group_train, fraction,
                                            min_reference = 1L) {
  groups <- unique(group_train)
  if (length(groups) < 2L || !is.finite(fraction) ||
      fraction <= 0 || fraction >= 1) {
    .ss_input_reject("INVALID_GROUPS",
                     "calibration requires at least two groups and a fraction between zero and one",
                     "fit_conformal_anomaly")
  }
  groups <- groups[sample.int(length(groups))]
  sizes <- vapply(groups, function(g) sum(group_train == g), integer(1L))
  n <- length(group_train)
  prefix <- cumsum(sizes)[-length(sizes)]
  counts <- c(prefix, n - prefix)
  candidates <- which(counts <= n - min_reference)
  if (!length(candidates)) {
    .ss_input_reject("INVALID_GROUPS",
                     "grouped calibration leaves fewer than k reference rows",
                     "fit_conformal_anomaly")
  }
  chosen <- candidates[which.min(abs(counts[candidates] - fraction * n))]
  take <- if (chosen <= length(prefix)) {
    seq_len(chosen)
  } else {
    seq.int(chosen - length(prefix) + 1L, length(groups))
  }
  which(group_train %in% groups[take])
}

# Keep every held-out environment intact and exclude its biological groups
# from the training complement, including groups observed in other environments.
.ss_grouped_environment_splits <- function(env, group_train) {
  groups <- .ss_validate_groups(group_train, length(env), "fit_sel_stablemate")
  out <- lapply(sort(unique(as.character(env))), function(e) {
    test <- which(env == e)
    train <- which(!groups %in% groups[test])
    if (!length(train)) {
      .ss_input_reject("INVALID_GROUPS",
                       "an environment holdout leaves no independent training groups",
                       "fit_sel_stablemate")
    }
    list(train = train, test = test)
  })
  out
}

# Prototypical episodes use whole groups on either the support or query side.
# Row-count targets can be exceeded by the last selected group.
.ss_grouped_episode_plan <- function(idx_case, idx_ctrl, hp, group_train) {
  n <- length(idx_case) + length(idx_ctrl)
  groups <- .ss_validate_groups(group_train, n, "fit_proto_net")
  y <- integer(n)
  y[idx_case] <- 1L
  take <- function(pool, target) {
    chosen <- integer()
    for (cl in c(0L, 1L)) {
      candidates <- unique(groups[pool[y[pool] == cl]])
      candidates <- candidates[sample.int(length(candidates))]
      count <- 0L
      for (g in candidates) {
        rows <- pool[groups[pool] == g]
        chosen <- union(chosen, rows)
        count <- count + sum(y[rows] == cl)
        if (count >= target) break
      }
    }
    sort(chosen)
  }
  lapply(seq_len(hp$epochs), function(ep) {
    seed <- as.integer((as.double(hp$seed) + ep) %% .Machine$integer.max)
    fold <- .ss_grouped_inner_folds(y, groups, k = 2L, seed = seed)
    support <- take(which(fold == 1L), hp$n_support)
    query <- take(which(fold == 2L), hp$n_query)
    list(support = support, query = query, query_label = y[query])
  })
}

# StableMate's weighted squared-error stability objective, with biological
# groups excluded from each environment's training complement. In-sample terms
# follow obj_pss; only the held-out training indices change.
.sel_stablemate_grouped_pss <- function(Y, X, env, pred_in, splits) {
  env_tab <- table(env)
  w <- as.numeric(1 / env_tab[as.character(env)])
  w <- w / sum(w) * length(Y)
  design <- cbind(1, X[, pred_in == 1L, drop = FALSE])
  fit_predict <- function(train, test) {
    if (sum(pred_in) == 0L) {
      return(rep(sum(w[train] * Y[train]) / sum(w[train]), length(test)))
    }
    mod <- stats::lm.wfit(design[train, , drop = FALSE], Y[train], w[train])
    beta <- mod$coefficients
    beta[is.na(beta)] <- 0
    as.vector(design[test, , drop = FALSE] %*% beta)
  }
  pss <- 0
  for (split in splits) {
    test <- if (length(env_tab) > 2L) split$test else seq_along(Y)
    pred <- fit_predict(split$train, test)
    pss <- pss + sum(w[test] * (pred - Y[test])^2)
    if (length(env_tab) > 2L) {
      test <- split$test
      pred <- fit_predict(test, test)
      pss <- pss + sum(w[test] * (pred - Y[test])^2)
    }
  }
  pss
}
