#' @title Modern Tabular Learners for OmicSelector 2.0
#'
#' @description
#' Optional mlr3 learner wrappers for modern tabular methods commonly used as
#' strong baselines on biomarker-style classification problems. The wrappers are
#' designed for small-\eqn{P}, medium-\eqn{N}, and class-imbalanced omics tasks.
#'
#' @details
#' The learners in this file are intentionally lazy about backend availability:
#' they can be constructed without the Python or CatBoost runtime being present,
#' but training will fail with an actionable error if the required dependency is
#' missing.
#'
#' Backends:
#' - `make_tabnet_learner()`: `pytorch_tabnet` through `reticulate`
#' - `make_tabm_learner()`: `tabm` or `rtdl_revisiting_models` through `reticulate`
#' - `make_fttransformer_learner()`: `rtdl_revisiting_models` through `reticulate`
#' - `make_catboost_learner()`: `catboost` R package
#' - `make_tabpfn_learner()`: `tabpfn` through `reticulate`
#'
#' GPU selection follows `torch.cuda.is_available()` when `device = "auto"`.
#'
#' @name TabularDL
NULL


#' @keywords internal
.tabdl_require_namespace <- function(pkg, message_text = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is.null(message_text)) {
      message_text <- sprintf(
        "Package '%s' is required for this learner. Install it before training.",
        pkg
      )
    }
    stop(message_text, call. = FALSE)
  }
}


#' @keywords internal
.tabdl_require_python_module <- function(module, hint = NULL) {
  .tabdl_require_namespace(
    "reticulate",
    "Package 'reticulate' is required for Python-backed tabular learners.\nInstall with: install.packages('reticulate')"
  )

  if (!reticulate::py_module_available(module)) {
    if (is.null(hint)) {
      hint <- sprintf(
        "Python module '%s' is not available in the active reticulate environment.",
        module
      )
    }
    stop(hint, call. = FALSE)
  }

  invisible(TRUE)
}


#' @keywords internal
.tabdl_torch_cuda_available <- function() {
  if (requireNamespace("torch", quietly = TRUE)) {
    cuda_r <- tryCatch(torch::cuda_is_available(), error = function(e) FALSE)
    if (isTRUE(cuda_r)) {
      return(TRUE)
    }
  }

  if (requireNamespace("reticulate", quietly = TRUE)) {
    py_cuda <- tryCatch({
      if (!reticulate::py_module_available("torch")) {
        return(FALSE)
      }
      torch_py <- reticulate::import("torch", delay_load = TRUE)
      isTRUE(torch_py$cuda$is_available())
    }, error = function(e) FALSE)
    if (isTRUE(py_cuda)) {
      return(TRUE)
    }
  }

  FALSE
}


#' @keywords internal
.tabdl_resolve_device <- function(device = c("auto", "cpu", "cuda", "gpu")) {
  device <- match.arg(device)

  if (identical(device, "auto")) {
    return(if (.tabdl_torch_cuda_available()) "cuda" else "cpu")
  }

  if (identical(device, "gpu")) {
    return(if (.tabdl_torch_cuda_available()) "cuda" else "cpu")
  }

  if (identical(device, "cuda") && !.tabdl_torch_cuda_available()) {
    warning("CUDA requested but not available; falling back to CPU.", call. = FALSE)
    return("cpu")
  }

  device
}


#' @keywords internal
.tabdl_feature_names <- function(task) {
  weight_cols <- character()
  if (!is.null(task$col_roles$weights_learner)) {
    weight_cols <- c(weight_cols, task$col_roles$weights_learner)
  }
  if (!is.null(task$col_roles$weights_measure)) {
    weight_cols <- c(weight_cols, task$col_roles$weights_measure)
  }
  setdiff(task$feature_names, unique(weight_cols))
}


#' @keywords internal
.tabdl_feature_matrix <- function(task, feature_names = NULL) {
  if (is.null(feature_names)) {
    feature_names <- .tabdl_feature_names(task)
  }

  x <- as.data.frame(task$data(cols = feature_names), stringsAsFactors = FALSE)
  bad <- names(x)[!vapply(x, function(col) {
    is.numeric(col) || is.integer(col) || is.logical(col)
  }, logical(1))]

  if (length(bad) > 0L) {
    stop(
      sprintf(
        "These learners support numeric/integer/logical features only. Unsupported columns: %s",
        paste(utils::head(bad, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  x_mat <- data.matrix(x)
  storage.mode(x_mat) <- "double"
  x_mat
}


#' @keywords internal
.tabdl_encode_outcome <- function(task) {
  y_raw <- as.character(task$truth())
  class_names <- task$class_names

  if (anyNA(y_raw)) {
    stop("Missing target labels are not supported by these learners.", call. = FALSE)
  }

  if (length(class_names) == 2L) {
    positive_class <- class_names[[1L]]
    negative_class <- class_names[[2L]]
    y_backend <- ifelse(y_raw == positive_class, 1L, 0L)
    backend_class_order <- c(negative_class, positive_class)

    return(list(
      y_raw = y_raw,
      y_backend = as.integer(y_backend),
      class_names = class_names,
      backend_class_order = backend_class_order,
      positive_class = positive_class,
      negative_class = negative_class,
      is_binary = TRUE
    ))
  }

  y_factor <- factor(y_raw, levels = class_names)
  if (anyNA(y_factor)) {
    stop("Could not align the target labels with task$class_names.", call. = FALSE)
  }

  list(
    y_raw = y_raw,
    y_backend = as.integer(y_factor) - 1L,
    class_names = class_names,
    backend_class_order = class_names,
    positive_class = NULL,
    negative_class = NULL,
    is_binary = FALSE
  )
}


#' @keywords internal
.tabdl_extract_task_weights <- function(task) {
  weights_dt <- tryCatch(task$weights_learner, error = function(e) NULL)
  if (is.null(weights_dt) || nrow(weights_dt) == 0L) {
    return(NULL)
  }

  idx <- match(task$row_ids, weights_dt$row_id)
  weights <- weights_dt$weight[idx]
  if (all(is.na(weights))) {
    return(NULL)
  }

  as.numeric(weights)
}


#' @keywords internal
.tabdl_class_weight_vector <- function(y_raw, class_names, class_weights) {
  if (is.null(class_weights) || identical(class_weights, "none")) {
    return(NULL)
  }

  if (identical(class_weights, "balanced")) {
    counts <- table(factor(y_raw, levels = class_names))
    n <- length(y_raw)
    k <- length(class_names)
    w <- n / (k * pmax(as.numeric(counts), 1))
    names(w) <- class_names
    return(w)
  }

  if (!is.numeric(class_weights)) {
    stop(
      "`class_weights` must be 'balanced', 'none', NULL, or a numeric vector.",
      call. = FALSE
    )
  }

  if (!is.null(names(class_weights))) {
    missing_names <- setdiff(class_names, names(class_weights))
    if (length(missing_names) > 0L) {
      stop(
        sprintf(
          "Missing class weights for classes: %s",
          paste(missing_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    weights <- as.numeric(class_weights[class_names])
    names(weights) <- class_names
    return(weights)
  }

  if (length(class_weights) != length(class_names)) {
    stop(
      sprintf(
        "Unnamed `class_weights` must have length %d.",
        length(class_names)
      ),
      call. = FALSE
    )
  }

  weights <- as.numeric(class_weights)
  names(weights) <- class_names
  weights
}


#' @keywords internal
.tabdl_backend_class_weights <- function(outcome, class_weights) {
  weights <- .tabdl_class_weight_vector(
    y_raw = outcome$y_raw,
    class_names = outcome$class_names,
    class_weights = class_weights
  )

  if (is.null(weights)) {
    return(NULL)
  }

  if (isTRUE(outcome$is_binary)) {
    return(c(
      unname(weights[[outcome$negative_class]]),
      unname(weights[[outcome$positive_class]])
    ))
  }

  unname(weights[outcome$class_names])
}


#' @keywords internal
.tabdl_sample_weights <- function(outcome, task_weights = NULL, class_weights = "balanced", indices = NULL) {
  if (is.null(indices)) {
    indices <- seq_along(outcome$y_backend)
  }

  weights <- rep(1, length(indices))

  if (!is.null(task_weights)) {
    weights <- weights * as.numeric(task_weights[indices])
  }

  class_weight_vec <- .tabdl_backend_class_weights(outcome, class_weights)
  if (!is.null(class_weight_vec)) {
    weights <- weights * class_weight_vec[outcome$y_backend[indices] + 1L]
  }

  if (all(abs(weights - 1) < .Machine$double.eps^0.5, na.rm = TRUE)) {
    return(NULL)
  }

  as.numeric(weights)
}


#' @keywords internal
.tabdl_validation_split <- function(y, validation_fraction = 0.15, seed = 1L) {
  idx <- seq_along(y)

  if (length(idx) < 10L || validation_fraction <= 0) {
    return(list(train = idx, valid = integer()))
  }

  old_seed <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  valid_idx <- integer()

  for (cls in sort(unique(y))) {
    cls_idx <- idx[y == cls]
    if (length(cls_idx) < 2L) {
      next
    }
    n_valid <- max(1L, floor(length(cls_idx) * validation_fraction))
    n_valid <- min(n_valid, length(cls_idx) - 1L)
    valid_idx <- c(valid_idx, sample(cls_idx, size = n_valid, replace = FALSE))
  }

  valid_idx <- sort(unique(valid_idx))
  train_idx <- setdiff(idx, valid_idx)

  if (length(valid_idx) == 0L) {
    return(list(train = idx, valid = integer()))
  }

  if (length(unique(y[train_idx])) < length(unique(y)) ||
      length(unique(y[valid_idx])) < length(unique(y))) {
    return(list(train = idx, valid = integer()))
  }

  list(train = train_idx, valid = valid_idx)
}


#' @keywords internal
.tabdl_training_bundle <- function(task) {
  feature_names <- .tabdl_feature_names(task)
  outcome <- .tabdl_encode_outcome(task)

  c(
    list(
      x = .tabdl_feature_matrix(task, feature_names = feature_names),
      feature_names = feature_names,
      task_weights = .tabdl_extract_task_weights(task)
    ),
    outcome
  )
}


#' @keywords internal
.tabdl_binary_prob_matrix <- function(prob_positive, class_names) {
  prob_positive <- pmin(pmax(as.numeric(prob_positive), 0), 1)
  prob <- cbind(
    prob_positive,
    1 - prob_positive
  )
  colnames(prob) <- class_names
  prob
}


#' @keywords internal
.tabdl_multiclass_prob_matrix <- function(prob, class_names) {
  prob <- as.matrix(prob)
  storage.mode(prob) <- "double"

  if (ncol(prob) != length(class_names)) {
    stop(
      sprintf(
        "Expected %d probability columns, got %d.",
        length(class_names),
        ncol(prob)
      ),
      call. = FALSE
    )
  }

  colnames(prob) <- class_names
  prob
}


#' @keywords internal
.tabdl_reorder_python_prob <- function(prob, class_names, backend_classes = NULL, is_binary = TRUE) {
  if (isTRUE(is_binary)) {
    if (is.null(dim(prob))) {
      return(.tabdl_binary_prob_matrix(prob, class_names))
    }

    prob <- as.matrix(prob)
    storage.mode(prob) <- "double"

    if (!is.null(backend_classes) && length(backend_classes) == 2L) {
      pos_col <- match(1L, as.integer(backend_classes))
      if (!is.na(pos_col)) {
        return(.tabdl_binary_prob_matrix(prob[, pos_col], class_names))
      }
    }

    if (ncol(prob) >= 2L) {
      return(.tabdl_binary_prob_matrix(prob[, 2L], class_names))
    }

    .tabdl_binary_prob_matrix(prob[, 1L], class_names)
  } else {
    prob <- as.matrix(prob)
    storage.mode(prob) <- "double"
    if (!is.null(backend_classes)) {
      order_idx <- match(0:(length(class_names) - 1L), as.integer(backend_classes))
      if (!anyNA(order_idx)) {
        prob <- prob[, order_idx, drop = FALSE]
      }
    }
    .tabdl_multiclass_prob_matrix(prob, class_names)
  }
}


#' @keywords internal
.tabdl_warn_ignored_weights <- function(learner_id) {
  warning(
    sprintf(
      "%s does not expose per-row sample weighting through this wrapper; task weights are ignored.",
      learner_id
    ),
    call. = FALSE
  )
}


#' @keywords internal
.tabdl_load_python_helpers <- function() {
  .tabdl_require_namespace(
    "reticulate",
    "Package 'reticulate' is required for Python-backed tabular learners.\nInstall with: install.packages('reticulate')"
  )

  if (reticulate::py_exists("omicselector_fit_rtdl")) {
    return(invisible(TRUE))
  }

  code <- paste(
    c(
      "import copy",
      "import random",
      "import numpy as np",
      "",
      "def omicselector_make_class_weight_dict(weights):",
      "    return {int(i): float(w) for i, w in enumerate(weights)}",
      "",
      "def _omicselector_seed(seed):",
      "    random.seed(int(seed))",
      "    np.random.seed(int(seed))",
      "    try:",
      "        import torch",
      "        torch.manual_seed(int(seed))",
      "        if torch.cuda.is_available():",
      "            torch.cuda.manual_seed_all(int(seed))",
      "        if hasattr(torch.backends, 'cudnn'):",
      "            torch.backends.cudnn.deterministic = True",
      "            torch.backends.cudnn.benchmark = False",
      "    except Exception:",
      "        pass",
      "",
      "def _omicselector_resolve_device(device='auto'):",
      "    import torch",
      "    if device is None or device == 'auto':",
      "        return 'cuda' if torch.cuda.is_available() else 'cpu'",
      "    if device == 'gpu':",
      "        return 'cuda' if torch.cuda.is_available() else 'cpu'",
      "    if device == 'cuda' and not torch.cuda.is_available():",
      "        return 'cpu'",
      "    return device",
      "",
      "def _omicselector_tabm_class():",
      "    try:",
      "        import tabm",
      "        if hasattr(tabm, 'TabM'):",
      "            return tabm.TabM",
      "    except Exception:",
      "        pass",
      "    import rtdl_revisiting_models as rm",
      "    if hasattr(rm, 'TabM'):",
      "        return rm.TabM",
      "    if hasattr(rm, 'tabm') and hasattr(rm.tabm, 'TabM'):",
      "        return rm.tabm.TabM",
      "    raise ImportError(\"TabM is unavailable. Install the official 'tabm' package or a build of rtdl_revisiting_models that exposes TabM.\")",
      "",
      "def _omicselector_build_model(model_name, n_features, n_classes, config):",
      "    import rtdl_revisiting_models as rm",
      "    d_out = 1 if n_classes == 2 else int(n_classes)",
      "    if model_name == 'fttransformer':",
      "        if hasattr(rm.FTTransformer, 'get_default_kwargs'):",
      "            kwargs = dict(rm.FTTransformer.get_default_kwargs(n_blocks=int(config.get('n_blocks', 3))))",
      "        else:",
      "            kwargs = {}",
      "        kwargs['d_out'] = d_out",
      "        if config.get('d_block') is not None:",
      "            kwargs['d_block'] = int(config['d_block'])",
      "        if config.get('attention_dropout') is not None:",
      "            kwargs['attention_dropout'] = float(config['attention_dropout'])",
      "        if config.get('ffn_dropout') is not None:",
      "            kwargs['ffn_dropout'] = float(config['ffn_dropout'])",
      "        if config.get('ffn_d_hidden') is not None:",
      "            kwargs['ffn_d_hidden'] = int(config['ffn_d_hidden'])",
      "        return rm.FTTransformer(n_cont_features=int(n_features), cat_cardinalities=[], **kwargs)",
      "    if model_name == 'tabm':",
      "        TabM = _omicselector_tabm_class()",
      "        kwargs = {",
      "            'n_num_features': int(n_features),",
      "            'cat_cardinalities': [],",
      "            'd_out': d_out,",
      "        }",
      "        for key in ('k', 'd_block', 'dropout', 'n_blocks', 'arch_type'):",
      "            if config.get(key) is not None:",
      "                kwargs[key] = config.get(key)",
      "        if hasattr(TabM, 'make'):",
      "            return TabM.make(**kwargs)",
      "        return TabM(**kwargs)",
      "    raise ValueError(f'Unknown RTDL model: {model_name}')",
      "",
      "def _omicselector_logits_to_probs(logits, n_classes, is_tabm, torch_module):",
      "    if is_tabm and hasattr(logits, 'ndim') and logits.ndim >= 3:",
      "        logits = logits.mean(dim=1)",
      "    if n_classes == 2:",
      "        logits = logits.squeeze(-1)",
      "        p1 = torch_module.sigmoid(logits)",
      "        return p1",
      "    return torch_module.softmax(logits, dim=-1)",
      "",
      "def _omicselector_loss(logits, yb, wb, n_classes, is_tabm, torch_module):",
      "    F = torch_module.nn.functional",
      "    if n_classes == 2:",
      "        if is_tabm and hasattr(logits, 'ndim') and logits.ndim >= 3:",
      "            targets = yb[:, None, None].expand_as(logits)",
      "            loss_raw = F.binary_cross_entropy_with_logits(logits, targets, reduction='none')",
      "            loss_raw = loss_raw.mean(dim=1).squeeze(-1)",
      "        else:",
      "            loss_raw = F.binary_cross_entropy_with_logits(logits.squeeze(-1), yb, reduction='none')",
      "    else:",
      "        if is_tabm and hasattr(logits, 'ndim') and logits.ndim >= 3:",
      "            b, k, c = logits.shape",
      "            loss_raw = F.cross_entropy(",
      "                logits.reshape(b * k, c),",
      "                yb[:, None].expand(b, k).reshape(-1),",
      "                reduction='none'",
      "            ).reshape(b, k).mean(dim=1)",
      "        else:",
      "            loss_raw = F.cross_entropy(logits, yb, reduction='none')",
      "    if wb is not None:",
      "        return (loss_raw * wb).sum() / torch_module.clamp(wb.sum(), min=1e-8)",
      "    return loss_raw.mean()",
      "",
      "def omicselector_fit_rtdl(model_name, X_train, y_train, X_valid=None, y_valid=None, sample_weight=None, config=None):",
      "    import torch",
      "    config = {} if config is None else dict(config)",
      "    X_train = np.asarray(X_train, dtype=np.float32)",
      "    y_train = np.asarray(y_train)",
      "    if X_valid is not None:",
      "        X_valid = np.asarray(X_valid, dtype=np.float32)",
      "        y_valid = np.asarray(y_valid)",
      "    if sample_weight is not None:",
      "        sample_weight = np.asarray(sample_weight, dtype=np.float32)",
      "    _omicselector_seed(int(config.get('seed', 1)))",
      "    device = _omicselector_resolve_device(config.get('device', 'auto'))",
      "    n_classes = int(len(np.unique(y_train)))",
      "    is_binary = n_classes == 2",
      "    means = X_train.mean(axis=0).astype(np.float32)",
      "    sds = X_train.std(axis=0).astype(np.float32)",
      "    sds[sds == 0] = 1.0",
      "    X_train_s = ((X_train - means) / sds).astype(np.float32)",
      "    X_valid_s = None if X_valid is None else ((X_valid - means) / sds).astype(np.float32)",
      "    model = _omicselector_build_model(model_name, X_train_s.shape[1], n_classes, config).to(device)",
      "    if hasattr(model, 'make_parameter_groups'):",
      "        params = model.make_parameter_groups()",
      "    else:",
      "        params = model.parameters()",
      "    optimizer = torch.optim.AdamW(",
      "        params,",
      "        lr=float(config.get('learning_rate', 1e-3)),",
      "        weight_decay=float(config.get('weight_decay', 1e-5))",
      "    )",
      "    X_tensor = torch.as_tensor(X_train_s, dtype=torch.float32)",
      "    y_dtype = torch.float32 if is_binary else torch.long",
      "    y_tensor = torch.as_tensor(y_train, dtype=y_dtype)",
      "    if sample_weight is None:",
      "        w_tensor = None",
      "    else:",
      "        w_tensor = torch.as_tensor(sample_weight, dtype=torch.float32)",
      "    dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor) if w_tensor is None else torch.utils.data.TensorDataset(X_tensor, y_tensor, w_tensor)",
      "    generator = torch.Generator()",
      "    generator.manual_seed(int(config.get('seed', 1)))",
      "    loader = torch.utils.data.DataLoader(",
      "        dataset,",
      "        batch_size=int(config.get('batch_size', 128)),",
      "        shuffle=True,",
      "        drop_last=False,",
      "        generator=generator",
      "    )",
      "    best_state = copy.deepcopy(model.state_dict())",
      "    best_metric = float('-inf')",
      "    best_epoch = 1",
      "    stale = 0",
      "    max_epochs = int(config.get('max_epochs', 200))",
      "    patience = int(config.get('patience', 20))",
      "    for epoch in range(1, max_epochs + 1):",
      "        model.train()",
      "        for batch in loader:",
      "            if w_tensor is None:",
      "                xb, yb = batch",
      "                wb = None",
      "            else:",
      "                xb, yb, wb = batch",
      "                wb = wb.to(device)",
      "            xb = xb.to(device)",
      "            yb = yb.to(device)",
      "            optimizer.zero_grad(set_to_none=True)",
      "            logits = model(xb, None)",
      "            loss = _omicselector_loss(logits, yb, wb, n_classes, model_name == 'tabm', torch)",
      "            loss.backward()",
      "            optimizer.step()",
      "        if X_valid_s is None:",
      "            best_epoch = epoch",
      "            best_state = copy.deepcopy(model.state_dict())",
      "            continue",
      "        model.eval()",
      "        with torch.no_grad():",
      "            logits = model(torch.as_tensor(X_valid_s, dtype=torch.float32, device=device), None)",
      "            probs = _omicselector_logits_to_probs(logits, n_classes, model_name == 'tabm', torch)",
      "            if is_binary:",
      "                metric = float(probs.detach().cpu().numpy().mean())",
      "                try:",
      "                    from sklearn.metrics import roc_auc_score",
      "                    metric = float(roc_auc_score(y_valid, probs.detach().cpu().numpy()))",
      "                except Exception:",
      "                    pass",
      "            else:",
      "                prob_np = probs.detach().cpu().numpy()",
      "                true_prob = prob_np[np.arange(len(y_valid)), np.asarray(y_valid, dtype=np.int64)]",
      "                metric = float(np.mean(true_prob))",
      "        if metric > best_metric + 1e-8:",
      "            best_metric = metric",
      "            best_epoch = epoch",
      "            best_state = copy.deepcopy(model.state_dict())",
      "            stale = 0",
      "        else:",
      "            stale += 1",
      "            if patience > 0 and stale >= patience:",
      "                break",
      "    model.load_state_dict(best_state)",
      "    return {",
      "        'model': model,",
      "        'means': means,",
      "        'sds': sds,",
      "        'device': device,",
      "        'n_classes': n_classes,",
      "        'is_binary': bool(is_binary),",
      "        'best_epoch': int(best_epoch),",
      "    }",
      "",
      "def omicselector_predict_rtdl(bundle, X):",
      "    import torch",
      "    X = np.asarray(X, dtype=np.float32)",
      "    means = np.asarray(bundle['means'], dtype=np.float32)",
      "    sds = np.asarray(bundle['sds'], dtype=np.float32)",
      "    sds[sds == 0] = 1.0",
      "    X_s = ((X - means) / sds).astype(np.float32)",
      "    model = bundle['model']",
      "    device = bundle['device']",
      "    n_classes = int(bundle['n_classes'])",
      "    is_binary = bool(bundle['is_binary'])",
      "    model_name = 'tabm' if model.__class__.__name__.lower().startswith('tabm') else 'fttransformer'",
      "    model.eval()",
      "    with torch.no_grad():",
      "        logits = model(torch.as_tensor(X_s, dtype=torch.float32, device=device), None)",
      "        probs = _omicselector_logits_to_probs(logits, n_classes, model_name == 'tabm', torch)",
      "    probs = probs.detach().cpu().numpy()",
      "    return probs.astype(np.float64)",
      ""
    ),
    collapse = "\n"
  )

  reticulate::py_run_string(code)
  invisible(TRUE)
}


#' @keywords internal
.tabdl_train_tabnet <- function(task, config, class_weights) {
  .tabdl_require_python_module(
    "pytorch_tabnet.tab_model",
    paste(
      "Python module 'pytorch_tabnet' is required for TabNet.",
      "Install in the active reticulate environment with: pip install pytorch-tabnet"
    )
  )

  data <- .tabdl_training_bundle(task)
  device <- .tabdl_resolve_device(config$device)
  split <- .tabdl_validation_split(
    y = data$y_backend,
    validation_fraction = config$validation_fraction,
    seed = config$seed
  )

  if (!is.null(data$task_weights)) {
    .tabdl_warn_ignored_weights("TabNet")
  }

  class_weight_backend <- .tabdl_backend_class_weights(data, class_weights)
  weights_arg <- if (identical(class_weights, "balanced")) {
    1L
  } else if (is.null(class_weight_backend)) {
    0L
  } else {
    .tabdl_load_python_helpers()
    reticulate::py$omicselector_make_class_weight_dict(as.numeric(class_weight_backend))
  }

  tabnet <- reticulate::import("pytorch_tabnet.tab_model")
  torch_py <- reticulate::import("torch")

  model <- tabnet$TabNetClassifier(
    n_d = as.integer(config$n_d),
    n_a = as.integer(config$n_a),
    n_steps = as.integer(config$n_steps),
    gamma = config$gamma,
    n_independent = as.integer(config$n_independent),
    n_shared = as.integer(config$n_shared),
    lambda_sparse = config$lambda_sparse,
    optimizer_fn = torch_py$optim$Adam,
    optimizer_params = list(
      lr = config$learning_rate,
      weight_decay = config$weight_decay
    ),
    mask_type = config$mask_type,
    device_name = if (identical(device, "cuda")) "gpu" else "cpu",
    seed = as.integer(config$seed),
    verbose = 0L
  )

  fit_args <- list(
    X_train = data$x[split$train, , drop = FALSE],
    y_train = data$y_backend[split$train],
    max_epochs = as.integer(config$max_epochs),
    patience = as.integer(config$patience),
    batch_size = as.integer(config$batch_size),
    virtual_batch_size = as.integer(config$virtual_batch_size),
    num_workers = 0L,
    drop_last = FALSE,
    weights = weights_arg,
    eval_metric = if (isTRUE(data$is_binary)) list("auc") else list("logloss")
  )

  if (length(split$valid) > 0L) {
    fit_args$eval_set <- list(list(
      data$x[split$valid, , drop = FALSE],
      data$y_backend[split$valid]
    ))
    fit_args$eval_name <- list("valid")
  }

  do.call(model$fit, fit_args)

  backend_classes <- tryCatch(reticulate::py_to_r(model$classes_), error = function(e) NULL)

  list(
    model = model,
    feature_names = data$feature_names,
    class_names = data$class_names,
    backend_classes = backend_classes,
    is_binary = data$is_binary
  )
}


#' @keywords internal
.tabdl_predict_tabnet <- function(model_info, task) {
  x <- .tabdl_feature_matrix(task, feature_names = model_info$feature_names)
  prob_raw <- model_info$model$predict_proba(x)
  prob <- .tabdl_reorder_python_prob(
    prob = prob_raw,
    class_names = model_info$class_names,
    backend_classes = model_info$backend_classes,
    is_binary = model_info$is_binary
  )

  mlr3::PredictionClassif$new(
    task = task,
    row_ids = task$row_ids,
    prob = prob
  )
}


#' @keywords internal
.tabdl_train_rtdl <- function(task, config, class_weights, model_name) {
  .tabdl_require_python_module(
    "rtdl_revisiting_models",
    paste(
      "Python module 'rtdl_revisiting_models' is required for",
      if (identical(model_name, "fttransformer")) "FT-Transformer." else "TabM.",
      "Install in the active reticulate environment with: pip install rtdl_revisiting_models",
      if (identical(model_name, "tabm")) "and, if needed, pip install tabm" else ""
    )
  )

  .tabdl_load_python_helpers()

  data <- .tabdl_training_bundle(task)
  split <- .tabdl_validation_split(
    y = data$y_backend,
    validation_fraction = config$validation_fraction,
    seed = config$seed
  )

  x_train <- data$x[split$train, , drop = FALSE]
  y_train <- data$y_backend[split$train]
  x_valid <- if (length(split$valid) > 0L) data$x[split$valid, , drop = FALSE] else NULL
  y_valid <- if (length(split$valid) > 0L) data$y_backend[split$valid] else NULL
  sample_weights <- .tabdl_sample_weights(
    outcome = data,
    task_weights = data$task_weights,
    class_weights = class_weights,
    indices = split$train
  )

  bundle <- reticulate::py$omicselector_fit_rtdl(
    model_name = model_name,
    X_train = x_train,
    y_train = y_train,
    X_valid = x_valid,
    y_valid = y_valid,
    sample_weight = sample_weights,
    config = list(
      seed = as.integer(config$seed),
      device = .tabdl_resolve_device(config$device),
      learning_rate = config$learning_rate,
      weight_decay = config$weight_decay,
      batch_size = as.integer(config$batch_size),
      max_epochs = as.integer(config$max_epochs),
      patience = as.integer(config$patience),
      n_blocks = as.integer(config$n_blocks),
      d_block = if (!is.null(config$d_block)) as.integer(config$d_block) else NULL,
      attention_dropout = config$attention_dropout,
      ffn_dropout = config$ffn_dropout,
      ffn_d_hidden = if (!is.null(config$ffn_d_hidden)) as.integer(config$ffn_d_hidden) else NULL,
      k = if (!is.null(config$k)) as.integer(config$k) else NULL,
      dropout = config$dropout
    )
  )

  list(
    bundle = bundle,
    feature_names = data$feature_names,
    class_names = data$class_names,
    is_binary = data$is_binary
  )
}


#' @keywords internal
.tabdl_predict_rtdl <- function(model_info, task) {
  .tabdl_load_python_helpers()
  x <- .tabdl_feature_matrix(task, feature_names = model_info$feature_names)
  prob_raw <- reticulate::py$omicselector_predict_rtdl(model_info$bundle, x)

  prob <- if (isTRUE(model_info$is_binary)) {
    .tabdl_binary_prob_matrix(prob_raw, model_info$class_names)
  } else {
    .tabdl_multiclass_prob_matrix(prob_raw, model_info$class_names)
  }

  mlr3::PredictionClassif$new(
    task = task,
    row_ids = task$row_ids,
    prob = prob
  )
}


#' @keywords internal
.tabdl_train_catboost <- function(task, config, class_weights) {
  .tabdl_require_namespace(
    "catboost",
    "Package 'catboost' is required for this learner.\nInstall it in R before training."
  )

  data <- .tabdl_training_bundle(task)
  split <- .tabdl_validation_split(
    y = data$y_backend,
    validation_fraction = config$validation_fraction,
    seed = config$seed
  )
  device <- .tabdl_resolve_device(config$device)

  x_train <- data$x[split$train, , drop = FALSE]
  y_train <- data$y_backend[split$train]
  train_weights <- .tabdl_sample_weights(
    outcome = data,
    task_weights = data$task_weights,
    class_weights = class_weights,
    indices = split$train
  )

  learn_pool <- catboost::catboost.load_pool(
    data = x_train,
    label = y_train,
    weight = train_weights
  )

  params <- list(
    iterations = as.integer(config$iterations),
    learning_rate = config$learning_rate,
    depth = as.integer(config$depth),
    l2_leaf_reg = config$l2_leaf_reg,
    border_count = as.integer(config$border_count),
    random_seed = as.integer(config$seed),
    allow_writing_files = FALSE,
    verbose = FALSE,
    task_type = if (identical(device, "cuda")) "GPU" else "CPU",
    loss_function = if (isTRUE(data$is_binary)) "Logloss" else "MultiClass",
    eval_metric = if (isTRUE(data$is_binary)) "AUC" else "MultiClass",
    od_type = "Iter",
    od_wait = as.integer(config$patience),
    use_best_model = TRUE
  )

  if (length(split$valid) > 0L) {
    x_valid <- data$x[split$valid, , drop = FALSE]
    y_valid <- data$y_backend[split$valid]
    valid_weights <- .tabdl_sample_weights(
      outcome = data,
      task_weights = data$task_weights,
      class_weights = class_weights,
      indices = split$valid
    )

    valid_pool <- catboost::catboost.load_pool(
      data = x_valid,
      label = y_valid,
      weight = valid_weights
    )
    model <- catboost::catboost.train(
      learn_pool = learn_pool,
      test_pool = valid_pool,
      params = params
    )
  } else {
    params$od_type <- NULL
    params$od_wait <- NULL
    params$use_best_model <- NULL
    params <- params[!vapply(params, is.null, logical(1))]
    model <- catboost::catboost.train(
      learn_pool = learn_pool,
      params = params
    )
  }

  list(
    model = model,
    feature_names = data$feature_names,
    class_names = data$class_names,
    is_binary = data$is_binary
  )
}


#' @keywords internal
.tabdl_predict_catboost <- function(model_info, task) {
  .tabdl_require_namespace(
    "catboost",
    "Package 'catboost' is required for this learner.\nInstall it in R before prediction."
  )

  x <- .tabdl_feature_matrix(task, feature_names = model_info$feature_names)
  pool <- catboost::catboost.load_pool(data = x)
  prob_raw <- catboost::catboost.predict(
    model_info$model,
    pool,
    prediction_type = "Probability"
  )

  prob <- if (isTRUE(model_info$is_binary)) {
    if (is.null(dim(prob_raw))) {
      .tabdl_binary_prob_matrix(prob_raw, model_info$class_names)
    } else {
      prob_raw <- as.matrix(prob_raw)
      .tabdl_binary_prob_matrix(prob_raw[, ncol(prob_raw)], model_info$class_names)
    }
  } else {
    .tabdl_multiclass_prob_matrix(prob_raw, model_info$class_names)
  }

  mlr3::PredictionClassif$new(
    task = task,
    row_ids = task$row_ids,
    prob = prob
  )
}


#' @keywords internal
.tabdl_train_tabpfn <- function(task, config, class_weights) {
  .tabdl_require_python_module(
    "tabpfn",
    paste(
      "Python module 'tabpfn' is required for TabPFN.",
      "Install in the active reticulate environment with: pip install tabpfn",
      "and accept the Prior Labs Hugging Face license if gated weights are used."
    )
  )

  data <- .tabdl_training_bundle(task)
  if (!is.null(data$task_weights) || !(is.null(class_weights) || identical(class_weights, "none"))) {
    warning(
      "TabPFN does not expose class/sample weighting through this wrapper; ignoring `class_weights` and task weights.",
      call. = FALSE
    )
  }

  tabpfn <- reticulate::import("tabpfn")
  device <- .tabdl_resolve_device(config$device)
  model <- NULL

  if (!is.null(config$model_version) &&
      reticulate::py_has_attr(tabpfn$TabPFNClassifier, "create_default_for_version")) {
    constants <- reticulate::import("tabpfn.constants", delay_load = TRUE)
    model_version <- constants$ModelVersion[[config$model_version]]
    model <- tabpfn$TabPFNClassifier$create_default_for_version(model_version)
  }

  if (is.null(model)) {
    model <- tryCatch(
      tabpfn$TabPFNClassifier(
        device = device,
        random_state = as.integer(config$seed)
      ),
      error = function(e) {
        tabpfn$TabPFNClassifier(device = device)
      }
    )
  }

  model$fit(data$x, data$y_backend)
  backend_classes <- tryCatch(reticulate::py_to_r(model$classes_), error = function(e) NULL)

  list(
    model = model,
    feature_names = data$feature_names,
    class_names = data$class_names,
    backend_classes = backend_classes,
    is_binary = data$is_binary
  )
}


#' @keywords internal
.tabdl_predict_tabpfn <- function(model_info, task) {
  x <- .tabdl_feature_matrix(task, feature_names = model_info$feature_names)
  prob_raw <- model_info$model$predict_proba(x)
  prob <- .tabdl_reorder_python_prob(
    prob = prob_raw,
    class_names = model_info$class_names,
    backend_classes = model_info$backend_classes,
    is_binary = model_info$is_binary
  )

  mlr3::PredictionClassif$new(
    task = task,
    row_ids = task$row_ids,
    prob = prob
  )
}


#' @keywords internal
LearnerClassifTabNet <- R6::R6Class(
  "LearnerClassifTabNet",
  inherit = mlr3::LearnerClassif,
  public = list(
    config = NULL,
    class_weights = NULL,
    initialize = function(config, class_weights) {
      self$config <- config
      self$class_weights <- class_weights
      super$initialize(
        id = "classif.tabnet",
        predict_types = c("response", "prob"),
        feature_types = c("logical", "integer", "numeric"),
        properties = c("twoclass", "multiclass"),
        packages = "reticulate",
        label = "TabNet Classifier"
      )
    }
  ),
  private = list(
    .train = function(task) {
      .tabdl_train_tabnet(task, config = self$config, class_weights = self$class_weights)
    },
    .predict = function(task) {
      .tabdl_predict_tabnet(self$model, task)
    }
  )
)


#' @keywords internal
LearnerClassifTabM <- R6::R6Class(
  "LearnerClassifTabM",
  inherit = mlr3::LearnerClassif,
  public = list(
    config = NULL,
    class_weights = NULL,
    initialize = function(config, class_weights) {
      self$config <- config
      self$class_weights <- class_weights
      super$initialize(
        id = "classif.tabm",
        predict_types = c("response", "prob"),
        feature_types = c("logical", "integer", "numeric"),
        properties = c("twoclass", "multiclass", "weights"),
        packages = "reticulate",
        label = "TabM Classifier"
      )
    }
  ),
  private = list(
    .train = function(task) {
      .tabdl_train_rtdl(task, config = self$config, class_weights = self$class_weights, model_name = "tabm")
    },
    .predict = function(task) {
      .tabdl_predict_rtdl(self$model, task)
    }
  )
)


#' @keywords internal
LearnerClassifFTTransformer <- R6::R6Class(
  "LearnerClassifFTTransformer",
  inherit = mlr3::LearnerClassif,
  public = list(
    config = NULL,
    class_weights = NULL,
    initialize = function(config, class_weights) {
      self$config <- config
      self$class_weights <- class_weights
      super$initialize(
        id = "classif.fttransformer",
        predict_types = c("response", "prob"),
        feature_types = c("logical", "integer", "numeric"),
        properties = c("twoclass", "multiclass", "weights"),
        packages = "reticulate",
        label = "FT-Transformer Classifier"
      )
    }
  ),
  private = list(
    .train = function(task) {
      .tabdl_train_rtdl(task, config = self$config, class_weights = self$class_weights, model_name = "fttransformer")
    },
    .predict = function(task) {
      .tabdl_predict_rtdl(self$model, task)
    }
  )
)


#' @keywords internal
LearnerClassifCatBoost <- R6::R6Class(
  "LearnerClassifCatBoost",
  inherit = mlr3::LearnerClassif,
  public = list(
    config = NULL,
    class_weights = NULL,
    initialize = function(config, class_weights) {
      self$config <- config
      self$class_weights <- class_weights
      super$initialize(
        id = "classif.catboost_omic",
        predict_types = c("response", "prob"),
        feature_types = c("logical", "integer", "numeric"),
        properties = c("twoclass", "multiclass", "weights"),
        packages = "catboost",
        label = "CatBoost Classifier"
      )
    }
  ),
  private = list(
    .train = function(task) {
      .tabdl_train_catboost(task, config = self$config, class_weights = self$class_weights)
    },
    .predict = function(task) {
      .tabdl_predict_catboost(self$model, task)
    }
  )
)


#' @keywords internal
LearnerClassifTabPFN <- R6::R6Class(
  "LearnerClassifTabPFN",
  inherit = mlr3::LearnerClassif,
  public = list(
    config = NULL,
    class_weights = NULL,
    initialize = function(config, class_weights) {
      self$config <- config
      self$class_weights <- class_weights
      super$initialize(
        id = "classif.tabpfn",
        predict_types = c("response", "prob"),
        feature_types = c("logical", "integer", "numeric"),
        properties = c("twoclass", "multiclass"),
        packages = "reticulate",
        label = "TabPFN Classifier"
      )
    }
  ),
  private = list(
    .train = function(task) {
      .tabdl_train_tabpfn(task, config = self$config, class_weights = self$class_weights)
    },
    .predict = function(task) {
      .tabdl_predict_tabpfn(self$model, task)
    }
  )
)


#' @title Create a TabNet Learner
#'
#' @description
#' Create an `mlr3` classification learner backed by
#' `pytorch_tabnet::TabNetClassifier` through `reticulate`.
#'
#' @param n_d Width of the decision layer. Default: `16`.
#' @param n_a Width of the attention layer. Default: `16`.
#' @param n_steps Number of sequential attention steps. Default: `4`.
#' @param gamma Feature reusage penalty. Default: `1.5`.
#' @param n_independent Number of independent GLU blocks. Default: `2`.
#' @param n_shared Number of shared GLU blocks. Default: `2`.
#' @param lambda_sparse Sparsity regularization strength. Default: `1e-4`.
#' @param batch_size Minibatch size. Default: `128`.
#' @param virtual_batch_size Ghost batch normalization size. Default: `32`.
#' @param max_epochs Maximum epochs. Default: `200`.
#' @param patience Early stopping patience. Default: `20`.
#' @param learning_rate Optimizer learning rate. Default: `0.02`.
#' @param weight_decay L2 regularization for Adam. Default: `1e-5`.
#' @param mask_type Mask function, `"sparsemax"` or `"entmax"`. Default: `"entmax"`.
#' @param validation_fraction Fraction of the training fold used for internal
#'   early stopping. Default: `0.15`.
#' @param class_weights Class weighting scheme. Use `"balanced"` (default),
#'   `"none"`, or a numeric vector with one weight per class. Named vectors are
#'   matched against the task class labels.
#' @param device `"auto"` (default), `"cpu"`, or `"cuda"`.
#' @param seed Random seed used inside the backend. Default: `20260331`.
#'
#' @return An `mlr3::LearnerClassif` with `predict_type = "prob"`.
#'
#' @references
#' Arik, S. O., & Pfister, T. (2021). TabNet: Attentive Interpretable Tabular
#' Learning.
#'
#' @export
make_tabnet_learner <- function(n_d = 16L,
                                n_a = 16L,
                                n_steps = 4L,
                                gamma = 1.5,
                                n_independent = 2L,
                                n_shared = 2L,
                                lambda_sparse = 1e-4,
                                batch_size = 128L,
                                virtual_batch_size = 32L,
                                max_epochs = 200L,
                                patience = 20L,
                                learning_rate = 0.02,
                                weight_decay = 1e-5,
                                mask_type = c("entmax", "sparsemax"),
                                validation_fraction = 0.15,
                                class_weights = "balanced",
                                device = c("auto", "cpu", "cuda"),
                                seed = 20260331L) {
  mask_type <- match.arg(mask_type)
  device <- match.arg(device)

  learner <- LearnerClassifTabNet$new(
    config = list(
      n_d = as.integer(n_d),
      n_a = as.integer(n_a),
      n_steps = as.integer(n_steps),
      gamma = gamma,
      n_independent = as.integer(n_independent),
      n_shared = as.integer(n_shared),
      lambda_sparse = lambda_sparse,
      batch_size = as.integer(batch_size),
      virtual_batch_size = as.integer(virtual_batch_size),
      max_epochs = as.integer(max_epochs),
      patience = as.integer(patience),
      learning_rate = learning_rate,
      weight_decay = weight_decay,
      mask_type = mask_type,
      validation_fraction = validation_fraction,
      device = device,
      seed = as.integer(seed)
    ),
    class_weights = class_weights
  )
  learner$predict_type <- "prob"
  learner
}


#' @title Create a TabM Learner
#'
#' @description
#' Create an `mlr3` classification learner for the modern MLP baseline TabM
#' through `reticulate`. The wrapper accepts either the official `tabm` package
#' or a compatible `rtdl_revisiting_models` build that exposes `TabM`.
#'
#' @param d_block Hidden block width. Default: `128`.
#' @param n_blocks Number of MLP blocks. Default: `3`.
#' @param dropout Model dropout. Default: `0.10`.
#' @param k Optional TabM ensemble width if supported by the backend. Default:
#'   `NULL`.
#' @param batch_size Minibatch size. Default: `128`.
#' @param max_epochs Maximum epochs. Default: `200`.
#' @param patience Early stopping patience. Default: `20`.
#' @param learning_rate Optimizer learning rate. Default: `1e-3`.
#' @param weight_decay AdamW weight decay. Default: `1e-5`.
#' @param validation_fraction Fraction of the training fold used for internal
#'   validation. Default: `0.15`.
#' @param class_weights Class weighting scheme. Use `"balanced"` (default),
#'   `"none"`, or a numeric vector with one weight per class.
#' @param device `"auto"` (default), `"cpu"`, or `"cuda"`.
#' @param seed Random seed used inside the backend. Default: `20260331`.
#'
#' @return An `mlr3::LearnerClassif` with `predict_type = "prob"`.
#'
#' @references
#' Gorishniy, Y., et al. (2024). TabM: Advancing Tabular Deep Learning with
#' Parameter-Efficient Ensembling.
#'
#' @export
make_tabm_learner <- function(d_block = 128L,
                              n_blocks = 3L,
                              dropout = 0.10,
                              k = NULL,
                              batch_size = 128L,
                              max_epochs = 200L,
                              patience = 20L,
                              learning_rate = 1e-3,
                              weight_decay = 1e-5,
                              validation_fraction = 0.15,
                              class_weights = "balanced",
                              device = c("auto", "cpu", "cuda"),
                              seed = 20260331L) {
  device <- match.arg(device)

  learner <- LearnerClassifTabM$new(
    config = list(
      d_block = as.integer(d_block),
      n_blocks = as.integer(n_blocks),
      dropout = dropout,
      k = if (is.null(k)) NULL else as.integer(k),
      batch_size = as.integer(batch_size),
      max_epochs = as.integer(max_epochs),
      patience = as.integer(patience),
      learning_rate = learning_rate,
      weight_decay = weight_decay,
      validation_fraction = validation_fraction,
      device = device,
      seed = as.integer(seed),
      attention_dropout = NULL,
      ffn_dropout = NULL,
      ffn_d_hidden = NULL
    ),
    class_weights = class_weights
  )
  learner$predict_type <- "prob"
  learner
}


#' @title Create an FT-Transformer Learner
#'
#' @description
#' Create an `mlr3` classification learner backed by
#' `rtdl_revisiting_models::FTTransformer` through `reticulate`.
#'
#' @param n_blocks Number of transformer blocks. Default: `3`.
#' @param d_block Transformer width. Default: `192`.
#' @param attention_dropout Attention dropout. Default: `0.20`.
#' @param ffn_dropout Feed-forward dropout. Default: `0.20`.
#' @param ffn_d_hidden Hidden width in the feed-forward network. Default:
#'   `NULL`, which leaves the backend default unchanged.
#' @param batch_size Minibatch size. Default: `128`.
#' @param max_epochs Maximum epochs. Default: `200`.
#' @param patience Early stopping patience. Default: `20`.
#' @param learning_rate Optimizer learning rate. Default: `1e-3`.
#' @param weight_decay AdamW weight decay. Default: `1e-5`.
#' @param validation_fraction Fraction of the training fold used for internal
#'   validation. Default: `0.15`.
#' @param class_weights Class weighting scheme. Use `"balanced"` (default),
#'   `"none"`, or a numeric vector with one weight per class.
#' @param device `"auto"` (default), `"cpu"`, or `"cuda"`.
#' @param seed Random seed used inside the backend. Default: `20260331`.
#'
#' @return An `mlr3::LearnerClassif` with `predict_type = "prob"`.
#'
#' @references
#' Gorishniy, Y., et al. (2021). Revisiting Deep Learning Models for Tabular Data.
#'
#' @export
make_fttransformer_learner <- function(n_blocks = 3L,
                                       d_block = 192L,
                                       attention_dropout = 0.20,
                                       ffn_dropout = 0.20,
                                       ffn_d_hidden = NULL,
                                       batch_size = 128L,
                                       max_epochs = 200L,
                                       patience = 20L,
                                       learning_rate = 1e-3,
                                       weight_decay = 1e-5,
                                       validation_fraction = 0.15,
                                       class_weights = "balanced",
                                       device = c("auto", "cpu", "cuda"),
                                       seed = 20260331L) {
  device <- match.arg(device)

  learner <- LearnerClassifFTTransformer$new(
    config = list(
      n_blocks = as.integer(n_blocks),
      d_block = as.integer(d_block),
      attention_dropout = attention_dropout,
      ffn_dropout = ffn_dropout,
      ffn_d_hidden = if (is.null(ffn_d_hidden)) NULL else as.integer(ffn_d_hidden),
      batch_size = as.integer(batch_size),
      max_epochs = as.integer(max_epochs),
      patience = as.integer(patience),
      learning_rate = learning_rate,
      weight_decay = weight_decay,
      validation_fraction = validation_fraction,
      device = device,
      seed = as.integer(seed),
      k = NULL,
      dropout = NULL
    ),
    class_weights = class_weights
  )
  learner$predict_type <- "prob"
  learner
}


#' @rdname make_fttransformer_learner
#' @export
make_tabtransformer_learner <- function(...) {
  make_fttransformer_learner(...)
}


#' @title Create a CatBoost Learner
#'
#' @description
#' Create an `mlr3` classification learner backed by the `catboost` R package.
#'
#' @param iterations Maximum boosting iterations. Default: `500`.
#' @param learning_rate Learning rate. Default: `0.03`.
#' @param depth Tree depth. Default: `6`.
#' @param l2_leaf_reg L2 regularization on leaf values. Default: `5`.
#' @param border_count Number of numerical split candidates. Default: `64`.
#' @param patience Overfitting detector patience for the internal validation
#'   split. Default: `40`.
#' @param validation_fraction Fraction of the training fold reserved for
#'   early-stopping. Default: `0.15`.
#' @param class_weights Class weighting scheme. Use `"balanced"` (default),
#'   `"none"`, or a numeric vector with one weight per class.
#' @param device `"auto"` (default), `"cpu"`, or `"cuda"`.
#' @param seed Random seed used inside the backend. Default: `20260331`.
#'
#' @return An `mlr3::LearnerClassif` with `predict_type = "prob"`.
#'
#' @references
#' Prokhorenkova, L., et al. (2018). CatBoost: unbiased boosting with
#' categorical features.
#'
#' @export
make_catboost_learner <- function(iterations = 500L,
                                  learning_rate = 0.03,
                                  depth = 6L,
                                  l2_leaf_reg = 5,
                                  border_count = 64L,
                                  patience = 40L,
                                  validation_fraction = 0.15,
                                  class_weights = "balanced",
                                  device = c("auto", "cpu", "cuda"),
                                  seed = 20260331L) {
  device <- match.arg(device)

  learner <- LearnerClassifCatBoost$new(
    config = list(
      iterations = as.integer(iterations),
      learning_rate = learning_rate,
      depth = as.integer(depth),
      l2_leaf_reg = l2_leaf_reg,
      border_count = as.integer(border_count),
      patience = as.integer(patience),
      validation_fraction = validation_fraction,
      device = device,
      seed = as.integer(seed)
    ),
    class_weights = class_weights
  )
  learner$predict_type <- "prob"
  learner
}


#' @title Create a TabPFN Learner
#'
#' @description
#' Create an `mlr3` classification learner backed by `tabpfn` through
#' `reticulate`. This is a zero-shot tabular foundation model suited to small
#' and medium tabular problems. TabPFN may require gated model access through
#' Hugging Face depending on the selected model version.
#'
#' @param model_version Optional TabPFN model version identifier exposed by
#'   `tabpfn.constants.ModelVersion`, for example `"V2_5"`. Default: `NULL`,
#'   which uses the package default.
#' @param class_weights Class weighting scheme. Included for interface
#'   consistency, but ignored by the current wrapper because TabPFN does not
#'   expose training-time class weighting here.
#' @param device `"auto"` (default), `"cpu"`, or `"cuda"`.
#' @param seed Random seed used inside the backend. Default: `20260331`.
#'
#' @return An `mlr3::LearnerClassif` with `predict_type = "prob"`.
#'
#' @references
#' Hollmann, N., et al. (2025). TabPFN-2.
#'
#' @export
make_tabpfn_learner <- function(model_version = NULL,
                                class_weights = "balanced",
                                device = c("auto", "cpu", "cuda"),
                                seed = 20260331L) {
  device <- match.arg(device)

  learner <- LearnerClassifTabPFN$new(
    config = list(
      model_version = model_version,
      device = device,
      seed = as.integer(seed)
    ),
    class_weights = class_weights
  )
  learner$predict_type <- "prob"
  learner
}
