# DANN encoder with separate cohort and kit adversarial heads.

if (!exists(".lka_as_matrix", mode = "function")) {
  stop("Learned kit-aware helper functions are unavailable.", call. = FALSE)
}

.lka_dann_module <- NULL

.lka_build_dann_module <- function() {
  if (!is.null(.lka_dann_module)) return(.lka_dann_module)
  .lka_dann_module <<- torch::nn_module(
    "LkaDannKitExtended",
    initialize = function(n_features, hidden_dim, latent_dim, n_cohorts, n_kits,
                          dropout = 0.10) {
      self$encoder <- torch::nn_sequential(
        torch::nn_linear(n_features, hidden_dim),
        torch::nn_relu(),
        torch::nn_dropout(dropout),
        torch::nn_linear(hidden_dim, latent_dim),
        torch::nn_relu()
      )
      self$decoder <- torch::nn_sequential(
        torch::nn_linear(latent_dim, hidden_dim),
        torch::nn_relu(),
        torch::nn_linear(hidden_dim, n_features)
      )
      self$disease <- torch::nn_linear(latent_dim, 1)
      self$cohort <- torch::nn_sequential(
        torch::nn_linear(latent_dim, max(4L, latent_dim)),
        torch::nn_relu(),
        torch::nn_linear(max(4L, latent_dim), n_cohorts)
      )
      self$kit <- torch::nn_sequential(
        torch::nn_linear(latent_dim, max(4L, latent_dim)),
        torch::nn_relu(),
        torch::nn_linear(max(4L, latent_dim), n_kits)
      )
    },
    encode = function(x) {
      self$encoder(x)
    },
    project = function(x) {
      self$decoder(self$encoder(x))
    },
    forward = function(x) {
      z <- self$encoder(x)
      list(recon = self$decoder(z),
           disease_logit = self$disease(z)$squeeze(2),
           cohort_logit = self$cohort(z),
           kit_logit = self$kit(z))
    }
  )
  .lka_dann_module
}

.lka_dann_project <- function(model, expr) {
  if (!identical(model$framework, "torch")) {
    z <- .lka_standardize_apply(expr, model$scaler)
    return(z)
  }
  x <- .lka_standardize_apply(expr, model$scaler)
  model$torch_model$eval()
  torch::with_no_grad({
    out <- model$torch_model$project(torch::torch_tensor(x, dtype = torch::torch_float()))
  })
  z <- torch::as_array(out)
  colnames(z) <- model$features
  rownames(z) <- rownames(expr)
  z
}

fit_dann_kit_extended <- function(train_expr, train_meta = NULL,
                                  train_disease_labels,
                                  train_kit_labels,
                                  train_cohort_labels,
                                  panel_features = NULL,
                                  epochs = 8L,
                                  batch_size = 128L,
                                  hidden_dim = NULL,
                                  latent_dim = 12L,
                                  lambda_cohort = 0.15,
                                  lambda_kit = 0.15,
                                  disease_weight = 1.0,
                                  learning_rate = 1e-3,
                                  seed = 42L,
                                  use_torch = TRUE) {
  x <- .lka_as_matrix(train_expr)
  y <- .lka_as_binary(train_disease_labels)
  if (length(y) != nrow(x)) stop("train_disease_labels length does not match train_expr rows.")
  if (is.null(panel_features)) panel_features <- colnames(x)
  panel_features <- intersect(panel_features, colnames(x))
  if (length(panel_features) == 0L) stop("panel_features do not overlap train_expr.")
  kit <- .lka_factor_state(train_kit_labels)
  cohort <- .lka_factor_state(train_cohort_labels)
  if (length(kit$index) != nrow(x) || length(cohort$index) != nrow(x)) {
    stop("kit/cohort label lengths must match train_expr rows.")
  }

  if (!use_torch || !.lka_torch_available()) {
    return(.lka_fit_linear_fallback(x, y, train_kit_labels, train_cohort_labels,
                                    panel_features, method = "dann_kit_extended"))
  }

  fit <- tryCatch({
    torch::torch_manual_seed(as.integer(seed))
    scaler <- .lka_standardize_fit(x)
    z <- .lka_standardize_apply(x, scaler)
    n_features <- ncol(z)
    if (is.null(hidden_dim)) hidden_dim <- min(96L, max(24L, ceiling(sqrt(n_features) * 3L)))
    latent_dim <- min(as.integer(latent_dim), max(4L, hidden_dim - 1L))
    mod <- .lka_build_dann_module()(n_features = n_features,
                                    hidden_dim = as.integer(hidden_dim),
                                    latent_dim = as.integer(latent_dim),
                                    n_cohorts = length(cohort$levels),
                                    n_kits = length(kit$levels))
    opt_heads <- torch::optim_adam(c(mod$cohort$parameters, mod$kit$parameters),
                                   lr = learning_rate)
    opt_main <- torch::optim_adam(c(mod$encoder$parameters, mod$decoder$parameters,
                                    mod$disease$parameters),
                                  lr = learning_rate, weight_decay = 1e-5)
    x_t <- torch::torch_tensor(z, dtype = torch::torch_float())
    y_t <- torch::torch_tensor(as.numeric(y), dtype = torch::torch_float())
    cohort_t <- torch::torch_tensor(as.integer(cohort$index), dtype = torch::torch_long())
    kit_t <- torch::torch_tensor(as.integer(kit$index), dtype = torch::torch_long())
    n <- nrow(z)
    batch_size <- as.integer(max(16L, min(batch_size, n)))
    for (ep in seq_len(max(1L, as.integer(epochs)))) {
      ord <- sample.int(n)
      grl_scale <- 2 / (1 + exp(-10 * (ep - 1) / max(1, epochs - 1))) - 1
      for (start in seq(1L, n, by = batch_size)) {
        idx <- ord[start:min(n, start + batch_size - 1L)]
        xb <- x_t[idx, ]
        yb <- y_t[idx]
        cb <- cohort_t[idx]
        kb <- kit_t[idx]

        z_det <- mod$encode(xb)$detach()
        loss_heads <- torch::nnf_cross_entropy(mod$cohort(z_det), cb) +
          torch::nnf_cross_entropy(mod$kit(z_det), kb)
        opt_heads$zero_grad()
        loss_heads$backward()
        opt_heads$step()

        out <- mod(xb)
        recon_loss <- torch::nnf_mse_loss(out$recon, xb)
        disease_loss <- torch::nnf_binary_cross_entropy_with_logits(out$disease_logit, yb)
        cohort_loss <- if (length(cohort$levels) > 1L)
          torch::nnf_cross_entropy(out$cohort_logit, cb) else torch::torch_tensor(0)
        kit_loss <- if (length(kit$levels) > 1L)
          torch::nnf_cross_entropy(out$kit_logit, kb) else torch::torch_tensor(0)
        loss <- recon_loss + disease_weight * disease_loss -
          (lambda_cohort * grl_scale) * cohort_loss -
          (lambda_kit * grl_scale) * kit_loss
        opt_main$zero_grad()
        loss$backward()
        opt_main$step()
      }
    }
    mod$eval()
    torch::with_no_grad({
      projected <- mod$project(x_t)
    })
    projected <- torch::as_array(projected)
    colnames(projected) <- colnames(x)
    weights <- .lka_feature_weights(projected, y, colnames(x))
    list(method = "dann_kit_extended",
         framework = "torch",
         framework_detail = "R torch alternating adversarial update; encoder step maximizes cohort and kit cross-entropy with separate lambda weights, equivalent to gradient reversal on encoder parameters",
         torch_model = mod,
         scaler = scaler,
         features = colnames(x),
         feature_weights = weights,
         panel_features = panel_features,
         kit_levels = kit$levels,
         cohort_levels = cohort$levels,
         epochs = as.integer(epochs),
         lambda_cohort = lambda_cohort,
         lambda_kit = lambda_kit,
         latent_dim = as.integer(latent_dim),
         hidden_dim = as.integer(hidden_dim))
  }, error = function(e) {
    warning("fit_dann_kit_extended torch path failed; using linear fallback: ",
            conditionMessage(e), call. = FALSE)
    .lka_fit_linear_fallback(x, y, train_kit_labels, train_cohort_labels,
                             panel_features, method = "dann_kit_extended")
  })
  fit
}

score_dann_kit_extended <- function(model, test_expr, panel_features = NULL) {
  if (!identical(model$framework, "torch")) {
    return(.lka_score_linear_fallback(model, test_expr, panel_features))
  }
  projected <- .lka_dann_project(model, test_expr)
  .lka_score_projection(projected, model$feature_weights,
                        panel_features %||% model$panel_features)
}

transform_dann_kit_extended <- function(model, expr) {
  .lka_dann_project(model, expr)
}
