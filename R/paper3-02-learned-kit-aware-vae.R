# Kit-conditional variational autoencoder for within-sample score projection.

if (!exists(".lka_as_matrix", mode = "function")) {
  stop("Learned kit-aware helper functions are unavailable.", call. = FALSE)
}

.lka_vae_module <- NULL

.lka_build_vae_module <- function() {
  if (!is.null(.lka_vae_module)) return(.lka_vae_module)
  .lka_vae_module <<- torch::nn_module(
    "LkaKitConditionalVae",
    initialize = function(n_features, n_kits, hidden_dim, latent_dim,
                          dropout = 0.10) {
      self$n_kits <- n_kits
      self$encoder <- torch::nn_sequential(
        torch::nn_linear(n_features + n_kits, hidden_dim),
        torch::nn_relu(),
        torch::nn_dropout(dropout)
      )
      self$mu <- torch::nn_linear(hidden_dim, latent_dim)
      self$logvar <- torch::nn_linear(hidden_dim, latent_dim)
      self$decoder <- torch::nn_sequential(
        torch::nn_linear(latent_dim + n_kits, hidden_dim),
        torch::nn_relu(),
        torch::nn_linear(hidden_dim, n_features)
      )
      self$disease <- torch::nn_linear(latent_dim, 1)
    },
    encode = function(x, kit_onehot) {
      h <- self$encoder(torch::torch_cat(list(x, kit_onehot), dim = 2))
      list(mu = self$mu(h), logvar = self$logvar(h))
    },
    reparameterize = function(mu, logvar) {
      eps <- torch::torch_randn_like(mu)
      mu + eps * torch::torch_exp(0.5 * logvar)
    },
    decode = function(z, kit_onehot) {
      self$decoder(torch::torch_cat(list(z, kit_onehot), dim = 2))
    },
    forward = function(x, kit_onehot) {
      enc <- self$encode(x, kit_onehot)
      z <- self$reparameterize(enc$mu, enc$logvar)
      list(recon = self$decode(z, kit_onehot),
           mu = enc$mu,
           logvar = enc$logvar,
           disease_logit = self$disease(enc$mu)$squeeze(2))
    }
  )
  .lka_vae_module
}

.lka_vae_project <- function(model, expr, test_kit_labels = NULL) {
  if (!identical(model$framework, "torch")) {
    return(.lka_standardize_apply(expr, model$scaler))
  }
  x <- .lka_standardize_apply(expr, model$scaler)
  x_t <- torch::torch_tensor(x, dtype = torch::torch_float())
  model$torch_model$eval()

  kit_indices <- NULL
  kit_weights <- NULL
  if (!is.null(test_kit_labels)) {
    kit_indices <- match(as.character(test_kit_labels), model$kit_levels)
    kit_indices[is.na(kit_indices)] <- model$default_kit_index
  } else {
    kit_indices <- seq_along(model$kit_levels)
    kit_weights <- model$kit_prior
  }

  torch::with_no_grad({
    if (!is.null(test_kit_labels)) {
      oh <- .lka_one_hot(kit_indices, length(model$kit_levels))
      oh_t <- torch::torch_tensor(oh, dtype = torch::torch_float())
      enc <- model$torch_model$encode(x_t, oh_t)
      out <- model$torch_model$decode(enc$mu, oh_t)
      projected <- torch::as_array(out)
    } else {
      accum <- matrix(0, nrow = nrow(x), ncol = ncol(x))
      for (kk in kit_indices) {
        oh <- .lka_one_hot(rep(kk, nrow(x)), length(model$kit_levels))
        oh_t <- torch::torch_tensor(oh, dtype = torch::torch_float())
        enc <- model$torch_model$encode(x_t, oh_t)
        out <- model$torch_model$decode(enc$mu, oh_t)
        accum <- accum + as.numeric(kit_weights[kk]) * torch::as_array(out)
      }
      projected <- accum
    }
  })
  colnames(projected) <- model$features
  rownames(projected) <- rownames(expr)
  projected
}

fit_kit_conditional_vae <- function(train_expr, train_meta = NULL,
                                    train_disease_labels,
                                    train_kit_labels,
                                    train_cohort_labels = NULL,
                                    panel_features = NULL,
                                    epochs = 8L,
                                    batch_size = 128L,
                                    hidden_dim = NULL,
                                    latent_dim = 12L,
                                    beta_kl = 0.03,
                                    disease_weight = 0.7,
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
  if (length(kit$index) != nrow(x)) stop("train_kit_labels length must match train_expr rows.")

  if (!use_torch || !.lka_torch_available()) {
    return(.lka_fit_linear_fallback(x, y, train_kit_labels, train_cohort_labels %||% "unknown",
                                    panel_features, method = "kit_conditional_vae"))
  }

  fit <- tryCatch({
    torch::torch_manual_seed(as.integer(seed))
    scaler <- .lka_standardize_fit(x)
    z <- .lka_standardize_apply(x, scaler)
    n_features <- ncol(z)
    if (is.null(hidden_dim)) hidden_dim <- min(96L, max(24L, ceiling(sqrt(n_features) * 3L)))
    latent_dim <- min(as.integer(latent_dim), max(4L, hidden_dim - 1L))
    mod <- .lka_build_vae_module()(n_features = n_features,
                                   n_kits = length(kit$levels),
                                   hidden_dim = as.integer(hidden_dim),
                                   latent_dim = as.integer(latent_dim))
    opt <- torch::optim_adam(mod$parameters, lr = learning_rate, weight_decay = 1e-5)
    x_t <- torch::torch_tensor(z, dtype = torch::torch_float())
    y_t <- torch::torch_tensor(as.numeric(y), dtype = torch::torch_float())
    kit_oh <- .lka_one_hot(kit$index, length(kit$levels))
    kit_t <- torch::torch_tensor(kit_oh, dtype = torch::torch_float())
    n <- nrow(z)
    batch_size <- as.integer(max(16L, min(batch_size, n)))
    for (ep in seq_len(max(1L, as.integer(epochs)))) {
      ord <- sample.int(n)
      for (start in seq(1L, n, by = batch_size)) {
        idx <- ord[start:min(n, start + batch_size - 1L)]
        out <- mod(x_t[idx, ], kit_t[idx, ])
        recon_loss <- torch::nnf_mse_loss(out$recon, x_t[idx, ])
        kld <- -0.5 * torch::torch_mean(1 + out$logvar -
                                          out$mu$pow(2) - out$logvar$exp())
        disease_loss <- torch::nnf_binary_cross_entropy_with_logits(out$disease_logit, y_t[idx])
        loss <- recon_loss + beta_kl * kld + disease_weight * disease_loss
        opt$zero_grad()
        loss$backward()
        opt$step()
      }
    }
    mod$eval()
    kit_prior <- as.numeric(table(factor(kit$index, levels = seq_along(kit$levels))))
    kit_prior <- kit_prior / sum(kit_prior)
    default_kit_index <- which.max(kit_prior)
    tmp_model <- list(framework = "torch", torch_model = mod, scaler = scaler,
                      kit_levels = kit$levels, kit_prior = kit_prior,
                      default_kit_index = default_kit_index,
                      features = colnames(x))
    projected <- .lka_vae_project(tmp_model, x)
    weights <- .lka_feature_weights(projected, y, colnames(x))
    c(tmp_model, list(method = "kit_conditional_vae",
                      framework_detail = "R torch conditional VAE; kit one-hot is concatenated to encoder and decoder, and test-time projection marginalizes over the training kit prior when test kit is not supplied",
                      feature_weights = weights,
                      panel_features = panel_features,
                      epochs = as.integer(epochs),
                      beta_kl = beta_kl,
                      latent_dim = as.integer(latent_dim),
                      hidden_dim = as.integer(hidden_dim)))
  }, error = function(e) {
    warning("fit_kit_conditional_vae torch path failed; using linear fallback: ",
            conditionMessage(e), call. = FALSE)
    .lka_fit_linear_fallback(x, y, train_kit_labels, train_cohort_labels %||% "unknown",
                             panel_features, method = "kit_conditional_vae")
  })
  fit
}

score_kit_conditional_vae <- function(model, test_expr, panel_features = NULL,
                                      test_kit_labels = NULL) {
  if (!identical(model$framework, "torch")) {
    return(.lka_score_linear_fallback(model, test_expr, panel_features))
  }
  projected <- .lka_vae_project(model, test_expr, test_kit_labels = test_kit_labels)
  .lka_score_projection(projected, model$feature_weights,
                        panel_features %||% model$panel_features)
}

transform_kit_conditional_vae <- function(model, expr, test_kit_labels = NULL) {
  .lka_vae_project(model, expr, test_kit_labels = test_kit_labels)
}
