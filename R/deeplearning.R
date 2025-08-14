#' Deep Learning Functions for OmicSelector
#'
#' This module provides comprehensive deep learning capabilities for miRNA biomarker
#' analysis using Keras/TensorFlow neural networks with hyperparameter optimization,
#' autoencoder support, and transfer learning.
#'
#' @name deeplearning
#' @docType package
NULL

#' Crea  codename <- sub(pattern = "(.*)[.].*$", replacement = "\\1", basename(output_file))e Keras Neural Network Model
#'
#' Creates a Keras neural network model based on hyperparameters for miRNA classification.
#' Supports multiple layers, dropout, regularization, and various activation functions.
#'
#' @param i Row index in hyperparameters dataframe
#' @param hyperparameters Dataframe with model configuration parameters
#' @param how_many_features Number of input features (miRNAs)
#'
#' @return Compiled Keras model ready for training
#'
#' @details
#' The hyperparameters dataframe should contain columns:
#' - layer1, layer2, layer3: Number of neurons in each layer
#' - activation_function_layer1/2/3: Activation functions
#' - dropout_layer1/2/3: Dropout rates
#' - layer1/2/3_regularizer: Whether to use L1 regularization
#' - optimizer: Optimization algorithm
#'
#' @examples
#' \dontrun{
#' hyperparams <- data.frame(
#'   layer1 = 64, layer2 = 32, layer3 = 0,
#'   activation_function_layer1 = "relu",
#'   activation_function_layer2 = "relu",
#'   activation_function_layer3 = "relu",
#'   dropout_layer1 = 0.2, dropout_layer2 = 0.1, dropout_layer3 = 0,
#'   layer1_regularizer = TRUE, layer2_regularizer = FALSE, layer3_regularizer = FALSE,
#'   optimizer = "adam"
#' )
#' model <- create_keras_model(1, hyperparams, 100)
#' }
#'
#' @export
create_keras_model <- function(i, hyperparameters, how_many_features) {
  # Check for keras availability
  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("Package 'keras' is required for deep learning functionality. Please install it with: install.packages('keras')")
  }
  
  # Validate inputs
  if (missing(i) || missing(hyperparameters) || missing(how_many_features)) {
    stop("All parameters (i, hyperparameters, how_many_features) are required")
  }
  
  if (i > nrow(hyperparameters) || i < 1) {
    stop("Index i must be within range of hyperparameters rows")
  }
  
  if (how_many_features <= 0) {
    stop("Number of features must be positive")
  }
  
  tryCatch({
    hyperparameters <- as.data.frame(hyperparameters)
    
    # Create sequential model
    model <- keras::keras_model_sequential()
    
    # First layer with optional regularization
    if (as.logical(hyperparameters[i, 10])) {
      keras::layer_dense(
        model,
        units = as.numeric(hyperparameters[i, 1]),
        kernel_regularizer = keras::regularizer_l1(l = 0.001),
        activation = as.character(hyperparameters[i, 4]),
        input_shape = c(how_many_features)
      )
    } else {
      keras::layer_dense(
        model,
        units = as.numeric(hyperparameters[i, 1]),
        activation = as.character(hyperparameters[i, 4]),
        input_shape = c(how_many_features)
      )
    }
    
    # First dropout layer
    if (as.numeric(hyperparameters[i, 7]) > 0) {
      keras::layer_dropout(model, rate = as.numeric(hyperparameters[i, 7]))
    }
    
    # Second layer (optional)
    if (as.numeric(hyperparameters[i, 2]) > 0) {
      if (as.logical(hyperparameters[i, 11])) {
        keras::layer_dense(
          model,
          units = as.numeric(hyperparameters[i, 2]),
          activation = as.character(hyperparameters[i, 5]),
          kernel_regularizer = keras::regularizer_l1(l = 0.001)
        )
      } else {
        keras::layer_dense(
          model,
          units = as.numeric(hyperparameters[i, 2]),
          activation = as.character(hyperparameters[i, 5])
        )
      }
    }
    
    # Second dropout layer
    if (as.numeric(hyperparameters[i, 2]) > 0 && as.numeric(hyperparameters[i, 8]) > 0) {
      keras::layer_dropout(model, rate = as.numeric(hyperparameters[i, 8]))
    }
    
    # Third layer (optional)
    if (as.numeric(hyperparameters[i, 3]) > 0) {
      if (as.logical(hyperparameters[i, 12])) {
        keras::layer_dense(
          model,
          units = as.numeric(hyperparameters[i, 3]),
          activation = as.character(hyperparameters[i, 6]),
          kernel_regularizer = keras::regularizer_l1(l = 0.001)
        )
      } else {
        keras::layer_dense(
          model,
          units = as.numeric(hyperparameters[i, 3]),
          activation = as.character(hyperparameters[i, 6])
        )
      }
    }
    
    # Third dropout layer
    if (as.numeric(hyperparameters[i, 3]) > 0 && as.numeric(hyperparameters[i, 9]) > 0) {
      keras::layer_dropout(model, rate = as.numeric(hyperparameters[i, 9]))
    }
    
    # Output layer for binary classification
    keras::layer_dense(model, units = 2, activation = 'softmax')
    
    # Compile model
    compiled_model <- keras::compile(
      model,
      optimizer = as.character(hyperparameters[i, 13]),
      loss = 'binary_crossentropy',
      metrics = 'accuracy'
    )
    
    return(compiled_model)
    
  }, error = function(e) {
    stop(paste("Error creating Keras model:", e$message))
  })
}

#' Deep Learning Training for miRNA Biomarker Selection
#'
#' Comprehensive deep learning training pipeline with hyperparameter optimization,
#' autoencoder support, and extensive model evaluation for miRNA biomarker analysis.
#'
#' @param selected_miRNAs Character vector of miRNA names, or '.' for all hsa miRNAs
#' @param wd Working directory containing mixed_train.csv, mixed_test.csv, mixed_validation.csv
#' @param SMOTE Logical, whether to use SMOTE for class balancing
#' @param keras_batch_size Integer, batch size for neural network training
#' @param clean_temp_files Logical, whether to clean temporary files after completion
#' @param save_threshold_trainacc Minimum training accuracy threshold to save model
#' @param save_threshold_testacc Minimum test accuracy threshold to save model
#' @param keras_epochae Maximum epochs for autoencoder training
#' @param keras_epoch Maximum epochs for classification model training
#' @param keras_patience Early stopping patience (epochs without improvement)
#' @param hyperparameters Dataframe with hyperparameter grid for model search
#' @param add_features_to_predictions Whether to include all features in prediction output
#' @param keras_threads Number of parallel threads for model training
#' @param start Starting index in hyperparameter grid
#' @param end Ending index in hyperparameter grid
#' @param output_file Output CSV filename for results
#' @param save_all_vars Whether to save all variables in model workspace
#' @param automatic_weight Whether to use automatic class weighting
#'
#' @return List containing:
#' \describe{
#'   \item{results}{Dataframe with performance metrics for each model}
#'   \item{best_models}{Information about best performing models}
#'   \item{model_paths}{Paths to saved model files}
#'   \item{status}{Completion status and summary}
#' }
#'
#' @details
#' This function implements a comprehensive deep learning pipeline that includes:
#' \itemize{
#'   \item Hyperparameter grid search across network architectures
#'   \item Optional autoencoder pre-training for feature reduction
#'   \item Parallel model training with early stopping
#'   \item Comprehensive model evaluation on train/test/validation sets
#'   \item Automatic model saving based on performance thresholds
#'   \item ROC analysis and optimal cutpoint determination
#'   \item Support for class balancing via SMOTE or automatic weighting
#' }
#'
#' The hyperparameters dataframe should contain columns for network architecture:
#' layer1, layer2, layer3 (neuron counts), activation functions, dropout rates,
#' regularization flags, optimizer choice, autoencoder settings, and scaling options.
#'
#' @examples
#' \dontrun{
#' # Basic deep learning with default hyperparameters
#' results <- train_deep_learning(
#'   selected_miRNAs = '.',  # Use all hsa_ miRNAs
#'   wd = 'data/',
#'   keras_epoch = 500,
#'   keras_threads = 2
#' )
#' 
#' # Advanced hyperparameter search
#' custom_grid <- expand.grid(
#'   layer1 = c(32, 64, 128),
#'   layer2 = c(0, 16, 32),
#'   layer3 = 0,
#'   activation_function_layer1 = c('relu', 'sigmoid'),
#'   activation_function_layer2 = 'relu',
#'   activation_function_layer3 = 'relu',
#'   dropout_layer1 = c(0, 0.1, 0.2),
#'   dropout_layer2 = 0,
#'   dropout_layer3 = 0,
#'   layer1_regularizer = c(TRUE, FALSE),
#'   layer2_regularizer = FALSE,
#'   layer3_regularizer = FALSE,
#'   optimizer = c('adam', 'rmsprop'),
#'   autoencoder = c(0, 7),  # 0 = no autoencoder, 7 = 7-neuron bottleneck
#'   balanced = FALSE,
#'   scaled = TRUE,
#'   stringsAsFactors = FALSE
#' )
#' 
#' results <- train_deep_learning(
#'   hyperparameters = custom_grid,
#'   save_threshold_trainacc = 0.90,
#'   save_threshold_testacc = 0.85
#' )
#' }
#'
#' @seealso \code{\link{create_keras_model}}, \code{\link{predict_deep_learning}}
#' @export
train_deep_learning <- function(selected_miRNAs = '.', 
                                wd = getwd(),
                                SMOTE = FALSE, 
                                keras_batch_size = 64,
                                clean_temp_files = TRUE,
                                save_threshold_trainacc = 0.85,
                                save_threshold_testacc = 0.8,
                                keras_epochae = 1000,
                                keras_epoch = 500,
                                keras_patience = 50,
                                hyperparameters = NULL,
                                add_features_to_predictions = FALSE,
                                keras_threads = ceiling(parallel::detectCores()/2),
                                start = 1,
                                end = NULL,
                                output_file = 'deeplearning_results.csv',
                                save_all_vars = FALSE,
                                automatic_weight = FALSE) {

  # Input validation and package checks
  if (!requireNamespace('keras', quietly = TRUE)) {
    stop('Package "keras" is required for deep learning functionality. Please install it with: install.packages("keras")')
  }
  
  required_packages <- c('foreach', 'doParallel', 'pROC', 'caret', 'cutpointr', 'ggplot2', 'gridExtra')
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop('Required packages missing: ', paste(missing_packages, collapse = ', '))
  }
  
  # Load required libraries
  library(keras)
  library(foreach)
  library(doParallel)
  library(pROC)
  library(caret)
  library(ggplot2)
  library(gridExtra)
  
  message('OmicSelector: Starting comprehensive deep learning pipeline...')
  
  # Create default hyperparameters if not provided
  if (is.null(hyperparameters)) {
    hyperparameters <- expand.grid(
      layer1 = seq(32, 128, by = 32),
      layer2 = c(0, seq(16, 64, by = 16)),
      layer3 = c(0, seq(8, 32, by = 8)),
      activation_function_layer1 = c('relu', 'sigmoid'),
      activation_function_layer2 = c('relu', 'sigmoid'),
      activation_function_layer3 = c('relu', 'sigmoid'),
      dropout_layer1 = c(0, 0.1, 0.2),
      dropout_layer2 = c(0, 0.1),
      dropout_layer3 = c(0),
      layer1_regularizer = c(TRUE, FALSE),
      layer2_regularizer = c(TRUE, FALSE),
      layer3_regularizer = c(FALSE),
      optimizer = c('adam', 'rmsprop', 'sgd'),
      autoencoder = c(0, 7, -7),  # 0=none, 7=regular, -7=sparse
      balanced = SMOTE,
      scaled = c(TRUE, FALSE),
      stringsAsFactors = FALSE
    )
    
    # Reduce grid size for reasonable computation time
    if (nrow(hyperparameters) > 100) {
      hyperparameters <- hyperparameters[sample(nrow(hyperparameters), 100), ]
      message('Reduced hyperparameter grid to 100 combinations for feasible computation')
    }
  }
  
  if (is.null(end)) end <- nrow(hyperparameters)
  
  # Setup directories and working environment
  old_wd <- getwd()
  setwd(wd)
  on.exit(setwd(old_wd))
  
  codename <- sub(pattern = '(.*)[.].*$', replacement = '\\1', basename(output_file))
  
  # Create necessary directories
  if (!dir.exists('temp')) dir.create('temp')
  if (!dir.exists('models')) dir.create('models')
  if (!dir.exists(paste0('models/', codename))) dir.create(paste0('models/', codename))
  
  # Determine temp directory
  if (dir.exists('/OmicSelector')) {
    temp_dir <- file.path(getwd(), 'temp-deeplearning')
    if (!dir.exists(temp_dir)) dir.create(temp_dir)
  } else {
    temp_dir <- tempdir()
  }
  
  message('Using temporary directory: ', temp_dir)
  
  # Save hyperparameters
  data.table::fwrite(hyperparameters, paste0('hyperparameters_', output_file))
  
  # Load datasets using OmicSelector data loading
  tryCatch({
    if (requireNamespace('OmicSelector', quietly = TRUE)) {
      dane <- OmicSelector::OmicSelector_load_datamix(wd = wd, replace_smote = FALSE, remove_zero_var = FALSE)
      train <- dane[[1]]
      test <- dane[[2]] 
      valid <- dane[[3]]
      train_smoted <- dane[[4]]
    } else {
      # Fallback data loading
      train <- data.table::fread('mixed_train.csv')
      test <- data.table::fread('mixed_test.csv')
      valid <- data.table::fread('mixed_validation.csv')
      train_smoted <- train  # Placeholder
    }
  }, error = function(e) {
    stop('Could not load datasets. Ensure mixed_train.csv, mixed_test.csv, and mixed_validation.csv exist in working directory.')
  })
  
  if (SMOTE) train <- train_smoted
  
  message('Data loaded successfully:')
  message('  Training samples: ', nrow(train))
  message('  Test samples: ', nrow(test))
  message('  Validation samples: ', nrow(valid))
  
  # Set up parallel processing
  message('Setting up parallel processing with ', keras_threads, ' cores...')
  cl <- makeCluster(keras_threads, outfile = file.path('temp', paste0(ceiling(as.numeric(Sys.time())), 'deeplearning_cluster.log')))
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  
  # Main parallel training loop
  message('Starting parallel hyperparameter search...')
  message('Processing hyperparameter sets ', start, ' to ', end, ' (', end-start+1, ' total)')
  
  set.seed(42)  # For reproducibility
  
  final_results <- foreach(i = start:end, 
                          .combine = rbind, 
                          .verbose = FALSE,
                          .inorder = FALSE, 
                          .errorhandling = 'remove',
                          .export = ls(envir = environment()),
                          .packages = c('keras', 'ggplot2', 'dplyr', 'data.table', 'pROC', 'caret', 'cutpointr', 'gridExtra')) %dopar% {
    
    # Individual model training function
    train_single_model(i, hyperparameters, train, test, valid, selected_miRNAs, 
                      temp_dir, keras_batch_size, keras_epochae, keras_epoch, 
                      keras_patience, save_threshold_trainacc, save_threshold_testacc,
                      add_features_to_predictions, save_all_vars, automatic_weight,
                      codename, old_wd, SMOTE)
  }
  
  # Process results
  if (!is.null(final_results) && nrow(final_results) > 0) {
    # Save comprehensive results
    data.table::fwrite(final_results, output_file)
    
    # Identify best models
    if ('test_Accuracy' %in% colnames(final_results)) {
      best_test_idx <- which.max(final_results$test_Accuracy)
      best_train_idx <- which.max(final_results$training_Accuracy)
      
      best_models <- list(
        best_test = final_results[best_test_idx, ],
        best_train = final_results[best_train_idx, ]
      )
    } else {
      best_models <- NULL
    }
    
    message('Deep learning training completed successfully!')
    message('Results saved to: ', output_file)
    message('Models trained: ', nrow(final_results))
    
    return(list(
      results = final_results,
      best_models = best_models,
      hyperparameters = hyperparameters,
      output_file = output_file,
      temp_dir = temp_dir,
      status = 'completed'
    ))
    
  } else {
    warning('No models completed successfully')
    return(list(
      status = 'failed',
      message = 'No models completed training',
      hyperparameters = hyperparameters
    ))
  }
}

#' Train Single Deep Learning Model (Internal Helper)
#'
#' Internal function to train a single model with given hyperparameters.
#' This function encapsulates the complete training pipeline for one model.
#'
#' @param i Index of hyperparameter set
#' @param hyperparameters Hyperparameter dataframe
#' @param train Training dataset
#' @param test Test dataset  
#' @param valid Validation dataset
#' @param selected_miRNAs Selected miRNA features
#' @param temp_dir Temporary directory for model files
#' @param keras_batch_size Batch size for training
#' @param keras_epochae Autoencoder epochs
#' @param keras_epoch Classification epochs
#' @param keras_patience Early stopping patience
#' @param save_threshold_trainacc Training accuracy threshold
#' @param save_threshold_testacc Test accuracy threshold
#' @param add_features_to_predictions Whether to include features in output
#' @param save_all_vars Whether to save all workspace variables
#' @param automatic_weight Whether to use automatic class weighting
#' @param codename Model codename for saving
#' @param old_wd Original working directory
#' @param SMOTE Whether SMOTE was used
#'
#' @return Dataframe row with model results
#' @keywords internal
train_single_model <- function(i, hyperparameters, train, test, valid, selected_miRNAs,
                              temp_dir, keras_batch_size, keras_epochae, keras_epoch,
                              keras_patience, save_threshold_trainacc, save_threshold_testacc,
                              add_features_to_predictions, save_all_vars, automatic_weight,
                              codename, old_wd, SMOTE) {
  
  library(keras)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(pROC)
  library(caret)
  library(cutpointr)
  library(gridExtra)
  
  set.seed(42)
  Sys.setenv(TF_FORCE_GPU_ALLOW_GROWTH = 'true')
  
  start_time <- Sys.time()
  
  # Generate model ID
  model_id <- paste0(format(i, scientific = FALSE), "-", ceiling(as.numeric(Sys.time())))
  if (SMOTE) model_id <- paste0(format(i, scientific = FALSE), "-SMOTE-", ceiling(as.numeric(Sys.time())))
  
  # Create model directory
  model_dir <- file.path(temp_dir, "models", paste0("keras", model_id))
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
  
  message("Training model: ", model_id)
  print(hyperparameters[i, ])
  
  # Initialize results
  temp_results <- data.frame(model_id = model_id)
  
  tryCatch({
    # Prepare data
    x_train <- train %>%
      {if (selected_miRNAs[1] != ".") dplyr::select(., all_of(selected_miRNAs)) else dplyr::select(., starts_with("hsa"))} %>%
      as.matrix()
    
    y_train <- train %>%
      dplyr::select("Class") %>%
      as.matrix()
    y_train[, 1] <- ifelse(y_train[, 1] == "Case", 1, 0)
    
    x_test <- test %>%
      {if (selected_miRNAs[1] != ".") dplyr::select(., all_of(selected_miRNAs)) else dplyr::select(., starts_with("hsa"))} %>%
      as.matrix()
    
    y_test <- test %>%
      dplyr::select("Class") %>%
      as.matrix()
    y_test[, 1] <- ifelse(y_test[, 1] == "Case", 1, 0)
    
    x_valid <- valid %>%
      {if (selected_miRNAs[1] != ".") dplyr::select(., all_of(selected_miRNAs)) else dplyr::select(., starts_with("hsa"))} %>%
      as.matrix()
    
    y_valid <- valid %>%
      dplyr::select("Class") %>%
      as.matrix()
    y_valid[, 1] <- ifelse(y_valid[, 1] == "Case", 1, 0)
    
    # Handle scaling
    if (hyperparameters[i, "scaled"]) {
      x_train_scale <- scale(x_train)
      col_mean_train <- attr(x_train_scale, "scaled:center")
      col_sd_train <- attr(x_train_scale, "scaled:scale")
      
      saveRDS(col_mean_train, file.path(model_dir, "col_mean_train.RDS"))
      saveRDS(col_sd_train, file.path(model_dir, "col_sd_train.RDS"))
      
      x_test_scale <- scale(x_test, center = col_mean_train, scale = col_sd_train)
      x_valid_scale <- scale(x_valid, center = col_mean_train, scale = col_sd_train)
    } else {
      x_train_scale <- x_train
      x_test_scale <- x_test
      x_valid_scale <- x_valid
    }
    
    # Handle automatic weighting
    if (automatic_weight) {
      counter <- table(y_train)
      majority <- max(counter)
      weights <- majority / counter
      class_weights <- as.list(weights)
      names(class_weights) <- names(counter)
    } else {
      class_weights <- NULL
    }
    
    # Autoencoder preprocessing if specified
    autoencoder_value <- hyperparameters[i, "autoencoder"]
    if (autoencoder_value != 0) {
      # Build and train autoencoder
      n1 <- abs(autoencoder_value)
      n3 <- ncol(x_train_scale)
      n2 <- max(7, ceiling(n3/2))
      
      input_layer <- layer_input(shape = c(n3))
      
      if (autoencoder_value > 0) {
        # Regular autoencoder
        encoder <- input_layer %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i, "activation_function_layer3"])) %>%
          layer_dense(units = n1, activation = "softmax")
        
        decoder <- encoder %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i, "activation_function_layer3"])) %>%
          layer_dense(units = n3, activation = as.character(hyperparameters[i, "activation_function_layer3"]))
      } else {
        # Sparse autoencoder  
        encoder <- input_layer %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i, "activation_function_layer3"]),
                     kernel_regularizer = regularizer_l1(l = 0.01)) %>%
          layer_dense(units = n1, activation = "softmax",
                     kernel_regularizer = regularizer_l1(l = 0.01))
        
        decoder <- encoder %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i, "activation_function_layer3"]),
                     kernel_regularizer = regularizer_l1(l = 0.01)) %>%
          layer_dense(units = n3, activation = as.character(hyperparameters[i, "activation_function_layer3"]))
      }
      
      ae_model <- keras_model(inputs = input_layer, outputs = decoder)
      ae_model %>% compile(
        loss = "mean_absolute_error",
        optimizer = as.character(hyperparameters[i, "optimizer"]),
        metrics = c("mean_squared_error")
      )
      
      # Callbacks
      ae_early_stop <- callback_early_stopping(monitor = "val_loss", mode = "min", patience = keras_patience)
      ae_checkpoint <- callback_model_checkpoint(
        filepath = file.path(model_dir, "autoencoderweights.hdf5"),
        save_best_only = TRUE,
        save_weights_only = TRUE,
        monitor = "val_loss",
        verbose = 0
      )
      
      # Train autoencoder
      ae_history <- fit(ae_model,
                       x = x_train_scale,
                       y = x_train_scale,
                       epochs = keras_epochae,
                       batch_size = keras_batch_size,
                       shuffle = TRUE,
                       verbose = 0,
                       validation_data = list(x_test_scale, x_test_scale),
                       callbacks = list(ae_checkpoint, ae_early_stop))
      
      saveRDS(ae_history, file.path(model_dir, "ae_history.RDS"))
      
      # Create encoder model and extract features
      encoder_model <- keras_model(inputs = input_layer, outputs = encoder)
      encoder_model %>% load_model_weights_hdf5(file.path(model_dir, "autoencoderweights.hdf5"),
                                               skip_mismatch = TRUE, by_name = TRUE)
      
      save_model_hdf5(encoder_model, file.path(model_dir, "autoencoder.hdf5"))
      
      # Generate deep features
      x_train_scale <- predict(encoder_model, x_train_scale)
      x_test_scale <- predict(encoder_model, x_test_scale)
      x_valid_scale <- predict(encoder_model, x_valid_scale)
      
      data.table::fwrite(x_train_scale, file.path(model_dir, "deepfeatures_train.csv"))
      data.table::fwrite(x_test_scale, file.path(model_dir, "deepfeatures_test.csv"))
      data.table::fwrite(x_valid_scale, file.path(model_dir, "deepfeatures_valid.csv"))
    }
    
    # Build classification model
    dnn_model <- create_keras_model(i, hyperparameters, ncol(x_train_scale))
    
    # Callbacks for classification training
    early_stop <- callback_early_stopping(monitor = "val_loss", mode = "min", patience = keras_patience)
    checkpoint <- callback_model_checkpoint(
      filepath = file.path(model_dir, "finalmodel.hdf5"),
      save_best_only = TRUE,
      monitor = "val_loss",
      verbose = 0
    )
    
    # Train classification model
    history <- fit(dnn_model,
                  x = x_train_scale,
                  y = to_categorical(y_train),
                  epochs = keras_epoch,
                  batch_size = keras_batch_size,
                  shuffle = TRUE,
                  verbose = 0,
                  validation_data = list(x_test_scale, to_categorical(y_test)),
                  callbacks = list(checkpoint, early_stop),
                  class_weight = class_weights)
    
    saveRDS(history, file.path(model_dir, "history.RDS"))
    
    # Load best model and make predictions
    dnn_model <- load_model_hdf5(file.path(model_dir, "finalmodel.hdf5"))
    
    y_train_pred <- predict(dnn_model, x_train_scale)
    y_test_pred <- predict(dnn_model, x_test_scale)
    y_valid_pred <- predict(dnn_model, x_valid_scale)
    
    # Determine optimal cutoff using training data
    pred_df <- data.frame(Class = train$Class, Pred = y_train_pred[, 2])
    cutoff_result <- cutpointr::cutpointr(pred_df, Pred, Class, pos_class = "Case", metric = youden)
    optimal_cutoff <- cutoff_result$optimal_cutpoint
    
    temp_results[1, "training_AUC"] <- cutoff_result$AUC
    temp_results[1, "cutoff"] <- optimal_cutoff
    
    # Training performance
    pred_df$PredClass <- ifelse(pred_df$Pred >= optimal_cutoff, "Case", "Control")
    pred_df$PredClass <- factor(pred_df$PredClass, levels = c("Control", "Case"))
    cm_train <- confusionMatrix(pred_df$PredClass, pred_df$Class, positive = "Case")
    
    temp_results[1, "training_Accuracy"] <- cm_train$overall[1]
    temp_results[1, "training_Sensitivity"] <- cm_train$byClass[1]
    temp_results[1, "training_Specificity"] <- cm_train$byClass[2]
    temp_results[1, "training_PPV"] <- cm_train$byClass[3]
    temp_results[1, "training_NPV"] <- cm_train$byClass[4]
    temp_results[1, "training_F1"] <- cm_train$byClass[7]
    
    # Test performance
    pred_df <- data.frame(Class = test$Class, Pred = y_test_pred[, 2])
    pred_df$PredClass <- ifelse(pred_df$Pred >= optimal_cutoff, "Case", "Control")
    pred_df$PredClass <- factor(pred_df$PredClass, levels = c("Control", "Case"))
    cm_test <- confusionMatrix(pred_df$PredClass, pred_df$Class, positive = "Case")
    
    temp_results[1, "test_Accuracy"] <- cm_test$overall[1]
    temp_results[1, "test_Sensitivity"] <- cm_test$byClass[1]
    temp_results[1, "test_Specificity"] <- cm_test$byClass[2]
    temp_results[1, "test_PPV"] <- cm_test$byClass[3]
    temp_results[1, "test_NPV"] <- cm_test$byClass[4]
    temp_results[1, "test_F1"] <- cm_test$byClass[7]
    
    # Validation performance
    pred_df <- data.frame(Class = valid$Class, Pred = y_valid_pred[, 2])
    pred_df$PredClass <- ifelse(pred_df$Pred >= optimal_cutoff, "Case", "Control")
    pred_df$PredClass <- factor(pred_df$PredClass, levels = c("Control", "Case"))
    cm_valid <- confusionMatrix(pred_df$PredClass, pred_df$Class, positive = "Case")
    
    temp_results[1, "valid_Accuracy"] <- cm_valid$overall[1]
    temp_results[1, "valid_Sensitivity"] <- cm_valid$byClass[1]
    temp_results[1, "valid_Specificity"] <- cm_valid$byClass[2]
    temp_results[1, "valid_PPV"] <- cm_valid$byClass[3]
    temp_results[1, "valid_NPV"] <- cm_valid$byClass[4]
    temp_results[1, "valid_F1"] <- cm_valid$byClass[7]
    
    # Save model files and results
    data.table::fwrite(cbind(hyperparameters[i, ], temp_results), file.path(model_dir, "wyniki.csv"))
    
    # Save weights and model
    model_weights <- get_weights(dnn_model)
    saveRDS(model_weights, file.path(model_dir, "finalmodel_weights.RDS"))
    save_model_weights_hdf5(dnn_model, file.path(model_dir, "finalmodel_weights.hdf5"))
    saveRDS(dnn_model, file.path(model_dir, "finalmodel.RDS"))
    
    # Create training plots
    train_loss <- history$metrics$loss
    val_loss <- history$metrics$val_loss
    train_acc <- history$metrics$accuracy
    val_acc <- history$metrics$val_accuracy
    
    if (!is.null(train_loss) && !is.null(val_loss)) {
      loss_df <- data.frame(
        epoch = seq_along(train_loss),
        training = train_loss,
        validation = val_loss
      ) %>%
        tidyr::gather(key = "type", value = "loss", -epoch)
      
      loss_plot <- ggplot(loss_df, aes(x = epoch, y = loss, color = type)) +
        geom_line() +
        labs(title = "Model Loss", x = "Epoch", y = "Loss") +
        theme_minimal()
    }
    
    if (!is.null(train_acc) && !is.null(val_acc)) {
      acc_df <- data.frame(
        epoch = seq_along(train_acc),
        training = train_acc,
        validation = val_acc
      ) %>%
        tidyr::gather(key = "type", value = "accuracy", -epoch)
      
      acc_plot <- ggplot(acc_df, aes(x = epoch, y = accuracy, color = type)) +
        geom_line() +
        labs(title = "Model Accuracy", x = "Epoch", y = "Accuracy") +
        theme_minimal()
      
      if (exists("loss_plot")) {
        combined_plot <- gridExtra::grid.arrange(loss_plot, acc_plot, nrow = 2)
        ggplot2::ggsave(file.path(model_dir, "training.png"), combined_plot)
      }
    }
    
    # Save model if it meets thresholds
    train_acc_val <- temp_results[1, "training_Accuracy"]
    test_acc_val <- temp_results[1, "test_Accuracy"]
    
    if (!is.na(train_acc_val) && !is.na(test_acc_val) &&
        train_acc_val > save_threshold_trainacc && test_acc_val > save_threshold_testacc) {
      
      # Create zip archive of model
      model_files <- list.files(model_dir, full.names = TRUE, recursive = TRUE)
      zip_path <- file.path(old_wd, "models", codename, paste0(codename, "_", model_id, ".zip"))
      
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zip(zip_path, model_files, mode = "cherry-pick")
      } else {
        # Fallback to system zip
        system(paste("cd", dirname(model_dir), "&& zip -r", zip_path, basename(model_dir)))
      }
      
      message("Model ", model_id, " saved: train_acc=", round(train_acc_val, 3), 
              " test_acc=", round(test_acc_val, 3))
    }
    
    # Return results combined with hyperparameters
    return(cbind(hyperparameters[i, ], temp_results))
    
  }, error = function(e) {
    message("Error training model ", model_id, ": ", e$message)
    return(data.frame(model_id = model_id, error = e$message))
  })
}

#' Predict Using Trained Deep Learning Model
#'
#' Makes predictions on new data using a trained deep learning model.
#' Handles scaling, autoencoder preprocessing, and provides comprehensive evaluation.
#'
#' @param model_path Path to ZIP file containing trained model
#' @param new_dataset New dataset for prediction
#' @param new_scaling Logical, perform new scaling on data
#' @param old_train_csv_to_restore_scaling Path to original training data for scaling
#' @param override_cutoff Numeric, custom classification cutoff
#' @param blinded Logical, whether validation is blinded (no true labels)
#'
#' @return List containing predictions, performance metrics, and model information
#'
#' @details
#' This function loads a trained model and makes predictions on new data.
#' It automatically handles:
#' - Feature matching between training and prediction data
#' - Data scaling using original parameters or new scaling
#' - Autoencoder preprocessing if used in training
#' - Classification using optimal cutoff or custom threshold
#' - Performance evaluation if true labels are available
#'
#' @examples
#' \dontrun{
#' predictions <- predict_deep_learning(
#'   model_path = "trained_model.zip",
#'   new_dataset = new_data,
#'   blinded = FALSE
#' )
#' print(predictions$confusion_matrix)
#' print(predictions$roc_auc)
#' }
#'
#' @export
predict_deep_learning <- function(model_path = "our_models/model5.zip",
                                 new_dataset = NULL,
                                 new_scaling = FALSE,
                                 old_train_csv_to_restore_scaling = NULL,
                                 override_cutoff = NULL,
                                 blinded = FALSE) {
  
  # Check required packages
  required_packages <- c("keras", "data.table", "pROC", "caret", "dplyr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop(paste("Required packages missing:", paste(missing_packages, collapse = ", ")))
  }
  
  # Load required libraries
  library(keras)
  library(data.table)
  library(pROC)
  library(caret)
  library(dplyr)
  
  # Validate inputs
  if (!file.exists(model_path)) {
    stop(paste("Model file does not exist:", model_path))
  }
  
  if (is.null(new_dataset)) {
    stop("New dataset is required for prediction")
  }
  
  message("Loading model and making predictions...")
  
  tryCatch({
    # Load model from ZIP
    model_files <- utils::unzip(model_path, list = TRUE)
    model_path_in_zip <- model_files[grepl("finalmodel.hdf5", model_files$Name), "Name"][1]
    
    if (is.na(model_path_in_zip)) {
      stop("No finalmodel.hdf5 found in model ZIP file")
    }
    
    utils::unzip(model_path, model_path_in_zip, exdir = tempdir())
    model_path_unzipped <- file.path(tempdir(), model_path_in_zip)
    
    if (!file.exists(model_path_unzipped)) {
      stop("Model file does not exist in provided ZIP. Is it correct?")
    }
    
    init_model <- keras::load_model_hdf5(model_path_unzipped)
    message("Model loaded successfully")
    print(init_model)
    
    # Load model configuration
    config_files <- model_files[grepl("wyniki.csv", model_files$Name), "Name"]
    if (length(config_files) == 0) {
      stop("Configuration file (wyniki.csv) not found in model package")
    }
    
    utils::unzip(model_path, config_files[1], exdir = tempdir())
    config_path_unzipped <- file.path(tempdir(), config_files[1])
    pre_conf <- data.table::fread(config_path_unzipped)
    
    # Extract network features from formula
    if (!"formula" %in% colnames(pre_conf)) {
      stop("Model configuration missing formula specification")
    }
    
    network_features <- strsplit(x = as.character(pre_conf[1, "formula"]), split = " + ", fixed = TRUE)[[1]]
    
    # Check feature availability in new dataset
    missing_features <- network_features[!network_features %in% colnames(new_dataset)]
    if (length(missing_features) > 0) {
      stop(paste0("The new dataset does not contain features: ", paste(missing_features, collapse = ", ")))
    }
    
    # Prepare feature matrix
    new_x <- new_dataset %>% 
      dplyr::select(all_of(network_features)) %>% 
      as.matrix()
    
    # Prepare labels if not blinded
    if (!blinded) {
      if (!"Class" %in% colnames(new_dataset)) {
        stop("Class column required for non-blinded validation")
      }
      
      new_y <- new_dataset %>%
        dplyr::select("Class") %>%
        as.matrix()
      new_y[, 1] <- ifelse(new_y[, 1] == "Case", 1, 0)
      new_dataset$Class <- factor(new_dataset$Class, levels = c("Control", "Case"))
    }
    
    # Handle scaling if model was trained with scaling
    col_mean_train <- NA
    col_sd_train <- NA
    
    if (as.character(pre_conf[1, "scaled"]) == "TRUE") {
      message("This network was trained with scaling")
      
      # Try to load scaling parameters from model
      mean_files <- model_files[grepl("col_mean_train.RDS", model_files$Name), "Name"]
      if (length(mean_files) > 0) {
        utils::unzip(model_path, mean_files[1], exdir = tempdir())
        mean_path_unzipped <- file.path(tempdir(), mean_files[1])
        if (file.exists(mean_path_unzipped)) {
          col_mean_train <- readRDS(mean_path_unzipped)
        } else {
          message("Scaling mean not saved in model")
        }
      }
      
      sd_files <- model_files[grepl("col_sd_train.RDS", model_files$Name), "Name"]
      if (length(sd_files) > 0) {
        utils::unzip(model_path, sd_files[1], exdir = tempdir())
        sd_path_unzipped <- file.path(tempdir(), sd_files[1])
        if (file.exists(sd_path_unzipped)) {
          col_sd_train <- readRDS(sd_path_unzipped)
        } else {
          message("Scaling SD not saved in model")
        }
      }
      
      # Apply scaling
      if (new_scaling) {
        message("Applying new scaling")
        new_x <- scale(new_x)
        new_x[is.nan(new_x)] <- 0
        message("New scaling performed with col_mean: ", paste(round(attr(new_x, "scaled:center"), 4), collapse = ", "))
        message("New col_sd: ", paste(round(attr(new_x, "scaled:scale"), 4), collapse = ", "))
      } else {
        if (is.null(old_train_csv_to_restore_scaling)) {
          # Use saved scaling parameters
          if (!any(is.na(col_mean_train)) && !any(is.na(col_sd_train))) {
            new_x <- scale(new_x, center = col_mean_train, scale = col_sd_train)
            new_x[is.nan(new_x)] <- 0
          } else {
            warning("Scaling parameters not available, proceeding without scaling")
          }
        } else {
          # Restore scaling from provided file
          if (!file.exists(old_train_csv_to_restore_scaling)) {
            stop("A file to restore scaling from does not exist")
          }
          
          message("Restoring scaling from provided file")
          temp_train <- data.table::fread(old_train_csv_to_restore_scaling)
          temp_train_x <- dplyr::select(temp_train, all_of(network_features))
          temp_train_x <- scale(temp_train_x)
          col_mean_train <- attr(temp_train_x, "scaled:center")
          col_sd_train <- attr(temp_train_x, "scaled:scale")
          
          new_x <- scale(new_x, center = col_mean_train, scale = col_sd_train)
          new_x[is.nan(new_x)] <- 0
          
          message("Scaling restored. Final col_mean: ", paste(round(col_mean_train, 4), collapse = ", "))
          message("Final col_sd: ", paste(round(col_sd_train, 4), collapse = ", "))
        }
      }
    }
    
    # Handle autoencoder if present
    model_autoencoder <- NULL
    if (as.character(pre_conf[1, "autoencoder"]) != "0") {
      message("Loading autoencoder model")
      
      ae_files <- model_files[grepl("autoencoder.hdf5", model_files$Name), "Name"]
      if (length(ae_files) > 0) {
        utils::unzip(model_path, ae_files[1], exdir = tempdir())
        ae_path_unzipped <- file.path(tempdir(), ae_files[1])
        
        if (!file.exists(ae_path_unzipped)) {
          stop("Autoencoder model file does not exist in provided file")
        }
        
        model_autoencoder <- keras::load_model_hdf5(ae_path_unzipped)
        print(model_autoencoder)
        
        # Transform features using autoencoder
        new_x_original <- new_x
        new_x <- predict(model_autoencoder, new_x) %>% as.matrix()
        message("Features transformed using autoencoder: ", ncol(new_x_original), " -> ", ncol(new_x))
      } else {
        warning("Autoencoder specified but not found in model files")
      }
    }
    
    # Get cutoff for prediction
    cutoff <- as.numeric(pre_conf[1, "cutoff"])
    if (!is.null(override_cutoff)) {
      cutoff <- as.numeric(override_cutoff)
      message("Using override cutoff: ", cutoff)
    }
    
    if (is.na(cutoff)) {
      stop("Cutoff is not numeric. Check your override_cutoff parameter or the value set in the model")
    }
    
    message("Using cutoff: ", cutoff)
    
    # Make predictions
    predictions <- predict(init_model, new_x)
    
    if (blinded) {
      # Blinded prediction (no true labels)
      pred <- data.frame(Pred = predictions[, 2])
      pred$Prediction <- ifelse(pred$Pred >= cutoff, "Case", "Control")
      confusion_matrix <- NA
      roc_result <- NA
      roc_auc <- NA
      
      message("Blinded predictions completed")
    } else {
      # Full evaluation with true labels
      pred <- data.frame(
        Class = factor(ifelse(new_y == 1, "Case", "Control"), levels = c("Control", "Case")),
        Pred = predictions[, 2]
      )
      pred$Prediction <- ifelse(pred$Pred >= cutoff, "Case", "Control")
      pred$Correctness <- ifelse(pred$Prediction == pred$Class, "Correct", "Incorrect")
      
      # Calculate confusion matrix
      confusion_matrix <- caret::confusionMatrix(
        as.factor(pred$Prediction),
        as.factor(pred$Class),
        positive = "Case"
      )
      
      # Calculate ROC and AUC
      roc_result <- pROC::roc(pred$Class ~ pred$Pred, quiet = TRUE)
      roc_auc <- pROC::ci.auc(roc_result)
      
      message("Performance evaluation completed:")
      message("  Accuracy: ", round(confusion_matrix$overall["Accuracy"], 3))
      message("  Sensitivity: ", round(confusion_matrix$byClass["Sensitivity"], 3))
      message("  Specificity: ", round(confusion_matrix$byClass["Specificity"], 3))
      message("  AUC: ", round(as.numeric(roc_result$auc), 3))
    }
    
    # Compile comprehensive results
    final_return <- list(
      predictions = pred,
      network_config = pre_conf,
      new_dataset = new_dataset,
      new_dataset_x = new_x,
      network_features = network_features,
      cutoff = cutoff,
      col_mean_train = col_mean_train,
      col_sd_train = col_sd_train,
      confusion_matrix = confusion_matrix,
      roc = roc_result,
      roc_auc = roc_auc,
      model = init_model,
      autoencoder = model_autoencoder
    )
    
    message("Deep learning prediction completed successfully")
    return(final_return)
    
  }, error = function(e) {
    stop(paste("Error in deep learning prediction:", e$message))
  })
}

#' Transfer Learning Neural Network
#'
#' Implements transfer learning using a pre-trained neural network model.
#' Loads a trained model, adapts it to new data, and performs fine-tuning
#' with optional layer freezing for domain adaptation.
#'
#' @param selected_miRNAs Vector of miRNA names or "." for all hsa_ features
#' @param new_scaling Logical, whether to perform new scaling on the data
#' @param model_path Path to ZIP file containing pre-trained model
#' @param save_scaling Logical, whether to save scaling parameters
#' @param old_train_csv_to_restore_scaling Path to original training CSV for scaling restoration
#' @param freeze_from Integer, starting layer for weight freezing (0 = no freezing)
#' @param freeze_to Integer, ending layer for weight freezing
#' @param train Training dataset for transfer learning
#' @param test Test dataset for evaluation
#' @param valid Validation dataset for evaluation
#' @param keras_epoch Maximum epochs for fine-tuning
#' @param keras_batch_size Batch size for training
#' @param keras_patience Early stopping patience
#'
#' @return List containing:
#' \describe{
#'   \item{initial_performance}{Performance metrics of pre-trained model on new data}
#'   \item{transfer_performance}{Performance metrics after transfer learning}
#'   \item{model_paths}{Paths to saved models}
#'   \item{scaling_info}{Scaling parameters and configuration}
#'   \item{cutoffs}{Optimal classification cutoffs for both scenarios}
#' }
#'
#' @details
#' This function implements comprehensive transfer learning for neural networks:
#' \itemize{
#'   \item Loads pre-trained model and evaluates on new dataset
#'   \item Handles feature scaling using either new data or original training data
#'   \item Performs optional layer freezing for controlled fine-tuning
#'   \item Re-calibrates model on new training data
#'   \item Compares performance before and after transfer learning
#'   \item Saves all intermediate models and evaluation metrics
#' }
#'
#' The transfer learning process includes:
#' 1. Loading pre-trained model and evaluating initial performance
#' 2. Setting up appropriate data scaling strategy
#' 3. Freezing specified layers if requested
#' 4. Fine-tuning on new training data with early stopping
#' 5. Comprehensive evaluation on train/test/validation sets
#' 6. Saving final transfer-learned model
#'
#' @examples
#' \dontrun{
#' # Basic transfer learning
#' results <- transfer_learning_neural_network(
#'   selected_miRNAs = ".",
#'   model_path = "models/pretrained_model.zip",
#'   train = train_data,
#'   test = test_data,
#'   valid = valid_data
#' )
#' 
#' # Transfer learning with layer freezing
#' results <- transfer_learning_neural_network(
#'   model_path = "models/pretrained_model.zip",
#'   freeze_from = 1,
#'   freeze_to = 2,
#'   new_scaling = FALSE,
#'   keras_epoch = 1000
#' )
#' }
#'
#' @seealso \code{\link{train_deep_learning}}, \code{\link{predict_deep_learning}}
#' @export
transfer_learning_neural_network <- function(selected_miRNAs = ".", 
                                           new_scaling = TRUE,
                                           model_path = "tcga_models/pancreatic_tcga_165592-1601307872.zip",
                                           save_scaling = TRUE,
                                           old_train_csv_to_restore_scaling = "../tcga/mixed_train.csv",
                                           freeze_from = 1,
                                           freeze_to = 2,
                                           train = NULL,
                                           test = NULL,
                                           valid = NULL,
                                           keras_epoch = 5000,
                                           keras_batch_size = 64,
                                           keras_patience = 200) {
  
  # Input validation and package checks
  if (!requireNamespace("keras", quietly = TRUE)) {
    stop("Package 'keras' is required for transfer learning functionality. Please install it with: install.packages('keras')")
  }
  
  required_packages <- c("pROC", "caret", "cutpointr", "ggplot2", "dplyr", "data.table")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "))
  }
  
  # Load required libraries
  library(keras)
  library(pROC)
  library(caret)
  library(cutpointr)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  
  message("OmicSelector: Starting transfer learning neural network...")
  
  # Load datasets if not provided
  if (is.null(train) || is.null(test) || is.null(valid)) {
    message("Loading default datasets...")
    tryCatch({
      train <- data.table::fread("circ_data/mixed_train.csv")
      test <- data.table::fread("circ_data/mixed_test.csv")
      valid <- data.table::fread("circ_data/mixed_valid.csv")
    }, error = function(e) {
      stop("Could not load datasets. Please provide train, test, and valid datasets or ensure data files exist.")
    })
  }
  
  # Initialize results storage
  temp_results <- data.frame()
  
  # Prepare feature matrices
  x_train <- train %>%
    {if (selected_miRNAs[1] != ".") dplyr::select(., all_of(selected_miRNAs)) else dplyr::select(., starts_with("hsa"))} %>%
    as.matrix()
  
  y_train <- train %>%
    dplyr::select("Class") %>%
    as.matrix()
  y_train[, 1] <- ifelse(y_train[, 1] == "Case", 1, 0)
  
  x_test <- test %>%
    {if (selected_miRNAs[1] != ".") dplyr::select(., all_of(selected_miRNAs)) else dplyr::select(., starts_with("hsa"))} %>%
    as.matrix()
  
  y_test <- test %>%
    dplyr::select("Class") %>%
    as.matrix()
  y_test[, 1] <- ifelse(y_test[, 1] == "Case", 1, 0)
  
  x_valid <- valid %>%
    {if (selected_miRNAs[1] != ".") dplyr::select(., all_of(selected_miRNAs)) else dplyr::select(., starts_with("hsa"))} %>%
    as.matrix()
  
  y_valid <- valid %>%
    dplyr::select("Class") %>%
    as.matrix()
  y_valid[, 1] <- ifelse(y_valid[, 1] == "Case", 1, 0)
  
  message("Data prepared: train=", nrow(train), " test=", nrow(test), " valid=", nrow(valid))
  
  # Load pre-trained model
  message("Loading pre-trained model from: ", model_path)
  
  tryCatch({
    # Extract model from ZIP
    model_files <- utils::unzip(model_path, list = TRUE)
    model_path_in_zip <- model_files[grepl("finalmodel.hdf5", model_files$Name), "Name"][1]
    
    if (is.na(model_path_in_zip)) {
      stop("No finalmodel.hdf5 found in model ZIP file")
    }
    
    utils::unzip(model_path, model_path_in_zip, exdir = tempdir())
    model_path_unzipped <- file.path(tempdir(), model_path_in_zip)
    init_model <- keras::load_model_hdf5(model_path_unzipped)
    
    message("Pre-trained model loaded successfully")
    
    # Load model configuration
    config_files <- model_files[grepl("wyniki.csv", model_files$Name), "Name"]
    if (length(config_files) > 0) {
      utils::unzip(model_path, config_files[1], exdir = tempdir())
      config_path_unzipped <- file.path(tempdir(), config_files[1])
      pre_conf <- data.table::fread(config_path_unzipped)
    } else {
      warning("Model configuration not found, using defaults")
      pre_conf <- data.frame(scaled = TRUE, optimizer = "adam")
    }
    
    # Store initial weights
    init_weights <- keras::get_weights(init_model)
    
  }, error = function(e) {
    stop("Failed to load pre-trained model: ", e$message)
  })
  
  # Handle data scaling
  if (as.logical(pre_conf[1, "scaled"])) {
    if (new_scaling) {
      message("Applying new scaling based on current training data...")
      x_train_scale <- scale(x_train)
      col_mean_train <- attr(x_train_scale, "scaled:center")
      col_sd_train <- attr(x_train_scale, "scaled:scale")
      
      if (save_scaling) {
        saveRDS(col_mean_train, "transfer_col_mean_train.RDS")
        saveRDS(col_sd_train, "transfer_col_sd_train.RDS")
      }
      
      x_test_scale <- scale(x_test, center = col_mean_train, scale = col_sd_train)
      x_valid_scale <- scale(x_valid, center = col_mean_train, scale = col_sd_train)
      
    } else {
      message("Using original scaling from pre-trained model...")
      # Try to load scaling from model first
      scaling_files <- model_files[grepl("col_mean_train.RDS", model_files$Name), "Name"]
      
      if (length(scaling_files) > 0) {
        utils::unzip(model_path, scaling_files, exdir = tempdir())
        col_mean_train <- readRDS(file.path(tempdir(), scaling_files[1]))
        
        scaling_files_sd <- model_files[grepl("col_sd_train.RDS", model_files$Name), "Name"]
        if (length(scaling_files_sd) > 0) {
          utils::unzip(model_path, scaling_files_sd, exdir = tempdir())
          col_sd_train <- readRDS(file.path(tempdir(), scaling_files_sd[1]))
        } else {
          stop("Scaling parameters incomplete in model file")
        }
      } else {
        # Fallback to regenerating from old training data
        if (file.exists(old_train_csv_to_restore_scaling)) {
          message("Regenerating scaling from original training data...")
          tcga_train <- data.table::fread(old_train_csv_to_restore_scaling)
          tcga_train_matrix <- tcga_train %>%
            {if (selected_miRNAs[1] != ".") dplyr::select(., all_of(selected_miRNAs)) else dplyr::select(., starts_with("hsa"))} %>%
            as.matrix() %>%
            scale()
          
          col_mean_train <- attr(tcga_train_matrix, "scaled:center")
          col_sd_train <- attr(tcga_train_matrix, "scaled:scale")
        } else {
          stop("Cannot restore original scaling: file not found and scaling parameters not in model")
        }
      }
      
      if (save_scaling) {
        saveRDS(col_mean_train, "transfer_col_mean_train.RDS")
        saveRDS(col_sd_train, "transfer_col_sd_train.RDS")
      }
      
      x_train_scale <- scale(x_train, center = col_mean_train, scale = col_sd_train)
      x_test_scale <- scale(x_test, center = col_mean_train, scale = col_sd_train)
      x_valid_scale <- scale(x_valid, center = col_mean_train, scale = col_sd_train)
    }
  } else {
    message("No scaling applied (model was trained without scaling)")
    x_train_scale <- x_train
    x_test_scale <- x_test
    x_valid_scale <- x_valid
  }
  
  # Evaluate initial model performance on new data
  message("Evaluating pre-trained model on new data...")
  
  y_train_pred_init <- predict(init_model, x_train_scale)
  y_test_pred_init <- predict(init_model, x_test_scale)
  y_valid_pred_init <- predict(init_model, x_valid_scale)
  
  # Determine optimal cutoff for initial model
  pred_df <- data.frame(
    Class = as.factor(ifelse(y_train == 1, "Case", "Control")), 
    Pred = y_train_pred_init[, 2]
  )
  
  cutoff_initial <- cutpointr::cutpointr(pred_df, Pred, Class, pos_class = "Case", metric = youden)
  optimal_cutoff_initial <- cutoff_initial$optimal_cutpoint
  
  # Save initial cutoff plot
  ggplot2::ggsave("initial_cutoff.png", plot(cutoff_initial))
  
  # Initial performance metrics
  temp_results[1, "initial_training_AUC"] <- cutoff_initial$AUC
  temp_results[1, "initial_cutoff"] <- optimal_cutoff_initial
  
  # Training performance (initial)
  pred_df$PredClass <- ifelse(pred_df$Pred >= optimal_cutoff_initial, "Case", "Control")
  pred_df$PredClass <- factor(pred_df$PredClass, levels = c("Control", "Case"))
  cm_train_init <- confusionMatrix(pred_df$PredClass, pred_df$Class, positive = "Case")
  
  temp_results[1, "initial_training_Accuracy"] <- cm_train_init$overall[1]
  temp_results[1, "initial_training_Sensitivity"] <- cm_train_init$byClass[1]
  temp_results[1, "initial_training_Specificity"] <- cm_train_init$byClass[2]
  temp_results[1, "initial_training_F1"] <- cm_train_init$byClass[7]
  
  # Test performance (initial)
  pred_df_test <- data.frame(Class = as.factor(test$Class), Pred = y_test_pred_init[, 2])
  pred_df_test$PredClass <- ifelse(pred_df_test$Pred >= optimal_cutoff_initial, "Case", "Control")
  pred_df_test$PredClass <- factor(pred_df_test$PredClass, levels = c("Control", "Case"))
  cm_test_init <- confusionMatrix(pred_df_test$PredClass, pred_df_test$Class, positive = "Case")
  
  temp_results[1, "initial_test_Accuracy"] <- cm_test_init$overall[1]
  temp_results[1, "initial_test_Sensitivity"] <- cm_test_init$byClass[1]
  temp_results[1, "initial_test_Specificity"] <- cm_test_init$byClass[2]
  temp_results[1, "initial_test_F1"] <- cm_test_init$byClass[7]
  
  # Validation performance (initial)
  pred_df_valid <- data.frame(Class = as.factor(valid$Class), Pred = y_valid_pred_init[, 2])
  pred_df_valid$PredClass <- ifelse(pred_df_valid$Pred >= optimal_cutoff_initial, "Case", "Control")
  pred_df_valid$PredClass <- factor(pred_df_valid$PredClass, levels = c("Control", "Case"))
  cm_valid_init <- confusionMatrix(pred_df_valid$PredClass, pred_df_valid$Class, positive = "Case")
  
  temp_results[1, "initial_valid_Accuracy"] <- cm_valid_init$overall[1]
  temp_results[1, "initial_valid_Sensitivity"] <- cm_valid_init$byClass[1]
  temp_results[1, "initial_valid_Specificity"] <- cm_valid_init$byClass[2]
  temp_results[1, "initial_valid_F1"] <- cm_valid_init$byClass[7]
  
  # Save initial evaluation results
  saveRDS(cm_train_init, "initial_cm_train.RDS")
  saveRDS(cm_test_init, "initial_cm_test.RDS") 
  saveRDS(cm_valid_init, "initial_cm_valid.RDS")
  
  message("Initial model performance evaluated")
  
  # Setup transfer learning
  message("Setting up transfer learning...")
  
  # Clone model for transfer learning
  trans_model <- clone_model(init_model)
  
  # Set up callbacks
  early_stop <- callback_early_stopping(monitor = "val_loss", mode = "min", patience = keras_patience)
  checkpoint <- callback_model_checkpoint(
    filepath = "transfer_model.hdf5",
    save_best_only = TRUE,
    monitor = "val_loss",
    verbose = 0
  )
  
  # Compile transfer model
  optimizer_name <- if ("optimizer" %in% colnames(pre_conf)) as.character(pre_conf$optimizer[1]) else "adam"
  compile(trans_model, 
          loss = 'binary_crossentropy',
          metrics = 'accuracy', 
          optimizer = optimizer_name)
  
  # Set initial weights
  keras::set_weights(trans_model, init_weights)
  
  # Apply layer freezing if specified
  if (freeze_from > 0) {
    message("Freezing layers ", freeze_from, " to ", freeze_to)
    freeze_weights(trans_model, from = freeze_from, to = freeze_to)
    compile(trans_model,
            loss = 'binary_crossentropy',
            metrics = 'accuracy',
            optimizer = optimizer_name)
  }
  
  # Perform transfer learning (fine-tuning)
  message("Starting transfer learning fine-tuning...")
  
  history <- fit(trans_model,
                x = x_train_scale,
                y = to_categorical(y_train),
                epochs = keras_epoch,
                batch_size = keras_batch_size,
                shuffle = TRUE,
                verbose = 0,
                validation_data = list(x_test_scale, to_categorical(y_test)),
                callbacks = list(checkpoint, early_stop))
  
  message("Transfer learning completed")
  
  # Save transfer learning history
  saveRDS(history, "transfer_history.RDS")
  
  # Save final transfer model
  save_model_hdf5(trans_model, "transfer_model_final.hdf5")
  
  # Evaluate transfer learning results
  message("Evaluating transfer learning performance...")
  
  y_train_pred_trans <- predict(trans_model, x_train_scale)
  y_test_pred_trans <- predict(trans_model, x_test_scale)
  y_valid_pred_trans <- predict(trans_model, x_valid_scale)
  
  # Determine optimal cutoff for transfer model
  pred_df_trans <- data.frame(
    Class = as.factor(ifelse(y_train == 1, "Case", "Control")),
    Pred = y_train_pred_trans[, 2]
  )
  
  cutoff_transfer <- cutpointr::cutpointr(pred_df_trans, Pred, Class, pos_class = "Case", metric = youden)
  optimal_cutoff_transfer <- cutoff_transfer$optimal_cutpoint
  
  # Save transfer cutoff plot
  ggplot2::ggsave("transfer_cutoff.png", plot(cutoff_transfer))
  
  # Transfer performance metrics
  temp_results[2, "transfer_training_AUC"] <- cutoff_transfer$AUC
  temp_results[2, "transfer_cutoff"] <- optimal_cutoff_transfer
  
  # Training performance (transfer)
  pred_df_trans$PredClass <- ifelse(pred_df_trans$Pred >= optimal_cutoff_transfer, "Case", "Control")
  pred_df_trans$PredClass <- factor(pred_df_trans$PredClass, levels = c("Control", "Case"))
  cm_train_trans <- confusionMatrix(pred_df_trans$PredClass, pred_df_trans$Class, positive = "Case")
  
  temp_results[2, "transfer_training_Accuracy"] <- cm_train_trans$overall[1]
  temp_results[2, "transfer_training_Sensitivity"] <- cm_train_trans$byClass[1]
  temp_results[2, "transfer_training_Specificity"] <- cm_train_trans$byClass[2]
  temp_results[2, "transfer_training_F1"] <- cm_train_trans$byClass[7]
  
  # Test performance (transfer)
  pred_df_test_trans <- data.frame(Class = as.factor(test$Class), Pred = y_test_pred_trans[, 2])
  pred_df_test_trans$PredClass <- ifelse(pred_df_test_trans$Pred >= optimal_cutoff_transfer, "Case", "Control")
  pred_df_test_trans$PredClass <- factor(pred_df_test_trans$PredClass, levels = c("Control", "Case"))
  cm_test_trans <- confusionMatrix(pred_df_test_trans$PredClass, pred_df_test_trans$Class, positive = "Case")
  
  temp_results[2, "transfer_test_Accuracy"] <- cm_test_trans$overall[1]
  temp_results[2, "transfer_test_Sensitivity"] <- cm_test_trans$byClass[1]
  temp_results[2, "transfer_test_Specificity"] <- cm_test_trans$byClass[2]
  temp_results[2, "transfer_test_F1"] <- cm_test_trans$byClass[7]
  
  # Validation performance (transfer)
  pred_df_valid_trans <- data.frame(Class = as.factor(valid$Class), Pred = y_valid_pred_trans[, 2])
  pred_df_valid_trans$PredClass <- ifelse(pred_df_valid_trans$Pred >= optimal_cutoff_transfer, "Case", "Control")
  pred_df_valid_trans$PredClass <- factor(pred_df_valid_trans$PredClass, levels = c("Control", "Case"))
  cm_valid_trans <- confusionMatrix(pred_df_valid_trans$PredClass, pred_df_valid_trans$Class, positive = "Case")
  
  temp_results[2, "transfer_valid_Accuracy"] <- cm_valid_trans$overall[1]
  temp_results[2, "transfer_valid_Sensitivity"] <- cm_valid_trans$byClass[1]
  temp_results[2, "transfer_valid_Specificity"] <- cm_valid_trans$byClass[2]
  temp_results[2, "transfer_valid_F1"] <- cm_valid_trans$byClass[7]
  
  # Save transfer evaluation results
  saveRDS(cm_train_trans, "transfer_cm_train.RDS")
  saveRDS(cm_test_trans, "transfer_cm_test.RDS")
  saveRDS(cm_valid_trans, "transfer_cm_valid.RDS")
  
  # Save comprehensive results
  data.table::fwrite(temp_results, "transfer_learning_results.csv")
  
  message("Transfer learning neural network completed successfully!")
  
  # Create summary comparison
  initial_perf <- temp_results[1, c("initial_test_Accuracy", "initial_test_Sensitivity", "initial_test_Specificity", "initial_test_F1")]
  transfer_perf <- temp_results[2, c("transfer_test_Accuracy", "transfer_test_Sensitivity", "transfer_test_Specificity", "transfer_test_F1")]
  
  message("Performance Comparison (Test Set):")
  message("  Initial Model - Accuracy: ", round(initial_perf$initial_test_Accuracy, 3), 
          " Sensitivity: ", round(initial_perf$initial_test_Sensitivity, 3),
          " F1: ", round(initial_perf$initial_test_F1, 3))
  message("  Transfer Model - Accuracy: ", round(transfer_perf$transfer_test_Accuracy, 3),
          " Sensitivity: ", round(transfer_perf$transfer_test_Sensitivity, 3), 
          " F1: ", round(transfer_perf$transfer_test_F1, 3))
  
  return(list(
    initial_performance = list(
      train = cm_train_init,
      test = cm_test_init,
      valid = cm_valid_init,
      cutoff = optimal_cutoff_initial,
      auc = cutoff_initial$AUC
    ),
    transfer_performance = list(
      train = cm_train_trans,
      test = cm_test_trans,
      valid = cm_valid_trans,
      cutoff = optimal_cutoff_transfer,
      auc = cutoff_transfer$AUC
    ),
    results_dataframe = temp_results,
    model_paths = list(
      initial_model = model_path_unzipped,
      transfer_model = "transfer_model_final.hdf5"
    ),
    scaling_info = list(
      used_scaling = as.logical(pre_conf[1, "scaled"]),
      new_scaling = new_scaling,
      col_mean = if (exists("col_mean_train")) col_mean_train else NULL,
      col_sd = if (exists("col_sd_train")) col_sd_train else NULL
    ),
    freeze_info = list(
      freeze_from = freeze_from,
      freeze_to = freeze_to,
      layers_frozen = freeze_from > 0
    ),
    status = "completed"
  ))
}

#' Predict Using Transfer Learning Model
#'
#' Makes predictions using a transfer learning model (fine-tuned from pre-trained model).
#' Specifically designed for models that have been adapted via transfer learning.
#'
#' @param model_path Path to ZIP file containing transfer learning model
#' @param new_dataset New dataset for prediction
#' @param new_scaling Logical, whether to perform new scaling on data
#' @param old_train_csv_to_restore_scaling Path to original training data for scaling
#' @param override_cutoff Numeric, custom classification cutoff
#' @param blinded Logical, whether validation is blinded (no true labels)
#'
#' @return List containing predictions, performance metrics, and model information
#'
#' @details
#' This function is specifically designed for models created through transfer learning.
#' It looks for the transfer learning model file (trans_model.hdf5) and uses
#' transfer learning-specific scaling parameters if available.
#'
#' The function handles:
#' \itemize{
#'   \item Loading transfer learning model from ZIP archive
#'   \item Applying appropriate scaling (transfer learning specific or original)
#'   \item Making predictions with optimal cutoff from transfer learning
#'   \item Comprehensive evaluation if labels are provided
#'   \item Support for blinded predictions
#' }
#'
#' @examples
#' \dontrun{
#' # Transfer learning prediction with evaluation
#' results <- predict_transfer_learning(
#'   model_path = "models/transfer_model.zip",
#'   new_dataset = new_data,
#'   blinded = FALSE
#' )
#' print(results$confusion_matrix)
#' 
#' # Blinded prediction
#' predictions <- predict_transfer_learning(
#'   model_path = "models/transfer_model.zip",
#'   new_dataset = new_data,
#'   blinded = TRUE
#' )
#' }
#'
#' @seealso \code{\link{transfer_learning_neural_network}}, \code{\link{predict_deep_learning}}
#' @export
predict_transfer_learning <- function(model_path = "our_models/transfer_model.zip",
                                     new_dataset = NULL,
                                     new_scaling = TRUE,
                                     old_train_csv_to_restore_scaling = NULL,
                                     override_cutoff = NULL,
                                     blinded = TRUE) {
  
  # Check required packages
  required_packages <- c("keras", "data.table", "pROC", "caret", "dplyr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop(paste("Required packages missing:", paste(missing_packages, collapse = ", ")))
  }
  
  # Load required libraries
  library(keras)
  library(data.table)
  library(pROC)
  library(caret)
  library(dplyr)
  
  # Validate inputs
  if (!file.exists(model_path)) {
    stop(paste("Model file does not exist:", model_path))
  }
  
  if (is.null(new_dataset)) {
    stop("New dataset is required for prediction")
  }
  
  message("Loading transfer learning model and making predictions...")
  
  tryCatch({
    # Load transfer learning model from ZIP
    model_files <- utils::unzip(model_path, list = TRUE)
    
    # Look for transfer learning model file
    transfer_model_files <- model_files[grepl("trans_model.hdf5|transfer_model.hdf5|transfer_model_final.hdf5", model_files$Name), "Name"]
    
    if (length(transfer_model_files) == 0) {
      # Fallback to regular model if transfer model not found
      warning("Transfer learning model not found, using regular model")
      transfer_model_files <- model_files[grepl("finalmodel.hdf5", model_files$Name), "Name"]
    }
    
    if (length(transfer_model_files) == 0) {
      stop("No model file found in ZIP archive")
    }
    
    model_path_in_zip <- transfer_model_files[1]
    utils::unzip(model_path, model_path_in_zip, exdir = tempdir())
    model_path_unzipped <- file.path(tempdir(), model_path_in_zip)
    
    if (!file.exists(model_path_unzipped)) {
      stop("Model file does not exist in provided ZIP. Is it correct?")
    }
    
    init_model <- keras::load_model_hdf5(model_path_unzipped)
    message("Transfer learning model loaded successfully")
    print(init_model)
    
    # Load model configuration
    config_files <- model_files[grepl("transfer_learning_results.csv|wyniki.csv", model_files$Name), "Name"]
    if (length(config_files) == 0) {
      stop("Configuration file not found in model package")
    }
    
    utils::unzip(model_path, config_files[1], exdir = tempdir())
    config_path_unzipped <- file.path(tempdir(), config_files[1])
    pre_conf <- data.table::fread(config_path_unzipped)
    
    # Extract network features
    if ("formula" %in% colnames(pre_conf)) {
      network_features <- strsplit(x = as.character(pre_conf[1, "formula"]), split = " + ", fixed = TRUE)[[1]]
    } else {
      # Try to infer features from dataset columns (hsa_ prefix)
      network_features <- colnames(new_dataset)[grepl("^hsa_", colnames(new_dataset))]
      if (length(network_features) == 0) {
        stop("Cannot determine network features from model configuration")
      }
      warning("Formula not found in configuration, using all hsa_ features")
    }
    
    # Check feature availability in new dataset
    missing_features <- network_features[!network_features %in% colnames(new_dataset)]
    if (length(missing_features) > 0) {
      stop(paste0("The new dataset does not contain features: ", paste(missing_features, collapse = ", ")))
    }
    
    # Prepare feature matrix
    new_x <- new_dataset %>% 
      dplyr::select(all_of(network_features)) %>% 
      as.matrix()
    
    # Prepare labels if not blinded
    if (!blinded) {
      if (!"Class" %in% colnames(new_dataset)) {
        stop("Class column required for non-blinded validation")
      }
      
      new_y <- new_dataset %>%
        dplyr::select("Class") %>%
        as.matrix()
      new_y[, 1] <- ifelse(new_y[, 1] == "Case", 1, 0)
      new_dataset$Class <- factor(new_dataset$Class, levels = c("Control", "Case"))
    }
    
    # Handle scaling for transfer learning model
    col_mean_train <- NA
    col_sd_train <- NA
    
    # Check if model was trained with scaling
    scaling_used <- any(c("scaled", "transfer_scaled") %in% colnames(pre_conf)) && 
                   any(grepl("TRUE|T", pre_conf[1, grepl("scaled", colnames(pre_conf))]))
    
    if (scaling_used) {
      message("This transfer learning network was trained with scaling")
      
      # Try to load transfer learning specific scaling parameters
      transfer_mean_files <- model_files[grepl("transfer_col_mean_train.RDS", model_files$Name), "Name"]
      if (length(transfer_mean_files) > 0) {
        utils::unzip(model_path, transfer_mean_files[1], exdir = tempdir())
        mean_path_unzipped <- file.path(tempdir(), transfer_mean_files[1])
        if (file.exists(mean_path_unzipped)) {
          col_mean_train <- readRDS(mean_path_unzipped)
          message("Transfer learning scaling parameters loaded")
        }
      }
      
      transfer_sd_files <- model_files[grepl("transfer_col_sd_train.RDS", model_files$Name), "Name"]
      if (length(transfer_sd_files) > 0) {
        utils::unzip(model_path, transfer_sd_files[1], exdir = tempdir())
        sd_path_unzipped <- file.path(tempdir(), transfer_sd_files[1])
        if (file.exists(sd_path_unzipped)) {
          col_sd_train <- readRDS(sd_path_unzipped)
        }
      }
      
      # Fallback to original scaling parameters if transfer parameters not found
      if (any(is.na(col_mean_train))) {
        mean_files <- model_files[grepl("col_mean_train.RDS", model_files$Name), "Name"]
        if (length(mean_files) > 0) {
          utils::unzip(model_path, mean_files[1], exdir = tempdir())
          mean_path_unzipped <- file.path(tempdir(), mean_files[1])
          if (file.exists(mean_path_unzipped)) {
            col_mean_train <- readRDS(mean_path_unzipped)
            message("Using original model scaling parameters")
          }
        }
      }
      
      if (any(is.na(col_sd_train))) {
        sd_files <- model_files[grepl("col_sd_train.RDS", model_files$Name), "Name"]
        if (length(sd_files) > 0) {
          utils::unzip(model_path, sd_files[1], exdir = tempdir())
          sd_path_unzipped <- file.path(tempdir(), sd_files[1])
          if (file.exists(sd_path_unzipped)) {
            col_sd_train <- readRDS(sd_path_unzipped)
          }
        }
      }
      
      # Apply scaling
      if (new_scaling) {
        message("Applying new scaling")
        new_x <- scale(new_x)
        new_x[is.nan(new_x)] <- 0
        message("New scaling performed")
      } else {
        if (is.null(old_train_csv_to_restore_scaling)) {
          # Use saved scaling parameters
          if (!any(is.na(col_mean_train)) && !any(is.na(col_sd_train))) {
            new_x <- scale(new_x, center = col_mean_train, scale = col_sd_train)
            new_x[is.nan(new_x)] <- 0
            message("Applied saved scaling parameters")
          } else {
            warning("Scaling parameters not available, proceeding without scaling")
          }
        } else {
          # Restore scaling from provided file
          if (!file.exists(old_train_csv_to_restore_scaling)) {
            stop("A file to restore scaling from does not exist")
          }
          
          message("Restoring scaling from provided file")
          temp_train <- data.table::fread(old_train_csv_to_restore_scaling)
          temp_train_x <- dplyr::select(temp_train, all_of(network_features))
          temp_train_x <- scale(temp_train_x)
          col_mean_train <- attr(temp_train_x, "scaled:center")
          col_sd_train <- attr(temp_train_x, "scaled:scale")
          
          new_x <- scale(new_x, center = col_mean_train, scale = col_sd_train)
          new_x[is.nan(new_x)] <- 0
          message("Scaling restored from provided file")
        }
      }
    }
    
    # Get cutoff for prediction
    # Try transfer learning specific cutoff first
    cutoff <- NA
    if ("transfer_cutoff" %in% colnames(pre_conf)) {
      cutoff <- as.numeric(pre_conf[nrow(pre_conf), "transfer_cutoff"])  # Use last row for transfer learning results
    } else if ("cutoff" %in% colnames(pre_conf)) {
      cutoff <- as.numeric(pre_conf[1, "cutoff"])
    }
    
    if (!is.null(override_cutoff)) {
      cutoff <- as.numeric(override_cutoff)
      message("Using override cutoff: ", cutoff)
    }
    
    if (is.na(cutoff)) {
      cutoff <- 0.5
      warning("Using default cutoff of 0.5")
    }
    
    message("Using cutoff: ", cutoff)
    
    # Make predictions
    predictions <- predict(init_model, new_x)
    
    if (blinded) {
      # Blinded prediction (no true labels)
      pred <- data.frame(Pred = predictions[, 2])
      pred$Prediction <- ifelse(pred$Pred >= cutoff, "Case", "Control")
      confusion_matrix <- NA
      roc_result <- NA
      roc_auc <- NA
      
      message("Blinded transfer learning predictions completed")
    } else {
      # Full evaluation with true labels
      pred <- data.frame(
        Class = factor(ifelse(new_y == 1, "Case", "Control"), levels = c("Control", "Case")),
        Pred = predictions[, 2]
      )
      pred$Prediction <- ifelse(pred$Pred >= cutoff, "Case", "Control")
      pred$Correctness <- ifelse(pred$Prediction == pred$Class, "Correct", "Incorrect")
      
      # Calculate confusion matrix
      confusion_matrix <- caret::confusionMatrix(
        as.factor(pred$Prediction),
        as.factor(pred$Class),
        positive = "Case"
      )
      
      # Calculate ROC and AUC
      roc_result <- pROC::roc(pred$Class ~ pred$Pred, quiet = TRUE)
      roc_auc <- pROC::ci.auc(roc_result)
      
      message("Transfer learning performance evaluation completed:")
      message("  Accuracy: ", round(confusion_matrix$overall["Accuracy"], 3))
      message("  Sensitivity: ", round(confusion_matrix$byClass["Sensitivity"], 3))
      message("  Specificity: ", round(confusion_matrix$byClass["Specificity"], 3))
      message("  AUC: ", round(as.numeric(roc_result$auc), 3))
    }
    
    # Compile comprehensive results
    final_return <- list(
      predictions = pred,
      network_config = pre_conf,
      new_dataset = new_dataset,
      new_dataset_x = new_x,
      network_features = network_features,
      cutoff = cutoff,
      col_mean_train = col_mean_train,
      col_sd_train = col_sd_train,
      confusion_matrix = confusion_matrix,
      roc = roc_result,
      roc_auc = roc_auc,
      model = init_model,
      model_type = "transfer_learning"
    )
    
    message("Transfer learning prediction completed successfully")
    return(final_return)
    
  }, error = function(e) {
    stop(paste("Error in transfer learning prediction:", e$message))
  })
}

# Backward compatibility aliases
#' @rdname create_keras_model
#' @export
OmicSelector_keras_create_model <- create_keras_model

#' @rdname train_deep_learning
#' @export
OmicSelector_deep_learning <- train_deep_learning

#' @rdname transfer_learning_neural_network
#' @export
OmicSelector_transfer_learning_neural_network <- transfer_learning_neural_network

#' @rdname predict_deep_learning
#' @export
OmicSelector_deep_learning_predict <- predict_deep_learning

#' @rdname predict_transfer_learning
#' @export
OmicSelector_deep_learning_predict_transfered <- predict_transfer_learning

#' Load Deep Learning Extension
#'
#' Loads the complete deep learning extension for full functionality.
#' This provides access to the original comprehensive implementation.
#'
#' @export
load_deeplearning_extension <- function() {
  extension_path <- system.file("extensions", "deeplearning.R", package = "OmicSelector")
  
  if (extension_path == "") {
    extension_path <- file.path(dirname(dirname(getwd())), "extensions", "deeplearning.R")
  }
  
  if (file.exists(extension_path)) {
    source(extension_path)
    message("Deep learning extension loaded successfully")
  } else {
    stop("Deep learning extension file not found")
  }
}
