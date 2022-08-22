# -*- coding: utf-8 -*-
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(caret))
suppressMessages(library(epiDisplay))
suppressMessages(library(pROC))
suppressMessages(library(ggplot2))
suppressMessages(library(DMwR))
suppressMessages(library(ROSE))
suppressMessages(library(gridExtra))
suppressMessages(library(gplots))
suppressMessages(library(devtools))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
library(OmicSelector)
suppressMessages(library(funModeling))
message("OmicSelector: DeepLearning extension loaded.")

if(!dir.exists(paste0("models"))) { dir.create(paste0("models")) }
if(!dir.exists(paste0("temp"))) { dir.create(paste0("temp")) }

OmicSelector_keras_create_model <- function(i, hyperparameters, how_many_features = ncol(x_train_scale)) {
  # tempmodel <- keras_model_sequential() %>%
  #   { if(hyperparameters[i,10]==T) { layer_dense(. , units = hyperparameters[i,1], kernel_regularizer = regularizer_l2(l = 0.001),
  #                                                activation = hyperparameters[i,4], input_shape = c(ncol(x_train_scale))) } else
  #                                  { layer_dense(. , units = hyperparameters[i,1], activation = hyperparameters[i,4],
  #                                                                input_shape = c(ncol(x_train_scale))) } } %>%
  #   { if(hyperparameters[i,7]>0) { layer_dropout(. , rate = hyperparameters[i,7]) } else { . } } %>%
  #   { if(hyperparameters[i,2]>0) {
  #   if(hyperparameters[i,11]==T) { layer_dense(. , units = hyperparameters[i,2], activation = hyperparameters[i,5],
  #               kernel_regularizer = regularizer_l2(l = 0.001)) } else {
  #                 layer_dense(units = hyperparameters[i,2], activation = hyperparameters[i,5]) } } }  %>%
  #   { if(hyperparameters[i,8]>0) { layer_dropout(rate = hyperparameters[i,8]) } else { . } } %>%
  #   { if(hyperparameters[i,3]>0) {
  #   if(hyperparameters[i,12]==T) { layer_dense(units = hyperparameters[i,3], activation = hyperparameters[i,6],
  #                                                kernel_regularizer = regularizer_l2(l = 0.001)) } else
  #                                   {layer_dense(units = hyperparameters[i,3], activation = hyperparameters[i,6])} } else { . } } %>%
  #   { if(hyperparameters[i,9]>0) { layer_dropout(rate = hyperparameters[i,9]) } else { . } } %>%
  #   layer_dense(units = 1, activation = 'sigmoid')
  library(keras)
  hyperparameters = as.data.frame(hyperparameters)
  tempmodel <- keras_model_sequential()
  if(hyperparameters[i,10]==T) { layer_dense(tempmodel , units = as.character(hyperparameters[i,1]), kernel_regularizer = regularizer_l1(l = 0.001),
                                             activation = as.character(hyperparameters[i,4]), input_shape = c(how_many_features)) } else
                                             { layer_dense(tempmodel , units = as.character(hyperparameters[i,1]), activation = as.character(hyperparameters[i,4]),
                                                           input_shape = c(how_many_features)) }
  if(hyperparameters[i,7]>0) { layer_dropout(tempmodel , rate = as.numeric(hyperparameters[i,7])) }
  if(hyperparameters[i,2]>0) {
    if(hyperparameters[i,11]==T) { layer_dense(tempmodel , units = as.character(hyperparameters[i,2]), activation = as.character(hyperparameters[i,5]),
                                               kernel_regularizer = regularizer_l1(l = 0.001)) } else
                                               {layer_dense(tempmodel, units = as.character(hyperparameters[i,2]), activation = as.character(hyperparameters[i,5])) } }

  if(hyperparameters[i,2]>0 & hyperparameters[i,8]>0) { layer_dropout(tempmodel, rate = as.numeric(hyperparameters[i,8])) }
  if(hyperparameters[i,3]>0) {
    if(hyperparameters[i,12]==T) { layer_dense(tempmodel, units = as.character(hyperparameters[i,3]), activation = as.character(hyperparameters[i,6]),
                                               kernel_regularizer = regularizer_l1(l = 0.001)) } else
                                               { layer_dense(tempmodel, units = as.character(hyperparameters[i,3]), activation = as.character(hyperparameters[i,6]))} }
  if(hyperparameters[i,3]>0 & hyperparameters[i,9]>0) { layer_dropout(rate = as.numeric(hyperparameters[i,9])) }
  layer_dense(tempmodel, units = 2, activation = 'softmax')

  print(tempmodel)


  dnn_class_model = keras::compile(tempmodel, optimizer = as.character(hyperparameters[i,13]),
                                   loss = 'binary_crossentropy',
                                   metrics = 'accuracy')

}

# autoencoder ma softmax na deep feature
# jesli zmienna autoencoder >0 budujemy autoencoder 3-wartwowy, jesli <0 budujemy autoencoder 3-warstwowy z regularyzacja (sparse autoencoder)
OmicSelector_deep_learning = function(selected_miRNAs = ".", wd = getwd(),
                            SMOTE = F, keras_batch_size = 64, clean_temp_files = T,
                            save_threshold_trainacc = 0.85, save_threshold_testacc = 0.8, keras_epochae = 5000,
                            keras_epoch = 2000, keras_patience = 50,
                            hyperparameters = expand.grid(layer1 = seq(3,11, by = 2), layer2 = c(0,seq(3,11, by = 2)), layer3 = c(0,seq(3,11, by = 2)),
                                                          activation_function_layer1 = c("relu","sigmoid"), activation_function_layer2 = c("relu","sigmoid"), activation_function_layer3 = c("relu","sigmoid"),
                                                          dropout_layer1 = c(0, 0.1), dropout_layer2 = c(0, 0.1), dropout_layer3 = c(0),
                                                          layer1_regularizer = c(T,F), layer2_regularizer = c(T,F), layer3_regularizer = c(T,F),
                                                          optimizer = c("adam","rmsprop","sgd"), autoencoder = c(0,7,-7), balanced = SMOTE, formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3], scaled = c(T,F),
                                                          stringsAsFactors = F), add_features_to_predictions = F,
                            keras_threads = ceiling(parallel::detectCores()/2), start = 1, end = nrow(hyperparameters), output_file = "deeplearning_results.csv", save_all_vars = F, automatic_weight = F)
{
  # library(OmicSelector)
  # OmicSelector_load_extension("deeplearning")
  codename = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(output_file))
  #options(warn=-1)
  oldwd = getwd()
  setwd = setwd(wd)
  set.seed(1)
if(dir.exists("/OmicSelector")) {
  if(!dir.exists(paste0(getwd(), "/temp-deeplearning"))) { dir.create(paste0(getwd(), "/temp-deeplearning")) }
    temp_dir = paste0(getwd(), "/temp-deeplearning")
   } else {
     temp_dir = tempdir()
   }
  if(!dir.exists("temp")) { dir.create("temp") }
  if(!dir.exists("models")) { dir.create("models") }
  options(bitmapType = 'cairo', device = 'png')
  library(plyr)
  library(dplyr)
  library(keras)
  library(foreach)
  library(doParallel)
  #library(doSNOW)
  library(data.table)
  fwrite(hyperparameters, paste0("hyperparameters_",output_file))

  dane = OmicSelector_load_datamix(wd = wd, replace_smote = F, remove_zero_var = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
  
  #train = data.table::fread("mixed_train.csv") %>% dplyr::select(starts_with("hsa"), Class)
  colnames(train)
  
  if (SMOTE == T) { train = train_smoted }
  message("Checkpoint passed: load lib and data")

  #cores=detectCores()
  cat(paste0("\nTemp dir: ", temp_dir, "\n"))
  cat("\nStarting preparing cluster..\n")
  #cl <- makePSOCKcluster(keras_threads) #not to overload your computer
  clusterlogfile = paste0("temp/", ceiling(as.numeric(Sys.time())), "deeplearning_cluster.log")
  cl = makeCluster(keras_threads, outfile=clusterlogfile)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  cat("\nCluster prepared..\n")
  #message("Checkpoint passed: cluster prepared")




  #args= names(mget(ls()))
  #export = export[!export %in% args]

  # tu musi isc iteracja
  cat(paste0("\nStarting parallel loop.. There are: ", end-start+1, " hyperparameter sets to be checked.\n"))
  final <- foreach(i=as.numeric(start):as.numeric(end), .combine=rbind, .verbose=F, .inorder=F, .errorhandling = 'remove', .export = ls()
                   #,.packages = loadedNamespaces()
  ) %dopar% {

    Sys.setenv(TF_FORCE_GPU_ALLOW_GROWTH = 'true')
    library(OmicSelector)
    OmicSelector_load_extension("deeplearning")

    if(!dir.exists(paste0(temp_dir,"/models"))) { dir.create(paste0(temp_dir,"/models")) }
    if(!dir.exists(paste0(temp_dir,"/temp"))) { dir.create(paste0(temp_dir,"/temp")) }
    start_time <- Sys.time()
    library(keras)
    library(ggplot2)
    library(dplyr)
    library(data.table)
    set.seed(1)

    # library(tensorflow)
    # gpu = tf$test$is_gpu_available()
    # if(gpu) {
    # gpux <- tf$config$experimental$get_visible_devices('GPU')[[1]]
    # tf$config$experimental$set_memory_growth(device = gpux, enable = TRUE)
    # py_run_string("gpus = tf.config.experimental.list_physical_devices('GPU')")
    # py_run_string("tf.config.experimental.set_virtual_device_configuration(gpus[0],[tf.config.experimental.VirtualDeviceConfiguration(memory_limit=2048)])")

    # }

    cat("\nStarting hyperparameters..\n")
    print(hyperparameters[i,])
    message(paste0("OmicSelector: Starting training network no ", i, "."))
    message(paste0(hyperparameters[i,], collapse = ", "))

    options(bitmapType = 'cairo', device = 'png')
    model_id = paste0(format(i, scientific = FALSE), "-", ceiling(as.numeric(Sys.time())))
    if(SMOTE == T) { model_id = paste0(format(i, scientific = FALSE), "-SMOTE-", ceiling(as.numeric(Sys.time()))) }
    tempwyniki = data.frame(model_id=model_id)
    tempwyniki[1, "model_id"] = model_id
    #message("Checkpoint passed: chunk 1")


    if(!dir.exists(paste0(temp_dir,"/models"))) { dir.create(paste0(temp_dir,"/models"))}
    if(!dir.exists(paste0(temp_dir,"/models/keras",model_id))) { dir.create(paste0(temp_dir,"/models/keras",model_id))}
    cat(paste0("\nTraining model: ",temp_dir,"/models/keras",model_id,"\n"))
    #message("Checkpoint passed: chunk 2")
    #pdf(paste0(temp_dir,"/models/keras",model_id,"/plots.pdf"), paper="a4")

    con <- file(paste0(temp_dir,"/models/keras",model_id,"/training.log"))
    sink(con, append=TRUE, split =TRUE)
    sink(con, append=TRUE, type="message")

    early_stop <- callback_early_stopping(monitor = "val_loss", mode="min", patience = keras_patience)
    cp_callback <- callback_model_checkpoint(
      filepath =  paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"),
      save_best_only = TRUE, period = 10, monitor = "val_loss",
      verbose = 0
    )
    ae_cp_callback <- callback_model_checkpoint(
      filepath =  paste0(temp_dir,"/models/keras",model_id,"/autoencoderweights.hdf5"),
      save_best_only = TRUE, save_weights_only = T, period = 10, monitor = "val_loss",
      verbose = 0
    )

    x_train <- train %>%
      { if (selected_miRNAs[1] != ".") { dplyr::select(train, selected_miRNAs) } else { dplyr::select(train, starts_with("hsa")) } } %>%
      as.matrix()
    y_train <- train %>%
      dplyr::select("Class") %>%
      as.matrix()
    y_train[,1] = ifelse(y_train[,1] == "Case",1,0)


    x_test <- test %>%
      { if (selected_miRNAs[1] != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>%
      as.matrix()
    y_test <- test %>%
      dplyr::select("Class") %>%
      as.matrix()
    y_test[,1] = ifelse(y_test[,1] == "Case",1,0)

    x_valid <- valid %>%
      { if (selected_miRNAs[1] != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>%
      as.matrix()
    y_valid <- valid %>%
      dplyr::select("Class") %>%
      as.matrix()
    y_valid[,1] = ifelse(y_valid[,1] == "Case",1,0)
    #message("Checkpoint passed: chunk 3")

    if(automatic_weight) {
      counter=funModeling::freq(to_categorical(y_train), plot=F) %>% select(var, frequency)
      majority=max(counter$frequency)
      counter$weight=ceil(majority/counter$frequency)
      l_weights=setNames(as.list(counter$weight), counter$var)
      message(l_weights)
    } else {
      l_weights = NULL
    }

    if(hyperparameters[i, 17] == T) {
      x_train_scale = x_train %>% scale()

      col_mean_train <- attr(x_train_scale, "scaled:center")
      saveRDS(col_mean_train, paste0(temp_dir,"/models/keras",model_id,"/col_mean_train.RDS"))
      col_sd_train <- attr(x_train_scale, "scaled:scale")
      saveRDS(col_sd_train, paste0(temp_dir,"/models/keras",model_id,"/col_sd_train.RDS"))

      x_test_scale <- x_test %>%
        scale(center = col_mean_train,
              scale = col_sd_train)

      x_valid_scale <- x_valid %>%
        scale(center = col_mean_train,
              scale = col_sd_train)
    } else {
      x_train_scale = x_train
      x_test_scale <- x_test
      x_valid_scale <- x_valid
    }

    input_layer <- layer_input(shape = c(ncol(x_train_scale)))
    #message("Checkpoint passed: chunk 4")

    #psych::describe(x_test_scale)

    # czy autoenkoder?
    if(hyperparameters[i,14] != 0) {


      n1 <- hyperparameters[i,14]
      n3 <- ncol(x_train_scale)

      n2 = 7
      if(ceiling(n3/2) > 7) { n2 = ceiling(n3/2) }

      if (hyperparameters[i,14]>0) {
        encoder <-
          input_layer %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i,6]))    %>%
          layer_dense(units = n1, activation = "softmax")  # dimensions of final encoding layer

        decoder <- encoder %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i,6])) %>%
          layer_dense(units = n3, as.character(hyperparameters[i,6]))  # dimension of original variable
      }
      else {
        n1 = -n1 # korekta dla ujemnej wartosci w hiperparametrach
        encoder <-
          input_layer %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i,6]), kernel_regularizer = regularizer_l1(l = 0.01))    %>%
          layer_dense(units = n1, activation = "softmax", kernel_regularizer = regularizer_l1(l = 0.01))  # dimensions of final encoding layer

        decoder <- encoder %>%
          layer_dense(units = n2, activation = as.character(hyperparameters[i,6]), kernel_regularizer = regularizer_l1(l = 0.01))    %>%
          layer_dense(units = n3, as.character(hyperparameters[i,6]))  # dimension of original variable
      }
      ae_model <- keras_model(inputs = input_layer, outputs = decoder)
      ae_model

      ae_model %>%
        keras::compile(loss = "mean_absolute_error",
                       optimizer = as.character(hyperparameters[i,13]),
                       metrics = c("mean_squared_error"))

      summary(ae_model)
      #message("Checkpoint passed: chunk 5")



      #ae_tempmodelfile = tempfile()
      ae_history <- fit(ae_model, x = x_train_scale,
                        y = x_train_scale,
                        epochs = keras_epochae, batch_size = keras_batch_size,
                        shuffle = T, verbose = 0,
                        view_metrics = FALSE,
                        validation_data = list(x_test_scale, x_test_scale),
                        callbacks = list(ae_cp_callback, early_stop))
      #file.copy(ae_tempmodelfile, paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"), overwrite = T)

      #ae_history_df <- as.data.frame(ae_history)
      #fwrite(ae_history_df, paste0(temp_dir,"/models/keras",model_id,"/autoencoder_history_df.csv.gz"))
      saveRDS(ae_history, paste0(temp_dir,"/models/keras",model_id,"/ae_history.RDS"))

      compare_cx <- data.frame(
        train_loss = ae_history$metrics$loss,
        test_loss = ae_history$metrics$val_loss
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)

      plot1 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("epoch") +
        ylab("loss") +
        theme_bw()

      #message("Checkpoint passed: chunk 6")
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training_autoencoder.png"), grid.arrange(plot1, nrow =1, top = "Training of autoencoder"))

      cat("\n- Reloading autoencoder to get weights...\n")
      encoder_model <- keras_model(inputs = input_layer, outputs = encoder)
      encoder_model %>% load_model_weights_hdf5(paste0(temp_dir,"/models/keras",model_id,"/autoencoderweights.hdf5"),
                                                skip_mismatch = T,
                                                by_name = T)

      cat("\n- Autoencoder saving...\n")
      save_model_hdf5(encoder_model, paste0(temp_dir,"/models/keras",model_id,"/autoencoder.hdf5"))


      cat("\n- Creating deep features...\n")
      #message("Checkpoint passed: chunk 7")

      ae_x_train_scale <- encoder_model %>%
        predict(x_train_scale) %>%
        as.matrix()
      fwrite(ae_x_train_scale, paste0(temp_dir,"/models/keras",model_id,"/deepfeatures_train.csv"))

      ae_x_test_scale <- encoder_model %>%
        predict(x_test_scale) %>%
        as.matrix()
      fwrite(ae_x_test_scale, paste0(temp_dir,"/models/keras",model_id,"/deepfeatures_test.csv"))

      ae_x_valid_scale <- encoder_model %>%
        predict(x_valid_scale) %>%
        as.matrix()
      fwrite(ae_x_valid_scale, paste0(temp_dir,"/models/keras",model_id,"/deepfeatures_valid.csv"))

      # # podmiana eby nie edytowa kodu
      x_train_scale = as.matrix(ae_x_train_scale)
      x_test_scale = as.matrix(ae_x_test_scale)
      x_valid_scale = as.matrix(ae_x_valid_scale)

      cat("\n- Training model based on deep features...\n")
      dnn_class_model <- OmicSelector_keras_create_model(i, hyperparameters = hyperparameters, how_many_features = ncol(x_train_scale))
      #message("Checkpoint passed: chunk 8")
      message("Starting training...")
      #tempmodelfile = tempfile()
      history <- fit(dnn_class_model, x = x_train_scale,
                     y = to_categorical(y_train),
                     epochs = keras_epoch,
                     validation_data = list(x_test_scale, to_categorical(y_test)),
                     callbacks = list(cp_callback, early_stop),
                     verbose = 0,
                     view_metrics = FALSE,
                     batch_size = keras_batch_size, shuffle = T, class_weight = l_weights)
      message(history)
      print(history)
      #message("Checkpoint passed: chunk 8")
      #plot(history, col="black")
      #history_df <- as.data.frame(history)
      # fwrite(history_df, paste0(temp_dir,"/models/keras",model_id,"/history_df.csv.gz"))
      saveRDS(history, paste0(temp_dir,"/models/keras",model_id,"/history.RDS"))
      #message("Checkpoint passed: chunk 9")

      cat("\n- Saving history and plots...\n")
      compare_cx <- data.frame(
        train_loss = history$metrics$loss,
        test_loss = history$metrics$val_loss
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)

      plot1 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("Epoch") +
        ylab("Loss") +
        theme_bw()
      #message("Checkpoint passed: chunk 10")

      compare_cx <- data.frame(
        train_accuracy = history$metrics$accuracy,
        test_accuracy = history$metrics$val_accuracy
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)

      plot2 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("Epoch") +
        ylab("Accuracy") +
        theme_bw()
      #message("Checkpoint passed: chunk 12")

      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training.png"), grid.arrange(plot1, plot2, nrow =2, top = "Training of final neural network"))

      # ewaluacja
      cat("\n- Saving final model...\n")
      dnn_class_model = load_model_hdf5(paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"))
      y_train_pred <- predict(object = dnn_class_model, x = x_train_scale)
      y_test_pred <- predict(object = dnn_class_model, x = x_test_scale)
      y_valid_pred <- predict(object = dnn_class_model, x = x_valid_scale)
      #message("Checkpoint passed: chunk 13")


      # wybranie odciecia
      pred = data.frame(`Class` = train$Class, `Pred` = y_train_pred)
      library(cutpointr)
      cutoff = cutpointr(pred, Pred.2, Class, pos_class = "Case", metric = youden)
      summary(cutoff)
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/cutoff.png"), plot(cutoff))
      wybrany_cutoff = cutoff$optimal_cutpoint
      #wyniki[i, "training_AUC"] = cutoff$AUC
      tempwyniki[1, "training_AUC"] = cutoff$AUC
      #wyniki[i, "cutoff"] = wybrany_cutoff
      tempwyniki[1, "cutoff"] = wybrany_cutoff
      cat(paste0("\n\n---- TRAINING AUC: ",cutoff$AUC," ----\n\n"))
      cat(paste0("\n\n---- OPTIMAL CUTOFF: ",wybrany_cutoff," ----\n\n"))
      cat(paste0("\n\n---- TRAINING PERFORMANCE ----\n\n"))
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
      cm_train = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
      print(cm_train)
      #message("Checkpoint passed: chunk 14")

      t1_roc = pROC::roc(Class ~ as.numeric(Pred.2), data=pred)
      tempwyniki[1, "training_AUC2"] = t1_roc$auc
      tempwyniki[1, "training_AUC_lower95CI"] = as.character(ci(t1_roc))[1]
      tempwyniki[1, "training_AUC_upper95CI"] = as.character(ci(t1_roc))[3]
      saveRDS(t1_roc, paste0(temp_dir,"/models/keras",model_id,"/training_ROC.RDS"))
      #message("Checkpoint passed: chunk 15")

      tempwyniki[1, "training_Accuracy"] = cm_train$overall[1]
      tempwyniki[1, "training_Sensitivity"] = cm_train$byClass[1]
      tempwyniki[1, "training_Specificity"] = cm_train$byClass[2]
      tempwyniki[1, "training_PPV"] = cm_train$byClass[3]
      tempwyniki[1, "training_NPV"] = cm_train$byClass[4]
      tempwyniki[1, "training_F1"] = cm_train$byClass[7]
      saveRDS(cm_train, paste0(temp_dir,"/models/keras",model_id,"/cm_train.RDS"))
      #message("Checkpoint passed: chunk 16")

      cat(paste0("\n\n---- TESTING PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = test$Class, `Pred` = y_test_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
      cm_test = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
      print(cm_test)
      tempwyniki[1, "test_Accuracy"] = cm_test$overall[1]
      tempwyniki[1, "test_Sensitivity"] = cm_test$byClass[1]
      tempwyniki[1, "test_Specificity"] = cm_test$byClass[2]
      tempwyniki[1, "test_PPV"] = cm_test$byClass[3]
      tempwyniki[1, "test_NPV"] = cm_test$byClass[4]
      tempwyniki[1, "test_F1"] = cm_test$byClass[7]
      saveRDS(cm_test, paste0(temp_dir,"/models/keras",model_id,"/cm_test.RDS"))
      #message("Checkpoint passed: chunk 17")

      cat(paste0("\n\n---- VALIDATION PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = valid$Class, `Pred` = y_valid_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
      cm_valid = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
      print(cm_valid)
      tempwyniki[1, "valid_Accuracy"] = cm_test$overall[1]
      tempwyniki[1, "valid_Sensitivity"] = cm_test$byClass[1]
      tempwyniki[1, "valid_Specificity"] = cm_test$byClass[2]
      tempwyniki[1, "valid_PPV"] = cm_test$byClass[3]
      tempwyniki[1, "valid_NPV"] = cm_test$byClass[4]
      tempwyniki[1, "valid_F1"] = cm_test$byClass[7]
      saveRDS(cm_valid, paste0(temp_dir,"/models/keras",model_id,"/cm_valid.RDS"))
      #message("Checkpoint passed: chunk 18")

      if(add_features_to_predictions) {
        mix = rbind(train,test,valid)
        mixx = rbind(x_train_scale, x_test_scale, x_valid_scale)
        y_mixx_pred <- predict(object = dnn_class_model, x = mixx)

        mix$Podzial = c(rep("Training",nrow(train)),rep("Test",nrow(test)),rep("Validation",nrow(valid)))
        mix$Pred = y_mixx_pred[,2]
        mix$PredClass = ifelse(mix$Pred >= wybrany_cutoff, "Case", "Control")
        fwrite(mix, paste0(temp_dir,"/models/keras",model_id,"/data_predictions.csv.gz")) } else {
          mix = rbind(train,test,valid)
          mixx = rbind(x_train_scale, x_test_scale, x_valid_scale)
          y_mixx_pred <- predict(object = dnn_class_model, x = mixx)

          mix2 = data.frame(
            `Podzial` = c(rep("Training",nrow(train)),rep("Test",nrow(test)),rep("Validation",nrow(valid))),
            `Pred` = y_mixx_pred[,2],
            `PredClass` = ifelse(y_mixx_pred[,2] >= wybrany_cutoff, "Case", "Control")
          )

          fwrite(mix2, paste0(temp_dir,"/models/keras",model_id,"/data_predictions.csv.gz"))

        }
      fwrite(cbind(hyperparameters[i,], tempwyniki), paste0(temp_dir,"/models/keras",model_id,"/wyniki.csv"))
      #message("Checkpoint passed: chunk 19")

      wagi = keras::get_weights(dnn_class_model)
      saveRDS(wagi, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.RDS"))
      save_model_weights_hdf5(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.hdf5"))
      saveRDS(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel.RDS"))
      #message("Checkpoint passed: chunk 20")

      # czy jest sens zapisywac?
      #sink()
      #sink(type="message")
      message(paste0("\n ",model_id, ": ", tempwyniki[1, "training_Accuracy"], " / ", tempwyniki[1, "test_Accuracy"], " ==> ", tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc))
      cat(paste0("\n ",model_id, ": ", tempwyniki[1, "training_Accuracy"], " / ", tempwyniki[1, "test_Accuracy"], " ==> ", tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc))
      if(tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc) {
        # zapisywanie modelu do waciwego katalogu
        if (save_all_vars) { save(list = ls(all=TRUE), file = paste0(temp_dir,"/models/keras",model_id,"/all.Rdata.gz"), compress = "gzip", compression_level = 9) }
        if(!dir.exists(paste0("models/",codename,"/"))) { dir.create(paste0("models/",codename,"/")) }
        if(dir.exists("/OmicSelector")) {
          system(paste0("/bin/bash -c 'screen -dmS save_network zip -9 -r ", paste0(oldwd,"/models/",codename,"/",codename, "_", model_id,".zip") ," ", paste0(temp_dir,"/models/keras",model_id), "/'"))
          # Debug:
          # OmicSelector_log(paste0("Running: /bin/bash -c 'screen -dmS save_network zip -9 -r ", paste0(oldwd,"/models/",codename,"/",codename, "_", model_id,".zip") ," ", paste0(temp_dir,"/models/keras",model_id), "/'"), "task.log")
        } else {
          zip(paste0(oldwd,"models/",codename,"/",codename, "_", model_id,".zip"),list.files(paste0(temp_dir,"/models/keras",model_id), full.names = T, recursive = T, include.dirs = T)) }
        # file.copy(list.files(paste0(temp_dir,"/models/keras",model_id), pattern = "_wyniki.csv$", full.names = T, recursive = T, include.dirs = T),paste0("temp/",codename,"_",model_id,"_deeplearningresults.csv"))
        if (!dir.exists(paste0("temp/",codename,"/"))) { dir.create(paste0("temp/",codename,"/")) }
        file.copy(list.files(paste0(temp_dir,"/models/keras",model_id), pattern = "_wyniki.csv$", full.names = T, recursive = T, include.dirs = T),paste0("temp/",codename,"/",model_id,"_deeplearningresults.csv"))
        #message("Checkpoint passed: chunk 21")
        #dev.off()
      }
    } else {
      x_train <- train %>%
        { if (selected_miRNAs[1] != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>%
        as.matrix()
      y_train <- train %>%
        dplyr::select("Class") %>%
        as.matrix()
      y_train[,1] = ifelse(y_train[,1] == "Case",1,0)
      #message("Checkpoint passed: chunk 22")


      x_test <- test %>%
        { if (selected_miRNAs[1] != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>%
        as.matrix()
      y_test <- test %>%
        dplyr::select("Class") %>%
        as.matrix()
      y_test[,1] = ifelse(y_test[,1] == "Case",1,0)

      x_valid <- valid %>%
        { if (selected_miRNAs[1] != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>%
        as.matrix()
      y_valid <- valid %>%
        dplyr::select("Class") %>%
        as.matrix()
      y_valid[,1] = ifelse(y_valid[,1] == "Case",1,0)
      #message("Checkpoint passed: chunk 23")

      if(hyperparameters[i, 17] == T) {
        x_train_scale = x_train %>% scale()

        col_mean_train <- attr(x_train_scale, "scaled:center")
        col_sd_train <- attr(x_train_scale, "scaled:scale")

        x_test_scale <- x_test %>%
          scale(center = col_mean_train,
                scale = col_sd_train)

        x_valid_scale <- x_valid %>%
          scale(center = col_mean_train,
                scale = col_sd_train)
      }
      else {
        x_train_scale = x_train
        x_test_scale <- x_test
        x_valid_scale <- x_valid
      }

      #message("Checkpoint passed: chunk 24")
      dnn_class_model <- OmicSelector_keras_create_model(i, hyperparameters = hyperparameters, how_many_features = ncol(x_train_scale))

      message("Starting training...")
      history <-  fit(dnn_class_model, x = x_train_scale,
                      y = to_categorical(y_train),
                      epochs = keras_epoch,
                      validation_data = list(x_test_scale, to_categorical(y_test)),
                      callbacks = list(
                        cp_callback,
                        #callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.1),
                        #callback_model_checkpoint(paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5")),
                        early_stop),
                      verbose = 0,
                      view_metrics = FALSE,
                      batch_size = keras_batch_size, shuffle = T, class_weight = l_weights)
      print(history)
      #message("Checkpoint passed: chunk 25")
      message(history)
      #plot(history, col="black")
      saveRDS(history, paste0(temp_dir,"/models/keras",model_id,"/history.RDS"))
      #message("Checkpoint passed: chunk 26")
      #fwrite(history_df, paste0(temp_dir,"/models/keras",model_id,"/history_df.csv.gz"))


      # pdf(paste0(temp_dir,"/models/keras",model_id,"/plots.pdf"))
      compare_cx <- data.frame(
        train_loss = history$metrics$loss,
        test_loss = history$metrics$val_loss
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)

      plot1 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("epoch") +
        ylab("loss") +
        theme_bw()
      #message("Checkpoint passed: chunk 27")

      compare_cx <- data.frame(
        train_accuracy = history$metrics$accuracy,
        test_accuracy = history$metrics$val_accuracy
      ) %>%
        tibble::rownames_to_column() %>%
        mutate(rowname = as.integer(rowname)) %>%
        tidyr::gather(key = "type", value = "value", -rowname)

      plot2 = ggplot(compare_cx, aes(x = rowname, y = value, color = type)) +
        geom_line() +
        xlab("epoch") +
        ylab("loss") +
        theme_bw()
      #message("Checkpoint passed: chunk 28")

      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training.png"), grid.arrange(plot1, plot2, nrow =2, top = "Training of final neural network"))

      # ewaluacja
      dnn_class_model = load_model_hdf5(paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5"))
      y_train_pred <- predict(object = dnn_class_model, x = x_train_scale)
      y_test_pred <- predict(object = dnn_class_model, x = x_test_scale)
      y_valid_pred <- predict(object = dnn_class_model, x = x_valid_scale)
      #message("Checkpoint passed: chunk 29")


      # wybranie odciecia
      pred = data.frame(`Class` = train$Class, `Pred` = y_train_pred)
      library(cutpointr)
      cutoff = cutpointr(pred, Pred.2, Class, pos_class = "Case", metric = youden)
      print(summary(cutoff))
      ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/cutoff.png"), plot(cutoff))
      wybrany_cutoff = cutoff$optimal_cutpoint
      #wyniki[i, "training_AUC"] = cutoff$AUC
      tempwyniki[1, "training_AUC"] = cutoff$AUC
      #wyniki[i, "cutoff"] = wybrany_cutoff
      tempwyniki[1, "cutoff"] = wybrany_cutoff
      #message("Checkpoint passed: chunk 30")

      cat(paste0("\n\n---- TRAINING PERFORMANCE ----\n\n"))
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
      cm_train = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
      print(cm_train)
      #message("Checkpoint passed: chunk 31")

      t1_roc = pROC::roc(Class ~ as.numeric(Pred.2), data=pred)
      tempwyniki[1, "training_AUC2"] = t1_roc$auc
      tempwyniki[1, "training_AUC_lower95CI"] = as.character(ci(t1_roc))[1]
      tempwyniki[1, "training_AUC_upper95CI"] = as.character(ci(t1_roc))[3]
      saveRDS(t1_roc, paste0(temp_dir,"/models/keras",model_id,"/training_ROC.RDS"))
      #ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training_ROC.png"), grid.arrange(plot(t1_roc), nrow =1, top = "Training ROC curve"))
      #message("Checkpoint passed: chunk 32")

      tempwyniki[1, "training_Accuracy"] = cm_train$overall[1]
      tempwyniki[1, "training_Sensitivity"] = cm_train$byClass[1]
      tempwyniki[1, "training_Specificity"] = cm_train$byClass[2]
      tempwyniki[1, "training_PPV"] = cm_train$byClass[3]
      tempwyniki[1, "training_NPV"] = cm_train$byClass[4]
      tempwyniki[1, "training_F1"] = cm_train$byClass[7]
      saveRDS(cm_train, paste0(temp_dir,"/models/keras",model_id,"/cm_train.RDS"))
      #message("Checkpoint passed: chunk 33")

      cat(paste0("\n\n---- TESTING PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = test$Class, `Pred` = y_test_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
      cm_test = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
      print(cm_test)
      tempwyniki[1, "test_Accuracy"] = cm_test$overall[1]
      tempwyniki[1, "test_Sensitivity"] = cm_test$byClass[1]
      tempwyniki[1, "test_Specificity"] = cm_test$byClass[2]
      tempwyniki[1, "test_PPV"] = cm_test$byClass[3]
      tempwyniki[1, "test_NPV"] = cm_test$byClass[4]
      tempwyniki[1, "test_F1"] = cm_test$byClass[7]
      saveRDS(cm_test, paste0(temp_dir,"/models/keras",model_id,"/cm_test.RDS"))
      #message("Checkpoint passed: chunk 34")

      cat(paste0("\n\n---- VALIDATION PERFORMANCE ----\n\n"))
      pred = data.frame(`Class` = valid$Class, `Pred` = y_valid_pred)
      pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
      pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
      cm_valid = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
      print(cm_valid)
      tempwyniki[1, "valid_Accuracy"] = cm_valid$overall[1]
      tempwyniki[1, "valid_Sensitivity"] = cm_valid$byClass[1]
      tempwyniki[1, "valid_Specificity"] = cm_valid$byClass[2]
      tempwyniki[1, "valid_PPV"] = cm_valid$byClass[3]
      tempwyniki[1, "valid_NPV"] = cm_valid$byClass[4]
      tempwyniki[1, "valid_F1"] = cm_valid$byClass[7]
      saveRDS(cm_valid, paste0(temp_dir,"/models/keras",model_id,"/cm_valid.RDS"))
      #message("Checkpoint passed: chunk 35")

      if(add_features_to_predictions) {
        mix = rbind(train,test,valid)
        mixx = rbind(x_train_scale, x_test_scale, x_valid_scale)
        y_mixx_pred <- predict(object = dnn_class_model, x = mixx)

        mix$Podzial = c(rep("Training",nrow(train)),rep("Test",nrow(test)),rep("Validation",nrow(valid)))
        mix$Pred = y_mixx_pred[,2]
        mix$PredClass = ifelse(mix$Pred >= wybrany_cutoff, "Case", "Control")
        fwrite(mix, paste0(temp_dir,"/models/keras",model_id,"/data_predictions.csv.gz")) } else {
          mix = rbind(train,test,valid)
          mixx = rbind(x_train_scale, x_test_scale, x_valid_scale)
          y_mixx_pred <- predict(object = dnn_class_model, x = mixx)

          mix2 = data.frame(
            `Podzial` = c(rep("Training",nrow(train)),rep("Test",nrow(test)),rep("Validation",nrow(valid))),
            `Pred` = y_mixx_pred[,2],
            `PredClass` = ifelse(y_mixx_pred[,2] >= wybrany_cutoff, "Case", "Control")
          )

          fwrite(mix2, paste0(temp_dir,"/models/keras",model_id,"/data_predictions.csv.gz"))

        }
      fwrite(cbind(hyperparameters[i,], tempwyniki), paste0(temp_dir,"/models/keras",model_id,"/wyniki.csv"))

      wagi = keras::get_weights(dnn_class_model)
      saveRDS(wagi, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.RDS"))
      save_model_weights_hdf5(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel_weights.hdf5"))
      saveRDS(dnn_class_model, paste0(temp_dir,"/models/keras",model_id,"/finalmodel.RDS"))
      #message("Checkpoint passed: chunk 36")


      # czy jest sens zapisywac?

      cat(paste0("\n ",model_id, ": ", tempwyniki[1, "training_Accuracy"], " / ", tempwyniki[1, "test_Accuracy"], " ==> ", tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc))
      if(tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc) {
        # zapisywanie modelu do waciwego katalogu
        #message("Checkpoint passed: chunk 37e")
        if (save_all_vars) { save(list = ls(all=TRUE), file = paste0(temp_dir,"/models/keras",model_id,"/all.Rdata.gz"), compress = "gzip", compression_level = 9) }
        #message("Checkpoint passed: chunk 37d")
        if(!dir.exists(paste0("models/",codename,"/"))) { dir.create(paste0("models/",codename,"/")) }
        #message("Checkpoint passed: chunk 37c")
        if(dir.exists("/OmicSelector")) {
          system(paste0("/bin/bash -c 'screen -dmS save_network zip -9 -r ", paste0(oldwd,"/models/",codename,"/",codename, "_", model_id,".zip") ," ", paste0(temp_dir,"/models/keras",model_id), "/'"))
          # Debug:
          # OmicSelector_log(paste0("Running: /bin/bash -c 'screen -dmS save_network zip -9 -r ", paste0(oldwd,"/models/",codename,"/",codename, "_", model_id,".zip") ," ", paste0(temp_dir,"/models/keras",model_id), "/'"), "task.log")
        } else {
          zip(paste0(oldwd,"/models/",codename,"/",codename, "_", model_id,".zip"),list.files(paste0(temp_dir,"/models/keras",model_id), full.names = T, recursive = T, include.dirs = T)) }
        #message("Checkpoint passed: chunk 37b")
        if (!dir.exists(paste0("temp/",codename,"/"))) { dir.create(paste0("temp/",codename,"/")) }
        #message("Checkpoint passed: chunk 37a")
        file.copy(list.files(paste0(temp_dir,"/models/keras",model_id), pattern = "_wyniki.csv$", full.names = T, recursive = T, include.dirs = T),paste0("models/",codename,"/",model_id,"_deeplearningresults.csv"))
        #message("Checkpoint passed: chunk 37")
      } }

    sink()
    sink(type="message")
    OmicSelector_log(logfile = "task.log",  message_to_log = paste0("OmicSelector: Finished training network id: ",model_id, " : training_acc: ", tempwyniki[1, "training_Accuracy"], ", testing_acc: ", tempwyniki[1, "test_Accuracy"], " ==> worth_saving: ", tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc))
    # message(paste0("OmicSelector: Finished training network id: ",model_id, " : training_acc: ", tempwyniki[1, "training_Accuracy"], ", testing_acc: ", tempwyniki[1, "test_Accuracy"], " ==> worth_saving: ", tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc))
    #dev.off()
    tempwyniki2 = cbind(hyperparameters[i,],tempwyniki)
    tempwyniki2[1,"name"] = paste0(codename,"_", model_id)
    tempwyniki2[1,"worth_saving"] = as.character(tempwyniki[1, "training_Accuracy"]>save_threshold_trainacc & tempwyniki[1, "test_Accuracy"]>save_threshold_testacc)
    end_time <- Sys.time()
    tempwyniki2[1,"training_time"] = as.character(end_time - start_time)
    tempwyniki2
  }

  if(nrow(final) == 0) { stop("Networks are not being trained. Something is wrong.") }

  saveRDS(final, paste0(output_file,".RDS"))
  try({ OmicSelector_log("All done!! Ending batch..", "task.log") })
  if (file.exists(output_file)) {
    tempfi = data.table::fread(output_file)
    final = rbind(tempfi, final) }
  data.table::fwrite(final, output_file)
  setwd(oldwd)
  #options(warn=0)
  # sprztanie
  if(clean_temp_files) {
    try({ OmicSelector_log(paste("Cleaning temp files... Here is cluster log for reference:") ,"task.log") })
    try({ OmicSelector_log(paste(readLines(clusterlogfile), collapse="\n") ,"task.log") })
    K <- backend()
    K$clear_session()
    unlink(paste0(normalizePath(temp_dir), "/", dir(temp_dir)), recursive = TRUE)
  }

  return(final)
}

#' OmicSelector_transfer_learning_neural_network()
#'
#' This function allows to perform simple transfer learning of the network created in the process above by freezing the weights on particular layers.
#'
#' @param selected_miRNAs Which features should be selected?
#' @param new_scaling If to generate new scaling for scaled models (it is recommended for recalibration procedures). If set to false you should set old_train_csv_to_restore_scaling - used to recreate initial scaling.
#' @param model_path Point zip file to be used as ininital model for transfer learning.
#' @param save_scaling Whould you like to save scaling parameters for further reference? Default: yes
#' @param freeze_from Which layers to freeze? Set from in https://keras.rstudio.com/reference/freeze_weights.html. If set to 0 or below - the algorithm will omit freezing.
#' @param freeze_to Which layers to freeze? Set to in https://keras.rstudio.com/reference/freeze_weights.html
#' @param train Provide train set.
#' @param test Provide test set.
#' @param valid Provide validation set.
#'
#' @return Results of uncalibrated and recalibrated models.
#'
#' @export
OmicSelector_transfer_learning_neural_network = function(selected_miRNAs = ".", new_scaling = TRUE, model_path = "tcga_models/pancreatic_tcga_165592-1601307872.zip", save_scaling = TRUE, 
old_train_csv_to_restore_scaling = "../tcga/mixed_train.csv" # used if save_scaling is F and no scaling properties were saved in the model; this regenerates scaling based on previous training set.
, freeze_from = 1
, freeze_to = 2,
train = fread("circ_data/mixed_train.csv"),
test = fread("circ_data/mixed_test.csv"),
valid = fread("circ_data/mixed_valid.csv")) {
tempwyniki = data.frame()
library(keras)
# library(OmicSelector)
# OmicSelector_load_extension("deeplearning")
library(data.table)


x_train <- train %>%
  { if (selected_miRNAs[1] != ".") { dplyr::select(., selected_miRNAs) } else { dplyr::select(., starts_with("hsa")) } } %>%
  as.matrix()
y_train <- train %>%
  dplyr::select("Class") %>%
  as.matrix()
y_train[,1] = ifelse(y_train[,1] == "Case",1,0)


x_test <- test %>%
  { if (selected_miRNAs[1] != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>%
  as.matrix()
y_test <- test %>%
  dplyr::select("Class") %>%
  as.matrix()
y_test[,1] = ifelse(y_test[,1] == "Case",1,0)

x_valid <- valid %>%
  { if (selected_miRNAs[1] != ".") { dplyr::select(.,selected_miRNAs) } else { dplyr::select(.,starts_with("hsa")) } } %>%
  as.matrix()
y_valid <- valid %>%
  dplyr::select("Class") %>%
  as.matrix()
y_valid[,1] = ifelse(y_valid[,1] == "Case",1,0)



#####
# INIT MODEL

model_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("finalmodel.hdf5",Name))[1,"Name"]
unzip(model_path, model_path_in_zip, exdir = tempdir())
model_path_unzipped = paste0(tempdir(), "/", model_path_in_zip)
init_model = keras::load_model_hdf5(model_path_unzipped)
init_model
wagi_wyjsciowe = keras::get_weights(init_model)

wyniki_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("wyniki.csv",Name))[1,"Name"]
unzip(model_path, wyniki_path_in_zip, exdir = tempdir())
wyniki_path_unzipped = paste0(tempdir(), "/", wyniki_path_in_zip)
pre_conf = data.table::fread(wyniki_path_unzipped)


#
# If scaled
if(as.logical(pre_conf[1,"scaled"])) {
  if(new_scaling) {
    x_train_scale = x_train %>% scale()
    
    col_mean_train <- attr(x_train_scale, "scaled:center")
    col_sd_train <- attr(x_train_scale, "scaled:scale")
    
    saveRDS(col_mean_train, "transfer_col_mean_train.RDS")
    saveRDS(col_sd_train, "transfer_col_sd_train.RDS")
    
    x_test_scale <- x_test %>%
      scale(center = col_mean_train,
            scale = col_sd_train)
    
    x_valid_scale <- x_valid %>%
      scale(center = col_mean_train,
            scale = col_sd_train)
  } else {
    tcga_train = fread(old_train_csv_to_restore_scaling)
    tcga_train <- tcga_train %>%
      { if (selected_miRNAs[1] != ".") { dplyr::select(., selected_miRNAs) } else { dplyr::select(., starts_with("hsa")) } } %>%
      as.matrix() %>% scale()
    col_mean_train <- attr(tcga_train, "scaled:center")
    col_sd_train <- attr(tcga_train, "scaled:scale")
    
    saveRDS(col_mean_train, "transfer_col_mean_train.RDS")
    saveRDS(col_sd_train, "transfer_col_sd_train.RDS")
    
    
    x_train_scale <- x_train %>%
      scale(center = col_mean_train,
            scale = col_sd_train)
    
    x_test_scale <- x_test %>%
      scale(center = col_mean_train,
            scale = col_sd_train)
    
    x_valid_scale <- x_valid %>%
      scale(center = col_mean_train,
            scale = col_sd_train)
  } 
} else {
  x_train_scale <- x_train
  x_test_scale <- x_test 
  x_valid_scale <- x_valid 
}


preds = predict(init_model, x_train_scale)[,2]
library(pROC)


# wybranie odciecia
pred = data.frame(`Class` = as.factor(ifelse(y_train==1,"Case","Control")), `Pred` = predict(init_model, x_train_scale))
library(cutpointr)
cutoff = cutpointr(pred, Pred.2, Class, pos_class = "Case", metric = youden)
print(summary(cutoff))
ggplot2::ggsave(file = paste0("cutoff.png"), plot(cutoff))
wybrany_cutoff = cutoff$optimal_cutpoint


y_train_pred <- predict(object = init_model, x = x_train_scale)
y_test_pred <- predict(object = init_model, x = x_test_scale)
y_valid_pred <- predict(object = init_model, x = x_valid_scale)
#message("Checkpoint passed: chunk 29")

#wyniki[i, "training_AUC"] = cutoff$AUC
tempwyniki[1, "new_training_AUC"] = cutoff$AUC
#wyniki[i, "cutoff"] = wybrany_cutoff
tempwyniki[1, "new_cutoff"] = wybrany_cutoff
#message("Checkpoint passed: chunk 30")

cat(paste0("\n\n---- TRAINING PERFORMANCE ----\n\n"))
pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
cm_train = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
print(cm_train)
#message("Checkpoint passed: chunk 31")

t1_roc = pROC::roc(Class ~ as.numeric(Pred.2), data=pred)
tempwyniki[1, "new_training_AUC2"] = t1_roc$auc
tempwyniki[1, "new_training_AUC_lower95CI"] = as.character(ci(t1_roc))[1]
tempwyniki[1, "new_training_AUC_upper95CI"] = as.character(ci(t1_roc))[3]
saveRDS(t1_roc, paste0("init_training_ROC.RDS"))
#ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training_ROC.png"), grid.arrange(plot(t1_roc), nrow =1, top = "Training ROC curve"))
#message("Checkpoint passed: chunk 32")

tempwyniki[1, "new_training_Accuracy"] = cm_train$overall[1]
tempwyniki[1, "new_training_Sensitivity"] = cm_train$byClass[1]
tempwyniki[1, "new_training_Specificity"] = cm_train$byClass[2]
tempwyniki[1, "new_training_PPV"] = cm_train$byClass[3]
tempwyniki[1, "new_training_NPV"] = cm_train$byClass[4]
tempwyniki[1, "new_training_F1"] = cm_train$byClass[7]
saveRDS(cm_train, paste0("init_cm_train.RDS"))
#message("Checkpoint passed: chunk 33")

cat(paste0("\n\n---- TESTING PERFORMANCE ----\n\n"))
pred = data.frame(`Class` = as.factor(test$Class), `Pred` = y_test_pred)
pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
cm_test = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
print(cm_test)
tempwyniki[1, "new_test_Accuracy"] = cm_test$overall[1]
tempwyniki[1, "new_test_Sensitivity"] = cm_test$byClass[1]
tempwyniki[1, "new_test_Specificity"] = cm_test$byClass[2]
tempwyniki[1, "new_test_PPV"] = cm_test$byClass[3]
tempwyniki[1, "new_test_NPV"] = cm_test$byClass[4]
tempwyniki[1, "new_test_F1"] = cm_test$byClass[7]
saveRDS(cm_test, paste0("init_cm_test.RDS"))
#message("Checkpoint passed: chunk 34")

cat(paste0("\n\n---- VALIDATION PERFORMANCE ----\n\n"))
pred = data.frame(`Class` = as.factor(valid$Class), `Pred` = y_valid_pred)
pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
cm_valid = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
print(cm_valid)
tempwyniki[1, "new_valid_Accuracy"] = cm_valid$overall[1]
tempwyniki[1, "new_valid_Sensitivity"] = cm_valid$byClass[1]
tempwyniki[1, "new_valid_Specificity"] = cm_valid$byClass[2]
tempwyniki[1, "new_valid_PPV"] = cm_valid$byClass[3]
tempwyniki[1, "new_valid_NPV"] = cm_valid$byClass[4]
tempwyniki[1, "new_valid_F1"] = cm_valid$byClass[7]
saveRDS(cm_valid, paste0("init_cm_valid.RDS"))
#message("Checkpoint passed: chunk 35")


#####
# RECALIBRATION
scenario_i = 2
scenario = "trans1"
trans_model = clone_model(init_model)


early_stop <- callback_early_stopping(monitor = "val_loss", mode="min", patience = 200)
cp_callback <- callback_model_checkpoint(
  filepath =  paste0(scenario,"_model.hdf5"),
  save_best_only = TRUE, period = 10, monitor = "val_loss",
  verbose = 0
)
keras_batch_size = 64



# Opcje transferu:
compile(trans_model, loss = 'binary_crossentropy',
        metrics = 'accuracy', optimizer = pre_conf$optimizer)

keras::set_weights(trans_model, keras::get_weights(init_model))
if(freeze_from>0) {
freeze_weights(trans_model, from = freeze_from, to = freeze_to)
compile(trans_model, loss = 'binary_crossentropy',
        metrics = 'accuracy', optimizer = pre_conf$optimizer) }

# Rekalibracja:
history <-  fit(trans_model, x = x_train_scale,
                y = to_categorical(y_train),
                epochs = 5000,
                validation_data = list(x_test_scale, to_categorical(y_test)),
                callbacks = list(
                  cp_callback,
                  #callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.1),
                  #callback_model_checkpoint(paste0(temp_dir,"/models/keras",model_id,"/finalmodel.hdf5")),
                  early_stop),
                verbose = 0,
                view_metrics = FALSE,
                batch_size = keras_batch_size, shuffle = T)
history

keras::get_weights(init_model)
keras::get_weights(trans_model)

trans_model %>% save_model_hdf5("trans_model.hdf5")

# Ocena po rekalibracji:
preds = predict(trans_model, x_train_scale)[,2]
library(pROC)


# wybranie odciecia
pred = data.frame(`Class` = as.factor(ifelse(y_train==1,"Case","Control")), `Pred` = predict(trans_model, x_train_scale))
library(cutpointr)
cutoff = cutpointr(pred, Pred.2, Class, pos_class = "Case", metric = youden)
print(summary(cutoff))
ggplot2::ggsave(file = paste0("trans_cutoff.png"), plot(cutoff))
wybrany_cutoff = cutoff$optimal_cutpoint


y_train_pred <- predict(object = trans_model, x = x_train_scale)
y_test_pred <- predict(object = trans_model, x = x_test_scale)
y_valid_pred <- predict(object = trans_model, x = x_valid_scale)
#message("Checkpoint passed: chunk 29")

#wyniki[i, "training_AUC"] = cutoff$AUC
tempwyniki[scenario_i, "new_training_AUC"] = cutoff$AUC
#wyniki[i, "cutoff"] = wybrany_cutoff
tempwyniki[scenario_i, "new_cutoff"] = wybrany_cutoff
#message("Checkpoint passed: chunk 30")

cat(paste0("\n\n---- TRAINING PERFORMANCE ----\n\n"))
pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
cm_train = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
print(cm_train)
#message("Checkpoint passed: chunk 31")

t1_roc = pROC::roc(Class ~ as.numeric(Pred.2), data=pred)
tempwyniki[scenario_i, "new_training_AUC2"] = t1_roc$auc
tempwyniki[scenario_i, "new_training_AUC_lower95CI"] = as.character(ci(t1_roc))[1]
tempwyniki[scenario_i, "new_training_AUC_upper95CI"] = as.character(ci(t1_roc))[3]
saveRDS(t1_roc, paste0("new_trans_training_ROC.RDS"))
#ggplot2::ggsave(file = paste0(temp_dir,"/models/keras",model_id,"/training_ROC.png"), grid.arrange(plot(t1_roc), nrow =1, top = "Training ROC curve"))
#message("Checkpoint passed: chunk 32")

tempwyniki[scenario_i, "new_training_Accuracy"] = cm_train$overall[1]
tempwyniki[scenario_i, "new_training_Sensitivity"] = cm_train$byClass[1]
tempwyniki[scenario_i, "new_training_Specificity"] = cm_train$byClass[2]
tempwyniki[scenario_i, "new_training_PPV"] = cm_train$byClass[3]
tempwyniki[scenario_i, "new_training_NPV"] = cm_train$byClass[4]
tempwyniki[scenario_i, "new_training_F1"] = cm_train$byClass[7]
saveRDS(cm_train, paste0("trans_cm_train.RDS"))
#message("Checkpoint passed: chunk 33")

cat(paste0("\n\n---- TESTING PERFORMANCE ----\n\n"))
pred = data.frame(`Class` = as.factor(test$Class), `Pred` = y_test_pred)
pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
cm_test = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
print(cm_test)
tempwyniki[scenario_i, "new_test_Accuracy"] = cm_test$overall[1]
tempwyniki[scenario_i, "new_test_Sensitivity"] = cm_test$byClass[1]
tempwyniki[scenario_i, "new_test_Specificity"] = cm_test$byClass[2]
tempwyniki[scenario_i, "new_test_PPV"] = cm_test$byClass[3]
tempwyniki[scenario_i, "new_test_NPV"] = cm_test$byClass[4]
tempwyniki[scenario_i, "new_test_F1"] = cm_test$byClass[7]
saveRDS(cm_test, paste0("trans_cm_test.RDS"))
#message("Checkpoint passed: chunk 34")

cat(paste0("\n\n---- VALIDATION PERFORMANCE ----\n\n"))
pred = data.frame(`Class` = as.factor(valid$Class), `Pred` = y_valid_pred)
pred$PredClass = ifelse(pred$Pred.2 >= wybrany_cutoff, "Case", "Control")
pred$PredClass = factor(pred$PredClass, levels = c("Control","Case"))
cm_valid = caret::confusionMatrix(pred$PredClass, pred$Class, positive = "Case")
print(cm_valid)
tempwyniki[scenario_i, "new_valid_Accuracy"] = cm_valid$overall[1]
tempwyniki[scenario_i, "new_valid_Sensitivity"] = cm_valid$byClass[1]
tempwyniki[scenario_i, "new_valid_Specificity"] = cm_valid$byClass[2]
tempwyniki[scenario_i, "new_valid_PPV"] = cm_valid$byClass[3]
tempwyniki[scenario_i, "new_valid_NPV"] = cm_valid$byClass[4]
tempwyniki[scenario_i, "new_valid_F1"] = cm_valid$byClass[7]
saveRDS(cm_valid, paste0("trans_cm_valid.RDS"))
#message("Checkpoint passed: chunk 35")


resu = cbind(pre_conf, tempwyniki)
resu$type = c("crude","recalibrated")
resu$based_on = basename(model_path)
resu$based_on_path = model_path
new_name = make.names(paste0("transfer-",as.numeric(Sys.time())))
resu$new_model_name = c("",new_name)
fwrite(resu, "resu.csv")

zip(zipfile = paste0(new_name,".zip"), files = c(list.files(".", pattern = ".hdf5"), list.files(".", pattern = ".csv"), list.files(".", pattern = ".RDS"), list.files(".", pattern = ".png")))
unlink(c(list.files(".", pattern = ".hdf5"), list.files(".", pattern = ".csv"), list.files(".", pattern = ".RDS"), list.files(".", pattern = ".png")))

resu
}

#' OmicSelector_deep_learning_predict()
#'
#' This function allows to predict using developed deep neural networks.
#'
#' @param model_path Path to model zip generated by OmicSelector.
#' @param new_dataset New dataset to perform prediction on.
#' @param new_scaling If the network is scaled, do you want to perform new scaling of your data?
#' @param old_train_csv_to_restore_scaling Default: NULL. Path to csv file the scaling can be restored from. This is optional and allows you to play with scaling. If you deleted col_mean_train and col_sd_train files from your model zip you can setup your own scaling factors by providing the file to this variable.
#' @param override_cutoff Default: NULL. Almost always you would like to use the cutoff develped in the process of the network training, but if you wish to use your own cutoff you can play with this parameter. Otherwise, keep it NULL.
#' @param blinded Is your validation blinded? This function can assess the predictive abilities by comparing the predictions to variable 'Class'. If your data doesn't have one or you are performing blinded validatin - set it to TRUE. Otherwise, keep it FALSE.
#'
#' @return Results of validation.
#'
#' @export
#'
OmicSelector_deep_learning_predict = function(model_path = "our_models/model5.zip",
                                              new_dataset = data.table::fread("Data/ks_data.csv"),
                                              new_scaling = F,
                                              old_train_csv_to_restore_scaling = NULL,
                                              override_cutoff = NULL,
                                              blinded = F
) {
  # Load model data
  library(keras)
  library(reticulate)
  model_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("finalmodel.hdf5",Name))[1,"Name"]
  unzip(model_path, model_path_in_zip, exdir = tempdir())
  model_path_unzipped = paste0(tempdir(), "/", model_path_in_zip)

  if(!file.exists(model_path_unzipped)) { stop("Model file does not exist in provided file. Is it the correct one?") }

  init_model = keras::load_model_hdf5(model_path_unzipped)
  print(init_model)


  wyniki_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("wyniki.csv",Name))[1,"Name"]
  unzip(model_path, wyniki_path_in_zip, exdir = tempdir())
  wyniki_path_unzipped = paste0(tempdir(), "/", wyniki_path_in_zip)

  if(!file.exists(wyniki_path_unzipped)) { stop("Configuration file does not exist in model package. Is it damaged?") }

  pre_conf = data.table::fread(wyniki_path_unzipped)
  #colnames(pre_conf)

  # Match variables
  library(dplyr)
  network_features = strsplit(x = as.character(pre_conf[1,"formula"]),split = " + ", fixed = T)[[1]]

  if(sum(is.na(match(network_features, colnames(new_dataset))))>0) {
    stop(paste0("The new dataset does not contain features: ", paste0(network_features[which(is.na(match(network_features, colnames(new_dataset))))], collapse = ", ")))
  }



  new_x <- new_dataset %>% dplyr::select(., all_of(network_features)) %>% as.matrix()
  if(blinded == FALSE) {
    new_y <- new_dataset %>%
      dplyr::select("Class") %>%
      as.matrix()
    new_y[,1] = ifelse(new_y[,1] == "Case",1,0)
    new_dataset$Class = factor(new_dataset$Class, levels = c("Control","Case")) }



  # If scaled
  col_mean_train = NA
  col_sd_train = NA
  if (pre_conf[1,"scaled"] == "TRUE") {
    cat("This network is scaled.")
    temp_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("col_mean_train.RDS",Name))[1,"Name"]
    unzip(model_path, temp_path_in_zip, exdir = tempdir())
    temp_path_unzipped = paste0(tempdir(), "/", temp_path_in_zip)
    if (file.exists(temp_path_unzipped)) { col_mean_train = readRDS(temp_path_unzipped) } else { cat(" Scaling mean not saved in model.") }

    temp_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("col_sd_train.RDS",Name))[1,"Name"]
    unzip(model_path, temp_path_in_zip, exdir = tempdir())
    temp_path_unzipped = paste0(tempdir(), "/", temp_path_in_zip)
    if (file.exists(temp_path_unzipped)) { col_sd_train = readRDS(temp_path_unzipped) } else { cat(" Scaling SD not saved in model.") }




    if(new_scaling) {
      new_x = new_x %>% scale()
      new_x[is.nan(new_x)] <- 0
      cat("\n\nNew scaling was performed with col_mean: ")
      print(attr(new_x, "scaled:center"))
      cat("\n\nNew col_sd: ")
      print(attr(new_x, "scaled:scale"))
    } else {
      if(is.null(old_train_csv_to_restore_scaling)) {
        new_x = new_x %>%
          scale(center = col_mean_train,
                scale = col_sd_train)
        new_x[is.nan(new_x)] <- 0
      } else {

        if(!file.exists(old_train_csv_to_restore_scaling)) { stop(" A file to restore scaling from does not exist.") }
        temp_train = data.table::fread(old_train_csv_to_restore_scaling)
        temp_train_x = dplyr::select(temp_train, all_of(network_features))
        temp_train_x = temp_train_x %>% scale()
        col_mean_train <- attr(temp_train_x, "scaled:center")
        col_sd_train <- attr(temp_train_x, "scaled:scale")
        new_x = new_x %>%
          scale(center = col_mean_train,
                scale = col_sd_train)
        new_x[is.nan(new_x)] <- 0
        cat("\n\nScaling was restored from provided file. Final col_mean: ")
        print(attr(new_x, "scaled:center"))
        cat("\n\nFinal col_sd: ")
        print(attr(new_x, "scaled:scale"))
      }}
  } # koniec jesli scaled

  model_autoencoder = NULL
  if (pre_conf[1,"autoencoder"] != "0") {
      model_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("autoencoder.hdf5",Name))[1,"Name"]
      unzip(model_path, model_path_in_zip, exdir = tempdir())
      model_path_unzipped = paste0(tempdir(), "/", model_path_in_zip)

      if(!file.exists(model_path_unzipped)) { stop("Autoencoder model file does not exist in provided file. Is it the correct one?") }

      model_autoencoder = keras::load_model_hdf5(model_path_unzipped)
      print(model_autoencoder)

      new_x_old = new_x
      new_x <- model_autoencoder %>%
        predict(new_x) %>%
        as.matrix()
  }


  # Przewidywanie
  library(pROC)
  cutoff = as.numeric(pre_conf[1,"cutoff"])
  if(!is.null(override_cutoff)) { cutoff = as.numeric(override_cutoff) }
  if(is.na(cutoff)) { stop("Why is cutoff not numeric? Check your override_cutoff parameter or the value set in the model.") }

  if(blinded){
    predictions = predict(init_model, new_x)
    pred = data.frame(`Pred` = predictions[,2])
    pred$Prediction = ifelse(pred$Pred >= cutoff, "Case", "Control")
    confusion_matrix = NA
    roc = NA
    roc_auc = NA

  } else {
    predictions = predict(init_model, new_x)
    pred = data.frame(`Class` = factor(ifelse(new_y==1,"Case","Control"), levels = c("Control","Case")), `Pred` = predictions[,2])
    pred$Prediction = ifelse(pred$Pred >= cutoff, "Case", "Control")
    pred$Correctness = ifelse(pred$Prediction == pred$Class, "Correct", "Incorrect")
    confusion_matrix = caret::confusionMatrix(as.factor(pred$Prediction), as.factor(pred$Class), positive = "Case")
    roc = pROC::roc(pred$Class ~ pred$Pred)
    roc_auc = pROC::ci.auc(roc(pred$Class ~ pred$Pred))
  }



  final_return = list(`predictions` = pred,
                      `network_config` = pre_conf,
                      `new_dataset` = new_dataset,
                      `new_dataset_x` = new_x,
                      `network_features` = network_features,
                      `cutoff` = cutoff,
                      `col_mean_train` = col_mean_train,
                      `col_sd_train` = col_sd_train,
                      `confusion_matrix` = confusion_matrix,
                      `roc` = roc,
                      `roc_auc` = roc_auc,
                      `model` = init_model,
                      `autoencoder` = model_autoencoder)
  final_return
}

OmicSelector_deep_learning_predict_transfered = function(model_path = "our_models/model5.zip",
                                                         new_dataset = data.table::fread("Data/ks_data.csv"),
                                                         new_scaling = T,
                                                         old_train_csv_to_restore_scaling = NULL,
                                                         override_cutoff = NULL,
                                                         blinded = T) {
  # Load model data
  library(keras)
  library(reticulate)
  model_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("trans_model.hdf5",Name))[1,"Name"]
  unzip(model_path, model_path_in_zip, exdir = tempdir())
  model_path_unzipped = paste0(tempdir(), "/", model_path_in_zip)
  
  if(!file.exists(model_path_unzipped)) { stop("Model file does not exist in provided file. Is it the correct one?") }
  
  init_model = keras::load_model_hdf5(model_path_unzipped)
  print(init_model)
  
  
  wyniki_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("resu.csv",Name))[1,"Name"]
  unzip(model_path, wyniki_path_in_zip, exdir = tempdir())
  wyniki_path_unzipped = paste0(tempdir(), "/", wyniki_path_in_zip)
  
  if(!file.exists(wyniki_path_unzipped)) { stop("Configuration file does not exist in model package. Is it damaged?") }
  
  pre_conf = data.table::fread(wyniki_path_unzipped)
  #colnames(pre_conf)
  
  # Match variables
  library(dplyr)
  network_features = strsplit(x = as.character(pre_conf[1,"formula"]),split = " + ", fixed = T)[[1]]
  
  if(sum(is.na(match(network_features, colnames(new_dataset))))>0) {
    stop(paste0("The new dataset does not contain features: ", paste0(network_features[which(is.na(match(network_features, colnames(new_dataset))))], collapse = ", ")))
  }
  
  
  
  new_x <- new_dataset %>% dplyr::select(., all_of(network_features)) %>% as.matrix()
  if(blinded == FALSE) {
    
    new_y <- new_dataset %>%
      dplyr::select("Class") %>%
      as.matrix()
    new_y[,1] = ifelse(new_y[,1] == "Case",1,0)
    new_dataset$Class = factor(new_dataset$Class, levels = c("Control","Case")) }
  
  
  
  # If scaled
  col_mean_train = NA
  col_sd_train = NA
  if (pre_conf[1,"scaled"] == "TRUE") {
    cat("This network is scaled.")
    temp_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("transfer_col_mean_train.RDS",Name))[1,"Name"]
    unzip(model_path, temp_path_in_zip, exdir = tempdir())
    temp_path_unzipped = paste0(tempdir(), "/", temp_path_in_zip)
    if (file.exists(temp_path_unzipped)) { col_mean_train = readRDS(temp_path_unzipped) } else { cat(" Scaling mean not saved in model.") }
    
    temp_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("transfer_col_sd_train.RDS",Name))[1,"Name"]
    unzip(model_path, temp_path_in_zip, exdir = tempdir())
    temp_path_unzipped = paste0(tempdir(), "/", temp_path_in_zip)
    if (file.exists(temp_path_unzipped)) { col_sd_train = readRDS(temp_path_unzipped) } else { cat(" Scaling SD not saved in model.") }
    
    
    
    
    if(new_scaling) {
      new_x = new_x %>% scale()
      new_x[is.nan(new_x)] <- 0
      cat("\n\nNew scaling was performed with col_mean: ")
      print(attr(new_x, "scaled:center"))
      cat("\n\nNew col_sd: ")
      print(attr(new_x, "scaled:scale"))
    } else {
      if(is.null(old_train_csv_to_restore_scaling)) {
        new_x = new_x %>%
          scale(center = col_mean_train,
                scale = col_sd_train)
        new_x[is.nan(new_x)] <- 0
      } else {
        
        if(!file.exists(old_train_csv_to_restore_scaling)) { stop(" A file to restore scaling from does not exist.") }
        temp_train = data.table::fread(old_train_csv_to_restore_scaling)
        temp_train_x = dplyr::select(temp_train, all_of(network_features))
        temp_train_x = temp_train_x %>% scale()
        col_mean_train <- attr(temp_train_x, "scaled:center")
        col_sd_train <- attr(temp_train_x, "scaled:scale")
        new_x = new_x %>%
          scale(center = col_mean_train,
                scale = col_sd_train)
        new_x[is.nan(new_x)] <- 0
        cat("\n\nScaling was restored from provided file. Final col_mean: ")
        print(attr(new_x, "scaled:center"))
        cat("\n\nFinal col_sd: ")
        print(attr(new_x, "scaled:scale"))
      }}
  } # koniec jesli scaled

    model_autoencoder = NULL
    if (pre_conf[1,"autoencoder"] != "0") {
      model_path_in_zip = dplyr::filter(unzip(model_path, list = T), grepl("autoencoder.hdf5",Name))[1,"Name"]
      unzip(model_path, model_path_in_zip, exdir = tempdir())
      model_path_unzipped = paste0(tempdir(), "/", model_path_in_zip)

      if(!file.exists(model_path_unzipped)) { stop("Autoencoder model file does not exist in provided file. Is it the correct one?") }

      model_autoencoder = keras::load_model_hdf5(model_path_unzipped)
      print(model_autoencoder)

      new_x_old = new_x
      new_x <- model_autoencoder %>%
        predict(new_x) %>%
        as.matrix()
    }
  
  # Przewidywanie
  library(pROC)
  cutoff = as.numeric(pre_conf$new_cutoff[2])
  if(!is.null(override_cutoff)) { cutoff = as.numeric(override_cutoff) }
  if(is.na(cutoff)) { stop("Why is cutoff not numeric? Check your override_cutoff parameter or the value set in the model.") }
  
  if(blinded){
    predictions = predict(init_model, new_x)
    pred = data.frame(`Pred` = predictions[,2])
    pred$Prediction = ifelse(pred$Pred >= cutoff, "Case", "Control")
    confusion_matrix = NA
    roc = NA
    roc_auc = NA
    
  } else {
    predictions = predict(init_model, new_x)
    pred = data.frame(`Class` = factor(ifelse(new_y==1,"Case","Control"), levels = c("Control","Case")), `Pred` = predictions[,2])
    pred$Prediction = ifelse(pred$Pred >= cutoff, "Case", "Control")
    pred$Correctness = ifelse(pred$Prediction == pred$Class, "Correct", "Incorrect")
    confusion_matrix = caret::confusionMatrix(as.factor(pred$Prediction), as.factor(pred$Class), positive = "Case")
    roc = pROC::roc(pred$Class ~ pred$Pred)
    roc_auc = pROC::ci.auc(roc(pred$Class ~ pred$Pred))
  }
  
  
  
  final_return = list(`predictions` = pred,
                      `network_config` = pre_conf,
                      `new_dataset` = new_dataset,
                      `new_dataset_x` = new_x,
                      `network_features` = network_features,
                      `cutoff` = cutoff,
                      `col_mean_train` = col_mean_train,
                      `col_sd_train` = col_sd_train,
                      `confusion_matrix` = confusion_matrix,
                      `roc` = roc,
                      `roc_auc` = roc_auc)
  final_return
}