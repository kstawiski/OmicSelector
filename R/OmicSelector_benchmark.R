#' OmicSelector_benchmark
#'
#' Second most important function in the package.
#' Using the formulas selected by `OmicSelector_OmicSelector` function, it test derived miRNA sets in a systematic manner using multiple model induction methods.
#' This function allows to benchmark miRNA sets in context of their potential for diagnostic test creation.
#' Hidden feature of this package is application of `mxnet`. Note that `mxnet` has to be installed and configured seperatly.
#'
#' @param wd Working directory here `OmicSelector_OmicSelector` was also working.
#' @param input_fomulas List of formulas as created by `OmicSelector_OmicSelector` or `OmicSelector_merge_formulas`. Those formulas will be check in benchmark.
#' @param algorithms Caret methods that will be checked in benchmark processing. By default the logistic regression is always included.
#' @param output_file Out csv file for the benchmark.
#' @param holdout Best set of hyperparameters can be selected using: (1) if TURE - using hold-out validation on test set, (2) if FALSE - using 10-fold cross-validation repeated 5 times.
#' @param stamp Character vector or timestamp to make the benchmark unique.
#' @param cores Number of cores using in parallel processing.
#' @param mxnet Whether to use mxnet. Default: F
#' @param search_iters_mxnet Number of iterations in mxnet-based neural network creation. Default: 5000
#' @param gpu Wheter to use GPU in mxnet and keras processing. Default: F
#' @param search_iters The number of random hyperparameters tested in the process of model induction.
#' @param keras_epochs Number of epochs used in keras-based methods, if keras methods are used. (e.g. "mlpKerasDropout", "mlpKerasDecay")
#' @param keras_threds This package supports training of keras networks in parallel. Here you can set the number of threads used. (e.g. "mlpKerasDropout", "mlpKerasDecay")
#'
#' @return Results of benchmark. Note that benchmark files are also saved in working directory (`wd`).
#'
#' @export
OmicSelector_benchmark = function(wd = getwd(), search_iters = 2000, keras_epochs = 5000, keras_threads = floor(parallel::detectCores()/2), search_iters_mxnet = 5000,
                        cores = detectCores()-1, input_formulas = readRDS("featureselection_formulas_final.RDS"),
                        output_file = "benchmark.csv", mxnet = F, gpu = F,
                        #algorithms = c("nnet","svmRadial", "svmLinear","rf","C5.0","mlp", "mlpML","xgbTree"),
                        algorithms = c("mlp", "mlpML", "svmRadial", "svmLinear","rf","C5.0", "rpart", "rpart2", "ctree"),
                        holdout = T, stamp = as.character(as.numeric(Sys.time()))) {
  suppressMessages(library(plyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(edgeR))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(rsq))
  suppressMessages(library(MASS))
  suppressMessages(library(Biocomb))
  suppressMessages(library(caret))
  suppressMessages(library(dplyr))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(pROC))
  suppressMessages(library(ggplot2))
  suppressMessages(library(DMwR))
  suppressMessages(library(stringr))
  suppressMessages(library(psych))
  suppressMessages(library(C50))
  suppressMessages(library(randomForest))
  suppressMessages(library(nnet))
  suppressMessages(library(reticulate))
  suppressMessages(library(stargazer))
  suppressMessages(library(plyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(edgeR))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(rsq))
  suppressMessages(library(MASS))
  suppressMessages(library(Biocomb))
  suppressMessages(library(caret))
  suppressMessages(library(dplyr))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(pROC))
  suppressMessages(library(ggplot2))
  suppressMessages(library(DMwR))
  suppressMessages(library(ROSE))
  suppressMessages(library(gridExtra))
  suppressMessages(library(gplots))
  suppressMessages(library(remotes))
  suppressMessages(library(stringr))
  suppressMessages(library(data.table))
  suppressMessages(library(tidyverse))


  if("mxnet" %in% rownames(installed.packages()) == FALSE && mxnet == T) {
    stop("Mxnet R package is not installed. Please set mxnet to FALSE or build and install mxnet R package. If you don't know how, just use our docker-based enviorment.")
  }

  if(!dir.exists("temp")) { dir.create("temp") }

  use_condaenv("tensorflow")
  zz <- file(paste0("temp/benchmark",stamp,".log"), open = "wt")
  pdf(paste0("temp/benchmark",stamp,".pdf"))
  sink(zz)
  sink(zz, type = "message")
  suppressMessages(library(doParallel))


  oldwd = getwd()
  setwd(wd)

  formulas = input_formulas

  wybrane = list()
  ile_miRNA = numeric()
  for (i in 1:length(formulas)) {
    wybrane[[i]] = all.vars(as.formula(formulas[[i]]))[-1]
    ile_miRNA[i] = length(all.vars(as.formula(formulas[[i]])))-1
  }

  hist(ile_miRNA, breaks = 150)
  psych::describe(ile_miRNA)

  # Wywalamy formuły z więcej niż 15 miRNA
  #formulas_old = formulas
  #formulas = formulas[which(ile_miRNA <= 15)]
  #print(formulas)

  dane = OmicSelector_load_datamix(wd = wd, replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  if(!dir.exists("models")) { dir.create("models") }

  wyniki = data.frame(method = names(formulas))
  wyniki$SMOTE = ifelse(grepl("SMOTE", names(formulas)),"Yes","No")


  if(mxnet == T) {
    suppressMessages(library(mxnet))
    #gpu = system("nvidia-smi", intern = T)
    if (gpu == T) {
      mxctx = mx.gpu()
      print("MXNET is using GPU not CPU :)")
    } else {
      mxctx = mx.cpu()
      print("MXNET is using CPU not GPU :(")
    }


    algorytmy = c("glm", "mxnetAdam", algorithms)
  } else {
    algorytmy = c("glm",algorithms)
  }





  for (ii in 1:length(algorytmy)) {
    algorytm = algorytmy[ii]
    print(algorytm)

    suppressMessages(library(doParallel))
    suppressMessages(library(doParallel))
    if(grepl("Keras",algorytm)) {
      cl <- makePSOCKcluster(useXDR = TRUE, keras_threads)
      registerDoParallel(cl)
      # on.exit(stopCluster(cl))
    } else {
      cl <- makePSOCKcluster(useXDR = TRUE, cores-1)
      registerDoParallel(cl)
      # on.exit(stopCluster(cl))
    }

    for (i in 1:length(formulas)) {
      print(paste0("Testing formula: ", as.character(formulas[[i]])))

      wyniki$miRy[i] = as.character(formulas[[i]])
      temptrain = train
      if (wyniki$SMOTE[i] == "Yes") { temptrain = train_smoted }

      # Hold-out czy CV?
      temptrainold = temptrain[complete.cases(temptrain), ] # Keep only the complete rows - useful for SMOTE on low number of samples
      if(holdout == T) {
        #fit_on = list(rs1 = 1:nrow(temptrain), rs2 = 1:nrow(temptrain))
        #pred_on = list(rs1 = (nrow(temptrain)+1):((nrow(temptrain)+1)+nrow(test)), rs2 = ((nrow(temptrain)+1)+nrow(test)+1):((nrow(temptrain)+1)+nrow(test)+1+nrow(valid)))
        #temptrain = rbind.fill(temptrain,test,valid)
        fit_on = list(rs1 = 1:nrow(temptrain))
        pred_on = list(rs1 = (nrow(temptrain)+1):((nrow(temptrain))+nrow(test)))
        temptrain = rbind.fill(temptrain,test)
      }

     #Debug:
     temptrain = temptrain[complete.cases(temptrain), ] # Keep only the complete rows - useful for SMOTE on low number of samples
     head(temptrain)



      # wyniki2 = tryCatch({
      if(algorytm == "mxnet") {
        hyperparameters = expand.grid(layer1 = unique(ceiling(seq(2, ncol(trainx)/2, length.out = 10))), layer2 = unique(ceiling(seq(0, ncol(trainx)/4, length.out = 5))), layer3 =unique(ceiling(seq(0, ncol(trainx)/8, length.out = 5))), activation = c('relu', 'sigmoid'),
                                      dropout = c(0,0.05),learning.rate= c(0.00001, 0.01), momentum = c(0, 0.8, 0.99))
        train_control <- trainControl(method="cv", repeats=5, number = 10, classProbs = TRUE,verboseIter = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), ctx = mxctx, optimizer = 'sgd',
                              #optimizer_params=(('learning_rate',0.1),('lr_scheduler',lr_sch)),
                              preProc = c("center", "scale"),
                              #epoch.end.callback = mx.callback.early.stop(5,10,NULL, maximize=TRUE, verbose=TRUE),
                              #eval.data = list(data=dplyr::select(test, starts_with("hsa")),label=dplyr::select(test, Class)),
                              epoch.end.callback=mx.callback.early.stop(30, 30),
                              #preProc = c("center", "scale"),
                              num.round = search_iters_mxnet, data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)
        print(model1$finalModel)
      } else if(algorytm == "mxnetAdam") {
        hyperparameters = expand.grid(layer1 = unique(ceiling(seq(2, ncol(trainx)/2, length.out = 10))), layer2 = unique(ceiling(seq(0, ncol(trainx)/4, length.out = 5))), layer3 =unique(ceiling(seq(0, ncol(trainx)/8, length.out = 5))), activation = c('relu', 'sigmoid'),
                                      dropout = c(0,0.05), beta1=0.9, beta2=0.999, learningrate= c(0.001))
        train_control <- trainControl(method="cv", repeats=5, number = 10, classProbs = TRUE,verboseIter = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), ctx = mxctx,
                              #optimizer_params=(('learning_rate',0.1),('lr_scheduler',lr_sch)),
                              preProc = c("center", "scale"),
                              #epoch.end.callback = mx.callback.early.stop(5,10,NULL, maximize=TRUE, verbose=TRUE),
                              #eval.data = list(data=dplyr::select(test, starts_with("hsa")),label=dplyr::select(test, Class)),
                              epoch.end.callback=mx.callback.early.stop(30, 30),
                              #preProc = c("center", "scale"),
                              num.round = search_iters_mxnet, data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)
        print(model1$finalModel)
      } else if(grepl("Keras",algorytm)) {
        train_control <- trainControl(method="cv", number = 5, search="random", classProbs = TRUE, verboseIter = TRUE,
                                      summaryFunction = twoClassSummary, savePredictions = TRUE, indexFinal = fit_on[[1]])
        #if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
        #                                                 classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, tuneLength = search_iters,
                              epochs = keras_epochs)
        print(model1$finalModel)
      } else if (algorytm == "glm") {
        #train_control <- trainControl(method="repeatedcv", repeats=5, number = 10, search="random", classProbs = TRUE, verboseIter = TRUE,
        #                              summaryFunction = twoClassSummary, savePredictions = TRUE)
        #if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
        #                                                 classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        train_control = trainControl(method= "none")
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method="glm", family="binomial")
        print(model1$finalModel)
      } else {
        train_control <- trainControl(method="repeatedcv", repeats=5, number = 10, search="random", classProbs = TRUE, verboseIter = TRUE,
                                      summaryFunction = twoClassSummary, savePredictions = TRUE)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, tuneLength = search_iters)
        print(model1$finalModel)
      }


      modelname = as.numeric(Sys.time())
      print(paste0("MODELID: ", modelname))
      saveRDS(model1, paste0("models/",modelname,".RDS"))
      wyniki[i,paste0(algorytm,"_modelname")] = modelname

      t1_roc = roc(temptrainold$Class ~ predict(model1, newdata = temptrainold , type = "prob")[,2])
      wyniki[i,paste0(algorytm,"_train_ROCAUC")] = t1_roc$auc
      wyniki[i,paste0(algorytm,"_train_ROCAUC_lower95CI")] = as.character(ci(t1_roc))[1]
      wyniki[i,paste0(algorytm,"_train_ROCAUC_upper95CI")] = as.character(ci(t1_roc))[3]

      t1 = caret::confusionMatrix(predict(model1, newdata = temptrainold), as.factor(temptrainold$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_train_Accuracy")] = t1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_train_Sensitivity")] = t1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_train_Specificity")] = t1$byClass["Specificity"]

      v1 = caret::confusionMatrix(predict(model1, newdata = test), as.factor(test$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_test_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_test_Sensitivity")]  = v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_test_Specificity")]  = v1$byClass["Specificity"]

      v1 = caret::confusionMatrix(predict(model1, newdata = valid), as.factor(valid$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_valid_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_valid_Sensitivity")]= v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_valid_Specificity")] = v1$byClass["Specificity"]
      # return(wyniki)
      # }
      # , warning = function(warning_condition) {
      #   print(paste("WARNING:  ",warning_condition))
      #   return(wyniki)
      # }, error = function(error_condition) {
      #   print(paste("WARNING:  ",error_condition))
      #   wyniki[i,paste0(algorytm,"_train_Accuracy")] = as.character(error_condition)
      #   wyniki[i,paste0(algorytm,"_train_Sensitivity")] = NA
      #   wyniki[i,paste0(algorytm,"_train_Specificity")] = NA
      #
      #   wyniki[i,paste0(algorytm,"_test_Accuracy")]  = NA
      #   wyniki[i,paste0(algorytm,"_test_Sensitivity")]  = NA
      #   wyniki[i,paste0(algorytm,"_test_Specificity")]  = NA
      #
      #   wyniki[i,paste0(algorytm,"_valid_Accuracy")]  = NA
      #   wyniki[i,paste0(algorytm,"_valid_Sensitivity")]= NA
      #   wyniki[i,paste0(algorytm,"_valid_Specificity")] = NA
      #   return(wyniki)
      # }, finally={
      #   stargazer(wyniki, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)
      # })
      # wyniki2 = wyniki
      # colnames(wyniki2) = make.names(colnames(wyniki2), unique = TRUE, allow_ = FALSE)
      # stargazer(wyniki2, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)
      write.csv(wyniki,paste0("temp/",output_file))
    }

    stopCluster(cl)

  }

# wyniki2 = wyniki
# colnames(wyniki2) = make.names(colnames(wyniki2), unique = TRUE, allow_ = FALSE)
# stargazer(wyniki2, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)

  write.csv(wyniki,output_file)
  setwd(oldwd)

  sink(type = "message")
  sink()
  dev.off()
  return(wyniki)
}
