#' OmicSelector_xgboost
#'
#' Train xgboost model with Bayesian optimalization.
#' Code based on http://www.mysmu.edu/faculty/jwwang/post/hyperparameters-tuning-for-xgboost-using-bayesian-optimization/
#'
#' @param features Vector of features to be used. If "all", all features starting with 'hsa' will be used.
#' @param train Training dataset with column Class ('Case' vs. 'Control') and features starting with 'hsa'.
#' @param test Testing dataset with column Class ('Case' vs. 'Control') and features starting with 'hsa'.
#' @param valid Testing dataset with column Class ('Case' vs. 'Control') and features starting with 'hsa'.
#' @param eta Bonderies of 'eta' parameter in XGBoost training, must be a vector of 2.
#' @param gamma Bonderies of 'gamma' parameter in XGBoost training, must be a vector of 2.
#' @param max_depth Bonderies of 'max_depth' parameter in XGBoost training, must be a vector of 2.
#' @param min_child_weight Bonderies of 'min_child_weight' parameter in XGBoost training, must be a vector of 2.
#' @param subsample Bonderies of 'subsample' parameter in XGBoost training, must be a vector of 2.
#' @param nfold Bonderies of 'nfold' parameter in XGBoost training, must be a vector of 2.
#' 
#' @return Xgboost model
#'
#' @export
OmicSelector_xgboost = function(features = "all", train = OmicSelector_load_datamix(use_smote_not_rose = T)[[1]], test = OmicSelector_load_datamix(use_smote_not_rose = T)[[2]], valid = OmicSelector_load_datamix(use_smote_not_rose = T)[[2]],
                                eta = c(0, 1),
                    gamma =c(0, 100),
                    max_depth = c(2L, 10L), # L means integers
                    min_child_weight = c(1, 25),
                    subsample = c(0.25, 1),
                    nfold = c(3L, 10L), initPoints = 8, iters.n=10){
    library(xgboost)
    library(ParBayesianOptimization)
    library(mlbench) 
    library(dplyr)
    library(recipes)
    library(resample)
    library(ROCR)
    library(ggplot2)

    tempwyniki = data.frame()
    tempwyniki[1, "method"] = "xgboost"
    tempwyniki[1, "features"] = paste0(features, collapse = " + ")
    tempwyniki[1, "eta"] = paste0(eta, collapse = " , ")
    tempwyniki[1, "gamma"] = paste0(gamma, collapse = " , ")
    tempwyniki[1, "max_depth"] = paste0(max_depth, collapse = " , ")
    tempwyniki[1, "min_child_weight"] = paste0(min_child_weight, collapse = " , ")
    tempwyniki[1, "subsample"] = paste0(subsample, collapse = " , ")
    tempwyniki[1, "nfold"] = paste0(nfold, collapse = " , ")
    tempwyniki[1, "initPoints"] = paste0(initPoints, collapse = " , ")
    tempwyniki[1, "iters.n"] = paste0(iters.n, collapse = " , ")

    if(features != "all")
    {
        library(dplyr)
        train = dplyr::select(train, Class, all_of(features))
        test = dplyr::select(test, Class, all_of(features))
        valid = dplyr::select(valid, Class, all_of(features))
    } else {
        train = dplyr::select(train, Class, starts_with('hsa'))
        test = dplyr::select(test, Class, starts_with('hsa'))
        valid = dplyr::select(valid, Class, starts_with('hsa'))
    }
    
    eq = as.formula("Class ~ . - 1")
    # Make matrices for training data
    #x <- model.matrix(eq, data = train)
    #y <- ifelse(train$Class == "Case", 1, 0)
    x <- model.matrix(eq, data = train)
    y <- ifelse(train$Class == "Case", 1, 0)
    assign("x", x, envir = .GlobalEnv)
    assign("y", y, envir = .GlobalEnv)

    if(sum(train$Class == "Control") / sum(train$Class == "Case")>1) {
    scale_pos_weight = c(0,sum(train$Class == "Control") / sum(train$Class == "Case")) }
    else { scale_pos_weight = c(sum(train$Class == "Control") / sum(train$Class == "Case"),5) }

    # Make matrices for testing data
    xvals <- model.matrix(eq, data = test)
    yvals <- ifelse(test$Class == "Case", 1, 0)

    # Make matrices for valid data
    xvals2 <- model.matrix(eq, data = valid)
    yvals2 <- ifelse(valid$Class == "Case", 1, 0)

    xgb_train = xgb.DMatrix(data = x, label = y)
    xgb_test = xgb.DMatrix(data = xvals, label = yvals)
    xgb_valid = xgb.DMatrix(data = xvals2, label = yvals2)

    #watchlist = list(train=xgb_train, test=xgb_test)
    #model = xgb.train(data = xgb_train, max.depth = 10, watchlist=watchlist, nrounds = 100, params = 
    #              list(booster = "gbtree", objective = "binary:logistic", eval_metric = "auc", scale_pos_weight = scale_pos))

    bounds <- list(
                    eta = eta,
                    gamma =gamma,
                    max_depth =max_depth, # L means integers
                    min_child_weight = min_child_weight,
                    subsample = subsample,
                    nfold = nfold,
                    scale_pos_weight = scale_pos_weight
                    )
    
    set.seed(2021)

    library(doParallel)
    no_cores <- detectCores() / 2 # use half my CPU cores to avoid crash  
    cl <- makeCluster(no_cores) # make a cluster
    registerDoParallel(cl) # register a parallel backend
    clusterExport(cl, c('x','y')) # import objects outside
    clusterEvalQ(cl,expr= { # launch library to be used in FUN
    library(xgboost)
    library(OmicSelector)
        
    OmicSelector_xgboost_scoring_function <- function(eta, gamma, max_depth, min_child_weight, subsample, nfold) {

    library(xgboost)
    library(ParBayesianOptimization)
    library(mlbench) 
    library(dplyr)
    library(recipes)
    library(resample)
    library(ROCR)
    library(ggplot2)
  dtrain <- xgb.DMatrix(x, label = y, missing = NA)

  pars <- list(
    eta = eta,
    gamma = gamma,
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    subsample = subsample,
    scale_pos_weight = scale_pos_weight,

    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "auc",
    verbosity = 0
  )
  
  xgbcv <- xgb.cv(
    params = pars,
    data = dtrain,
    
    nfold = nfold,
    # scale_pos_weight = scale_pos,
    nrounds = 1000,
    prediction = TRUE,
    showsd = TRUE,
    early_stopping_rounds = 50,
    maximize = TRUE,
    stratified = TRUE
  )
  
  # required by the package, the output must be a list
  # with at least one element of "Score", the measure to optimize
  # Score must start with capital S
  # For this case, we also report the best num of iteration
  return(
    list(
      Score = max(xgbcv$evaluation_log$test_auc_mean),
      nrounds = xgbcv$best_iteration
    )
  )
}
    
    })

    time_withparallel <- system.time(
    opt_obj <- bayesOpt(
    FUN = OmicSelector_xgboost_scoring_function,
    bounds = bounds,
    initPoints = initPoints,
    iters.n = iters.n,
    parallel = TRUE
    ))

    stopCluster(cl) # stop the cluster
    registerDoSEQ() # back to serial computing
    
    # take the optimal parameters for xgboost()
    params <- list(eta = getBestPars(opt_obj)[1],
                gamma = getBestPars(opt_obj)[2],
                max_depth = getBestPars(opt_obj)[3],
                min_child_weight = getBestPars(opt_obj)[4],
                subsample = getBestPars(opt_obj)[5],
                nfold = getBestPars(opt_obj)[6],
                objective = "binary:logistic",
                scale_pos_weight = getBestPars(opt_obj)[7])

    # the numrounds which gives the max Score (auc)
    numrounds <- opt_obj$scoreSummary$nrounds[
    which(opt_obj$scoreSummary$Score
            == max(opt_obj$scoreSummary$Score))]

    fit_tuned <- xgboost(params = params,
                        data = x,
                        label = y,
                        nrounds = numrounds,
                        eval_metric = "auc")


    y_train_pred <- predict(fit_tuned, x, type = "response")
    y_test_pred <- predict(fit_tuned, xvals, type = "response")
    y_valid_pred <- predict(fit_tuned, xvals2, type = "response")

    pred = data.frame(`Class` = train$Class, `Pred.2` = y_train_pred)
    library(cutpointr)
    cutoff = cutpointr(pred, Pred.2, Class, pos_class = "Case", metric = youden)
    wybrany_cutoff = cutoff$optimal_cutpoint

    tempwyniki[1, "cutoff"] = wybrany_cutoff

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

    t1_roc = pROC::roc(Class ~ as.numeric(Pred.2), data=pred)
    tempwyniki[1, "training_AUC2"] = t1_roc$auc
    tempwyniki[1, "training_AUC_lower95CI"] = as.character(ci(t1_roc))[1]
    tempwyniki[1, "training_AUC_upper95CI"] = as.character(ci(t1_roc))[3]
    #message("Checkpoint passed: chunk 15")

    tempwyniki[1, "training_Accuracy"] = cm_train$overall[1]
    tempwyniki[1, "training_Sensitivity"] = cm_train$byClass[1]
    tempwyniki[1, "training_Specificity"] = cm_train$byClass[2]
    tempwyniki[1, "training_PPV"] = cm_train$byClass[3]
    tempwyniki[1, "training_NPV"] = cm_train$byClass[4]
    tempwyniki[1, "training_F1"] = cm_train$byClass[7]
    #message("Checkpoint passed: chunk 16")

    cat(paste0("\n\n---- TESTING PERFORMANCE ----\n\n"))
    pred = data.frame(`Class` = test$Class, `Pred.2` = y_test_pred)
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
    #message("Checkpoint passed: chunk 17")

    cat(paste0("\n\n---- VALIDATION PERFORMANCE ----\n\n"))
    pred = data.frame(`Class` = valid$Class, `Pred.2` = y_valid_pred)
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

    wyniki_modelu = list("model" = fit_tuned, "features" = features, "params" = params, "numround" = numrounds, "cutoff" = cutoff, "selected_cutoff" = wybrany_cutoff, "metric" = tempwyniki, "cm_train" = cm_train, "training_roc" = t1_roc, "cm_test" = cm_test, "cm_valid" = cm_valid, "train" = train, "test" = test, "valid" = valid)
    return(wyniki_modelu)
}


OmicSelector_xgboost_scoring_function <- function(eta, gamma, max_depth, min_child_weight, subsample, nfold, scale_pos_weight) {

  dtrain <- xgb.DMatrix(x, label = y, missing = NA)

  pars <- list(
    eta = eta,
    gamma = gamma,
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    subsample = subsample,
    scale_pos_weight = scale_pos_weight,

    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "auc",
    verbosity = 0
  )
  
  xgbcv <- xgb.cv(
    params = pars,
    data = dtrain,
    
    nfold = nfold,
    # scale_pos_weight = scale_pos,
    nrounds = 1000,
    prediction = TRUE,
    showsd = TRUE,
    early_stopping_rounds = 50,
    maximize = TRUE,
    stratified = TRUE
  )
  
  # required by the package, the output must be a list
  # with at least one element of "Score", the measure to optimize
  # Score must start with capital S
  # For this case, we also report the best num of iteration
  return(
    list(
      Score = max(xgbcv$evaluation_log$test_auc_mean),
      nrounds = xgbcv$best_iteration
    )
  )
}