#' OmicSelector_iteratedRFE
#'
#' Helper for Marcin Kaszkowiak propritary method.
#' Very longs.. needs optimizing.
#'
#' @export
OmicSelector_iteratedRFE <- function(trainSet, testSet = NULL, initFeatures = colnames(trainSet), classLab, checkNFeatures = 25, votingIterations = 1000, useCV = F, nfolds = 10, initRandomState = 42 ) {

  set.seed(initRandomState)
  initFeatures <- initFeatures[initFeatures != classLab]
  #prepare output data structures
  resAcc <- c(rep(0, checkNFeatures))
  resVotes <- data.frame(matrix(0, nrow = length(initFeatures), ncol = checkNFeatures), row.names = initFeatures)
  for(i in 1:checkNFeatures) colnames(resVotes)[i] <- toString(i)
  resTop <- list()


  for (i in 1:votingIterations) {

    print(paste0("iteration ", i))

    if(useCV == F) {
      params <- rfeControl(functions = rfFuncs, saveDetails = T)
      iter <- rfeIter(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), testX = testSet[, initFeatures], testY = as.factor(testSet[, classLab]), sizes = 1:checkNFeatures,
                      metric = "Accuracy", rfeControl = params)

      for(j in 1:checkNFeatures) {
        tmp <- iter$pred[iter$pred$Variables == j, ]

        acc <- length(which(tmp$pred == tmp$obs)) / nrow(tmp) #calculate and add accuracy
        resAcc[j] <- resAcc[j] + acc

        selected <- iter$finalVariables[[j+1]]$var
        numb <- iter$finalVariables[[j+1 ]]$Variables[1]

        resVotes[selected, numb] <- resVotes[selected, numb] + 1
      }


    }
    else {

      seeds <- vector(mode = "list", length = nfolds + 1) # add random seeds for cross validation
      for(i in 1:nfolds) seeds[[i]] <- sample.int(1000000000, checkNFeatures + 1)
      seeds[nfolds + 1] <- sample.int(1000000000, 1)

      params <- rfeControl(functions = rfFuncs, number = nfolds, saveDetails = T)
      iter <- rfe(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), sizes = 1:checkNFeatures, rfeControl = params)

      for(j in 1:checkNFeatures) {
        tmp <- iter$variables[iter$variables$Variables == j, ]
        for(k in tmp$var) resVotes[k, j] <- resVotes[k, j] + 1 # increase a voting score for each fold

        resAcc[j] <- resAcc[j] + iter$results[iter$results$Variables == j, "Accuracy"]

      }
    }
  }

  resAcc <- resAcc / votingIterations #make average accuracy

  for(i in 1:ncol(resVotes)) resTop[[i]] <- rownames(resVotes[order(-resVotes[, i])[1:i], ])

  returning <- list(data.frame(resAcc), resVotes, resTop)
  names(returning) <- c("accuracyPerNFeatures", "votesPerN", "topFeaturesPerN")
  return(returning)

}
