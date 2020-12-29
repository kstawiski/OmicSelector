#Copyright (C) 2018 Bing Zhu
#================================================================
# BalanceCascade
#========================================================================
# Reference
# Xu-Ying, L., W. Jianxin, et al. (2009). "Exploratory # Undersampling 
# for Class-Imbalance Learning."  Systems, Man, and Cybernetics, Part B: 
# Cybernetics, IEEE Transactions on 39(2): 539-550.
#========================================================================
BalanceCascade <-
    function(x, ...)
        UseMethod("BalanceCascade")


BalanceCascade.data.frame  <-
    function (x, y, iter = 4)
    {
        # Input:
        #        x: A data frame of the predictors from training data
        #        y: A vector of response variable from training data
        #     iter: Iterations to train base classifiers
        # allowParallel: A logical number to control the parallel computing. If allowParallel =TRUE, the function is run using parallel techniques
        
        source(system.file("extdata", "IRIC/R/Ensemble-based level/BalanceBoost.R", package = "OmicSelector"))
        funcCall <- match.call(expand.dots = FALSE)
        data <- data.frame(x, y)
        tgt <- length(data)
        classTable   <- table(data[, tgt])
        classTable   <- sort(classTable, decreasing = TRUE)
        classLabels  <- names(classTable)
        indexMaj <- which(data[, tgt] == classLabels[1])
        indexMin <- which(data[, tgt] == classLabels[2])
        numMin <- length(indexMin)
        numMaj <- length(indexMaj)
        FP <- (numMin/numMaj)^(1/(iter-1))
        
        #initialization
        x.nam <- names(x)
        form <- as.formula(paste("y ~ ", paste(x.nam, collapse = "+")))
        H      <- list()
        thresh <- rep(NA, iter)
        
        for (i in seq(iter)){
            if (length(indexMaj) < numMin)
                numMin  <- length(indexMaj)
            indexMajSampling <- sample(indexMaj, numMin)
            dataCurrent <- data[c(indexMin, indexMajSampling),]      
            H[[i]] <- bboost.data.frame(dataCurrent[, -tgt], dataCurrent[,tgt], type = "AdaBoost")
            pred   <- predict(H[[i]], data[c(indexMaj), -tgt], data[c(indexMaj), tgt], type ="probability") 
            sortIndex   <- order(pred[, 2], decreasing = TRUE)
            numkeep     <- round(length(indexMaj)*FP)
            thresh[i]   <- pred[sortIndex[numkeep],2]*sum(H[[i]]$alpha)   
            indexMaj    <- indexMaj[sortIndex[1:numkeep]]    
        }
        
        iter   <- sum(sapply(H,"[[", 5))
        fits   <- unlist(lapply(H,"[[", 6), recursive = FALSE) 
        alphas <- unlist(lapply(H,"[[", 7))
        
        structure(
            list( call        = funcCall   ,
                  iter        = iter       ,
                  classLabels = classLabels, 
                  base        = H[[1]]$base,
                  alphas      = alphas      ,
                  fits        = fits       ,
                  thresh      = sum(thresh))  ,
            class = "BalanceCascade")
        
    }

predict.BalanceCascade<-    
    function(obj, x,  type = "class")
    {
        
        #  input 
        #     obj: Output from BalanceCascade.data.frame
        #       x: A data frame of the predictors from testing data
      
        if(is.null(x)) stop("please provide predictors for prediction")  
        if (!type %in% c("class", "probability"))
            stop("wrong setting with type")
        data <- x
        classLabels <- obj$classLabels
        numClass    <- length(classLabels)
        numIns      <- dim(data)[1]
        weight      <- obj$alphas
        btPred      <- sapply(obj$fits, obj$base$pred, data = data, type ="class")    
        classfinal  <- matrix(0, ncol = numClass, nrow = numIns)
        colnames(classfinal) <- classLabels
        for (i in 1:numClass){
            classfinal[, i] <- matrix(as.numeric(btPred == classLabels[i]), nrow = numIns)%*%weight
        }
        if (type == "class")
        {
            classfinal <- classfinal - obj$thresh
            out <- factor(classLabels[apply(classfinal, 1, which.max)], levels = classLabels)
            
        } else {
            out <- data.frame(classfinal/rowSums(classfinal))  
            
        }
        out
    }








