
# =======================================================
#  RandomSampling: random oversampling and undersampling
# =======================================================

RandomSampling <-
    function(x, y, percOver = 0, percUnder = 6.8)
        # Inputs
        #       x: A data frame of the predictors from training data
        #       y: A vector of response variable from training data
        #    percOver: Oversampling percentage
        #    percUnder: Undersampling percentage
    {
        if (percUnder > 100)
            stop("percUnder must be less than 100")
        
        # the column where the target variable is
        data <- data.frame(x, y)
        tgt <- length(data)
        classTable <- table(data[, tgt])
        
        # find the minirty and majority class
        minCl <- names(which.min(classTable))
        majCl <- names(which.max(classTable))
        
        # get the cases of the minority and majority class
        indexMin <- which(data[, tgt] == minCl)
        numMin  <- length(indexMin)
        indexMaj <- which(data[, tgt] == majCl)
        numMaj  <- length(indexMaj)
        
        # get the number of instances after sampling
        numMajUnder <- round(percUnder*numMaj/100)
        numMinOver  <- round(percOver*numMin/100)
        if (numMinOver + numMin > numMajUnder)
            warning("More minority instances than majority instances ")
        
        indexMajUnder <- sample(indexMaj, numMajUnder)
        indexMinOver  <- sample(indexMin, numMinOver, replace = TRUE) 
        newIndex  <- c(indexMajUnder, indexMinOver, indexMin)
        data[newIndex, ]
    }






