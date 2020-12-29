# =============================================================================
#  CLUS: clustering-based undersampling method
# =============================================================================
# Yen, S.-J. and Y.-S. Lee (2009). 
# "Cluster-based under-sampling approaches for imbalanced data distributions."
# Expert Systems with Applications 36(3): 5718-5727.
#-----------------------------------------------------------------------------

CLUS <-
    function(x, y, k = 3, m = 1.5)
        # Inputs:
        #       x : A data frame of the predictors from training data
        #       y : A vector of response variable from training data
        #       k : Number of clusters
        #       m : Imbalanced ratio in output dataset
    {
        # find the majority and minority instances
        data <- data.frame(x, y)
        tgt <- length(data)
        classTable <- table(data[, tgt])
        numRow <- dim(data)[1]
        
        # find the minirty and majority
        minCl    <- names(which.min(classTable))
        indexMin <- which(data[, tgt] == minCl)
        numMin   <- length(indexMin)
        majCl    <- names(which.max(classTable))
        #indexMaj <- which(data[, tgt] == majCl)
        numMajFinal <- m*numMin
        
        source(system.file("extdata", "IRIC/R/Data level/Numeralize.R", package = "OmicSelector"))
        dataX            <- Numeralize(data[, -tgt])
        clusteringMoldel <- kmeans(dataX, k)
        mebership        <- clusteringMoldel$cluster
        indexGrouping    <- split(1:numRow, mebership)
        ratio         <- rep(0, k)
        MajCluster   <- list()
        numMajcluster <- rep(0, k)
        
        for (i in 1:k)
        {
            indexMajLocal    <- which(data[indexGrouping[[i]], tgt] == majCl) 
            indexMajGlobal   <- indexGrouping[[i]][indexMajLocal]
            numMincluster    <- sum(data[indexGrouping[[i]], tgt] == minCl)
            MajCluster[[i]]  <- indexMajGlobal
            numMajcluster[i] <- length(indexMajGlobal)
            if (numMincluster == 0) 
                numMincluster <- 1
            ratio[i] <- numMajcluster[i]/numMincluster
        }
        
        #control the imbalance ratio in the output datasets
        ratio <- ratio/sum(ratio)
        indexMajFinal <- c()
        for (i in 1:k)
        {
            if  (ratio[i]!=0) 
            {
                numMajUnder   <- round(numMajFinal*ratio[i])
                indexSelection    <- sample(1:numMajcluster[i], numMajUnder, replace = TRUE)
                indexMajSelection <- MajCluster[[i]][indexSelection]
                indexMajFinal     <- c(indexMajFinal, indexMajSelection)
            }
        }
        
        newData   <- rbind(data[indexMin, ], data[indexMajFinal, ])
        return(newData)
        
    }






