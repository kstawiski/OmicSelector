# ======================================================================
#  ADASYN: adaptive synthetic sampling
# ======================================================================
# Reference:
# He, H., Y. Bai, et al. (2008). 
# "ADASYN: Adaptive synthetic sampling approach for imbalanced learning."
# IEEE Neural Networks(IJCNN)  
# -----------------------------------------------------------------------

ADASYN <-
    function (x, y, beta = 0.65, k = 5)
        # Inputs:
        #       x : A data frame of the predictors from training data
        #       y : A vector of response variable from training data
        #    beta : Balance level (0, 1], when beta=1, the dataset is fully balanced
        #       k : Number of nearest neighbors
    {
        # find the column of class attribute 
        data <- data.frame(x, y)
        tgt <- length(data)
        classTable <- table(data[, tgt])
        
        # find the minority and majority class labels
        minCl  <- names(which.min(classTable))
        majCl  <- names(which.max(classTable))
        indexMin <- which(data[, tgt] == minCl)
        indexMaj <- which(data[, tgt] == majCl)
        numMin  <- length(indexMin)
        numMaj  <- length(indexMaj)
        numGen   <- (numMaj - numMin)*beta # number of instances to be generated
        
        numRow <- dim(data)[1]
        numCol <- dim(data)[2]
        ratio <- rep(0, numMin)
        
        # move the class attribute to the last column
        if (tgt < numCol)
        {
            cols <- 1:numCol
            cols[c(tgt, numCol)] <- cols[c(numCol, tgt)]
            data <- data[, cols]
        }
        
        # transform factor into integer
        nomatr <- c() 
        dataTransformed <- matrix(nrow = numRow, ncol = numCol-1)
        for (col in 1:(numCol - 1))      
        {
            if (class(data[, col]) == "factor")
            {
                nomatr <- c(nomatr, col)
                dataTransformed[, col] <- as.integer(data[, col])
            }
            else dataTransformed[, col] <- data[, col] 
        } 
        
        # transform dataset into numeric matrix
        source(system.file("extdata", "IRIC/R/Data level/Numeralize.R", package = "OmicSelector"))
        numerMatrix <- Numeralize(data[, -numCol]) 
        
        library("RANN")
        indexOrder <- nn2(numerMatrix, numerMatrix[indexMin, ], k+1)$nn.idx
        
        for (i in 1:numMin)
        {
            kNNs     <- indexOrder[i, 2:(k+1) ]
            numMajNN  <- sum(data[kNNs, numCol] == majCl)
            ratio[i] <- numMajNN/k
        }
        
        if (sum(ratio) == 0)
        {
            stop("no instance will be generated")
        } else {
            ratio <- ratio/sum(ratio)
        }
        newExs  <- matrix(nrow = round(numGen), ncol = numCol-1) #new instance set gernerated
        numCumExs <- 0
        indexOrderMin <- nn2(numerMatrix[indexMin, ], numerMatrix[indexMin, ], k+1)$nn.idx
        
        source(system.file("extdata", "IRIC/R/Data level/SmoteExs.R", package = "OmicSelector"))
        for (i in 1:numMin)
        { 
            numExs <- floor(numGen*ratio[i])
            if (numExs != 0)
            {
                kNNsMin  <- indexOrderMin[i, 2:(k+1)] 
                datakNNsMin <- dataTransformed[indexMin[kNNsMin], ]
                newExs[(numCumExs + 1):(numCumExs + numExs), ] <- InsExs(dataTransformed[indexMin[i], ], datakNNsMin, numExs, nomatr)
                numCumExs <- numCumExs + numExs
            }
        }  
        newExs <- newExs[1:numCumExs, ]
        newExs <- data.frame(newExs)
        
        for (i in nomatr)
        {
            newExs[, i] <- factor(newExs[, i], levels=1:nlevels(data[, i]), labels = levels(data[, i]))
        }
        newExs[, numCol]   <- factor(rep(minCl, numCumExs), levels = levels(data[, numCol]))
        colnames(newExs)  <- colnames(data)
        
        # rearrange of variable structure
        if (tgt < numCol) 
        {
            newExs <- newExs[, cols]
            data   <- data[, cols]
        }
        
        # produce final data set 
        newData <- rbind(data, newExs)
        return(newData)
    }