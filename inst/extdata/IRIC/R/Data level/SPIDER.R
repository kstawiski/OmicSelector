#Copyright (C) 2018 Bing Zhu
# ============================================================
#  SPIDER: selective pre-processing
# ============================================================
# Reference:
# Stefanowski, J. and S. Wilk (2008). 
# "Selective pre-processing of imbalanZced data for 
# improving classification performance."
# Data Warehousing and Knowledge Discovery, Springer: 283-292.
# ------------------------------------------------------------

SPIDER <-
    function(x, y, method = "weak", allowParallel = TRUE)
        # Inputs:
        #       X   : A data frame of the predictors from training data
        #       y   : A vector of response variable from training data
        #    method : Type of modification of the minority class in the second phase, including ??weak??, ??relabel??, ??strong??
    {
        data <- data.frame(x, y)
        numRow <- dim(data)[1]
        numCol <- dim(data)[2]
        
        # find class attribute
        tgt <- length(data)
        classTable <- table(data[, tgt])
        
        # find the minority instances
        minCl    <- names(which.min(classTable))
        indexMin <- which(data[, tgt] == minCl)
        numMin    <- length(indexMin)
        
        # find the majority instances
        majCl    <- names(which.max(classTable))
        indexMaj <- which(data[, tgt] == majCl)
        numMaj   <- length(indexMaj)
        #dataMaj <- data[indexMaj, ]
        indexSafe  <- rep(FALSE, numRow)
        dataCopy <- data.frame(matrix(ncol = numCol, nrow = numMin*5))
        numCopy    <- 0
        
        source(system.file("extdata", "IRIC/R/Data level/Numeralize.R", package = "OmicSelector"))
        dataTransformed <- Numeralize(data[, -tgt])
        
        library("RANN")
        indexOrder <- nn2(dataTransformed, dataTransformed, 6)$nn.idx
        indexOrderMin <- indexOrder[indexMin, ] 
        nomatr <- which(sapply(data, class) == "factor")
        # find the neigborhood and safty flag
        if(allowParallel) {
            library("parallel")
            cl <- makeCluster(2)
            classMatrix <- matrix(data[indexOrder[, 2:4], tgt], nrow = numRow)
            classCount  <- parApply (cl, classMatrix, 1, table)
            modeColumn  <- parLapply(cl, classCount, which.max)
            classMode   <- parSapply(cl, modeColumn, names)
            indexSafe[data[, tgt]== classMode] <- TRUE
            stopCluster(cl)
            indexMajLeft <- indexMaj[indexSafe[indexMaj]]
        } else {
            for (i in 1:numRow)      
            {    
                classCount <- table(data[indexOrder[i, (2:4)], tgt])
                classMode <- names(which.max(classtable))
                if (data[i, tgt] == classMode)
                    indexSafe[i] <- TRUE
            }
            indexMajLeft <- indexMaj[indexSafe[indexMaj]]
        }
        
        if (method == "relabel")
        {
            dataRelabel <- data.frame(matrix(ncol = numCol, nrow = numMaj))
            colnames(dataRelabel) <- colnames(data)
            numRelabel <- 0
        }
        newCopy <- NULL
        for (i in 1:numMin)
        { 
            # weak method
            if (method == "weak" && indexSafe[indexMin[i]] == FALSE)
            {     
                newCopy <- Copying(data[indexMin[i],  ], data[indexOrderMin[i, 2:4], ], indexSafe[indexOrderMin[i, 2:4]], majCl, tgt)
            }
            
            # relabel method
            if (method == "relabel" && indexSafe[indexMin[i]] == FALSE)
            {
                newCopy   <- Copying(data[indexMin[i], ], data[indexOrderMin[i, 2:4], ], indexSafe[indexOrderMin[i, 2:4]], majCl, tgt)
                kNNsMin    <- indexOrderMin[i, 2:4]
                kNNsMaj    <- which(data[kNNsMin, tgt] == majCl)
                kNNsNoise  <- which(indexSafe[kNNsMin] == FALSE)
                kNNsMajNoise <- intersect(kNNsMaj, kNNsNoise)
                numRelaCurrent <- length(kNNsMajNoise)
                if (numRelaCurrent > 0)
                {
                    dataRelabel[(numRelabel+1):(numRelabel+ numRelaCurrent), ] <- data[kNNsMin[kNNsMajNoise], ]
                    numRelabel    <- numRelabel + numRelaCurrent
                }
            }
            
            # strong method
            if (method == "strong" && indexSafe[indexMin[i]] == TRUE)
            {     
                newCopy <- Copying(data[indexMin[i], ], data[indexOrderMin[i, 2:4], ], indexSafe[indexOrderMin[i, 2:4]], majCl, tgt)
            }
            if (method == "strong" && indexSafe[indexMin[i]] == FALSE)
            {              
                classCount <- table(data[indexOrderMin[i, 2:6], tgt])
                classMode <- names(which.max(classCount))
                if (classMode == minCl)
                {
                    newCopy <- Copying(data[indexMin[i], ], data[indexOrderMin[i, 2:4], ], indexSafe[indexOrderMin[i, 2:4]], majCl, tgt)
                } else {
                    newCopy <- Copying(data[indexMin[i], ], data[indexOrderMin[i, 2:6], ], indexSafe[indexOrderMin[i, 2:6]], majCl, tgt)
                }
            }
            
            # add new copies
            if (!is.null(newCopy))
            { 
                numNewCopy <- nrow(newCopy)
                dataCopy[(numCopy + 1):(numCopy + numNewCopy),] <- newCopy
                numCopy <- numCopy + numNewCopy
            }
        }
        
        if (numCopy > 0)
        {
            dataCopy <- dataCopy[1:numCopy, ]
            colnames(dataCopy) <- colnames(data)
        } else {
            dataCopy <- NULL
        }
        
        if (method == "relabel")
        {    
            if (numRelabel > 0)
                dataRelabel <- dataRelabel[1:numRelabel, , drop = FALSE ]
        }
        
        for (i in nomatr)
        {
            dataCopy[, i] <- factor(dataCopy[, i], levels = 1:nlevels(data[, i]), labels = levels(data[, i]))
            if (method == "relabel")
                dataRelabel[, i] <- factor(dataRelabel[, i], levels = 1:nlevels(data[, i]), labels = levels(data[, i]))
        }
        
        if (method == "relabel")
        {
            dataRelabel[, tgt] <- data[indexMin[1], tgt]
            newData <- rbind(data[indexMajLeft, ],  data[indexMin, ], dataCopy, dataRelabel)
        } else {
            newData <- rbind(data[indexMajLeft, ],  data[indexMin, ], dataCopy)
        }
        return(newData) 
    }

# ===================================
# Copying: Duplicate(copy) instances 
# ===================================

Copying<-
    function(instance, datakNNs, safe, majCl, tgt)
    { 
        kNNsMaj     <- which(datakNNs[, tgt] == majCl)
        kNNsSafe    <- which(safe == TRUE)
        kNNsMajSafe  <- intersect(kNNsMaj, kNNsSafe)
        if (length(kNNsMajSafe) > 0)
        {
            index   <- rep(1, length(kNNsMajSafe))
            newCopy <- instance[index, ]
        } else {
            newCopy <- NULL
        }
        return(newCopy)
    }



