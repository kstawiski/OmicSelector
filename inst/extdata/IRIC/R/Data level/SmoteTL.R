#Copyright (C) 2018 Bing Zhu
# ==============================================
#  SmoteTL: Smote sampling+TomekLinks
# ==============================================

SmoteTL <-
    function(x, y, percOver = 1400, k = 5)
        # Inputs
        #      x    : A data frame of the predictors from training data
        #      y    : A vector of response variable from training data
        #   per_over: Number of new instance generated for each minority instance
        #   k       : Number of nearest neighbors used in Smote
    { 
        source(system.file("extdata", "IRIC/R/Data level/SMOTE.R", package = "OmicSelector"))
        newData <- SMOTE(x, y, percOver, k)
        tgt <- length(newData)
        indexTL <- TomekLink(tgt, newData)
        newDataRemoved <- newData[!indexTL, ]
        return(newDataRemoved)
    }


# ==========================================
#  TomekLink: find the TomekLink
# ==========================================

TomekLink <-
    function(tgt, data)
        # Inputs:
        #   form: model formula
        #   data: dataset 
        # Output:
        #   logical vector indicating whether a instance is in TomekLinks
    {
        
        indexTomek <- rep(FALSE, nrow(data))
        
        # find the column of class variable
        classTable <- table(data[, tgt])
        
        # seperate the group
        majCl <- names(which.max(classTable))
        minCl <- names(which.min(classTable))
        
        # get the instances of the larger group
        indexMin <- which(data[, tgt] == minCl)
        #numMin  <- length(indexMin)
        
        
        # convert dataset in numeric matrix
        source(system.file("extdata", "IRIC/R/Data level/Numeralize.R", package = "OmicSelector"))
        dataTransformed <- Numeralize(data[, -tgt])
        
        # generate indicator matrix 
        require("RANN")
        indexOrder1  <- nn2(dataTransformed, dataTransformed[indexMin, ], k = 2)$nn.idx
        indexTomekCa <- data[indexOrder1[, 2], tgt] == majCl
        if (sum(indexTomekCa) > 0)
        {
            TomekCa <- cbind(indexMin[indexTomekCa],indexOrder1[indexTomekCa, 2])
            
            # find nearest neighbour of potential majority instance
            indexOrder2 <- nn2(dataTransformed, dataTransformed[TomekCa[, 2], ], k = 2)$nn.idx
            indexPaired <- indexOrder2[ ,2] == TomekCa[, 1]
            if (sum(indexPaired) > 0)
            {
                indexTomek[TomekCa[indexPaired, 1]] <- TRUE  
                indexTomek[TomekCa[indexPaired, 2]] <- TRUE  
            }
        }
        return(indexTomek)
    }







