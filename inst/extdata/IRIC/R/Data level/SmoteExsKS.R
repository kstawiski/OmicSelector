#=========================================================
# SmoteExs: obtain Smote instances for minority instances
#=========================================================

SmoteExsKS<-
  function(data, howManyNew, k)
    # Input:
    #     data      : dataset of the minority instances
    #     howManyNew   : how many new cases - oversampling 
    #     k         : number of nearest neighours
    
  {
    # transform factors into integer
    nomAtt  <- c()
    numRow  <- dim(data)[1]
    numCol  <- dim(data)[2]
    dataX   <- data[ ,-numCol]
    dataTransformed <- matrix(nrow = numRow, ncol = numCol-1)
    for (col in 1:(numCol-1))
    { 
      if (is.factor(data[, col])) 
      {
        dataTransformed[, col] <- as.integer(data[, col])
        nomAtt <- c(nomAtt , col)
      } else {
        dataTransformed[, col] <- data[, col]
      }
    }
    numExs  <-  howManyNew # this is the number of artificial instances generated
    newExs  <-  matrix(ncol = numCol-1, nrow = numRow*numExs)
    
    indexDiff <- sapply(dataX, function(x) length(unique(x)) > 1)
    source(system.file("extdata", "IRIC/R/Data level/Numeralize.R", package = "OmicSelector"))
    #system.file("extdata", "IRIC/R/Data level/Numeralize.R")
    numerMatrix <- Numeralize(dataX[ ,indexDiff])
    require("RANN")
    id_order <- nn2(numerMatrix, numerMatrix, k+1)$nn.idx
    for(i in 1:numRow) 
    {
      kNNs   <- id_order[i, 2:(k+1)]
      newIns <- InsExs(dataTransformed[i, ], dataTransformed[kNNs, ], numExs, nomAtt)
      newExs[((i-1)*numExs+1):(i*numExs), ] <- newIns
    }
    
    # get factors as in the original data.
    newExs <- data.frame(newExs)
    for(i in nomAtt)
    {
      newExs[, i] <- factor(newExs[, i], levels = 1:nlevels(data[, i]), labels = levels(data[, i]))
    }
    newExs[, numCol] <- factor(rep(data[1, numCol], nrow(newExs)), levels=levels(data[, numCol]))
    colnames(newExs) <- colnames(data)
    return(newExs)
  }

#=================================================================
# InsExs: generate Synthetic instances from nearest neighborhood
#=================================================================

InsExs <-
  function(instance, dataknns, numExs, nomAtt)
    # Input:
    #    instance : selected instance
    #    dataknns : nearest instance set
    #    numExs   : number of new intances generated for each instance
    #    nomAtt   : indicators of factor variables 
  {
    numRow  <- dim(dataknns)[1]
    numCol  <- dim(dataknns)[2]
    newIns <- matrix (nrow = numExs, ncol = numCol)
    neig   <- sample(1:numRow, size = numExs, replace = TRUE)
    
    # generated  attribute values 
    insRep  <- matrix(rep(instance, numExs), nrow = numExs, byrow = TRUE)
    diffs   <- dataknns[neig,] - insRep
    newIns  <- insRep + runif(1)*diffs
    # randomly change nominal attribute
    for (j in nomAtt)
    {
      newIns[, j]   <- dataknns[neig, j]
      indexChange   <- runif(numExs) < 0.5
      newIns[indexChange, j] <- insRep[indexChange, j]
    }
    return(newIns) 
  }
