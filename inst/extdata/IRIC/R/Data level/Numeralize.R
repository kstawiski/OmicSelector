# ======================================================
# Numeralize: convert dataset into numeric matrix
# ======================================================

Numeralize <-
  function(data, form = NULL)
  { 
    if (!is.null(form))
    {
      tgt    <- which(names(data) == as.character(form[[2]]))
      dataY <- data[drop = FALSE,, tgt]
      dataX <- data[, -tgt]
    } else {
      dataX <- data
    }
    numRow      <- dim(dataX)[1]
    #numCol      <- dim(dataX)[2]
    indexOrder      <- sapply(dataX, is.ordered)
    indexMultiValue <- sapply(dataX, nlevels)>2
    indexNominal    <- !indexOrder & indexMultiValue
    numerMatrixNames<- NULL
    if (all(indexNominal))
    {
      numerMatrix   <- NULL
    } else {
      numerMatrix      <- dataX[drop = FALSE, ,!indexNominal]
      numerMatrixNames <- colnames(numerMatrix)
      numerMatrix      <- data.matrix(numerMatrix)
      Min              <- apply(numerMatrix, 2, min) 
      range            <- apply(numerMatrix, 2, max)-Min
      numerMatrix      <- scale(numerMatrix, Min, range)[, ]
    }
    
    if (any(indexNominal))
    {
      
      BiNames     <- NULL
      dataNominal <- dataX[drop = FALSE, ,indexNominal]
      numNominal  <- sum(indexNominal)
      if (numNominal>1)
      {
        dimEx <- sum(sapply(dataX[,indexNominal], nlevels))
      } else {
        dimEx <- nlevels(dataX[, indexNominal])
      }
      dataBinary  <- matrix(nrow = numRow, ncol = dimEx )
      cl <- 0
      for (i in 1:numNominal)     
      { 
        numCat <- nlevels(dataNominal[, i])
        for (j in 1:numCat)
        {     
          value <- levels(dataNominal[, i])[j]
          ind  <- (dataNominal[,i] == value)
          dataBinary[, cl+1] <- as.integer(ind)
          BiNames[cl+1]   <- paste(names(dataNominal)[i], "_", value, sep="")
          cl <- cl+1
        }
      }
      numerMatrix  <- cbind(numerMatrix, dataBinary)
      colnames(numerMatrix) <- c(numerMatrixNames, BiNames)
    }  
    
    if (!is.null(form))
    {
      numerMatrix <- data.frame(numerMatrix)
      numerMatrix <- cbind(numerMatrix, dataY)
    }
    return(numerMatrix)
  }
