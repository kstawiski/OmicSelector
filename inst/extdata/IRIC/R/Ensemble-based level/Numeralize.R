# ======================================================
# Numeralize: convert dataset into numeric matrix
# ======================================================

Numeralize <-
  function(data, form = NULL)
  { 
    if (!is.null(form))
    {
      tgt    <- which(names(data) == as.character(form[[2]]))
      data_Y <- data[drop = FALSE,, tgt]
      data_X <- data[, -tgt]
    } else {
      data_X <- data
    }
    n_row      <- dim(data_X)[1]
    n_col      <- dim(data_X)[2]
    id_order      <- sapply(data_X, is.ordered)
    id_multiValue <- sapply(data_X, nlevels)>2
    id_nominal    <- !id_order & id_multiValue
    numerMatrixnames<- NULL
    if (all(id_nominal))
    {
      numerMatrix   <- NULL
    } else {
      numerMatrix      <- data_X[drop = FALSE, ,!id_nominal]
      numerMatrixnames <- colnames(numerMatrix)
      numerMatrix      <- data.matrix(numerMatrix)
      Min              <- apply(numerMatrix, 2, min) 
      range            <- apply(numerMatrix, 2, max)-Min
      numerMatrix      <- scale(numerMatrix, Min, range)[, ]
    }
    
    if (any(id_nominal))
    {
      
      BiNames      <- NULL
      data_nominal <- data_X[drop = FALSE, ,id_nominal]
      n_nominal    <- sum(id_nominal)
      if (n_nominal>1)
      {
        dimEx <- sum(sapply(data_X[,id_nominal], nlevels))
      } else {
        dimEx <- nlevels(data_X[, id_nominal])
      }
      data_binary  <- matrix(nrow = n_row, ncol = dimEx )
      cl <- 0
      for (i in 1:n_nominal)     
      { 
        n_cat <- nlevels(data_nominal[, i])
        for (j in 1:n_cat)
        {     
          value <- levels(data_nominal[, i])[j]
          id_i  <- (data_nominal[,i] == value)
          data_binary[, cl+1] <- as.integer(id_i)
          BiNames[cl+1]   <- paste(names(data_nominal)[i], "_", value, sep="")
          cl <- cl+1
        }
        
      }
      numerMatrix  <- cbind(numerMatrix, data_binary)
      colnames(numerMatrix) <- c(numerMatrixnames, BiNames)
    }  
    
    if (!is.null(form))
    {
      numerMatrix <- data.frame(numerMatrix)
      numerMatrix <- cbind(numerMatrix, data_Y)
    }
    return(numerMatrix)
  }
