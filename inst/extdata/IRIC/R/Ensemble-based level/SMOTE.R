# =====================================================
#  SMOTE sampling
# =====================================================

SMOTE <- 
  function(form, data, perOver = 500, k = 5)
    # INPUTS:
    #    form: model formula
    #    data: original  dataset 
    #    perOver/100: number of new instance generated for each minority instance 
    #    k: number of nearest  neighbours 
  {
    
    # find the class variable
    tgt <- which(names(data) == as.character(form[[2]]))
    classTable<- table(data[, tgt])
    numCol <- dim(data)[2]
    
    # find the minority and majority instances
    minClass  <- names(which.min(classTable))
    indexMin  <- which(data[, tgt] == minClass)
    numMin    <- length(indexMin)
    majClass  <- names(which.max(classTable))
    indexMaj  <- which(data[, tgt] == majClass)
    numMaj    <- length(indexMaj)
    
    # move the class variable to the last column

    if (tgt < numCol)
    {
      cols <- 1:numCol
      cols[c(tgt, numCol)] <- cols[c(numCol, tgt)]
      data <- data[, cols]
    }
    # generate synthetic minority instances
    source("SmoteExs.R")
    if (perOver < 100)
    {
      indexMinSelect <- sample(1:numMin, round(numMin*perOver/100))
      dataMinSelect  <- data[indexMin[indexMinSelect], ]
      perOver <- 100
    } else {
      dataMinSelect <- data[indexMin, ]
    }
    
    newExs <- SmoteExs(dataMinSelect, perOver, k)
    
    # move the class variable back to original position
    if (tgt < numCol) 
    {
      newExs <- newExs[, cols]
      data   <- data[, cols]
    }
    
    # unsample for the majority intances
    newData <- rbind(data, newExs)

    return(newData)
  }

