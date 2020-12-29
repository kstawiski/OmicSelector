FindLabel <-
  function(label){
        out <- list()
        classTable  <- table(label)
        classTable  <- sort(classTable, decreasing = TRUE)
        classLabels <- names(classTable)
        negLabel  <- classLabels[1]
        posLabel  <- classLabels[2]
        out$neg <-  negLabel
        out$pos <-  posLabel
        out
      }
