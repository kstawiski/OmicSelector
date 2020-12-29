#Copyright (C) 2018 Bing Zhu
# ===================================================================
#  Cost-sensitive C4.5 decision tree for binary classification
# ===================================================================

CSC45 <-  
    function(x, ...) 
        UseMethod("CSC45")

# default method
CSC45.data.frame <- 
    function(x, y, pruning = TRUE, MDL = TRUE, minIns = 2, costRatio  = 11/56)
    {           
        # Input:
        #     x         :  A data frame of the predictors from training data
        #     y         :  A vector of response variable from training data
        #     pruning   :  A logical number to determine whether to prune the tree. If pruning=TRUE, do the pruning process
        #     minIns    :  Minimum number of instances for split
        #     costRatio :  CostRatio between Majority class and Minority class

        # input check 
        Call <- match.call()
        indx <- match(c("x", "y"), names(Call), nomatch = 0L)
        if (indx[1] == 0L | indx[2] == 0L) stop("'predictors' and 'label' arguments are required")
        # initialization
        
        numRow      <- dim(x)[1]
        numColumn   <- dim(x)[2]
        classTable  <- table(y)
        classTable  <- sort(classTable, decreasing = TRUE)
        classLabels <- names(classTable)
        indexMajorIns <- y == classLabels[1] 
        attIndex    <- c(1 : numColumn)
        insIndex    <- c(1 : numRow)
        costAtt  <- rep(NA, numRow)
        costAtt[indexMajorIns]  <- numRow/((sum(indexMajorIns) + sum(!indexMajorIns)*costRatio))
        costAtt[!indexMajorIns] <- costRatio*numRow/(sum(indexMajorIns) + sum(!indexMajorIns)*costRatio) 
        
        
        # find the categorical and continous attribute
        isOrdered     <- sapply(x, is.ordered)
        isFactor      <- sapply(x, is.factor)
        isCategorical <- !isOrdered & isFactor  
        categoricalIndex <- which(isCategorical == TRUE)
        continousIndex   <- which(isCategorical == FALSE)
        
        # function for calculating entropy
        Entropy  <- function(x, cost){
            subset <- split(cost, x)
            prop   <- sapply(subset, sum)/sum(cost)
            H     <- -sum(ifelse(prop > 0, prop*log2(prop), 0))
            return(H)
        }
        
        
        # calculate entropy of categorical attribute
        CategoricalEntropy <- function(attIndex, insIndex)
        { 
            label <- y[insIndex]
            att   <- x[insIndex, attIndex] 
            att   <- factor(att) # delete empty level   
            cost  <- costAtt[insIndex]
            labelSubset   <- split(label, att)
            costSubset    <- split(cost, att)
            lengthSet     <- sapply(labelSubset, length)
            weightSet     <- sapply(costSubset, sum)
            minCheck      <- sum(lengthSet >= minIns)
            if (minCheck >= 2){
                labelEntropy  <- Entropy(label, cost) 
                attEntropy    <- Entropy(att, cost) 
                subsetEntropy <- mapply(Entropy, labelSubset, costSubset) 
                weight        <- weightSet/sum(cost)
                splitEntropy  <- crossprod(weight, subsetEntropy)
                gainRatio     <- (labelEntropy - splitEntropy)/attEntropy
            }
            else
                gainRatio <- -Inf 
            return (gainRatio) 
        } 
        
        # function for calculating entropy of continous attribute  
        ContinousEntropy <- function(attIndex, insIndex)
        {
            gainRatio  <- -Inf
            splitPoint <- -Inf
            att <- x[insIndex, attIndex]
            numValue   <- length(unique(att))
            if (numValue != 1){ 
                numIns     <- length(insIndex)
                label <- y[insIndex]
                sortIndex   <- order(att, label)
                attSorted   <- att[sortIndex]
                labelSorted <- label[sortIndex]
                splitCandidate  <- which(attSorted[-1] != attSorted[-numIns])
                indexNext   <- splitCandidate + 1    
                indexWithin <- c(1, splitCandidate + 1) 
                isEqualNext <- labelSorted[splitCandidate] == labelSorted[indexNext]
                isEqualWithin  <- labelSorted[c(splitCandidate,numIns)] == labelSorted[indexWithin]
                unuseCandidate <- isEqualNext &  isEqualWithin[-1]  & isEqualWithin[-c(numValue)] 
                splitCandidate <- splitCandidate[!unuseCandidate] 
                minCheck    <- splitCandidate  >= minIns &  splitCandidate <= numIns - minIns
                if (sum(minCheck) != 0) { 
                    cost  <- costAtt[insIndex]
                    costSorted  <- cost[sortIndex]
                    labelEntropy <- Entropy(labelSorted, costSorted) 
                    splitCandidate  <- splitCandidate[minCheck]
                    costSum <- sum(cost)
                    if (length(splitCandidate) == 1) {
                        entropyLeft <- Entropy(labelSorted[1:splitCandidate], costSorted[1:splitCandidate])
                        entropyRight<- Entropy(labelSorted[(splitCandidate+1):numIns], costSorted[(splitCandidate+1):numIns])  
                        weight      <- c(sum(costSorted[1:splitCandidate])/costSum, 1 - sum(costSorted[1:splitCandidate])/costSum)    
                        splitEntropy<- crossprod(weight, c(entropyLeft, entropyRight))
                        attEntropy  <- -sum(weight*log2(weight))
                        splitPoint  <- attSorted[splitCandidate]
                    }
                    
                    if (length(splitCandidate) > 1) {
                        
                        entropyLeft <- sapply(splitCandidate, function(x) Entropy(labelSorted[1:x], costSorted[1:x]))
                        entropyRight<- sapply(splitCandidate, function(x) Entropy(labelSorted[(x+1):numIns], costSorted[(x+1):numIns]))
                        weightLeft  <- sapply(splitCandidate, function(x) sum(costSorted[1:x]))   
                        weightRight <- sapply(splitCandidate, function(x) sum(costSorted[(x+1):numIns])) 
                        weightCandidate <- cbind(weightLeft, weightRight)/costSum
                        splitEntropyCandicate <- rowSums(weightCandidate*cbind(entropyLeft, entropyRight))
                        opt          <- which.min(splitEntropyCandicate)  
                        splitEntropy <- splitEntropyCandicate[opt]
                        weight       <- weightCandidate[opt,] 
                        attEntropy  <- -sum(weight*log2(weight))
                        splitPoint   <- attSorted[splitCandidate[opt]] 
                    }
                    if (MDL)
                        gainRatio   <- (labelEntropy - splitEntropy - log2(numValue-1)/numIns)/attEntropy  # MDL correction 
                    else
                        gainRatio   <- (labelEntropy - splitEntropy)/attEntropy
                }
            }
            optimalSplit <- c(gainRatio = gainRatio, splitPoint = splitPoint)
            return (optimalSplit)
        }
        
        # create new node
        CreateNode <- function (nodeCode = NULL, isLeaf = TRUE, insIndex = NULL, splitInfo = "NO" )
        {
            numIns <- length(insIndex)
            label  <- y[insIndex]
            cost   <- costAtt[insIndex]
            #classTable <- table(label)
            #classTable <- sort(classTable, decreasing = TRUE)
            indexMaj <- label == classLabels[1]
            indexMin <- label == classLabels[2]
            
            if (sum(cost[indexMaj]) > sum(cost[indexMin])){
                predClass <- classLabels[1]
                accuracy <- sum(cost[indexMaj])/sum(cost)
            } else {
                predClass <- classLabels[2]
                accuracy <- sum(cost[indexMin])/sum(cost)
            }
            error <- 1 - accuracy
            upperError <- numIns*(error + 0.69*0.69/(2*numIns)+0.69*sqrt(error/numIns+error*error/numIns+0.69*0.69/(4*numIns*numIns)))/(1+0.69*0.69/numIns)
            predProb   <- accuracy
            newNode <- list(
                nodeCode    = nodeCode    ,
                isLeaf      = isLeaf      ,
                insIndex    = insIndex    ,
                splitInfo   = splitInfo   ,
                predClass   = predClass   ,
                predProb    = predProb    ,
                upperError  = upperError
            ) 
            class(newNode) <- "node" 
            newNode
        }
        
        ReshapeTree <- function(tree){ 
            newList  <- vector(mode = "list", length = 10)
            newTree  <- vector(mode = "list", length = 10)
            ind <- 0
            for (i in sequence(length(tree))){   
                if (class(tree[[i]]) == "node")
                    numElement <- 1
                else 
                    numElement <- length(tree[[i]])  
                if ((ind + numElement) > length(newTree))
                    newTree <- append(newTree, newList)
                if (numElement == 1)
                    newTree[[ind+1]] <- tree[[i]]
                else
                    newTree[(ind+1):(ind + numElement)] <- unlist(tree[i], recursive = FALSE)
                ind <- ind + numElement
            }
            newTree <- newTree[1:ind]
            return(newTree)
        }
        
        
        
        # function for dividing and conquer partition
        Partition <- function(insIndex, nodeCode, attIndex) 
        {
            classDistribution <- table(y[insIndex])
            isLeaf<- any(classDistribution == length(insIndex)) | length(insIndex) < 2*minIns 
            
            if (!isLeaf) 
            {  
                
                # calculate information gain for categorical attribute
                maxCategoricalEntropy <- -Inf
                if (sum(isCategorical[attIndex]) != 0){  
                    currentCategorical <- intersect(categoricalIndex, attIndex)
                    entropyCategorical <- sapply(currentCategorical, CategoricalEntropy, insIndex = insIndex)
                    maxCategoricalEntropy <- max(entropyCategorical) 
                }   
                
                # calculate information gain for continous attribute
                maxContinousEntropy <- -Inf
                if (sum(!isCategorical[attIndex]) != 0){ 
                    currentContinous  <- intersect(continousIndex, attIndex)
                    entropyContinous  <- sapply(currentContinous, ContinousEntropy, insIndex = insIndex)
                    maxContinousEntropy <- max(entropyContinous[1, ])
                }
                isLeaf <- TRUE 
                if (maxCategoricalEntropy > maxContinousEntropy){
                    isLeaf      <- FALSE
                    splitAttIndex    <- which.max(entropyCategorical)
                    splitRule   <- x[insIndex, currentCategorical[splitAttIndex]]
                    splitRule   <- factor(splitRule)
                    splitInfo   <- list(splitAtt = currentCategorical[splitAttIndex], gainRatio = entropyCategorical[splitAttIndex], splitType = "Discrete",  splitPoint = levels(splitRule))
                    attNext     <- setdiff(attIndex, currentCategorical[splitAttIndex])
                    nodeCodeSet <- paste0(nodeCode, 1:nlevels(splitRule))
                    splitSubset <- split(insIndex, splitRule)
                }
                if (maxContinousEntropy > maxCategoricalEntropy){
                    isLeaf    <- FALSE
                    splitAttIndex  <- which.max(entropyContinous[1,])     
                    splitInfo <- list(splitAtt = currentContinous[splitAttIndex], gainRatio = entropyContinous[1, splitAttIndex],  splitType = "Continous",  splitPoint = entropyContinous[2, splitAttIndex] )
                    splitRule <- x[insIndex, currentContinous[splitAttIndex]] > splitInfo$splitPoint
                    attNext   <- attIndex
                    nodeCodeSet <- paste0(nodeCode, 1:2)
                    splitSubset <- split(insIndex, splitRule)
                }
            }
            
            if (!isLeaf){
                childNodes <- mapply(Partition, splitSubset, nodeCodeSet , MoreArgs = list(attIndex = attNext), SIMPLIFY= FALSE)
                names(childNodes) <- NULL
                isAllLeaf <- sum(sapply(childNodes, class) == "node") == length(childNodes)
                if (!isAllLeaf)
                    childNodes <- ReshapeTree(childNodes)
                interNode <- CreateNode(nodeCode = nodeCode, isLeaf = FALSE, insIndex = insIndex,  splitInfo = splitInfo)
                subsTree  <- append(childNodes, list(interNode))
            } else {
                leafNode <- CreateNode( nodeCode =  nodeCode, isLeaf = TRUE, insIndex = insIndex)   
                subsTree <- leafNode
            }
            return (subsTree)
        }
        
        rootCode <- "0"
        tree <- Partition(insIndex, rootCode, attIndex)
        
        # function for pruning
        prune <- function(node){
            if (node$isLeaf == TRUE){
                leafError <- node$upperError
                leafInfo  <- list(leafError = node$upperError, leafSet= node$nodeCode, pruneSet = NULL, newleaf = NULL) 
                return(leafInfo)
            } else {
                childIndex <- FindChild(node, nodeCodeList) 
                branches   <- lapply(tree[childIndex], prune)
                sumLeafError <- sum(sapply(branches,"[[", 1))
                parentError  <- node$upperError
                preLeafSet   <- unlist(lapply(branches,"[[", 2))
                leafSet <- c(preLeafSet, node$nodeCode)
                if (sumLeafError > parentError){
                    leafInfo <- list(leafError = sumLeafError, leafSet = leafSet, pruneSet = preLeafSet, newLeafset = node$nodeCode)
                } else {
                    pruneSet  <- unlist(lapply(branches,"[[", 3))
                    newLeafset   <- unlist(lapply(branches,"[[", 4))
                    leafInfo <- list(leafError = sumLeafError, leafSet = leafSet, pruneSet = pruneSet, newLeafset = newLeafset)
                }
                return(leafInfo)
            }
        }
        # prune tree 
        rootIndex <- length(tree)
        nodeCodeList <- sapply(tree, '[[', 1)
        
        if (pruning){
            pruneResult <- prune(tree[[rootIndex]])
            prunedNodeInex <- nodeCodeList %in% pruneResult$pruneSet
            newLeafIndex   <- nodeCodeList %in% pruneResult$newLeafset
            CreateNewLeaf <- function(node){
                node$isLeaf <- TRUE 
                node$splitInfo <- "NO"
                node
            }
            if (any(newLeafIndex))
                tree[newLeafIndex] <- lapply(tree[newLeafIndex], CreateNewLeaf)
            tree <- tree[!prunedNodeInex]
        }
        
        # define the C4.5 class 
        object <- list(
            pruning      = pruning    ,
            tree         = tree       ,
            classLabels  = classLabels,
            costRatio    = costRatio  ,
            #formula      = form       , 
            MDL          = MDL
        )
        class(object) <- "CSC45"
        object
    }

predict.CSC45 <- function(object, x.test, type = "class")
{
    if (!type %in% c("class", "prob"))
        stop("type must be class or prob")
    data       <- x.test
    tree       <- object$tree
    numRow     <- dim(data)[1]
    insIndex   <- c(1 : numRow)
    rootIndex  <- length(tree)
    nodeCodeList <- sapply(tree, '[[', 1)
        
    
    GetPredict <- function(node, insIndex)
    {
        if (node$isLeaf == TRUE){
            numIns  <- length(insIndex)
            predictSubset <- data.frame(classPredict = rep(node$predClass, numIns), probPredict = rep(node$predProb, numIns), index = insIndex)
        } else {
            numUnknown <- 0
            att <- data[insIndex, node$splitInfo$splitAtt]
            if (node$splitInfo$splitType == "Discrete")
            {
                attValue    <- unique(att)
                indexKnown  <- att %in% node$splitInfo$splitPoint
                numUnknown  <- sum(!indexKnown)
                insIndexKnown <- insIndex[indexKnown]
                if (length(insIndexKnown) > 0){
                    splitRule   <- factor(att[indexKnown])
                    childSelect <- which (node$splitInfo$splitPoint %in% attValue)
                    childSelectCode   <- paste0(node$nodeCode,  childSelect)     
                }
                
            } else {
                insIndexKnown <- insIndex
                splitRule   <- att > node$splitInfo$splitPoint
                if (sum(splitRule) == 0)
                    childSelectCode <- paste0(node$nodeCode, "1")
                
                if (sum(!splitRule) == 0)
                    childSelectCode <- paste0(node$nodeCode, "2")
                
                if ((sum(splitRule) > 0) & (sum(!splitRule) > 0))
                    childSelectCode <- paste0(node$nodeCode, 1:2)
            }
            predictSubset <- data.frame()
            if (length(insIndexKnown) > 0){
                childSelectIndex  <- which(nodeCodeList %in%  childSelectCode)
                if (length(childSelectCode) == 1){
                    predictSubset <- GetPredict(tree[[childSelectIndex]], insIndexKnown)
                } else {
                    splitSubset <- split(insIndexKnown, splitRule)
                    predictList <- mapply(GetPredict, tree[childSelectIndex], splitSubset, SIMPLIFY = FALSE)
                    predictSubset <- do.call(rbind, predictList)
                }  
            }
            if  (numUnknown > 0){
                predictUnkown <- data.frame(classPredict = rep(node$predClass, numUnknown), probPredict = rep(node$predProb, numUnknown), index = insIndex[!indexKnown])
                predictSubset <- rbind(predictSubset, predictUnkown)
                
            }
        }
        return(predictSubset)
    }
    
    prediction <- GetPredict(tree[[rootIndex]], insIndex)
    sortIndex  <- order(prediction[,3])
    prediction <- prediction[sortIndex,]
    if (type == "prob"){
        output <- matrix(rep(NA, 2*length(insIndex)), ncol = 2)
        index <- prediction[, 1] == object$classLabel[2]
        output[index,  2]  <- prediction[index, 2]
        output[!index, 2]  <- 1- prediction[!index, 2]
        output[, 1]        <- 1- output[, 2]
        colnames(output) <- object$classLabel 
    } else {
        output <- prediction[, 1]
    }
    return(output)
}


# find child
FindChild <- function(node, nodeCodeList)
{
    if (node$splitInfo$splitType == "Discrete")
        childCodeSet <- paste0(node$nodeCode, 1:length(node$splitInfo$splitPoit))
    else 
        childCodeSet <- paste0(node$nodeCode, 1:2)
    childIndex <- nodeCodeList %in%  childCodeSet
    return(childIndex)
}

