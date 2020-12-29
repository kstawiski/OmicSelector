# IRIC： Integrated R Library for Imbalanced Classification
IRIC is an R library for imbalanced classification,  which will bring convenience to users by integrating a wide set of solutions into one library.



- [Introduction](#Introduction)
- [Installation](#Installation)
- [Examples](#Examples)

---

## Introduction	
The current version of IRIC (v1.1) provides a set of 19 approaches for imbalanced classification, in which 8 approaches are new implementations in R. All these approaches can be classified into 3 strategies: data level, algorithm level and ensemble-based strategy. In addition, we provide parallel implementations of Bagging-based solution to improve the efficiency of model building. All approaches in IRIC are presented in the table below.


<table border=0 cellpadding=0 cellspacing=0 width=843 style='border-collapse:
 collapse;table-layout:fixed;width:632pt'>
 <col width=189 style='mso-width-source:userset;mso-width-alt:6048;width:142pt'>
 <col width=171 style='mso-width-source:userset;mso-width-alt:5472;width:128pt'>
 <col width=411 style='mso-width-source:userset;mso-width-alt:13152;width:308pt'>
 <col width=72 style='width:54pt'>
 <tr height=22 style='height:16.5pt'>
  <th height=22 class=xl70 width=189 style='height:16.5pt;width:142pt'>Strategy</th>
  <th class=xl70 width=171 style='width:128pt'>Submodule</th>
  <th class=xl70 width=411 style='width:308pt'>Method</th>
 </tr>
 <tr height=22 style='height:16.5pt'>
  <td height=22 class=xl70 style='height:16.5pt;border-top:none'>Algorithm
  level</td>
  <td class=xl70 style='border-top:none'>Cost-sensitive learning</td>
  <td class=xl70 style='border-top:none'>CSC4.5</td>
 </tr>
 <tr height=22 style='height:16.5pt'>
  <td rowspan=3 height=66 class=xl71 style='border-bottom:.5pt solid black;
  height:49.5pt;border-top:none'>Data level</td>
  <td class=xl66 style='border-top:none'>Oversampling</td>
  <td class=xl66 style='border-top:none'>ADASYN, MWMOTE, Random
  Oversampling, SMOTE</td>
 </tr>
 <tr height=22 style='height:16.5pt'>
  <td height=22 class=xl65 style='height:16.5pt'>Undersampling</td>
  <td class=xl65>CLUS，Random Undersampling</td>
 </tr>
 <tr height=22 style='height:16.5pt'>
  <td height=22 class=xl69 style='height:16.5pt'>Hybrid Sampling</td>
  <td class=xl69>SmoteENN, SmoteTL, SPIDER</td>
 </tr>
 <tr height=22 style='height:16.5pt'>
  <td rowspan=3 height=66 class=xl67 style='border-bottom:.5pt solid black;
  height:49.5pt'>Ensemble-based learning</td>
  <td class=xl65>BalanceBagging</td>
  <td class=xl65>RBBagging, ROSBagging,RUSBagging,SMOTEBagging</td>
 </tr>
 <tr height=22 style='height:16.5pt'>
  <td height=22 class=xl65 style='height:16.5pt'>BalanceBoost</td>
  <td class=xl65> AdaC2, RUSBoost, SMOTEBoost</td>
 </tr>
 <tr height=22 style='height:16.5pt'>
  <td height=22 class=xl69 style='height:16.5pt'>Hybrid Ensemble</td>
  <td class=xl69>BalanceCascade, EasyEnsemble </td>
 </tr>
</table>




## Installation
Download the code from GitHub repository before and then apply the techniques.  R version >= 3.1.
## Examples
SMOTE(Data level), CSC4.5 (Algorithm level) and RBBagging (Ensemble-based level) are presented as examples of IRIC's usage.
- [SMOTE](#SMOTE)
- [CSC4.5](#CSC4.5)
- [RBBagging](#RBBagging)
#### SMOTE
```
#Example of SMOTE
#Load the package caret for data partitioning
library(caret)
#Load data set
load("Korean.RDa")
#Run the script file of SMOTE
source("SMOTE.R")
#data split
sub <- createDataPartition(Korean$Churn,p=0.75,list=FALSE)
trainset <- Korean[sub,]
testset <- Korean[-sub,]
x <- trainset[, -11]
y <- trainset[, 11]
#call the SMOTE
newData<- SMOTE(x, y)
```
#### CSC4.5
```
#Example of CSC4.5 
#load CSC4.5
source("CSC45.R")
library(caret)
#Load data set
load("Korean.RDa")
#Data split
sub <- createDataPartition(Korean$Churn,p=0.75,list=FALSE)
trainset <- Korean[sub,]
testset <- Korean[-sub,]
x <- trainset[, -11]
y <- trainset[, 11]
#training model
model <- CSC45(x, y, pruning = TRUE)
#Prediction
output <- predict (model, x) 
```
#### RBBagging
```
#Example of RBBagging 
#Load the package caret for data partitioning
library (caret) 
#Load data set 
load(”Korean.RDa”) 
#Run the script file of RBBagging 
source(”BalanceBagging.R”)
#Data split
sub <- createDataPartition(Korean$Churn, p=0.75,list=FALSE)
trainset <- Korean[sub,]
testset <- Korean[-sub,]
x <- trainset[, -11]
y <- trainset[, 11]
#call the RBBaging for model training train 
model <- bbagging(x, y, type=”RBBagging", allowParallel=TRUE)
#prediction
output <- predict (model, x)
```


