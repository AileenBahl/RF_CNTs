require(tools)
require(xlsx)
require(class)
require(randomForest)
require(caret)
library(imbalance)

set.seed(2509)

myPath <- "I:/CaseStudy_OmicsMWCNT/"
myFile <- "NM_GSEA_metaInfo_final.xlsx"
group1 <- "LongTangled"
group2 <- "LongRigid"


### read and select data

myPathwaysTable <- read.xlsx(paste0(myPath, myFile),1) 

myPathwaysTable_groupsSelected <- myPathwaysTable[which(myPathwaysTable$group == group1| myPathwaysTable$group == group2),]
myDescriptorsAndTox <- myPathwaysTable_groupsSelected[,which(grepl("HALLMARK",colnames(myPathwaysTable_groupsSelected)) | grepl("LungInflammationKeyEvent",colnames(myPathwaysTable_groupsSelected)) | grepl("group",colnames(myPathwaysTable_groupsSelected)))]


### oversampling

numberAdditionalSamples <- max(table(myDescriptorsAndTox$group))-min(table(myDescriptorsAndTox$group))
myDescriptorsAndTox_oversampled <- mwmote(myDescriptorsAndTox, numInstances = numberAdditionalSamples, classAttr = "group")
myDescriptorsAndTox_new <- rbind(myDescriptorsAndTox,myDescriptorsAndTox_oversampled)

myDescriptors <- myDescriptorsAndTox_new[,-which(colnames(myDescriptorsAndTox_new)== "group")]
myDescriptors_numeric <- as.data.frame(apply(myDescriptors,2,as.numeric))
myTox <- myDescriptorsAndTox_new[,which(colnames(myDescriptorsAndTox_new)== "group")]


### perform rf

predictRF <- function(i){
  myRF <- randomForest(myDescriptors_numeric[-i,],factor(myTox[-i], levels=c(group1, group2)), myDescriptors_numeric[i,],factor(myTox[i], levels=c(group1, group2)), ntree=5000, importance=T)
  myImportance <- round(importance(myRF)[,3], 2)
  myPredictedClass <- myRF$test$predicted
  myList=list(prediction=myPredictedClass,importance=myImportance)
  return(myList)
}

set.seed(2509)
flds <- createFolds(myTox, k = 10, list = TRUE, returnTrain = FALSE)

myPredictions <- sapply(flds,predictRF)


myComparisonOfClasses <- data.frame(PredictedClass=unlist(myPredictions[1,]),TrueClass=myTox[unlist(flds)])
rownames(myComparisonOfClasses) <- gsub(".*\\.","",names(unlist(myPredictions[1,])))
print("RF - full model:")

myImportanceMatrix <- do.call(cbind, myPredictions[2,])
myParameterImportance <- sort(rowMeans(myImportanceMatrix),decreasing=T)

print(myParameterImportance)
print(myComparisonOfClasses)

myCorrectPrediction <- sum(myComparisonOfClasses$PredictedClass==myComparisonOfClasses$TrueClass)
print("Number of correct predictions:")
print(myCorrectPrediction)

mySensitivity <- sum(myComparisonOfClasses$PredictedClass == group1 & myComparisonOfClasses$TrueClass == group1) / sum(myComparisonOfClasses$PredictedClass == group1)
print("Sensitivity:")
print(mySensitivity)

mySpecificity <- sum(myComparisonOfClasses$PredictedClass == group2 & myComparisonOfClasses$TrueClass == group2) / sum(myComparisonOfClasses$PredictedClass == group2)
print("Specificity:")
print(mySpecificity)

myBalancedAccuracy <- (mySensitivity + mySpecificity)/2 
print("Balanced accuracy:")
print(myBalancedAccuracy)
cat("\n")


### Prediction for short fibres

myShort <- myPathwaysTable[which(myPathwaysTable$group == "ShortFibre"),]

myDescriptors_short_numeric <- as.data.frame(apply(myShort[,c(2:44)],2,as.numeric))

myRF <- randomForest(myDescriptors_numeric,factor(myTox, levels=c(group1,group2)), myDescriptors_short_numeric, ntree=5000, importance=T)
myImportance <- round(importance(myRF)[,3], 2)

myPredictedClass <- myRF$test$predicted

print(data.frame(Predicted=myPredictedClass, Sample=myShort$condition))


