###############################################################
# Credit: Gaurav Sharma                            #
# Email : ergauravrandev@gmail.com                                    #
###############################################################

cat("\nSTART\n")
startTime = proc.time()[3]
startTime



##Include Library

#install.packages("FSelector")
#install.packages("caret", dependencies = TRUE)
library(randomForest)
library(FSelector)
library(caret)
library(pROC)





##Variable Declaration and Seed Setting

modelName <- "RandomForest"
modelName
InputDataFileName="Kmeans4.csv"
InputDataFileName

seed=sample(99999:999999, 1)
seed
set.seed(seed)



##Data Division for train/test

training=70      # define percentage



##Load the data from file

dataset <- read.csv(InputDataFileName)    # Read the datafile
head(dataset)   # Show Top 6 records
nrow(dataset)   # Show number of records
names(dataset)  # Show fields names or columns names
ncol(dataset)


##Count total number of observations/rows.

totalDataset <- nrow(dataset)
totalDataset



##Choose Target variable

target  <- names(dataset)[1]   # i.e. cluster
target



##Choose inputs Variables

inputs <- setdiff(names(dataset),target)
inputs
length(inputs)

selectedInputs <- inputs


##Select Training Data Set ( Select random indices )

trainSample <- sample(totalDataset, totalDataset * training/100)
head(trainSample)     # Show Top 6 indecies
length(trainSample)   # Show number of indecies
ncol(trainSample)

testSample <- setdiff(seq_len(nrow(dataset)), trainSample)
head(testSample)      # Show Top 6 indecies
length(testSample) 
ncol(testSample)

####Feature Weighting (4 techniques, uncomment to run)

#weights = information.gain(cluster~., dataset, unit = "log2")

#weights <- gain.ratio(cluster~., dataset)

#weights = chi.squared(cluster~., dataset)

weights = random.forest.importance(cluster~., dataset, importance.type = 1)


subset = cutoff.k(weights, 5)

evaluator = function(subset) {
  k = 5  
  splits <- runif(nrow(dataset))
  results = sapply(1:k, function(i) {
    test.idx <- (splits >= (i - 1) / k) & (splits < i / k)
    train.idx <- !test.idx
    test <- dataset[test.idx, , drop=FALSE]
    train <- dataset[train.idx, , drop=FALSE]
    tree  = randomForest(as.simple.formula(subset, "cluster"), train, ntree=500,mtry=5)
    error.rate = sum(test$cluster != round(predict(tree, test),0)) / nrow(test)
    return(1 - error.rate)
  })
  return(mean(results))
}


####Optimal Feature subset finding (4 techniques, uncomment to run)

#attr.subset = backward.search(names(dataset)[-1], evaluator)

attr.subset = forward.search(names(dataset)[-1], evaluator)

#attr.subset = best.first.search(names(dataset)[-1], evaluator, max.backtracks = 10)

#attr.subset = hill.climbing.search(names(dataset)[-1], evaluator)


##Feature Selection based on Optimal Feature subset finding technique

selectedInputs <- attr.subset

##Select Training Data Set ( Select random indices )

trainDataset <- dataset[trainSample,c(selectedInputs, target)]
head(trainDataset)    # Show Top 6 records
nrow(trainDataset)

##Select Testing Data Set

testDataset <- dataset[testSample,c(selectedInputs, target)]
head(testDataset)
nrow(testDataset)

##Model Building: Training

formula <- as.formula(paste(target,"~",paste(c(selectedInputs),collapse = "+")))
formula
print(formula)

model   <- randomForest(formula, trainDataset, ntree=500,mtry=5)
print(model)

##Prediction (Testing)

Predicted <- round(predict(model, testDataset, type = 'class'))
Predicted

##Extracting Actual

Actual <- as.double(unlist(testDataset[target]))
head(Actual)

##Model Evaluation

accuracy <- round(mean(Actual==Predicted)*100,2)
accuracy

cm <- confusionMatrix(Actual, Predicted)    #Generates confusion Matrix

Pred <- (predict(model, testDataset, type = 'response' ))
roc.multi <- multiclass.roc(testDataset$cluster, Pred)
roc.multi                                                 #Generates AUC 


## K-fold cross validation (10-fold)

eval = function(attr.subset) {
  k = 10 
  splits <- runif(nrow(trainDataset))
  results = sapply(1:k, function(i) {
    test.id <- (splits >= (i - 1) / k) & (splits < i / k)
    train.id <- !test.id
    test <- dataset[test.id, , drop=FALSE]
    train <- dataset[train.id, , drop=FALSE]
    tree  = randomForest(formula, trainDataset, ntree=500,mtry=5)
    error.rate = sum(test$cluster != round(predict(tree, test),0)) / nrow(test)
    #error.rate = sum(test$cluster != predict(tree, test, type="class")) / nrow(test)
    return(1 - error.rate)
  })
  return (results)
}

Kfold <- eval()

totalTime = proc.time()[3] - startTime
print (totalTime)

##Save evaluation resut

result <- data.frame(modelName,accuracy, totalTime)[1:1,]
print(result)

##Writing to file

write.csv(result, file=paste(modelName,"-Evaluation-Result.csv",sep=''), row.names=FALSE)

write.csv(data.frame(Actual,Predicted), file=paste(modelName,"-ActualPredicted-Result.csv",sep=''), row.names=FALSE)

write.csv(cm$table, file=paste(modelName,"-cm.csv",sep=''), row.names=TRUE)

write.csv(cm$overall, file=paste(modelName,"-cm-overall.csv",sep=''), row.names=TRUE)

write.csv(cm$byClass, file=paste(modelName,"-cm-class.csv",sep=''), row.names=TRUE)

write.csv(Kfold, file=paste(modelName,"-CV.csv",sep=''), row.names=FALSE)

save.image(file=paste(modelName,"-Model.RData",sep=''))

