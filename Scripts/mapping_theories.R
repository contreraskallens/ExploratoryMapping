#### libraries ####
setwd(dir = 'GitHub/exploratory-mapping-draft/ExploratoryMapping/Scripts/')

require("MASS") #used to calculate the projection of new data in old SVD space.

# plotting #

require("reshape2")
require("ggplot2")
require("ggthemes")
require("scales")

# dendrogram tree cutting from https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/BranchCutting/ #
require("dynamicTreeCut")
require("bio3d")
require("moduleColor")

# paralellization of prediction machine and UI #
require("log4r")
require("foreach")
require("doParallel")

####FUNCTIONS####

cosineGen <- function(matrix){
  lengthVec <- sqrt(rowSums(matrix * matrix))
  tcrossprod(matrix) / (lengthVec %o% lengthVec)
}  #Function that generates the cosine between each row of a matrix.

calcEnt <- function(matrix){
  workMatrix <- t(matrix) #transposes for division
  a <- workMatrix / rowSums(workMatrix) #generates a probability matrix
  b <- 1 + ((rowSums(a * log2(a), na.rm = T)) / log2(dim(matrix)[1])) #calculate entropy (1 + (sum of probability times log2 probability, divided by the total number of documents)).
  workMatrix <- log(workMatrix[which(b > 0, arr.ind = T), ] + 1) #log normalizes frequency matrix and deletes 0 entropy ones.
  workMatrix <- workMatrix * b[b > 0] #weight log normalized matrix by multiplying terms with entropy higher than 0 by its entropy.
  return(t(workMatrix)) #returns original, non-transposed matrix.
} #calculates the entropy of each term of a matrix. Uses formula in Martin & Berry, 2007, "Mathematical foundations behind latent semantic analysis".

compareTheories <- function(matrix, cat){
  ##doc: takes matrix as the result of an svd(u). populates a pre-allocated list of every topic in external topicList (not generalized yet) with matrices of an "index" of similarity using each dimension (column) in matrix for each topic.
  resultList <- lapply(lapply(1:8, matrix, data = 0, nrow = 8, ncol = dim(matrix)[2]), function(x){row.names(x) <- topicList; return(x)}) #pre-allocation of list
  cosineAvgsList <- lapply(1:8, matrix, data = 0, nrow = 8, ncol = dim(matrix)[1]) #pre-allocation of second list
  names(resultList) <- topicList #name for easy accessing theories
  indexMatrix <- matrix(FALSE, nrow = dim(matrix)[1], ncol = 8, dimnames = list(1:dim(matrix)[1], topicList)) #pre-allocates logical matrix
  for(topic in topicList){indexMatrix[,topic] <- cat[,2] == topic} #populates logical matrix with a logical mask reflecting catalog$id
  docsByTop <- colSums(indexMatrix) #number of documents for each topic
  n <- 1 #counter for each dimension
  while(n <= dim(matrix)[2]){ #loops through dimensions
    database <- matrix[, 1:n] #slices dimensions
    if(n == 1){database <- cbind(database, 0)} #if it has only one dimension, then add a column of 0s to make cosineGen work
    database <- cosineGen(database) #produces a x by x matrix of cosines between each paper.
    database[is.na(database)] <- 0 #replaces NA with 0.
    meanMatrix <- crossprod(indexMatrix, database) #produces a matrix with the sum of cosines of each paper with each of the topics
    meanMatrix <- meanMatrix / docsByTop #produces a matrix with the mean cosine of each paper with each of the topics
    cosineAvgsList[[n]] <- meanMatrix #stores the matrix of means in a list with n as index for dimensions used.
    meanMatrix <- meanMatrix %*% indexMatrix #produces a vector with the sum of mean cosines for each topic against each topic in dimension n
    meanMatrix <- t(meanMatrix) / docsByTop #produces a vector of the means of sums of mean cosines for each topic against each topic in dimension n.
    for(topic in topicList){ #loops through topics to populate results of all cosines.
      resultList[[topic]][, n] <- meanMatrix[topic,]
    }
    n = n + 1
  }
  returnList = list(resultList,cosineAvgsList) #makes list of lists with results and mean cosines
  return(returnList) #returns everything
} #function for calculating cosines of the whole matrix. "matrix" is the result of an SVD; "cat" is the catalog to obtain topic information (in this case, catalog$id). Returns a list of two lists: [[1]] is all cosines by paper, [[2]] is a list of matrices of mean distance of each paper with each of the other topics. [[1]][n] and [[2]][n] are the different dimensions resulting from SVD. 

plotTopicDiff <- function(topic, resultsList){
  workTable <- as.data.frame(melt(resultsList[[1]][[topic]], varnames = c("topic", "dimension"), value.name = "cosine")) #long form for ggplot
  plot <- ggplot(data = workTable, aes(x = dimension, y = cosine, color = topic, group = topic)) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) + labs(title = paste("Mean cosine of", topic, "papers with other theories and itself across dimensions")) + geom_line() + scale_colour_solarized("red") + geom_point(size = 0.7, shape = 3) + guides(colour = guide_legend(override.aes = list(size=3)))
  print(plot)
  } #function for plotting the mean distance of every topic with all other topics. "topic" is one of the topics of topicList; "resultsList" is the object that compareTheories() returns.

####DATA AND BASIC CLEANUP####

##ORIGINAL##

#loads files and catalog for original#

freqMatrix <- as.matrix(read.table('Data/Original/document_by_term.txt', sep='\t', header = T))[, -1] #loads the DBT minus one column, the identifier in the text file. 
row.names(freqMatrix) <- c(1:nrow(freqMatrix)) #row names with docID
freqMatrix <- freqMatrix[, apply(freqMatrix, 2, function(x){sum(x==0) < (dim(freqMatrix)[1] - 5)})] #removes columns with words that appear in fewer than 5 documents.

freqMatrix <- freqMatrix[which(rowSums(freqMatrix) > 0), ] #eliminates documents with 0 terms after cleanup of terminology
catalog <- read.table('catalog.txt', stringsAsFactors = F, sep = '\t', fill = T, quote = "") #loads catalog
catalog <- catalog[row.names(freqMatrix), ] #catalog also has row names as docID
colnames(catalog) = c('id','topic','year','authors','title','journal','abstract') #variable names for catalog
topicList <- unique(catalog$topic) #list of theories for analysis.

##REPLICATION##

#the same, but for replication documents#

repFreqMatrix <- as.matrix(read.table('rep_document_by_term.txt', sep='\t', header = T, quote = ""))[, -1] 
row.names(repFreqMatrix) <- c(1:nrow(repFreqMatrix))
repFreqMatrix <- repFreqMatrix[, apply(repFreqMatrix, 2, function(x){sum(x==0) < (dim(repFreqMatrix)[1] - 5)})] #removes columns with words that appear in fewer than 5 documents.
repFreqMatrix <- repFreqMatrix[which(rowSums(repFreqMatrix) > 0), ]
repCatalog <- read.table('rep_catalog.txt', stringsAsFactors = F, sep = '\t', fill = T, quote = "")
repCatalog <- repCatalog[row.names(repFreqMatrix), ]
colnames(repCatalog) = c('id','topic','year','authors','title','journal','abstract')
topicList <- unique(repCatalog$topic)

#NULL HYPOTHESIS#

#If you want to test the null hypothesis, change the parameter to T. It randomizes the theory of the papers in the databases.#

nullHyp <- F
if(nullHyp == T ){
  catalog$topic <- sample(catalog$topic, length(catalog$topic))
  repCatalog$topic <- sample(repCatalog$topic, length(repCatalog$topic))
}

####DATA PROCESSING (Latent Semantic Analysis)####
      
##ENTROPY##

#this uses the entropy function in calcEnt to weight the matrices with a log-entropy function#

#ORIGINAL#

cleanData <- calcEnt(freqMatrix) #entropy
cleanData[is.na(cleanData)] <- 0 #replace NA with 0.

#REPLICATION#
repCleanData <- calcEnt(repFreqMatrix) 
repCleanData[is.na(repCleanData)] <- 0

##SVD##

#ORIGINAL#

#dimensionality reduction#

wholeMacaroni = svd(cleanData) #partial SVD of 150 dimensions.

#REPLICATION#

repofWholeMacaroni = svd(repCleanData)

#ALLOCATE LOADING MATRICES#

documentLoadings <- wholeMacaroni$u %*% diag(wholeMacaroni$d)
termLoadings <- wholeMacaroni$v %*% diag(wholeMacaroni$d)
row.names(documentLoadings) <- row.names(cleanData)
row.names(termLoadings) <- colnames(cleanData)
termLoadings <- termLoadings[,1:15] 
documentLoadings <- documentLoadings[, 1:15]

#REPLICATION#

repDocumentLoadings <- repofWholeMacaroni$u %*% diag(repofWholeMacaroni$d)
repTermLoadings <- repofWholeMacaroni$v %*% diag(repofWholeMacaroni$d)
row.names(repDocumentLoadings) <- row.names(repCleanData)
row.names(repTermLoadings) <- colnames(repCleanData)
repTermLoadings <- repTermLoadings[,1:80] 
repDocumentLoadings <- repDocumentLoadings[, 1:80]

####GLM Models####

# This part of the script has the GLM models of each theory that attempt to predict the theory belonging of each paper #

## PARAMETERS ##

# Set parameters for the prediction. Min and max number of dimensions are used to control the number of dimensions that are to be used in the construction of the models. The procedure loops through the dimensions resulting from the SVD. Starts at minNumberOfDimensions (default: 3), stops at maxNumberOfDimensions (default: 50). "Method" refers to the data used to build and train the models. With "free", dimensions are selected for how well they predict theory belonging against every theory. With "cluster", the training is stratified to the "most similar" theories; e.g. "computational" is built using the dimensions that best predict computational papers when compared to 'bayesian' and 'connectionist'. "Repeats" is the number of iterations of the predicting process. "Source" controls which data set is to be used: "original" uses the original dataset, "replication" uses the replication data, and "cross" uses the projection of the replication data into the SVD space of the original dataset to predict their theories with the models built with the original dataset. #

minNumberOfDimensions <- 2 #lower boundary of D. does not work if lower than 2.
maxNumberOfDimensions <- 15 # upper boundary of D
repeats <- 100 # how many repetitions of prediction should be averaged?
method <- "free" # "cluster" or "free".
source <- "original" #"original" or "replication", "cross".

##OBJECTS##

# Loads the different objects depending on the parameter "source" on line 155.

if(source == "original"){
  my_catalog <- catalog
  my_svd <- documentLoadings
}
if(source == "replication"){
  my_catalog <- repCatalog
  my_svd <- repDocumentLoadings
}
if(source == "cross"){
  my_catalog <- catalog
  cross_catalog <- repCatalog
  my_svd <- documentLoadings

  replicationProjection <- matrix(0, nrow = nrow(repFreqMatrix), ncol = ncol(freqMatrix))   #create matrix for projection and populate#
  row.names(replicationProjection) <- row.names(repFreqMatrix) #set row names as wholeMacaroni
  colnames(replicationProjection) <- colnames(freqMatrix)
  sharedWords <- colnames(repFreqMatrix)[which(colnames(repFreqMatrix) %in% colnames(freqMatrix))] # obtain all the shared words between original data and replication data.
  replicationProjection[,sharedWords] <- repFreqMatrix[,sharedWords] #build a new DbT that has the documents in replication as rows and the shared words between both datasets as columns.
  
  cross_svd <- replicationProjection %*% ginv(termLoadings) #projects the replication data into the SV  space of the original dataset..
}

### TOPIC BEST PREDICTORS ###

## At first, the dimensions that differentiate the papers of a theory when compared to the other ones are selected. The first 80 of them are stored. #

##FREE FOR ALL##

#for each topic, fill the "best predictors" matrix with the dimensions that most differentiate that theory from the other 8 theories.

bestPredictors <- c()
predictorRatings <- c() #store ratings to compare positive versus negative in term loading inspection
for (topic in topicList) { 
  glmOutput = glm(my_catalog$topic==topic~.,data=data.frame(my_svd[,1:maxNumberOfDimensions]),family=binomial)
  bestPredictors = rbind(bestPredictors,data.frame(t(sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),ind=T,decreasing=T)$ix[1:maxNumberOfDimensions]))) #TODO: fix this. 2 is here to skip intercept.
  predictorRatings <- rbind(predictorRatings, glmOutput$coefficients[2:length(glmOutput$coefficients)][sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]), ind = T, decreasing = T)$ix])
}
row.names(bestPredictors) <- topicList
row.names(predictorRatings) <- topicList

##STRATIFIED SAMPLING##

#for each topic, fill the "best predictors" matrix with the dimensions that most differentiate that theory from the other theories in the cluster. Cluster 1 is "Classic" with computational, bayesian and connectionist. Cluster 2 is "alt" with ecological, embodied, dynamical, distributed, enactive.

topicListClassic <- c(topicList[1], topicList[2], topicList[8]) # topic list of classical cluster.
topicListAlt <- setdiff(topicList, topicListClassic) #topic list of alt cluster.
catalogClassic <- my_catalog[which(my_catalog$topic %in% topicListClassic),] #catalog of classic papers
catalogAlt <- my_catalog[which(my_catalog$topic %in% topicListAlt),] #catalog of alt papers

#then, same procedure as above but using these clusters#

bestPredictorsClusters <- c()
predictorRatingsClusters <- c()

for (topic in topicListClassic) {
  glmOutput = glm(catalogClassic$topic==topic~.,data=data.frame(my_svd[which(my_catalog$topic %in% topicListClassic, arr.ind = T),1:maxNumberOfDimensions]),family=binomial)
  bestPredictorsClusters = rbind(bestPredictorsClusters,data.frame(t(sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),ind=T,decreasing=T)$ix[1:maxNumberOfDimensions]))) #TODO: fix this. 2 is here to skip intercept.
  predictorRatingsClusters <- rbind(predictorRatingsClusters, glmOutput$coefficients[2:length(glmOutput$coefficients)][sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]), ind = T, decreasing = T)$ix])
  }

for (topic in topicListAlt) {
  glmOutput = glm(catalogAlt$topic==topic~.,data=data.frame(my_svd[which(my_catalog$topic %in% topicListAlt, arr.ind = T),1:maxNumberOfDimensions]),family=binomial)
  bestPredictorsClusters = rbind(bestPredictorsClusters,data.frame(t(sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),ind=T,decreasing=T)$ix[1:maxNumberOfDimensions]))) #TODO: fix this. 2 is here to skip intercept.
  predictorRatingsClusters <- rbind(predictorRatingsClusters, glmOutput$coefficients[2:length(glmOutput$coefficients)][sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]), ind = T, decreasing = T)$ix])
}
row.names(bestPredictorsClusters) <- c(topicListClassic, topicListAlt)
row.names(predictorRatingsClusters) <- c(topicListClassic, topicListAlt)

### EXECUTION OF MODEL ###

# Each theory has a model built and trained using the dimensions specified in minNumberOfDimensions and maxNumberOfDimensions. Each model is of the theories is a GLM. In each iteration, a random training set of 600 papers is selected i (random theory belonging if "free", 60 per theory if "cluster"). The remaining papers are presented to each GLM model and the probability returned by the model that the paper belongs to that theory is collected. The highest prediction value is selected as the "predicted" theory and stored. The results of this prediction are collected in each iteration. The number of iterations is specified in parameters above. The function is parallelized using the "foreach" and "doparallel" packages. "cl" controls the number of parallel processes; change to fit the number of cores in CPU. Each parallel iteration is one of the dimensions between minNumberOfDimensions and maxNumberOfDimensions # 

dimensionVec <- c(minNumberOfDimensions:maxNumberOfDimensions)

# text file to monitor the parallel foreach #

logger = create.logger()
logfile(logger) = 'monitor.log'
level(logger) = 'INFO'


cl <- makeCluster(7)
registerDoParallel(cl)


listResults <- foreach(dimension=dimensionVec, .verbose = T, .packages = "log4r") %dopar% { #parallelized foreach loop with each dimension. Stored in a list object containing the aggregate matrices of iterations controlled in repeat, for each number of dimensions used.
  
  resultsListModel <- lapply(c(1:repeats), matrix, nrow = 8, ncol = 8) #pre allocate the result list with the number of iterations selected in parameters.
  s = 1 #controller for the number of iterations.
  
  while(s <= repeats){ # repeats the process of training-prediction as specified in parameters.
    if(method == "free") {
      trainingSet = sample(1:nrow(my_svd),600) #not controlled training
      predictors <- bestPredictors
    }
    if(method == "cluster"){
      trainingSet <- c() # controlled training set for equal representation of each topic. 
      for(topic in topicList){
        trainingSet <- c(trainingSet, sample(which(my_catalog$topic==topic), round(length(which(my_catalog$topic==topic)) * 0.7)))
      }
      predictors <- bestPredictorsClusters    }
    if(source == "original" | source == "replication"){ #if the procedure is either predicting original for predicting original, or replication for predicting replication, set the trainingset as the rest of the papers.
    testSet = setdiff(1:nrow(my_svd),trainingSet)
    }
    if(source == "cross"){ #if using original to predict replication, trainingset is original dataset, and testset is replication set.
      trainingSet <- c(1:(nrow(my_svd)))
      testSet <- c(1:(nrow(cross_svd)))
    }

    predictionResults = c()
    
    for (topic in topicList) { #loop through the models of each theory
      trainingdata <- data.frame(my_svd[trainingSet,unlist(predictors[topic,1:dimension])]) # prepare training data of model by using the dimensions selected as the best predictors for each topic and the documents selected to be training.
      glmTopic = glm(my_catalog$topic[trainingSet]==topic~., data=trainingdata, family=binomial) #build the model of the topic.

      if(source == "original" | source == "replication"){ #if predicting inside the dataset
      testdata = data.frame(my_svd[testSet,unlist(predictors[topic,1:dimension])]) #prepare the data to be predicted
      predicted = predict.glm(glmTopic,newdata=testdata,type="response") #store the probability that the paper belongs to the theory being tested
      predictionResults = cbind(predictionResults,scale(predicted)) #add to a matrix and scale
      }
      if(source == "cross"){ #modify procedure above to account for original data being training and replication being test
        predicted = predict.glm(glmTopic,newdata=data.frame(cross_svd[testSet,unlist(predictors[topic,1:dimension])]),type="response")
        predictionResults = cbind(predictionResults,scale(predicted))
      }
    }
    # the predictions of each theory are aggregated in a matrix. each row of matrix is a document, each column is the probability that it belongs to that theory.
    
    predictionResults = data.frame(predictionResults)
    colnames(predictionResults) = topicList
    
    if(source == "original" | source == "replication"){
    predictionResults$topic = my_catalog$topic[testSet] #add a column with the correct theory of each of the papers in the testset
    predictionResults$predicted_topic = topicList[max.col(predictionResults[,1:8])] #add a column with the highest prediction of the models for that paper.
    resultTable <- (table(predictionResults$topic,predictionResults$predicted_topic) / as.vector(table(my_catalog$topic[testSet])))*100 #generates a frequency table of how many times each topic was predicted as each other topic. then, transforms into percentages.
    
    }
    
    if(source == "cross") { #same procedure, but for cross prediction.
      predictionResults$topic = cross_catalog$topic[testSet]
      predictionResults$predicted_topic = topicList[max.col(predictionResults[,1:8])]
      resultTable <- (table(predictionResults$topic,predictionResults$predicted_topic) / as.vector(table(cross_catalog$topic[testSet])))*100
    }
    
    resultsListModel[[s]] <- resultTable # store this iteration for final aggregation in list.
    s = s + 1
    if(s %% (round(repeats * 0.2
                   )) == 0){
      info(logger, paste("dimension ", dimension, ", ", "iteration number", s))
    }
  }
  
  #aggregate the results of the iterations#
  
  finalPredictionTable <- matrix(0, nrow = 8, ncol = 8) #pre-allocate final table
  
  #iterate through list of predictions to aggregate the results#
  
  for(matrix in resultsListModel) { # sums every result table resulting from the iterations.
    finalPredictionTable <- finalPredictionTable + matrix 
  }
  finalPredictionTable <- finalPredictionTable / length(resultsListModel) #divide by total number of iterations to aggregate
  return(finalPredictionTable) #return this table to be added to the list that foreach is constructing,
}

stopCluster(cl)

names(listResults) <- dimensionVec #name the objects in the list with the dimensions used in calculating them

dimEvMat <- matrix (0, nrow = length(dimensionVec), ncol = 9) # pre-allocate matrix for the evaluation of the effectiveness of models for each dimension
row.names(dimEvMat) <- as.character(dimensionVec) # name rows by dimension
colnames(dimEvMat) <- c(topicList, "mean") # name columns by each theory and designate a column to store the mean effectiveness
for(dimension in dimensionVec) { dimEvMat[as.character(dimension),] <- c(diag(listResults[[as.character(dimension)]]), mean(diag(listResults[[as.character(dimension)]])))} #fill list with the values of the diagonal of confusability matrices and its mean effectiveness.

  ###Plotting###

# scripts for producing the prediction related plots using ggplot2 #

## confusability matrix by dimension ##

dimension = 20 # parameter for choosing the number of dimensions to be used in the plot

topicMatrix <- listResults[[as.character(dimension)]] #extract the matrix of the chosen value of D
meltedResults <- melt(topicMatrix, varnames = c("Topic1", "Topic2"), value.name = "Percentage.Predicted")
heatmap <- ggplot(meltedResults, aes(y = Topic1, x = ordered(Topic2, levels = rev(sort(unique(Topic2)))))) + geom_tile(aes(fill = Percentage.Predicted)) + coord_equal() + scale_fill_gradient(limits = c(0, 100), low="white", high="seagreen", guide =  guide_colorbar(title = paste("% Predicted", "\n"))) + xlab("") + ylab("") + theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle=330, hjust=0.4, vjust = 0.7, size = 14)) + geom_text(aes(label = paste(round(Percentage.Predicted, 1), "%", sep = "")), colour = "gray25", size = 5)
print(heatmap)

## horizontal heatmap of effectiveness of each dimension by topic and mean #

meltedDimEv <- melt(dimEvMat[, 1:9], varnames = c("D", "topic"), value.name = "Effectiveness") #melted effectiveness for ggplot2
heatmap <- ggplot(meltedDimEv, aes(y = topic, x = ordered(D))) + geom_tile(aes(fill = Effectiveness), color = "white") + coord_equal() + scale_fill_gradient(limits = c(0, 100), low="white", high="seagreen") + xlab("") + ylab("") + theme(axis.text = element_text(size = 12)) + geom_text(aes(label = paste(round(Effectiveness, 0))), size = 4, colour = "gray25")
print(heatmap)


 # panel of effectiveness for each theory including mean #

meanResultsTable <- meltedDimEv[which(meltedDimEv$topic != "mean"),] #get mean out of data to generate panel of 8 theories
meanResultsTable <- as.data.frame(meanResultsTable) #long form for ggplot
plot <- ggplot(data = meanResultsTable, aes(x = D, y = Effectiveness, color = topic, group = topic)) + ylim (0, 100) + theme_gray() + geom_line(size = 1) + ylab("Mean Performance") + scale_colour_brewer(type = "qual", palette = "Paired", guide = F) + geom_point(size = 1, shape = 3) + facet_wrap(~topic, ncol = 2, nrow = 5)
print(plot)

# line plot of only mean effectiveness #

meanPerf <- data.frame("Mean.Effectiveness" = rowMeans(dimEvMat), "D" = dimensionVec) #mean performance melted dataframe
plot <- ggplot(data = meanPerf, aes(y = Mean.Effectiveness, x = D)) + theme_gray() + geom_line(size = 1.5, color = "seagreen") + geom_point(size = 1.5, shape = 3, color = "seagreen") + scale_y_continuous(limits = c(0, 100)) + labs(y = "Mean Performance (%)")
print(plot)

####LINEAR MODELS WITH COSINE MATRICES####

## this section generates a comparison of the similarity between and within theories by using the cosine generated by pairs of abstracts in as a vector space representation in the space of the SVD. Similarity data is then used to measure self similarity and other-similarity with linear models. ##

## Cosine matrices ##

originalCosines <- compareTheories(documentLoadings, catalog) # has two elements: [[1]] list of a 8 matrices, one for each theory, with mean distance of that theory with each other theory. (2) list of one matrix for each dimension with theories as rows and all papers as column. cells show the mean distance of that theory with that paper.
originalAverages <- originalCosines[[2]] #extracts the average distance of theory by paper.
replicationCosines <- compareTheories(repDocumentLoadings, repCatalog) #same procedure for replication
replicationAverages <- replicationCosines[[2]]

#plots of mean distance of theories with other theories for dimension

plotTopicDiff("bayesian", originalCosines)
plotTopicDiff("symbolic", originalCosines)
plotTopicDiff("connectionism", originalCosines)
plotTopicDiff("embodied", originalCosines)
plotTopicDiff("distributed", originalCosines)
plotTopicDiff("enactive", originalCosines)
plotTopicDiff("dynamical", originalCosines)
plotTopicDiff("ecological", originalCosines)

##PARAMETERS##

# these parameters control the models of self-similarity and other-similarity. "testingSelf" defines the source to be used (original data or replication data). "dimensionToTest" controls the number of dimensions used in the analysis (D) #

testingSelf <- "original" #original or replication
dimension <- 5  #specify the dimension to be tested

##OBJECTS##

#loads the data based on the parameters provided#

if(testingSelf == "original"){
  avgList <- originalAverages
  my_cat <- catalog
}
if(testingSelf == "replication"){
  avgList <- replicationAverages
  my_cat <- repCatalog
}

### SELF-SIMILARITY ###

# measure of self-similarity of each theory. takes the average distance of each theory with each paper of that same theory and then models that relation with linear model #

allDat = c()

for (topic in topicList) { #populate the mean cosine of theory with papers of that theory data
  dat = avgList[[dimension]][which(topic==topicList),my_cat$topic==topic]
  allDat = rbind(allDat,data.frame(topic=topic,cosine=dat))
}

#build a linear model of mean distance of theory with papers of that theory and the theory#

lmObject = lm(cosine~topic,data=allDat)
summary(lmObject)

#visualize with boxplot#

boxplot <- ggplot(allDat, aes(x = topic, y = cosine)) + geom_boxplot(fill = "#fdf6e3", colour = "#2aa198", outlier.color = "#2aa198") + scale_x_discrete(name = "Theory") + scale_y_continuous(name = "Mean self-cosine", breaks = seq(-0.1, 1, .10), limits = c(-0.1, 1)) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) + labs(title = paste("Self-similarity of different theories of dimension number", dimension), subtitle = "Cosines of members of a theory with other members of that theory as predicted by theory membership")
boxplot

##EXTERNAL DIFFERENCE##

# measure of the similarity between the theory and the difference of similarity between the mean distance with its own papers, and the mean distance with papers of the most similar theory other than itself ("closest neighbor"). higher values indicate more difference between the theory and the most similar theory to it.#

allDatMax = c() # data of similarity with closest neighbor

for (topic in topicList) { 
  selfdat = avgList[[dimension]][which(topic==topicList),my_cat$topic==topic] #collect the data of similarity with its own papers
  otherdat = avgList[[dimension]][which(topic != topicList),my_cat$topic==topic] #collect the mean distance of other theories with each of the papers of the theory being looped in the for loop
  maxdat = apply(otherdat,2,function(x) { return(max(x))}) #collect the highest value
  maxdat = selfdat - maxdat #store the difference between the theory the paper belongs to and the most similar other theory
  allDatMax = rbind(allDatMax,data.frame(topic=topic,cosine=maxdat)) #matrix with the differences for each theory.
}

lmObjectMax = lm(cosine~topic,data=allDatMax) #linear model for other-similarity
summary(lmObjectMax) #summary

#visualize with boxplot#

boxplotMax <- ggplot(allDatMax, aes(x = topic, y = cosine)) + geom_boxplot(fill = "#fdf6e3", colour = "#2aa198", outlier.color = "#2aa198") + scale_x_discrete(name = "Theory") + scale_y_continuous(name = "Mean distance") + labs(title = paste("Distance from nearest theory of different theories", testingSelf), subtitle = paste("Difference between mean distance for own theory and mean distance from the theory with the highest cosine predicted by theory of", testingSelf)) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) 
print(boxplotMax)

####CLUSTER ANALYSIS USING COSINE DATA####

#this section uses the similarity data of cosines gathered in the previous section to visualize the clustering of each theories with a heatmap of the similarity matrix and a hierarchical cluster analysis#

##PARAMETERS##

#parameters that control the value of D ("dimension") and the source of the data ("dendrogramMode")#

dimension <- 20 #change dimension being considered
dendrogramMode <- "original" #original or replication

##OBJECTS##

#loads objects based on parameters#

if(dendrogramMode == "original"){
  avgMatrix <- originalCosines[[2]][[dimension]]
  my_cat <- catalog
  my_svd <- documentLoadings
}

if(dendrogramMode == "replication"){
  avgMatrix <- replicationCosines[[2]][[dimension]]
  my_cat <- repCatalog
  my_svd <- repDocumentLoadings
}

## similarity matrix of theories. higher value is less distance. ##

#generates a similarity matrix of cosine data to visualize in a heatmap #

topicListReordered = c("bayesian", "connectionism", "symbolic", "distributed", "dynamical", "enactive", "ecological", "embodied") #reorder topics to better visualize if intuitive cluster shows

simMatrix = matrix(0,nrow=8,ncol=8) #allocates similarity matrix 
row.names(simMatrix) = topicListReordered #name dimensions of similarity matrix (theory x theory)
colnames(simMatrix) = topicListReordered

for (topic in topicListReordered) { #loops through first topic (row)
  for(i in 1:8){ #loops through second topic (column)
    simMatrix[topic,topicListReordered[i]] = mean(avgMatrix[topic, which(my_cat$topic==topicListReordered[i])]) #gather the mean similarity between theory 1 and the papers of theory 2
  }
}


# visualize similarity matrix with a heatmap #

meltedDistMatrix <- melt(simMatrix, varnames = c("Topic1", "Topic2"), value.name = "Closeness") #long form of simMatrix for ggplot
heatmap <- ggplot(meltedDistMatrix, aes(x = Topic1, y = ordered(Topic2, levels = rev(sort(unique(Topic2)))))) + geom_tile(aes(fill = Closeness), colour = "white") + coord_equal() +  scale_fill_gradient(low="white", high="seagreen", limits = c(0, 1), guide =  guide_colorbar(title = paste("Cosine", "\n"))) + xlab("") + ylab("") + theme(axis.text = element_text(size = 12), axis.text.x = element_text(angle = 330, hjust = 0.4, vjust = 0.7))
print(heatmap)

# hierarchical cluster analysis of the distance matrix #

topic_hclust <- hclust(dist(simMatrix, upper = T), method = "average" ) #applies hclust algorithm to the similarity matrix

#visualize topic_hclust with a dendrogram and mark clusters using cutreedynamic package (https://www.rdocumentation.org/packages/dynamicTreeCut/versions/1.63-1/topics/cutreeDynamic)#

par(lwd = 2, cex.axis = 1.2, las = 1) #graphical parameters for the axes and lines of dendrogram
print(hclustplot(topic_hclust, colors = labels2colors(cutreeDynamic(topic_hclust, minClusterSize = 1, method = "hybrid", deepSplit = 0, distM = as.matrix(dist(simMatrix)), pamStage = T), colorSeq = c("#2E8B57", "#FF2052", "#8b2e62")), fillbox = T, las = 0, cex = 1.2, mar = c(3, 3, 2, 0.5), font = 2)) #produces a dendrogram and marks the clusters with dynamictreecut. the minimum cluster size is set to 1, the parameter controlling the stringency of the clustering is set to the default (0). we used the "hybrid" method which takes both the dendrogram and a distance matrix to generate clusters. PAM was not used. We provided the colors (colorSeq) of the fillbox marking the clusters for up to 3 clusters.   

#rotation and dimension inspection#

##building a matrix with dimensions selected based on prediction performance##

termVarimax <- varimax(termLoadings)
termLoadingsVarimax <- unclass(termVarimax$loadings)
documentLoadingsVarimax <- documentLoadings %*% termVarimax$rotmat

bestPredictorsVarimax <- c()
predictorRatingsVarimax <- c() #store ratings to compare positive versus negative in term loading inspection
for (topic in topicList) { 
  glmOutput = glm(my_catalog$topic==topic~.,data=data.frame(documentLoadingsVarimax[,1:maxNumberOfDimensions]),family=binomial)
  bestPredictorsVarimax = rbind(bestPredictorsVarimax,data.frame(t(sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),ind=T,decreasing=T)$ix[1:maxNumberOfDimensions]))) #TODO: fix this. 2 is here to skip intercept.
  predictorRatingsVarimax <- rbind(predictorRatingsVarimax, glmOutput$coefficients[2:length(glmOutput$coefficients)][sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]), ind = T, decreasing = T)$ix])
}
row.names(bestPredictorsVarimax) <- topicList
row.names(predictorRatingsVarimax) <- topicList

wordList <- list()

for(n in 1:length(topicList)){ # generate a data frame for each theory
  topic <- topicList[n]
  predictors <- unlist(bestPredictorsVarimax[topic,])
  ratings <- unlist(predictorRatingsVarimax[topic,])
  ratingsValence <- ratings > 0 # boolean for sign of the predictor. true if > 0, false if < 0
  words <- c()
  i = 1
  for(i in c(1:length(predictors))){
    dimension <- predictors[i]
    words <- c(words, toString(head(names(sort(termLoadingsVarimax[,dimension], decreasing = T)), 150)))
  }
  nameOfDF <- paste(topic, "Meanings", sep = "")
  assign(nameOfDF, data.frame(predictors, ratings, ratingsValence, words))
  wordList[[n]] <- get(x = nameOfDF)
}

names(wordList) <- topicList

#wordclouds#

for(topic in topicList){
  positiveWords <- c()
  negativeWords <- c()
  meanings <- wordList[[topic]]
  for(i in 1:10){
    if(meanings$ratingsValence[i]){
    positiveWords <- c(positiveWords, unlist(strsplit(as.character(meanings$words[i]), ',')))
    } else{
      negativeWords <- c(negativeWords, unlist(strsplit(as.character(meanings$words[i]), ',')))
      }
    }
  positivedf <- data.frame(names(sort(table(positiveWords), decreasing = T)), as.numeric(sort(table(positiveWords), decreasing = T)))
  negativedf <- data.frame( names(sort(table(negativeWords), decreasing = T)), as.numeric(sort(table(negativeWords), decreasing = T)))
  colnames(positivedf) <- c("word","freq")
  colnames(negativedf) <- c("word", "freq")
  png(filename = paste(topic, "PositiveCloud.png", sep = ""), units = "px", width = 3000, height = 2000)
  wordcloud(positivedf$word, positivedf$freq, scale = c(7, .5), random.order = F, rot.per = 0, use.r.layout = F, min.freq = 2)
  dev.off()
  png(filename = paste(topic, "NegativeCloud.png", sep = ""), units = "px", width = 3000, height = 2000)
  wordcloud(negativedf$word, negativedf$freq, scale = c(7, .5), random.order = F, rot.per = 0, use.r.layout = F, min.freq = 2)
  dev.off()
}

