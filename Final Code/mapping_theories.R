#### libraries ####

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

# paralellization of prediction machine #
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
    if(n == 1){dt <- cbind(database, 0)} #if it has only one dimension, then add a column of 0s to make cosineGen work
    database <- cosineGen(dt) #produces a x by x matrix of cosines between each paper.
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
  plot <- ggplot(data = workTable, aes(x = dimension, y = cosine, color = topic, group = topic)) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) + labs(title = paste("Mean cosine of", x, "papers with other theories and itself across dimensions")) + geom_line() + scale_colour_solarized("red") + geom_point(size = 0.7, shape = 3) + guides(colour = guide_legend(override.aes = list(size=3)))
  print(plot)
  } #function for plotting the mean distance of every topic with all other topics. "topic" is one of the topics of topicList; "resultsList" is the object that compareTheories() returns.

####RAW DATA AND WORD FREQUENCY####

##ORIGINAL##

#loads files and catalog for original#

freqMatrix <- as.matrix(read.table('document_by_term.txt', sep='\t', header = T))[, -1] #loads the DBT minus one column. 
row.names(freqMatrix) <- c(1:nrow(freqMatrix)) #row names with docID
freqMatrix <- freqMatrix[, apply(freqMatrix, 2, function(x){sum(x==0) < 1047})] #removes columns with words that appear in fewer than 5 documents.

freqMatrix <- freqMatrix[which(rowSums(freqMatrix) > 0), ] #eliminates documents with 0 terms after cleanup of terminology
catalog <- read.table('catalog.txt', stringsAsFactors = F, sep = '\t', fill = T) #loads catalog
catalog <- catalog[row.names(freqMatrix), ] #catalog also has row names as docID
colnames(catalog) = c('id','topic','year','authors','title','journal','abstract') #variable names for catalog
topicList <- unique(catalog$topic) #list of theories for analysis.

##REPLICATION##

#the same, but for replication documents#

repFreqMatrix <- as.matrix(read.table('rep_document_by_term.txt', sep='\t', header = T, quote = ""))[, -1] 
row.names(repFreqMatrix) <- c(1:nrow(repFreqMatrix))
repFreqMatrix <- repFreqMatrix[, apply(repFreqMatrix, 2, function(x){sum(x==0) < 1009})] #removes columns with words that appear in fewer than 5 documents.
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

wholeMacaroni = svd(cleanData, nu = 150, nv = 150) #partial SVD of 150 dimensions.
row.names(wholeMacaroni$u) <- row.names(cleanData) #puts docIDS in the matrices resulting from SVD.
row.names(wholeMacaroni$v) <- colnames(cleanData)

#REPLICATION#

repofWholeMacaroni = svd(repCleanData, nu = 150, nv = 150)
row.names(repofWholeMacaroni$u) <- row.names(repCleanData)
row.names(repofWholeMacaroni$v) <- colnames(repCleanData)

  
####GLM Models####

# This part of the script has the GLM models of each theory that attempt to predict the theory belonging of each paper #

## PARAMETERS ##

# Set parameters for the prediction. Min and max number of dimensions are used to control the number of dimensions that are to be used in the construction of the models. The procedure loops through the dimensions resulting from the SVD. Starts at minNumberOfDimensions (default: 3), stops at maxNumberOfDimensions (default: 50). "Method" refers to the data used to build and train the models. With "free", dimensions are selected for how well they predict theory belonging against every theory. With "cluster", the training is stratified to the "most similar" theories; e.g. "computational" is built using the dimensions that best predict computational papers when compared to 'bayesian' and 'connectionist'. "Repeats" is the number of iterations of the predicting process. "Source" controls which data set is to be used: "original" uses the original dataset, "replication" uses the replication data, and "cross" uses the projection of the replication data into the SVD space of the original dataset to predict their theories with the models built with the original dataset. #

minNumberOfDimensions <- 3 #lower boundary of D
maxNumberOfDimensions <- 50 # upper boundary of D
repeats <- 1000 # how many repetitions of prediction should be averaged?
method <- "cluster" # "cluster" or "free".
source <- "original" #"original" or "replication", "cross".

##OBJECTS##

if(source == "original"){
  my_catalog <- catalog
  my_svd <- wholeMacaroni$u
}
if(source == "replication"){
  my_catalog <- repCatalog
  my_svd <- repofWholeMacaroni$u
}
if(source == "cross"){
  my_catalog <- catalog
  cross_catalog <- repCatalog
  my_svd <- wholeMacaroni$u
  #create matrix for projection and populate#
  replicationProjection <- matrix(0, nrow = 964, ncol = 3611)
  row.names(replicationProjection) <- row.names(repFreqMatrix)
  colnames(replicationProjection) <- colnames(freqMatrix)
  sharedWords <- colnames(repFreqMatrix)[which(colnames(repFreqMatrix) %in% colnames(freqMatrix))]
  replicationProjection[,sharedWords] <- repFreqMatrix[,sharedWords]
  #project replica matrix onto previous svd space#
  cross_svd <- replicationProjection %*% ginv(diag(wholeMacaroni$d[1:150]) %*% t(wholeMacaroni$v))
}


### TOPIC BEST PREDICTORS ###

##FREE FOR ALL##

bestPredictors <- c()
for (topic in topicList) {
  glmOutput = glm(my_catalog$topic==topic~.,data=data.frame(my_svd[,1:80]),family=binomial)
  bestPredictors = rbind(bestPredictors,data.frame(topic,(t(sort(glmOutput$coefficients,ind=T,decreasing=T)$ix[1:maxNumberOfDimensions]))))
}
row.names(bestPredictors) <- topicList

##STRATIFIED SAMPLING##

topicListClassic <- c(topicList[1], topicList[2], topicList[8]) # topic list of traditional approaches (comp, bayes, conn)
topicListAlt <- setdiff(topicList, topicListClassic)
catalogClassic <- my_catalog[which(my_catalog$topic %in% topicListClassic),]

catalogAlt <- my_catalog[which(my_catalog$topic %in% topicListAlt),]

bestPredictorsClusters <- c()

for (topic in topicListClassic) {
  glmOutput = glm(catalogClassic$topic==topic~.,data=data.frame(my_svd[which(my_catalog$topic %in% topicListClassic, arr.ind = T),1:80]),family=binomial)
  bestPredictorsClusters = rbind(bestPredictorsClusters,data.frame(topic,(t(sort(glmOutput$coefficients,ind=T,decreasing=T)$ix[1:maxNumberOfDimensions]))))
}

for (topic in topicListAlt) {
  glmOutput = glm(catalogAlt$topic==topic~.,data=data.frame(my_svd[which(my_catalog$topic %in% topicListAlt, arr.ind = T),1:80]),family=binomial)
  bestPredictorsClusters = rbind(bestPredictorsClusters,data.frame(topic,(t(sort(glmOutput$coefficients,ind=T,decreasing=T)$ix[1:maxNumberOfDimensions]))))
}

row.names(bestPredictorsClusters) <- c(topicListClassic, topicListAlt)


### EXECUTION OF MODEL ###

dimensionVec <- c(minNumberOfDimensions:maxNumberOfDimensions)

cl <- makeCluster(8)
registerDoParallel(cl)

listResults <- foreach(dimension=dimensionVec, .verbose = T) %dopar% { #dirty for loop for evaluating self identification by number of dimensionsx
  resultsListModel <- lapply(c(1:repeats), matrix, nrow = 8, ncol = 8)
  s = 1
  while(s <= repeats){
    #training and test set#
    if(method == "free") {
      trainingSet = sample(1:nrow(my_svd),600) #not controlled training
      predictors <- bestPredictors
    }
    if(method == "cluster"){
      trainingSet <- c() # controlled training set for equal representation of each topic. Stratified sampling.
      for(topic in topicList){
        trainingSet <- c(trainingSet, sample(which(my_catalog$topic==topic), 60)) #takes indices of each topic
      }
      predictors <- bestPredictorsClusters    }
    if(source == "original" | source == "replication"){
    testSet = setdiff(1:nrow(my_svd),trainingSet)
    }
    if(source == "cross"){
      trainingSet <- c(1:(nrow(my_svd)))
      testSet <- c(1:(nrow(cross_svd)))
    }
    
    #GLM model with predictors, and prediction of every paper in testdata by each model of each topic#
    
    predictionResults = c()
    for (topic in topicList) {
      trainingdata <- data.frame(my_svd[trainingSet,unlist(predictors[topic,2:dimension])])
      glmTopic = glm(my_catalog$topic[trainingSet]==topic~., data=trainingdata, family=binomial)
      #glmTopic = glm(my_catalog$topic[trainingSet]==topic~., data=data.frame(my_svd[trainingSet,unlist(predictors[topic,2:dimension])]), family=binomial)
      if(source == "original" | source == "replication"){
      testdata = data.frame(my_svd[testSet,unlist(predictors[topic,2:dimension])])
      
      predicted = predict.glm(glmTopic,newdata=testdata,type="response")
      #predicted = predict.glm(glmTopic,newdata=data.frame(my_svd[testSet,unlist(predictors[topic,2:dimension])]),type="response")
      predictionResults = cbind(predictionResults,scale(predicted))
      }
      if(source == "cross"){
        predicted = predict.glm(glmTopic,newdata=data.frame(cross_svd[testSet,unlist(predictors[topic,2:dimension])]),type="response")
        predictionResults = cbind(predictionResults,scale(predicted))
      }
    }
    
    #aggregate predictions for each topic into a table#
    
    predictionResults = data.frame(predictionResults)
    colnames(predictionResults) = topicList
    if(source == "original" | source == "replication"){
    predictionResults$topic = my_catalog$topic[testSet]
    predictionResults$predicted_topic = topicList[max.col(predictionResults[,1:8])]
    resultTable <- t(t(table(predictionResults$topic,predictionResults$predicted_topic) / as.vector(table(my_catalog$topic[testSet])))*100)
    }
    if(source == "cross") {
      predictionResults$topic = cross_catalog$topic[testSet]
      predictionResults$predicted_topic = topicList[max.col(predictionResults[,1:8])]
      resultTable <- t(t(table(predictionResults$topic,predictionResults$predicted_topic) / as.vector(table(cross_catalog$topic[testSet])))*100)
    }
    
    #allocate aggregate data of prediction (mean effectiveness by topic)#
    
    finalRunTable <- matrix(0, nrow = 8, ncol = 8) ## to avoid dimension errors when summing, when topic is not predicted
    row.names(finalRunTable) <- topicList
    colnames(finalRunTable) <- topicList
    for(topic in colnames(resultTable)){finalRunTable[, topic] = resultTable[, topic]}
    resultsListModel[[s]] <- finalRunTable
    if(s %% 100 == 0) {print(paste("repetition number", s))}
    s = s + 1
    }
  finalPredictionTable <- matrix(0, nrow = 8, ncol = 8)
  
  #iterate through list of predictions to aggregate the results#
  
  for(matrix in resultsListModel) {
    finalPredictionTable <- finalPredictionTable + matrix
    }
    finalPredictionTable <- finalPredictionTable / length(resultsListModel)
    #for(topic in topicList){#populate results looping through topics
      #predResultsxDim[topic,dimension - 2] <- finalPredictionTable[topic,topic]
    #}
  return(finalPredictionTable)
} #for loop for evaluating multiple dimensions etc. needs to be cleaned up.

stopCluster(cl)

names(listResults) <- dimensionVec

###LOADING SAVED RESULTS###

#LOAD ONE AT A TIME. LOADS AS OBJECT "LISTRESULTS", WITH A LIST OF 50 DIMENSIONS WITH MATRIX AGGREGATING 1000 PASSES FOR EACH.#

#load(file = "PredictionResultsOriginal50D.RData")
#load(file = "PredictionResultsReplication50D.RData")
#load(file = "PredictionResultsCross50D.RData")
load(file = "PredictionResultsOriginal50DNullHypothesis.RData")

dimEvMat <- matrix (0, nrow = length(dimensionVec), ncol = 9) # matrix for the evaluation of mean effectiveness of each dimension
row.names(dimEvMat) <- as.character(dimensionVec)
colnames(dimEvMat) <- c(topicList, "mean")

pdf(file = "confmatrixD20NH.pdf", width = 10, height = 8)

for(dimension in dimensionVec){
  topicMatrix <- listResults[[as.character(dimension)]]
  dimEvMat[as.character(dimension),] <- c(diag(topicMatrix), mean(diag(topicMatrix)))
  #plot confusability matrix each dimension#
  meltedResults <- melt(topicMatrix, varnames = c("Topic1", "Topic2"), value.name = "Percentage.Predicted")
  heatmap <- ggplot(meltedResults, aes(y = Topic1, x = ordered(Topic2, levels = rev(sort(unique(Topic2)))))) + geom_tile(aes(fill = Percentage.Predicted)) + coord_equal() + scale_fill_gradient(limits = c(0, 100), low="white", high="seagreen", guide =  guide_colorbar(title = paste("% Predicted", "\n"))) + xlab("") + ylab("") + theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle=330, hjust=0.4, vjust = 0.7, size = 14)) + geom_text(aes(label = paste(round(Percentage.Predicted, 1), "%", sep = "")), colour = "gray25", size = 5)
  print(heatmap)
}
dev.off()

meltedDimEv <- melt(dimEvMat[, 1:8], varnames = c("D", "topic"), value.name = "Effectiveness")
heatmap <- ggplot(meltedDimEv, aes(y = topic, x = ordered(D))) + geom_tile(aes(fill = Effectiveness), color = "white") + coord_equal() + scale_fill_gradient(limits = c(0, 100), low="white", high="seagreen") + xlab("") + ylab("") + theme(axis.text = element_text(size = 12)) + geom_text(aes(label = paste(round(Effectiveness, 0))), size = 4, colour = "gray25")
print(heatmap)
ggsave(filename = "effectivenessheatmap50DOGNH.pdf", device = "pdf", width = 35, height = 7.5, units = "in", dpi = 1200)

meanResultsTable <- as.data.frame(meltedDimEv) #long form for ggplot
plot <- ggplot(data = meanResultsTable, aes(x = D, y = Effectiveness, color = topic, group = topic)) + ylim (0, 50) + theme_gray() + geom_line(size = 1) + scale_colour_brewer(type = "qual", palette = "Set2", guide = F) + geom_point(size = 0.5, shape = 3) + facet_wrap(~topic, ncol = 2, nrow = 5)
print(plot)
ggsave(filename = "effectivenesslineplotallD-OG-NullHyp.pdf", device = "pdf", width = 174, height = 174, units = "mm", dpi = 1200)

meanPerf <- data.frame("Mean.Effectiveness" = rowMeans(dimEvMat), "D" = dimensionVec) #mean performance
plot <- ggplot(data = meanPerf, aes(y = Mean.Effectiveness, x = D)) + theme_gray() + geom_line(size = 1.5, color = "seagreen") + geom_point(size = 1.5, shape = 3, color = "seagreen") + scale_y_continuous(limits = c(0, 100)) + labs(y = "Mean Effectiveness (%)")
print(plot)
ggsave(filename = "meaneffectiveness50D-OG-NullHyp.pdf", device = "pdf", width = 174, height = 70, units = "mm", dpi = 1200)


####linear models with cosine matrix single dimension####

##ALL OBJECTS##

originalCosines <- compareTheories(wholeMacaroni$u, catalog) # has two elements: all theories (1), averages (2)
originalAverages <- originalCosines[[2]]
replicationCosines <- compareTheoriesVect(repofWholeMacaroni$u, repCatalog)  # has two elements: all theories (1), averages (2)
replicationAverages <- replicationCosines[[2]]

plotTopicDiff("bayesian", originalCosines)
plotTopicDiff("symbolic", originalCosines)

##PARAMETERS##

testingSelf <- "original" #original or replication
dimensionToTest <- 10 #specify the dimension to be tested

##OBJECTS##

if(testingSelf == "original"){
  avgList <- originalAverages
  my_cat <- catalog
}
if(testingSelf == "replication"){
  avgList <- replicationAverages
  my_cat <- repCatalog
}

### SELF-SIMILARITY ###

allDat = c()
#data collection from averages#
for (topic in topicList) {
  dat = avgList[[dimensionToTest]][which(topic==topicList),my_cat$topic==topic]
  allDat = rbind(allDat,data.frame(topic=topic,cosine=dat))
}
#model#
lmObject = lm(cosine~topic,data=allDat)
summary(lmObject)
#plot#
boxplot <- ggplot(allDat, aes(x = topic, y = cosine)) + geom_boxplot(fill = "#fdf6e3", colour = "#2aa198", outlier.color = "#2aa198") + scale_x_discrete(name = "Theory") + scale_y_continuous(name = "Mean self-cosine", breaks = seq(-0.1, 1, .10), limits = c(-0.1, 1)) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) + labs(title = paste("Self-similarity of different theories of dimension number", dimensionToTest), subtitle = "Cosines of members of a theory with other members of that theory as predicted by theory membership")
boxplot
ggsave(filename = "BoxplotSelfSimilarityD10NullHyp.pdf", device = "pdf", width = 174, height = 70, units = "mm", dpi = 1200)

##EXTERNAL DIFFERENCE##

allDatMax = c() # closest neighbor
allDatMean = c() # mean distance

for (topic in topicList) {
  selfdat = avgList[[10]][which(topic==topicList),my_cat$topic==topic]
  # the AVERAGE across category similarity for this topic
  otherdat = avgList[[10]][which(topic != topicList),my_cat$topic==topic]
  maxdat = apply(otherdat,2,function(x) { return(max(x))})
  maxdat = selfdat - maxdat
  meandat = apply(otherdat,2,function(x) { return(mean(x))})
  meandat = selfdat - meandat
  allDatMax = rbind(allDatMax,data.frame(topic=topic,cosine=maxdat)) #matrix of what the maximum distance for each paper of each topic is from the most similar theory (another measure of self similarity). If higher, topic is farther from nearest theory.
  allDatMean = rbind(allDatMean, data.frame(topic=topic, cosine=meandat))
}
lmObjectMax = lm(cosine~topic,data=allDatMax) #linear model for max
summary(lmObjectMax) #summary
lmObjectMean = lm(cosine~topic,data=allDatMean) #linear model for mean
summary(lmObjectMean)

#boxplot of max
boxplotMax <- ggplot(allDatMax, aes(x = topic, y = cosine)) + geom_boxplot(fill = "#fdf6e3", colour = "#2aa198", outlier.color = "#2aa198") + scale_x_discrete(name = "Theory") + scale_y_continuous(name = "Mean distance") + labs(title = paste("Distance from nearest theory of different theories", testingSelf), subtitle = paste("Difference between mean distance for own theory and mean distance from the theory with the highest cosine predicted by theory of", testingSelf)) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) 
print(boxplotMax)
ggsave(filename = "BoxplotMaxOtherSimD10NullHyp.pdf", device = "pdf", width = 174, height = 70, units = "mm", dpi = 1200)


#boxplot of mean
boxplotMean <- ggplot(allDatMean, aes(x = topic, y = cosine)) + geom_boxplot(fill = "#fdf6e3", colour = "#2aa198", outlier.color = "#2aa198") + scale_x_discrete(name = "Theory") + scale_y_continuous(name = "Mean distance", labels = comma) + labs(title = paste("Mean distance from different theories", testingSelf), subtitle = paste("Difference between mean distance for own theory and mean distance from mean distance from other theories predicted by theory of"), testingSelf) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) 
print(boxplotMean)
ggsave(filename = "BoxplotMeanOtherD10NullHyp.pdf", device = "pdf", width = 174, height = 70, units = "mm", dpi = 1200)



####DENDROGRAMS####

##PARAMETERS##
dimension <- 10 #change dimension being considered
dendrogramMode <- "original" #original or replication

##OBJECTS##

if(dendrogramMode == "original"){
  avgMatrix <- originalCosines[[2]][[dimension]]
  my_cat <- catalog
  my_svd <- wholeMacaroni$u
}

if(dendrogramMode == "replication"){
  avgMatrix <- replicationCosines[[2]][[dimension]]
  my_cat <- repCatalog
  my_svd <- repofWholeMacaroni$u
}

##TOPICS##
#uses matrices of AVERAGE cosine of each topic with respective paper.

topicListReordered = c("bayesian", "connectionism", "symbolic", "distributed", "dynamical", "enactive", "ecological", "embodied")

dists = matrix(0,nrow=8,ncol=8) #allocates dist matrix
row.names(dists) = topicListReordered
colnames(dists) = topicListReordered

for (topic in topicListReordered) { #loops through topics. topic1 in comments
  for(i in 1:8){ #second loop through topics. topic2 in comments
    dists[topic,topicListReordered[i]] = mean(avgMatrix[topic, which(my_cat$topic==topicListReordered[i])]) #mean avg of topic 1 with every topic2 spits out avg cosine topicxtopic matrix.
  }
}

topic_hclust <- hclust(dist(dists, upper = T), method = "average" ) #hclust algorithm application.

#topicxtopic heatmap#
meltedDists <- melt(dists, varnames = c("Topic1", "Topic2"), value.name = "Closeness")
heatmap <- ggplot(meltedDists, aes(x = Topic1, y = ordered(Topic2, levels = rev(sort(unique(Topic2)))))) + geom_tile(aes(fill = Closeness), colour = "white") + coord_equal() +  scale_fill_gradient(low="white", high="seagreen", limits = c(0, 1), guide =  guide_colorbar(title = paste("Cosine", "\n"))) + xlab("") + ylab("") + theme(axis.text = element_text(size = 12), axis.text.x = element_text(angle = 330, hjust = 0.4, vjust = 0.7))
print(heatmap)
ggsave(filename = "similarity-heatmap.pdf", device = "pdf", width = 174, height = 130, units = "mm", dpi = 1200)


#hclustplot(topic_hclust, colors = labels2colors(cutreeDynamic(topic_hclust, minClusterSize = 1, method = "hybrid", deepSplit = 0, distM = as.matrix(dist(dists)), pamStage = F)), fillbox = T, main = paste("Dendrogram of topics for dimensions 1 to", dimension))
pdf(file = "dendrogramD10NH.pdf", width = 10, height = 7)
par(lwd = 2, cex.axis = 1.2, las = 1)
print(hclustplot(topic_hclust, colors = labels2colors(cutreeDynamic(topic_hclust, minClusterSize = 1, method = "hybrid", deepSplit = 0, distM = as.matrix(dist(dists)), pamStage = F), colorSeq = c("#2E8B57", "#FF2052", 	"#8b2e62")), fillbox = T, las = 0, cex = 1.2, mar = c(3, 3, 2, 0.5), font = 2))

dev.off()
