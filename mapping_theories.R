# Libraries ---------------------------------------------------------------

# Plotting
require("reshape2")
require("ggplot2")
require("ggthemes")
require("scales")

# Pendrogram tree cutting from https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/BranchCutting/
require("dynamicTreeCut")
require("bio3d")
require("moduleColor")

# Paralellization of prediction machine and logging of silent parallel cores.
require("log4r")
require("foreach")
require("doParallel")

# Linear algebra (inverse matrices and orthonormalization)
require('MASS')
require('far')

# Word clouds
require('wordcloud')

# Functions ---------------------------------------------------------------

generateCosineMatrix <- function(matrix) {
  lengthOfVector <- sqrt(rowSums(matrix * matrix))
  tcrossprod(matrix) / (lengthOfVector %o% lengthOfVector)
}  #Function that generates the cosine between each row of a matrix.

calculateEntropy <- function(matrix) {
  temporaryMatrix <- t(matrix) #transposes for division
  a <-
    temporaryMatrix / rowSums(temporaryMatrix) #generates a probability matrix
  b <-
    1 + ((rowSums(a * log2(a), na.rm = T)) / log2(dim(matrix)[1])) #calculate entropy (1 + (sum of probability times log2 probability, divided by the total number of documents)).
  temporaryMatrix <-
    log(temporaryMatrix[which(b > 0, arr.ind = T), ] + 1) #log normalizes frequency matrix and deletes 0 entropy ones.
  temporaryMatrix <-
    temporaryMatrix * b[b > 0] #weight log normalized matrix by multiplying terms with entropy higher than 0 by its entropy.
  return(t(temporaryMatrix)) #returns original, non-transposed matrix.
} #calculates the entropy of each term of a matrix. Uses formula in Martin & Berry, 2007, "Mathematical foundations behind latent semantic analysis".

compareTheories <- function(matrix, cat) {
  ##doc: takes matrix as the result of an svd(u). populates a pre-allocated list of every theory in external theoryList (not generalized yet) with matrices of an "index" of similarity using each dimension (column) in matrix for each theory.
  listOfResults <-
    lapply(lapply(
      1:8,
      matrix,
      data = 0,
      nrow = 8,
      ncol = dim(matrix)[2]
    ), function(x) {
      row.names(x) <- theoryList
      return(x)
    }) #pre-allocation of list
  cosineAvgsList <-
    lapply(
      1:8,
      matrix,
      data = 0,
      nrow = 8,
      ncol = dim(matrix)[1]
    ) #pre-allocation of second list
  names(listOfResults) <-
    theoryList #name for easy accessing theories
  logicalMatrixTheories <-
    matrix(
      FALSE,
      nrow = dim(matrix)[1],
      ncol = 8,
      dimnames = list(1:dim(matrix)[1], theoryList)
    ) #pre-allocates logical matrix
  for (theory in theoryList) {
    logicalMatrixTheories[, theory] <-
      cat[, 2] == theory
  } #populates logical matrix with a logical mask reflecting catalog$id
  docsByTop <-
    colSums(logicalMatrixTheories) #number of documents for each theory
  n <- 1 #counter for each dimension
  while (n <= dim(matrix)[2]) {
    #loops through dimensions
    database <- matrix[, 1:n] #slices dimensions
    if (n == 1) {
      database <-
        cbind(database, 0)
    } #if it has only one dimension, then add a column of 0s to make cosineGen work
    database <-
      generateCosineMatrix(database) #produces a x by x matrix of cosines between each paper.
    database[is.na(database)] <- 0 #replaces NA with 0.
    meanMatrix <-
      crossprod(logicalMatrixTheories, database) #produces a matrix with the sum of cosines of each paper with each of the theorys
    meanMatrix <-
      meanMatrix / docsByTop #produces a matrix with the mean cosine of each paper with each of the theorys
    cosineAvgsList[[n]] <-
      meanMatrix #stores the matrix of means in a list with n as index for dimensions used.
    meanMatrix <-
      meanMatrix %*% logicalMatrixTheories #produces a vector with the sum of mean cosines for each theory against each theory in dimension n
    meanMatrix <-
      t(meanMatrix) / docsByTop #produces a vector of the means of sums of mean cosines for each theory against each theory in dimension n.
    for (theory in theoryList) {
      #loops through theorys to populate results of all cosines.
      listOfResults[[theory]][, n] <- meanMatrix[theory,]
    }
    n = n + 1
  }
  returnList = list(listOfResults, cosineAvgsList) #makes list of lists with results and mean cosines
  return(returnList) #returns everything
} #function for calculating cosines of the whole matrix. "matrix" is the result of an SVD; "cat" is the catalog to obtain theory information (in this case, catalog$id). Returns a list of two lists: [[1]] is all cosines by paper, [[2]] is a list of matrices of mean distance of each paper with each of the other theorys. [[1]][n] and [[2]][n] are the different dimensions resulting from SVD.

plotTheoryDistances <- function(theory, resultsList) {
  temporaryMatrix <-
    as.data.frame(melt(
      resultsList[[1]][[theory]],
      varnames = c("theory", "dimension"),
      value.name = "cosine"
    )) #long form for ggplot
  plot <-
    ggplot(data = temporaryMatrix,
           aes(
             x = dimension,
             y = cosine,
             color = theory,
             group = theory
           )) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) + labs(
             title = paste(
               "Mean cosine of",
               theory,
               "papers with other theories and itself across dimensions"
             )
           ) + geom_line() + scale_colour_solarized("red") + geom_point(size = 0.7, shape = 3) + guides(colour = guide_legend(override.aes = list(size =
                                                                                                                                                    3)))
  print(plot)
} #function for plotting the mean distance of every theory with all other theorys. "theory" is one of the theorys of theoryList; "resultsList" is the object that compareTheories() returns.


# Loading and cleaning data -----------------------------------------------

## Original Data

# Load document by term matrix and catalog for original data. Assigns ID.

frequencyMatrix <-
  as.matrix(read.table('document_by_term.txt', sep = '\t', header = T))[, -1] # Minus one column, the identifier in the text file.
row.names(frequencyMatrix) <-
  c(1:nrow(frequencyMatrix)) # DocID assignment.

# Remove columns with words that appear in fewer than 5 documents.
frequencyMatrix <-
  frequencyMatrix[, apply(frequencyMatrix, 2, function(x) {
    sum(x == 0) < (dim(frequencyMatrix)[1] - 5)
  })]

# Eliminate hand-coded list of words encoded in allWords.csv
cleanWords <-
  read.csv(file = 'allWords.csv', header = T, sep = ',')
frequencyMatrix <- frequencyMatrix[, cleanWords$Include.]

# Eliminates documents with 0 terms after cleanup of terminology
frequencyMatrix <-
  frequencyMatrix[which(rowSums(frequencyMatrix) > 0), ]

# Load catalog, assign column names to type of data, extract the different theories being tested.
catalog <-
  read.table(
    'catalog.txt',
    stringsAsFactors = F,
    sep = '\t',
    fill = T,
    quote = ""
  )
catalog <-
  catalog[row.names(frequencyMatrix), ] #limit catalog to documents in frequencyMatrix
colnames(catalog) = c('id',
                      'theory',
                      'year',
                      'authors',
                      'title',
                      'journal',
                      'abstract')
theoryList <- unique(catalog$theory) #list of theories for analysis.

## Replication Data

# Same procedure, but for the replication data.
# Loads frequency matrix, modifies terminology and loads catalog

repFrequencyMatrix <-
  as.matrix(read.table(
    'rep_document_by_term.txt',
    sep = '\t',
    header = T,
    quote = ""
  ))[, -1]
row.names(repFrequencyMatrix) <- c(1:nrow(repFrequencyMatrix))

repFrequencyMatrix <-
  repFrequencyMatrix[, apply(repFrequencyMatrix, 2, function(x) {
    sum(x == 0) < (dim(repFrequencyMatrix)[1] - 5)
  })]
repFrequencyMatrix <-
  repFrequencyMatrix[which(rowSums(repFrequencyMatrix) > 0), ]

repCatalog <-
  read.table(
    'rep_catalog.txt',
    stringsAsFactors = F,
    sep = '\t',
    fill = T,
    quote = ""
  )
repCatalog <- repCatalog[row.names(repFrequencyMatrix), ]
colnames(repCatalog) = c('id',
                         'theory',
                         'year',
                         'authors',
                         'title',
                         'journal',
                         'abstract')

## Null Hypothesis

# If nullHypothesis = T, catalog$topic is randomized.

nullHypothesis <- T
if (nullHypothesis == T) {
  catalog$theory <- sample(catalog$theory, length(catalog$theory))
  repCatalog$theory <-
    sample(repCatalog$theory, length(repCatalog$theory))
}

# Examine quality of data for cosine similarity through average lengths of
# abstracts for each theory

# Original Data
lengthsOfTheories <-
  data.frame(
    AvgLength = double(length = 8),
    StandardDevLength = double(length = 8),
    theory = theoryList,
    row.names = theoryList
  )

for (theory in theoryList) {
  matrixOfTheory <- frequencyMatrix[which(catalog$theory == theory), ]
  lengthsOfTheory <- rowSums(matrixOfTheory)
  lengthsOfTheories[theory, "AvgLength"] <- mean(lengthsOfTheory)
  lengthsOfTheories[theory, "StandardDevLength"] <-
    sd(lengthsOfTheory)
}

# Examine significance of relation between theory and abstract length through
# AOV and LM
allLengthsMatrix <-
  data.frame(theory = factor(catalog$theory),
             length = rowSums(frequencyMatrix))
AOVOfLength <- aov(length ~ theory, data = allLengthsMatrix)
summary(AOVOfLength)
LMOfLength <- lm(length ~ theory, data = allLengthsMatrix)
summary(LMOfLength)

# Replication Data
lengthsOfTheoriesRep <-
  data.frame(
    AvgLength = double(length = 8),
    StandardDevLength = double(length = 8),
    theory = theoryList,
    row.names = theoryList
  )

for (theory in theoryList) {
  matrixOfTheory <-
    repFrequencyMatrix[which(repCatalog$theory == theory), ]
  lengthsOfTheory <- rowSums(matrixOfTheory)
  lengthsOfTheoriesRep[theory, "AvgLength"] <- mean(lengthsOfTheory)
  lengthsOfTheoriesRep[theory, "StandardDevLength"] <-
    sd(lengthsOfTheory)
}

allLengthsMatrixRep <-
  data.frame(theory = factor(repCatalog$theory),
             length = rowSums(repFrequencyMatrix))
AOVOfLengthRep <- aov(length ~ theory, data = allLengthsMatrixRep)
summary(AOVOfLengthRep)
LMOfLengthRep <- lm(length ~ theory, data = allLengthsMatrixRep)
summary(LMOfLengthRep)

# Data weighting and dimensionality reduction -----------------------------

## Weighting

# Use entropy function calculateEntropy() to weight matrices using log-entropy.

# Original Data
weightedFrequencyMatrix <-
  calculateEntropy(frequencyMatrix) #entropy
weightedFrequencyMatrix[is.na(weightedFrequencyMatrix)] <-
  0 #replace NA with 0.

# Replication
repWeightedFrequencyMatrix <- calculateEntropy(repFrequencyMatrix)
repWeightedFrequencyMatrix[is.na(repWeightedFrequencyMatrix)] <- 0

## Singular Value Decomposition

# Original Data
reducedMatrixSVD = svd(weightedFrequencyMatrix)

# Replication Data
repReducedMatrixSVD = svd(repWeightedFrequencyMatrix)

## Generate document and term loadings using SVD matrices

# Multiply matrices and then assign row/column names to coincide with original matrices.

# Original Data
documentLoadings <- reducedMatrixSVD$u %*% diag(reducedMatrixSVD$d)
termLoadings <- reducedMatrixSVD$v %*% diag(reducedMatrixSVD$d)
row.names(documentLoadings) <- row.names(weightedFrequencyMatrix)
row.names(termLoadings) <- colnames(weightedFrequencyMatrix)

# Replication Data
repDocumentLoadings <-
  repReducedMatrixSVD$u %*% diag(repReducedMatrixSVD$d)
repTermLoadings <-
  repReducedMatrixSVD$v %*% diag(repReducedMatrixSVD$d)
row.names(repDocumentLoadings) <-
  row.names(repWeightedFrequencyMatrix)
row.names(repTermLoadings) <- colnames(repWeightedFrequencyMatrix)

# Reduce the dimensionality of the document loading matrices
termLoadings <- termLoadings[, 1:200]
documentLoadings <- documentLoadings[, 1:200]
repTermLoadings <- repTermLoadings[, 1:200]
repDocumentLoadings <- repDocumentLoadings[, 1:200]

# Predict theory using Generalized Linear Models --------------------------

## Parameters of the prediction

# Min and max number of dimensions are used to control the number of
# dimensions that are to be used in the construction of the models.
# The procedure loops through the dimensions resulting from the SVD.
# Starts at minNumberOfDimensions, stops at maxNumberOfDimensions.
minNumberOfDimensions <- 2 # does not work if lower than 2.
maxNumberOfDimensions <- 50

# Controls the maximum number of dimensions to be considered as
# "best predictors" for each theory.
dimensionsForBestPredictors <- 80

# "Repeats" is the number of iterations of the predicting process to be aggregated.
iterations <- 10000

# "Source" controls which data set is to be used: "original" uses the original dataset,
# "replication" uses the replication data, and "cross" uses the projection of the
# replication data into the semantic space of the original dataset to predict their
# theories with the models built with the original dataset.\
source <- "cross"

## Load objects based on parameters

if (source == "original") {
  theCatalog <- catalog
  theSVD <- documentLoadings
}
if (source == "replication") {
  theCatalog <- repCatalog
  theSVD <- repDocumentLoadings
}
if (source == "cross") {
  theCatalog <- catalog
  crossCatalog <- repCatalog
  theSVD <- documentLoadings
  
  # Project the replication data on the space generated by the original data.
  replicationProjection <-
    matrix(0,
           nrow = nrow(repFrequencyMatrix),
           ncol = ncol(frequencyMatrix))
  
  # Set names of columns and rows to access through previous matrices
  row.names(replicationProjection) <- row.names(repFrequencyMatrix)
  colnames(replicationProjection) <- colnames(frequencyMatrix)
  
  # Determine the words that appear in both original data and replication data.
  sharedWords <-
    colnames(repFrequencyMatrix)[which(colnames(repFrequencyMatrix) %in% colnames(frequencyMatrix))]
  
  # Build a new DbT matrix that has the documents in replication as rows and the
  # shared words between both datasets as columns.
  replicationProjection[, sharedWords] <-
    repFrequencyMatrix[, sharedWords]
  
  # Projects the replication data into the SV space of the original dataset.
  crossSVD <-
    replicationProjection %*% t(ginv(termLoadings)) # Moore-Penrose Generalized Inverse
}

## Best predicting dimensions for each theory

# For each theory, fill the "best predictors" matrix with the dimensions that most
# differentiate that theory from the other 7 theories using GLMS.

bestPredictors <- c()
predictorRatings <-
  c() # Store the ratings of each dimension to use in evaluating dimension meaning afterwards.
for (theory in theoryList) {
  glmOutput = glm(theCatalog$theory == theory ~ .,
                  data = data.frame(theSVD[, 1:dimensionsForBestPredictors]),
                  family = binomial)
  # Sort the dimensions according to the absolute value of their rating for each theory and store.
  # Absolute value allows to use both positive and negative prediction for each theory.
  bestPredictors = rbind(bestPredictors, data.frame(t(sort(
    abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),
    # 2 is here to skip intercept.
    ind = T,
    decreasing = T
  )$ix[1:maxNumberOfDimensions])))
  
  # Also store the sorted ratings
  predictorRatings <-
    rbind(predictorRatings, glmOutput$coefficients[2:length(glmOutput$coefficients)][sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),
                                                                                          ind = T,
                                                                                          decreasing = T)$ix])
}
row.names(bestPredictors) <- theoryList
row.names(predictorRatings) <- theoryList

## Prediction

# The prediction is parallelized using the "foreach" and "doparallel" packages.
# In each iteration, a random training set of 70% of the papers of each theory is selected.
# The remaining papers are presented to each GLM model and the probability returned
# by the model that the paper belongs to that theory is collected.
# The highest prediction value is selected as the "predicted" theory and stored.
# The number of iterations to be aggregated for each dimension is controlled
# by parameter iterations. Iterations are averaged.

# Each parallel process uses a different dimension from the set of dimensions between
# minNumberOfDimensions and maxNumberOfDimensions.
dimensionsForTesting <-
  c(minNumberOfDimensions:maxNumberOfDimensions)

# Parallel process monitored by a logger object (log4r).
logger = create.logger()
logfile(logger) = 'monitor.log'
level(logger) = 'INFO'

# clusters controls the number of parallel processes; change to fit the number of cores in CPU.
clusters <-
  makeCluster(detectCores() - 1) # default: number of cores in CPU - 1.
registerDoParallel(clusters)

# The results of this prediction are collected for each dimension in listOfResults.
listOfResults <-
  foreach(dimension = dimensionsForTesting,
          .packages = "log4r") %dopar% {
            # Pre-allocate the result list for this dimension.
            resultsListModel <-
              lapply(c(1:iterations), matrix, nrow = 8, ncol = 8)
            
            # Repeats the process of training-prediction a number of times specified by parameter "iterations".
            iteration <- 1
            while (iteration <= iterations) {
              trainingSet <- c()
              for (theory in theoryList) {
                trainingSet <-
                  c(trainingSet, sample(which(theCatalog$theory == theory), round(length(
                    which(theCatalog$theory == theory)
                  ) * 0.7)))
              }
              
              # If the procedure is either predicting original for predicting original, or replication for predicting replication, set the trainingset as the rest of the papers.
              if (source == "original" |
                  source == "replication") {
                testSet = setdiff(1:nrow(theSVD), trainingSet)
              }
              
              # If using original to predict replication, trainingset is original dataset, and testset is replication set.
              if (source == "cross") {
                trainingSet <- c(1:(nrow(theSVD)))
                testSet <- c(1:(nrow(crossSVD)))
              }
              predictionResults = c()
              
              # Loop through the models of each theory
              for (theory in theoryList) {
                # Take only the dimensions in bestPredictors from the document vector.
                trainingdata <-
                  data.frame(theSVD[trainingSet, unlist(bestPredictors[theory, 1:dimension])]) # prepare training data of model by using the dimensions selected as the best predictors for each theory and the documents selected to be training.
                
                # Build the GLM of the theory trained on trainingSet
                glmOfTheory = glm(theCatalog$theory[trainingSet] == theory ~ .,
                                  data = trainingdata,
                                  family = binomial)
                
                # If predicting inside the dataset
                if (source == "original" |
                    source == "replication") {
                  testdata = data.frame(theSVD[testSet, unlist(bestPredictors[theory, 1:dimension])])
                  
                  # Store the probability that the paper belongs to the theory being tested
                  predicted = predict.glm(glmOfTheory, newdata = testdata, type = "response")
                  predictionResults = cbind(predictionResults, scale(predicted)) # Add to a matrix and scale
                }
                
                # If cross-predicting, original data is training and replication is test
                if (source == "cross") {
                  predicted = predict.glm(glmOfTheory,
                                          newdata = data.frame(crossSVD[testSet, unlist(bestPredictors[theory, 1:dimension])]),
                                          type = "response")
                  predictionResults = cbind(predictionResults, scale(predicted))
                }
              }
              
              # Aggregate predictions for each theory in a dataframe.
              # Rows are documents, columns are the probability that doc belongs to that theory.
              predictionResults = data.frame(predictionResults)
              colnames(predictionResults) = theoryList
              
              # Evalute if the predicted theory is correct.
              if (source == "original" | source == "replication") {
                predictionResults$theory = theCatalog$theory[testSet]  # Add a column with the correct theory
                predictionResults$predicted_theory = theoryList[max.col(predictionResults[, 1:8])] # Add a column with the predicted theory
                
                # Determine how many times each theory was predicted as each other theory in percentages.
                resultTable <-
                  (
                    table(
                      predictionResults$theory,
                      predictionResults$predicted_theory
                    ) / as.vector(table(theCatalog$theory[testSet]))
                  ) * 100
                
              }
              
              # Same procedure, modified for cross prediction.
              if (source == "cross") {
                predictionResults$theory = crossCatalog$theory[testSet]
                predictionResults$predicted_theory = theoryList[max.col(predictionResults[, 1:8])]
                resultTable <-
                  (
                    table(
                      predictionResults$theory,
                      predictionResults$predicted_theory
                    ) / as.vector(table(crossCatalog$theory[testSet]))
                  ) * 100
              }
              
              # Store this iteration for final aggregation.
              resultsListModel[[iteration]] <-
                resultTable
              iteration = iteration + 1
              
              # Periodically output to log file to monitor process.
              if (iteration %% (round(iterations * 0.2)) == 0) {
                info(logger,
                     paste(
                       "dimension ",
                       dimension,
                       ", ",
                       "iteration number",
                       iteration
                     ))
              }
            }
            
            # Average the results of the tables generated by each iteration.
            finalPredictionTable <-
              matrix(0, nrow = 8, ncol = 8)
            for (matrix in resultsListModel) {
              finalPredictionTable <- finalPredictionTable + matrix
            }
            finalPredictionTable <-
              finalPredictionTable / length(resultsListModel)
            
            # Return averaged table. This table is then returned in a list by foreach.
            return(finalPredictionTable)
          }

stopCluster(clusters) # Stop the parallel clusters

# Name the objects in the list with the dimensions used
names(listOfResults) <-
  dimensionsForTesting

# Allocate results to evaluate the performance of each dimension
performanceMatrix <-
  matrix (0, nrow = length(dimensionsForTesting), ncol = 9)
row.names(performanceMatrix) <-
  as.character(dimensionsForTesting)
colnames(performanceMatrix) <-
  c(theoryList, "mean") # Add a column for average performance
for (dimension in dimensionsForTesting) {
  performanceMatrix[as.character(dimension),] <-
    c(diag(listOfResults[[as.character(dimension)]]), mean(diag(listOfResults[[as.character(dimension)]])))
} # Fill list with the percentage of correct predictions for each theory

## Plotting prediction results

# Confusion matrix for each dimension
dimension = 5 # Specify the dimension to use for producing the confusion matrix
predictionMatrix <-
  listOfResults[[as.character(dimension)]]
meltedResults <-
  melt(
    predictionMatrix,
    varnames = c("theory1", "theory2"),
    value.name = "Percentage.Predicted"
  )

# Confusion matrix for chosen dimension as a heatmap
heatmap <-
  ggplot(meltedResults, aes(y = theory1, x = ordered(theory2, levels = rev(sort(
    unique(theory2)
  ))))) + geom_tile(aes(fill = Percentage.Predicted)) + coord_equal() + scale_fill_gradient(
    limits = c(0, 100),
    low = "white",
    high = "seagreen",
    guide =  guide_colorbar(title = paste("Performance (%)", "\n"))
  ) + xlab("") + ylab("") + theme(
    axis.text = element_text(size = 14),
    axis.text.x = element_text(
      angle = 330,
      hjust = 0.4,
      vjust = 0.7,
      size = 14
    )
  ) + geom_text(aes(label = paste(round(
    Percentage.Predicted, 1
  ), "%", sep = "")), colour = "gray25", size = 5)
print(heatmap)

# Barplot of performance for each theory at the chosen dimension.
# Includes solid line for mean performance, dotdashed line for chance.
dimensionDataFrame <- performanceMatrix[as.character(dimension),]
meanPerfromanceDimension <- dimensionDataFrame["mean"]
dimensionDataFrame <-
  dimensionDataFrame[which(names(dimensionDataFrame) != "mean")]
dimensionDataFrame <-
  cbind(melt(dimensionDataFrame), names(dimensionDataFrame))
colnames(dimensionDataFrame) <- c("Performance", "Theory")
barplot <-
  ggplot(data = dimensionDataFrame, aes(
    x = Theory,
    y = Performance,
    fill = Theory,
    guide = F
  )) +
  geom_bar(stat = "identity") +
  scale_colour_brewer(type = "qual",
                      palette = "Paired",
                      guide = F) +
  geom_hline(yintercept = meanPerfromanceDimension,
             linetype = "solid",
             size = 1.2) +
  geom_hline(yintercept = 12.5,
             linetype = "dotdash",
             size = 1.2) +
  ylim(c(0, 100)) + guides(fill = F) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
print(barplot)

# Performance of each theory for each dimension used, including mean. As heatmap.
meltedResults <-
  melt(performanceMatrix[, 1:9],
       varnames = c("D", "theory"),
       value.name = "Performance")
heatmap <-
  ggplot(meltedResults, aes(y = theory, x = ordered(D))) +
  geom_tile(aes(fill = Performance), color = "white") + coord_equal() +
  scale_fill_gradient(limits = c(0, 100),
                      low = "white",
                      high = "seagreen") +
  xlab("") + ylab("") + theme(axis.text = element_text(size = 12)) +
  geom_text(aes(label = paste(round(Performance, 0))), size = 4, colour = "gray25")
print(heatmap)


# Panel of line plots for each theory across all dimensions used
meanResultsTable <-
  as.data.frame(meltedResults[which(meltedResults$theory != "mean"),]) # Don't show mean on this plot
plot <-
  ggplot(data = meanResultsTable, aes(
    x = D,
    y = Performance,
    color = theory,
    group = theory
  )) +
  ylim (0, 100) + theme_gray() + geom_line(size = 1) + ylab("Performance") +  xlab("Dimensions") +
  scale_colour_brewer(type = "qual",
                      palette = "Paired",
                      guide = F) +
  geom_point(size = 1, shape = 3) + theme(text = element_text(size = 12)) + facet_wrap(~ theory, ncol = 2, nrow = 5)
print(plot)

# Line plot of mean performance

meanPerformance <-
  data.frame("Mean.Performance" = rowMeans(performanceMatrix),
             "D" = dimensionsForTesting)
plot <-
  ggplot(data = meanPerformance, aes(y = Mean.Performance, x = D)) +
  theme_gray() + geom_line(size = 1.5, color = "seagreen") +
  geom_point(size = 1.5,
             shape = 3,
             color = "seagreen") + scale_y_continuous(limits = c(0, 100)) +
  labs(y = "Mean Performance (%)") + theme(text = element_text(size = 15))
print(plot)


# Cosine similarity among theories ----------------------------------------

## This section generates a comparison of the similarity between and within
## theories by using the cosine generated by pairs of abstracts in as a vector
## space representation in the space of the SVD. Similarity data is then used to
## measure self similarity models and clustering.

## Generate cosine matrices for each dimension

originalCosines <-
  compareTheories(documentLoadings, catalog) # Has two elements: [[1]] list of a 8 matrices, one for each theory, with mean distance of that theory with each other theory. (2) list of one matrix for each dimension with theories as rows and all papers as column. cells show the mean distance of that theory with that paper.
originalAverages <-
  originalCosines[[2]] # Extracts the average distance of theory by paper.
replicationCosines <-
  compareTheories(repDocumentLoadings, repCatalog)
replicationAverages <- replicationCosines[[2]]

# Plots of mean distance of theories with other theories for dimension
plotTheoryDistances("bayesian", originalCosines)
plotTheoryDistances("symbolic", originalCosines)
plotTheoryDistances("connectionism", originalCosines)
plotTheoryDistances("embodied", originalCosines)
plotTheoryDistances("distributed", originalCosines)
plotTheoryDistances("enactive", originalCosines)
plotTheoryDistances("dynamical", originalCosines)
plotTheoryDistances("ecological", originalCosines)

## Linear model predicting self-similarity

# Parameters

# "source defines the dataset to be used: "original" or "replication".
# "dimension" controls the number of dimensions used in the analysis (D).
source <- "original"
dimension <- 10

if (source == "original") {
  listOfAverages <- originalAverages
  paperCatalog <- catalog
}
if (source == "replication") {
  listOfAverages <- replicationAverages
  paperCatalog <- repCatalog
}

# Take the average distance of each theory with each paper of that same theory
# and then model that relation with linear model

allData = c()

for (theory in theoryList) {
  # Populate the mean cosine of theory with papers of that theory data
  data = listOfAverages[[dimension]][which(theory == theoryList), paperCatalog$theory == theory]
  allData = rbind(allData, data.frame(theory = theory, cosine = data))
}

# Build a linear model of mean distance of theory with papers of that theory and the theory#
selfSimilarityModel = lm(cosine ~ theory, data = allData)
summary(selfSimilarityModel)

#visualize with boxplot#

boxplot <-
  ggplot(allData, aes(x = theory, y = cosine)) + geom_boxplot(fill = "#fdf6e3",
                                                              colour = "#2aa198",
                                                              outlier.color = "#2aa198") + scale_x_discrete(name = "Theory") + scale_y_continuous(
                                                                name = "Mean self-cosine",
                                                                breaks = seq(-0.1, 1, .10),
                                                                limits = c(-0.1, 1)
                                                              ) + theme_solarized(base_size = 14) + theme(axis.text = element_text(colour = "#586e75")) + labs(
                                                                title = paste(
                                                                  "Self-similarity of different theories of dimension number",
                                                                  dimension
                                                                ),
                                                                subtitle = "Cosines of members of a theory with other members of that theory as predicted by theory membership"
                                                              )

## Cluster analysis

# This section uses the similarity data of cosines gathered previously

# Parameters
# value of D ("dimension") and source of the data ("dendrogramModesource")
dimension <- 10 #change dimension being considered
source <- "original" #original or replication

if (source == "original") {
  matrixOfAverages <- originalCosines[[2]][[dimension]]
  paperCatalog <- catalog
  documents <- documentLoadings
}

if (source == "replication") {
  matrixOfAverages <- replicationCosines[[2]][[dimension]]
  paperCatalog <- repCatalog
  documents <- repDocumentLoadings
}

# Generates a similarity matrix of cosine data to visualize in a heatmap
theoryListReordered = c(
  "bayesian",
  "connectionism",
  "symbolic",
  "distributed",
  "dynamical",
  "enactive",
  "ecological",
  "embodied"
) # Reorder theories to better visualize if intuitive cluster shows

similarityMatrix = matrix(0, nrow = 8, ncol = 8)
row.names(similarityMatrix) = theoryListReordered
colnames(similarityMatrix) = theoryListReordered

for (theory in theoryListReordered) {
  #loops through first theory (row)
  for (i in 1:8) {
    #loops through second theory (column)
    similarityMatrix[theory, theoryListReordered[i]] = mean(matrixOfAverages[theory, which(paperCatalog$theory ==
                                                                                             theoryListReordered[i])])
  }
}

# Visualize similarity matrix with a heatmap
meltedSimilarityMatrix <-
  melt(
    similarityMatrix,
    varnames = c("theory1", "theory2"),
    value.name = "Similarity"
  ) #long form of simMatrix for ggplot
heatmap <-
  ggplot(meltedSimilarityMatrix, aes(x = theory1, y = ordered(theory2, levels = rev(sort(
    unique(theory2)
  ))))) + geom_tile(aes(fill = Similarity), colour = "white") + coord_equal() +  scale_fill_gradient(
    low = "white",
    high = "seagreen",
    limits = c(0, 1),
    guide =  guide_colorbar(title = paste("Cosine", "\n"))
  ) + xlab("") + ylab("") + theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 330,
      hjust = 0.4,
      vjust = 0.7
    )
  ) + geom_text(aes(label = paste(round(Similarity, 2))), size = 5, colour = "gray25")
print(heatmap)

# Hierarchical cluster analysis of the distance matrix
theory_hclust <-
  hclust(dist(similarityMatrix, upper = T), method = "average") #applies hclust algorithm to the similarity matrix

# Visualize theory_hclust with a dendrogram and mark clusters using
# cutreedynamic package The minimum cluster size is set to 1, the parameter
# controlling the stringency of the clustering is set to the default (0). We
# used the "hybrid" method which takes both the dendrogram and a distance matrix
# to generate clusters. PAM was not used. We provided the colors (colorSeq) of
# the fillbox marking the clusters for up to 3 clusters.
par(lwd = 2,
    cex.axis = 1.2,
    las = 1)
print(
  hclustplot(
    theory_hclust,
    colors = labels2colors(
      cutreeDynamic(
        theory_hclust,
        minClusterSize = 1,
        method = "hybrid",
        deepSplit = 0,
        distM = as.matrix(dist(similarityMatrix)),
        pamStage = T
      ),
      colorSeq = c("#2E8B57", "#FF2052", "#8b2e62")
    ),
    fillbox = T,
    las = 0,
    cex = 1.2,
    mar = c(3, 3, 2, 0.5),
    font = 2
  )
)


# Exploration of dimension meanings using a varimax rotated matrix --------

# Varimax rotation of the term matrix and subsequent rotation of document
# loadings
termVarimax <- varimax(termLoadings)
termLoadingsVarimax <- unclass(termVarimax$loadings)
documentLoadingsVarimax <- documentLoadings %*% termVarimax$rotmat

# Use same procedure as before to rank the dimensions according to how well
# they predict theory.
bestPredictorsVarimax <- c()
predictorRatingsVarimax <-
  c() # Store ratings to compare positive versus negative in term loading inspection
for (theory in theoryList) {
  glmOutput = glm(
    paperCatalog$theory == theory ~ .,
    data = data.frame(documentLoadingsVarimax[, 1:dim(documentLoadingsVarimax)[2]]),
    family = binomial
  )
  bestPredictorsVarimax = rbind(bestPredictorsVarimax, data.frame(t(sort(
    abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),
    ind = T,
    decreasing = T
  )$ix[1:dim(documentLoadingsVarimax)[2]])))
  predictorRatingsVarimax <-
    rbind(predictorRatingsVarimax, glmOutput$coefficients[2:length(glmOutput$coefficients)][sort(abs(glmOutput$coefficients[2:length(glmOutput$coefficients)]),
                                                                                                 ind = T,
                                                                                                 decreasing = T)$ix])
}
row.names(bestPredictorsVarimax) <- theoryList
row.names(predictorRatingsVarimax) <- theoryList

# For each theory, build a data frame storing the dimensions that best predict
# them, their coefficient, the sign of the coefficient, and the top 100 terms.
listOfDimensionTerms <- list()
for (n in 1:length(theoryList)) {
  theory <- theoryList[n]
  predictors <- unlist(bestPredictorsVarimax[theory,])
  ratings <- unlist(predictorRatingsVarimax[theory,])
  ratingsValence <-
    ratings > 0 # Boolean for sign of the predictor. true if > 0, false if < 0
  words <- c()
  i = 1
  for (i in c(1:length(predictors))) {
    dimension <- predictors[i]
    words <-
      c(words, toString(head(names(
        sort(termLoadingsVarimax[, dimension], decreasing = T)
      ), 100))) # Extract 100 top words
  }
  nameOfDF <-
    paste(theory, "Meanings", sep = "") # Generate the name of the object
  assign(nameOfDF,
         data.frame(predictors, ratings, ratingsValence, words))
  listOfDimensionTerms[[n]] <- get(x = nameOfDF)
}

names(listOfDimensionTerms) <- theoryList

# Extract the words for each theory and weight them using coefficient of
# dimension * loading of the term in that dimension.
# The code saves one wordcloud for each theory, either positive or negatives.
# To produce positive wordclouds, comment line 968 and 1005 through 1010.
# To produce negative wordclouds, comment line 967 and 999 through 1004.
for (theory in theoryList) {
  wordsAndRating <- data.frame()
  for (i in 1:dim(wordList[[theory]][1])) {
    if (wordList[[theory]]$ratingsValence[i]) {
      # if (!wordList[[theory]]$ratingsValence[i]){
      words <-
        unlist(strsplit(as.character(wordList[[theory]]$words[i]), ', '))
      dimRating <-
        rep(abs(wordList[[theory]]$ratings[i]), each = length(words))
      weights <- c()
      for (s in 1:length(words)) {
        word <- words[s]
        dimensionRating <- dimRating[s]
        dimension <- wordList[[theory]]$predictors[i]
        wordWeight <-
          dimensionRating * termLoadingsVarimax[word, dimension]
        if (wordWeight < 0) {
          wordWeight <- 0
        }
        weights <- c(weights, wordWeight)
      }
      wordsAndRating <- rbind(wordsAndRating, cbind(words, weights))
    }
  }
  colnames(wordsAndRating) <- c("words", "Weight")
  wordsAndRating$words <- as.character(wordsAndRating$words)
  wordsAndRating$Weight <-
    as.numeric(levels(wordsAndRating$Weight))[wordsAndRating$Weight]
  
  wordsAndRatingAggregated <-
    aggregate(data = wordsAndRating, Weight ~ words, FUN = sum)
  wordsAndRatingForWordCloud <-
    wordsAndRatingAggregated[which(wordsAndRatingAggregated$Weight > median(wordsAndRatingAggregated$Weight)), ]
  wordsAndRatingForWordCloud$Weight <-
    round(wordsAndRatingForWordCloud$Weight)
  pdf(
    file = paste(theory, "Positive.pdf", sep = ""),
    width = 20,
    height = 20,
    onefile = F
  )
  #pdf(
  #  file = paste(theory, "Negative.pdf", sep = ""),
  #  width = 20,
  #  height = 20,
  #  onefile = F
  #)
  wordcloud(
    wordsAndRatingForWordCloud$words,
    wordsAndRatingForWordCloud$Weight,
    scale = c(10, 3),
    random.order = F,
    rot.per = 0,
    use.r.layout = F,
    max.words = 50
  )
  dev.off()
}