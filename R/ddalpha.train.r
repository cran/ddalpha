################################################################################
# File:             ddalpha.train.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     28.02.2013
# 
# Contains the training function of the DDalpha-classifier.
# 
# For a description of the algorithm, see:
#   Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric 
#     classification based on data depth. Statistical Papers.
#   Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world 
#     data with the DDalpha-procedure. Mimeo.
################################################################################

ddalpha.train <- function(data, 
                          depth = "randomTukey", 
                          separator = "alpha", 
                          outsider.methods = "LDA", 
                          outsider.settings = NULL, 
                          aggregation.method = "majority", 
                          knnrange = NULL, 
                          num.chunks = 10, 
                          num.directions = 1000, 
                          use.convex = FALSE, 
                          max.degree = 3, 
                          mah.estimate = "moment", 
                          mah.parMcd = 0.75, 
                          mah.priors = NULL){
  # Check for data consistency #1
  if (!is.matrix(data) 
      || !is.numeric(data)){
    warning("Argument data has unacceptable format. Classifier can not be trained!!!")
    return (NULL)
  }
  
  # Raw data processing
  ddalpha <- .ddalpha.create.structure(data)
  
  # Check for data consistency #2
  if (ddalpha$numPatterns < 2){
    warning("Number of patterns is < 2. Classifier can not be trained!!!")
    return (NULL)
  }
  # TODO ddalpha$numPatterns should be restricted from above as well
  if (ddalpha$dimension < 2){
    warning("Data dimension is < 2. Classifier can not be trained!!!")
    return (NULL)
  }
  for (i in 1:ddalpha$numPatterns){
    if (ddalpha$patterns[[i]]$cardinality < ddalpha$dimension + 1){
      warning("At least in one patern number of the points < (dimension + 1). Classifier can not be trained!!!")
      return (NULL)
    }
  }
  
  # Reassigning the properties
  if (!is.character(depth)
      || length(depth) != 1
      || !(depth %in% c("zonoid", "randomTukey", "Mahalanobis", "projectionRandom", "projectionLinearize", "spatial"))){
    ddalpha$methodDepth <- "randomTukey"
    warning("Argument \"depth\" not specified correctly. \"randomTukey\" is used as a default value")
  }else{
    ddalpha$methodDepth <- depth
  }
  if (!is.character(separator)
      || length(separator) != 1
      || !(separator %in% c("alpha", "knnlm"))){
    ddalpha$methodSeparator <- "alpha"
    warning("Argument \"separator\" not specified correctly. \"alpha\" is used as a default value")
  }else{
    ddalpha$methodSeparator <- separator
  }
  if (!is.character(aggregation.method)
      || length(aggregation.method) != 1
      || !(aggregation.method %in% c("majority", "sequent"))){
    ddalpha$methodAggregation <- "majority"
    warning("Argument \"aggregation.method\" not specified correctly. \"majority\" is used as a default value")
  }else{
    ddalpha$methodAggregation <- aggregation.method
  }
  if (ddalpha$methodAggregation == "majority"){
    maxChunks <- ddalpha$patterns[[ddalpha$numPatterns]]$cardinality + ddalpha$patterns[[ddalpha$numPatterns - 1]]$cardinality
  }else{
    maxChunks <- ddalpha$numPoints
  }
  if (!is.numeric(num.chunks) 
      || is.na(num.chunks) 
      || length(num.chunks) != 1 
      || !.is.wholenumber(num.chunks) 
      || !(num.chunks > 0 && num.chunks <= maxChunks)){
    ddalpha$numChunks <- maxChunks
    warning("Argument \"num.chunks\" not specified correctly. ", maxChunks, " is used instead")
  }else{
    ddalpha$numChunks <- num.chunks
  }
  if (is.null(knnrange)
      || !is.numeric(knnrange) 
      || is.na(knnrange) 
      || length(knnrange) != 1 
      || !.is.wholenumber(knnrange) 
      || !(knnrange >=2 && knnrange <= ceiling(ddalpha$numPoints/2))){
    
    if (is.null(knnrange))      knnrange <- 10*( (ddalpha$numPoints)^(1/ddalpha$numPatterns) ) + 1
    
    knnrange <- min(knnrange, ceiling(ddalpha$numPoints/2))
    knnrange <- max(knnrange, 2)
    
    ddalpha$knnrange <- knnrange
    if (!is.null(knnrange)) warning("Argument \"knnrange\" not specified correctly. ", knnrange, " is used instead")
  }else{
    ddalpha$knnrange <- knnrange
  }
  if (ddalpha$methodDepth == "randomTukey" 
      && (!is.numeric(num.directions) 
          || is.na(num.directions) 
          || length(num.directions) != 1 
          || !.is.wholenumber(num.directions) 
          || !(num.directions > 1 && num.directions < 10000000))
      ){
    ddalpha$numDirections <- 1000
    warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
  }else{
    ddalpha$numDirections <- num.directions
  }
  if (!is.logical(use.convex) 
      || length(use.convex) != 1 
      || !(use.convex %in% c(TRUE, FALSE))){
    ddalpha$useConvex <- FALSE
    warning("Argument \"use.convex\" not specified correctly. FALSE is used as a default value")
  }else{
    ddalpha$useConvex <- use.convex
  }
  if (!is.numeric(max.degree) 
      || is.na(max.degree) 
      || length(max.degree) != 1 
      || !.is.wholenumber(max.degree) 
      || !(max.degree %in% 1:10)){
    ddalpha$maxDegree <- 3
    warning("Argument \"max.degree\" not specified correctly. 3 is used as a default value")
  }else{
    ddalpha$maxDegree <- max.degree
  }
  if (!is.character(mah.estimate) 
      || length(mah.estimate) != 1 
      || !(mah.estimate %in% c("moment", "MCD"))){
    warning("In treatment number ", i, ": Argument \"mah.estimate\" not specified correctly. \"moment\" is used as a default value")
    ddalpha$mahEstimate <- "moment"
  }else{
    ddalpha$mahEstimate <- mah.estimate
  }
  if (!is.vector(mah.priors, mode = "double") 
      || is.na(min(mah.priors)) 
      || length(mah.priors) != ddalpha$numPatterns 
      || min(mah.priors) <= 0 
      || max(mah.priors) <= 0){
    if (!is.null(mah.priors)){
      warning("In treatment number ", i, ": Argument \"mah.priors\" not specified correctly. Defaults in the form of class portions are applied")
    }
    ddalpha$mahPriors <- NULL
  }else{
    ddalpha$mahPriors <- mah.priors/sum(mah.priors)
  }
  if (!is.vector(mah.parMcd, mode = "double") 
      || is.na(min(mah.parMcd)) 
      || length(mah.parMcd) != 1 
      || mah.parMcd < 0.5 
      || mah.parMcd > 1){
    if (ddalpha$mahEstimate == "MCD"){
      warning("In treatment number ", i, ": Argument \"mah.parMcd\" not specified correctly. 0.75 is used as a default value")
    }
    ddalpha$mahParMcd <- 0.75
  }else{
    ddalpha$mahParMcd <- mah.parMcd
  }
  
  # Calculate depths
  ddalpha <- .ddalpha.learn.depth(ddalpha)
  
  # Learn classification machine
  if (ddalpha$methodSeparator == "alpha"){
    ddalpha <- .ddalpha.learn.alpha(ddalpha)
  }
  if (ddalpha$methodSeparator == "knnlm"){
    ddalpha <- .ddalpha.learn.knnlm(ddalpha)
  }
  
  # Learn outsider treatments if needed
  if (!(ddalpha$methodDepth %in% 
          c("Mahalanobis", "projectionRandom", "projectionLinearize", "spatial"))){
    ddalpha <- .ddalpha.learn.outsiders(ddalpha = ddalpha, 
                                        methodsOutsider = outsider.methods, 
                                        settingsOutsider = outsider.settings)
  }
  class(ddalpha) <- "ddalpha"
  
  return (ddalpha)
}
