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
                          pretransform = NULL,
                          mah.parMcd = 0.75,
                          use.convex = FALSE,
                          
                          ...
      
#                           # knn
#                           knnrange = NULL, 
#                           # alpha
#                           num.chunks = 10, 
#                           max.degree = 3, 
#                          
#                           # randomTukey depth
#                           num.directions = 1000,                          
#                           # Mahalanobis depth
#                           mah.estimate = "moment", 
#                           mah.parMcd = 0.75, 
#                           mah.priors = NULL
                          
                          ){
  # Check for data consistency #1
  if (!(is.matrix(data) && is.numeric(data)
      || is.data.frame(data) && prod(sapply(data[,-ncol(data)], is.numeric)))){
    stop("Argument data has unacceptable format. Classifier can not be trained!!!")
  }

  # Raw data processing
  ddalpha <- .ddalpha.create.structure(data)
  
  # Check for data consistency #2
  if (ddalpha$numPatterns < 2){
    stop("Number of patterns is < 2. Classifier can not be trained!!!")
  }
  # TODO ddalpha$numPatterns should be restricted from above as well
  if (ddalpha$dimension < 2){
    stop("Data dimension is < 2. Classifier can not be trained!!!")
  }
  for (i in 1:ddalpha$numPatterns){
    if (ddalpha$patterns[[i]]$cardinality < ddalpha$dimension + 1){
      stop("At least in one patern number of the points < (dimension + 1). Classifier can not be trained!!!")
    }
  }
    
  ## Validating the properties
  depthsThatNeedNoScaling = c("zonoid", "randomTukey", "Mahalanobis", "projectionRandom", "projectionLinearize", "spatial")
  supportedDepths = c(depthsThatNeedNoScaling, "potential (not implemented)")
  
  if (!is.character(depth)
      || length(depth) != 1
      || !(depth %in% supportedDepths)){
    ddalpha$methodDepth <- "randomTukey"
    warning("Argument \"depth\" not specified correctly. \"randomTukey\" is used as a default value")
  }else{
    ddalpha$methodDepth <- depth
  }
  if (!is.character(separator)
      || length(separator) != 1
      || !(separator %in% c("alpha", "polynomial", "knnlm", "maxD"))){
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
  
  ddalpha$needtransform = 0
  if (!is.null(pretransform))
    if (ddalpha$methodDepth %in% depthsThatNeedNoScaling){
      warning("The used depth method is affine-invariant and pretransform doesn't influence the result. The data won't be transformed.")
    } else if (pretransform == "1MahMom" || pretransform == "1MahMCD"){
      ddalpha$needtransform = 1
      
      if (pretransform == "1MahMom")
        mm <- mah.moment(data[,1:(ncol(data)-1)])
      else   # "1MahMCD"
        mm <- mah.mcd(data[,1:(ncol(data)-1)], mah.parMcd)              
      
      for (i in 1:ddalpha$numPatterns){
        ddalpha$patterns[[i]]$transformer <- MahMomentTransformer(mm$mu, mm$b)
        ddalpha$patterns[[i]]$points <- ddalpha$patterns[[i]]$transformer(ddalpha$patterns[[i]]$points)
      }
    } else if (pretransform == "NMahMom" || pretransform == "NMahMCD"){
      ddalpha$needtransform = 2
      
      for (i in 1:ddalpha$numPatterns){
        if (pretransform == "NMahMom")
          mm <- mah.moment(ddalpha$patterns[[i]]$points)
        else   # "NMahMCD"
          mm <- mah.mcd(ddalpha$patterns[[i]]$points, mah.parMcd)
        
        ddalpha$patterns[[i]]$transformer <- MahMomentTransformer(mm$mu, mm$b)
      }
    } else {
      warning("Argument pretransform has unacceptable format. The data won't be transformed.")
    }  
  
  # appends ddalpha with the values from the given list (lst)
  ddalpha.append <- function(lst){ if(is.list(lst)) for(k in names(lst))  ddalpha[[k]] <<- lst[[k]] }
  # calls the validation method for the selected separator || depth
  # NO errors or warnings if the function doesn't exist!!
  # ddalpha is READONLY inside the validators
  validate <- function(method_name){
    f <- try(match.fun(paste0(method_name, ".validate")), silent = T)
    if (is.function(f)){
      lst  <- f(ddalpha, ...)
      ddalpha.append(lst) 
    }
  }
   
  ## Separator parameters validation
 
  validate(ddalpha$methodSeparator)
  
  ## Depth parameters validation
 
  if (!is.logical(use.convex) 
      || length(use.convex) != 1 
      || !(use.convex %in% c(TRUE, FALSE))){
    ddalpha$useConvex <- FALSE
    warning("Argument \"use.convex\" not specified correctly. FALSE is used as a default value")
  }else{
    ddalpha$useConvex <- use.convex
  }

  validate(ddalpha$methodDepth)
  
  ## The learning procedures

  # Calculate depths
  ddalpha <- .ddalpha.learn.depth(ddalpha)
  
  # Learn classification machine
  if (ddalpha$methodSeparator == "alpha"){
    ddalpha <- .ddalpha.learn.alpha(ddalpha)
  }
  if (ddalpha$methodSeparator == "polynomial"){
    ddalpha <- .ddalpha.learn.polynomial(ddalpha)
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

################################################################################
# Validation functions
################################################################################

alpha.validate  <- function(ddalpha, num.chunks = 10, max.degree = 3,...){

  if (ddalpha$methodAggregation == "majority"){
    maxChunks <- ddalpha$patterns[[ddalpha$numPatterns]]$cardinality + ddalpha$patterns[[ddalpha$numPatterns - 1]]$cardinality
  }else{
    maxChunks <- ddalpha$numPoints
  }
  
  if (is.character(num.chunks) && toupper(num.chunks)=="MAX")
    num.chunks <- maxChunks
  else if (!is.numeric(num.chunks) 
      || is.na(num.chunks) 
      || length(num.chunks) != 1 
      || !.is.wholenumber(num.chunks) 
      || !(num.chunks > 0 && num.chunks <= maxChunks)){
    num.chunks <- maxChunks
    warning("Argument \"num.chunks\" not specified correctly. ", maxChunks, " is used instead")
  }
  
  if(!is.numeric(max.degree) 
     || is.na(max.degree) 
     || length(max.degree) != 1 
     || !.is.wholenumber(max.degree) 
     || !(max.degree %in% 1:10)){
    max.degree <- 3
    warning("Argument \"max.degree\" not specified correctly. 3 is used as a default value")
  }
  
  return (list(numChunks = num.chunks, maxDegree = max.degree))  
}

polynomial.validate  <- alpha.validate  # the same  

knnlm.validate  <- function(ddalpha, knnrange = 10*( (ddalpha$numPoints)^(1/ddalpha$numPatterns) ) + 1,...){
  isnull = missing(knnrange) || is.null(knnrange)

  if (is.character(knnrange) && toupper(knnrange)=="MAX")
    knnrange = ceiling(ddalpha$numPoints/2)
  else if(is.null(knnrange)
     || !is.numeric(knnrange) 
     || is.na(knnrange) 
     || length(knnrange) != 1 
     || !.is.wholenumber(knnrange) 
     || !(knnrange >=2 && knnrange <= ceiling(ddalpha$numPoints/2))){    
    knnrange <- 10*( (ddalpha$numPoints)^(1/ddalpha$numPatterns) ) + 1   # Default
    
    knnrange <- min(knnrange, ceiling(ddalpha$numPoints/2))
    knnrange <- max(knnrange, 2)
    
    if (!isnull) warning("Argument \"knnrange\" not specified correctly. ", knnrange, " is used instead")
  }
  return (list(knnrange = knnrange))
}

randomTukey.validate  <- function(ddalpha, num.directions = 1000,...){
  if(!is.numeric(num.directions) 
     || is.na(num.directions) 
     || length(num.directions) != 1 
     || !.is.wholenumber(num.directions) 
     || !(num.directions > 1 && num.directions < 10000000) ){
    num.directions <- 1000
    warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
  }
  return (list(numDirections = num.directions))
}

Mahalanobis.validate  <- function(ddalpha, mah.estimate = "moment", mah.priors = NULL, mah.parMcd = 0.75, ...){ 
  if (!is.character(mah.estimate) 
      || length(mah.estimate) != 1 
      || !(mah.estimate %in% c("moment", "MCD"))){
    mah.estimate <- "moment"
    warning("Argument \"mah.estimate\" not specified correctly. \"moment\" is used as a default value")
  }
  
  if (!is.vector(mah.priors, mode = "double") 
      || is.na(min(mah.priors)) 
      || length(mah.priors) != ddalpha$numPatterns 
      || min(mah.priors) <= 0 
      || max(mah.priors) <= 0){
    if (!is.null(mah.priors)){
      warning("Argument \"mah.priors\" not specified correctly. Defaults in the form of class portions are applied")
    }
    mah.priors <- NULL
  }else{
    mah.priors <- mah.priors/sum(mah.priors)
  }
  
  ret <- list(mahEstimate = mah.estimate, mahPriors = mah.priors)
  
  if (mah.estimate == "MCD"){
    if (is.null(mah.parMcd) 
        || !is.vector(mah.parMcd, mode = "double") 
        || is.na(min(mah.parMcd)) 
        || length(mah.parMcd) != 1 
        || mah.parMcd < 0.5 
        || mah.parMcd > 1){
      mah.parMcd <- 0.75
      warning("Argument \"mah.parMcd\" not specified correctly. 0.75 is used as a default value")
    }
    ret$mahParMcd = mah.parMcd
  }
  return (ret)
}

