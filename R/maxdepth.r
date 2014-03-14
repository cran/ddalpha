maxdepth.createStructure <- function(data){

  # Elemantary statistics
  dimension <- ncol(data) - 1
  numOfPoints <- nrow(data)
  classNames <- unique(data[,dimension + 1])
  numOfClasses <- length(classNames)
  
  # Ordering patterns according to their cardinalities
  classCardinalities <- rep(0, numOfClasses)
  for (i in 1:numOfClasses){
    classCardinalities[i] <- nrow(data[data[,dimension + 1] == classNames[i],])
  }
  
  # Creating pattern templates
  patterns <- as.list("")
  for (i in 1:numOfClasses){
    maxCarIndex <- which.max(classCardinalities)
    # Creating a single template
    pattern.index       <- i
    pattern.points      <- data[data[,dimension + 1] == 
                                  classNames[maxCarIndex],1:dimension]
    pattern.name        <- classNames[maxCarIndex]
    pattern.cardinality <- classCardinalities[maxCarIndex]
    pattern.depths      <- matrix(
      rep(0, numOfClasses*classCardinalities[maxCarIndex]), 
      nrow = classCardinalities[maxCarIndex], ncol = numOfClasses)
    pattern.center      <- NULL
    pattern.cov         <- NULL
    pattern.sigma       <- NULL
    pattern.centerMcd   <- NULL
    pattern.covMcd      <- NULL
    pattern.sigmaMcd    <- NULL
    pattern.votes       <- 0
    pattern <- structure(
      list(index = pattern.index, 
        points = pattern.points, 
        name = pattern.name, 
        cardinality = pattern.cardinality, 
        depths = pattern.depths, 
        center = pattern.center, 
        cov = pattern.cov, 
        sigma = pattern.sigma, 
        centerMcd = pattern.centerMcd, 
        covMcd = pattern.covMcd, 
        sigmaMcd = pattern.sigmaMcd, 
        votes = pattern.votes), 
      .Names = c("index", "points", "name", "cardinality", "depths", "center", 
                 "cov", "sigma", "centerMcd", "covMcd", "sigmaMcd", 
                 "votes"))
    # Adding pattern template to the list of patterns
    patterns[[i]] <- pattern
    # Deleting processed pattern
    classCardinalities[maxCarIndex] <- -1    
  }
  # Creating overall structure
  maxdepth <- structure(
    list(dimension = dimension,          
         numPatterns = numOfClasses, 
         numPoints = numOfPoints, 
         patterns = patterns, 
         classifiers = as.list(""), 
         numClassifiers = 0, 
         methodDepth = "majority", 
         methodOutsider = "RandProp", 
         maxDepth = "Mahalanobis", 
         proirs = NULL, 
         mcdEstimate = "moment", 
         mcdAlpha = 0.75, 
         numDirections = 1000), 
    .Names = c("dimension", "numPatterns", "numPoints", "patterns", 
               "classifiers", "numClassifiers", "methodDepth", 
               "methodOutsider", "maxDepth", "priors", "mcdEstimate", 
               "mcdAlpha", "numDirections"))
  
  return (maxdepth)
}

maxdepth.fillStructure <- function(maxdepth){

  if (is.null(maxdepth$priors)){
    maxdepth$priors <- c()
    for (i in 1:maxdepth$numPatterns){
      maxdepth$priors[i] <- maxdepth$patterns[[i]]$cardinality/maxdepth$numPoints
    }
  }
  if (maxdepth$maxDepth == "Mahalanobis"){
    for (i in 1:maxdepth$numPatterns){
      # Calculating raw and robust estimates (package "robustbase" needed)
      if (maxdepth$mcdEstimate == "moment"){
        maxdepth$patterns[[i]]$center <- colMeans(maxdepth$patterns[[i]]$points)
        maxdepth$patterns[[i]]$cov    <- cov(maxdepth$patterns[[i]]$points)
        try(
          maxdepth$patterns[[i]]$sigma  <- solve(maxdepth$patterns[[i]]$cov)
        )
      }
      if (maxdepth$mcdEstimate == "MCD"){
        try(
          estimate <- covMcd(maxdepth$patterns[[i]]$points, maxdepth$mcdAlpha)
        )
        try(
          maxdepth$patterns[[i]]$centerMcd <- estimate$center
        )
        try(
          maxdepth$patterns[[i]]$covMcd    <- estimate$cov
        )
        try(
          maxdepth$patterns[[i]]$sigmaMcd  <- solve(estimate$cov)
        )
      }
    }
  }

  return (maxdepth)
}

maxdepth.train <- function(data, depth = "Mahalanobis", 
                           outsider.method = "RandProp", 
                           priors = NULL, 
                           mah.estimate = "moment", 
                           mah.parMcd = 0.75, 
                           num.directions = 1000){

  maxdepth <- maxdepth.createStructure(data)
  maxdepth$maxDepth <- depth
  maxdepth$methodOutsider <- outsider.method
  # Checks
#   if (!is.character(depth)
#       || length(depth) != 1
#       || !(depth %in% c("Mahalanobis", "spatial", "Tukey", "Liu", "Oja", "projectionRandom"))){
#     maxdepth$maxDepth <- "randomTukey"
#     warning("Argument \"depth\" not specified correctly. \"Mahalanobis\" is used as a default value")
#   }else{
     maxdepth$maxDepth <- depth
#   }
#   if (!is.character(mah.estimate) 
#       || length(mah.estimate) != 1 
#       || !(mah.estimate %in% c("moment", "MCD"))){
#     warning("Argument \"mah.estimate\" not specified correctly. \"moment\" is used as a default value")
#     maxdepth$mcdEstimate <- "moment"
#   }else{
     maxdepth$mcdEstimate <- mah.estimate
#   }
  if (!is.vector(priors, mode = "double") 
      || is.na(min(priors)) 
      || length(priors) != maxdepth$numPatterns 
      || min(priors) <= 0 
      || max(priors) <= 0){
    if (!is.null(priors)){
      warning("Argument \"priors\" not specified correctly. Defaults in the form of class portions are applied")
    }
    maxdepth$priors <- NULL
  }else{
    maxdepth$priors <- priors/sum(priors)
  }
#   if (!is.vector(mah.parMcd, mode = "double") 
#       || is.na(min(mah.parMcd)) 
#       || length(mah.parMcd) != 1 
#       || mah.parMcd < 0.5 
#       || mah.parMcd > 1){
#     if (tmp.mah.estimate == "MCD"){
#       warning("Argument \"mah.parMcd\" not specified correctly. 0.75 is used as a default value")
#     }
#     maxdepth$mcdAlpha <- 0.75
#   }else{
     maxdepth$mcdAlpha <- mah.parMcd
#   }
#   if ((maxdepth$methodDepth == "randomTukey"
#        || maxdepth$methodDepth == "projectionRandom")
#       && (!is.numeric(num.directions) 
#           || is.na(num.directions) 
#           || length(num.directions) != 1 
#           || !.is.wholenumber(num.directions) 
#           || !(num.directions > 1 && num.directions < 10000000))
#   ){
#     maxdepth$numDirections <- 1000
#     warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
#   }else{
     maxdepth$numDirections <- num.directions
#   }
  # Preparing data
  maxdepth <- maxdepth.fillStructure(maxdepth)

  return (maxdepth)
}

maxdepth.classify <- function(objects, maxdepth, outsider.method = "RandProp"){

  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  output <- as.list("")
  
  # Calculating depths of the objects according to each pattern
  depths <- matrix(rep(0, maxdepth$numPatterns*nrow(objects)), nrow=nrow(objects), ncol=maxdepth$numPatterns)
  for (i in 1:maxdepth$numPatterns){
    depths[,i] <- maxdepth.depth(objects, maxdepth, i)
  }
  for (i in 1:nrow(objects)){
    # Preparing voting system
    votes <- rep(0, maxdepth$numPatterns)
    # Checking if the object is an outsider (lies outside the convex hull of all patterns)
    if (length(depths[i,depths[i,] <= 0]) == maxdepth$numPatterns){
      # Outsider detected! Classify
      if (maxdepth$methodOutsider == "Ignore"){
        output[[i]] <- "Ignored"
        next
      }
      if (maxdepth$methodOutsider == "Right"){
        output[[i]] <- "Right"
        next
      }
      if (maxdepth$methodOutsider == "RandProp"){
        curMin <- 0
        rndVal <- runif(1) * maxdepth$numPoints
        for (j in 1:maxdepth$numPatterns){
          tmpVal <- rndVal - curMin
          if (tmpVal <= maxdepth$patterns[[j]]$cardinality || j == maxdepth$numPatterns){
            output[[i]] <- maxdepth$patterns[[j]]$name
          }else{
            curMin <- curMin + maxdepth$patterns[[j]]$cardinality
          }
        }      
      }
      if (maxdepth$methodOutsider == "RandEqual"){
        curMin <- 0
        rndVal <- sample(1:maxdepth$numPatterns, 1)
        output[[i]] <- maxdepth$patterns[[rndVal]]$name
      }
    }else{
      # Object is inside the convex hull at least of one of the patterns. Maximum depth classifier can be applied    
      for (j in 1:maxdepth$numPatterns){
        votes[j] <- depths[i,j]*maxdepth$priors[j]
      }
    }
    output[[i]] <- maxdepth$patterns[[which.max(votes)]]$name
  }
  
#   for (i in 1:nrow(objects)){
#     output[[i]] <- maxdepth.classifyOne(objects[i,], maxdepth)
#   }

  return (output)
}

maxdepth.numOutsiders <- function(objects, maxdepth){

  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  outsiders <- as.vector(TRUE)
  for (i in 1:nrow(objects)){
    outsiders[i] <- maxdepth.isOutsider(objects[i,], maxdepth)
  }

  return (length(outsiders[outsiders == TRUE]))
}

maxdepth.isOutsider <- function(object, maxdepth){

  # Calculating depths of the object w.r.t. each pattern
  depths <- rep(0, maxdepth$numPatterns)
  for (i in 1:maxdepth$numPatterns){
    depths[i] <- maxdepth.depth(object, maxdepth, i)
  }
  
  # Checking if the object is an outsider (lies outside the convex hull of all patterns)
  if (length(depths[depths <= 0]) == maxdepth$numPatterns){
    return (TRUE)
  }else{
    return (FALSE)
  }
}

maxdepth.classifyOne <- function(object, maxdepth){
  # Classifying a single object
  
  # Calculating depths of the object according to each pattern
  depths <- rep(0, maxdepth$numPatterns)
  for (i in 1:maxdepth$numPatterns){
    depths[i] <- maxdepth.depth(object, maxdepth, i)
  }
  # Preparing voting system
  votes <- rep(0, maxdepth$numPatterns)
  # Checking if the object is an outsider (lies outside the convex hull of all patterns)
  if (length(depths[depths <= 0]) == maxdepth$numPatterns){
    # Outsider detected! Classify
    if (maxdepth$methodOutsider == "Ignore"){
      return ("Ignored")
    }
    if (maxdepth$methodOutsider == "Right"){
      return ("Right")
    }
    if (maxdepth$methodOutsider == "RandProp"){
      curMin <- 0
      rndVal <- runif(1) * maxdepth$numPoints
      for (i in 1:maxdepth$numPatterns){
        tmpVal <- rndVal - curMin
        if (tmpVal <= maxdepth$patterns[[i]]$cardinality || i == maxdepth$numPatterns){
          return (maxdepth$patterns[[i]]$name)
        }else{
          curMin <- curMin + maxdepth$patterns[[i]]$cardinality
        }
      }      
    }
    if (maxdepth$methodOutsider == "RandEqual"){
      curMin <- 0
      rndVal <- sample(1:maxdepth$numPatterns, 1)
      return (maxdepth$patterns[[rndVal]]$name)
    }
  }else{
    # Object is inside the convex hull at least of one of the patterns. Maximum depth classifier can be applied    
    for (i in 1:maxdepth$numPatterns){
      votes[i] <- depths[i]*maxdepth$priors[i]
    }
  }

  return(maxdepth$patterns[[which.max(votes)]]$name)
}

################################################################################
# Functions for intermediate calculations are presented below
################################################################################

maxdepth.depth <- function(object, maxdepth, patternIndex){

  depth <- -1
  if (maxdepth$maxDepth == "Mahalanobis"){
    if (maxdepth$mcdEstimate == "moment"){
      depth <- .Mahalanobis_depth(object, 
                                  maxdepth$patterns[[patternIndex]]$center, 
                                  maxdepth$patterns[[patternIndex]]$sigma)
    }
    if (maxdepth$mcdEstimate == "MCD"){
      depth <- .Mahalanobis_depth(object, 
                                  maxdepth$patterns[[patternIndex]]$centerMcd, 
                                  maxdepth$patterns[[patternIndex]]$sigmaMcd)
    }
  }
  if (maxdepth$maxDepth == "Tukey"){
    depth <- depth(object, maxdepth$patterns[[patternIndex]]$points, method = "Tukey")
  }
  if (maxdepth$maxDepth == "Liu"){
    depth <- depth(object, maxdepth$patterns[[patternIndex]]$points, method = "Liu")
  }
  if (maxdepth$maxDepth == "Oja"){
    depth <- depth(object, maxdepth$patterns[[patternIndex]]$points, method = "Oja")
  }
  if (maxdepth$maxDepth == "spatial"){
    depth <- depth.spatial(object, maxdepth$patterns[[patternIndex]]$points)
  }
  if (maxdepth$maxDepth == "projectionRandom"){
    depth <- depth.projection(maxdepth$patterns[[patternIndex]]$points, object, num.directions = maxdepth$numDirections)
  }
  return (depth)
}
