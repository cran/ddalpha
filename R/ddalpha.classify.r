ddalpha.classify <- function(objects, 
                             ddalpha, 
                             outsider.method = "LDA", 
                             use.convex = NULL){
  # Checks
  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  if (!is.numeric(objects)){
    warning("Argument \"objects\" has unacceptable format. Classification can not be performed!!!")
    return (NULL)
  } 
  if (ncol(objects) != ddalpha$dimension){
    warning("Dimension of the objects to be classified does not correspond to the dimension of the trained classifier. Classification can not be performed!!!")
    return (NULL)
  }
  if (!is.character(outsider.method) 
      || length(outsider.method) != 1){
    warning("Argument \"outsidet.method\" not specified correctly. Outsiders will be ignored!!!")
    outsider.method <- NULL
  }
  if (is.null(use.convex)){
    use.convex <- ddalpha$useConvex
  }
  depths <- matrix(nrow=0, ncol=ddalpha$numPatterns) #?
  freePoints <- matrix(nrow=0, ncol=ncol(objects)) #?
  
  # Define points that can be classified by the DD-Alpha and the outsiders
  if (use.convex){
    points <- ddalpha$patterns[[1]]$points
    cardinalities <- c(ddalpha$patterns[[1]]$cardinality)
    for (i in 2:ddalpha$numPatterns){
      points <- rbind(points, ddalpha$patterns[[i]]$points)
      cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
    }
    classifiable <- .are_classifiable(objects, points, cardinalities)
    classifiableIndices <- which(classifiable == 1)
    if (length(classifiableIndices) == 0){
      depths <- matrix(nrow=0, ncol=ddalpha$numPatterns)
      freePoints <- objects
    }else{
      if (ddalpha$methodDepth == "randomTukey"){
        depths <- .halfspace_depths(ddalpha, objects[classifiableIndices,])
      }
      if (ddalpha$methodDepth == "zonoid"){
        depths <- .zonoid_depths(ddalpha, objects[classifiableIndices,], 0)
      }
      freePoints <- objects[-classifiableIndices,]
    }
  }else{
    if (ddalpha$methodDepth == "randomTukey"){
      depths <- .halfspace_depths(ddalpha, objects)
    }
    if (ddalpha$methodDepth == "zonoid"){
      depths <- .zonoid_depths(ddalpha, objects, 0)
    }
    classifiableIndices <- c()
    for (i in 1:nrow(depths)){
      if (sum(depths[i,]) > 0){
        classifiableIndices <- c(classifiableIndices, i)
      }
    }
    if (length(classifiableIndices) == 0){
      depths <- matrix(nrow=0, ncol=ddalpha$numPatterns)
      freePoints <- objects
    }else{
      depths <- depths[classifiableIndices,]
      freePoints <- objects[-classifiableIndices,]
    }
  }
  
  # Classify with the pure DD-ALpha
  resultsDepths <- list()
  votes <- matrix(rep(0, nrow(depths)*ddalpha$numPatterns), nrow=nrow(depths))
  toClassify <- as.double(as.vector(t(depths)))
  m <- as.integer(nrow(depths))
  q <- as.integer(ncol(depths))
  for (i in 1:ddalpha$numClassifiers){
    result <- .C("AlphaClassify", toClassify, m, q, as.integer(ddalpha$classifiers[[i]]$degree), as.double(ddalpha$classifiers[[i]]$hyperplane), output=integer(m))$output
    for (j in 1:m){
      if (result[j] > 0){
        votes[j,ddalpha$classifiers[[i]]$index0] <- votes[j,ddalpha$classifiers[[i]]$index0] + 1
      }else{
        votes[j,ddalpha$classifiers[[i]]$index1] <- votes[j,ddalpha$classifiers[[i]]$index1] + 1
      }
    }
  }
  for (i in 1:m){
    resultsDepths[[i]] <- ddalpha$patterns[[which.max(votes[i,])]]$name
  }
  
  # Classify Outsiders
  resultsOutsiders <- as.list(rep("Ignored", nrow(freePoints)))
  if (!is.null(outsider.method)){
    for (i in 1:length(ddalpha$methodsOutsider)){
      if (toupper(ddalpha$methodsOutsider[[i]]$name) == toupper(outsider.method)){
        resultsOutsiders <- .ddalpha.classify.outsiders(freePoints, ddalpha, ddalpha$methodsOutsider[[i]])
        break
      }
    }
  }
  
  # Merge classifiable and outsiders
  results <- list()
  counterDepths <- 1
  counterOutsiders <- 1
  for (i in 1:nrow(objects)){
    if (i %in% classifiableIndices){
      results[[i]] <- resultsDepths[[counterDepths]]
      counterDepths <- counterDepths + 1
    }else{
      results[[i]] <- resultsOutsiders[[counterOutsiders]]
      counterOutsiders <- counterOutsiders + 1
    }
  }
  return (results)
}
