################################################################################
# File:             ddalpha.classify.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     15.05.2013
# 
# Contains the classification function of the DDalpha-classifier.
# 
# For a description of the algorithm, see:
#   Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric 
#     classification based on data depth. Statistical Papers.
#   Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world 
#     data with the DDalpha-procedure. Mimeo.
################################################################################

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
      if (ddalpha$methodDepth == "Mahalanobis"){
        depths <- .Mahalanobis_depths(ddalpha, objects[classifiableIndices,])
      }
      if (   ddalpha$methodDepth == "projectionRandom" 
          || ddalpha$methodDepth == "projectionLinearize"){
        depths <- .projection_depths(ddalpha, objects[classifiableIndices,])
      }
      if (ddalpha$methodDepth == "spatial"){
        depths <- .spatial_depths(ddalpha, objects[classifiableIndices,])
      }
      freePoints <- matrix(objects[-classifiableIndices,], 
                           nrow=nrow(objects)-length(classifiableIndices))
    }
  }else{
    if (ddalpha$methodDepth == "randomTukey"){
      depths <- .halfspace_depths(ddalpha, objects)
    }
    if (ddalpha$methodDepth == "zonoid"){
      depths <- .zonoid_depths(ddalpha, objects, 0)
    }
    if (ddalpha$methodDepth == "Mahalanobis"){
      depths <- .Mahalanobis_depths(ddalpha, objects)
    }
    if (   ddalpha$methodDepth == "projectionRandom" 
           || ddalpha$methodDepth == "projectionLinearize"){
      depths <- .projection_depths(ddalpha, objects)
    }
    if (ddalpha$methodDepth == "spatial"){
      depths <- .spatial_depths(ddalpha, objects)
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
      depths <- matrix(depths[classifiableIndices,], 
                       nrow=length(classifiableIndices), ncol=ddalpha$numPatterns)
      freePoints <- matrix(objects[-classifiableIndices,], 
                           nrow=nrow(objects)-length(classifiableIndices), 
                           ncol=ncol(objects))
    }
  }
  
  # Classify with the pure DD-ALpha
  resultsDepths <- list()
  if (nrow(depths) > 0){
    if (ddalpha$methodSeparator == "alpha"){
      votes <- matrix(rep(0, nrow(depths)*ddalpha$numPatterns), nrow=nrow(depths), ncol=ddalpha$numPatterns)
      toClassify <- as.double(as.vector(t(depths)))
      m <- as.integer(nrow(depths))
      q <- as.integer(ncol(depths))
      for (i in 1:ddalpha$numClassifiers){
        result <- .C("AlphaClassify", 
                     toClassify, 
                     m, 
                     q, 
                     as.integer(ddalpha$classifiers[[i]]$degree), 
                     as.double(ddalpha$classifiers[[i]]$hyperplane), 
                     output=integer(m))$output
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
    }
    if (ddalpha$methodSeparator == "knnlm"){
      z <- as.vector(t(depths))
      output <- .C("KnnClassify", 
                   as.double(z), 
                   as.integer(nrow(depths)), 
                   as.double(ddalpha$knnX), 
                   as.integer(ddalpha$knnY), 
                   as.integer(ddalpha$numPoints), 
                   as.integer(ddalpha$numPatterns), 
                   as.integer(ddalpha$knnK), 
                   as.integer(2), 
                   output=integer(nrow(depths)))$output
      for (i in 1:nrow(depths)){
        resultsDepths[[i]] <- ddalpha$patterns[[output[i] + 1]]$name
      }
      
#       dists <- ddalpha$knnD
# #      print(ddalpha$knnY)
#       for (i in 1:nrow(depths)){
#         newDists <- c()
# #        print(ddalpha$knnX)
# #        print(depths)
#         for (j in 1:nrow(ddalpha$knnX)){
#           newDists[j] <- knn.dist(rbind(ddalpha$knnX[j,], depths[i,]), dist.meth="maximum")[1,2]
#         }
# #        print(newDists)
#         dists <- cbind(ddalpha$knnD, newDists)
#         dists <- rbind(dists, c(newDists, 0))
# #        print(dists)
#         resultsDepths[[i]] <- as.numeric(knn.predict(1:length(ddalpha$knnY), 
#                                                      length(ddalpha$knnY) + 1, 
#                                                      ddalpha$knnY, 
#                                                      dists, k=ddalpha$knnK, 
#                                                      agg.meth="majority", 
#                                                      ties.meth="first"))
# #        print(resultsDepths[[i]])
#       }
    }
  }
  
  # Classify Outsiders
  resultsOutsiders <- as.list(rep("Ignored", nrow(freePoints)))
  if (length(resultsOutsiders) > 0 && !is.null(outsider.method)){
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
