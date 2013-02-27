.ddalpha.create.structure <- function(data){
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
    pattern.points      <- data[data[,dimension + 1] == classNames[maxCarIndex],1:dimension]
    pattern.name        <- classNames[maxCarIndex]
    pattern.cardinality <- classCardinalities[maxCarIndex]
    pattern.depths      <- matrix(rep(0, numOfClasses*classCardinalities[maxCarIndex]), nrow = classCardinalities[maxCarIndex], ncol = numOfClasses)
    pattern.votes       <- 0
    pattern <- structure(
      list(index = pattern.index, 
        points = pattern.points, 
        name = pattern.name, 
        cardinality = pattern.cardinality, 
        depths = pattern.depths, 
        votes = pattern.votes), 
      .Names = c("index", "points", "name", "cardinality", "depths", "votes"))
    # Adding pattern template to the list of patterns
    class(pattern)<-"ddalpha.pattern"
    patterns[[i]] <- pattern
    # Deleting processed pattern
    classCardinalities[maxCarIndex] <- -1    
  }
  # Creating overall structure
  ddalpha <- structure(
    list(raw <- data, 
         dimension = dimension, 
         numPatterns = numOfClasses, 
         numPoints = numOfPoints, 
         patterns = patterns, 
         classifiers = list(), 
         numClassifiers = 0, 
         methodDepth <- "randomTukey",
         methodAggregation = "majority", 
         methodsOutsider = NULL, 
         numDirections = 1000, 
         directions <- NULL, 
         projections <- NULL, 
         sameDirections <- TRUE, 
         useConvex <- FALSE, 
         maxDegree = 3, 
         numChunks = numOfPoints, 
         treatments <- c("LDA", "KNNAFF", "KNN", "DEPTH.MAHALANOBIS", "RANDEQUAL", "RANDPROP", "IGNORE")), 
    .Names = c("raw", "dimension", "numPatterns", "numPoints", "patterns", 
               "classifiers", "numClassifiers", "methodDepth", "methodAggregation", "methodsOutsider", 
               "numDirections", "directions", "projections", "sameDirections", "useConvex", 
               "maxDegree", "numChunks", "treatments"))

  return (ddalpha)
}

.ddalpha.learn.depth <- function(ddalpha){
  # If it's the random Tukey depth, compute it first
  if (ddalpha$methodDepth == "randomTukey"){
    dSpaceStructure <- .halfspace_space(ddalpha)
    ddalpha$directions <- dSpaceStructure$directions
    ddalpha$projections <- dSpaceStructure$projections
    tmpDSpace <- dSpaceStructure$dspace
  }
  # Calculating depths in each pattern
  for (i in 1:ddalpha$numPatterns){
    if (ddalpha$methodDepth == "randomTukey"){
      # The random Tukey depth is already calculated, just distribute
      ddalpha$patterns[[i]]$depths <- tmpDSpace[1:ddalpha$patterns[[i]]$cardinality,]
      tmpDSpace <- tmpDSpace[-(1:ddalpha$patterns[[i]]$cardinality),]
    }
    if (ddalpha$methodDepth == "zonoid"){
      # Calculate depths for the class w.r.t all classes, saying to which of the classes the chunk belongs
      ddalpha$patterns[[i]]$depths <- .zonoid_depths(ddalpha, ddalpha$patterns[[i]]$points, i)
    }
  }

  return (ddalpha)
}

.ddalpha.learn.alpha <- function(ddalpha){
  # Separating (calculating extensions and normals)
  counter <- 1
  # Determining multi-class behaviour
  if (ddalpha$methodAggregation == "majority"){
    for (i in 1:(ddalpha$numPatterns - 1)){
      for (j in (i + 1):ddalpha$numPatterns){
        # Creating a classifier
        hyperplane <- .alpha_learn(ddalpha$maxDegree, rbind(ddalpha$patterns[[i]]$depths, ddalpha$patterns[[j]]$depths), ddalpha$patterns[[i]]$cardinality, ddalpha$patterns[[j]]$cardinality, ddalpha$numChunks)
        classifier.index          <- counter
        classifier.index0         <- i
        classifier.index1         <- j
        classifier.hyperplane     <- hyperplane$por
        classifier.degree         <- hyperplane$deg
        classifier.dimProperties  <- length(hyperplane$por) - 1
        classifier.dimFeatures    <- hyperplane$dim
        # Adding the classifier to the list of classifiers
        ddalpha$classifiers[[counter]] <- structure(
          list(index = classifier.index, 
               index0 = classifier.index0, 
               index1 = classifier.index1, 
               hyperplane = classifier.hyperplane, 
               degree = classifier.degree, 
               dimProperties = classifier.dimProperties, 
               dimFeatures = classifier.dimFeatures), 
          .Names = c("index", "index0", "index1", "hyperplane", "degree", "dimProperties", "dimFeatures"))
        counter <- counter + 1
      }
    }
    ddalpha$numClassifiers <- counter - 1
  }
  if (ddalpha$methodAggregation == "sequent"){
  # !!!DEBUG!!! - Do not forget about the same structure in the knnAff
    for (i in 1:ddalpha$numPatterns){
      anotherClass <- NULL
      for (j in 1:ddalpha$numPatterns){
        if (j != i){
          anotherClass <- rbind(anotherClass, ddalpha$patterns[[j]]$depths)
        }
      }
      hyperplane <- .alpha_learn(ddalpha$maxDegree, rbind(ddalpha$patterns[[i]]$depths, anotherClass), ddalpha$patterns[[i]]$cardinality, nrow(anotherClass), ddalpha$numChunks)
      classifier.index          <- counter
      classifier.index0         <- i
      classifier.index1         <- -1
      classifier.hyperplane     <- hyperplane$por
      classifier.degree         <- hyperplane$deg
      classifier.dimProperties  <- length(hyperplane$por) - 1
      classifier.dimFeatures    <- hyperplane$dim
      # Adding the classifier to the list of classifiers
      ddalpha$classifiers[[i]] <- structure(
        list(index = classifier.index, 
             index0 = classifier.index0, 
             index1 = classifier.index1, 
             hyperplane = classifier.hyperplane, 
             degree = classifier.degree, 
             dimProperties = classifier.dimProperties, 
             dimFeatures = classifier.dimFeatures), 
        .Names = c("index", "index0", "index1", "hyperplane", "degree", "dimProperties", "dimFeatures"))
    }
    ddalpha$numClassifiers <- ddalpha$numPatterns
  }

  return (ddalpha)
}

.ddalpha.learn.outsiders <- function(ddalpha, methodsOutsider = "LDA", settingsOutsider = NULL){
  # Refine treatments
  if (is.null(settingsOutsider)){
    ddalpha$methodsOutsider <- .parse.methods(methodsOutsider)
  }else{
    ddalpha$methodsOutsider <- .parse.settings(ddalpha, settingsOutsider)
  }
  # Train treatments
  for (i in 1:length(ddalpha$methodsOutsider)){
    if (ddalpha$methodsOutsider[[i]]$method == "LDA"){
      ddalpha$methodsOutsider[[i]] <- .lda_learn(ddalpha, ddalpha$methodsOutsider[[i]])
      next
    }
    if (ddalpha$methodsOutsider[[i]]$method == "KNN"){
      ddalpha$methodsOutsider[[i]] <- .knn_learn(ddalpha, ddalpha$methodsOutsider[[i]])
      next
    }
    if (ddalpha$methodsOutsider[[i]]$method == "KNNAff"){
      ddalpha$methodsOutsider[[i]] <- .knnAff_learn(ddalpha, ddalpha$methodsOutsider[[i]])
      next
    }
    if (ddalpha$methodsOutsider[[i]]$method == "depth.Mahalanobis"){
      ddalpha$methodsOutsider[[i]] <- .mah_learn(ddalpha, ddalpha$methodsOutsider[[i]])
      next
    }
  }
  
  return(ddalpha)
}

.ddalpha.classify.outsiders<- function (objects, ddalpha, settings){
  if (settings$method == "LDA"){
    return (.lda_classify(objects, ddalpha, settings))
  }
  if (settings$method == "KNN"){
    return (.knn_classify(objects, ddalpha, settings))
  }
  if (settings$method == "KNNAff"){
    return (.knnAff_classify(objects, ddalpha, settings))
  }
  if (settings$method == "depth.Mahalanobis"){
    return (.mah_classify(objects, ddalpha, settings))
  }
  if (settings$method == "Ignore"){
    return (.ignore_classify(nrow(objects)))
  }
  if (settings$method == "Mark"){
    return (.right_classify(nrow(objects)))
  }
  if (settings$method == "RandEqual"){
    return (.randequal_classify(nrow(objects), ddalpha))
  }
  if (settings$method == "RandProp"){
    return (.randprop_classify(nrow(objects), ddalpha, settings))
  }
}

################################################################################
# Functions used for intermediate calculations and checks are presented below
################################################################################

.are_classifiable <- function(objects, points, cardinalities){
  x <- as.vector(t(points))
  dimension <- ncol(points)
  numClasses <- length(cardinalities)
  o <- as.vector(t(objects))
  numObjects <- nrow(objects)
  result <- .C("IsInConvexes", 
               as.double(x), 
               as.integer(dimension), 
               as.integer(cardinalities), 
               as.integer(numClasses), 
               as.double(o), 
               as.integer(numObjects), 
               isInConvexes=integer(numObjects))$isInConvexes
  return (result)
}

.halfspace_space <- function(ddalpha){
  points <- NULL
  cardinalities <- NULL
  for (i in 1:ddalpha$numPatterns){
    points <- rbind(points, ddalpha$patterns[[i]]$points)
    cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
  }
  x <- as.vector(t(points))
  c <- as.vector(cardinalities)
  k <- ddalpha$numDirections
#  sameDirs <- 1 #?
  if (ddalpha$sameDirections){
    rez <- .C("HDSpace", as.double(x), as.integer(ncol(points)), as.integer(c), as.integer(ddalpha$numPatterns), as.integer(k), as.integer(1), dspc=double(nrow(points)*ddalpha$numPatterns), dirs=double(k*ncol(points)), prjs=double(k*nrow(points)))
    return (list(dspace=matrix(rez$dspc, nrow=nrow(points), ncol=ddalpha$numPatterns, byrow=TRUE), directions=rez$dirs, projections=rez$prjs))
  }else{
    rez <- .C("HDSpace", as.double(x), as.integer(ncol(points)), as.integer(c), as.integer(ddalpha$numPatterns), as.integer(k), as.integer(0), dspc=double(nrow(points)*ddalpha$numPatterns), dirs=double(1), prjs=double(1))
    return (list(dspace=matrix(rez$dspc, nrow=nrow(points), ncol=ddalpha$numPatterns, byrow=TRUE), directions=0, projections=0))
  }
}

.halfspace_depths <- function(ddalpha, objects){
  points <- NULL
  cardinalities <- NULL
  for (i in 1:ddalpha$numPatterns){
    points <- rbind(points, ddalpha$patterns[[i]]$points)
    cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
  }
  x <- as.vector(t(points))
  y <- as.vector(t(objects))
  c <- as.vector(cardinalities)
  k <- ddalpha$numDirections
  if (ddalpha$sameDirections){
    result <- .C("HDepth", 
                 as.double(x), 
                 as.double(y), 
                 as.integer(nrow(objects)), 
                 as.integer(ncol(points)), 
                 as.integer(c), 
                 as.integer(ddalpha$numPatterns), 
                 as.double(ddalpha$directions), 
                 as.double(ddalpha$projections), 
                 as.integer(k), 
                 as.integer(1), 
                 depths=double(ddalpha$numPatterns*nrow(objects)))
  }else{
    result <- .C("HDepth", 
                 as.double(x), 
                 as.double(y), 
                 as.integer(nrow(objects)), 
                 as.integer(ncol(points)), 
                 as.integer(c), 
                 as.integer(ddalpha$numPatterns), 
                 as.double(0), 
                 as.double(0), 
                 as.integer(k), 
                 as.integer(0), 
                 depths=double(ddalpha$numPatterns*nrow(objects)))
  }
  return (matrix(result$depths, nrow=nrow(objects), ncol=ddalpha$numPatterns, byrow=TRUE))
}

.zonoid_depths <- function(ddalpha, objects, ownPattern = 0){
  depths <- NULL
  for (i in 1:ddalpha$numPatterns){
    pattern <- ddalpha$patterns[[i]]$points
    x <- as.vector(t(pattern))
    y <- as.vector(t(objects))
    ds <- .C("ZDepth", as.double(x), as.double(y), as.integer(nrow(pattern)), as.integer(nrow(objects)), as.integer(ncol(pattern)), depths=double(nrow(objects)))$depths
    if (i == ownPattern){
      ds <- replace(ds, which(ds < 1/nrow(pattern)), 1/nrow(pattern))
    }else{
      ds <- replace(ds, which(ds < 1/nrow(pattern)), 0)
    }
    depths <- cbind(depths, ds)
  }
  return (depths)
}

.mahalanobis_depth <- function(point, center, sigma){
  d <- 1/(1 + (point - center) %*% sigma %*% (point - center))
  return (d)
}

.alpha_learn <- function(maxDegree, data, numClass1, numClass2, numChunks){
  points <- as.vector(t(data))
  numPoints <- numClass1 + numClass2
  dimension <- ncol(data)
  cardinalities <- c(numClass1, numClass2)
  upToPower <- maxDegree
  minFeatures <- 2
  maxExtDimension <- (factorial(dimension + maxDegree) / (factorial(dimension)*factorial(maxDegree))) - 1;
  
  p <- .C("AlphaLearnCV", as.double(points), as.integer(numPoints), as.integer(dimension), as.integer(cardinalities),  as.integer(upToPower), as.integer(numChunks), as.integer(minFeatures), portrait=double(maxExtDimension + 1))$portrait
  
  degree <- p[1];
  extDimension <- (factorial(dimension + degree) / (factorial(dimension)*factorial(degree))) - 1;
  d <- 0
  for (i in 2:(extDimension + 1)){
    if (p[i] != 0){
      d <- d + 1
    }
  }
  return(list(por=p,dim=d,deg=degree))
}

.parse.methods <- function(methods){
  if (!is.vector(methods) || length(methods) <= 0){return(list())}
  methods.refined <- unique(toupper(methods))
  
  treatments.settings <- list()
  counter <- 1
  for (i in 1:length(methods.refined)){
    supported <- FALSE
    if (methods.refined[i] == "LDA"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "LDA", 
             method = "LDA", 
             priors = NULL, 
             lda = NULL), 
        .Names = c("name", "method", "priors", "lda"))
    }
    if (methods.refined[i] == "KNNAFF"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "KNNAff", 
             method = "KNNAff", 
             knnAff.methodAggregation = "majority", 
             knnAff.range = -1, 
             knnAff.k = -1, 
             knnAff.classifiers = NULL), 
        .Names = c("name", "method", "knnAff.methodAggregation", "knnAff.range", "knnAff.k", "knnAff.classifiers"))
    }
    if (methods.refined[i] == "KNN"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "KNN", 
             method = "KNN", 
             knn.range = -1, 
             knn.k = -1, 
             knn.train = NULL, 
             knn.cl = NULL), 
        .Names = c("name", "method", "knn.range", "knn.k", "knn.train", "knn.cl"))
    }
    if (methods.refined[i] == "DEPTH.MAHALANOBIS"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "depth.Mahalanobis", 
             method = "depth.Mahalanobis", 
             mah.estimate = "moment", 
             priors = NULL, 
             mah.classes = NULL, 
             mcd.alpha = 0.5), 
        .Names = c("name", "method", "mah.estimate", "priors", "mah.classes", "mcd.alpha"))
    }
    if (methods.refined[i] == "RANDEQUAL"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "RandEqual", 
             method = "RandEqual"), 
        .Names = c("name", "method"))
    }
    if (methods.refined[i] == "RANDPROP"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "RandProp", 
             method = "RandProp", 
             priors = NULL), 
        .Names = c("name", "method", "priors"))
    }
    if (methods.refined[i] == "IGNORE"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "Ignore", 
             method = "Ignore"), 
        .Names = c("name", "method"))
    }
    if (supported){
      treatments.settings[[counter]] <- treatment.settings
      counter <- counter + 1
    }
  }
  return(treatments.settings)
}

.parse.settings <- function(ddalpha, settings){
  if (!is.list(settings) || length(settings) <= 0){return(list())}
  treatments.names <- c()
  
  treatments.settings <- list()
  counter <- 1
  for (i in 1:length(settings)){
    supported <- FALSE
    if (!is.list(settings[[i]]) 
        || is.null(settings[[i]]$name) 
        || is.null(settings[[i]]$method)){
      warning("In treatment number ", i, ": The treatment has unacceptable format. The treatment will be ignored")
      next
    }
    if (!is.character(settings[[i]]$method) 
        || length(settings[[i]]$method) != 1 
        || !(toupper(settings[[i]]$method) %in% ddalpha$treatments)){
      warning("In treatment number ", i, ": The method name of the treatment is not acceptable. The treatment will be ignored")
      next
    }
    if (!is.character(settings[[i]]$name) 
        || length(settings[[i]]$name) != 1){
      warning("In treatment number ", i, ": The name of the treatment is not acceptable. The treatment will be ignored")
      next
    }
    if (settings[[i]]$name %in% treatments.names){
      warning("In treatment number ", i, ": Treatment with the name ", settings[[i]]$name, " already exists. The treatment will be ignored")
      next
    }else{
      treatments.names <- c(treatments.names, settings[[i]]$name)
    }
    if (toupper(settings[[i]]$method) == "LDA"){
      tmp.name <- settings[[i]]$name
      if (!is.vector(settings[[i]]$priors, mode = "double") 
          || is.na(min(settings[[i]]$priors)) 
          || length(settings[[i]]$priors) != ddalpha$numPatterns 
          || min(settings[[i]]$priors) <= 0 
          || max(settings[[i]]$priors) <= 0){
        warning("In treatment number ", i, ": Argument \"priors\" not specified correctly. Defaults in the form of class portions are applied")
        tmp.priors <- NULL
      }else{
        tmp.priors <- settings[[i]]$priors/sum(settings[[i]]$priors)
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "LDA", 
             priors = tmp.priors, 
             lda = NULL), 
        .Names = c("name", "method", "priors", "lda"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "KNNAFF"){
      tmp.name <- settings[[i]]$name
      if (!is.character(settings[[i]]$knnAff.methodAggregation) 
          || length(settings[[i]]$knnAff.methodAggregation) != 1 
          || !(settings[[i]]$knnAff.methodAggregation %in% c("majority", "sequent"))){
        warning("In treatment number ", i, ": Argument \"knnAff.methodAggregation\" not specified correctly. \"majority\" is used as a default value")
        tmp.knnAff.methodAggregation <- "majority"
      }else{
        tmp.knnAff.methodAggregation <- settings[[i]]$knnAff.methodAggregation
      }
      if (!is.numeric(settings[[i]]$knnAff.range) 
          || is.na(settings[[i]]$knnAff.range) 
          || length(settings[[i]]$knnAff.range) != 1 
          || !.is.wholenumber(settings[[i]]$knnAff.range) 
          || !(settings[[i]]$knnAff.range >= 2 
               && settings[[i]]$knnAff.range <= (ddalpha$patterns[[ddalpha$numPatterns]]$cardinality + ddalpha$patterns[[ddalpha$numPatterns - 1]]$cardinality - 1) 
               || settings[[i]]$knnAff.range == -1)){
        warning("In treatment number ", i, ": Argument \"knnAff.range\" not specified correctly. Defaults are applied")
        tmp.knnAff.range <- -1
      }else{
        tmp.knnAff.range <- settings[[i]]$knnAff.range
      }
      if (!is.numeric(settings[[i]]$knnAff.k) 
          || is.na(settings[[i]]$knnAff.k) 
          || length(settings[[i]]$knnAff.k) != 1 
          || !.is.wholenumber(settings[[i]]$knnAff.k) 
          || !(settings[[i]]$knnAff.k >= 1 
               && settings[[i]]$knnAff.k <= (ddalpha$patterns[[ddalpha$numPatterns]]$cardinality + ddalpha$patterns[[ddalpha$numPatterns - 1]]$cardinality) 
               || settings[[i]]$knnAff.k == -1)){
        warning("In treatment number ", i, ": Argument \"knnAff.k\" not specified correctly. Defaults are applied")
        tmp.knnAff.k <- -1
      }else{
        tmp.knnAff.k <- settings[[i]]$knnAff.k
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "KNNAff", 
             knnAff.methodAggregation = tmp.knnAff.methodAggregation, 
             knnAff.range = tmp.knnAff.range, 
             knnAff.k = tmp.knnAff.k, 
             knnAff.classifiers = NULL), 
        .Names = c("name", "method", "knnAff.methodAggregation", "knnAff.range", "knnAff.k", "knnAff.classifiers"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "KNN"){
      tmp.name <- settings[[i]]$name
      if (!is.numeric(settings[[i]]$knn.range) 
          || is.na(settings[[i]]$knn.range) 
          || length(settings[[i]]$knn.range) != 1 
          || !.is.wholenumber(settings[[i]]$knn.range) 
          || !(settings[[i]]$knn.range >= 2 
               && settings[[i]]$knn.range <= (ddalpha$numPoints - 1) 
               || settings[[i]]$knn.range == -1)){
        warning("In treatment number ", i, ": Argument \"knn.range\" not specified correctly. Defaults are applied")
        tmp.knn.range <- -1
      }else{
        tmp.knn.range <- settings[[i]]$knn.range
      }
      if (!is.numeric(settings[[i]]$knn.k) 
          || is.na(settings[[i]]$knn.k) 
          || length(settings[[i]]$knn.k) != 1 
          || !.is.wholenumber(settings[[i]]$knn.k) 
          || !(settings[[i]]$knn.k >= 1 
               && settings[[i]]$knn.k <= (ddalpha$numPoints) 
               || settings[[i]]$knn.k == -1)){
        warning("In treatment number ", i, ": Argument \"knn.k\" not specified correctly. Defaults are applied")
        tmp.knn.k <- -1
      }else{
        tmp.knn.k <- settings[[i]]$knn.k
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "KNN", 
             knn.range = tmp.knn.range, 
             knn.k = tmp.knn.k, 
             knn.train = NULL, 
             knn.cl = NULL), 
        .Names = c("name", "method", "knn.range", "knn.k", "knn.train", "knn.cl"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "DEPTH.MAHALANOBIS"){
      tmp.name <- settings[[i]]$name
      if (!is.character(settings[[i]]$mah.estimate) 
          || length(settings[[i]]$mah.estimate) != 1 
          || !(settings[[i]]$mah.estimate %in% c("moment", "MCD"))){
        warning("In treatment number ", i, ": Argument \"mah.estimate\" not specified correctly. \"moment\" is used as a default value")
        tmp.mah.estimate <- "moment"
      }else{
        tmp.mah.estimate <- settings[[i]]$mah.estimate
      }
      if (!is.vector(settings[[i]]$priors, mode = "double") 
          || is.na(min(settings[[i]]$priors)) 
          || length(settings[[i]]$priors) != ddalpha$numPatterns 
          || min(settings[[i]]$priors) <= 0 
          || max(settings[[i]]$priors) <= 0){
        warning("In treatment number ", i, ": Argument \"priors\" not specified correctly. Defaults in the form of class portions are applied")
        tmp.priors <- NULL
      }else{
        tmp.priors <- settings[[i]]$priors/sum(settings[[i]]$priors)
      }
      if (!is.vector(settings[[i]]$mcd.alpha, mode = "double") 
          || is.na(min(settings[[i]]$mcd.alpha)) 
          || length(settings[[i]]$mcd.alpha) != 1 
          || settings[[i]]$mcd.alpha < 0.5 
          || settings[[i]]$mcd.alpha > 1){
        if (tmp.mah.estimate == "MCD"){
          warning("In treatment number ", i, ": Argument \"mcd.alpha\" not specified correctly. 0.75 is used as a default value")
        }
        tmp.mcd.alpha <- 0.75
      }else{
        tmp.mcd.alpha <- settings[[i]]$mcd.alpha
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "depth.Mahalanobis", 
             mah.estimate = tmp.mah.estimate, 
             priors = tmp.priors, 
             mah.classes = NULL, 
             mcd.alpha = tmp.mcd.alpha), 
        .Names = c("name", "method", "mah.estimate", "priors", "mah.classes", "mcd.alpha"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "RANDEQUAL"){
      tmp.name <- settings[[i]]$name
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "RandEqual"), 
        .Names = c("name", "method"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "RANDPROP"){
      tmp.name <- settings[[i]]$name
      if (!is.vector(settings[[i]]$priors, mode = "double") 
          || is.na(min(settings[[i]]$priors)) 
          || length(settings[[i]]$priors) != ddalpha$numPatterns 
          || min(settings[[i]]$priors) <= 0 
          || max(settings[[i]]$priors) <= 0){
        warning("In treatment number ", i, ": Argument \"priors\" not specified correctly. Defaults in the form of class portions are applied")
        tmp.priors <- NULL
      }else{
        tmp.priors <- settings[[i]]$priors/sum(settings[[i]]$priors)
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "RandProp", 
             priors = NULL), 
        .Names = c("name", "method", "priors"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "IGNORE"){
      tmp.name <- settings[[i]]$name
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "Ignore"), 
        .Names = c("name", "method"))
      supported <- TRUE
    }
    if (supported){
      treatments.settings[[counter]] <- treatment.settings
      counter <- counter + 1
    }
  }

  return (treatments.settings)
}

.lda_learn <- function(ddalpha, settings){
  new.frm <- data.frame(ddalpha$raw)
  settings$lda <- MASS::lda(formula=as.formula(paste("X", ncol(ddalpha$raw), " ~ .", sep="")), data=data.frame(ddalpha$raw), priors=settings$priors)
  settings$priors <- settings$lda$prior
  return (settings)
}

.lda_classify <- function(objects, ddalpha, settings){
  return (as.list(predict(settings$lda, data.frame(objects))$class))
}

.knnAff_learn <- function(ddalpha, settings){
  counter <- 1
  # Determining multi-class behaviour
  if (settings$knnAff.methodAggregation == "majority"){
    for (i in 1:(ddalpha$numPatterns - 1)){
      for (j in (i + 1):ddalpha$numPatterns){
        # Creating a classifier
        classifier.index          <- counter
        classifier.index0         <- i
        classifier.index1         <- j
        classifier.points         <- as.double(t(rbind(ddalpha$patterns[[i]]$points, ddalpha$patterns[[j]]$points)))
        classifier.cardinalities  <- as.integer(c(ddalpha$patterns[[i]]$cardinality, ddalpha$patterns[[j]]$cardinality))
        if (settings$knnAff.k < 1 || settings$knnAff.k > (ddalpha$patterns[[i]]$cardinality + ddalpha$patterns[[j]]$cardinality - 1))
        {
          if (settings$knnAff.range < 2 || settings$knnAff.range > (ddalpha$patterns[[i]]$cardinality + ddalpha$patterns[[j]]$cardinality - 1)){
            maxk <- 10*( (ddalpha$numPoints)^(1/ddalpha$dimension) ) + 1
          }else{
            maxk <- settings$knnAff.range
          }
          maxk <- min(maxk, ddalpha$patterns[[i]]$cardinality + ddalpha$patterns[[j]]$cardinality - 1)
          maxk <- max(maxk, 2)
          classifier.range <- maxk
          classifier.k <- as.integer(.C("KnnAffInvLearnJK", 
                             classifier.points, 
                             as.integer(ddalpha$dimension), 
                             classifier.cardinalities, 
                             as.integer(maxk), 
                             k=integer(1))$k)
        }else{
          classifier.range <- settings$knnAff.range
          classifier.k <- as.integer(settings$knnAff.k)
        }
        # Adding the classifier to the list of classifiers
        settings$knnAff.classifiers[[counter]] <- structure(
          list(index = classifier.index, 
               index0 = classifier.index0, 
               index1 = classifier.index1, 
               points = classifier.points, 
               cardinalities = classifier.cardinalities, 
               k = classifier.k, 
               range = classifier.range), 
          .Names = c("index", "index0", "index1", "points", "cardinalities", "k", "range"))
        counter <- counter + 1
      }
    }
  }
  if (settings$knnAff.methodAggregation == "sequent"){
    for (i in 1:ddalpha$numPatterns){
      anotherClass <- NULL
      for (j in 1:ddalpha$numPatterns){
        if (j != i){
          anotherClass <- rbind(anotherClass, ddalpha$patterns[[j]]$points)
        }
      }
      classifier.index          <- counter
      classifier.index0         <- i
      classifier.index1         <- -1
      classifier.points         <- as.double(t(rbind(ddalpha$patterns[[i]]$points, anotherClass)))
      classifier.cardinalities  <- as.integer(c(ddalpha$patterns[[i]]$cardinality, nrow(anotherClass)))
      if (settings$knnAff.k < 1 || settings$knnAff.k > ddalpha$numPoints)
      {
        if (settings$knnAff.range < 2 || settings$knnAff.range > (ddalpha$numPoints - 1)){
          maxk <- 10*( (ddalpha$numPoints)^(1/ddalpha$dimension) ) + 1
        }else{
          maxk <- settings$knnAff.range
        }
        maxk <- min(maxk, ddalpha$numPoints - 1)
        maxk <- max(maxk, 2)
        classifier.range <- maxk
        classifier.k <- as.integer(.C("KnnAffInvLearnJK", 
                           classifier.points, 
                           as.integer(ddalpha$dimension), 
                           classifier.cardinalities, 
                           as.integer(maxk), 
                           k=integer(1))$k)
      }else{
        classifier.range <- settings$knnAff.range
        classifier.k <- as.integer(settings$knnAff.k)
      }
      # Adding the classifier to the list of classifiers
      settings$knnAff.classifiers[[counter]] <- structure(
        list(index = classifier.index, 
             index0 = classifier.index0, 
             index1 = classifier.index1, 
             points = classifier.points, 
             cardinalities = classifier.cardinalities, 
             k = classifier.k, 
             range = classifier.range), 
        .Names = c("index", "index0", "index1", "points", "cardinalities", "k", "range"))
      counter <- counter + 1
    }
  }
  return (settings)
}

.knnAff_classify <- function(objects, ddalpha, settings){
  # Correct input data
  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  # Initialization of the vote array
  votes <- matrix(rep(0, nrow(objects)*ddalpha$numPatterns), nrow=nrow(objects))
  for (i in 1:length(settings$knnAff.classifiers)){
    res <- .C("KnnAffInvClassify", 
              as.double(t(objects)), 
              as.integer(nrow(objects)), 
              settings$knnAff.classifiers[[i]]$points, 
              as.integer(ddalpha$dimension), 
              settings$knnAff.classifiers[[i]]$cardinalities, 
              settings$knnAff.classifiers[[i]]$k, 
              output=integer(nrow(objects)))$output
    for (j in 1:nrow(objects)){
      if (res[j] == 0){
        votes[j,settings$knnAff.classifiers[[i]]$index0] <- votes[j,settings$knnAff.classifiers[[i]]$index0] + 1
      }else{
        votes[j,settings$knnAff.classifiers[[i]]$index1] <- votes[j,settings$knnAff.classifiers[[i]]$index1] + 1
      }
    }
  }
  # Collect results
  results <- list()
  for (i in 1:nrow(objects)){
    results[[i]] <- ddalpha$patterns[[which.max(votes[i,])]]$name
  }
  return (results)
}

.knn_learn <- function(ddalpha, settings){
  settings$knn.train <- ddalpha$raw[,1:ddalpha$dimension]
  settings$knn.cl <- ddalpha$raw[,ddalpha$dimension + 1]
  if (settings$knn.k < 1 || settings$knn.k > ddalpha$numPoints){
    if (settings$knn.range < 1 || settings$knn.range > ddalpha$numPoints - 1){
      settings$knn.range <- 10*( (ddalpha$numPoints)^(1/ddalpha$dimension) ) + 1
      settings$knn.range <- min(settings$knn.range, ddalpha$numPoints - 1)
      settings$knn.range <- max(settings$knn.range, 2)
    }
    cv.err <- c()
    ks <- 1:settings$knn.range
    for (i in ks){
      newpre <- as.vector(class::knn.cv(settings$knn.train, settings$knn.cl, k = i))
      cv.err <- c(cv.err, sum(settings$knn.cl != newpre))
    }
    settings$knn.k <- ks[which.min(cv.err)] 
  }
  return (settings)
}

.knn_classify <- function(objects, ddalpha, settings){
  return (class::knn(settings$knn.train, objects, settings$knn.cl, settings$knn.k))
}

.mah_learn <- function(ddalpha, settings){
  settings$mah.classes <- list()
  if (is.null(settings$priors)){
    settings$priors <- c()
    for (i in 1:ddalpha$numPatterns){
      settings$priors[i] <- ddalpha$patterns[[i]]$cardinality/ddalpha$numPoints
    }
  }
  for (i in 1:ddalpha$numPatterns){
    class.prior <- NULL
    class.mean <- NULL
    class.cov <- NULL
    class.sigma <- NULL
    class.prior <- settings$priors[i]
    if (settings$mah.estimate == "moment"){
      class.mean <- colMeans(ddalpha$patterns[[i]]$points)
      class.cov <- cov(ddalpha$patterns[[i]]$points)
      class.sigma <- solve(class.cov)
    }
    if (settings$mah.estimate == "MCD"){
      estimate <- robustbase::covMcd(ddalpha$patterns[[i]]$points, alpha=settings$mcd.alpha)
      class.mean <- estimate$center
      class.cov <- estimate$cov
      class.sigma <- solve(class.cov)
    }
    settings$mah.classes[[i]] <- structure(
      list(index = i, 
           prior = class.prior, 
           mean = class.mean, 
           cov = class.cov, 
           sigma = class.sigma), 
      .Names = c("index", "prior", "mean", "cov", "sigma"))
  }
  return (settings)
}

.mah_classify <- function(objects, ddalpha, settings){
  # Correct input data
  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  # Initialization of the vote array
  votes <- matrix(rep(0, nrow(objects)*ddalpha$numPatterns), nrow=nrow(objects))
  for (i in 1:nrow(objects)){
    for (j in 1:length(settings$mah.classes)){
      votes[i,j] <- settings$mah.classes[[j]]$prior*.mahalanobis_depth(objects[i,], settings$mah.classes[[j]]$mean, settings$mah.classes[[j]]$sigma)
    }
  }
  # Collect results
  results <- list()
  for (i in 1:nrow(objects)){
    results[[i]] <- ddalpha$patterns[[which.max(votes[i,])]]$name
  }
  return (results)
}

.ignore_classify <- function(nobjects){
  return (as.list(rep("Ignored", nobjects)))
}

.right_classify <- function(nobjects){
  return (as.list(rep("Outsider", nobjects)))
}

.randequal_classify <- function(nobjects, ddalpha){
  results <- list()
  for (i in 1:nobjects){
    results[[i]] <- ddalpha$patterns[[sample(1:ddalpha$numPatterns,1)]]$name
  }
  return (results)
}

.randprop_classify <- function(nobjects, ddalpha, settings){
  priors <- settings$priors
  if (is.null(priors)){
    priors <- c()
    for (i in 1:ddalpha$numPatterns){
      priors[i] <- ddalpha$patterns[[i]]$cardinality/ddalpha$numPoints
    }
  }
  results <- list()
  for (i in 1:nobjects){
    results[[i]] <- ddalpha$patterns[[sample(1:ddalpha$numPatterns,1,prob = priors)]]$name
  }
  return (results)
}

# Function is taken from the R-documentation, "Examples" to the function "is.integer"
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

print.ddalpha <- function(x, ...){
  cat("ddalpha:\n")
  cat("\t num.points = ", x$numPoints, 
      ", dimension = ", x$dimension, 
      ", num.patterns = ", x$numPatterns, "\n", sep="")
  cat("\t depth = \"", x$methodDepth, "\"\n", sep="")
  cat("\t aggregation.method = \"", x$methodAggregation, "\"\n", sep="")
  cat("\t num.chunks =", x$numChunks, "\n")
  if (x$methodDepth == "randomTukey"){
    cat("\t num.directions =", x$numDirections, "\n")
  }
  cat("\t use.convex =", x$useConvex, "\n")
  cat("\t max.degree =", x$maxDegree, "\n")
  cat("patterns:\n")
  for (i in 1:length(x$patterns)){
    cat("\t ");print(x$patterns[[i]])
  }
  cat("num.classifiers =", x$numClassifiers, "\n")
  cat("outsider.methods:\n")
  for (i in 1:length(x$methodsOutsider)){
    cat ("\t \"",x$methodsOutsider[[i]]$name, "\":\n", sep="")
    cat("\t\t method = \"", x$methodsOutsider[[i]]$method, "\"\n", sep="")
    if (x$methodsOutsider[[i]]$method == "LDA" 
     || x$methodsOutsider[[i]]$method == "RandProp" 
     || x$methodsOutsider[[i]]$method == "depth.Mahalanobis"){
      cat("\t\t priors =", x$methodsOutsider[[i]]$priors, "\n")
    }
    if (x$methodsOutsider[[i]]$method == "KNNAff"){
      cat("\t\t aggregation.method = \"", 
          x$methodsOutsider[[i]]$knnAff.methodAggregation, "\"\n", sep="")
      for (j in 1:length(x$methodsOutsider[[i]]$knnAff.classifiers)){
        cat("\t\t k.range = ", format(paste("1:", 
          x$methodsOutsider[[i]]$knnAff.classifiers[[j]]$range, sep=""), 
          justify="right", width=10), sep="")
      }
      cat("\n")
      for (j in 1:length(x$methodsOutsider[[i]]$knnAff.classifiers)){
        cat("\t\t k       = ", format( 
          x$methodsOutsider[[i]]$knnAff.classifiers[[j]]$k, 
          justify="right", width=10), sep="")
      }
      cat("\n")
    }
    if (x$methodsOutsider[[i]]$method == "KNN"){
      cat("\t\t k.range = 1:", x$methodsOutsider[[i]]$knn.range, "\n", sep="")
      cat("\t\t k =", x$methodsOutsider[[i]]$knn.k, "\n")
    }
    if (x$methodsOutsider[[i]]$method == "depth.Mahalanobis"){
      cat("\t\t estimate = \"", x$methodsOutsider[[i]]$mah.estimate, "\"\n", sep="")
      if (x$methodsOutsider[[i]]$mah.estimate == "MCD"){
        cat("\t\t mcd.alpha = ", x$methodsOutsider[[i]]$mcd.alpha, "\n")
      }
    }
  }
}

print.ddalpha.pattern <- function(x, ...){
  cat("pattern[", x$index, "]:", sep="")
  cat("\t ", x$cardinality, " points, label = \"", x$name, "\"\n", sep="")
}
