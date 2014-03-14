compclassf.train <- function(dataf, labels, 
                            to.equalize = TRUE, 
                            to.reduce = FALSE, 
                            classifier.type = "ddalpha", 
                            ...){
  # Trains the functional componentwise classifier
  # Args:
  #   dataf:  list containing lists (functions) of two vectors of equal length, 
  #           named "args" and "vals": arguments sorted in ascending order and 
  #           corresponding them values respectively
  #   labels: output labels of the functinal observations
  #   other arguments: TODO
  # Returns:
  #   Functional componentwise clasifier
  
  # Check "dataf"
  # TODO
  
  # Check "labels"
  # TODO
  
  # Check classifier.type
  
  # Bring to finite dimension
  
  # Pointize
  points <- GetPointsDHB12(dataf, labels, to.equalize, to.reduce)
  # CV
  arg.indices <- getBestSpaceDHB12(points$data, classifier.type, num.chunks=10, ...) #, ...
  data <- points$data[,c(arg.indices,ncol(points$data))]
  # Apply chosen classifier to train the data
  if (classifier.type == "ddalpha"){
    classifier <- ddalpha.train(data, ...)
  }
  if (classifier.type == "maxdepth"){
    classifier <- maxdepth.train(data, ...)
  }
  if (classifier.type == "knnaff"){
    classifier <- knnaff.train(data, i = 0, ...)
  }
  if (classifier.type == "lda"){
    classifier <- lda.train(data, ...)
  }
  if (classifier.type == "qda"){
    classifier <- qda.train(data, ...)
  }
  # Create the eventual output structure
  compclassf <- structure(
    list(dataf = points$dataf, 
      labels = points$labels, 
      adc.method = "equalCover", 
      adc.args = list(instance = "val", numFcn = ncol(points$data) - 1, numDer = 0), 
      adc.transmat = points$transmat, 
      the.args = arg.indices, 
      data = points$data, 
      classifier.type = classifier.type, 
      classifier = classifier), 
    .Names = c("dataf", "labels", "adc.method", "adc.args", "adc.transmat", 
               "the.args",  "data", "classifier.type", "classifier"))
  class(compclassf) <- "compclassf"
  
  return (compclassf)
}

compclassf.classify <- function(objectsf, compclassf, ...){
  # Classifies functions
  # Args:
  #   objectsf: sample to classify, a list containing lists (functions) of 
  #             two vectors of equal length, named "args" and "vals": 
  #             arguments sorted in ascending order and corresponding them 
  #             values respectively
  #   compclassf: functional DDalpha-classifier
  # Returns:
  #   List of labels assigned to the functions from "objectsf"

  # Check "objectsf"
  # TODO

  # Prepare to multivariate classification
  objectsf.equalized <- equalize(objectsf)
  if (compclassf$adc.method == "equalCover"){
    if (compclassf$adc.args$instance == "val"){
      input <- getValGrid(objectsf.equalized, 
                          compclassf$adc.args$numFcn, compclassf$adc.args$numDer)
    }
    if (compclassf$adc.args$instance == "avr"){
      input <- getAvrGrid(objectsf.equalized, 
                          compclassf$adc.args$numFcn, compclassf$adc.args$numDer)
    }
    if (!is.null(compclassf$adc.transmat)){
      input <- input%*%compclassf$adc.transmat
    }
  }
  input <- input[,compclassf$the.args]
  # Classify and assign class labels
  if (compclassf$classifier.type == "ddalpha"){
    output <- ddalpha.classify(objects = input, compclassf$classifier, ...)
  }
  if (compclassf$classifier.type == "maxdepth"){
    output <- maxdepth.classify(objects = input, compclassf$classifier, ...)
  }
  if (compclassf$classifier.type == "knnaff"){
    output <- knnaff.classify(objects = input, compclassf$classifier, ...)
  }
  if (compclassf$classifier.type == "lda"){
    output <- lda.classify(objects = input, compclassf$classifier, ...)
  }
  if (compclassf$classifier.type == "qda"){
    output <- qda.classify(objects = input, compclassf$classifier, ...)
  }
  classes <- list()
  for (i in 1:length(output)){
#    if (is.numeric(output[[i]])){
      classes[[i]] <- compclassf$labels[[ output[[i]] ]]
#    }else{
#      classes[[i]] <- output[[i]]
#    }
  }

  return (classes)
}

is.in.convexf <- function(objectsf, dataf, cardinalities, 
                          adc.method = "equalCover", 
                          adc.args = list(instance = "val", 
                                          numFcn = 5, 
                                          numDer = 5)){
  # Checks if the function(s) lie(s) inside convex hulls of the 
  # functions from the sample in the projection space
  # Args:
  #   objectsf:      list containing lists (functions) of two vectors of equal 
  #                  length, named "args" and "vals": arguments sorted in 
  #                  ascending order and corresponding them values 
  #                  respectively. These functions are supposed to be checked 
  #                  for 'outsiderness'
  #   dataf:         list containing lists (functions) of two vectors of equal 
  #                  length, named "args" and "vals": arguments sorted in 
  #                  ascending order and corresponding them values respectively
  #   cardinalities: cardinalities of the classes in "dataf"
  #   other arguments: TODO
  # Returns:
  #   Vector with 1s for those lying inside of at least one of the convex hulls 
  #   of the classes and 0s for those lying beyond them
  
  # Data-consistency checks
  # TODO
  
  # Project "objectsf" into a multivariate space
  objectsf.equalized <- equalize(objectsf)
  if (adc.method == "equalCover"){
    if (adc.args$instance == "val"){
      objects <- getValGrid(objectsf.equalized, adc.args$numFcn, adc.args$numDer)
    }
    if (adc.args$instance == "avr"){
      objects <- getAvrGrid(objectsf.equalized, adc.args$numFcn, adc.args$numDer)
    }
  }
  # Project "dataf" into a multivariate space
  dataf.equalized <- equalize(dataf)
  if (adc.method == "equalCover"){
    if (adc.args$instance == "val"){
      data <- getValGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
    }
    if (adc.args$instance == "avr"){
      data <- getAvrGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
    }    
  }
  in.convex <- is.in.convex(objects, data, cardinalities)
  
  return (in.convex)
}

################################################################################
# Functions below are used for intermediate computations                       #
################################################################################

getVapnikBound <- function(points, dim = NULL){
  n <- nrow(points)
  d <- ncol(points) - 1
  lda <- lda.train(points)
  result <- lda.classify(points[,1:d], lda)
  empRisk <- sum(result != points[,d + 1])/n
  
  if (is.null(dim)){
    dim <- d
  }
  nu <- 1/n
  if (empRisk == 0){
    epsilon <- 1 - n^(-(dim + 1)/n)
  }else{
    epsilon <- sqrt( ( (dim + 1)*log(n) )/(2*n) )
  }
  
  return (empRisk + epsilon)
}

GetPointsDHB12 <- function(dataf, labels, to.equalize=T, to.reduce=F){
  # Numerize labels
  names <- unique(labels)
  output <- rep(0, length(labels))
  for (i in 1:length(labels)){
    for (j in 1:length(names)){
      if (labels[[i]] == names[[j]]){
        output[i] = j
        break
      }
    }
  }
  # Prepare data
  if (to.equalize){
    num.times = length(dataf[[1]]$args)
    dataf.equalized <- equalize(dataf)
    adc.args = list(instance = "val", 
                    numFcn = num.times, 
                    numDer = 0)
    input <- getValGrid(dataf.equalized, adc.args$numFcn, adc.args$numDer)
  }else{
    input <- NULL
    for (i in 1:length(dataf)){
      input <- rbind(input, dataf[[i]]$vals)
    }
  }
  transmat <- NULL
  if (to.reduce){# Reduce dimension if needed
    princomps <- NULL
    newDim <- ncol(input)
    for (i in 1:length(names)){
      classi <- input[output == i,1:ncol(input)]
      princompsi <- prcomp(x=classi, tol=sqrt(.Machine$double.eps))
      newDimi <- sum(princompsi$sdev > sqrt(.Machine$double.eps))
      if (newDimi < newDim){
        newDim <- newDimi
        princomps <- princompsi
      }
    }
    transmat <- NULL
    if (newDim < ncol(input)){
      transmat <- matrix(as.vector(princomps$rotation[,1:newDim]), ncol=newDim)
      input <- input%*%transmat
    }
  }
  # Combine data
  data <- cbind(input, output, deparse.level=0)
  return (list(data = data, dataf = dataf.equalized, labels = names, transmat = transmat))
}

getBestSpaceDHB12 <- function(data, 
                              classifier.type = "ddalpha", 
                              num.chunks = 10, 
                              ...){
  indices.num <- ncol(data) - 1
  indices.avlbl <- rep(TRUE, indices.num)
  indices.best <- c()
  error.last <- nrow(data) + 1
  r <- 0
  while (sum(indices.avlbl) > 0){
    # If this is the first iteration search through all possible pairs
    if (r == 0){
      # Generate all combinations with smallest distance 2
      combinations <- combn((1:indices.num)[indices.avlbl], 2)
      tmp.cmb <- rbind(combinations[-1,], rep(-1000000, ncol(combinations)))
      tmp.cmb <- (tmp.cmb - combinations)==T
      combinations <- combinations[,apply(tmp.cmb, 2, sum)==0]
      # Choose the best combination
      errors <- c()
      for (i in 1:ncol(combinations)){
        cat("r = ", r, ": ", i, "/", ncol(combinations), ".\n", sep="")
        errors <- c(errors, 0)
        # Actually CV
        num.points <- nrow(data)
        indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
        for (j in 1:num.chunks){
          # Determine points to be taken off
          take.off <- (indices.off + j)[(indices.off + j) <= num.points]
          # Apply chosen classifier
          if (classifier.type == "ddalpha"){
            classifier <- ddalpha.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- ddalpha.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "maxdepth"){
            classifier <- maxdepth.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- maxdepth.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "knnaff"){
            classifier <- knnaff.train(data[-take.off,c(combinations[,i], indices.num + 1)], i = i, ...)
            results <- knnaff.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "lda"){
            classifier <- lda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- lda.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "qda"){
            classifier <- qda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- qda.classify(data[take.off,combinations[,i]], classifier)
          }
          # Collect errors
          errors[i] <- errors[i] + sum(unlist(results) != data[take.off,indices.num + 1])
        }
      }
      # Collect results
      error.last <- min(errors)
      indices.best <- combinations[,which.min(errors)]
      indices.avlbl <- rep(TRUE, indices.num)
      indices.to.dsbl <- unique(c(indices.best, indices.best - 1, indices.best + 1))
      indices.to.dsbl <- indices.to.dsbl[indices.to.dsbl >= 1 && indices.to.dsbl <= indices.num]
      indices.avlbl[indices.to.dsbl] <- FALSE
      r <- 2
      next
    }
    # First, sequential approach
    errors <- c()
    variants <- c()
    for (i in 1:indices.num){
      if (indices.avlbl[i]){
        errors <- c(errors, 0)
        variants <- c(variants, i)
        indices.cur <- c(indices.best, i)
        # Actually CV
        num.points <- nrow(data)
        indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
        for (j in 1:num.chunks){
          # Determine points to be taken off
          take.off <- (indices.off + j)[(indices.off + j) <= num.points]
          # Apply chosen classifier
          if (classifier.type == "ddalpha"){
            classifier <- ddalpha.train(data[-take.off,c(indices.cur, indices.num + 1)], ...)
            results <- ddalpha.classify(data[take.off,indices.cur], classifier)
          }
          if (classifier.type == "maxdepth"){
            classifier <- maxdepth.train(data[-take.off,c(indices.cur, indices.num + 1)], ...)
            results <- maxdepth.classify(data[take.off,indices.cur], classifier)
          }
          if (classifier.type == "knnaff"){
            classifier <- knnaff.train(data[-take.off,c(indices.cur, indices.num + 1)], i = i, ...)
            results <- knnaff.classify(data[take.off,indices.cur], classifier)
          }
          if (classifier.type == "lda"){
            classifier <- lda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- lda.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "qda"){
            classifier <- qda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- qda.classify(data[take.off,combinations[,i]], classifier)
          }
          # Collect errors
          errors[i] <- errors[i] + sum(unlist(results) != data[take.off,indices.num + 1])
        }
      }
    }
    error.best <- min(errors)
    best.i <- variants[which.min(errors)[1]]
    indices.new <- c(indices.best, best.i)
    # Refinements for r=2, 3 and 4
    if (r %in% 2:3){
      # Define the grid
      if (r == 2){step <- 10}
      if (r == 3){step <- 5}
      grid.one <- c(-(1:step*2), 0, 1:step*2)
      grid <- c()
      for (i in 1:length(indices.new)){
        grid <- c(grid, indices.new[i] + grid.one)
      }
      grid <- unique(grid)
      grid <- sort(grid[(grid >= 1) & (grid <= indices.num)])
      # Generate all combinations with smallest distance 2
      combinations <- combn(grid, r + 1)
      tmp.cmb <- rbind(combinations[-1,], rep(-1000000, ncol(combinations)))
      tmp.cmb <- (tmp.cmb - combinations)==T
      combinations <- combinations[,apply(tmp.cmb, 2, sum)==0]
      # Choose the best combination
      #indices.grid <- (1:indices.num)[indices.avlbl & ((1:induces.num) %in% grid)]
      # Go through the combinations
      errors <- c()
      #combinations <- combn(indices.grid, r + 1)
      for (i in 1:ncol(combinations)){
        cat("r = ", r, ": ", i, "/", ncol(combinations), ".\n", sep="")
        errors <- c(errors, 0)
        # Actually CV
        num.points <- nrow(data)
        indices.off <- num.chunks*(0:(ceiling(num.points/num.chunks) - 1))
        for (j in 1:num.chunks){
          # Determine points to be taken off
          take.off <- (indices.off + j)[(indices.off + j) <= num.points]
          # Apply chosen classifier
          if (classifier.type == "ddalpha"){
            classifier <- ddalpha.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- ddalpha.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "maxdepth"){
            classifier <- maxdepth.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- maxdepth.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "knnaff"){
            classifier <- knnaff.train(data[-take.off,c(combinations[,i], indices.num + 1)], i = i, ...)
            results <- knnaff.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "lda"){
            classifier <- lda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- lda.classify(data[take.off,combinations[,i]], classifier)
          }
          if (classifier.type == "qda"){
            classifier <- qda.train(data[-take.off,c(combinations[,i], indices.num + 1)], ...)
            results <- qda.classify(data[take.off,combinations[,i]], classifier)
          }
          # Collect errors
          errors[i] <- errors[i] + sum(unlist(results) != data[take.off,indices.num + 1])
        }
      }
      error.best <- min(errors)
      indices.cur <- combinations[,which.min(errors)]
    }else{
      indices.cur <- indices.new
    }
    if (error.best < error.last){
      indices.best <- indices.cur
      error.last <- error.best
      indices.avlbl <- rep(TRUE, indices.num)
      indices.to.dsbl <- unique(c(indices.best, indices.best - 1, indices.best + 1))
      indices.to.dsbl <- indices.to.dsbl[indices.to.dsbl >= 1 && indices.to.dsbl <= indices.num]
      indices.avlbl[indices.to.dsbl] <- FALSE
      r <- r + 1
    }else{
      break
    }
  }
  return (indices.best)
}

# equalize <- function(dataf){
#   # 1. Adjusts the data to have equal (the largest) argument interval
#   # 2. Calclates - numerically - derivative
#   # Args:
#   #   dataf: list containing lists (functions) of two vectors 
#   #          of equal length, named "args" and "vals": arguments 
#   #          sorted in ascending order and corresponding them 
#   #          values respectively
#   # Returns:
#   #   The list of lists of the same structure, 'equalized', 
#   #   contating derivatives as 3rd vector named "der1"
#   
#   # Check whether every function contains fields "args" and "vals", 
#   # whether they are numerical and of equal length, have no NAs or 
#   # ties and are sorted in ascending order
#   # TODO
#   
#   # 1.
#   # Get argument bounds
#   min <- Inf
#   max <- -Inf
#   for (i in 1:length(dataf)){
#     if (dataf[[i]]$args[1] < min){min <- dataf[[i]]$args[1]}
#     if (dataf[[i]]$args[length(dataf[[i]]$args)] > max){
#       max <- dataf[[i]]$args[length(dataf[[i]]$args)]
#     }
#   }
#   # and apply them to equalize functions timely
#   for (i in 1:length(dataf)){
#     if (dataf[[i]]$args[1] > min){
#       dataf[[i]]$args <- c(min, dataf[[i]]$args)
#       dataf[[i]]$vals <- c(dataf[[i]]$vals[1], dataf[[i]]$vals)
#     }
#     if (dataf[[i]]$args[length(dataf[[i]]$args)] < max){
#       dataf[[i]]$args <- c(dataf[[i]]$args, max)
#       dataf[[i]]$vals <- c(dataf[[i]]$vals, 
#                            dataf[[i]]$vals[length(dataf[[i]]$vals)])
#     }
#     # Computational trick - add "-1" to the "left"
#     dataf[[i]]$args <- c(min - 1, dataf[[i]]$args)
#     dataf[[i]]$vals <- c(dataf[[i]]$vals[1], dataf[[i]]$vals)
#   }
#   
#   # 2.
#   for (i in 1:length(dataf)){
#     dataf[[i]] <- derive(dataf[[i]])
#   }
#   
#   return (dataf)
# }
# 
# derive <- function(fcn){
#   # Adds 1st derivative to the function: a vector named "der1"
#   # Args:
#   #   fcn: function, a list of two vectors of equal length, named "args" 
#   #   and "vals": arguments sorted in ascending order and corresponding 
#   #   them values respectively
#   # Returns:
#   #   The list of of the same structure contating derivative as 3rd 
#   #   vector named "der1"
#   
#   fcn$der1 <- rep(0, length(fcn$args))
#   fcn$der1[1] = -Inf
#   fcn$der1[2] = 0
#   for (i in 3:length(fcn$der1)){
#     fcn$der1[i] <- (fcn$vals[i] - fcn$vals[i - 1])/
#       (fcn$args[i] - fcn$args[i - 1])
#   }
#   
#   return (fcn)
# }
# 
# getValue <- function(fcn, arg, fcnInstance){
#   # Gets the value of the function or its derivative for the given argument 
#   # value
#   # Args:
#   #   fcn:         function, a list of vectors of equal length, named "args" 
#   #                (arguments), "vals" (function values) [and it's 
#   #                derivatives of order "I", named derI]; arguments (and 
#   #                corresponding values) sorted in ascending order
#   #   arg:         argument value at which the function (derivative) value 
#   #                is to be taken
#   #   fcnInstance: inctance to be evaluated; "vals" for the function values, 
#   #                "der1" for the first derivative
#   # Returns:
#   #   Value of the function (derivative)
#   
#   # Check "arg", evt. "fcnInstance"
#   # TODO
#   
#   # Find corresponding interval
#   index <- 2
#   while (arg > fcn$args[index]){
#     index <- index + 1
#   }
#   # Get the function value(s) by linear approximation
#   if (fcnInstance == "vals"){
#     value <- fcn$vals[index - 1] + 
#       (fcn$vals[index] - fcn$vals[index - 1])*
#       ((arg - fcn$args[index - 1])/(fcn$args[index] - fcn$args[index - 1]))
#   }
#   if (fcnInstance == "der1"){
#     value <- fcn$der1[index]
#   }
#   
#   return (value)
# }
# 
# getAvrValue <- function(fcn, argFrom, argTo, fcnInstance){
#   # Gets the average value of the function or its derivative on the given 
#   # interval
#   # Args:
#   #   fcn:         function, a list of vectors of equal length, named "args" 
#   #                (arguments), "vals" (function values) [and it's 
#   #                derivatives of order "I", named derI]; arguments (and 
#   #                corresponding values) sorted in ascending order
#   #   argFrom:     argument value from which the function (derivative) value 
#   #                is to be averaged
#   #   argTo:       argument value to which the function (derivative) value 
#   #                is to be averaged
#   #   fcnInstance: inctance to be evaluated; "vals" for the function values, 
#   #                "der1" for the first derivative
#   # Returns:
#   #   Average value of the function (derivative) on the interval 
#   #   (argFrom, argTo)
#   
#   # Check "argFrom" and "argTo", evt. "fcnInstance"
#   # TODO
#   
#   # Define 'from' and 'to' interval
#   indexFrom <- 2
#   while (argFrom > fcn$args[indexFrom]){
#     indexFrom <- indexFrom + 1
#   }
#   indexTo <- 2
#   while (argTo > fcn$args[indexTo]){
#     indexTo <- indexTo + 1
#   }
#   average <- 0
#   valTo <- getValue(fcn, argTo, fcnInstance)
#   # Integrate
#   curArgFrom <- argFrom
#   if (fcnInstance == "vals"){
#     valFrom <- getValue(fcn, curArgFrom, "vals")
#     while(indexFrom < indexTo){
#       average <- average + (valFrom + fcn$vals[indexFrom])*
#         (fcn$args[indexFrom] - curArgFrom)/2
#       valFrom <- fcn$vals[indexFrom]
#       curArgFrom <- fcn$args[indexFrom]
#       indexFrom <- indexFrom + 1
#     }
#     average <- average + (valFrom + valTo)*(argTo - curArgFrom)/2
#   }
#   if (fcnInstance == "der1"){
#     while(indexFrom < indexTo){
#       average <- average + (fcn$der1[indexFrom])*
#         (fcn$args[indexFrom] - curArgFrom)
#       curArgFrom <- fcn$args[indexFrom]
#       indexFrom <- indexFrom + 1
#     }
#     average <- average + valTo*(argTo - curArgFrom)
#   }
#   average <- average/(argTo - argFrom)
#   
#   return (average)
# }
# 
# getValGrid <- function(dataf, numFcn, numDer){
#   # Represents a function sample as a multidimensional (d="numFcn"+"numDer") 
#   # one averaging for that each function and it derivative on "numFcn" 
#   # (resp. "numDer") equal nonoverlapping covering intervals
#   # Args:
#   #   dataf:  list containing lists (functions) of vectors of equal length, 
#   #           first two named "args" and "vals" are arguments sorted in 
#   #           ascending order and having same bounds for all functions and 
#   #           corresponding them values respectively
#   #   numFcn: number of function intervals
#   #   numDer: number of first-derivative intervals
#   # Returns:
#   #   Matrix - a multidimensional presentation of the functional sample
#   
#   # Get argument bounds ("dataf" is equalized)
#   min <- dataf[[1]]$args[1]
#   max <- dataf[[1]]$args[length(dataf[[1]]$args)]
#   # Get argument grid
#   args <- dataf[[1]]$args
#   argsFcn <- min + 0:numFcn*(max - min)/(numFcn - 1)
#   argsDer <- min + 0:numDer*(max - min)/(numDer - 1)
#   # Get function/derivative grid
#   fcnGrid <- matrix(nrow = length(dataf), ncol = numFcn)
#   derGrid <- matrix(nrow = length(dataf), ncol = numDer)
#   if (numFcn > 0){
#     for (i in 1:length(dataf)){
#       # Set merging algorithm (Novikoff, Piter)
#       cArgs <- 1
#       cArgsFcn <- 1
#       fcnGrid[i,1] <- dataf[[i]]$vals[1]
#       while (cArgsFcn != numFcn){
#         #                print(argsFcn)
#         #                print(fcnGrid[i,])
#         #                cat(cArgs, " and ", cArgsFcn, "\n")
#         #                cat(args[cArgs + 1], " and ", argsFcn[cArgsFcn + 1], "\n")
#         if (args[cArgs + 1] < argsFcn[cArgsFcn + 1]){
#           cArgs <- cArgs + 1
#         }else{
#           nextArg <- argsFcn[cArgsFcn + 1]
#           fcnGrid[i,cArgsFcn + 1] <- dataf[[i]]$vals[cArgs] + (nextArg - args[cArgs])*dataf[[i]]$der1[cArgs + 1]
#           if (args[cArgs + 1] == argsFcn[cArgsFcn + 1]){
#             cArgs <- cArgs + 1
#           }
#           cArgsFcn <- cArgsFcn + 1
#         }
#       }
#     }
#   }
#   if (numDer > 0){
#     for (i in 1:length(dataf)){
#       # Again, set merging algorithm (Novikoff, Piter)
#       cArgs <- 1
#       cArgsDer <- 1
#       derGrid[1] <- dataf[[i]]$ders[2]
#       while (cArgsDer != numDer){
#         #        print(argsDer)
#         #        print(derGrid[i,])
#         #        cat(cArgs, " and ", cArgsDer, "\n")
#         #        cat(args[cArgs + 1], " and ", argsDer[cArgsDer + 1], "\n")
#         if (args[cArgs + 1] < argsDer[cArgsDer + 1]){
#           cArgs <- cArgs + 1
#         }else{
#           derGrid[i,cArgsDer + 1] <- dataf[[i]]$der1[cArgs + 1]
#           if (args[cArgs + 1] == argsDer[cArgsDer + 1]){
#             cArgs <- cArgs + 1
#           }
#           cArgsDer <- cArgsDer + 1
#         }
#       }
#     }
#   }
#   mvX <- cbind(fcnGrid, derGrid)
#   return (mvX)
# }
# 
# getAvrGrid <- function(dataf, numFcn, numDer){
#   # Represents a function sample as a multidimensional (d="numFcn"+"numDer") 
#   # one averaging for that each function and it derivative on "numFcn" 
#   # (resp. "numDer") equal nonoverlapping covering intervals
#   # Args:
#   #   dataf:  list containing lists (functions) of vectors of equal length, 
#   #           first two named "args" and "vals" are arguments sorted in 
#   #           ascending order and having same bounds for all functions and 
#   #           corresponding them values respectively
#   #   numFcn: number of function intervals
#   #   numDer: number of first-derivative intervals
#   # Returns:
#   #   Matrix - a multidimensional presentation of the functional sample
#   
#   # Get argument bounds ("dataf" is equalized)
#   min <- dataf[[1]]$args[1] + 1
#   max <- dataf[[1]]$args[length(dataf[[1]]$args)]
#   # Get argument grid
#   argsFcn <- min + 0:numFcn*(max - min)/numFcn
#   argsDer <- min + 0:numDer*(max - min)/numDer
#   # Get function/derivative grid
#   fcnGrid <- matrix(nrow = length(dataf), ncol = numFcn)
#   derGrid <- matrix(nrow = length(dataf), ncol = numDer)
#   if (numFcn > 0){
#     for (i in 1:length(dataf)){
#       for (j in 1:numFcn){
#         fcnGrid[i,j] <- getAvrValue(dataf[[i]], argsFcn[j], argsFcn[j + 1], 
#                                     "vals")
#       }
#     }
#   }
#   if (numDer > 0){
#     for (i in 1:length(dataf)){
#       for (j in 1:numDer){
#         derGrid[i,j] <- getAvrValue(dataf[[i]], argsDer[j], argsDer[j + 1], 
#                                     "der1")
#       }
#     }
#   }
#   mvX <- cbind(fcnGrid, derGrid)
#   
#   return (mvX)
# }
