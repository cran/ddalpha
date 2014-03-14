################################################################################
# File:             depth.space.randomTukey.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     28.02.2013
# 
# Computation of depth space based on the random Tukey data depth.
################################################################################

depth.space.randomTukey <- function(data, cardinalities, num.directions = 1000){
  if (!is.numeric(data)
      || !is.matrix(data) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.vector(cardinalities, mode = "numeric") 
      || is.na(min(cardinalities)) 
      || sum(.is.wholenumber(cardinalities)) != length(cardinalities) 
      || min(cardinalities) <= 0 
      || sum(cardinalities) != nrow(data)){
    stop("Argument \"cardinalities\" should be a vector of cardinalities of the classes in \"data\" ")
  }
  if (sum(cardinalities < ncol(data) + 1) != 0){
    stop("Not in all classes sufficiently enough objetcs")
  }
  if (!is.numeric(num.directions) 
      || is.na(num.directions) 
      || length(num.directions) != 1 
      || !.is.wholenumber(num.directions) 
      || !(num.directions > 1 && num.directions < 10000000)){
    numDirections <- 1000
    warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
  }else{
    numDirections <- num.directions
  }
  
  x <- as.vector(t(data))
  c <- as.vector(cardinalities)
  k <- numDirections
  rez <- .C("HDSpace", 
            as.double(x), 
            as.integer(ncol(data)), 
            as.integer(c), 
            as.integer(length(cardinalities)), 
            as.integer(k), as.integer(1), 
            dspc=double(nrow(data)*length(cardinalities)), 
            dirs=double(k*ncol(data)), 
            prjs=double(k*nrow(data)))
  depth.space <- matrix(rez$dspc, nrow=nrow(data), ncol=length(cardinalities), byrow=TRUE)
  
  return (depth.space)
}
