################################################################################
# File:             depth.randomTukey.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     28.02.2013
# 
# Computation of the random Tukey data depth.
################################################################################

depth.randomTukey <- function(x, data, num.directions = 1000){
  if (!is.matrix(x) 
      && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.numeric(x)){
    stop("Argument \"x\" should be numeric")
  }
  if (ncol(x) != ncol(data)){
    stop("Dimensions of the arguments \"x\" and \"data\" should coincide")
  }
  if (ncol(data) + 1 > nrow(data)){ #?
    stop("To few data points")
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

  points <- as.vector(t(data))
  objects <- as.vector(t(x))
  c <- as.vector(nrow(data))
  k <- numDirections
  ds <- .C("HDepth", 
           as.double(points), 
           as.double(objects), 
           as.integer(nrow(x)), 
           as.integer(ncol(data)), 
           as.integer(c), 
           as.integer(1), 
           as.double(0), 
           as.double(0), 
           as.integer(k), 
           as.integer(1), # use the same directions and projections
           depths=double(nrow(x)))$depths
  
  return (ds)
}
