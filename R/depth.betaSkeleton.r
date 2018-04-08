################################################################################
# File:             depth.betaSkeleton.r
# Created by:       Pavlo Mozharovskyi
# First published:  08.03.2018
# Last revised:     08.03.2018
# 
# Computation of the beta-skeleton depth.
################################################################################

depth.betaSkeleton <- function(x, data, beta = 2, distance = "Lp", Lp.p = 2, 
                               mah.estimate = "moment", mah.parMcd = 0.75){
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
  if (nrow(data) < 2){
    stop("Too few data points")
  }
  if (!is.numeric(beta) || length(beta) > 1 || beta < 1){
    stop("Argument \"beta\" should be a real value >= 1")
  }
  if (!is.numeric(Lp.p) || length(Lp.p) > 1 || Lp.p < 0){
    stop("Argument \"Lp.p\" should be a real value >= 0")
  }
  if (toupper(distance) == "MAHALANOBIS"){
    distance.code <- 5L
    if(toupper(mah.estimate) == "NONE"){
      sigma = diag(ncol(data))
    } else {
      if(toupper(mah.estimate) == "MOMENT"){
        cov <- cov(data)
      } else if(toupper(mah.estimate) == "MCD"){
        cov <- covMcd(data, mah.parMcd)$cov
      } else {stop("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")}
      if(sum(is.na(cov)) == 0){
        sigma <- solve(cov)
      } else{
        sigma = diag(ncol(data))
        warning("Covariance estimate not found, no affine-invariance-adjustment")
      }
    }
  }else{
    sigma = 0;
    if (toupper(distance) == "LP"){
      distance.code <- 4L
      if (Lp.p == 1){
        distance.code <- 1L
      }
      if (Lp.p == 2){
        distance.code <- 2L
      }
      if (is.infinite(Lp.p) && Lp.p > 0){
        distance.code <- 3L
      }
    } else {stop("Argument \"distance\" should be either \"Lp\" or \"Mahalanobis\"")}
  }
  
  points <- as.vector(t(data))
  objects <- as.vector(t(x))
  ds <- .C("BetaSkeletonDepth", 
           as.double(points), 
           as.double(objects), 
           as.integer(nrow(data)), 
           as.integer(nrow(x)), 
           as.integer(ncol(data)), 
           as.double(beta), 
           as.integer(distance.code), 
           as.double(Lp.p), 
           as.double(sigma), 
           depths=double(nrow(x)))$depths
  
  return (ds)
}
