################################################################################
# File:             depth.L2.r
# Created by:       Pavlo Mozharovskyi
# First published:  08.03.2018
# Last revised:     08.03.2018
# 
# Computation of the L2-depth.
################################################################################

depth.L2 <- function(x, data, mah.estimate = "moment", mah.parMcd = 0.75){
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if(is.data.frame(data))
    data = data.matrix(data)
  if (!is.matrix(x)){
    if(is.vector(x))
      x <- matrix(x, nrow=1)
    if(is.data.frame(x))
      x = data.matrix(x)
  }
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

  depths <- rep(-1, nrow(x))
  for (i in 1:nrow(x)){
    tmp1 <- t(x[i,] - t(data))
    tmp2 <- tmp1 %*% sigma
    depths[i] <- 1/(1 + mean(sqrt(rowSums(tmp2 * tmp1))))
  }
  return (depths)
}
