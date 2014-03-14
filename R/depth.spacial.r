depth.spatial <- function(x, data){
  if (!is.matrix(x) 
      && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  mean <- colMeans(data)
  cov <- cov(data)
  cov.eig <- eigen(cov)
  B <- cov.eig$vectors %*% diag(sqrt(cov.eig$values))
  lambda <- solve(B)
  depths <- rep(-1, nrow(x))
  for (i in 1:nrow(x)){
    tmp1 <- t(lambda %*% (x[i,] - t(data)))
    tmp1 <- tmp1[which(rowSums(tmp1) != 0),]
    tmp2 <- 1/sqrt(rowSums(tmp1^2))
    depths[i] <- 1 - sqrt(sum((colSums(tmp2*tmp1)/nrow(data))^2))
  }
  return (depths)
}
