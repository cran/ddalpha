\name{depth.projection}
\alias{depth.projection}
\title{
  Calculate the Projection Depth
}
\description{
  Calculates the projection depth of points w.r.t. a multivariate data set.
}
\usage{
  depth.projection(x, data, method = "random", num.directions = 1000, seed = 0)
}
\arguments{
  \item{x}{
    Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
  }
  \item{data}{
    Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
  }
  \item{method}{
  to be used in calculations.
  
  \code{"random"} Here the depth is determined as the minimum univariate depth of the data projected on lines in several directions. The directions are distributed uniformly on the \eqn{(d-1)}-sphere; the same direction set is used for all points.
  
  \code{"linearize"} The Nelder-Mead method for function minimization, taken from Olsson, Journal of Quality Technology, 1974, 6, 56. R-codes of this function were written by Subhajit Dutta.
  }
  \item{num.directions}{
Number of random directions to be generated for \code{method = "random"}.
}
  \item{seed}{
the random seed. The dafault value \code{seed=0} makes no changes.
}
}
\details{
  Calculates projection depth. Projection depth, similar to Mahalanobis depth, is based on a measure of outlyingness,
used by Stahel (1981) and Donoho (1982), and has been first formulated by Liu (1992). The
worst case outlyingness is obtained by maximizing an outlyingness measure over all univariate
projections. In practice most
often \emph{median}, and \emph{median absolute deviation from the mediam }(MAD), are used as they are robust measures.
}
\value{
  Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
  
Donoho, D.L. (1982). Breakdown properties of multivariate location estimators. Ph.D. qualifying paper. Department of Statistics, Harvard University.

Liu, R.Y. (1992). Data depth and multivariate rank tests. In: Dodge, Y. (ed.), L1-Statistics and Related Methods, North-Holland, Amsterdam, 279--294.

Liu, X. and Zuo, Y. (2014b). Computing projection depth and its associated estimators.
\emph{Statistics and Computing}  \bold{24} 51--63.

Stahel, W.A. (1981). Robust estimation: infinitesimal optimality and covariance matrix
estimators. Ph.D. thesis (in German). Eidgenossische Technische Hochschule Zurich.

Zuo, Y.J. and Lai, S.Y. (2011). Exact computation of bivariate projection depth and the
Stahel-Donoho estimator. \emph{Computational Statistics and Data Analysis}  \bold{55} 1173--1179.

}
\seealso{
\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth. 

\code{\link{depth.randomTukey}} for calculation of the random Tukey depth.

\code{\link{depth.spatial}} for calculation of spatial depth. 

\code{\link{depth.zonoid}} for calculation of zonoid depth.

}
\examples{
  # 5-dimensional normal distribution
  data <- mvrnorm(1000, rep(0, 5), 
                  matrix(c(1, 0, 0, 0, 0, 
                           0, 2, 0, 0, 0, 
                           0, 0, 3, 0, 0, 
                           0, 0, 0, 2, 0, 
                           0, 0, 0, 0, 1),
                         nrow = 5))
  x <- mvrnorm(10, rep(1, 5), 
               matrix(c(1, 0, 0, 0, 0, 
                        0, 1, 0, 0, 0, 
                        0, 0, 1, 0, 0, 
                        0, 0, 0, 1, 0, 
                        0, 0, 0, 0, 1),
                      nrow = 5))
  
  depths <- depth.projection(x, data, method = "random", num.directions = 1000)
  cat("Depths random: ", depths, "\n")
  depths <- depth.projection(x, data, method = "linearize")
  cat("Depths linearize: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
