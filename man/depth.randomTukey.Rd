\name{depth.randomTukey}
\alias{depth.randomTukey}
\title{
Calculate the Random Tukey Depth
}
\description{
Calculates the random Tukey depth of points w.r.t. a multivariate data set.
}
\usage{
depth.randomTukey(x, data, num.directions = 1000, seed = 0)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{num.directions}{
Number of random directions to be generated.
}
  \item{seed}{
the random seed. The dafault value \code{seed=0} makes no changes.
}
}
\details{
Calculates Tukey (= halfspace, location) depth (Tukey, 1975) approximately using the random Tukey depth method proposed by Cuesta-Albertos and Nieto-Reyes (2008). Here the depth is determined as the minimum univariate Tukey depth of the - on lines in several directions - projected data. The directions are distributed uniformly on the \eqn{(d-1)}-sphere; the same direction set is used for all points.
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Cuesta-Albertos, J.A. and Nieto-Reyes, A. (2008), The random Tukey depth, \emph{Computational Statistics and Data Analysis}, \bold{52}, 4979--4988.

Lange, T., Mosler, K. and Mozharovskyi, P. (2014), DD\eqn{\alpha}-classification of asymmetric and fat-tailed data, In: Spiliopoulou, M., Schmidt-Thieme, L., Janning, R. (eds), \emph{Data Analysis, Machine Learning and Knowledge Discovery}, Springer (Berlin), 71--78.

Tukey, J.W. (1974), Mathematics and the picturing of data, In: \emph{Proceeding of the International Congress of Mathematicians}, Vancouver, 523--531.
}
\seealso{
\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth. 

\code{\link{depth.projection}} for calculation of projection depth. 

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
                
depths <- depth.randomTukey(x, data)
cat("Depths: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
