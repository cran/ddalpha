\name{depth.zonoid}
\alias{depth.zonoid}
\title{
Calculate Zonoid Depth
}
\description{
Calculates the zonoid depth of points w.r.t. a multivariate data set.
}
\usage{
depth.zonoid(x, data)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
}
\details{
Calculates zonoid depth (Koshevoy and Mosler, 1997; Mosler, 2002) exactly based on the algorithm of Dyckerhoff, Koshevoy and Mosler (1996), implemented in C++ (and provided) by Rainer Dyckerhoff.
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Dyckerhoff, R., Koshevoy, G. and Mosler, K. (1996), Zonoid data depth: theory and computation. In: Prat A. (ed), \emph{COMPSTAT 1996. Proceedings in computational statistics}, Physica-Verlag (Heidelberg), 235--240.

Koshevoy, G. and Mosler, K. (1997), Zonoid trimming for multivariate distributions, \emph{Annals of Statistics}, \bold{25}, 1998--2017.

Mosler, K. (2002), \emph{Multivariate dispersion, central regions and depth: the lift zonoid approach}, Springer (New York).
}
\author{
The algorithm for computation of zonoid depth (Dyckerhoff, Koshevoy and Mosler, 1996) has been implemented in C++ by Rainer Dyckerhoff.
}
\seealso{
\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth. 

\code{\link{depth.projection}} for calculation of projection depth. 

\code{\link{depth.randomTukey}} for calculation of the random Tukey depth.

\code{\link{depth.spatial}} for calculation of spatial depth. 

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
                
depths <- depth.zonoid(x, data)
cat("Depths: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
