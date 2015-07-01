\name{depth.space.simplicialVolume}
\alias{depth.space.simplicialVolume}
\title{
Calculate Depth Space using simplicial volume Depth
}
\description{
Calculates the representation of the training classes in depth space using simplicial volume depth.
}
\usage{
depth.space.simplicialVolume(data, cardinalities, exact = F, k = 0.05, seed = 0)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{exact}{
\code{exact=F} (by default) implies the approximative algorithm, considering \code{k} simplices, \code{exact=T} implies the exact algorithm.
}
  \item{k}{
Number (\eqn{k>1}) or portion (if \eqn{0<k<1}) of simplices that are considered if \code{exact=F}. If \eqn{k>1}, then the algorithmic complexity is polynomial in \eqn{d} but is independent of the number of observations in \code{data}, given \eqn{k}. If \eqn{0<k<1}, then the algorithmic complexity is exponential in the number of observations in \code{data}, but the calculation precision stays approximately the same.
}
  \item{seed}{
The random seed. The dafault value \code{seed=0} makes no changes.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.simplicialVolume}}, see References below for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{

Oja, H. (1983). Descriptive statistics for multivariate distributions. \emph{Statistics & Probability Letters}, \bold{1(6)}, 327--332. 

Zuo, Y.J. and Serfling, R. (2000). General notions of statistical depth function. \emph{The Annals of Statistics} \bold{28} 461--482.

}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.simplicialVolume}} for calculation of simplicial depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using Oja depth
depth.space.simplicialVolume(data, c(10, 10))
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
