\name{depth.space.zonoid}
\alias{depth.space.zonoid}
\title{
Calculate Depth Space using Zonoid Depth
}
\description{
Calculates the depth space of the training sample using zonoid depth.
}
\usage{
depth.space.zonoid(data, cardinalities)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects, representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.zonoid}}, see 'References' for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Dyckerhoff, R., Koshevoy, G. and Mosler, K. (1996), Zonoid data depth: theory and computation. In: \emph{Prat A. (ed) COMPSTAT 1996 Proceedings in computational statistics}, Physica-Verlag, Heidelberg, pp 235--240.

Koshevoy, G. and Mosler, K. (1997), Zonoid trimming for multivariate distributions, \emph{Annals of Statistics}, \bold{25}, 1998--2017.

Lange, T., Mosler, K. and Mozharovskyi, P. (2012), Fast nonparametric classification based on data depth, \emph{Statistical Papers}.

Mosler, K. (2002), \emph{Multivariate dispersion, central regions and depth: the lift zonoid approach}, Springer, New York.
}
\author{
The algorithm for computation of zonoid depth (Dyckerhoff, Koshevoy and Mosler, 1996) has been implemented in C++ by Rainer Dyckerhoff.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.zonoid}} for calculation of zonoid depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using zonoid depth
depth.space.zonoid(data, c(10, 10))
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
