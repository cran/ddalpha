\name{depth.space.*}
\alias{depth.space.*}

\title{
Calculate Depth Space using the given Depth
}
\description{
Calculates the representation of the training classes in depth space.

The detailed descriptions are found in the corresponding topics.
}
\usage{
## Mahalanobis depth
# depth.space.Mahalanobis(data, cardinalities, mah.estimate = "moment", mah.parMcd = 0.75)

## projection depth
# depth.space.projection(data, cardinalities, method = "random", num.directions = 1000)

## Random Tukey depth
# depth.space.randomTukey(data, cardinalities, num.directions = 1000)

## spatial depth
# depth.space.spatial(data, cardinalities)

## zonoid depth
# depth.space.zonoid(data, cardinalities)

}

\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
}

\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}

\seealso{

\code{\link{depth.space.Mahalanobis}}

\code{\link{depth.space.projection}}

\code{\link{depth.space.randomTukey}}

\code{\link{depth.space.spatial}}

\code{\link{depth.space.zonoid}}

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
