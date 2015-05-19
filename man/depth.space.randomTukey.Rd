\name{depth.space.randomTukey}
\alias{depth.space.randomTukey}
\title{
Calculate Depth Space using the Random Tukey Depth
}
\description{
Calculates the representation of the training classes in depth space using the random Tukey depth.
}
\usage{
depth.space.randomTukey(data, cardinalities, num.directions = 1000, seed = 0)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{num.directions}{
Number of random directions to be generated.
}
  \item{seed}{
the random seed. The dafault value \code{seed=0} makes no changes.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.randomTukey}}, see 'References' for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Cuesta-Albertos, J.A. and Nieto-Reyes, A. (2008), The randomTukey depth, \emph{Computational Statistics and Data Analysis}, \bold{52}, 4979--4988.

Lange, T., Mosler, K. and Mozharovskyi, P. (2014), DD\eqn{\alpha}-classification of asymmetric and fat-tailed data, In: Spiliopoulou, M., Schmidt-Thieme, L., Janning, R. (eds), \emph{Data Analysis, Machine Learning and Knowledge Discovery}, Springer (Berlin), 71--78.

Lange, T., Mosler, K. and Mozharovskyi, P. (2014), Fast nonparametric classification based on data depth, \emph{Statistical Papers}, \bold{55}, 49--69.

Tukey, J.W. (1974), Mathematics and the picturing of data, In: \emph{Proceeding of the International Congress of Mathematicians}, Vancouver, 523--531.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.randomTukey}} for calculation of the random Tukey depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using the random Tukey depth
depth.space.zonoid(data, c(10, 10))
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
