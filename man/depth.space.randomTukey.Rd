\name{depth.space.randomTukey}
\alias{depth.space.randomTukey}
\title{
Calculate Depth Space using the Random Tukey Depth
}
\description{
Calculates the depth space of the training sample using the random Tukey depth.
}
\usage{
depth.space.randomTukey(data, cardinalities, num.directions = 1000)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects, representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{num.directions}{
Number of random directions to be generated.
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

Lange, T., Mosler, K. and Mozharovskyi, P. (2012), Fast nonparametric classification based on data depth, \emph{Statistical Papers}.

Lange, T., Mosler, K. and Mozharovskyi, P. (2013), DD\eqn{\alpha}-classification of asymmetric and fat-tailed data, Proceedings of the \emph{36th Annual Conference of the German Classification Society} on "Data Analysis, Machine Learning and Knowledge Discovery", Hildesheim, 1-3 August 2012.

Tukey, J.W. (1974), Mathematics and the picturing of data, In: \emph{Proceeding of the International Congress of Mathematicians}, Vancouver, pp 523--531.
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
