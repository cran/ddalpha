\name{ddalpha.classify}
\alias{ddalpha.classify}
\title{
Classify using DD\eqn{\alpha}-Classifier
}
\description{
Classifies data using the DD\eqn{\alpha}-classifier and a specified outsider treatment.
}
\usage{
ddalpha.classify(objects, ddalpha, outsider.method = "LDA", use.convex = NULL)
}
\arguments{
  \item{objects}{
Matrix containing objects to be classified; each row is one \eqn{d}-dimensional object.
}
  \item{ddalpha}{
DD\eqn{\alpha}-classifier (obtained by \code{\link{ddalpha.train}}).
}
  \item{outsider.method}{
Character string, name of a treatment to be used for outsiders; one of those trained by \code{\link{ddalpha.train}}. If the treatment was specified using the argument \code{outsider.methods} then use the name of the method.
}
  \item{use.convex}{
Logical variable indicating whether outsiders should be determined as the points not contained in any of the convex hulls of the classes from the training sample (\code{TRUE}) or those having zero depth w.r.t. each class from the training sample (\code{FALSE}). For \code{depth =} \code{"zonoid"} both values give the same result. If \code{NULL} the value specified in DD\eqn{\alpha}-classifier (in \code{\link{ddalpha.train}}) is used.
}
}
\details{
Only one outsider treatment can be specified.

See Lange, Mosler and Mozharovskyi (2014) for details and additional information.
}
\value{
List containing class labels, or character string "Ignored" for the outsiders if "Ignore" was specified as the outsider treating method.
}
\references{
Dyckerhoff, R., Koshevoy, G. and Mosler, K. (1996), Zonoid data depth: theory and computation. In: Prat A. (ed), \emph{COMPSTAT 1996. Proceedings in computational statistics}, Physica-Verlag (Heidelberg), 235--240.

Lange, T., Mosler, K. and Mozharovskyi, P. (2014), Fast nonparametric classification based on data depth, \emph{Statistical Papers}, \bold{55}, 49--69.
}
\author{
The algorithm for computation of zonoid depth (Dyckerhoff, Koshevoy and Mosler, 1996) has been implemented in C++ by Rainer Dyckerhoff.
}
\seealso{
\code{\link{ddalpha.train}} to train the DD\eqn{\alpha}-classifier.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 200 training objects and 200 to test with
class1 <- mvrnorm(200, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(200, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
trainIndices <- c(1:100)
testIndices <- c(101:200)
propertyVars <- c(1:2)
classVar <- 3
trainData <- rbind(cbind(class1[trainIndices,], rep(1, 100)), 
                   cbind(class2[trainIndices,], rep(2, 100)))
testData <- rbind(cbind(class1[testIndices,], rep(1, 100)), 
                  cbind(class2[testIndices,], rep(2, 100)))
data <- list(train = trainData, test = testData)

# Train the DDalpha-Classifier (zonoid depth, maximum Mahalanobis depth 
# classifier with defaults as outsider treatment)
ddalpha <- ddalpha.train(data$train, 
                         depth = "zonoid", 
                         outsider.methods = "depth.Mahalanobis")
# Get the classification error rate
classes <- ddalpha.classify(data$test[,propertyVars], ddalpha, 
                            outsider.method = "depth.Mahalanobis")
cat("Classification error rate: ", 
    sum(unlist(classes) != data$test[,classVar])/200, ".\n", sep="")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }
