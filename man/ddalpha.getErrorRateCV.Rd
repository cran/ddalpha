\name{ddalpha.getErrorRateCV}
\alias{ddalpha.getErrorRateCV}
\title{
Test DD-Classifier
}
\description{
Performs a cross-validation procedure over the given data. 
On each of \code{times} steps \code{size} observations are removed from the data, the DD-Classifier is trained on these data and tested on the removed observations.
}
\usage{
ddalpha.getErrorRateCV(data, size = 0.3, times = 10,  ...)
}
\arguments{
  \item{data}{
Matrix containing training sample where each of \eqn{n} rows is one object of the training sample where first \eqn{d} entries are inputs and the last entry is output (class label).
}
  \item{size}{
  the excluded sequences size. Either an integer between \eqn{1} and \eqn{n}, or a fraction of data between \eqn{0} and \eqn{1}.
}
  \item{times}{
  the number of times the classifier is trained.
}
  \item{\dots}{
additional parameters passed to \code{\link{ddalpha.train}}
}
}

\value{

  \item{errors}{
  the part of incorrectly classified data
  }
  \item{time}{
  the mean training time
  }
  \item{time_sd}{
  the standart seviation training time
  }

}


\seealso{
\code{\link{ddalpha.train}} to train the DD\eqn{\alpha}-classifier, 
\code{\link{ddalpha.classify}} for classification using DD\eqn{\alpha}-classifier, 
\code{\link{ddalpha.test}} to test the DD-classifier on particular learning and testing data,
\code{\link{ddalpha.getErrorRateJK}} to get error rate of the DD-classifier on particular data.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 200 training objects and 200 to test with
class1 <- mvrnorm(150, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(150, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
propertyVars <- c(1:2)
classVar <- 3
data <- rbind(cbind(class1, rep(1, 150)), cbind(class2, rep(2, 150)))

# Train 1st DDalpha-classifier (default settings) 
# and get the classification error rate
stat <- ddalpha.getErrorRateCV(data, size = 10, times = 10)
cat("1. Classification error rate (defaults): ", 
    stat$error, ".\n", sep = "")

# Train 2nd DDalpha-classifier (zonoid depth, maximum Mahalanobis 
# depth classifier with defaults as outsider treatment) 
# and get the classification error rate
stat2 <- ddalpha.getErrorRateCV(data, depth = "zonoid", 
                                outsider.methods = "depth.Mahalanobis", size = 0.2, times = 10)
cat("2. Classification error rate (depth.Mahalanobis): ", 
    stat2$error, ".\n", sep = "")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ benchmark }
