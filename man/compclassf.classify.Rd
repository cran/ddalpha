\name{compclassf.classify}
\alias{compclassf.classify}
\title{
Classify using Functional componentwise classifier
}
\description{
Classifies data using the functional componentwise classifier.
}
\usage{
compclassf.classify(objectsf, compclassf, ...)
}
\arguments{
  \item{objectsf}{list containing lists (functions) of two vectors of equal length, named "args" and "vals": arguments sorted in ascending order and corresponding them values respectively
}
  \item{compclassf}{
Functional componentwise classifier (obtained by \code{\link{compclassf.train}}).
}
  \item{\dots}{
additional parameters, passed to the classifier, selected with parameter \code{classifier.type} in \code{\link{compclassf.train}}.
}
}



\value{
List containing class labels.
}
\references{
Delaigle, A., Hall, P. and Bathia, N. (2012). Componentwise classification and clustering of functional data. Biometrika 99 299-313.
}
\seealso{
\code{\link{compclassf.train}} to train the functional componentwise classifier.
}
\examples{

\dontrun{
## load the Growth dataset
dataf = dataf.growth()

learn = c(head(dataf$dataf, 49), tail(dataf$dataf, 34))
labels =c(head(dataf$labels, 49), tail(dataf$labels, 34)) 
test = tail(head(dataf$dataf, 59), 10)    # elements 50:59. 5 girls, 5 boys

c = compclassf.train (learn, labels, classifier.type = "ddalpha")

classified = compclassf.classify(test, c)

print(unlist(classified))

}
}

\keyword{ functional }
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }