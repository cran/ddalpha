\name{ddalphaf.classify}
\alias{ddalphaf.classify}
\title{
Classify using Functional DD\eqn{\alpha}-Classifier
}
\description{
Classifies data using the functional DD\eqn{\alpha}-classifier.
}
\usage{
ddalphaf.classify(objectsf, ddalphaf, ...)
}
\arguments{
  \item{objectsf}{list containing lists (functions) of two vectors of equal length, named "args" and "vals": arguments sorted in ascending order and corresponding them values respectively
}
  \item{ddalphaf}{
Functional DD\eqn{\alpha}-classifier (obtained by \code{\link{ddalphaf.train}}).
}
  \item{\dots}{
additional parameters, passed to the classifier, selected with parameter \code{classifier.type} in \code{\link{ddalphaf.train}}.
}
}



\value{
List containing class labels.
}
\references{
Mozharovskyi, P. (2015), \emph{Contributions to Depth-based Classification and Computation of the Tukey Depth}, Verlag Dr. Kovac (Hamburg).

Mosler, K., & Mozharovskyi, P. (2014). Fast DD-classification of functional data. arXiv preprint arXiv:1403.1158.
}
\seealso{
\code{\link{ddalphaf.train}} to train the functional DD\eqn{\alpha}-classifier.
}
\examples{

\dontrun{
## load the Growth dataset
dataf = dataf.growth()

learn = c(head(dataf$dataf, 49), tail(dataf$dataf, 34))
labels= c(head(dataf$labels, 49), tail(dataf$labels, 34)) 
test  = tail(head(dataf$dataf, 59), 10)    # elements 50:59. 5 girls, 5 boys

c = ddalphaf.train (learn, labels, classifier.type = "ddalpha")

classified = ddalphaf.classify(test, c)

print(unlist(classified))

}
}

\keyword{ functional }
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }