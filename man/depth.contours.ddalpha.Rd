\name{depth.contours.ddalpha}
\alias{depth.contours.ddalpha}
\title{
Depth contours
}
\description{
Builds the data depth contours for multiclass 2-dimensional data using the trained classifier.
}
\usage{
depth.contours.ddalpha(ddalpha,
              main = "", xlab="", ylab = "", 
              drawplot = T, frequency=100)
}

\arguments{
  \item{ddalpha}{
DD\eqn{\alpha}-classifier (obtained by \code{\link{ddalpha.train}}).
}
  \item{main}{
an overall title for the plot: see \code{\link{title}}
}
  \item{xlab, ylab}{
labels of the axes
}
  \item{drawplot}{
if set to false, the contours are built on the existing plot. 
}
  \item{frequency}{
number of points on each direction, x and y. Impacts the smoothness of the contures. 
}
}

\seealso{
\code{\link{depth.*}}

\code{\link{depth.contours}}

\code{\link{depth.graph}}
}
\examples{

\dontrun{

par(mfrow = c(2,2))
data(hemophilia)

ddalpha = ddalpha.train(hemophilia, depth = "none")
depth.contours.ddalpha(ddalpha, main = "data")

for (depth in c("zonoid", "Mahalanobis", "projection", "spatial")){
  ddalpha = ddalpha.train(hemophilia, depth = depth)
  depth.contours.ddalpha(ddalpha, main = depth)
}

for (depth in c("halfspace", "simplicial", "simplicialVolume")){
  ddalpha = ddalpha.train(hemophilia, depth = depth, exact = T)
  depth.contours.ddalpha(ddalpha, main = depth)
}


}
}

\keyword{ visualization }
