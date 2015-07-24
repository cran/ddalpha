\name{draw.ddplot}
\alias{draw.ddplot}
\title{
Draw \emph{DD}-plot
}
\description{
The function draws the \emph{DD}-plot either of the existing DD\eqn{\alpha}-classifier of the depth space.
}
\usage{
draw.ddplot(ddalpha, depth.space, cardinalities, 
            main = "DD plot", xlab = "C1", ylab = "C2", 
            classes = c(1, 2), colors = c("red", "blue", "green"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddalpha}{
DD\eqn{\alpha}-classifier (obtained by \code{\link{ddalpha.train}}).
}
  \item{depth.space}{
The ready depth space obtained by \code{\link{depth.space.*}}}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{main}{
an overall title for the plot: see \code{\link{title}}
}
  \item{xlab, ylab}{
class labels
}
  \item{classes}{
vector of numbers of two classes used for depth calculation
}
  \item{colors}{
vector of the classes' colors
}
}

\seealso{
\code{\link{ddalpha.train}}

\code{\link{depth.space.*}}

}
\examples{
  data = getdata("kidney")
  
  #1. using the existing ddalpha classifier
  ddalpha = ddalpha.train(data, depth = "spatial")
  draw.ddplot(ddalpha, main = "DD-plot")
  
  #2. using depth.space.*
  # Sort the data w.r.t. classes
  data = rbind(data[data$C == 1,], data[data$C == 2,])
  cardinalities = c(sum(data$C == 1), sum(data$C == 2))
  
  dspace = depth.space.spatial(data[,-6], cardinalities = cardinalities)
  draw.ddplot(depth.space = dspace, cardinalities = cardinalities, 
              main = "DD-plot", xlab = 1, ylab = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualization }
