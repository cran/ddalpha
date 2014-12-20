\name{plotf}
\alias{plotf}
\title{
Plot the functional data
}
\description{
Plots the functional data given in the form which is described in the topic \code{\link{dataf.*}}.
}
\usage{
plotf(dataf, 
      main = "Functional data", xlab = "args", ylab = "vals", 
      colors = c("red", "blue", "green", "black", "orange", "pink"))
}

\arguments{
  \item{dataf}{
The functional data as in the topic \code{\link{dataf.*}}.
}

  \item{main}{
an overall title for the plot: see \code{\link{title}}
}
  \item{xlab}{
a title for the x axis: see \code{\link{title}}
}
  \item{ylab}{
a title for the y axis: see \code{\link{title}}
}
  \item{colors}{
the colors for the classes of the data. The colors are applied to the classes sorted in alphabetical order.
}
}

\seealso{
\code{\link{dataf.*}} for functional data description
}
\examples{

\dontrun{
  ## load the Growth dataset
  dataf = dataf.growth()

  labels = unlist(dataf$labels)
  plotf(dataf, 
    main = paste("Growth: girls red (", sum(labels == "girl"), "),", 
                      " boys blue (", sum(labels == "boy"), ")", sep=""),
    xlab="Year", ylab="Height, cm",
    colors = c("blue", "red")   # in alphabetical order of class labels   
    )
}

}
\keyword{ visualization }
\keyword{ functional }
