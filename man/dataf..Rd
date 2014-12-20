\name{dataf.*}
\alias{dataf.*}
\docType{data}
\title{
Functional data sets
}
\description{
The functions generate data sets of functional two-dimensional data of two or more classes.
}
\usage{
# dataf.[name]()
}
\format{
  The functional data as a data structure.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates 
      \describe{
        \item{\code{args}}{The arguments vector}
        \item{\code{vals}}{The values vector}
       }
    }
    \item{\code{labels}}{The classes of the objects}
  }
}
\details{
More details about the datasets in the topics:

\code{\link{dataf.geneexp}}

\code{\link{dataf.growth}}

\code{\link{dataf.medflies}}

\code{\link{dataf.tecator}}

The following datasets provide simulated data:

\code{\link{dataf.sim.1.CFF07}}

\code{\link{dataf.sim.2.CFF07}}

}
\seealso{
\code{\link{plotf}} for building plots of functional data
}
\examples{
## load the Growth dataset
dataf = dataf.growth()

## view the classes
unique(dataf$labels)

## access the 5th point of the 2nd object
dataf$dataf[[2]]$args[5]
dataf$dataf[[2]]$vals[5]

\dontrun{plotf(dataf)}
}
\keyword{datasets}
\keyword{functional}
