\name{ddalpha-package}
\alias{ddalpha-package}
\alias{ddalpha}
\docType{package}
\title{
Depth-Based Classification and Calculation of Data Depth
}
\description{
The package provides many procedures for calculating the depth of points in an empirical distribution for many notions of data depth. 
Further it provides implementations for depth-based classification, for multivariate and functional data.

The package implements the DD\eqn{\alpha}-classifier (Lange, Mosler and Mozharovskyi, 2014), a nonparametric procedure for supervised binary classification with \eqn{q\ge 2} classes. In the training step, the sample is first transformed into a \eqn{q}-dimensional cube of depth vectors, then a linear separation rule in its polynomial extension is constructed with the \eqn{\alpha}-procedure. The classification step involves alternative treatments of 'outsiders'.
}
\details{
\tabular{ll}{
Package: \tab ddalpha\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2016-07-25\cr
License: \tab GPL-2\cr
}
Use \code{\link{ddalpha.train}} to train the DD-classifier and \code{\link{ddalpha.classify}} to classify with it. 
Load sample classification problems using \code{\link{getdata}}. The package contains 50 classification problems built of 33 sets of real data.

The list of the implemented depths is found in topic \code{\link{depth.*}}. Corresponding to them depth representations of the data are obtained with \code{\link{depth.space.*}}. Functions \code{\link{depth.contours}} and \code{\link{depth.contours.ddalpha}} build depth contours, and \code{\link{depth.graph}} builds depth graphs for two-dimensional data. Function \code{\link{draw.ddplot}} draws DD-plot for the existing DD-classifier, or for pre-calculated depth space.

The package supports user-defined depths and classifiers, see topic \code{\link{Custom Methods}}. A pre-calculated DD-plot may also be used as \code{data}, see topic \code{\link{ddalpha.train}}.


\code{\link{is.in.convex}} shows whether an object is no 'outsider', i.e. can be classified by its depth values. 
Outsiders are alternatively classified by LDA, kNN and maximum Mahalanobis depth as well as by random assignment.

Use \code{\link{compclassf.train}} and \code{\link{ddalphaf.train}} to train the functional DD-classifiers and \code{\link{compclassf.classify}} \code{\link{ddalpha.classify}} to classify with them. Load sample functional classification problems with \code{\link{dataf.*}}. The package contains 4 functional data sets and 2 data set generators. The functional data are visualized with \code{\link{plotf}}.


}
\author{
Oleksii Pokotylo, <alexey.pokotylo at gmail.com>

Pavlo Mozharovskyi, <mozharovskyi at statistik.uni-koeln.de>

Rainer Dyckerhoff, <rainer.dyckerhoff at statistik.uni-koeln.de>
}
\references{
Mozharovskyi, P. (2015), \emph{Contributions to Depth-based Classification and Computation of the Tukey Depth}, Verlag Dr. Kovac (Hamburg).

Lange, T., Mosler, K. and Mozharovskyi, P. (2014), DD\eqn{\alpha}-classification of asymmetric and fat-tailed data, In: Spiliopoulou, M., Schmidt-Thieme, L., Janning, R. (eds), \emph{Data Analysis, Machine Learning and Knowledge Discovery}, Springer (Berlin), 71--78.

Lange, T., Mosler, K. and Mozharovskyi, P. (2014), Fast nonparametric classification based on data depth, \emph{Statistical Papers}, \bold{55}, 49--69.

Mozharovskyi, P., Mosler, K. and Lange, T. (2013), Classifying real-world data with the DD\eqn{\alpha}-procedure, \emph{Mimeo}.
}
\keyword{ package }
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }
\seealso{
\code{\link{ddalpha.train}}, \code{\link{ddalpha.classify}}, 

\code{\link{ddalphaf.train}}, \code{\link{ddalphaf.classify}}, \code{\link{compclassf.train}}, \code{\link{compclassf.classify}}

\code{\link{depth.*}}, \code{\link{depth.space.*}}, 

\code{\link{is.in.convex}},

\code{\link{getdata}}, \code{\link{dataf.*}},

\code{\link{plotf}}, \code{\link{depth.graph}}.
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

# Train the DDalpha-classifier
ddalpha <- ddalpha.train(data$train)

# Classify by means of DDalpha-classifier
classes <- ddalpha.classify(data$test[,propertyVars], ddalpha)
cat("Classification error rate:", 
    sum(unlist(classes) != data$test[,classVar])/200, "\n")

# Calculate zonoid depth of top 10 testing objects w.r.t. 1st class
depths.zonoid <- depth.zonoid(data$test[1:10,propertyVars], 
                              data$train[trainIndices,propertyVars])
cat("Zonoid depths:", depths.zonoid, "\n")

# Calculate the random Tukey depth of top 10 testing objects w.r.t. 1st class
depths.halfspace <- depth.halfspace(data$test[1:10,propertyVars], 
                                        data$train[trainIndices,propertyVars])
cat("Random Tukey depths:", depths.halfspace, "\n")

# Calculate depth space with zonoid depth
dspace.zonoid <- depth.space.zonoid(data$train[,propertyVars], c(100, 100))

# Calculate depth space with the exact Tukey depth
dspace.halfspace <- depth.space.halfspace(data$train[,propertyVars], c(100, 100), exact = TRUE)

# Count outsiders
numOutsiders = sum(rowSums(is.in.convex(data$test[,propertyVars], 
                                data$train[,propertyVars], c(100, 100))) == 0)
cat(numOutsiders, "outsiders found in the testing sample.\n")
}
