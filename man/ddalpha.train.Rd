\name{ddalpha.train}
\alias{ddalpha.train}
\title{
Train DD\eqn{\alpha}-Classifier
}
\description{
Trains the DD\eqn{\alpha}-classifier (Lange, Mosler and Mozharovskyi, 2014; Mozharovskyi, Mosler and Lange, 2013) using a training sample according to given parameters. The DD\eqn{\alpha}-classifier is a nonparametric procedure that first transforms the training sample into the depth space calculating for that depth of each point w.r.t each class (dimension of this space equals the number of classes in the training sample then), and then constructs a linear separating rule in the polynomial extension of the depth space with the \eqn{\alpha}-procedure (Vasil'ev, 2003); maximum degree of the polynomial products is determined via cross-validation (in the depth space). If in the classification phase an object does not belong to the convex hull of at least one class (we mention such an object as an 'outsider'), it is mapped into the origin of the depth space and hence cannot be classified in the depth space. For these objects, after 'outsiderness' has been assured, an outsider treatment, i.e. a classification procedure functioning outside convex hulls of the classes is applied; it has to be trained first too.

The current realization of the DD\eqn{\alpha}-classifier allows for several alternative outsider treatments; they involve different traditional classification methods, see 'Details' and 'Arguments' for parameters needed. 

The function allows for classification with \eqn{q\ge 2} classes, see \code{aggregation.method} in 'Arguments'.
}
\usage{
ddalpha.train(data, 
              depth = "randomTukey", 
              separator = "alpha", 
              outsider.methods = "LDA", 
              outsider.settings = NULL, 
              aggregation.method = "majority", 
              knnrange = NULL, 
              num.chunks = 10, 
              num.directions = 1000, 
              use.convex = FALSE, 
              max.degree = 3, 
              mah.estimate = "moment", 
              mah.parMcd = 0.75, 
              mah.priors = NULL)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is one object of the training sample where first \eqn{d} entries are inputs and the last entry is output (class label).
}
  \item{depth}{
Character string determining which depth notion to use; can be \code{"randomTukey"} (the default) or \code{"zonoid"}.
}
  \item{separator}{
The method used for separation on the DD-plot; can be \code{"alpha"} (the default) or \code{"knnlm"}.
}
  \item{outsider.methods}{
Vector of character strings each being a name of a basic outsider method for eventual classification; possible names are: \code{"LDA"} (the default), \code{"kNN"}, \code{"kNNAff"}, \code{"depth.Mahalanobis"}, \code{"RandProp"}, \code{"RandEqual"} and \code{"Ignore"}. Each method can be specified only once, replications are ignored. By specifying treatments in such a way only a basic treatment method can be chosen (by the name), and the default settings for each of the methods are applied, see 'Details'.
}
  \item{outsider.settings}{
List containing outsider treatments each described by a list of parameters including a name, see 'Details' and 'Examples'. Each method can be used multiply with (not necessarily) different parameters, just the name should be unique, entries with the repeating names are ignored.
}
  \item{aggregation.method}{
Character string determining which method to apply to aggregate binary classification results during multiclass classification; can be \code{"majority"} (the default) or \code{"sequent"}. If \code{"majority"}, \eqn{q(q-1)/2} (with \eqn{q} being the number of classes in the training sample) binary classifiers are trained, the classification results are aggregated using the majority voting, where classes with larger proportions in the training sample (eventually with the earlier entries in the \code{data}) are preferred when tied. If \code{"sequent"}, \eqn{q} binary 'one against all'-classifiers are trained and ties during the classification are resolved as before.
}
  \item{knnrange}{
The maximal number of neighbours for kNN separation.
}
  \item{num.chunks}{
Number of chunks to split data into when cross-validating the \eqn{\alpha}-procedure; should be \eqn{>0}, and smaller than the total number of points in the two smallest classes when \code{aggregation.method =} \code{"majority"} and smaller than the total number of points in the training sample when \code{aggregation.method =} \code{"sequent"}.
}
  \item{num.directions}{
Number of directions to use when calculating the random Tukey depth (i.e. when \code{depth =} \code{"randomTukey"}); should be \eqn{>1}.
}
  \item{use.convex}{
Logical variable indicating whether outsiders should be determined exactly, i.e. as the points not contained in any of the convex hulls of the classes from the training sample (\code{TRUE}), or those having zero depth w.r.t. each class from the training sample (\code{FALSE}). For \code{depth =} \code{"zonoid"} both values give the same result.
}
  \item{max.degree}{
Maximum of the range of degrees of the polynomial depth space extension over which the \eqn{\alpha}-procedure is to be cross-validated; can be 1, 2 or 3.
}
  \item{mah.estimate}{
%%     ~~Describe \code{num.directions} here~~
}
  \item{mah.parMcd}{
%%     ~~Describe \code{use.convex} here~~
}
  \item{mah.priors}{
%%     ~~Describe \code{max.degree} here~~
}
}
\details{
An outsider treatment is a supplementary classifier for data that lie outside the convex hulls of all \eqn{q} training classes.
Available methods are: Linear Discriminant Analysis (referred to as "LDA"), see \code{\link{lda}}; \eqn{k}-Nearest-Neighbor Classifier ("kNN"), see \code{\link{knn}}, \code{\link{knn.cv}}; Affine-Invariant kNN ("kNNAff"), an affine-invariant version of the kNN, suited only for binary classification (some aggregation is used with multiple classes) and not accounting for ties (at all), but very fast by that; Maximum Mahalanobis Depth Classifier ("depth.Mahalanobis"), the outsider is referred to a class w.r.t. which it has the highest depth value scaled by (approximated) priors; Proportional Randomization ("RandProp"), the outsider is referred to a class randomly with probability equal to it (approximated) prior; Equal Randomization ("RandEqual"), the outsider is referred to a class randomly, chances for each class are equal; Ignoring ("Ignore"), the outsider is not classified, the string "Ignored" is returned instead.

An outsider treatment is specified by a list containing a name and parameters:

\code{name} is a character string, name of the outsider treatment to be freely specified; should be unique; is obligatory.

\code{method} is a character string, name of the method to use, can be \code{"LDA"}, \code{"kNN"}, \code{"kNNAff"}, \code{"depth.Mahalanobis"}, \code{"RandProp"}, \code{"RandEqual"} and \code{"Ignore"}; is obligatory.

\code{priors} is a numerical vector specifying prior probabilities of classes; class portions in the training sample are used by the default. \code{priors} is used in methods "LDA", "depth.Mahalanobis" and "RandProp".

\code{knn.k} is the number of the nearest neighbors taken into account; can be between \eqn{1} and the number of points in the training sample. Set to \eqn{-1} (the default) to be determined by the leave-one-out cross-validation. \code{knn.k} is used in method "kNN".

\code{knn.range} is the upper bound on the range over which the leave-one-out cross-validation is performed (the lower bound is \eqn{1}); can be between \eqn{2} and the number of points in the training sample \eqn{-1}. Set to \eqn{-1} (the default) to be calculated automatically accounting for number of points and dimension. \code{knn.range} is used in method "kNN".

\code{knnAff.methodAggregation} is a character string specifying the aggregation technique for method "kNNAff"; works in the same way as the function argument \code{aggregation.method}. \code{knnAff.methodAggregation} is used in method "kNNAff".

\code{knnAff.k} is the number of the nearest neighbors taken into account; should be at least \eqn{1} and up to the number of points in the training sample when \code{knnAff.methodAggregation =} \code{"sequent"}, and up to the total number of points in the training sample when \code{knnAff.methodAggregation =} \code{"majority"}. Set to \eqn{-1} (the default) to be determined by the leave-one-out cross-validation. \code{knnAff.k} is used in method "kNNAff".

\code{knnAff.range} is the upper bound on the range over which the leave-one-out cross-validation is performed (the lower bound is \eqn{1}); should be \eqn{>1} and smaller than the total number of points in the two smallest classes when \code{knnAff.methodAggregation =} \code{"majority"}, and \eqn{>1} and smaller than the total number of points in the training sample when \code{knnAff.methodAggregation =} \code{"sequent"}. Set to \eqn{-1} to be calculated automatically accounting for number of points and dimension. \code{knnAff.range} is used in method "kNNAff".

\code{mah.estimate} is a character string specifying which estimates to use when calculating the Mahalanobis depth; can be \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. \code{mah.estimate} is used in method "depth.Mahalanobis".

\code{mcd.alpha} is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used in method "depth.Mahalanobis" when \code{mah.estimate =} \code{"MCD"}.
}
\value{
Trained DD\eqn{\alpha}-classifier containing following - rather informative - fields:
\item{num.points}{Total number of points in the training sample.}
\item{dimension}{Dimension of the original space.}
\item{depth}{Character string determining which depth notion to use.}
\item{methodAggregation}{Character string determining which method to apply to aggregate binary classification results.}
\item{num.chunks}{Number of chunks data has been split into when cross-validating the \eqn{\alpha}-procedure.}
\item{num.directions}{Number of directions used for approximating the Tukey depth (when it is used).}
\item{use.convex}{Logical variable indicating whether outsiders should be determined exactly when classifying.}
\item{max.degree}{Maximum of the range of degrees of the polynomial depth space extension over which the \eqn{\alpha}-procedure has been cross-validated.}
\item{patterns}{Classes of the training sample.}
\item{num.classifiers}{Number of binary classifiers trained.}
\item{outsider.methods}{Treatments to be used to classify outsiders.}
}
\references{
Dyckerhoff, R., Koshevoy, G. and Mosler, K. (1996), Zonoid data depth: theory and computation. In: Prat A. (ed), \emph{COMPSTAT 1996. Proceedings in computational statistics}, Physica-Verlag (Heidelberg), 235--240.

Lange, T., Mosler, K. and Mozharovskyi, P. (2014), Fast nonparametric classification based on data depth, \emph{Statistical Papers}, \bold{55}, 49--69.

Mozharovskyi, P., Mosler, K. and Lange, T. (2013), Classifying real-world data with the DD\eqn{\alpha}-procedure, \emph{Mimeo}.

Vasil'ev, V.I. (2003), The reduction principle in problems of revealing regularities I, \emph{Cybernetics and Systems Analysis}, \bold{39}, 686--694.
}
\author{
The algorithm for computation of zonoid depth (Dyckerhoff, Koshevoy and Mosler, 1996) has been implemented in C++ by Rainer Dyckerhoff.
}
\seealso{
\code{\link{ddalpha.classify}} for classification using DD\eqn{\alpha}-classifier, \code{\link{depth.space.zonoid}} and \code{\link{depth.space.randomTukey}} for calculation of depth spaces, \code{\link{is.in.convex}} to check whether a point is not an outsider.
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

# Train 1st DDalpha-classifier (default settings) 
# and get the classification error rate
ddalpha1 <- ddalpha.train(data$train)
classes1 <- ddalpha.classify(data$test[,propertyVars], ddalpha1)
cat("1. Classification error rate (defaults): ", 
    sum(unlist(classes1) != data$test[,classVar])/200, ".\n", sep = "")

# Train 2nd DDalpha-classifier (zonoid depth, maximum Mahalanobis 
# depth classifier with defaults as outsider treatment) 
# and get the classification error rate
ddalpha2 <- ddalpha.train(data$train, depth = "zonoid", 
                          outsider.methods = "depth.Mahalanobis")
classes2 <- ddalpha.classify(data$test[,propertyVars], ddalpha2, 
                               outsider.method = "depth.Mahalanobis")
cat("2. Classification error rate (depth.Mahalanobis): ", 
    sum(unlist(classes2) != data$test[,classVar])/200, ".\n", sep = "")

# Train 3rd DDalpha-classifier (100 random directions for the Tukey depth, 
# adjusted maximum Mahalanobis depth classifier 
# and equal randomization as outsider treatments) 
# and get the classification error rates
treatments <- list(list(name = "mahd1", method = "depth.Mahalanobis", 
                        mah.estimate = "MCD", mcd.alpha = 0.75, priors = c(1, 1)/2), 
                   list(name = "rand1", method = "RandEqual"))
ddalpha3 <- ddalpha.train(data$train, outsider.settings = treatments, 
                          num.direction = 100)
classes31 <- ddalpha.classify(data$test[,propertyVars], ddalpha3, 
                              outsider.method = "mahd1")
classes32 <- ddalpha.classify(data$test[,propertyVars], ddalpha3, 
                              outsider.method = "rand1")
cat("3. Classification error rate (by treatments):\n")
cat("   Error (mahd1): ", 
    sum(unlist(classes31) != data$test[,classVar])/200, ".\n", sep = "")
cat("   Error (rand1): ", 
    sum(unlist(classes32) != data$test[,classVar])/200, ".\n", sep = "")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }
