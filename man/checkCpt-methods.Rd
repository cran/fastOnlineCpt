% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastOnlineCpt.r
\docType{methods}
\name{checkCpt}
\alias{checkCpt}
\alias{checkCpt,fastOnlineCpt-method}
\title{Test if a changepoint has occurred.}
\usage{
checkCpt(obj, screenOutput = TRUE)

\S4method{checkCpt}{fastOnlineCpt}(obj, screenOutput = TRUE)
}
\arguments{
\item{obj}{An object of the class "fastOnlineCpt".}

\item{screenOutput}{Boolean variable to indicate if the test result for a changepoint should be printed on the screen (default TRUE).}
}
\value{
An object of the class "fastOnlineCpt".
}
\description{
Test if a changepoint has occurred.
}
\examples{
library(fastOnlineCpt)
alpha <- 0.01
halfspent <- 100
spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
obj <- fastOnlineCpt(spending_sequence)
p <- 10
n <- 50
data <- matrix(rnorm(p*n,mean=0),ncol=n)
obj <- addData(obj,data)
obj <- checkCpt(obj)
data <- matrix(rnorm(p*n,mean=1),ncol=n)
obj <- addData(obj,data)
obj <- checkCpt(obj)

}
\references{
Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
}
