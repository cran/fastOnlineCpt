% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastOnlineCpt.r
\docType{class}
\name{fastOnlineCpt-class}
\alias{fastOnlineCpt-class}
\title{S4 class providing functionality to detect multivariate changepoints in an online setting.}
\description{
Provides method "addData" to add new incoming data for one or more time points, "checkCpt" to test for a changepoint, "lastCptTest" to query the last test result of the function "checkCpt", and "resetAlgorithm" to reset the algorithm in order to detect a new changepoint.
}
\section{Slots}{

\describe{
\item{\code{spending_sequence}}{A function handle which returns a testing level used in multiple testing.}

\item{\code{data}}{Environment variable to store incoming data as a matrix.}

\item{\code{T}}{The current time point.}

\item{\code{S}}{Internal variable of the algorithm (modified cusum statistic).}

\item{\code{s}}{Internal variable of the algorithm (modified cusum statistic).}

\item{\code{nTest}}{Internal variable of the algorithm (counter for the multiple testing correction).}

\item{\code{lastCptTest}}{Internal variable to store the last test result, which can be queried with a member function.}
}}

\examples{
library(fastOnlineCpt)
alpha <- 0.01
halfspent <- 100
spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
obj <- fastOnlineCpt(spending_sequence)

}
\references{
Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
}
