# *********************************************************************************************************************************************
# Methodology for the article "Online multivariate changepoint detection with type I error control and constant time/memory updates per series"
# 
# The package provides functionality to detect multiple changepoints with O(p) effort per time step, where p is the number of time series.
# The main functions provided in the package are:
# 
# 1) A constructor "fastOnlineCpt" to create a new object of the class "fastOnlineCpt".
# 2) A function "addData" to add new incoming data for one or more time points.
# 3) A function "checkCpt" to perform a test for a changepoint.
# 4) A function "lastCptTest" to query the last result of the changepoint test.
# 5) A function "resetAlgorithm" to reset the algorithm in order to detect a new changepoint. The algorithm can be reset at any point in time.
# 
# *********************************************************************************************************************************************



#' S4 class providing functionality to detect multivariate changepoints in an online setting.
#' 
#' @description Provides method "addData" to add new incoming data for one or more time points, "checkCpt" to test for a changepoint, "lastCptTest" to query the last test result of the function "checkCpt", and "resetAlgorithm" to reset the algorithm in order to detect a new changepoint.
#' 
#' @slot spending_sequence A function handle which returns a testing level used in multiple testing.
#' @slot data Environment variable to store incoming data as a matrix.
#' @slot T The current time point.
#' @slot S Internal variable of the algorithm (modified cusum statistic).
#' @slot s Internal variable of the algorithm (modified cusum statistic).
#' @slot nTest Internal variable of the algorithm (counter for the multiple testing correction).
#' @slot lastCptTest Internal variable to store the last test result, which can be queried with a member function.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
#' 
#' @examples
#' library(fastOnlineCpt)
#' alpha <- 0.01
#' halfspent <- 100
#' spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
#' obj <- fastOnlineCpt(spending_sequence)
#' 
#' @export
setClass("fastOnlineCpt", representation=representation(spending_sequence="function", data="environment", T="numeric", S="numeric", s="numeric", nTest="numeric", lastCptTest="environment"))



#' Initialize a new object of the class "fastOnlineCpt". This object allows one to add data in an online fashion and test for a changepoint.
#' 
#' @param spending_sequence A function \eqn{f(n)} of one argument which for every \eqn{n} returns a testing threshold such that \eqn{\sum_{n=1}^\infty f(n) = \alpha}, where \eqn{\alpha} is the desired level of the test over the (possibly infinite) time horizon.
#' 
#' @return A new object of the class "fastOnlineCpt".
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
#' 
#' @examples
#' library(fastOnlineCpt)
#' alpha <- 0.01
#' halfspent <- 100
#' spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
#' obj <- fastOnlineCpt(spending_sequence)
#' 
#' @export
fastOnlineCpt <- function(spending_sequence) {
	obj <- new("fastOnlineCpt")
	obj@spending_sequence <- spending_sequence
	obj@data$x <- NULL
	obj@T <- 0
	obj@S <- 0
	obj@s <- 0
	obj@nTest <- 0
	obj@lastCptTest$x <- NA
	return(obj)
}



#' Add new \eqn{p}-dimensional data point, where \eqn{p} is the number of time series being monitored.
#' 
#' @param obj An object of the class "fastOnlineCpt".
#' @param data The new data of dimension \eqn{n \times p} to be added, where \eqn{n} is the number of new time points being added.
#' 
#' @return An object of the class "fastOnlineCpt".
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
#' 
#' @examples
#' library(fastOnlineCpt)
#' alpha <- 0.01
#' halfspent <- 100
#' spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
#' obj <- fastOnlineCpt(spending_sequence)
#' p <- 10
#' n <- 50
#' data <- matrix(rnorm(p*n,mean=0),ncol=n)
#' obj <- addData(obj,data)
#' 
#' @rdname addData-methods
#' @docType methods
#' @export
setGeneric("addData", function(obj,data) standardGeneric("addData"))
#' @rdname addData-methods
setMethod("addData", signature(obj="fastOnlineCpt"),
	function(obj,data) {
		# add data
		obj@data$x <- cbind(obj@data$x,data)
		# update T, S, and s
		while(ncol(obj@data$x)>1) {
			obj@T <- obj@T + 2
			obj@S <- obj@S + obj@data$x[,1] + obj@data$x[,2]
			obj@s <- obj@s + obj@S/obj@T
			obj@data$x <- obj@data$x[,-(1:2)]
		}
		return(obj)
	}
)



#' Test if a changepoint has occurred.
#' 
#' @param obj An object of the class "fastOnlineCpt".
#' @param screenOutput Boolean variable to indicate if the test result for a changepoint should be printed on the screen (default TRUE).
#' 
#' @return An object of the class "fastOnlineCpt".
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
#' 
#' @examples
#' library(fastOnlineCpt)
#' alpha <- 0.01
#' halfspent <- 100
#' spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
#' obj <- fastOnlineCpt(spending_sequence)
#' p <- 10
#' n <- 50
#' data <- matrix(rnorm(p*n,mean=0),ncol=n)
#' obj <- addData(obj,data)
#' obj <- checkCpt(obj)
#' data <- matrix(rnorm(p*n,mean=1),ncol=n)
#' obj <- addData(obj,data)
#' obj <- checkCpt(obj)
#' 
#' @rdname checkCpt-methods
#' @docType methods
#' @export
setGeneric("checkCpt", function(obj,screenOutput=TRUE) standardGeneric("checkCpt"))
#' @rdname checkCpt-methods
setMethod("checkCpt", signature(obj="fastOnlineCpt"),
	function(obj,screenOutput=TRUE) {
		# compute V and Z
		V2 <- obj@T * (obj@s/(obj@T/2) - obj@S/obj@T)**2
		Z <- sum(V2)
		p <- length(V2)
		pval <- exp( -0.5 * ((sqrt(2*Z)-sqrt(2*p))**2 + log(2*pi)) )
		obj@nTest <- obj@nTest + 1
		testing_level <- obj@spending_sequence(obj@nTest)
		location <- ifelse(pval < testing_level, obj@T/2, NA)
		obj@lastCptTest$x <- c(Z,pval,testing_level,location)
		if(screenOutput) {
			cat("Testing for a changepoint\n")
			cat("*************************\n")
			cat("Time point:          ",obj@T,"\n")
			cat("Number of the test:  ",obj@nTest,"\n")
			cat("Test statistic:      ",Z,"\n")
			cat("P-value:             ",pval,"\n")
			cat("Testing level:       ",testing_level,"\n")
			cat("Changepoint location:",location,"\n")
		}
		return(obj)
	}
)



#' Return the last result of the changepoint test performed with the function "checkCpt" as a vector.
#' 
#' @param obj An object of the class "fastOnlineCpt".
#' 
#' @return A 5-dimensional vector containing the number of the test, the value of the Z-statistic, the p-value, the available testing level, and the changepoint location if a changepoint has been detected or NA otherwise. If no previous test has been performed, NA is returned.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
#' 
#' @examples
#' library(fastOnlineCpt)
#' alpha <- 0.01
#' halfspent <- 100
#' spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
#' obj <- fastOnlineCpt(spending_sequence)
#' p <- 10
#' n <- 50
#' data <- matrix(rnorm(p*n,mean=0),ncol=n)
#' obj <- addData(obj,data)
#' obj <- checkCpt(obj)
#' print(lastCptTest(obj))
#' data <- matrix(rnorm(p*n,mean=1),ncol=n)
#' obj <- addData(obj,data)
#' obj <- checkCpt(obj)
#' print(lastCptTest(obj))
#' 
#' @rdname lastCptTest-methods
#' @docType methods
#' @export
setGeneric("lastCptTest", function(obj) standardGeneric("lastCptTest"))
#' @rdname lastCptTest-methods
setMethod("lastCptTest", signature(obj="fastOnlineCpt"),
	function(obj) {
		if(length(obj@lastCptTest$x)==1) {
			return(NA)
		}
		else {
			return(c(obj@nTest,obj@lastCptTest$x))
		}
	}
)



#' Reset the algorithm in order to detect a new changepoint. The algorithm can be reset at any point in time. To ensure valid multiple testing corrections, the time horizon is not reset to zero.
#' 
#' @param obj An object of the class "fastOnlineCpt".
#' 
#' @return An object of the class "fastOnlineCpt".
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G. (2021). Online multivariate changepoint detection with type I error control and constant time/memory updates per series. Under review.
#' 
#' @examples
#' library(fastOnlineCpt)
#' alpha <- 0.01
#' halfspent <- 100
#' spending_sequence <- function(n) { (n/(n+halfspent) - (n-1)/(n-1+halfspent)) * alpha }
#' obj <- fastOnlineCpt(spending_sequence)
#' p <- 10
#' n <- 50
#' data <- matrix(rnorm(p*n,mean=0),ncol=n)
#' obj <- addData(obj,data)
#' obj <- resetAlgorithm(obj)
#' 
#' @rdname resetAlgorithm-methods
#' @docType methods
#' @export
setGeneric("resetAlgorithm", function(obj) standardGeneric("resetAlgorithm"))
#' @rdname resetAlgorithm-methods
setMethod("resetAlgorithm", signature(obj="fastOnlineCpt"),
	function(obj) {
		obj@data$x <- NULL
		obj@T <- 0
		obj@S <- 0
		obj@s <- 0
		obj@lastCptTest$x <- NA
		return(obj)
	}
)
