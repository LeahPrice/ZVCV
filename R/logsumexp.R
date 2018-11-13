#' A stable function for computing log(sum(exp(x)))
#'
#' @param x      The values for which you want to compute log(sum(exp(x)))

#' @return 				The stable result of log(sum(exp(x)))
#'
#' @name helper_functions
logsumexp <- function(x){
	myMax <- max(x)
	x <- x - myMax
	return (log(sum(exp(x))) + myMax)
}