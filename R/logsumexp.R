#' Useful helper functions
#'
#' The function \code{logsumexp} is used for stable computation of log(sum(exp(x))), which is useful when summing weights for example.
#'
#' @param x      The values for which you want to compute log(sum(exp(x)))

#' @return 				The stable result of log(sum(exp(x)))
#'
#' @name helper_functions
#' @seealso See \link{ZVCV} for more package details.
logsumexp <- function(x){
	myMax <- max(x)
	x <- x - myMax
	return (log(sum(exp(x))) + myMax)
}