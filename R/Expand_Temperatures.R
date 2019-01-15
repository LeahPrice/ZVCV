#' Adjusting the temperature schedule
#'
#' @description This function is used to adjust the temperature schedule so that it is more (or less) strict than the original.
#'
#' @param temperatures      A vector of length \eqn{T} temperatures for the power posterior temperatures.
#' @param loglike       An \eqn{N} by \eqn{T} matrix of log likelihood values corresponding to \code{samples}.
#' @param rho         The tolerance for the new temperatures, which are selected so that the CESS at each temperature is \eqn{\rho*N} where \eqn{N} is the population size.
#' @param bisec_tol     The tolerance for the bisection method used in selecting temperatures. The default is \code{.Machine$double.eps^0.25}

#' @return 				A list is returned, containing the following components:
#' \itemize{
#' \item \code{temperatures_all}: The new set of temperatures of length \eqn{tau}.
#' \item \code{relevant_samples}: A vector of length \eqn{tau} containinng indices to show which particle sets the new temperatures are based on.
#' \item \code{logw}: An \eqn{N} by \eqn{tau} matrix of log normalised weights of the particles
#' }
#'
#' @references
#' South, L. F., Oates, C. J., Mira, A., & Drovandi, C. (2019). Regularised zero variance control variates for high-dimensional variance reduction.
#'
#' @author Leah F. South
#' @seealso See \code{\link{evidence}} for functions to estimate the evidence, \code{\link{VDP}} for an example and \link{ZVCV} for more package details.
Expand_Temperatures  <- function(temperatures, loglike, rho, bisec_tol = .Machine$double.eps^0.25){
	# This function can be used to add extra temperatures based on a more
	# strict ESS schedule than the move.
	
	N <- nrow(loglike)
	
	tt <- 1 # the index of the new set of temperatures
	t_existing <- 1 # the index of the existing set of temperatures
	T_existing <- length(temperatures)
	logw <- matrix(-log(N),nrow = N, ncol = 1)
	ESS <- N
	CESS <- N
	
	relevant_samples <- t_existing
	temperatures_all <- temperatures[tt]
	
	while (temperatures_all[tt]!=1){
		tt <- tt+1
		
		#Testing temperatures=1
		ESSdiff <- compute_cess_diff(logw[,tt-1],1,temperatures_all[tt-1],loglike[,t_existing],rho)
		
		#Choosing next temperature
		if (ESSdiff >= 0){
			temperatures_all <- c(temperatures_all,1)
		} else {
			temperatures_new <- uniroot(compute_cess_diff, lower=temperatures_all[tt-1], upper=1, old_logw = logw[,tt-1], temp_old = temperatures_all[tt-1], log_incr_w = loglike[,t_existing], rho = rho, tol = bisec_tol)$root
			temperatures_all <- c(temperatures_all,temperatures_new)
		}
		
		CESS <- c(CESS,compute_cess_diff(logw[,tt-1],temperatures_all[tt],temperatures_all[tt-1],loglike[,t_existing],rho) + N*rho)
		
		if (temperatures_all[tt]>=temperatures[min(t_existing+1,T_existing)]){ # if greater than or equal to the next scheduled temperature
			#t_existing <- t_existing + 1
			t_existing <- max(which(temperatures<=temperatures_all[tt]))
		}
		
		if (temperatures_all[tt]==1){
			t_existing <- T_existing
		}
		
		logw <- cbind(logw,(temperatures_all[tt] - temperatures[t_existing])*loglike[,t_existing])
		logw[,tt] <- logw[,tt] - logsumexp(logw[,tt])
		
		relevant_samples <- c(relevant_samples,t_existing)
		
		ESS <- c(ESS,exp(-logsumexp(2*(logw[,tt]))))
	}
	
	return (list(temperatures_all = temperatures_all, relevant_samples = relevant_samples, logw = logw))
}

# Computes the conditional effective sample size - rho*N. Used within the Expand_Temperatures function only.
compute_cess_diff <- function(old_logw,temp_new,temp_old,log_incr_w,rho){
	N <- length(log_incr_w)
	prep1 <- old_logw + (temp_new-temp_old)*log_incr_w
	numer <- log(N) + 2*logsumexp(prep1)
	prep2 <- old_logw + 2*(temp_new-temp_old)*log_incr_w
	denom <- logsumexp(prep2)
	return(exp(numer-denom) - rho*N)
}

