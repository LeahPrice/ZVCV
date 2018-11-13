#' Adjusting the temperature schedule
#'
#' @description This function is used to adjust the temperature schedule so that it is more (or less) strict than the original.
#'
#' @param gammavar      A vector of length \eqn{T} temperatures for the power posterior temperatures.
#' @param loglike       An \eqn{N} by \eqn{T} matrix of log likelihood values corresponding to \code{samples}.
#' @param alpha         The tolerance for the new temperatures, which are selected so that the CESS at each temperature is \eqn{\alpha*N} where \eqn{N} is the population size.
#' @param bisec_tol     The tolerance for the bisection method used in selecting temperatures. The default is \code{.Machine$double.eps^0.25}

#' @return 				A list is returned, containing the following components:
#' \itemize{
#' \item \code{gammavar_all}: The new set of temperatures of length \eqn{tau}.
#' \item \code{relevant_samples}: A vector of length \eqn{tau} containinng indices to show which particle sets the new temperatures are based on.
#' \item \code{logw}: An \eqn{N} by \eqn{tau} matrix of log normalised weights of the particles
#' }
#'
#' @references
#' South, L. F., Mira, A., & Drovandi, C. (2018). Regularised zero variance control variates.
#'
#' @author Leah F. South
Expand_Temperatures  <- function(gammavar, loglike, alpha, bisec_tol = .Machine$double.eps^0.25){
# This function can be used to add extra temperatures based on a more
# strict ESS schedule than the move.

N <- nrow(loglike)

tt <- 1 # the index of the new set of temperatures
t_existing <- 1 # the index of the existing set of temperatures
T_existing <- length(gammavar)
logw <- matrix(-log(N),nrow = N, ncol = 1)
ESS <- N
CESS <- N

relevant_samples <- t_existing
gammavar_all <- gammavar[tt]

while (gammavar_all[tt]!=1){
    tt <- tt+1
    
    #Testing gammavar=1
    ESSdiff <- compute_cess_diff(logw[,tt-1],1,gammavar_all[tt-1],loglike[,t_existing],alpha)
    
    #Choosing next temperature
    if (ESSdiff >= 0){
        gammavar_all <- c(gammavar_all,1)
    } else {
        gammavar_new <- uniroot(compute_cess_diff, lower=gammavar_all[tt-1], upper=1, old_logw = logw[,tt-1], temp_old = gammavar_all[tt-1], log_w = loglike[,t_existing], alpha = alpha, tol = bisec_tol)$root
        gammavar_all <- c(gammavar_all,gammavar_new)
    }
    
    CESS <- c(CESS,compute_cess_diff(logw[,tt-1],gammavar_all[tt],gammavar_all[tt-1],loglike[,t_existing],alpha) + N*alpha)
    
    if (gammavar_all[tt]>=gammavar[min(t_existing+1,T_existing)]){ # if greater than or equal to the next scheduled temperature
        #t_existing <- t_existing + 1
        t_existing <- max(which(gammavar<=gammavar_all[tt]))
    }
    
    if (gammavar_all[tt]==1){
        t_existing <- T_existing
    }
    
    logw <- cbind(logw,(gammavar_all[tt] - gammavar[t_existing])*loglike[,t_existing])
    logw[,tt] <- logw[,tt] - logsumexp(logw[,tt])
    
    relevant_samples <- c(relevant_samples,t_existing)
    
    ESS <- c(ESS,exp(-logsumexp(2*(logw[,tt]))))
}

return (list(gammavar_all = gammavar_all, relevant_samples = relevant_samples, logw = logw))
}
