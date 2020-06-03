#' Evidence estimation with ZV-CV
#'
#' The function \code{evidence_CTI} uses ZV-CV on the controlled thermodynamic integration estimator for the normalising constant.
#'
#' @param samples       An \eqn{N} by \eqn{d} by \eqn{T} matrix of samples from the \eqn{T} power posteriors, where \eqn{N} is the number of samples and \eqn{d} is the dimension of the target distribution
#' @param loglike       An \eqn{N} by \eqn{T} matrix of log likelihood values corresponding to \code{samples}
#' @param der_loglike   An \eqn{N} by \eqn{d} by \eqn{T} matrix of the derivatives of the log likelihood with respect to the parameters, with parameter values corresponding to \code{samples}
#' @param der_logprior   An \eqn{N} by \eqn{d} by \eqn{T} matrix of the derivatives of the log prior with respect to the parameters, with parameter values corresponding to \code{samples}
#' @param temperatures      A vector of length \eqn{T} of temperatures for the power posterior temperatures
#' @param temperatures_all  An adjusted vector of length \eqn{tau} of temperatures. Better performance should be obtained with a more conservative temperature schedule. See \code{\link{Expand_Temperatures}} for a function to adjust the temperatures.
#' @param most_recent   A vector of length \eqn{tau} which gives the indices in the original temperatures that the new temperatures correspond to.
#' @param obs_estim_choose See \code{\link{zvcv}}.
#' @param obs_estim     See \code{\link{zvcv}}.
#' @param options       See \code{\link{zvcv}}.
#' 
#' @return 				The function \code{evidence_CTI}  returns a list, containing the following components:
#' \itemize{
#' \item \code{log_evidence_PS1}: The 1st order quadrature estimate for the log normalising constant
#' \item \code{log_evidence_PS2}: The 2nd order quadrature estimate for the log normalising constant
#' \item \code{regression_LL}: The set of \eqn{tau} \code{zvcv} type returns for the 1st order quadrature expectations
#' \item \code{regression_vLL}: The set of \eqn{tau} \code{zvcv} type returns for the 2nd order quadrature expectations
#' }
#'
#' @references
#' Mira, A., Solgi, R., & Imparato, D. (2013). Zero variance Markov chain Monte Carlo for Bayesian estimators. Statistics and Computing, 23(5), 653-662.
#'
#' South, L. F., Oates, C. J., Mira, A., & Drovandi, C. (2019). Regularised zero variance control variates for high-dimensional variance reduction. \url{https://arxiv.org/abs/1811.05073}
#' 
#' @author  Leah F. South
#' @seealso See \code{\link{Expand_Temperatures}} for a function that can be used to find stricter (or less stricter) temperature schedules based on the conditional effective sample size. See an example at \code{\link{VDP}} and see \link{ZVCV} for more package details.
#' 
#' @name evidence
evidence_CTI_CF <- function(samples, loglike, der_loglike, der_logprior, temperatures, temperatures_all, most_recent, steinOrder, kernel_function, sigma_list, folds = 5){
	# Stepping stone identity for evidence.
	
	N <- NROW(samples)
	d <- NCOL(samples)
	TT <- length(temperatures_all)
	
	log_weights <- matrix(,nrow=N,ncol=TT)
	log_weights[,1] <- -log(N)*rep(1,N)
	for (tt in 2:TT){
		log_weights[,tt] <- (temperatures_all[tt] - temperatures[most_recent[tt]])*loglike[,most_recent[tt]]
		log_weights[,tt] <- log_weights[,tt] - logsumexp(log_weights[,tt])
	}
	
	samples <- array(samples[,,most_recent],c(N,d,TT))
	loglike <- loglike[,most_recent]
	der_loglike <- array(der_loglike[,,most_recent],c(N,d,TT))
	der_logprior <- array(der_logprior[,,most_recent],c(N,d,TT))
	
	regression_LL <- list()
	regression_vLL <- list()
	expectation_LL <- selected_LL_CF <- rep(0,TT)
	expectation_vLL <- selected_vLL_CF <- rep(0,TT)
	
	for (tt in 1:TT){
		log_weights_curr <- log_weights[,tt]
		
		uniqs <- !duplicated(samples[,1,tt])
		if (length(uniqs!=N)){
			# adjusting the weights
			inds <- order(samples[uniqs,1,tt]) # find ordering of unique values (with have to undo ordering after below command which automatically reorders)
			num_dups <- data.frame(a=samples[,1,tt]) %>% group_by(a) %>% group_size() # find the number of duplicates using dplyr
			log_weights_curr <- log_weights_curr[uniqs] # adjusting weights
			log_weights_curr <- log_weights_curr + log(num_dups[order(inds)]) # adjusting weights 
		}
		
		if (d==1){
			samples_curr <- matrix(samples[uniqs,,tt],ncol=1)
			derivatives <- matrix(temperatures_all[tt]*der_loglike[uniqs,,tt] + der_logprior[uniqs,,tt],ncol=1)
		} else{
			samples_curr <- samples[uniqs,,tt]
			derivatives <- temperatures_all[tt]*der_loglike[uniqs,,tt] + der_logprior[uniqs,,tt]
		}
		
		integrand <- loglike[uniqs,tt]
		
		Z <- squareNorm(samples_curr)
		K0_list <- list()
		for (k in 1:length(sigma_list)){
			K0_list[[k]] <- K0_fn(samples_curr,derivatives,sigma_list[[k]],steinOrder,kernel_function, Z = Z)
		}
		
		regression_LL[[tt]] <- CF_crossval(integrand, samples_curr, derivatives, K0_list = K0_list, folds = folds, input_weights = exp(log_weights_curr))
		selected_LL_CF[tt] <- regression_LL[[tt]]$optinds
		expectation_LL[tt] <- regression_LL[[tt]]$expectation
		
		if (is.na(expectation_LL[tt])==FALSE){
			integrand <- (loglike[uniqs,tt] - expectation_LL[tt])^2
			regression_vLL[[tt]] <- CF_crossval(integrand, samples_curr, derivatives, K0_list = K0_list, folds = folds, input_weights = exp(log_weights_curr))
			selected_vLL_CF[tt] <- regression_vLL[[tt]]$optinds
			expectation_vLL[tt] <- regression_vLL[[tt]]$expectation
		} else {
			expectation_vLL[tt] <- NaN
		}
	}
	
	log_evidence_PS1 <- sum( quad1(temperatures_all,expectation_LL) ) # first order approximation
	log_evidence_PS2 <- log_evidence_PS1 - sum( quad2(temperatures_all,expectation_vLL) ) # full second order approximation
	
	return(list(log_evidence_PS1 = log_evidence_PS1, log_evidence_PS2 = log_evidence_PS2, 
							regression_LL = regression_LL, regression_vLL = regression_vLL, selected_LL_CF = selected_LL_CF, selected_vLL_CF = selected_vLL_CF))
	
}


quad1 <- function(temperatures,means){
	# Quadrature first order approximation
	
	TT <- length(temperatures)
	temp <- (temperatures[2:TT]-temperatures[1:(TT-1)])/2*(means[1:(TT-1)]+means[2:TT])
	return(temp)
}


quad2 <- function(temperatures,vars){
	# Quadrature second order approximation
	
	TT <- length(temperatures)
	temp <- (temperatures[2:TT]-temperatures[1:(TT-1)])^2/12*(vars[2:TT]-vars[1:(TT-1)])
	return(temp)
}


