#' Evidence estimation with ZV-CV
#'
#' The functions \code{evidence_CTI} and \code{evidence_CTI_CF} can be used to improve upon the thermodynamic integration (TI) estimate of the normalising constant with ZV-CV and CF, respectively. The functions \code{evidence_SMC} and \code{evidence_SMC_CF} do the same thing for the sequential Monte Carlo (SMC) normalising constant identity.
#'
#' @param samples       An \eqn{N} by \eqn{d} by \eqn{T} matrix of samples from the \eqn{T} power posteriors, where \eqn{N} is the number of samples and \eqn{d} is the dimension of the target distribution
#' @param loglike       An \eqn{N} by \eqn{T} matrix of log likelihood values corresponding to \code{samples}
#' @param der_loglike   An \eqn{N} by \eqn{d} by \eqn{T} matrix of the derivatives of the log likelihood with respect to the parameters, with parameter values corresponding to \code{samples}
#' @param der_logprior   An \eqn{N} by \eqn{d} by \eqn{T} matrix of the derivatives of the log prior with respect to the parameters, with parameter values corresponding to \code{samples}
#' @param temperatures      A vector of length \eqn{T} of temperatures for the power posterior temperatures
#' @param temperatures_all  An adjusted vector of length \eqn{tau} of temperatures. Better performance should be obtained with a more conservative temperature schedule. See \code{\link{Expand_Temperatures}} for a function to adjust the temperatures.
#' @param most_recent   A vector of length \eqn{tau} which gives the indices in the original temperatures that the new temperatures correspond to.
#' @param est_inds     (optional) A vector of indices for the estimation-only samples. The default when \code{est_inds} is missing or \code{NULL} is to perform both estimation of the control variates and evaluation of the integral using all samples. Otherwise, the samples from \code{est_inds} are used in estimating the control variates and the remainder are used in evaluating the integral. Splitting the indices in this way can be used to reduce bias from adaption and to make computation feasible for very large sample sizes (small \code{est_inds} is faster), but in general in will increase the variance of the estimator.
#' @param options       A list of control variate specifications for ZV-CV. This can be a single list containing the elements below (the defaults are used for elements which are not specified). Alternatively, it can be a list of lists containing any or all of the elements below. Where the latter is used, the function \code{zvcv} automatically selects the best performing option based on cross-validation. 
#' @param folds  The number of folds used in k-fold cross-validation for selecting the optimal control variate. For ZV-CV, this may include selection of the optimal polynomial order, regression type and subset of parameters depending on \code{options}. For CF, this includes the selection of the optimal tuning parameters in \code{sigma_list}. The default is five.
#' 
#' @return 				The function \code{evidence_CTI}  returns a list, containing the following components:
#' \itemize{
#' \item \code{log_evidence_PS1}: The 1st order quadrature estimate for the log normalising constant
#' \item \code{log_evidence_PS2}: The 2nd order quadrature estimate for the log normalising constant
#' \item \code{regression_LL}: The set of \eqn{tau} \code{zvcv} type returns for the 1st order quadrature expectations
#' \item \code{regression_vLL}: The set of \eqn{tau} \code{zvcv} type returns for the 2nd order quadrature expectations
#' }
#'
#' @inheritSection K0_fn On the choice of \eqn{\sigma}, the kernel and the Stein order
#'
#' @references
#' Mira, A., Solgi, R., & Imparato, D. (2013). Zero variance Markov chain Monte Carlo for Bayesian estimators. Statistics and Computing, 23(5), 653-662.
#'
#' South, L. F., Oates, C. J., Mira, A., & Drovandi, C. (2019). Regularised zero variance control variates for high-dimensional variance reduction. \url{https://arxiv.org/abs/1811.05073}
#' 
#' @author  Leah F. South
#' @seealso See an example at \code{\link{VDP}} and see \link{ZVCV} for more package details. See \code{\link{Expand_Temperatures}} for a function that can be used to find stricter (or less stricter) temperature schedules based on the conditional effective sample size. 
#' 
#' @name evidence
evidence_CTI <- function(samples, loglike, der_loglike, der_logprior, temperatures, temperatures_all, most_recent, est_inds, options, folds = 5){# =  = list(polyorder = 2, regul_reg = TRUE, alpha_elnet = 1, nfolds = 10, apriori, intercept = TRUE)){
	# Stepping stone identity for evidence.
	
	N <- NROW(samples)
	d <- NCOL(samples)
	TT <- length(temperatures_all)
	
	options <- clean_options(options,N,d)
	
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
	
	#First, second and third order polynomial control variates
	regression_LL <- list()
	regression_vLL <- list()
	expectation_LL <- rep(0,TT)
	expectation_vLL <- rep(0,TT)
	
	for (tt in 1:TT){
		derivatives <- temperatures_all[tt]*der_loglike[,,tt] + der_logprior[,,tt]
		
		integrand <- loglike[,tt]
		
		regression_LL[[tt]] <- zvcv(integrand, samples[,,tt], derivatives, log_weights[,tt], est_inds = est_inds, options = options, folds = folds)
		expectation_LL[tt] <- regression_LL[[tt]]$expectation
		
		if (is.na(expectation_LL[tt])==FALSE){
			integrand <- (loglike[,tt] - expectation_LL[tt])^2
			regression_vLL[[tt]] <- zvcv(integrand, samples[,,tt], derivatives, log_weights[,tt], est_inds = est_inds, options = options, folds = folds)
			expectation_vLL[tt] <- regression_vLL[[tt]]$expectation
		} else {
			expectation_vLL[tt] <- NaN
		}
	}
	
	log_evidence_PS1 <- sum( quad1(temperatures_all,expectation_LL) ) # first order approximation
	log_evidence_PS2 <- log_evidence_PS1 - sum( quad2(temperatures_all,expectation_vLL) ) # full second order approximation
	
	return(list(log_evidence_PS1 = log_evidence_PS1, log_evidence_PS2 = log_evidence_PS2, 
							regression_LL = regression_LL, regression_vLL = regression_vLL))
	
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


