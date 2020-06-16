# The function \code{evidence_SMC} uses ZV-CV on the SMC estimator for the normalising constant.
#'
#' @return 				The function \code{evidence_SMC}  returns a list, containing the following components:
#' \itemize{
#' \item \code{log_evidence}: The logged SMC estimate for the normalising constant
#' \item \code{regression_SMC}: The set of \eqn{tau} \code{zvcv} type returns for the expectations
#' }
#' 
#' @rdname evidence
evidence_SMC <- function(samples, loglike, der_loglike, der_logprior, temperatures, temperatures_all, most_recent, est_inds, options, folds = 5){# = list(polyorder = 2, regul_reg = TRUE, alpha_elnet = 1, nfolds = 10, apriori, intercept = TRUE)){
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
	expectation_SMC <- rep(0,TT-1)
	regression_SMC <- list()	
	for (tt in 1:(TT-1)){
		derivatives <- temperatures_all[tt]*der_loglike[,,tt] + der_logprior[,,tt]
		
		#log_integrand <- log_weights[,tt] + loglike[,tt]*(temperatures_all[tt+1]-temperatures_all[tt])
		log_integrand <- loglike[,tt]*(temperatures_all[tt+1]-temperatures_all[tt])
		
		regression_SMC[[tt]] <- zvcv(log_integrand, samples[,,tt], derivatives, log_weights[,tt], integrand_logged = TRUE, est_inds = est_inds, options = options, folds = folds)
		expectation_SMC[tt] <- regression_SMC[[tt]]$expectation
	}
	
	log_evidence <- sum(expectation_SMC)
	
	return(list(log_evidence = log_evidence, regression_SMC = regression_SMC))
	
}

