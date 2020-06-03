#' The function \code{evidence_SMC} uses ZV-CV on the SMC estimator for the normalising constant.
#'
#' @return 				The function \code{evidence_SMC}  returns a list, containing the following components:
#' \itemize{
#' \item \code{log_evidence}: The logged SMC estimate for the normalising constant
#' \item \code{regression_SMC}: The set of \eqn{tau} \code{zvcv} type returns for the expectations
#' }
#' 
#' @rdname evidence
evidence_SMC_CF <- function(samples, loglike, der_loglike, der_logprior, temperatures, temperatures_all, most_recent, steinOrder, kernel_function, sigma_list, folds = 5){
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
	
	#First, second and third order polynomial control variates
	expectation_SMC <- selected_CF <- rep(0,TT-1)
	regression_SMC <- list()	
	for (tt in 1:(TT-1)){
		
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
		
		log_integrand <- loglike[uniqs,tt]*(temperatures_all[tt+1]-temperatures_all[tt])
		max_integrand <- max(log_integrand)
		
		Z <- squareNorm(samples_curr)
		K0_list <- list()
		for (k in 1:length(sigma_list)){
			K0_list[[k]] <- K0_fn(samples_curr,derivatives,sigma_list[[k]],steinOrder,kernel_function, Z = Z)
		}
		
		regression_SMC[[tt]] <- CF_crossval(exp(log_integrand-max_integrand), samples_curr, derivatives, K0_list = K0_list, folds = folds, input_weights = exp(log_weights_curr))
		selected_CF[tt] <- regression_SMC[[tt]]$optinds
		regression_SMC[[tt]]$expectation <- log(regression_SMC[[tt]]$expectation) + max_integrand
		expectation_SMC[tt] <- regression_SMC[[tt]]$expectation
		if (is.na(expectation_SMC[tt])){
			regression_SMC[[tt]]$expectation <- log(mean(exp(log_integrand - max_integrand))) + max_integrand
			expectation_SMC[tt] <- regression_SMC[[tt]]$expectation
		}
		
	}
	
	log_evidence <- sum(expectation_SMC)
	
	return(list(log_evidence = log_evidence, regression_SMC = regression_SMC, selected_CF = selected_CF))
	
}

