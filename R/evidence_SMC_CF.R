#' The function \code{evidence_SMC_CF} uses CF on the SMC estimator for the normalising constant.
#'
#' @return 				The function \code{evidence_SMC_CF}  returns a list, containing the following components:
#' \itemize{
#' \item \code{log_evidence}: The logged SMC estimate for the normalising constant
#' \item \code{regression_SMC}: The set of \eqn{tau} \code{CF_crossval} type returns for the expectations
#' \item \code{selected_CF}: The set of \eqn{tau} selected tuning parameters from \code{sigma_list} for the expectations
#' }
#' 
#' @rdname evidence
evidence_SMC_CF <- function(samples, loglike, der_loglike, der_logprior, temperatures, temperatures_all, most_recent, est_inds, steinOrder, kernel_function, sigma_list, folds = 5){
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
    
    
    
    if (is.null(est_inds)){
      # est_inds is not specified: Remove all duplicates and adjust weights accordingly.
      inds_unique <- !duplicated(samples)
      if (sum(inds_unique)!=N){
        # adjusting the weights
        inds <- order(samples[inds_unique,1]) # find ordering of unique values (with have to undo ordering after below command which automatically reorders)
        num_dups <- data.frame(a=samples[,1]) %>% group_by(a) %>% group_size() # find the number of duplicates using dplyr
        log_weights_curr <- log_weights_curr[inds_unique] # adjusting weights
        log_weights_curr <- log_weights_curr + log(num_dups[order(inds)]) # adjusting weights
      }
      
      if (d==1){
        samples_curr <- matrix(samples[inds_unique,,tt],ncol=1)
        derivatives <- matrix(temperatures_all[tt]*der_loglike[inds_unique,,tt] + der_logprior[inds_unique,,tt],ncol=1)
      } else{
        samples_curr <- samples[inds_unique,,tt]
        derivatives <- temperatures_all[tt]*der_loglike[inds_unique,,tt] + der_logprior[inds_unique,,tt]
      }
      
      log_integrand <- loglike[inds_unique,tt]*(temperatures_all[tt+1]-temperatures_all[tt])
    } else{
      # est_inds is specified: Remove duplicates in estimation indices. Duplicates are fine in the evaluation set so no need to remove those or adjust weights.
      inds_all <- 1:N
      # can have duplicated in eval_inds but not est_inds
      to_remove <- est_inds[duplicated(samples[est_inds,,drop=FALSE])] 
      N_new <- N - length(to_remove)
      inds_all <- inds_all[-to_remove]
      inds_all[!to_remove] <- 1:N_new
      est_inds <- inds_all[est_inds]
      
      if (d==1){
        samples_curr <- matrix(samples[-to_remove,,tt],ncol=1)
        derivatives <- matrix(temperatures_all[tt]*der_loglike[-to_remove,,tt] + der_logprior[-to_remove,,tt],ncol=1)
      } else{
        samples_curr <- samples[-to_remove,,tt]
        derivatives <- temperatures_all[tt]*der_loglike[-to_remove,,tt] + der_logprior[-to_remove,,tt]
      }
      
      log_integrand <- loglike[-to_remove,tt]*(temperatures_all[tt+1]-temperatures_all[tt])
    }
    
    max_integrand <- max(log_integrand)
    
    Z <- squareNorm(samples_curr)
    K0_list <- list()
    for (k in 1:length(sigma_list)){
      K0_list[[k]] <- K0_fn(samples_curr,derivatives,sigma_list[[k]],steinOrder,kernel_function, Z = Z)
    }
    
    regression_SMC[[tt]] <- CF_crossval(exp(log_integrand-max_integrand), samples_curr, derivatives, K0_list = K0_list, folds = folds, est_inds = est_inds, log_weights = log_weights_curr)
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

