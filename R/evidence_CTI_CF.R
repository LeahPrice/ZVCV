# The function \code{evidence_CTI_CF} uses CF on the thermodynamic integration estimator for the normalising constant.
#'
#' @return 				The function \code{evidence_CTI_CF}  returns a list, containing the following components:
#' \itemize{
#' \item \code{log_evidence_PS1}: The 1st order quadrature estimate for the log normalising constant
#' \item \code{log_evidence_PS2}: The 2nd order quadrature estimate for the log normalising constant
#' \item \code{regression_LL}: The set of \eqn{tau} \code{CF_crossval} type returns for the 1st order quadrature expectations
#' \item \code{regression_vLL}: The set of \eqn{tau} \code{CF_crossval} type returns for the 2nd order quadrature expectations
#' \item \code{selected_LL_CF}: The set of \eqn{tau} selected tuning parameters from \code{sigma_list} for the 1st order quadrature expectations.
#' \item \code{selected_vLL_CF}: The set of \eqn{tau} selected tuning parameters from \code{sigma_list} for the 2nd order quadrature expectations.
#' }
#'
#' @inheritParams CF_crossval
#' @rdname evidence
evidence_CTI_CF <- function(samples, loglike, der_loglike, der_logprior, temperatures, temperatures_all, most_recent, est_inds, steinOrder, kernel_function, sigma_list, folds = 5){
  # Stepping stone identity for evidence.
  
  N <- NROW(samples)
  d <- NCOL(samples)
  TT <- length(temperatures_all)
  
  if (missing(est_inds)) {
    est_inds <- NULL
  }
  
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
      
      integrand <- loglike[inds_unique,tt]
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
      
      integrand <- loglike[-to_remove,tt]
    }
    
    Z <- squareNorm(samples_curr)
    K0_list <- list()
    for (k in 1:length(sigma_list)){
      K0_list[[k]] <- K0_fn(samples_curr,derivatives,sigma_list[[k]],steinOrder,kernel_function, Z = Z)
    }
    
    regression_LL[[tt]] <- CF_crossval(integrand, samples_curr, derivatives, K0_list = K0_list, folds = folds, est_inds = est_inds, log_weights = log_weights_curr)
    selected_LL_CF[tt] <- regression_LL[[tt]]$optinds
    expectation_LL[tt] <- regression_LL[[tt]]$expectation
    
    if (is.na(expectation_LL[tt])==FALSE){
      if (is.null(est_inds)){
        integrand <- (loglike[inds_unique,tt] - expectation_LL[tt])^2
      } else{
        integrand <- (loglike[-to_remove,tt] - expectation_LL[tt])^2
      }
      regression_vLL[[tt]] <- CF_crossval(integrand, samples_curr, derivatives, K0_list = K0_list, folds = folds, est_inds = est_inds, log_weights = log_weights_curr)
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


