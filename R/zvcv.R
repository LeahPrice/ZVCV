#' ZV-CV for general expectations
#'
#' @description The function \code{zvcv} is used to perform (regularised) ZV-CV given a set of samples, derivatives and function evaluations.
#'
#' @param integrand		The integrand, i.e. a set of \eqn{N} evaluations of the function of interest.
#' @param samples       An \eqn{N} by \eqn{d} matrix of samples from the target
#' @param derivatives	An \eqn{N} by \eqn{d} matrix of derivatives of the log target with respect to the parameters
#' @param log_weight    (optional) A vector of length \eqn{N} containing log weights of the samples. The default is equal weights.
#' @param integrand_logged (optional) Sometimes it is better to input the integrand on the logged scale for stability. If the actual integrand is the exponential of \code{integrand}, then \code{integrand_logged = TRUE}. Otherwise, the default of \code{integrand_logged = FALSE} should be used.
#' @param obs_estim_choose (optional) The indices of the samples to be used in estimating the MSE for cross validation. This is not relevant if \code{options} is a single list of specifications. If \code{options} is a list of multiple control variate specifications, then the default is to use 10 percent of the full data set for estimating the MSE.
#' @param obs_estim     (optional) The indices of the samples to be used in evaluating the integrand. If this is missing or \code{NULL} then the full data set is used for both estimating the polynomial and evaluating the integrand. Otherwise, the samples from \code{obs_estim} are used in evaluating the integral and the remainder are used in estimating the polynomial. 
#' @param options       A list of control variate specifications. This can be a single list containing the elements below (the defaults are used for elements which are not specified). Alternatively, it can be a list of lists containing any or all of the elements below. Where the latter is used, the function \code{zvcv} automatically selects the best performing option based on cross-validation. 
#' \itemize{
#' \item \code{polyorder}:   The order of the polynomial, with a default of 2. Fast implementations of polynomial orders 1, 2, 3 and 4 are available but higher order polynomials are also possible.
#' \item \code{regul_reg}:   A flag for whether regularised regression is to be used. The default is TRUE, i.e. regularised regression is used.
#' \item \code{alpha_elnet}:   The alpha parameter for elastic net. The default is 1, which correponds to LASSO. A value of 0 would correspond to ridge regression.
#' \item \code{nfolds}:   The number of folds used in cross-validation to select lambda for LASSO or elastic net. The default is 10.
#' \item \code{apriori}:   The indices of the parameters to use in the polynomial. The default is to use all parameters \eqn{1:d} where \eqn{d} is the dimension of the target. If only the first and third derivatives should be used, then this would be specified by using \code{apriori} = c(1,3) (alternatively, this can be done by only using the relevant columns in \code{samples} and \code{derivatives}).
#' \item \code{intercept}:   A flag for whether the intercept should be estimated or fixed to 0. The default is to include an intercept (\code{intercept = TRUE}) as this tends to lead to better variance reductions. Note that an \code{intercept = TRUE} flag may be changed to \code{intercept = FALSE} within the function if \code{integrand_logged = TRUE} and a \code{NaN} is encountered. See South et al. (2018) for further details.
#' }
#' 
#' @return 				A list is returned, containing the following components:
#' \itemize{
#' \item \code{expectation}: The estimate of the expectation.
#' \item \code{num_select}: The number of non-zero coefficients in the polynomial.
#' \item \code{mse}: The mean square error for the evaluation set.
#' \item \code{integrand_logged}: The \code{integrand_logged} input stored for reference.
#' \item \code{obs_estim}: The \code{obs_estim} input stored for reference.
#' \item \code{polyorder}: The \code{polyorder} value used in the final estimate.
#' \item \code{regul_reg}: The \code{regul_reg} flag used in the final estimate.
#' \item \code{alpha_elnet}: The \code{alpha_elnet} value used in the final estimate.
#' \item \code{nfolds}: The \code{nfolds} value used in the final estimate.
#' \item \code{apriori}:  The \code{apriori} value used in the final estimate.
#' \item \code{intercept}: The \code{intercept} flag used in the final estimate.
#' }
#'
#' @examples
#' # Estimating the mean of theta1 when theta is bivariate normally distributed with:
#' mymean <- c(1,2)
#' mycov <- matrix(c(1,0.5,0.5,2),nrow=2)
#' 
#' # Perfect draws from the target distribution (could be replaced with approximate draws from e.g. MCMC or SMC)
#' N <- 50
#' require(mvtnorm)
#' set.seed(1)
#' samples <- rmvnorm(N, mean = mymean, sigma = mycov)
#' integrand <- samples[,1]
#' derivatives <- t( apply(samples,1,function(x) -solve(mycov)%*%(x - mymean)) ) # derivatives of Gaussian wrt x
#' 
#' # Estimates without ZV-CV (i.e. vanilla Monte Carlo integration)
#' mean(integrand) # Without the ZVCV package
#' zvcv(integrand,samples,derivatives,options = list(polyorder = 0))$expectation # With the ZVCV package
#' 
#' # ZV-CV with fixed specifications
#' sprintf("%.15f",zvcv(integrand,samples,derivatives)$expectation) # 2nd order polynomial with LASSO
#' sprintf("%.15f",zvcv(integrand,samples,derivatives,options = list(polyorder = 2, regul_reg = FALSE))$expectation) # 2nd order polynomial with OLS
#' 
#' # ZV-CV with cross validation
#' myopts <- list(list(polyorder = Inf, regul_reg = FALSE),list(polyorder = Inf)) # Choose between OLS and LASSO, with the order selected using cross validation
#' temp <- zvcv(integrand,samples,derivatives,options = myopts) 
#' temp$polyorder # The chosen control variate order
#' sprintf("%.15f",temp$expectation) # The expectation based on the minimum MSE control variate
#' temp$regul_reg # Flag for if the chosen control variate uses regularisation
#' 
#' @references
#' Mira, A., Solgi, R., & Imparato, D. (2013). Zero variance Markov chain Monte Carlo for Bayesian estimators. Statistics and Computing, 23(5), 653-662.
#'
#' South, L. F., Mira, A., & Drovandi, C. (2018). Regularised zero variance control variates.
#'
#' @author 								Leah F. South
#' @seealso 							See \code{\link{evidence}} for functions which use \code{zvcv} to estimate the normalising constant of the posterior
#' @export zvcv
zvcv <- function(integrand, samples, derivatives, log_weight, integrand_logged = FALSE, obs_estim_choose, obs_estim, options = list(polyorder = 2, regul_reg = TRUE, alpha_elnet = 1, nfolds = 10, apriori = (1:NCOL(samples)), intercept = TRUE)) {
	
	if (any(c("polyorder","regul_reg","alpha_elnet","nfolds","apriori","intercept") %in% names(options))){
		options <- rep(list(options),1)
	}
	num_options <- length(options)
	N <- NROW(integrand)
	d <- NCOL(samples)
	
	if (missing(log_weight)) { log_weight <- rep(0,N) } # weights are normalised in another function
	if (missing(obs_estim_choose) & num_options>1) { obs_estim_choose = 1:which.min(abs(cumsum(exp(log_weight - logsumexp(log_weight)))-0.1)) } # gets 10% of the samples (taking into account weighting)
	if (missing(obs_estim)) { obs_estim <- NULL }
	for (i in 1:num_options){
		if (!("polyorder" %in% names(options[[i]]))) { options[[i]]$polyorder <- 2 }
		if (!("regul_reg" %in% names(options[[i]]))) { options[[i]]$regul_reg <- TRUE }
		if (!("alpha_elnet" %in% names(options[[i]]))) { options[[i]]$alpha_elnet <- 1 }
		if (!("nfolds" %in% names(options[[i]]))) { options[[i]]$nfolds <- 10 }
		if (!("apriori" %in% names(options[[i]]))) { options[[i]]$apriori <- 1:d }
		if (!("intercept" %in% names(options[[i]]))) { options[[i]]$intercept <- TRUE}
	}
	
	
	if (num_options>1){
		res <- list()
		mse <- rep(0,num_options)
		for (i in 1:num_options){
			# If polyorder is infinity then keep increasing polynomial order until mse stops decreasing
			if (is.infinite(options[[i]]$polyorder)){
				polyorder_curr <- -1
				ret_curr <- list(mse = Inf)
				mse_diff <- -Inf
				while (mse_diff<0){
					ret_old <- ret_curr
					polyorder_curr <- polyorder_curr + 1
					ret_curr <- zv(integrand, samples, derivatives, log_weight, integrand_logged, obs_estim = obs_estim_choose, polyorder_curr, options[[i]]$regul_reg, options[[i]]$alpha_elnet, options[[i]]$nfolds, options[[i]]$apriori, options[[i]]$intercept)
					mse_diff <- ret_curr$mse - ret_old$mse
					if (is.na(mse_diff)){
						mse_diff <- Inf # stop if get a NA fit
					}
				}
				res[[i]] <- ret_old
				options[[i]]$polyorder <- polyorder_curr - 1
			} else{
				res[[i]] <- zv(integrand, samples, derivatives, log_weight, integrand_logged, obs_estim = obs_estim_choose, options[[i]]$polyorder, options[[i]]$regul_reg, options[[i]]$alpha_elnet, options[[i]]$nfolds, options[[i]]$apriori, options[[i]]$intercept)
			}
			mse[i] <- res[[i]]$mse
		}
		
		i <- which.min(mse)
		
	} else{
		i <- 1
	}
	
	res_final <- zv(integrand, samples, derivatives, log_weight, integrand_logged, obs_estim = obs_estim, options[[i]]$polyorder, options[[i]]$regul_reg, options[[i]]$alpha_elnet, options[[i]]$nfolds, options[[i]]$apriori, options[[i]]$intercept)
	
	
	return(res_final)
}
