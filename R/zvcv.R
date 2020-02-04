#' ZV-CV for general expectations
#'
#' @description The function \code{zvcv} is used to perform (regularised) ZV-CV given a set of samples, derivatives and function evaluations.
#'
#' @param integrand		An \eqn{N} by \eqn{k} matrix of integrands (evaluations of the function of interest)
#' @param samples       An \eqn{N} by \eqn{d} matrix of samples from the target
#' @param derivatives	An \eqn{N} by \eqn{d} matrix of derivatives of the log target with respect to the parameters
#' @param log_weight    (optional) A vector of length \eqn{N} containing log weights of the samples. The default is equal weights.
#' @param integrand_logged (optional) Sometimes it is better to input the integrand on the logged scale for stability. If the actual integrand is the exponential of \code{integrand}, then \code{integrand_logged = TRUE}. Otherwise, the default of \code{integrand_logged = FALSE} should be used.
#' @param obs_estim_choose (optional) The indices of the samples to be used in estimating the MSE for cross validation. This is not relevant if \code{options} is a single list of specifications. If \code{options} is a list of multiple control variate specifications, then the default is to use 10 percent of the full data set for estimating the MSE.
#' @param obs_estim     (optional) The indices of the samples to be used in evaluating the integrand. If this is missing or \code{NULL} then the full data set is used for both estimating the polynomial and evaluating the integrand. Otherwise, the samples from \code{obs_estim} are used in evaluating the integral and the remainder are used in estimating the polynomial. 
#' @param options       A list of control variate specifications. This can be a single list containing the elements below (the defaults are used for elements which are not specified). Alternatively, it can be a list of lists containing any or all of the elements below. Where the latter is used, the function \code{zvcv} automatically selects the best performing option based on cross-validation. 
#' \itemize{
#' \item \code{polyorder}:   The order of the polynomial, with a default of 2. A value of \code{Inf} will get the cross-validation method to choose between orders.
#' \item \code{regul_reg}:   A flag for whether regularised regression is to be used. The default is TRUE, i.e. regularised regression is used.
#' \item \code{alpha_elnet}:   The alpha parameter for elastic net. The default is 1, which correponds to LASSO. A value of 0 would correspond to ridge regression.
#' \item \code{nfolds}:   The number of folds used in cross-validation to select lambda for LASSO or elastic net. The default is 10.
#' \item \code{apriori}:   The indices of the parameters to use in the polynomial. The default is to use all parameters \eqn{1:d} where \eqn{d} is the dimension of the target. If only the first and third derivatives should be used, then this would be specified by using \code{apriori} = c(1,3) (alternatively, this can be done by only using the relevant columns in \code{samples} and \code{derivatives}).
#' \item \code{intercept}:   A flag for whether the intercept should be estimated or fixed to the empirical mean of the integrand in the estimation set. The default is to include an intercept (\code{intercept = TRUE}) as this tends to lead to better variance reductions. Note that an \code{intercept = TRUE} flag may be changed to \code{intercept = FALSE} within the function if \code{integrand_logged = TRUE} and a \code{NaN} is encountered. See South et al. (2018) for further details.
#' \item \code{polyorder_max}:   The maximum allowable polynomial order. This may be used to prevent memory issues in the case that the polynomial order is selected automatically. You will be prompted to see if you wish to continue if the selected polynomial order results in a design matrix with more than ten million elements. A default maximum polynomial order will be selected if the \code{polyorder} is infinite and in this case a warning will be given. Recall that setting your default R settings to \code{options(warn=1)} will ensure that you receive these warnings in real time. Optimal polynomial order selection may go to at most this maximum value, or it may stop earlier. 
#' }
#' @param folds_choose  The number of folds used in k-fold cross-validation for selecting the optimal control variate. Depending on the \code{options}, this may include selection of the optimal polynomial order, regression type and subset of parameters in the polynomial. The default is five.
#' 
#' @return 				A list is returned, containing the following components:
#' \itemize{
#' \item \code{expectation}: The estimates of the expectations.
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
#' # Estimating some expectations when theta is bivariate normally distributed with:
#' mymean <- c(1,2)
#' mycov <- matrix(c(1,0.5,0.5,2),nrow=2)
#' 
#' # Perfect draws from the target distribution (could be replaced with
#' # approximate draws from e.g. MCMC or SMC)
#' N <- 100
#' require(mvtnorm)
#' set.seed(1)
#' samples <- rmvnorm(N, mean = mymean, sigma = mycov)
#' derivatives <- t( apply(samples,1,function(x) -solve(mycov)%*%(x - mymean)) ) # derivatives of Gaussian wrt x
#' 
#' # The integrands are the marginal posterior means of theta, the variances and the covariance (true values are c(1,2,1,2,0.5))
#' integrand <- cbind(samples,(samples[,1] - 1)^2,(samples[,2] - 2)^2, (samples[,1] - 1)*(samples[,2] - 2))
#' 
#' # Estimates without ZV-CV (i.e. vanilla Monte Carlo integration)
#' # Without the ZVCV package
#' colMeans(integrand)
#' # With the ZVCV package
#' zvcv(integrand,samples,derivatives,options = list(polyorder = 0))$expectation
#' 
#' # ZV-CV with fixed specifications
#' # 2nd order polynomial with LASSO
#' sprintf("%.15f",zvcv(integrand, samples, derivatives)$expectation)
#' # 2nd order polynomial with OLS
#' sprintf("%.15f",zvcv(integrand, samples, derivatives,
#' 										 options = list(polyorder = 2, regul_reg = FALSE))$expectation) 
#' 
#' # ZV-CV with cross validation
#' # Choose between OLS and LASSO, with the order selected using cross validation
#' # For this example, polyorder = 1 with OLS is exact for the first two integrands and polyorder = 2 with OLS is exact for the last three integrands
#' myopts <- list(list(polyorder = Inf, regul_reg = FALSE),list(polyorder = Inf)) 
#' temp <- zvcv(integrand,samples,derivatives,options = myopts) 
#' temp$polyorder # The chosen control variate order
#' temp$regul_reg # Flag for if the chosen control variate uses regularisation
#' sprintf("%.15f",temp$expectation) # The expectation based on the minimum MSE control variate
#' 
#' 
#' @references
#' Mira, A., Solgi, R., & Imparato, D. (2013). Zero variance Markov chain Monte Carlo for Bayesian estimators. Statistics and Computing, 23(5), 653-662.
#'
#' South, L. F., Oates, C. J., Mira, A., & Drovandi, C. (2019). Regularised zero variance control variates for high-dimensional variance reduction. \url{https://arxiv.org/abs/1811.05073}
#'
#' @author 								Leah F. South
#' @seealso 							See \code{\link{evidence}} for functions which use \code{zvcv} to estimate the normalising constant of the posterior. See an example at \code{\link{VDP}} and see \link{ZVCV} for more package details.
#' @export zvcv
zvcv <- function(integrand, samples, derivatives, log_weight, integrand_logged = FALSE, obs_estim_choose, obs_estim, options = list(polyorder = 2, regul_reg = TRUE, alpha_elnet = 1, nfolds = 10, apriori = (1:NCOL(samples)), intercept = TRUE, polyorder_max = Inf), folds_choose = 5) {
	
	N <- NROW(integrand)
	N_expectations <- NCOL(integrand)
	d <- NCOL(samples)
	
	options <- clean_options(options,N,d)
	num_options <- length(options)
	
	if (missing(log_weight)) { log_weight <- rep(0,N) } # weights are normalised in another function
	if (missing(obs_estim_choose) & (num_options>1 | is.infinite(options[[1]]$polyorder))) {  obs_estim_choose <- split(sample(1:N),rep(1:folds_choose, ceiling(N/folds_choose),length.out = N))  }
	if (missing(obs_estim)) { obs_estim <- NULL }
	
	if (is.null(ncol(integrand))){
		integrand <- matrix(integrand,ncol=1)
	}
	
	if (num_options>1 | is.infinite(options[[1]]$polyorder)){
		mse <- matrix(NaN,nrow=num_options,ncol=N_expectations)
		for (i in 1:num_options){
			# If polyorder is infinity then keep increasing polynomial order until mse stops decreasing
			if (is.infinite(options[[i]]$polyorder)){
				options[[i]]$polyorder <- rep(NaN,N_expectations)
				polyorder_track <- 0
				polyorder_curr <- rep(0,N_expectations)
				mse_old <- rep(0,N_expectations)
				mse_new <- rep(Inf,N_expectations)
				inds_to_continue <- 1:N_expectations
				while (length(inds_to_continue)!=0){
					mse_old[inds_to_continue] <- mse_new[inds_to_continue]
					polyorder_track <- polyorder_track + 1
					polyorder_curr[inds_to_continue] <- polyorder_curr[inds_to_continue] + 1
					if  (any(polyorder_curr > options[[i]]$polyorder_max)){
						break
					}
					mse_temp <- matrix(NaN,nrow=folds_choose,ncol=length(inds_to_continue))
					for (k in 1:folds_choose){
						res_curr <- zv(integrand[,inds_to_continue], samples, derivatives, log_weight, integrand_logged, obs_estim = obs_estim_choose[[k]], polyorder_track, options[[i]]$regul_reg, options[[i]]$alpha_elnet, options[[i]]$nfolds, options[[i]]$apriori, options[[i]]$intercept)
						mse_temp[k,] <- res_curr$mse
						mse_temp[k,which(is.na(mse_temp[k,]))] <- Inf # stop if get a NA fit
					}
					mse_new[inds_to_continue] <- colMeans(mse_temp)
					inds_to_continue <- which(((mse_new - mse_old) < 0))
				}
				options[[i]]$polyorder <- polyorder_curr - 1
				mse[i,] <- mse_old
			} else{
				mse_temp <- matrix(NaN,nrow=folds_choose,ncol=N_expectations)
				for (k in 1:folds_choose){
					res <- zv(integrand, samples, derivatives, log_weight, integrand_logged, obs_estim = obs_estim_choose[[k]], options[[i]]$polyorder, options[[i]]$regul_reg, options[[i]]$alpha_elnet, options[[i]]$nfolds, options[[i]]$apriori, options[[i]]$intercept)
					mse_temp[k,] <- res$mse
				}
				mse[i,] <- colMeans(mse_temp)
				options[[i]]$polyorder <- rep(options[[i]]$polyorder,N_expectations)
			}
		}
		ind_select <- apply(mse,2,which.min)
		
		r_expectation <- r_num_select <- r_mse <-
			r_integrand_logged <- r_polyorder <- r_regul_reg <- r_alpha_elnet<-
			r_nfolds <- r_intercept <- rep(NaN,N_expectations)
		r_apriori <- list()
		
		for (j in unique(ind_select)){
			for (k in unique(options[[j]]$polyorder)){
				inds <- which(ind_select==j)
				inds <- inds[options[[j]]$polyorder[inds]==k]
				temp <- zv(matrix(integrand[,inds],ncol=length(inds)), samples, derivatives, log_weight, integrand_logged, obs_estim = obs_estim, k, options[[j]]$regul_reg, options[[j]]$alpha_elnet, options[[j]]$nfolds, options[[j]]$apriori, options[[j]]$intercept)
				r_expectation[inds] <- temp$expectation
				r_num_select[inds] <- temp$num_select
				r_mse[inds] <- temp$mse
				r_integrand_logged[inds] <- temp$integrand_logged
				r_polyorder[inds] <- temp$polyorder
				r_regul_reg[inds] <- temp$regul_reg
				r_alpha_elnet[inds] <- temp$alpha_elnet
				r_nfolds[inds] <- temp$nfolds
				r_intercept[inds] <- temp$intercept
				for (zz in inds){ r_apriori[[zz]] <- temp$apriori }
			}
		}
		res_final <- list(expectation = r_expectation, num_select = r_num_select, mse = r_mse,
								integrand_logged = r_integrand_logged, obs_estim = obs_estim, polyorder = r_polyorder, regul_reg = r_regul_reg, alpha_elnet = r_alpha_elnet,
								nfolds = r_nfolds, apriori = r_apriori, intercept = r_intercept)
		
	} else{
		ind_select <- rep(1,N_expectations)
		res_final <- zv(integrand, samples, derivatives, log_weight, integrand_logged, obs_estim = obs_estim, options[[1]]$polyorder, options[[1]]$regul_reg, options[[1]]$alpha_elnet, options[[1]]$nfolds, options[[1]]$apriori, options[[1]]$intercept)
	}
	
	return(res_final)
}
