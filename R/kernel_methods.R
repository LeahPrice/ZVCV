#' Approximate semi-exact control functionals (aSECF)
#' 
#' This function performs approximate semi-exact control functionals as described in South et al (2020). It uses a nystrom approximation and conjugate gradient to speed up SECF.
#' This is faster than \code{\link{SECF}} for large \eqn{N}. If you would like to choose
#' between different kernels using cross-validation, then you can use \code{\link{aSECF_crossval}}.
#'
#' @param integrand		An \eqn{N} by \eqn{k} matrix of integrands (evaluations of the function of interest)
#' @param samples		An \eqn{N} by \eqn{d} matrix of samples from the target
#' @param derivatives	An \eqn{N} by \eqn{d} matrix of derivatives of the log target with respect to the parameters
#' @param polyorder (optional)		The order of the polynomial to be used in the parametric component, with a default of \eqn{1}. We recommend keeping this value low (e.g. only 1-2).
#' @param steinOrder (optional)	This is the order of the Stein operator. The default is \code{1} in the control functionals paper (Oates et al, 2017) and \code{2} in the semi-exact control functionals paper (South et al, 2020).  The following values are currently available: \code{1} for all kernels and \code{2} for "gaussian", "matern" and "RQ". See below for further details.
#' @param kernel_function (optional)		Choose between "gaussian", "matern", "RQ", "product" or "prodsim". See below for further details.
#' @param sigma (optional)			The tuning parameters of the specified kernel. This involves a single length-scale parameter in "gaussian" and "RQ", a length-scale and a smoothness parameter in "matern" and two parameters in "product" and "prodsim". See below for further details.
#' @param K0 (optional) The kernel matrix. One can specify either this or all of \code{sigma}, \code{steinOrder} and \code{kernel_function}. The former involves pre-computing the kernel matrix using \code{\link{K0_fn}} and is more efficient when using multiple estimators out of \code{\link{CF}}, \code{\link{SECF}} and  \code{\link{aSECF}} or when using the cross-validation functions.
#' @param nystrom_inds (optional) The sample indices to be used in the Nystrom approximation.
#' @param est_inds (optional) The default is to perform estimation and evaluation using the full set of samples. This argument lets you specify that the \code{est_inds} indices are used for estimation only and the remaining samples are used for evaluation only. This can be used to reduce bias from adaption and also to make computation feasible for very large sample sizes (small \code{est_inds} is faster). 
#' @param subset (optional) The subset of parameters to be used in the polynomial. Typically this argument would only be used if the dimension of the problem is very large.
#' @param conjugate_gradient (optional) A flag for whether to perform conjugate gradient to further speed up the nystrom approximation (the default is true).
#' @param reltol (optional) The relative tolerance for choosing when the stop conjugate gradient iterations (the default is 1e-02).
#' using \code{\link{squareNorm}}, as long as the \code{nystrom_inds} are \code{NULL}.
#' @param diagnostics (optional) A flag for whether to return the necessary outputs for plotting or estimating using the fitted model. The default is \code{false} since this requires some additional computation when \code{est_inds} is \code{NULL}.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{expectation}: The estimate(s) of the (\eqn{k}) expectations(s).
#' \item \code{cond_no}: (Only if \code{conjugate_gradient} = \code{TRUE}) The condition number of the matrix being solved using conjugate gradient.
#' \item \code{iter}: (Only if \code{conjugate_gradient} = \code{TRUE}) The number of conjugate gradient iterations
#' \item \code{f_true}: (Only if \code{est_inds} is not \code{NULL}) The integrands for the evaluation set. This should be the same as integrands[setdiff(1:N,est_inds),].
#' \item \code{f_hat}: (Only if \code{est_inds} is not \code{NULL}) The fitted values for the integrands in the evaluation set. This can be used to help assess the performance of the Gaussian process model.
#' \item \code{a}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{a} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{b}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{b} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{ny_inds}: (Only if \code{diagnostics} = \code{TRUE}) The indices of the samples used in the nystrom approximation (this will match nystrom_inds if this argument was not \code{NULL}).
#' } 
#'
#' @inheritSection K0_fn On the choice of \eqn{\sigma}, the kernel and the Stein order
#'
#' @references
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' @author Leah F. South
#' @seealso \code{\link{aSECF_crossval}} for a function to choose between different kernels for this estimator.
aSECF <- function(integrands,samples,derivatives, polyorder = NULL, steinOrder = NULL, kernel_function = NULL, sigma = NULL, K0 = NULL,nystrom_inds = NULL, est_inds = NULL, subset = NULL, conjugate_gradient = TRUE, reltol = 1e-02, diagnostics = FALSE){
	
	N <- NROW(samples)
	N_expectations <- NCOL(integrands)
	
	## convert integrand vector to an Nx1 matrix if necessary
	if (is.null(ncol(integrands))){
		integrands <- matrix(integrands,nrow=N,ncol=1)
	}
	
	if (is.null(ncol(samples))){
		samples <- matrix(samples,nrow=N,ncol=1)
		derivatives <- matrix(derivatives,nrow=N,ncol=1)
	}
	
	d <- NCOL(samples)
	if (!is.null(est_inds)){
		N <- length(est_inds)
	}
	if (!is.null(polyorder)){
		if (choose(d+polyorder,d) >= N){
			stop("The polyorder is too high for this sample size.")
		}
	} else if ((d >= N) && is.null(subset)){
		stop("The dimension is too large for this sample size. Consider increasing the sample size or using the subset argument.")
	} else if (length(subset) >= N){
		stop("The dimension is too large for this sample size. Consider reducing the number of terms in the subset argument.")
	}
	
	if (is.null(est_inds)){
		
		temp <- aSECF_cpp_prep(integrands, samples, derivatives, getX = getX, polyorder, steinOrder, kernel_function, sigma, K0, subset, nystrom_inds, conjugate_gradient)
		
		A <- temp$A
		b <- temp$b
		B2 <- temp$B2
		cond_no <- temp$cond_no
		m0 <- temp$m0
		Q <- NCOL(temp$phi)
		
		if (diagnostics){
			B1 <- temp$B1
			ny_inds <- temp$ny_inds
			
			a <- matrix(NaN,nrow=m0,ncol=N_expectations)
			beta <- matrix(NaN,nrow=Q,ncol=N_expectations)
		}
		
		
		expectation <- rep(NaN, nrow=N_expectations)
		iter <- rep(NaN,N_expectations)
		
		for (i in 1:N_expectations){
			
			if (conjugate_gradient){
				B2_inv <- solve(B2)
				xinit <- c(rep(0,m0),B2_inv[,1]*mean(integrands[,i]))
				ab_tilde <- lsolve.cg(A, b[,i], xinit = xinit, reltol = reltol, preconditioner = diag(ncol(A)), adjsym = TRUE, verbose = FALSE) # maxiter = 10, 
				expectation[i] <- matrix(B2[1,],nrow=1)%*%ab_tilde$x[(m0+1):(m0+Q)]
				iter[i] <- ab_tilde$iter
			} else{
				ab_tilde <- solve(nearPD(A),b[,i]) # Replace this with coordinate descent
				expectation[i] <- ab_tilde[m0+1]
				iter[i] <- NaN
			}
			
			if (diagnostics){
				a[,i] <-  B1%*%ab_tilde$x[1:m0]
				beta[,i] <-  B2%*%ab_tilde$x[(m0+1):(m0+Q)]
			}
		}
		if (conjugate_gradient){
			
			if (diagnostics){
				res <- list(expectation = expectation, cond_no=cond_no, iter = iter, a = a, b = beta, ny_inds = ny_inds)
			} else{
				res <- list(expectation = expectation, cond_no=cond_no, iter = iter)
			}
		} else{
			if (diagnostics){
				res <- list(expectation = expectation, a = a, b = beta, ny_inds = ny_inds)
			} else{
				res <- list(expectation = expectation)
			}
		}
		
		
	} else{
		res <- aSECF_unbiased_cpp_prep(integrands, samples, derivatives, est_inds, getX = getX, aSECF_mse_linsolve = aSECF_mse_linsolve, polyorder, steinOrder, kernel_function, sigma, K0, subset, nystrom_inds, conjugate_gradient, reltol, diagnostics)
	}
	
	return(res)
}


#' Control functionals (CF)
#' 
#' This function performs control functionals as described in Oates et al (2017).
#' To choose between different kernels using cross-validation, use \code{\link{CF_crossval}}.
#'
#' @inheritParams aSECF
#' @param one_in_denom (optional) Whether or not to include a \eqn{1 + } in the denominator of the control functionals estimator, as in equation 2 on p703 of Oates et al (2017). The \eqn{1 +} in the denominator is an arbitrary choice so we set it to zero by default.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{expectation}: The estimate(s) of the (\eqn{k}) expectation(s).
#' \item \code{f_true}: (Only if \code{est_inds} is not \code{NULL}) The integrands for the evaluation set. This should be the same as integrands[setdiff(1:N,est_inds),].
#' \item \code{f_hat}: (Only if \code{est_inds} is not \code{NULL}) The fitted values for the integrands in the evaluation set. This can be used to help assess the performance of the Gaussian process model.
#' \item \code{a}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{a} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + 1*b} for heldout K0 and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b}.
#' \item \code{b}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{b} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + 1*b} for heldout K0 and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b}.
#' \item \code{ksd}: (Only if \code{diagnostics} = \code{TRUE}) An estimated kernel Stein discrepancy based on the fitted model that can be used for diagnostic purposes. See South et al (2020) for further details.
#' \item \code{bound_const}: (Only if \code{diagnostics} = \code{TRUE} and \code{est_inds}=\code{NULL}) This is such that the absolute error for the estimator should be less than \eqn{ksd \times bound_const}.
#' } 
#'
#' @inheritSection K0_fn On the choice of \eqn{\sigma}, the kernel and the Stein order
#'
#' @section Warning:
#' Solving the linear system in CF has \eqn{O(N^3)} complexity and is therefore not suited to large \eqn{N}. Using \eqn{est_inds} will instead have an \eqn{O(N_0^3)} cost in solving the linear system and an \eqn{O((N-N_0)^2)} cost in handling the remaining samples, where \eqn{N_0} is the length of \eqn{est_inds}. This can be much cheaper for large \eqn{N}.
#'
#' @references
#' Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.
#'
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' @author Leah F. South
#' @seealso \code{\link{CF_crossval}} for a function to choose between different kernels for this estimator.
CF <- function(integrands, samples, derivatives, steinOrder = NULL, kernel_function = NULL, sigma = NULL, K0 = NULL, est_inds = NULL, one_in_denom = FALSE, diagnostics = FALSE){
	
	N <- NROW(samples)
	d <- NCOL(samples)
	
	## converting any vectors to matrices as required.
	if (is.null(ncol(integrands))){
		integrands <- matrix(integrands,nrow=N,ncol=1)
	}
	
	if (is.null(ncol(samples))){
		samples <- matrix(samples,nrow=N,ncol=1)
		derivatives <- matrix(derivatives,nrow=N,ncol=1)
	}
	
	if (is.null(est_inds)){
		temp <- CF_cpp(integrands, samples, derivatives, steinOrder, kernel_function, sigma, K0, one_in_denom, diagnostics)
	} else{
		temp <- CF_unbiased_cpp(integrands, samples, derivatives, est_inds, steinOrder, kernel_function, sigma, K0, one_in_denom, diagnostics)
	}
	return (temp)
}

#' Control functionals (CF) with cross-validation
#' 
#' This function chooses between a list of kernel tuning parameters (\code{sigma_list}) or a list of K0 matrices (\code{K0_list}) for
#' the control functionals method described in Oates et al (2017). The latter requires
#' calculating and storing kernel matrices using \code{\link{K0_fn}} but it is more flexible
#' because it can be used to choose the Stein operator order and the kernel function, in addition
#' to its parameters. It is also faster to pre-specify \code{\link{K0_fn}}.
#' For estimation with fixed kernel parameters, use \code{\link{CF}}.
#'
#' @inheritParams aSECF
#' @param sigma_list (optional between this and \code{K0_list})			A list of tuning parameters for the specified kernel. This involves a list of single length-scale parameter in "gaussian" and "RQ", a list of vectors containing length-scale and smoothness parameters in "matern" and a list of vectors of the two parameters in "product" and "prodsim". See below for further details. When \code{sigma_list} is specified and not \code{K0_list}, the \eqn{K0} matrix is computed twice for each selected tuning parameter.
#' @param K0_list (optional between this and \code{sigma_list}) A list of kernel matrices, which can be calculated using \code{\link{K0_fn}}.
#' @param one_in_denom (optional) Whether or not to include a \eqn{1 + } in the denominator of the control functionals estimator, as in equation 2 on p703 of Oates et al (2017). The \eqn{1 +} in the denominator is an arbitrary choice so we set it to zero by default.
#' @param folds (optional) The number of folds for cross-validation. The default is five.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{expectation}: The estimate(s) of the (\eqn{k}) expectation(s).
#' \item \code{mse}: A matrix of the cross-validation mean square prediction errors. The number of columns is the number of tuning options given and the number of rows is \eqn{k}, the number of integrands of interest.
#' \item \code{optinds}: The optimal indices from the list for each expectation.
#' \item \code{f_true}: (Only if \code{est_inds} is not \code{NULL}) The integrands for the evaluation set. This should be the same as integrands[setdiff(1:N,est_inds),].
#' \item \code{f_hat}: (Only if \code{est_inds} is not \code{NULL}) The fitted values for the integrands in the evaluation set. This can be used to help assess the performance of the Gaussian process model.
#' \item \code{a}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{a} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + 1*b} for heldout K0 and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b}.
#' \item \code{b}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{b} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + 1*b} for heldout K0 and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b}.
#' \item \code{ksd}: (Only if \code{diagnostics} = \code{TRUE}) An estimated kernel Stein discrepancy based on the fitted model that can be used for diagnostic purposes. See South et al (2020) for further details.
#' \item \code{bound_const}: (Only if \code{diagnostics} = \code{TRUE} and \code{est_inds}=\code{NULL}) This is such that the absolute error for the estimator should be less than \eqn{ksd \times bound_const}.
#' } 
#'
#' @inheritSection K0_fn On the choice of \eqn{\sigma}, the kernel and the Stein order
#'
#' @section Warning:
#' Solving the linear system in CF has \eqn{O(N^3)} complexity and is therefore not suited to large \eqn{N}. Using \eqn{est_inds} will instead have an \eqn{O(N_0^3)} cost in solving the linear system and an \eqn{O((N-N_0)^2)} cost in handling the remaining samples, where \eqn{N_0} is the length of \eqn{est_inds}. This can be much cheaper for large \eqn{N}.
#'
#' @references
#' Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.
#'
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' @author Leah F. South
#' @seealso \code{\link{CF}} for a function to perform control functionals with fixed kernel specifications.
CF_crossval <- function(integrands, samples, derivatives, steinOrder = NULL, kernel_function = NULL, sigma_list = NULL, K0_list = NULL, est_inds = NULL, one_in_denom = FALSE, folds = NULL, diagnostics = FALSE){
	
	N <- NROW(samples)
	d <- NCOL(samples)
	
	## converting any vectors to matrices as required.
	if (is.null(ncol(integrands))){
		integrands <- matrix(integrands,nrow=N,ncol=1)
	}
	
	if (is.null(ncol(samples))){
		samples <- matrix(samples,nrow=N,ncol=1)
		derivatives <- matrix(derivatives,nrow=N,ncol=1)
	}
	
	temp <- CF_crossval_cpp(integrands, samples, derivatives, steinOrder, kernel_function, sigma_list, K0_list, folds, est_inds, one_in_denom, diagnostics)
	
	return (temp)
}

#' Semi-exact control functionals (SECF)
#' 
#' This function performs semi-exact control functionals as described in South et al (2020).
#' To choose between different kernels using cross-validation, use \code{\link{SECF_crossval}}.
#'
#' @inheritParams aSECF
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{expectation}: The estimate(s) of the (\eqn{k}) expectation(s).
#' \item \code{f_true}: (Only if \code{est_inds} is not \code{NULL}) The integrands for the evaluation set. This should be the same as integrands[setdiff(1:N,est_inds),].
#' \item \code{f_hat}: (Only if \code{est_inds} is not \code{NULL}) The fitted values for the integrands in the evaluation set. This can be used to help assess the performance of the Gaussian process model.
#' \item \code{a}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{a} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{b}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{b} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{ksd}: (Only if \code{diagnostics} = \code{TRUE}) An estimated kernel Stein discrepancy based on the fitted model that can be used for diagnostic purposes. See South et al (2020) for further details.
#' \item \code{bound_const}: (Only if \code{diagnostics} = \code{TRUE} and \code{est_inds}=\code{NULL}) This is such that the absolute error for the estimator should be less than \eqn{ksd \times bound_const}.
#' }
#'
#' @inheritSection K0_fn On the choice of \eqn{\sigma}, the kernel and the Stein order 
#'
#' @section Warning:
#' Solving the linear system in SECF has \eqn{O(N^3+Q^3)} complexity where \eqn{N} is the sample size and \eqn{Q} is the number of terms in the polynomial.
#' Standard SECF is therefore not suited to large \eqn{N}. The method aSECF is designed for larger \eqn{N} and details can be found at \code{\link{aSECF}} and in South et al (2020).
#' An alternative would be to use \eqn{est_inds} which has \eqn{O(N_0^3 + Q^3)} complexity in solving the linear system and \eqn{O((N-N_0)^2)} complexity in
#' handling the remaining samples, where \eqn{N_0} is the length of \eqn{est_inds}. This can be much cheaper for small \eqn{N_0} but the estimation of the
#' Gaussian process model is only done using \eqn{N_0} samples and the evaluation of the integral only uses \eqn{N-N_0} samples.
#'
#' @references
#' Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.
#'
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' @author Leah F. South
#' @seealso \code{\link{SECF_crossval}} for a function to choose between different kernels for this estimator.
SECF <- function(integrands,samples,derivatives, polyorder = NULL, steinOrder = NULL, kernel_function = NULL, sigma = NULL, K0 = NULL, est_inds = NULL,subset = NULL, diagnostics = FALSE){
	
	N <- NROW(samples)
	N_expectations <- NCOL(integrands)
	
	## convert integrand vector to an Nx1 matrix if necessary
	if (is.null(ncol(integrands))){
		integrands <- matrix(integrands,nrow=N,ncol=1)
	}
	
	if (is.null(ncol(samples))){
		samples <- matrix(samples,nrow=N,ncol=1)
		derivatives <- matrix(derivatives,nrow=N,ncol=1)
	}
	
	d <- NCOL(samples)
	if (!is.null(est_inds)){
		N <- length(est_inds)
	}
	if (!is.null(polyorder)){
		if (choose(d+polyorder,d) >= N){
			stop("The polyorder is too high for this sample size.")
		}
	} else if ((d >= N) && is.null(subset)){
		stop("The dimension is too large for this sample size. Consider increasing the sample size or using the subset argument.")
	} else if (length(subset) >= N){
		stop("The dimension is too large for this sample size. Consider reducing the number of terms in the subset argument.")
	}
	
	if (is.null(est_inds)){
		temp <- SECF_cpp(integrands, samples, derivatives, getX = getX, polyorder, steinOrder, kernel_function, sigma, K0, subset, diagnostics)
	} else{
		temp <- SECF_unbiased_cpp(integrands, samples, derivatives, est_inds, getX = getX, polyorder, steinOrder, kernel_function, sigma, K0, subset, diagnostics)
	}
	
	return(temp)
}


#' Semi-exact control functionals (SECF) with cross-validation
#' 
#' This function chooses between a list of kernel tuning parameters (\code{sigma_list}) or a list of K0 matrices (\code{K0_list}) for
#' the semi-exact control functionals method described in South et al (2020). The latter requires
#' calculating and storing kernel matrices using \code{\link{K0_fn}} but it is more flexible
#' because it can be used to choose the Stein operator order and the kernel function, in addition
#' to its parameters. It is also faster to pre-specify \code{\link{K0_fn}}.
#' For estimation with fixed kernel parameters, use \code{\link{SECF}}.
#'
#' @inheritParams aSECF
#' @param sigma_list (optional between this and \code{K0_list})			A list of tuning parameters for the specified kernel. This involves a list of single length-scale parameter in "gaussian" and "RQ", a list of vectors containing length-scale and smoothness parameters in "matern" and a list of vectors of the two parameters in "product" and "prodsim". See below for further details. When \code{sigma_list} is specified and not \code{K0_list}, the \eqn{K0} matrix is computed twice for each selected tuning parameter.
#' @param K0_list (optional between this and \code{sigma_list}) A list of kernel matrices, which can be calculated using \code{\link{K0_fn}}.
#' @param folds (optional) The number of folds for cross-validation. The default is five.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{expectation}: The estimate(s) of the (\eqn{k}) expectation(s).
#' \item \code{mse}: A matrix of the cross-validation mean square prediction errors. The number of columns is the number of tuning options given and the number of rows is \eqn{k}, the number of integrands of interest.
#' \item \code{optinds}: The optimal indices from the list for each expectation.
#' \item \code{f_true}: (Only if \code{est_inds} is not \code{NULL}) The integrands for the evaluation set. This should be the same as integrands[setdiff(1:N,est_inds),].
#' \item \code{f_hat}: (Only if \code{est_inds} is not \code{NULL}) The fitted values for the integrands in the evaluation set. This can be used to help assess the performance of the Gaussian process model.
#' \item \code{a}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{a} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{b}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{b} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{ksd}: (Only if \code{diagnostics} = \code{TRUE}) An estimated kernel Stein discrepancy based on the fitted model that can be used for diagnostic purposes. See South et al (2020) for further details.
#' \item \code{bound_const}: (Only if \code{diagnostics} = \code{TRUE} and \code{est_inds}=\code{NULL}) This is such that the absolute error for the estimator should be less than \eqn{ksd \times bound_const}.
#' } 
#'
#' @inheritSection K0_fn On the choice of \eqn{\sigma}, the kernel and the Stein order
#'
#' @section Warning:
#' Solving the linear system in SECF has \eqn{O(N^3+Q^3)} complexity where \eqn{N} is the sample size and \eqn{Q} is the number of terms in the polynomial.
#' Standard SECF is therefore not suited to large \eqn{N}. The method aSECF is designed for larger \eqn{N} and details can be found at \code{\link{aSECF}} and in South et al (2020).
#' An alternative would be to use \eqn{est_inds} which has \eqn{O(N_0^3 + Q^3)} complexity in solving the linear system and \eqn{O((N-N_0)^2)} complexity in
#' handling the remaining samples, where \eqn{N_0} is the length of \eqn{est_inds}. This can be much cheaper for large \eqn{N} but the estimation of the
#' Gaussian process model is only done using \eqn{N_0} samples and the evaluation of the integral only uses \eqn{N-N_0} samples.
#'
#' @references
#' Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.
#'
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' @author Leah F. South
#' @seealso \code{\link{SECF}} for a function to perform semi-exact control functionals with fixed kernel specifications.
SECF_crossval <- function(integrands,samples,derivatives, polyorder = NULL, steinOrder = NULL, kernel_function = NULL, sigma_list = NULL, K0_list = NULL, est_inds = NULL, subset = NULL, folds = NULL, diagnostics = FALSE){
	
	N <- NROW(samples)
	N_expectations <- NCOL(integrands)
	
	## convert integrand vector to an Nx1 matrix if necessary
	if (is.null(ncol(integrands))){
		integrands <- matrix(integrands,nrow=N,ncol=1)
	}
	
	if (is.null(ncol(samples))){
		samples <- matrix(samples,nrow=N,ncol=1)
		derivatives <- matrix(derivatives,nrow=N,ncol=1)
	}
	
	d <- NCOL(samples)
	if (!is.null(est_inds)){
		N <- length(est_inds)
	}
	if (is.null(folds)){
		N_perfit <- floor(0.8*N)
	} else{
		N_perfit <- floor((folds-1)/folds*N)
	}
	if (!is.null(polyorder)){
		if (choose(d+polyorder,d) >= N_perfit){
			stop("The polyorder is too high for this sample size and number of folds.")
		}
	} else if ((d >= N_perfit) && is.null(subset)){
		stop("The dimension is too large for this sample size and number of folds. Consider increasing the sample size, reducing the number of cross-validation folds or using the subset argument.")
	} else if (length(subset) >= N_perfit){
		stop("The dimension is too large for this sample size and number of folds. Consider reducing the number of cross-validation folds or reducing the number of terms in the subset argument.")
	}
	
	temp <- SECF_crossval_cpp(integrands, samples, derivatives, getX = getX, polyorder, steinOrder, kernel_function, sigma_list, K0_list, subset, folds, est_inds, diagnostics)
	
	return(temp)
}


# An internal function used for performing conjugate gradient in aSECF. This function is used for cross-validation in aSECF and for aSECF when samples are split for estimation and evaluation. The conjugate gradient step is done in R so that an R library can be used.
aSECF_mse_linsolve <- function(integrands,samples,derivatives, polyorder = NULL, steinOrder = NULL, kernel_function = NULL, sigma = NULL, K0 = NULL, subset = NULL, nystrom_inds = NULL, conjugate_gradient = TRUE, reltol = 1e-02){
	
	N <- NROW(samples)
	N_expectations <- NCOL(integrands)
	
	## convert integrand vector to an Nx1 matrix if necessary
	if (is.null(ncol(integrands))){
		integrands <- matrix(integrands,nrow=N,ncol=1)
	}
	
	if (is.null(ncol(samples))){
		samples <- matrix(samples,nrow=N,ncol=1)
		derivatives <- matrix(derivatives,nrow=N,ncol=1)
	}
	
	temp <- aSECF_cpp_prep(integrands, samples, derivatives, getX = getX, polyorder, steinOrder, kernel_function, sigma, K0, subset, nystrom_inds, conjugate_gradient)
	
	A <- temp$A
	b <- temp$b
	B1 <- temp$B1
	B2 <- temp$B2
	cond_no <- temp$cond_no
	m0 <- temp$m0
	
	Q <- NCOL(temp$phi)
	
	expectation <- rep(NaN, nrow=N_expectations)
	
	ab_tilde <- list() # approximate solution, including weights for kernel part and then beta.
	
	for (i in 1:N_expectations){
		ab_tilde[[i]] <- list()
		if (conjugate_gradient){
			B2_inv <- solve(B2)
			xinit <- c(rep(0,m0),B2_inv[,1]*mean(integrands[,i]))
			# ab_tilde[[i]] <- lsolve.cg(A, b[,i], xinit = xinit, reltol = reltol, preconditioner = diag(ncol(A)), adjsym = TRUE, verbose = FALSE)$x #, maxiter = 10
			temp <- lsolve.cg(A, b[,i], xinit = xinit, reltol = reltol, preconditioner = diag(ncol(A)), adjsym = TRUE, verbose = FALSE) #, maxiter = 10
			temp$x[1:m0] <-  B1%*%temp$x[1:m0]
			temp$x[(m0+1):(m0+Q)] <-  B2%*%temp$x[(m0+1):(m0+Q)]
			ab_tilde[[i]]$sol <- temp$x #, maxiter = 10
			ab_tilde[[i]]$iter <- temp$iter
			ab_tilde[[i]]$cond_no <- cond_no
		} else{
			ab_tilde[[i]]$sol <- solve(nearPD(A),b)
		}
	}
	
	return(ab_tilde)
}



#' Approximate semi-exact control functionals (aSECF) with cross-validation
#' 
#' This function chooses between a list of kernel tuning parameters (\code{sigma_list}) or a list of K0 matrices (\code{K0_list}) for
#' the approximate semi-exact control functionals method described in South et al (2020). The latter requires
#' calculating and storing kernel matrices using \code{\link{K0_fn}} but it is more flexible
#' because it can be used to choose the Stein operator order and the kernel function, in addition
#' to its parameters. It is also faster to pre-specify \code{\link{K0_fn}}.
#' For estimation with fixed kernel parameters, use \code{\link{aSECF}}.
#'
#' @inheritParams aSECF
#' @param sigma_list (optional between this and \code{K0_list})			A list of tuning parameters for the specified kernel. This involves a list of single length-scale parameter in "gaussian" and "RQ", a list of vectors containing length-scale and smoothness parameters in "matern" and a list of vectors of the two parameters in "product" and "prodsim". See below for further details. When \code{sigma_list} is specified and not \code{K0_list}, the \eqn{K0} matrix is computed twice for each selected tuning parameter.
#' @param K0_list (optional between this and \code{sigma_list}) A list of kernel matrices, which can be calculated using \code{\link{K0_fn}}.
#' @param num_nystrom (optional) The number of samples to use in the Nystrom approximation, with a default of ceiling(sqrt(N)). The nystrom indices cannot be passed in here because of the way the cross-validation has been set up.
#' @param folds (optional) The number of folds for cross-validation. The default is five.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{expectation}: The estimate(s) of the (\eqn{k}) expectations(s).
#' \item \code{mse}: A matrix of the cross-validation mean square prediction errors. The number of columns is the number of tuning options given and the number of rows is \eqn{k}, the number of integrands of interest.
#' \item \code{optinds}: The optimal indices from the list for each expectation.
#' \item \code{cond_no}: (Only if \code{conjugate_gradient} = \code{TRUE}) The condition number of the matrix being solved using conjugate gradient.
#' \item \code{iter}: (Only if \code{conjugate_gradient} = \code{TRUE}) The number of conjugate gradient iterations
#' \item \code{f_true}: (Only if \code{est_inds} is not \code{NULL}) The integrands for the evaluation set. This should be the same as integrands[setdiff(1:N,est_inds),].
#' \item \code{f_hat}: (Only if \code{est_inds} is not \code{NULL}) The fitted values for the integrands in the evaluation set. This can be used to help assess the performance of the Gaussian process model.
#' \item \code{a}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{a} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{b}: (Only if \code{diagnostics} = \code{TRUE}) The value of \eqn{b} as described in South et al (2020), where predictions are of the form \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices and estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]}.
#' \item \code{ny_inds}: (Only if \code{diagnostics} = \code{TRUE}) The indices of the samples used in the nystrom approximation (this will match nystrom_inds if this argument was not \code{NULL}).
#' } 
#'
#' @inheritSection K0_fn On the choice of \eqn{\sigma}, the kernel and the Stein order
#'
#' @references
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' @author Leah F. South
#' @seealso \code{\link{aSECF_crossval}} for a function to choose between different kernels for this estimator.
aSECF_crossval <- function(integrands,samples,derivatives, polyorder = NULL, steinOrder = NULL, kernel_function = NULL, sigma_list = NULL, est_inds = NULL, subset = NULL, num_nystrom = NULL, conjugate_gradient = TRUE, reltol = 1e-02, folds = NULL, diagnostics = FALSE){
	
	N <- NROW(samples)
	N_expectations <- NCOL(integrands)
	
	## convert integrand vector to an Nx1 matrix if necessary
	if (is.null(ncol(integrands))){
		integrands <- matrix(integrands,nrow=N,ncol=1)
	}
	
	if (is.null(ncol(samples))){
		samples <- matrix(samples,nrow=N,ncol=1)
		derivatives <- matrix(derivatives,nrow=N,ncol=1)
	}
	
	d <- NCOL(samples)
	if (!is.null(est_inds)){
		N <- length(est_inds)
	}
	if (is.null(folds)){
		N_perfit <- floor(0.8*N)
	} else{
		N_perfit <- floor((folds-1)/folds*N)
	}
	if (!is.null(polyorder)){
		if (choose(d+polyorder,d) >= N_perfit){
			stop("The polyorder is too high for this sample size and number of folds.")
		}
	} else if ((d >= N_perfit) && is.null(subset)){
		stop("The dimension is too large for this sample size and number of folds. Consider increasing the sample ize, reducing the number of cross-validation folds or using the subset argument.")
	} else if (length(subset) >= N_perfit){
		stop("The dimension is too large for this sample size and number of folds. Consider reducing the number of cross-validation folds or reducing the number of terms in the subset argument.")
	}
	
	if (is.null(num_nystrom)){
		num_nystrom <- ceiling(sqrt(N))
	}
	
	mse <- aSECF_crossval_cpp(integrands, samples, derivatives, getX = getX, aSECF_mse_linsolve = aSECF_mse_linsolve, num_nystrom = num_nystrom, polyorder, steinOrder, kernel_function, sigma_list, subset, folds, conjugate_gradient, reltol = reltol, est_inds = est_inds)
	
	opt_indices <- apply(mse,1,which.min)
	
	expectation <- rep(NaN,N_expectations)
	if (conjugate_gradient){
		cond_no <- iter <- rep(NaN,N_expectations)
	}
	if (!is.null(est_inds)){
		f_true <- f_hat <- matrix(NaN,nrow=N-length(est_inds),ncol=N_expectations)
	}
	if (diagnostics){
		a <- matrix(NaN,nrow=ceiling(sqrt(N)),ncol=N_expectations)
		ny_inds <- matrix(NaN,nrow=num_nystrom,ncol=N_expectations)
		
	}
	for (j in unique(opt_indices)){
		inds <- which(opt_indices==j)
		nystrom_inds <- sample(1:N,num_nystrom) # Randomly select the final nystrom_inds
		temp <- aSECF(matrix(integrands[,inds],ncol=length(inds)), samples, derivatives, polyorder, steinOrder, kernel_function, sigma_list[[j]], NULL, nystrom_inds, est_inds, subset, conjugate_gradient, reltol, diagnostics)
		
		if(diagnostics & (j==(unique(opt_indices))[1])){
			b <- matrix(NaN,nrow=length(temp$b),ncol=N_expectations)
		}
		expectation[inds] <- temp$expectation
		if (!is.null(est_inds)){
			f_true[,inds] <- temp$f_true
			f_hat[,inds] <- temp$f_hat
		}
		if (conjugate_gradient){
			cond_no[inds] <- temp$cond_no
			iter[inds] <- temp$iter
		}
		if (diagnostics){
			a[,inds] <- temp$a
			b[,inds] <- temp$b
			if (length(inds)==1){
				ny_inds[,inds] <- nystrom_inds
			} else{
				for (zz in 1:length(inds)){
					ny_inds[,inds[zz]] <- nystrom_inds
				}
			}
		}
	}
	
	if (!diagnostics){
		if (conjugate_gradient & !is.null(est_inds)){
			return(list(expectation=expectation,f_true=f_true,f_hat=f_hat,iter=iter,cond_no=cond_no,mse=mse,optinds=opt_indices))
		} else if (conjugate_gradient){
			return(list(expectation=expectation,iter=iter,cond_no=cond_no,mse=mse,optinds=opt_indices))
		} else{
			return(list(expectation=expectation,mse=mse,optinds=opt_indices))
		}
	} else {
		if (conjugate_gradient & !is.null(est_inds)){
			return(list(expectation=expectation,f_true=f_true,f_hat=f_hat,iter=iter,cond_no=cond_no,mse=mse,optinds=opt_indices,
									a=a,b=b,ny_inds=ny_inds))
		} else if (conjugate_gradient){
			return(list(expectation=expectation,iter=iter,cond_no=cond_no,mse=mse,optinds=opt_indices,
									a=a,b=b,ny_inds=ny_inds))
		} else{
			return(list(expectation=expectation,mse=mse,optinds=opt_indices,
									a=a,b=b,ny_inds=ny_inds))
		}
	}
	
	
}

#' Phi matrix calculation
#' 
#' This function calculates the \eqn{\Phi} matrix, which is a second order Stein operator applied
#' to a polynomial. See South et al (2020) for further details. This function is not required for
#' estimation but may be useful when evaluation samples are not initially available since
#' estimators using heldout samples are of the form \eqn{mean(f - f_hat) + b[1]} where \eqn{f_hat = K0*a + Phi*b} for heldout K0 and Phi matrices.
#' 
#' @inheritParams aSECF
#'
#' @return An \eqn{N} by \eqn{Q} matrix (where Q is determined by the polynomial order and the subset).
#'
#' @references
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' @author Leah F. South
Phi_fn <- function(samples,derivatives,polyorder=NULL,subset=NULL){
	return(Phi_fn_cpp(samples,derivatives, getX = getX, polyorder,subset))
}
