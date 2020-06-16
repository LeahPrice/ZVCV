#' Zero-Variance Control Variates
#' 
#' @description 
#' This package can be used to perform post-hoc variance reduction of Monte Carlo estimators when the derivatives of the log target are available.
#' The main functionality is available through the following functions. 
#' All of these use a set of \eqn{N} \eqn{d}-dimensional samples along with the associated derivatives of the log target. 
#' You can evaluate posterior expectations of \eqn{k} functions.
#'
#' \itemize{
#'      \item \code{\link{zvcv}}: For estimating expectations using (regularised) zero-variance control variates (ZV-CV, Mira et al, 2013; South et al, 2018).
#'      This function can also be used to choose between various versions of ZV-CV using cross-validation.
#'      \item \code{\link{CF}}: For estimating expectations using control functionals (CF, Oates et al, 2017). 
#'      \item \code{\link{SECF}}: For estimating expectations using semi-exact control functionals (SECF, South et al, 2020).
#'      \item \code{\link{aSECF}}: For estimating expectations using approximate semi-exact control functionals (aSECF, South et al, 2020). 
#'      \item \code{\link{CF_crossval}}: CF with cross-validation tuning.
#'      \item \code{\link{SECF_crossval}}: SECF with cross-validation tuning.
#'      \item \code{\link{aSECF_crossval}}: aSECF with cross-validation tuning.
#'}
#'
#' ZV-CV is exact for polynomials of order at most \code{polyorder} under Gaussian targets and is fast for large \eqn{N} (although
#' setting a limit on \code{polyorder} through \code{polyorder_max} is recommended for large \eqn{N}).
#' CF is a non-parametric approach that offers better than the standard Monte Carlo convergence rates. 
#' SECF has both a parametric and a non-parametric component and it offers the advantages of both for an additional computational cost. The cost of
#' SECF is reduced in aSECF using nystrom approximations and conjugate gradient.
#'
#' @section Helper functions:
#' \itemize{
#'      \item \code{\link{getX}}: Calculates the design matrix for ZV-CV (without the column of 1's for the intercept)
#'      \item \code{\link{medianTune}}: Calculates the median heuristic for use in e.g. the Gaussian, Matern and rational quadratic kernels. Using the median heuristic is an alternative to cross-validation.
#'      \item \code{\link{K0_fn}}: Calculates the \eqn{K_0} matrix. The output of this function can be used as an argument to \code{\link{CF}}, \code{\link{CF_crossval}},
#'      \code{\link{SECF}}, \code{\link{SECF_crossval}}, \code{\link{aSECF}} and \code{\link{aSECF_crossval}}.
#'      The kernel matrix is automatically computed in all of the above methods, but it is faster to calculate
#'      in advance when using more than one of the above functions and when using any of the crossval functions.
#'      \item \code{\link{Phi_fn}}: Calculates the Phi matrix for SECF and aSECF (similar to \code{getX} but with different arguments and it includes the column of 1's)
#'      \item \code{\link{squareNorm}}: Gets the matrix of square norms which is needed for all kernels.
#'      Calculating this can help to save time if you are also interested in calculating the median heuristic, handling multiple tuning parameters or trying other kernels.
#'      \item \code{\link{nearPD}}: Finds the nearest symmetric positive definite matrix to the given matrix, for handling numerical issues.
#'      \item \code{\link{logsumexp}}: Performs stable computation of the log sum of exponential (useful when handling the sum of weights)
#'      }
#' 
#' @section Evidence estimation:
#' The following functions are used to estimate the evidence (the normalisiing constant of the posterior) as described in South et al (2018). They are relevant when
#' sequential Monte Carlo with an annealing schedule has been used to collect the samples, and therefore are not of interest to those who are interested in
#' variance reduction based on vanilla MCMC.
#' \itemize{
#'      \item \code{\link{evidence_CTI}} and \code{\link{evidence_CTI_CF}}: Functions to estimate the evidence using thermodynamic integration (TI) with ZV-CV and CF, respectively
#'      \item \code{\link{evidence_SMC}} and \code{\link{evidence_SMC_CF}}: Function to estimate the evidence using the SMC evidence identity with ZV-CV and CF, respectively.
#'}
#' 
#' The function \code{\link{Expand_Temperatures}} can be used to adjust the temperature schedule so that it is more (or less) strict than the original schedule of \eqn{T} temperatures.
#'
#'@examples
#' # A real data example using ZV-CV is available at \link{VDP}. This involves estimating posterior expectations and the evidence from SMC samples.
#' 
#' # The remainder of this section is duplicating (albeit with a different random seed) Figure 2a of South et al. (2020).
#' 
#' N_repeats <- 10 # For speed, the actual code uses 100 
#' N_all <- c(25,50,100) # For speed, the actual code uses c(10,25,50,100,250,500,1000) 
#' sigma_list <- list(10^(-1.5),10^(-1),10^(-0.5),1,10^(0.5),10)
#' folds <- 5
#' d <- 4
#' 
#' integrand_fn <- function(x){ return (1 + x[,2] + 0.1*x[,1]*x[,2]*x[,3] + sin(x[,1])*exp(-(x[,2]*x[,3])^2)) }
#' 
#' results <- data.frame()
#' for (N in N_all){
#' 
#'   # identify the largest polynomial order that can be fit without regularisation for auto ZV-CV
#'   max_r <- 0
#'   while (choose(d + max_r + 1,d)<(0.5*N)){ # 0.5*N because of the 2-fold cross-validation
#'   	max_r <- max_r + 1
#'   }
#' 
#'   MC <- ZV1 <- ZV2 <- ZVchoose <- CF <- SECF1 <- aSECF1 <- SECF2 <- aSECF2 <- rep(NaN,N_repeats)
#'   CF_medHeur <- SECF1_medHeur <- aSECF1_medHeur <- SECF2_medHeur <- aSECF2_medHeur <- rep(NaN,N_repeats)
#'   for (i in 1:N_repeats){     
#'     x <- matrix(rnorm(N*d),ncol=d)
#'     u <- -x
#'     f <- integrand_fn(x)
#'     
#'     MC[i] <- mean(f)
#'     ZV1[i] <- zvcv(f,x,u,options=list(polyorder=1,regul_reg=FALSE))$expectation
#'     if (N > choose(d+2,d)){  # Checking if the sample size is large enough to accommodation a second order polynomial
#'       ZV2[i] <- zvcv(f,x,u,options=list(polyorder=2,regul_reg=FALSE))$expectation
#'     }
#'     myopts <- list(list(polyorder=Inf,regul_reg=FALSE,polyorder_max=max_r),list(polyorder=Inf,nfolds=4))
#'     ZVchoose[i] <- zvcv(f,x,u,options=myopts,folds = 2)$expectation
#'     
#'     # Calculating the kernel matrix in advance for CF and SECF
#'     K0_list <- list()
#'     for (j in 1:length(sigma_list)){
#'       K0_list[[j]] <- K0_fn(x,u,sigma_list[[j]],steinOrder=2,kernel_function="RQ")
#'     }
#'     
#'     CF[i] <- CF_crossval(f,x,u,K0_list=K0_list,folds = 2)$expectation
#'     SECF1[i] <- SECF_crossval(f,x,u,K0_list=K0_list,folds = 2)$expectation
#'     aSECF1[i] <- aSECF_crossval(f,x,u,steinOrder=2,kernel_function="RQ",sigma_list=sigma_list,reltol=1e-05,folds = 2)$expectation
#'     if (max_r>=2){
#'       SECF2[i] <- SECF_crossval(f,x,u,polyorder=2,K0_list=K0_list,folds = 2)$expectation
#'       aSECF2[i] <- aSECF_crossval(f,x,u,polyorder=2,steinOrder=2,kernel_function="RQ",sigma_list=sigma_list,reltol=1e-05,folds = 2)$expectation
#'     }
#'
#'     medHeur <- medianTune(x)
#'     K0_medHeur <- K0_fn(x,u,medHeur,steinOrder=2,kernel_function="RQ")
#'     CF_medHeur[i] <- CF(f,x,u,K0=K0_medHeur)$expectation
#'     SECF1_medHeur[i] <- SECF(f,x,u,K0=K0_medHeur)$expectation
#'     aSECF1_medHeur[i] <- aSECF(f,x,u,steinOrder=2,kernel_function="RQ",sigma=medHeur,reltol=1e-05)$expectation
#'     if (max_r>=2){
#'       SECF2_medHeur[i] <- SECF(f,x,u,polyorder=2,K0=K0_medHeur)$expectation
#'       aSECF2_medHeur[i] <- aSECF(f,x,u,polyorder=2,steinOrder=2,kernel_function="RQ",sigma=medHeur,reltol=1e-05)$expectation
#'     }
#'     
#'     print(sprintf("--%d",i))
#'   }
#'   # Adding the results to a data frame
#'   MSE_crude <- mean((MC - 1)^2)
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = 1, type = "MC")) 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((ZV1 - 1)^2), type = "ZV")) 
#'   results <- rbind(results,data.frame(N=N, order = "2", efficiency = MSE_crude/mean((ZV2 - 1)^2), type = "ZV")) 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((ZVchoose - 1)^2), type = "ZVchoose")) 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((CF - 1)^2), type = "CF")) 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((SECF1 - 1)^2), type = "SECF")) 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((aSECF1 - 1)^2), type = "aSECF")) 
#'   if ((0.5*N) > choose(d+2,d)){
#'     results <- rbind(results,data.frame(N=N, order = "2", efficiency = MSE_crude/mean((SECF2 - 1)^2), type = "SECF")) 
#'     results <- rbind(results,data.frame(N=N, order = "2", efficiency = MSE_crude/mean((aSECF2 - 1)^2), type = "aSECF")) 
#'   }
#' 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((CF_medHeur - 1)^2), type = "CF_medHeur")) 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((SECF1_medHeur - 1)^2), type = "SECF_medHeur")) 
#'   results <- rbind(results,data.frame(N=N, order = "1 or NA", efficiency = MSE_crude/mean((aSECF1_medHeur - 1)^2), type = "aSECF_medHeur")) 
#'   if ((0.5*N) > choose(d+2,d)){
#'     results <- rbind(results,data.frame(N=N, order = "2", efficiency = MSE_crude/mean((SECF2_medHeur - 1)^2), type = "SECF_medHeur")) 
#'     results <- rbind(results,data.frame(N=N, order = "2", efficiency = MSE_crude/mean((aSECF2_medHeur - 1)^2), type = "aSECF_medHeur")) 
#'   }
#'   print(N)
#' }
#' 
#' 
#' # Plotting results where cross-validation is used for kernel methods
#' require(ggplot2)
#' require(ggthemes)
#' a <- ggplot(data=subset(results,!(type %in% c("CF_medHeur","SECF_medHeur","aSECF_medHeur","SECF_medHeur","aSECF_medHeur"))),
#'             aes(x=N,y=efficiency,col=type,linetype=order)) + scale_color_pander() + 
#'   ggtitle("") + geom_line(size=1.5) + scale_x_log10() + scale_y_log10() + 
#'   annotation_logticks(base=10) + labs(x="N",y="Efficiency",color="Method",linetype="Polynomial Order") + theme_minimal(base_size = 15) +
#'   theme(legend.key.size = unit(0.5, "cm"),legend.key.width =  unit(1, "cm")) +
#'   guides(linetype = guide_legend(override.aes = list(size=1),title.position = "top"), color = guide_legend(override.aes = list(size=1),title.position = "top"))
#' print(a)
#' 
#' 
#' # Plotting results where the median heuristic is used for kernel methods
#' b <- ggplot(data=subset(results,!(type %in% c("CF","SECF","aSECF","SECF","aSECF"))),
#'             aes(x=N,y=efficiency,col=type,linetype=order)) + scale_color_pander() + 
#'   ggtitle("") + geom_line(size=1.5) + scale_x_log10() + scale_y_log10() + 
#'   annotation_logticks(base=10) + labs(x="N",y="Efficiency",color="Method",linetype="Polynomial Order") + theme_minimal(base_size = 15) +
#'   theme(legend.key.size = unit(0.5, "cm"),legend.key.width =  unit(1, "cm")) +
#'   guides(linetype = guide_legend(override.aes = list(size=1),title.position = "top"), color = guide_legend(override.aes = list(size=1),title.position = "top"))
#' print(b)
#' 
#'
#' @references
#' Mira, A., Solgi, R., & Imparato, D. (2013). Zero variance Markov chain Monte Carlo for Bayesian estimators. Statistics and Computing, 23(5), 653-662.
#' 
#' South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
#'
#' South, L. F., Oates, C. J., Mira, A., & Drovandi, C. (2018). Regularised zero-variance control variates for high-dimensional variance reduction. \url{https://arxiv.org/abs/1811.05073}
#'
#' @author Leah F. South
#' @name ZVCV_package
"_PACKAGE"
#> [1] "_PACKAGE"