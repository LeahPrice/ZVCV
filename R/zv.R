zv <- function(integrand, samples, derivatives, log_weight, integrand_logged, obs_estim, polyorder, regul_reg, alpha_elnet, nfolds, apriori, intercept) {
	nan_flag <- TRUE
	N <- NROW(integrand)
	N_expectations <- NCOL(integrand)
	d <- NCOL(samples)
	
	if (d==1){
		samples <- array(samples,c(N,d))
		derivatives <- array(derivatives,c(N,d))
	}
	
	if (is.null(ncol(integrand))){
		integrand <- matrix(integrand,ncol=1)
	}
	
	if (integrand_logged){
		max_integrand <- apply(integrand,2,max)
		integrand <- exp(test - matrix(1,nrow=N)%*%max_integrand)
	}
	
	# Setting default weight
	if (missing(log_weight)){
		log_weight <- -log(N)*rep(1,N)
	}
	
	# Here onwards depends on options
	
	if (is.null(obs_estim)){
		obs_estim <- 1:N
		obs_fit <- 1:N
	} else {
		obs_fit <- setdiff(1:N,obs_estim) 
	}
	
	# Using only relevant variables
	samples <- array(samples[,apriori],c(N,length(apriori)))
	derivatives <- array(derivatives[,apriori],c(N,length(apriori)))
	
	# Splitting the data for fitting and estimation
	log_weight_estim <- log_weight[obs_estim] - logsumexp(log_weight[obs_estim])
	weight_estim <- exp(log_weight_estim)
	y_estim <- integrand[obs_estim,]
	
	# Get the data for fitting
	log_weight_fit <- log_weight[obs_fit] - logsumexp(log_weight[obs_fit])
	weight_fit <- exp(log_weight_fit)
	y_fit <- integrand[obs_fit,]
	
	if (N_expectations==1){
		y_estim <- matrix(y_estim,ncol=1)
		y_fit <- matrix(y_fit,ncol=1)
	}
	
	# Check if the polynomial order is 0.
	if (polyorder == 0){
		if (!integrand_logged){
			expectation <- matrix(weight_fit,nrow=1)%*%y_fit
		} else {
			expectation <- apply(y_fit,2,function(x) logsumexp(log(x) + log_weight_fit)) + max_integrand
		}
		mse <-  colMeans((y_estim - matrix(1,nrow=NROW(y_estim))%*%expectation)^2)
		return(list(expectation = expectation, num_select = NaN, aic = NaN, bic = NaN, R2 = NaN, adjR2 = NaN, mse = mse,
								integrand_logged = integrand_logged, obs_estim = obs_estim, polyorder = polyorder, regul_reg = regul_reg, alpha_elnet = alpha_elnet, nfolds = nfolds, apriori = apriori, intercept = intercept))
	} else{
		X_estim <- getX(matrix(samples[obs_estim,],nrow=length(obs_estim),ncol=length(apriori)), matrix(derivatives[obs_estim,],nrow=length(obs_estim),ncol=length(apriori)), polyorder)
		X_fit <- getX(samples[obs_fit,], derivatives[obs_fit,], polyorder)
	}
	
	
	
	# may need to rerun without an intercept if there are issues
	while (nan_flag){
		
		if (!intercept){
			y_fit <- y_fit - matrix(1,nrow=NROW(y_fit))%*%weight_fit%*%y_fit
		}
		
		
		if (regul_reg){
			num_select <- rep(NaN,NCOL(y_fit))
			if (intercept){
				coefs <- matrix(NaN,nrow=NCOL(X_fit)+1,ncol=NCOL(y_fit))
				for (j in 1:NCOL(y_fit)){
					fit <- glmnet(X_fit, y_fit[,j], family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, nlambda = 3)
					mymax <- max(fit$lambda)
					myseq <- c(0,exp(seq(log(mymax),min(-12,log(min(fit$lambda))),length.out = 99)))
					fit <- cv.glmnet(X_fit, y_fit[,j], family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, lambda=myseq, nfolds = nfolds)
					
					coefs[,j] <- coef(fit,s = "lambda.min")[,1]
					num_select[j] <- unname(fit$nzero[which(fit$lambda == fit$lambda.min)])
				}
			} else{
				coefs <- matrix(NaN,nrow=NCOL(X_fit),ncol=NCOL(y_fit))
				for (j in 1:NCOL(y_fit)){
					fit <- glmnet(X_fit, y_fit[,j], intercept = FALSE, family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, nlambda = 3)
					mymax <- max(fit$lambda)
					myseq <- c(0,exp(seq(log(mymax),min(-12,log(min(fit$lambda))),length.out = 99)))
					fit <- cv.glmnet(X_fit, y_fit[,j], intercept = FALSE, family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, lambda=myseq, nfolds = nfolds)
					
					coefs[,j] <- coef(fit,s = "lambda.min")
					num_select[j] <- unname(fit$nzero[which(fit$lambda == fit$lambda.min)])
				}
			}
			
			
			
		} else {
			if (intercept){
				fit <- lm(y_fit ~ X_fit, weights = weight_fit)
			} else {
				fit <- lm(y_fit ~ X_fit - 1, weights = weight_fit)
			}
			coefs <- coef(fit) # columns are different integrands
			num_select <- rep(NCOL(X_fit),NCOL(y_fit))
		}
		
		if (is.null(ncol(coefs))){
			coefs <- matrix(coefs,ncol=1)
		}
		
		if (intercept | regul_reg){
			fitteds_estim <- X_estim%*%coefs[2:NROW(coefs),]
			mse <- colMeans((y_estim - fitteds_estim - matrix(1,nrow=NROW(y_estim))%*%coefs[1,])^2)
		} else {
			fitteds_estim <- X_estim%*%coefs
			mse <- colMeans((y_estim - fitteds_estim - matrix(1,nrow=NROW(y_fit))%*%weight_fit%*%y_fit)^2)
		}
		
		
		if (any(is.na(fitteds_estim))){
			fitteds_estim <- matrix(0,NROW(y_estim),NCOL(y_estim))
		}
		
		integrand_new <- y_estim - fitteds_estim
		
		if (!integrand_logged){
			expectation <- weight_estim%*%integrand_new
		} else if(any(integrand_new<0)){
			expectation <- log(weight_estim%*%integrand_new) + max_integrand
		} else {
			expectation <- apply(integrand_new,2,function(x) logsumexp(log(x) + log_weight_estim)) + max_integrand
		}
		
		if (any(is.nan(expectation)) & intercept){
			intercept <- FALSE
			print("NaN return. Rerunning now with a fixed intercept.")
		}
		else if (any(is.nan(expectation))){
			nan_flag <- FALSE
			print("Error - intercept was fixed and yet the expectation was NaN")
		} else {
			nan_flag <- FALSE
		}
	}
	
	
	return(list(expectation = expectation, num_select = num_select, mse = mse,
							integrand_logged = integrand_logged, obs_estim = obs_estim, polyorder = polyorder, regul_reg = regul_reg, alpha_elnet = alpha_elnet, nfolds = nfolds, apriori = apriori, intercept = intercept))
}
