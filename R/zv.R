zv <- function(integrand, samples, derivatives, log_weight, integrand_logged, obs_estim, polyorder, regul_reg, alpha_elnet, nfolds, apriori, intercept) {
	nan_flag <- TRUE
	N <- NROW(integrand)
	d <- NCOL(samples)
	
	if (d==1){
		samples <- array(samples,c(N,d))
		derivatives <- array(derivatives,c(N,d))
	}
	
	if (integrand_logged){
		max_integrand <- max(integrand)
		integrand <- exp(integrand - max_integrand)
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
	y_estim <- integrand[obs_estim]
	
	# Get the data for fitting
	log_weight_fit <- log_weight[obs_fit] - logsumexp(log_weight[obs_fit])
	weight_fit <- exp(log_weight_fit)
	y_fit <- integrand[obs_fit]
	
	# Check if the polynomial order is 0.
	if (polyorder == 0){
		if (!integrand_logged){
			expectation <- sum(weight_fit*y_fit)
		} else {
			expectation <- logsumexp(log(y_fit) + log_weight_fit) + max_integrand
		}
		mse <-  mean((y_estim - expectation)^2)
		return(list(expectation = expectation, num_select = NaN, aic = NaN, bic = NaN, R2 = NaN, adjR2 = NaN, mse = mse,
								integrand_logged = integrand_logged, obs_estim = obs_estim, polyorder = polyorder, regul_reg = regul_reg, alpha_elnet = alpha_elnet, nfolds = nfolds, apriori = apriori, intercept = intercept))
	} else{
		X_estim <- getX(matrix(samples[obs_estim,],nrow=length(obs_estim),ncol=length(apriori)), matrix(derivatives[obs_estim,],nrow=length(obs_estim),ncol=length(apriori)), polyorder)
		X_fit <- getX(samples[obs_fit,], derivatives[obs_fit,], polyorder)
	}
	
	
	
	# may need to rerun without an intercept if there are issues
	while (nan_flag){
		
		if (!intercept){
			y_fit <- y_fit - sum(y_fit*weight_fit)
		}
		
		
		if (regul_reg){
			if (intercept){
				fit <- glmnet(X_fit, y_fit, family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, nlambda = 3)
				mymax <- max(fit$lambda)
				myseq <- c(0,exp(seq(log(mymax),min(-12,log(min(fit$lambda))),length.out = 99)))
				fit <- cv.glmnet(X_fit, y_fit, family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, lambda=myseq, nfolds = nfolds)
			} else{
				fit <- glmnet(X_fit, y_fit, intercept = FALSE, family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, nlambda = 3)
				mymax <- max(fit$lambda)
				myseq <- c(0,exp(seq(log(mymax),min(-12,log(min(fit$lambda))),length.out = 99)))
				fit <- cv.glmnet(X_fit, y_fit, intercept = FALSE, family = 'gaussian', alpha = alpha_elnet, weights = weight_fit, lambda=myseq, nfolds = nfolds)
			}
			
			R2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == fit$lambda.min)]
			adjR2 <- 1 - (1-R2)*(length(y_fit) - 1)/(length(y_fit) - NCOL(X_fit) - 1)
			
			coefs <- coef(fit,s = "lambda.min")
			
			num_select <- unname(fit$nzero[which(fit$lambda == fit$lambda.min)])
			
		} else {
			if (intercept){
				fit <- lm(y_fit ~ X_fit, weights = weight_fit)
			} else {
				fit <- lm(y_fit ~ X_fit - 1, weights = weight_fit)
			}
			coefs <- coef(fit)
			num_select <- NCOL(X_fit)
		}
		
		if (intercept | regul_reg){
			fitteds_estim <- X_estim%*%coefs[2:length(coefs)]
			mse <- mean((y_estim - fitteds_estim - coefs[1])^2)
		} else {
			fitteds_estim <- X_estim%*%coefs
			mse <- mean((y_estim - fitteds_estim - sum(y_fit*weight_fit))^2)
		}
		
		
		if (sum(is.na(fitteds_estim))>0){
			fitteds_estim <- rep(0,NROW(y_estim))
		}
		
		integrand_new <- y_estim - fitteds_estim
		
		if (!integrand_logged){
			expectation <- sum(weight_estim*integrand_new)
		} else if(sum(integrand_new<0) != 0){
			expectation <- log(sum(integrand_new*weight_estim)) + max_integrand
		} else {
			expectation <- logsumexp(log(integrand_new) + log_weight_estim) + max_integrand
		}
		
		if (is.nan(expectation) & intercept){
			intercept <- FALSE
			print("NaN return. Rerunning now with a fixed intercept.")
		}
		else if (is.nan(expectation)){
			nan_flag <- FALSE
			print("Error - intercept was fixed and yet the expectation was NaN")
		} else {
			nan_flag <- FALSE
		}
	}
	
	
	return(list(expectation = expectation, num_select = num_select, mse = mse,
    integrand_logged = integrand_logged, obs_estim = obs_estim, polyorder = polyorder, regul_reg = regul_reg, alpha_elnet = alpha_elnet, nfolds = nfolds, apriori = apriori, intercept = intercept))
}
