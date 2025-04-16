zv <- function(integrand, samples, derivatives, log_weights, integrand_logged, est_inds, polyorder, regul_reg, alpha_elnet, nfolds, apriori, intercept) {
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
		integrand <- exp(integrand - matrix(1,nrow=N)%*%max_integrand)
	}
	
	# Setting default weight
	if (missing(log_weights)){
		log_weights <- -log(N)*rep(1,N)
	}
	
	# Here onwards depends on options
	
	if (is.null(est_inds)){
		est_inds <- 1:N
		eval_inds <- 1:N
	} else {
		eval_inds <- setdiff(1:N,est_inds) 
	}
	
	# Using only relevant variables
	samples <- array(samples[,apriori],c(N,length(apriori)))
	derivatives <- array(derivatives[,apriori],c(N,length(apriori)))
	
	# Splitting the data for fitting and estimation
	log_weights_eval <- log_weights[eval_inds] - logsumexp(log_weights[eval_inds])
	weight_eval <- exp(log_weights_eval)
	y_eval <- integrand[eval_inds,]
	
	# Get the data for fitting
	log_weights_estim <- log_weights[est_inds] - logsumexp(log_weights[est_inds])
	weight_estim <- exp(log_weights_estim)
	y_estim <- integrand[est_inds,]
	
	if (N_expectations==1){
		y_eval <- matrix(y_eval,ncol=1)
		y_estim <- matrix(y_estim,ncol=1)
	}
	
	# Check if the polynomial order is 0.
	if (polyorder == 0){
		if (!integrand_logged){
			expectation <- matrix(weight_estim,nrow=1)%*%y_estim
		} else {
			expectation <- apply(y_estim,2,function(x) logsumexp(log(x) + log_weights_estim)) + max_integrand
		}
		mse <-  colMeans((y_eval - matrix(1,nrow=NROW(y_eval))%*%expectation)^2)
		return(list(expectation = expectation, num_select = NaN, aic = NaN, bic = NaN, R2 = NaN, adjR2 = NaN, mse = mse,
								integrand_logged = integrand_logged, est_inds = est_inds, polyorder = polyorder, regul_reg = regul_reg, alpha_elnet = alpha_elnet, nfolds = nfolds, apriori = apriori, intercept = intercept))
	} else{
		X_eval <- getX(matrix(samples[eval_inds,],nrow=length(eval_inds),ncol=length(apriori)), matrix(derivatives[eval_inds,],nrow=length(eval_inds),ncol=length(apriori)), polyorder)
		X_estim <- getX(samples[est_inds,], derivatives[est_inds,], polyorder)
	}
	
	
	
	# may need to rerun without an intercept if there are issues
	while (nan_flag){
		
		if (!intercept){
			y_estim <- y_estim - matrix(1,nrow=NROW(y_estim))%*%weight_estim%*%y_estim
		}
		
	  if (regul_reg){
			num_select <- rep(NaN,NCOL(y_estim))
			if (intercept){
				coefs <- matrix(NaN,nrow=NCOL(X_estim)+1,ncol=NCOL(y_estim))
				for (j in 1:NCOL(y_estim)){
					fit <- glmnet(X_estim, y_estim[,j], family = 'gaussian', alpha = alpha_elnet, weights = weight_estim, nlambda = 3)
					mymax <- max(fit$lambda)
					myseq <- c(0,exp(seq(log(mymax),min(-12,log(min(fit$lambda))),length.out = 99)))
					fit <- cv.glmnet(X_estim, y_estim[,j], family = 'gaussian', alpha = alpha_elnet, weights = weight_estim, lambda=myseq, nfolds = nfolds)
					
					coefs[,j] <- coef(fit,s = "lambda.min")[,1]
					num_select[j] <- unname(fit$nzero[which(fit$lambda == fit$lambda.min)])
				}
			} else{
				coefs <- matrix(NaN,nrow=NCOL(X_estim),ncol=NCOL(y_estim))
				for (j in 1:NCOL(y_estim)){
					fit <- glmnet(X_estim, y_estim[,j], intercept = FALSE, family = 'gaussian', alpha = alpha_elnet, weights = weight_estim, nlambda = 3)
					mymax <- max(fit$lambda)
					myseq <- c(0,exp(seq(log(mymax),min(-12,log(min(fit$lambda))),length.out = 99)))
					fit <- cv.glmnet(X_estim, y_estim[,j], intercept = FALSE, family = 'gaussian', alpha = alpha_elnet, weights = weight_estim, lambda=myseq, nfolds = nfolds)
					
					coefs[,j] <- coef(fit,s = "lambda.min")
					num_select[j] <- unname(fit$nzero[which(fit$lambda == fit$lambda.min)])
				}
			}
			
			
			
		} else {
			if (intercept){
				fit <- lm(y_estim ~ X_estim, weights = weight_estim, tol = .Machine$double.xmin)
			} else {
				fit <- lm(y_estim ~ X_estim - 1, weights = weight_estim, tol = .Machine$double.xmin)
			}
			coefs <- coef(fit) # columns are different integrands
			num_select <- rep(NCOL(X_estim),NCOL(y_estim))
		}
		
		if (is.null(ncol(coefs))){
			coefs <- matrix(coefs,ncol=1)
		}
		
		if (intercept | regul_reg){
			fitteds_eval <- X_eval%*%coefs[2:NROW(coefs),]
			mse <- colMeans((y_eval - fitteds_eval - matrix(1,nrow=NROW(y_eval))%*%coefs[1,])^2)
		} else {
			fitteds_eval <- X_eval%*%coefs
			mse <- colMeans((y_eval - fitteds_eval - matrix(1,nrow=NROW(y_eval))%*%weight_estim%*%y_estim)^2)
		}
		
		
		if (any(is.na(fitteds_eval))){
			fitteds_eval <- matrix(0,NROW(y_eval),NCOL(y_eval))
		}
		
		integrand_new <- y_eval - fitteds_eval
		
		if (!integrand_logged){
			expectation <- weight_eval%*%integrand_new
		} else if(any(integrand_new<0)){
			expectation <- log(weight_eval%*%integrand_new) + max_integrand
		} else {
			expectation <- apply(integrand_new,2,function(x) logsumexp(log(x) + log_weights_eval)) + max_integrand
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
	
	
	return(list(expectation = expectation, num_select = num_select, mse = mse, coefs = coefs,
							integrand_logged = integrand_logged, est_inds = est_inds,
							polyorder = polyorder, regul_reg = regul_reg, alpha_elnet = alpha_elnet,
							nfolds = nfolds, apriori = apriori, intercept = intercept))
}
