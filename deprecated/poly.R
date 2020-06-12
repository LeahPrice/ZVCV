# #' The function getX is used to get the matrix of covariates for the regression based on a specified polynomial order
# #' @rdname zvcv
getX <- function(samples, derivatives, polyorder){
	
	N <- NROW(samples)
	d <- NCOL(samples)
	
	if (d==1 | (polyorder %in% c(1,2,3,4) & polyorder <= d)){
		return (getX_fast(samples,derivatives,polyorder))
	}
	
	Z <- -0.5*derivatives
	
	if (d==1){
		samples <- array(samples,c(N,d))
		Z <- array(Z,c(N,d))
	}
	
	X <- matrix(,nrow=N,ncol=0)
	for (mypol in 1:polyorder){
		alpha <- compositions(mypol,d) # from library(partitions)
		num_monos <- NCOL(alpha)
		for (k in 1:num_monos){
			term <- 0
			# First order derivatives wrt parameter j
			for (j in 1:d){
				# my_set_diff <- setdiff(c(1:d),j)
				# term <- term + max(alpha[j,k],0)*samples[,j]^(alpha[j,k] - 1)*Z[,j]  * apply(samples[,my_set_diff]^alpha[my_set_diff,k],1,prod)
				myprod <- max(alpha[j,k],0)*samples[,j]^(alpha[j,k] - 1)*Z[,j]
				for (jj in setdiff(c(1:d),j)){
					myprod <- myprod*samples[,jj]^alpha[jj,k]
				}
				term <- term + myprod
			}
			# Second order derivatives wrt parameter j
			for (j in 1:d){
				# my_set_diff <- setdiff(c(1:d),j)
				# term <- term + -0.5*max((alpha[j,k])*(alpha[j,k]-1),0)*samples[,j]^(alpha[j,k] - 2) * apply(samples[,my_set_diff]^alpha[my_set_diff,k],1,prod)
				myprod <- max((alpha[j,k])*(alpha[j,k]-1),0)*samples[,j]^(alpha[j,k] - 2)
				for (jj in setdiff(c(1:d),j)){
					myprod <- myprod*samples[,jj]^alpha[jj,k]
				}
				term <- term + -0.5*myprod
			}
			X <- cbind(X,term)
		}
	}
	
	return (X)
	
}

getX_fast <- function(samples, derivatives, polyorder){
	
    N <- NROW(samples)
	n_theta <- NCOL(samples)
	
	Z <- -0.5*derivatives
    
    if (n_theta==1){
        X <- matrix(,nrow=N,ncol=0)
        for (i in 1:polyorder){
            X <- cbind(X,i*samples^(i-1)*Z-0.5*i*(i-1)*samples^(i-2))
        }	
        return (X)
    }
    
	if (polyorder==1){
		return (Z)
	}
	
	
	squared <- 2*samples * Z - 1
	
	twoway <- matrix(, nrow = N, ncol = choose(n_theta,2))
	pos <- 0
	for(j in 1:(n_theta-1))
	{
		for(ii in (j+1):n_theta)
		{
			pos <- pos + 1
			twoway[,pos] <- samples[,ii]*Z[,j] + samples[,j]*Z[,ii]
		}
	}
	
	poly2 <- cbind(Z, squared, twoway)
	
	if (polyorder==2){
		return (poly2)
	}
	
	poly3 <- poly2
	
	cubed <- matrix(, nrow = N, ncol = n_theta)
	pos <- 0
	for(j in 1:n_theta)
	{
		pos <- pos + 1
		cubed[,pos] <- -3*samples[,j] + 3*samples[,j]^2*Z[,j]
	}
	
	threeway <- matrix(,nrow = N, ncol = choose(n_theta,3) + 2*choose(n_theta,2))
	
	pos <- 0
	# Three unique
	for(j in 1:(n_theta-2))
	{
		for(ii in (j+1):(n_theta-1))
		{
			for(kk in (ii+1):n_theta)
			{
				pos <- pos + 1
				threeway[,pos] <- samples[,j]*samples[,ii]*Z[,kk] + samples[,j]*samples[,kk]*Z[,ii] + samples[,ii]*samples[,kk]*Z[,j]
			}
		}
	}
	
	# A double up in the three
	for(j in 1:n_theta)
	{
		for(ii in 1:n_theta)
		{
			if (ii!=j){
				pos <- pos + 1
				threeway[,pos] <- -samples[,ii] + 2*samples[,j]*samples[,ii]*Z[,j] + samples[,j]^2*Z[,ii]
			}
		}
	}
	
	poly3 <- cbind(poly3, cubed)
	poly3 <- cbind(poly3, threeway)
	
	if (polyorder==3){
		return (poly3)
	}
	
	poly4 <- poly3
	
	fourway <- matrix(,nrow = N,ncol = choose(n_theta+4,n_theta) - choose(n_theta+3,n_theta))
	
	pos <- 0
	# Four unique
	for(j in 1:(n_theta-3))
	{
		for(ii in (j+1):(n_theta-2))
		{
			for(kk in (ii+1):(n_theta-1))
			{
				for(ll in (kk+1):n_theta)
				{
					pos <- pos + 1
					fourway[,pos] <- samples[,j]*samples[,ii]*samples[,kk]*Z[,ll] + samples[,j]*samples[,ii]*samples[,ll]*Z[,kk] + samples[,j]*samples[,kk]*samples[,ll]*Z[,ii] + samples[,ii]*samples[,kk]*samples[,ll]*Z[,j]
				}
			}
		}
	}
	
	# double for one part
	for(j in 1:n_theta)
	{
		for(ii in 1:(n_theta-1))
		{
			for(kk in (ii+1):n_theta)
			{
				if ( (ii!=j) && (ii!=kk) && (j!=kk) )
				{
					pos <- pos + 1
					fourway[,pos] <- -samples[,ii]*samples[,kk] + 2*samples[,j]*samples[,ii]*samples[,kk]*Z[,j] + samples[,j]^2*samples[,kk]*Z[,ii] + samples[,j]^2*samples[,ii]*Z[,kk]
				}
			}
		}
	}
	
	# double for two parts
	for(j in 1:(n_theta-1))
	{
		for(ii in (j+1):n_theta)
		{
			pos <- pos + 1
			fourway[,pos] <- -samples[,ii]^2 - samples[,j]^2 + 2*samples[,j]*samples[,ii]^2*Z[,j] + 2*samples[,j]^2*samples[,ii]*Z[,ii]
		}
	}
	
	# triple
	for(j in 1:n_theta)
	{
		for(ii in 1:n_theta)
		{
			if (ii!=j){
				pos <- pos + 1
				fourway[,pos] <- -3*samples[,j]*samples[,ii] + 3*samples[,j]^2*samples[,ii]*Z[,j] + samples[,j]^3*Z[,ii]
			}
		}
	}
	
	# to power four
	for(j in 1:n_theta)
	{
		pos <- pos + 1
		fourway[,pos] <- -6*samples[,j]^2 + 4*samples[,j]^3*Z[,j]
	}
	
	poly4 <- cbind(poly3,fourway)
	
	return (poly4)
	
}
