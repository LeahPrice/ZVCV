#' ZV-CV design matrix
#'
#' The function \code{getX} is used to get the matrix of covariates for the regression based on a specified polynomial order.
#' 
#' @param samples       An \eqn{N} by \eqn{d} matrix of samples from the target
#' @param derivatives	An \eqn{N} by \eqn{d} matrix of derivatives of the log target with respect to the parameters
#' @param polyorder     The order of the polynomial.
#' 
#' @return 			    The design matrix for the regression (except for the column of 1's for the intercept).
#' @seealso				\code{\link{Phi_fn}} for a very similar function for use in semi-exact control functionals. The function \code{\link{Phi_fn}} essentially gets the same matrix but with a column of ones added.
getX <- function(samples, derivatives, polyorder){
	
	N <- NROW(samples)
	d <- NCOL(samples)
	
	# There is a fast (hard-coded) version of this for d=1, polyorder<=4<=d in getX_fast.
	if (d==1 | (polyorder %in% c(1,2,3,4) & polyorder <= d)){
		return (getX_fast(samples,derivatives,polyorder))
	}
	
	if (d==1){
		samples <- array(samples,c(N,d))
		Z <- array(Z,c(N,d))
	}
	
	# If the partitions function is available, then using the compositions function is the fastest approach that I've come across.
	if (requireNamespace("partitions", quietly = TRUE)){
		alpha <- list()
		for (i in 1:polyorder){
			alpha[[i]] <- partitions::compositions(i,d)
		}
		
		return (getPoly_withpackage(samples,derivatives,alpha))
	}
	
	# However, the function partitions appears not to be Linux friendly, so the last resort is below.
	mymat <- matrix(0,nrow=polyorder,ncol=d)
	mymat[,1] <- 1:polyorder
	curr_inds <- 1:polyorder
	while (1){
		for (inds in curr_inds){
			for (j in 2:d){
				for (k in 1:polyorder){
					temp = mymat[inds,]
					if (sum(temp)<k){
						temp[j] <- k - sum(temp)
					}
					if ((temp[j]<=temp[j-1]) && sum(apply(mymat, 1, function(x) identical(temp, x)))==0){
						mymat <- rbind(mymat,temp)
					}
				}
			}
		}
		if (NROW(mymat)==curr_inds[length(curr_inds)]){
			break;
		}
		curr_inds <- (curr_inds[length(curr_inds)] + 1):NROW(mymat)
	}
	mymat <- mymat[order(rowSums(mymat),decreasing=TRUE),]
	
	final_mat_dup <- get_all_combins(mymat, polyorder)
	
	return(getPoly_withoutpackage(samples,derivatives,final_mat_dup))
	
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
