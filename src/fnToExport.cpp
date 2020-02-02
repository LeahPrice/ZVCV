#include <kernels.h>
#include <fnToExport.h>
#include <ctime>

using namespace std;

//' Squared norm matrix calculation
//' 
//' This function gets the matrix of square norms which is needed for the Gaussian, Matern and rational quadratic kernels.
//' Calculating this can help to save time if you are also interested in calculating the median heuristic, handling multiple tuning parameters
//' or trying other kernels in this group.
//'
//' @param samples An \eqn{N} by \eqn{d} matrix of samples from the target
//' @param nystrom_inds The (optional) sample indices to be used in the Nystrom approximation (for when using approximate SECF).
//'
//' @return An \eqn{N} by \eqn{N} matrix of squared norms between samples (or \eqn{N} by \eqn{k} where \eqn{k} is the length of \code{nystrom_inds}).
//'
//' @author Leah F. South
//' @seealso See \code{\link{medianTune}} and \code{\link{K0_fn}} for functions which use this.
// [[Rcpp::export]]
arma::mat squareNorm(const arma::mat & samples, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue){
  unsigned int N = samples.n_rows;
    
  unsigned int j;
  if (nystrom_inds.isNull()){ //full kernel matrix
    arma::mat z(N,N,arma::fill::zeros);
    for (unsigned int i = 0; i<N; i++){
      for (j = i; j<N; j++){
        z(i,j) = pow(norm(samples.row(i).t()-samples.row(j).t(),2),2);
        z(j,i) = z(i,j);
      }
    }
    return( z );
  } else{ //subset-based approach for Nystrom approximation
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat z(N,m0,arma::fill::zeros);
    for (unsigned int i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        z(i,j) = pow(norm(samples.row(i).t()-samples.row(inds(j)).t(),2),2);
      }
    }
    return( z );
  }
  
  return ( std::numeric_limits<arma::mat>::quiet_NaN() );
}


//' Median heuristic
//' 
//' This function calculates the median heuristic for use in e.g. the Gaussian, Matern and rational quadratic kernels.
//'
//' @param samples An \eqn{N} by \eqn{d} matrix of samples from the target
//' @param Z (optional) An NxN matrix of square norms, which can be calculated
//' using \code{\link{squareNorm}}, as long as the \code{nystrom_inds} are \code{NULL}.
//'
//' @return The median heuristic, which can then be used as the length-scale parameter in the Gaussian, Matern and rational quadratic kernels
//'
//' @references
//' Garreau, D., Jitkrittum, W. and Kanagawa, M. (2017). Large sample analysis of the median heuristic.  \url{https://arxiv.org/abs/1707.07269}
//'
//' @author Leah F. South
//' @seealso See \code{\link{medianTune}} and \code{\link{K0_fn}} for functions which use this.
// [[Rcpp::export]]
double medianTune(const arma::mat & samples, const Rcpp::Nullable<Rcpp::NumericMatrix> & Z = R_NilValue){
  
  unsigned int N = samples.n_rows;
  
  arma::mat dist_mat;
  if (Z.isNotNull()){
    dist_mat = Rcpp::as<arma::mat>(Z);
  } else {
    dist_mat = squareNorm(samples);
  }
  
  // Taking vector of dist_mat[i,j] where i<j (makes a slight difference compared to using the full matrix)
  unsigned int len = static_cast<unsigned int>(N*(N-1)/2);
  arma::vec dists(len);
  unsigned int k = 0;
  unsigned int i, j;
  for (i=0; i<(N-1); i++){
    if ((i+1)<=N){
      for (j=(i+1); j<N; j++){
        dists(k) = dist_mat(i,j);
        k++;
      }
    }
  }
  
  return( sqrt(arma::median(dists)/2) );
}



//' Kernel matrix calculation
//' 
//' This function calculates the full \eqn{K_0} matrix, which is a first or second order Stein operator applied to
//' a standard kernel. 
//' The output of this function can be used as an argument to \code{\link{CF}}, \code{\link{CF_crossval}},
//' \code{\link{SECF}}, \code{\link{SECF_crossval}}, \code{\link{aSECF}} and \code{\link{aSECF_crossval}}.
//' The kernel matrix is automatically computed in all of the above methods, but it is faster to calculate
//' in advance when using more than one of the above functions and when using any of the crossval functions.
//' 
//' @param samples An \eqn{N} by \eqn{d} matrix of samples from the target
//' @param derivatives	An \eqn{N} by \eqn{d} matrix of derivatives of the log target with respect to the parameters
//' @param sigma			The tuning parameters of the specified kernel. This involves a single length-scale parameter in "gaussian" and "RQ", a length-scale and a smoothness parameter in "matern" and two parameters in "product" and "prodsim". See below for further details.
//' @param steinOrder	This is the order of the Stein operator. The default is \code{1} in the control functionals paper (Oates et al, 2017) and \code{2} in the semi-exact control functionals paper (South et al, 2020).  The following values are currently available: \code{1} for all kernels and \code{2} for "gaussian", "matern" and "RQ". See below for further details.
//' @param kernel_function		Choose between "gaussian", "matern", "RQ", "product" or "prodsim". See below for further details.
//' @param Z (optional) An \eqn{N} by \eqn{N} (or \eqn{N} by \eqn{k} where \eqn{k} is the length of \code{nystrom_inds}). This can be calculated using \code{\link{squareNorm}}.
//' @param nystrom_inds (optional) The sample indices to be used in the Nystrom approximation (for when using approximate semi-exact control functionals).
//'
//' @return An \eqn{N} by \eqn{N} kernel matrix (or \eqn{N} by \eqn{k} where \eqn{k} is the length of \code{nystrom_inds}).
//'
//' @section On the choice of \eqn{\sigma}, the kernel and the Stein order:
//' The kernel in Stein-based kernel methods is \eqn{L_x L_y k(x,y)} where \eqn{L_x} is a first or second order Stein operator in \eqn{x} and \eqn{k(x,y)} is some generic kernel to be specified.
//'
//' The Stein operators for distribution \eqn{p(x)} are defined as:
//' \itemize{
//' \item \strong{\code{steinOrder=1}}: \eqn{L_x g(x) = \nabla_x^T g(x) + \nabla_x \log p(x)^T g(x)} (see e.g. Oates el al (2017))
//' \item \strong{\code{steinOrder=2}}: \eqn{L_x g(x) = \Delta_x g(x) + \nabla_x log p(x)^T \nabla_x g(x)} (see e.g. South el al (2020))
//' }
//' Here \eqn{\nabla_x} is the first order derivative wrt \eqn{x} and \eqn{\Delta_x = \nabla_x^T \nabla_x} is the Laplacian operator.
//'
//' The generic kernels which are implemented in this package are listed below.  Note that the input parameter \strong{\code{sigma}} defines the kernel parameters \eqn{\sigma}. 
//' \itemize{
//' \item \strong{\code{"gaussian"}}: A Gaussian kernel,
//' \deqn{k(x,y) = exp(-z(x,y)/\sigma^2)}
//' \item \strong{{\code{"matern"}}}: A Matern kernel with \eqn{\sigma = (\lambda,\nu)},
//' \deqn{k(x,y) = bc^{\nu}z(x,y)^{\nu/2}K_{\nu}(c z(x,y)^{0.5})} where \eqn{b=2^{1-\nu}(\Gamma(\nu))^{-1}}, \eqn{c=(2\nu)^{0.5}\lambda^{-1}} and \eqn{K_{\nu}(x)} is the modified Bessel function of the second kind. Note that \eqn{\lambda} is the length-scale parameter and \eqn{\nu} is the smoothness parameter (which defaults to 2.5 for \eqn{steinOrder=1} and 4.5 for \eqn{steinOrder=2}).
//' \item \strong{\code{"RQ"}}: A rational quadratic kernel,
//' \deqn{k(x,y) = (1+\sigma^{-2}z(x,y))^{-1}}
//' \item \strong{\code{"product"}}: The product kernel that appears in Oates et al (2017) with \eqn{\sigma = (a,b)}
//' \deqn{k(x,y) = (1+a z(x) + a z(y))^{-1} exp(-0.5 b^{-2} z(x,y)) }
//' \item \strong{\code{"prodsim"}}: A slightly different product kernel with \eqn{\sigma = (a,b)} (see e.g. \url{https://www.imperial.ac.uk/inference-group/projects/monte-carlo-methods/control-functionals/}),
//' \deqn{k(x,y) = (1+a z(x))^{-1}(1 + a z(y))^{-1} exp(-0.5 b^{-2} z(x,y)) }
//' }
//' In the above equations, \eqn{z(x) = \sum_j x[j]^2} and \eqn{z(x,y) = \sum_j (x[j] - y[j])^2}. For the last two kernels, \code{steinOrder} must be \code{1}. Each combination of \code{steinOrder} and \code{kernel_function} above is currently hard-coded but it may be possible to extend this to other kernels in future versions using autodiff. The calculations for the first three kernels above are detailed in South et al (2020).
//'
//' @references
//' Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.
//'
//' South, L. F., Karvonen, T., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/}
//'
//' @author Leah F. South
// [[Rcpp::export]]
arma::mat K0_fn(const arma::mat & samples, const arma::mat & derivatives, arma::vec sigma, unsigned int steinOrder, std::string kernel_function, const Rcpp::Nullable<Rcpp::NumericMatrix> & Z = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue){
  unsigned int len_sigma = sigma.n_rows;
  
  arma::mat z = getSqNorm(samples,nystrom_inds,Z); // getting squared norms.
  
  if ((steinOrder!=1) & (steinOrder!=2)){
    throw(Rcpp::exception("The steinOrder should be either 1 or 2."));
  }
  
  if (kernel_function=="gaussian"){//(strcmp(kernel_function,"gaussian")==0){
    if (len_sigma!=1){
      throw(Rcpp::exception("The gaussian kernel should have a single parameter."));
    }
    double sigma_new = arma::as_scalar(sigma);
    return (gaussian_k(steinOrder,samples,derivatives,sigma_new,kernel_function,z,nystrom_inds));
    
  } else if (kernel_function=="matern"){
    z = z + pow(10.0,-10.0);
    if (len_sigma==1){
      arma::vec sigma_new(2);
      sigma_new(0) = sigma(0);
      if (steinOrder==2){
        sigma_new(1) = 4.5;
        Rcpp::warning("Using a default nu of 4.5 for the matern kernel.");
      } else{
        sigma_new(1) = 2.5;
        Rcpp::warning("Using a default nu of 2.5 for the matern kernel.");
      }
      return (matern_k(steinOrder,samples,derivatives,sigma_new,kernel_function,z,nystrom_inds));
    } else if (len_sigma==2){
      if ((sigma(1)<4.5) & (steinOrder==2)){
        throw(Rcpp::exception("The nu parameter for the matern kernel should be at least 4.5 for steinOrder=2 in this implementation."));
      } else if ((sigma(1)<2.5) & (steinOrder==1)){
        throw(Rcpp::exception("The nu parameter for the matern kernel should be at least 2.5 for steinOrder=1 in this implementation."));
      }
      return (matern_k(steinOrder,samples,derivatives,sigma,kernel_function,z,nystrom_inds));
    }
    throw(Rcpp::exception("The matern kernel should have one or two parameters."));
    
  } else if (kernel_function=="RQ"){
    if (len_sigma!=1){
      throw(Rcpp::exception("The rational quadratic kernel should a single parameter."));
    }
    double sigma_new = arma::as_scalar(sigma);
    return (RQ_k(steinOrder,samples,derivatives,sigma_new,kernel_function,z,nystrom_inds));
    
  } else if (kernel_function=="product"){
    if ((steinOrder==2) or (len_sigma!=2)){
      throw(Rcpp::exception("The product kernel is only implemented for steinOrder=1 and two parameters."));
    }
    return (product_k(1,samples,derivatives,sigma,kernel_function,z,nystrom_inds));
    
  } else if (kernel_function=="prodsim"){
    if ((steinOrder==2) or (len_sigma!=2)){
      throw(Rcpp::exception("The prodsim kernel is only implemented for steinOrder=1 and two parameters."));
    }
    return (prodsim_k(1,samples,derivatives,sigma,kernel_function,z,nystrom_inds));
    
  }
  
  throw(Rcpp::exception("Enter a valid kernel name."));
  return ( std::numeric_limits<arma::mat>::quiet_NaN() );
  
}




//' Nearest symmetric positive definite matrix
//' 
//' This function finds the nearest symmetric positive definite matrix to the given matrix.
//' It is used throughout the package to handle numerical issues in matrix inverses
//' and cholesky decompositions.
//'
//' @param K0 A square matrix
//'
//' @return The closest symmetric positive definite matrix to K0.
//'
//' @references
//' Higham, N. J. (1988). Computing a nearest symmetric positive semidefinite matrix. Linear Algebra and its Applications, 103, 103-118.
//'
//' D'Errico, J. (2013). nearestSPD Matlab function. \url{https://uk.mathworks.com/matlabcentral/fileexchange/42885-nearestspd}.
//'
//' @author Adapted from Matlab code by John D'Errico
// [[Rcpp::export]]
arma::mat nearPD(arma::mat K0){
	// Based off of Matlab code by John D'Errico
	
	if (K0.is_sympd() == false){
		unsigned int N = K0.n_rows;
		unsigned int N1 = K0.n_cols;
		
		if (N != N1){
			throw(Rcpp::exception("K0 must be a square matrix"));
		}
		
		// symmetry
		K0 = (K0 + K0.t())/2;
		
		// svd
		arma::mat U;
		arma::vec s;
		arma::mat V;
		
		arma::svd(U,s,V,K0);
		arma::mat H = V*arma::diagmat(s)*V.t();
		
		K0 = (K0+H)/2;
		
		// symmetry
		K0 = (K0 + K0.t())/2;
		
		bool SPDcheck = false;
		double k = 0.0;
		double mineig;
		while (SPDcheck == false){
			SPDcheck = K0.is_sympd();
			k = k + 1.0;
			if (SPDcheck == false){
				mineig = arma::as_scalar(arma::min(arma::eig_sym(K0)));
				K0 = K0 + (-mineig*pow(k,2.0) + arma::eps(mineig))*arma::eye<arma::mat>(N,N);
			}
		}
	}
	
	return(K0);
}

// The internal function for getting the Phi matrix in SECF. There is an R wrapper for this function.
// [[Rcpp::export]]
arma::mat Phi_fn_cpp(const arma::mat & samples, const arma::mat & derivatives, Rcpp::Function getX, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue){
  unsigned int N = samples.n_rows;
  
  unsigned int poly_q;
  if (polyorder.isNotNull()) {
    poly_q = Rcpp::as<unsigned int>(polyorder);
  } else{
    poly_q = 1;
  }
  
  // phi matrix of 2nd order stein operator acting on functions at each point
  // phi is based on a given order polynomial
  // subs is the subset of variables used
  arma::mat phi;
  if (subset.isNotNull()){
    arma::uvec subs = Rcpp::as<arma::uvec>(subset) -1;
    arma::mat samples_sub = samples.cols(subs);
    arma::mat derivatives_sub = derivatives.cols(subs);
    phi = Rcpp::as<arma::mat>( getX(samples_sub,derivatives_sub,poly_q) );
  } else{
    phi = Rcpp::as<arma::mat>( getX(samples,derivatives,poly_q) );
  }
  phi.insert_cols(0,arma::ones<arma::colvec>(N));
  
  return(phi);
}

// The internal function for performing control functionals when samples are NOT split for estimation and evaluation. There is an R wrapper that combines this function and CF_unbiased_cpp.
// [[Rcpp::export]]
Rcpp::List CF_cpp(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<arma::vec> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, bool one_in_denom = false, bool diagnostics = false){
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  arma::mat K0mat = getK0(samples,derivatives,K0,steinOrder,kernel_function,sigma);
  
  K0mat = nearPD(K0mat);
  
  if (!diagnostics){
    arma::vec K0inv_ones = arma::solve(K0mat,arma::ones<arma::vec>(N));
    arma::vec expectation;
    if (one_in_denom){
      expectation = integrands.t()*K0inv_ones / arma::as_scalar(1.0 + arma::ones<arma::rowvec>(N) * K0inv_ones);
    } else{
      expectation = integrands.t()*K0inv_ones / arma::as_scalar(arma::ones<arma::rowvec>(N) * K0inv_ones);
    }
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation) );
  }
  
  arma::mat K0inv = K0mat.i();
  arma::vec K0inv_ones = K0inv*arma::ones<arma::vec>(N);
  arma::vec weights;
  if (one_in_denom){
    weights = K0inv_ones/ arma::as_scalar(1.0 + arma::ones<arma::rowvec>(N) * K0inv_ones);
  } else{
    weights = K0inv_ones/ arma::as_scalar(arma::ones<arma::rowvec>(N) * K0inv_ones);
  }
  arma::vec expectation = integrands.t()*weights;
  arma::mat diff = integrands - arma::ones<arma::vec>(N)*expectation.t();
  arma::mat a = K0inv*diff;
  double ksd = pow(arma::as_scalar( weights.t() * K0mat * weights ), 0.5);
  arma::vec bound_const(N_expectations);
  for (unsigned int i = 0; i<N_expectations; i++){
    bound_const(i) = pow(arma::as_scalar( a.col(i).t() * diff.col(i) ), 0.5); // a.t()*K0*a
  }
  
  return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation,
                              Rcpp::Named("a") = a, Rcpp::Named("b") = expectation.t(), Rcpp::Named("ksd") = ksd, Rcpp::Named("bound_const") = bound_const) );
}

// The internal function for performing control functionals when samples ARE split for estimation and evaluation. There is an R wrapper that combines this function and CF_cpp.
// [[Rcpp::export]]
Rcpp::List CF_unbiased_cpp(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, arma::uvec est_inds, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<arma::vec> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, bool one_in_denom = false, bool diagnostics = false){
  unsigned int N = samples.n_rows;
  
  arma::mat K0mat = getK0(samples,derivatives,K0,steinOrder,kernel_function,sigma);
  
  arma::uvec all_indx = arma::linspace<arma::uvec >(0, N-1,N);
  
  est_inds = est_inds - 1;
  unsigned int N_est = est_inds.n_rows;
  
  arma::uvec eval_inds = all_indx;
  // eval_inds.shed_rows(est_inds);
  arma::uvec est_inds2 = sort(est_inds,"descend");
  for (unsigned int zz = 0; zz<N_est; zz++){
    eval_inds.shed_row(est_inds2(zz));
  }
  
  arma::mat K0mat_est = nearPD(K0mat.submat(est_inds,est_inds));
  arma::mat K0mat_eval = K0mat.submat(eval_inds,est_inds);
  
  arma::mat f_est = integrands.rows(est_inds);
  arma::mat f_eval = integrands.rows(eval_inds);
  
  // Estimation
  arma::mat K0inv = arma::inv(K0mat_est);
  arma::vec K0inv_ones = K0inv*arma::ones<arma::vec>(N_est);
  //arma::vec K0inv_ones = arma::solve(K0mat_est,arma::ones<arma::vec>(N_est));
  arma::vec weights;
  if (one_in_denom){
    weights = K0inv_ones/ arma::as_scalar(1.0 + arma::ones<arma::rowvec>(N_est) * K0inv_ones);
  } else{
    weights = K0inv_ones/ arma::as_scalar(arma::ones<arma::rowvec>(N_est) * K0inv_ones);
  }
  arma::rowvec beta = weights.t()*f_est;
  arma::mat diff = f_est - arma::ones<arma::vec>(N_est)*beta;
  arma::mat a = K0inv*diff;
  arma::mat f_hat = arma::ones<arma::vec>(N - N_est)*beta + K0mat_eval*a;
  arma::vec expectation = 1/(N - N_est)*arma::sum(f_eval - f_hat,0).t() + beta.t();
  
  if (!diagnostics){
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat));
  }
  
  double ksd = pow(arma::as_scalar( weights.t() * K0mat_est * weights ), 0.5);
  
  return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat,
                              Rcpp::Named("a") = a, Rcpp::Named("b") = beta, Rcpp::Named("ksd") = ksd) );
  
}

// The internal function for performing cross-validation to select from a list of sigma (kernel tuning parameter) values or
// a list of K0 matrices (which can correspond to different sigma, steinOrder and kernel_function values) in control functionals. There is an R wrapper for this function.
// If the K0 matrix is not given, this function recalculates the K0 matrix once for each selected tuning parameter (there can be multiple tuning parameters selected if there are multiple integrands).
// This double-up on computation can be avoided by using the K0 argument.
// [[Rcpp::export]]
Rcpp::List CF_crossval_cpp(arma::mat integrands, arma::mat samples, arma::mat derivatives, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<Rcpp::List> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::List> & K0 = R_NilValue, Rcpp::Nullable<unsigned int> folds = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds = R_NilValue, bool one_in_denom = false, bool diagnostics = false){
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  unsigned int n_sigma = 0;
  unsigned int n_K0 = 0;
  bool sigma_given = false;
  bool K0_given = false;
  
  Rcpp::List sig_list;
  Rcpp::List K0_list;
  
  if (sigma.isNotNull()) {
    sig_list = Rcpp::as<Rcpp::List>(sigma);
    sigma_given = true;
    n_sigma = sig_list.length();
  }
  
  if (K0.isNotNull()) {
    K0_list = Rcpp::as<Rcpp::List>(K0);
    K0_given = true;
    n_K0 = K0_list.length();
  }
  
  if ( (K0_given == false) & (sigma_given == false)){
    throw(Rcpp::exception("Must give sigma and/or K0 to do cross validation."));
  }
  if (sigma_given & K0_given & (n_K0!=n_sigma)){
    throw(Rcpp::exception("If K0 and sigma are both given then their sizes need to match."));
  }
  if (sigma_given==false){
    n_sigma = n_K0;
  }
  
  arma::mat mse(N_expectations,n_sigma);
  
  for (unsigned int j=0; j<n_sigma; j++){
    //Rcpp::Rcout << "Iteration " << j << std::endl;
    if (K0_given){
      Rcpp::NumericMatrix K0curr = K0_list[j];
      Rcpp::Nullable<Rcpp::NumericMatrix> Kcurr(K0curr);
      mse.col(j) = CF_mse_cpp(integrands, samples, derivatives, steinOrder, kernel_function, R_NilValue, Kcurr, folds, est_inds, one_in_denom);
    } else {
      Rcpp::NumericVector sigs = sig_list[j];
      Rcpp::Nullable<Rcpp::NumericVector> sigs1(sigs);
      mse.col(j) = CF_mse_cpp(integrands, samples, derivatives, steinOrder, kernel_function, sigs1, R_NilValue, folds, est_inds, one_in_denom);
    }
  }
  
  arma::ucolvec opt_indices = arma::index_min(mse,1);
  
  arma::ucolvec inds_unique = unique(opt_indices);
  unsigned int length_unique = inds_unique.n_rows;
  Rcpp::List temp;
  
  arma::vec expectation(N_expectations);
  unsigned int N_eval;
  if (est_inds.isNotNull()){
    arma::uvec inds_est = Rcpp::as<arma::uvec>(est_inds);
    N_eval = N - inds_est.n_rows;
  } else{
    N_eval = 0;
  }
  arma::mat f_true(N_eval,N_expectations);
  arma::mat f_hat(N_eval,N_expectations);
  arma::mat a(N - N_eval,N_expectations);
  arma::rowvec beta(N_expectations);
  arma::vec ksd(N_expectations);
  arma::vec bound_const(N_expectations);
  
  Rcpp::NumericVector MyNumVec;
  Rcpp::NumericMatrix test;
  for (unsigned int i = 0; i<length_unique; i++){
    arma::uvec q1 = find(opt_indices == inds_unique(i));
    arma::mat integrands_temp(N,q1.n_rows);
    for (unsigned int k = 0; k<(q1.n_rows); k++){
      integrands_temp.col(k) = integrands.col(q1(k));
    }
    
    if (K0_given){
      Rcpp::NumericMatrix K0curr_opt = K0_list[inds_unique(i)];
      Rcpp::Nullable<Rcpp::NumericMatrix> Kcurr_opt(K0curr_opt);
      if (est_inds.isNull()){
        temp = CF_cpp(integrands_temp, samples, derivatives, steinOrder, kernel_function, R_NilValue, Kcurr_opt, one_in_denom, diagnostics);
      } else{
        arma::uvec inds_est = Rcpp::as<arma::uvec>(est_inds);
        temp = CF_unbiased_cpp(integrands_temp, samples, derivatives, inds_est, steinOrder, kernel_function, R_NilValue, Kcurr_opt, one_in_denom, diagnostics);
      }
    } else {
      Rcpp::NumericVector sigs_opt = sig_list[inds_unique(i)];
      //arma::vec sigs_opt_vec = Rcpp::as<arma::vec>(wrap(sigs_opt));
      Rcpp::Nullable<arma::vec> sigs1_opt(wrap(sigs_opt));
      if (est_inds.isNull()){
        temp = CF_cpp(integrands_temp, samples, derivatives, steinOrder, kernel_function, sigs1_opt, R_NilValue, one_in_denom, diagnostics);
      } else{
        arma::uvec inds_est = Rcpp::as<arma::uvec>(est_inds);
        temp = CF_unbiased_cpp(integrands_temp, samples, derivatives, inds_est, steinOrder, kernel_function, sigs1_opt, R_NilValue, one_in_denom, diagnostics);
      }
    }
    for (unsigned int kk = 0; kk<q1.n_rows; kk++){
      MyNumVec = temp["expectation"];
      expectation(q1(kk)) = MyNumVec[kk];
      if (est_inds.isNotNull()){
        test = Rcpp::as<Rcpp::NumericMatrix>(wrap(temp["f_true"]));
        f_true.col(q1(kk)) = Rcpp::as<arma::vec>(wrap(test.column(kk)));
        test = Rcpp::as<Rcpp::NumericMatrix>(wrap(temp["f_hat"]));
        f_hat.col(q1(kk)) = Rcpp::as<arma::vec>(wrap(test.column(kk)));
      }
      if (diagnostics){
        test = Rcpp::as<Rcpp::NumericMatrix>(wrap(temp["a"]));
        a.col(q1(kk)) = Rcpp::as<arma::vec>(wrap(test.column(kk)));
        MyNumVec = temp["b"];
        beta(q1(kk)) = MyNumVec[kk];
        MyNumVec = temp["ksd"];
        ksd(q1(kk)) = temp["ksd"];
        if (est_inds.isNull()){
          MyNumVec = temp["bound_const"];
          bound_const(q1(kk)) = MyNumVec[kk];
        }
      }
    }
  }  
  
  if (!diagnostics){
    if (est_inds.isNull()){
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1) );
    } else{
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_true, Rcpp::Named("f_hat") = f_hat, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1) );
    }
  }
  
  if (est_inds.isNull()){
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1,
                                Rcpp::Named("a") = a, Rcpp::Named("b") = beta, Rcpp::Named("ksd") = ksd, Rcpp::Named("bound_const") = bound_const) );
  } else{
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_true, Rcpp::Named("f_hat") = f_hat, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1,
                                Rcpp::Named("a") = a, Rcpp::Named("b") = beta, Rcpp::Named("ksd") = ksd) );
  }
  
}

// The internal function for performing semi-exact control functionals when samples are NOT split for estimation and evaluation. There is an R wrapper that combines this function and SECF_unbiased_cpp.
// [[Rcpp::export]]
Rcpp::List SECF_cpp(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, Rcpp::Function getX, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<arma::vec> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, bool diagnostics = false){
  unsigned int N_expectations = integrands.n_cols;
  
  arma::mat phi = Phi_fn_cpp(samples,derivatives,getX,polyorder,subset);
  
  arma::mat K0mat = getK0(samples,derivatives,K0,steinOrder,kernel_function,sigma);
  
  K0mat = nearPD(K0mat);
  
  if (!diagnostics){
    arma::mat K0inv_phi = arma::solve(K0mat,phi);
    arma::mat beta = arma::solve( phi.t() * K0inv_phi, K0inv_phi.t() )*integrands;
    arma::vec expectation = beta.row(0).t();
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation) );
  }
  
  arma::mat K0inv = K0mat.i();
  arma::mat K0inv_phi = K0inv*phi;
  arma::mat weights_prep = arma::solve( phi.t() * K0inv_phi, K0inv_phi.t() );
  arma::vec weights = weights_prep.row(0).t();
  arma::mat beta = weights_prep*integrands;
  arma::vec expectation = beta.row(0).t();
  arma::mat diff = integrands - phi * beta;
  arma::mat a = K0inv*diff;
  arma::vec bound_const(N_expectations);
  for (unsigned int i = 0; i<N_expectations; i++){
    bound_const(i) = pow(arma::as_scalar( a.col(i).t() * diff.col(i) ), 0.5);
  }
  double ksd = pow(arma::as_scalar( weights.t() * K0mat * weights ), 0.5);
  return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation,
                              Rcpp::Named("a") = a, Rcpp::Named("b") = beta, Rcpp::Named("ksd") = ksd, Rcpp::Named("bound_const") = bound_const) );
  
}

// The internal function for performing semi-exact control functionals when samples ARE split for estimation and evaluation. There is an R wrapper that combines this function and SECF_cpp.
// [[Rcpp::export]]
Rcpp::List SECF_unbiased_cpp(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, arma::uvec est_inds, Rcpp::Function getX, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<arma::vec> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, bool diagnostics = false){
  
  unsigned int N = samples.n_rows;
  
  arma::mat phi = Phi_fn_cpp(samples,derivatives,getX,polyorder,subset);
  
  arma::mat K0mat = getK0(samples,derivatives,K0,steinOrder,kernel_function,sigma);
  
  
  arma::uvec all_indx = arma::linspace<arma::uvec >(0, N-1,N);
  est_inds = est_inds - 1;
  unsigned int N_est = est_inds.n_rows;
  
  arma::uvec eval_inds = all_indx;
  // eval_inds.shed_rows(est_inds);
  arma::uvec est_inds2 = sort(est_inds,"descend");
  for (unsigned int zz = 0; zz<N_est; zz++){
    eval_inds.shed_row(est_inds2(zz));
  }
  
  arma::mat K0mat_est = nearPD(K0mat.submat(est_inds,est_inds));
  arma::mat K0mat_eval = K0mat.submat(eval_inds,est_inds);
  arma::mat phi_est = phi.rows(est_inds);
  arma::mat phi_eval = phi.rows(eval_inds);
  
  
  arma::mat f_est = integrands.rows(est_inds);
  arma::mat f_eval = integrands.rows(eval_inds);
  
  arma::mat K0inv = arma::inv(K0mat_est);
  arma::mat K0inv_phi = K0inv*phi_est;
  //arma::mat K0inv_phi = arma::solve(K0mat_est,phi_est);
  
  arma::mat weights_prep = arma::solve( phi_est.t() * K0inv_phi, K0inv_phi.t() );
  arma::mat beta = weights_prep * f_est;
  arma::mat diff = f_est-phi_est * beta;
  arma::mat a = K0inv*diff;
  arma::mat f_hat = phi_eval * beta + K0mat_eval * a;
  arma::vec expectation =  1/(N - N_est)*arma::sum(f_eval - f_hat,0).t() + beta.row(0).t();
  
  if (!diagnostics){
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat) );
  }
  
  double ksd = pow(arma::as_scalar( weights_prep.row(0) * K0mat_est * weights_prep.row(0).t() ), 0.5);
  return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat,
                              Rcpp::Named("a") = a, Rcpp::Named("b") = beta, Rcpp::Named("ksd") = ksd) );
  
}

// The internal function for performing cross-validation to select from a list of sigma (kernel tuning parameter) values or
// a list of K0 matrices (which can correspond to different sigma, steinOrder and kernel_function values) in semi-exact control functionals. There is an R wrapper for this function.
// If the K0 matrix is not given, this function recalculates the K0 matrix once for each selected tuning parameter (there can be multiple tuning parameters selected if there are multiple integrands).
// This double-up on computation can be avoided by using the K0 argument.
// [[Rcpp::export]]
Rcpp::List SECF_crossval_cpp(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, Rcpp::Function getX, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<Rcpp::List> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::List> & K0 = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, Rcpp::Nullable<unsigned int> folds = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds = R_NilValue, bool diagnostics = false){
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  unsigned int n_sigma = 0;
  unsigned int n_K0 = 0;
  bool sigma_given = false;
  bool K0_given = false;
  
  Rcpp::List sig_list;
  Rcpp::List K0_list;
  
  if (sigma.isNotNull()) {
    sig_list = Rcpp::as<Rcpp::List>(sigma);
    sigma_given = true;
    n_sigma = sig_list.length();
  }
  
  std::string kernel_fun;
  if (K0.isNotNull()) {
    K0_list = Rcpp::as<Rcpp::List>(K0);
    K0_given = true;
    n_K0 = K0_list.length();
  }
  
  if ( (K0_given == false) & (sigma_given == false)){
    throw(Rcpp::exception("Must give sigma and/or K0 to do cross validation."));
  }
  if (sigma_given & K0_given & (n_K0!=n_sigma)){
    throw(Rcpp::exception("If K0 and sigma are both given then their sizes need to match."));
  }
  if (sigma_given==false){
    n_sigma = n_K0;
  }
  
  arma::mat mse(N_expectations,n_sigma);
  
  for (unsigned int j=0; j<n_sigma; j++){
    //Rcpp::Rcout << "Iteration " << j << std::endl;
    if (K0_given){
      Rcpp::NumericMatrix K0curr = K0_list[j];
      Rcpp::Nullable<Rcpp::NumericMatrix> Kcurr(K0curr);
      mse.col(j) = SECF_mse_cpp(integrands, samples, derivatives, getX, polyorder, steinOrder, kernel_function, R_NilValue, Kcurr, subset, folds, est_inds);
    } else {
      Rcpp::NumericVector sigs = sig_list[j];
      Rcpp::Nullable<Rcpp::NumericVector> sigs1(sigs);
      mse.col(j) = SECF_mse_cpp(integrands, samples, derivatives, getX, polyorder, steinOrder, kernel_function, sigs1, R_NilValue, subset, folds, est_inds);
    }
  }
  
  arma::ucolvec opt_indices = arma::index_min(mse,1);
  arma::ucolvec inds_unique = unique(opt_indices);
  unsigned int length_unique = inds_unique.n_rows;
  Rcpp::List temp;
  
  arma::vec expectation(N_expectations);
  unsigned int N_eval, N_estim;
  if (est_inds.isNotNull()){
    arma::uvec inds_est = Rcpp::as<arma::uvec>(est_inds);
    N_estim = inds_est.n_rows;
    N_eval = N - N_estim;
  } else{
    N_estim = N;
    N_eval = 0;
  }
  
  unsigned int Q = Phi_fn_cpp(samples.row(0), derivatives.row(0), getX, polyorder).n_cols;
  
  arma::mat f_true(N_eval,N_expectations);
  arma::mat f_hat(N_eval,N_expectations);
  arma::mat a(N_estim,N_expectations);
  arma::mat beta(Q,N_expectations);
  arma::vec ksd(N_expectations);
  arma::vec bound_const(N_expectations);
  
  Rcpp::NumericVector MyNumVec;
  Rcpp::NumericMatrix test;
  for (unsigned int i = 0; i<length_unique; i++){
    arma::uvec q1 = find(opt_indices == inds_unique(i));
    arma::mat integrands_temp(N,q1.n_rows);
    for (unsigned int k = 0; k<(q1.n_rows); k++){
      integrands_temp.col(k) = integrands.col(q1(k));
    }
    
    if (K0_given){
      Rcpp::NumericMatrix K0curr_opt = K0_list[inds_unique(i)];
      Rcpp::Nullable<Rcpp::NumericMatrix> Kcurr_opt(K0curr_opt);
      if (est_inds.isNull()){
        temp = SECF_cpp(integrands_temp, samples, derivatives, getX, polyorder, steinOrder, kernel_function, R_NilValue, Kcurr_opt, subset, diagnostics);
      } else{
        arma::uvec inds_est = Rcpp::as<arma::uvec>(est_inds);
        temp = SECF_unbiased_cpp(integrands_temp, samples, derivatives, inds_est, getX, polyorder, steinOrder, kernel_function, R_NilValue, Kcurr_opt, subset, diagnostics);
      }
    } else {
      Rcpp::NumericVector sigs_opt = sig_list[inds_unique(i)];
      //arma::vec sigs_opt_vec = Rcpp::as<arma::vec>(wrap(sigs_opt));
      Rcpp::Nullable<arma::vec> sigs1_opt(wrap(sigs_opt));
      if (est_inds.isNull()){
        temp = SECF_cpp(integrands_temp, samples, derivatives, getX, polyorder, steinOrder, kernel_function, sigs1_opt, R_NilValue, subset, diagnostics);
      } else{
        arma::uvec inds_est = Rcpp::as<arma::uvec>(est_inds);
        temp = SECF_unbiased_cpp(integrands_temp, samples, derivatives, inds_est, getX, polyorder, steinOrder, kernel_function, sigs1_opt, R_NilValue, subset, diagnostics);
      }
    }
    for (unsigned int kk = 0; kk<q1.n_rows; kk++){
      MyNumVec = temp["expectation"];
      expectation(q1(kk)) = MyNumVec[kk];
      if (est_inds.isNotNull()){
        test = Rcpp::as<Rcpp::NumericMatrix>(wrap(temp["f_true"]));
        f_true.col(q1(kk)) = Rcpp::as<arma::vec>(wrap(test.column(kk)));
        test = Rcpp::as<Rcpp::NumericMatrix>(wrap(temp["f_hat"]));
        f_hat.col(q1(kk)) = Rcpp::as<arma::vec>(wrap(test.column(kk)));
      }
      if (diagnostics){
        test = Rcpp::as<Rcpp::NumericMatrix>(wrap(temp["a"]));
        a.col(q1(kk)) = Rcpp::as<arma::vec>(wrap(test.column(kk)));
        test = Rcpp::as<Rcpp::NumericMatrix>(wrap(temp["b"]));
        beta.col(q1(kk)) = Rcpp::as<arma::vec>(wrap(test.column(kk)));
        ksd(q1(kk)) = temp["ksd"];
        if (est_inds.isNull()){
          MyNumVec = temp["bound_const"];
          bound_const(q1(kk)) = MyNumVec[kk];
        }
      }
    }
  }
  
  if (!diagnostics){
    if (est_inds.isNull()){
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1) );
    } else{
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_true, Rcpp::Named("f_hat") = f_hat, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1) );
    }
  }
  
  if (est_inds.isNull()){
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1,
                                Rcpp::Named("a") = a, Rcpp::Named("b") = beta, Rcpp::Named("ksd") = ksd, Rcpp::Named("bound_const") = bound_const) );
  } else{
    return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_true, Rcpp::Named("f_hat") = f_hat, Rcpp::Named("mse") = mse, Rcpp::Named("optinds") = opt_indices + 1,
                                Rcpp::Named("a") = a, Rcpp::Named("b") = beta, Rcpp::Named("ksd") = ksd) );
  }
  
}

// This internal function is used to get the linear system that requires solving in approximate semi-exact control functionals. The conjugate gradient step is done in R so that an R library can be used.
// [[Rcpp::export]]
Rcpp::List aSECF_cpp_prep(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, Rcpp::Function getX, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<arma::vec> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue, bool conjugate_gradient = true){
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  unsigned int m0;
  arma::uvec inds;
  if (nystrom_inds.isNotNull()){
    inds = Rcpp::as<arma::uvec>(nystrom_inds) -1;
    m0 = inds.n_rows;
  } else{
    m0 = std::ceil( std::sqrt( static_cast<double>(N) ) );
    inds = Rcpp::RcppArmadillo::sample(arma::linspace<arma::uvec>(0, N-1, N), m0, false);
  }
  
  Rcpp::IntegerVector temp_inds = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(inds + 1));
  Rcpp::Nullable<Rcpp::IntegerVector> inds_nullable(temp_inds);
  
  arma::mat phi = Phi_fn_cpp(samples,derivatives,getX,polyorder,subset);
  unsigned int q = phi.n_cols;
  arma::mat phi_sub = phi.rows(inds);
  
  arma::mat K0mat;
  std::string kernel_fun;
  arma::vec sigmas;
  if (K0.isNotNull()) {
    K0mat = Rcpp::as<arma::mat>(K0);
    if ((K0mat.n_rows == N) & (K0mat.n_cols == N)){
      K0mat = K0mat.cols(inds);
    } else if (K0mat.n_cols != m0){
      throw(Rcpp::exception("If K0 is given, then nystrom_inds must also be given. The number of columns in K0 match either the number of samples or the length of the indices."));
    }
  } else {
    K0mat = getK0(samples,derivatives,K0,steinOrder,kernel_function,sigma,inds_nullable);
  }
  
  arma::mat A;
  arma::mat b(m0+q,N_expectations);
  arma::mat B2 = arma::eye(q,q);
  arma::mat B1 = arma::eye(m0,m0);
  // double cond_no, cond_no_1, cond_no_2, cond_no_3;
  
  if (conjugate_gradient){
    arma::mat temp, first, second, third;
    temp =  static_cast<double>(N)/static_cast<double>(m0)*K0mat.rows(inds)*K0mat.rows(inds) + phi_sub*phi_sub.t() ;
    B1 = (arma::chol(nearPD(arma::inv(nearPD(temp))))).t();
    arma::mat phit_phi = phi.t() * phi;
    B2 = (arma::chol(nearPD(arma::inv(nearPD(phit_phi))))).t();
    
    arma::mat B1t_K0t = B1.t() * K0mat.t();
    arma::mat B1t_phisub = B1.t() * phi_sub;
    arma::mat phi_B2 = phi * B2;
    
    third = arma::eye(q,q); // B2.t() * phit_phi * B2;
    first = B1t_K0t*B1t_K0t.t() + B1t_phisub*B1t_phisub.t(); // B1.t() * (K0mat.t() * K0mat + phi_sub * phi_sub.t()) * B1;
    second = B1t_K0t * phi_B2; // B1.t() * K0mat.t()  * phi * B2;
    
    A = arma::join_horiz(first, second);
    A = arma::join_vert(A,arma::join_horiz(second.t(),third));
    
    arma::vec f;
    
    for (unsigned int i = 0; i<N_expectations; i++){
      f = integrands.col(i);
      b.col(i) = arma::join_vert( B1t_K0t*f, phi_B2.t()*f);// arma::join_vert( B1.t()*K0mat.t()*f, B2.t()*phi.t()*f);
    }
    
  } else{
    arma::mat first, second, third;
    
    third = phi.t() * phi;
    first = K0mat.t() * K0mat + phi_sub * phi_sub.t();
    second = K0mat.t()  * phi;
    
    A = arma::join_horiz(first, second);
    A = arma::join_vert(A,arma::join_horiz(second.t(),third));
    
    arma::vec f;
    
    for (unsigned int i = 0; i<N_expectations; i++){
      f = integrands.col(i);
      b.col(i) = arma::join_vert(K0mat.t()*f, phi.t()*f);
    }
  }
  
  double cond_no = arma::cond(A);
  
  return ( Rcpp::List::create(Rcpp::Named("A") = A, Rcpp::Named("b") = b, Rcpp::Named("B2") = B2, Rcpp::Named("B1") = B1, Rcpp::Named("m0") = m0, Rcpp::Named("cond_no") = cond_no,
                                          Rcpp::Named("phi") = phi, Rcpp::Named("K0") = K0mat, Rcpp::Named("ny_inds") = inds + 1) );
}



// This internal function is used to performing approximate semi-exact control functionals when samples ARE split for estimation and evaluation.
// There is an R wrapper that combines this function and the version of aSECF for when samples are not split.
// [[Rcpp::export]]
Rcpp::List aSECF_unbiased_cpp_prep(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, arma::uvec est_inds, Rcpp::Function getX, Rcpp::Function aSECF_mse_linsolve, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<arma::vec> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue, bool conjugate_gradient = true, double reltol = 0.01, bool diagnostics = false){
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  arma::uvec all_indx = arma::linspace<arma::uvec >(0, N-1,N);
  est_inds = est_inds - 1;
  unsigned int N_est = est_inds.n_rows;
  unsigned int m0;
  arma::uvec inds_ny, inds_ny_withinfull;
  if (nystrom_inds.isNull()){
    m0 = std::ceil( std::sqrt( static_cast<double>(N) ) );
    if (m0>N_est){
      throw(Rcpp::exception("The (defaulted) number of nystrom_inds is larger than the estimation sample size."));
    }
    inds_ny = Rcpp::RcppArmadillo::sample(arma::linspace<arma::uvec>(0, N_est-1, N_est), m0, false);
    inds_ny_withinfull = est_inds.rows(inds_ny);
  } else{
    inds_ny_withinfull = Rcpp::as<arma::uvec>(nystrom_inds) -1;
    m0 = inds_ny_withinfull.n_rows;
    inds_ny = arma::uvec(m0);
    arma::uvec overlap = arma::intersect(est_inds,inds_ny_withinfull);
    if (overlap.n_rows!=m0){
      throw(Rcpp::exception("The nystrom_inds are not a unique subset of the estimation set est_inds."));
    }
    arma::uvec temp;
    for (unsigned int kk=0; kk<m0; kk++){
      inds_ny(kk) = arma::as_scalar(find(est_inds == inds_ny_withinfull(kk)));
    }
  }
  
  unsigned int N_eval = N - N_est;
  
  arma::uvec eval_inds = all_indx;
  // eval_inds.shed_rows(est_inds);
  arma::uvec est_inds2 = sort(est_inds,"descend");
  for (unsigned int zz = 0; zz<N_est; zz++){
    eval_inds.shed_row(est_inds2(zz));
  }
  
  Rcpp::IntegerVector temp_inds = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(inds_ny + 1));
  Rcpp::Nullable<Rcpp::IntegerVector> inds_nullable(temp_inds);
  
  Rcpp::IntegerVector temp_inds2 = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(inds_ny_withinfull + 1));
  Rcpp::Nullable<Rcpp::IntegerVector> inds_nullable_withinfull(temp_inds2);
  
  arma::mat K0mat = getK0(samples, derivatives, R_NilValue, steinOrder, kernel_function, sigma, inds_nullable_withinfull, R_NilValue);
  
  arma::mat K0_est = K0mat.rows(est_inds);
  arma::mat K0_eval = K0mat.rows(eval_inds);
  
  arma::mat integrands_e = integrands.rows(est_inds);
  arma::mat samples_e = samples.rows(est_inds);
  arma::mat derivatives_e = derivatives.rows(est_inds);
  
  Rcpp::List ab_tilde = aSECF_mse_linsolve(integrands_e,samples_e,derivatives_e,polyorder, R_NilValue, R_NilValue, R_NilValue, K0_est, subset, inds_nullable, conjugate_gradient, reltol);
  
  arma::mat phi_eval = Phi_fn_cpp(samples.rows(eval_inds),derivatives.rows(eval_inds),getX,polyorder,subset);
  
  arma::uvec beta_rows = arma::linspace<arma::uvec>(m0, m0 + phi_eval.n_cols - 1, phi_eval.n_cols);
  arma::uvec a_rows = arma::linspace<arma::uvec>(0, m0-1, m0);
  
  Rcpp::List ab_tilde_curr;
  Rcpp::NumericVector res;
  arma::vec expectation(N_expectations);
  arma::vec res2;
  arma::mat beta(phi_eval.n_cols,N_expectations);
  arma::mat a(m0,N_expectations);
  arma::mat f_hat(N_eval,N_expectations);
  arma::mat f_eval = integrands.rows(eval_inds);
  for (unsigned int i = 0; i<N_expectations; i++){
    ab_tilde_curr = ab_tilde[i];
    res = ab_tilde_curr["sol"];
    res2 = Rcpp::as<arma::vec>(res);
    
    beta.col(i) = res2.rows(beta_rows);
    a.col(i) = res2.rows(a_rows);
    
    f_hat.col(i) = phi_eval * beta.col(i) + K0_eval * a.col(i);
    
    expectation(i) =  1/(N - N_est)*arma::sum(f_eval.col(i) - f_hat.col(i)) + arma::as_scalar(beta(0,i));
  }
  
  
  if (conjugate_gradient){
    arma::uvec iter(N_expectations);
    double cond_no = ab_tilde_curr["cond_no"]; // The same across iterations
    for (unsigned int i = 0; i<N_expectations; i++){
      ab_tilde_curr = ab_tilde[i];
      iter(i) = Rcpp::as<unsigned int>(ab_tilde_curr["iter"]);
    }
    if (!diagnostics){
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat,
                                  Rcpp::Named("iter") = iter, Rcpp::Named("cond_no") = cond_no) );
    } else{
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat,
                                  Rcpp::Named("iter") = iter, Rcpp::Named("cond_no") = cond_no,
                                  Rcpp::Named("a") = a, Rcpp::Named("b") = beta,
                                  Rcpp::Named("ny_inds") = inds_ny_withinfull + 1) );
    }
  } else {
    if (!diagnostics){
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat) );
    } else{
      return ( Rcpp::List::create(Rcpp::Named("expectation") = expectation, Rcpp::Named("f_true") = f_eval, Rcpp::Named("f_hat") = f_hat,
                                  Rcpp::Named("a") = a, Rcpp::Named("b") = beta,
                                  Rcpp::Named("ny_inds") = inds_ny_withinfull + 1) );
    }
  }
  
}

// The internal function for performing cross-validation to determine the mean square predictive error using a list of sigma (kernel tuning parameter) values or
// a list of K0 matrices (which can correspond to different sigma, steinOrder and kernel_function values) in semi-exact control functionals.
// The final estimation is performed in an R wrapper for this function.
// If the K0 matrix is not given, this procedure involves recalculating the K0 matrix once for each selected tuning parameter (there can be multiple tuning parameters selected if there are multiple integrands).
// This double-up on computation can be avoided by using the K0 argument.
// [[Rcpp::export]]
arma::mat aSECF_crossval_cpp(const arma::mat & integrands, const arma::mat & samples, const arma::mat & derivatives, Rcpp::Function getX, Rcpp::Function aSECF_mse_linsolve, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<Rcpp::List> sigma = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, Rcpp::Nullable<unsigned int> folds = R_NilValue, bool conjugate_gradient = true, double reltol = 0.01, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds = R_NilValue){
  
  unsigned int N_expectations = integrands.n_cols;
  
  unsigned int n_sigma;
  
  Rcpp::List sig_list;
  Rcpp::List K0_list;
  
  if (sigma.isNotNull()) {
    sig_list = Rcpp::as<Rcpp::List>(sigma);
    n_sigma = sig_list.length();
  } else{
    throw(Rcpp::exception("Must give sigma."));
  }
  
  arma::mat mse(N_expectations,n_sigma);
  
  for (unsigned int j=0; j<n_sigma; j++){
    //Rcpp::Rcout << "Iteration " << j << std::endl;
    Rcpp::NumericVector sigs = sig_list[j];
    Rcpp::Nullable<Rcpp::NumericVector> sigs1(sigs);
    mse.col(j) = aSECF_mse_cpp(integrands, samples, derivatives, getX, aSECF_mse_linsolve, polyorder, steinOrder, kernel_function, sigs1, subset, folds, conjugate_gradient, reltol, est_inds);
  }
  return ( mse );
}


// An internal function used by CF_crossval_cpp. Given a single kernel and tuning parameter, this function uses (folds)-fold cross-validation to get the approximate mean square predictive error using the fitted gaussian process models from different estimation sets.
arma::vec CF_mse_cpp(arma::mat integrands, arma::mat samples, arma::mat derivatives, Rcpp::Nullable<unsigned int> steinOrder, Rcpp::Nullable<Rcpp::String> kernel_function, Rcpp::Nullable<Rcpp::NumericVector> sigma, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0, Rcpp::Nullable<unsigned int> folds, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds, bool one_in_denom){
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  unsigned int n_fold;
  if (folds.isNotNull()) {
    n_fold = Rcpp::as<unsigned int>(folds);
  } else{
    n_fold = 5;
  }
  
  arma::uvec inds_est;
  if (est_inds.isNotNull()) {
    inds_est = Rcpp::as<arma::uvec>(est_inds) - 1;
    integrands = integrands.rows(inds_est);
    samples = samples.rows(inds_est);
    derivatives = derivatives.rows(inds_est);
    N = inds_est.n_rows;
  }
  
  
  arma::mat K0mat;
  if (K0.isNotNull()) {
    K0mat = Rcpp::as<arma::mat>(K0);
    if (est_inds.isNotNull()) {
      K0mat = K0mat.submat(inds_est,inds_est);
    }
  } else {
    Rcpp::NumericVector sigs = Rcpp::as<Rcpp::NumericVector>(sigma);
    Rcpp::Nullable<arma::vec> sigs1(wrap(sigs));
    K0mat = getK0(samples, derivatives, R_NilValue, steinOrder, kernel_function, sigs1);
  }
  
  K0mat = nearPD(K0mat);
  
  arma::uvec resampled = Rcpp::RcppArmadillo::sample(arma::linspace<arma::uvec>(0, N-1, N), N, false);
  
  arma::mat crossval_mse(N_expectations,n_fold);
  
  unsigned int num_each = std::ceil( static_cast<double>(N)/ static_cast<double>(n_fold) );
  unsigned int num_fewer = boost::math::iround(static_cast<double>(num_each)*static_cast<double>(n_fold) - static_cast<double>(N));
  
  
  arma::uvec holdout_indx, holdout_indx2, keep_indx;
  arma::uvec all_indx = arma::linspace<arma::uvec >(0, N-1,N);
  arma::mat K0inv, K0_curr, K0_holdout, f_curr_prep, f_holdout_prep;
  arma::vec K0inv_ones, weights, f_curr, f_holdout, K0inv_f, f_hat;
  double temp_1, mu;
  unsigned int N_curr, N_holdout;
  unsigned int ii_sum = 0;
  for (unsigned int ii = 0; ii<n_fold; ii++){
    if (static_cast<double>(ii)<=(static_cast<double>(num_fewer)-1.0)){
      holdout_indx = arma::linspace<arma::uvec>( ii_sum, ii_sum + num_each - 2, num_each-1);
      ii_sum = ii_sum + num_each - 1;
    } else{
      holdout_indx = arma::linspace<arma::uvec>( ii_sum, ii_sum + num_each - 1, num_each);
      ii_sum = ii_sum + num_each;
    }
    holdout_indx = resampled.rows(holdout_indx.rows(find(holdout_indx<=(N-1))));
    N_holdout = holdout_indx.n_rows;
    
    keep_indx = all_indx;
    // keep_indx.shed_rows(holdout_indx);
    holdout_indx2 = sort(holdout_indx,"descend");
    for (unsigned int zz = 0; zz<N_holdout; zz++){
      keep_indx.shed_row(holdout_indx2(zz));
    }
    
    K0_curr = K0mat.submat(keep_indx,keep_indx);
    K0_holdout = K0mat.submat(holdout_indx,keep_indx);
    
    N_curr = keep_indx.n_rows;
    
    K0_curr = nearPD(K0_curr);
    K0inv = arma::inv(K0_curr);
    K0inv_ones = K0inv*arma::ones<arma::vec>(N_curr);
    if (one_in_denom){
      temp_1 = arma::as_scalar(1.0 + arma::ones<arma::rowvec>(N_curr) * K0inv_ones);
    } else{
      temp_1 = arma::as_scalar(arma::ones<arma::rowvec>(N_curr) * K0inv_ones);
    }
    
    for (unsigned int i = 0; i<N_expectations; i++){      
      f_curr_prep = integrands.rows(keep_indx);
      f_curr = f_curr_prep.col(i);
      f_holdout_prep = integrands.rows(holdout_indx);
      f_holdout = f_holdout_prep.col(i);
      
      K0inv_f = K0inv*f_curr;
      mu = arma::as_scalar(arma::ones<arma::rowvec>(N_curr) *K0inv_f)/temp_1;
      f_hat = K0_holdout * K0inv_f + (1.0 - K0_holdout * K0inv_ones)*mu;
      crossval_mse(i,ii) = pow(norm(f_hat-f_holdout,2),2);
    }
  }
  
  arma::vec mse = arma::mean(crossval_mse,1);
  
  // return ( Rcpp::List::create(Rcpp::Named("mse") = mse) );
  return( mse );
}



// An internal function used by SECF_crossval_cpp. Given a single kernel and tuning parameter, this function uses (folds)-fold cross-validation to get the approximate mean square predictive error using the fitted gaussian process models from different estimation sets.
arma::vec SECF_mse_cpp(arma::mat integrands, arma::mat samples, arma::mat derivatives, Rcpp::Function getX, Rcpp::Nullable<unsigned int> polyorder, Rcpp::Nullable<unsigned int> steinOrder, Rcpp::Nullable<Rcpp::String> kernel_function, Rcpp::Nullable<Rcpp::NumericVector> sigma, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0, Rcpp::Nullable<Rcpp::IntegerVector> subset, Rcpp::Nullable<unsigned int> folds, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds){
  
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  unsigned int n_fold;
  if (folds.isNotNull()) {
    n_fold = Rcpp::as<unsigned int>(folds);
  } else{
    n_fold = 5;
  }
  
  arma::uvec inds_est;
  if (est_inds.isNotNull()) {
    inds_est = Rcpp::as<arma::uvec>(est_inds) - 1;
    integrands = integrands.rows(inds_est);
    samples = samples.rows(inds_est);
    derivatives = derivatives.rows(inds_est);
    N = inds_est.n_rows;
  }
  
  arma::mat phi = Phi_fn_cpp(samples,derivatives,getX,polyorder,subset);
  
  arma::mat K0mat;
  if (K0.isNotNull()) {
    K0mat = Rcpp::as<arma::mat>(K0);
    if (est_inds.isNotNull()) {
      K0mat = K0mat.submat(inds_est,inds_est);
    }
  } else {
    Rcpp::NumericVector sigs = Rcpp::as<Rcpp::NumericVector>(sigma);
    Rcpp::Nullable<arma::vec> sigs1(wrap(sigs));
    K0mat = getK0(samples, derivatives, R_NilValue, steinOrder, kernel_function, sigs1);
  }
  
  K0mat = nearPD(K0mat);
  
  arma::uvec resampled = Rcpp::RcppArmadillo::sample(arma::linspace<arma::uvec>(0, N-1, N), N, false);
  
  arma::mat crossval_mse(N_expectations,n_fold);
  
  unsigned int num_each = std::ceil( static_cast<double>(N)/ static_cast<double>(n_fold) );
  unsigned int num_fewer = boost::math::iround(static_cast<double>(num_each)*static_cast<double>(n_fold) - static_cast<double>(N));
  
  arma::uvec holdout_indx, holdout_indx2, keep_indx;
  arma::uvec all_indx = arma::linspace<arma::uvec >(0, N-1,N);
  arma::mat K0inv, K0_curr, K0_holdout, phi_curr, phi_holdout, K0inv_phi, f_curr_prep, f_holdout_prep;
  arma::vec beta, a, weights, f_curr, f_holdout, K0inv_f, f_hat, mu;
  unsigned int N_holdout;
  unsigned int ii_sum = 0;
  for (unsigned int ii = 0; ii<n_fold; ii++){
    if (static_cast<double>(ii)<=(static_cast<double>(num_fewer)-1.0)){
      holdout_indx = arma::linspace<arma::uvec>( ii_sum, ii_sum + num_each - 2, num_each-1);
      ii_sum = ii_sum + num_each - 1;
    } else{
      holdout_indx = arma::linspace<arma::uvec>( ii_sum, ii_sum + num_each - 1, num_each);
      ii_sum = ii_sum + num_each;
    }
    holdout_indx = resampled.rows(holdout_indx.rows(find(holdout_indx<=(N-1))));
    
    N_holdout = holdout_indx.n_rows;
    
    keep_indx = all_indx;
    holdout_indx2 = sort(holdout_indx,"descend");
    for (unsigned int zz = 0; zz<N_holdout; zz++){
      keep_indx.shed_row(holdout_indx2(zz));
    }
    
    K0_curr = K0mat.submat(keep_indx,keep_indx);
    K0_holdout = K0mat.submat(holdout_indx,keep_indx);
    
    K0_curr = nearPD(K0_curr);
    phi_curr = phi.rows(keep_indx);
    phi_holdout = phi.rows(holdout_indx);
    
    K0inv = arma::inv(K0_curr);
    K0inv_phi = K0inv*phi_curr;
    
    arma::mat weights_prep = arma::solve( phi_curr.t() * K0inv_phi, K0inv_phi.t() );
    
    for (unsigned int i = 0; i<N_expectations; i++){      
      f_curr_prep = integrands.rows(keep_indx);
      f_curr = f_curr_prep.col(i);
      f_holdout_prep = integrands.rows(holdout_indx);
      f_holdout = f_holdout_prep.col(i);
      
      beta = weights_prep * f_curr;
      
      mu = phi_curr * beta;
      a = K0inv*(f_curr-mu);
      f_hat = phi_holdout * beta + K0_holdout * a;
      crossval_mse(i,ii) = pow(norm(f_hat-f_holdout,2),2);
    }
  }
  
  arma::vec mse = arma::mean(crossval_mse,1);
  
  // return ( Rcpp::List::create(Rcpp::Named("mse") = mse) );
  return( mse );
}



// An internal function used by aSECF_crossval_cpp. Given a single kernel and tuning parameter, this function uses (folds)-fold cross-validation to get the approximate mean square predictive error using the fitted gaussian process models from different estimation sets.
arma::vec aSECF_mse_cpp(arma::mat integrands, arma::mat samples, arma::mat derivatives, Rcpp::Function getX, Rcpp::Function aSECF_mse_linsolve, Rcpp::Nullable<unsigned int> polyorder, Rcpp::Nullable<unsigned int> steinOrder, Rcpp::Nullable<Rcpp::String> kernel_function, Rcpp::Nullable<Rcpp::NumericVector> sigma, Rcpp::Nullable<Rcpp::IntegerVector> subset, Rcpp::Nullable<unsigned int> folds, bool conjugate_gradient, double reltol, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds){
  
  unsigned int N = samples.n_rows;
  unsigned int N_expectations = integrands.n_cols;
  
  unsigned int n_fold;
  if (folds.isNotNull()) {
    n_fold = Rcpp::as<unsigned int>(folds);
  } else{
    n_fold = 5;
  }
  
  
  arma::uvec inds_est;
  if (est_inds.isNotNull()) {
    inds_est = Rcpp::as<arma::uvec>(est_inds) - 1;
    integrands = integrands.rows(inds_est);
    samples = samples.rows(inds_est);
    derivatives = derivatives.rows(inds_est);
    N = inds_est.n_rows;
  }
  
  unsigned int m0 = std::ceil( std::sqrt( static_cast<double>(N) ) );
  arma::uvec inds_ny, inds_ny_withinfull; // Not giving inds_nystrom as can't really do cross-validation with a set set of inds. Don't get benefit of doing K0 once either.
  
  arma::uvec resampled = Rcpp::RcppArmadillo::sample(arma::linspace<arma::uvec>(0, N-1, N), N, false);
  
  arma::mat crossval_mse(N_expectations,n_fold);
  
  unsigned int num_each = std::ceil( static_cast<double>(N)/ static_cast<double>(n_fold) );
  unsigned int num_fewer = boost::math::iround(static_cast<double>(num_each)*static_cast<double>(n_fold) - static_cast<double>(N));
  
  arma::uvec holdout_indx, holdout_indx2, keep_indx;
  arma::uvec all_indx = arma::linspace<arma::uvec >(0, N-1,N);
  arma::mat K0_curr, K0_holdout, phi_holdout, K0inv_phi, f_curr_prep, f_holdout_prep;
  arma::vec beta, a, weights, f_curr, f_holdout, K0inv_f, f_hat, mu;
  unsigned int N_holdout,N_estim;
  unsigned int ii_sum = 0;
  for (unsigned int ii = 0; ii<n_fold; ii++){
    if (static_cast<double>(ii)<=(static_cast<double>(num_fewer)-1.0)){
      holdout_indx = arma::linspace<arma::uvec>( ii_sum, ii_sum + num_each - 2, num_each-1);
      ii_sum = ii_sum + num_each - 1;
    } else{
      holdout_indx = arma::linspace<arma::uvec>( ii_sum, ii_sum + num_each - 1, num_each);
      ii_sum = ii_sum + num_each;
    }
    holdout_indx = resampled.rows(holdout_indx.rows(find(holdout_indx<=(N-1))));
    
    N_holdout = holdout_indx.n_rows;
    N_estim = N - N_holdout;
    
    keep_indx = all_indx;
    holdout_indx2 = sort(holdout_indx,"descend");
    for (unsigned int zz = 0; zz<N_holdout; zz++){
      keep_indx.shed_row(holdout_indx2(zz));
    }
    
    inds_ny = Rcpp::RcppArmadillo::sample(arma::linspace<arma::uvec>(0, N_estim-1, N_estim), m0, false);
    inds_ny_withinfull = keep_indx.rows(inds_ny);
    
    Rcpp::IntegerVector temp_inds = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(inds_ny + 1));
    Rcpp::Nullable<Rcpp::IntegerVector> inds_nullable(temp_inds);
    
    Rcpp::IntegerVector temp_inds2 = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(inds_ny_withinfull + 1));
    Rcpp::Nullable<Rcpp::IntegerVector> inds_nullable_withinfull(temp_inds2);
    
    Rcpp::NumericVector sigs = Rcpp::as<Rcpp::NumericVector>(sigma);
    Rcpp::Nullable<arma::vec> sigs1(wrap(sigs));
    arma::mat K0mat = getK0(samples, derivatives, R_NilValue, steinOrder, kernel_function, sigs1, inds_nullable_withinfull, R_NilValue);
    
    arma::mat K0_curr = K0mat.rows(keep_indx);
    arma::mat K0_holdout = K0mat.rows(holdout_indx);
    
    arma::mat integrands_e = integrands.rows(keep_indx);
    arma::mat samples_e = samples.rows(keep_indx);
    arma::mat derivatives_e = derivatives.rows(keep_indx);
    
    Rcpp::List ab_tilde = aSECF_mse_linsolve(integrands_e,samples_e,derivatives_e,polyorder, R_NilValue, R_NilValue, R_NilValue, K0_curr, subset, inds_nullable, conjugate_gradient, reltol);
    
    arma::mat phi_holdout = Phi_fn_cpp(samples.rows(holdout_indx),derivatives.rows(holdout_indx),getX,polyorder,subset);
    
    arma::uvec beta_rows = arma::linspace<arma::uvec>(m0, m0 + phi_holdout.n_cols - 1, phi_holdout.n_cols);
    arma::uvec a_rows = arma::linspace<arma::uvec>(0, m0-1, m0);
    
    Rcpp::List ab_tilde_curr;
    Rcpp::NumericVector res;
    arma::vec res2;
    f_holdout_prep = integrands.rows(holdout_indx);
    for (unsigned int i = 0; i<N_expectations; i++){
      ab_tilde_curr = ab_tilde[i];
      res = ab_tilde_curr["sol"];
      res2 = Rcpp::as<arma::vec>(res);
      
      f_holdout = f_holdout_prep.col(i);
      
      arma::vec beta = res2.rows(beta_rows);
      arma::vec a = res2.rows(a_rows);
      
      f_hat = phi_holdout * beta + K0_holdout * a;
      crossval_mse(i,ii) = pow(norm(f_hat-f_holdout,2),2);
    }
  }
  
  arma::vec mse = arma::mean(crossval_mse,1);
  
  return( mse );
}


// An internal function for getting the matrix of square norms (which may be given or may require calculating).
arma::mat getSqNorm(const arma::mat & samples, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds, const Rcpp::Nullable<Rcpp::NumericMatrix> & Z){
  unsigned int N = samples.n_rows;
  
  // Checking if this has already been calculated (e.g. in cross-validation)
  if (Z.isNotNull()){
    if (nystrom_inds.isNull()){
      return ( Rcpp::as<arma::mat>(Z) );
    } else {
      arma::mat z = Rcpp::as<arma::mat>(Z);
      
      arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) -1;
      unsigned int m0 = inds.n_rows;
      
      if ((z.n_rows == N) & (z.n_cols == N)){
        z = z.cols(inds);
      } else if (z.n_cols != m0){
        throw(Rcpp::exception("If z is given, then nystrom_inds must also be given. The number of columns in z match either the number of samples or the length of the indices."));
      }
      return (z);
    }
  }
  
  arma::mat z = squareNorm(samples, nystrom_inds);
  return(z);
  
}

// An internal function for getting the K0 matrix (which may be given or may require calculating). This is where default kernel specifications are set.
arma::mat getK0(const arma::mat & samples, const arma::mat & derivatives, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0, Rcpp::Nullable<unsigned int> steinOrder, Rcpp::Nullable<Rcpp::String> kernel_function, Rcpp::Nullable<arma::vec> sigma, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds, const Rcpp::Nullable<Rcpp::NumericMatrix> & Z){
  arma::mat K0mat;
  std::string kernel_fun;
  arma::vec sigmas;
  unsigned int steinOrd;
  if (K0.isNotNull()) {
    K0mat = Rcpp::as<arma::mat>(K0);
  } else {
    if (kernel_function.isNotNull()) {
      kernel_fun = Rcpp::as<string>(kernel_function);
    } else{
      kernel_fun = "gaussian";
    }
    if (sigma.isNotNull()){
      sigmas = Rcpp::as<arma::vec>(sigma);
    } else{
      sigmas.ones(1);
      sigmas[0] = medianTune(samples);
    }
    if (steinOrder.isNotNull()){
      steinOrd = Rcpp::as<unsigned int>(steinOrder);
    } else{
      steinOrd = 1;
    }
    K0mat = K0_fn(samples,derivatives,sigmas,steinOrd,kernel_fun,Z,nystrom_inds);
  }
  return(K0mat);
}




