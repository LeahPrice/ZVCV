#include <kernels.h>

using namespace std;

arma::mat gaussian_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, double kernel_params, std::string kernel_function, const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds){
  
  unsigned int N = samples.n_rows;
  double d = static_cast<double>(samples.n_cols);
  
  // arma::mat z = getSqNorm(samples,nystrom_inds); // getting squared norms.
  arma::vec phi_z(4), x, y, ux, uy; // initialising derivatives wrt z
  
  double myvar = pow(kernel_params,2.0);
  
  if (nystrom_inds.isNull() & (steinOrder==2)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::mat K0(N,N,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = pow(-1.0/myvar,kk)*exp(-z(i,j)/myvar);
        }
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        K0(i,j) = 16*pow(z(i,j),2)*phi_z(3) + 16*(2+d)*z(i,j)*phi_z(2) + 4*(2+d)*d*phi_z(1) +
          4*(2*z(i,j)*phi_z(2) + (2+d)*phi_z(1))*arma::dot(ux-uy,x-y) - 4*phi_z(1)*arma::dot(ux,x-y)*arma::dot(uy,x-y) -
          2*phi_z(0)*arma::dot(ux,uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNull() & (steinOrder==1)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::mat K0(N,N,arma::fill::zeros);
    double kn;
    arma::vec dx1_k, dy1_k;
    arma::mat dx1dy1_k;
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = pow(-1.0/myvar,kk)*exp(-z(i,j)/myvar);
        }
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        kn = exp(-z(i,j)/myvar);
        
        dx1_k = 2*phi_z(0)*(x-y);
        dy1_k = -2*phi_z(0)*(x-y);
        dx1dy1_k = -4*phi_z(1)*(x-y)*(x-y).t() - 2*phi_z(0)*arma::ones<arma::mat>(d,d);
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
    
  } else if (nystrom_inds.isNotNull() & (steinOrder==2)){ //subset-based approach for Nystrom approximation
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = pow(-1.0/myvar,kk)*exp(-z(i,j)/myvar);
        }
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        
        K0(i,j) = 16*pow(z(i,j),2)*phi_z(3) + 16*(2+d)*z(i,j)*phi_z(2) + 4*(2+d)*d*phi_z(1) +
          4*(2*z(i,j)*phi_z(2) + (2+d)*phi_z(1))*arma::dot(ux-uy,x-y) - 4*phi_z(1)*arma::dot(ux,x-y)*arma::dot(uy,x-y) -
          2*phi_z(0)*arma::dot(ux,uy);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNotNull() & (steinOrder==1)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    double kn;
    arma::vec dx1_k, dy1_k;
    arma::mat dx1dy1_k;
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = pow(-1.0/myvar,kk)*exp(-z(i,j)/myvar);
        }
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        kn = exp(-z(i,j)/myvar);
        
        dx1_k = 2*phi_z(0)*(x-y);
        dy1_k = -2*phi_z(0)*(x-y);
        dx1dy1_k = -4*phi_z(1)*(x-y)*(x-y).t() - 2*phi_z(0)*arma::ones<arma::mat>(d,d);
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
      }
    }
    return(K0);
  }
  
  Rcpp::Rcout << "Error in gaussian_k" << std::endl;
  return ( std::numeric_limits<arma::mat>::quiet_NaN() );
  
}


arma::mat matern_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, arma::vec kernel_params, std::string kernel_function, const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds){
  
  unsigned int N = samples.n_rows;
  double d = static_cast<double>(samples.n_cols);
  
  // arma::mat z = getSqNorm(samples,nystrom_inds); // getting squared norms.
  //z = z + pow(10.0,-10.0);
  arma::vec phi_z(4), x, y, ux, uy; // initialising derivatives wrt z
  
  
  double lambda = kernel_params(0);
  double nu = kernel_params(1);
  
  double b = pow(2.0,1.0-nu)/exp(lgamma(nu));
  double c = sqrt(2.0*nu)/lambda;
  
  double kn_const = b*pow(c,nu);
  arma::vec phi_z_part1(4);
  
  double kk = 0.0;
  for (unsigned int k = 0; k<4; k++){
    kk += 1.0;
    phi_z_part1(k) = pow(-1.0/2.0,kk)*b*pow(c,nu+kk);
  }
  
  
  if (nystrom_inds.isNull() & (steinOrder==2)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::mat K0(N,N,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) =  phi_z_part1(k)*pow(z(i,j),(nu-kk)/2.0)*boost::math::cyl_bessel_k(nu-kk,c*sqrt(z(i,j)));
        }
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        K0(i,j) = 16*pow(z(i,j),2)*phi_z(3) + 16*(2+d)*z(i,j)*phi_z(2) + 4*(2+d)*d*phi_z(1) +
          4*(2*z(i,j)*phi_z(2) + (2+d)*phi_z(1))*arma::dot(ux-uy,x-y) - 4*phi_z(1)*arma::dot(ux,x-y)*arma::dot(uy,x-y) -
          2*phi_z(0)*arma::dot(ux,uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNull() & (steinOrder==1)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::mat K0(N,N,arma::fill::zeros);
    double kn;
    arma::vec dx1_k, dy1_k;
    arma::mat dx1dy1_k;
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = phi_z_part1(k)*pow(z(i,j),(nu-kk)/2.0)*boost::math::cyl_bessel_k(nu-kk,c*sqrt(z(i,j)));
        }
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        kn = kn_const*pow(z(i,j),nu/2.0)*boost::math::cyl_bessel_k(nu,c*sqrt(z(i,j)));
        
        dx1_k = 2*phi_z(0)*(x-y);
        dy1_k = -2*phi_z(0)*(x-y);
        dx1dy1_k = -4*phi_z(1)*(x-y)*(x-y).t() - 2*phi_z(0)*arma::ones<arma::mat>(d,d);
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
    
  } else if (nystrom_inds.isNotNull() & (steinOrder==2)){ //subset-based approach for Nystrom approximation
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = phi_z_part1(k)*pow(z(i,j),(nu-kk)/2.0)*boost::math::cyl_bessel_k(nu-kk,c*sqrt(z(i,j)));
        }
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        
        K0(i,j) = 16*pow(z(i,j),2)*phi_z(3) + 16*(2+d)*z(i,j)*phi_z(2) + 4*(2+d)*d*phi_z(1) +
          4*(2*z(i,j)*phi_z(2) + (2+d)*phi_z(1))*arma::dot(ux-uy,x-y) - 4*phi_z(1)*arma::dot(ux,x-y)*arma::dot(uy,x-y) -
          2*phi_z(0)*arma::dot(ux,uy);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNotNull() & (steinOrder==1)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    double kn;
    arma::vec dx1_k, dy1_k;
    arma::mat dx1dy1_k;
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = phi_z_part1(k)*pow(z(i,j),(nu-kk)/2.0)*boost::math::cyl_bessel_k(nu-kk,c*sqrt(z(i,j)));
        }
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        kn = kn_const*pow(z(i,j),nu/2.0)*boost::math::cyl_bessel_k(nu,c*sqrt(z(i,j)));
        
        dx1_k = 2*phi_z(0)*(x-y);
        dy1_k = -2*phi_z(0)*(x-y);
        dx1dy1_k = -4*phi_z(1)*(x-y)*(x-y).t() - 2*phi_z(0)*arma::ones<arma::mat>(d,d);
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
      }
    }
    return(K0);
  }
  
  Rcpp::Rcout << "Error in matern_k" << std::endl;
  return ( std::numeric_limits<arma::mat>::quiet_NaN() );
  
}


arma::mat RQ_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, double kernel_params, std::string kernel_function, const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds){
  
  unsigned int N = samples.n_rows;
  double d = static_cast<double>(samples.n_cols);
  
  //arma::mat z = getSqNorm(samples,nystrom_inds); // getting squared norms.
  arma::vec phi_z(4), x, y, ux, uy; // initialising derivatives wrt z
  
  double myprecision = pow(kernel_params,-2.0);
  
  arma::mat c = 1.0 + myprecision*z;
  
  arma::vec phi_z_part1(4);
  
  double kk = 0.0;
  double factorial = 1.0;
  for (unsigned int k = 0; k<4; k++){
    kk += 1.0;
    factorial *= kk;
    phi_z_part1(k) = pow(-1.0,kk)*factorial*pow(myprecision,kk);
  }
  
  
  if (nystrom_inds.isNull() & (steinOrder==2)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::mat K0(N,N,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = phi_z_part1(k)*pow(c(i,j),-kk-1.0);
        }
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        K0(i,j) = 16*pow(z(i,j),2)*phi_z(3) + 16*(2+d)*z(i,j)*phi_z(2) + 4*(2+d)*d*phi_z(1) +
          4*(2*z(i,j)*phi_z(2) + (2+d)*phi_z(1))*arma::dot(ux-uy,x-y) - 4*phi_z(1)*arma::dot(ux,x-y)*arma::dot(uy,x-y) -
          2*phi_z(0)*arma::dot(ux,uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNull() & (steinOrder==1)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::mat K0(N,N,arma::fill::zeros);
    double kn;
    arma::vec dx1_k, dy1_k;
    arma::mat dx1dy1_k;
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = phi_z_part1(k)*pow(c(i,j),-kk-1.0);
        }
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        kn = 1.0/c(i,j);
        
        dx1_k = 2*phi_z(0)*(x-y);
        dy1_k = -2*phi_z(0)*(x-y);
        dx1dy1_k = -4*phi_z(1)*(x-y)*(x-y).t() - 2*phi_z(0)*arma::ones<arma::mat>(d,d);
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
    
  } else if (nystrom_inds.isNotNull() & (steinOrder==2)){ //subset-based approach for Nystrom approximation
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = phi_z_part1(k)*pow(c(i,j),-kk-1.0);
        }
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        
        K0(i,j) = 16*pow(z(i,j),2)*phi_z(3) + 16*(2+d)*z(i,j)*phi_z(2) + 4*(2+d)*d*phi_z(1) +
          4*(2*z(i,j)*phi_z(2) + (2+d)*phi_z(1))*arma::dot(ux-uy,x-y) - 4*phi_z(1)*arma::dot(ux,x-y)*arma::dot(uy,x-y) -
          2*phi_z(0)*arma::dot(ux,uy);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNotNull() & (steinOrder==1)){
    unsigned int i, j; // indices for the loop over N
    unsigned int k; // for use in the derivatives
    double kk;
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    double kn;
    arma::vec dx1_k, dy1_k;
    arma::mat dx1dy1_k;
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        kk = 0.0;
        for (k = 0; k<4; k++){
          kk += 1.0;
          phi_z(k) = phi_z_part1(k)*pow(c(i,j),-kk-1.0);
        }
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        kn = 1.0/c(i,j);
        
        dx1_k = 2*phi_z(0)*(x-y);
        dy1_k = -2*phi_z(0)*(x-y);
        dx1dy1_k = -4*phi_z(1)*(x-y)*(x-y).t() - 2*phi_z(0)*arma::ones<arma::mat>(d,d);
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
      }
    }
    return(K0);
  }
  
  Rcpp::Rcout << "Error in RQ_k" << std::endl;
  return ( std::numeric_limits<arma::mat>::quiet_NaN() );
  
}


arma::mat product_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, arma::vec kernel_params, std::string kernel_function, const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds){
  
  unsigned int N = samples.n_rows;
  double d = static_cast<double>(samples.n_cols);
  
  //arma::mat z = getSqNorm(samples,nystrom_inds); // getting squared norms.
  arma::vec x, y, ux, uy; // initialising derivatives wrt z
  
  double a = kernel_params(0);
  double beta = -0.5*pow(kernel_params(1),-2.0);
  
  if (nystrom_inds.isNull() & (steinOrder==2)){
    double f, dx2_f, dy2_f, dx2dy2_f, g, dx2_g, dy2_g, dx2dy2_g, dx2dy2_k;
    arma::vec dx1_f, dy1_f, dx2dy1_f, dx1dy2_f, dx1_g, dy1_g, dx2dy1_g, dx1dy2_g, dx1dy2_k, dx2dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::mat K0(N,N,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        // First part of product, (1 + alpha_1*norm(x)^2 + alpha_1*norm(y)^2)^(-1)
        f = 1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)); // pow(1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)),-1.0);
        //Rcpp::Rcout << "f is " << f << std::endl;
        dx1_f = -2*a*x*pow(f,-2.0);
        dy1_f = -2*a*y*pow(f,-2.0);
        dx1dy1_f = 8*pow(a,2.0)*x*y.t()*pow(f,-3.0);
        dx2_f = -2*d*a*pow(f,-2.0) + 8*pow(a,2.0)*pow(norm(x,2),2.0)*pow(f,-3.0);
        dy2_f = -2*d*a*pow(f,-2.0) + 8*pow(a,2.0)*pow(norm(y,2),2.0)*pow(f,-3.0);
        dx2dy1_f = 8*d*pow(a,2.0)*y*pow(f,-3.0) - 48*pow(a,3.0)*y*pow(norm(x,2),2.0)*pow(f,-4.0);
        dx1dy2_f = 8*d*pow(a,2.0)*x*pow(f,-3.0) - 48*pow(a,3.0)*x*pow(norm(y,2),2.0)*pow(f,-4.0);
        dx2dy2_f = 8*pow(d,2.0)*pow(a,2.0)*pow(f,-3.0) -
          48*d*pow(a,3.0)*(pow(norm(x,2),2.0)+pow(norm(y,2),2.0))*pow(f,-4.0) +
          384*pow(a,4.0)*pow(norm(x,2),2.0)*pow(norm(y,2),2.0)*pow(f,-5.0);
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        dx2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dy2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dx2dy1_g = -8*z(i,j)*pow(beta,3.0)*g*(x-y) - 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx1dy2_g = 8*z(i,j)*pow(beta,3.0)*g*(x-y) + 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx2dy2_g = 16*pow(z(i,j),2.0)*pow(beta,4.0)*g + 16*(2+d)*z(i,j)*pow(beta,3.0)*g + 4*(2+d)*d*pow(beta,2.0)*g;
        
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        dx1dy2_k = 2*(dx1dy1_f * dy1_g) + 2*(dx1dy1_g * dy1_f) + (dx1_f * dy2_g) + (dx1_g * dy2_f) + pow(f,-1.0)*(dx1dy2_g) + g*(dx1dy2_f);
        dx2dy1_k = 2*(dx1dy1_f * dx1_g) + 2*(dx1dy1_g * dx1_f) + (dx2_f * dy1_g) + (dx2_g * dy1_f) + pow(f,-1.0)*(dx2dy1_g) + g*(dx2dy1_f);
        dx2dy2_k = 2*arma::dot(dx2dy1_f, dy1_g) +
          2*arma::dot(dx2dy1_g, dy1_f) +
          2*arma::dot(dx1_f, dx1dy2_g) +
          2*arma::dot(dx1_g, dx1dy2_f) +
          f*dx2dy2_g +
          g*dx2dy2_f +
          4*arma::dot(dx1dy1_f, dx1dy1_g) +
          2*(dx2_f * dy2_g);
        
        K0(i,j) = arma::as_scalar(ux.t()*dx1dy1_k*uy + ux.t()*dx1dy2_k + uy.t()*dx2dy1_k + dx2dy2_k);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNull() & (steinOrder==1)){
    double f, g, kn;
    arma::vec dx1_f, dy1_f, dx1_g, dy1_g, dx1_k, dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::mat K0(N,N,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        f = 1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)); // pow(1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)),-1.0);
        dx1_f = -2*a*x*pow(f,-2.0);
        dy1_f = -2*a*y*pow(f,-2.0);
        dx1dy1_f = 8*pow(a,2.0)*x*y.t()*pow(f,-3.0);
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        
        kn = pow(f,-1.0)*g;
        dx1_k = pow(f,-1.0)*dx1_g + g*dx1_f;
        dy1_k = pow(f,-1.0)*dy1_g + g*dy1_f;
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
    
  } else if (nystrom_inds.isNotNull() & (steinOrder==2)){ //subset-based approach for Nystrom approximation
    double f, dx2_f, dy2_f, dx2dy2_f, g, dx2_g, dy2_g, dx2dy2_g, dx2dy2_k;
    arma::vec dx1_f, dy1_f, dx2dy1_f, dx1dy2_f, dx1_g, dy1_g, dx2dy1_g, dx1dy2_g, dx1dy2_k, dx2dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        
        // First part of product, (1 + alpha_1*norm(x)^2 + alpha_1*norm(y)^2)^(-1)
        f = 1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)); // pow(1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)),-1.0);
        //Rcpp::Rcout << "f is " << f << std::endl;
        dx1_f = -2*a*x*pow(f,-2.0);
        dy1_f = -2*a*y*pow(f,-2.0);
        dx1dy1_f = 8*pow(a,2.0)*x*y.t()*pow(f,-3.0);
        dx2_f = -2*d*a*pow(f,-2.0) + 8*pow(a,2.0)*pow(norm(x,2),2.0)*pow(f,-3.0);
        dy2_f = -2*d*a*pow(f,-2.0) + 8*pow(a,2.0)*pow(norm(y,2),2.0)*pow(f,-3.0);
        dx2dy1_f = 8*d*pow(a,2.0)*y*pow(f,-3.0) - 48*pow(a,3.0)*y*pow(norm(x,2),2.0)*pow(f,-4.0);
        dx1dy2_f = 8*d*pow(a,2.0)*x*pow(f,-3.0) - 48*pow(a,3.0)*x*pow(norm(y,2),2.0)*pow(f,-4.0);
        dx2dy2_f = 8*pow(d,2.0)*pow(a,2.0)*pow(f,-3.0) -
          48*d*pow(a,3.0)*(pow(norm(x,2),2.0)+pow(norm(y,2),2.0))*pow(f,-4.0) +
          384*pow(a,4.0)*pow(norm(x,2),2.0)*pow(norm(y,2),2.0)*pow(f,-5.0);
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        dx2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dy2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dx2dy1_g = -8*z(i,j)*pow(beta,3.0)*g*(x-y) - 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx1dy2_g = 8*z(i,j)*pow(beta,3.0)*g*(x-y) + 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx2dy2_g = 16*pow(z(i,j),2.0)*pow(beta,4.0)*g + 16*(2+d)*z(i,j)*pow(beta,3.0)*g + 4*(2+d)*d*pow(beta,2.0)*g;
        
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        dx1dy2_k = 2*(dx1dy1_f * dy1_g) + 2*(dx1dy1_g * dy1_f) + (dx1_f * dy2_g) + (dx1_g * dy2_f) + pow(f,-1.0)*(dx1dy2_g) + g*(dx1dy2_f);
        dx2dy1_k = 2*(dx1dy1_f * dx1_g) + 2*(dx1dy1_g * dx1_f) + (dx2_f * dy1_g) + (dx2_g * dy1_f) + pow(f,-1.0)*(dx2dy1_g) + g*(dx2dy1_f);
        dx2dy2_k = 2*arma::dot(dx2dy1_f, dy1_g) +
          2*arma::dot(dx2dy1_g, dy1_f) +
          2*arma::dot(dx1_f, dx1dy2_g) +
          2*arma::dot(dx1_g, dx1dy2_f) +
          f*dx2dy2_g +
          g*dx2dy2_f +
          4*arma::dot(dx1dy1_f, dx1dy1_g) +
          2*(dx2_f * dy2_g);
        
        K0(i,j) = arma::as_scalar(ux.t()*dx1dy1_k*uy + ux.t()*dx1dy2_k + uy.t()*dx2dy1_k + dx2dy2_k);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNotNull() & (steinOrder==1)){
    double f, g, kn;
    arma::vec dx1_f, dy1_f, dx1_g, dy1_g, dx1_k, dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        
        f = 1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)); // pow(1+a*(pow(norm(x,2),2.0) + pow(norm(y,2),2.0)),-1.0);
        dx1_f = -2*a*x*pow(f,-2.0);
        dy1_f = -2*a*y*pow(f,-2.0);
        dx1dy1_f = 8*pow(a,2.0)*x*y.t()*pow(f,-3.0);
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        
        kn = pow(f,-1.0)*g;
        dx1_k = pow(f,-1.0)*dx1_g + g*dx1_f;
        dy1_k = pow(f,-1.0)*dy1_g + g*dy1_f;
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
      }
    }
    return(K0);
  }
  
  Rcpp::Rcout << "Error in product_k" << std::endl;
  return ( std::numeric_limits<arma::mat>::quiet_NaN() );
  
}


arma::mat prodsim_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, arma::vec kernel_params, std::string kernel_function, const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds){
  
  unsigned int N = samples.n_rows;
  double d = static_cast<double>(samples.n_cols);
  
  double a = kernel_params(0);
  double beta = -0.5*pow(kernel_params(1),-2.0);
  
  //arma::mat z = getSqNorm(samples,nystrom_inds); // getting squared norms.
  arma::vec x, y, ux, uy; // initialising derivatives wrt z
  
  if (nystrom_inds.isNull() & (steinOrder==2)){
    double f1, f2, f, dx2_f1, dx2_f, dy2_f2, dy2_f, dx2dy2_f, g, dx2_g, dy2_g, dx2dy2_g, dx2dy2_k;
    arma::vec dx1_f1, dx1_f, dy1_f2, dy1_f, dx2dy1_f, dx1dy2_f, dx1_g, dy1_g, dx2dy1_g, dx1dy2_g, dx1dy2_k, dx2dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::mat K0(N,N,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        // First part of product, (1 + alpha_1*norm(x)^2 + alpha_1*norm(y)^2)^(-1)
        f1 = 1+a*pow(norm(x,2),2.0);
        f2 = 1+a*pow(norm(y,2),2.0);
        f = f1*f2;
        dx1_f1 = -2*a*x*pow(f1,-2.0);
        dx1_f = dx1_f1*pow(f2,-1.0);
        dy1_f2 = -2*a*y*pow(f2,-2.0);
        dy1_f = dy1_f2*pow(f1,-1.0);
        dx1dy1_f = dx1_f*dy1_f.t();
        dx2_f1 = -2*d*a*pow(f1,-2.0) + 8*pow(a,2.0)*pow(norm(x,2),2.0)*pow(f1,-3.0);
        dx2_f = dx2_f1*pow(f2,-1.0);
        dy2_f2 = -2*d*a*pow(f2,-2.0) + 8*pow(a,2.0)*pow(norm(y,2),2.0)*pow(f2,-3.0);
        dy2_f = dy2_f2*pow(f1,-1.0);
        dx2dy1_f = dx2_f1*dy1_f2;
        dx1dy2_f = dx1_f1*dy2_f2;
        dx2dy2_f = dx2_f1*dy2_f2;
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        dx2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dy2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dx2dy1_g = -8*z(i,j)*pow(beta,3.0)*g*(x-y) - 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx1dy2_g = 8*z(i,j)*pow(beta,3.0)*g*(x-y) + 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx2dy2_g = 16*pow(z(i,j),2.0)*pow(beta,4.0)*g + 16*(2+d)*z(i,j)*pow(beta,3.0)*g + 4*(2+d)*d*pow(beta,2.0)*g;
        
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        dx1dy2_k = 2*(dx1dy1_f * dy1_g) + 2*(dx1dy1_g * dy1_f) + (dx1_f * dy2_g) + (dx1_g * dy2_f) + pow(f,-1.0)*(dx1dy2_g) + g*(dx1dy2_f);
        dx2dy1_k = 2*(dx1dy1_f * dx1_g) + 2*(dx1dy1_g * dx1_f) + (dx2_f * dy1_g) + (dx2_g * dy1_f) + pow(f,-1.0)*(dx2dy1_g) + g*(dx2dy1_f);
        dx2dy2_k = 2*arma::dot(dx2dy1_f, dy1_g) +
          2*arma::dot(dx2dy1_g, dy1_f) +
          2*arma::dot(dx1_f, dx1dy2_g) +
          2*arma::dot(dx1_g, dx1dy2_f) +
          f*dx2dy2_g +
          g*dx2dy2_f +
          4*arma::dot(dx1dy1_f, dx1dy1_g) +
          2*(dx2_f * dy2_g);
        
        K0(i,j) = arma::as_scalar(ux.t()*dx1dy1_k*uy + ux.t()*dx1dy2_k + uy.t()*dx2dy1_k + dx2dy2_k);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNull() & (steinOrder==1)){
    double f1, f2, f, g, kn;
    arma::vec dx1_f1, dx1_f, dy1_f2, dy1_f, dx1_g, dy1_g, dx1_k, dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::mat K0(N,N,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = i; j<N; j++){
        x = samples.row(i).t();
        y = samples.row(j).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(j).t();
        
        // First part of product, (1 + alpha_1*norm(x)^2 + alpha_1*norm(y)^2)^(-1)
        f1 = 1+a*pow(norm(x,2),2.0);
        f2 = 1+a*pow(norm(y,2),2.0);
        f = f1*f2;
        dx1_f1 = -2*a*x*pow(f1,-2.0);
        dx1_f = dx1_f1*pow(f2,-1.0);
        dy1_f2 = -2*a*y*pow(f2,-2.0);
        dy1_f = dy1_f2*pow(f1,-1.0);
        dx1dy1_f = dx1_f1*dy1_f2.t();
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        
        kn = pow(f,-1.0)*g;
        dx1_k = pow(f,-1.0)*dx1_g + g*dx1_f;
        dy1_k = pow(f,-1.0)*dy1_g + g*dy1_f;
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
        K0(j,i) = K0(i,j);
      }
    }
    return(K0);
    
  } else if (nystrom_inds.isNotNull() & (steinOrder==2)){ //subset-based approach for Nystrom approximation
    double f1, f2, f, dx2_f1, dx2_f, dy2_f2, dy2_f, dx2dy2_f, g, dx2_g, dy2_g, dx2dy2_g, dx2dy2_k;
    arma::vec dx1_f1, dx1_f, dy1_f2, dy1_f, dx2dy1_f, dx1dy2_f, dx1_g, dy1_g, dx2dy1_g, dx1dy2_g, dx1dy2_k, dx2dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        
        // First part of product, (1 + alpha_1*norm(x)^2 + alpha_1*norm(y)^2)^(-1)
        f1 = 1+a*pow(norm(x,2),2.0);
        f2 = 1+a*pow(norm(y,2),2.0);
        f = f1*f2;
        dx1_f1 = -2*a*x*pow(f1,-2.0);
        dx1_f = dx1_f1*pow(f2,-1.0);
        dy1_f2 = -2*a*y*pow(f2,-2.0);
        dy1_f = dy1_f2*pow(f1,-1.0);
        dx1dy1_f = dx1_f*dy1_f.t();
        dx2_f1 = -2*d*a*pow(f1,-2.0) + 8*pow(a,2.0)*pow(norm(x,2),2.0)*pow(f1,-3.0);
        dx2_f = dx2_f1*pow(f2,-1.0);
        dy2_f2 = -2*d*a*pow(f2,-2.0) + 8*pow(a,2.0)*pow(norm(y,2),2.0)*pow(f2,-3.0);
        dy2_f = dy2_f2*pow(f1,-1.0);
        dx2dy1_f = dx2_f1*dy1_f2;
        dx1dy2_f = dx1_f1*dy2_f2;
        dx2dy2_f = dx2_f1*dy2_f2;
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        dx2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dy2_g = 4*z(i,j)*pow(beta,2.0)*g + 2*d*beta*g;
        dx2dy1_g = -8*z(i,j)*pow(beta,3.0)*g*(x-y) - 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx1dy2_g = 8*z(i,j)*pow(beta,3.0)*g*(x-y) + 4*(2+d)*pow(beta,2.0)*g*(x-y);
        dx2dy2_g = 16*pow(z(i,j),2.0)*pow(beta,4.0)*g + 16*(2+d)*z(i,j)*pow(beta,3.0)*g + 4*(2+d)*d*pow(beta,2.0)*g;
        
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        dx1dy2_k = 2*(dx1dy1_f * dy1_g) + 2*(dx1dy1_g * dy1_f) + (dx1_f * dy2_g) + (dx1_g * dy2_f) + pow(f,-1.0)*(dx1dy2_g) + g*(dx1dy2_f);
        dx2dy1_k = 2*(dx1dy1_f * dx1_g) + 2*(dx1dy1_g * dx1_f) + (dx2_f * dy1_g) + (dx2_g * dy1_f) + pow(f,-1.0)*(dx2dy1_g) + g*(dx2dy1_f);
        dx2dy2_k = 2*arma::dot(dx2dy1_f, dy1_g) +
          2*arma::dot(dx2dy1_g, dy1_f) +
          2*arma::dot(dx1_f, dx1dy2_g) +
          2*arma::dot(dx1_g, dx1dy2_f) +
          f*dx2dy2_g +
          g*dx2dy2_f +
          4*arma::dot(dx1dy1_f, dx1dy1_g) +
          2*(dx2_f * dy2_g);
        
        K0(i,j) = arma::as_scalar(ux.t()*dx1dy1_k*uy + ux.t()*dx1dy2_k + uy.t()*dx2dy1_k + dx2dy2_k);
      }
    }
    return(K0);
  } else if (nystrom_inds.isNotNull() & (steinOrder==1)){
    double f1, f2, f, g, kn;
    arma::vec dx1_f1, dx1_f, dy1_f2, dy1_f, dx1_g, dy1_g, dx1_k, dy1_k;
    arma::mat dx1dy1_f, dx1dy1_g, dx1dy1_k;
    
    unsigned int i, j; // indices for the loop over N
    arma::uvec inds = Rcpp::as<arma::uvec>(nystrom_inds) - 1;
    unsigned int m0 = inds.n_rows;
    arma::mat K0(N,m0,arma::fill::zeros);
    for (i = 0; i<N; i++){
      for (j = 0; j<m0; j++){
        x = samples.row(i).t();
        y = samples.row(inds(j)).t();
        ux = derivatives.row(i).t();
        uy = derivatives.row(inds(j)).t();
        
        // First part of product, (1 + alpha_1*norm(x)^2 + alpha_1*norm(y)^2)^(-1)
        f1 = 1+a*pow(norm(x,2),2.0);
        f2 = 1+a*pow(norm(y,2),2.0);
        f = f1*f2;
        dx1_f1 = -2*a*x*pow(f1,-2.0);
        dx1_f = dx1_f1*pow(f2,-1.0);
        dy1_f2 = -2*a*y*pow(f2,-2.0);
        dy1_f = dy1_f2*pow(f1,-1.0);
        dx1dy1_f = dx1_f1*dy1_f2.t();
        
        // Second part of product, exp(beta * norm(x-y)^2)
        g = exp(beta*z(i,j));
        dx1_g = 2*beta*g*(x-y);
        dy1_g = -2*beta*g*(x-y);
        dx1dy1_g = -4*pow(beta,2.0)*g*(x-y)*(x-y).t() - 2*beta*g;
        
        kn = pow(f,-1.0)*g;
        dx1_k = pow(f,-1.0)*dx1_g + g*dx1_f;
        dy1_k = pow(f,-1.0)*dy1_g + g*dy1_f;
        dx1dy1_k = dy1_f * dx1_g.t() + dx1_f * dy1_g.t() + pow(f,-1.0)*dx1dy1_g + g*dx1dy1_f;
        
        K0(i,j) = arma::as_scalar(arma::trace(dx1dy1_k) + ux.t()*dy1_k + uy.t()*dx1_k + ux.t()*kn*uy);
      }
    }
    return(K0);
  }
  
  Rcpp::Rcout << "Error in prodsim_k" << std::endl;
  return ( std::numeric_limits<arma::mat>::quiet_NaN() );
  
}
