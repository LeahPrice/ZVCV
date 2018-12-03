#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <cmath>

using namespace std;

// getPoly() function callable from R via Rcpp::
// [[Rcpp::export]]
arma::mat getPoly(arma::mat samples, arma::mat derivatives, Rcpp::List combinations){//, unsigned long polyorder) {  
	
	unsigned long N = samples.n_rows;
	unsigned long d = samples.n_cols;
	
	derivatives = -0.5*derivatives;
	
	unsigned long polyorder = combinations.size();
	
	unsigned int mypol, k, j, jj;
	
	
	unsigned long num_monos, total_monos;
	arma::vec term;
	
	total_monos = 0;
	for (mypol = 0; mypol<polyorder; mypol++){
		arma::Mat<int> alpha = combinations[mypol];
		total_monos += alpha.n_cols;
	}
	
	arma::mat X(N,total_monos);
	unsigned long ii = 0;
	
	for (mypol = 0; mypol<polyorder; mypol++){
		arma::Mat<int> alpha = combinations[mypol];
		num_monos = alpha.n_cols;
		for (k = 0; k<num_monos; k++){
			term = arma::zeros<arma::vec>(N);
			// First order derivatives wrt parameter j
			for (j = 0; j<d; j++){
				arma::vec myprod = std::max(alpha(j,k),0)*pow(samples.col(j),alpha(j,k) - 1) % derivatives.col(j);
				for (jj=0; jj<d; jj++){
					if (jj!=j){myprod = myprod % pow(samples.col(jj),alpha(jj,k));}
				}
				term = term + myprod;
			}
			// Second order derivatives wrt parameter j
			for (j = 0; j<d; j++){
				arma::vec myprod = std::max(alpha(j,k)*(alpha(j,k)-1),0) * pow(samples.col(j),alpha(j,k) - 2);
				for (jj=0; jj<d; jj++){
					if (jj!=j){myprod = myprod % pow(samples.col(jj),alpha(jj,k));}
				}
				term = term + -0.5*myprod;
			}
			X.col(ii) = term;
			ii++;
		}
	}
	
	return X;
}

