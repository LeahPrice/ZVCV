#define BOOST_DISABLE_ASSERTS

#ifndef ZVCV_KERNELS_H
#define ZVCV_KERNELS_H 1.0

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <cmath>

#include <boost/math/special_functions/bessel.hpp>

#include <RcppArmadillo.h>

//#include <algorithm>

// The stein operator applied to each kernel (k_0 for the specified kernel)

arma::mat gaussian_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, double kernel_params, std::string kernel_function,  const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue);

arma::mat matern_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, arma::vec kernel_params, std::string kernel_function,  const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue);

arma::mat RQ_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, double kernel_params, std::string kernel_function,  const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue);

arma::mat product_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, arma::vec kernel_params, std::string kernel_function,  const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue);

arma::mat prodsim_k(unsigned int steinOrder, const arma::mat & samples, const arma::mat & derivatives, arma::vec kernel_params, std::string kernel_function,  const arma::mat & z, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue);

#endif
