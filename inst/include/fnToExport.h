#ifndef ZVCV_FNTOEXPORT_H
#define ZVCV_FNTOEXPORT_H 1.0

#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <cmath>

#include <boost/math/special_functions/round.hpp>
#include <RcppArmadilloExtensions/sample.h>

//#include <algorithm>


// An internal function for getting the K0 matrix (which may be given or may require calculating). This is where default kernel specifications are set.
arma::mat getK0(const arma::mat & samples, const arma::mat & derivatives, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<arma::vec> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & Z = R_NilValue);

// Finds the nearest positive definite matrix to handle numerical issues in matrix inverses.
arma::mat nearPD(arma::mat K0);	// Based off of Matlab code by John D'Errico

// An internal function for getting the matrix of square norms (which may be given or may require calculating).
arma::mat getSqNorm(const arma::mat & samples, const Rcpp::Nullable<Rcpp::IntegerVector> & nystrom_inds = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & Z = R_NilValue);

// An internal function used by CF_crossval_cpp. Given a single kernel and tuning parameter, this function uses (folds)-fold cross-validation to get the approximate mean square predictive error using the fitted gaussian process models from different estimation sets.
arma::vec CF_mse_cpp(arma::mat integrands, arma::mat samples, arma::mat derivatives, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, Rcpp::Nullable<unsigned int> folds = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds = R_NilValue, bool one_in_denom = false);

// An internal function used by SECF_crossval_cpp. Given a single kernel and tuning parameter, this function uses (folds)-fold cross-validation to get the approximate mean square predictive error using the fitted gaussian process models from different estimation sets.
arma::vec SECF_mse_cpp(arma::mat integrands, arma::mat samples, arma::mat derivatives, Rcpp::Function getX, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> sigma = R_NilValue, const Rcpp::Nullable<Rcpp::NumericMatrix> & K0 = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, Rcpp::Nullable<unsigned int> folds = R_NilValue, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds = R_NilValue);

// An internal function used by aSECF_crossval_cpp. Given a single kernel and tuning parameter, this function uses (folds)-fold cross-validation to get the approximate mean square predictive error using the fitted gaussian process models from different estimation sets.
arma::vec aSECF_mse_cpp(arma::mat integrands, arma::mat samples, arma::mat derivatives, Rcpp::Function getX, Rcpp::Function aSECF_mse_linsolve, Rcpp::Nullable<unsigned int> polyorder = R_NilValue, Rcpp::Nullable<unsigned int> steinOrder = R_NilValue, Rcpp::Nullable<Rcpp::String> kernel_function = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> sigma = R_NilValue, Rcpp::Nullable<Rcpp::IntegerVector> subset = R_NilValue, Rcpp::Nullable<unsigned int> folds = R_NilValue, bool conjugate_gradient = true, double reltol = 0.01, const Rcpp::Nullable<Rcpp::IntegerVector> & est_inds = R_NilValue);

#endif
