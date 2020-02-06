#ifndef ZVCV_GETPOLY_H
#define ZVCV_GETPOLY_H 1.0

#include <RcppArmadillo.h>
#include <math.h>
#include <cmath>

#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <cmath>

using namespace std;


// Gets all potential combinations of the values in mymat
arma::umat get_all_combins(arma::umat mymat, unsigned int polyorder);

// Turns the matrix of all combinations of polynomial orders into the matrix of covariates using the samples and derivatives. Does this WITHOUT using the R package partitions which cannot be used on Linux.
arma::mat getPoly_withoutpackage(arma::mat samples, arma::mat derivatives, arma::mat combinations);

// Turns the matrix of all combinations of polynomial orders into the matrix of covariates using the samples and derivatives. Does this WITH the R package partitions which cannot be used on Linux.
arma::mat getPoly_withpackage(arma::mat samples, arma::mat derivatives, Rcpp::List combinations);

// The number of ways you can choose k values from a list of n without taking order into account (i.e. the number of combinations)
unsigned int choose(unsigned int n, unsigned int k);

// Gets the matrix of combinations of m values from the vector x (which is length > m)
// Converted from base R's combn function
arma::umat combn(arma::uvec x, unsigned int m);

// Gets the vector of elements that are in x but not y
arma::uvec setdiff(arma::uvec x, arma::uvec y);

#endif
