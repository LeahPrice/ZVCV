# ZVCV 0.1.0

* Initial release
* Includes 1st - 4th order zero variance control variates
* Regularisation is based on lasso or more generally elastic net where lambda is chosen using cross-validation and the elastic net parameter is specified

# ZVCV 0.1.1

* Allows for any order polynomial (with fast implementations available for polynomial orders Q = 1-4 or dimension d = 1)
* A cross validation method which chooses the polynomial order starting at 1 and going to infinity is now implemented (this new method could allow for super-root-N convergence)

# ZVCV 0.1.2

* Speeding up the higher order polynomial matrix getter using C++

# ZVCV 0.1.3

* Adding k-fold cross validation to select the polynomial order and doing some documentation improvements

# ZVCV 1.1.0

* Added control functionals, semi-exact control functionals and approximate semi-exact control functionals

# ZVCV 2.1.0

* Making Linux friendly
* Adding more checks of input arguments
* Removing duplicates in kernel methods
* Returning the estimated coefficients in zvcv
* Changing some input arguments for zvcv:
- log_weight --> log_weights
- folds_choose --> folds
- obs_estim --> est_inds and is used to specify the estimation/fitting only samples (with the remainder being used for evaluation of the integrand
- REMOVED obs_estim_choose, the option to specify the samples for each cross-validation fold. I think this level of flexibility would rarely be required and it cause confusion when est_inds is specified.

# ZVCV 2.1.1

* Fixing a bug in the sample pre-processing for the special case of no duplicates + split estimation
* Other small bug fixes