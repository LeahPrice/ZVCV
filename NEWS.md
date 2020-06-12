# ZVCV 0.1.0

* Initial release
* Includes 1st - 4th order zero variance control variates
* Regularisation is based on lasso or more generally elastic net where lambda is chosen using cross-validation and the elastic net parameter is specified

# ZVCV 0.1.1

* Allows for any order polynomial (with fast implementations available for polynomial orders Q = 1-4 or dimension d = 1)
* A cross validation method which chooses the polynomial order starting at 0 and going to infinity is now implemented (this new method could allow for super-root-N convergence)

# ZVCV 0.1.2

* Speeding up the higher order polynomial matrix getter using C++

# ZVCV 0.1.3

* Adding k-fold cross validation to select the polynomial order and doing some documentation improvements

# ZVCV 1.1.0

* Added control functionals, semi-exact control functionals and approximate semi-exact control functionals
