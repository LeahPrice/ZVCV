# ZVCV 0.1.0

* Initial release
* Includes 1st - 4th order zero variance control variates
* Regularisation is based on lasso or more generally elastic net where lambda is chosen using cross-validation and the elastic net parameter is specified

# ZVCV 0.1.1

* Allows for any order polynomial (with fast implementations available for polynomial orders Q = 1-4 or dimension d = 1)
* A cross validation method which chooses the polynomial order starting at 0 and going to infinity is now implemented (this new method could allow for super-root-N convergence)
