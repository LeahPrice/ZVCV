# ZVCV
R package for derivative-based variance reduction

To download this package, use the following R code:
```{r}
library(devtools)

install_github("LeahPrice/ZVCV") 
library(ZVCV)

?ZVCV # use this to navigate through the documentation. ?zvcv has some example code
```


Details of the methods in this package can be found at:

* Mira, A., Solgi, R., & Imparato, D. (2013). *Zero variance Markov chain Monte Carlo for Bayesian estimators*. Statistics and Computing, 23(5), 653-662.
* Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.
* South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  https://arxiv.org/abs/2002.00033
* South, L. F., Oates, C. J., Mira, A., & Drovandi, C. (2018). *Regularised zero variance control variates*. https://arxiv.org/abs/1811.05073
