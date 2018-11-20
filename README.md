# ZVCV
R package for Regularised Zero-Variance Control Variates

To download this package, use the following R code:
```{r}
library(devtools)

install_github("LeahPrice/ZVCV") 
library(ZVCV)

?ZVCV # use this to navigate through the documentation. ?zvcv has some example code
```

I plan to add more details about the ZV-CV method to the documentation. I also need to add more details about the evidence estimators and give an example of evidence estimation with the package.

Details of the ZV-CV method can be found in:

Mira, A., Solgi, R., & Imparato, D. (2013). *Zero variance Markov chain Monte Carlo for Bayesian estimators*. Statistics and Computing, 23(5), 653-662.

The regularised ZV-CV method is described in the following paper which will appear on arXiv within the next few days:

South, L. F., Mira, A., & Drovandi, C. (2018). *Regularised zero variance control variates*.
