\docType{data}
\name{VDP}
\alias{VDP}
\title{Example of estimation using SMC}
\format{A list containing the following :
\describe{
 \item{N}{The size of the SMC population}
 \item{rho}{The tolerance for the new temperatures, which are selected so that the CESS at each temperature is \eqn{\rho*N} where \eqn{N} is the population size.}
 \item{temperatures}{A vector of length \eqn{T} of inverse power posterior temperatures}
 \item{samples}{An \eqn{N} by \eqn{d} by \eqn{T} matrix of samples from the \eqn{T} power posteriors, where \eqn{d} is the dimension of the target distribution. The samples are transformed to be on the log scale and all derivatives are with respect to log samples.}
 \item{loglike}{An \eqn{N} by \eqn{T} matrix of log likelihood values corresponding to \code{samples}}
 \item{logprior}{An \eqn{N} by \eqn{T} matrix of log prior values corresponding to \code{samples}}
 \item{der_loglike}{An \eqn{N} by \eqn{d} by \eqn{T} matrix of the derivatives of the log likelihood with respect to the parameters, with parameter values corresponding to \code{samples}}
 \item{der_logprior}{ An \eqn{N} by \eqn{d} by \eqn{T} matrix of the derivatives of the log prior with respect to the parameters, with parameter values corresponding to \code{samples}}
}}
\usage{
data(VDP)
}
\description{
This example illustrates how ZV-CV can be used for post-processing of results from SMC. In particular, we use ZV-CV to estimate posterior expectations and the evidence for a single SMC run of this example based on the Van der Pol oscillatory differential equations (Van der Pol, 1926). Further details about this example and applications to ZV-CV can be found in Oates et al. (2017) and South et al. (2019).

Given that the focus of this R package is on ZV-CV, we assume that samples have already been obtained from SMC and put into the correct format. One could use the R package \code{RcppSMC} or implement their own sampler in order to obtain results like this. The key is to make sure the derivatives of the log likelihood and log prior are stored, along with the inverse temperatures.

}
\examples{
set.seed(1)

# Load the SMC results
data(VDP) 

# Set up the list of control variates to choose from 
options <- list()
# Vanilla Monte Carlo
options[[1]] <- list(polyorder = 0)
# Standard ZV-CV with polynomial order selected through cross-validation
options[[2]] <- list(polyorder = Inf, regul_reg = FALSE)

##############################
# Posterior expectation - The true expectation is 0.9852 to 4 decimal places
##############################

# Note the exp() because samples and derivatives were stored on the log scale
# but we are interested in the expectation on the original scale
posterior <- zvcv(exp(VDP$samples[,,8]), VDP$samples[,,8],
VDP$der_loglike[,,8] + VDP$der_logprior[,,8], options = options) 
posterior$expectation # The posterior expectation estimate
posterior$polyorder # The selected polynomial order

##############################
# Evidence estimation - The true logged evidence is 10.36 to 2 decimal places
##############################

\dontrun{
# Getting additional temperatures based on maintaing a CESS of 0.99N rather than 0.9N
temp <- Expand_Temperatures(VDP$temperatures, VDP$loglike, 0.99)
VDP$temperatures_new <- temp$temperatures_all # the new temperatures
VDP$most_recent <- temp$relevant_samples # the samples associated with the new temperatures
}

\dontshow{
# Getting additional temperatures based on maintaing a CESS of 0.91N rather than 0.9N
temp <- Expand_Temperatures(VDP$temperatures, VDP$loglike, 0.91)
VDP$temperatures_new <- temp$temperatures_all # the new temperatures
VDP$most_recent <- temp$relevant_samples # the samples associated with the new temperatures
}

# Evidence estimation using the SMC identity
Z_SMC <- evidence_SMC(VDP$samples, VDP$loglike, VDP$der_loglike, VDP$der_logprior,
VDP$temperatures, VDP$temperatures_new, VDP$most_recent, options = options)
Z_SMC$log_evidence

# Evidence estimation using the CTI identity
Z_CTI <- evidence_CTI(VDP$samples, VDP$loglike, VDP$der_loglike, VDP$der_logprior,
VDP$temperatures, VDP$temperatures_new, VDP$most_recent, options = options)
Z_CTI$log_evidence_PS2
}
\references{ 
Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.

South, L. F., Oates, C. J., Mira, A., & Drovandi, C. (2019). Regularised zero-variance control variates for high-dimensional variance reduction.

Van der Pol, B. (1926). On relaxation-oscillations. The London, Edinburgh and Dublin Philosophical Magazine and Journal of Science, 2(11), 978-992.
}
\seealso{
See \link{ZVCV} for more package details.
}
