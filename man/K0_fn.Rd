% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{K0_fn}
\alias{K0_fn}
\title{Kernel matrix calculation}
\usage{
K0_fn(
  samples,
  derivatives,
  sigma,
  steinOrder,
  kernel_function,
  Z = NULL,
  nystrom_inds = NULL
)
}
\arguments{
\item{samples}{An \eqn{N} by \eqn{d} matrix of samples from the target}

\item{derivatives}{An \eqn{N} by \eqn{d} matrix of derivatives of the log target with respect to the parameters}

\item{sigma}{The tuning parameters of the specified kernel. This involves a single length-scale parameter in "gaussian" and "RQ", a length-scale and a smoothness parameter in "matern" and two parameters in "product" and "prodsim". See below for further details.}

\item{steinOrder}{This is the order of the Stein operator. The default is \code{1} in the control functionals paper (Oates et al, 2017) and \code{2} in the semi-exact control functionals paper (South et al, 2020).  The following values are currently available: \code{1} for all kernels and \code{2} for "gaussian", "matern" and "RQ". See below for further details.}

\item{kernel_function}{Choose between "gaussian", "matern", "RQ", "product" or "prodsim". See below for further details.}

\item{Z}{(optional) An \eqn{N} by \eqn{N} (or \eqn{N} by \eqn{m} where \eqn{m} is the length of \code{nystrom_inds}). This can be calculated using \code{\link{squareNorm}}.}

\item{nystrom_inds}{(optional) The sample indices to be used in the Nystrom approximation (for when using aSECF).}
}
\value{
An \eqn{N} by \eqn{N} kernel matrix (or \eqn{N} by \eqn{m} where \eqn{m} is the length of \code{nystrom_inds}).
}
\description{
This function calculates the full \eqn{K_0} matrix, which is a first or second order Stein operator applied to
a standard kernel. 
The output of this function can be used as an argument to \code{\link{CF}}, \code{\link{CF_crossval}},
\code{\link{SECF}}, \code{\link{SECF_crossval}}, \code{\link{aSECF}} and \code{\link{aSECF_crossval}}.
The kernel matrix is automatically computed in all of the above methods, but it is faster to calculate
in advance when using more than one of the above functions and when using any of the crossval functions.
}
\section{On the choice of \eqn{\sigma}, the kernel and the Stein order}{

The kernel in Stein-based kernel methods is \eqn{L_x L_y k(x,y)} where \eqn{L_x} is a first or second order Stein operator in \eqn{x} and \eqn{k(x,y)} is some generic kernel to be specified.

The Stein operators for distribution \eqn{p(x)} are defined as:
\itemize{
\item \strong{\code{steinOrder=1}}: \eqn{L_x g(x) = \nabla_x^T g(x) + \nabla_x \log p(x)^T g(x)} (see e.g. Oates el al (2017))
\item \strong{\code{steinOrder=2}}: \eqn{L_x g(x) = \Delta_x g(x) + \nabla_x log p(x)^T \nabla_x g(x)} (see e.g. South el al (2020))
}
Here \eqn{\nabla_x} is the first order derivative wrt \eqn{x} and \eqn{\Delta_x = \nabla_x^T \nabla_x} is the Laplacian operator.

The generic kernels which are implemented in this package are listed below.  Note that the input parameter \strong{\code{sigma}} defines the kernel parameters \eqn{\sigma}. 
\itemize{
\item \strong{\code{"gaussian"}}: A Gaussian kernel,
\deqn{k(x,y) = exp(-z(x,y)/\sigma^2)}
\item \strong{{\code{"matern"}}}: A Matern kernel with \eqn{\sigma = (\lambda,\nu)},
\deqn{k(x,y) = bc^{\nu}z(x,y)^{\nu/2}K_{\nu}(c z(x,y)^{0.5})} where \eqn{b=2^{1-\nu}(\Gamma(\nu))^{-1}}, \eqn{c=(2\nu)^{0.5}\lambda^{-1}} and \eqn{K_{\nu}(x)} is the modified Bessel function of the second kind. Note that \eqn{\lambda} is the length-scale parameter and \eqn{\nu} is the smoothness parameter (which defaults to 2.5 for \eqn{steinOrder=1} and 4.5 for \eqn{steinOrder=2}).
\item \strong{\code{"RQ"}}: A rational quadratic kernel,
\deqn{k(x,y) = (1+\sigma^{-2}z(x,y))^{-1}}
\item \strong{\code{"product"}}: The product kernel that appears in Oates et al (2017) with \eqn{\sigma = (a,b)}
\deqn{k(x,y) = (1+a z(x) + a z(y))^{-1} exp(-0.5 b^{-2} z(x,y)) }
\item \strong{\code{"prodsim"}}: A slightly different product kernel with \eqn{\sigma = (a,b)} (see e.g. \url{https://www.imperial.ac.uk/inference-group/projects/monte-carlo-methods/control-functionals/}),
\deqn{k(x,y) = (1+a z(x))^{-1}(1 + a z(y))^{-1} exp(-0.5 b^{-2} z(x,y)) }
}
In the above equations, \eqn{z(x) = \sum_j x[j]^2} and \eqn{z(x,y) = \sum_j (x[j] - y[j])^2}. For the last two kernels, the code only has implementations for \code{steinOrder}=\code{1}. Each combination of \code{steinOrder} and \code{kernel_function} above is currently hard-coded but it may be possible to extend this to other kernels in future versions using autodiff. The calculations for the first three kernels above are detailed in South et al (2020).
}

\references{
Oates, C. J., Girolami, M. & Chopin, N. (2017). Control functionals for Monte Carlo integration. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79(3), 695-718.

South, L. F., Karvonen, T., Nemeth, C., Girolami, M. and Oates, C. J. (2020). Semi-Exact Control Functionals From Sard's Method.  \url{https://arxiv.org/abs/2002.00033}
}
\author{
Leah F. South
}
