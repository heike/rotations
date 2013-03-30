\name{Fisher}
\alias{dfisher}
\alias{Fisher}
\alias{pfisher}
\alias{rfisher}
\title{The Matrix Fisher Distribution}
\usage{
  dfisher(r, kappa = 1, nu = NULL, Haar = TRUE)

  pfisher(q, kappa = 1, nu = NULL, lower.tail = TRUE)

  rfisher(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{r,q}{vector of quantiles}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required}

  \item{kappa}{concentration paramter}

  \item{nu}{circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}
}
\value{
  \code{dfisher} gives the density, \code{pfisher} gives
  the distribution function, \code{rfisher} generates
  random deviates
}
\description{
  Density and random generation for the matrix Fisher
  distribution with concentration kappa (\eqn{\kappa})
}
\details{
  The matrix Fisher distribution with concentration kappa
  (or circular variance nu) has density
  \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}{C(r|\kappa)=exp[2\kappa
  cos(r)][1-cos(r)]/(2\pi[I0(2\kappa)-I1(2\kappa)])} where
  \eqn{\mathrm{I_p}(\cdot)}{Ip()} denotes the Bessel
  function of order \eqn{p} defined as
  \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos
  r}dr}{Ip(\kappa)} is the modified Bessel function with
  parameters \eqn{p} and \eqn{kappa}.
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package
}

