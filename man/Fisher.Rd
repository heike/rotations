\name{Fisher}
\alias{Fisher}
\alias{dfisher}
\alias{rfisher}
\title{The Matrix Fisher Distribution}
\usage{
  dfisher(r, kappa = 1, nu = NULL, Haar = TRUE)

  rfisher(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{r}{vector of quantiles}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required}

  \item{kappa}{concentration paramter}

  \item{nu}{circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}
}
\value{
  \code{dfisher} gives the density, \code{rfisher}
  generates random deviates
}
\description{
  Density and random generation for the matrix Fisher
  distribution with concentration kappa
}
\details{
  The matrix Fisher distribution with concentration kappa
  (or circular variance nu) has density
  \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}
  where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel
  function of order \eqn{p} defined as
  \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos
  r}dr}.
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package
}

