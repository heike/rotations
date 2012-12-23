\name{Mises}
\alias{Mises}
\alias{dvmises}
\alias{rvmises}
\title{The circular-von Mises distribution}
\usage{
  dvmises(r, kappa = 1, nu = NULL, Haar = T)

  rvmises(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{r}{vector of quantiles}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required}

  \item{kappa}{concentration paramter}

  \item{nu}{The circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}
}
\value{
  \code{dvmises} gives the density, \code{rvmises}
  generates random deviates
}
\description{
  Density for the the circular von Mises-based distribution
  with concentration kappa
}
\details{
  The circular von Mises-based distribution has the density
  \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi
  \mathrm{I_0}(\kappa)}e^{\kappa cos(r)}.}
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package
}

