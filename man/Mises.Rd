\name{Mises}
\alias{dvmises}
\alias{Mises}
\alias{rvmises}
\title{The circular-von Mises distribution}
\usage{
  dvmises(r, kappa = 1, nu = NULL, Haar = T,
    lower.tail = TRUE)

  rvmises(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{r}{vector of quantiles}

  \item{n}{number of observations}

  \item{kappa}{concentration paramter}

  \item{nu}{The circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}

  \item{lower.tail}{logical; if TRUE probabilites are
  \eqn{P(X\le x)}}
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

