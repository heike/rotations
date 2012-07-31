\name{dvmises}
\alias{dvmises}
\title{Density function for circular von Mises distribution}
\usage{
  dvmises(r, kappa = 1, Haar = F)
}
\arguments{
  \item{r}{value at which to evaluate the distribution
  function}

  \item{kappa}{concentration paramter}

  \item{Haar}{logical, if density is evaluated with respect
  to Haar measure or Lebesgue}
}
\value{
  value of circular-von Mises distribution with
  concentration \eqn{\kappa} evaluated at r
}
\description{
  The circular von Mises-based distribution has the density
  \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi
  \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.  This function
  allows the use to evaluate \eqn{C_\mathrm{M}(r|\kappa)}
  at angle \eqn{r} given a concentration parameter
  \eqn{\kappa}.
}
\seealso{
  \code{\link{rvmises}},
  \code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dcayley}}
}

