\name{dvmises}
\alias{dvmises}
\title{The circular-von Mises distribution}
\usage{
  dvmises(r, kappa = 1, nu = NULL, Haar = T,
    lower.tail = TRUE)
}
\arguments{
  \item{r}{vector of quantiles}

  \item{kappa}{concentration paramter}

  \item{nu}{The circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}

  \item{lower.tail}{logical; if TRUE probabilites are
  \eqn{P(X\leq x)}}
}
\value{
  value of circular-von Mises distribution with
  concentration \eqn{\kappa} evaluated at r
}
\description{
  Density for the the circular von Mises-based distribution
  with concentration kappa
}
\details{
  The circular von Mises-based distribution has the density
  \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi
  \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.
}
\seealso{
  \code{\link{rvmises}},
  \code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dcayley}}
}

