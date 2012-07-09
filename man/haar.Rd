\name{haar}
\alias{cayley}
\alias{fisher}
\alias{haar}
\alias{mises}
\title{Angle distributions}
\usage{
  haar <- function(x) fisher <- function(x, kappa) cayley
    <- function(x, kappa) mises <- function(x, kappa)
}
\arguments{
  \item{x}{vector of angles (between -pi and pi)}

  \item{kappa}{concentration parameter, kappa > 0}
}
\description{
  Densities of the most commonly used angle distributions:
  Haar, Fisher, Cayley, and von Mises.
}
\details{
  The Haar density for angles is given as
  \eqn{(1-\cos(x))/(2*\pi)}, where \eqn{x \in (-\pi, \pi)}.
  When the fix axis \eqn{U} of rotation \eqn{R} is randomly
  picked, rotations with angles from the Haar measure are
  uniformly distributed on the space of all 3d rotations,
  \eqn{SO(3)}.
}
\examples{
x <- seq(-pi, pi, by=0.01)
plot(x, haar(x), type="l")
lines(x, fisher(x, kappa=1), col=2)
lines(x, cayley(x, kappa=1), col=3)

plot(x, fisher(x, kappa=1)/haar(x), type="l")
lines(x, cayley(x, kappa=1)/haar(x), col=3)
# Cayley with the same concentration parameter is heavier tailed than Fisher.
}

