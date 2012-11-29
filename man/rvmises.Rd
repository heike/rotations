\name{rvmises}
\alias{rvmises}
\title{The circular-von Mises distribution}
\usage{
  rvmises(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{n}{The number of angles desired}

  \item{kappa}{The concentration parameter of the
  distribution}

  \item{nu}{An alternative to kappa; circular variance}
}
\value{
  S3 \code{rvmises} object; a vector of n angles following
  the von Mises Circular distribution with concentration
  kappa and mean/mode 0
}
\description{
  Random generation for the the circular von Mises-based
  distribution with concentration kappa
}
\details{
  The circular von Mises-based distribution has the density
  \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi
  \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.
}
\examples{
r<-rvmises(20,0.01)
}

