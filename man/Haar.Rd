\name{Haar}
\alias{dhaar}
\alias{Haar}
\alias{phaar}
\alias{rhaar}
\title{Haar Measure}
\usage{
  dhaar(r)

  phaar(q, lower.tail = TRUE)

  rhaar(n)
}
\arguments{
  \item{r,q}{vector of quantiles}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required}
}
\value{
  \code{dhaar} gives the density, \code{phaar} gives the
  distribution function, \code{rhaar} generates random
  deviates
}
\description{
  Uniform density on the circle
}
\details{
  The uniform density on the circle (also referred to as
  Haar measure) has the density
  \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}.}{C(r)=1-cos(r)/2\pi.}
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package
}

