\name{Haar}
\alias{dhaar}
\alias{Haar}
\alias{rhaar}
\title{Haar Measure}
\usage{
  dhaar(r, lower.tail = TRUE)

  rhaar(n)
}
\arguments{
  \item{r}{Where the density is being evaluated}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required}

  \item{lower.tail}{logical; if TRUE probabilites are
  \eqn{P(X\le x)}}
}
\value{
  \code{dhaar} gives the density, \code{rhaar} generates
  random deviates
}
\description{
  Uniform density on the circle
}
\details{
  The uniform density on the circle (also referred to as
  Haar measure) has the density
  \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}.}
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package
}

