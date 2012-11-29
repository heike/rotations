\name{dhaar}
\alias{dhaar}
\title{Haar Measure}
\usage{
  dhaar(r, lower.tail = TRUE)
}
\arguments{
  \item{r}{Where the density is being evaluated}

  \item{lower.tail}{logical; if TRUE probabilites are
  \eqn{P(X\le x)}}
}
\value{
  the probability density evaluated at r
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
  \code{\link{rhaar}},
  \code{\link{dfisher}},\code{\link{dvmises}},\code{\link{dcayley}}
}

