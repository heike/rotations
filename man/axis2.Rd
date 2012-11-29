\name{axis2}
\alias{axis2}
\alias{axis2.EA}
\alias{axis2.Q4}
\alias{axis2.SO3}
\title{Misorientation Axis}
\usage{
  axis2(R)

  \method{axis2}{SO3} (R)

  \method{axis2}{Q4} (q)

  \method{axis2}{EA} (eur)
}
\arguments{
  \item{R}{3-by-3 matrix in SO3}
}
\value{
  axis in form of three dimensional vector of length one.

  \code{NULL}

  \code{NULL}

  \code{NULL}
}
\description{
  Find the misorientation axis of a rotation
}
\details{
  Every rotation can be thought of as some reference
  coordinate system rotated about an axis through an angle.
  These quantites are referred to as the misorientation
  axis and misorientation angle, respectively, in the
  material sciences literature. This function returns the
  misorentation axis associated with a rotation assuming
  the reference coordinate system is the identity.
}
\seealso{
  \code{\link{angle}}
}

