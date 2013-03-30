\name{angle}
\alias{angle}
\alias{angle.Q4}
\alias{angle.SO3}
\title{Misorientation Angle}
\usage{
  angle(Rs)

  \method{angle}{SO3} (Rs)

  \method{angle}{Q4} (Qs)
}
\arguments{
  \item{Rs}{rotation matrix}
}
\value{
  angle of rotation
}
\description{
  Find the misorientation angle of a rotation
}
\details{
  Every rotation can be thought of as some reference
  coordinate system rotated about an axis through an angle.
  These quantites are referred to as the misorientation
  axis and misorientation angle, respectively, in the
  material sciences literature. This function returns the
  misorentation angle associated with a rotation assuming
  the reference coordinate system is the identity.
}
\seealso{
  \code{\link{axis}}
}

