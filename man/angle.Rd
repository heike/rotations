\name{angle}
\alias{angle}
\title{Find the angle of rotation R}
\usage{
  angle(Rs)
}
\arguments{
  \item{R}{rotation matrix in form of a 3-by-3 matrix in
  SO3}
}
\value{
  angle of rotation
}
\description{
  Extract angle from rotation matrix R. The extraction is
  based on the Rodrigues' Angle Axis representation, which
  leads us directy to the result that for any rotation
  matrix, the angle to the identity matrix is given as
}
\seealso{
  \code{\link{axis}}
}

