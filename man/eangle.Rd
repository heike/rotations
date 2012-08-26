\name{eangle}
\alias{eangle}
\title{Find the angle of rotation R}
\usage{
  eangle(Rs)
}
\arguments{
  \item{R}{3-by-3 matrix in SO3 whos angle of rotation is
  needed}
}
\value{
  angle of rotation
}
\description{
  This function will find the angle of rotation matrix R.
  Used mostly to reparameterize between quaternions and
  rotations.
}
\seealso{
  \code{\link{eaxis}}
}

