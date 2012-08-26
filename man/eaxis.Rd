\name{eaxis}
\alias{eaxis}
\title{Find the axis of rotation R}
\usage{
  eaxis(R)
}
\arguments{
  \item{R}{3-by-3 matrix in SO3 whos axis is needed}
}
\value{
  three dimensional axis of length one
}
\description{
  This function will find the axis of rotation matrix R.
  The simple calculation is based on Rodrigues formula and
  noticing that R - t(R) can be simplified greatly.  Used
  mostly to reparameterize between quaternions and
  rotations.
}
\seealso{
  \code{\link{eangle}}
}

