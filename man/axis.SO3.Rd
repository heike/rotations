\name{axis.SO3}
\alias{axis.SO3}
\title{Find the axis of rotation R}
\usage{
  axis.SO3(R)
}
\arguments{
  \item{R}{3-by-3 matrix in SO3}
}
\value{
  axis in form of three dimensional vector of length one.
}
\description{
  This function will find the axis of rotation matrix R.
  The simple calculation is based on Rodrigues formula and
  noticing that R - t(R) can be simplified greatly.
}
\seealso{
  \code{\link{angle}}
}

