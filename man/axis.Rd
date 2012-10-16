\name{axis}
\alias{axis}
\alias{axis.EA}
\alias{axis.Q4}
\alias{axis.SO3}
\title{Find the axis of rotation R}
\usage{
  axis(R)

  \method{axis}{SO3} (R)

  \method{axis}{Q4} (q)

  \method{axis}{EA} (eur)
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
  This function will find the axis of rotation matrix R.
  The simple calculation is based on Rodrigues formula and
  noticing that R - t(R) can be simplified greatly.
}
\seealso{
  \code{\link{angle}}
}

