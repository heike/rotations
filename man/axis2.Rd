\name{axis2}
\alias{axis2}
\alias{axis2.EA}
\alias{axis2.Q4}
\alias{axis2.SO3}
\title{Find the axis of rotation R}
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
  This function will find the axis of rotation matrix R.
  The simple calculation is based on Rodrigues formula and
  noticing that R - t(R) can be simplified greatly.
}
\seealso{
  \code{\link{angle}}
}

