\name{axis.Q4}
\alias{axis.Q4}
\title{Find the axis of rotation R}
\usage{
  axis.Q4(q)
}
\arguments{
  \item{q}{unit quaternion}
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

