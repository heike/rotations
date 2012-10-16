\name{axis.EA}
\alias{axis.EA}
\title{Find the axis of rotation R}
\usage{
  axis.EA(eur)
}
\arguments{
  \item{eur}{Euler angle}
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

