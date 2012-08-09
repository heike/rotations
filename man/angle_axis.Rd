\name{angle_axis}
\alias{angle_axis}
\title{Method for creating rotations using the angle axis representation}
\usage{
  angle_axis(U, theta)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{theta}{angle between -pi and pi}
}
\value{
  Used in \code{\link{eyeBall}} to orient the data properly
}
\description{
  Angle-axis representation based on the Rodrigues formula.
}

