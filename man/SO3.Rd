\name{SO3}
\alias{SO3}
\title{Method for creating a rotation using the angle axis representation}
\usage{
  SO3(U, theta)
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

