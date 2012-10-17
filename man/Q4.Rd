\name{Q4}
\alias{Q4}
\alias{Q4.default}
\alias{Q4.SO3}
\title{Form a unit quaterion}
\usage{
  Q4(U, ...)

  \method{Q4}{default} (U, theta)

  \method{Q4}{SO3} (eur)

  \method{Q4}{EA} (eur)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{theta}{angle between -pi and pi}
}
\value{
  unit quaternion of class "Q4"

  \code{NULL}

  \code{NULL}

  \code{NULL}
}
\description{
  Create a unit quaternion representing the rotation of the
  identity matrix about the axis U throught the angle theta
}

