\name{Q4}
\alias{Q4}
\alias{Q4.default}
\alias{Q4.EA}
\alias{Q4.SO3}
\title{Quaternions}
\usage{
  Q4(U, ...)

  \method{Q4}{default} (U, theta)

  \method{Q4}{SO3} (R)

  \method{Q4}{EA} (eur)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{...}{additional arguments}
}
\value{
  unit quaternion of class "Q4"

  \code{NULL}

  \code{NULL}

  \code{NULL}
}
\description{
  Create a unit quaternion
}
\details{
  Create quaternion representing the rotation of the
  identity matrix about the axis U throught the angle
  theta.  This can be accomplished by providing the axis
  and angle explicitly or by providing the rotation in some
  other form, e.g. a matrix in SO(3) or Euler angles.
}
\seealso{
  Other Q4: \code{\link{as.Q4}}, \code{\link{id.Q4}},
  \code{\link{Q4-class}}

  Other Q4: \code{\link{as.Q4}}, \code{\link{id.Q4}},
  \code{\link{Q4-class}}

  Other Q4: \code{\link{as.Q4}}, \code{\link{id.Q4}},
  \code{\link{Q4-class}}

  Other Q4: \code{\link{as.Q4}}, \code{\link{id.Q4}},
  \code{\link{Q4-class}}
}

