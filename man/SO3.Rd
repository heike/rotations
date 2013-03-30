\name{SO3}
\alias{SO3}
\alias{SO3.default}
\alias{SO3.Q4}
\title{Matrix in SO(3)}
\usage{
  SO3(U, ...)

  \method{SO3}{default} (U, theta = NULL)

  \method{SO3}{Q4} (q)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{...}{additional arguments}
}
\value{
  matrix of rotations in SO3 format

  \code{NULL}

  \code{NULL}
}
\description{
  Create a rotation matrix
}
\details{
  Create matrix in SO(3) representing the rotation of the
  identity matrix about the axis U throught the angle
  theta.  This can be accomplished by providing the axis
  and angle explicitly or by providing the rotation in some
  other form, e.g. a vector of Euler angles or unit
  quaternion.
}
\seealso{
  Other SO3: \code{\link{as.SO3}}, \code{\link{id.SO3}},
  \code{\link{is.Q4}}, \code{\link{is.SO3}},
  \code{\link{SO3-class}}

  Other SO3: \code{\link{as.SO3}}, \code{\link{id.SO3}},
  \code{\link{is.Q4}}, \code{\link{is.SO3}},
  \code{\link{SO3-class}}

  Other SO3: \code{\link{as.SO3}}, \code{\link{id.SO3}},
  \code{\link{is.Q4}}, \code{\link{is.SO3}},
  \code{\link{SO3-class}}
}

