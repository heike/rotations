\name{EA}
\alias{EA}
\alias{EA.default}
\alias{EA.Q4}
\alias{EA.SO3}
\title{Euler Angles}
\usage{
  EA(U, ...)

  \method{EA}{default} (U, theta)

  \method{EA}{SO3} (R)

  \method{EA}{Q4} (Qs)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{...}{additional arguments}
}
\value{
  Euler angle of class "EA"

  \code{NULL}

  \code{NULL}

  \code{NULL}
}
\description{
  Create a Euler angles
}
\details{
  Create Euler angles representing the rotation of the
  identity matrix about the axis U throught the angle
  theta.  This can be accomplished by providing the axis
  and angle explicitly or by providing the rotation in some
  other form, e.g. a matrix in SO(3) or unit quaternion.
}
\seealso{
  Other EA: \code{\link{as.EA}}, \code{\link{EA-class}},
  \code{\link{id.EA}}

  Other EA: \code{\link{as.EA}}, \code{\link{EA-class}},
  \code{\link{id.EA}}

  Other EA: \code{\link{as.EA}}, \code{\link{EA-class}},
  \code{\link{id.EA}}
}

