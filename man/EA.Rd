\name{EA}
\alias{EA}
\alias{EA.default}
\alias{EA.Q4}
\alias{EA.SO3}
\title{Create Euler angles}
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
  Create Euler angles representing the rotation of the
  identity matrix about the axis U throught the angle theta
}

