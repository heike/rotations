\name{SO3}
\alias{SO3}
\alias{SO3.default}
\alias{SO3.EA}
\alias{SO3.Q4}
\title{Method for creating a rotation in SO3 format using the angle axis representation}
\usage{
  SO3(U, ...)

  \method{SO3}{default} (U, theta)

  \method{SO3}{EA} (eur)

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

  \code{NULL}
}
\description{
  Angle-axis representation based on the Rodrigues formula.
}

