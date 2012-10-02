\name{SO3}
\alias{SO3}
\title{Method for creating a rotation in SO3 format using the angle axis representation}
\usage{
  SO3(U, theta)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{theta}{angle between -pi and pi}
}
\value{
  matrix of rotations in SO3 format
}
\description{
  Angle-axis representation based on the Rodrigues formula.
}

