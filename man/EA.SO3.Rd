\name{EA.SO3}
\alias{EA.SO3}
\title{Translate from matrix to Euler angle format}
\usage{
  EA.SO3(rot)
}
\arguments{
  \item{rot}{a rotation matrix in the form of a 3x3 matrix}
}
\value{
  a vector of Euler angles
}
\description{
  Based on the Z-X-Z definition of Euler angles, this will
  take a 3x3 rotation matrix and return the Euler angle
  reprentation of that rotation.
}

