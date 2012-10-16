\name{EA}
\alias{EA}
\title{Create Euler angles}
\usage{
  EA(U, theta)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{theta}{angle between -pi and pi}
}
\value{
  Euler angle of class "EA"
}
\description{
  Create Euler angles representing the rotation of the
  identity matrix about the axis U throught the angle theta
}

