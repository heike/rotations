\name{Q4}
\alias{Q4}
\title{Form a unit quaterion}
\usage{
  Q4(U, theta)
}
\arguments{
  \item{U}{three-dimensional vector describing the fix axis
  of the rotation}

  \item{theta}{angle between -pi and pi}
}
\value{
  unit quaternion of class "Q4"
}
\description{
  Create a unit quaternion representing the rotation of the
  identity matrix about the axis U throught the angle theta
}

