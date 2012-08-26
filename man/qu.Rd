\name{qu}
\alias{qu}
\title{Reparameterize a rotation matrix as a unit quaternion}
\usage{
  qu(R)
}
\arguments{
  \item{R}{a rotation matrix in SO3}
}
\value{
  a unit quaternion of class Q4
}
\description{
  The rotation R is taken, the angle and axis of rotation
  are found and a unit quaternion is formed from the
  results.
}

