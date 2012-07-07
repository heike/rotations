\name{rotation}
\alias{rotation}
\title{Class of rotation matrices.}
\description{
  An S4 class that stores a matrix describing a 3 by 3
  rotational matrix.
}
\examples{
a <- new("rotation") ## a is identity
b <- new("rotation", R=matrix(runif(9), nrow=3))
## b is not a rotation matrix
}

