\name{riedist.SO3}
\alias{riedist.SO3}
\title{Riemannian Distance Between Two Random Rotations in matrix format}
\usage{
  riedist.SO3(R, S = diag(1, 3, 3))
}
\arguments{
  \item{R}{The estimate of the central direction}

  \item{S}{The true central direction}
}
\value{
  S3 \code{riedist} object; a number between 0 and pi that
  is the shortest geodesic curve connecting two matrices,
  i.e., the Riemannian distance
}
\description{
  This function will calculate the riemannian distance
  between an estimate of the central direction (in matrix
  or vector form) and the central direction.  By default
  the central direction is taken to be the identity matrix,
  but any matrix in SO(3) will work.  This distance is
  equivalent to the angle of rotation of R'S
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
Sp<-mean(Rs)
riedist.SO3(Sp,diag(1,3,3))
}

