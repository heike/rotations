\name{riedist}
\alias{riedist}
\title{Riemannian Distance Between Two Random Rotations}
\usage{
  riedist(R, S = diag(1, 3, 3))
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
  but any matrix in SO(3) will work.  It calls the matrix
  log and matrix exponential functions also given here.
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
Sp<-arith.mean(Rs)
riedist(Sp,diag(1,3,3))
}

