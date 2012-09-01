\name{dist}
\alias{dist}
\alias{dist.Q4}
\alias{dist.SO3}
\title{Riemannian Distance Between Two Rotations}
\usage{
  dist(Rs, method = "euclidean", p = 1)

  \method{dist}{SO3} (Rs, method = "euclidean", p = 1)

  \method{dist}{Q4} (Qs, method = "euclidean", p = 1)
}
\arguments{
  \item{Rs}{The estimate of the central direction}

  \item{S}{The true central direction}
}
\value{
  S3 \code{dist} object; a number between 0 and pi that is
  the shortest geodesic curve connecting two matrices,
  i.e., the Riemannian distance

  \code{NULL}

  \code{NULL}
}
\description{
  This function will calculate the riemannian or Euclidean
  distance between two rotations.  If only one rotation is
  specified the other will be set to the identity and the
  distance between the two is returned.
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
Sp<-mean(Rs)
dist.SO3(Sp)
}

