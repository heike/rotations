\name{dist.SO3}
\alias{dist.SO3}
\title{Distance Between Two Rotations}
\usage{
  dist.SO3(R, S = id.SO3, method = "euclidean", p = 1)
}
\arguments{
  \item{R}{rotation in SO3 representation}

  \item{S}{rotation in SO3 representation}

  \item{method}{calculate riemannian or euclidean distance}

  \item{p}{the power of the respective distance}
}
\value{
  the pth power of the euclidean or riemannian distance
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

