\name{dist.SO3}
\alias{dist.SO3}
\title{Distance Between Two Rotations}
\usage{
  \method{dist}{SO3} (R1, R2 = id.SO3,
    method = "projected", p = 1)
}
\arguments{
  \item{R1}{rotation in SO3 representation}

  \item{R2}{rotation in SO3 representation}

  \item{method}{calculate intrinsic or projected distance}

  \item{p}{the power of the respective distance}
}
\value{
  the pth power of the intrinsic or projected distance
}
\description{
  This function will calculate the intrinsic (Riemannian)
  or projected (Euclidean) distance between two rotations.
  If only one rotation is specified the other will be set
  to the identity and the distance between the two is
  returned.
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
Sp<-mean(Rs)
dist.SO3(Sp)
}

