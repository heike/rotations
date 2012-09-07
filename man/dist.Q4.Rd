\name{dist.Q4}
\alias{dist.Q4}
\title{Distance Between Two Rotations}
\usage{
  dist.Q4(Q1, Q2 = id.Q4, method = "projected", p = 1)
}
\arguments{
  \item{Q1}{First rotation in quaternion form}

  \item{Q2}{Second rotation in quaternion form, identity by
  default}

  \item{method}{calculate intrinsic or projected distance}

  \item{p}{the power of the respective distance}
}
\value{
  the pth power of the intrinsic or projected distance
  between Q1 and Q2
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
Qs<-genR(r,space="Q4")
Qp<-mean(Qs)
dist(Qp)
}

