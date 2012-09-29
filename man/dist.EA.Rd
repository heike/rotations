\name{dist.EA}
\alias{dist.EA}
\title{Distance Between Two Rotations}
\usage{
  dist.EA(EA1, EA2 = id.EA, method = "projected", p = 1)
}
\arguments{
  \item{s}{The estimate of the central direction}

  \item{method}{calculate intrinsic or projected distance}

  \item{p}{the power of the respective distance}
}
\value{
  the pth power of the projected or intrinsic distance
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
EAs<-genR(r,space="EA")
EAp<-mean(EAs)
dist(EAp)
}

