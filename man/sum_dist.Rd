\name{sum_dist}
\alias{sum_dist}
\alias{sum_dist.Q4}
\alias{sum_dist.SO3}
\title{Sample Distance}
\usage{
  sum_dist(Rs, S = genR(0, space = class(Rs)),
    method = "projected", p = 1)

  \method{sum_dist}{SO3} (Rs, S = id.SO3,
    method = "projected", p = 1)

  \method{sum_dist}{Q4} (Qs, S = id.Q4,
    method = "projected", p = 1)
}
\arguments{
  \item{Rs}{a matrix of rotation observations, one row per
  observation}

  \item{S}{the individual matrix of interest, usually an
  estimate of the mean}

  \item{method}{type of distance used method in 'projected'
  or 'intrinsic'}

  \item{p}{the order of the distances to compute}
}
\value{
  the sum of the pth order distance between each sample in
  Rs and S

  \code{NULL}

  \code{NULL}
}
\description{
  Compute the sum of the \eqn{p^{\text{th}}}{pth} order
  distances between Rs and S
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
Sp<-mean(Rs)
sum_dist(Rs,S=Sp,p=2)
}

