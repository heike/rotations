\name{SumDist}
\alias{SumDist}
\title{Compute the sum of the \eqn{p^{\text{th}}} order distances between Rs and S}
\usage{
  SumDist(Rs, S = diag(1, 3, 3), p)
}
\arguments{
  \item{Rs}{numeric matrix with sample size n rows and m
  columns}

  \item{S}{the matrix to compute the sum of distances
  between each row of Rs with}

  \item{p}{the order of the distances to compute}
}
\value{
  list of size two \item{Rieman}{the sum of
  \eqn{p^{\text{th}}} order Riemannian distances}
  \item{Euclid}{the sum of \eqn{p^{\text{th}}} order
  Euclidean distances}
}
\description{
  Compute the sum of the \eqn{p^{\text{th}}} order
  distances between Rs and S
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
Sp<-arith_mean(Rs)
SumDist(Rs,S=Sp,2)
}

