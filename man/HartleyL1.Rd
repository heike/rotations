\name{HartleyL1}
\alias{HartleyL1}
\title{Compute the geometric median of a sample of random rotations}
\usage{
  HartleyL1(Rs, epsilon = 1e-05, maxIter = 1000)
}
\arguments{
  \item{Rs}{the sample \eqn{n \times 9} matrix with rows
  corresponding to observations}

  \item{epsilon}{the stopping rule for the iterative
  algorithm}

  \item{maxIter}{integer, the maximum number of iterations
  allowed}
}
\value{
  list \item{S}{the element in SO(3) minimizing the sum of
  first order Riemannian distances for sample Rs}
  \item{iter}{the number of iterations needed to converge
  or not}
}
\description{
  This function uses the algorithm published by Hartley to
  estimate the principle direction of a sample of random
  rotaions with the point in \eqn{SO(3)} that minimizes the
  sum of first order Riemannian distances, aka the
  geometric median and denoted \eqn{\widetilde{\bm S}_G}.
  More explicitly \deqn{\widetilde{\bm
  S}_G=\widetilde{\bm{S}}_{G}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_G(\bm{R}_i,\bm{S})}.
}
\examples{
r<-rvmises(20,0.1)
Rs<-genR(r)
HartleyL1(Rs)
}
\seealso{
  \code{\link{MantonL2}}, \code{\link{arith.mean}},
  \code{\link{rmedian}}
}

