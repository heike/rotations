\name{rmedian}
\alias{rmedian}
\title{Compute the minimizer of the first order Euclidean distances.}
\usage{
  rmedian(Rs, epsilon = 1e-05, maxIter = 2000)
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
  a list \item{S}{the element in SO(3) minimizing the sum
  of first order Euclidean distances for sample Rs}
  \item{iter}{the number of iterations needed to converge
  or not}
}
\description{
  The embeded median type estimator we call the projected
  median and is given by
  \deqn{\widetilde{\bm{S}}_P=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_E(\bm{R}_i,\bm{S})}. The algorithm
  used is a modified Weiszfeld algorithm and is similar to
  the algorithm proposed by Hartley to compute the
  geometric median \eqn{\widetilde{\bm S}_G}.
}
\examples{
r<-rcayley(50,1)
 Rs<-genR(r)
 rmedian(Rs)
}
\seealso{
  \code{\link{MantonL2}}, \code{\link{HartleyL1}},
  \code{\link{arith.mean}}
}

