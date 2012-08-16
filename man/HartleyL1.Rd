\name{HartleyL1}
\alias{HartleyL1}
\title{Compute the geometric median of a sample of random rotations.  To be deleted}
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
  THIS HAS BEEN REPLACED WITH 'median.SO3'!!  This function
  uses the algorithm developed in \cite{Hartley et al.
  (2011)} to estimate the principle direction of a sample
  of random rotaions with the point in \eqn{SO(3)} that
  minimizes the sum of first order Riemannian distances,
  aka the geometric median and denoted \eqn{\widetilde{\bm
  S}_G}.  More explicitly \deqn{\widetilde{\bm
  S}_G=\widetilde{\bm{S}}_{G}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_G(\bm{R}_i,\bm{S})}.
}
\examples{
r<-rvmises(20,0.1)
Rs<-genR(r)
HartleyL1(Rs)
}
\references{
  Hartley R, Aftab K and Trumpf J (2011). "L1 rotation
  averaging using the Weiszfeld algorithm." In _2011 IEEE
  Conference on Computer Vision and Pattern Recognition
  (CVPR)_, pp. 3041-3048. IEEE.
}
\seealso{
  \code{\link{MantonL2}}, \code{\link{arith_mean}},
  \code{\link{rmedian}}
}

