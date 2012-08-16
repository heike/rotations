\name{SO3.median}
\alias{SO3.median}
\title{Compute the projected or intrinsic median estimate of the central direction}
\usage{
  SO3.median(Rs, type = "projected", epsilon = 1e-05,
    maxIter = 2000)
}
\arguments{
  \item{Rs}{the sample \eqn{n \times 9} matrix with rows
  corresponding to observations}

  \item{type}{String indicating 'projeted' or 'intrinsic'
  type mean estimator}

  \item{epsilon}{the stopping rule for the iterative
  algorithm}

  \item{maxIter}{integer, the maximum number of iterations
  allowed}
}
\value{
  S the element in SO(3) minimizing the sum of first order
  Euclidean or Riemannian distances for sample Rs
}
\description{
  The median-type estimators are defined as
  \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}.  If the choice
  of distance metrid, \eqn{d_D}, is Riemannian then the
  estimator is called the intrinsic, and if the distance
  metric in Euclidean then it projected. The algorithm used
  in the intrinsic case is discussed in \cite{Hartley et
  al. (2011)} and the projected case was written by the
  authors.
}
\examples{
r<-rcayley(50,1)
Rs<-genR(r)
median(Rs)
}
\references{
  Hartley R, Aftab K and Trumpf J (2011). "L1 rotation
  averaging using the Weiszfeld algorithm." In _2011 IEEE
  Conference on Computer Vision and Pattern Recognition
  (CVPR)_, pp. 3041-3048. IEEE.
}
\seealso{
  \code{\link{SO3.mean}}
}

