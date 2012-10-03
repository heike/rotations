\name{median.SO3}
\alias{median.SO3}
\title{Compute the projected or intrinsic median estimate of the central direction}
\usage{
  median.SO3(x, type = "projected", epsilon = 1e-05,
    maxIter = 2000, na.rm = FALSE)
}
\arguments{
  \item{x}{A \eqn{n\times 9} matrix where each row
  corresponds to a random rotation in matrix form}

  \item{type}{String indicating 'projected' or 'intrinsic'
  type mean estimator}

  \item{epsilon}{the stopping rule for the iterative
  algorithm}

  \item{maxIter}{integer, the maximum number of iterations
  allowed}

  \item{na.rm}{not used, only there for inconsistency with
  median.default}
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
  \code{\link{mean.SO3}}
}

