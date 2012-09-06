\name{median.EA}
\alias{median.EA}
\title{Compute the projected or intrinsic median estimate of the central direction}
\usage{
  median.EA(EAs, type = "projected", epsilon = 1e-05,
    maxIter = 2000)
}
\arguments{
  \item{EAs}{A \eqn{n\times 3} matrix where each row
  corresponds to a random rotation in Euler angle form}

  \item{type}{String indicating 'projeted' or 'intrinsic'
  type mean estimator}

  \item{epsilon}{the stopping rule for the iterative
  algorithm}

  \item{maxIter}{integer, the maximum number of iterations
  allowed}
}
\value{
  S the vector of Euler angles minimizing the sum of first
  order Euclidean or Riemannian distances for sample Rs
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
EA<-genR(r,space="R3")
median.R3(EA)
}
\references{
  Hartley R, Aftab K and Trumpf J (2011). "L1 rotation
  averaging using the Weiszfeld algorithm." In _2011 IEEE
  Conference on Computer Vision and Pattern Recognition
  (CVPR)_, pp. 3041-3048. IEEE.
}
\seealso{
  \code{\link{median.SO3}}
}

