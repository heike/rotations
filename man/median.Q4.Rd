\name{median.Q4}
\alias{median.Q4}
\title{Compute the projected or intrinsic median estimate of the central direction}
\usage{
  median.Q4(Qs, type = "projected", epsilon = 1e-05,
    maxIter = 2000)
}
\arguments{
  \item{Qs}{A \eqn{n\times 4} matrix where each row
  corresponds to a random rotation in unit quaternion form}

  \item{type}{String indicating 'projeted' or 'intrinsic'
  type mean estimator}

  \item{epsilon}{the stopping rule for the iterative
  algorithm}

  \item{maxIter}{integer, the maximum number of iterations
  allowed}
}
\value{
  S the unit quaternion minimizing the sum of first order
  Euclidean or Riemannian distances for sample Rs
}
\description{
  The median-type estimators are defined as
  \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}.  If the choice
  of distance metrid, \eqn{d_D}, is Riemannian then the
  estimator is called the intrinsic, and if the distance
  metric in Euclidean then it projected. The algorithm used
  in the intrinsic case is discussed in \cite{hartley11}
  and the projected case was written by the authors.
}
\examples{
r<-rcayley(50,1)
Qs<-genR(r,space="Q4")
median(Qs)
}
\seealso{
  \code{\link{median.SO3}}
}

