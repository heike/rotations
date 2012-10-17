\name{median.SO3}
\alias{median.SO3}
\title{Compute the projected or intrinsic median estimate of the central direction}
\usage{
  median.SO3(Rs, type = "projected", epsilon = 1e-05,
    maxIter = 2000, na.rm = FALSE)
}
\arguments{
  \item{Rs}{A \eqn{n\times 9} matrix where each row
  corresponds to a random rotation in matrix form}

  \item{type}{choice of geometry: intrinsic or projected}

  \item{epsilon}{the stopping rule for the iterative
  process}

  \item{maxIter}{maximum number of iterations to allow}

  \item{na.rm}{how to handle NAs}
}
\value{
  the median
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
\references{
  Hartley R, Aftab K and Trumpf J (2011). "L1 rotation
  averaging using the Weiszfeld algorithm." In _2011 IEEE
  Conference on Computer Vision and Pattern Recognition
  (CVPR)_, pp. 3041-3048. IEEE.
}
\seealso{
  \code{\link{mean.SO3}}
}

