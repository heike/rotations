\name{median}
\alias{median}
\alias{median.EA}
\alias{median.Q4}
\alias{median.SO3}
\title{Compute the projected or intrinsic median estimate of the central direction}
\usage{
  \method{median}{SO3} (x, ...)

  \method{median}{SO3} (Rs, type = "projected",
    epsilon = 1e-05, maxIter = 2000, na.rm = FALSE)

  \method{median}{Q4} (x, type = "projected",
    epsilon = 1e-05, maxIter = 2000, na.rm = FALSE)

  \method{median}{EA} (x, type = "projected",
    epsilon = 1e-05, maxIter = 2000, na.rm = FALSE)
}
\arguments{
  \item{x}{A \eqn{n\times 9} matrix where each row
  corresponds to a random rotation in matrix form}

  \item{...}{additional arguments}
}
\value{
  the median

  \code{NULL}

  \code{NULL}

  \code{NULL}
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

