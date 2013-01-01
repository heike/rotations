\name{weighted.mean.SO3}
\alias{weighted.mean.SO3}
\title{Weighted Mean Rotation}
\usage{
  \method{weighted.mean}{SO3} (Rs, w, type = "projected",
    epsilon = 1e-05, maxIter = 2000, ...)
}
\arguments{
  \item{Rs}{A n-by-9 matrix where each row corresponds to a
  random rotation in matrix form}

  \item{w}{a numerical vector of weights the same length as
  Rs giving the weights to use for elements of Rs}

  \item{type}{String indicating 'projeted' or 'intrinsic'
  type mean estimator}

  \item{epsilon}{Stopping rule for the intrinsic method}

  \item{maxIter}{The maximum number of iterations allowed
  before returning most recent estimate}

  \item{...}{only used for consistency with mean.default}
}
\value{
  weighted projected mean of the sample
}
\description{
  Compute the weighted intrinsic or projected mean of a
  sample of rotations
}
\details{
  This function takes a sample of \eqn{3\times 3} rotations
  (in the form of a n-by-9 matrix where n>1 is the sample
  size) and returns the weighted projected arithmetic mean
  denoted \eqn{\widehat{\bm S}_P} or intrinsic mean
  \eqn{\widehat{\bm S}_G} according to the \code{type}
  option. For a sample of \eqn{n} random rotations
  \eqn{\bm{R}_i\in SO(3)$, $i=1,2,\dots,n}, the mean-type
  estimator is defined as
  \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})} where
  \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i} and
  the distance metric \eqn{d_D} is the Riemannian or
  Euclidean.  For more on the projected mean see
  \cite{Moakher (2002)} and for the intrinsic mean see
  \cite{Manton (2004)}.
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
wt<-abs(1/r)
weighted.mean(Rs,wt)
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.

  Manton J (2004). "A globally convergent numerical
  algorithm for computing the centre of mass on compact Lie
  groups." In _8th Conference on Control, Automation,
  Robotics and Vision, (ICARCV) _, volume 3, pp. 2211-2216.
  IEEE.
}
\seealso{
  \code{\link{median.SO3}} \code{\link{mean.SO3}}
}
