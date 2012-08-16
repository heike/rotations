\name{SO3.mean}
\alias{SO3.mean}
\title{Compute the projected or intrinsic mean estimate of the central direction}
\usage{
  SO3.mean(Rs, type = "projected", epsilon = 1e-05,
    maxIter = 2000)
}
\arguments{
  \item{Rs}{A sample of n \eqn{3\times 3} random rotations}

  \item{type}{String indicating 'projeted' or 'intrinsic'
  type mean estimator}

  \item{epsilon}{Stopping rule for the intrinsic method}

  \item{maxIter}{The maximum number of iterations allowed
  before returning most recent estimate}
}
\value{
  projected or intrinsic mean of the sample
}
\description{
  This function takes a sample of \eqn{3\times 3} rotations
  (in the form of a \eqn{n\times 9} matrix where n is the
  sample size) and returns the projected arithmetic mean
  denoted \eqn{\widehat{\bm S}_P} or intrinsic mean
  \eqn{\widehat{\bm S}_G} according to the \code{type}
  option. For a sample of \eqn{n} random rotations
  \eqn{\bm{R}_i\in SO(3)$, $i=1,2,\dots,n}, the mean-type
  estimator is defined as
  \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})} where
  \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i} and
  the distance metric \edq{d_D} is the Riemannian or
  Euclidean.  For more on the projected mean see
  \cite{Moakher (2002)} and for the intrinsic mean see
  \cite{Manton (2004)}.
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
mean(Rs)
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
  \code{\link{SO3.median}}
}

