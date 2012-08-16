\name{arith_mean}
\alias{arith_mean}
\title{Projected Arithmetic Mean \eqn{\hat{S}_P}}
\usage{
  arith_mean(Rs)
}
\arguments{
  \item{Rs}{A sample of n \eqn{3\times 3} random rotations}
}
\value{
  A \eqn{3\times 3} matrix in SO(3) called the Projected
  arithmetic mean
}
\description{
  THIS HAS BEEN REPLACED WITH 'mean.SO3'!!  This function
  takes a sample of \eqn{3\times 3} rotations (in the form
  of a \eqn{n\times 9} matrix where n is the sample size)
  and returns the projected arithmetic mean denoted
  \eqn{\widehat{\bm S}_P}. For a sample of $n$ random
  rotations \eqn{\bm{R}_i\in SO(3)$, $i=1,2,\dots,n}, this
  mean-type estimator is defined as \deqn{\[
  \widehat{\bm{S}}_P=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_E^2(\bm{R}_i,\bm{S})=\argmax_{\bm{S}\in
  SO(3)}\tr(\bm{S}^{\top}\bar{\bm{R}}) \]} where
  \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i}. First
  the mean of each element is calculated then that matrix
  is projected to SO(3) in accordance with the procedure
  presented in Moahker's 2003 paper
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
arith_mean(Rs)
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}
\seealso{
  \code{\link{MantonL2}}, \code{\link{HartleyL1}},
  \code{\link{rmedian}}
}

