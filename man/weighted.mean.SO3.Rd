\name{weighted.mean.SO3}
\alias{weighted.mean.SO3}
\title{Weighted Mean Rotation}
\usage{
  \method{weighted.mean}{SO3} (Rs, w, ...)
}
\arguments{
  \item{Rs}{A n-by-9 matrix where each row corresponds to a
  random rotation in matrix form}

  \item{w}{a numerical vector of weights the same length as
  Rs giving the weights to use for elements of Rs}

  \item{...}{only used for consistency with mean.default}
}
\value{
  weighted projected mean of the sample
}
\description{
  Compute the weighted projected mean of a sample of
  rotations
}
\details{
  This function takes a sample of \eqn{3\times 3} rotations
  (in the form of a n-by-9 matrix where n>1 is the sample
  size) and returns the weighted projected arithmetic mean
  denoted \eqn{\widehat{\bm S}_P}. For a sample of \eqn{n}
  random rotations \eqn{\bm{R}_i\in SO(3)$,
  $i=1,2,\dots,n}, the mean-type estimator is defined as
  \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})} where
  \eqn{\bar{\bm{R}}=\sum_{i=1}^nw_i\bm{R}_i} and the
  distance metric \eqn{d_D} is the Euclidean.  For more on
  the projected mean see \cite{Moakher (2002)}.
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r)
wt<-1:20
weighted.mean(Rs,wt)
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}
\seealso{
  \code{\link{median.SO3}} \code{\link{mean.SO3}}
}

