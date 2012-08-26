\name{CIradius.Q4}
\alias{CIradius.Q4}
\title{Find confidence interval radius for central orientation estimate}
\usage{
  CIradius.Q4(Qs, fun = "mean", B = 1000, m = n, q = 0.95,
    ...)
}
\arguments{
  \item{Qs}{sample of size n rotations in matrix format}

  \item{fun}{The method of central orientation estimation;
  e.g. mean or median}

  \item{B}{bootstrap iterations}

  \item{m}{bootstrap sample size}

  \item{q}{level of confidence}

  \item{...}{additional arguments passed to 'fun' such as
  'type' and 'epsilon'}
}
\value{
  the radius of the confidence cone
}
\description{
  Works the same was as 'CIradius.SO3' except with the
  quaternion representation of rotations.
}
\examples{
Qs<-genR(rvmises(20),space='Q4')
 CIradius.SO3(Qs,mean.Q4,type='projected')
}
\references{
  bingham10
}
\seealso{
  \code{\link{CIradius.SO3}}
}

