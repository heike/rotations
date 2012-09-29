\name{CIradius}
\alias{CIradius}
\alias{CIradius.EA}
\alias{CIradius.Q4}
\alias{CIradius.SO3}
\title{Find confidence interval radius for central orientation estimate}
\usage{
  CIradius(Rs, main = mean, B = 1000, m = nrow(Rs),
    alpha = 0.95, ...)

  \method{CIradius}{SO3} (Rs, main = mean, B = 1000,
    m = nrow(Rs), alpha = 0.95, ...)

  \method{CIradius}{Q4} (Qs, main = mean, B = 1000,
    m = nrow(Qs), alpha = 0.95, ...)

  \method{CIradius}{EA} (EAs, main = mean, B = 1000,
    m = nrow(EAs), alpha = 0.95, ...)
}
\arguments{
  \item{Rs}{sample of size n rotations in matrix format}

  \item{main}{The method of central orientation estimation;
  e.g. mean or median, other methods can be used, but will
  result in a warning}

  \item{B}{bootstrap iterations}

  \item{m}{bootstrap sample size}

  \item{alpha}{level of confidence}

  \item{...}{additional arguments passed to 'fun' such as
  'type' and 'epsilon'}
}
\value{
  the radius of the confidence cone

  \code{NULL}

  \code{NULL}

  \code{NULL}
}
\description{
  Given a sample or rotations and estimator for the central
  orientation, a confidence interval radius will be
  estimated based on a bootstrap procedure.  The level of
  confidence, bootstrap iterations and bootstrap sample
  sizes are all options for the user to set.  The procedure
  is as follows: 0) Estimate Shat from (R_1,...,R_n) in Rs
  1) Sample R_1,...,R_m from Rs with replacement 2)
  Estimate ShatStar from bootstrap sample 3) Compute
  That=d_r(Shat,ShatStar) 4) Repeat steps 1-3 B times 5)
  Report q% percentile of That to be the radius of the
  confidence 'cone'
}
\examples{
Rs<-genR(rvmises(20))
 CIradius(Rs,mean,type='projected')  ## calls CIradius.SO3
 Qs<-genR(rvmises(20),space='Q4')
 CIradius(Qs,mean,type='projected')  ## calls CIradius.Q4
}
\references{
  bingham10
}

