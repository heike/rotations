\name{mean.Q4}
\alias{mean.Q4}
\title{Compute the projected or intrinsic mean estimate of the central direction}
\usage{
  mean.Q4(Qs, type = "projected", epsilon = 1e-05,
    maxIter = 2000)
}
\arguments{
  \item{Qs}{A \eqn{n\times 4} matrix where each row
  corresponds to a random rotation in unit quaternion}

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
  This function takes a sample of n unit quaternions and
  approximates the mean rotation.  If the projected mean is
  called for then the quaternions are turned
  reparameterized to matrices and mean.SO3 is called.  If
  the intrinsic mean is called then according to
  \cite{gramkow01} a better approximation is achieved by
  taking average quaternion and normalizing.  Our
  simulations don't match this claim.
}
\examples{
r<-rvmises(20,0.01)
Qs<-genR(r,space="Q4")
mean(Qs,type='intrinsic')
}
\seealso{
  \code{\link{mean.SO3}}
}

