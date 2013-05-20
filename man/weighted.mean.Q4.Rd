\name{weighted.mean.Q4}
\alias{weighted.mean.Q4}
\title{Weighted Rotation Median}
\usage{
  \method{weighted.mean}{Q4} (Qs, w, type = "projected",
    epsilon = 1e-05, maxIter = 2000)
}
\arguments{
  \item{Qs}{A \eqn{n\times 4} matrix where each row
  corresponds to a random rotation in unit quaternion}

  \item{w}{a numerical vector of weights the same length as
  Rs giving the weights to use for elements of Rs}

  \item{type}{String indicating 'projeted' or 'intrinsic'
  type mean estimator}

  \item{epsilon}{Stopping rule for the intrinsic method}

  \item{maxIter}{The maximum number of iterations allowed
  before returning most recent estimate}
}
\value{
  weighted projected or intrinsic mean of the sample
}
\description{
  Compute the weighted projected or intrinsic mean of a
  sample of rotations
}
\details{
  This function takes a sample of n unit quaternions and
  approximates the mean rotation.  If the projected mean is
  called for then the quaternions are turned
  reparameterized to matrices and mean.SO3 is called.  If
  the intrinsic mean is called then according to
  \cite{Gramkow (2001)} a better approximation is achieved
  by taking average quaternion and normalizing.  Our
  simulations don't match this claim.
}
\examples{
r<-rvmises(20,0.01)
wt<-abs(1/r)
Qs<-genR(r,space="Q4")
weighted.mean(Qs,wt)
}
\references{
  Gramkow C (2001). "On averaging rotations." _Journal of
  Mathematical Imaging and Vision_, *15*(1), pp. 7-16.

  Manton J (2004). "A globally convergent numerical
  algorithm for computing the centre of mass on compact Lie
  groups." In _8th Conference on Control, Automation,
  Robotics and Vision, (ICARCV) _, volume 3, pp. 2211-2216.
  IEEE.

  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}
\seealso{
  \code{\link{mean.SO3}}
}

