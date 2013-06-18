\name{mean.Q4}
\alias{mean.Q4}
\title{Mean Rotation}
\usage{
  \method{mean}{Q4} (Qs, type = "projected",
    epsilon = 1e-05, maxIter = 2000)
}
\arguments{
  \item{Qs}{A \eqn{n\times 4}{n-by-4} matrix where each row
  corresponds to a random rotation in unit quaternion}

  \item{type}{String indicating 'projeted' or 'geometric'
  type mean estimator}

  \item{epsilon}{Stopping rule for the geometric method}

  \item{maxIter}{The maximum number of iterations allowed
  before returning most recent estimate}
}
\value{
  projected or geometric mean of the sample
}
\description{
  Compute the projected or geometric mean of a sample of
  rotations
}
\details{
  This function takes a sample of \eqn{n} unit quaternions
  and approximates the mean rotation.  If the projected
  mean is called for then the according to \cite{Tyler
  (1981)} an estimate of the mean is the eigenvector
  corresponding to the largest eigen value of
  \eqn{\frac{1}{n}\sum_{i=1}^nq_i^\top q_i}{Q`Q/n}.  If the
  geometric mean is called then the quaternions are
  transformed into \eqn{3\times 3}{3-by-3} matrices and the
  \code{mean.SO3} function is called.
}
\examples{
Qs<-ruars(20,rcayley,space="Q4")
mean(Qs,type='geometric')
}
\references{
  Manton J (2004). "A globally convergent numerical
  algorithm for computing the centre of mass on compact Lie
  groups." In _8th Conference on Control, Automation,
  Robotics and Vision, (ICARCV) _, volume 3, pp. 2211-2216.
  IEEE.

  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.

  Tyler DE (1981). "Asymptotic inference for eigenvectors."
  _The Annals of Statistics_, pp. 725-736.
}
\seealso{
  \code{\link{mean.SO3}}
}

