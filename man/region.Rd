\name{region}
\alias{region}
\alias{region.Q4}
\alias{region.SO3}
\title{Confidence Region for Mean Rotation}
\usage{
  region(Qs, method, alpha, ...)

  \method{region}{Q4} (Qs, method, alpha, ...)

  \method{region}{SO3} (Rs, method, alpha, ...)
}
\arguments{
  \item{Rs,Qs}{A \eqn{n\times p}{n-by-p} matrix where each
  row corresponds to a random rotation in matrix (p=9) or
  quaternion form (p=4)}

  \item{method}{Character string specifying which type of
  interval is required}

  \item{alpha}{The alpha level desired, e.g. 0.95 or 0.90}

  \item{...}{Additional arguments}
}
\value{
  radius of the confidence region centered at the projected
  mean
}
\description{
  Find the radius of a \eqn{100\alpha%} confidence region
  for the projected mean
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Qs,method='rancourt',alpha=0.9)
}
\references{
  Rancourt D, Rivest L and Asselin J (2000). "Using
  orientation statistics to investigate variations in human
  kinematics." _Journal of the Royal Statistical Society:
  Series C (Applied Statistics)_, *49*(1), pp. 81-94.
}

