\name{fisherCR}
\alias{fisherCR}
\alias{fisherCR.Q4}
\alias{fisherCR.SO3}
\title{Fisher Mean Polax Axis CR Method}
\usage{
  fisherCR(Qs, a, boot, m)

  \method{fisherCR}{Q4} (Qs, a, boot = T, m = 300)

  \method{fisherCR}{SO3} (Rs, a, boot = T, m = 300)
}
\arguments{
  \item{Qs}{A n-by-4 matrix where each row corresponds to a
  random rotation in matrix form}

  \item{a}{The alpha level desired}

  \item{boot}{Should the bootstrap or normal theory
  critical value be used}

  \item{m}{number of bootstrap replicates to use to
  estimate critical value}
}
\value{
  radius of the confidence region centered at the projected
  mean

  \code{NULL}

  \code{NULL}
}
\description{
  Find the radius of a 100(1-a)% confidence region for the
  projected mean \cite{Fisher et al. (1996)}
}
\details{
  This works in the same way as done in \cite{Bingham et
  al. (2009)} which assumes rotational symmetry and is
  therefore conservative.
}
\examples{
Qs<-ruars(20,rcayley,kappa=100,space='Q4')
region(Qs,method='fisher',alpha=0.9)
}
\references{
  Fisher N, Hall P, Jing B and Wood A (1996). "Improved
  pivotal methods for constructing confidence regions with
  directional data." _Journal of the American Statistical
  Association_, *91*(435), pp. 1062-1070.

  Bingham M, Nordman D and Vardeman S (2009). "Modeling and
  Inference for Measured Crystal Orientations and a
  Tractable Class of Symmetric Distributions for Rotations
  in three Dimensions." _Journal of the American
  Statistical Association_, *104*(488), pp. 1385-1397.
}

