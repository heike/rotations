\name{rancourtCR.Q4}
\alias{rancourtCR.Q4}
\title{Rancourt CR Method}
\usage{
  \method{rancourtCR}{Q4} (Qs, a)
}
\arguments{
  \item{Qs}{A n-by-4 matrix where each row corresponds to a
  random rotation in matrix form}

  \item{method}{Character string specifying which type of
  interval is required}

  \item{a}{The alhpa level desired}
}
\value{
  radius of the confidence region centered at the projected
  mean
}
\description{
  Find the radius of a 100(1-a)% confidence region for the
  projected mean \cite{Rancourt et al. (2000)}
}
\details{
  This works in the same way as done in \cite{Bingham et
  al. (2009)} which assumes rotational symmetry and is
  therefore conservative.
}
\examples{
Qs<-ruars(20,rcayley,kappa=100,space='Q4')
region(Qs,method='rancourt',alpha=0.9)
}
\references{
  Rancourt D, Rivest L and Asselin J (2000). "Using
  orientation statistics to investigate variations in human
  kinematics." _Journal of the Royal Statistical Society:
  Series C (Applied Statistics)_, *49*(1), pp. 81-94.

  Bingham M, Nordman D and Vardeman S (2009). "Modeling and
  Inference for Measured Crystal Orientations and a
  Tractable Class of Symmetric Distributions for Rotations
  in three Dimensions." _Journal of the American
  Statistical Association_, *104*(488), pp. 1385-1397.
}

