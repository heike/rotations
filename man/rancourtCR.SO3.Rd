\name{rancourtCR.SO3}
\alias{rancourtCR.SO3}
\title{Rancourt CR Method}
\usage{
  \method{rancourtCR}{SO3} (Rs, a)
}
\arguments{
  \item{Rs}{A n-by-9 matrix where each row corresponds to a
  random rotation in matrix form}

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
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='rancourt',a=0.9)
}
\references{
  Rancourt D, Rivest L and Asselin J (2000). "Using
  orientation statistics to investigate variations in human
  kinematics." _Journal of the Royal Statistical Society:
  Series C (Applied Statistics)_, *49*(1), pp. 81-94.
}

