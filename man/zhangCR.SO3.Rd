\name{zhangCR.SO3}
\alias{zhangCR.SO3}
\title{Zhang CR Method}
\usage{
  \method{zhangCR}{SO3} (Rs, alpha, m = 300, pivot = T)
}
\arguments{
  \item{Rs}{A n-by-9 matrix where each row corresponds to a
  random rotation in matrix form}

  \item{m}{Number of replicates to use to estiamte cut
  point}

  \item{alpha}{The alhpa level desired}

  \item{pivot}{should the pivotal (T) or non-pivotal (F)
  method be used}
}
\value{
  radius of the confidence region centered at the projected
  mean
}
\description{
  Find the radius of a 100(1-a)% confidence region for the
  projected mean
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='zhang',alpha=0.9)
}

