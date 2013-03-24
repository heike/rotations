\name{zhangCR}
\alias{zhangCR}
\alias{zhangCR.Q4}
\alias{zhangCR.SO3}
\title{Zhang CR Method}
\usage{
  zhangCR(Qs, a, m, pivot, estimator)

  \method{zhangCR}{SO3} (Rs, a, m = 300, pivot = T,
    estimator = "mean")

  \method{zhangCR}{Q4} (Qs, alpha, m = 300, pivot = T,
    estimator = "mean")
}
\arguments{
  \item{Qs}{A n-by-p matrix where each row corresponds to a
  random rotation in matrix or quaternion form}

  \item{a}{The alhpa level desired}

  \item{m}{Number of replicates to use to estiamte cut
  point}

  \item{pivot}{should the pivotal (T) or non-pivotal (F)
  method be used}

  \item{estimator}{Mean or median}
}
\value{
  radius of the confidence region centered at the projected
  mean

  \code{NULL}

  \code{NULL}
}
\description{
  Find the radius of a 100(1-a)% confidence region for the
  projected mean
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='zhang',alpha=0.9)
}

