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
  \item{Rs,Qs}{A \eqn{n\times p}{n-by-p} matrix where each
  row corresponds to a random rotation in matrix (p=9) or
  quaternion form (p=4)}

  \item{a}{The alpha level desired, e.g. 0.95 or 0.90}

  \item{m}{Number of replicates to use to estiamte cut
  point}

  \item{pivot}{should the pivotal (T) or non-pivotal (F)
  method be used}

  \item{estimator}{Mean or median}
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
region(Rs,method='zhang',alpha=0.9)
}

