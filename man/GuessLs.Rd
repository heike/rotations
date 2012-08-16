\name{GuessLs}
\alias{GuessLs}
\title{\eqn{d_R^p} based estimators of the Central Direction}
\usage{
  GuessLs(Rs, maxe = 0.001, p)
}
\arguments{
  \item{Rs}{The sample of random rotations in the form of
  an nx9 matrix}

  \item{maxe}{The stopping rule of the optimization
  process}

  \item{p}{The order of the riemannian distance to
  minimize}
}
\value{
  S3 \item{Shatp}{The estimated random rotation
  \eqn{\widehat{\bm S}_P}}
}
\description{
  THIS NEEDS TO BE REMOVED.  This functions is slow but it
  computes the element of SO(3) that minimizes the sum of
  the pth order Riemannian distances. It returns the the
  random roatation Shat.  It calls the function SumDistR
  which calculates the sum of the pth order Riemannian
  distances between the sample Rs and S.  This needs to be
  trashed, most likely.
}

