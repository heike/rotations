\name{dist}
\alias{dist}
\alias{dist.Q4}
\alias{dist.SO3}
\title{Rotational Distance}
\usage{
  dist(x, ...)

  \method{dist}{SO3} (R1, R2 = id.SO3,
    method = "projected", p = 1)

  \method{dist}{Q4} (Q1, Q2 = id.Q4, method = "projected",
    p = 1)
}
\arguments{
  \item{R1}{(or Q1) a rotation in matrix or quaternion
  representation}

  \item{R2}{(or Q2) the second rotation in the same
  parameterization as R1}

  \item{method}{String indicating 'projected' or
  'intrinsic' method of distance}

  \item{p}{the order of the distance}

  \item{...}{Additional arguments}
}
\value{
  the rotational distance between R1 and R2
}
\description{
  Calculate the Euclidean or Riemannian distance between
  two rotations
}
\details{
  This function will calculate the intrinsic (Riemannian)
  or projected (Euclidean) distance between two rotations.
  If only one rotation is specified the other will be set
  to the identity and the distance between the two is
  returned.  For rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
  both in \eqn{SO(3)}, the Euclidean distance between them
  is \deqn{||R_1-R_2||_F}{||R1-R2||} where
  \eqn{||\cdot||_F}{|| ||} is the Frobenius norm. The
  intrinsic distance is defined as \deqn{||Log(R_1^\top
  R_2)||_F}{||Log(R1'R2)||} where \eqn{Log} is the matrix
  logarithm, and it corresponds to the misorientation angle
  of \eqn{R_1^\top R_2}{R1'R2}.
}

