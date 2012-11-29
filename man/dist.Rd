\name{dist}
\alias{dist}
\alias{dist.EA}
\alias{dist.Q4}
\alias{dist.SO3}
\title{Rotational Distance}
\usage{
  dist(x, ...)

  \method{dist}{SO3} (R1, R2 = id.SO3,
    method = "projected", p = 1)

  \method{dist}{Q4} (R1, R2 = id.Q4, method = "projected",
    p = 1)

  \method{dist}{EA} (R1, R2 = id.EA, method = "projected",
    p = 1)
}
\arguments{
  \item{x}{rotation in SO3 representation}

  \item{...}{Additional arguments}
}
\value{
  the distance between x and something else

  \code{NULL}

  \code{NULL}

  \code{NULL}
}
\description{
  This function will calculate the intrinsic (Riemannian)
  or projected (Euclidean) distance between two rotations.
  If only one rotation is specified the other will be set
  to the identity and the distance between the two is
  returned.
}

