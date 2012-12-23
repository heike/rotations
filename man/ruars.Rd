\name{ruars}
\alias{ruars}
\title{UARS random deviates}
\usage{
  ruars(n, S = diag(3), kappa = 1, rangle, space = "SO3",
    ...)
}
\arguments{
  \item{n}{number of observations. If \code{length(n)>1},
  the length is taken to be n}

  \item{S}{principal direction of the distribution}

  \item{kappa}{concentration of the distribution}

  \item{rangle}{The function from which to simulate angles:
  e.g. rcayley, rvmises, rhaar, rfisher}

  \item{space}{Indicates the desired representation: matrix
  (SO3), quaternion (Q4) or Euler angles (EA)}

  \item{...}{additional arguments passed to the angular
  function}
}
\value{
  random deviates from the specified UARS distribution
}
\description{
  Produce random deviates from a chosen UARS distribution.
}

