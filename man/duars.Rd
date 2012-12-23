\name{duars}
\alias{duars}
\title{UARS density function}
\usage{
  duars(o, S = diag(3), kappa = 1, dangle, ...)
}
\arguments{
  \item{o}{Value at which to evaluate the UARS density}

  \item{S}{principal direction of the distribution}

  \item{kappa}{concentration of the distribution}

  \item{dangle}{The function to evaulate the angles from:
  e.g. dcayley, dvmises, dfisher, dhaar}

  \item{...}{additional arguments passed to the angular
  distribution}
}
\value{
  density value at o
}
\description{
  Evaluate the UARS density with a given angular
  distribution.
}

