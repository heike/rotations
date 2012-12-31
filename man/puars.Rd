\name{puars}
\alias{puars}
\title{UARS distribution function}
\usage{
  puars(os, S = diag(3), kappa = 1, pangle = NULL, ...)
}
\arguments{
  \item{os}{Value at which to evaluate the UARS density}

  \item{S}{principal direction of the distribution}

  \item{kappa}{concentration of the distribution}

  \item{dangle}{The function to evaulate the angles from:
  e.g. dcayley, dvmises, dfisher, dhaar.  If left at NULL,
  the empirical CDF is used}

  \item{...}{additional arguments passed to the angular
  distribution}
}
\value{
  cdf evaulated at each os value
}
\description{
  Evaluate the UARS distributions fuction with a given
  angular distribution
}

