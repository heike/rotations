\name{duars}
\alias{duars}
\title{UARS density function}
\usage{
  duars(o, S = diag(3), kappa, c, ...)
}
\arguments{
  \item{o}{Value at which to evaluate the UARS density}

  \item{S}{principal direction of the distribution}

  \item{kappa}{concentration of the distribution}

  \item{c}{Circular distribution, choices: 'dcayley',
  'dfisher', 'dhaar' and 'dvmises'}
}
\value{
  density value at o
}
\description{
  Evaluate the UARS density with a given angular
  distribution.
}

