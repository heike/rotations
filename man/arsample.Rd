\name{arsample}
\alias{arsample}
\title{Accept/reject algorithm written by Dr. Hofmann}
\usage{
  arsample(f, g, M, kappa, ...)
}
\arguments{
  \item{f}{target density}

  \item{g}{sampling density}

  \item{M}{maximum in uniform density}

  \item{kappa}{second parameter in the target density}

  \item{...}{additional arguments passed to samping
  density, g}
}
\value{
  a random observation from target density
}
\description{
  Accept/reject algorithm written by Dr. Hofmann
}
\author{
  Heike Hofmann
}

