\name{rar}
\alias{rar}
\title{Sample of size n from target density f}
\usage{
  rar(n, f, g, M, ...)
}
\arguments{
  \item{n}{number of sample wanted}

  \item{f}{target density}

  \item{g}{sampling distribution}

  \item{M}{maximum number in uniform proposal density}

  \item{...}{additional arguments sent to arsample}
}
\value{
  a vector of size n of observations from target density
}
\description{
  Sample of size n from target density f
}
\examples{
# sample from haar distribution
x <- rar(10000, haar, runif, 1/pi, min=-pi, max=pi)

kappa=0.5
M <- max(fisher(seq(-pi, pi, length=1000), kappa))
x.fisher <- rar(10000, fisher, runif, M, min=-pi, max=pi, kappa=kappa)
}
\author{
  Heike Hofmann
}

