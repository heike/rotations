\name{rar}
\alias{rar}
\title{Random simulation from angle distributions}
\usage{
  rar(n, f, g, M, ...)
}
\arguments{
  \item{n}{sample size}

  \item{f}{density to be sampled from}

  \item{g}{density, enveloping f, i.e. f < g}

  \item{M}{real valued constant}

  \item{...}{parameters for densities f and g}
}
\description{
  acceptance-rejection random sampling from angle
  distributions
}
\examples{
# sample from haar distribution
x <- rar(10000, haar, runif, 1/pi, min=-pi, max=pi)

kappa=0.5
M <- max(fisher(seq(-pi, pi, length=1000), kappa))
x.fisher <- rar(10000, fisher, runif, M, min=-pi, max=pi, kappa=kappa)
}

