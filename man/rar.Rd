\name{rar}
\alias{rar}
\title{Call arsample 'n' times to get a sample of size n from target density f}
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
  Call arsample 'n' times to get a sample of size n from
  target density f
}
\author{
  Heike Hofmann
}

