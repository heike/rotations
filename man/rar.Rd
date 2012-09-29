\name{rar}
\alias{rar}
\title{Sample of size n from target density f}
\usage{
  rar(n, f, M, ...)
}
\arguments{
  \item{n}{number of sample wanted}

  \item{f}{target density}

  \item{M}{maximum number in uniform proposal density}

  \item{...}{additional arguments sent to arsample}
}
\value{
  a vector of size n of observations from target density
}
\description{
  Sample of size n from target density f
}
\author{
  Heike Hofmann
}

