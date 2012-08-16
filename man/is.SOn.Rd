\name{is.SOn}
\alias{is.SOn}
\title{A function to determine if a given matrix is in \eqn{SO(n)} or not.  Should be deleted}
\usage{
  is.SOn(x)
}
\arguments{
  \item{x}{numeric \eqn{n \times n} matrix or vector of
  length \eqn{n^2}}
}
\value{
  logical T if the matrix is in SO(n) and false otherwise
}
\description{
  A function to determine if a given matrix is in
  \eqn{SO(n)} or not.  Should be deleted
}
\examples{
is.SOn(diag(1,3,3))
is.SOn(1:9)
}

