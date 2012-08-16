\name{is.SO3}
\alias{is.SO3}
\title{A function to determine if a given matrix is in \eqn{SO(3)} or not.}
\usage{
  is.SO3(x)
}
\arguments{
  \item{x}{numeric \eqn{n \times n} matrix or vector of
  length \eqn{n^2}}
}
\value{
  logical T if the matrix is in SO(3) and false otherwise
}
\description{
  A function to determine if a given matrix is in
  \eqn{SO(3)} or not.
}
\examples{
is.SO3(diag(1,3,3))
is.SO3(1:9)
}

