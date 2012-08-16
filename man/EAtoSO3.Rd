\name{EAtoSO3}
\alias{EAtoSO3}
\title{A function that will take in a Euler angle and return a rotation matrix in vector format}
\usage{
  EAtoSO3(eur)
}
\arguments{
  \item{eur}{numeric Euler angle representation of an
  element in SO(3)}
}
\value{
  numeric \eqn{9\times 1} vector of a matrix in SO(3)
}
\description{
  A function that will take in a Euler angle and return a
  rotation matrix in vector format
}
\examples{
eaExample<-c(pi/2,3*pi/4,0)
SO3Dat<-EAtoSO3(eaExample)
is.SO3(SO3Dat)
}
\seealso{
  \code{\link{is.SO3}} can be used to check the output of
  this function
}

