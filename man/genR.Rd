\name{genR}
\alias{genR}
\title{Generate rotation matrix given misorientation angle, r}
\usage{
  genR(r, S = diag(1, 3, 3))
}
\arguments{
  \item{r}{The angle through which all three dimensions are
  rotated after the axis was picked uniformly on the unit
  sphere}

  \item{S}{the principle direction}
}
\value{
  a \eqn{n\times 9} matrix in SO(3) with misorientation
  angle r and principal direction S
}
\description{
  A function that generates a random rotation in
  \eqn{SO(3)} following a Uniform-Axis random roation
  distribution with central direction S The exact form of
  the UARS distribution depends upon the distribution of
  the roation r
}
\examples{
r<-rvmises(20,0.01)
genR(r)
}

