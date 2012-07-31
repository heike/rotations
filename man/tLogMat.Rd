\name{tLogMat}
\alias{tLogMat}
\title{Log of a matrix times some center S}
\usage{
  tLogMat(x, S)
}
\arguments{
  \item{x}{vector of length 9}

  \item{S}{\eqn{3\times 3} matrix}
}
\value{
  skew-symmetic matrix \eqn{\log(\bm{S}^\top\bm{x})}
}
\description{
  Used to speed up the Riemannian based estimators
}
\examples{
rs<-rvmises(20,1)
Rs<-genR(rs)
apply(Rs,1,tLogMat,S=diag(1,3,3))
}
\seealso{
  \code{\link{MantonL2}}, \code{\link{HartleyL1}}
}

