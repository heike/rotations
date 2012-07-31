\name{projMatrix}
\alias{projMatrix}
\title{The projection of an arbitrary \eqn{3\times 3} matrix into \eqn{SO(3)}}
\usage{
  projMatrix(M)
}
\arguments{
  \item{M}{\eqn{3\times 3} matrix to project}
}
\value{
  projection of \eqn{\bm M} into \eqn{SO(3)}
}
\description{
  This function uses the process given in Moakher 2002 to
  project an arbitrary \eqn{3\times 3} matrix into
  \eqn{SO(3)}.
}
\examples{
M<-matrix(rnorm(9),3,3)
projMatrix(M)
}
\seealso{
  \code{\link{arith.mean}}, \code{\link{rmedian}}
}

