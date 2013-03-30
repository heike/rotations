\name{project.SO3}
\alias{project.SO3}
\title{Projection Procedure}
\usage{
  project.SO3(M)
}
\arguments{
  \item{M}{\eqn{3\times 3}{3-by-3} matrix to project}
}
\value{
  projection of \eqn{\bm M}{M} into \eqn{SO(3)}
}
\description{
  Project an arbitrary \eqn{3\times 3}{3-by-3} matrix into
  SO(3)
}
\details{
  This function uses the process given in \cite{Moakher
  (2002)} to project an arbitrary \eqn{3\times 3}{3-by-3}
  matrix into \eqn{SO(3)}.
}
\examples{
M<-matrix(rnorm(9),3,3)
project.SO3(M)
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}
\seealso{
  \code{\link{mean.SO3}}, \code{\link{median.SO3}}
}

