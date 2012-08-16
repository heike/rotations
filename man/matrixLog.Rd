\name{matrixLog}
\alias{matrixLog}
\title{This fuction will compute the natural logarithm of a matrix in SO(n).  It uses the special case of the Taylor expansion for SO(n) matrices.}
\usage{
  matrixLog(R)
}
\arguments{
  \item{R}{numeric matrix in \eqn{SO(n)}}
}
\value{
  mlog numeric matrix \eqn{\log(R)}
}
\description{
  For details see \cite{Moakher (2002)}
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}

