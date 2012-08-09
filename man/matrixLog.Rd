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
  For details see \cite{moakher02}
}

