\name{matrixExp}
\alias{matrixExp}
\title{This fuction will compute the natural exponential of skew-symmetric matrix.}
\usage{
  matrixExp(A)
}
\arguments{
  \item{A}{3-dimensional skew-symmetric matrix, i.e.,
  \eqn{\bm A=-\bm A^\top}}
}
\value{
  numeric matrix \eqn{e^{\bm A}}
}
\description{
  See \cite{Moakher (2002)}
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}

