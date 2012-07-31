\name{matrixExp}
\alias{matrixExp}
\title{This fuction will compute the natural exponential of skew-symmetric matrix.  It uses the special case of the Taylor expansion for SO(n) matrices.}
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
  This fuction will compute the natural exponential of
  skew-symmetric matrix.  It uses the special case of the
  Taylor expansion for SO(n) matrices.
}

