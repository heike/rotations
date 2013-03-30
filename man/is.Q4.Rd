\name{is.Q4}
\alias{is.Q4}
\title{A function to determine if a given matrix is in unit quaternion or not.}
\usage{
  is.Q4(x)
}
\arguments{
  \item{x}{numeric 1-by-4 vector of length 9}
}
\value{
  logical T if the vector is a unit quaternion and false
  otherwise
}
\description{
  A function to determine if a given matrix is in unit
  quaternion or not.
}
\seealso{
  Other SO3: \code{\link{as.SO3}}, \code{\link{id.SO3}},
  \code{\link{is.SO3}}, \code{\link{SO3}},
  \code{\link{SO3-class}}, \code{\link{SO3.default}},
  \code{\link{SO3.Q4}}
}

